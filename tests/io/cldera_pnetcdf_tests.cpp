#include "io/cldera_pnetcdf.hpp"

#include <ekat/ekat_assert.hpp>

#include <catch2/catch.hpp>

#include <pnetcdf.h>
#include <iostream>

void write ()
{
  using namespace cldera;

  // Get comm specs
  ekat::Comm comm(MPI_COMM_WORLD);
  const int rank = comm.rank();
  const int size = comm.size();

  // Create file
  auto file = open_file ("test.nc",comm,IOMode::Write);

  // Add dims
  add_dim (*file,"lat",12);
  add_dim (*file,"lon",36);
  add_dim (*file,"dim2",2);

  REQUIRE_THROWS (add_dim(*file,"lat",10)); // Dim already added

  // Add vars
  add_var(*file,"T","double",{"lat","lon"},true);
  add_var(*file,"V","float",{"lat","lon","dim2"},false);
  add_var(*file,"I","long long",{"lat","lon"},true);

  REQUIRE_THROWS (add_var(*file,"V","float",{"lat","lon","dim2"},false)); // Var already added
  REQUIRE_THROWS (add_var(*file,"W","floatsies",{"lat","lon","dim2"},false)); // Unknown dtype
  REQUIRE_THROWS (add_var(*file,"W","float",{"lat","lon","dim3"},false)); // Unknown dim

  // Partition along lat dimension
  std::vector<int> my_lat;
  for (int i=0; i<12; ++i) {
    if (i%size == rank) {
      my_lat.push_back(i);
    }
  }
  add_decomp (*file,"lat",my_lat);

  REQUIRE_THROWS (add_decomp (*file,"lat",my_lat)); // Decomp for "lat" already added
  REQUIRE_THROWS (add_var(*file,"W","float",{"lon","lat","dim2"},false)); // Decomp dim is not first

  // Close define phase
  enddef (*file);

  REQUIRE_THROWS (add_dim(*file,"dim3",3)); // Not in def mode

  // Create data in an rank-independent way
  int nlats = my_lat.size();
  std::vector<double>     tdata(nlats*36);
  std::vector<float>      vdata(nlats*36*2);
  std::vector<long long>  idata(nlats*36);
  for (int i=0; i<nlats; ++i) {
    int idx = my_lat[i];
    std::iota(tdata.begin()+i*36,  tdata.begin()+(i+1)*36,idx*36);
    std::iota(vdata.begin()+i*36*2,vdata.begin()+(i+1)*36*2,idx*36*2);
    std::iota(idata.begin()+i*36,  idata.begin()+(i+1)*36,idx*36);
  }

  write_var(*file,"T",tdata.data());
  write_var(*file,"V",vdata.data());
  write_var(*file,"I",idata.data());

  // Check write exceptions
  REQUIRE_THROWS (write_var(*file,"W",vdata.data())); // Not a valid var
  REQUIRE_THROWS (write_var(*file,"I",vdata.data())); // Wrong data type ptr

  close_file(*file);
}

void read ()
{
  using namespace cldera;

  ekat::Comm comm(MPI_COMM_WORLD);
  const int rank = comm.rank();
  const int size = comm.size();

  auto file = open_file ("test.nc",comm,IOMode::Read);

  // Partition along lat dimension, but differently from how it was done during write
  // Divide lat entries evenly, and do some extra work to account for remainders
  int nglat = file->dims.at("lat")->len;
  int nglon = file->dims.at("lon")->len;
  std::vector<int> my_lat;
  int num_local = nglat / size;
  int rem = nglat % size;
  if (rem>0 && rank<rem) { ++num_local; }

  EKAT_REQUIRE_MSG (num_local>=1,
      "Error! This test cannot be run with a number of ranks larger than the lat size.\n");

  std::vector<int> offsets (size+1,0);
  comm.all_gather(&num_local,&offsets[1],1);
  for (int r=0; r<size; ++r) {
    offsets[r+1] += offsets[r];
  }
  REQUIRE (offsets[size]==nglat);

  for (int i=0; i<num_local; ++i) {
    my_lat.push_back(i+offsets[rank]);
  }
  add_decomp (*file,"lat",my_lat);

  // Check dims
  REQUIRE (file->dims.size()==4); // "time" is automatically added
  REQUIRE (file->dims.at("lat")->len==nglat);
  REQUIRE (file->dims.at("lon")->len==nglon);
  REQUIRE (file->dims.at("dim2")->len==2);

  // Check vars
  REQUIRE (file->vars.size()==3);

  auto T = file->vars.at("T");
  REQUIRE (T->dims.size()==3);
  REQUIRE (T->dims[0]->name=="time");
  REQUIRE (T->dims[1]->name=="lat");
  REQUIRE (T->dims[2]->name=="lon");
  REQUIRE (T->dtype=="double");

  auto V = file->vars.at("V");
  REQUIRE (V->dims.size()==3);
  REQUIRE (V->dims[0]->name=="lat");
  REQUIRE (V->dims[1]->name=="lon");
  REQUIRE (V->dims[2]->name=="dim2");
  REQUIRE (V->dtype=="float");

  auto I = file->vars.at("I");
  REQUIRE (I->dims.size()==3);
  REQUIRE (I->dims[0]->name=="time");
  REQUIRE (I->dims[1]->name=="lat");
  REQUIRE (I->dims[2]->name=="lon");
  REQUIRE (I->dtype=="long long");

  // Check data
  int nlats = my_lat.size();
  std::vector<double>     tdata(nlats*nglon);
  std::vector<float>      vdata(nlats*nglon*2);
  std::vector<long long>  idata(nlats*nglon);

  read_var(*file,"T",tdata.data());
  read_var(*file,"V",vdata.data());
  read_var(*file,"I",idata.data());

  for (int i=0; i<nlats; ++i) {
    int idx = my_lat[i];
    for (int j=0; j<nglon; ++j) {
      REQUIRE(tdata[i*nglon+j]==(idx*nglon+j));
      REQUIRE(vdata[i*nglon*2+2*j]==(idx*nglon*2+2*j));
      REQUIRE(vdata[i*nglon*2+2*j+1]==(idx*nglon*2+2*j+1));
      REQUIRE(idata[i*nglon+j]==(idx*nglon+j));
    }
  }

  // Check read exceptions
  REQUIRE_THROWS (read_var(*file,"V",vdata.data(),1)); // Not a time-dep var
  REQUIRE_THROWS (read_var(*file,"I",idata.data(),2)); // Not a valid record
  REQUIRE_THROWS (read_var(*file,"I",vdata.data()));   // Wrong data type ptr

  close_file(*file);
}

TEST_CASE ("pnetcdf_io") {
  write ();
  read  ();
}

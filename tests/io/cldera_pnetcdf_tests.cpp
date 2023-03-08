#include "io/cldera_pnetcdf.hpp"

#include <ekat/ekat_assert.hpp>

#include <catch2/catch.hpp>

#include <iostream>

void write_bad ()
{
  using namespace cldera;
  using namespace cldera::io::pnetcdf;

  // Get comm specs
  ekat::Comm comm(MPI_COMM_WORLD);
  const int rank = comm.rank();
  const int size = comm.size();

  EKAT_REQUIRE_MSG (size<=4,
      "Error! This test is meant for up to 4 MPI ranks.\n"
      "  Number MPI ranks found: " + std::to_string(size) + "\n");

  const int nglat = 12;
  const int nglon = 12;
  const int nlat = nglat / size;
  const int nlon = nglon / size;

  // Create file
  auto file = open_file ("test_bad_np" + std::to_string(size) +".nc",comm,IOMode::Write);

  // Add dims
  add_dim (*file,"lat",nlat,true);
  add_dim (*file,"lon",nlon,true);

  REQUIRE_THROWS (add_dim(*file,"lat",2*nlat)); // Dim already added with different length

  // Add vars
  add_var(*file,"T","double",{"lat","lon"},true);
  REQUIRE_THROWS (add_var(*file,"W","floatsies",{"lat","lon","dim2"},false)); // Unknown dtype
  REQUIRE_THROWS (add_var(*file,"W","float",{"lat","lon","dim3"},false)); // Unknown dim
  REQUIRE_THROWS (add_var(*file,"T","double",{"lat","lon"},false)); // Different time dep
  REQUIRE_THROWS (add_var(*file,"T","double",{"lat"},true)); // Different layout
  REQUIRE_THROWS (add_var(*file,"T","int",{"lat","lon"},true)); // Different data type

  REQUIRE_THROWS (set_att(*file,"blah","XYZ",10)); // Not a valid var name

  REQUIRE_THROWS(add_decomp (*file,"lat",{})); // Can't add decomp while in def mode

  // Close define phase
  enddef (*file);

  REQUIRE_THROWS (add_dim(*file,"dim3",3)); // Not in def mode

  // Check write exceptions
  double ddata;
  int idata;
  double* ptr;
  REQUIRE_THROWS (write_var(*file,"W",&ddata)); // Not a valid var
  REQUIRE_THROWS (write_var(*file,"T",&idata)); // Wrong data type ptr
  REQUIRE_THROWS (write_var(*file,"T",ptr));    // Invalid ptr

  // Add decomposition
  REQUIRE_THROWS(add_decomp (*file,"lat",{})); // entries don't add up

  // Partition along lat dimension
  std::vector<int> my_lat;
  for (int i=0; i<nlat; ++i) {
    my_lat.push_back(i*size+rank);
  }
  REQUIRE_THROWS(add_decomp (*file,"lat",my_lat)); // 'T' is decomposed along two dimensions

  close_file(*file);
}

void write ()
{
  using namespace cldera;
  using namespace cldera::io::pnetcdf;

  // Get comm specs
  ekat::Comm comm(MPI_COMM_WORLD);
  const int rank = comm.rank();
  const int size = comm.size();

  EKAT_REQUIRE_MSG (size<=4,
      "Error! This test is meant for up to 4 MPI ranks.\n"
      "  Number MPI ranks found: " + std::to_string(size) + "\n");

  const int nglat = 12;
  const int nglon = 12;
  const int nlat = nglat / size;

  // Create file
  auto file = open_file ("test_np" + std::to_string(size) +".nc",comm,IOMode::Write);

  // Add dims
  add_dim (*file,"lat",nlat,true);
  add_dim (*file,"lon",nglon);
  add_dim (*file,"dim2",2);

  // Add vars
  add_var(*file,"T","double",{"lat","lon"},true);
  add_var(*file,"V","float",{"lat","lon","dim2"},false);
  add_var(*file,"I","long long",{"lat","lon"},true);
  add_var(*file,"V","float",{"lat","lon","dim2"},false); // Var already added, but with same specs, so should be fine

  // Partition along lat dimension
  std::vector<int> my_lat;
  for (int i=0; i<nglat; ++i) {
    if (i%size == rank) {
      my_lat.push_back(i);
    }
  }

  set_att (*file,"my_int","NC_GLOBAL",static_cast<int>(42));
  set_att (*file,"my_dbl","NC_GLOBAL",static_cast<double>(42));
  set_att (*file,"units","V",std::string("m/s"));

  // Close define phase
  enddef (*file);

  // Add decomposition
  add_decomp (*file,"lat",my_lat);

  // Create data in an rank-independent way
  int nlats = my_lat.size();
  std::vector<double>     tdata(nlats*nglon);
  std::vector<float>      vdata(nlats*nglon*2);
  std::vector<long long>  idata(nlats*nglon);
  for (int i=0; i<nlats; ++i) {
    int idx = my_lat[i];
    std::iota(tdata.begin()+i*nglon,  tdata.begin()+(i+1)*nglon,idx*nglon);
    std::iota(vdata.begin()+i*nglon*2,vdata.begin()+(i+1)*nglon*2,idx*nglon*2);
    std::iota(idata.begin()+i*nglon,  idata.begin()+(i+1)*nglon,idx*nglon);
  }

  write_var(*file,"T",tdata.data());
  write_var(*file,"V",vdata.data());
  write_var(*file,"I",idata.data());

  close_file(*file);
}

void read ()
{
  using namespace cldera;
  using namespace cldera::io::pnetcdf;

  ekat::Comm comm(MPI_COMM_WORLD);
  const int rank = comm.rank();
  const int size = comm.size();

  auto file = open_file ("test_np" + std::to_string(size) +".nc",comm,IOMode::Read);

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
  REQUIRE (file->dims.at("lat")->glen==nglat);
  REQUIRE (file->dims.at("lat")->len==num_local);
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
      if(not(tdata[i*nglon+j]==(idx*nglon+j))) {
        std::cout << "tdata:\n";
        for (int i=0; i<nlats; ++i) {
          for (int j=0; j<nglon; ++j) {
            std::cout << " " << tdata[i*nglon+j];
          }
          std::cout << "\n";
        }
      }
      REQUIRE(tdata[i*nglon+j]==(idx*nglon+j));
      REQUIRE(vdata[i*nglon*2+2*j]==(idx*nglon*2+2*j));
      REQUIRE(vdata[i*nglon*2+2*j+1]==(idx*nglon*2+2*j+1));
      REQUIRE(idata[i*nglon+j]==(idx*nglon+j));
    }
  }
  int my_int;
  double my_dbl;
  std::string units;
  get_att (*file,"my_int","NC_GLOBAL",my_int);
  get_att (*file,"my_dbl","NC_GLOBAL",my_dbl);
  get_att (*file,"units","V",units);

  REQUIRE (my_int==42);
  REQUIRE (my_dbl==42);
  REQUIRE (units=="m/s");

  // Check read exceptions
  REQUIRE_THROWS (read_var(*file,"V",vdata.data(),1)); // Not a time-dep var
  REQUIRE_THROWS (read_var(*file,"I",idata.data(),2)); // Not a valid record
  REQUIRE_THROWS (read_var(*file,"I",vdata.data()));   // Wrong data type ptr

  close_file(*file);
}

TEST_CASE ("pnetcdf_io") {
  write_bad ();
  write ();
  read  ();
}

#include "profiling/stats/cldera_field_stat_factory.hpp"
#include "profiling/stats/cldera_field_bounding_box.hpp"
#include "profiling/stats/cldera_field_pnetcdf_reference.hpp"
#include "profiling/stats/cldera_field_zonal_mean.hpp"

#include <ekat/mpi/ekat_comm.hpp>

#include <catch2/catch.hpp>

#include <map>

TEST_CASE ("stats - pnetcdf") {
  using namespace cldera;

  ekat::Comm comm(MPI_COMM_WORLD);

  // open the file
  std::string pnetcdf_filename = "../../data/sample_database.nc";
  std::shared_ptr<io::pnetcdf::NCFile> file = io::pnetcdf::open_file(pnetcdf_filename,comm,io::pnetcdf::IOMode::Read);

  // grab relevant dims (time x lev x col)
  auto timedim = file->dims.at("time");
  auto levdim = file->dims.at("lev");
  auto coldim = file->dims.at("ncol");

  // set "realistic" field params using the dims
  const int ntime = timedim->len;
  const int nlev = levdim->len;
  const int ncol = coldim->len;

  // print them for my own sanity
  std::cout << "timedim->len=" << timedim->len << std::endl;
  std::cout << "levdim->len=" << levdim->len << std::endl;
  std::cout << "coldim->len=" << coldim->len << std::endl;

  // set a decomp along cols like E3SM (all on one proc for now)
  std::vector<int> my_cols;
  for (int i=0; i<coldim->len; ++i) {
    my_cols.push_back(i);
  }
  io::pnetcdf::add_decomp(*file,"ncol",my_cols);

  // grab relevant vars
  auto latvar = file->vars.at("lat"); // defined on cols
  auto lonvar = file->vars.at("lon"); // defined on cols
  auto areavar = file->vars.at("area"); // defined on cols
  auto timevar = file->vars.at("time"); // defined on time

  // print them for my own sanity
  std::cout << "latvar->dims=" << latvar->dims.size() << std::endl;
  std::cout << "latvar->dims[0]->name=" << latvar->dims[0]->name << std::endl;
  std::cout << "latvar->size=" << latvar->size << std::endl;
  std::cout << "latvar->dtype=" << latvar->dtype << std::endl;
  std::cout << "lonvar->dims=" << lonvar->dims.size() << std::endl;
  std::cout << "lonvar->dims[0]->name=" << lonvar->dims[0]->name << std::endl;
  std::cout << "lonvar->size=" << lonvar->size << std::endl;
  std::cout << "lonvar->dtype=" << lonvar->dtype << std::endl;
  std::cout << "areavar->dims=" << areavar->dims.size() << std::endl;
  std::cout << "areavar->dims[0]->name=" << areavar->dims[0]->name << std::endl;
  std::cout << "areavar->size=" << areavar->size << std::endl;
  std::cout << "areavar->dtype=" << areavar->dtype << std::endl;
  std::cout << "timevar->dims=" << timevar->dims.size() << std::endl;
  std::cout << "timevar->dims[0]->name=" << timevar->dims[0]->name << std::endl;
  std::cout << "timevar->size=" << timevar->size << std::endl;
  std::cout << "timevar->dtype=" << timevar->dtype << std::endl;

  // dump lat,lon into data
  std::vector<double>     latdata(ncol);
  std::vector<double>     londata(ncol);
  std::vector<double>     areadata(ncol);
  std::vector<double>     timedata(ntime);

  read_var(*file,"lat",latdata.data());
  read_var(*file,"lon",londata.data());
  read_var(*file,"area",areadata.data());
  // 
  for(int itime=0; itime<ntime; ++itime)
    read_var(*file,"time",&timedata[itime],itime);

  // // print lats for sanity (looks good)
  // std::cout << "lat=[";
  // for(auto l : latdata)
  //   std::cout << l << ", ";
  // std::cout << "]" << std::endl;

  // // print lons for sanity (looks good)
  // std::cout << "lon=[";
  // for(auto l : londata)
  //   std::cout << l << ", ";
  // std::cout << "]" << std::endl;

  // // print areas for sanity (looks good)
  // std::cout << "area=[";
  // for(auto a : areadata)
  //   std::cout << a << ", ";
  // std::cout << "]" << std::endl;

  // // print time for sanity (looks good)
  // std::cout << "time=[";
  // for(auto t : timedata)
  //   std::cout << t*86400 << ", ";
  // std::cout << "]" << std::endl;

  // now grab the T variable depending on the timestep
  auto Tvar = file->vars.at("T"); // defined on lev x cols, for each time
  std::vector<double> tdata(nlev*ncol);
  read_var(*file,"T",tdata.data(),0);

  // // print T for sanity (seems correct)
  // std::cout << "T=[";
  // for(auto t : tdata)
  //   std::cout << t << ", ";
  // std::cout << "]" << std::endl;

  // Allocate a very simple field
  int nparts = 1;
  int part_size = ncol / nparts;
  int part_dim = 2;
  Field foo("foo", {ntime,nlev,ncol}, {"lev", "dim", "ncol"}, nparts, part_dim);
  std::vector<std::vector<Real>> foo_data (nparts,std::vector<Real>(ntime*nlev*part_size));
  for (int i = 0; i < nparts; ++i) {
    std::iota(foo_data[i].begin(), foo_data[i].end(), i*ntime*nlev*part_size);
    foo.set_part_size(i, part_size);
    foo.set_part_data(i, foo_data[i].data());
  }
  foo.commit();

  // Allocate lat/lon
  std::shared_ptr<Field> lat(new Field("lat", {ncol}, {"ncol"}, nparts, 0));
  std::shared_ptr<Field> lon(new Field("lon", {ncol}, {"ncol"}, nparts, 0));
  for (int i = 0; i < nparts; ++i) {
    lat->set_part_size(i, part_size);
    lat->set_part_data(i, latdata.data());
    lon->set_part_size(i, part_size);
    lon->set_part_data(i, londata.data());
  }
  lat->commit();
  lon->commit();

  // make a dumb column partition for now
  std::vector<double>     colgiddata(ncol);
  std::iota(colgiddata.begin(), colgiddata.end(), ncol);
  std::shared_ptr<Field> col_gids(new Field("col_gids", {ncol}, {"ncol"}, nparts, 0));
  for (int i = 0; i < nparts; ++i) {
    col_gids->set_part_size(i, part_size);
    col_gids->set_part_data(i, colgiddata.data());
  }
  col_gids->commit();

  // Set pnetcdf stat parameters
  // auto pnetcdf_reference_pl = ekat::ParameterList("pnetcdf_reference");
  // pnetcdf_reference_pl.set<std::vector<Real>>("Latitude Bounds", {0.0, 1.0}); // latitude is measured in radians
  // pnetcdf_reference_pl.set<std::vector<Real>>("Longitude Bounds", {0.0, 1.0}); // longitude is measured in radians
  // pnetcdf_reference_pl.set<std::string>("Pnetcdf Filename", pnetcdf_filename); // use the same filename that we pulled earlier
  // pnetcdf_reference_pl.set<std::string>("Reference Field Name", "T"); // use the reference temperature from file in the computation
  // pnetcdf_reference_pl.set<std::string>("Reference Deviation Field Name", "T"); // TODO: technically the wrong thing to do, but works

  // Get the stat from the stat factory
  //auto stat = create_stat(pnetcdf_reference_pl, comm);
  //auto pnetcdf_reference_stat = dynamic_cast<FieldPnetcdfReference *>(stat.get());
  //REQUIRE_THROWS(pnetcdf_reference_stat->compute(foo)); // initialize() required
  //pnetcdf_reference_stat->initialize(lat, lon, col_gids); // give it the lat and lon
  //const auto pnetcdf_comparison_field = pnetcdf_reference_stat->compute(foo);

}

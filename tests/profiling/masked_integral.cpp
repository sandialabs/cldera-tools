#include "profiling/stats/cldera_field_stat.hpp"
#include "profiling/stats/cldera_register_stats.hpp"
#include "io/cldera_pnetcdf.hpp"

#include <ekat/mpi/ekat_comm.hpp>

#include <catch2/catch.hpp>

#include <map>

TEST_CASE ("vertical_contraction") {
  using namespace cldera;

  // REQUIRE_THROWS causes some start_timer calls to not be matched
  // by a corresponding stop_timer. When that happens, a subsequent
  // call to start_timer for that timer would throw. For this test,
  // we can just disable timings
  timing::TimingSession::instance().toggle_session(false);

  constexpr auto cp  = DataAccess::Copy;

  register_stats ();

  ekat::Comm comm(MPI_COMM_WORLD);

  // open the sample file
  std::string pnetcdf_filename = "../../data/ipcc_mask_ne4pg2.nc";
  auto file = io::pnetcdf::open_file(pnetcdf_filename,comm,io::pnetcdf::IOMode::Read);

  // Grab ncol from file
  const int ncol = file->dims.at("ncol")->len;
  const int nlevs = 10;
  const std::string lev_name = "lev";

  int ymd = 20220915;
  int tod = 43000;
  TimeStamp time(ymd,tod);

  // Read mask field from file
  Field mask("mask", {ncol}, {"ncol"}, DataAccess::Copy, DataType::IntType);
  mask.commit();
  auto mask_data = mask.data_nonconst<int>();
  io::pnetcdf::read_var(*file,"mask",mask_data);

  // Count how many cols are in each region
  std::map<int,int> col_count;
  for (int i=0; i<mask.layout().size(); ++i) {
    ++col_count[mask_data[i]];
  }
  const int num_regions = col_count.size();
  const auto stat_dim_name = "dim" + std::to_string(num_regions);

  // Helper lambda, to create/init a masked integral stat
  auto create_stat = [&] (const Field& f, const Field* w = nullptr) {
    ekat::ParameterList pl("masked_integral");
    pl.set<std::string>("type","masked_integral");
    pl.set<std::string>("mask_field","mask");

    std::map<std::string,Field> aux_fields;
    aux_fields["mask"] = mask;

    if (w!=nullptr) {
      pl.set<std::string>("weight_field",w->name());
      aux_fields[w->name()] = *w;
    }

    auto stat = StatFactory::instance().create("masked_integral",comm,pl);

    stat->set_field(f);
    stat->set_aux_fields(aux_fields);
    stat->create_stat_field ();
    return stat;
  };

  auto int2real = [](const Field& f, const std::string& name) {
    Field fr(name,f.layout(),DataAccess::Copy,DataType::RealType);
    fr.commit();
    auto in  = f.data<int>();
    auto out = fr.data_nonconst<Real>();
    for (int i=0; i<f.layout().size(); ++i) {
      out[i] = in[i];
    }
    return fr;
  };
  auto mask_real = int2real(mask,"w");

  auto extrude_mask = [&](const std::string& name, const int nlevs, const bool lev_dim_first) {
    const auto& ml = mask.layout();
    auto names = ml.names();
    auto dims  = ml.dims();

    if (lev_dim_first) {
      names.insert(names.begin(),lev_name);
      dims.insert(dims.begin(),nlevs);
    } else {
      names.push_back(lev_name);
      dims.push_back(nlevs);
    }
    FieldLayout fl(dims,names);
    Field f(name,fl,DataAccess::Copy,DataType::IntType);
    f.commit();

    auto fview = f.nd_view_nonconst<int,2>();
    auto mview = mask.nd_view_nonconst<int,1>();
    for (int icol=0; icol<ncol; ++icol) {
      for (int ilev=0; ilev<nlevs; ++ilev) {
        if (lev_dim_first) {
          fview(ilev,icol) = ilev*mview(icol);
        } else {
          fview(icol,ilev) = ilev*mview(icol);
        }
      }
    }
    return f;
  };

  // Also write results to file, so we can use it for downstream unit tests
  auto ofile = io::pnetcdf::open_file("results.nc",comm,io::pnetcdf::IOMode::Write);
  io::pnetcdf::add_dim(*ofile,"lev",nlevs,false);
  io::pnetcdf::add_dim(*ofile,"ncol",ncol,false);
  io::pnetcdf::add_dim(*ofile,stat_dim_name,num_regions,false);

  io::pnetcdf::add_time(*ofile,"double");
  io::pnetcdf::add_var(*ofile,"s2d_no_w","double",{stat_dim_name},true);
  io::pnetcdf::add_var(*ofile,"s2d_w","double",{stat_dim_name},true);
  io::pnetcdf::add_var(*ofile,"s3d_no_w","double",{lev_name,stat_dim_name},true);
  io::pnetcdf::add_var(*ofile,"s3d_w","double",{lev_name,stat_dim_name},true);
  io::pnetcdf::add_var(*ofile,"mask","double",{"ncol"},false);
  io::pnetcdf::add_var(*ofile,"mask_ids","int",{stat_dim_name},false);

  io::pnetcdf::enddef(*ofile);

  // Write mask, and mask values
  io::pnetcdf::write_var(*ofile,"mask",mask_real.data<Real>());
  std::vector<int> mask_ids;
  for (int i=0; i<ncol; ++i) {
    if (std::find(mask_ids.begin(),mask_ids.end(),mask_data[i])==mask_ids.end()) {
      mask_ids.push_back(mask_data[i]);
    }
  }
  io::pnetcdf::write_var(*ofile,"mask_ids",mask_ids.data());

  {
    // Integrate mask field itself
    auto stat = create_stat(mask);
    auto out = stat->compute(time);

    // Check layout
    REQUIRE (out.layout().rank()==1);
    REQUIRE (out.layout().dims()[0]==num_regions);
    REQUIRE (out.layout().names()[0]==stat_dim_name);

    // Check data. Unweighted integral of C is C*ncol_in_region
    auto out_data = out.nd_view<Real,1>();
    int idx = 0;
    for (auto it : col_count) {
      auto mid = it.first;
      auto n   = it.second;
      REQUIRE (out_data(idx) == mid*n);
      ++idx;
    }
    io::pnetcdf::write_var(*ofile,"s2d_no_w",out_data.data());
  }

  {
    // Integrate mask field itself, using mask itself as a weight
    auto stat = create_stat(mask,&mask_real);
    auto out = stat->compute(time);

    // Check layout
    REQUIRE (out.layout().rank()==1);
    REQUIRE (out.layout().dims()[0]==num_regions);
    REQUIRE (out.layout().names()[0]==stat_dim_name);

    // Check data. Integral of C with C as weight is C*C*ncol_in_region
    auto out_data = out.nd_view<Real,1>();
    int idx = 0;
    for (auto it : col_count) {
      auto mid = it.first;
      auto n   = it.second;
      REQUIRE (out_data(idx) == mid*mid*n);
      ++idx;
    }
    io::pnetcdf::write_var(*ofile,"s2d_w",out_data.data());
  }

  {
    // Integrate a field where f(lev,col) = lev*mask(col)
    for (bool lev_dim_first : {true, false}) { 
      auto f = extrude_mask ("f3d",nlevs,lev_dim_first);

      auto stat = create_stat(f);
      auto out = stat->compute(time);
      const auto& ol = out.layout();

      // Check layout
      REQUIRE (ol.rank()==2);

      REQUIRE (ol.dims()[0]==(lev_dim_first ? nlevs : num_regions));
      REQUIRE (ol.names()[0]==(lev_dim_first ? lev_name : stat_dim_name));

      REQUIRE (ol.dims()[1]==(lev_dim_first ? num_regions : nlevs));
      REQUIRE (ol.names()[1]==(lev_dim_first ? stat_dim_name : lev_name));

      // Check data. Integral of ilev*C is ilev*C*ncol_in_region
      auto out_data = out.nd_view<Real,2>();
      int idx = 0;
      for (auto it : col_count) {
        auto mid = it.first;
        auto n   = it.second;
        for (int ilev=0; ilev<nlevs; ++ilev) {
          if (lev_dim_first) {
            REQUIRE (out_data(ilev,idx) == ilev*mid*n);
          } else {
            REQUIRE (out_data(idx,ilev) == ilev*mid*n);
          }
        }
        ++idx;
      }
      if (lev_dim_first) {
        // Do this only for 'normal' e3sm layout
          io::pnetcdf::write_var(*ofile,"s3d_no_w",out_data.data());
      }
    }
  }

  {
    // Integrate a field where f(lev,col) = lev*mask(col), using mask itself as a weight
    for (bool lev_dim_first : {true, false}) { 
      auto f = extrude_mask ("f3d",nlevs,lev_dim_first);

      auto stat = create_stat(f,&mask_real);
      auto out = stat->compute(time);
      const auto& ol = out.layout();

      // Check layout
      REQUIRE (ol.rank()==2);

      REQUIRE (ol.dims()[0]==(lev_dim_first ? nlevs : num_regions));
      REQUIRE (ol.names()[0]==(lev_dim_first ? lev_name : stat_dim_name));

      REQUIRE (ol.dims()[1]==(lev_dim_first ? num_regions : nlevs));
      REQUIRE (ol.names()[1]==(lev_dim_first ? stat_dim_name : lev_name));

      // Check data. Integral of ilev*C with C as weight is ilev*C*C*ncol_in_region
      auto out_data = out.nd_view<Real,2>();
      int idx = 0;
      for (auto it : col_count) {
        auto mid = it.first;
        auto n   = it.second;
        for (int ilev=0; ilev<nlevs; ++ilev) {
          if (lev_dim_first) {
            REQUIRE (out_data(ilev,idx) == ilev*mid*mid*n);
          } else {
            REQUIRE (out_data(idx,ilev) == ilev*mid*mid*n);
          }
        }
        ++idx;
      }
      if (lev_dim_first) {
        // Do this only for 'normal' e3sm layout
          io::pnetcdf::write_var(*ofile,"s3d_w",out_data.data());
      }
    }
  }
  io::pnetcdf::update_time(*ofile,1.0);
  io::pnetcdf::close_file(*ofile);
}

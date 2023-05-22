#include "profiling/stats/cldera_field_vertical_contraction.hpp"
#include "profiling/stats/cldera_register_stats.hpp"

#include <ekat/mpi/ekat_comm.hpp>

#include <catch2/catch.hpp>

#include <map>

TEST_CASE ("vertical_contraction") {
  using namespace cldera;

  register_stats ();
  auto& factory = StatFactory::instance();

  // REQUIRE_THROWS causes some start_timer calls to not be matched
  // by a corresponding stop_timer. When that happens, a subsequent
  // call to start_timer for that timer would throw. For this test,
  // we can just disable timings
  timing::TimingSession::instance().toggle_session(false);

  ekat::Comm comm(MPI_COMM_WORLD);

  int ymd = 20220915;
  int tod = 43000;
  TimeStamp time(ymd,tod);

  // Field dimensions
  constexpr int nlevs = 4;
  constexpr int ndims = 3;
  constexpr int ncols = 4;
  constexpr int nparts = 2;
  constexpr int part_size = ncols / nparts;
  constexpr int lev_beg = 1;
  constexpr int lev_end = 2;
  constexpr int nlev_integral = lev_end - lev_beg + 1;
  constexpr double w_val = 2.0;
  constexpr auto cp = DataAccess::Copy;

  // Create input fields
  Field s2d("scalar2d",{ncols},{"ncol"},nparts, 0,cp);
  Field s3d_fwd("scalar3d",{nlevs,      ncols},{"lev",      "ncol"},nparts, 1,cp);
  Field v3d_fwd("vector3d",{nlevs,ndims,ncols},{"lev","dim","ncol"},nparts, 2,cp);
  Field s3d_bwd("scalar3d",{ncols,      nlevs},{"ncol",      "lev"},nparts, 0,cp);
  Field v3d_bwd("vector3d",{ncols,ndims,nlevs},{"ncol","dim","lev"},nparts, 0,cp);
  Field w2d_fwd ("w2d",{nlevs,ncols},{"lev","ncol"},nparts,1,cp);
  Field w2d_bwd ("w2d",{ncols,nlevs},{"ncol","lev"},nparts,0,cp);
  Field w1d ("w1d",FieldLayout({nlevs},{"lev"}),cp);
  for (int i = 0; i < nparts; ++i) {
    s2d.set_part_size(i, part_size);
    s3d_fwd.set_part_size(i, part_size);
    v3d_fwd.set_part_size(i, part_size);
    w2d_fwd.set_part_size(i, part_size);
    s3d_bwd.set_part_size(i, part_size);
    v3d_bwd.set_part_size(i, part_size);
    w2d_bwd.set_part_size(i, part_size);
  }
  s2d.commit();
  s3d_fwd.commit();
  v3d_fwd.commit();
  w2d_fwd.commit();
  s3d_bwd.commit();
  v3d_bwd.commit();
  w2d_bwd.commit();
  w1d.commit();
  w1d.deep_copy(w_val);
  w2d_fwd.deep_copy(w_val);
  w2d_bwd.deep_copy(w_val);

  for (int p = 0; p < nparts; ++p) {
    auto s3d_fwd_pv = s3d_fwd.part_nd_view_nonconst<Real,2>(p);
    auto v3d_fwd_pv = v3d_fwd.part_nd_view_nonconst<Real,3>(p);
    auto s3d_bwd_pv = s3d_bwd.part_nd_view_nonconst<Real,2>(p);
    auto v3d_bwd_pv = v3d_bwd.part_nd_view_nonconst<Real,3>(p);

    for (int col=0; col<part_size; ++col) {
      int gcol = p*part_size+col;
      for (int lev=0; lev<nlevs; ++lev) {
        s3d_fwd_pv(lev,col) = gcol*nlevs + lev;
        s3d_bwd_pv(col,lev) = gcol*nlevs + lev;
        for (int dim=0; dim<ndims; ++dim) {
          v3d_fwd_pv(lev,dim,col) = gcol*ndims*nlevs + dim*nlevs + lev;
          v3d_bwd_pv(col,dim,lev) = gcol*ndims*nlevs + dim*nlevs + lev;
        }
      }
    }
  }

  std::map<std::string,Field> aux_fields =
  {
    {"w1d", w1d}
  };

  // We test correctness by comparing global_sum(f)==vert_sum(horiz_sum)
  std::vector<int> lev_bounds{lev_beg,lev_end};
  for (const std::string& weight : {"NONE","w1d","w2d"}) {
    std::cout << " -> weight field: " << weight << "\n";
    const bool weighted = weight!="NONE";
    std::list<bool> w_const_vals = {true};
    if (weighted) {
      w_const_vals.push_back(false);
    }

    auto w = weighted ? w_val : 1.0;

    for (bool w_const : w_const_vals) {
      std::cout << "   -> const weight: " << (w_const ? "yes" : "no") << "\n";
      for (bool average : {false,true}) {
        std::cout << "     -> average: " << (average ? "yes" : "no") << "\n";
        // Create stat
        ekat::ParameterList pl("vert_int");
        pl.set("level_bounds",lev_bounds);
        pl.set("average",average);
        pl.set("weight_field",weight);
        pl.set("constant_weight",w_const);

        auto stat = factory.create("vertical_contraction",comm,pl);

        // Input field does not have vertical dim
        REQUIRE_THROWS (stat->set_field(s2d));

        std::cout << "        - 3d scalar field (nlevs,ncols) ..." << std::endl;
        stat->set_field(s3d_fwd);
        if (weighted) {
          aux_fields["w2d"] = w2d_fwd;
          stat->set_aux_fields(aux_fields);
        }
        auto s3d_fwd_vi = stat->compute(time);

        REQUIRE (s3d_fwd_vi.layout().rank()==1);
        REQUIRE (s3d_fwd_vi.layout().size()==ncols);

        auto s3d_fwd_vi_view = s3d_fwd_vi.nd_view<Real,1>();
        for (int col=0; col<ncols; ++col) {
          Real tgt = 0;
          for (int lev=lev_bounds[0];lev<=lev_bounds[1]; ++lev) {
            tgt += (col*nlevs + lev) * w;
          }

          if (average) {
            tgt /= nlev_integral * w;
          }
          REQUIRE (s3d_fwd_vi_view(col)==tgt);
        }
        std::cout << "          3d scalar field (nlevs,ncols) ... OK!" << std::endl;

        // Check s3d_bwd against s3d_fwd
        std::cout << "        - 3d scalar field (ncols,nlevs) ..." << std::endl;
        stat->set_field(s3d_bwd);
        if (weighted) {
          aux_fields["w2d"] = w2d_bwd;
          stat->set_aux_fields(aux_fields);
        }
        auto s3d_bwd_vi = stat->compute(time);

        REQUIRE (s3d_bwd_vi.layout().rank()==1);
        REQUIRE (s3d_bwd_vi.layout().size()==ncols);

        auto s3d_bwd_vi_view = s3d_bwd_vi.nd_view<Real,1>();
        for (int col=0; col<ncols; ++col) {
          REQUIRE (s3d_bwd_vi_view(col)==s3d_fwd_vi_view(col));
        }
        std::cout << "        - 3d scalar field (ncols,nlevs) ... OK!" << std::endl;

        // Test v3d
        std::cout << "        - 3d vector field (nlevs,ndims,ncols) ..." << std::endl;
        stat->set_field(v3d_fwd);
        if (weighted) {
          aux_fields["w2d"] = w2d_fwd;
          stat->set_aux_fields(aux_fields);
        }

        auto v3d_fwd_vi = stat->compute(time);

        REQUIRE (v3d_fwd_vi.layout().rank()==2);
        REQUIRE (v3d_fwd_vi.layout().dims()==std::vector<int>{ndims,ncols});

        auto v3d_fwd_vi_view = v3d_fwd_vi.nd_view<Real,2>();
        for (int col=0; col<ncols; ++col) {
          for (int dim=0; dim<ndims; ++dim) {
            Real tgt = 0;
            for (int lev=lev_bounds[0];lev<=lev_bounds[1]; ++lev) {
              tgt += (col*ndims*nlevs + dim*nlevs + lev) * w;
            }

            if (average) {
              tgt /= nlev_integral * w;
            }
            REQUIRE (v3d_fwd_vi_view(dim,col)==tgt);
          }
        }
        std::cout << "        - 3d vector field (nlevs,ndims,ncols) ... OK!" << std::endl;

        // Check v3d_bwd against v3d_fwd
        std::cout << "        - 3d vector field (ncols,ndims,nlevs) ..." << std::endl;
        stat->set_field(v3d_bwd);
        if (weighted) {
          aux_fields["w2d"] = w2d_bwd;
          stat->set_aux_fields(aux_fields);
        }

        auto v3d_bwd_vi = stat->compute(time);

        REQUIRE (v3d_bwd_vi.layout().rank()==2);
        REQUIRE (v3d_bwd_vi.layout().dims()==std::vector<int>{ncols,ndims});

        auto v3d_bwd_vi_view = v3d_bwd_vi.nd_view<Real,2>();
        for (int col=0; col<ncols; ++col) {
          for (int dim=0; dim<ndims; ++dim) {
            REQUIRE (v3d_bwd_vi_view(col,dim)==v3d_fwd_vi_view(dim,col));
          }
        }
        std::cout << "        - 3d vector field (ncols,ndims,nlevs) ... OK!" << std::endl;
      }
    }
  }
}

#include "profiling/stats/cldera_field_stat.hpp"
#include "profiling/stats/cldera_field_identity.hpp"
#include "profiling/stats/cldera_register_stats.hpp"

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

  register_stats ();

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
  constexpr auto cp = DataAccess::Copy;

  const int mask_value = 42;
  const Real w_value = 100.0;

  // Create input fields, mask field, and masked input field (f=0 if mask!=mask_value)
  Field s2d("scalar2d", {            ncols}, {              "ncol"}, nparts, 0,cp);
  Field s3d("scalar3d", {nlevs,      ncols}, {"lev",        "ncol"}, nparts, 1,cp);
  Field v3d("vector3d", {nlevs,ndims,ncols}, {"lev", "dim", "ncol"}, nparts, 2,cp);
  
  Field mask("mask", {ncols}, {"ncol"}, nparts, 0,cp,DataType::IntType);
  Field weight("weight", {ncols}, {"ncol"}, nparts, 0,cp);

  // Create copy of input fields, where f=0 where mask!=mask_val
  Field s2d_masked("masked_scalar2d", {            ncols}, {              "ncol"}, nparts, 0,cp);
  Field s3d_masked("masked_scalar3d", {nlevs,      ncols}, {"lev",        "ncol"}, nparts, 1,cp);
  Field v3d_masked("masked_vector3d", {nlevs,ndims,ncols}, {"lev", "dim", "ncol"}, nparts, 2,cp);

  for (int p = 0; p < nparts; ++p) {
    s2d.set_part_size(p, part_size);
    s3d.set_part_size(p, part_size);
    v3d.set_part_size(p, part_size);

    mask.set_part_size(p, part_size);
    weight.set_part_size(p, part_size);

    s2d_masked.set_part_size(p, part_size);
    s3d_masked.set_part_size(p, part_size);
    v3d_masked.set_part_size(p, part_size);
  }
  s2d.commit();
  s3d.commit();
  v3d.commit();
  mask.commit();
  weight.commit();
  s2d_masked.commit();
  s3d_masked.commit();
  v3d_masked.commit();

  for (int p = 0; p < nparts; ++p) {
    std::vector<Real> data;
    // s2d
    data.resize(part_size);
    std::iota(data.begin(), data.end(), p*part_size);
    s2d.copy_part_data(p, data.data());

    // s3d
    data.resize(part_size*nlevs);
    for (int lev=0; lev<nlevs; ++lev) {
      auto beg = data.begin()+lev*part_size;
      auto end = beg + part_size;
      std::iota(beg,end,lev*ncols + p*part_size);
    }
    s3d.copy_part_data(p, data.data());

    // v3d
    data.resize(part_size*nlevs*ndims);
    for (int lev=0; lev<nlevs; ++lev) {
      for (int dim=0; dim<ndims; ++dim) {
        auto beg = data.begin()+lev*ndims*part_size + dim*part_size;
        auto end = beg + part_size;
        std::iota(beg,end,lev*ndims*ncols + dim*ncols + p*part_size);
      }
    }
    v3d.copy_part_data(p, data.data());

    // Set m=mask_value if the col index is even
    auto m = mask.part_view_nonconst<int>(p);
    auto w = weight.part_view_nonconst<Real>(p);
    for (int i=0; i<part_size; ++i) {
      int icol = i+part_size*p;
      m(i) = icol%2==0 ? mask_value : 1;
      w(i) = w_value;
    }

    // Set masked field to 0 if m!=mask_value
    auto s2d_pv = s2d.part_nd_view<Real,1>(p);
    auto s3d_pv = s3d.part_nd_view<Real,2>(p);
    auto v3d_pv = v3d.part_nd_view<Real,3>(p);
    auto ms2d_pv = s2d_masked.part_nd_view_nonconst<Real,1>(p);
    auto ms3d_pv = s3d_masked.part_nd_view_nonconst<Real,2>(p);
    auto mv3d_pv = v3d_masked.part_nd_view_nonconst<Real,3>(p);
    for (int i=0; i<part_size; ++i) {
      ms2d_pv(i) = m(i)==mask_value ? s2d_pv(i) : 0;
      for (int k=0; k<nlevs; ++k) {
        ms3d_pv(k,i) = m(i)==mask_value ? s3d_pv(k,i) : 0;
        for (int j=0; j<ndims; ++j) {
          mv3d_pv(k,j,i) = m(i)==mask_value ? v3d_pv(k,j,i) : 0;
        }
      }
    }
  }

  auto& s = StatFactory::instance();
  auto sum_stat = s.create("sum_along_columns",comm,ekat::ParameterList(""));

  // Compute target values (unweighted), by doing a global_sum of
  // s2d_masked, s3d_masked, v3d_masked
  std::map<std::string,Field> tgt_stats;
  sum_stat->set_field(s2d_masked);
  tgt_stats[s2d.name()] = sum_stat->compute(time);
  sum_stat->set_field(s3d_masked);
  tgt_stats[s3d.name()] = sum_stat->compute(time);
  sum_stat->set_field(v3d_masked);
  tgt_stats[v3d.name()] = sum_stat->compute(time);

  std::map<std::string,Field> aux_fields =
  {
    { "weight", weight },
    { "mask"  , mask   }
  };

  for (auto weighted : {false, true}) {
    std::cout << " -> use weight: " << (weighted ? "yes\n" : "no\n");
    // Create masked integral stat
    ekat::ParameterList pl("masked_integral");
    pl.set<std::string>("type","masked_integral");
    pl.set<std::string>("mask_dim_name","ncol");
    pl.set<std::string>("mask_field","mask");
    pl.set("mask_value",mask_value);

    if (weighted) {
      pl.set<std::string>("weight_field","weight");
      pl.set("weight_field_constant",true);
    }

    auto stat = s.create("masked_integral",comm,pl);

    std::map<std::string,Field> stat_value;
    for (auto f : { s2d, s3d, v3d} ) {
      std::cout << "   -> field: " << f.name() << " ...\n";

      stat->set_field(f);
      stat->set_aux_fields(aux_fields);
      auto stat_f = stat->compute(time);
      const auto& tgt_f = tgt_stats.at(f.name());

      REQUIRE (stat_f.layout()==tgt_f.layout());

      auto data = stat_f.data<Real>();
      auto tgt_data = tgt_f.data<Real>();

      int sz = tgt_f.layout().size();
      for (int i=0; i<sz; ++i) {
        auto factor = weighted ? w_value : 1;
        REQUIRE (data[i]==tgt_data[i]*factor);
      }
      std::cout << "   -> field: " << f.name() << " ... OK!\n";
    }
  }
}

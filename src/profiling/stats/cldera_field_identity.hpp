#ifndef CLDERA_FIELD_IDENTITY_HPP
#define CLDERA_FIELD_IDENTITY_HPP

#include "cldera_field_stat.hpp"

namespace cldera
{

class FieldIdentity : public FieldStat
{
public:
  FieldIdentity (const ekat::Comm& comm,
                 const ekat::ParameterList& pl)
   : FieldStat(comm,pl)
  { /* Nothing to do here */ }

  std::string type () const { return "identity"; }

  // Given a field, return the layout that the computed stat will have
  FieldLayout stat_layout (const FieldLayout& field_layout) const {
    return field_layout;
  }

  void compute_impl (const Field& f, Field& stat) const {
    const auto dt = f.data_type();
    if (dt==IntType) {
      do_compute_impl<int>(f,stat);
    } else if (dt==RealType) {
      do_compute_impl<Real>(f,stat);
    } else {
      EKAT_ERROR_MSG ("[FieldIdentity] Unrecognized/unsupported data type (" + e2str(dt) + ")\n");
    }
  }
private:
  template<typename T, int N>
  using view_Nd_host = typename KokkosTypesHost::template view_ND<T,N>;

  template<typename T>
  void do_compute_impl (const Field& f, Field& stat) const {
    using pair_t = Kokkos::pair<int,int>;
    constexpr auto ALL = Kokkos::ALL();

    auto stat_data = stat.data_nonconst<T>();
    const int rank = f.layout().rank();
    const auto& fdims = f.layout().dims();

    const int nparts = f.nparts();
    const int part_dim = f.part_dim();
    int part_start = 0;
    for (int p=0; p<nparts; ++p) {
      const auto& part_layout = f.part_layout(p);
      const int part_size = part_layout.dims()[part_dim];
      const int part_end = part_start + part_size;
      auto fpart_data = f.part_data<const T>(p);
      switch (rank) {
        case 1:
        {
          view_Nd_host<T,1> stat_view (stat_data,fdims[0]);
          view_Nd_host<const T,1> fpart_view (fpart_data,part_layout.kokkos_layout());
          auto stat_subview = Kokkos::subview(stat_view,pair_t(part_start,part_end));
          Kokkos::deep_copy(stat_subview,fpart_view);
          break;
        }
        case 2:
        {
          view_Nd_host<T,2> stat_view (stat_data,fdims[0],fdims[1]);
          view_Nd_host<const T,2> fpart_view (fpart_data,part_layout.kokkos_layout());
          if (part_dim==0) {
            auto stat_subview = Kokkos::subview(stat_view,pair_t(part_start,part_end),ALL);
            Kokkos::deep_copy(stat_subview,fpart_view);
          } else {
            auto stat_subview = Kokkos::subview(stat_view,ALL,pair_t(part_start,part_end));
            Kokkos::deep_copy(stat_subview,fpart_view);
          }
          break;
        }
        case 3:
        {
          view_Nd_host<T,3> stat_view (stat_data,fdims[0],fdims[1],fdims[2]);
          view_Nd_host<const T,3> fpart_view (fpart_data,part_layout.kokkos_layout());
          if (part_dim==0) {
            auto stat_subview = Kokkos::subview(stat_view,pair_t(part_start,part_end),ALL,ALL);
            Kokkos::deep_copy(stat_subview,fpart_view);
          } else if (part_dim==1) {
            auto stat_subview = Kokkos::subview(stat_view,ALL,pair_t(part_start,part_end),ALL);
            Kokkos::deep_copy(stat_subview,fpart_view);
          } else {
            auto stat_subview = Kokkos::subview(stat_view,ALL,ALL,pair_t(part_start,part_end));
            Kokkos::deep_copy(stat_subview,fpart_view);
          }
          break;
        }
        default:
          EKAT_ERROR_MSG ("[FieldIdentity] Unsupported field rank\n"
              "  - field name: " + f.name() + "\n"
              "  - field rank: " + std::to_string(rank) + "\n");
      }
      part_start = part_end;
    }
  }
};

} // namespace cldera

#endif // CLDERA_FIELD_IDENTITY_HPP

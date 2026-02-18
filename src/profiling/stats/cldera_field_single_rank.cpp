#include "cldera_field_single_rank.hpp"
#include "profiling/utils/cldera_subview_utils.hpp"
#include "profiling/cldera_mpi_timing_wrappers.hpp"

namespace cldera
{

void FieldSingleRank::
create_stat_field () {

  // initialize to zero
  m_local_size = 0;
  m_global_size = 0;
  m_part_sizes = std::vector<int>();

  // first, find out how many columns we need in total
  const int nparts = m_field.nparts();
  const int part_dim = m_field.part_dim();
  auto non_part_layout = m_field.layout().strip_dim(part_dim);

  // loop over parts and grab part dims
  for (int p=0; p<nparts; ++p) {
    const auto& part_layout = m_field.part_layout(p);
    const int part_size = part_layout.dims()[part_dim];
    m_part_sizes.push_back(part_size);
    m_local_size += part_size;
  }

  // count up global columns
  track_mpi_all_reduce(m_comm,&m_local_size,&m_global_size,1,MPI_SUM, name());

  // create a new layout which is size m_global_size on rank 0 and size 0 on all other ranks
  auto m_stat_layout = m_field.layout();
  // if (m_comm.am_i_root()) {
  //   m_stat_layout[m_stat_layout.dim_idx("ncol")] = m_global_size;
  // } else {
  //   m_stat_layout[m_stat_layout.dim_idx("ncol")] = 0;
  // }

  m_stat_field = Field(name(), m_stat_layout, DataAccess::Copy, stat_data_type());
}

void FieldSingleRank::
compute_impl () {
  const auto dt = m_field.data_type();
  const int rank = m_field.layout().rank();

  // assert we're a single rank because otherwise we'll have to think harder
  EKAT_REQUIRE_MSG (rank<=1,
      "[FieldSingleRank] Unsupported field rank.\n"
      " - field name: " + m_field.name() + "\n"
      " - field rank: " + std::to_string(rank) + "\n");

  if (dt==IntType) {
    if (rank==1) {
      do_compute_impl<int,1>();
    }
  } else if (dt==RealType) {
    if (rank==1) {
      do_compute_impl<Real,1>();
    }
  } else {
    EKAT_ERROR_MSG (
        "[FieldSingleRank] Unrecognized/unsupported data type\n"
        " - field name: " + m_field.name() + "\n"
        " - data type : " + e2str(dt) + "\n");
  }
}

template<typename T, int N>
void FieldSingleRank::
do_compute_impl ()
{
  auto stat_view = m_stat_field.nd_view_nonconst<T,N>();

  // second, fill a buffer with our local data to send to the root rank
  T* local_data;

  // first, find out how many columns we need in total
  const int nparts = m_field.nparts();
  const int part_dim = m_field.part_dim();
  auto non_part_layout = m_field.layout().strip_dim(part_dim);

  // loop over parts and grab part dims
  for (int p=0; p<nparts; ++p) {
    const auto& part_layout = m_field.part_layout(p);
    const int part_size = part_layout.dims()[part_dim];
    const int part_offset = m_field.part_offset(p);
    auto fpart_view = m_field.part_nd_view<const T,N>(p);
    
    
    // TODO: copy to local_data
  }

  // finally, send the buffer to the root rank and store it in stat_field
  // TODO: do this

  // TODO: delete this
  // for (int p=0; p<nparts; ++p) {
  //   const auto& part_layout = m_field.part_layout(p);
  //   const int part_size = part_layout.dims()[part_dim];
  //   const int part_offset = m_field.part_offset(p);
  //   auto fpart_view = m_field.part_nd_view<const T,N>(p);
  //   for (int i=0; i<part_size; ++i) {
  //     if constexpr (N==1) {
  //       stat_view(i+part_offset) = fpart_view(i);
  //     } else {
  //       auto stat_slice = slice(stat_view,part_dim,i);
  //       auto f_slice    = slice(fpart_view,part_dim,i);
  //       for (int j=0; j<non_part_layout.extent(0); ++j) {
  //         if constexpr (N==2) {
  //           stat_slice(j) = f_slice(j);
  //         } else {
  //           for (int k=0; k<non_part_layout.extent(0); ++k) {
  //             stat_slice(j,k) = f_slice(j,k);
  //           }
  //         }
  //       }
  //     }
  //   }
  // }
}

} // namespace cldera

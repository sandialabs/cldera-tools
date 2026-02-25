#include "cldera_field_single_rank.hpp"
#include "profiling/utils/cldera_subview_utils.hpp"
#include "profiling/cldera_mpi_timing_wrappers.hpp"

#include <mpi.h>

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

  // loop over parts and grab part dims, counting the local size
  for (int p=0; p<nparts; ++p) {
    const auto& part_layout = m_field.part_layout(p);
    const int part_size = part_layout.dims()[part_dim];
    m_part_sizes.push_back(part_size);
    m_local_size += part_size;
  }

  // sum up local columns for global columns
  track_mpi_all_reduce(m_comm,&m_local_size,&m_global_size,1,MPI_SUM, name());

  // if we're root, we have a m_global_size columns, otherwise there's no data
  // TODO: if the data has other dimensions (e.g. lev, cmp) we need to make sure we allocate the right size
  if (m_comm.am_i_root()) {
    auto m_stat_layout = m_field.layout();
    //m_stat_layout[m_stat_layout.dim_idx("ncol")] = m_global_size; // TODO: resize this layout correctly
    m_stat_field = Field(name(),m_stat_layout,DataAccess::Copy,stat_data_type());
  } else {
    m_stat_field = Field(name(),{},{},DataAccess::Copy);
  }

  // commit the data, which allocates everything and finishes setup
  m_stat_field.commit();

}

void FieldSingleRank::
compute_impl () {
  const auto dt = m_field.data_type();
  const int rank = m_field.layout().rank();

  // TODO: check this for the rank>=2 case
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

  // fill a buffer with all part data combined to send to the root rank
  T* local_data; // TODO: there is likely a smarter way to do this with hijacking pointers to view data
  local_data = (T*)malloc(m_local_size*sizeof(T));

  T* root_data;
  if (m_comm.am_i_root()) {
    root_data = (T*)malloc(m_global_size*sizeof(T));
  }

  // first, find out how many columns we need in total
  const int nparts = m_field.nparts();
  const int part_dim = m_field.part_dim();
  auto non_part_layout = m_field.layout().strip_dim(part_dim);

  int offset = 0;

  // loop over parts and grab part dims
  for (int p=0; p<nparts; ++p) {
    const auto& part_layout = m_field.part_layout(p);
    const int part_size = part_layout.dims()[part_dim];
    const int part_offset = m_field.part_offset(p);
    auto fpart_view = m_field.part_nd_view<const T,N>(p);
    
    // TODO: copy each part to local_view
    for (int i=0; i<fpart_view.extent(0); ++i) {
      local_data[offset+i] = fpart_view(part_offset+i);
    }
    offset += part_size;
  }

  // very standard MPI setup for the gatherv since ekat doesn't have a wrapper for it
  // send all sizes to the root rank, compute offsets on the root rank based on the sizes
  int* recvcounts;
  int* displs;
  if (m_comm.am_i_root()) {
    recvcounts = (int*)malloc(m_comm.size() * sizeof(int));
    displs = (int*)malloc(m_comm.size() * sizeof(int));
  }
  MPI_Gather(&m_local_size, 1, MPI_INT, recvcounts, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (m_comm.am_i_root()) {
    displs[0] = 0;
    for (int i=1; i<m_comm.size(); i++) {
      displs[i] = displs[i-1] + recvcounts[i-1];
    }
  }

  // send the buffer to the root rank and store it in stat_field
  MPI_Datatype mpi_type = MPI_DOUBLE; // TODO: allow any type
  MPI_Gatherv(local_data, m_local_size, mpi_type, root_data, recvcounts, displs, mpi_type, 0, MPI_COMM_WORLD);

  // finally, copy that information in the stat field


  // TODO: fix memory leaks
  if (m_comm.am_i_root()) {
    free(recvcounts);
    free(displs);
    free(local_data);
  }
}

} // namespace cldera

#ifndef CLDERA_FIELD_HPP
#define CLDERA_FIELD_HPP

#include "cldera_profiling_types.hpp"
#include <ekat/kokkos/ekat_kokkos_types.hpp>

#include <vector>
#include <string>

namespace cldera
{

/*
 * A small struct holding extents of a field
 *
 * It is basically a std::vector<int>, with a few utilities
 */

class FieldLayout {
public:
  FieldLayout (const std::vector<int>& dims) {
    EKAT_REQUIRE_MSG (dims.size()>0, "Error! Invalid rank.\n");
    for (auto d : dims) {
      EKAT_REQUIRE_MSG (d>0, "Error! Invalid dimension (" + std::to_string(d) + "\n");
    }
    m_dims = dims;
  }

  const std::vector<int>& dims () const { return m_dims; }

  int operator[]  (const int i) const { return m_dims[i]; }

  int rank () const { return dims().size(); }
  long long size () const {
    long long s = 1;
    for (auto d : dims()) {
      s *= d;
    }
    return s;
  }
private:
  std::vector<int> m_dims;
};

/*
 * A small struct holding data and metadata of a field
 *
 * The basic metadata is given by the 'name' and 'dims' members.
 *
 * The Field is assumed to be potentially partitioned in memory.
 * The partitioning is described by
 *  - nparts: number of partitions
 *  - part_dim: index of dimension along which field is partitioned
 *  - part_sizes: the size of each partition
 *
 * Notice that each partition is stored contiguously in memory.
 * The method part_layout allows to get the layout of the partition,
 * so that one can safely iterate over the partition.
 *
 * Note: a contiguous/non-partitioned field simply has nparts=1.
 */

class Field {
public:
  using KT_h = ekat::KokkosTypes<ekat::HostDevice>;
  template<typename T, typename MT = Kokkos::MemoryManaged>
  using hview_1d = typename KT_h::template view_1d<T,MT>;

  Field (const std::string& n, const std::vector<int>& d,
         const int nparts, const int part_dim);

  // Shortcuts for single-partition field
  Field (const std::string& n, const std::vector<int>& d);
  Field (const std::string& n, const std::vector<int>& d, const Real* data);

  const std::string& name () const { return m_name; }

  // Get rank-global (not partitioned) and part layouts
  const FieldLayout& layout () const { return m_layout; }
        FieldLayout part_layout (const int ipart) const;

  // Set part specs
  void set_part_size  (const int ipart, const int part_size);
  void set_part_data  (const int ipart, const Real* data);
  void set_part_size_and_data  (const int ipart, const int part_size, const Real* data);
  void set_data (const Real* data);

  void commit ();

  // Get data as view
  hview_1d<const Real> get_part_view (const int ipart) const;
  hview_1d<const Real> get_view () const;

  // Get raw data
  const Real* get_part_data (const int ipart) const;
  const Real* get_data () const;

  // Query status
  int nparts () const { return m_nparts; }
  bool committed () const { return m_committed; }

private:

  void check_part_idx (const int i) const;

  std::string                       m_name;
  FieldLayout                       m_layout;
  int                               m_nparts;
  int                               m_part_dim;
  std::vector<long long>            m_part_sizes;
  std::vector<hview_1d<const Real>> m_data;
  bool                              m_committed = false;
};

inline const Real*
Field::get_part_data (const int ipart) const {
  return get_part_view(ipart).data();
}

inline const Real*
Field::get_data () const {
  return get_view().data();
}

inline auto
Field::get_part_view (const int ipart) const
 -> hview_1d<const Real>
{
  EKAT_REQUIRE_MSG (m_committed,
      "Error! Field '" + m_name + "' was not committed yet.\n");

  check_part_idx(ipart);

  return m_data[ipart];
}

inline auto
Field::get_view () const
 -> hview_1d<const Real>
{
  EKAT_REQUIRE_MSG (m_nparts==1,
      "Error! Field::get_view is only available for non-partitioned fields.\n");

  return get_part_view(0);
}

} // namespace cldera

#endif // CLDERA_FIELD_HPP

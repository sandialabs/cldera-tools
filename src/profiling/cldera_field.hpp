#ifndef CLDERA_FIELD_HPP
#define CLDERA_FIELD_HPP

#include "cldera_field_layout.hpp"
#include "cldera_profiling_types.hpp"

#include <ekat/ekat_assert.hpp>

#include <vector>
#include <string>

namespace cldera
{

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

enum class DataAccess {
  Copy,
  View
};

class Field {
public:
  // Not really needed for actual work, but it allows e.g. to use op[] on std::map
  Field () = default;

  Field (const std::string& n, const FieldLayout& fl,
         const int nparts, const int part_dim,
         const DataAccess cv = DataAccess::View);
  Field (const std::string& n,
         const std::vector<int>& dims,
         const std::vector<std::string>& dimnames,
         const int nparts, const int part_dim,
         const DataAccess cv = DataAccess::View);

  // Shortcuts for single-partition field
  Field (const std::string& n, const FieldLayout& fl,
         const DataAccess cv = DataAccess::View);
  Field (const std::string& n,
         const std::vector<int>& dims,
         const std::vector<std::string>& dimnames,
         const DataAccess cv = DataAccess::View);
  Field (const std::string& n, const FieldLayout& fl, const Real* data);
  Field (const std::string& n,
         const std::vector<int>& dims,
         const std::vector<std::string>& dimnames,
         const Real* data);

  const std::string& name () const { return m_name; }

  // Get rank-global (not partitioned) and part layouts
  const FieldLayout& layout () const { return m_layout; }
        FieldLayout part_layout (const int ipart) const;

  // Set part specs
  void set_part_size  (const int ipart, const int part_size);
  void set_part_data  (const int ipart, const Real* data);
  void set_data (const Real* data);

  void commit ();

  // Copy into managed views, if m_data_access=Copy
  void copy_part_data (const int ipart, const Real* data);
  void copy_data (const Real* data);

  // Get data as view
  view_1d_host<const Real> part_view (const int ipart) const;
  view_1d_host<      Real> part_view_nonconst (const int ipart);
  view_1d_host<const Real> view () const;
  view_1d_host<      Real> view_nonconst ();

  // Get raw data
  const Real* part_data (const int ipart) const;
        Real* part_data_nonconst (const int ipart);
  const Real* data () const;
        Real* data_nonconst ();

  // Query status
  int nparts () const { return m_nparts; }
  bool committed () const { return m_committed; }

private:

  void check_part_idx (const int i) const;

  std::string             m_name;
  FieldLayout             m_layout;
  int                     m_nparts   = -1;  // Set to something invalid for default ctor
  int                     m_part_dim = -1;
  std::vector<long long>  m_part_sizes;
  bool                    m_committed = false;
  DataAccess              m_data_access;

  std::vector<view_1d_host<const Real>>   m_data;
  std::vector<view_1d_host<      Real>>   m_data_nonconst;
};

// =================== IMPLEMENTATION =================== //

inline const Real*
Field::part_data (const int ipart) const {
  return part_view(ipart).data();
}

inline Real*
Field::part_data_nonconst (const int ipart) {
  return part_view_nonconst(ipart).data();
}

inline const Real*
Field::data () const {
  return view().data();
}

inline Real*
Field::data_nonconst () {
  return view_nonconst().data();
}

inline  view_1d_host<const Real>
Field::part_view (const int ipart) const
{
  EKAT_REQUIRE_MSG (m_committed,
      "Error! Field '" + m_name + "' was not committed yet.\n");

  check_part_idx(ipart);

  return m_data[ipart];
}

inline  view_1d_host<Real>
Field::part_view_nonconst (const int ipart)
{
  EKAT_REQUIRE_MSG (m_committed,
      "Error! Field '" + m_name + "' was not committed yet.\n");

  check_part_idx(ipart);

  return m_data_nonconst[ipart];
}

inline view_1d_host<const Real>
Field::view () const
{
  EKAT_REQUIRE_MSG (m_nparts==1,
      "Error! Field::view is only available for non-partitioned fields.\n");

  return part_view(0);
}

inline view_1d_host<Real>
Field::view_nonconst ()
{
  EKAT_REQUIRE_MSG (m_nparts==1,
      "Error! Field::view is only available for non-partitioned fields.\n");

  return part_view_nonconst(0);
}

} // namespace cldera

#endif // CLDERA_FIELD_HPP

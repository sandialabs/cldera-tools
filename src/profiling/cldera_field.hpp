#ifndef CLDERA_FIELD_HPP
#define CLDERA_FIELD_HPP

#include "cldera_profiling_types.hpp"

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
 *  - part_beg: where each partition starts along the $part_dim dimension
 *
 * Notice that each partition is stored contiguously in memory.
 * The array part_beg is such that partition i corresponds to indices
 * [part_beg[i], part_beg[i+1]) in the unpartitioned field.
 * The method part_dims allows to get the dimensions of the partition,
 * so that one can safely iterate over the partition.
 *
 * Note: a contiguous/non-partitioned field simply has nparts=1.
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

  int rank () const { return dims().size(); }
  int size () const {
    int s = 1;
    for (auto d : dims()) {
      s *= d;
    }
    return s;
  }
private:
  std::vector<int> m_dims;
};

class Field {
public:
  Field (const std::string& n, const std::vector<int>& d,
         const int nparts, const int part_dim);

  // Shortcuts for single-partition field
  Field (const std::string& n, const std::vector<int>& d);
  Field (const std::string& n, const std::vector<int>& d, const Real* data);

  const std::string& name () const { return m_name; }

  const FieldLayout& layout () const { return m_layout; }

  FieldLayout part_layout (const int p) const {
    EKAT_REQUIRE_MSG (m_committed,
        "Error! Field '" + m_name + "' was not yet committed.\n");
    check_part_idx(p);

    std::vector<int> d = m_layout.dims();
    d[m_part_dim] = m_part_beg[p+1] - m_part_beg[p];
    return FieldLayout{d};
  }

  void set_part_data (const int i, const int beg, const Real* data);
  // Shortcut for single-partition field only
  void set_data (const Real* data);

  int nparts () const { return m_nparts; }
  void commit ();

  const Real* get_part_data (const int i) const;
  // Shortcut for single-partition field only
  const Real* get_data () const;

  bool committed () const { return m_committed; }

private:
  void check_part_idx (const int i) const;

  std::string               m_name;
  FieldLayout               m_layout;
  int                       m_nparts;
  int                       m_part_dim;
  std::vector<int>          m_part_beg;
  std::vector<const Real*>  m_data;
  bool                      m_committed = false;
};

inline const Real*
Field::get_part_data (const int i) const {
  EKAT_REQUIRE_MSG (m_committed,
      "Error! Field '" + m_name + "' was not committed yet.\n");

  check_part_idx(i);

  return m_data[i];
}

inline const Real*
Field::get_data () const {
  EKAT_REQUIRE_MSG (m_nparts==1,
      "Error! Field::get_data is only available for non-partitioned fields.\n");

  return get_part_data(0);
}

} // namespace cldera

#endif // CLDERA_FIELD_HPP

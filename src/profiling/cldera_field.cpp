#include "cldera_field.hpp"

#include <ekat/util/ekat_string_utils.hpp>

#include <numeric>
#include <algorithm>

namespace cldera
{

Field::
Field(const std::string& n, const std::vector<int>& d)
 : Field(n,d,1,0)
{
  // We can go ahead and set sizes
  set_part_size(m_part_dim,m_layout[m_part_dim]);
}

Field::
Field(const std::string& n, const std::vector<int>& d, const Real* data)
 : Field(n,d)
{
  set_data(data);
}

Field::
Field (const std::string& n, const std::vector<int>& d,
       const int nparts, const int part_dim)
 : m_name(n)
 , m_layout(d)
{
  EKAT_REQUIRE_MSG (part_dim>=0 && part_dim<m_layout.rank(),
      "Error! Invalid partition dimension.\n"
      "  - Field name: " + m_name + "\n"
      "  - Field rank: " + std::to_string(m_layout.rank()) + "\n"
      "  - Part dim  : " + std::to_string(part_dim) + "\n");
  EKAT_REQUIRE_MSG (nparts>=1 && nparts<=d[part_dim],
      "Error! Invalid number of parts.\n"
      "  - Field name : " + m_name + "\n"
      "  - Num parts  : " + std::to_string(nparts) + "\n"
      "  - Part dim   : " + std::to_string(part_dim) + "\n"
      "  - Field dims : [" + ekat::join(m_layout.dims(),",") + "]\n");

  m_nparts = nparts;
  m_part_dim = part_dim;
  m_data.resize(nparts);
  m_part_sizes.resize(nparts,-1);
}

FieldLayout Field::
part_layout (const int ipart) const {
  EKAT_REQUIRE_MSG (m_committed,
      "Error! Field '" + m_name + "' was not yet committed.\n");
  check_part_idx(ipart);

  std::vector<int> d = m_layout.dims();
  d[m_part_dim] = m_part_sizes[ipart];
  return FieldLayout{d};
}

void Field::
set_part_size (const int ipart, const int part_size)
{
  check_part_idx(ipart);
  EKAT_REQUIRE_MSG (part_size>=0 && part_size<=m_layout[m_part_dim],
      "[Field::set_part_size]\n"
      "  Error! Invalid part size.\n"
      "    - Field name: " + m_name + "\n"
      "    - Part index: " + std::to_string(ipart) + "\n"
      "    - Part size : " + std::to_string(part_size) + "\n");
  EKAT_REQUIRE_MSG (m_part_sizes[ipart]==-1 || m_part_sizes[ipart]==part_size,
      "[Field::set_part_data]\n"
      "  Error! Part size was already set to a different value.\n"
      "    - Field name    : " + m_name + "\n"
      "    - Part index    : " + std::to_string(ipart) + "\n"
      "    - Curr part size: " + std::to_string(m_part_sizes[ipart]) + "\n"
      "    - New part size : " + std::to_string(part_size) + "\n");

  m_part_sizes[ipart] = part_size;
}

void Field::
set_part_size_and_data (const int ipart, const int part_size, const Real* data)
{
  set_part_size(ipart,part_size);
  set_part_data(ipart,data);
}

void Field::
set_data (const Real* data)
{
  EKAT_REQUIRE_MSG (m_nparts==1,
      "Error! Field::set_data is only available for non-partitioned fields.\n");
  set_part_data(0,data);

  // Since there's only one part, we can commit here
  commit();
}

void Field::commit () {
  if (m_committed) {
    // Should we error out instead?
    return;
  }

  // Checks
  int part_sizes_sum = 0;
  for (int ipart=0; ipart<m_nparts; ++ipart) {
    EKAT_REQUIRE_MSG (m_data[ipart].data()!=nullptr,
        "[Field::commit]\n"
        "  Error! Part data was not set for at least one partition.\n"
        "    - Field name: " + m_name + "\n"
        "    - Part index: " + std::to_string(ipart) + "\n");
    for (int jpart=0; jpart<ipart; ++jpart) {
      EKAT_REQUIRE_MSG (m_data[ipart].data()!=m_data[jpart].data(),
        "[Field::commit]\n"
        "  Error! Found two partitions that have the same data pointer.\n"
        "    - Field name: " + m_name + "\n"
        "    - Part 1 index: " + std::to_string(ipart) + "\n"
        "    - Part 2 index: " + std::to_string(jpart) + "\n");
    }
    part_sizes_sum += m_part_sizes[ipart];
  }
  EKAT_REQUIRE_MSG (part_sizes_sum==m_layout[m_part_dim],
      "[Field::commit]\n"
      "  Error! Partition sizes do not add up to layout dimension.\n"
      "    - Field name     : " + m_name + "\n"
      "    - Field dims     : [" + ekat::join(m_layout.dims(),",") + "]\n"
      "    - Part dim       : " + std::to_string(m_part_dim) + "\n"
      "    - Dim[Part dim]  : " + std::to_string(m_layout[m_part_dim]) + "\n"
      "    - Sum(part sizes): " + std::to_string(part_sizes_sum) + "\n");

  m_committed = true;
}

void Field::
check_part_idx (const int ipart) const {
  EKAT_REQUIRE_MSG (ipart>=0 && ipart<m_nparts,
      "Error! Invalid part index.\n"
      "  - Field name: " + m_name + "\n"
      "  - Part index: " + std::to_string(ipart) + "\n"
      "  - Num parts : " + std::to_string(m_nparts) + "\n");
}

void Field::
set_part_data (const int ipart, const Real* data)
{
  EKAT_REQUIRE_MSG (m_data[ipart].data()==nullptr,
      "[Field::set_part_data]\n"
      "  Error! Part data was already set.\n"
      "    - Field name: " + m_name + "\n"
      "    - Part index: " + std::to_string(ipart) + "\n");
  EKAT_REQUIRE_MSG (data!=nullptr,
      "[Field::set_part_data]\n"
      "  Error! Invalid part data pointer.\n"
      "    - Field name: " + m_name + "\n"
      "    - Part index: " + std::to_string(ipart) + "\n");

  m_data[ipart] = hview_1d<const Real> (data,part_layout(ipart).size());
}

} // namespace cldera

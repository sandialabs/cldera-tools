#include "cldera_field.hpp"

#include <ekat/util/ekat_string_utils.hpp>

#include <numeric>
#include <algorithm>

namespace cldera
{

Field::
Field(const std::string& n, const FieldLayout& fl,
      const DataAccess cv, const DataType dt)
 : Field(n,fl,1,0,cv,dt)
{
  // We can go ahead and set sizes
  if (m_layout.rank()>0) {
    set_part_size(m_part_dim,m_layout.extent(m_part_dim));
  } else {
    set_part_size(m_part_dim,1);
  }
}

Field::
Field (const std::string& n, const std::vector<int>& dims,
       const std::vector<std::string>& dimnames,
       const DataAccess cv, const DataType dt)
 : Field(n,FieldLayout(dims,dimnames),cv, dt)
{
  // Nothing to do here
}

Field::
Field (const std::string& n, const FieldLayout& fl,
       const int nparts, const int part_dim,
       const DataAccess cv, const DataType dt)
 : m_name(n)
 , m_layout(fl)
 , m_data_access(cv)
 , m_data_type(dt)
{
  EKAT_REQUIRE_MSG (is_valid(m_data_type),
      "Error! Invalid input data type: " + e2str(dt) + "\n");
  EKAT_REQUIRE_MSG (part_dim>=0 && part_dim<(m_layout.rank()==0 ? 1 : m_layout.rank()),
      "Error! Invalid partition dimension.\n"
      "  - Field name: " + m_name + "\n"
      "  - Field rank: " + std::to_string(m_layout.rank()) + "\n"
      "  - Part dim  : " + std::to_string(part_dim) + "\n");
  EKAT_REQUIRE_MSG (nparts>=1 && nparts<=(m_layout.rank()==0 ? 1 : m_layout.extent(part_dim)),
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

Field::
Field (const std::string& n,
       const std::vector<int>& dims,
       const std::vector<std::string>& dimnames,
       const int nparts, const int part_dim,
       const DataAccess cv, const DataType dt)
 : Field(n,FieldLayout(dims,dimnames),nparts,part_dim,cv,dt)
{
  // Nothing to do here
}

FieldLayout Field::
part_layout (const int ipart) const {
  check_part_idx(ipart);
  EKAT_REQUIRE_MSG (m_part_sizes[ipart]!=-1,
      "[Field::part_layout]\n"
      "  Error! Part size was not yet set.\n"
      "    - Field name: " + m_name + "\n"
      "    - Part index: " + std::to_string(ipart) + "\n");

  if (m_layout.rank()==0) {
    return FieldLayout();
  } else {
    std::vector<int> d = m_layout.dims();
    d[m_part_dim] = m_part_sizes[ipart];
    return FieldLayout{d,m_layout.names()};
  }
}

void Field::
set_part_size (const int ipart, const int part_size)
{
  check_part_idx(ipart);
  EKAT_REQUIRE_MSG (part_size>=0 && part_size<=(m_layout.rank()==0 ? 1 : m_layout.extent(m_part_dim)),
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

void Field::commit () {
  if (m_committed) {
    // Should we error out instead?
    return;
  }

  if (m_data_access==DataAccess::Copy) {
    // Allocate the views
    m_data_nonconst.resize(m_nparts);
    for (int i=0; i<m_nparts; ++i) {
      auto iname = m_name + "_" + std::to_string(i);
      m_data_nonconst[i] = view_1d_host<char> (iname,size_of(m_data_type)*part_layout(i).size());
      m_data[i] = m_data_nonconst[i];
    }
  }

  // Checks
  int part_sizes_sum = 0;
  for (int ipart=0; ipart<m_nparts; ++ipart) {
    EKAT_REQUIRE_MSG (m_data[ipart].data()!=nullptr,
        "[Field::commit]\n"
        "  Error! Part data was not set for at least one partition.\n"
        "    - Field name: " + m_name + "\n"
        "    - Part index: " + std::to_string(ipart) + "\n");

    // Note: the following should be automatically true for Copy fields
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
  EKAT_REQUIRE_MSG (part_sizes_sum==(m_layout.rank()==0 ? 1 : m_layout.extent(m_part_dim)),
      "[Field::commit]\n"
      "  Error! Partition sizes do not add up to layout dimension.\n"
      "    - Field name     : " + m_name + "\n"
      "    - Field dims     : [" + ekat::join(m_layout.dims(),",") + "]\n"
      "    - Part dim       : " + std::to_string(m_part_dim) + "\n"
      "    - Dim[Part dim]  : " + std::to_string(m_layout.extent(m_part_dim)) + "\n"
      "    - Sum(part sizes): " + std::to_string(part_sizes_sum) + "\n");

  m_committed = true;
}

void Field::
check_single_part (const std::string& method_name) const {
  EKAT_REQUIRE_MSG (m_nparts==1,
      "Error! Method Field::" + method_name + " only available for single-part fields.\n");
}

void Field::
check_part_idx (const int ipart) const {
  EKAT_REQUIRE_MSG (ipart>=0 && ipart<m_nparts,
      "Error! Invalid part index.\n"
      "  - Field name: " + m_name + "\n"
      "  - Part index: " + std::to_string(ipart) + "\n"
      "  - Num parts : " + std::to_string(m_nparts) + "\n");
}

} // namespace cldera

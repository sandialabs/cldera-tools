#include "cldera_field.hpp"

#include <ekat/util/ekat_string_utils.hpp>

#include <numeric>
#include <algorithm>

namespace cldera
{

Field::
Field(const std::string& n, const FieldLayout& fl, const DataAccess cv)
 : Field(n,fl,1,0,cv)
{
  // We can go ahead and set sizes
  set_part_size(m_part_dim,m_layout.extent(m_part_dim));
}

Field::
Field (const std::string& n,
       const std::vector<int>& dims,
       const std::vector<std::string>& dimnames,
       const Real* data)
 : Field(n,FieldLayout(dims,dimnames),data)
{
  // Nothing to do here
}

Field::
Field(const std::string& n, const FieldLayout& fl, const Real* data)
 : Field(n,fl,DataAccess::View)
{
  set_data(data);
}

Field::
Field (const std::string& n,
       const std::vector<int>& dims,
       const std::vector<std::string>& dimnames,
       const DataAccess cv)
 : Field(n,FieldLayout(dims,dimnames),cv)
{
  // Nothing to do here
}

Field::
Field (const std::string& n, const FieldLayout& fl,
       const int nparts, const int part_dim, const DataAccess cv)
 : m_name(n)
 , m_layout(fl)
 , m_data_access(cv)
{
  EKAT_REQUIRE_MSG (part_dim>=0 && part_dim<m_layout.rank(),
      "Error! Invalid partition dimension.\n"
      "  - Field name: " + m_name + "\n"
      "  - Field rank: " + std::to_string(m_layout.rank()) + "\n"
      "  - Part dim  : " + std::to_string(part_dim) + "\n");
  EKAT_REQUIRE_MSG (nparts>=1 && nparts<=m_layout.extent(part_dim),
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
       const DataAccess cv)
 : Field(n,FieldLayout(dims,dimnames),nparts,part_dim,cv)
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

  std::vector<int> d = m_layout.dims();
  d[m_part_dim] = m_part_sizes[ipart];
  return FieldLayout{d,m_layout.names()};
}

void Field::
set_part_size (const int ipart, const int part_size)
{
  check_part_idx(ipart);
  EKAT_REQUIRE_MSG (part_size>=0 && part_size<=m_layout.extent(m_part_dim),
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
copy_part_data (const int ipart, const Real* data)
{
  EKAT_REQUIRE_MSG (m_data_access==DataAccess::Copy,
      "[Field::copy_part_data]\n"
      "  Error! Attempt to copy data, but field data access is not 'Copy'.\n"
      "    - Field name: " + m_name + "\n");

  check_part_idx(ipart);
  EKAT_REQUIRE_MSG (data!=nullptr,
      "[Field::copy_part_data]\n"
      "  Error! Invalid part data pointer.\n"
      "    - Field name: " + m_name + "\n"
      "    - Part index: " + std::to_string(ipart) + "\n");

  view_1d_host<const Real,Kokkos::MemoryUnmanaged> v(data,part_layout(ipart).size());
  Kokkos::deep_copy(m_data_nonconst[ipart],v);
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

void Field::
copy_data (const Real* data)
{
  copy_part_data(0,data);
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
      m_data_nonconst[i] = view_1d_host<Real> (iname,part_layout(i).size());
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
  EKAT_REQUIRE_MSG (part_sizes_sum==m_layout.extent(m_part_dim),
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
  EKAT_REQUIRE_MSG (m_data_access==DataAccess::View,
      "[Field::set_part_data]\n"
      "  Error! Attempt to set data pointer, but data access is 'Copy'.\n"
      "    - Field name: " + m_name + "\n");

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

  m_data[ipart] = view_1d_host<const Real> (data,part_layout(ipart).size());
}

} // namespace cldera

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
  // We can go ahead and set the part extent
  if (m_layout.rank()>0) {
    set_part_extent(m_part_dim,m_layout.extent(m_part_dim));
  } else {
    set_part_extent(m_part_dim,1);
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
       const DataAccess cv, const DataType dt,
       int part_dim_alloc_size)
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

  // EKAT_REQUIRE_MSG (cv==DataAccess::View || part_dim_alloc_size==-1,
  //   "Error! If DataAccess==Copy, do not attempt to prescribe part_dim_alloc_size.\n"
  //   "  - Field name         : " + m_name + "\n"
  //   "  - Part dim           : " + std::to_string(part_dim) + "\n"
  //   "  - Field dims         : [" + ekat::join(m_layout.dims(),",") + "]\n"
  //   "  - part dim alloc size: [" + std::to_string(part_dim_alloc_size) + "]\n");

  // If part_dim_alloc_size==-1, it means alloc_size = part_extent
  m_part_dim_alloc_size = part_dim_alloc_size;

  m_nparts = nparts;
  m_part_dim = part_dim;
  m_data.resize(nparts);
  m_part_extents.resize(nparts,-1);
}

Field::
Field (const std::string& n,
       const std::vector<int>& dims,
       const std::vector<std::string>& dimnames,
       const int nparts, const int part_dim,
       const DataAccess cv, const DataType dt,
       const int part_dim_alloc_size)
 : Field(n,FieldLayout(dims,dimnames),nparts,part_dim,cv,dt,part_dim_alloc_size)
{
  // Nothing to do here
}

FieldLayout Field::
part_layout (const int ipart) const {
  check_part_idx(ipart);
  EKAT_REQUIRE_MSG (m_part_extents[ipart]!=-1,
      "[Field::part_layout]\n"
      "  Error! Part extent was not yet set.\n"
      "    - Field name: " + m_name + "\n"
      "    - Part index: " + std::to_string(ipart) + "\n");

  if (m_layout.rank()==0) {
    return FieldLayout();
  } else {
    std::vector<int> dims = m_layout.dims();
    dims[m_part_dim] = m_part_extents[ipart];
    if (m_part_dim_alloc_size==-1) {
      return FieldLayout(dims,m_layout.names());
    } else {
      std::vector<int> alloc_dims = m_layout.dims();
      alloc_dims[m_part_dim] = m_part_dim_alloc_size;
      return FieldLayout(dims,alloc_dims,m_layout.names());
    }
  }
}

void Field::
set_part_extent (const int ipart, const int part_extent)
{
  check_part_idx(ipart);
  EKAT_REQUIRE_MSG (part_extent>=0 && part_extent<=(m_layout.rank()==0 ? 1 : m_layout.extent(m_part_dim)),
      "[Field::set_part_extent]\n"
      "  Error! Invalid part extent.\n"
      "    - Field name: " + m_name + "\n"
      "    - Part index: " + std::to_string(ipart) + "\n"
      "    - Part extent : " + std::to_string(part_extent) + "\n");
  EKAT_REQUIRE_MSG (m_part_extents[ipart]==-1 || m_part_extents[ipart]==part_extent,
      "[Field::set_part_data]\n"
      "  Error! Part extent was already set to a different value.\n"
      "    - Field name    : " + m_name + "\n"
      "    - Part index    : " + std::to_string(ipart) + "\n"
      "    - Curr part extent: " + std::to_string(m_part_extents[ipart]) + "\n"
      "    - New part extent : " + std::to_string(part_extent) + "\n");
  EKAT_REQUIRE_MSG (m_part_dim_alloc_size==-1 || m_part_dim_alloc_size>=part_extent,
      "[Field::set_part_extent]\n"
      "  Error! Part extent exceeds alloc size along partitioned dim.\n"
      "    - Field name     : " + m_name + "\n"
      "    - Part index     : " + std::to_string(ipart) + "\n"
      "    - Part extent    : " + std::to_string(part_extent) + "\n"
      "    - Part alloc size: " + std::to_string(m_part_dim_alloc_size) + "\n");

  m_part_extents[ipart] = part_extent;
}

void Field::commit () {
  if (m_committed) {
    // It's better to allow this, in case we commit one field, and later
    // call commit_all_fields in the profiling interface.
    return;
  }

  if (m_data_access==DataAccess::Copy) {
    // Allocate the views
    m_data_nonconst.resize(m_nparts);
    for (int i=0; i<m_nparts; ++i) {
      auto pl = part_layout(i);
      auto iname = m_name + "_" + std::to_string(i);
      m_data_nonconst[i] = view_1d_host<char> (iname,size_of(m_data_type)*pl.alloc_size());
      m_data[i] = m_data_nonconst[i];
    }
  }

  // Checks
  int part_extents_sum = 0;
  m_part_offsets.resize(m_nparts);
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
    m_part_offsets[ipart] = part_extents_sum;
    part_extents_sum += m_part_extents[ipart];
  }
  EKAT_REQUIRE_MSG (part_extents_sum==(m_layout.rank()==0 ? 1 : m_layout.extent(m_part_dim)),
      "[Field::commit]\n"
      "  Error! Partition extents do not add up to layout dimension.\n"
      "    - Field name       : " + m_name + "\n"
      "    - Field dims       : [" + ekat::join(m_layout.dims(),",") + "]\n"
      "    - Part dim         : " + std::to_string(m_part_dim) + "\n"
      "    - Dim[Part dim]    : " + std::to_string(m_layout.extent(m_part_dim)) + "\n"
      "    - Sum(part extents): " + std::to_string(part_extents_sum) + "\n");

  m_committed = true;
}

Field Field::clone() const {
  check_committed(true,"Field::clone");
  Field f(m_name,m_layout,m_nparts,m_part_dim,DataAccess::Copy,m_data_type);

  // Recall, public/private is based on type, not instance, so we can
  // access f's members
  f.m_part_dim_alloc_size = m_part_dim_alloc_size;
  f.m_part_extents = m_part_extents;
  f.m_committed = true;
  f.m_data_nonconst.resize(m_nparts);
  f.m_data.resize(m_nparts);

  for (int i=0; i<m_nparts; ++i) {
    auto iname = m_name + "_" + std::to_string(i);
    f.m_data_nonconst[i] = view_1d_host<char> (iname,size_of(m_data_type)*part_layout(i).size());
    f.m_data[i] = f.m_data_nonconst[i];
    Kokkos::deep_copy(f.m_data_nonconst[i],m_data[i]);
  }

  return f;
}

Field Field::read_only() const {
  Field f (*this);
  // Ensure copy_part_data cannot be called on this clone
  f.m_data_access = DataAccess::View;
  return f;
}

void Field::deep_copy (const Field& src) {
  check_committed(true,"Field::deep_copy");
  EKAT_REQUIRE_MSG (src.m_committed,
      "Error! Cannot deep copy if src field is not yet committed.\n"
      " Source field: " + src.m_name + "\n"
      " Target field: " + m_name + "\n");
  EKAT_REQUIRE_MSG (m_layout==src.m_layout,
      "Error! Cannot deep copy from a field with different layout.\n"
      " Source field: " + src.m_name + "\n"
      " Target field: " + m_name + "\n");
  EKAT_REQUIRE_MSG (m_nparts==src.m_nparts,
      "Error! Cannot deep copy from a field with different number of parts.\n"
      " Source field: " + src.m_name + "\n"
      " Target field: " + m_name + "\n");
  for (int p=0; p<m_nparts; ++p) {
    EKAT_REQUIRE_MSG (part_layout(p)==src.part_layout(p),
        "Error! Cannot deep copy from a field with different part layouts.\n"
        " Source field: " + src.m_name + "\n"
        " Target field: " + m_name + "\n"
        " Incompatible part index: " + std::to_string(p) + "\n"
        " Source part layout: " + src.part_layout(p).to_string() + "\n"
        " Target part layout: " + part_layout(p).to_string() + "\n"
        );
    auto src_pv = src.m_data[p];
    auto tgt_pv = m_data_nonconst[p];
    Kokkos::deep_copy(tgt_pv,src_pv);
  }
}

int Field::
part_offset (const int ipart) const {
  check_part_idx(ipart);
  return m_part_offsets[ipart];
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

void Field::
check_committed (const bool expected, const std::string& context) const {
  EKAT_REQUIRE_MSG (m_committed==expected,
      "Error! Field was not in the expected committed state.\n"
      "  - Field Name: " + m_name + "\n"
      "  - Expected  : " + std::string(expected ? "committed" : "not committed") + "\n"
      "  - Actual    : " + std::string(m_committed ? "committed" : "not committed") + "\n"
      "  - Context   : " + context + "\n");
}

void Field::
check_rank (const int N) const {
  EKAT_REQUIRE_MSG (N==layout().rank(),
      "Error! Wrong value for field rank.\n"
      " - field name    : " + m_name + "\n"
      " - field rank    : " + std::to_string(layout().rank()) + "\n"
      " - requested rank: " + std::to_string(N) + "\n");
}

} // namespace cldera

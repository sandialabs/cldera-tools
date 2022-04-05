#include "cldera_field.hpp"

#include <numeric>
#include <algorithm>

namespace cldera
{

Field::
Field(const std::string& n, const std::vector<int>& d)
 : Field(n,d,1,0)
{
  // Nothing to do
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
      "  - Field name   : " + m_name + "\n"
      "  - Num parts    : " + std::to_string(nparts) + "\n"
      "  - Part dim     : " + std::to_string(part_dim) + "\n"
      "  - Dim[Part dim]: " + std::to_string(d[part_dim]) + "\n");

  m_nparts = nparts;
  m_part_dim = part_dim;
  m_data.resize(nparts,nullptr);
  m_part_beg.resize(nparts+1,-1);
}

void Field::
set_part_data (const int i, const int beg, const Real* data)
{
  check_part_idx(i);
  EKAT_REQUIRE_MSG (m_data[i]==nullptr,
      "Error! Part was already set.\n"
      "  - Field name: " + m_name + "\n"
      "  - Part index: " + std::to_string(i) + "\n");

  m_data[i] = data;
  m_part_beg[i] = beg;
}

void Field::
set_data (const Real* data)
{
  EKAT_REQUIRE_MSG (m_nparts==1,
      "Error! Field::set_data is only available for non-partitioned fields.\n");
  set_part_data(0,0,data);

  // Since there's only one part, we can commit here
  commit();
}

void Field::commit () {
  if (m_committed) {
    // Should we error out instead?
    return;
  }

  // First, rearrange partitions, to ensure they are ordered along partition_dim
  // We need this to be able to use consecutive entries in m_part_beg as bounds
  // This trick allows to get the permutation that orders partitions
  std::vector<int> idx(m_nparts);
  std::iota(idx.begin(),idx.end(),0);
  std::sort(idx.begin(),idx.end(),
      [&](const int& lhs, const int& rhs) -> bool {
          return m_part_beg[lhs]<m_part_beg[rhs];
  });

  auto data = m_data;
  auto pbeg = m_part_beg;
  for (int i=0; i<m_nparts; ++i) {
    m_part_beg[i] = pbeg[idx[i]];
    m_data[i]     = data[idx[i]];
  }
  m_part_beg[m_nparts] = m_layout.dims()[m_part_dim];

  // Checks
  auto it1 = std::adjacent_find(m_part_beg.begin(),m_part_beg.end());
  EKAT_REQUIRE_MSG (it1==m_part_beg.end(),
      "[Field::commit]\n"
      "  Error! Found two partitions that overlap.\n");

  auto tos = [](Real* p) -> std::string {
    std::stringstream ss;
    ss << (void*)p;
    return ss.str();
  };

  std::sort(data.begin(),data.end());
  auto it2 = std::adjacent_find(data.begin(),data.end());
  EKAT_REQUIRE_MSG (it2==data.end(),
      "[ProfilingArchive::commit_partitioned_field]\n"
      "  Error! Found two partitions that have the same data pointer.\n");

  m_committed = true;
}

void Field::
check_part_idx (const int i) const {
  EKAT_REQUIRE_MSG (i>=0 && i<m_nparts,
      "Error! Invalid part index.\n"
      "  - Field name: " + m_name + "\n"
      "  - Part index: " + std::to_string(i) + "\n"
      "  - Num parts : " + std::to_string(m_nparts) + "\n");
}


} // namespace cldera

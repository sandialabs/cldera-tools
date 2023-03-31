#include "cldera_field_layout.hpp"

#include <ekat/util/ekat_string_utils.hpp>
#include <ekat/ekat_assert.hpp>

namespace cldera
{

FieldLayout::
FieldLayout (const std::vector<int>& dims,
             const std::vector<int>& alloc_dims,
             const std::vector<std::string>& names)
{
  for (auto d : dims) {
    EKAT_REQUIRE_MSG (d>0, "Error! Invalid dimension (" + std::to_string(d) + "\n");
  }
  m_dims = dims;

  EKAT_REQUIRE_MSG (names.size()==m_dims.size(),
      "Error! Size of names and dims array must match.\n");
  m_names = names;
  for (int i=0; i<rank(); ++i) {
    m_kokkos_layout.dimension[i] = alloc_dims[i];
  }
}

long long FieldLayout::
size () const {
  long long s = 1;
  for (auto d : dims()) {
    s *= d;
  }
  return s;
}

int FieldLayout::
dim_idx (const std::string& name) const {
  EKAT_REQUIRE_MSG (has_dim_name (name),
      "Error! Cannot get dimension index in layout: dimension not found.\n"
      " - stored names: " + ekat::join(m_names,",") + "\n"
      " - input name  : " + name + "\n");

  return std::distance(m_names.begin(),std::find(m_names.begin(),m_names.end(),name));
}

bool FieldLayout::
has_dim_name (const std::string& name) const {
  return std::find(m_names.begin(),m_names.end(),name)!=m_names.end();
}

FieldLayout FieldLayout::
strip_dim (const std::string& name) const {
  auto d = dims();
  auto n = names();
  auto pos = dim_idx(name);

  d.erase(d.begin()+pos);
  n.erase(n.begin()+pos);
  return FieldLayout(d,n);
}

int FieldLayout::
extent (const std::string& name) {
  auto it = std::find(m_names.begin(),m_names.end(),name);
  EKAT_REQUIRE_MSG (it!=m_names.end(),
      "Error! Dimension name not found.\n"
      " - stored names: " + ekat::join(m_names,",") + "\n"
      " - input name  : " + name + "\n");

  EKAT_REQUIRE_MSG (std::find(std::next(it),m_names.end(),name)==m_names.end(),
      "Error! Muyltiple dimensions have the requested name.\n"
      " - stored names: " + ekat::join(m_names,",") + "\n"
      " - input name  : " + name + "\n");
  return m_dims[std::distance(m_names.begin(),it)];
}

std::string FieldLayout::to_string () const {
  std::string names = "<" + ekat::join(m_names,",") + ">";
  std::string dims  = "(" + ekat::join(m_dims,",") + ")";

  return names + " " + dims;
}

} // namespace cldera

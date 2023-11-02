#ifndef CLDERA_FIELD_LAYOUT_HPP
#define CLDERA_FIELD_LAYOUT_HPP

#include "cldera_profiling_types.hpp"

#include <ekat/util/ekat_string_utils.hpp>
#include <ekat/ekat_assert.hpp>

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
  // It's the layout of a scalar wrapped in a Field
  FieldLayout () = default;

  FieldLayout (const std::vector<int>& dims,
               const std::vector<std::string>& names)
    : FieldLayout (dims,dims,names)
  {
    // Nothing to do here
  }

  FieldLayout (const std::vector<int>& dims,
               const std::vector<int>& alloc_dims,
               const std::vector<std::string>& names);

  const std::vector<int>& dims () const { return m_dims; }
  const std::vector<std::string>& names () const { return m_names; }

  bool has_dim (const std::string& name) const {
    auto it = std::find(m_names.begin(),m_names.end(),name);
    return it!=m_names.end();
  }
  int idim (const std::string& name) const {
    auto it = std::find(m_names.begin(),m_names.end(),name);
    EKAT_REQUIRE_MSG(it!=m_names.end(),
        "Error! Input dim name not found in this layout.\n"
        " - dim names : [" + ekat::join(m_names,",") + "]\n"
        " - input name: " + name + "\n");
    return std::distance(m_names.begin(),it);
  }
  int extent (const std::string& name) const {
    return m_dims[idim(name)];
  }

  int extent (const int i) const { return m_dims[i]; }
  std::string name (const int i) const { return m_names[i]; }

  int dim_idx (const std::string& name) const;
  int rank () const { return dims().size(); }
  long long size () const;

  bool has_dim_name (const std::string& name) const;

  FieldLayout strip_dim (const std::string& name) const;

  int extent (const std::string& name);

  std::string to_string () const;

  long long alloc_size () const {
    long long s = 1;
    for (int i=0; i<rank(); ++i) {
      s *= m_kokkos_layout.dimension[i];
    }
    return s;
  }

  const Kokkos::LayoutRight& kokkos_layout () const { return m_kokkos_layout; }

  friend bool operator== (const FieldLayout& lhs, const FieldLayout& rhs);
private:
  std::vector<int>          m_dims;
  std::vector<std::string>  m_names;
  Kokkos::LayoutRight       m_kokkos_layout;
};

inline bool operator== (const FieldLayout& lhs, const FieldLayout& rhs) {
  return lhs.m_dims  == rhs.m_dims &&
         lhs.m_kokkos_layout == rhs.m_kokkos_layout &&
         lhs.m_names == rhs.m_names;
}

inline bool operator!= (const FieldLayout& lhs, const FieldLayout& rhs) {
  return not (lhs==rhs);
}

} // namespace cldera

#endif // CLDERA_FIELD_LAYOUT_HPP

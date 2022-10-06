#ifndef CLDERA_FIELD_LAYOUT_HPP
#define CLDERA_FIELD_LAYOUT_HPP

#include "cldera_profiling_types.hpp"

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
  {
    for (auto d : dims) {
      EKAT_REQUIRE_MSG (d>0, "Error! Invalid dimension (" + std::to_string(d) + "\n");
    }
    m_dims = dims;

    EKAT_REQUIRE_MSG (names.size()==m_dims.size(),
        "Error! Size of names and dims array must match.\n");
    m_names = names;
  }

  const std::vector<int>& dims () const { return m_dims; }
  const std::vector<std::string>& names () const { return m_names; }

  int extent (const int i) const { return m_dims[i]; }
  std::string name (const int i) const { return m_names[i]; }

  int rank () const { return dims().size(); }
  long long size () const {
    long long s = 1;
    for (auto d : dims()) {
      s *= d;
    }
    return s;
  }

  friend bool operator== (const FieldLayout& lhs, const FieldLayout& rhs);
private:
  std::vector<int>          m_dims;
  std::vector<std::string>  m_names;
};

inline bool operator== (const FieldLayout& lhs, const FieldLayout& rhs) {
  return lhs.m_dims  == rhs.m_dims &&
         lhs.m_names == rhs.m_names;
}

inline bool operator!= (const FieldLayout& lhs, const FieldLayout& rhs) {
  return not (lhs==rhs);
}

} // namespace cldera

#endif // CLDERA_FIELD_LAYOUT_HPP

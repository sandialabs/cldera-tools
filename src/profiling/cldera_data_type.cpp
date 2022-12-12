#include "cldera_data_type.hpp"
#include "cldera_profiling_types.hpp"

namespace cldera {

std::string e2str (const DataType dt) {
  std::string s;
  switch (dt) {
    case RealType:
      s = "real"; break;
    case IntType:
      s = "int";  break;
    default:
      s = "invalid";
  }
  return s;
}

bool is_valid (const DataType dt) {
  return dt==IntType || dt==RealType;
}

template<>
DataType get_data_type<int> () {
  return IntType;
}

template<>
DataType get_data_type<Real> () {
  return RealType;
}

size_t size_of (const DataType dt) {
  switch (dt) {
    case RealType:
      return sizeof(Real);
    case IntType:
      return sizeof(int);
    default:
      return 0;
  }
}

} // namespace cldera

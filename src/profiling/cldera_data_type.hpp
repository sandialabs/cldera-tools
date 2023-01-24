#ifndef CLDERA_DATA_TYPE_HPP
#define CLDERA_DATA_TYPE_HPP

#include <string>

namespace cldera {

enum DataType {
  Invalid  = 0,
  RealType = 1,
  IntType  = 2
};

std::string e2str (const DataType dt);

DataType str2data_type (const std::string& dt);

bool is_valid (const DataType dt);

template<typename T>
typename std::enable_if<not std::is_const<T>::value,DataType>::type
get_data_type ();
template<typename T>
typename std::enable_if<std::is_const<T>::value,DataType>::type
get_data_type () {
  return get_data_type<typename std::remove_const<T>::type>();
}

size_t size_of (const DataType dt);

} // namespace cldera

#endif // CLDERA_DATA_TYPE_HPP

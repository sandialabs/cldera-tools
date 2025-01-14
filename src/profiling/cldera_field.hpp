#ifndef CLDERA_FIELD_HPP
#define CLDERA_FIELD_HPP

#include "cldera_field_layout.hpp"
#include "cldera_data_type.hpp"

#include <ekat/ekat_assert.hpp>

#include <vector>
#include <string>

namespace cldera
{

/*
 * A small struct holding data and metadata of a field
 *
 * The basic metadata is given by the 'name' and 'dims' members.
 *
 * The Field is assumed to be potentially partitioned in memory.
 * The partitioning is described by
 *  - nparts: number of partitions
 *  - part_dim: index of dimension along which field is partitioned
 *  - part_extent: the extent along part_dim of each partition
 *
 * Notice that each partition is stored contiguously in memory.
 * The method part_layout allows to get the layout of the partition,
 * so that one can safely iterate over the partition.
 *
 * Note: a contiguous/non-partitioned field simply has nparts=1.
 */

enum class DataAccess {
  Copy,
  View
};

class Field {
public:
  // Not really needed for actual work, but it allows e.g. to use op[] on std::map
  Field () = default;

  Field (const std::string& n, const FieldLayout& fl,
         const int nparts, const int part_dim,
         const DataAccess cv = DataAccess::View,
         const DataType dt = RealType,
         const int part_dim_alloc_size = -1);
  Field (const std::string& n,
         const std::vector<int>& dims,
         const std::vector<std::string>& dimnames,
         const int nparts, const int part_dim,
         const DataAccess cv = DataAccess::View,
         const DataType dt = RealType,
         const int part_dim_alloc_size = -1);

  // Shortcuts for single-partition field
  Field (const std::string& n, const FieldLayout& fl,
         const DataAccess cv = DataAccess::View,
         const DataType dt = RealType);
  Field (const std::string& n,
         const std::vector<int>& dims,
         const std::vector<std::string>& dimnames,
         const DataAccess cv = DataAccess::View,
         const DataType dt = RealType);

  template<typename T>
  Field (const std::string& n, const FieldLayout& fl, const T* data);

  template<typename T>
  Field (const std::string& n,
         const std::vector<int>& dims,
         const std::vector<std::string>& dimnames,
         const T* data);

  Field& operator= (const Field& src) = default;

  const std::string& name () const { return m_name; }

  // Get rank-global (not partitioned) and part layouts
  const FieldLayout& layout () const { return m_layout; }
        FieldLayout part_layout (const int ipart) const;

  // Set part specs
  void set_part_extent  (const int ipart, const int part_extent);
  template<typename T>
  void set_part_data  (const int ipart, const T* data);
  template<typename T>
  void set_data (const T* data);

  void commit ();

  void rename (const std::string& name) { m_name = name; }

  // Copy into managed views, if m_data_access=Copy
  template<typename T>
  void copy_part_data (const int ipart, const T* data);
  template<typename T>
  void copy_data (const T* data);

  // Get data as view
  template<typename T>
  view_1d_host<const T> part_view (const int ipart) const;
  template<typename T>
  view_1d_host<      T> part_view_nonconst (const int ipart);
  template<typename T>
  view_1d_host<const T> view () const;
  template<typename T>
  view_1d_host<      T> view_nonconst ();

  template<typename T,int N>
  view_Nd_host<const T,N> part_nd_view (const int ipart) const;
  template<typename T,int N>
  view_Nd_host<      T,N> part_nd_view_nonconst (const int ipart);
  template<typename T,int N>
  view_Nd_host<const T,N> nd_view () const;
  template<typename T,int N>
  view_Nd_host<      T,N> nd_view_nonconst ();

  // Get raw data
  template<typename T>
  const T* part_data (const int ipart) const;
  template<typename T>
        T* part_data_nonconst (const int ipart);
  template<typename T>
  const T* data () const;
  template<typename T>
        T* data_nonconst ();

  // Query status
  int nparts () const { return m_nparts; }
  int part_dim () const { return m_part_dim; }
  int part_offset (const int ipart) const;
  bool committed () const { return m_committed; }
  DataAccess data_access () const { return m_data_access; }
  DataType data_type () const { return m_data_type; }

  Field clone () const;
  Field read_only () const;

  // Perform this = beta*this + alpha*x
  template<typename T>
  void update (const Field& x, const T alpha, const T beta);

  template<typename T>
  void scale (const T alpha);

  template<typename T>
  void deep_copy (const T val);

  void deep_copy (const Field& src);
private:

  // Methods to go to and from the internal char* storage
  template<typename T>
  char* ptr2char (T* p) const { return reinterpret_cast<char*>(p); }
  template<typename T>
  const char* ptr2char (const T* p) const { return reinterpret_cast<const char*>(p); }

  template<typename T>
  T* char2ptr (char* p) const { return reinterpret_cast<T*>(p); }
  template<typename T>
  const T* char2ptr (const char* p) const { return reinterpret_cast<const T*>(p); }

  template<typename T>
  void update_impl (const Field& x, const T alpha, const T beta);
  template<typename T>
  void scale_impl (const T alpha);

  // Check methods are called with the right inputs for this field
  template<typename T>
  void check_data_type () const;
  void check_single_part (const std::string& method_name) const;
  void check_part_idx (const int i) const;
  void check_committed (const bool expected, const std::string& context) const;
  void check_rank (const int N) const;

  std::string       m_name;
  FieldLayout       m_layout;
  int               m_nparts   = -1;  // Set to something invalid for default ctor
  int               m_part_dim = -1;
  std::vector<int>  m_part_extents;
  std::vector<int>  m_part_offsets;
  int               m_part_dim_alloc_size;
  bool              m_committed = false;
  DataAccess        m_data_access;
  DataType          m_data_type;

  // Store data as char
  std::vector<view_1d_host<const char>>   m_data;
  std::vector<view_1d_host<      char>>   m_data_nonconst;
};

// =================== IMPLEMENTATION =================== //

template<typename T>
Field::
Field (const std::string& n,
       const std::vector<int>& dims,
       const std::vector<std::string>& dimnames,
       const T* data)
 : Field(n,FieldLayout(dims,dimnames),data)
{
  // Nothing to do here
}

template<typename T>
Field::
Field(const std::string& n, const FieldLayout& fl, const T* data)
 : Field(n,fl,DataAccess::View)
{
  set_data(data);
}

template<typename T>
void Field::
set_data (const T* data)
{
  check_single_part("set_data");
  set_part_data(0,data);

  // Since there's only one part, we can commit here
  commit();
}

template<typename T>
void Field::
set_part_data (const int ipart, const T* data)
{
  check_committed(false,"Field::set_part_data");
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
  check_data_type<T>();

  const auto alloc_size = size_of(m_data_type)*part_layout(ipart).alloc_size();
  m_data[ipart] = view_1d_host<const char> (ptr2char(data),alloc_size);
}

template<typename T>
void Field::
copy_data (const T* data)
{
  check_single_part("copy_data");
  copy_part_data(0,data);
}

template<typename T>
void Field::
copy_part_data (const int ipart, const T* data)
{
  check_committed(true,"Field::copy_part_data");
  EKAT_REQUIRE_MSG (m_data_access==DataAccess::Copy,
      "[Field::copy_part_data]\n"
      "  Error! Attempt to copy data, but field data access is not 'Copy'.\n"
      "    - Field name: " + m_name + "\n");
  EKAT_REQUIRE_MSG (data!=nullptr,
      "[Field::copy_part_data]\n"
      "  Error! Invalid part data pointer.\n"
      "    - Field name: " + m_name + "\n"
      "    - Part index: " + std::to_string(ipart) + "\n");

  check_data_type<T>();
  check_part_idx(ipart);

  view_1d_host<const T,Kokkos::MemoryUnmanaged> v(data,part_layout(ipart).alloc_size());
  Kokkos::deep_copy(part_view_nonconst<T>(ipart),v);
}

template<typename T>
const T*
Field::part_data (const int ipart) const {
  return part_view<const T>(ipart).data();
}

template<typename T>
T*
Field::part_data_nonconst (const int ipart) {
  return part_view_nonconst<T>(ipart).data();
}

template<typename T>
const T*
Field::data () const {
  return view<const T>().data();
}

template<typename T>
T*
Field::data_nonconst () {
  return view_nonconst<T>().data();
}

template<typename T>
view_1d_host<const T>
Field::part_view (const int ipart) const
{
  check_committed(true,"Field::part_view");
  check_data_type<T>();
  check_part_idx(ipart);

  const auto data = char2ptr<const T>(m_data[ipart].data());
  const auto size = part_layout(ipart).alloc_size();
  return view_1d_host<const T>(data,size);
}

template<typename T>
view_1d_host<T>
Field::part_view_nonconst (const int ipart)
{
  check_committed(true,"Field::part_view_nonconst");
  check_data_type<T>();
  check_part_idx(ipart);

  const auto data = char2ptr<T>(m_data_nonconst[ipart].data());
  const auto size = part_layout(ipart).alloc_size();
  return view_1d_host<T>(data,size);
}

template<typename T>
view_1d_host<const T>
Field::view () const
{
  check_single_part("view");
  return part_view<const T>(0);
}

template<typename T>
view_1d_host<T>
Field::view_nonconst ()
{
  check_single_part("view_nonconst");
  return part_view_nonconst<T>(0);
}

template<typename T, int N>
view_Nd_host<const T,N>
Field::part_nd_view (const int ipart) const
{
  check_committed(true,"Field::part_nd_view");
  check_rank(N);
  check_data_type<T>();
  check_part_idx(ipart);

  const auto data = char2ptr<const T>(m_data[ipart].data());
  return view_Nd_host<const T,N>(data,part_layout(ipart).kokkos_layout());
}

template<typename T, int N>
view_Nd_host<T,N>
Field::part_nd_view_nonconst (const int ipart)
{
  check_committed(true,"Field::part_nd_view_nonconst");
  check_rank(N);
  check_data_type<T>();
  check_part_idx(ipart);

  const auto data = char2ptr<T>(m_data_nonconst[ipart].data());
  return view_Nd_host<T,N>(data,part_layout(ipart).kokkos_layout());
}

template<typename T, int N>
view_Nd_host<const T,N>
Field::nd_view () const
{
  check_single_part("nd_view");
  return part_nd_view<const T,N>(0);
}

template<typename T, int N>
view_Nd_host<T,N>
Field::nd_view_nonconst ()
{
  check_single_part("nd_view_nonconst");
  return part_nd_view_nonconst<T,N>(0);
}

template<typename T>
void Field::update (const Field& x, const T alpha, const T beta)
{
  check_committed(true,"Field::update");
  x.check_committed(true,"Field::update");
  EKAT_REQUIRE_MSG (layout()==x.layout(),
      "Error! Field::update called with fields of different layout."
      "  - tgt field name: " + name() + "\n"
      "  - src field name: " + x.name() + "\n"
      "  - tgt field layout: " + layout().to_string() + "\n"
      "  - src field layout: " + x.layout().to_string() + "\n");
  EKAT_REQUIRE_MSG (m_data_type==DataType::RealType or get_data_type<T>()==DataType::IntType,
      "Error! Field::update called with bad coefficients type.\n"
      "  - field name: " + name() + "\n"
      "  - field data type: " + e2str(m_data_type) + "\n"
      "  - coeff data type: " + e2str(get_data_type<T>()) + "\n");

  EKAT_REQUIRE_MSG (x.nparts()==nparts(),
      "Error! Field::update called with fields of different number of parts."
      "  - tgt field name: " + name() + "\n"
      "  - src field name: " + x.name() + "\n"
      "  - tgt field nparts: " + std::to_string(nparts()) + "\n"
      "  - src field nparts: " + std::to_string(x.nparts()) + "\n");

  if (m_data_type==DataType::IntType) {
    update_impl<int>(x,static_cast<int>(alpha),static_cast<int>(beta));
  } else if (m_data_type==DataType::RealType) {
    update_impl<Real>(x,static_cast<Real>(alpha),static_cast<Real>(beta));
  } else {
    EKAT_ERROR_MSG (
        "[Field::update] Error! Invalid/unsupported data type.\n"
        "  - field name: " + name() + "\n"
        "  - field data type: " + e2str(m_data_type) + "\n");
  }
}

template<typename T>
void Field::scale (const T alpha)
{
  check_committed(true,"Field::scale");

  EKAT_REQUIRE_MSG (m_data_type==DataType::RealType or get_data_type<T>()==DataType::IntType,
      "Error! Field::scale called with bad coefficient type.\n"
      "  - field name: " + name() + "\n"
      "  - field data type: " + e2str(m_data_type) + "\n"
      "  - coeff data type: " + e2str(get_data_type<T>()) + "\n");

  if (m_data_type==DataType::IntType) {
    scale_impl<int>(static_cast<int>(alpha));
  } else if (m_data_type==DataType::RealType) {
    scale_impl<Real>(static_cast<Real>(alpha));
  } else {
    EKAT_ERROR_MSG (
        "[Field::scale] Error! Invalid/unsupported data type.\n"
        "  - field name: " + name() + "\n"
        "  - field data type: " + e2str(m_data_type) + "\n");
  }
}

template<typename T>
void Field::
deep_copy (const T val)
{
  auto run = [&](const auto value) {
    using value_t = typename std::remove_cv<decltype(value)>::type;
    for (int p=0; p<m_nparts; ++p) {
      auto pv = part_view_nonconst<value_t>(p);
      Kokkos::deep_copy(pv,value);
    }
  };

  if (m_data_type==DataType::IntType) {
    run(static_cast<int>(val));
  } else if (m_data_type==DataType::RealType) {
    run(static_cast<Real>(val));
  } else {
    EKAT_ERROR_MSG (
        "[Field::update] Error! Invalid/unsupported data type.\n"
        "  - field name: " + name() + "\n"
        "  - field data type: " + e2str(m_data_type) + "\n");
  }
}

template<typename T>
void Field::update_impl (const Field& x, const T alpha, const T beta)
{
  for (int p=0; p<m_nparts; ++p) {
    auto yv = part_view_nonconst<T>(p);
    auto xv = x.part_view<T>(p);
    auto pl = part_layout(p);
    for (int i=0; i<pl.size(); ++i) {
      yv(i) = beta*yv(i) + alpha*xv(i);
    }
  }
}

template<typename T>
void Field::scale_impl (const T alpha)
{
  for (int p=0; p<m_nparts; ++p) {
    auto yv = part_view_nonconst<T>(p);
    auto pl = part_layout(p);
    for (int i=0; i<pl.size(); ++i) {
      yv(i) *= alpha;
    }
  }
}

template<typename T>
void Field::check_data_type () const
{
  EKAT_REQUIRE_MSG (m_data_type==get_data_type<T>(),
      "Error! Attempt to use wrong data type for this field:\n"
      "  - field name: " + m_name + "\n"
      "  - stored data type: " + e2str(m_data_type) + "\n"
      "  - input  data type: " + e2str(get_data_type<T>()) + "\n");
}

} // namespace cldera

#endif // CLDERA_FIELD_HPP

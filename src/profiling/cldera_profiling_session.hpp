#ifndef CLDERA_PROFILING_SESSION_HPP
#define CLDERA_PROFILING_SESSION_HPP

#include <ekat/mpi/ekat_comm.hpp>

namespace cldera {

class ProfilingSession {
public:
  static ProfilingSession& instance () {
    static ProfilingSession s;
    return s;
  };

  void init (const ekat::Comm& comm);
  void clean_up ();

private:
  ProfilingSession () = default;

  ekat::Comm  m_comm;
  bool m_inited = false;
};

} // namespace cldera

#endif // CLDERA_PROFILING_SESSION_HPP

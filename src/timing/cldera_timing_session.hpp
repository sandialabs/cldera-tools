#ifndef CLDERA_TIMING_SESSION_HPP
#define CLDERA_TIMING_SESSION_HPP

#include "cldera_timer_history.hpp"

#include <ekat/mpi/ekat_comm.hpp>

#include <map>
#include <ostream>
#include <string>

namespace cldera {
namespace timing {

// This class stores timers, so that it can later dump
// all stats to file. The class follows the singleton
// pattern, so the same data can be accessed from anywhere
// in the host app.

struct TimingSession
{
public:
  template<typename T>
  using strmap_t = std::map<std::string,T>;

  static TimingSession& instance () {
    static TimingSession ts;
    return ts;
  }

  // Dump all timer history stats to file
  void dump (      std::ostream& out,
             const ekat::Comm& comm) const;

  // Start/stop a given timer
  void start_timer (const std::string& timer_name);
  void stop_timer (const std::string& timer_name);

  // Toggle on/off actual timing
  void toggle_session (const bool on);

  // Clean up the class
  void clean_up ();
private:

  TimingSession () = default;

  // map[timer_name] = timer_history
  strmap_t<TimerHist>   timers;

  bool session_active = true;
};

} // namespace timing
} // namespace cldera

#endif // CLDERA_TIMING_SESSION_HPP

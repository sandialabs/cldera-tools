#include "cldera_timing_session.hpp"

#include <ekat/ekat_assert.hpp>

#include <fstream>
#include <iomanip>
#include <utility>

// A simple struct holding a double and an int
namespace {
struct DblInt{
  double d;
  int    i;
};
}

// Specialize ekat::Comm::all_reduce for DblInt
namespace ekat {
template<>
void Comm::all_reduce<DblInt> (const DblInt* my_vals, DblInt* result, const int count, const MPI_Op op) const
{
  check_mpi_inited();
  MPI_Allreduce(my_vals,result,count,MPI_DOUBLE_INT,op,m_mpi_comm);
}

} // namespace ekat

namespace cldera {
namespace timing {

void TimingSession::
dump (std::ostream& out, const ekat::Comm& comm) const
{
  out << "+---------------------------------------------------------------------------------------+\n";
  out << "|                                  CLDERA TIMING STATS                                  |\n";
  out << "+---------------------------------------------------------------------------------------+\n";
  out << "|              Timer Name                 |   N   |     Max (rank)   |     Min (rank)   |\n";
  out << "+---------------------------------------------------------------------------------------+\n";

  // Sanity check: all timers must be present on all ranks
  int my_ntimers = timers.size();
  int max_ntimers,min_ntimers;
  comm.all_reduce(&my_ntimers,&max_ntimers,1,MPI_MAX);
  comm.all_reduce(&my_ntimers,&min_ntimers,1,MPI_MIN);
  EKAT_REQUIRE_MSG (max_ntimers==min_ntimers,
      "Error! MPI ranks do not all store the same number of timers.\n");

  // std::ostringstream right_float_fmt;
  auto right_float_fmt = [](std::ostream& out) -> std::ostream& {
    out << std::right << std::setfill(' ') << std::setw(12) << std::setprecision(3) << std::fixed;
    return out;
  };
  for (const auto& it : timers) {
    const auto& n = it.first;
    const auto& t = it.second;

    out << "| " << std::left  << std::setfill(' ') << std::setw(40) << n;

    out << "| " << std::right << std::setfill(' ') << std::setw(5) << t.num_hist() << " ";

    DblInt max = {t.max().count(), comm.rank()};
    DblInt min = {t.min().count(), comm.rank()};

    DblInt gmax,gmin;
    comm.all_reduce(&min,&gmin,1,MPI_MINLOC);
    comm.all_reduce(&max,&gmax,1,MPI_MAXLOC);

    out << "| " << right_float_fmt << gmax.d << " (" << gmax.i << ") ";
    out << "| " << right_float_fmt << gmin.d << " (" << gmin.i << ") ";
    out << "|\n";
  }

  out << "+---------------------------------------------------------------------------------------+\n";
}

void TimingSession::
start_timer (const std::string& timer_name)
{
  if (session_active) {
    auto& timer = timers[timer_name];
    EKAT_REQUIRE_MSG (not timer.is_running(),
        "Error! Clock was already running!\n");
    timer.start();
  }
}

void TimingSession::
stop_timer (const std::string& timer_name)
{
  if (session_active) {
    auto& timer = timers[timer_name];
    EKAT_REQUIRE_MSG (timer.is_running(),
        "Error! Clock was not running!\n");
    timer.stop();
  }
}

void TimingSession::
toggle_session (const bool on)
{
  session_active = on;
}

void TimingSession::
clean_up ()
{
  TimingSession ts;
  std::swap(ts,*this);
}

} // namespace timing
} // namespace cldera

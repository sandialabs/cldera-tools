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
  int rank_width = 1;
  int comm_size = comm.size();
  while (comm_size>=10) {
    ++rank_width;
    comm_size /= 10;
  }

  int rw_left = rank_width / 2;
  int rw_right = rank_width - rw_left;

  std::string dashes(rank_width,'-');
  std::string dashes_beg = "+" + dashes;
  std::string dashes_end = dashes + "+\n";
  std::string space (rank_width,' ');
  std::string space_bar = space + "|\n";
  std::string bar_space = "|" + space;
  std::string sp_l (rw_left,' ');
  std::string sp_r (rw_right,' ');

  out << dashes_beg + "---------------------------------------------------------------------------------" + dashes_end;
  out << bar_space + "                              CLDERA TIMING STATS                                " + space_bar;
  out << dashes_beg + "---------------------------------------------------------------------------------" + dashes_end;
  out << "|              Timer Name                 | Count |" + sp_l + "   Max (rank)  " + sp_r + "|" + sp_l + "   Min (rank)  " + sp_r + "|\n";
  out << dashes_beg + "---------------------------------------------------------------------------------" + dashes_end;

  // Sanity check: all timers must be present on all ranks
  int my_ntimers = timers.size();
  int max_ntimers,min_ntimers;
  comm.all_reduce(&my_ntimers,&max_ntimers,1,MPI_MAX);
  comm.all_reduce(&my_ntimers,&min_ntimers,1,MPI_MIN);
  EKAT_REQUIRE_MSG (max_ntimers==min_ntimers,
      "Error! MPI ranks do not all store the same number of timers.\n");

  // std::ostringstream right_float_fmt;
  auto right_float_fmt = [](std::ostream& out) -> std::ostream& {
    out << std::right << std::setfill(' ') << std::setw(10) << std::setprecision(3) << std::fixed;
    return out;
  };
  for (const auto& it : timers) {
    const auto& n = it.first;
    const auto& t = it.second;

    out << "| " << std::left  << std::setfill(' ') << std::setw(40) << n;

    out << "| " << std::right << std::setfill(' ') << std::setw(5) << t.count() << " ";

    DblInt local = {t.elapsed().count(), comm.rank()};
    DblInt gmin,gmax;
    comm.all_reduce(&local,&gmin,1,MPI_MINLOC);
    comm.all_reduce(&local,&gmax,1,MPI_MAXLOC);

    out << "| " << right_float_fmt << gmax.d << " (" << std::setw(rank_width) << gmax.i << ") ";
    out << "| " << right_float_fmt << gmin.d << " (" << std::setw(rank_width) << gmin.i << ") ";
    out << "|\n";
  }

  out << dashes_beg + "---------------------------------------------------------------------------------" + dashes_end;
}

void TimingSession::
start_timer (const std::string& timer_name)
{
  if (session_active) {
    auto& timer = timers[timer_name];
    timer.start();
  }
}

void TimingSession::
stop_timer (const std::string& timer_name)
{
  if (session_active) {
    auto& timer = timers[timer_name];
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

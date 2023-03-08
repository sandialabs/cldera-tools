#ifndef CLDERA_TIMER_HISTORY_HPP
#define CLDERA_TIMER_HISTORY_HPP

#include <chrono>

namespace cldera {
namespace timing {

struct TimerHist {
  using duration_type = std::chrono::duration<double>;
  using clock_type    = std::chrono::high_resolution_clock;
  using time_type     = std::chrono::time_point<clock_type>;

  TimerHist () {
    d_min = duration_type::max();
    d_max = duration_type(0);
    nhist = 0;
  }

  void start () {
    running = true;
    t = clock_type::now();
  }

  void stop () {
    auto dt = clock_type::now() - t;
    running = false;
    update(dt);
  }

  duration_type max () const { return d_max; }
  duration_type min () const { return d_min; }
  int num_hist () const { return nhist; }
  bool is_running () const { return running; }

private:
  void update (const duration_type& dt) {
    d_max = d_max<dt ? dt : d_max;
    d_min = d_min>dt ? dt : d_min;
    ++nhist;
  }

  duration_type d_max;
  duration_type d_min;

  int nhist;
  bool running = false;

  time_type t;
};

} // namespace timing
} // namespace cldera

#endif // CLDERA_TIMER_HISTORY_HPP


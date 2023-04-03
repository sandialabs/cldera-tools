#ifndef CLDERA_TIMER_HPP
#define CLDERA_TIMER_HPP

#include <ekat/ekat_assert.hpp>

#include <chrono>

namespace cldera {
namespace timing {

struct Timer {
  using duration_type = std::chrono::duration<double>;
  using clock_type    = std::chrono::high_resolution_clock;
  using time_type     = std::chrono::time_point<clock_type>;

  Timer () {
    m_elapsed = duration_type(0);
    m_count = 0;
  }

  void start () {
    EKAT_REQUIRE_MSG (not m_running,
        "Error! Clock was already running!\n");
    m_running = true;
    t = clock_type::now();
  }

  void stop () {
    auto dt = clock_type::now() - t;
    EKAT_REQUIRE_MSG (m_running,
        "Error! Timer was not running!\n");
    m_running = false;
    update(dt);
  }

  duration_type elapsed () const { return m_elapsed; }
  int count () const { return m_count; }
  bool running () const { return m_running; }

private:
  void update (const duration_type& dt) {
    m_elapsed += dt;
    ++m_count;
  }

  duration_type m_elapsed;

  int m_count;
  bool m_running = false;

  time_type t;
};

} // namespace timing
} // namespace cldera

#endif // CLDERA_TIMER_HPP

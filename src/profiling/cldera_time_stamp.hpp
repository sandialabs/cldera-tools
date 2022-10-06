#ifndef CLDERA_TIME_STAMP_HPP
#define CLDERA_TIME_STAMP_HPP

#include <ekat/ekat_assert.hpp>

#include <iostream>
#include <iomanip>

namespace cldera {

inline bool is_leap (const int y) {
  if (y%4 != 0) {
    return false;
  }

  return y%100 !=0 or (y/100) % 4 == 0;
}

inline int days_in_month (const int mm, const int yy) {
  static constexpr int dom[] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
  static constexpr int leap_dom[] = {31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
  if (is_leap(yy)) {
    return leap_dom[mm-1];
  } else {
    return dom[mm-1];
  }
}

inline int days_in_year (const int yy) {
  return is_leap(yy) ? 366 : 365;
}

// An std::pair would work too, but ymd/tod convey
// more meaning than first/second.
// NOTE: ymd is in YYYYMMDD format, while tod is second of the day
struct TimeStamp {
  static constexpr int spd = 86400;

  TimeStamp (const int ymd, const int tod)
   : m_ymd(ymd), m_tod(tod)
  {
    EKAT_REQUIRE_MSG (year()>=1 && month()<=3000,
        "Error! Invalid year in YYYYMMDD input.\n"
        " - YYYYMMDD: " + std::to_string(m_ymd) + "\n"
        " - year   : " + std::to_string(year()) + "\n");
    EKAT_REQUIRE_MSG (month()>=1 && month()<=12,
        "Error! Invalid month in YYYYMMDD input.\n"
        " - YYYYMMDD: " + std::to_string(m_ymd) + "\n"
        " - month   : " + std::to_string(month()) + "\n");
    EKAT_REQUIRE_MSG (day()>=1 && day()<=days_in_month(month(),year()),
        "Error! Invalid month in YYYYMMDD input.\n"
        " - YYYYMMDD: " + std::to_string(m_ymd) + "\n"
        " - day     : " + std::to_string(day()) + "\n");

    EKAT_REQUIRE_MSG (tod>=0 && tod<86400,
        "Error! Invalid TOD input: " + std::to_string(m_tod) + "\n");
  }

  int year  () const { return (m_ymd / 100) / 100; }
  int month () const { return (m_ymd / 100) % 100; }
  int day   () const { return  m_ymd % 100;        }

  int hour () const { return (m_tod / 60) / 60; }
  int min  () const { return (m_tod / 60) % 60; }
  int sec  () const { return  m_tod % 60;       }

  double frac_of_year () const {
    double d = day()-1 + m_tod/double(spd);
    for (int m=1; m<month(); ++m) {
      d += days_in_month(m,year());
    }
    return d/days_in_year(year());
  }

  double frac_of_day () const {
    return (hour()*3600 + min()*60 + sec() ) / double(spd);
  }

  int ymd () const { return m_ymd; }
  int tod () const { return m_tod; }

  std::string to_string () const {
    std::stringstream ss;
    ss << *this;
    return ss.str();
  }

  TimeStamp& operator+= (const int secs);

  friend double operator-(const TimeStamp&, const TimeStamp&);
  friend bool operator==(const TimeStamp&, const TimeStamp&);
  friend bool operator<(const TimeStamp&, const TimeStamp&);
  friend std::ostream& operator<<(std::ostream&, const TimeStamp&);

private:
  int m_ymd;
  int m_tod;
};

inline TimeStamp& TimeStamp::operator+= (const int secs)
{
  EKAT_REQUIRE_MSG (secs>=0, "Error! Cannot unwind time, sorry.\n");
  m_tod += secs;
  int tod   = m_tod % spd;
  int carry = m_tod / spd;
  if (carry==0) {
    return *this;
  }
  m_tod = tod;
  int yy = year();
  int mm = month();
  int dd = day();
  while (carry>0) {
    ++dd;
    --carry;
    if (dd > days_in_month(mm,yy)) {
      dd = 1;
      ++mm;
      if (mm > 12) {
        mm = 1;
        ++yy;
      }
    }
  }
  m_ymd = yy*10000 + mm*100 + dd;
  return *this;
}

inline double operator- (const TimeStamp& lhs, const TimeStamp& rhs) {
  if (lhs.m_ymd<rhs.m_ymd ) {
    return rhs-lhs;
  }
  if (lhs.m_ymd==rhs.m_ymd) {
    return lhs.frac_of_day() - rhs.frac_of_day();
  }
  if (lhs.year()==rhs.year()) {
    return days_in_year(lhs.year())*(lhs.frac_of_year() - rhs.frac_of_year());
  }
  
  double diff = 0;

  // Remainder of year in rhs
  diff += days_in_year(rhs.year()) * (1-rhs.frac_of_year());

  for (int y=rhs.year()+1; y<lhs.year(); ++y) {
    // Full years in between lhs and rhs
    diff += days_in_year(y);
  }

  // Chunk of year in lhs
  diff += lhs.frac_of_year()*days_in_year(lhs.year());

  return diff;
}

inline bool operator< (const TimeStamp& lhs, const TimeStamp& rhs) {
  return lhs.m_ymd<rhs.m_ymd || (lhs.m_ymd==rhs.m_ymd && lhs.m_tod<rhs.m_tod);
}

inline bool operator== (const TimeStamp& lhs, const TimeStamp& rhs) {
  return lhs.m_ymd==rhs.m_ymd && lhs.m_tod==rhs.m_tod;
}

inline std::ostream& operator<< (std::ostream& out, const TimeStamp& t) {
  out << std::setw(4) << std::setfill('0') << t.year()
      << std::setw(1) << "-"
      << std::setw(2) << std::setfill('0') << t.month()
      << std::setw(1) << "-"
      << std::setw(2) << std::setfill('0') << t.day()
      << std::setw(1) << "-"
      << std::setw(5) << std::setfill('0') << t.m_tod;
  return out;
}

} // namespace cldera

#endif // CLDERA_TIME_STAMP_HPP

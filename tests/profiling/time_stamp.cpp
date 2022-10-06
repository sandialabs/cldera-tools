#include <catch2/catch.hpp>

#include "profiling/cldera_time_stamp.hpp"

#include <ekat/mpi/ekat_comm.hpp>

TEST_CASE ("test_manager") {
  using namespace cldera;

  const ekat::Comm comm(MPI_COMM_WORLD);

  // Free functions
  REQUIRE (not is_leap(1999));
  REQUIRE (is_leap(2000));
  REQUIRE (is_leap(2004));
  REQUIRE (not is_leap(2100));

  REQUIRE (days_in_year(1999)==365);
  REQUIRE (days_in_year(2000)==366);

  REQUIRE (days_in_month(2,2000)==29);
  REQUIRE (days_in_month(2,2001)==28);

  // Bad ctors
  REQUIRE_THROWS (TimeStamp(50001301,0));       // bad year (5000)
  REQUIRE_THROWS (TimeStamp(20001301,0));       // bad month
  REQUIRE_THROWS (TimeStamp(20001131,0));       // bad day
  REQUIRE_THROWS (TimeStamp(20000101,90000));   // bad TOD

  // Accessors
  {
    TimeStamp t1(20000330,43261);

    REQUIRE (t1.year()==2000);
    REQUIRE (t1.month()==3);
    REQUIRE (t1.day()==30);

    REQUIRE (t1.hour()==12);
    REQUIRE (t1.min()==1);
    REQUIRE (t1.sec()==1);

    REQUIRE (t1.frac_of_day()==43261/86400.0);
    REQUIRE (t1.frac_of_year()== (t1.frac_of_day() + 31+29+29)/366);
  }

  // Comparison
  {
    TimeStamp t1(20001130,43200);
    TimeStamp t2(20001130,43200);
    TimeStamp t3(20001030,43200);

    REQUIRE (t1==t2);
    REQUIRE (t3<t2);
  }

  // Difference
  {
    TimeStamp t1(20021130,43200);
    TimeStamp t2(20011130,43200);
    TimeStamp t3(20011231,43200);
    TimeStamp t4(20010101,43200);

    REQUIRE ( (t1-t2) == 365 );
    REQUIRE ( (t3-t4) == Approx(364) ); // Approx b/c of roundoffs in double calculations
  }

  // Streaming
  {
    TimeStamp t1(20000102,60);

    std::ostringstream ss;
    ss << t1;
    REQUIRE (ss.str()=="2000-01-02-00060");
  }

  // Update
  {
    TimeStamp t1(20000101,0);
    t1 += 86400;
    REQUIRE (t1.day()==2);
    t1 += 86400*30;
    REQUIRE (t1.month()==2);
    t1 += 86400*(29+31+30+31+30+31+31+30+31+30+31);
    REQUIRE (t1.year()==2001);
  }
}

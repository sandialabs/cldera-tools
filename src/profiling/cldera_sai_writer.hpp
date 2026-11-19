#ifndef CLDERA_SAI_WRITER_HPP
#define CLDERA_SAI_WRITER_HPP

#include "cldera_field.hpp"
#include "cldera_time_stamp.hpp"

#include <ekat/ekat_parameter_list.hpp>
#include <ekat/mpi/ekat_comm.hpp>

#include <fstream>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace cldera {

class ProfilingArchive;

/*
 * SAIWriter: writes SO2 injection rates into a *writable* host-app field.
 *
 * The host app (EAM) registers its 'sai' forcing array as a writable View
 * field (see cldera_set_field_part_data_nonconst_c). At every call to
 * cldera_compute_stats, after the controller (if any) has run, this class:
 *
 *   1. zeroes the whole field;
 *   2. for each injection site, finds every column inside the site's lat/lon
 *      box and, in that column, the model LAYER containing the target height;
 *   3. computes the global injection volume of each site (MPI all-reduce);
 *   4. converts the per-site rate [Tg/yr] (from the controller) into
 *      molec/cm3/s and writes it into the field;
 *   5. recomputes the written mass from exactly the same cells/thicknesses
 *      (mass-closure check), and logs what was injected and where;
 *   6. (BMW, G6) integrates the requested and written mass over each "period"
 *      (from one set_rates call to the next, e.g. a control interval such as MAMJ),
 *      and writes one line per site and period to the summary log.
 *
 * Box: a column belongs to a site if |lat-lat_site| < 'Box Half Width Lat Deg' and the
 *   longitude distance < 'Box Half Width Lon Deg' (both default to 'Box Half Width Deg').
 *   With 'Require Columns Per Site: true' (default), a site with no columns is an error.
 *
 * Vertical storage ordering (IMPORTANT):
 *   Model level index k=0 is the model TOP (zi[0] is the top interface).
 *   If the forcing file has an 'altitude' dimension (EAM file%alt_data), EAM
 *   stores the data BOTTOM-UP, and extfrc_set later flips it with
 *   data(:,pver:1:-1,:). In that case we must write at storage level
 *   nlev-1-k so that, after the flip, the rate ends up at model level k.
 *   This is controlled by the 'storage_bottom_up' flag passed by the host app.
 *
 * Required fields in the archive (all registered by EAM):
 *   lat, lon [deg] (ncol), area [sr] (ncol), zi [m above sfc] (ncol,ilev),
 *   plus the writable target field itself (ncol,lev).
 */
class SAIWriter {
public:
  struct Site {
    std::string name;
    double lat;       // [deg]
    double lon;       // [deg, 0-360]
    double alt_m;     // target height above surface [m]
  };

  SAIWriter (const ekat::Comm& comm,
             const ekat::ParameterList& params,
             const std::vector<Site>& sites);

  // Set by the host app when registering the target field. MUST be called
  // (apply() errors out otherwise), so we never silently assume an ordering.
  void set_storage_bottom_up (const bool bottom_up) { m_storage_bottom_up = bottom_up; m_storage_set = true; }

  const std::string& target_field_name () const { return m_field_name; }
  int num_sites () const { return m_sites.size(); }
  const std::vector<Site>& sites () const { return m_sites; }

  // Set the per-site injection rates [Tg/yr]. 'period' is a label for the mass budget
  // (e.g. the control interval name); a call to set_rates closes the current budget period
  // and opens a new one.
  void set_rates (const std::vector<double>& rates_tgyr, const TimeStamp& t,
                  const std::string& period = "");
  const std::vector<double>& rates () const { return m_rates_tgyr; }
  bool has_rates () const { return m_has_rates; }

  // Zero the target field and write the current rates into it.
  void apply (ProfilingArchive& archive, const TimeStamp& t);

  // End of run: close the current budget period (adds the last step) and flush logs.
  void finalize ();

private:
  struct ColEntry {
    int ipart;
    int icol;        // column index within the part
    int site;        // site index
    int lev_model;   // model layer index (0 = top)
    int lev_store;   // storage layer index (accounts for bottom-up storage)
    double dz_m;     // layer thickness [m]
    double area_sr;  // cell area [sr]
    double lat, lon; // [deg]
  };

  void build_entries (ProfilingArchive& archive);
  void write_map_log (const TimeStamp& t);
  void write_inj_log (const TimeStamp& t,
                      const std::vector<double>& molec,
                      const std::vector<double>& written_tgyr);

  // Mass budget (BMW): integrate requested and written mass over each period
  void budget_accumulate_until (const TimeStamp& t);
  void budget_close_period (const TimeStamp& t_end, const bool run_end);

  ekat::Comm          m_comm;
  std::string         m_field_name;
  std::vector<Site>   m_sites;
  double              m_d_lat;              // box half width in latitude [deg]
  double              m_d_lon;              // box half width in longitude [deg]
  bool                m_require_columns;    // error out if a site has no columns
  int                 m_log_every;          // log injection every N calls (root only)
  bool                m_storage_bottom_up = false;
  bool                m_storage_set = false;

  std::vector<double> m_rates_tgyr;
  bool                m_has_rates = false;
  double              m_default_rate_tgyr;  // used until controller provides rates

  std::vector<ColEntry> m_entries;          // local columns hit by some site
  std::vector<double>   m_glb_volume;       // per-site global volume [m*sr]
  std::vector<int>      m_glb_ncols;        // per-site global number of columns
  int                   m_num_calls = 0;

  // Mass budget state (identical on all ranks; only root writes)
  std::string           m_period_label;
  TimeStamp             m_period_start;
  bool                  m_period_open = false;
  double                m_period_seconds = 0;
  std::vector<double>   m_period_req_tg;    // integral of requested rate [Tg]
  std::vector<double>   m_period_wrt_tg;    // integral of written (closure) rate [Tg]
  std::vector<double>   m_period_rate;      // rate set at the period start [Tg/yr]
  bool                  m_has_prev = false;
  TimeStamp             m_prev_t;
  std::vector<double>   m_prev_written;     // written rate of the last apply [Tg/yr]
  std::vector<double>   m_prev_requested;   // requested rate of the last apply [Tg/yr]
  double                m_dt_min = -1;      // smallest spacing between apply calls [s]

  // NOTE: shared_ptr (not a plain ofstream) so that SAIWriter is copyable, which is
  //       required to store it in the ProfilingContext (ekat::any).
  std::shared_ptr<std::ofstream> m_inj_log;  // root only
  std::shared_ptr<std::ofstream> m_sum_log;  // root only
  std::string         m_map_log_name;
};

} // namespace cldera

#endif // CLDERA_SAI_WRITER_HPP

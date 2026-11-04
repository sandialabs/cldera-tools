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
 *      (mass-closure check), and logs what was injected and where.
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

  // Set the per-site injection rates [Tg/yr]
  void set_rates (const std::vector<double>& rates_tgyr, const TimeStamp& t);
  const std::vector<double>& rates () const { return m_rates_tgyr; }
  bool has_rates () const { return m_has_rates; }

  // Zero the target field and write the current rates into it.
  void apply (ProfilingArchive& archive, const TimeStamp& t);

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

  ekat::Comm          m_comm;
  std::string         m_field_name;
  std::vector<Site>   m_sites;
  double              m_d_latlon;           // box half width [deg]
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

  // NOTE: shared_ptr (not a plain ofstream) so that SAIWriter is copyable, which is
  //       required to store it in the ProfilingContext (ekat::any).
  std::shared_ptr<std::ofstream> m_inj_log;  // root only
  std::string         m_map_log_name;
};

} // namespace cldera

#endif // CLDERA_SAI_WRITER_HPP

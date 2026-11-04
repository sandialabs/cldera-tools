#include "cldera_sai_writer.hpp"

#include "cldera_profiling_archive.hpp"

#include <ekat/ekat_assert.hpp>

#include <cmath>
#include <cstdio>
#include <iomanip>
#include <set>
#include <sstream>

namespace cldera {

namespace {

// Physical constants. NOTE: EAM uses a 365-day (noleap) calendar.
constexpr double kAvogadro    = 6.02214076e23;  // molec/mol
constexpr double kMolWeightSO2 = 64.066;        // g/mol
constexpr double kSecPerYear  = 365.0*86400.0;  // s/yr
constexpr double kGramPerTg   = 1.0e12;         // g/Tg
constexpr double kEarthRadius = 6.37122e6;      // m (EAM shr_const_rearth)
constexpr double kPi          = 3.14159265358979323846;

inline double rad2deg (const double r) { return r*180.0/kPi; }

// Distance in longitude, accounting for wrap-around at 0/360
inline double lon_dist (const double a, const double b) {
  double d = std::fabs(a-b);
  d = std::fmod(d,360.0);
  return d>180.0 ? 360.0-d : d;
}

// Tg/yr -> molecules/s
inline double tgyr_to_molec_per_s (const double tgyr) {
  return tgyr * kGramPerTg / kSecPerYear / kMolWeightSO2 * kAvogadro;
}

} // anonymous namespace

SAIWriter::
SAIWriter (const ekat::Comm& comm,
           const ekat::ParameterList& params_in,
           const std::vector<Site>& sites)
 : m_comm (comm)
 , m_sites (sites)
{
  // Non-const copy: ParameterList::get with a default value is non-const
  ekat::ParameterList params = params_in;
  m_field_name        = params.get<std::string>("Field","forcing_sai");
  m_d_latlon          = params.get<double>("Box Half Width Deg");
  m_log_every         = params.get<int>("Log Frequency",1);
  m_default_rate_tgyr = params.get<double>("Default Rate Tg/yr",0.0);
  const auto inj_log_name = params.get<std::string>("Injection Log File","sai_injection_log.csv");
  m_map_log_name      = params.get<std::string>("Column Map File","sai_injection_columns.csv");

  EKAT_REQUIRE_MSG (m_sites.size()>0,
      "[SAIWriter] Error! No injection sites defined (check 'Injections' in the controller yaml).\n");
  EKAT_REQUIRE_MSG (m_d_latlon>0,
      "[SAIWriter] Error! 'Box Half Width Deg' must be positive.\n");
  EKAT_REQUIRE_MSG (m_log_every>=1,
      "[SAIWriter] Error! 'Log Frequency' must be >= 1.\n");

  m_rates_tgyr.assign(m_sites.size(),m_default_rate_tgyr);

  if (m_comm.am_i_root()) {
    m_inj_log = std::make_shared<std::ofstream>(inj_log_name);
    EKAT_REQUIRE_MSG (m_inj_log->good(),
        "[SAIWriter] Error! Could not open injection log file '" + inj_log_name + "'.\n");
    *m_inj_log << "# SAI injection log written by cldera-tools SAIWriter\n"
              << "# rate_tgyr: requested rate (from controller); written_tgyr: mass written into the\n"
              << "#   forcing field, recomputed from the same cells/layer thicknesses;\n"
              << "# molec_cm3_s: uniform rate written in each cell of the site's box.\n"
              << "# A 'TOTAL' row per call gives the total as EAM's SO2_CLXF will see it (re-reads the\n"
              << "#   whole field, applies the storage flip, and integrates over all levels).\n"
              << "time,site,lat,lon,alt_m,ncols,rate_tgyr,volume_m3,molec_cm3_s,written_tgyr\n";
    m_inj_log->flush();

    printf(" [CLDERA] SAIWriter: target field '%s', %d sites, box half width %.3f deg\n",
           m_field_name.c_str(), (int)m_sites.size(), m_d_latlon);
    for (const auto& s : m_sites) {
      printf(" [CLDERA]   site %-10s lat=%8.3f lon=%8.3f alt=%8.1f m\n",
             s.name.c_str(), s.lat, s.lon, s.alt_m);
    }
  }
}

void SAIWriter::
set_rates (const std::vector<double>& rates_tgyr, const TimeStamp& t)
{
  EKAT_REQUIRE_MSG (rates_tgyr.size()==m_sites.size(),
      "[SAIWriter] Error! Number of rates does not match number of sites.\n"
      " - num rates: " + std::to_string(rates_tgyr.size()) + "\n"
      " - num sites: " + std::to_string(m_sites.size()) + "\n");
  for (auto r : rates_tgyr) {
    EKAT_REQUIRE_MSG (std::isfinite(r) and r>=0,
        "[SAIWriter] Error! Injection rates must be finite and non-negative.\n");
  }
  m_rates_tgyr = rates_tgyr;
  m_has_rates = true;

  if (m_comm.am_i_root()) {
    printf(" [CLDERA] SAIWriter: new injection rates at %s [Tg/yr]:", t.to_string().c_str());
    for (auto r : m_rates_tgyr) printf(" %g", r);
    printf("\n");
  }
}

void SAIWriter::
build_entries (ProfilingArchive& archive)
{
  // Columns that fall in each site's box. Lat/lon/area never change, so we do this only once.
  // The target layer (and its thickness) is recomputed at every call in apply(), since zi changes.
  const auto& lat  = archive.get_field("lat");
  const auto& lon  = archive.get_field("lon");
  const auto& area = archive.get_field("area");
  const auto& tgt  = archive.get_field(m_field_name);

  EKAT_REQUIRE_MSG (lat.nparts()==tgt.nparts() and lon.nparts()==tgt.nparts() and area.nparts()==tgt.nparts(),
      "[SAIWriter] Error! lat/lon/area and target field have different number of parts.\n");

  m_entries.clear();
  for (int p=0; p<tgt.nparts(); ++p) {
    auto lat_v  = lat.part_view<const Real>(p);
    auto lon_v  = lon.part_view<const Real>(p);
    auto area_v = area.part_view<const Real>(p);
    const int ncols = lat.part_layout(p).extent("ncol");
    for (int i=0; i<ncols; ++i) {
      // EAM state lat/lon are in radians
      const double la = rad2deg(lat_v(i));
      const double lo = rad2deg(lon_v(i));
      for (int s=0; s<num_sites(); ++s) {
        if (std::fabs(la-m_sites[s].lat)<m_d_latlon and lon_dist(lo,m_sites[s].lon)<m_d_latlon) {
          ColEntry e;
          e.ipart = p;
          e.icol = i;
          e.site = s;
          e.lev_model = -1;
          e.lev_store = -1;
          e.dz_m = 0;
          e.area_sr = area_v(i);
          e.lat = la;
          e.lon = lo;
          m_entries.push_back(e);
        }
      }
    }
  }
}

void SAIWriter::
apply (ProfilingArchive& archive, const TimeStamp& t)
{
  EKAT_REQUIRE_MSG (archive.has_field(m_field_name),
      "[SAIWriter] Error! Target field '" + m_field_name + "' not found in the archive.\n"
      "  Did the host app register it?\n");
  auto& tgt = archive.get_field(m_field_name);
  EKAT_REQUIRE_MSG (m_storage_set,
      "[SAIWriter] Error! The host app did not specify the vertical storage order of the target field.\n"
      "  Call cldera_set_sai_storage_bottom_up when registering it.\n");
  EKAT_REQUIRE_MSG (tgt.committed() and not tgt.is_read_only(),
      "[SAIWriter] Error! Target field '" + m_field_name + "' must be committed and writable.\n");
  EKAT_REQUIRE_MSG (tgt.layout().rank()==2 and tgt.layout().names()[0]=="lev" and tgt.layout().names()[1]=="ncol",
      "[SAIWriter] Error! Target field must have layout (lev,ncol) [C ordering].\n"
      " - layout: " + tgt.layout().to_string() + "\n");

  const auto& zi = archive.get_field("zi");
  const int nlev = tgt.layout().extent("lev");
  EKAT_REQUIRE_MSG (zi.layout().extent("ilev")==nlev+1,
      "[SAIWriter] Error! zi must have nlev+1 interfaces.\n");

  if (m_num_calls==0) {
    build_entries(archive);
  }

  // 1. Zero the whole field: the controller is the only SAI source.
  tgt.deep_copy(0.0);

  // 2. Find target layer and its thickness for each entry; accumulate local volumes
  std::vector<double> loc_vol(num_sites(),0.0);
  for (auto& e : m_entries) {
    auto zi_v = zi.part_nd_view<const Real,2>(e.ipart);  // (ilev, col)
    const double z = m_sites[e.site].alt_m;
    // zi decreases with index (0=model top). The first interface below z
    // is ilev; the layer containing z is k=ilev-1, bounded by zi(k) and zi(k+1).
    int k = -1;
    for (int ilev=1; ilev<=nlev; ++ilev) {
      if (zi_v(ilev,e.icol) < z) {
        k = ilev-1;
        break;
      }
    }
    if (k<0) {
      // Target below the lowest interface (should not happen). Skip column.
      e.lev_model = e.lev_store = -1;
      e.dz_m = 0;
      continue;
    }
    e.lev_model = k;
    // EAM alt_data files are stored bottom-up; extfrc_set flips them back with
    // data(:,pver:1:-1,:). Write at the mirrored storage index so that, after the
    // flip, the value lands at model layer k.
    e.lev_store = m_storage_bottom_up ? (nlev-1-k) : k;
    e.dz_m = zi_v(k,e.icol) - zi_v(k+1,e.icol);
    loc_vol[e.site] += e.area_sr*e.dz_m;
  }

  // 3. Global volume of each site [m*sr]
  m_glb_volume.assign(num_sites(),0.0);
  m_comm.all_reduce(loc_vol.data(),m_glb_volume.data(),num_sites(),MPI_SUM);

  // 4. Convert Tg/yr -> molec/cm3/s, and write
  std::vector<double> molec(num_sites(),0.0);
  for (int s=0; s<num_sites(); ++s) {
    const double vol_cm3 = m_glb_volume[s]*kEarthRadius*kEarthRadius*1.0e6;
    if (vol_cm3>0) {
      molec[s] = tgyr_to_molec_per_s(m_rates_tgyr[s]) / vol_cm3;
    } else if (m_num_calls==0 and m_comm.am_i_root()) {
      printf(" [CLDERA] SAIWriter WARNING: site '%s' has zero injection volume (no columns in box?).\n"
             "          Nothing will be injected at this site.\n", m_sites[s].name.c_str());
    }
  }

  std::vector<double> loc_written(num_sites(),0.0);
  for (const auto& e : m_entries) {
    if (e.lev_store<0) continue;
    auto f_v = tgt.part_nd_view_nonconst<Real,2>(e.ipart); // (lev, col)
    // Additive, in case boxes of different sites overlap
    f_v(e.lev_store,e.icol) += molec[e.site];
    loc_written[e.site] += molec[e.site]*e.area_sr*e.dz_m;
  }

  // 5a. Mass closure per site [Tg/yr]
  std::vector<double> written(num_sites(),0.0);
  m_comm.all_reduce(loc_written.data(),written.data(),num_sites(),MPI_SUM);
  const double to_tgyr = kEarthRadius*kEarthRadius*1.0e6*kSecPerYear/kAvogadro*kMolWeightSO2/kGramPerTg;
  for (auto& w : written) w *= to_tgyr;

  // 5b. Independent check, mimicking EAM's SO2_CLXF: re-read the WHOLE field in all
  //     affected columns, map storage level -> model level (undoing the flip), and
  //     integrate with the model layer thickness. If the flip or the indexing were wrong,
  //     this would differ from the requested total.
  double loc_clxf = 0;
  {
    std::set<std::pair<int,int>> visited;
    for (const auto& e : m_entries) {
      if (not visited.insert({e.ipart,e.icol}).second) continue;
      auto f_v  = tgt.part_nd_view<const Real,2>(e.ipart);
      auto zi_v = zi.part_nd_view<const Real,2>(e.ipart);
      for (int j=0; j<nlev; ++j) {
        const int k = m_storage_bottom_up ? (nlev-1-j) : j;
        loc_clxf += f_v(j,e.icol)*(zi_v(k,e.icol)-zi_v(k+1,e.icol))*e.area_sr;
      }
    }
  }
  double clxf = 0;
  m_comm.all_reduce(&loc_clxf,&clxf,1,MPI_SUM);
  clxf *= to_tgyr;

  // 6. Logs
  if (m_num_calls==0) {
    write_map_log(t);
  }
  if (m_comm.am_i_root() and m_num_calls%m_log_every==0) {
    write_inj_log(t,molec,written);
    double req = 0;
    for (auto r : m_rates_tgyr) req += r;
    *m_inj_log << t.to_string() << ",TOTAL,,,,," << req << ",,," << clxf << "\n";
    m_inj_log->flush();
    if (req>0 and std::fabs(clxf-req)>1e-6*req) {
      printf(" [CLDERA] SAIWriter WARNING at %s: CLXF-equivalent total %.6g Tg/yr != requested %.6g Tg/yr\n",
             t.to_string().c_str(),clxf,req);
    }
  }

  ++m_num_calls;
}

void SAIWriter::
write_inj_log (const TimeStamp& t,
               const std::vector<double>& molec,
               const std::vector<double>& written_tgyr)
{
  // Per-site numbers. Column counts were computed (globally) in write_map_log.
  for (int s=0; s<num_sites(); ++s) {
    const auto& site = m_sites[s];
    *m_inj_log << t.to_string() << "," << site.name << ","
              << site.lat << "," << site.lon << "," << site.alt_m << ","
              << (s<(int)m_glb_ncols.size() ? m_glb_ncols[s] : -1) << ","
              << std::setprecision(10)
              << m_rates_tgyr[s] << ","
              << m_glb_volume[s]*kEarthRadius*kEarthRadius << ","
              << molec[s] << ","
              << written_tgyr[s] << "\n";
  }
}

void SAIWriter::
write_map_log (const TimeStamp& t)
{
  // Gather all (local) column entries to root, and write one row per (site,column).
  constexpr int nvals = 7; // site, lat, lon, lev_model, lev_store, dz, area
  std::vector<double> loc;
  loc.reserve(nvals*m_entries.size());
  for (const auto& e : m_entries) {
    loc.push_back(e.site);
    loc.push_back(e.lat);
    loc.push_back(e.lon);
    loc.push_back(e.lev_model);
    loc.push_back(e.lev_store);
    loc.push_back(e.dz_m);
    loc.push_back(e.area_sr);
  }

  const MPI_Comm mpicomm = m_comm.mpi_comm();
  const int root = m_comm.root_rank();
  int nloc = loc.size();
  std::vector<int> counts, displs;
  if (m_comm.am_i_root()) {
    counts.resize(m_comm.size());
  }
  MPI_Gather(&nloc,1,MPI_INT,counts.data(),1,MPI_INT,root,mpicomm);
  std::vector<double> all;
  int tot = 0;
  if (m_comm.am_i_root()) {
    displs.resize(m_comm.size(),0);
    for (int r=0; r<m_comm.size(); ++r) {
      displs[r] = tot;
      tot += counts[r];
    }
    all.resize(std::max(tot,1));
  }
  MPI_Gatherv(loc.data(),nloc,MPI_DOUBLE,all.data(),counts.data(),displs.data(),MPI_DOUBLE,root,mpicomm);

  // Per-site global column count (needed by the injection log on root)
  std::vector<int> loc_n(num_sites(),0);
  for (const auto& e : m_entries) ++loc_n[e.site];
  m_glb_ncols.assign(num_sites(),0);
  m_comm.all_reduce(loc_n.data(),m_glb_ncols.data(),num_sites(),MPI_SUM);

  if (m_comm.am_i_root()) {
    std::ofstream out(m_map_log_name);
    out << "# Columns and layers receiving SAI injection (written once, at " << t.to_string() << ")\n"
        << "# lev_model: model layer index (0-based, 0=model top); lev_store: index used in the\n"
        << "#   forcing array (differs from lev_model if storage is bottom-up, i.e. EAM alt_data);\n"
        << "# Fortran (1-based) model level = lev_model+1. dz_m/area_sr at the time of writing.\n"
        << "# storage_bottom_up = " << (m_storage_bottom_up ? 1 : 0) << "\n"
        << "site,site_lat,site_lon,site_alt_m,col_lat,col_lon,lev_model,lev_store,dz_m,area_sr\n";
    for (int i=0; i+nvals<=tot; i+=nvals) {
      const int s = static_cast<int>(all[i]);
      out << m_sites[s].name << "," << m_sites[s].lat << "," << m_sites[s].lon << "," << m_sites[s].alt_m << ","
          << std::setprecision(8) << all[i+1] << "," << all[i+2] << ","
          << static_cast<int>(all[i+3]) << "," << static_cast<int>(all[i+4]) << ","
          << all[i+5] << "," << all[i+6] << "\n";
    }
    out.close();

    printf(" [CLDERA] SAIWriter: injection columns per site:");
    for (int s=0; s<num_sites(); ++s) printf(" %s=%d", m_sites[s].name.c_str(), m_glb_ncols[s]);
    printf("\n [CLDERA] SAIWriter: column map written to '%s'\n", m_map_log_name.c_str());
  }
}

} // namespace cldera

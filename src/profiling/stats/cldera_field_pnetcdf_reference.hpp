#ifndef SRC_PROFILING_STATS_CLDERA_FIELD_PNETCDF_REFERENCE_HPP_
#define SRC_PROFILING_STATS_CLDERA_FIELD_PNETCDF_REFERENCE_HPP_

#include "profiling/stats/cldera_field_stat.hpp"
#include "profiling/stats/cldera_field_stat_utils.hpp"
#include "profiling/stats/cldera_field_stat_factory.hpp"
#include "profiling/stats/cldera_field_identity.hpp"

#include "io/cldera_pnetcdf.hpp"

#include <ekat/ekat_parameter_list.hpp>
#include <ekat/mpi/ekat_comm.hpp>

#include <memory>

namespace cldera {

class FieldPnetcdfReference : public FieldSinglePartStat
{
public:
  FieldPnetcdfReference (const ekat::Comm& comm,
                         const ekat::ParameterList& pl)
   : FieldSinglePartStat (comm,pl)
   , m_lat_bounds(pl.get<std::vector<Real>>("Latitude Bounds"))
   , m_lon_bounds(pl.get<std::vector<Real>>("Longitude Bounds"))
   , m_pnetcdf_filename(pl.get<std::string>("Pnetcdf Filename"))
   , m_ref_field_name(pl.get<std::string>("Reference Field Name"))
   , m_ref_deviation_name(pl.get<std::string>("Reference Deviation Field Name"))
   , m_time_step_ratio(pl.get<int>("Time Step Ratio"))
   , m_mask_val(m_params.get("Mask Value",0.0))
  { /* Nothing to do here */ }

  ~FieldPnetcdfReference() {
    io::pnetcdf::close_file(*m_pnetcdf_file);
  }

  std::string type() const override { return "pnetcdf_reference"; }

  void reset () {
    m_lat = m_lon = m_colgids = nullptr;
    m_timeindex = m_pnetcdf_timeindex = 0;;
    if (m_pnetcdf_file) {
      io::pnetcdf::close_file(*m_pnetcdf_file);
    }
    m_ref_var_dims.clear();
    m_refvar_data.clear();
    m_refvardev_data.clear();

    m_inited = false;
  }

  void initialize(const std::shared_ptr<const Field>& lat, const std::shared_ptr<const Field>& lon,
                  const std::shared_ptr<const Field>& col_gids) {
    if (m_inited) {
      return;
    }

    EKAT_REQUIRE_MSG(lat->name() == "lat" && lon->name() == "lon",
        "Error! Field names are not lat and lon!\n");
    m_lat = lat;
    m_lon = lon;
    m_colgids = col_gids;
    m_timeindex = 0;
    m_pnetcdf_timeindex = 0;

    // open the pnetcdf file
    m_pnetcdf_file = io::pnetcdf::open_file(m_pnetcdf_filename,m_comm,io::pnetcdf::IOMode::Read);
    
    // grab relevant dims (time x lev x col)
    auto timedim = m_pnetcdf_file->dims.at("time");
    auto levdim = m_pnetcdf_file->dims.at("lev");
    auto coldim = m_pnetcdf_file->dims.at("ncol");

    // set "realistic" field params using the dims
    int m_ntime = timedim->len;
    int m_nlev = levdim->len;
    int m_ncol = coldim->len;
    
    // we need a decomp to read the file - use the column GIDs from the archive
    // first, make an identity stat to get them on one part
    auto identity_pl = ekat::ParameterList("identity");
    auto column_stat = create_stat(identity_pl, m_comm);
    auto column_id_stat = dynamic_cast<FieldIdentity *>(column_stat.get());
    auto singlepartgids = column_id_stat->compute(*m_colgids);

    // second, use the single-part stat to make a decomp
    auto ngids = singlepartgids.layout().size();
    std::vector<int> my_cols(ngids);
    auto gids = singlepartgids.data<int>();
    std::copy(gids,gids+ngids,my_cols.begin());
    for(int i=0; i<my_cols.size(); ++i) {
      my_cols[i] = my_cols[i] - 1; // GH: offset column map by 1 since E3SM orders starting at 1, but we need to start at 0
    }
    io::pnetcdf::add_decomp(*m_pnetcdf_file,"ncol",my_cols);

    // TODO: use something like this to handle timestepping sync
    // m_times = std::vector<double>(m_ntime,0.0);
    // for(unsigned int i=0; i<m_ntime; ++i) {
    //   read_var(*m_pnetcdf_file,"time",&(m_times[i]),i);
    // }

    // grab variables from the file now that we've decomposed it
    // we want to copy the data we need and then close the file
    std::cout << "Looking for ref_var " + m_ref_field_name << std::endl;
    auto ref_var = m_pnetcdf_file->vars.at(m_ref_field_name);
    std::cout << "Looking for ref_dev_var " + m_ref_deviation_name << std::endl;
    auto ref_dev_var = m_pnetcdf_file->vars.at(m_ref_deviation_name);

    // store layout, and check that the fields have the same layouts and sizes
    m_ref_var_dims = std::vector<int>(ref_var->dims.size(),0);
    EKAT_REQUIRE_MSG(ref_var->dims.size() == ref_dev_var->dims.size(),
        "Error! The pnetcdf reference field and reference deviation field do not have the dims size!\n");
    for(unsigned int i=0; i<ref_var->dims.size(); ++i) {
      m_ref_var_dims[i] = ref_var->dims[i]->len;
      EKAT_REQUIRE_MSG(ref_var->dims[i]->len == ref_dev_var->dims[i]->len,
          "Error! The reference field and reference deviation field have unequal len on dim " + std::to_string(i) + "\n");
      EKAT_REQUIRE_MSG(ref_var->dims[i]->name == ref_dev_var->dims[i]->name,
          "Error! The reference field and reference deviation field have unequal name on dim " + std::to_string(i) + "\n");
    }

    // load the correct slice of the variable of interest
    m_refvar_data = std::vector<float>(m_nlev*m_ncol,0.0); // reference variable data
    m_refvardev_data = std::vector<float>(m_nlev*m_ncol,0.0); // reference variable deviation data
    read_var(*m_pnetcdf_file,m_ref_field_name,m_refvar_data.data(), m_pnetcdf_timeindex);
    read_var(*m_pnetcdf_file,m_ref_deviation_name,m_refvardev_data.data(),m_pnetcdf_timeindex);

    m_inited = true;
  }

protected:
  void compute_impl(const Field& f, Field& stat) const override {
    EKAT_REQUIRE_MSG(m_lat != nullptr && m_lon != nullptr && m_colgids != nullptr,
        "Error! lat/lon/col_gids fields not initialized!\n");

    const int nparts = f.nparts();
    EKAT_REQUIRE_MSG(nparts == m_lat->nparts() && nparts == m_lon->nparts() && nparts == m_colgids->nparts(),
        "Error! Field " + f.name() + " should have the same number of parts as lat/lon/col_gids!\n");

    const auto dt = f.data_type();
    if(dt==IntType) {
      do_compute_impl<int>(f,stat);
    } else if(dt==RealType) {
      do_compute_impl<Real>(f,stat);
    } else {
      EKAT_ERROR_MSG("[FieldPnetcdfReference] Unrecognized/unsupported data type (" + e2str(dt) + ")\n");
    }
  }

  template <typename T>
  void do_compute_impl(const Field& f, Field& stat) const {

    // this increments m_timestep and only updates the data every time_step_ratio steps
    update_refvar_data();

    // use the provided strides to make things easier
    // assuming the same strides for the field and the pnetcdf file are applicable
    const auto& stat_strides = compute_stat_strides(f.layout());
    auto pnetcdf_reference_field = stat.view_nonconst<T>();
    const int field_rank = f.layout().rank();
    const int field_part_dim = f.part_dim();
    for (int ipart = 0; ipart < f.nparts(); ++ipart) {
      const auto& field_part_data = f.part_data<const T>(ipart);
      const auto& lat_part_data = m_lat->part_data<const Real>(ipart);
      const auto& lon_part_data = m_lon->part_data<const Real>(ipart);
      const auto& field_part_layout = f.part_layout(ipart);
      const auto& field_part_dims = field_part_layout.dims();
      for (int field_part_index = 0; field_part_index < field_part_layout.size(); ++field_part_index) {
        const int stat_index = compute_stat_index(
            ipart, field_part_index, field_rank, field_part_dim, field_part_dims, stat_strides);
        const int geo_part_index = compute_geo_part_index(field_part_index, field_rank, field_part_dim, field_part_dims);
        const Real lat_val = lat_part_data[geo_part_index];
        const Real lon_val = lon_part_data[geo_part_index];
        if (lat_val > m_lat_bounds.min && lat_val < m_lat_bounds.max &&
            lon_val > m_lon_bounds.min && lon_val < m_lon_bounds.max) {
          double tmp = (m_refvar_data[field_part_index] - field_part_data[field_part_index])/(m_refvardev_data[field_part_index]);
          pnetcdf_reference_field(stat_index) = tmp*tmp;
        }
        else
          pnetcdf_reference_field(stat_index) = m_mask_val;
      }
    }
  }

  void update_refvar_data() const {
    // TODO: handle the time stepping sync between E3SM and pnetcdf better
    m_timeindex++;
    m_pnetcdf_timeindex = m_timeindex/m_time_step_ratio;

    if(m_timeindex % m_time_step_ratio == 0) {
      // grab variables from the file now that we've decomposed it
      // we want to copy the data we need and then close the file
      auto ref_var = m_pnetcdf_file->vars.at(m_ref_field_name);
      auto ref_dev_var = m_pnetcdf_file->vars.at(m_ref_deviation_name);

      // load the correct slice of the variable of interest
      m_refvar_data = std::vector<float>(m_nlev*m_ncol,0.0); // reference variable data
      m_refvardev_data = std::vector<float>(m_nlev*m_ncol,0.0); // reference variable deviation data
      read_var(*m_pnetcdf_file,m_ref_field_name,m_refvar_data.data(),m_pnetcdf_timeindex);
      read_var(*m_pnetcdf_file,m_ref_deviation_name,m_refvardev_data.data(),m_pnetcdf_timeindex);
    }
  }

  /// the communicator
  const ekat::Comm m_comm;
  /// the fields from E3SM
  std::shared_ptr<const Field> m_lat, m_lon, m_colgids;
  
  /// bounds for masking latitude and longitude (radians)
  const Bounds m_lat_bounds, m_lon_bounds;
  /// mask value (default: 0.0)
  const Real m_mask_val;
  
  /// filename for the pnetcdf reference data
  const std::string m_pnetcdf_filename;
  /// pointer to the pnetcdf file
  std::shared_ptr<io::pnetcdf::NCFile> m_pnetcdf_file;
  
  /// number of levels in the reference data
  int m_nlev;
  /// number of columns in the reference data
  int m_ncol;
  
  /// dimensions of reference data, may be beneficial later
  std::vector<int> m_ref_var_dims;
  /// the number of times this stat has been called (assumes once per time step)
  mutable int m_timeindex;
  /// the time index for the pnetcdf file
  mutable int m_pnetcdf_timeindex;
  /// the time step ratio such that m_timeindex = m_pnetcdf_timeindex / m_time_step_ratio
  const int m_time_step_ratio;
  
  /// name of the reference field
  const std::string m_ref_field_name;
  /// values of reference data
  mutable std::vector<float> m_refvar_data;
  
  /// name of the reference deviation field
  const std::string m_ref_deviation_name;
  /// values of reference deviation data
  mutable std::vector<float> m_refvardev_data;

  bool m_inited = false;
};

} // namespace cldera

#endif /* SRC_PROFILING_STATS_CLDERA_FIELD_PNETCDF_REFERENCE_HPP_ */

#ifndef CLDERA_FIELD_STAT_PIPE_HPP
#define CLDERA_FIELD_STAT_PIPE_HPP

#include "cldera_field_stat.hpp"

namespace cldera
{

class FieldStatPipe : public FieldStat
{
public:
  FieldStatPipe (const ekat::Comm& comm,
                 const ekat::ParameterList& pl)
   : FieldStat(comm,pl)
  {
    auto& outer_pl = m_params.sublist("outer");
    auto& inner_pl = m_params.sublist("inner");

    // This ensures that all pipe stats' inner/outer stats have different names.
    // If names clash, the start/stop timers may run into problems like calling
    // start_timer twice in a row without a corresponding stop_timer in between.
    if (not outer_pl.isParameter("name")) {
      outer_pl.set("name",m_name + "_outer");
    }
    if (not inner_pl.isParameter("name")) {
      inner_pl.set("name",m_name + "_inner");
    }

    auto& f = StatFactory::instance();
    m_outer = f.create(outer_pl.get<std::string>("type"),comm,outer_pl);
    m_inner = f.create(inner_pl.get<std::string>("type"),comm,inner_pl);
  }

  std::string type () const { return "pipe"; }

  FieldLayout stat_layout (const FieldLayout& fl) const {
    return m_outer->stat_layout(m_inner->stat_layout(fl));
  }

  void create_stat_field () {
    m_outer->create_stat_field();
    m_stat_field = m_outer->get_stat_field();
    m_stat_field.rename(m_name);
  }

  std::vector<std::string> get_aux_fields_names () const {
    std::vector<std::string> aux_fnames;
    for (const auto& it : m_outer->get_aux_fields_names()) {
      aux_fnames.push_back(it);
    }
    for (const auto& it : m_inner->get_aux_fields_names()) {
      aux_fnames.push_back(it);
    }

    // sort and remove duplicates
    std::sort(aux_fnames.begin(),aux_fnames.end());
    auto it = std::unique(aux_fnames.begin(),aux_fnames.end());
    aux_fnames.erase(it,aux_fnames.end());
    return aux_fnames;
  }

protected:

  void set_field_impl (const Field& f) {
    m_inner->set_field(f);

    // If no aux fields are needed by the inner stat,
    // we can create the inner stat field, and set it as input
    // to the outer stat
    if (m_inner->get_aux_fields_names().size()==0) {
      m_inner->create_stat_field();
      m_outer->set_field(m_inner->get_stat_field());
    }
  }

  void set_aux_fields_impl () {
    m_inner->set_aux_fields(m_aux_fields);

    // We can let m_inner create the stat field, so we can finally start
    // to setup the outer stat
    m_inner->create_stat_field ();

    // Inner stat may have created some aux fields, so grab them
    for (const auto& it : m_inner->get_aux_fields()) {
      m_aux_fields[it.first] = it.second;
    }
    
    m_outer->set_field(m_inner->get_stat_field());
    m_outer->set_aux_fields(m_aux_fields);

    // Outer stat may have created some aux fields, so grab them
    for (const auto& it : m_outer->get_aux_fields()) {
      m_aux_fields[it.first] = it.second;
    }
  }

  void compute_impl () {
    m_inner->compute(m_timestamp);
    m_outer->compute(m_timestamp);
  }

  std::shared_ptr<FieldStat> m_inner;
  std::shared_ptr<FieldStat> m_outer;
};

} // namespace cldera

#endif // CLDERA_FIELD_STAT_PIPE_HPP

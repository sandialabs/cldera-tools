#include "cldera_graph_factory.hpp"

namespace cldera
{

std::shared_ptr<Graph> GraphFactory::build_graph(std::ostream& out) const
{
  using vos_t = std::vector<std::string>;
  
  if(m_verbose)
    out << "GraphFactory::buildGraph()\n" << "  Loading input file " << m_filename << "\n";

  std::map<std::string, std::shared_ptr<GraphVertex> > vertices;
  std::map<std::string, vos_t> edges;

  if(m_verbose) 
    out << "  Grabbing DAG sublist...\n";

  // TODO: no longer checks filename to see if it's valid. will still probably crash in the case of an invalid filename, but 
  // it may be more cryptic than previous

  const auto& dag_list = m_params.sublist("DAG");
  for (auto nameptr = dag_list.sublists_names_cbegin(); nameptr!=dag_list.sublists_names_cend(); ++nameptr)
  {
    std::string name = *nameptr;
    if(m_verbose)
      out << "  Processing species " << name << "...\n";

    const auto& species_list = dag_list.sublist(name);

    GraphVertex vertex(name);
    vertices[name] = std::make_shared<GraphVertex>(vertex);

    if(species_list.isParameter("Affects"))
    {
      // Two routes: Affects is a list, or Affects is a string
      if(species_list.isType<vos_t>("Affects"))
      {
        vos_t affected_species = species_list.get<vos_t>("Affects");
        edges[name] = affected_species;
        for (std::string& affected_name : affected_species)
          if(m_verbose)
            out << "  Species " << name << " affects " << affected_name << "\n";
      }
      else if(species_list.isType<std::string>("Affects"))
      {
        vos_t affected_species;
        affected_species.push_back(species_list.get<std::string>("Affects"));

        edges[name] = affected_species;
        for (std::string& affected_name : affected_species)
          if(m_verbose)
            out << "  Species " << name << " affects " << affected_name << "\n";
      }
    }
  }

  Graph graph(vertices, edges);
  return std::make_shared<Graph>(graph);
}

} // namespace cldera

#include "cldera_graph_factory.hpp"

namespace cldera
{

std::shared_ptr<Graph> GraphFactory::build_graph(std::ostream& out) const
{
  using vos_t = std::vector<std::string>;
  
  if(m_verbose)
    out << "GraphFactory::buildGraph()\n" << "  Loading input file " << m_filename << "\n";

  EKAT_REQUIRE_MSG (std::ifstream(m_filename).good(), "Error! GraphFactory: Filename of " + m_filename + " is invalid!\n");
  ekat::ParameterList params(ekat::parse_yaml_file(m_filename));

  std::map<std::string, std::shared_ptr<GraphVertex> > vertices;
  std::map<std::string, vos_t> edges;

  if(m_verbose) 
    out << "  Grabbing DAG sublist...\n";

  const auto& dag_list = params.sublist("DAG");
  for (auto nameptr = dag_list.sublists_names_cbegin(); nameptr!=dag_list.sublists_names_cend(); ++nameptr)
  {
    std::string name = *nameptr;
    if(m_verbose)
      out << "  Processing species " << name << "...\n";

    const auto& species_list = dag_list.sublist(name);

    double activation_threshold = 0.0;
    if(species_list.isParameter("Activation Threshold"))
      activation_threshold = species_list.get<double>("Activation Threshold");

    double deactivation_threshold = 0.0;
    if(species_list.isParameter("Deactivation Threshold")) 
      deactivation_threshold = species_list.get<double>("Deactivation Threshold");

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

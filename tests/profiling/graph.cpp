#include <catch2/catch.hpp>

#include "profiling/cldera_profiling_archive.hpp"
#include "profiling/cldera_graph.hpp"
#include "profiling/cldera_graph_vertex.hpp"

#include <ekat/ekat_parameter_list.hpp>
#include <ekat/io/ekat_yaml.hpp>
#include <ekat/util/ekat_string_utils.hpp>
#include <ekat/ekat_assert.hpp>

#include <fstream>
#include <memory>
#include <vector>

TEST_CASE ("graph") {

  using namespace cldera;

  using requests_t = std::map<std::string,std::vector<cldera::StatType>>;
  using vos_t = std::vector<std::string>;

  std::string filename = "./cldera_graph_input.yaml";
  if (std::ifstream(filename).good()) {
    ekat::ParameterList params(ekat::parse_yaml_file(filename));

    std::map<std::string, std::shared_ptr<cldera::GraphVertex> > vertices;
    std::map<std::string, std::vector<std::string> > edges;

    std::cout << "Grabbing DAG sublist..." << std::endl;
    const auto& dag_list = params.sublist("DAG");
    for (auto nameptr = dag_list.sublists_names_cbegin(); nameptr!=dag_list.sublists_names_cend(); ++nameptr)
    {
      std::string name = *nameptr;
      std::cout << "Processing species " << name << "..." << std::endl;

      const auto& species_list = dag_list.sublist(name);

      double activation_threshold = 0.0;
      if(species_list.isParameter("Activation Threshold"))
        activation_threshold = species_list.get<double>("Activation Threshold");

      double deactivation_threshold = 0.0;
      if(species_list.isParameter("Deactivation Threshold")) 
        deactivation_threshold = species_list.get<double>("Deactivation Threshold");

      cldera::GraphVertex vertex(name,activation_threshold,deactivation_threshold);
      vertices[name] = std::make_shared<cldera::GraphVertex>(vertex);

      if(species_list.isParameter("Affects"))
      {
        // Two routes: Affects is a list, or Affects is a string
        if(species_list.isType<std::vector<std::string> >("Affects"))
        {
          std::vector<std::string> affected_species = species_list.get<std::vector<std::string> >("Affects");
          edges[name] = affected_species;
          for (std::string& affected_name : affected_species)
            std::cout << "Species " << name << " affects " << affected_name << std::endl;
        }
        else if(species_list.isType<std::string>("Affects"))
        {
          std::vector<std::string> affected_species;
          affected_species.push_back(species_list.get<std::string>("Affects"));

          edges[name] = affected_species;
          for (std::string& affected_name : affected_species)
            std::cout << "Species " << name << " affects " << affected_name << std::endl;
        }
      }
    }

    cldera::Graph graph(vertices, edges);
    graph.generate_dot_graph();
    REQUIRE(graph.get_num_vertices() == 4);
    REQUIRE(graph.get_num_edges() == 3);
    
  } else {
    std::cout << " -> WARNING: no 'cldera_graph_input.yaml' file found." << std::endl;
    REQUIRE(false);
  }
}

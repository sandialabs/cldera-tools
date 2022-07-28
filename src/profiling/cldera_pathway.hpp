#ifndef CLDERA_PATHWAY_HPP
#define CLDERA_PATHWAY_HPP

#include "cldera_graph.hpp"
#include "cldera_field_test.hpp"

#include <iostream>
#include <vector>
#include <string>

namespace cldera
{

enum DetectionStatus {
  IGNORE = 0,
  WATCH = 1,
  DETECTED = 2
};

class Pathway {
public:
  // Take a Graph and a map from vertex names to a vector of FieldTests
  Pathway(const Graph& graph, const std::map<std::string, std::vector<std::shared_ptr<FieldTest> > >& tests);
  //Pathway(std::shared_ptr<Graph> graph, std::map<std::string, std::vector<std::shared_ptr<FieldTest> > >& tests);

  // Take a Graph and a map from vertex names to a vector of FieldTests
  //Pathway(Graph& graph, std::map<std::string, std::vector<std::shared_ptr<FieldTest> > >& tests, );

  // Return true if all tests pass. Return false and record values when any fail
  bool run_pathway_tests(const ekat::Comm& comm, const TimeStamp& time);

  // Draw the current pathway status to a DAG
  void generate_dot_graph(std::ostream& out);

  // Dump the pathway data from all tests to a yaml
  void dump_test_history_to_yaml(const std::string& file_name);

private:

  void initialize_detection_status();

  // The graph to run pathway tests on
  Graph m_graph;

  // A map from vertex names to lists of tests to be run by the pathway
  std::map<std::string, std::vector<std::shared_ptr<FieldTest> > > m_tests;

  // A map from vertex names to a DetectionStatus in the pathway
  std::map<std::string, DetectionStatus> m_detection_status;

  // A map from vertex names to detection history
  std::map<std::string, History> m_detection_data;
};

} // namespace cldera

#endif // CLDERA_PATHWAY_HPP

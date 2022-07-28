#include "cldera_pathway.hpp"

#include <ekat/io/ekat_yaml.hpp>
#include <ekat/ekat_parameter_list.hpp>
#include <ekat/ekat_assert.hpp>

namespace cldera
{

Pathway::Pathway(const Graph& graph, const std::map<std::string, std::vector<std::shared_ptr<FieldTest> > >& tests)
 : m_graph(graph),
   m_tests(tests)
{
  // Check that graph and the map are consistent with each other; i.e. check keys of "tests" are all vertices of "graph"
  for(auto test_pair : m_tests) {
    std::string key = test_pair.first;
    EKAT_REQUIRE_MSG (m_graph.has_vertex(key), "Error! Pathway constructor: Test vertex " + key + " does not exist on associated graph!\n");
  }

  initialize_detection_status();
}

void Pathway::initialize_detection_status()
{
  // TODO: allow more fine-grained control over the various pathway functionalities in this interface.

  // CASE 1) initialize all vertices with WATCH status
  for(auto test_pair : m_tests) {
    std::string key = test_pair.first;
    m_detection_status[key] = DetectionStatus::WATCH;
  }

  // CASE 2) initialize parentless nodes with WATCH status and others with IGNORE
  for(auto test_pair : m_tests) {
    std::string key = test_pair.first;

    std::vector<std::string> parents = m_graph.get_parents(key);
    if(parents.size()==0)
      m_detection_status[key] = DetectionStatus::WATCH;
    else
      m_detection_status[key] = DetectionStatus::IGNORE;
  }
}

bool Pathway::run_pathway_tests(const ekat::Comm& comm, const TimeStamp& time)
{
  bool tests_passed = true;

  std::vector<std::string> failed_tests;

  // loop over all key, vector<ptr<FieldTest> > pairs in the pathway
  for(auto test_pair : m_tests)
  {
    std::string key = test_pair.first;

    // if the variable in question is at an elevated status, grab the tests
    if(m_detection_status[key] >= DetectionStatus::WATCH)
    {
      std::vector<std::shared_ptr<FieldTest> > tests = test_pair.second;

      // loop over and run all the tests
      for(auto test : tests)
      {
        bool result = test->test(comm,{time.ymd,time.tod});

        // if it fails mark the result as failure, detect the node, and add children
        if(result == false)
        {
          tests_passed = false;
          m_detection_status[key] = DetectionStatus::DETECTED;
          failed_tests.push_back(key);
        }
      }
    }
  }

  // batch update to mark children of failed tests with "WATCH"
  for(std::string key : failed_tests)
    for(std::string child : m_graph.get_children(key))
      m_detection_status[child] = DetectionStatus::WATCH;

  return tests_passed;
}

void Pathway::generate_dot_graph(std::ostream& out)
{
  out << "digraph graphname {\n";

  // for each vertex in the graph, create a node with its detection status
  for(auto& pair : m_detection_status)
  {
    std::string name = pair.first;
    std::string node_color = "#000000"; // dot supports #RGB and #RGBA formats
    std::string font_color = "#000000";
    out << "  " << name << " [color=\"" << node_color << "\",fontcolor=\"" << font_color << "\"]\n";
  }

  // for each edge from first to second, create an edge in dot format: "first -> second"
  for(auto& pair : m_graph.get_edges())
  {
    for (const auto& dst : pair.second)
    {
      std::string arrow_color = "#000000";
      out << "  " << pair.first << " -> " << dst << " [color=\"" << arrow_color << "\"]\n";
    }
  }

  // output the data
  out << "  labelloc=\"t\"\n"
      << "  label=\"Pathway Diagram\\ntime=19910701.101000\"\n"
      << "}\n";
}

void Pathway::dump_test_history_to_yaml(const std::string& file_name)
{
  ekat::ParameterList pathway_list("Pathway tests");
  
  for(const auto& test_pair : m_tests)
  {
    std::string vertex_name = test_pair.first;
    auto& vertex_list = pathway_list.sublist(vertex_name);
    for(const auto& field_test : test_pair.second)
    {
      History test_history = field_test->get_test_history();
      vertex_list.set("Values",test_history.values());

      auto& times_vector = vertex_list.get<std::vector<std::string>>("Timestamps",{});
      const std::vector<TimeStamp> times = test_history.times();
      for (unsigned int i=0; i<times.size(); ++i)
      {
        times_vector.push_back(std::to_string(times[i].ymd)+"."+std::to_string(times[i].tod));
      }
    }
  }

  ekat::write_yaml_file(file_name,pathway_list);
}

} // namespace cldera

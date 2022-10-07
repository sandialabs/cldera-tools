#include "cldera_pathway_factory.hpp"

namespace cldera
{

std::shared_ptr<Pathway> PathwayFactory::build_pathway(std::ostream& out) const
{
  using vos_t = std::vector<std::string>;
  
  // Strategy:
  // 1. Build a graph based on the DAG in the input file
  // 2. Build field tests based on the parameters specified in the input file
  // 3. Link field tests to DAG vertices by reading the corresponding test names
  // 4. Call pathway(graph, tests)

  // Sanity checking and file parsing
  if(m_verbose)
    out << "PathwayFactory::buildPathway()\n" << "  Loading input file " << m_filename << "\n";
  EKAT_REQUIRE_MSG (std::ifstream(m_filename).good(), "Error! PathwayFactory: Filename of " + m_filename + " is invalid!\n");
  ekat::ParameterList params(ekat::parse_yaml_file(m_filename));

  // 1. Build the graph
  ekat::ParameterList pathway_sublist = params.sublist("Pathway");
  if(m_verbose)
    out << "PathwayFactory::build_pathway()\n" << "  Creating graph via GraphFactory\n";
  cldera::GraphFactory factory(pathway_sublist, m_verbose);
  std::shared_ptr<cldera::Graph> graph = factory.build_graph();

  // 2. Build the field tests
  if(m_verbose)
    out << "PathwayFactory::build_pathway()\n" << "  Creating field tests via FieldTestFactory\n";
  FieldTestFactory field_test_factory(m_filename, m_fields, m_comm);
  std::map<std::string, std::shared_ptr<FieldTest> > tests = field_test_factory.build_field_tests();

  // 3. Link field tests to DAG vertices (fill vertex_tests using entries in tests)
  if(m_verbose)
    out << "PathwayFactory::build_pathway()\n" << "  Linking test names to vertices\n";
  std::map<std::string, std::vector<std::shared_ptr<FieldTest> > > vertex_tests;
  ekat::ParameterList dag_sublist = pathway_sublist.sublist("DAG");
  std::vector<std::string> graph_vertices = graph->get_vertex_names();
  for(std::string &vertex_name : graph_vertices)
  {
    if(m_verbose)
      out << "Processing vertex " << vertex_name << "...\n";

    ekat::ParameterList vertex_sublist = dag_sublist.sublist(vertex_name);

    if(vertex_sublist.isParameter("Activation Tests"))
    {
      if(m_verbose)
        out << "Adding tests for " << vertex_name << "...\n";
      
      std::vector<std::string> activation_test_names = vertex_sublist.get<std::vector<std::string>>("Activation Tests");
      for(std::string &test_name : activation_test_names)
      {
        if(m_verbose)
          out << "Linking test " << test_name << " to vertex " << vertex_name << "...\n";

        vertex_tests[vertex_name].push_back(tests[test_name]);
      }
    }
  }

  // 4. Build the pathway
  if(m_verbose)
    out << "PathwayFactory::build_pathway()\n" << "  Building pathway\n";
  Pathway pathway(*graph, vertex_tests);
  return std::make_shared<Pathway>(pathway);
}

} // namespace cldera

#include "cldera_field_test_factory.hpp"

namespace cldera
{

std::map<std::string, std::shared_ptr<FieldTest> > FieldTestFactory::build_field_tests(std::ostream& out)
{
  using vos_t = std::vector<std::string>;

  std::map<std::string, std::shared_ptr<FieldTest> > tests;

  if(m_verbose)
    out << "FieldTestFactory::buildFieldTest()\n" << "  Loading input file " << m_filename << "\n";

  EKAT_REQUIRE_MSG (std::ifstream(m_filename).good(), "Error! FieldTestFactory: Filename of " + m_filename + " is invalid!\n");
  ekat::ParameterList params(ekat::parse_yaml_file(m_filename));

  if(m_verbose) 
    out << "  Grabbing Test sublist...\n";

  const auto& test_list = params.sublist("Tests");
  for (auto nameptr = test_list.sublists_names_cbegin(); nameptr!=test_list.sublists_names_cend(); ++nameptr)
  {
    std::string name = *nameptr; // test name
    
    // Get the sublist we're on
    ekat::ParameterList test_params = test_list.sublist(name);

    EKAT_REQUIRE_MSG(test_params.isParameter("Type"), "Error! FieldTestFactory: Test named " + name + " does not have a valid Type!\n");

    std::string test_type = test_params.get<std::string>("Type");
    std::string field_name = test_params.get<std::string>("Field");
    const auto& field = *m_fields[field_name];

    if(m_verbose)
      out << "  Processing test " << name << " of type " << test_type << " on field " << field_name << "...\n";

    auto& test = tests[name];

    // when adding new FieldTest children, add an entry here as well so the profiling tools can find it
    if(test_type == "Bounds")
    {
      std::vector<Real> params = test_params.get<std::vector<Real> >("Params");
      if(m_verbose)
        out << "    Min=" << params[0] << ", Max=" << params[1] << "\n";
      Bounds<Real> bounds(params);
      test = std::make_shared<BoundsFieldTest>(name,field,bounds,m_comm);
    }
    else if(test_type == "Min")
    {
      std::vector<Real> params = test_params.get<std::vector<Real> >("Params");
      if(m_verbose)
        out << "    Min=" << params[0] << "\n";
      test = std::make_shared<MinFieldTest>(name,field,params[0],m_comm);
    }
    else if(test_type == "Max")
    {
      std::vector<Real> params = test_params.get<std::vector<Real> >("Params");
      if(m_verbose)
        out << "    Max=" << params[0] << "\n";
      test = std::make_shared<MaxFieldTest>(name,field,params[0],m_comm);
    }
    else
    {
      EKAT_ERROR_MSG("FieldTestFactory: Invalid test type " + test_type + " encountered!");
    }
    // We save history for all tests
    test->set_save_history(true);
  }

  return tests;
}

} // namespace cldera

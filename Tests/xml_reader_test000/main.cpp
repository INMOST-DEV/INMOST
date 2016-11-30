
#include <cstdlib>
#include <inmost.h>
#include <istream>
#include <string>
#include <cassert>

using namespace INMOST;

const std::string xml_string = "<SolverParameters>\n"
        "            <!-- Links to files for external solvers, prefix processed inside -->\n"
        "            <PETSc File=\"petsc_options.txt\"/>\n"
        "            <Trilinos_Ifpack file=\"trilinos_ifpack.xml\"/>\n"
        "            <Trilinos_Belos file=\"trilinos_belos.xml\"/>\n"
        "            <Trilinos_ML file=\"trilinos_ml.xml\"/>\n"
        "            <Trilinos_AztecOO file=\"trilinos_aztec.xml\"/>\n"
        "            <!-- internal solvers parameters -->\n"
        "            <inner_ilu2>\n"
        "                        <!-- for prefix=\"pressure_solver\" -->\n"
        "                        <pressure_solver>\n"
        "                                   <absolute_tolerance value=\"1.0e-8\"/>\n"
        "                                   <relative_tolerance value=\"1.0e-4\"/>\n"
        "                                   <drop_tolerance value=\"1.0e-2\"/>\n"
        "                                   <reuse_tolerance value=\"1.0e-3\"/>\n"
        "                        </pressure_solver>\n"
        "                        <!-- for prefix=\"diffusion_solver\" -->\n"
        "                        <diffusion_solver> \n"
        "                                   <absolute_tolerance value=\"1.0e-6\"/>\n"
        "                                   <relative_tolerance value=\"1.0e-8\"/>\n"
        "                                   <drop_tolerance value=\"1.0e-2\"/>\n"
        "                                   <reuse_tolerance value=\"1.0e-3\"/>\n"
        "                        </diffusion_solver>\n"
        "                        <!-- for prefix=\"advection_solver\" -->\n"
        "            </inner_ilu2>\n"
        "            <inner_mptiluc2>\n"
        "                        <!-- link external files -->\n"
        "                        <pressure_solver Include=\"mptiluc2_pressure_solver.xml\"/>\n"
        "            </inner_mptiluc2>\n"
        "            <!-- link external files -->\n"
        "            <inner_mptilu2 Include=\"mptilu2.xml\"/>\n"
        "</SolverParameters>\n";

int main(int argc, char **argv) {
    std::stringstream input(xml_string);

    assert(!input.fail());

    XMLReader reader("", input);
    XMLReader::XMLTree tree = reader.ReadXML();

    //Check root
    assert(tree.tag.name == "SolverParameters");
    assert(tree.tag.attributes.size() == 0);
    assert(tree.children.size() == 8);

    //Check PETSc child
    XMLReader::XMLTree petsc_tree = tree.children[0];
    assert(petsc_tree.tag.name == "PETSc");
    assert(petsc_tree.tag.attributes.size() == 1);
    assert(petsc_tree.tag.attributes[0].name == "File");
    assert(petsc_tree.tag.attributes[0].value == "petsc_options.txt");

    //Check Trilinos_Ifpack child
    XMLReader::XMLTree Trilinos_Ifpack_tree = tree.children[1];
    assert(Trilinos_Ifpack_tree.tag.name == "Trilinos_Ifpack");
    assert(Trilinos_Ifpack_tree.tag.attributes.size() == 1);
    assert(Trilinos_Ifpack_tree.tag.attributes[0].name == "file");
    assert(Trilinos_Ifpack_tree.tag.attributes[0].value == "trilinos_ifpack.xml");

    //Check Trilinos_Belos child
    XMLReader::XMLTree Trilinos_Belos_tree = tree.children[2];
    assert(Trilinos_Belos_tree.tag.name == "Trilinos_Belos");
    assert(Trilinos_Belos_tree.tag.attributes.size() == 1);
    assert(Trilinos_Belos_tree.tag.attributes[0].name == "file");
    assert(Trilinos_Belos_tree.tag.attributes[0].value == "trilinos_belos.xml");

    //Check Trilinos_Belos child
    XMLReader::XMLTree Trilinos_ML_tree = tree.children[3];
    assert(Trilinos_ML_tree.tag.name == "Trilinos_ML");
    assert(Trilinos_ML_tree.tag.attributes.size() == 1);
    assert(Trilinos_ML_tree.tag.attributes[0].name == "file");
    assert(Trilinos_ML_tree.tag.attributes[0].value == "trilinos_ml.xml");

    //Check Trilinos_Belos child
    XMLReader::XMLTree Trilinos_AztecOO_tree = tree.children[4];
    assert(Trilinos_AztecOO_tree.tag.name == "Trilinos_AztecOO");
    assert(Trilinos_AztecOO_tree.tag.attributes.size() == 1);
    assert(Trilinos_AztecOO_tree.tag.attributes[0].name == "file");
    assert(Trilinos_AztecOO_tree.tag.attributes[0].value == "trilinos_aztec.xml");

    //Check inner_ilu2 child
    XMLReader::XMLTree inner_ilu2_tree = tree.children[5];
    assert(inner_ilu2_tree.tag.name == "inner_ilu2");
    assert(inner_ilu2_tree.tag.attributes.size() == 0);
    assert(inner_ilu2_tree.children.size() == 2);

    //Check inner_ilu2 pressure_solver
    XMLReader::XMLTree inner_ilu2_pressure_solver = inner_ilu2_tree.children[0];
    assert(inner_ilu2_pressure_solver.tag.name == "pressure_solver");
    assert(inner_ilu2_pressure_solver.tag.attributes.size() == 0);
    assert(inner_ilu2_pressure_solver.children.size() == 4);
    assert(inner_ilu2_pressure_solver.children[0].tag.name == "absolute_tolerance");
    assert(inner_ilu2_pressure_solver.children[0].tag.attributes.size() == 1);
    assert(inner_ilu2_pressure_solver.children[0].tag.attributes[0].name == "value");
    assert(inner_ilu2_pressure_solver.children[0].tag.attributes[0].value == "1.0e-8");
    assert(inner_ilu2_pressure_solver.children[1].tag.name == "relative_tolerance");
    assert(inner_ilu2_pressure_solver.children[1].tag.attributes.size() == 1);
    assert(inner_ilu2_pressure_solver.children[1].tag.attributes[0].name == "value");
    assert(inner_ilu2_pressure_solver.children[1].tag.attributes[0].value == "1.0e-4");
    assert(inner_ilu2_pressure_solver.children[2].tag.name == "drop_tolerance");
    assert(inner_ilu2_pressure_solver.children[2].tag.attributes.size() == 1);
    assert(inner_ilu2_pressure_solver.children[2].tag.attributes[0].name == "value");
    assert(inner_ilu2_pressure_solver.children[2].tag.attributes[0].value == "1.0e-2");
    assert(inner_ilu2_pressure_solver.children[3].tag.name == "reuse_tolerance");
    assert(inner_ilu2_pressure_solver.children[3].tag.attributes.size() == 1);
    assert(inner_ilu2_pressure_solver.children[3].tag.attributes[0].name == "value");
    assert(inner_ilu2_pressure_solver.children[3].tag.attributes[0].value == "1.0e-3");

    //Check inner_ilu2 diffusion_solver
    XMLReader::XMLTree inner_ilu2_diffusion_solver = inner_ilu2_tree.children[1];
    assert(inner_ilu2_diffusion_solver.tag.name == "diffusion_solver");
    assert(inner_ilu2_diffusion_solver.tag.attributes.size() == 0);
    assert(inner_ilu2_diffusion_solver.children.size() == 4);
    assert(inner_ilu2_diffusion_solver.children[0].tag.name == "absolute_tolerance");
    assert(inner_ilu2_diffusion_solver.children[0].tag.attributes.size() == 1);
    assert(inner_ilu2_diffusion_solver.children[0].tag.attributes[0].name == "value");
    assert(inner_ilu2_diffusion_solver.children[0].tag.attributes[0].value == "1.0e-6");
    assert(inner_ilu2_diffusion_solver.children[1].tag.name == "relative_tolerance");
    assert(inner_ilu2_diffusion_solver.children[1].tag.attributes.size() == 1);
    assert(inner_ilu2_diffusion_solver.children[1].tag.attributes[0].name == "value");
    assert(inner_ilu2_diffusion_solver.children[1].tag.attributes[0].value == "1.0e-8");
    assert(inner_ilu2_diffusion_solver.children[2].tag.name == "drop_tolerance");
    assert(inner_ilu2_diffusion_solver.children[2].tag.attributes.size() == 1);
    assert(inner_ilu2_diffusion_solver.children[2].tag.attributes[0].name == "value");
    assert(inner_ilu2_diffusion_solver.children[2].tag.attributes[0].value == "1.0e-2");
    assert(inner_ilu2_diffusion_solver.children[3].tag.name == "reuse_tolerance");
    assert(inner_ilu2_diffusion_solver.children[3].tag.attributes.size() == 1);
    assert(inner_ilu2_diffusion_solver.children[3].tag.attributes[0].name == "value");
    assert(inner_ilu2_diffusion_solver.children[3].tag.attributes[0].value == "1.0e-3");

    //Check inner_mptiluc2 child
    XMLReader::XMLTree inner_mptiluc2_tree = tree.children[6];
    assert(inner_mptiluc2_tree.tag.name == "inner_mptiluc2");
    assert(inner_mptiluc2_tree.tag.attributes.size() == 0);
    assert(inner_mptiluc2_tree.children.size() == 1);
    assert(inner_mptiluc2_tree.children[0].tag.name == "pressure_solver");
    assert(inner_mptiluc2_tree.children[0].tag.attributes.size() == 0);

    //Check inner_mptilu2 child
    XMLReader::XMLTree inner_mptilu2_tree = tree.children[7];
    assert(inner_mptilu2_tree.tag.name == "inner_mptilu2");
    assert(inner_mptilu2_tree.tag.attributes.size() == 0);
    assert(inner_mptilu2_tree.children.size() == 0);

    return 0;
} 
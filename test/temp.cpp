//#include "../watersystem.h"
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/parsers.hpp>
#include <iostream>

using namespace std;
namespace po = boost::program_options;

int main (int argc, char **argv) {

// Declare the supported options.
po::options_description desc("Allowed options");
desc.add_options()
    ("help", "produce help message")
    ("file", po::value<char*>(), "choose a file")
;

po::variables_map vm;
po::store(po::parse_command_line(argc, argv, desc), vm);
po::notify(vm);

if (vm.count("help")) {
    cout << desc << "\n";
    return 1;
}

if (vm.count("file")) {
    cout << "File loaded is "
 << vm["file"].as<char*>() << ".\n";
} else {
    cout << "File was not give\n";
}

return 0;
}

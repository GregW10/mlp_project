#include <iostream>
#include <fstream>
#include "../glib/misc/gregparse.hpp"

int main(int argc, char **argv) {
    gtd::parser parser{argc, argv};
    const char *arg = parser.get_arg(std::regex{R"(^r-p?v?a?c?e?$)"});
    std::cout << (arg ? arg : "nullptr") << std::endl;
    return 0;
}

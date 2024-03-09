#include "../glib/ml/gregffnn.hpp"
#include "../glib/misc/gregparse.hpp"
#include "nnsup.hpp"
#include <algorithm>

#define DEF_NUM_PASSES ((uint64_t) 1'000)

int main(int argc, char **argv) {
    gtd::parser parser{argc, argv};
    const char *data_dir = parser.get_arg("--data_dir");
    if (!data_dir)
        data_dir = ".";
    bool pp = parser.get_arg("--preprocess", false);
    uint64_t num_passes = parser.get_arg("--num_passes", DEF_NUM_PASSES);
    if (!num_passes) {
        std::cerr << "Error: number of passes over entire dataset cannot be zero.\n";
        return 1;
    }
    std::vector<std::string> files;
    struct dirent *entry;
    DIR *dir;
    if (!(dir = opendir(data_dir))) {
        std::cerr << "Error: could not open directory \"" << data_dir << "\".\n";
        return 1;
    }
    while ((entry = readdir(dir)))
        if (gtd::endswith(entry->d_name, ".3bod"))
            files.emplace_back(entry->d_name);
    closedir(dir);
    if (files.empty()) {
        std::cerr << "Error: no .3bod files found in specified directory \"" << data_dir << "\".\n";
        return 1;
    }
    uint64_t pass = 0;
    std::mt19937_64 rng{std::random_device{}()};
    while (pass++ < num_passes) {
        std::shuffle(files.begin(), files.end(), rng); // shuffle .3bod files so pass over them isn't in the same order
    }
    return 0;
}

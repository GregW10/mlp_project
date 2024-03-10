#include "../glib/ml/gregffnn.hpp"
#include "../glib/misc/gregparse.hpp"
#include "nnsup.hpp"
#include <algorithm>

#define DEF_NUM_PASSES ((uint64_t) 1'000)

[[noreturn]] void pp_error() {
    std::cerr << "Error: only one of \"--just_preprocess\", \"--preprocess\" or \"--not_preprocessed\" can be "
                 "specified.\n";
    exit(1);
}

template <typename T> requires (std::is_floating_point_v<T>)
std::pair<void (*)(T&), void (*)(T&)> afunc_selector(const char *arg) {
    /* Function to select appropriate activation function and its derivative from command-line argument `arg`. */
    if (!arg)
        throw std::invalid_argument{"Error: activation function description cannot be nullptr.\n"};
    if (gtd::str_eq(arg, "sigmoid"))
        return {gml::activations::sigmoid<T>, gml::activations::sigmoid_d<T>};
    if (gtd::str_eq(arg, "softsign"))
        return {gml::activations::softsign<T>, gml::activations::softsign_d<T>};
    std::cerr << "Error: invalid activation function \"" << arg << "\" passed.\n";
    exit(1);
}

template <typename dtype, typename wtype> requires (std::is_floating_point_v<dtype> && std::is_floating_point_v<wtype>)
void run_sim() {
    /* Runs the entire simulation. `dtype` is the data type of the data and `wtype` is the type to use in the NN. */
}

int main(int argc, char **argv) {
    gtd::parser parser{argc, argv};
    const char *data_dir = parser.get_arg("--data_dir");
    if (!data_dir)
        data_dir = ".";
    /* There are 4 possibilities for running this program, selection of which is controlled by the 3 booleans below.
     * If:
     *     - All 3 booleans are false, then the program interprets the data in `data_dir` as being preprocessed and runs
     *       the neural net directly on these (thus looks for .3bodpp files)
     *     - `jpp` is true, then the program just performs preprocessing on the allegedly unprocessed data in `data_dir`
     *       and saves it all to a new folder
     *     - `pp` is true, then the allegedly unprocessed data in `data_dir` is preprocessed and written all to a new
     *       folder, and the neural network is run on this freshly preprocessed data
     *     - `npp` is true, then the neural network runs directly on the allegedly unprocessed data in `data_dir` (and
     *       thus looks for .3bod files)
     * Any other combination is invalid. */
    bool jpp = parser.get_arg("--just_preprocess", false); // only preprocess data in `data_dir`, don't run neural net
    bool pp = parser.get_arg("--preprocess", false); // preprocess data in `data_dir` and run neural net
    bool npp = parser.get_arg("--not_preprocessed", false); // data in `data_dir` isn't preprocessed and shouldn't be
    if (jpp) {
        if (pp || npp)
            pp_error();

    }
    const char *afunc_arg = parser.get_arg("--activation");
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

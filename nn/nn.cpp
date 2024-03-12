#include "../glib/ml/gregffnn.hpp"
#include "../glib/misc/gregparse.hpp"
#include "../glib/misc/gregmisc.hpp"
#include "nnsup.hpp"
#include <algorithm>

#define NORM_RGX std::regex{R"(^--norm=(m|l)?t?p?v?$)"} // improve this regex
#define LAYERS_RGX std::regex{R"(^--layers=\d{1,17}-\d{1,17}:\d{1,8}(,\d{1,17}-\d{1,17}:\d{1,8})*$)"}

#define DEF_NUM_PASSES ((uint64_t) 1'000)

template <typename dtype, typename wtype> requires (std::is_floating_point_v<dtype> && std::is_floating_point_v<wtype>)
void run_sim() {
    /* Runs the entire simulation. `dtype` is the data type of the data and `wtype` is the type to use in the NN. */
}

template <typename T> requires (std::is_floating_point_v<T>)
gtd::normaliser<T> *get_normaliser(const char *norm_arg) {
    // logic...
    if (!norm_arg)
        return nullptr;
    norm_arg += 7;
    if (!*norm_arg)
        return nullptr;
    bool _m = false;
    bool _l = false;
    bool _t = false;
    bool _p = false;
    bool _v = false;
    do {
        if (*norm_arg == 'm')
            _m = true;
        if (*norm_arg == 'l')
            _l = true;
        if (*norm_arg == 't')
            _t = true;
        if (*norm_arg == 'p')
            _p = true;
        if (*norm_arg == 'v')
            _v = true;
    } while (*++norm_arg);
    if (_m) {
        if (_t) {
            if (_p) {
                if (_v)
                    return new gtd::gen_normaliser<T, gtd::mass_normaliser<T>, gtd::time_normaliser<T>, gtd::pos_normaliser<T>, gtd::vel_normaliser<T>>{};
                else
                    return new gtd::gen_normaliser<T, gtd::mass_normaliser<T>, gtd::time_normaliser<T>, gtd::pos_normaliser<T>>{};
            }
            else {
                if (_v)
                    return new gtd::gen_normaliser<T, gtd::mass_normaliser<T>, gtd::time_normaliser<T>, gtd::vel_normaliser<T>>{};
                else
                    return new gtd::gen_normaliser<T, gtd::mass_normaliser<T>, gtd::time_normaliser<T>>{};
            }
        }
        else {
            if (_p) {
                if (_v)
                    return new gtd::gen_normaliser<T, gtd::mass_normaliser<T>, gtd::pos_normaliser<T>, gtd::vel_normaliser<T>>{};
                else
                    return new gtd::gen_normaliser<T, gtd::mass_normaliser<T>, gtd::pos_normaliser<T>>{};
            }
            else {
                if (_v)
                    return new gtd::gen_normaliser<T, gtd::mass_normaliser<T>, gtd::vel_normaliser<T>>{};
                else
                    return new gtd::gen_normaliser<T, gtd::mass_normaliser<T>>{};
            }
        }
    } else {
        if (_l) { // `_m` and `_l` are mutually exclusive (this was already checked by the regex)
            if (_t) {
                if (_p) {
                    if (_v)
                        return new gtd::gen_normaliser<T, gtd::log_mass_normaliser<T>, gtd::time_normaliser<T>, gtd::pos_normaliser<T>, gtd::vel_normaliser<T>>{};
                    else
                        return new gtd::gen_normaliser<T, gtd::log_mass_normaliser<T>, gtd::time_normaliser<T>, gtd::pos_normaliser<T>>{};
                }
                else {
                    if (_v)
                        return new gtd::gen_normaliser<T, gtd::log_mass_normaliser<T>, gtd::time_normaliser<T>, gtd::vel_normaliser<T>>{};
                    else
                        return new gtd::gen_normaliser<T, gtd::log_mass_normaliser<T>, gtd::time_normaliser<T>>{};
                }
            }
            else {
                if (_p) {
                    if (_v)
                        return new gtd::gen_normaliser<T, gtd::log_mass_normaliser<T>, gtd::pos_normaliser<T>, gtd::vel_normaliser<T>>{};
                    else
                        return new gtd::gen_normaliser<T, gtd::log_mass_normaliser<T>, gtd::pos_normaliser<T>>{};
                }
                else {
                    if (_v)
                        return new gtd::gen_normaliser<T, gtd::log_mass_normaliser<T>, gtd::vel_normaliser<T>>{};
                    else
                        return new gtd::gen_normaliser<T, gtd::log_mass_normaliser<T>>{};
                }
            }
        }
        else {
            if (_t) {
                if (_p) {
                    if (_v)
                        return new gtd::gen_normaliser<T, gtd::time_normaliser<T>, gtd::pos_normaliser<T>, gtd::vel_normaliser<T>>{};
                    else
                        return new gtd::gen_normaliser<T, gtd::time_normaliser<T>, gtd::pos_normaliser<T>>{};
                }
                else {
                    if (_v)
                        return new gtd::gen_normaliser<T, gtd::time_normaliser<T>, gtd::vel_normaliser<T>>{};
                    else
                        return new gtd::gen_normaliser<T, gtd::time_normaliser<T>>{};
                }
            }
            else {
                if (_p) {
                    if (_v)
                        return new gtd::gen_normaliser<T, gtd::pos_normaliser<T>, gtd::vel_normaliser<T>>{};
                    else
                        return new gtd::gen_normaliser<T, gtd::pos_normaliser<T>>{};
                }
                else {
                    // if (_v)
                        return new gtd::gen_normaliser<T, gtd::vel_normaliser<T>>{};
                    // else
                    //     return new gtd::gen_normaliser<T>{}; // this is not an option, regex anyway assures this
                }
            }
        }
    }
    // return new gtd::mass_normaliser<T>{}; // placeholder until logic above is implemented
}

template <typename T> requires (std::unsigned_integral<T>)
T to_unsigned(const char **str) {
    if (!*str || !**str)
        return -1;
    T val = 0;
    const char *original = *str;
    while (**str) {
        if (**str < 48 || **str > 57) {
            if (original == *str)
                return -1;
            return val;
        }
        val *= 10;
        val += *(*str)++ - 48;
    }
    return val;
}

template <typename W>
void build_layers(gml::ffnn<W> &_nn, const char *layers) {
    if (!layers)
        return;
    layers += 9; // move pointer to point to first digit
    uint64_t idim;
    uint64_t odim;
    uint32_t fid;
    std::pair<void (*)(T&), void (*)(T&)> funcs;
    while (true) {
        idim = to_unsigned<uint64_t>(&layers);
        ++layers; // left pointing to hyphen so move forward by one to next digit
        odim = to_unsigned<uint64_t>(&layers);
        ++layers; // left pointing to colon so move forward by one to next digit
        fid = to_unsigned<uint32_t>(&layers);
        funcs = gml::activations::get_func_by_id<W>(fid);
        _nn.emplace_back(idim, odim, funcs.first, funcs.second); // later add in support for different W initialisations
        if (*layers++ == ',')
            continue;
        else
            break;
    }
}

template <typename T>
void entry_to_vector(const typename gtd::f3bodr<T>::entry_type &_e, gml::vector<T> *_v) {
    // figure out how to know at which offset to start (input vs output vector)
}

template <typename D, typename W>
void build_and_run_nn(const std::vector<std::string> &files, const char *layer_str, uint64_t num_passes,
                      uint64_t pairs_per_file, uint64_t batch_size, long double _lr,
                      const char *from_model, const std::vector<std::string> &val_files, bool alloc = true) {
    if (files.empty())
        return;
    typename std::vector<std::string>::size_type train_size = files.size();
    if (batch_size > train_size)
        throw std::invalid_argument{"Error: batch size cannot be greater than number of examples in training set.\n"};
    if (!batch_size)
        throw std::invalid_argument{"Error: batch size cannot be zero.\n"};
    if (layer_str && from_model)
        throw std::invalid_argument{"Error: either the model can be built from scratch, or loaded from a .nnw file, "
                                    "but not both!\n"};
    uint64_t tcounter = 0;
    uint64_t bcounter;
    uint64_t pcounter;
    std::mt19937_64 rng{std::random_device{}()};
    std::uniform_int_distribution<uint64_t> dist{};
    const typename gtd::f3bodr<D>::hdr_t *header{};
    uint64_t idx_lo;
    uint64_t idx_hi;
    uint64_t idx_temp;
    gml::vector<D> _input{22};
    gml::vector<D> *_y;
    gml::vector<W> _f{18};
    W _loss;
    if (alloc) {
        std::vector<gtd::f3bodr<D>> readers{};
        readers.reserve(files.size());
        for (const std::string &_str : files)
            readers.emplace_back(_str.c_str()); // exception will occur here in case of T mismatch
        gml::ffnn<W> net;
        if (layer_str)
            build_layers(net, layer_str);
        else
            net.load_model(from_model);
        while (tcounter++ < num_passes) {
            std::shuffle(readers.begin(), readers.end(), rng); // to make sure data is never presented in same order
            bcounter = 0;
            for (const gtd::f3bodr<D> &_rdr : readers) {
                header = &_rdr.header();
                dist.param({0, header->N});
                _input[0] = header->masses[0];
                _input[1] = header->masses[1];
                _input[2] = header->masses[2];
                pcounter = 0;
                while (pcounter++ < pairs_per_file) {
                    rand_idx:
                    idx_lo = dist(rng);
                    idx_hi = dist(rng);
                    if (idx_lo == idx_hi)
                        goto rand_idx;
                    if (idx_lo > idx_hi) {
                        idx_temp = idx_lo;
                        idx_lo = idx_hi;
                        idx_hi = idx_temp;
                    }
                    // forward and backwards passes
                    // entry_to_vector(_rdr.entry_at(idx_lo), &_input);
                    gml::gen::copy(_input.begin() + 3, _rdr.entry_at(idx_lo).positions, 18);
                    _input[21] = (idx_hi - idx_lo)*header->dt; // time-step between input and output
                    _f = &net.forward_pass(_input);
                    // entry_to_vector(_rdr.entry_at(idx_hi), &_y);
                    gml::gen::copy(_y.begin(), _rdr.entry_at(idx_hi).positions, 18);
                    _loss = net.backward_pass(_y);
                    if (++bcounter == batch_size) {
                        net.update_params(_lr);
                        bcounter = 0;
                    }
                }
            }
        }
        /* while (tcounter++ < train_size) {
            bcounter = 0;
            while (bcounter++ < batch_size && tcounter++ < train_size) {

            }
        } */
    }
}

template <typename T> requires (std::is_floating_point_v<T>)
std::vector<std::string> *pproc(const char *norm_str, gtd::normaliser<T> **ptr, bool jpp) {
    if (!ptr)
        throw std::invalid_argument{"Error: pointer to gtd::normaliser<T> pointer cannot be nullptr.\n"};
    *ptr = norm_str ? get_normaliser<T>(norm_str) :
            new gtd::gen_normaliser<T, gtd::mass_normaliser<T>, gtd::time_normaliser<T>>{};
    if (!*ptr)
        throw std::invalid_argument{"Empty normalisation string provided.\n"};
    std::pair<std::string, std::vector<std::string>*> _pair = gtd::preprocess(data_dir, **ptr, jpp);
    delete *ptr;
    std::cout << "Normalisation values written to: \"" << _pair.first << "\"\n";
    return _pair.second;
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
    if (jpp + pp + npp > 1) {
        std::cerr << "Error: only one of \"--just_preprocess\", \"--preprocess\" or \"--not_preprocessed\" can be "
                     "specified.\n";
        return 1;
    }
    const char *ftype = parser.get_arg("--fp_type"); // floating point type of weights and biases in NN
    const char *norm_str = parser.get_arg(NORM_RGX);
    const char *layers = parser.get_arg(LAYERS_RGX);
    if (!layers)
        layers = "--layers=22-128:1,128-64:1,64-18:0";
    if (!parser.empty()) {
        std::cerr << "Error: unrecognised arguments:\n";
        for (const char &*arg : parser) {
            std::cerr << arg << '\n';
        }
        return 1;
    }
    if (jpp) {
        if (ftype) {
            std::cerr << "Error: floating point type can only be specified for when the NN will be run, which is not "
                         "the case if the \"--just_preprocess\" flag is passed.\n";
            return 1;
        }
        gtd::normaliser<long double> *nptr_ld{};
        try {
            pproc<long double>(norm_str, &nptr_ld, true);
        } catch (const std::invalid_argument &e) {
            std::cerr << "Error: " << e.what() << '\n';
            return 1;
        } catch (const std::exception&) {
            delete [] *nptr_ld;
            gtd::normaliser<double> *nptr_d{};
            try {
                pproc<double>(norm_str, &nptr_d, true);
            } catch (const std::exception&) {
                delete [] *nptr_d;
                gtd::normaliser<float> *nptr_f{};
                try {
                    pproc<float>(norm_str, &nptr_f, true);
                } catch (const std::exception &e) {
                    delete [] *nptr_f;
                    std::cerr << "Error: " << e.what() << '\n';
                    return 1;
                }
            }
        }
        return 0;
    } // after this block it is certain the NN will be run
    bool nn_ld = false;
    bool nn_d = false;
    bool nn_f = false;
    if (!ftype)
        nn_ld = ftype = "long double";
    else {
        if (gtd::str_eq(ftype, "long double"))
            nn_ld = true;
        else if (gtd::str_eq(ftype, "double"))
            nn_d = true;
        else if (gtd::str_eq(ftype, "float"))
            nn_f = true;
        else {
            std::cerr << "Error: unrecognised floating point type specified \"" << ftype << "\".\n";
            return 1;
        }
    }
    std::vector<std::string> *_files{};
    bool d_ld = false;
    bool d_d = false;
    bool d_f = false;
    if (pp) {
        gtd::normaliser<long double> *nptr_ld{};
        try {
            _files = pproc<long double>(norm_str, &nptr_ld, false);
            d_ld = true;
        } catch (const std::invalid_argument &e) {
            std::cerr << "Error: " << e.what() << '\n';
            return 1;
        } catch (const std::exception&) {
            delete [] *nptr_ld;
            gtd::normaliser<double> *nptr_d{};
            try {
                _files = pproc<double>(norm_str, &nptr_d, false);
                d_d = true;
            } catch (const std::exception&) {
                delete [] *nptr_d;
                gtd::normaliser<float> *nptr_f{};
                try {
                    _files = pproc<float>(norm_str, &nptr_f, false);
                    d_f = true;
                } catch (const std::exception &e) {
                    delete [] *nptr_f;
                    std::cerr << "Error: " << e.what() << '\n';
                    return 1;
                }
            }
        }
    }
    else if (npp) {
        _files = gtd::find_files(data_dir, ".3bod");
        if (!_files) {
            std::cerr << "Error: could not find any .3bod files in directory \"" << data_dir << "\".\n";
            return 1;
        }
    } else {
        _files = gtd::find_files(data_dir, ".3bodpp");
        if (!_files) {
            std::cerr << "Error: could not find any .3bodpp files in directory \"" << data_dir << "\".\n";
            return 1;
        }
    }
    /* uint64_t num_passes = parser.get_arg("--num_passes", DEF_NUM_PASSES);
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
    } */
    delete [] _files;
    return 0;
}

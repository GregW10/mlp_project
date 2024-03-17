#include <algorithm>
#include "../glib/ml/gregffnn.hpp"
#include "../glib/misc/gregparse.hpp"
#include "../glib/misc/gregmisc.hpp"
#include "nnsup.hpp"
#include <csignal>
#include <atomic>

#define NORM_RGX std::regex{R"(^--norm=(m|l)?t?p?v?$)"} // improve this regex
#define LAYERS_RGX std::regex{R"(^--layers=\d{1,17}-\d{1,17}:\d{1,8}(,\d{1,17}-\d{1,17}:\d{1,8})*$)"}

#define DEF_NUM_PASSES ((uint64_t) 1'000)
#define DEF_LR (1/powl(2, 10))
#define DEF_PPF ((uint64_t) 1'000)
#define DEF_REC_FREQ ((uint64_t) 1)

// std::atomic<bool> received_signal{false};
std::atomic<int> signal_number{};

void signal_handler(int signal) {
    // received_signal.store(true);
    signal_number.store(signal);
}

uint64_t lossrf;

std::ofstream logger{}; // global so I don't have to worry about passing it to functions

template <typename T> requires (std::is_floating_point_v<T>)
gtd::normaliser<T> *get_normaliser(const char *norm_arg) {
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
    std::pair<void (*)(W&), void (*)(W&)> funcs;
    while (true) {
        idim = to_unsigned<uint64_t>(&layers);
        ++layers; // left pointing to hyphen so move forward by one to next digit
        odim = to_unsigned<uint64_t>(&layers);
        ++layers; // left pointing to colon so move forward by one to next digit
        fid = to_unsigned<uint32_t>(&layers);
        funcs = gml::activations::get_func_by_id<W>(fid);
        // Must add in support later for different initialisations:
        _nn.emplace_back(idim, odim, funcs.first, funcs.second, gml::GLOROT_UNIFORM);
        if (*layers++ == ',')
            continue;
        else
            break;
    }
}

// namespace gtd {
//     template <>
//     void swap(gtd::f3bodr<long double> &f1, gtd::f3bodr<long double> &f2) {
//         gml::gen::memswap(&f1, &f2);
//     }
//     template <>
//     void swap(gtd::f3bodr<double> &f1, gtd::f3bodr<double> &f2) {
//         gml::gen::memswap(&f1, &f2);
//     }
//     template <>
//     void swap(gtd::f3bodr<float> &f1, gtd::f3bodr<float> &f2) {
//         gml::gen::memswap(&f1, &f2);
//     }
// }

template <typename D, typename W>
void build_and_run_nn(std::unique_ptr<std::vector<std::string>> &train_files,
                      std::unique_ptr<std::vector<std::string>> &val_files,
                      const char *layer_str,
                      const char *from_model,
                      uint64_t num_passes,
                      uint64_t pairs_per_file,
                      uint64_t batch_size,
                      W _lr,
                      const char *output_nnw,
                      bool alloc = true) {
    if (!train_files)
        throw std::invalid_argument{"Error: pointer to training files std::vector cannot be nullptr.\n"};
    if (train_files->empty())
        throw std::invalid_argument{"Error: no training files found.\n"};
    typename std::vector<std::string>::size_type num_train_files = train_files->size();
    // if (batch_size > train_size)
    //     throw std::invalid_argument{"Error: batch size cannot be greater than number of examples in training set.\n"};
    uint64_t tot_examples;
    if (batch_size > (tot_examples = num_passes*num_train_files*pairs_per_file))
        throw std::invalid_argument{"Error: batch size cannot be greater than the total number of training examples "
                                    "that will be seen.\n"};
    if (!batch_size)
        throw std::invalid_argument{"Error: batch size cannot be zero.\n"};
    if (layer_str && from_model)
        throw std::invalid_argument{"Error: either the model can be built from scratch, or loaded from a .nnw file, "
                                    "but not both!\n"};
    if (lossrf > tot_examples/batch_size) // check this
        throw std::logic_error{"Error: the loss recording frequency cannot be greater than the total number of weight "
                               "updates foreseen to occur.\n"};
    std::chrono::time_point<std::chrono::system_clock> _start;
    std::chrono::time_point<std::chrono::system_clock> _end;
    uint64_t tcounter = 0;
    uint64_t bcounter = 0;
    uint64_t pcounter;
    uint64_t fcounter;
    std::mt19937_64 rng{std::random_device{}()};
    std::uniform_int_distribution<uint64_t> dist{};
    const typename gtd::f3bodr<D>::hdr_t *header{};
    uint64_t idx_lo;
    uint64_t idx_hi;
    uint64_t idx_temp;
    gml::vector<D> _input{(uint64_t) 22};
    gml::vector<D> _y{(uint64_t) 18};
    const gml::vector<W> *_f;
    W _tloss; // loss per training example
    W _mloss; // mean loss across batch
    gml::ffnn<W> net;
    if (layer_str)
        build_layers(net, layer_str);
    else
        net.load_model(from_model);
    logger << ".nnw file size = " << net.nnw_fsize() << " bytes\nNN f.p. type: ";
    if constexpr (std::same_as<W, long double>)
        logger << "\"long double\"\n";
    else if constexpr (std::same_as<W, double>)
        logger << "\"double\"\n";
    else
        logger << "\"float\"\n";
    logger << "Size of NN f.p. type = " << sizeof(W) << " bytes\nData f.p. type: ";
    if constexpr (std::same_as<D, long double>)
        logger << "\"long double\"\n";
    else if constexpr (std::same_as<D, double>)
        logger << "\"double\"\n";
    else
        logger << "\"float\"\n";
    logger << "Size of data f.p. type = " << sizeof(D) << " bytes\n";
    char *_ltime{};
    uint64_t rf_counter = 0;
    std::ofstream tlosses{"train_losses.csv", std::ios_base::out | std::ios_base::trunc};
    if (!tlosses.good())
        throw std::ios_base::failure{"Error: could not open \"train_losses.csv\" file.\n"};
    tlosses << "epoch,file_idx,loss\r\n";
    // make separate branch for val if `val_files` and only there open file
    if (alloc) {
        std::vector<gtd::f3bodr<D>> readers{};
        readers.reserve(num_train_files);
        for (const std::string &_str : *train_files)
            readers.emplace_back(_str.c_str()); // exception will occur here in case of T mismatch
        // delete [] train_files;
        train_files.reset();
        _ltime = gml::gen::now_str();
        logger << "Training start time: " << _ltime << '\n' << std::flush;
        delete [] _ltime;
        _start = std::chrono::system_clock::now();
        // think about recording epoch 0 loss (have to add stuff in `gregffnn.hpp`)
        while (tcounter++ < num_passes) {
            std::shuffle(readers.begin(), readers.end(), rng); // to make sure data is never presented in same order
            // bcounter = 0;
            fcounter = 0;
            for (const gtd::f3bodr<D> &_rdr : readers) {
                header = &_rdr.header();
                if (header->N <= 1) // case for empty or single-entry .3bod/.3bodpp file
                    continue;
                dist.param(std::uniform_int_distribution<uint64_t>::param_type{0, header->N - 1});
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
                    gml::gen::copy(_input.begin() + 3, _rdr.entry_at(idx_lo).positions, 18);
                    _input[21] = (idx_hi - idx_lo)*header->dt; // time-step between input and output
                    _f = &net.forward_pass(_input);
                    gml::gen::copy(_y.begin(), _rdr.entry_at(idx_hi).positions, 18);
                    _tloss = net.backward_pass(_y);
                    // std::cerr << "T loss: " << _tloss << std::endl;
                    ++fcounter;
                    if (++bcounter == batch_size) {
                        _mloss = net.update_params(_lr);
                        bcounter = 0;
                        std::cout << "Loss: " << _mloss << std::endl;
                        if (++rf_counter == lossrf) {
                            tlosses << tcounter << ',' << fcounter << ',' << _mloss << "\r\n";
                        }
                    }
                    if (signal_number.load()) {
                        std::cout << "Received signal ";
                        switch (signal_number.load()) {
                            case SIGINT:
                                std::cout << "SIGINT";
                                break;
                            case SIGTERM:
                                std::cout << "SIGTERM";
                                break;
                            case SIGHUP:
                                std::cout << "SIGHUP";
                                break;
                            case SIGQUIT:
                                std::cout << "SIGQUIT";
                        }
                        std::cout << ".\nStopping training and saving weights...\n";
                        goto _out;
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
    if (bcounter) {
        _mloss = net.update_params(_lr);
        std::cout << "Loss: " << _mloss << std::endl;
    }
    _out:
    _end = std::chrono::system_clock::now();
    _ltime = gml::gen::now_str();
    logger << "Training end time: " << _ltime << '\n';
    delete [] _ltime;
    tlosses.close();
    logger << "Actual num. passes over training set = " << tcounter << '\n';
    try {
        if (output_nnw) {
            std::cout << "Neural Network weights written to \"" << net.to_nnw(output_nnw) << "\"." << std::endl;
            logger << "NN weights filepath: \"" << output_nnw << "\"\n";
        }
        else {
            char *_ptr;
            std::cout << "Neural Network weights written to \"" << (_ptr = net.to_nnw()) << "\"." << std::endl;
            logger << "NN weights filepath: \"" << _ptr << "\"\n";
            /* Does not result in unfreed memory in case of an exception, as I have used `std::unique_ptr within my
             * `.to_nnw()` function, such that its destructor would be called during stack unwinding, as the exception
             * would be caught below. If no exception occurs, the `std::unique_ptr` releases the pointer. */
            delete [] _ptr;
        }
    } catch (const std::exception &_e) {
        std::cerr << "Error: weights were not written due to an exception occurring.\nwhat(): " << _e.what() << '\n';
    }
    std::cout << "Total training time took "
              << std::chrono::duration_cast<std::chrono::nanoseconds>(_end - _start).count()/BILLION << " seconds."
              << std::endl;
    logger.close();
}

int main(int argc, char **argv) {
    if (signal(SIGINT, signal_handler) == SIG_ERR) {
        std::cerr << "Error: could not register signal_handler() with signal() for SIGINT.\n";
        return 1;
    }
    if (signal(SIGTERM, signal_handler) == SIG_ERR) {
        std::cerr << "Error: could not register signal_handler() with signal() for SIGTERM.\n";
        return 1;
    }
    if (signal(SIGHUP, signal_handler) == SIG_ERR) {
        std::cerr << "Error: could not register signal_handler() with signal() for SIGHUP.\n";
        return 1;
    }
    if (signal(SIGQUIT, signal_handler) == SIG_ERR) {
        std::cerr << "Error: could not register signal_handler() with signal() for SIGQUIT.\n";
        return 1;
    }
    printf("PID: %d\n----------------\n", getpid());
    gtd::parser parser{argc, argv};
    bool _print = parser.get_arg("--print", false);
    if (_print) {
        printf("Experiment started with arguments:\n");
        for (const auto &[_, arg] : parser)
            printf("\t%s\n", arg.c_str());
    }
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
    const char *norm_str = parser.get_arg(NORM_RGX);
    if (norm_str && (npp || (!jpp && !pp))) {
        std::cerr << "Error: normalisation type cannot be specified when normalisation will not take place.\n";
        return 1;
    }
    const char *val_dir = parser.get_arg("--val_dir");
    // char *normv_tpath{};
    // char *normv_vpath{};
    std::unique_ptr<std::vector<std::string>> dppf{};
    std::unique_ptr<std::vector<std::string>> vppf{};
    if (jpp) {
        if (!parser.empty()) {
            std::cerr << "Error: unrecognised arguments:\n";
            for (const auto &[_, arg] : parser)
                fprintf(stderr, "\t%s\n", arg.c_str());
            fprintf(stderr, "These arguments could be valid for when the NN will be run, but the \"--just_preprocess\" "
                            "means the NN will not be run.\n");
            return 1;
        }
        try {
            try { // my first time using `std::unique_ptr<>` !!
#define PP_BLOCK(type, bval) \
                std::unique_ptr<gtd::normaliser<type>> \
                normaliser{norm_str ? get_normaliser<type>(norm_str) : \
                new gtd::gen_normaliser<type, \
                gtd::mass_normaliser<type>, gtd::time_normaliser<type>>}; \
                if (!normaliser) { \
                    std::cerr << "Error: invalid normaliser argument provided (empty).\n"; \
                    return 1; \
                } \
                gtd::preprocess(data_dir, val_dir, dppf, vppf, *normaliser, bval); \
                normaliser.reset(); // free the normaliser after preprocessing as it's no longer needed
                PP_BLOCK(long double, true)
            } catch (const gtd::invalid_3bod_ld_size &_e) {
                std::cerr << _e.what() << '\n';
                return 1;
            } catch (const gtd::invalid_3bod_format&) {
                try {
                    PP_BLOCK(double, true)
                } catch (const gtd::invalid_3bod_format&) {
                    try {
                        PP_BLOCK(float, true)
                    } catch (const gtd::invalid_3bod_format &_inv) {
                        std::cerr << _inv.what() << '\n';
                        return 1;
                    }
                }
            }
        } catch (const std::exception &_e) {
            std::cerr << _e.what() << '\n';
            return 1;
        }
        // delete [] normv_tpath;
        // delete [] normv_vpath;
        return 0;
    } // after this block it is certain the NN will be run
    const char *ftype = parser.get_arg("--fp_type"); // floating point type of weights and biases in NN
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
    long double learning_rate = parser.get_arg("--learning_rate", std::numeric_limits<long double>::quiet_NaN());
    if (learning_rate != learning_rate)
        learning_rate = DEF_LR;
    else {
        if (learning_rate <= 0) {
            std::cerr << "Error: learning rate must be positive.\n";
            return 1;
        }
    }
    uint64_t num_passes = parser.get_arg("--num_passes", DEF_NUM_PASSES);
    if (!num_passes) {
        std::cerr << "Error: number of passes over entire dataset cannot be zero.\n";
        return 1;
    }
    uint64_t pairs_per_file = parser.get_arg("--pairs_per_file", DEF_PPF);
    if (!pairs_per_file) {
        std::cerr << "Error: number of pairs per file must be positive.\n";
        return 1;
    }
    uint64_t batch_size = parser.get_arg("--batch_size", DEF_PPF);
    if (!batch_size) {
        std::cerr << "Error: batch size must be positive.\n";
        return 1;
    } // Doesn't include the last update if the below were not a whole number:
    // uint64_t tot_num_updates = (num_passes*num_files*pairs_per_file)/batch_size; // floor division is good here
    // if (batch_size > num_passes*pairs_per_file) {
    //     std::cerr << "Error: the batch size cannot be greater than the total number of training examples.\n";
    //     return 1;
    // }
    bool allocate_dataset = parser.get_arg("--allocate_dataset", true);
    const char *model_str = parser.get_arg("--from_model");
    const char *layers = parser.get_arg(LAYERS_RGX);
    const char *weights_path = parser.get_arg("-n");
    const char *edir = parser.get_arg("-o");
    lossrf = parser.get_arg("--loss_rec_freq", DEF_REC_FREQ);
    if (!parser.empty()) {
        std::cerr << "Error: unrecognised arguments:\n";
        for (const auto &[_, arg] : parser)
            fprintf(stderr, "\t%s\n", arg.c_str());
        return 1;
    }
    if (!lossrf) {
        std::cerr << "Error: loss rec. freq. cannot be zero.\n";
        return 1;
    }
    char from_model[PATH_MAX];
    if (model_str) {
        if (layers) {
            std::cerr << "Error: both the model from which to load and the layout of the NN cannot be specified.\n";
            return 1;
        }
        if (!realpath(model_str, from_model)) {
            std::cerr << "Error: could not obtain absolute resolved pathname for \"" << model_str << "\".\n";
            return 1;
        }
        model_str = from_model;
    } else {
        if (!layers)
            layers = "--layers=22-128:1,128-64:1,64-18:0";
    }
    char cwd[PATH_MAX]; // for log file, in case `data_dir` is "."
    if (*data_dir == '.' && !*(data_dir + 1)) {
        if (!getcwd(cwd, PATH_MAX)) {
            std::cerr << "Error: could not obtain current working directory.\nReason: " << strerror(errno) << '\n';
            return 1;
        }
        data_dir = cwd;
    }
    std::unique_ptr<char[]> exp_dir{};
#define CREATE_DIR \
    if (!edir) { \
        exp_dir.reset(gml::gen::now_str("NN_experiment_", "")); \
        edir = exp_dir.get(); \
    } \
    if (mkdir(edir, S_IRWXU | S_IRWXG | S_IRWXO) == -1) { \
        std::cerr << "Error: could not create experiment directory \"" << edir << "\".\n"; \
        return 1; \
    } \
    if (chdir(edir) == -1) { \
        std::cerr << "Error: could not change directory to experiment directory \"" << edir << "\".\n"; \
        return 1; \
    } \
    std::cout << "Created and changed directory to experiment directory \"" << edir << "\".\n";
    char *exp_info = new char[gtd::strlen_c(edir) + 5];
    gtd::strcpy_c(exp_info, edir);
    gtd::strcat_c(exp_info, ".txt");
    logger.open(exp_info, std::ios_base::out | std::ios_base::trunc);
    if (!logger.good()) {
        std::cerr << "Error: could not open log file \"" << exp_info << "\".\n";
        delete [] exp_info;
        return 1;
    }
    delete [] exp_info;
    logger << "Experiment name: \"" << edir << "\"\nTraining data directory: \"" << data_dir
           << "\"\nValidation data directory: \"" << (val_dir ? val_dir : "(None)") << "\nNormalisation type: "
           << (norm_str ? norm_str + 7 : "N/A") << "\nLearning rate = " << learning_rate
           << "\nRequested num. passes over training set = " << num_passes
           << "\nNumber of pairs to generate per file = " << pairs_per_file << "\nBatch size = " << batch_size
           << "\nAllocate dataset: " << std::boolalpha << allocate_dataset;
    if (model_str)
        logger << "\nModel pre-loaded from: \"" << model_str << "\"\n";
    else
        logger << "\nModel layers: " << (layers + 9) << '\n';
    /* if (weights_path)
        logger << "Weights output path: \"" << weights_path << "\"\n";
    else
        logger << "Weights output path: (to be computed at exit)\n"; */
    if (pp) {
        try {
            try {
                PP_BLOCK(long double, false)
#define NN_BLOCK(type) \
                if (nn_ld) \
                    build_and_run_nn<type, long double>(dppf, vppf, layers, model_str, num_passes, \
                                                        pairs_per_file, batch_size, learning_rate, weights_path, \
                                                        allocate_dataset); \
                else if (nn_d) \
                    build_and_run_nn<type, double>(dppf, vppf, layers, model_str, num_passes, \
                                                   pairs_per_file, batch_size, learning_rate, weights_path, \
                                                   allocate_dataset); \
                else \
                    build_and_run_nn<type, float>(dppf, vppf, layers, model_str, num_passes, \
                                                  pairs_per_file, batch_size, learning_rate, weights_path, \
                                                  allocate_dataset);
                CREATE_DIR
                NN_BLOCK(long double)
            } catch (const gtd::invalid_3bod_ld_size &_e) {
                std::cerr << _e.what() << '\n';
                return 1;
            } catch (const gtd::invalid_3bod_format&) {
                try {
                    PP_BLOCK(double, false)
                    CREATE_DIR
                    NN_BLOCK(double)
                } catch (const gtd::invalid_3bod_format&) {
                    try {
                        PP_BLOCK(float, false)
                        CREATE_DIR
                        NN_BLOCK(float)
                    } catch (const gtd::invalid_3bod_format &_inv) {
                        std::cerr << _inv.what() << '\n';
                        return 1;
                    }
                }
            }
        } catch (const std::exception &_e) {
            std::cerr << _e.what() << '\n';
            return 1;
        }
        return 0;
    }
    if (npp) {
        dppf.reset(gtd::find_files(data_dir, ".3bod", true));
        if (!dppf) {
            std::cerr << "Error: could not find any training .3bod files in directory \"" << data_dir << "\".\n";
            return 1;
        }
        if (val_dir) {
            vppf.reset(gtd::find_files(val_dir, ".3bod", true));
            if (!vppf) {
                std::cerr << "Error: could not find any validation .3bod files in directory \"" << val_dir << "\".\n";
                return 1;
            }
        }
    } else {
        dppf.reset(gtd::find_files(data_dir, ".3bodpp", true));
        if (!dppf) {
            std::cerr << "Error: could not find any training .3bodpp files in directory \"" << data_dir << "\".\n";
            return 1;
        }
        if (val_dir) {
            vppf.reset(gtd::find_files(val_dir, ".3bodpp", true));
            if (!vppf) {
                std::cerr << "Error: could not find any validation .3bodpp files in directory \"" << val_dir << "\".\n";
                return 1;
            }
        }
    }
    CREATE_DIR
    try {
        try {
            NN_BLOCK(long double)
        } catch (const gtd::invalid_3bod_format&) {
            try {
                NN_BLOCK(double)
            } catch (const gtd::invalid_3bod_format&) {
                try {
                    NN_BLOCK(float)
                } catch (const gtd::invalid_3bod_format &_inv) {
                    std::cerr << "Error: " << _inv.what() << '\n';
                    return 1;
                }
            }
        }
    } catch (const std::exception &_e) {
        std::cerr << "Error: " << _e.what() << '\n';
        return 1;
    }
    return 0;
}

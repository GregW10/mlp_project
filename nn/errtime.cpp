#include "../data_management/datsup.hpp"
#include "../glib/ml/gregffnn.hpp"
#include "../glib/misc/gregparse.hpp"
#include "nnsup.hpp"
#include <map>

const char *tdir{};
const char *from_model{};
const char *outf{};
uint64_t ppf;
uint64_t npasses;
std::streamsize prec;

std::map<long double, std::pair<long double, uint64_t>> final_tot_losses;

std::mutex mutex;

uint64_t unique_pairs(uint64_t _num) {
    return (_num*(_num - 1))/2;
}

template <typename D, typename W> requires (std::is_floating_point_v<D> && std::is_floating_point_v<W>)
void run_nn(gml::ffnn<W> *net, const gtd::f3bodr<D> *readers, uint64_t _size) {
    if (!net || !readers || !_size)
        return;
    // gml::ffnn<W> net{from_model}; // this is what might throw
    std::mt19937_64 rng{std::random_device{}()};
    std::uniform_int_distribution<uint64_t> dist;
    uint64_t pcounter = 0;
    uint64_t ncounter = 0;
    const gtd::f3bodr<D> *_rdr = readers;
    std::vector<uint64_t> max_pairs;
    max_pairs.reserve(_size);
    while (pcounter++ < _size)
        max_pairs.push_back(unique_pairs(_rdr++->header().N));
    const typename gtd::f3bodr<D>::hdr_t *header;
    std::map<long double, std::pair<long double, uint64_t>> losses;
    // uint64_t *max_ppf;
    uint64_t idx_temp;
    uint64_t idx_lo;
    uint64_t idx_hi;
    gml::vector<D> _input{(uint64_t) 22};
    gml::vector<D> _y{(uint64_t) 18};
    const gml::vector<W> *_f{};
    uint64_t _nm1;
    std::pair<long double, uint64_t> *ptr{};
    // std::cout << "Got before main loop" << std::endl;
    while (ncounter++ < npasses) {
        _rdr = readers;
        for (const uint64_t &_mppf : max_pairs) {
            header = &_rdr->header();
            if (header->N <= 1) { // case for empty or single-entry .3bodpp file
                ++_rdr;
                continue;
            }
            // std::cout << "ppf: " << ppf << ", N: " << header->N << ", _mppf: " << _mppf << std::endl;
            _input[0] = header->masses[0];
            _input[1] = header->masses[1];
            _input[2] = header->masses[2];
            // std::cout << "Got before ppf." << std::endl;
            if (_mppf > ppf) {
                dist.param(std::uniform_int_distribution<uint64_t>::param_type{0, header->N - 1});
                pcounter = 0;
                while (pcounter++ < ppf) {
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
                    gml::gen::memcopy(_input.begin() + 3, _rdr->entry_at(idx_lo).positions, sizeof(D), 18);
                    _input[21] = header->dt*(idx_hi - idx_lo);
                    gml::gen::memcopy(_y.begin(), _rdr->entry_at(idx_hi).positions, sizeof(D), 18);
                    ptr = &losses[(long double) _input[21]];
                    ptr->first += net->fpass_loss(_input, _y);
                    ++(ptr->second);
                }
            } else {
                _nm1 = header->N - 1;
                for (idx_lo = 0; idx_lo < _nm1; ++idx_lo) {
                    gml::gen::memcopy(_input.begin() + 3, _rdr->entry_at(idx_lo).positions, sizeof(D), 18);
                    for (idx_hi = idx_lo + 1; idx_hi < header->N; ++idx_hi) {
                        _input[21] = header->dt*(idx_hi - idx_lo);
                        gml::gen::memcopy(_y.begin(), _rdr->entry_at(idx_hi).positions, sizeof(D), 18);
                        ptr = &losses[(long double) _input[21]];
                        ptr->first += net->fpass_loss(_input, _y);
                        ++(ptr->second);
                    }
                }
            }
        }
        ++_rdr;
    }
    std::lock_guard<std::mutex> guard{mutex};
    for (const auto &[_dt, _pair] : losses) {
        ptr = &final_tot_losses[_dt];
        ptr->first += _pair.first;
        ptr->second += _pair.second;
        // ++(ptr->second);
    }
}

template <typename D, typename W> requires (std::is_floating_point_v<D> && std::is_floating_point_v<W>)
void run_thread(gml::ffnn<W> *net, const gtd::f3bodr<D> *readers, uint64_t _size) {
    try {
        run_nn<D, W>(net, readers, _size);
    } catch (const std::exception &_e) {
        std::cerr << "An exception occurred.\nwhat(): " << _e.what() << '\n';
    }
}

template <typename D, typename W> requires (std::is_floating_point_v<D> && std::is_floating_point_v<W>)
void dispatch_threads(const std::vector<gtd::f3bodr<D>> &readers, unsigned int maxt) {
    std::vector<gml::ffnn<W>> nns;
    unsigned int counter = maxt;
    while (counter --> 0)
        nns.emplace_back(from_model);
    std::vector<std::thread> threads;
    typename std::vector<gtd::f3bodr<D>>::size_type _size = readers.size();
    const gtd::f3bodr<D> *ptr = readers.data();
    gml::ffnn<W> *nptr = nns.data();
    if (maxt >= _size) {
        while (maxt --> 0)
            threads.emplace_back(&run_thread<D, W>, nptr++, ptr++, 1);
    } else {
        counter = 1;
        long double num_per_thread = ((long double) _size)/maxt;
        uint64_t next_offset;
        uint64_t curr_offset = 0;
        while (counter <= maxt) {
            next_offset = (uint64_t) roundl(counter++*num_per_thread);
            // std::cout << "Counter: " << counter << ", next_offset: " << next_offset << ", curr_offset: " << curr_offset << std::endl;
            threads.emplace_back(&run_thread<D, W>, nptr++, ptr + curr_offset, next_offset - curr_offset);
            curr_offset = next_offset;
        }
    }
    for (std::thread &_t : threads)
        _t.join();
}

template <typename D> requires (std::is_floating_point_v<D>)
void write_results(std::ostream &os) {
    std::unique_ptr<std::vector<std::string>> ptr{gtd::find_files(tdir, ".normv", true)};
    if (ptr && ptr->size() == 1) {
        gtd::time_normaliser<D> normaliser;
        try {
            normaliser.load_normv(ptr->operator[](0).c_str());
        } catch (...) {
            ptr.reset();
            goto _rest;
        }
        ptr.reset();
        os << "dt,error\r\n";
        for (const auto &[_dt, _pair] : final_tot_losses)
            os << _dt*normaliser.max_time() << ',' << _pair.first/_pair.second << "\r\n";
        return;
    }
    ptr.reset();
    _rest:
    os << "normalised_dt,error\r\n";
    for (const auto &[_dt, _pair] : final_tot_losses)
        os << _dt << ',' << _pair.first/_pair.second << "\r\n";
}

template <typename D> requires (std::is_floating_point_v<D>)
void run(const std::vector<std::string> &files) {
    std::vector<gtd::f3bodr<D>> readers;
    typename std::vector<gtd::f3bodr<D>>::size_type _size = files.size();
    readers.reserve(_size);
    for (const std::string &_f : files)
        readers.emplace_back(_f.c_str());
    std::shuffle(readers.begin(), readers.end(), std::mt19937_64{std::random_device{}()});
    unsigned int maxt = std::thread::hardware_concurrency();
    try {
        dispatch_threads<D, long double>(readers, maxt);
    } catch (const gml::exceptions::invalid_nnw_format&) {
        try {
            dispatch_threads<D, double>(readers, maxt);
        } catch (const gml::exceptions::invalid_nnw_format&) {
            dispatch_threads<D, float>(readers, maxt); // let it get caught in `main` if it fails
        }
    }
    std::ofstream out{outf, std::ios_base::out | std::ios_base::trunc};
    if (!out.good()) {
        std::cerr << "Error: could not open \"" << outf << "\" - writing to standard output instead.\n";
        std::cout.precision(prec);
        write_results<D>(std::cout);
        return;
    }
    out.precision(prec);
    write_results<D>(out);
    out.close();
}

int main(int argc, char **argv) {
    gtd::parser parser{argc, argv};
    from_model = parser.get_arg("--from_model");
    if (!from_model) {
        std::cerr << "Error: path to .nnw model weights must be passed.\n";
        return 1;
    }
    tdir = parser.get_arg("--test_dir");
    if (!tdir)
        tdir = ".";
    outf = parser.get_arg("-o");
    if (!outf)
        outf = "test_err_vs_time.csv";
    npasses = parser.get_arg("--num_passes", (uint64_t) 1);
    if (!npasses) {
        std::cerr << "Error: number of passes cannot be zero.\n";
        return 1;
    }
    ppf = parser.get_arg("--ppf", (uint64_t) 1'000);
    if (!ppf) {
        std::cerr << "Error: number of pairs per file must be positive.\n";
        return 1;
    }
    std::unique_ptr<std::vector<std::string>> files{gtd::find_files(tdir, ".3bodpp", true)};
    if (!files) {
        std::cerr << "Error: no .3bodpp files found within \"" << tdir << "\".\n";
        return 1;
    }
    prec = parser.get_arg("--precision", (std::streamsize) 20);
    try {
        try {
            run<long double>(*files);
        } catch (const gtd::invalid_3bod_ld_size &_e) {
            std::cerr << "what(): " << _e.what() << '\n';
            return 1;
        } catch (const gtd::invalid_3bod_format&) {
            try {
                run<double>(*files);
            } catch (const gtd::invalid_3bod_format&) {
                try {
                    run<float>(*files);
                } catch (const gtd::invalid_3bod_format &_e) {
                    std::cerr << "what(): " << _e.what() << '\n';
                    return 1;
                }
            }
        }
    } catch (const std::exception &_e) {
        std::cerr << "An exception occurred.\nwhat(): " << _e.what() << '\n';
        return 1;
    }
    return 0;
}

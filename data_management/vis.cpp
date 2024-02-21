#include "datsup.hpp"
#include <dirent.h>

#ifdef __linux__ // necessary since Linux and MacOS have diff. names for the `struct timespec` fields in `struct stat`
#define BUFF_SEC buff.st_mtim.tv_sec
#define BUFF_NSEC buff.st_mtim.tv_nsec
#elif defined(__APPLE__)
#define BUFF_SEC buff.st_mtimespec.tv_sec
#define BUFF_NSEC buff.st_mtimespec.tv_nsec
#endif

#define VIEW_ALL R"(^view$)"
#define VIEW_HEADER R"(^header$)"
#define VIEW_RANGE R"(^r:\d{0,18}-\d{0,18}$)"
#define VIEW_ENTRY R"(^e:\d{1,18}$)"
#define DISTANCE R"(^dist$)"
#define DISTANCE_RANGE R"(^dist:\d{1,18}-\d{1,18}$)"
#define MAX_SEP
#define MAX_VEL
#define MAX_FORCE
#define MAX_ACC
#define MIN_SEP
#define MIN_VEL
#define MEAN_SEP // with SD?
// include COM calculations
// SHOULD COMBINE LAST 7 ABOVE OPTIONS INTO ONE AS THEY INVOLVE SIMILAR CALCULATIONS (would be a waste to sep. them)
// perhaps actually last 9, and can be a range as well (maybe call it "advanced analysis" or something)

std::regex va_rgx{VIEW_ALL};
std::regex vh_rgx{VIEW_HEADER};
std::regex vr_rgx{VIEW_RANGE};
std::regex ve_rgx{VIEW_ENTRY};
std::regex d_rgx{DISTANCE};
std::regex dr_rgx{DISTANCE_RANGE};

template <typename T>
using etype = typename gtd::f3bodr<T>::entry_type;

template <typename T>
using htype = typename gtd::f3bodr<T>::hdr_t;

[[noreturn]] void log_no_readings() {
    std::cerr << BOLD_TXT(MAGENTA_TXT("No evolutions recorded in .3bod file.\n"));
    exit(1);
}

#define BOLD "\033[1m" // bold text
#define UDL "\033[4m" // underlined text
#define UDLRST "\033[24m" // reset underlined text
#define RST "\033[0m" // reset text flags
#define ERRCOL "\033[38;5;214m" // "Error: " colour
#define ERRTCOL "\033[38;5;226m" // colour of text coming after "Error: "
#define ETCOL "\033[38;5;174m" // entry text colour
#define ENCOL "\033[38;5;199m" // entry number colour
#define BTCOL "\033[96m" // body text colour
#define BNCOL "\033[95m" // body number colour
#define PVCOL "\033[92m" // pos/vel text colour
#define VLCOL "\033[39m" // vector letter colour
#define PMCOL "\033[38;5;93m" // plus/minus colour
#define UCOL "\033[38;5;15m" // units colour
#define NNCOL "\033[38;5;118m" // numeric number colour
#define LNCOL "\033[38;5;2m" // logging number colour
#define HCOL "\033[38;5;123m" // hyphen colour
#define SHCOL "\033[38;5;145m" // sub-hyphen colour
#define HDFCOL "\033[38;5;21m" // header field colour
#define ICHCOL "\033[38;5;123m" // info char colour

template <typename T> requires (std::is_floating_point_v<T>)
void print_entry(const etype<T> &e) {
    T *pptr = (T*) &e;
    T *vptr = ((T*) &e) + 9; // offset of 9 to reach velocities
    unsigned char counter = 0;
    while (counter++ < 3) {
        std::cout << BOLD BTCOL "Body " BNCOL << +counter << SHCOL "\n------\n" RST PVCOL UDL "Position:" UDLRST " "
        BOLD NNCOL << *pptr++;
        if (*pptr < 0)
            std::cout << VLCOL "i " PMCOL "- " BOLD NNCOL << -*pptr++;
        else
            std::cout << VLCOL "i " PMCOL "+ " BOLD NNCOL << *pptr++;
        if (*pptr < 0)
            std::cout << VLCOL "j " PMCOL "- " BOLD NNCOL << -*pptr++;
        else
            std::cout << VLCOL "j " PMCOL "+ " BOLD NNCOL << *pptr++;
        std::cout << VLCOL "k " UCOL "m" RST "\n" PVCOL UDL "Velocity:" UDLRST " " BOLD NNCOL << *vptr++;
        if (*vptr < 0)
            std::cout << VLCOL "i " PMCOL "- " BOLD NNCOL << -*vptr++;
        else
            std::cout << VLCOL "i " PMCOL "+ " BOLD NNCOL << *vptr++;
        if (*vptr < 0)
            std::cout << VLCOL "j " PMCOL "- " BOLD NNCOL << -*vptr++;
        else
            std::cout << VLCOL "j " PMCOL "+ " BOLD NNCOL << *vptr++;
        std::cout << VLCOL "k " UCOL "m/s\n" /* SHCOL "------\n" */ RST;
    }
}

char *get_hyphens(uint64_t max_idx) {
    uint64_t num_digits = (uint64_t) (log10l((long double) max_idx) + 1);//will automatically be rounded down if decimal
    char *hyphens = new char[num_digits + 9]{};
    char *cptr = hyphens + num_digits + 7;
    *cptr-- = '\n';
    while (cptr > hyphens)
        *cptr-- = '-';
    *cptr = '\n';
    return hyphens;
}

template <typename T> requires (std::is_floating_point_v<T>)
void view_all(const gtd::f3bodr<T> &reader) {
    const htype<T> &hdr = reader.header();
    if (!hdr.N)
        log_no_readings();
    char *hyphens = get_hyphens(hdr.N);
    uint64_t index = 0;
    unsigned char counter;
    for (const auto &e : reader) {
        std::cout << ETCOL "Entry " BOLD LNCOL << index++ << HCOL << hyphens;
        print_entry<T>(e);
        std::cout.put('\n');
    }
    delete [] hyphens;
}

template <typename T> requires (std::is_floating_point_v<T>)
void view_header(const gtd::f3bodr<T> &reader) {
    const htype<T> &hdr = reader.header();
    std::cout << BOLD HDFCOL UDL "Head:" UDLRST << ' ' << ICHCOL << hdr.hd[0] << HDFCOL ", " ICHCOL  << hdr.hd[1] <<
    HDFCOL UDL "\nsizeof(T):" UDLRST NNCOL " " << +hdr.flt_size
    << UCOL " bytes\n" HDFCOL UDL "Collisions:" UDLRST " " ICHCOL
    << std::boolalpha << ((bool) (hdr.lds_coll & 0b00000001))
    << HDFCOL UDL "\nsizeof(long double):" UDLRST " " NNCOL
    << (hdr.lds_coll >> 1) << UCOL " bytes\n" HDFCOL UDL "Time-step:" UDLRST " " NNCOL
    << hdr.dt << UCOL " seconds\n" HDFCOL UDL "Num. Epochs:" UDLRST " " NNCOL << hdr.N
    << HDFCOL UDL "\nMasses:\n" UDLRST
    << BTCOL "\tBody " BNCOL "1" BTCOL ": " NNCOL << hdr.masses[0]
    << UCOL " kg\n\t" BTCOL "Body " BNCOL "2" BTCOL ": " NNCOL << hdr.masses[1]
    << UCOL " kg\n\t" BTCOL "Body " BNCOL "3" BTCOL ": " NNCOL << hdr.masses[2]
    << UCOL " kg\n" HDFCOL UDL "Radii:\n" UDLRST
    << BTCOL "\tBody " BNCOL "1" BTCOL ": " NNCOL
    << hdr.radii[0] << UCOL " m\n\t" BTCOL "Body " BNCOL "2" BTCOL ": " NNCOL
    << hdr.radii[1] << UCOL " m\n\t" BTCOL "Body " BNCOL "3" BTCOL ": " NNCOL << hdr.radii[2]
    << UCOL " m\n" HDFCOL UDL "nCORs:\n" UDLRST << BTCOL "\tBody " BNCOL "1" BTCOL ": " NNCOL << hdr.nCORs[0]
    << BTCOL "\n\tBody " BNCOL "2" BTCOL ": " NNCOL << hdr.nCORs[1] << BTCOL "\n\tBody " BNCOL "3" BTCOL ": " NNCOL
    << hdr.nCORs[2] << std::endl;
}

uint64_t to_uint64(const char *str, const char **endptr = nullptr, char terminator = 0) {
    if (!str || !*str)
        return -1;
    uint64_t val = 0;
    goto start;
    while (*str != terminator) {
        val *= 10;
        start:
        val += *str++ - 48;
    }
    if (endptr)
        *endptr = str;
    return val;
}

template <typename T> requires (std::is_floating_point_v<T>)
void view_range(const gtd::f3bodr<T> &reader, const char *action) { // over the CLOSED interval [lower,higher]
    const char *ptr = action + 2; // move past the "r:" to point to the first digit of the starting index
    uint64_t idx1 = to_uint64(ptr, &ptr, '-');
    uint64_t idx2 = to_uint64(ptr + 1); // point to the first digit of the end index
    if (idx1 > idx2) {
        std::cerr << BOLD ERRCOL UDL "Error:" UDLRST ERRTCOL
        " starting index in range specified is higher than end index.\n";
        exit(1);
    }
    const htype<T> &hdr = reader.header();
    if (idx1 >= hdr.N) {
        std::cerr << BOLD UDL ERRCOL "Error:" UDLRST ERRTCOL
        " starting index in range is out-of-bounds for a .3bod file with " NNCOL << hdr.N << ERRTCOL " entries.\n";
        exit(1);
    }
    if (idx2 >= hdr.N) {
        std::cerr << BOLD UDL ERRCOL "Error:" UDLRST ERRTCOL
        " end index in range is out-of-bounds for a .3bod file with " NNCOL << hdr.N << ERRTCOL " entries.\n";
        exit(1);
    }
    char *hyphens = get_hyphens(idx2);
    std::cout << "\033[38;5;6mRange: " BOLD "\033[38;5;7m[" LNCOL << idx1 << "\033[38;5;11m, " LNCOL
              << idx2 << "\033[38;5;7m]\n\n";
    auto it = reader.begin() + idx1;
    while (idx1 <= idx2) {
        std::cout << ETCOL "Entry " BOLD LNCOL << idx1++ << HCOL << hyphens;
        print_entry<T>(*it++);
        std::cout.put('\n');
    }
    delete [] hyphens;
}

template <typename T> requires (std::is_floating_point_v<T>)
void view_entry(const gtd::f3bodr<T> &reader, const char *action) {
    uint64_t idx = to_uint64(action + 2);
    const htype<T> &hdr = reader.header();
    if (idx >= hdr.N) {
        std::cerr << BOLD UDL ERRCOL "Error:" UDLRST ERRTCOL " specified index is out-of-range for a .3bod file with "
        NNCOL << hdr.N << ERRTCOL " entries.\n";
        exit(1);
    }
    char *hyphens = get_hyphens(idx);
    std::cout << ETCOL "Entry " BOLD LNCOL << idx << HCOL << hyphens;
    print_entry<T>(reader.entry_at(idx));
    std::cout.put('\n');
    delete [] hyphens;
}

template <typename T> requires (std::is_floating_point_v<T>)
void distance(const gtd::f3bodr<T> &reader) {

}

template <typename T> requires (std::is_floating_point_v<T>)
void distance_range(const gtd::f3bodr<T> &reader, const char *action) {

}

template <typename T> requires (std::is_floating_point_v<T>)
[[noreturn]] void perform_action(const gtd::f3bodr<T> &reader, gtd::parser &parser) {
    std::streamsize prec = parser.get_arg("--precision", (std::streamsize) std::cout.precision());
    if (prec < 0) {
        std::cerr << BOLD ERRCOL UDL "Error:" UDLRST ERRTCOL " precision cannot be negative.\n";
        exit(1);
    }
    std::cout.precision(prec);
    const char *action = parser.get_arg("--action");
    if (!action || std::regex_match(action, va_rgx))
        view_all(reader);
    else if (std::regex_match(action, vh_rgx))
        view_header(reader);
    else if (std::regex_match(action, vr_rgx))
        view_range(reader, action);
    else if (std::regex_match(action, ve_rgx))
        view_entry(reader, action);
    else if (std::regex_match(action, d_rgx))
        distance(reader);
    else if (std::regex_match(action, dr_rgx))
        distance_range(reader, action);
    else {
        std::cerr << BOLD UDL ERRCOL "Error:" UDLRST ERRTCOL " invalid action specified.\n";
        exit(1);
    }
    exit(0);
}

int main(int argc, char **argv) {
    gtd::parser parser{argc, argv};
    char *fpath = const_cast<char*>(parser.get_arg("--file"));
    bool free_path = false;
    if (!fpath)
        fpath = const_cast<char*>(parser.get_arg("-f"));
    if (!fpath) { // if no path to a .3bod file is given, I take the latest modified one found using directory entries
        DIR *dir = opendir(".");
        if (!dir) {
            std::cerr << BOLD UDL ERRCOL "Error:" UDLRST ERRTCOL " no path to .3bod file provided and could not open "
                         "current directory to search for .3bod files.\n";
            return 1;
        }
        struct dirent *entry;
        struct stat buff{};
        time_t latest_s = 0;
        decltype(BUFF_NSEC) latest_ns = 0;
        long max_len = pathconf(".", _PC_PATH_MAX);
        if (max_len == -1)
            max_len = sizeof(std::declval<struct dirent>().d_name);
        fpath = new char[max_len + 1]{}; // need to allocate memory as `entry` points to a changing static `dirent`
        while ((entry = readdir(dir))) {
            if (gtd::endswith(entry->d_name, ".3bod")) {
                if (stat(entry->d_name, &buff) != -1) { // ignore error and continue onto next file
                    if (BUFF_SEC > latest_s || (BUFF_SEC == latest_s && BUFF_NSEC > latest_ns)) {
                        gtd::strcpy_c(fpath, entry->d_name);
                        latest_s = BUFF_SEC;
                        latest_ns = BUFF_NSEC;
                    }
                }
            }
        }
        closedir(dir);
        if (!*fpath) {
            std::cerr << BOLD UDL ERRCOL "Error:" UDLRST ERRTCOL " no .3bod files found in current directory.\n";
            delete [] fpath;
            return 1;
        }
        free_path = true;
    }
    std::cout << "\033[1m\033[37mFile: \033[0m\033[32m\"\033[33m" << fpath << "\033[32m\"\033[0m\n\n";
    // Maybe manually open file here first to check its integrity
    try {
        try {
            gtd::f3bodr<long double> reader{fpath};
            perform_action(reader, parser);
        }
        catch (const gtd::invalid_3bod_ld_size&) {
            std::cerr << BOLD ERRCOL UDL "Error:" UDLRST ERRTCOL " \"" HDFCOL "sizeof(long double)" ERRTCOL
            "\" does not match that reported in the .3bod file.\n";
            return 1;
        }
        catch (const gtd::invalid_3bod_T_size&) {
            try {
                gtd::f3bodr<double> reader{fpath};
                perform_action(reader, parser);
            } catch (const gtd::invalid_3bod_T_size&) {
                try {
                    gtd::f3bodr<float> reader{fpath};
                    perform_action(reader, parser);
                } catch (const gtd::invalid_3bod_T_size&) {
                    std::cerr << BOLD_TXT(MAGENTA_TXT("Error: "))
                    YELLOW_TXT("reported floating point data type size does not match the size of a ")
                    GREEN_TXT("\"") BLUE_TXT("long double") GREEN_TXT("\"") YELLOW_TXT(", ") GREEN_TXT("\"")
                    BLUE_TXT("double") GREEN_TXT("\"") YELLOW_TXT(" or ") GREEN_TXT("\"") BLUE_TXT("float")
                    GREEN_TXT("\"") YELLOW_TXT(".\n");
                    return 1;
                }
            }
        }
    }
    catch (const gtd::invalid_3bod_format &e) {
        std::cerr << BOLD ERRCOL UDL "Error:" UDLRST ERRTCOL " invalid .3bod format.\n";
        return 1;
    }
    if (free_path)
        delete [] fpath;
    return 0;
}

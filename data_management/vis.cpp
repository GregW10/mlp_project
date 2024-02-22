#include "datsup.hpp"
#include <dirent.h>

#ifdef __linux__ // necessary since Linux and MacOS have diff. names for the `struct timespec` fields in `struct stat`
#define BUFF_SEC buff.st_mtim.tv_sec
#define BUFF_NSEC buff.st_mtim.tv_nsec
#elif defined(__APPLE__)
#define BUFF_SEC buff.st_mtimespec.tv_sec
#define BUFF_NSEC buff.st_mtimespec.tv_nsec
#endif

#define VIEW_HEADER R"(^h$)"
#define VIEW_ALL R"(^v-p?v?a?c?e?$)"
#define VIEW_RANGE R"(^v-p?v?a?c?e?:\d{1,18}-\d{1,18}$)"
#define VIEW_ENTRY R"(^v-p?v?a?c?e?:\d{1,18}$)"
#define ANALYSIS R"(^a$)"
#define ANALYSIS_RANGE R"(^a:\d{1,18}-\d{1,18}$)"
// #define DISTANCE R"(^dist$)"
// #define DISTANCE_RANGE R"(^dist:\d{1,18}-\d{1,18}$)"
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
std::regex d_rgx{ANALYSIS};
std::regex dr_rgx{ANALYSIS_RANGE};

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
#define COMPCOL "\033[38;5;134m" // COM pos. desc. colour
#define COMVCOL "\033[38;5;128m" // COM vel. desc. colour
#define ETCOL "\033[38;5;213m" // energy text color
#define ENCOL "\033[38;5;234m" // energy number color
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

void *masses{}; // pointer to masses (only used for analysis calculations)
char *path{};

std::vector<void*> to_free = {masses, path};

void free_all() {
    for (void* &vptr : to_free)
        delete [] vptr;
}

#pragma pack(push, 1) // not really necessary, since it's guaranteed that `sizeof(T)` >= 4
template <typename T> requires (std::is_floating_point_v<T>)
struct vector { // I am not using my `gtd::vector3D<T>` class in this program
    T x{}, y{}, z{};
};
#pragma pack(pop)

template <typename T> requires (std::is_floating_point_v<T>)
inline vector<T> operator+(const vector<T> &v1, const vector<T> &v2) {
    return {v1.x + v2.x, v1.y + v2.y, v1.z + v2.z};
}

template <typename T> requires (std::is_floating_point_v<T>)
vector<T> operator-(const vector<T> &v1, const vector<T> &v2) {
    return {v1.x - v2.x, v1.y - v2.y, v1.z - v2.z};
}

template <typename T> requires (std::is_floating_point_v<T>)
T operator*(const vector<T> &v1, const vector<T> &v2) { // dot product
    return v1.x*v2.x + v1.y*v2.y + v1.z*v2.z;
}

template <typename T> requires (std::is_floating_point_v<T>)
vector<T> operator*(T scalar, const vector<T> &vec) {
    return {scalar*vec.x, scalar*vec.y, scalar*vec.z};
}

template <typename T> requires (std::is_floating_point_v<T>)
vector<T> operator*(const vector<T> &vec, T scalar) {
    return {scalar*vec.x, scalar*vec.y, scalar*vec.z};
}

template <typename T> requires (std::is_floating_point_v<T>)
vector<T> operator/(const vector<T> &vec, T scalar) {
    return {vec.x/scalar, vec.y/scalar, vec.z/scalar};
}

template <typename T> requires (std::is_floating_point_v<T>)
vector<T> &operator+=(vector<T> &vec, const vector<T> &other) {
    vec.x += other.x;
    vec.y += other.y;
    vec.z += other.z;
    return vec;
}

template <typename T> requires (std::is_floating_point_v<T>)
vector<T> &operator-=(vector<T> &vec, const vector<T> &other) {
    vec.x -= other.x;
    vec.y -= other.y;
    vec.z -= other.z;
    return vec;
}

template <typename T> requires (std::is_floating_point_v<T>)
T mag_sq(const vector<T> &vec) {
    return vec.x**vec.x + vec.y**vec.y + vec.z**vec.z;
}

template <typename T> requires (std::is_floating_point_v<T>)
long double mag(const vector<T> &vec) {
    return sqrtl(vec.x**vec.x + vec.y**vec.y + vec.z**vec.z);
}

template <typename T> requires (std::is_floating_point_v<T>)
long double distance(const vector<T> &v1, const vector<T> &v2) {
    return mag(v1 - v2);
}

template <typename T> requires (std::is_floating_point_v<T>)
vector<T> com(const vector<T> &p1, const vector<T> &p2, const vector<T> &p3) { //hard-coded for 3 as only need that here
    return (*((T*) masses)*p1 + *(((T*) masses) + 1)*p2 + *(((T*) masses) + 2)*p3)/
           (*((T*) masses)    + *(((T*) masses) + 1)    + *(((T*) masses) + 2));
}

template <typename T> requires (std::is_floating_point_v<T>)
std::ostream &operator<<(std::ostream &os, const vector<T> &vec) {
    os << BOLD NNCOL << vec.x << VLCOL;
    if (vec.y < 0)
        std::cout << VLCOL "i " PMCOL "- " BOLD NNCOL << -vec.y;
    else
        std::cout << VLCOL "i " PMCOL "+ " BOLD NNCOL << vec.y;
    if (vec.z < 0)
        std::cout << VLCOL "j " PMCOL "- " BOLD NNCOL << -vec.z;
    else
        std::cout << VLCOL "j " PMCOL "+ " BOLD NNCOL << vec.z;
    return os << VLCOL "k" RST;
}

template <typename T> requires (std::is_floating_point_v<T>)
const vector<T> *calc_acc(const etype<T> &entry) {
    static vector<T> acc_vals[3];
    vector<T> r12 = *(((vector<T>*) &entry) + 1) - *((vector<T>*) &entry);
    vector<T> r13 = *(((vector<T>*) &entry) + 2) - *((vector<T>*) &entry);
    vector<T> r23 = *(((vector<T>*) &entry) + 2) - *(((vector<T>*) &entry) + 1);
    T r12_magsq = r12*r12;
    T r13_magsq = r13*r13;
    T r23_magsq = r23*r23;
    acc_vals[0] =  gtd::sys::G_SI*(r12*(*(((T*) masses) + 1)/(r12_magsq)) + r13*(*(((T*) masses) + 2)/(r13_magsq)));
    acc_vals[1] = -gtd::sys::G_SI*(r12*(*(((T*) masses))/(r12_magsq)) - r23*(*(((T*) masses) + 2)/(r23_magsq)));
    acc_vals[2] = -gtd::sys::G_SI*(r13*(*(((T*) masses))/(r13_magsq)) + r23*(*(((T*) masses) + 1)/(r23_magsq)));
    return acc_vals;
}

template <typename T> requires (std::is_floating_point_v<T>)
std::pair<T*, vector<T>> calc_e(const etype<T> &entry) { // returns energies, COM vel.
    static T e_vals[3]; // energy values, in order: kinetic energy, potential energy, total energy
    e_vals[1] = -((gtd::sys::G_SI*(*(((T*) masses) + 0))*(*(((T*) masses) + 1)))/
                        mag(*(((vector<T>*) &entry) + 1) - *(((vector<T>*) &entry) + 0)) +
                  (gtd::sys::G_SI*(*(((T*) masses) + 0))*(*(((T*) masses) + 2)))/
                        mag(*(((vector<T>*) &entry) + 2) - *(((vector<T>*) &entry) + 0)) +
                  (gtd::sys::G_SI*(*(((T*) masses) + 1))*(*(((T*) masses) + 2)))/
                        mag(*(((vector<T>*) &entry) + 2) - *(((vector<T>*) &entry) + 1)));
    vector<T> com_vel = com(*(((vector<T>*) &entry) + 3), *(((vector<T>*) &entry) + 4), *(((vector<T>*) &entry) + 5));
    e_vals[0] = 0.5l(*(((T*) masses) + 0)*(*(((vector<T>*) &entry) + 3) - com_vel) +
                     *(((T*) masses) + 1)*(*(((vector<T>*) &entry) + 4) - com_vel) +
                     *(((T*) masses) + 2)*(*(((vector<T>*) &entry) + 5) - com_vel));
    e_vals[2] = e_vals[0] + e_vals[1];
    return {e_vals, com_vel}; // have to calc. COM velocity, so might as well return it
}

template <typename T> requires (std::is_floating_point_v<T>)
std::tuple<vector<T>*, T*, vector<T>> calc_acc_e(const etype<T> &entry) { // returns accelerations, energies, COM vel.
    static vector<T> acc_vals[3];
    static T e_vals[3]; // energy values, in order: kinetic energy, potential energy, total energy
    vector<T> r12 = *(((vector<T>*) &entry) + 1) - *(((vector<T>*) &entry) + 0);
    vector<T> r13 = *(((vector<T>*) &entry) + 2) - *(((vector<T>*) &entry) + 0);
    vector<T> r23 = *(((vector<T>*) &entry) + 2) - *(((vector<T>*) &entry) + 1);
    T r12_magsq = r12*r12;
    T r13_magsq = r13*r13;
    T r23_magsq = r23*r23;
    acc_vals[0] =  gtd::sys::G_SI*(r12*(*(((T*) masses) + 1)/(r12_magsq)) + r13*(*(((T*) masses) + 2)/(r13_magsq)));
    acc_vals[1] = -gtd::sys::G_SI*(r12*(*(((T*) masses) + 0)/(r12_magsq)) - r23*(*(((T*) masses) + 2)/(r23_magsq)));
    acc_vals[2] = -gtd::sys::G_SI*(r13*(*(((T*) masses) + 0)/(r13_magsq)) + r23*(*(((T*) masses) + 1)/(r23_magsq)));
    e_vals[1] = -((gtd::sys::G_SI*(*(((T*) masses) + 0))*(*(((T*) masses) + 1)))/sqrtl(r12_magsq) +
                  (gtd::sys::G_SI*(*(((T*) masses) + 0))*(*(((T*) masses) + 2)))/sqrtl(r13_magsq) +
                  (gtd::sys::G_SI*(*(((T*) masses) + 1))*(*(((T*) masses) + 2)))/sqrtl(r23_magsq));
    vector<T> com_vel = com(*(((vector<T>*) &entry) + 3), *(((vector<T>*) &entry) + 4), *(((vector<T>*) &entry) + 5));
    e_vals[0] = 0.5l(*(((T*) masses) + 0)*(*(((vector<T>*) &entry) + 3) - com_vel) +
                     *(((T*) masses) + 1)*(*(((vector<T>*) &entry) + 4) - com_vel) +
                     *(((T*) masses) + 2)*(*(((vector<T>*) &entry) + 5) - com_vel));
    e_vals[2] = e_vals[0] + e_vals[1];
    return {acc_vals, e_vals, com_vel};
}

template <typename T> requires (std::is_floating_point_v<T>)
void print_pos(const vector<T> *pos) {
    std::cout << PVCOL UDL "Position:" UDLRST " " BOLD NNCOL << *pos << BOLD UCOL " m\n" RST;
}

template <typename T> requires (std::is_floating_point_v<T>)
void print_vel(const vector<T> *vel) {
    std::cout << PVCOL UDL "Velocity:" UDLRST " " BOLD NNCOL << *vel << BOLD UCOL " m/s\n" RST;
}

template <typename T> requires (std::is_floating_point_v<T>)
void print_acc(const vector<T> *acc) {
    std::cout << PVCOL UDL "Acceleration:" UDLRST " " BOLD NNCOL << *acc << BOLD UCOL " m/s^2\n" RST;
}

template <typename T, bool POS, bool VEL, bool ACC, bool COM, bool ENERGY> requires (std::is_floating_point_v<T>)
void print_entry(const etype<T> &entry) {
    static vector<T> *aptr;
    static T *eptr;
    unsigned char counter = 0;
    unsigned char idx = 0;
    if constexpr (ACC && ENERGY) {
        auto [acc, _e, _com_vel] = calc_acc_e<T>(entry);
        while (counter++ < 3) {
            std::cout << BOLD BTCOL "Body " BNCOL << +counter << SHCOL "\n------\n" RST;
            if constexpr (pos)
                print_pos(((vector<T>*) &entry) + idx);
            if constexpr (vel)
                print_vel(((vector<T>*) &entry) + 3 + idx);
            print_acc(*acc++);
            ++idx;
        }
        std::cout << HCOL << "----------\n" << ETCOL "KE: " ENCOL BOLD << *_e++ << RST UCOL " joules\n"
                << ETCOL "PE: " ENCOL BOLD << *_e++ << RST UCOL " joules\n"
                << ETCOL "E: " ENCOL BOLD << *_e << RST UCOL " joules\n" << RST;
        if constexpr (COM) {
            std::cout << COMPCOL "COM position: " BOLD NNCOL <<
                      com(*(((vector<T>*) &entry) + 0), *(((vector<T>*) &entry) + 1), *(((vector<T>*) &entry) + 2)) <<
                      UCOL " m\n" RST;
            std::cout << COMPCOL "COM velocity: " BOLD NNCOL << _com_vel << UCOL " m\n" RST;
        }
        return;
    }
    if constexpr (ACC) {
        auto acc = calc_acc<T>(entry);
        while (counter++ < 3) {
            std::cout << BOLD BTCOL "Body " BNCOL << +counter << SHCOL "\n------\n" RST;
            if constexpr (pos)
                print_pos(((vector<T>*) &entry) + idx);
            if constexpr (vel)
                print_vel(((vector<T>*) &entry) + 3 + idx);
            print_acc(*acc++);
            ++idx;
        }
        if constexpr (COM) {
            std::cout << COMPCOL "COM position: " BOLD NNCOL <<
                      com(*(((vector<T>*) &entry) + 0), *(((vector<T>*) &entry) + 1), *(((vector<T>*) &entry) + 2)) <<
                      UCOL " m\n" RST;
            std::cout << COMPCOL "COM velocity: " BOLD NNCOL <<
                      com(*(((vector<T>*) &entry) + 3), *(((vector<T>*) &entry) + 4), *(((vector<T>*) &entry) + 5)) <<
                      UCOL " m\n" RST;
        }
        return;
    }
    if constexpr (ENERGY) {
        auto [_e, _com_vel] = calc_e<T>(entry);
        while (counter++ < 3) {
            std::cout << BOLD BTCOL "Body " BNCOL << +counter << SHCOL "\n------\n" RST;
            if constexpr (pos)
                print_pos(((vector<T>*) &entry) + idx);
            if constexpr (vel)
                print_vel(((vector<T>*) &entry) + 3 + idx);
            ++idx;
        }
        std::cout << HCOL << "----------\n" << ETCOL "KE: " ENCOL BOLD << *_e++ << RST UCOL " joules\n"
                  << ETCOL "PE: " ENCOL BOLD << *_e++ << RST UCOL " joules\n"
                  << ETCOL "E: " ENCOL BOLD << *_e << RST UCOL " joules\n" << RST;
        if constexpr (COM) {
            std::cout << COMPCOL "COM position: " BOLD NNCOL <<
                      com(*(((vector<T>*) &entry) + 0), *(((vector<T>*) &entry) + 1), *(((vector<T>*) &entry) + 2)) <<
                      UCOL " m\n" RST;
            std::cout << COMPCOL "COM velocity: " BOLD NNCOL << _com_vel << UCOL " m\n" RST;
        }
        return;
    }
    while (counter++ < 3) {
        std::cout << BOLD BTCOL "Body " BNCOL << +counter << SHCOL "\n------\n" RST;
        if constexpr (pos)
            print_pos(((vector<T>*) &entry) + idx);
        if constexpr (vel)
            print_vel(((vector<T>*) &entry) + 3 + idx);
        ++idx;
    }
    if constexpr (COM) {
        std::cout << COMPCOL "COM position: " BOLD NNCOL <<
                  com(*(((vector<T>*) &entry) + 0), *(((vector<T>*) &entry) + 1), *(((vector<T>*) &entry) + 2)) <<
                  UCOL " m\n" RST;
        std::cout << COMPCOL "COM velocity: " BOLD NNCOL <<
                  com(*(((vector<T>*) &entry) + 3), *(((vector<T>*) &entry) + 4), *(((vector<T>*) &entry) + 5)) <<
                  UCOL " m\n" RST;
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

template <typename T, bool COM = false, bool ENERGY = false> requires (std::is_floating_point_v<T>)
[[noreturn]] void view_all(const gtd::f3bodr<T> &reader) {
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
    exit(0);
}

template <typename T> requires (std::is_floating_point_v<T>)
[[noreturn]] void view_header(const gtd::f3bodr<T> &reader) {
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
    exit(0);
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

template <typename T, bool COM = false, bool ENERGY = false> requires (std::is_floating_point_v<T>)
[[noreturn]] void view_range(const gtd::f3bodr<T> &reader, const char *action) { // over the CLOSED interval [lo,hi]
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
    exit(0);
}

template <typename T, bool COM = false, bool ENERGY = false> requires (std::is_floating_point_v<T>)
[[noreturn]] void view_entry(const gtd::f3bodr<T> &reader, const char *action) {
    uint64_t idx = to_uint64(action + 2);
    const htype<T> &hdr = reader.header();
    if (idx >= hdr.N) {
        std::cerr << BOLD UDL ERRCOL "Error:" UDLRST ERRTCOL " specified index is out-of-range for a .3bod file with "
        NNCOL << hdr.N << ERRTCOL " entries.\n";
        exit(1);
    }
    char *hyphens = get_hyphens(idx);
    std::cout << ETCOL "Entry " BOLD LNCOL << idx << HCOL << hyphens;
    print_entry<T, COM, ENERGY>(reader.entry_at(idx));
    std::cout.put('\n');
    delete [] hyphens;
    exit(0);
}

template <typename T> requires (std::is_floating_point_v<T>)
[[noreturn]] void distance(const gtd::f3bodr<T> &reader) {
    exit(0);
}

template <typename T> requires (std::is_floating_point_v<T>)
[[noreturn]] void distance_range(const gtd::f3bodr<T> &reader, const char *action) {
    exit(0);
}

template <typename T> requires (std::is_floating_point_v<T>)
[[noreturn]] analysis_range(const gtd::f3bodr<T> &reader, const char *action) {

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
    bool calc_com = parser.get_arg("--com", false);
    bool calc_e = parser.get_arg("--energy", false);
    if (action) {
        if (calc_com || calc_e) {
            masses = new T[3]{};
            gtd::copy(masses, reader.header().masses, 3*sizeof(T));
            if (std::regex_match(action, vh_rgx)) {
                std::cerr << BOLD ERRCOL UDL "Error:" UDLRST ERRTCOL
                             " COM and energy calculations are only performed when viewing entries.\n";
                exit(1);
            } else goto after; // this avoids the overhead of re-checking with `std::regex_match()` for the header
        }
    } else {
        if (calc_com && calc_e)
            view_all<T, true, true>(reader); // note that all these functions don't return
        if (calc_com && !calc_e)
            view_all<T, true, false>(reader);
        if (!calc_com && calc_e)
            view_all<T, false, true>(reader);
        view_all<T, false, false>(reader);
    }
    // bool centre = parser.get_arg("--centre", false);
    if (std::regex_match(action, vh_rgx))
        view_header(reader);
    after:
    if (calc_com && calc_e) {
        if (std::regex_match(action, va_rgx))
            view_all<T, true, true>(reader);
        if (std::regex_match(action, vr_rgx))
            view_range<T, true, true>(reader, action);
        if (std::regex_match(action, ve_rgx))
            view_entry<T, true, true>(reader, action);
    }
    /* if (std::regex_match(action, d_rgx))
        distance(reader);
    if (std::regex_match(action, dr_rgx))
        distance_range(reader, action); */
    std::cerr << BOLD UDL ERRCOL "Error:" UDLRST ERRTCOL " invalid action specified.\n";
    exit(1);
}

int main(int argc, char **argv) {
    if (atexit(free_all)) {
        std::cerr << BOLD ERRCOL UDL "Error:" UDLRST ERRTCOL " \"atexit\" error.\n";
        return 1;
    }
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
        path = fpath;
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
            return 1;
        }
        free_path = true;
    }
    std::cout << "\033[1m\033[37mFile: \033[0m\033[32m\"\033[33m" << fpath << "\033[32m\"\033[0m\n\n";
    // Maybe manually open file here first to check its integrity
    try {
        try {
            gtd::f3bodr<long double> reader{fpath};
            if (free_path) {
                delete [] fpath; // delete `fpath` and set `path` to `nullptr` so that my registered function does not
                path = nullptr; // attempt to free it
            }
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
                if (free_path) {
                    delete [] fpath;
                    path = nullptr;
                }
                perform_action(reader, parser);
            } catch (const gtd::invalid_3bod_T_size&) {
                try {
                    gtd::f3bodr<float> reader{fpath};
                    if (free_path) {
                        delete [] fpath;
                        path = nullptr;
                    }
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
    return 0;
}

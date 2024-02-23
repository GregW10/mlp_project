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
std::regex an_rgx{ANALYSIS};
std::regex anr_rgx{ANALYSIS_RANGE};

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
#define ENTRYTCOL "\033[38;5;174m" // entry text colour
#define ENTRYNCOL "\033[38;5;199m" // entry number colour
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

void *masses{}; // pointer to masses (only used for acc/e/com calculations)
char *path{};

// std::vector<void*> to_free = {masses, path};

template <typename T> requires (std::is_floating_point_v<T>)
void free_masses() {
    delete [] ((T*) masses);
}

void free_path() {
    delete [] path;
}

#pragma pack(push, 1) // not really necessary, since it's guaranteed that `sizeof(T)` >= 4
template <typename T> requires (std::is_floating_point_v<T>)
struct vector { // I am not using my `gtd::vector3D<T>` class in this program as it is not a POD-type
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

template <typename T, typename U> requires (std::is_floating_point_v<T> && std::is_floating_point_v<U>)
vector<T> operator*(U scalar, const vector<T> &vec) {
    return {scalar*vec.x, scalar*vec.y, scalar*vec.z};
}

template <typename T, typename U> requires (std::is_floating_point_v<T> && std::is_floating_point_v<U>)
vector<T> operator*(const vector<T> &vec, U scalar) {
    return {scalar*vec.x, scalar*vec.y, scalar*vec.z};
}

template <typename T, typename U> requires (std::is_floating_point_v<T> && std::is_floating_point_v<U>)
vector<T> operator/(const vector<T> &vec, U scalar) {
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
    return vec.x*vec.x + vec.y*vec.y + vec.z*vec.z;
}

template <typename T> requires (std::is_floating_point_v<T>)
long double mag(const vector<T> &vec) {
    return sqrtl(vec.x*vec.x + vec.y*vec.y + vec.z*vec.z);
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
const vector<T> *calc_acc(const etype<T> &entry, bool ret_dist = false) {
    static vector<T> acc_vals[6];
    vector<T> r12 = *(((vector<T>*) &entry) + 1) - *((vector<T>*) &entry);
    vector<T> r13 = *(((vector<T>*) &entry) + 2) - *((vector<T>*) &entry);
    vector<T> r23 = *(((vector<T>*) &entry) + 2) - *(((vector<T>*) &entry) + 1);
    T r12_magsq = r12*r12;
    T r13_magsq = r13*r13;
    T r23_magsq = r23*r23;
    acc_vals[0] =  gtd::sys::G_SI*(r12*(*(((T*) masses) + 1)/(r12_magsq)) + r13*(*(((T*) masses) + 2)/(r13_magsq)));
    acc_vals[1] = -gtd::sys::G_SI*(r12*(*(((T*) masses))/(r12_magsq)) - r23*(*(((T*) masses) + 2)/(r23_magsq)));
    acc_vals[2] = -gtd::sys::G_SI*(r13*(*(((T*) masses))/(r13_magsq)) + r23*(*(((T*) masses) + 1)/(r23_magsq)));
    if (ret_dist) {
        acc_vals[3] = r12;
        acc_vals[4] = r13;
        acc_vals[5] = r23;
    }
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
    e_vals[0] = 0.5l*(*(((T*) masses) + 0)*mag_sq(*(((vector<T>*) &entry) + 3) - com_vel) +
                      *(((T*) masses) + 1)*mag_sq(*(((vector<T>*) &entry) + 4) - com_vel) +
                      *(((T*) masses) + 2)*mag_sq(*(((vector<T>*) &entry) + 5) - com_vel));
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
            if constexpr (POS)
                print_pos(((vector<T>*) &entry) + idx);
            if constexpr (VEL)
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
            std::cout << COMVCOL "COM velocity: " BOLD NNCOL << _com_vel << UCOL " m\n" RST;
        }
        return;
    }
    if constexpr (ACC) {
        auto acc = calc_acc<T>(entry);
        while (counter++ < 3) {
            std::cout << BOLD BTCOL "Body " BNCOL << +counter << SHCOL "\n------\n" RST;
            if constexpr (POS)
                print_pos(((vector<T>*) &entry) + idx);
            if constexpr (VEL)
                print_vel(((vector<T>*) &entry) + 3 + idx);
            print_acc(*acc++);
            ++idx;
        }
        if constexpr (COM) {
            std::cout << COMPCOL "COM position: " BOLD NNCOL <<
                      com(*(((vector<T>*) &entry) + 0), *(((vector<T>*) &entry) + 1), *(((vector<T>*) &entry) + 2)) <<
                      UCOL " m\n" RST;
            std::cout << COMVCOL "COM velocity: " BOLD NNCOL <<
                      com(*(((vector<T>*) &entry) + 3), *(((vector<T>*) &entry) + 4), *(((vector<T>*) &entry) + 5)) <<
                      UCOL " m\n" RST;
        }
        return;
    }
    if constexpr (ENERGY) {
        auto [_e, _com_vel] = calc_e<T>(entry);
        while (counter++ < 3) {
            std::cout << BOLD BTCOL "Body " BNCOL << +counter << SHCOL "\n------\n" RST;
            if constexpr (POS)
                print_pos(((vector<T>*) &entry) + idx);
            if constexpr (VEL)
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
            std::cout << COMVCOL "COM velocity: " BOLD NNCOL << _com_vel << UCOL " m\n" RST;
        }
        return;
    }
    while (counter++ < 3) {
        std::cout << BOLD BTCOL "Body " BNCOL << +counter << SHCOL "\n------\n" RST;
        if constexpr (POS)
            print_pos(((vector<T>*) &entry) + idx);
        if constexpr (VEL)
            print_vel(((vector<T>*) &entry) + 3 + idx);
        ++idx;
    }
    if constexpr (COM) {
        std::cout << COMPCOL "COM position: " BOLD NNCOL <<
                  com(*(((vector<T>*) &entry) + 0), *(((vector<T>*) &entry) + 1), *(((vector<T>*) &entry) + 2)) <<
                  UCOL " m\n" RST;
        std::cout << COMVCOL "COM velocity: " BOLD NNCOL <<
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

template <typename T, bool POS, bool VEL, bool ACC, bool COM, bool ENERGY> requires (std::is_floating_point_v<T>)
[[noreturn]] void view_all(const gtd::f3bodr<T> &reader) {
    const htype<T> &hdr = reader.header();
    if (!hdr.N)
        log_no_readings();
    char *hyphens = get_hyphens(hdr.N);
    uint64_t index = 0;
    unsigned char counter;
    for (const auto &e : reader) {
        std::cout << ENTRYTCOL "Entry " BOLD LNCOL << index++ << HCOL << hyphens;
        print_entry<T, POS, VEL, ACC, COM, ENERGY>(e);
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
    << hdr.nCORs[2] << "\n" HDFCOL UDL "Total sim. time:" UDLRST " " NNCOL << hdr.dt*(hdr.N ? hdr.N - 1 : 0)
    << UCOL " s\n" RST;
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

template <typename T> requires (std::is_floating_point_v<T>)
std::pair<uint64_t, uint64_t> get_indices(const char *str, const htype<T> &hdr) {
    uint64_t idx1 = to_uint64(str, &str, '-');
    uint64_t idx2 = to_uint64(str + 1); // point to the first digit of the end index
    if (idx1 > idx2) {
        std::cerr << BOLD ERRCOL UDL "Error:" UDLRST ERRTCOL " starting index" NNCOL << idx1 << ERRTCOL
        " in range specified is higher than end index" NNCOL << idx2 << ERRTCOL ".\n" RST;
        exit(1);
    }
    if (idx1 >= hdr.N) {
        std::cerr << BOLD UDL ERRCOL "Error:" UDLRST ERRTCOL " starting index" NNCOL << idx1 << ERRTCOL
        " in range is out-of-bounds for a .3bod file with " NNCOL << hdr.N << ERRTCOL " entries.\n";
        exit(1);
    }
    if (idx2 >= hdr.N) {
        std::cerr << BOLD UDL ERRCOL "Error:" UDLRST ERRTCOL " end index" NNCOL << idx2 << ERRTCOL
        " in range is out-of-bounds for a .3bod file with " NNCOL << hdr.N << ERRTCOL " entries.\n";
        exit(1);
    }
    return {idx1, idx2};
}

template <typename T, bool POS, bool VEL, bool ACC, bool COM, bool ENERGY> requires (std::is_floating_point_v<T>)
[[noreturn]] void view_range(const gtd::f3bodr<T> &reader, const char *action) { // over the CLOSED interval [lo,hi]
    auto [idx1, idx2] = get_indices(action + 2, reader.header());//move past "r:" to point to first digit of start. idx
    char *hyphens = get_hyphens(idx2);
    std::cout << "\033[38;5;6mRange: " BOLD "\033[38;5;7m[" LNCOL << idx1 << "\033[38;5;11m, " LNCOL
              << idx2 << "\033[38;5;7m]\n\n";
    auto it = reader.begin() + idx1;
    while (idx1 <= idx2) {
        std::cout << ENTRYTCOL "Entry " BOLD LNCOL << idx1++ << HCOL << hyphens;
        print_entry<T, POS, VEL, ACC, COM, ENERGY>(*it++);
        std::cout.put('\n');
    }
    delete [] hyphens;
    exit(0);
}

template <typename T, bool POS, bool VEL, bool ACC, bool COM, bool ENERGY> requires (std::is_floating_point_v<T>)
[[noreturn]] void view_entry(const gtd::f3bodr<T> &reader, const char *action) {
    uint64_t idx = to_uint64(action + 2);
    const htype<T> &hdr = reader.header();
    if (idx >= hdr.N) {
        std::cerr << BOLD UDL ERRCOL "Error:" UDLRST ERRTCOL " specified index is out-of-range for a .3bod file with "
        NNCOL << hdr.N << ERRTCOL " entries.\n";
        exit(1);
    }
    char *hyphens = get_hyphens(idx);
    std::cout << ENTRYTCOL "Entry " BOLD LNCOL << idx << HCOL << hyphens;
    print_entry<T, POS, VEL, ACC, COM, ENERGY>(reader.entry_at(idx));
    std::cout.put('\n');
    delete [] hyphens;
    exit(0);
}

template <typename T> requires (std::is_floating_point_v<T>)
[[noreturn]] void perform_analysis(const gtd::f3bodr<T> &reader, uint64_t idx_lo, uint64_t idx_hi) {
    std::cout << "\033[38;5;6mRange: " BOLD "\033[38;5;7m[" LNCOL << idx_lo << "\033[38;5;11m, " LNCOL
              << idx_hi << "\033[38;5;7m]\n\n";
    uint64_t _num = idx_hi - idx_lo + 1;
    std::pair<std::pair<unsigned char, unsigned char>, T> max_sep;
    T max_vels[3]{};
    T max_acc[3]{}; // max_forces will follow from this
    std::pair<std::pair<unsigned char, unsigned char>, T> min_sep = {{}, std::numeric_limits<T>::infinity()};
    T min_vels[3] = {std::numeric_limits<T>::infinity()};
    T mean_seps[3]{};
    T mean_sq_seps[3]{}; // mean of the squares of the separation - used for SD calculation
    // T sep_sd = 0;
    T distances[3]{};
    auto it = reader.begin() + idx_lo;
    uint64_t counter = 0;
    vector<T> *p1, *p2, *p3, *v1, *v2, *v3, pp1{*((vector<T>*) it.operator->() + 0)},
               pp2{*((vector<T>*) it.operator->() + 3)}, pp3{*((vector<T>*) it.operator->() + 6)};
    const vector<T> *acc_sep{};
    T val1, val2, val3;
    while (counter++ < _num) { /*
        p1 = (vector<T>*) it->positions[0];
        p2 = (vector<T>*) it->positions[3];
        p3 = (vector<T>*) it->positions[6];
        v1 = (vector<T>*) it->velocities[0];
        v2 = (vector<T>*) it->velocities[3];
        v3 = (vector<T>*) it->velocities[6]; */
        acc_sep = calc_acc<T>(*it, true);
        if ((val1 = mag(*acc_sep)) > max_acc[0])
            max_acc[0] = val1;
        if ((val2 = mag(*++acc_sep)) > max_acc[1])
            max_acc[1] = val2;
        if ((val3 = mag(*++acc_sep)) > max_acc[2])
            max_acc[2] = val3;
        val1 = mag(*++acc_sep); // |r12|
        val2 = mag(*++acc_sep); // |r13|
        val3 = mag(*++acc_sep); // |r23|
        // mean_sep += val1 + val2 + val3;
        // mean_sq_sep += val1*val1 + val2*val2 + val3*val3;
        mean_seps[0] += val1;
        mean_seps[1] += val2;
        mean_seps[2] += val3;
        mean_sq_seps[0] += val1*val1;
        mean_sq_seps[1] += val2*val2;
        mean_sq_seps[2] += val3*val3;
        distances[0] += mag(*((vector<T>*) (it.operator->() + 0)) - pp1); // add dist between currpos and prevpos to tot.
        distances[1] += mag(*((vector<T>*) (it.operator->() + 3)) - pp2); // redundant zero-addition in first iteration
        distances[2] += mag(*((vector<T>*) (it.operator->() + 6)) - pp3);
        pp1 = *((vector<T>*) (it.operator->() + 0));
        pp2 = *((vector<T>*) (it.operator->() + 3));
        pp3 = *((vector<T>*) (it.operator->() + 6));
        if (val1 > max_sep.second) {
            if (val2 > max_sep.second)
                goto max2;
            if (val3 > max_sep.second)
                goto max3;
            max_sep.first.first = 1;
            max_sep.first.second = 2;
            max_sep.second = val1;
        }
        else if (val2 > max_sep.second) {
            max2:
            if (val3 > max_sep.second)
                goto max3;
            max_sep.first.first = 1;
            max_sep.first.second = 3;
            max_sep.second = val2;
        }
        else if (val3 > max_sep.second) {
            max3:
            max_sep.first.first = 2;
            max_sep.first.second = 3;
            max_sep.second = val3;
        }
        if (val1 < min_sep.second) {
            if (val2 < min_sep.second)
                goto min2;
            if (val3 < min_sep.second)
                goto min3;
            min_sep.first.first = 1;
            min_sep.first.second = 2;
            min_sep.second = val1;
        }
        else if (val2 < min_sep.second) {
            min2:
            if (val3 < min_sep.second)
                goto min3;
            min_sep.first.first = 1;
            min_sep.first.second = 3;
            min_sep.second = val2;
        }
        else if (val3 < min_sep.second) {
            min3:
            min_sep.first.first = 2;
            min_sep.first.second = 3;
            min_sep.second = val3;
        }
        if ((val1 = mag(*((vector<T>*) (it.operator->() + 9)))) > max_vels[0])
            max_vels[0] = val1;
        if ((val2 = mag(*((vector<T>*) (it.operator->() + 12)))) > max_vels[1])
            max_vels[1] = val2;
        if ((val3 = mag(*((vector<T>*) (it.operator->() + 15)))) > max_vels[2])
            max_vels[2] = val3;
        if (val1 < min_vels[0])
            min_vels[0] = val1;
        if (val2 < min_vels[1])
            min_vels[1] = val2;
        if (val3 < min_vels[2])
            min_vels[2] = val3;
        ++it;
        std::cout << "\rProcessed epoch " << counter << '/' << _num;
        std::cout.flush();
    }
    std::cout << "\n----------\n----------\nMax. separation = " << max_sep.second << ", between bodies "
              << +max_sep.first.first << " and " << +max_sep.first.second << "\n----------\n";
    std::cout << "\n----------\nMin. separation = " << min_sep.second << ", between bodies "
              << +min_sep.first.first << " and " << +min_sep.first.second << "\n----------\n";
    char c2;
    for (counter = 0, c2 = 1; counter < 3; ++counter, ++c2)
        std::cout << "Max. velocity of Body " << +c2 << " = " << max_vels[counter] << '\n';
    std::cout << "----------\n";
    for (counter = 0, c2 = 1; counter < 3; ++counter, ++c2)
        std::cout << "Min. velocity of Body " << +c2 << " = " << min_vels[counter] << '\n';
    std::cout << "----------\n";
    for (counter = 0, c2 = 1; counter < 3; ++counter, ++c2)
        std::cout << "Max. acceleration of Body " << +c2 << " = " << max_acc[counter] <<
        "\nMax. force on Body " << +c2 << " = " << max_acc[counter]*(*(((T*) masses) + counter)) << '\n';
    mean_seps[0] /= _num;
    mean_seps[1] /= _num;
    mean_seps[2] /= _num;
    mean_sq_seps[0] /= _num;
    mean_sq_seps[1] /= _num;
    mean_sq_seps[2] /= _num;
    T vars[3];
    vars[0] = mean_sq_seps[0] - mean_seps[0]*mean_seps[0];
    vars[1] = mean_sq_seps[1] - mean_seps[1]*mean_seps[1];
    vars[2] = mean_sq_seps[2] - mean_seps[2]*mean_seps[2];
    std::cout << "----------\nMean sep. between bodies 1 and 2 = " << mean_seps[0] << " m +/- " << sqrtl(vars[0]);
    std::cout << " m\nMean sep. between bodies 1 and 3 = " << mean_seps[1] << " m +/- " << sqrtl(vars[1]);
    std::cout << " m\nMean sep. between bodies 2 and 3 = " << mean_seps[2] << " m +/- " << sqrtl(vars[2]);
    std::cout << " m\nMean sep. = " << (mean_seps[0] + mean_seps[1] + mean_seps[2])/3 << " m +/- "
              << sqrtl(vars[0] + vars[1] + vars[2])/3
              << "\n----------\n";
    for (counter = 0, c2 = 1; counter < 3; ++counter, ++c2)
        std::cout << "Total distance covered by Body " << +c2 << " = " << distances[counter] << '\n';
    std::cout << "Total distance covered by all 3 bodies = " << distances[0] + distances[1] + distances[2] << std::endl;
    exit(0);
}

template <typename T> requires (std::is_floating_point_v<T>)
[[noreturn]] void analysis_range(const gtd::f3bodr<T> &reader, const char *action) {
    while (*action++ != ':');
    auto [idx1, idx2] = get_indices(action, reader.header());
    perform_analysis<T>(reader, idx1, idx2);
}

template <typename T> requires (std::is_floating_point_v<T>)
[[noreturn]] void analysis(const gtd::f3bodr<T> &reader) {
    perform_analysis<T>(reader, 0, reader.header().N - 1);
}

template <typename T> requires (std::is_floating_point_v<T>)
[[noreturn]] void perform_action(const gtd::f3bodr<T> &reader, gtd::parser &parser) {
    std::streamsize prec = parser.get_arg("--precision", (std::streamsize) std::cout.precision());
    if (prec < 0) {
        std::cerr << BOLD ERRCOL UDL "Error:" UDLRST ERRTCOL " precision cannot be negative.\n";
        exit(1);
    }
    std::cout.precision(prec);
#define FUNC_SELECT(func, ...) \
    if (_p && _v && _a && _c && _e) { \
        func<T, true, true, true, true, true>(__VA_ARGS__); \
    } \
    if (_p && _v && _a && _c && !_e) { \
        func<T, true, true, true, true, false>(__VA_ARGS__); \
    } \
    if (_p && _v && _a && !_c && _e) { \
        func<T, true, true, true, false, true>(__VA_ARGS__); \
    } \
    if (_p && _v && _a && !_c && !_e) { \
        func<T, true, true, true, false, false>(__VA_ARGS__); \
    } \
    if (_p && _v && !_a && _c && _e) { \
        func<T, true, true, false, true, true>(__VA_ARGS__); \
    } \
    if (_p && _v && !_a && _c && !_e) { \
        func<T, true, true, false, true, false>(__VA_ARGS__); \
    } \
    if (_p && _v && !_a && !_c && _e) { \
        func<T, true, true, false, false, true>(__VA_ARGS__); \
    } \
    if (_p && _v && !_a && !_c && !_e) { \
        func<T, true, true, false, false, false>(__VA_ARGS__); \
    } \
    if (_p && !_v && _a && _c && _e) { \
        func<T, true, false, true, true, true>(__VA_ARGS__); \
    } \
    if (_p && !_v && _a && _c && !_e) { \
        func<T, true, false, true, true, false>(__VA_ARGS__); \
    } \
    if (_p && !_v && _a && !_c && _e) { \
        func<T, true, false, true, false, true>(__VA_ARGS__); \
    } \
    if (_p && !_v && _a && !_c && !_e) { \
        func<T, true, false, true, false, false>(__VA_ARGS__); \
    } \
    if (_p && !_v && !_a && _c && _e) { \
        func<T, true, false, false, true, true>(__VA_ARGS__); \
    } \
    if (_p && !_v && !_a && _c && !_e) { \
        func<T, true, false, false, true, false>(__VA_ARGS__); \
    } \
    if (_p && !_v && !_a && !_c && _e) { \
        func<T, true, false, false, false, true>(__VA_ARGS__); \
    } \
    if (_p && !_v && !_a && !_c && !_e) { \
        func<T, true, false, false, false, false>(__VA_ARGS__); \
    } \
    if (!_p && _v && _a && _c && _e) { \
        func<T, false, true, true, true, true>(__VA_ARGS__); \
    } \
    if (!_p && _v && _a && _c && !_e) { \
        func<T, false, true, true, true, false>(__VA_ARGS__); \
    } \
    if (!_p && _v && _a && !_c && _e) { \
        func<T, false, true, true, false, true>(__VA_ARGS__); \
    } \
    if (!_p && _v && _a && !_c && !_e) { \
        func<T, false, true, true, false, false>(__VA_ARGS__); \
    } \
    if (!_p && _v && !_a && _c && _e) { \
        func<T, false, true, false, true, true>(__VA_ARGS__); \
    } \
    if (!_p && _v && !_a && _c && !_e) { \
        func<T, false, true, false, true, false>(__VA_ARGS__); \
    } \
    if (!_p && _v && !_a && !_c && _e) { \
        func<T, false, true, false, false, true>(__VA_ARGS__); \
    } \
    if (!_p && _v && !_a && !_c && !_e) { \
        func<T, false, true, false, false, false>(__VA_ARGS__); \
    } \
    if (!_p && !_v && _a && _c && _e) { \
        func<T, false, false, true, true, true>(__VA_ARGS__); \
    } \
    if (!_p && !_v && _a && _c && !_e) { \
        func<T, false, false, true, true, false>(__VA_ARGS__); \
    } \
    if (!_p && !_v && _a && !_c && _e) { \
        func<T, false, false, true, false, true>(__VA_ARGS__); \
    } \
    if (!_p && !_v && _a && !_c && !_e) { \
        func<T, false, false, true, false, false>(__VA_ARGS__); \
    } \
    if (!_p && !_v && !_a && _c && _e) { \
        func<T, false, false, false, true, true>(__VA_ARGS__); \
    } \
    if (!_p && !_v && !_a && _c && !_e) { \
        func<T, false, false, false, true, false>(__VA_ARGS__); \
    } \
    if (!_p && !_v && !_a && !_c && _e) { \
        func<T, false, false, false, false, true>(__VA_ARGS__); \
    } \
    if (!_p && !_v && !_a && !_c && !_e) { \
        func<T, false, false, false, false, false>(__VA_ARGS__); \
    }
    bool _p, _v, _a, _c, _e;
    const char *action = nullptr;
    if ((action = parser.get_arg(vh_rgx)))
        view_header(reader); // function doesn't return
    if ((action = parser.get_arg(va_rgx))) {
#define WITH_ACTION \
        if (!parser.empty()) {\
            std::cerr << BOLD ERRCOL UDL "Error:" UDLRST ERRTCOL " only one action can be specified.\n";\
            exit(1);\
        }\
        action += 2; \
        _p = gtd::contains(action, 'p'); \
        _v = gtd::contains(action, 'v'); \
        _a = gtd::contains(action, 'a'); \
        _c = gtd::contains(action, 'c'); \
        _e = gtd::contains(action, 'e'); \
        if (_a || _c || _e) { \
            masses = new T[3]{}; \
            gtd::copy(masses, reader.header().masses, 3*sizeof(T)); \
        }
        WITH_ACTION
        FUNC_SELECT(view_all, reader) // function doesn't return
    }
    if ((action = parser.get_arg(vr_rgx))) {
        WITH_ACTION
        FUNC_SELECT(view_range, reader, action)
    }
    if ((action = parser.get_arg(ve_rgx))) {
        WITH_ACTION
        FUNC_SELECT(view_entry, reader, action)
    }
    if ((action = parser.get_arg(an_rgx))) {
        // WITH_ACTION
        // FUNC_SELECT(analysis, reader)
        analysis<T>(reader);
    }
    if ((action = parser.get_arg(anr_rgx))) {
        // WITH_ACTION
        // FUNC_SELECT(analysis_range, reader, action)
        analysis<T>(reader, action);
    }
    std::cerr << BOLD UDL ERRCOL "Error:" UDLRST ERRTCOL " invalid action specified.\n";
    exit(1);
}

int main(int argc, char **argv) {
    if (atexit(free_path)) {
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
            if (atexit(free_masses<long double>)) {
                std::cerr << BOLD ERRCOL UDL "Error:" UDLRST ERRTCOL " \"atexit\" error.\n";
                return 1;
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
                if (atexit(free_masses<double>)) {
                    std::cerr << BOLD ERRCOL UDL "Error:" UDLRST ERRTCOL " \"atexit\" error.\n";
                    return 1;
                }
                perform_action(reader, parser);
            } catch (const gtd::invalid_3bod_T_size&) {
                try {
                    gtd::f3bodr<float> reader{fpath};
                    if (free_path) {
                        delete [] fpath;
                        path = nullptr;
                    }
                    if (atexit(free_masses<float>)) {
                        std::cerr << BOLD ERRCOL UDL "Error:" UDLRST ERRTCOL " \"atexit\" error.\n";
                        return 1;
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

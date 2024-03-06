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
#define VIEW_ALL R"(^v-(?:([pvace])(?!.*\1)){0,5}$)" // this regex took me AGES to craft!! ... simple as it looks :(
#define VIEW_RANGE R"(^v-(?:([pvace])(?!.*\1)){0,5}:\d{1,18}-\d{1,18}$)"
#define VIEW_ENTRY R"(^v-(?:([pvace])(?!.*\1)){0,5}:\d{1,18}$)"
#define ANALYSIS R"(^a$)"
#define ANALYSIS_RANGE R"(^a:\d{1,18}-\d{1,18}$)"

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
#define ENTRYTCOL "\033[38;5;10m" // entry text colour
#define ENTRYNCOL "\033[38;5;226m" // entry number colour
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
#define ANFCOL "\033[38;5;27m" // analysis field colour
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
/* NOTE: all the `+ 0` additions below are simply for visual clarity, as the reduction in performance is minimal. */
template <typename T> requires (std::is_floating_point_v<T>)
const vector<T> *calc_acc(const etype<T> &entry, bool ret_dist = false) {
    static vector<T> acc_vals[6];
    vector<T> r12 = *(((vector<T>*) &entry) + 1) - *((vector<T>*) &entry);
    vector<T> r13 = *(((vector<T>*) &entry) + 2) - *((vector<T>*) &entry);
    vector<T> r23 = *(((vector<T>*) &entry) + 2) - *(((vector<T>*) &entry) + 1);
    T r12_magsq = r12*r12;
    T r13_magsq = r13*r13;
    T r23_magsq = r23*r23;
    acc_vals[0] =  gtd::sys::G_SI*(r12*(*(((T*) masses) + 1)/powl(r12_magsq, 1.5l)) + r13*(*(((T*) masses) + 2)/powl(r13_magsq, 1.5l)));
    acc_vals[1] = -gtd::sys::G_SI*(r12*(*(((T*) masses) + 0)/powl(r12_magsq, 1.5l)) - r23*(*(((T*) masses) + 2)/powl(r23_magsq, 1.5l)));
    acc_vals[2] = -gtd::sys::G_SI*(r13*(*(((T*) masses) + 0)/powl(r13_magsq, 1.5l)) + r23*(*(((T*) masses) + 1)/powl(r23_magsq, 1.5l)));
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
    e_vals[0] = 0.5l*(*(((T*) masses) + 0)*mag_sq(*(((vector<T>*) &entry) + 3) - com_vel) +
                      *(((T*) masses) + 1)*mag_sq(*(((vector<T>*) &entry) + 4) - com_vel) +
                      *(((T*) masses) + 2)*mag_sq(*(((vector<T>*) &entry) + 5) - com_vel));
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
    acc_vals[0] =  gtd::sys::G_SI*(r12*(*(((T*) masses) + 1)/powl(r12_magsq, 1.5l)) + r13*(*(((T*) masses) + 2)/powl(r13_magsq, 1.5l)));
    acc_vals[1] = -gtd::sys::G_SI*(r12*(*(((T*) masses) + 0)/powl(r12_magsq, 1.5l)) - r23*(*(((T*) masses) + 2)/powl(r23_magsq, 1.5l)));
    acc_vals[2] = -gtd::sys::G_SI*(r13*(*(((T*) masses) + 0)/powl(r13_magsq, 1.5l)) + r23*(*(((T*) masses) + 1)/powl(r23_magsq, 1.5l)));
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
            print_acc(acc++);
            ++idx;
        }
        if constexpr (COM) {
            std::cout << COMPCOL "COM position: " BOLD NNCOL <<
                      com(*(((vector<T>*) &entry) + 0), *(((vector<T>*) &entry) + 1), *(((vector<T>*) &entry) + 2)) <<
                      UCOL BOLD " m\n" RST;
            std::cout << COMVCOL "COM velocity: " BOLD NNCOL << _com_vel << UCOL BOLD " m/s\n" RST;
        }
        std::cout << HCOL << "----------\n" UDL CYAN_TXT_START "In COM frame:\n" UDLRST << ETCOL "\tKE: " ENCOL BOLD
                  << *_e++ << UCOL " joules\n" RST
                  << ETCOL "\tPE: " ENCOL BOLD << *_e++ << UCOL " joules\n" RST
                  << ETCOL "\tE: " ENCOL BOLD << *_e << UCOL " joules\n" << RST;
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
            print_acc(acc++);
            ++idx;
        }
        if constexpr (COM) {
            std::cout << COMPCOL "COM position: " BOLD NNCOL <<
                      com(*(((vector<T>*) &entry) + 0), *(((vector<T>*) &entry) + 1), *(((vector<T>*) &entry) + 2)) <<
                      UCOL BOLD " m\n" RST;
            std::cout << COMVCOL "COM velocity: " BOLD NNCOL <<
                      com(*(((vector<T>*) &entry) + 3), *(((vector<T>*) &entry) + 4), *(((vector<T>*) &entry) + 5)) <<
                      UCOL BOLD " m/s\n" RST;
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
        std::cout << HCOL << "----------\n" UDL CYAN_TXT_START "In COM frame:\n" UDLRST << ETCOL "\tKE: " ENCOL BOLD
                  << *_e++ << UCOL " joules\n" RST
                  << ETCOL "\tPE: " ENCOL BOLD << *_e++ << UCOL " joules\n" RST
                  << ETCOL "\tE: " ENCOL BOLD << *_e << UCOL " joules\n" << RST;
        if constexpr (COM) {
            std::cout << COMPCOL "COM position: " BOLD NNCOL <<
                      com(*(((vector<T>*) &entry) + 0), *(((vector<T>*) &entry) + 1), *(((vector<T>*) &entry) + 2)) <<
                      UCOL BOLD " m\n" RST;
            std::cout << COMVCOL "COM velocity: " BOLD NNCOL << _com_vel << UCOL BOLD " m/s\n" RST;
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
                  UCOL BOLD " m\n" RST;
        std::cout << COMVCOL "COM velocity: " BOLD NNCOL <<
                  com(*(((vector<T>*) &entry) + 3), *(((vector<T>*) &entry) + 4), *(((vector<T>*) &entry) + 5)) <<
                  UCOL BOLD " m/s\n" RST;
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
    while (*action++ != ':');
    auto [idx1, idx2] = get_indices<T>(action, reader.header());//move past ':' to point to first digit of start. index
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
    while (*action++ != ':');
    uint64_t idx = to_uint64(action);
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
    T inf = std::numeric_limits<T>::infinity();
    std::cout << "\033[38;5;6mRange: " BOLD "\033[38;5;7m[" LNCOL << idx_lo << "\033[38;5;11m, " LNCOL
              << idx_hi << "\033[38;5;7m]\n\n" RST;
    uint64_t _num = idx_hi - idx_lo + 1;
    std::tuple<T, std::pair<unsigned char, unsigned char>, uint64_t> max_sep; // max. sep, body numbers, epoch number
    std::pair<T, uint64_t> max_vels[3] = {{0, (uint64_t) -1}, {0, (uint64_t) -1}, {0, (uint64_t) -1}};
    std::pair<T, uint64_t> max_acc[3] = {{0, (uint64_t) -1}, {0, (uint64_t) -1}, {0, (uint64_t) -1}};//max. forces after
    std::tuple<T, std::pair<unsigned char, unsigned char>, uint64_t> min_sep = {inf, {}, (uint64_t) -1};
    std::pair<T, uint64_t> min_vels[3] = {{inf, (uint64_t) -1}, {inf, (uint64_t) -1}, {inf, (uint64_t) -1}};
    T mean_seps[3]{};
    T mean_sq_seps[3]{}; // mean of the squares of the separation - used for SD calculation
    T distances[3]{};
    auto it = reader.begin() + idx_lo;
    uint64_t counter = idx_lo;
    vector<T> pp1{*(((vector<T>*) it->positions) + 0)}, pp2{*(((vector<T>*) it->positions) + 1)},
              pp3{*(((vector<T>*) it->positions) + 2)};
    const vector<T> *acc_sep{};
    T val1, val2, val3;
    while (counter <= idx_hi) {
        acc_sep = calc_acc<T>(*it, true);
        if ((val1 = mag(*acc_sep)) > max_acc[0].first) {
            max_acc[0].first = val1;
            max_acc[0].second = counter;
        }
        if ((val2 = mag(*++acc_sep)) > max_acc[1].first) {
            max_acc[1].first = val2;
            max_acc[1].second = counter;
        }
        if ((val3 = mag(*++acc_sep)) > max_acc[2].first) {
            max_acc[2].first = val3;
            max_acc[2].second = counter;
        }
        val1 = mag(*++acc_sep); // |r12|
        val2 = mag(*++acc_sep); // |r13|
        val3 = mag(*++acc_sep); // |r23|
        mean_seps[0] += val1;
        mean_seps[1] += val2;
        mean_seps[2] += val3;
        mean_sq_seps[0] += val1*val1;
        mean_sq_seps[1] += val2*val2;
        mean_sq_seps[2] += val3*val3;
        distances[0] += mag(*(((vector<T>*) it->positions) + 0) - pp1); // add dist between currpos and prevpos to tot.
        distances[1] += mag(*(((vector<T>*) it->positions) + 1) - pp2); // redundant zero-addition in first iteration
        distances[2] += mag(*(((vector<T>*) it->positions) + 2) - pp3);
        pp1 = *(((vector<T>*) it->positions) + 0);
        pp2 = *(((vector<T>*) it->positions) + 1);
        pp3 = *(((vector<T>*) it->positions) + 2);
        if (val1 > std::get<0>(max_sep)) {
            if (val2 > std::get<0>(max_sep))
                goto max2;
            if (val3 > std::get<0>(max_sep))
                goto max3;
            std::get<1>(max_sep).first = 1;
            std::get<1>(max_sep).second = 2;
            std::get<0>(max_sep) = val1;
            std::get<2>(max_sep) = counter;
        }
        else if (val2 > std::get<0>(max_sep)) {
            max2:
            if (val3 > std::get<0>(max_sep))
                goto max3;
            std::get<1>(max_sep).first = 1;
            std::get<1>(max_sep).second = 3;
            std::get<0>(max_sep) = val2;
            std::get<2>(max_sep) = counter;
        }
        else if (val3 > std::get<0>(max_sep)) {
            max3:
            std::get<1>(max_sep).first = 2;
            std::get<1>(max_sep).second = 3;
            std::get<0>(max_sep) = val3;
            std::get<2>(max_sep) = counter;
        }
        if (val1 < std::get<0>(min_sep)) {
            if (val2 < std::get<0>(min_sep))
                goto min2;
            if (val3 < std::get<0>(min_sep))
                goto min3;
            std::get<1>(min_sep).first = 1;
            std::get<1>(min_sep).second = 2;
            std::get<0>(min_sep) = val1;
            std::get<2>(min_sep) = counter;
        }
        else if (val2 < std::get<0>(min_sep)) {
            min2:
            if (val3 < std::get<0>(min_sep))
                goto min3;
            std::get<1>(min_sep).first = 1;
            std::get<1>(min_sep).second = 3;
            std::get<0>(min_sep) = val2;
            std::get<2>(min_sep) = counter;
        }
        else if (val3 < std::get<0>(min_sep)) {
            min3:
            std::get<1>(min_sep).first = 2;
            std::get<1>(min_sep).second = 3;
            std::get<0>(min_sep) = val3;
            std::get<2>(min_sep) = counter;
        }
        if ((val1 = mag(*(((vector<T>*) it->velocities) + 0))) > max_vels[0].first) {
            max_vels[0].first = val1;
            max_vels[0].second = counter;
        }
        if ((val2 = mag(*(((vector<T>*) it->velocities) + 1))) > max_vels[1].first) {
            max_vels[1].first = val2;
            max_vels[1].second = counter;
        }
        if ((val3 = mag(*(((vector<T>*) it->velocities) + 2))) > max_vels[2].first) {
            max_vels[2].first = val3;
            max_vels[2].second = counter;
        }
        if (val1 < min_vels[0].first) {
            min_vels[0].first = val1;
            min_vels[0].second = counter;
        }
        if (val2 < min_vels[1].first) {
            min_vels[1].first = val2;
            min_vels[1].second = counter;
        }
        if (val3 < min_vels[2].first) {
            min_vels[2].first = val3;
            min_vels[2].second = counter;
        }
        ++it;
        std::cout << RST ICHCOL "\rProcessed epoch \033[38;5;214m" BOLD << counter++ << "\033[38;5;226m/\033[38;5;10m" << idx_hi;
        std::cout.flush();
    }
#define EQUALS "\033[38;5;51m = " NNCOL
    std::cout << HCOL "\n----------\n" ANFCOL "Max. separation" EQUALS << std::get<0>(max_sep)
              << UCOL " m" ANFCOL ", between bodies " BNCOL
              << +std::get<1>(max_sep).first << ANFCOL " and " BNCOL
              << +std::get<1>(max_sep).second << ANFCOL " at " ENTRYTCOL "epoch " ENTRYNCOL << std::get<2>(max_sep)
              << HCOL "\n----------\n" ANFCOL "Min. separation" EQUALS
              << std::get<0>(min_sep) << UCOL " m" ANFCOL ", between bodies " BNCOL
              << +std::get<1>(min_sep).first << ANFCOL " and " BNCOL << +std::get<1>(min_sep).second
              << ANFCOL " at " ENTRYTCOL "epoch " ENTRYNCOL << std::get<2>(min_sep) << HCOL "\n----------";
    char c2;
    for (counter = 0, c2 = 1; counter < 3; ++counter, ++c2)
        std::cout << ANFCOL "\nMax. velocity of Body " BNCOL << +c2 << EQUALS << max_vels[counter].first
                  << UCOL " m/s " ANFCOL "at " ENTRYTCOL "epoch " ENTRYNCOL << max_vels[counter].second;
    std::cout << HCOL "\n----------";
    for (counter = 0, c2 = 1; counter < 3; ++counter, ++c2)
        std::cout << ANFCOL "\nMin. velocity of Body " BNCOL << +c2 << EQUALS << min_vels[counter].first
                  << UCOL " m/s " ANFCOL "at" ENTRYTCOL " epoch " ENTRYNCOL << min_vels[counter].second;
    std::cout << HCOL "\n----------";
    for (counter = 0, c2 = 1; counter < 3; ++counter, ++c2)
        std::cout << ANFCOL "\nMax. acceleration of Body " BNCOL << +c2 << EQUALS << max_acc[counter].first
                  << UCOL " m/s^2" ANFCOL "\nMax. force on Body " BNCOL << +c2 << EQUALS
                  << max_acc[counter].first*(*(((T*) masses) + counter))
                  << UCOL " N " ANFCOL "at" ENTRYTCOL " epoch " ENTRYNCOL << max_acc[counter].second << HCOL "\n----";
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
    std::cout << "------\n" ANFCOL "Mean sep. between bodies " BNCOL "1" ANFCOL " and " BNCOL "2" EQUALS << mean_seps[0]
              << UCOL " m " PMCOL "+/- " NNCOL << sqrtl(vars[0])
              << UCOL " m" ANFCOL "\nMean sep. between bodies " BNCOL "1" ANFCOL " and " BNCOL "3" EQUALS
              << mean_seps[1] << UCOL " m " PMCOL "+/- " NNCOL << sqrtl(vars[1])
              << UCOL " m" ANFCOL "\nMean sep. between bodies " BNCOL "2" ANFCOL " and " BNCOL "3" EQUALS
              << mean_seps[2] << UCOL " m " PMCOL "+/- " NNCOL << sqrtl(vars[2])
              << UCOL " m" HCOL "\n----\n" ANFCOL "Mean sep." EQUALS << (mean_seps[0] + mean_seps[1] + mean_seps[2])/3
              << UCOL " m " PMCOL "+/- " NNCOL << sqrtl(vars[0] + vars[1] + vars[2])/3
              << UCOL " m" HCOL "\n----------\n";
    for (counter = 0, c2 = 1; counter < 3; ++counter, ++c2)
        std::cout << ANFCOL "Total distance covered by Body " BNCOL << +c2 << EQUALS << distances[counter]
                  << UCOL " m\n";
    std::cout << HCOL "----\n" ANFCOL "Total distance covered by all 3 bodies" EQUALS
              << distances[0] + distances[1] + distances[2] << UCOL " m" HCOL "\n----------" << std::endl;
    exit(0);
}

template <typename T> requires (std::is_floating_point_v<T>)
[[noreturn]] void analysis_range(const gtd::f3bodr<T> &reader, const char *action) {
    while (*action++ != ':');
    auto [idx1, idx2] = get_indices<T>(action, reader.header());
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
    if (parser.empty())
        view_all<T, true, true, false, false, false>(reader);
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
        masses = new T[3]{};
        gtd::copy(masses, reader.header().masses, 3*sizeof(T));
        analysis<T>(reader);
    }
    if ((action = parser.get_arg(anr_rgx))) {
        // WITH_ACTION
        // FUNC_SELECT(analysis_range, reader, action)
        masses = new T[3]{};
        gtd::copy(masses, reader.header().masses, 3*sizeof(T));
        analysis_range<T>(reader, action);
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

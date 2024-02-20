#include "datsup.hpp"
#include <dirent.h>

#ifdef __linux__ // necessary since Linux and MacOS have diff. names for the `struct timespec` fields in `struct stat`
#define BUFF_SEC buff.st_mtim.tv_sec
#define BUFF_NSEC buff.st_mtim.tv_nsec
#elif defined(__APPLE__)
#define BUFF_SEC buff.st_mtimespec.tv_sec
#define BUFF_NSEC buff.st_mtimespec.tv_nsec
#endif

int main(int argc, char **argv) {
    gtd::parser parser{argc, argv};
    char *fpath = const_cast<char*>(parser.get_arg("--file"));
    bool free_path = false;
    if (!fpath)
        fpath = const_cast<char*>(parser.get_arg("-f"));
    if (!fpath) { // if no path to a .3bod file is given, I take the latest modified one found using directory entries
        DIR *dir = opendir(".");
        if (!dir) {
            std::cerr << "\033[1m\033[35mError: \033[0m\033[32mno path to .3bod file provided and could not open "
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
            std::cerr << "\033[1m\033[34mError: \033[0m\033[36m no .3bod files found in current directory.\n";
            delete [] fpath;
            return 1;
        }
        free_path = true;
    }
    std::cout << "\033[1m\033[37mFile: \033[0m\033[32m\"\033[33m" << fpath << "\033[32m\"" << std::endl;
    // Maybe manually open file here first to check its integrity
    try {
        try {
            gtd::f3bodr<long double> reader{fpath};
            // print entries
        }
        catch (const gtd::invalid_3bod_ld_size&) {
            std::cerr << BOLD_TXT(MAGENTA_TXT("Error: ")) YELLOW_TXT(" \"sizeof(long double)\" does not match that "
                                                                     "reported in the .3bod file.\n");
            return 1;
        }
        catch (const gtd::invalid_3bod_T_size&) {
            try {
                gtd::f3bodr<double> reader{fpath};
            } catch (const gtd::invalid_3bod_T_size&) {
                try {
                    gtd::f3bodr<float> reader{fpath};
                } catch (const gtd::invalid_3bod_T_size&) {
                    std::cerr << "Error: reported floating point data type size does not match the size of a "
                                 "\"long double\", \"double\" or \"float\".\n";
                    return 1;
                }
            }
        }
    }
    catch (const gtd::invalid_3bod_format &e) {
        std::cerr << BOLD_TXT(MAGENTA_TXT("Error: ")) YELLOW_TXT("invalid .3bod format.\n");
        return 1;
    }
    if (free_path)
        delete [] fpath;
    return 0;
}

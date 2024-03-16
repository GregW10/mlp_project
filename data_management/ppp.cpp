#include "../glib/misc/gregmisc.hpp"
#include "../glib/misc/gregparse.hpp"
#include "../glib/ml/greggen.hpp"
#include <iostream>
#include <filesystem>
#include <memory>

#define RST "\033[0m"
#define BOLD "\033[1m"
#define UDL "\033[4m"
#define UDLRST "\033[24m"
#define BLINK "\033[5m"
#define BLKRST "\033[25m"
#define ERRCOL "\033[38;5;9m"
#define ERRTCOL "\033[38;5;11m"
#define QUOTECOL "\033[38;5;10m"
#define PATHCOL "\033[38;5;21m"
#define B3F "\033[1m\033[38;5;87m.3bod\033[0m"
#define INFCOL "\033[38;5;33m"
#define NCOL "\033[38;5;3m"
#define BYTES "\033[38;5;204mbytes\033[0m"
#define URCARGCOL "\033[38;5;1m"

int main(int argc, char **argv) {
    gtd::parser parser{argc, argv};
    const char *ddir = parser.get_arg("--data_dir");
    if (!ddir)
        ddir = ".";
    const char *odir = parser.get_arg("-o");
    if (!parser.empty()) {
        std::cerr << ERRTCOL "Unrecognised arguments:\n" BOLD URCARGCOL;
        for (const auto &[_, arg] : parser)
            fprintf(stderr, "\t%s\n", arg.c_str());
        std::cerr << RST;
        return 1;
    }
    std::unique_ptr<char[]> optr;
    if (!odir) {
        optr.reset(gml::gen::now_str("gathered3bods_", ""));
        odir = optr.get();
    }
    if (mkdir(odir, S_IRWXU | S_IRWXG | S_IRWXO) == -1) {
        std::cerr << BOLD ERRCOL "Error: " RST ERRTCOL "could not create directory " QUOTECOL "\"" PATHCOL BOLD << odir
                  << RST QUOTECOL "\"" ERRTCOL ".\n" RST;
        return 1;
    }
    std::unique_ptr<std::vector<std::string>> files{};
    try {
        files.reset(gtd::find_files(ddir, ".3bod"));
    } catch (const std::exception &_e) {
        std::cerr << BOLD ERRTCOL << _e.what() << RST << '\n';
        return 1;
    }
    if (!files || files->empty()) {
        std::cerr << BOLD ERRCOL "Error: " RST ERRTCOL "could find any " B3F ERRTCOL " files within directory "
        QUOTECOL "\"" PATHCOL BOLD << ddir << RST QUOTECOL "\"" ERRTCOL ".\n" RST;
        return 1;
    }
    std::string opath = odir;
    optr.reset();
    if (opath.back() != '/')
        opath.push_back('/');
    struct stat buff{};
    uint64_t tot_size = 0;
    std::string::size_type odir_len = opath.size();
    uint64_t counter = 0;
    uint64_t num_files = files->size();
    for (const std::string &file : *files) {
        if (stat(file.c_str(), &buff) == -1) {
            std::cerr << BOLD ERRCOL "Error: " RST ERRTCOL "could not obtain file info for " QUOTECOL "\"" PATHCOL BOLD
                      << file << RST QUOTECOL "\"" ERRTCOL ".\n" RST;
            return 1;
        }
        if (!S_ISREG(buff.st_mode)) {
            std::cerr << BOLD ERRCOL "Error: " RST ERRTCOL "file " QUOTECOL "\"" PATHCOL BOLD
                      << file << RST QUOTECOL "\"" ERRTCOL " is not a regular file.\n" RST;
            return 1;
        }
        opath += std::to_string(counter++);
        opath.push_back('_');
        opath.append(file.begin() + file.rfind('/') + 1, file.end());
        std::filesystem::copy_file(file, opath);
        opath.erase(odir_len);
        tot_size += buff.st_size;
        std::cout << INFCOL "\rCopied " BOLD NCOL << counter << RST INFCOL "/" BOLD NCOL << num_files << RST INFCOL
        " files" << std::flush;
    }
    std::cout << INFCOL "\nCopied a total of " BOLD NCOL BLINK << num_files << RST INFCOL " files for a total of "
    BOLD NCOL BLINK << tot_size << BLKRST " " BYTES INFCOL " to directory " QUOTECOL "\"" PATHCOL BOLD << opath <<
    RST QUOTECOL "\"" INFCOL ".\n";
    return 0;
}

#include "../glib/ml/gregffnn.hpp"
#include "../glib/misc/gregparse.hpp"
#include <memory>

int main(int argc, char **argv) {
    gtd::parser parser{argc, argv};
    const char *fname = parser.get_arg("-f");
    if (!fname) {
        std::cerr << "Error: filename must be specified.\n";
        return 1;
    }
    std::unique_ptr<char[]> out_name{};
    const char *fout = parser.get_arg("-o");
    if (!fout) {
        char *ptr;
        out_name.reset((ptr = new char[9 + gtd::strlen_c(fname) + 1]));
        // fout = out_name.get();
        gtd::strcpy_c(ptr, "id_added_");
        gtd::strcpy_c(ptr + 9, fname);
        fout = ptr;
    }
    uint32_t id = parser.get_arg("--id", (uint32_t) 1);
    try {
        try {
            gml::ffnn<long double> ffnn{fname};
        } catch (const gml::exceptions::invalid_nnw_format&) {
            try {
                gml::ffnn<double> ffnn{fname};
            } catch (const gml::exceptions::invalid_nnw_format&) {
                try {
                    gml::ffnn<float> ffnn{fname};
                } catch (const gml::exceptions::invalid_nnw_format &e) {
                    struct stat buff{};
                    if (stat(fname, &buff) == -1) {
                        std::cerr << "Error: could not obtain file info for \"" << fname << "\".\n";
                        return 1;
                    }
                    if (!S_ISREG(buff.st_mode)) {
                        std::cerr << "Error: file \"" << fname << "\" is not a regular file.\n";
                        return 1;
                    }
                    if (buff.st_size < 4) {
                        std::cerr << "Error: the size of file \"" << fname << "\" is under 4 bytes.\n";
                        return 1;
                    }
                    std::ifstream in{fname, std::ios_base::in | std::ios_base::binary};
                    if (!in.good()) {
                        std::cerr << "Error: could not open input file \"" << fname << "\".\n";
                        return 1;
                    }
                    std::ofstream out{fout, std::ios_base::out | std::ios_base::binary};
                    if (!out.good()) {
                        std::cerr << "Error: could not open output file \"" << fout << "\".\n";
                        return 1;
                    }
                    char *buffer = new char[buff.st_size];
                    in.read(buffer, buff.st_size);
                    in.close();
                    out.write(buffer, 4);
                    out.write((char *) &id, sizeof(uint32_t));
                    out.write(buffer + 4, buff.st_size - 4);
                    delete [] buffer;
                    out.close();
                    std::cout << "Output file \"" << fout << "\" written.\n";
                    return 0;
                }
            }
        }
    } catch (const std::exception &_e) {
        std::cerr << "Error: some other exception occurred.\nwhat(): " << _e.what();
        return 1;
    }
    std::cout << ".nnw file \"" << fname << "\" appears to already be in order.\n";
    return 0;
}

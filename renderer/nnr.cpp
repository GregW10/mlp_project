#include "nnrsup.hpp"

template <typename T>
using nvf = typename gtd::normaliser<T>::normv_file;

int main(int argc, char **argv) {
    gtd::parser parser{argc, argv};
    const char *nnw = parser.get_arg("--from_model");
    if (!nnw) {
        std::cerr << "Error: model weights must be passed.\n";
        return 1;
    }
    const char *normv = parser.get_arg("--normv_file");
    if (!normv) {
        std::cerr << "Error: .normv file must be provided.\n";
        return 1;
    }
    struct stat buff{};
    if (stat(normv, &buff) == -1) {
        std::cerr << "Error: could not obtain info on file \"" << normv << "\".\n";
        return 1;
    }
    if (buff.st_size < sizeof(nvf<float>)) {
        std::cerr << "Error: file size of \"" << normv
                  << "\" (" << buff.st_size << " bytes) is smaller than size of smallest possible .normv file ("
                  << sizeof(nvf<float>) << " bytes).\n";
        return 1;
    }
    int fd = open(normv, O_RDONLY);
    if (fd == -1) {
        std::cerr << "Error: could not open \"" << normv << "\" for reading.\n";
        return 1;
    }
    char *ptr{};
    std::unique_ptr<char[]> buffer{(ptr = new char[buff.st_size])};
    read(fd, ptr, buff.st_size);
    close(fd);
    if (*ptr != 'N' || *(ptr + 1) != 'O' || *(ptr + 2) != 'R' || *(ptr + 3) != 'M') {
        std::cerr << "Error: invalid file header.\n";
        return 1;
    }
    uint32_t _st;
    uint32_t _sl;
    gtd::memcopy(&_st, ptr + 4, sizeof(uint32_t));
    gtd::memcopy(&_sl, ptr + 8, sizeof(uint32_t));
    return 0;
}

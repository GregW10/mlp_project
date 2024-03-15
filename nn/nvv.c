#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>
#include <dirent.h>
#include <sys/stat.h>

typedef enum boolean {
    false, true
} bool;

#pragma pack(push, 1)
typedef struct normv_file_ld {
    char header[4];
    uint32_t _sizeof_T;
    uint32_t _sizeof_ld;
    // uint32_t _sizeof_u64;
    uint64_t _norm_id;
    long double _maxm;
    long double _minm;
    long double _maxt;
    long double _maxp;
    long double _maxv;
} nvf_ld;

typedef struct normv_file_d {
    char header[4];
    uint32_t _sizeof_T;
    uint32_t _sizeof_ld;
    // uint32_t _sizeof_u64;
    uint64_t _norm_id;
    double _maxm;
    double _minm;
    long double _maxt;
    double _maxp;
    double _maxv;
} nvf_d;

typedef struct normv_file_f {
    char header[4];
    uint32_t _sizeof_T;
    uint32_t _sizeof_ld;
    // uint32_t _sizeof_u64;
    uint64_t _norm_id;
    float _maxm;
    float _minm;
    long double _maxt;
    float _maxp;
    float _maxv;
} nvf_f;
#pragma pack(pop)

uint64_t strlen_c(const char *str) {
    if (!str)
        return -1;
    uint64_t _len = 0;
    while (*str++) ++_len;
    return _len;
}

bool endswith(const char *str, const char *with) {
    if (!str || !with)
        abort();
    uint64_t _slen = strlen_c(str);
    uint64_t _wlen = strlen_c(with);
    if (_wlen > _slen)
        return false;
    str += _slen - _wlen;
    while (*str)
        if (*str++ != *with++)
            return false;
    return true;
}

void strcpy_c(char *dst, const char *src) {
    if (!dst || !src)
        return;
    while (*src)
        *dst++ = *src++;
    *dst = 0;
}

char *find_normv(void) {
    struct dirent *entry;
    DIR *dir;
    if (!(dir = opendir("."))) {
        fprintf(stderr, "Error: no .normv file provided and could not open current directory stream.\n");
        exit(1);
    }
    while ((entry = readdir(dir))) {
        if (endswith(entry->d_name, ".normv")) {
            char *_str = malloc(sizeof(char)*(strlen_c(entry->d_name) + 1));
            if (!_str) {
                fprintf(stderr, "Error: malloc error.\n");
                exit(1);
            }
            strcpy_c(_str, entry->d_name);
            return _str;
        }
    }
    return NULL;
}

int process_normv(const char *path) {
    if (!path)
        return 1;
    printf("File: \"%s\"\n", path);
    struct stat buff = {0};
    if (stat(path, &buff) == -1) {
        fprintf(stderr, "Error: could not obtain file info.\n");
        return 1;
    }
    if (!S_ISREG(buff.st_mode)) {
        fprintf(stderr, "Error: file is not a regular file.\n");
        return 1;
    }
    if (buff.st_size == sizeof(nvf_ld)) {
        FILE *fp = fopen(path, "rb");
        if (!fp) {
            fprintf(stderr, "Error: could not open file for reading.\n");
            return 1;
        }
        nvf_ld _info;
        fread(&_info, sizeof(nvf_ld), 1, fp);
        fclose(fp);
        printf("Header: %c, %c, %c, %c\n"
               "sizeof(T): %" PRIu32 " bytes\n"
               "sizeof(long double): %" PRIu32 " bytes\n""Normaliser ID: %" PRIu64 "\n"
               "Max. mass: %Lf\nMin. mass: %Lf\nMax. time: %Lf\nMax. pos:  %Lf\nMax. vel:  %Lf\n",
               _info.header[0], _info.header[1], _info.header[2],
               _info.header[3], _info._sizeof_T, _info._sizeof_ld, _info._norm_id, _info._maxm, _info._minm,
               _info._maxt, _info._maxp, _info._maxv);
        return 0;
    }
    if (buff.st_size == sizeof(nvf_d)) {
        FILE *fp = fopen(path, "rb");
        if (!fp) {
            fprintf(stderr, "Error: could not open file for reading.\n");
            return 1;
        }
        nvf_d _info;
        fread(&_info, sizeof(nvf_d), 1, fp);
        fclose(fp);
        printf("Header: %c, %c, %c, %c\n"
               "sizeof(T): %" PRIu32 " bytes\n"
               "sizeof(long double): %" PRIu32 " bytes\n""Normaliser ID: %" PRIu64 "\n"
               "Max. mass: %lf\nMin. mass: %lf\nMax. time: %Lf\nMax. pos:  %lf\nMax. vel:  %lf\n",
               _info.header[0], _info.header[1], _info.header[2],
               _info.header[3], _info._sizeof_T, _info._sizeof_ld, _info._norm_id, _info._maxm, _info._minm,
               _info._maxt, _info._maxp, _info._maxv);
        return 0;
    }
    if (buff.st_size == sizeof(nvf_f)) {
        FILE *fp = fopen(path, "rb");
        if (!fp) {
            fprintf(stderr, "Error: could not open file for reading.\n");
            return 1;
        }
        nvf_f _info;
        fread(&_info, sizeof(nvf_f), 1, fp);
        fclose(fp);
        printf("Header: %c, %c, %c, %c\n"
               "sizeof(T): %" PRIu32 " bytes\n"
               "sizeof(long double): %" PRIu32 " bytes\n""Normaliser ID: %" PRIu64 "\n"
               "Max. mass: %f\nMin. mass: %f\nMax. time: %Lf\nMax. pos:  %f\nMax. vel:  %f\n",
               _info.header[0], _info.header[1], _info.header[2],
               _info.header[3], _info._sizeof_T, _info._sizeof_ld, _info._norm_id, _info._maxm, _info._minm,
               _info._maxt, _info._maxp, _info._maxv);
        return 0;
    }
    fprintf(stderr, "Error: invalid file size.\n");
    return 1;
}

int main(int argc, char **argv) {
    if (argc > 2) {
        fprintf(stderr, "Error: invalid number of command-line arguments.\nUsage: nnv [opt:filename]\n");
        return 1;
    }
    if (argc == 2)
        return process_normv(*(argv + 1));
    char *path = find_normv();
    if (!path) {
        fprintf(stderr, "Error: no .normv file provided and could not find any in current directory.\n");
        return 1;
    }
    int _ret = process_normv(path);
    free(path);
    return _ret;
}

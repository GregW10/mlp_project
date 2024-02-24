#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>

typedef enum boolean {
	false, true
} bool;

uint64_t to_uint64(const char *str) {
	if (!str)
		exit(1);
	uint64_t tot = 0;
	goto start;
	while (*str) {
		tot *= 10;
		start:
		if (*str < 48 || *str > 57) {
			fprintf(stderr, "\033[1m\033[31mError: \033[32mcharacter '%c' is non-numeric.\n", *str);
			exit(1);
		}
		tot += *str++ - 48;
	}
	return tot;
}

bool equal(const char *s1, const char *s2) {
	if (!s1 || !s2)
		return false;
	while (*s1 || *s2)
		if (*s1++ != *s2++)
			return false;
	return true;
}

int main(int argc, char **argv) {
	if (argc != 3) {
		fprintf(stderr, "\033[1m\033[31mError: \033[38;5;11minvalid number of command-line arguments.\n");
		return 1;
	}
	const char *flt = *(argv + 1);
	uint64_t Nt = to_uint64(*(argv + 2));
	size_t flt_s;
	if (equal(flt, "long double"))
		flt_s = sizeof(long double);
	else if (equal(flt, "double"))
		flt_s = sizeof(double);
	else if (equal(flt, "float"))
		flt_s = sizeof(float);
	else {
		fprintf(stderr, "\033[1m\033[38;5;9mError: \033[32munrecognised floating-point type.\n");
		return 1;
	}
	printf("\033[38;5;93m.3bod file size \033[38;5;226m=\033[1m\033[38;5;10m %" PRIu64 "\033[0m \033[38;5;14mbytes\n",
           12 + 4*sizeof(long double) + 6*flt_s*(1 + 3*Nt));
	return 0;
}

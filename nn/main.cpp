#include "../glib/ml/gregmtx.hpp"
#include "../data_management/datsup.hpp"
#include <chrono>

#define NUMBER 262144

struct vecc {
    long double x, y, z;
    unsigned long long s;
    char c;
    char a;
    short b[3];
};

int main(int argc, char **argv) {
    gml::matrix<long double> mat{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
    std::cout << "mat.tsr bytes written: " << mat.to_tsr("mat.tsr") << std::endl;
    std::cout << mat << std::endl;
    gml::matrix<long double> m2 = {{5, 8, 3, 4, 1},
                                   {1, 9, 0, 4, 2},
                                   {7, 1, 4, 9, 9},
                                   {1, 5, 3, 7, 0},
                                   {4, 9, 7, 6, 2}};
    std::cout << m2 << std::endl;
    gml::tensor<long double> t{"mat.tsr"};
    std::cout << "t:\n" << t << std::endl;
    std::cout << "Success?\n";
    std::cout << std::boolalpha << gml::gen::endswith("I dunno bro!", "!") << std::endl;
    std::cout << "mat.mtsr bytes written: " << m2.to_tsr("mat.mtsr") << std::endl;
    gml::matrix<long double> m3{"mat.mtsr"};
    std::cout << "m3:\n" << m3 << std::endl;
    std::cout << "m3[2][2]: " << m3(2, 2) << ", " << m3.at(2, 2) << std::endl;
    std::vector<long double> d1;
    d1.reserve(NUMBER);
    std::vector<long double> d2;
    d2.reserve(NUMBER);
    std::mt19937_64 rng{std::random_device{}()};
    std::uniform_real_distribution<long double> dist{-1, 1};
    uint64_t counter = NUMBER;
    while (counter --> 0) {
        d1.push_back(dist(rng));
        d2.push_back(dist(rng));
    }
    gml::matrix<long double> r{d1, 512, 512};
    long double *p1 = new long double[NUMBER];
    long double *p2 = new long double[NUMBER];
    long double *org1 = p1;
    long double *org2 = p2;
    int i = 0, j = 0;
    std::chrono::time_point<std::chrono::system_clock> start1 = std::chrono::system_clock::now();
    for (; i < 512; ++i) {
        for (; j < 512; ++j) {
            *p1++ = r.at(i, j);
        }
    }
    std::chrono::time_point<std::chrono::system_clock> end1 = std::chrono::system_clock::now();
    std::chrono::time_point<std::chrono::system_clock> start2 = std::chrono::system_clock::now();
    i = 0; j = 0;
    for (; i < 512; ++i) {
        for (; j < 512; ++j) {
            *p2++ = r(i, j);
        }
    }
    std::chrono::time_point<std::chrono::system_clock> end2 = std::chrono::system_clock::now();
    std::chrono::nanoseconds duration1 = std::chrono::duration_cast<std::chrono::nanoseconds>(end1 - start1);
    std::chrono::nanoseconds duration2 = std::chrono::duration_cast<std::chrono::nanoseconds>(end2 - start2);
    std::cout << "Duration of .at(): " << duration1.count() << std::endl;
    std::cout << "Duration of .operator(): " << duration2.count() << std::endl;
    r.to_tsr("big.tsr");
    r.to_mtsr("big.mtsr");
    delete [] org1;
    delete [] org2;
    gml::matrix<long double> good{"big.tsr"};
    gml::matrix<long double> bad{"big.mtsr"};
    std::cout << std::boolalpha << (r == good && r == bad && good == bad) << std::endl;
    std::cout << "Long double: " << typeid(long double).name() << std::endl;
    std::cout << "Common between d and ld: " << typeid(std::common_type<double, long double>::type).name() << std::endl;
    std::cout << "Double: " << typeid(double).name() << std::endl;
    std::cout << "Common between d and f: " << typeid(std::common_type<double, float>::type).name() << std::endl;
    return 0;
}

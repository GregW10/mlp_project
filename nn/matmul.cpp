#include "../glib/ml/gregvct.hpp"
#include <random>
#include <chrono>

#define NUMBER 262144

#define TYPE long double

int main() {
    std::vector<TYPE> d1;
    d1.reserve(NUMBER);
    std::vector<TYPE> d2;
    d2.reserve(NUMBER);
    std::mt19937_64 rng{std::random_device{}()};
    std::uniform_real_distribution<TYPE> dist{-1, 1};
    // std::uniform_int_distribution<TYPE> dist{0, 10};
    uint64_t counter = NUMBER;
    while (counter --> 0) {
        d1.push_back(dist(rng));
        d2.push_back(dist(rng));
    }
    gml::matrix<TYPE> m1{d1.data(), 512, 512};
    gml::matrix<TYPE> m2{d2.data(), 512, 512};
    // std::cout << "m1:\n" << m1 << std::endl;
    // std::cout << "m2:\n" << m2 << std::endl;
    std::chrono::time_point<std::chrono::system_clock> start = std::chrono::system_clock::now();
    gml::matrix<TYPE> m3 = m1*m2;
    std::chrono::time_point<std::chrono::system_clock> end = std::chrono::system_clock::now();
    std::chrono::nanoseconds duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
    // std::cout << "m3:\n" << m3 << std::endl;
    std::cout.precision(10);
    std::cout << "Done! " << duration.count()/BILLION << " s taken to perform multiplication." << std::endl;
    m1.to_mtsr("m1.mtsr");
    m2.to_mtsr("m2.mtsr");
    m3.to_mtsr("m3.mtsr");
    gml::matrix<long double> test;
    std::cout << test << std::endl;
    std::cout << test.shape() << std::endl;
    std::cout << test.shape().volume() << std::endl;
    gml::matrix<long double> t2{3, 6};
    std::cout << t2 << std::endl;
    return 0;
}

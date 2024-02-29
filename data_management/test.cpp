#include <iostream>
#include <fstream>
#include "../glib/misc/gregparse.hpp"
#include "../glib/nbod/gregstr.hpp"

#pragma pack(push, 1)
template <typename T>
struct vec {
    T x{}, y{}, z{};
};
#pragma pack(pop)

class A {
    int x{}, y{}, z{};
public:
    A() = default;
    A(int _x, int _y, int _z) : x{_x}, y{_y}, z{_z} {}
    virtual A &negate() {
        x = -x;
        y = -y;
        z = -z;
        return *this;
    }
};

template <bool b>
class B : public A {
public:
    using A::A;
    A &negate() override {
        A &ref = *this;
        B<!b> &rb = dynamic_cast<B<!b>&>(ref);
        return rb;
    }
};

int main(int argc, char **argv) {
    gtd::parser parser{argc, argv};
    const char *arg = parser.get_arg(std::regex{R"(^r-p?v?a?c?e?$)"});
    std::cout << (arg ? arg : "nullptr") << std::endl;
    auto var = [](auto ...args){
        return (args + ...);
    };
    std::cout << "Sum: " << var(7, 6, 5, 4, 3, 2, 1) << std::endl;
    std::cout << std::boolalpha << gtd::is_integral(argc > 1 ? *(argv + 1) : "2382929") << std::endl;
    std::cout << std::boolalpha << gtd::contains(argc > 1 ? *(argv + 1) : "2382929", '1') << std::endl;
    vec<long double> v;
    std::cout << "x: " << v.x << ", y: " << v.y << ", z: " << v.z << std::endl;
    std::cout << sizeof(vec<char>) << std::endl;
    // if (setenv("TERM", "dumb", 1) == -1) // doesn't work :( or rather setting TERM to "dumb" doesn't stop coloured outp.
    //     std::cerr << "setenv error\n";
    // std::cout << "\033[1m\033[32mTesting...\n";
    B<true> b{3, 2, 1};
    b.negate();
    return 0;
}

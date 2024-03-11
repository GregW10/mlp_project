#include <iostream>

class A {
    char *a{};
public:
    A() : a{new char(50)} {}
    ~A() {delete a; std::cout << "A's dtor called." << std::endl;}
    void speak() {std::cout << "speaking" << std::endl;}
};

[[noreturn]] void bad(const A &a) {
    std::cout << "bad called. This function doesn't return." << std::endl;
    exit(0);
}

void good(const A &a) {
    std::cout << "good called. This function throws instead of using \"exit\"" << std::endl;
    A b;
    b.speak();
    throw 42;
}

int main(int argc, char **argv) {
    A a;
    if (argc == 2)
        bad(a);
    if (argc > 2)
        good(a);
    return 0;
}

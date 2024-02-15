#include <iostream>
#include <fstream>

namespace gml {
    template <typename T>
    concept Numeric = requires (T value) {
        T{};
        // lots of stuff...
    };
    template <Numeric T, bool alloc = true>
    class tensor {
        // lots more stuff...
        class shape {
            int n;
        public:
            shape(int num) requires (alloc) : n{num} {}
            shape(int num) requires (!alloc) : n{-num} {}
            void print() const noexcept {
                if constexpr (alloc)
                    std::cout << "alloc true: " << n << '\n';
                else
                    std::cout << "alloc false: " << n << '\n';
            }
            int get() const noexcept {
                return n;
            }
        };
        shape s{4};
    public:
        tensor() {
            std::cout << "tensor ctor.\n";
        }
        void shit() noexcept(alloc) {
            if constexpr (std::same_as<T, long double>)
                if (s.get() == 4)
                    throw std::invalid_argument{"Piss off.\n"};
            s.print();
        }
        template <Numeric U, Numeric V>
        friend bool operator==(const tensor<U>&, const tensor<V>&);
    };
    template <Numeric U, Numeric V>
    bool operator==(const tensor<U> &t1, const tensor<V> &t2) {
        return true;
    }
}

int main(int argc, char **argv) {
    if (argc != 2)
        return 1; /*
    gml::tensor<long double> t;
    gml::tensor<long double, false> t3;
    // std::cout << std::boolalpha << "t == t3: " << (t == t3) << std::endl;
    t.shit();
    t3.shit(); */
    std::cout << sizeof(std::ifstream) << std::endl;
    std::ifstream in{*(argv + 1), std::ios_base::in | std::ios_base::binary};
    if (!in.good())
        return 1;
    uint64_t val; /*
    in.seekg(4 + sizeof(long double));
    in.read((char *) &val, sizeof(uint64_t)); */
    in.read((char *) &val, 8);
    in.get();
    std::cout << std::boolalpha << "in.good(): " << in.good() << ", in.bad(): " << in.bad() << ", in.fail(): "
              << in.fail() << ", in.eof(): " << in.eof() << std::endl;
    in.close();
    std::cout << "Number of iterations: " << val << std::endl;
    return 0;
}

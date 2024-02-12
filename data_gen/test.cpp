#include <iostream>

namespace gml {
    template <typename T>
    concept Numeric = requires (T value) {
        T{};
        // lots of stuff...
    };
    template <Numeric T>
    class tensor {
        // lots more stuff...
    public:
        tensor() = default; // constructs an empty tensor
        template <Numeric U, Numeric V>
        friend bool operator==(const tensor<U>&, const tensor<V>&);
    };
    template <Numeric U, Numeric V>
    bool operator==(const tensor<U> &t1, const tensor<V> &t2) {
        return true;
    }
}

int main() {
    gml::tensor<long double> t;
    gml::tensor<long double> t3;
    std::cout << std::boolalpha << "t == t3: " << (t == t3) << std::endl;
    return 0;
}

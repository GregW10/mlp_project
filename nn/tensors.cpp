#include "../glib/ml/gregmtx.hpp"

int main() {
    gml::tensor<long double> tens;
    std::cout << "Tensor: " << tens << std::endl;
    std::cout << "Shape: " << tens.shape() << std::endl;
    std::cout << "Volume: " << tens.shape().volume() << std::endl;
    gml::tensor<long double> t2(4.3l);
    std::cout << "Tensor: " << t2 << std::endl;
    std::cout << "Shape: " << t2.shape() << std::endl;
    std::cout << "Volume: " << t2.shape().volume() << std::endl;
    tens.to_tsr("empty.tsr");
    t2.to_tsr("scalar.tsr");
    gml::tensor<long double> empty{"empty.tsr"};
    gml::tensor<long double> scalar{"scalar.tsr"};
    std::cout << "Empty: " << empty << std::endl;
    std::cout << "Scalar: " << scalar << std::endl;
    gml::tensor<long double> squashed{gml::tensor_shape{1, 2, 3, 4, 0, 0, 9, 12, 14, 3, 2, 0}};
    std::cout << squashed << std::endl;
    squashed.to_tsr("squished.tsr");
    gml::tensor<long double> squished{"squished.tsr"};
    std::cout << squished << std::endl;
    std::cout << squished.shape() << std::endl;
    std::cout << squished.shape().volume() << std::endl;
    return 0;
}

/* Test program to ensure my .3bod reader is working correctly. */

#include "datsup.hpp"

using etype = typename gtd::f3bodr<long double>::entry_type;

void add_entry(std::vector<etype> &entries, const gtd::sys &sys){
    static etype entry;
    long double *ptr = (long double *) &entry;
    char counter = 1;
    for (const gtd::bod_0f &b : sys) {
        std::cout << "x_" << +counter << ": " << (*ptr++ = b.pos()[0]) << " m\n"
                  << "y_" << +counter << ": " << (*ptr++ = b.pos()[1]) << " m\n"
                  << "z_" << +(counter++) << ": " << (*ptr++ = b.pos()[2]) << " m\n";
    }
    counter = 1;
    for (const gtd::bod_0f &b : sys) {
        std::cout << "vx_" << +counter << ": " << (*ptr++ = b.vel()[0]) << " m/s\n"
                  << "vy_" << +counter << ": " << (*ptr++ = b.vel()[1]) << " m/s\n"
                  << "vz_" << +(counter++) << ": " << (*ptr++ = b.vel()[2]) << " m/s\n";
    }
    entries.push_back(entry);
}

void print_entry(const etype &e) {
    long double *ptr = (long double *) &e;
    char counter = 0;
    while (counter++ < 3) {
        std::cerr << "x_" << +counter << ": " << *ptr++ << " m\n"
                  << "y_" << +counter << ": " << *ptr++ << " m\n"
                  << "z_" << +counter << ": " << *ptr++ << " m\n";
    }
    counter = 0;
    while (counter++ < 3) {
        std::cerr << "vx_" << +counter << ": " << *ptr++ << " m/s\n"
                  << "vy_" << +counter << ": " << *ptr++ << " m/s\n"
                  << "vz_" << +counter << ": " << *ptr++ << " m/s\n";
    }
}

int main(int argc, char **argv) {
    gtd::parser parser{argc, argv};
    // uint64_t evols = parser.get_arg("--evolutions", 1'000ull);
    uint64_t epochs = parser.get_arg("--epochs", 1'000ull);
    uint64_t iters = parser.get_arg("--iterations", 1'000ull);
    gtd::sys sys;
    sys.emplace_body(true, 1'000'000.0l, 1.0l, gtd::vec3{0, 5, 5}, gtd::vec3{});
    sys.emplace_body(true, 2'000'000.0l, 1.0l, gtd::vec3{0, -5, 5}, gtd::vec3{});
    sys.emplace_body(true, 4'000'000.0l, 1.0l, gtd::vec3{-5, -5, -5}, gtd::vec3{});
    sys.iters(iters);
    std::vector<etype> entries;
    entries.reserve(epochs + 1);
    std::cout << "----------\nEpoch 0\n----------\n";
    add_entry(entries, sys);
    putchar('\n');
    gtd::f3bodw<long double> writer{"log.3bod", &sys, epochs  + 1};
    writer.add_entry();
    uint64_t counter = 1;
    while (counter <= epochs) {
        sys.evolve();
        std::cout << "----------\nEpoch " << counter++ << "\n----------\n";
        add_entry(entries, sys);
        putchar('\n');
        writer.add_entry();
    }
    gtd::f3bodr<long double, false, false> reader{"log.3bod"};
    counter = 0;
    for (const auto &entry : reader) {
        std::cerr << "----------\nEpoch " << counter++ << "\n----------\n";
        print_entry(entry);
        fputc('\n', stderr);
    }
    std::cout << "Epoch 1000:\n";
    print_entry(reader.entry_at(1000));
    reader.close();
    return 0;
}

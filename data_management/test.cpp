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
    pid_t pid = fork();
    if (pid == 0) {
        // execlp("sleep", "sleep", "30", (char *) nullptr);
        const char *const args[] = {"sleep", "30", nullptr};
        execvp("sleep", const_cast<char *const *>(args));
    }
    std::cout << "Waiting for child process with PID " << pid << std::endl;
    wait(nullptr);
    std::cout << "Child process exited." << std::endl;
    return 0;
}

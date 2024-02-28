#include <iostream>
#include <unistd.h>
#include <sys/stat.h>
#include <string>
#include <set>
#include <fstream>
#include <regex>
#include <sstream>
#include <cstring>

int main(int argc, char **argv) {
    if (argc != 2) {
        std::cerr << "Error: filename containing commands expected.\n";
        return 1;
    }
    const char *path = *(argv + 1);
    struct stat buff{};
    if (stat(path, &buff) == -1) {
        std::cerr << "Error: could not obtain file information.\n";
        return 1;
    }
    if (!S_ISREG(buff.st_mode)) {
        std::cerr << "Error: file is not a regular file.\n";
        return 1;
    }
    std::vector<std::string> commands;
    std::string command;
    pid_t pid;
    auto print = [](const std::string &s){
        for (const char &c : s)
            std::cout << +c << std::endl;
    };
    std::ifstream in{path, std::ios_base::in};
    while (std::getline(in, command))
        if (!std::regex_match(command, std::regex{R"(^(#.*|\s*)$)"}))
            commands.push_back(command);
    in.close();
    for (std::string &cmd : commands) {
        if ((pid = fork()) == -1) {
            std::cerr << "Error: fork() error.\n";
            return 1;
        }
        if (pid == 0) { // child process
            std::vector<char*> args;
            char *arg_cstr;
            std::istringstream iss{cmd};
            std::string arg;
            std::string long_arg;
            while (iss >> arg) {
                if (arg.starts_with('\'')) {
                    arg.erase(arg.begin());
                    long_arg = arg;
                    while (iss >> arg)
                        long_arg += " " + arg;
                    if (!long_arg.ends_with('\'')) {
                        std::cerr << "Error: argument contains unterminated single quote.\n";
                        return 1;
                    }
                    long_arg.erase(long_arg.end() - 1);
                    arg_cstr = new char[long_arg.size() + 1];
                    strcpy(arg_cstr, long_arg.c_str());
                    long_arg.clear();
                    args.push_back(arg_cstr);
                } else if (arg.starts_with("\"")) {
                    arg.erase(arg.begin());
                    long_arg = arg;
                    while (iss >> arg)
                        long_arg += " " + arg;
                    if (!long_arg.ends_with('"')) {
                        std::cerr << "Error: argument contains unterminated double quote.\n";
                        return 1;
                    }
                    long_arg.erase(long_arg.end() - 1);
                    arg_cstr = new char[long_arg.size() + 1];
                    strcpy(arg_cstr, long_arg.c_str());
                    long_arg.clear();
                    args.push_back(arg_cstr);
                } else {
                    arg_cstr = new char[arg.size() + 1];
                    strcpy(arg_cstr, arg.c_str());
                    args.push_back(arg_cstr);
                }
            }
            args.push_back(nullptr);
            if (execvp(args[0], args.data()) == -1) {
                std::cerr << "Error: execvp() error.\n";
                return 1;
            }
        }
        if (wait(nullptr) != pid) {
            std::cerr << "Error: wait() error.\n";
            return 1;
        }
    }
    return 0;
}

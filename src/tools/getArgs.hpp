#ifndef GETARGS_HPP
#define GETARGS_HPP

#include <string>

#ifdef _WIN32
#include "unistdWindows.hpp"
#else
#include <unistd.h>
#endif

#include <cstring>
#include <tuple>
#include "type.hpp"

class argsController {
    int argc;
    char ** argv;
public:
    argsController(int a, char ** b) {
        argc = a;
        argv = b;
    }

    std::string get(const std::string & s) {
        for(int i = 0; i < argc; i++) {
            if(s == argv[i]) {
                return i + 1 == argc ? "" : argv[i+1];
            }
        }
        return "";
    }

    bool exist(const std::string & s) {
        for(int i = 0; i < argc; i++) {
            if(s == argv[i]) return true;
        }
        return false;
    }
};

#endif
#include "./configparser.h"

#include <string>
#include <iostream>
#include <fstream>
#include "rapidjson/filereadstream.h"

rapidjson::Document read_config(const std::string& fname, int verbose) {
    std::ifstream f(fname);
    if (!f.good()) {
        std::cout << "config file \""<<fname<<"\" failed to open with error ("<<errno<<"): " << std::strerror(errno);
    }
    std::string buf;
    // reserve just enough space for file stringbuffer
    f.seekg(0, std::ios::end);
    buf.reserve(f.tellg());
    f.seekg(0, std::ios::beg);
    buf.assign((std::istreambuf_iterator<char>(f)),
                std::istreambuf_iterator<char>());
    // convert to rapidjson stream
    rapidjson::StringStream stream(buf.c_str());
    rapidjson::Document doc;
    doc.ParseStream(stream);

    if (verbose) {
        std::cout << "CONFIG PARSER:" << std::endl << "--------------" << std::endl;
        std::cout << "Found " << doc.MemberCount() << " members in \""<<fname<<"\"" << std::endl;
    }

    return doc;
}

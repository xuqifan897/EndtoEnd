#include "logging.h"

#include <iostream>
#include <sstream>

void Logger::print_head(const std::string& title) {
    m_n = title.length() + 1;
    std::cout << std::string(m_n, '-') << std::endl;
    std::cout << title << ":" << std::endl;
    std::cout << std::string(m_n, '-') << std::endl;
}
void Logger::print_tail() {
    std::cout << std::string(m_n, '=') << std::endl;
}

std::string set_color(COLOR foreground, COLOR background) {
    return set_color(static_cast<unsigned int>(foreground), static_cast<unsigned int>(background));
}
std::string set_color(unsigned int foreground, unsigned int background) {
    std::ostringstream s;
    s << "\033[";

    if (!foreground && ! background) { s << "0"; } // reset colors if no params

    if (foreground) {
        s << (29 + foreground);

        if (background) { s << ";"; }
    }

    if (background) { s << (39 + background); }
    s << "m";

    return s.str();
}

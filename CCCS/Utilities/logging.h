#ifndef __LOGGING_H__
#define __LOGGING_H__

#include <string>

// class providing standardized logging routines for stdout printing used in cli programs
class Logger {
    private:
        unsigned int m_n = 0; // length of divider for previously called print_head()
    public:
        void print_head(const std::string& title);
        void print_tail();
};

enum COLOR {
    NONE = 0,
    BLACK, RED, GREEN,
    YELLOW, BLUE, MAGENTA,
    CYAN, WHITE
};
std::string set_color(unsigned int foreground, unsigned int background);
std::string set_color(COLOR foreground=COLOR::NONE, COLOR background=COLOR::NONE);

#endif // __LOGGING_H__


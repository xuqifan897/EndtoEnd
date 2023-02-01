#ifndef __IO_HELPERS_H__
#define __IO_HELPERS_H__

#include <string>
#include <array>
#include <vector>

// TODO: Add safe "open()" method that standardizes error handling when files cannot be found/opened
namespace dcio {
    // Returns filesize in bytes
    int filesize(const std::string& filename);

    std::string get_dirname(std::string path);
    std::string get_basename(std::string path);
    bool file_exists(const std::string& path);
    int copy_file(const std::string& oldpath, const std::string& newpath);
    int move_file(const std::string& oldpath, const std::string& newpath);
    std::array<std::string, 2> splitext(std::string fname);

    // get path type
    bool is_relative(const std::string& path);
    bool is_absolute(const std::string& path);
    bool is_unqualified(const std::string& path);

    // create directory with full path "dir" and permissions "perms", ignore EEXIST error if exist_ok==True
    int create_directory(const std::string& dir, const bool exist_ok=true, const mode_t perms=0775);
    bool dir_exists(const std::string& dir);
    std::vector<std::string> get_file_list(const std::string& dir, bool verbose=false);
    int remove_file(const std::string& fname);
    int remove_files(const std::vector<std::string>& files);
    int remove_directory(const std::string& dname, bool keepdir=false, bool recursive=false);
    int remove_directories(const std::vector<std::string>& dirs, bool keepdir=false, bool recursive=false);

    // test for comment line in file, given the comment_char
    bool is_comment_string(const std::string& str, const char comment_char);

    // split string into vector of tokens using char delimiters specified in string: "delim"
    void tokenize_string(const std::string& str, std::vector<std::string>& tokens, const std::string& delims=" ");

    // convert string to lowercase in place
    std::string tolower(std::string retval);
    char* tolower(char*);

    // test if string is a valid floating point number
    bool is_number(const std::string& test);

    std::string get_username();
};

#endif // __IO_HELPERS_H__

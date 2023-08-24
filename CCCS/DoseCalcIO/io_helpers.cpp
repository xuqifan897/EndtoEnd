#include "io_helpers.h"

#include <cstdlib>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/sendfile.h>
#include <pwd.h>
#include <fcntl.h>
#include <dirent.h>
#include <cerrno>
#include <string>
#include <cstring>
#include <cctype>
#include <cstdio>
#include <fstream>
#include <algorithm>
#include <libgen.h>
#include <cerrno>
#include <iostream>

namespace dcio {
    ///////////////////////////////////////////////////
    //////////////// FILESYSTEM ///////////////////////
    ///////////////////////////////////////////////////
    int filesize(const std::string& filename) {
        // struct stat filestatus;
        // stat(filename.c_str(), &filestatus);
        // return filestatus.st_size;

        std::ifstream infile(filename, std::ios::ate | std::ios::binary);
        return static_cast<int>(infile.tellg());
    }

    std::string get_dirname(std::string path) {
        return std::string(::dirname((char*)path.c_str()));
    }
    std::string get_basename(std::string path) {
        return std::string(::basename((char*)path.c_str()));
    }
    std::array<std::string, 2> splitext(std::string fname) {
        // get basename
        std::string bname = get_basename(fname);

        // find last '.'
        size_t pos;
        if ( std::string::npos != (pos = bname.find_last_of('.'))) {
            return { {bname.substr(0, pos), bname.substr(pos, std::string::npos)} };
        } else {
            // no dot found
            return { {bname, ""} };
        }


    }

    bool file_exists(const std::string& path) {
        struct stat buffer;
        return (stat(path.c_str(), &buffer) == 0);
    }
    int copy_file(const std::string& oldpath, const std::string& newpath) {
        int status;

        int read_fd;
        int write_fd;
        struct stat stat_buf;
        off_t offset = 0;

        read_fd = open(oldpath.c_str(), O_RDONLY);
        fstat (read_fd, &stat_buf);
        write_fd = open(newpath.c_str(), O_WRONLY | O_CREAT, stat_buf.st_mode);
        status = sendfile(write_fd, read_fd, &offset, stat_buf.st_size);
        close (read_fd); close (write_fd);

        if (status < 0) {
            printf("copy_file() failed with error (%d): %s\n", errno, std::strerror(errno));
        }
        return (status >= 0);
    }
    int move_file(const std::string& oldpath, const std::string& newpath) {
        int status;
        if (0 != (status = rename(oldpath.c_str(), newpath.c_str())) ) {
            // TODO: use copy_file && delete(oldpath) if error is from different partitions
            printf("move_file() failed with error (%d): %s\n", errno, std::strerror(errno));
        }
        return (status == 0);
    }


    bool is_relative(const std::string& path) {
        // check for "." or ".." in front
        if (path.empty()) {return false;}

        // check for relative path without dot qualification
        size_t idx;
        if (std::string::npos != (idx = path.find_first_of("/"))) {
            // check that not escaped
            if (path[idx-1] != '\\') { return true; }
        }

        // check for dot qualification
        if (path[0] == '.') {
            if (path[1] == '.' || path[1] == '/') { return true; }
        }
        return false;
    }
    bool is_absolute(const std::string& path) {
        // check for "/" in front
        if (path.empty()) {return false;}

        // check for root base
        if (path.front() == '/') { return true; }

        // check for home expansion
        if (path[0] == '~' && (path[1] == '/')) { return true; }

        return false;
    }
    bool is_unqualified(const std::string& path) { return (!is_relative(path) && !is_absolute(path)); }

    int create_directory(const std::string& dir, const bool exist_ok, const mode_t perms) {
#if defined(_WIN32)
        int _fail = _mkdir(dir.c_str());
#else
        int _fail = mkdir(dir.c_str(), perms);
#endif
        if (_fail && (errno!=EEXIST || !exist_ok)) {
            printf("create_directory() failed with error (%d): %s\n", errno, std::strerror(errno));
        }
        return errno;
    }
    bool dir_exists(const std::string& dir) {
        if (DIR* FD = opendir(dir.c_str())) {
            closedir(FD);
            return true;
        } else if (errno == ENOENT) {
            return false;
        } else {
            printf("dir_exists() failed with error (%d): %s\n", errno, std::strerror(errno));
            return false;
        }
    }
    std::vector<std::string> get_file_list(const std::string& dir, bool verbose) {
        DIR* FD;
        struct dirent* in_file;
        std::vector<std::string> files;

        /* Scanning the in directory */
        if (NULL == (FD = opendir(dir.c_str()))) {
            if (verbose) { fprintf(stderr, "Error : Failed to open input directory \"%s\" - %s\n", dir.c_str(), strerror(errno)); }
            return std::vector<std::string>();
        }
        while ((in_file = readdir(FD))) {
            /* On linux/Unix we don't want current and parent directories
             * On windows machine too
             */
            if (!strcmp (in_file->d_name, ".")) { continue; }
            if (!strcmp (in_file->d_name, "..")) { continue; }
            files.push_back(in_file->d_name);
        }
        return files;
    }
    bool is_regular_file(const std::string& path) {
        struct stat path_stat;
        stat(path.c_str(), &path_stat);
        return S_ISREG(path_stat.st_mode);
    }
    bool is_symlink(const std::string& path) {
        struct stat path_stat;
        stat(path.c_str(), &path_stat);
        return S_ISLNK(path_stat.st_mode);
    }
    bool is_directory(const std::string& path) {
        struct stat path_stat;
        stat(path.c_str(), &path_stat);
        return S_ISDIR(path_stat.st_mode);
    }
    int remove_file(const std::string& fname) {
        unlink(fname.c_str());
        return true;
    }
    int remove_files(const std::vector<std::string>& files) {
        for (const auto& fname : files) {
            if (!remove_file(fname)) { return false; }
        }
        return true;
    }
    int remove_directory(const std::string& dname, bool keepdir, bool recursive) {
        for (const auto& fname : get_file_list(dname)) {
            // check if directory
            if (recursive && is_directory(dname+"/"+fname)) {
                remove_directory(dname+"/"+fname+"/", keepdir, recursive);
            } else if (is_regular_file(dname+"/"+fname) || is_symlink(dname+"/"+fname)) {
                unlink((dname + "/" + fname).c_str());
            }
        }
        if (!keepdir) { rmdir(dname.c_str()); }
        return true;
    }
    int remove_directories(const std::vector<std::string>& dirs, bool keepdir, bool recursive) {
        for (const auto& dname : dirs) {
            if (!remove_directory(dname, keepdir, recursive)) { return false; }
        }
        return true;
    }

    bool is_comment_string(const std::string& str, const char comment_char) {
        size_t firstchar = str.find_first_not_of(" \t");
        if (firstchar != std::string::npos && str[firstchar] == comment_char) {
            return true;
        }
        return false;
    }

    void tokenize_string(const std::string& str, std::vector<std::string>& tokens, const std::string& delims) {
        // skip delims at beginning
        size_t lastpos = str.find_first_not_of(delims, 0);
        // Find one after end of first token (next delim position)
        size_t nextpos = str.find_first_of(delims, lastpos);
        while (nextpos != std::string::npos || lastpos != std::string::npos) {
            // add token to vector
            tokens.push_back(str.substr(lastpos, nextpos-lastpos));
            // update positions
            lastpos = str.find_first_not_of(delims, nextpos);
            nextpos = str.find_first_of(delims, lastpos);
        }
    }
    std::string tolower(std::string retval) {
        std::transform(retval.begin(), retval.end(), retval.begin(), static_cast<int(*)(int)>(std::tolower));
        return retval;
    }
    char* tolower(char* retval) {
        for(uint i=0; retval[i]; ++i) {
            retval[i] = std::tolower(retval[i]);
        }
        return retval;
    }
    bool is_number(const std::string& test) {
        char *endptr;
        errno = 0;

        strtod(test.c_str(),&endptr);
        if (errno) {
            // printf("conversion failed");
            return false;
        }
        if (test.length()==0) {
            // printf("empty string\n");
            return false;
        } else {
            // printf("f=%f\n", f);
            if (*endptr == 0) {
                // printf("entire string valid\n");
                return true;
            } else {
                // printf("extra characters: %s\n", endptr);
                return false;
            }
        }
    }

    std::string get_username() {
        struct passwd* pw;
        uid_t uid;
        uid = geteuid();
        pw = getpwuid(uid);
        if (pw) {
            return pw->pw_name;
        } else {
            char buf[40];
            sprintf(buf, "Cannot obtain username for UID \"%u\"", (unsigned int)uid);
            throw std::runtime_error("Cannot obtain username for UID \"%u\"");
        }
    }
};


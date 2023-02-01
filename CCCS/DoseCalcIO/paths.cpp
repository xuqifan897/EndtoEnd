#include "paths.h"
#include "./io_helpers.h"
#include <boost/filesystem.hpp>

Paths* Paths::m_pInstance = NULL;
void Paths::Initialize() {
    // check env vars
    if (const char* env_d = std::getenv("DOSECALC_DATA")) {
        m_data_dir = std::string(env_d);
    } else {
        m_data_dir = "./data";
    }
    std::string user;
    try {
        user = dcio::get_username();
    } catch (std::runtime_error) {
        user = "unknown";
    }
    m_temp_dir = "/tmp/dosecalc/"+user;
    boost::filesystem::create_directories(m_temp_dir);
}

Paths* Paths::Instance() {
    if (!m_pInstance) {
        m_pInstance = new Paths;
        m_pInstance->Initialize();
    }
    return m_pInstance;
}


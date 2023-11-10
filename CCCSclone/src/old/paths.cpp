#include "paths.h"
#include "configure.h"
#include "argparse.h"

#include <string>
#include <boost/filesystem.hpp>

namespace fs = boost::filesystem;

old::Paths* old::Paths::m_pInstance = nullptr;
void old::Paths::Initialize()
{
    this->m_data_dir = fs::path(DATA_FOLDER);
    this->m_temp_dir = fs::path(dev::getarg<std::string>("dataFolder"));
    this->m_result_dir = fs::path(dev::getarg<std::string>("resultFolder"));
}

old::Paths* old::Paths::Instance()
{
    if (!m_pInstance)
    {
        m_pInstance = new Paths;
        m_pInstance->Initialize();
    }
    return m_pInstance;
}
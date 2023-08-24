#include "./logconfig.h"

#include "dcmtk/config/osconfig.h"
#include "dcmtk/oflog/oflog.h"

LogConfig* LogConfig::m_pInstance = nullptr;
LogConfig* LogConfig::Instance() {
    // set logging prefs
    OFLog::configure( OFLogger::LogLevel::ERROR_LOG_LEVEL );

    if (m_pInstance == nullptr) { m_pInstance = new LogConfig; }
    return m_pInstance;
}

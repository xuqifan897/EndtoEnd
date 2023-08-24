#ifndef __LOGCONFIG_H__
#define __LOGCONFIG_H__

class LogConfig {
    public:
        static LogConfig* Instance();

    private:
        LogConfig() {}
        LogConfig(LogConfig const&) = delete;
        LogConfig& operator=(LogConfig const&) = delete;
        static LogConfig* m_pInstance;
};

#endif // __LOGCONFIG_H__

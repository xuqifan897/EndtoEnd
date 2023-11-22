#ifndef __PATHS_H__
#define __PATHS_H__

#include "argparse.h"
#include <string>
#include <boost/filesystem.hpp>

namespace old
{
    class Paths
    {
    public:
        static Paths* Instance();

        // Getters
        std::string temp_dir() {return m_temp_dir.string();}
        std::string data_dir() {return m_data_dir.string();}
        std::string result_dir() {return m_result_dir.string();}
        std::string debug_dir() {return m_debug_dir.string();}
        std::string kernel_dir() {return (m_data_dir / boost::filesystem::path("dsa")).string();}
        std::string spectra_dir() {return (m_data_dir / boost::filesystem::path("spectra")).string();}

        std::string radius_file() {return (boost::filesystem::path(kernel_dir()) / 
            boost::filesystem::path("radii.dat")).string();}
        std::string polar_file() {return boost::filesystem::path(kernel_dir() / 
            boost::filesystem::path("polar.dat")).string();}
        std::string conv_theta_file() {return (m_temp_dir / 
            boost::filesystem::path("convolution_theta_angles.raw")).string();}
        std::string conv_phi_file() {return (m_temp_dir / 
            boost::filesystem::path("convolution_phi_angles.raw")).string();}
        std::string cum_kernel_file() {return (m_temp_dir / 
            boost::filesystem::path("cumulative_kernel.h5")).string();}

    protected:
        boost::filesystem::path m_temp_dir;
        boost::filesystem::path m_data_dir;
        boost::filesystem::path m_result_dir;
        boost::filesystem::path m_debug_dir;
    
    private:
        Paths() {};
        static Paths* m_pInstance;
        void Initialize();
    };
}

#endif
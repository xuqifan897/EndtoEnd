#ifndef __PATHS_H__
#define __PATHS_H__

#include <string>

/* Singleton class that gets path specfications from environment vars or uses defaults */
class Paths {
    public:
        static Paths* Instance();

        // Getters
        std::string temp_dir() { return m_temp_dir; }
        std::string data_dir() { return m_data_dir; }
        std::string kernel_dir() { return data_dir() + "/dsa"; }
        std::string spectra_dir() { return data_dir() + "/spectra"; }

        std::string radius_file() { return kernel_dir() + "/radii.dat"; }
        std::string polar_file() { return data_dir() + "/dsa/polar.dat"; }
        std::string conv_theta_file() { return "convolution_theta_angles"; }
        std::string conv_phi_file() { return "convolution_phi_angles"; }
        std::string cum_kernel_file() { return temp_dir() + "/cumulative_kernel.h5"; }

    protected:
        std::string m_temp_dir;
        std::string m_data_dir;

    private:
        Paths() {};
        static Paths* m_pInstance;
        void Initialize();
};

#endif

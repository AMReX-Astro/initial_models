#ifndef MODEL_UTIL_H
#define MODEL_UTIL_H

#include <string>

#include <AMReX_REAL.H>

std::string num_to_unitstring(const amrex::Real value) {

    std::string numstr{};
    std::stringstream ss;

    if (value > 1.e5) {
        // work in km
        ss << std::setprecision(2) << std::fixed << value / 1.e5;
        numstr = ss.str() + "km";
    } else {
        // work in cm
        ss << std::setprecision(2) << std::fixed << value;
        numstr = ss.str() + "cm";
    }

    return numstr;
}



#endif

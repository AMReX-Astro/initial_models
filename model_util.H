#ifndef MODEL_UTIL_H
#define MODEL_UTIL_H

std::string num_to_unitstring(const Real value) {

    std::string numstr{};
    std::stringstream ss;

    if (value > 1.e5) {
        // work in km
        ss << std::setprecision(2) << std::fixed << value / 1.e5;
        numstr = ss.str() + "km";
        std::cout << "here!!!!!!!!!!!!!!!!!!!" << numstr << std::endl;
    } else {
        // work in cm
        ss << std::setprecision(2) << std::fixed << value;
        numstr = ss.str() + "cm";
    }

    return numstr;
}



#endif
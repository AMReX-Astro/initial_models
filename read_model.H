#ifndef READ_MODEL_H
#define READ_MODEL_H

#include <AMReX_Print.H>

const int MAX_VARNAME_LENGTH = 80;

// define convenient indices for the scalars

namespace model {

    constexpr int nvar = 5 + NumSpec;

    constexpr int idens = 0;
    constexpr int itemp = 1;
    constexpr int ipres = 2;
    constexpr int ientr = 3;
    constexpr int iyef  = 4;  // this is ye -- we add "f" for file to not clash with the network
    constexpr int ispec = 5;
}

struct initial_model_t {
    int npts;
    amrex::Array2D<amrex::Real, 0, NPTS_MODEL-1, 0, model::nvar-1> state;
    amrex::Array1D<amrex::Real, 0, NPTS_MODEL-1> r;
};

namespace model_string
{
    inline std::string& ltrim(std::string& s)
    {
	auto it = std::find_if(s.begin(), s.end(),
                               [](int c) {
                                   return !std::isspace(c);
                               });
	s.erase(s.begin(), it);
	return s;
    }

    inline std::string& rtrim(std::string& s)
    {
	auto it = std::find_if(s.rbegin(), s.rend(),
                               [](int c) {
                                   return !std::isspace(c);
                               });
	s.erase(it.base(), s.end());
	return s;
    }
}

AMREX_INLINE void
read_file(const std::string filename, initial_model_t& initial_model) {

    std::ifstream mf;

    mf.open(filename, std::ios::in);
    if (!mf.is_open()) {
        amrex::Error("Error opening the initial model");
    }

    std::string line;

    // first the header line -- this tells us the number of points

    getline(mf, line);
    std::string npts_string = line.substr(line.find("=")+1, line.length());
    initial_model.npts = std::stoi(npts_string);

    if (initial_model.npts > NPTS_MODEL) {
        amrex::Error("Error: model has more than NPTS_MODEL points,");
    }

    // next line tells use the number of variables

    getline(mf, line);
    std::string num_vars_string = line.substr(line.find("=")+1, line.length());
    int nvars_model_file = std::stoi(num_vars_string);

    // now read in the names of the variables

    std::vector<std::string> varnames_stored;
    for (int n = 0; n < nvars_model_file; n++) {
        getline(mf, line);
        std::string var_string = line.substr(line.find("#")+1, line.length());
        varnames_stored.push_back(model_string::ltrim(model_string::rtrim(var_string)));
    }


    // start reading in the data

    amrex::Vector<Real> vars_stored;
    vars_stored.resize(nvars_model_file);

    for (int i = 0; i < initial_model.npts; i++) {
        mf >> initial_model.r(i);

        for (int j = 0; j < nvars_model_file; j++) {
            mf >> vars_stored[j];
        }

        for (int j = 0; j < model::nvar; j++) {
            initial_model.state(i,j) = 0.0_rt;
        }


        for (int j = 0; j < nvars_model_file; j++) {

            // keep track of whether the current variable from the model
            // file is one that we care about
            bool found_model = false;

            if (varnames_stored[j] == "density") {
                initial_model.state(i,model::idens) = vars_stored[j];
                found_model = true;

            } else if (varnames_stored[j] == "temperature") {
                initial_model.state(i,model::itemp) = vars_stored[j];
                found_model = true;

            } else if (varnames_stored[j] == "pressure") {
                initial_model.state(i,model::ipres) = vars_stored[j];
                found_model = true;

            } else if (varnames_stored[j] == "ye" || varnames_stored[j] == "Ye") {
                initial_model.state(i,model::iyef) = vars_stored[j];
                found_model = true;

            } else {
                for (int comp = 0; comp < NumSpec; comp++) {
                    if (varnames_stored[j] == spec_names_cxx[comp] || varnames_stored[j] == short_spec_names_cxx[comp]) {
                        initial_model.state(i,model::ispec+comp) = vars_stored[j];
                        found_model = true;
                        break;
                    }
                }
            }

            // yell if we didn't find the current variable

            if (!found_model && i == 0) {
                amrex::Print() << Font::Bold << FGColor::Yellow << "[WARNING] variable not found: " << varnames_stored[j] << ResetDisplay << std::endl;
            }

        } // end loop over nvars_model_file

    }
}

#endif

#ifndef READ_MODEL_H
#define READ_MODEL_H

#include <AMReX_Print.H>

#include <model_util.H>
#include <coord_info.H>

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


AMREX_INLINE Real cfmt(const Real x) {

    Real safe_x;

    safe_x = x;
    if (std::abs(x) < 1.0e-98_rt && x != 0.0) {
        safe_x = 1.0e-98_rt;
    }

    return safe_x;
}

AMREX_INLINE void
write_model(std::string model_name,
            const Array1D<Real, 0, NPTS_MODEL-1>& xzn_hse,
            const Array2D<Real, 0, NPTS_MODEL-1, 0, model::nvar-1>& model_hse,
            const bool write_ye=false) {

    // Write data stored in `model_state` array to file

    int npts = problem_rp::nx;
    if (problem_rp::use_irreg_grid) {
        npts = get_irreg_nr();
    }

    std::string outfile{};

    // if we are basing the model off an model file from a stellar
    // evolution code, then we use that model file name as the root of
    // the output

    if (! problem_rp::model_file.empty()) {
        int ipos = problem_rp::model_file.find(".dat");
        if (ipos < 0) {
            ipos = problem_rp::model_file.find(".txt");
        }
        if (ipos < 0) {
            ipos = problem_rp::model_file.find(".raw");
        }

        outfile += problem_rp::model_file.substr(0, ipos) + ".";
    }

    outfile += model_name;

    std::string dx_str{};

    if (problem_rp::use_irreg_grid) {
        dx_str = "irreg";
    } else {
        Real dx = (xzn_hse(1) - xzn_hse(0));
        dx_str = "dx" + num_to_unitstring(dx);
    }

    outfile += "." + dx_str;

    std::cout << "writing " << model_name << " model to " << outfile << std::endl;

    std::ofstream of;
    of.open(outfile);

    int num_out = write_ye ? 4+NumSpec : 3+NumSpec;

    of << "# npts = " << npts << std::endl;
    of << "# num of variables = " << num_out << std::endl;
    of << "# density" << std::endl;
    of << "# temperature" << std::endl;
    of << "# pressure" << std::endl;
    if (write_ye) {
        of << "# Ye" << std::endl;
    }

    for (int n = 0; n < NumSpec; ++n) {
        of << "# " << spec_names_cxx[n] << std::endl;
    }

    for (int i = 0; i < npts; ++i) {
        of << std::setprecision(12) << std::setw(20) << xzn_hse(i) << " ";
        of << std::setprecision(12) << std::setw(20) << model_hse(i, model::idens) << " ";
        of << std::setprecision(12) << std::setw(20) << model_hse(i, model::itemp) << " ";
        of << std::setprecision(12) << std::setw(20) << model_hse(i, model::ipres) << " ";
        if (write_ye) {
            of << std::setprecision(12) << std::setw(20) << model_hse(i, model::iyef) << " ";
        }
        for (int n = 0; n < NumSpec; ++n) {
            of << std::setprecision(12) << std::setw(20) << cfmt(model_hse(i, model::ispec+n)) << " ";
        }
        of << std::endl;
    }

    of.close();

}


AMREX_INLINE void
read_file(const std::string filename, initial_model_t& initial_model) {

    std::ifstream mf;

    mf.open(filename, std::ios::in);
    if (!mf.is_open()) {
        amrex::Error("Error opening the initial model");
    }

    std::string line;

    // first the header line -- this will determine the format.  If it
    // tells us the number of points (npts), then we are the old style.

    getline(mf, line);

    int nvars_model_file;
    std::vector<std::string> varnames_stored;

    if (line.find("npts") != std::string::npos) {
        // old format. example:
        /*
        # npts = 3
        # number of variables = 2
        # radius 
        # density 
        # temperature
                 0.1           3         4.0
                 0.2           4         3.9
                 0.3           4         3.7
        */
        
        // ignore the first line, and detect the number of points
        // automatically below

        // next line tells use the number of variables

        getline(mf, line);

        if (line.find("num of variables") == std::string::npos && line.find("number of variables") == std::string::npos) {
            amrex::Error("Error: missing the number of variables header");
        }

        std::string num_vars_string = line.substr(line.find("=")+1, line.length());
        nvars_model_file = std::stoi(num_vars_string);

        // now read in the names of the variables

        for (int n = 0; n < nvars_model_file; n++) {
            getline(mf, line);
            std::string var_string = line.substr(line.find("#")+1, line.length());
            varnames_stored.push_back(model_string::ltrim(model_string::rtrim(var_string)));
        }
    } else {
        // new format. example:
        /*
        # radius density temperature
                 0.1           3         4.0
                 0.2           4         3.9
                 0.3           4         3.7
        */

        // remove the leading "#"
        if (!line.empty() && line[0] == '#') {
            line = line.substr(1);
        }

        std::istringstream iss(line);

        // in this format, the first column name is the coordinate --
        // we can skip that

        std::string word;
        iss >> word;

        while (iss >> word) {
            varnames_stored.push_back(word);
        }

        nvars_model_file = static_cast<int>(varnames_stored.size());
    }


    // start reading in the data

    amrex::Vector<Real> vars_stored;
    vars_stored.resize(nvars_model_file);

    int i{0};
    while (getline(mf, line)) {
        if (i > NPTS_MODEL) {
            amrex::Error("Error: model has more than NPTS_MODEL points. Increase MAX_NPTS_MODEL");
        }

        std::istringstream iss(line);
        iss >> initial_model.r(i);

        for (int j = 0; j < nvars_model_file; j++) {
            iss >> vars_stored[j];
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

        } // end loop over variables in line

        i++;
    }  // end of loop over lines in model file

    initial_model.npts = i;

    mf.close();
}

#endif

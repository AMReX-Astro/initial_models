#ifndef INIT_1D_H
#define INIT_1D_H

//  The model is placed into HSE by the following differencing:
//
//   (1/dr) [ <P>_i - <P>_{i-1} ] = (1/2) [ <rho>_i + <rho>_{i-1} ] g
//
//  This will be iterated over in tandem with the EOS call,
//  P(i-1) = P_eos(rho(i-1), T(i-1), X(i-1)
//
#include <AMReX_Array.H>

#include <sstream>
#include <string>

#include <network.H>

#include <extern_parameters.H>
#include <fundamental_constants.H>
#include <interpolate.H>
#include <coord_info.H>
#include <read_model.H>
#include <model_util.H>

using namespace amrex;

constexpr Real TOL = 1.e-10_rt;

constexpr int MAX_ITER = 250;


AMREX_INLINE void init_1d( ) {

    //Creates the type initial_model_t, parse it into read_model
    //and returns with all the file information.

    initial_model_t lagrangian_planar;
    read_file(problem_rp::model_file, lagrangian_planar);

    //Apply the shift

    for (int i = 0; i < lagrangian_planar.npts; ++i) {
        lagrangian_planar.r(i) = lagrangian_planar.r(i) - problem_rp::model_shift + problem_rp::layer_size;
    }

    // Create a new uniform grid with problem_rp::nx cells.

    if (problem_rp::nx > NPTS_MODEL) {
        Error("too many zones requested -- increase NPTS_MODEL");
    }

    Array1D<Real, 0, NPTS_MODEL-1> xzn_hse;
    Array1D<Real, 0, NPTS_MODEL-1> xznl_hse;
    Array1D<Real, 0, NPTS_MODEL-1> xznr_hse;
    Array2D<Real, 0, NPTS_MODEL-1, 0, model::nvar-1> model_hse;

    fill_coord_arrays(xzn_hse, xznl_hse, xznr_hse);

    eos_t eos_state;

    for (int i = 0; i < problem_rp::nx; ++i) {
        for (int n = 0; n < model::nvar; ++n) {

            if (n == model::itemp) {
                model_hse(i,n) = amrex::max(problem_rp::temp_cutoff, interpolate(xzn_hse(i), n, lagrangian_planar, true));
            } else {
                model_hse(i,n) = interpolate(xzn_hse(i), n,  lagrangian_planar, true);
            }

        }


        //Now we have to make the thermodynamics of our model
        //consistent with the defined e.o.s for each site i.

        eos_state.rho = model_hse(i,model::idens);
        eos_state.T = model_hse(i, model::itemp);

        for (int n = 0; n < NumSpec; ++n){
            eos_state.xn[n] = model_hse(i, model::ispec+n);
        }

        eos(eos_input_rt, eos_state);

        model_hse(i, model::ipres) = eos_state.p;
    }

    write_model("uniform", xzn_hse, model_hse);

    // Once, the thermodynamical consistent pressure variable is
    // defined in our model, we pick the biggest T, r[i] state.

    int index_base = -1;
    Real max_T = -1.0_rt;

    for (int i = 0; i < problem_rp::nx; ++i) {
        if (model_hse(i, model::itemp) > max_T) {
            index_base = i;
            max_T = model_hse(i, model::itemp);
        }
    }

    if (index_base == -1) {
        Error("ERROR: invalid base height");
    }

    //and make sure again that all thermodyamics is consistent for
    //the initial choice.

    eos_state.rho = model_hse(index_base, model::idens);
    eos_state.T = model_hse(index_base, model::itemp);

    for (int n = 0; n < NumSpec; ++n) {
        eos_state.xn[n] = model_hse(index_base, model::ispec+n);
    }

    eos(eos_input_rt, eos_state);

    model_hse(index_base, model::ipres) = eos_state.p;

    //test1: ~5.0e7
    //std::cout << model_hse(index_base, model::itemp) << std::endl;

    //Now, we have to integrate up:

    Real xn[NumSpec];
    Real dens_zone;
    Real temp_zone;
    Real pres_zone;
    bool converged_hse;
    Real dpd;
    Real dpdr;
    Real g_zone;
    Real delx;
    Real pwant;
    Real drho;

    // Integrating Up
    for (int i = index_base+1; i < problem_rp::nx; ++i) {

        delx = xzn_hse(i) - xzn_hse(i-1);

        // Compute the gravitational acceleration
        if (problem_rp::do_invsq_grav == 1) {
            g_zone = -C::Gconst*problem_rp::M_enclosed/std::pow(xznl_hse(i), 2);
        } else {
            g_zone = problem_rp::g_const;
        }

        // Now, we have to setup our initial guess of rho, T, and
        // composition variables, for each i-cell.

        dens_zone = model_hse(i, model::idens);
        temp_zone = amrex::max(problem_rp::temp_cutoff, model_hse(i, model::itemp));

        for (int n = 0; n < NumSpec; ++n) {
            xn[n] = model_hse(i, model::ispec+n);
        }

        // Now, for each i-index point we force the thermodynamical variables
        // to enter in HSE

        converged_hse = false;

        for (int iter=0; iter < MAX_ITER; ++iter) {
            pwant = model_hse(i-1, model::ipres) + delx*0.5_rt*(dens_zone +
                model_hse(i-1, model::idens))*g_zone;

            eos_state.T = temp_zone;
            eos_state.rho = dens_zone;

            for (int n = 0; n < NumSpec; ++n) {
                eos_state.xn[n] = xn[n];
            }

            eos(eos_input_rt, eos_state);

            pres_zone = eos_state.p;

            dpd = eos_state.dpdr;
            drho = (pwant - pres_zone)/(dpd - 0.5_rt*delx*g_zone);

            dens_zone = amrex::max(0.9_rt*dens_zone,
                    amrex::min(dens_zone+drho, 1.1_rt*dens_zone));

            if (std::abs(drho) < TOL*dens_zone) {
                converged_hse = true;
                break;
            }

            if (dens_zone < problem_rp::low_density_cutoff) {
                dens_zone = problem_rp::low_density_cutoff;
                temp_zone = problem_rp::temp_cutoff;
                converged_hse = true;
                break;
            }

        }

        if (!converged_hse) {
            std::cout << "ERROR: Zone " << i << " did not converged in init_1d" <<std::endl;
            std::cout << "Integrate up" << std::endl;
            std::cout << dens_zone << " " << temp_zone << std::endl;
            std::cout << pwant << std::endl;
            std::cout << drho << std::endl;
            Error("ERROR: HSE non convergence");
        }


        //Now, we need to make the thermodynamics consistent again:
        eos_state.T = temp_zone;
        eos_state.rho = dens_zone;

        for (int n = 0; n < NumSpec; ++n) {
            eos_state.xn[n] = xn[n];
        }

        eos(eos_input_rt, eos_state);

        pres_zone = eos_state.p;

        model_hse(i, model::idens) = dens_zone;
        model_hse(i, model::itemp) = temp_zone;
        model_hse(i, model::ipres) = pres_zone;
    }

    converged_hse = false;

    //Now we have to integrate down:
    for (int i = index_base - 1; i >= 0; --i) {

        delx = xzn_hse(i+1) - xzn_hse(i);

        //Compute the gravitational acceleration
        if (problem_rp::do_invsq_grav == 1) {
            g_zone = -C::Gconst*problem_rp::M_enclosed/std::pow(xznr_hse(i), 2);
        } else {
            g_zone = problem_rp::g_const;
        }

        dens_zone = model_hse(i, model::idens);
        temp_zone = amrex::max(problem_rp::temp_cutoff, model_hse(i, model::itemp));

        for (int n = 0; n < NumSpec; n++) {
            xn[n] = model_hse(i, model::ispec+n);
        }

        //We start the Newton's loop
        for (int iter = 0; iter < MAX_ITER; ++iter) {

            pwant = model_hse(i+1, model::ipres) -
                delx*0.5*(dens_zone + model_hse(i+1, model::idens))*g_zone;

            eos_state.T = temp_zone;
            eos_state.rho = dens_zone;

            for (int n = 0; n < NumSpec; ++n) {
                eos_state.xn[n] = xn[n];
            }

            eos(eos_input_rt, eos_state);

            pres_zone = eos_state.p;

            dpd = eos_state.dpdr;

            drho = (pwant - pres_zone)/(dpd + 0.5_rt*delx*g_zone);

            dens_zone = amrex::max(0.9_rt*dens_zone, amrex::min(dens_zone + drho, 1.1_rt*dens_zone));

            if (std::abs(drho) < TOL*dens_zone) {
                converged_hse = true;
                break;
            }
        }

        if (!converged_hse) {
            std::cout << "ERROR: Zone " << i << " did not converged in init_1d" << std::endl;
            std::cout << "Integrate down" << std::endl;
            std::cout << dens_zone << " " << temp_zone << std::endl;
            std::cout << pwant << std::endl;
            std::cout << drho  << std::endl;
            Error("ERROR: HSE non convergence");
        }


        //Now, we enforce one more time the eos:

        eos_state.T = temp_zone;
        eos_state.rho = dens_zone;

        for (int n=0; n<NumSpec; ++n) {
            eos_state.xn[n] = xn[n];
        }

        eos(eos_input_rt, eos_state);

        pres_zone = eos_state.p;

        model_hse(i, model::idens) = dens_zone;
        model_hse(i, model::itemp) = temp_zone;
        model_hse(i, model::ipres) = pres_zone;
    }


    write_model("hse", xzn_hse, model_hse);

    Real max_hse_error = -1.e30_rt;
    Real rhog;

    for (int i = 1; i < problem_rp::nx-1; ++i) {

        if (problem_rp::do_invsq_grav == 1) {
            g_zone == -C::Gconst * problem_rp::M_enclosed/ std::pow(xznr_hse(i-1), 2);
        } else {
            g_zone = problem_rp::g_const;
        }

        dpdr = (model_hse(i, model::ipres) - model_hse(i-1, model::ipres)) / delx;
        rhog = 0.5_rt * (model_hse(i, model::idens) + model_hse(i-1, model::idens)) * g_zone;

        if (dpdr != 0.0_rt && model_hse(i+1, model::idens) > problem_rp::low_density_cutoff) {
            max_hse_error = max(max_hse_error, std::abs(dpdr - rhog) / std::abs(dpdr));
        }
    }

    std::cout << "maximum HSE error = " << max_hse_error << std::endl;
}
#endif

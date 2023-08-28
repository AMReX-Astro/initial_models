#ifndef INIT_1D_H
#define INIT_1D_H

#include <AMReX_Array.H>


#include <sstream>

#include <extern_parameters.H>
#include <fundamental_constants.H>

#include <coord_info.H>
#include <model_util.H>

// we use this only for the indices
#include <read_model.H>

using namespace amrex;


// generate an initial model for spherical geometry with a uniform
// composition.  Here we take a central density and temperature and
// use HSE and constant entropy to generate the model.  We specific
// the mass of the convective zone, and for radii beyond that, we
// use constant T.

AMREX_INLINE void init_1d() {

    // TOL_HSE is the tolerance used when iterating over a zone to
    // force it into HSE by adjusting the current density (and
    // possibly temperature).  TOL_HSE should be very small (~
    // 1.e-10).

    const Real TOL_HSE = 1.e-10_rt;

    const int MAX_ITER = 250;

    const Real smallx = 1.e-10_rt;

    // get the species indices

    int ic12  = network_spec_index("carbon-12");
    int io16  = network_spec_index("oxygen-16");

    if (ic12 < 0 || io16 < 0) {
        amrex::Error("ERROR: species not defined");
    }

    if (problem_rp::cfrac < 0.0_rt || problem_rp::cfrac > 1.0_rt) {
        amrex::Error("ERROR: cfrac must be between 0 and 1");
    }

    Real xn_base[NumSpec] = {problem_rp::smallx};
    xn_base[ic12] = problem_rp::cfrac;
    xn_base[io16] = 1.0_rt - problem_rp::cfrac;

    // Create a 1-d uniform grid that is identical to the mesh that we
    // are mapping onto, and then we want to force it into HSE on that
    // mesh.

    if (problem_rp::nx > NPTS_MODEL) {
        amrex::Error("too many zones requested -- increase NPTS_MODEL");
    }

    Array1D<Real, 0, NPTS_MODEL-1> xzn_hse;
    Array1D<Real, 0, NPTS_MODEL-1> xznl;
    Array1D<Real, 0, NPTS_MODEL-1> xznr;

    Array2D<Real, 0, NPTS_MODEL-1, 0, model::nvar-1> model_hse;

    Array1D<Real, 0, NPTS_MODEL-1> M_enclosed;
    Array1D<Real, 0, NPTS_MODEL-1> entropy_want;

    // compute the coordinates of the new gridded function

    fill_coord_arrays(xzn_hse, xznl, xznr);

    bool fluff{false};

    // call the EOS one more time for this zone and then go on to the next

    eos_t eos_state;
    eos_state.T = problem_rp::temp_base;
    eos_state.rho = problem_rp::dens_base;
    for (int n = 0; n < NumSpec; ++n) {
        eos_state.xn[n] = xn_base[n];
    }

    // (t, rho) -> (p, s)

    eos(eos_input_rt, eos_state);

    // make the initial guess be completely uniform

    for (int i = 0; i < problem_rp::nx; ++i) {
        model_hse(i, model::idens) = eos_state.rho;
        model_hse(i, model::itemp) = eos_state.T;
        model_hse(i, model::ipres) = eos_state.p;

        entropy_want(i) = eos_state.s;

        for (int n = 0; n < NumSpec; ++n) {
            model_hse(i, model::ispec+n) = eos_state.xn[n];
        }
    }

    // keep track of the mass enclosed below the current zone

    M_enclosed(0) = (4.0_rt / 3.0_rt) * M_PI *
            (std::pow(xznr(0), 3) - std::pow(xznl(0), 3)) * model_hse(0, model::idens);


    // HSE + entropy solve

    bool isentropic{true};

    Real dens_zone;
    Real temp_zone;
    Real pres_zone;
    Real entropy;
    Real xn[NumSpec];

    int i_conv{0};

    Real delx = xzn_hse(1) - xzn_hse(0);

    for (int i = 1; i < problem_rp::nx; ++i) {

        // as the initial guess for the temperature and density, use
        // the previous zone

        dens_zone = model_hse(i-1, model::idens);
        temp_zone = model_hse(i-1, model::itemp);
        for (int n = 0; n < NumSpec; ++n) {
            xn[n] = model_hse(i, model::ispec+n);
        }

        Real g_zone = -C::Gconst * M_enclosed(i-1) / (xznl(i) * xznl(i));


        // iteration loop

        // start off the Newton loop by saying that the zone has not converged
        bool converged_hse{false};

        if (! fluff) {

            Real p_want;
            Real drho;
            Real dtemp;
            Real entropy_base;

            for (int iter = 0; iter < MAX_ITER; ++iter) {

                if (isentropic) {

                    p_want = model_hse(i-1, model::ipres) +
                        delx * 0.5_rt * (dens_zone + model_hse(i-1, model::idens)) * g_zone;

                    eos_state.T = temp_zone;  // initial guess
                    eos_state.rho = dens_zone;  // initial guess
                    eos_state.p = p_want;
                    eos_state.s = entropy_want(i);

                    for (int n = 0; n < NumSpec; ++n) {
                        eos_state.xn[n] = xn[n];
                    }

                    eos(eos_input_ps, eos_state);

                    drho = eos_state.rho - dens_zone;
                    dens_zone = eos_state.rho;

                    dtemp = eos_state.T - temp_zone;
                    temp_zone = eos_state.T;

                    // check if the density falls below our minimum
                    // cut-off -- if so, floor it

                    if (dens_zone < problem_rp::dens_conv_zone) {
                        i_conv = i;
                        isentropic = false;
                    }

                    if (std::abs(drho) < TOL_HSE * dens_zone &&
                        std::abs(dtemp) < TOL_HSE * temp_zone) {
                        converged_hse = true;
                        break;
                    }

                } else {

                    // do isothermal

                    p_want = model_hse(i-1, model::ipres) +
                        delx * 0.5_rt * (dens_zone + model_hse(i-1, model::idens)) * g_zone;

                    temp_zone = model_hse(i-1, model::itemp);

                    eos_state.T = temp_zone;
                    eos_state.rho = dens_zone; // initial guess
                    eos_state.p = p_want;

                    for (int n = 0; n < NumSpec; ++n) {
                        eos_state.xn[n] = xn[n];
                    }

                    eos(eos_input_tp, eos_state);

                    drho = eos_state.rho - dens_zone;
                    dens_zone = eos_state.rho;

                    if (std::abs(drho) < TOL_HSE * dens_zone) {
                        converged_hse = true;
                        break;
                    }

                }

                if (dens_zone < problem_rp::low_density_cutoff) {
                    dens_zone = problem_rp::low_density_cutoff;
                    temp_zone = problem_rp::temp_fluff;
                    converged_hse = true;
                    fluff = true;
                    break;
                }

            }  // thermo iteration loop

            if (! converged_hse) {

                std::cout << "Error zone " << i <<  " did not converge in init_1d" << std::endl;
                std::cout << dens_zone << " " << temp_zone << std::endl;
                std::cout << p_want;
                std::cout << drho;
                amrex::Error("Error: HSE non-convergence");
            }

            if (temp_zone < problem_rp::temp_fluff) {
                temp_zone = problem_rp::temp_fluff;
                isentropic = false;
            }

        } else {
            dens_zone = problem_rp::low_density_cutoff;
            temp_zone = problem_rp::temp_fluff;
        }


        // call the EOS one more time for this zone and then go on to
        // the next

        eos_state.T = temp_zone;
        eos_state.rho = dens_zone;
        for (int n = 0; n < NumSpec; ++n) {
            eos_state.xn[n] = xn[n];
        }

        // (t, rho) -> (p, s)

        eos(eos_input_rt, eos_state);

        pres_zone = eos_state.p;

        // update the thermodynamics in this zone

        model_hse(i, model::idens) = dens_zone;
        model_hse(i, model::itemp) = temp_zone;
        model_hse(i, model::ipres) = pres_zone;

        M_enclosed(i) = M_enclosed(i-1) +
            (4.0_rt / 3.0_rt) * M_PI * (xznr(i) - xznl(i)) *
            (std::pow(xznr(i), 2) + xznl(i) * xznr(i) + std::pow(xznl(i), 2)) * model_hse(i, model::idens);

        if (M_enclosed(i) > problem_rp::M_conv_zone * C::M_solar && isentropic) {
            i_conv = i;
            isentropic = false;
        }

    } // end loop over zones


    std::cout << "mass = " << M_enclosed(problem_rp::nx-1) / C::M_solar << std::endl;;
    std::cout << "mass of convective zone = " << M_enclosed(i_conv) / C::M_solar << std::endl;;

    std::string outfile = problem_rp::prefix + ".hse";

    write_model(outfile, xzn_hse, model_hse);

    // compute the maximum HSE error

    Real max_hse_error = -1.e30_rt;

    for (int i = 1; i < problem_rp::nx-1; ++i) {
        Real g_zone = -C::Gconst * M_enclosed(i-1) / (xznr(i-1) * xznr(i-1));
        Real dpdr = (model_hse(i, model::ipres) - model_hse(i-1, model::ipres)) / delx;
        Real rhog = 0.5_rt * (model_hse(i, model::idens) + model_hse(i-1, model::idens)) * g_zone;

        if (dpdr != 0.0_rt && model_hse(i+1, model::idens) > problem_rp::low_density_cutoff) {
            max_hse_error = amrex::max(max_hse_error, std::abs(dpdr - rhog) / std::abs(dpdr));
        }
    }

    std::cout << "maximum HSE error = " << max_hse_error << std::endl;
}
#endif


#ifndef INIT_1D_H
#define INIT_1D_H

//  Create a 1-d hydrostatic, isoentropic atmosphere given the temperature,
//  composition, and base density.  This version allows for an entropy
//  decrease below a height base_height, to make the lower region
//  convectively stable.
//
//  The model is placed into HSE by the following differencing:
//
//   (1/dr) [ <P>_i - <P>_{i-1} ] = (1/2) [ <rho>_i + <rho>_{i-1} ] g
//
//   We can take <P>_base and <rho>_base as given, and use the above
//   to find all others.
//
//   This will be iterated over in tandem with the EOS call,
//   P(i-1) = P_eos(rho(i-1), T(i-1), X(i-1)
//

#include <AMReX_Array.H>


#include <sstream>

#include <extern_parameters.H>
#include <fundamental_constants.H>

#include <read_model.H>

using namespace amrex;


AMREX_INLINE void init_1d() {

    // TOL is the tolerance used when iterating over a zone to force
    // it into HSE by adjusting the current density (and possibly
    // temperature).  TOL should be very small (~ 1.e-10).

    const Real TOL = 1.e-10_rt;

    const int MAX_ITER = 250;

    const Real smallx = 1.e-10_rt;

    // get the species indices

    int ic12 = network_spec_index("carbon-12");
    int io16 = network_spec_index("oxygen-16");

    // ash can be either Mg24 or 'ash', depending on the network

    int iash = network_spec_index("magnesium-24");
    if (iash < 0) {
        iash = network_spec_index("ash");
    }

    if (ic12 < 0 || io16 < 0 || iash < 0) {
        amrex::Error("ERROR: species not defined");
    }

    Real xn_base[NumSpec] = {smallx};

    xn_base[ic12] = 0.3_rt;
    xn_base[io16] = 0.7_rt;


    // Create a 1-d uniform grid that is identical to the mesh that we
    // are mapping onto, and then we want to force it into HSE on that
    // mesh.

    if (problem_rp::nx > NPTS_MODEL) {
        amrex::Error("too many zones requested -- increase NPTS_MODEL");
    }

    Array1D<Real, 0, NPTS_MODEL-1> xzn_hse;

    Array2D<Real, 0, NPTS_MODEL-1, 0, model::nvar-1> model_hse;

    Array1D<Real, 0, NPTS_MODEL-1> entropy_want;
    Array1D<Real, 0, NPTS_MODEL-1> entropy_store;

    // compute the coordinates of the new gridded function

    Real dCoord = (problem_rp::xmax - problem_rp::xmin) / static_cast<Real>(problem_rp::nx);

    for (int i = 0; i < problem_rp::nx; ++i) {
        xzn_hse(i) = problem_rp::xmin + (static_cast<Real>(i) + 0.5_rt) * dCoord;
    }


    int index_base = -1;

    // find the index of the base height

    for (int i = 0; i < problem_rp::nx; ++i) {
        if (xzn_hse(i) >= problem_rp::base_height) {
            index_base = i;
            break;
        }
    }

    if (index_base == -1) {
        amrex::Error("ERROR: invalid base_height");
    }


    // put the model onto our new uniform grid

    bool fluff = false;

    // all the material above base_height will have constant entropy

    eos_t eos_state;
    eos_state.T = problem_rp::temp_base;
    eos_state.rho = problem_rp::dens_base;
    for (int n = 0; n < NumSpec; ++n) {
        eos_state.xn[n] = xn_base[n];
    }

    eos(eos_input_rt, eos_state);

    Real entropy = eos_state.s;

    // make the initial guess be completely uniform

    for (int i = 0; i < problem_rp::nx; ++i) {
        model_hse(i, model::idens) = eos_state.rho;
        model_hse(i, model::itemp) = eos_state.T;
        model_hse(i, model::ipres) = eos_state.p;

        for (int n = 0; n < NumSpec; ++n) {
            model_hse(i, model::ispec+n) = eos_state.xn[n];
        }
    }


    // set the desired entropy profile

    for (int i = index_base; i < problem_rp::nx; ++i) {
        // entropy is constant
        entropy_want(i) = entropy;
        entropy_store(i) = entropy_want(i);
    }


    Real slope = (entropy / problem_rp::entropy_jump - entropy / problem_rp::entropy_end) /
        (xzn_hse(index_base) - xzn_hse(0));

    for (int i = index_base-1; i >= 0; --i) {

        // entropy gradient
        entropy_want(i) = slope * (xzn_hse(i) - xzn_hse(index_base-1)) + entropy / problem_rp::entropy_jump;
        entropy_store(i) = entropy_want(i);
    }


    // HSE + entropy solve

    // the HSE state will be done putting creating an isentropic state
    // until the temperature goes below problem_rp::temp_fluff -- then we will do
    // isothermal.  also, once the density goes below
    // problem_rp::low_density_cutoff, we stop HSE

    bool isentropic = true;

    // integrate up

    Real dens_zone;
    Real temp_zone;
    Real pres_zone;
    Real xn[NumSpec];

    Real delx = xzn_hse(1) - xzn_hse(0);

    Real p_want;
    Real dtemp;
    Real drho;

    for (int i = index_base+1; i < problem_rp::nx; ++i) {

        // as the initial guess for the temperature and density, use the previous
        // zone

        dens_zone = model_hse(i-1, model::idens);
        temp_zone = model_hse(i-1, model::itemp);
        for (int n = 0; n < NumSpec; ++n) {
            xn[n] = model_hse(i, model::ispec+n);
        }

        // iteration loop

        // start off the Newton loop by saying that the zone has not converged

        bool converged_hse = false;

        if (! fluff) {

            for (int iter = 0; iter < MAX_ITER; ++iter) {

                if (isentropic) {

                    // get the pressure we want from the HSE equation,
                    // just the zone below the current.  Note, we are
                    // using an average of the density of the two
                    // zones as an approximation of the interface
                    // value -- this means that we need to iterate for
                    // find the density and pressure that are
                    // consistent

                    // furthermore, we need to get the entropy that we
                    // need, which will come from adjusting the
                    // temperature in addition to the density.

                    // HSE differencing

                    p_want = model_hse(i-1, model::ipres) +
                        delx * 0.5_rt * (dens_zone + model_hse(i-1, model::idens)) * problem_rp::g_zone;

                    // now we have two functions to zero:
                    //   A = p_want - p(rho,T)
                    //   B = entropy_want - s(rho,T)
                    // We use a two dimensional Taylor expansion and
                    // find the deltas for both density and
                    // temperature

                    // (t, rho) -> (p, s)

                    eos_state.T = temp_zone;
                    eos_state.rho = dens_zone;
                    for (int n = 0; n < NumSpec; ++n) {
                        eos_state.xn[n] = xn[n];
                    }

                    eos(eos_input_rt, eos_state);

                    entropy = eos_state.s;
                    pres_zone = eos_state.p;

                    Real dpT = eos_state.dpdT;
                    Real dpd = eos_state.dpdr;
                    Real dsT = eos_state.dsdT;
                    Real dsd = eos_state.dsdr;

                    Real A = p_want - pres_zone;
                    Real B = entropy_want(i) - entropy;

                    dtemp = ((dsd / (dpd - 0.5_rt * delx * problem_rp::g_zone)) * A - B) /
                        (dsd * dpT / (dpd -0.5_rt * delx * problem_rp::g_zone) - dsT);

                    drho = (A - dpT * dtemp) / (dpd - 0.5_rt * delx * problem_rp::g_zone);

                    dens_zone =
                        amrex::max(0.9_rt * dens_zone,
                                   amrex::min(dens_zone + drho, 1.1_rt * dens_zone));

                    temp_zone =
                        amrex::max(0.9_rt * temp_zone,
                                   amrex::min(temp_zone + dtemp, 1.1_rt * temp_zone));


                    // check if the density falls below our minimum cut-off --
                    // if so, floor it

                    if (dens_zone < problem_rp::low_density_cutoff) {
                        dens_zone = problem_rp::low_density_cutoff;
                        temp_zone = problem_rp::temp_fluff;
                        converged_hse = true;
                        fluff = true;
                        break;
                    }

                    if (std::abs(drho) < TOL * dens_zone && std::abs(dtemp) < TOL * temp_zone) {
                        converged_hse = true;
                        break;
                    }

                } else {

                    // do isothermal

                    Real p_want = model_hse(i-1, model::ipres) +
                        delx * 0.5_rt * (dens_zone + model_hse(i-1, model::idens)) * problem_rp::g_zone;

                    temp_zone = problem_rp::temp_fluff;

                    // (t, rho) -> (p)

                    eos_state.T = temp_zone;
                    eos_state.rho = dens_zone;
                    for (int n = 0; n < NumSpec; ++n) {
                        eos_state.xn[n] = xn[n];
                    }

                    eos(eos_input_rt, eos_state);

                    entropy = eos_state.s;
                    pres_zone = eos_state.p;

                    Real dpd = eos_state.dpdr;

                    Real drho = (p_want - pres_zone) / (dpd - 0.5_rt * delx * problem_rp::g_zone);

                    dens_zone =
                        amrex::max(0.9_rt * dens_zone,
                                   amrex::min(dens_zone + drho, 1.1_rt * dens_zone));

                    if (std::abs(drho) < TOL * dens_zone) {
                        converged_hse = true;
                        break;
                    }

                    if (dens_zone < problem_rp::low_density_cutoff) {
                        dens_zone = problem_rp::low_density_cutoff;
                        temp_zone = problem_rp::temp_fluff;
                        converged_hse = true;
                        fluff = true;
                        break;
                    }

                } // isentropic

            }

            if (! converged_hse) {
                std::cout << "Error zone " << i <<  " did not converge in init_1d" << std::endl;
                std::cout << "integrate up" << std::endl;
                std::cout << dens_zone << " " << temp_zone << std::endl;
                std::cout << p_want << " " << entropy_want(i) << " " << entropy << std::endl;
                std::cout << drho << " " << dtemp << std::endl;
                amrex:Error("Error: HSE non-convergence");
            }

            if (temp_zone < problem_rp::temp_fluff) {
                temp_zone = problem_rp::temp_fluff;
                isentropic = false;
            }

        } else {
            dens_zone = problem_rp::low_density_cutoff;
            temp_zone = problem_rp::temp_fluff;
        }


        // call the EOS one more time for this zone and then go on to the next

        // (t, rho) -> (p)

        eos_state.T = temp_zone;
        eos_state.rho = dens_zone;
        for (int n = 0; n < NumSpec; ++n) {
            eos_state.xn[n] = xn[n];
        }

        eos(eos_input_rt, eos_state);

        pres_zone = eos_state.p;

        // update the thermodynamics in this zone

        model_hse(i, model::idens) = dens_zone;
        model_hse(i, model::itemp) = temp_zone;
        model_hse(i, model::ipres) = pres_zone;
        entropy_store(i) = entropy;

        // to make this process converge faster, set the density in the next zone to
        // the density in this zone

        // model_hse(i+1, idens) = dens_zone;

    }


    // integrate down

    for (int i = index_base-1; i >= 0; --i) {

        // as the initial guess for the temperature and density, use
        // the previous zone

        dens_zone = model_hse(i+1, model::idens);
        temp_zone = model_hse(i+1, model::itemp);
        for (int n = 0; n < NumSpec; ++n) {
            xn[n] = model_hse(i, model::ispec+n);
        }


        // iteration loop

        // start off the Newton loop by saying that the zone has not converged

        bool converged_hse = false;

        for (int iter = 0; iter < MAX_ITER; ++iter) {

            // get the pressure we want from the HSE equation, just
            // the zone below the current.  Note, we are using an
            // average of the density of the two zones as an
            // approximation of the interface value -- this means that
            // we need to iterate for find the density and pressure
            // that are consistent

            // furthermore, we need to get the entropy that we need,
            // which will come from adjusting the temperature in
            // addition to the density.

            // HSE differencing

            p_want = model_hse(i+1, model::ipres) -
                delx * 0.5_rt * (dens_zone + model_hse(i+1, model::idens)) * problem_rp::g_zone;

            // now we have two functions to zero:
            //   A = p_want - p(rho,T)
            //   B = entropy_want - s(rho,T)
            // We use a two dimensional Taylor expansion and find the
            // deltas for both density and temperature

            // (t, rho) -> (p, s)

            eos_state.T = temp_zone;
            eos_state.rho = dens_zone;
            for (int n = 0; n < NumSpec; ++n) {
                eos_state.xn[n] = xn[n];
            }

            eos(eos_input_rt, eos_state);

            entropy = eos_state.s;
            pres_zone = eos_state.p;

            Real dpT = eos_state.dpdT;
            Real dpd = eos_state.dpdr;
            Real dsT = eos_state.dsdT;
            Real dsd = eos_state.dsdr;

            Real A = p_want - pres_zone;
            Real B = entropy_want(i) - entropy;

            dtemp = ((dsd / (dpd + 0.5_rt * delx * problem_rp::g_zone)) * A - B) /
                (dsd * dpT / (dpd + 0.5_rt * delx * problem_rp::g_zone) - dsT);

            drho = (A - dpT * dtemp) / (dpd + 0.5_rt * delx * problem_rp::g_zone);

            dens_zone =
                amrex::max(0.9_rt * dens_zone,
                           amrex::min(dens_zone + drho, 1.1_rt * dens_zone));

            temp_zone =
                amrex::max(0.9_rt * temp_zone,
                           amrex::min(temp_zone + dtemp, 1.1_rt * temp_zone));

            if (std::abs(drho) < TOL * dens_zone && std::abs(dtemp) < TOL * temp_zone) {
                converged_hse = true;
                break;
            }

        }

        if (! converged_hse) {
            std::cout << "Error zone " << i << " did not converge in init_1d" << std::endl;
            std::cout << "integrate down" << std::endl;
            std::cout << dens_zone << " " << temp_zone << std::endl;
            std::cout << p_want << " " << entropy_want(i) << " " << entropy << std::endl;
            std::cout << drho << " " << dtemp << std::endl;
            amrex::Error("Error: HSE non-convergence");
        }

        // call the EOS one more time for this zone and then go on to the next

        // (t, rho) -> (p)

        eos_state.T = temp_zone;
        eos_state.rho = dens_zone;
        for (int n = 0; n < NumSpec; ++n) {
            eos_state.xn[n] = xn[n];
        }

        eos(eos_input_rt, eos_state);

        pres_zone = eos_state.p;

        // update the thermodynamics in this zone

        model_hse(i, model::idens) = dens_zone;
        model_hse(i, model::itemp) = temp_zone;
        model_hse(i, model::ipres) = pres_zone;
        entropy_store(i) = entropy;

    }

    write_model("model.hse", xzn_hse, model_hse);

    // compute the maximum HSE error

    Real max_hse_error = -1.e30_rt;

    for (int i = 1; i < problem_rp::nx-1; ++i) {
        Real dpdr = (model_hse(i, model::ipres) - model_hse(i-1, model::ipres)) / delx;
        Real rhog = 0.5_rt * (model_hse(i, model::idens) + model_hse(i-1, model::idens)) * problem_rp::g_zone;

        if (dpdr != 0.0_rt && model_hse(i+1, model::idens) > problem_rp::low_density_cutoff) {
            max_hse_error = amrex::max(max_hse_error, std::abs(dpdr - rhog) / std::abs(dpdr));
        }
    }

    std::cout << "maximum HSE error = " << max_hse_error << std::endl;
}
#endif

#ifndef INIT_1D_H
#define INIT_1D_H

#include <AMReX_Array.H>

#include <sstream>

#include <extern_parameters.H>
#include <fundamental_constants.H>

#include <coord_info.H>
#include <model_util.H>

// we only use the model namespace from here
#include <read_model.H>

using namespace amrex;


// generate an initial model for spherical geometry with a uniform
// composition.  Here we take a central density and temperature and
// use HSE and constant entropy to generate the model.  We specific
// the mass of the convective zone, and for radii beyond that, we
// use constant T.

// we allow for a non-uniform grid spacing, following the ideas in
// the MAESTROeX paper for mapping a spherical domain into a 1-d radial
// array.

AMREX_INLINE void init_1d() {

    // TOL_HSE is the tolerance used when iterating over a zone to
    // force it into HSE by adjusting the current density (and
    // possibly temperature).  TOL_HSE should be very small (~
    // 1.e-10).

    const Real TOL_HSE = 1.e-10_rt;

    const int MAX_ITER = 250;

    const Real TOL_MASS = 1.e-6_rt;

    // get the species indices

    const int ih1  = network_spec_index("hydrogen-1");
    const int ihe4  = network_spec_index("helium-4");

    if (ih1 < 0 || ihe4 < 0) {
        amrex::Error("ERROR: species not defined");
    }

    if (problem_rp::hefrac < 0.0_rt || problem_rp::hefrac > 1.0_rt) {
        amrex::Error("ERROR: hefrac must be between 0 and 1");
    }

    Real xn_core[NumSpec] = {0.0};
    xn_core[ihe4] = problem_rp::hefrac;
    xn_core[ih1] = 1.0_rt - problem_rp::hefrac;

    // Create a 1-d uniform grid that is identical to the mesh that we
    // are mapping onto, and then we want to force it into HSE on that
    // mesh.

    int nr = get_irreg_nr();

    if (nr > NPTS_MODEL) {
        amrex::Error("too many zones requested -- increase NPTS_MODEL");
    }

    Array1D<Real, 0, NPTS_MODEL-1> xzn_hse;
    Array1D<Real, 0, NPTS_MODEL-1> xznl;
    Array1D<Real, 0, NPTS_MODEL-1> xznr;

    Array2D<Real, 0, NPTS_MODEL-1, 0, model::nvar-1> model_hse;

    Array1D<Real, 0, NPTS_MODEL-1> M_enclosed;
    Array1D<Real, 0, NPTS_MODEL-1> entropy_want;

    // these are only needed for an irregular grid

    Array1D<Real, 0, NPTS_MODEL-1> delrl;
    Array1D<Real, 0, NPTS_MODEL-1> delrr;

    fill_coord_arrays_irreg(nr, xzn_hse, xznl, xznr, delrl, delrr);

    // We don't know what central density will give the desired total
    // mass, so we need to iterate over central density

    // we will do a secant iteration.  rho_c_old is the 'old' guess for
    // the central density and rho_c is the current guess.  After 2
    // loops, we can start estimating the density required to yield our
    // desired mass

    Real rho_c_old{-1.0_rt};
    Real rho_c{1.e3_rt};   // 1.e3 is a reasonable starting low mass star density
    Real mass_star{-1.0};
    Real mass_star_old{-1.0};

    bool mass_converged{false};

    for (int iter_mass = 0; iter_mass < MAX_ITER; ++iter_mass) {

        //std::cout << "mass iter = " << iter_mass << " " << rho_c << " " << problem_rp::temp_core << std::endl;

        bool fluff{false};


        // call the EOS one more time for this zone and then go on to the next

        eos_t eos_state;
        eos_state.T = problem_rp::temp_core;
        eos_state.rho = rho_c;
        for (int n = 0; n < NumSpec; ++n) {
            eos_state.xn[n] = xn_core[n];
        }

        // (t, rho) -> (p, s)

        eos(eos_input_rt, eos_state);

        // make the initial guess be completely uniform

        for (int i = 0; i < nr; ++i) {

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

        for (int i = 1; i < nr; ++i) {

            Real delx{0};
            Real rfrac{0};
            if (problem_rp::use_irreg_grid) {
                delx = delrr(i-1) + delrl(i);
                rfrac = delrr(i-1) / delx;
            } else {
                delx = xzn_hse(i) - xzn_hse(i-1);
                rfrac = 0.5;
            }

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

                for (int iter = 0; iter < MAX_ITER; ++iter) {

                    if (isentropic) {

                        p_want = model_hse(i-1, model::ipres) +
                            delx * ((1.0_rt - rfrac) * dens_zone + rfrac * model_hse(i-1, model::idens)) * g_zone;

                        // now we have two functions to zero:
                        //   A = p_want - p(rho,T)
                        //   B = entropy_want - s(rho,T)
                        // We use a two dimensional Taylor expansion
                        // and find the deltas for both density and
                        // temperature

                        eos_state.T = temp_zone;
                        eos_state.rho = dens_zone;
                        for (int n = 0; n < NumSpec; ++n) {
                            eos_state.xn[n] = xn[n];
                        }

                        // (t, rho) -> (p, s)
                        eos(eos_input_rt, eos_state);

                        entropy = eos_state.s;
                        pres_zone = eos_state.p;

                        Real dpT = eos_state.dpdT;
                        Real dpd = eos_state.dpdr;
                        Real dsT = eos_state.dsdT;
                        Real dsd = eos_state.dsdr;

                        Real A = p_want - pres_zone;
                        Real B = entropy_want(i) - entropy;

                        Real dAdT = -dpT;
                        Real dAdrho = (1.0_rt - rfrac) * delx * g_zone - dpd;
                        Real dBdT = -dsT;
                        Real dBdrho = -dsd;

                        dtemp = (B - (dBdrho / dAdrho) * A) /
                            ((dBdrho / dAdrho) * dAdT - dBdT);

                        drho = -(A + dAdT * dtemp) / dAdrho;

                        dens_zone =
                            amrex::max(0.9_rt * dens_zone,
                                       amrex::min(dens_zone + drho, 1.1_rt * dens_zone));

                        temp_zone =
                            amrex::max(0.9_rt * temp_zone,
                                       amrex::min(temp_zone + dtemp, 1.1_rt * temp_zone));

                        // check if the density falls below our minimum
                        // cut-off -- if so, floor it

                        if (dens_zone < problem_rp::low_density_cutoff) {

                            dens_zone = problem_rp::low_density_cutoff;
                            temp_zone = problem_rp::temp_fluff;
                            converged_hse = true;
                            fluff = true;
                            break;
                        }

                        if (std::abs(drho) < TOL_HSE * dens_zone &&
                            std::abs(dtemp) < TOL_HSE * temp_zone) {
                            converged_hse = true;
                            break;
                        }

                    } else {

                        // do isothermal

                        p_want = model_hse(i-1, model::ipres) +
                            delx * ((1.0_rt - rfrac) * dens_zone + rfrac * model_hse(i-1, model::idens)) * g_zone;

                        temp_zone = model_hse(i-1, model::itemp);

                        eos_state.T = temp_zone;
                        eos_state.rho = dens_zone;
                        for (int n = 0; n < NumSpec; ++n) {
                            eos_state.xn[n] = xn[n];
                        }

                        // (t, rho) -> (p, s)

                        eos(eos_input_rt, eos_state);

                        entropy = eos_state.s;
                        pres_zone = eos_state.p;

                        Real dpd = eos_state.dpdr;

                        drho = (p_want - pres_zone) / (dpd - (1.0_rt - rfrac) * delx * g_zone);

                        dens_zone =
                            amrex::max(0.9_rt * dens_zone,
                                       amrex::min(dens_zone + drho, 1.1_rt * dens_zone));

                        if (std::abs(drho) < TOL_HSE * dens_zone) {
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

        } // end loop over zones

        mass_star = M_enclosed(nr-1);

        std::cout << "mass = " << mass_star / C::M_solar << " central density = " << rho_c << std::endl;;

        if (rho_c_old < 0.0_rt) {
            // not enough iterations yet -- store the old central density and
            // mass and pick a new value
            rho_c_old = rho_c;
            mass_star_old = mass_star;

            rho_c = 0.5 * rho_c_old;

        } else {
            // have we converged
            if (std::abs(mass_star - problem_rp::M_tot * C::M_solar) /
                (problem_rp::M_tot * C::M_solar) < TOL_MASS) {
                mass_converged = true;
                break;
            }

            // do a secant iteration:
            // M_tot = M(rho_c) + dM/drho |_rho_c x drho + ...
            Real drho_c = (problem_rp::M_tot * C::M_solar - mass_star) /
                ((mass_star - mass_star_old) / (rho_c - rho_c_old));

            rho_c_old = rho_c;
            mass_star_old = mass_star;

            rho_c = std::min(1.5_rt * rho_c_old,
                             std::max((rho_c + drho_c), 0.5_rt * rho_c_old));

        }

    }  // end mass constraint loop

    if (! mass_converged) {
        amrex::Error("ERROR: mass did not converge");
    }


    std::cout << "mass = " << M_enclosed(nr-1) / C::M_solar << std::endl;;

    // output

    std::string out_base = problem_rp::prefix + ".hse";

    write_model(out_base, xzn_hse, model_hse);

    // compute the maximum HSE error

    Real max_hse_error = -1.e30_rt;

    for (int i = 1; i < nr; ++i) {
        Real g_zone = -C::Gconst * M_enclosed(i-1) / (xznr(i-1) * xznr(i-1));

        Real delx;
        Real rfrac;
        if (problem_rp::use_irreg_grid) {
            delx = delrl(i) + delrr(i-1);
            rfrac = delrr(i-1) / delx;
        } else {
            delx = xznr(i) - xznl(i);
            rfrac = 0.5;
        }

        Real dpdr = (model_hse(i, model::ipres) - model_hse(i-1, model::ipres)) / delx;
        Real rhog = ((1.0 - rfrac) * model_hse(i, model::idens) + rfrac * model_hse(i-1, model::idens)) * g_zone;

        if (dpdr != 0.0_rt && model_hse(i+1, model::idens) > problem_rp::low_density_cutoff) {
            max_hse_error = amrex::max(max_hse_error, std::abs(dpdr - rhog) / std::abs(dpdr));
        }
    }

    std::cout << "maximum HSE error = " << max_hse_error << std::endl;
}
#endif


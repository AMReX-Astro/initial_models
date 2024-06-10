#ifndef INIT_1D_H
#define INIT_1D_H

#include <AMReX_Array.H>


#include <sstream>

#include <extern_parameters.H>
#include <fundamental_constants.H>

#include <coord_info.H>
#include <read_model.H>
#include <model_util.H>

using namespace amrex;


// generate an initial model for an arbitrary-mass, isothermal C WD
// with an isentropic He envelope on the surface.

AMREX_INLINE void init_1d() {

    // TOL_WD_MASS is tolerance used for getting the total WD mass
    // equal to M_tot (defined below).  It can be reasonably small,
    // since there will always be a central density value that can
    // give the desired WD mass on the grid we use

    const Real TOL_WD_MASS = 1.e-3_rt;

    // TOL_HE_MASS is the tolerance used for getting the mass of the
    // He envelope.  This should not be too small, since the values of
    // the He envelope mass we can achieve will depend on our grid
    // spacing.

    const Real TOL_HE_MASS = 1.e-3_rt;

    const int MAX_ITER = 1000;


    // convert the envelope and WD mass into CGS

    problem_rp::M_tot *= C::M_solar;
    problem_rp::M_He *= C::M_solar;

    // get the species indices

    int ihe4 = network_spec_index("helium-4");
    int ic12 = network_spec_index("carbon-12");
    int io16 = network_spec_index("oxygen-16");

    int in14 = network_spec_index("nitrogen-14");

    if (ihe4 < 0 || ic12 < 0 || io16 < 0) {
        amrex::Error("ERROR: species not defined");
    }

    if (in14 < 0 && problem_rp::X_N14 > 0.0_rt) {
        amrex::Error("ERROR: N14 not defined");
    }

    Real xn_core[NumSpec] = {network_rp::small_x};
    Real xn_he[NumSpec] = {network_rp::small_x};;

    if (problem_rp::mixed_co_wd) {
        xn_core[ic12] = 0.5_rt - 0.5_rt * (NumSpec - 1) * network_rp::small_x;
        xn_core[io16] = 0.5_rt - 0.5_rt * (NumSpec - 1) * network_rp::small_x;
    } else {
        xn_core[ic12] = 1.0_rt - (NumSpec - 1) * network_rp::small_x;
    }

    if (problem_rp::X_N14 > 0.0_rt) {
        xn_he[in14] = problem_rp::X_N14;
    }

    if (problem_rp::X_C12 > 0.0_rt) {
        xn_he[ic12] = problem_rp::X_C12;
    }

    if (problem_rp::X_O16 > 0.0_rt) {
        xn_he[io16] = problem_rp::X_O16;
    }

    xn_he[ihe4] = 1.0_rt - xn_he[ic12] - xn_he[in14] - xn_he[io16] - (NumSpec-4) * network_rp::small_x;

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
    Array1D<Real, 0, NPTS_MODEL-1> cs_hse;
    Array1D<Real, 0, NPTS_MODEL-1> s_hse;

    // compute the coordinates of the new gridded function

    fill_coord_arrays(xzn_hse, xznl, xznr);

    // We don't know what WD central density will give the desired
    // total mass, so we need to iterate over central density

    // we will do a secant iteration.  rho_c_old is the 'old' guess
    // for the central density and rho_c is the current guess.  After
    // 2 loops, we can start estimating the density required to yield
    // our desired mass

    Real rho_c_old = -1.0_rt;
    Real rho_c = 1.e9_rt;  // 1.e9 is a reasonable starting WD central density

    // rho_he_old is the old guess for the density to transition to He,
    // where we will be isentropic, and rho_he is the current guess.

    Real rho_he_old = -1.0_rt;
    Real rho_he = 0.5_rt * rho_c;

    bool mass_converged = false;

    Real mass_wd, mass_wd_old;
    Real mass_he, mass_he_old;

    int ihe_layer;
    int ihe_entropy;
    int icutoff;

    Real delx = xzn_hse(1) - xzn_hse(0);

    for (int iter_mass = 0; iter_mass < MAX_ITER; ++iter_mass) {

        std::cout << "mass iter = " << iter_mass << " " << rho_c << " " << problem_rp::temp_core << std::endl;

        bool fluff = false;

        // we start at the center of the WD and integrate outward.
        // Initialize the central conditions.

        eos_t eos_state;
        eos_state.T = problem_rp::temp_core;
        eos_state.rho = rho_c;
        for (int n = 0; n < NumSpec; ++n) {
            eos_state.xn[n] = xn_core[n];
        }

        // (t, rho) -> (p, s)

        eos(eos_input_rt, eos_state);

        // make the initial guess be completely uniform

        for (int i = 0; i < problem_rp::nx; ++i) {
            model_hse(i, model::idens) = eos_state.rho;
            model_hse(i, model::itemp) = eos_state.T;
            model_hse(i, model::ipres) = eos_state.p;

            cs_hse(i) = eos_state.cs;
            s_hse(i) = eos_state.s;

            for (int n = 0; n < NumSpec; ++n) {
                model_hse(i, model::ispec+n) = eos_state.xn[n];
            }
        }


        // keep track of the mass enclosed below the current zone

        M_enclosed(0) = (4.0_rt / 3.0_rt) * M_PI *
            (std::pow(xznr(0), 3) - std::pow(xznl(0), 3)) * model_hse(0, model::idens);

        ihe_layer = -1;
        ihe_entropy = -1;
        icutoff = -1;

        // HSE + entropy solve

        bool isentropic;

        Real entropy_base;

        Real dens_zone;
        Real temp_zone;
        Real pres_zone;
        Real entropy;
        Real xn[NumSpec];

        for (int i = 1; i < problem_rp::nx; ++i) {

            // as the initial guess for the density, use the previous zone

            dens_zone = model_hse(i-1, model::idens);

            if (dens_zone > rho_he) {
                temp_zone = problem_rp::temp_core;
                for (int n = 0; n < NumSpec; ++n) {
                    xn[n] = xn_core[n];
                }

                isentropic = false;

            } else {

                if (ihe_layer == -1) {
                    ihe_layer = i;
                }


                // determine whether we are starting the ramp up.  We
                // will use a tanh profile, centered at
                // (xzn_hse(ihe_layer) + FOUR*delta).  The "+
                // FOUR*delta" enables us to capture the leading edge
                // of the profile.  Since rho_he is computed by
                // considering the integral of He on the grid,
                // shifting the profile by FOUR*delta doesn't affect
                // the overall mass.

                Real test = 0.5_rt * (1.0_rt + std::tanh((xzn_hse(i) - xzn_hse(ihe_layer) - 4.0_rt * problem_rp::delta) / problem_rp::delta));

                if (test < 0.999_rt) {

                    // small tanh ramp up regime

                    for (int n = 0; n < NumSpec; ++n) {
                        xn[n] = xn_core[n] + 0.5_rt * (xn_he[n] - xn_core[n]) *
                            (1.0_rt + std::tanh((xzn_hse(i) - xzn_hse(ihe_layer) - 4.0_rt * problem_rp::delta) / problem_rp::delta));
                    }

                    temp_zone = problem_rp::temp_core + 0.5_rt * (problem_rp::temp_base - problem_rp::temp_core)*
                        (1.0_rt + std::tanh((xzn_hse(i) - xzn_hse(ihe_layer) - 4.0_rt * problem_rp::delta) / problem_rp::delta));

                    isentropic = false;

                } else {

                    if (problem_rp::isothermal_layer) {
                        // isothermal He layer no matter what
                        temp_zone = problem_rp::temp_base;
                        isentropic = false;
                    } else {
                        // fully isentropic

                        if (ihe_entropy == -1) {
                            ihe_entropy = i;
                            temp_zone = problem_rp::temp_base;
                            isentropic = false;
                        } else {
                            temp_zone = model_hse(i-1, model::itemp);
                            isentropic = true;
                        }
                    }

                    for (int n = 0; n < NumSpec; ++n) {
                        xn[n] = xn_he[n];
                    }

                }

            }

            Real g_zone = -C::Gconst * M_enclosed(i-1) / (xznl(i) * xznl(i));


            // thermodynamic state iteration loop

            // start off the Newton loop by saying that the zone has not converged
            bool converged_hse = false;

            if (! fluff) {

                Real p_want;
                Real drho;
                Real dtemp;

                for (int iter = 0; iter < MAX_ITER; ++iter) {

                    if (isentropic) {

                        p_want = model_hse(i-1, model::ipres) +
                            delx * 0.5_rt * (dens_zone + model_hse(i-1, model::idens)) * g_zone;

                        // now we have two functions to zero:
                        //   A = p_want - p(rho,T)
                        //   B = entropy_base - s(rho,T)
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
                        Real B = entropy_base - entropy;

                        Real dAdT = -dpT;
                        Real dAdrho = 0.5_rt * delx * g_zone - dpd;
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

                        if (std::abs(drho) < problem_rp::tol_hse * dens_zone &&
                            std::abs(dtemp) < problem_rp::tol_hse * temp_zone) {
                            converged_hse = true;
                            break;
                        }

                    } else {

                        // the core is isothermal, so we just need to
                        // constrain the density and pressure to agree
                        // with the EOS and HSE

                        // We difference HSE about the interface
                        // between the current zone and the one just
                        // inside.

                        p_want = model_hse(i-1, model::ipres) +
                            delx * 0.5_rt * (dens_zone + model_hse(i-1, model::idens)) * g_zone;

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

                        drho = (p_want - pres_zone) / (dpd - 0.5_rt *delx * g_zone);

                        dens_zone =
                            amrex::max(0.9_rt * dens_zone,
                                       amrex::min(dens_zone + drho, 1.1_rt * dens_zone));

                        if (std::abs(drho) < problem_rp::tol_hse * dens_zone) {
                            converged_hse = true;
                            break;
                        }

                        if (dens_zone < problem_rp::low_density_cutoff) {

                            icutoff = i;
                            dens_zone = problem_rp::low_density_cutoff;
                            temp_zone = problem_rp::temp_fluff;
                            converged_hse = true;
                            fluff = true;
                            break;
                        }
                    }

                    if (temp_zone < problem_rp::temp_fluff && isentropic) {
                        temp_zone = problem_rp::temp_fluff;
                        isentropic = false;
                    }

                }  // thermo iteration loop


                if (! converged_hse) {

                    std::cout << "Error zone " << i <<  " did not converge in init_1d" << std::endl;
                    std::cout << dens_zone << " " << temp_zone << std::endl;
                    std::cout << p_want << std::endl;
                    std::cout << drho << std::endl;
                    amrex::Error("Error: HSE non-convergence");
                }

            } else {
                // fluff region
                dens_zone = problem_rp::low_density_cutoff;
                temp_zone = problem_rp::temp_fluff;
            }

            // call the EOS one more time for this zone and then go on
            // to the next

            eos_state.T = temp_zone;
            eos_state.rho = dens_zone;
            for (int n = 0; n < NumSpec; ++n) {
                eos_state.xn[n] = xn[n];
            }

            // (t, rho) -> (p, s)

            eos(eos_input_rt, eos_state);

            pres_zone = eos_state.p;

            // determine the entropy that we want to constrain to, if
            // this is the first zone of the He layer

            if (i == ihe_entropy) {
                entropy_base = entropy;
            }

            // update the thermodynamics in this zone

            model_hse(i, model::idens) = dens_zone;
            model_hse(i, model::itemp) = temp_zone;
            model_hse(i, model::ipres) = pres_zone;

            for (int n = 0; n < NumSpec; ++n) {
                model_hse(i, model::ispec+n) = xn[n];
            }

            M_enclosed(i) = M_enclosed(i-1) +
                (4.0_rt / 3.0_rt) * M_PI * (xznr(i) - xznl(i)) *
                (std::pow(xznr(i), 2) + xznl(i) * xznr(i) + std::pow(xznl(i), 2)) * model_hse(i, model::idens);

            cs_hse(i) = eos_state.cs;
            s_hse(i) = eos_state.s;

        } // end loop over zones

        mass_he = 0.0;
        mass_wd = 0.0;

        // it might be that we never reach the cutoff density in our
        // domain.  This is especially the case if we do an isothermal
        // model.  Make sure we integrate over everything in that
        // case.
        int max_index = icutoff == -1 ? problem_rp::nx : icutoff;

        for (int i = 0; i < max_index; ++i) {

            Real vol{0.0};
            if (i == 0) {
                vol = (4.0_rt / 3.0_rt) * M_PI * (std::pow(xznr(0), 3) - std::pow(xznl(0), 3));
            } else {
                vol = (4.0_rt / 3.0_rt) * M_PI *
                (xznr(i) - xznl(i)) * (std::pow(xznr(i), 2) + xznl(i) * xznr(i) + std::pow(xznl(i), 2));
            }

            // compute the total mass of the He layer and C/O WD
            // note: only count the C and O if we are in the He layer

            Real zone_X = model_hse(i, model::ispec+ihe4);
            if (problem_rp::X_N14 > 0.0) {
                zone_X += model_hse(i, model::ispec+in14);
            }

            if (problem_rp::X_C12 > 0.0 && zone_X > 10.0 * network_rp::small_x) {
                zone_X += model_hse(i, model::ispec+ic12);
            }

            if (problem_rp::X_O16 > 0.0 && zone_X > 10.0 * network_rp::small_x) {
                zone_X += model_hse(i, model::ispec+io16);
            }

            Real core_X{0.0};
            if (problem_rp::X_C12 == 0.0 && problem_rp::X_O16 == 0.0) {
                core_X += model_hse(i, model::ispec+ic12) + model_hse(i, model::ispec+io16);
            } else if (model_hse(i, model::ispec+ihe4) <= network_rp::small_x) {
                core_X += model_hse(i, model::ispec+ic12) + model_hse(i, model::ispec+io16);
            }

            mass_he += vol * model_hse(i, model::idens) * zone_X;
            mass_wd += vol * model_hse(i, model::idens) * core_X;
        }

        if (mass_wd == 0.0) {
            std::string err_file = "zero_mass";
            write_model(err_file, xzn_hse, model_hse);
            amrex::Error("zero mass");
        }

        if (rho_c_old < 0.0_rt) {
            // not enough iterations yet -- store the old central
            // density and mass and pick a new value

            rho_c_old = rho_c;
            mass_wd_old = mass_wd;

            rho_he_old = rho_he;
            mass_he_old = mass_he;

            rho_c = 0.5_rt * rho_c_old;
            rho_he = 0.5_rt * rho_he_old;

        } else {
            // have we converged

            if (std::abs(mass_wd - problem_rp::M_tot) / problem_rp::M_tot < TOL_WD_MASS &&
                std::abs(mass_he - problem_rp::M_He) / problem_rp::M_He < TOL_HE_MASS) {
                mass_converged = true;
                break;
            }

            // do a secant iteration:

            // M_tot = M(rho_c) + dM/drho |_rho_c x drho + ...

            Real drho_c = (problem_rp::M_tot - mass_wd) /
                ((mass_wd  - mass_wd_old) / (rho_c - rho_c_old));

            rho_c_old = rho_c;
            mass_wd_old = mass_wd;

            rho_c = amrex::min(1.1_rt * rho_c_old,
                               amrex::max((rho_c + drho_c), 0.9_rt * rho_c_old));


            Real drho_he = (problem_rp::M_He - mass_he) /
                ((mass_he  - mass_he_old) / (rho_he - rho_he_old));

            rho_he_old = rho_he;
            mass_he_old = mass_he;

            rho_he = amrex::min(1.5_rt * rho_he_old,
                                amrex::max((rho_he + drho_he), 0.75_rt * rho_he_old));

        }

    } // end mass constraint loop

    if (! mass_converged) {
        amrex::Error("ERROR: mass did not converge");
    }

    std::cout << "converged" << std::endl;
    std::cout << "central density of WD: " << rho_c << std::endl;
    std::cout << "density at base of He layer: " << rho_he << std::endl;


    std::cout << "final masses: " << std::endl;
    std::cout << " mass WD: " << mass_wd / C::M_solar << std::endl;
    std::cout << " mass He: " << mass_he / C::M_solar << std::endl;
    std::cout << ihe_layer << std::endl;

    // store the model

    Real dCoord = xzn_hse(1) - xzn_hse(0);
    std::string dxstr = "dx" + num_to_unitstring(dCoord);

    std::string delta_str = num_to_unitstring(problem_rp::delta);

    std::stringstream ss;
    ss << std::setprecision(2) << std::fixed << mass_wd / C::M_solar;
    std::string mass_wd_str = ss.str();

    std::stringstream ss2;
    ss2 << std::setprecision(3) << std::fixed << mass_he / C::M_solar;
    std::string mass_he_str = ss2.str();

    std::stringstream ss3;
    ss3 << problem_rp::temp_base;
    std::string temp_base_str = ss3.str();

    std::string outfile_base = problem_rp::prefix + ".M_WD-" + mass_wd_str + ".M_He-" + mass_he_str  + ".delta" + delta_str + ".temp_base-" + temp_base_str;


    std::string outfile = outfile_base + ".hse";

    if (problem_rp::mixed_co_wd) {
        outfile += ".CO";
    } else {
        outfile += ".C";
    }

    if (problem_rp::X_C12 > 0.0_rt) {
        outfile += ".C12";
    }

    if (problem_rp::X_N14 > 0.0_rt) {
        outfile += ".N14";
    }

    if (problem_rp::X_O16 > 0.0_rt) {
        outfile += ".O16";
    }

    write_model(outfile, xzn_hse, model_hse);

    // extra info

    outfile = outfile_base + ".extras";

    if (problem_rp::mixed_co_wd) {
        outfile += ".CO";
    } else {
        outfile += ".C";
    }

    if (problem_rp::X_C12 > 0.0_rt) {
        outfile += ".C12";
    }
    if (problem_rp::X_N14 > 0.0_rt) {
        outfile += ".N14";
    }
    if (problem_rp::X_O16 > 0.0_rt) {
        outfile += ".O16";
    }

    outfile += "." + dxstr;

    std::ofstream ef;
    ef.open(outfile);

    ef << "# npts = " << problem_rp::nx << std::endl;
    ef << "# cs" << std::endl;
    ef << "# entropy" << std::endl;

    for (int i = 0; i < problem_rp::nx; ++i) { 
        ef << std::setprecision(12) << std::setw(20) << xzn_hse(i) << " ";
        ef << std::setprecision(12) << std::setw(20) << cs_hse(i) << " ";
        ef << std::setprecision(12) << std::setw(20) << s_hse(i) << std::endl;
    }

    ef.close();

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

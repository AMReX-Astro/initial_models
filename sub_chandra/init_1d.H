#ifndef INIT_1D_H
#define INIT_1D_H

#include <AMReX_Array.H>


#include <sstream>

#include <extern_parameters.H>
#include <fundamental_constants.H>

using namespace amrex;



AMREX_INLINE std::string
num_to_unitstring(const Real value) {

    std::stringstream ss;
    std::string result;

    if (value > 1.e5_rt) {

        // work in km
        ss << std::setprecision(3) << value / 1.e5_rt;
        result = ss.str() + "km";
    } else {

        // work in cm

        ss << std::setprecision(3) << value;
        result = ss.str() + "cm";
    }

    return result;
}

// generate an initial model for an arbitrary-mass, isothermal C WD
// with an isentropic He envelope on the surface.

AMREX_INLINE void init_1d() {

    // define convenient indices for the scalars

    const int nvar = 3 + NumSpec;
    const int idens = 0;
    const int itemp = 1;
    const int ipres = 2;
    const int ispec = 3;

    // TOL_HSE is the tolerance used when iterating over a zone to
    // force it into HSE by adjusting the current density (and
    // possibly temperature).  TOL_HSE should be very small (~
    // 1.e-10).

    const Real TOL_HSE = 1.e-10_rt;

    // TOL_WD_MASS is tolerance used for getting the total WD mass
    // equal to M_tot (defined below).  It can be reasonably small,
    // since there will always be a central density value that can
    // give the desired WD mass on the grid we use

    const Real TOL_WD_MASS = 1.e-4_rt;

    // TOL_HE_MASS is the tolerance used for getting the mass of the
    // He envelope.  This should not be too small, since the values of
    // the He envelope mass we can achieve will depend on our grid
    // spacing.

    const Real TOL_HE_MASS = 2.e-2;

    const int MAX_ITER = 250;

    const Real smallx = 1.e-10_rt;



    // convert the envelope and WD mass into CGS

    M_tot *= C::M_solar;
    M_He *= C::M_solar;

    // get the species indices

    int ihe4 = network_spec_index("helium-4");
    int ic12 = network_spec_index("carbon-12");
    int io16 = network_spec_index("oxygen-16");

    int in14 = network_spec_index("nitrogen-14");

    if (ihe4 < 0 || ic12 < 0 || io16 < 0) {
        amrex::Error("ERROR: species not defined");
    }

    if (in14 < 0 && X_N14 > 0.0_rt) {
        amrex::Error("ERROR: N14 not defined");
    }

    Real xn_core[NumSpec] = {smallx};
    Real xn_he[NumSpec] = {smallx};;

    if (mixed_co_wd) {
        xn_core[ic12] = 0.5_rt - 0.5_rt * (NumSpec - 1) * smallx;
        xn_core[io16] = 0.5_rt - 0.5_rt * (NumSpec - 1) * smallx;
    } else {
        xn_core[ic12] = 1.0_rt - (NumSpec - 1) * smallx;
    }

    if (X_N14 > 0.0_rt) {
        xn_he[in14] = X_N14;
        xn_he[ihe4] = 1.0_rt - (NumSpec - 2) * smallx - X_N14;
    } else {
        xn_he[ihe4] = 1.0_rt - (NumSpec - 1) * smallx;
    }

    // Create a 1-d uniform grid that is identical to the mesh that we
    // are mapping onto, and then we want to force it into HSE on that
    // mesh.

    Array1D<Real, 0, NPTS_MODEL-1> xzn_hse;
    Array1D<Real, 0, NPTS_MODEL-1> xznl;
    Array1D<Real, 0, NPTS_MODEL-1> xznr;

    Array2D<Real, 0, NPTS_MODEL-1, 0, nvar-1> model_hse;

    Array1D<Real, 0, NPTS_MODEL-1> M_enclosed;
    Array1D<Real, 0, NPTS_MODEL-1> cs_hse;
    Array1D<Real, 0, NPTS_MODEL-1> s_hse;

    // compute the coordinates of the new gridded function

    Real dCoord = (xmax - xmin) / static_cast<Real>(nx);

    for (int i = 0; i < nx; ++i) {
        xznl(i) = xmin + static_cast<Real>(i) * dCoord;
        xznr(i) = xmin + (static_cast<Real>(i) + 1.0_rt) * dCoord;
        xzn_hse(i) = 0.5_rt * (xznl(i) + xznr(i));
    }

    // We don't know what WD central density will give the desired
    // total mass, so we need to iterate over central density

    // we will do a secant iteration.  rho_c_old is the 'old' guess
    // for the central density and rho_c is the current guess.  After
    // 2 loops, we can start estimating the density required to yield
    // our desired mass

    Real rho_c_old = -1.0_rt;
    Real rho_c = 1.e9_rt;  // 1.e9 is a reasonable starting WD central density

    // rho_he_old is the old guess for the density to transition to He,
    // where we will be isentropic, and rho_he is the currrent guess.

    Real rho_he_old = -1.0_rt;
    Real rho_he = 0.5_rt * rho_c;

    bool mass_converged = false;

    Real mass_wd, mass_wd_old;
    Real mass_he, mass_he_old;

    int ihe_layer;
    int ihe_entropy;
    int icutoff;

    for (int iter_mass = 0; iter_mass < MAX_ITER; ++iter_mass) {

        std::cout << "mass iter = " << iter_mass << " " << rho_c << " " << temp_core << std::endl;

        bool fluff = false;

        // we start at the center of the WD and integrate outward.
        // Initialize the central conditions.

        eos_t eos_state;
        eos_state.T = temp_core;
        eos_state.rho = rho_c;
        for (int n = 0; n < NumSpec; ++n) {
            eos_state.xn[n] = xn_core[n];
        }

        // (t, rho) -> (p, s)

        eos(eos_input_rt, eos_state);

        // make the initial guess be completely uniform

        for (int i = 0; i < nx; ++i) {
            model_hse(i,idens) = eos_state.rho;
            model_hse(i,itemp) = eos_state.T;
            model_hse(i,ipres) = eos_state.p;


            for (int n = 0; n < NumSpec; ++n) {
                model_hse(i,ispec+n) = eos_state.xn[n];
            }
        }


        // keep track of the mass enclosed below the current zone

        M_enclosed(0) = (4.0_rt / 3.0_rt) * M_PI *
            (std::pow(xznr(0), 3) - std::pow(xznl(0), 3)) * model_hse(0,idens);

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

        for (int i = 1; i < nx; ++i) {

            Real delx = xzn_hse(i) - xzn_hse(i-1);

            // as the initial guess for the density, use the previous zone

            dens_zone = model_hse(i-1, idens);

            if (dens_zone > rho_he) {
                temp_zone = temp_core;
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

                Real test = 0.5_rt * (1.0_rt + std::tanh((xzn_hse(i) - xzn_hse(ihe_layer) - 4.0_rt * delta) / delta));

                if (test < 0.999_rt) {

                    // small tanh ramp up regime

                    for (int n = 0; n < NumSpec; ++n) {
                        xn[n] = xn_core[n] + 0.5_rt * (xn_he[n] - xn_core[n]) *
                            (1.0_rt + std::tanh((xzn_hse(i) - xzn_hse(ihe_layer) - 4.0_rt * delta) / delta));
                    }

                    temp_zone = temp_core + 0.5_rt * (temp_base - temp_core)*
                        (1.0_rt + std::tanh((xzn_hse(i) - xzn_hse(ihe_layer) - 4.0_rt * delta) / delta));

                    isentropic = false;

                } else {

                    // fully isentropic

                    if (ihe_entropy == -1) {
                        ihe_entropy = i;
                        temp_zone = temp_base;
                        isentropic = false;
                    } else {
                        temp_zone = model_hse(i-1, itemp);
                        isentropic = true;
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

                        p_want = model_hse(i-1, ipres) +
                            delx * 0.5_rt * (dens_zone + model_hse(i-1, idens)) * g_zone;

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

                        if (dens_zone < low_density_cutoff) {

                            dens_zone = low_density_cutoff;
                            temp_zone = temp_fluff;
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

                        // the core is isothermal, so we just need to
                        // constrain the density and pressure to agree
                        // with the EOS and HSE

                        // We difference HSE about the interface
                        // between the current zone and the one just
                        // inside.

                        p_want = model_hse(i-1, ipres) +
                            delx * 0.5_rt * (dens_zone + model_hse(i-1, idens)) * g_zone;

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

                        if (std::abs(drho) < TOL_HSE * dens_zone) {
                            converged_hse = true;
                            break;
                        }

                        if (dens_zone < low_density_cutoff) {

                            icutoff = i;
                            dens_zone = low_density_cutoff;
                            temp_zone = temp_fluff;
                            converged_hse = true;
                            fluff = true;
                            break;
                        }
                    }

                    if (temp_zone < temp_fluff && isentropic) {
                        temp_zone = temp_fluff;
                        isentropic = false;
                    }

                }  // thermo iteration loop


                if (! converged_hse) {

                    std::cout << "Error zone " << i <<  " did not converge in init_1d" << std::endl;
                    std::cout << dens_zone << " " << temp_zone << std::endl;
                    std::cout << p_want;
                    std::cout << drho;
                    amrex::Error("Error: HSE non-convergence");
                }

            } else {
                // fluff region
                dens_zone = low_density_cutoff;
                temp_zone = temp_fluff;
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

            model_hse(i, idens) = dens_zone;
            model_hse(i, itemp) = temp_zone;
            model_hse(i, ipres) = pres_zone;

            for (int n = 0; n < NumSpec; ++n) {
                model_hse(i, ispec+n) = xn[n];
            }

            M_enclosed(i) = M_enclosed(i-1) +
                (4.0_rt / 3.0_rt) * M_PI * (xznr(i) - xznl(i)) *
                (std::pow(xznr(i), 2) + xznl(i) * xznr(i) + std::pow(xznl(i), 2)) * model_hse(i, idens);

            cs_hse(i) = eos_state.cs;
            s_hse(i) = eos_state.s;

        } // end loop over zones


        // compute the total mass of the He layer and C/O WD

        mass_he = (4.0_rt / 3.0_rt) * M_PI *
            (std::pow(xznr(0), 3) - std::pow(xznl(0), 3)) *
            model_hse(0, idens) * model_hse(0, ispec+ihe4);

        mass_wd = (4.0_rt / 3.0_rt) * M_PI *
            (std::pow(xznr(0), 3) - std::pow(xznl(0), 3)) *
            model_hse(0, idens) * (model_hse(0, ispec+ic12) + model_hse(0, ispec+io16));

        for (int i = 1; i < icutoff; ++i) {
            mass_he += (4.0_rt / 3.0_rt) * M_PI *
                (xznr(i) - xznl(i)) *
                (std::pow(xznr(i), 2) + xznl(i) * xznr(i) + std::pow(xznl(i), 2)) *
                model_hse(i,idens) * model_hse(i,ispec+ihe4);

            mass_wd += (4.0_rt / 3.0_rt) * M_PI *
                (xznr(i) - xznl(i)) *
                (std::pow(xznr(i), 2) + xznl(i) * xznr(i) + std::pow(xznl(i), 2)) *
                model_hse(i,idens) * (model_hse(i,ispec+ic12) + model_hse(i,ispec+io16));
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

            if (std::abs(mass_wd - M_tot) / M_tot < TOL_WD_MASS &&
                std::abs(mass_he - M_He) / M_He < TOL_HE_MASS) {
                mass_converged = true;
                break;
            }

            // do a secant iteration:

            // M_tot = M(rho_c) + dM/drho |_rho_c x drho + ...

            Real drho_c = (M_tot - mass_wd) /
                ((mass_wd  - mass_wd_old) / (rho_c - rho_c_old));

            rho_c_old = rho_c;
            mass_wd_old = mass_wd;

            rho_c = amrex::min(1.1_rt * rho_c_old,
                               amrex::max((rho_c + drho_c), 0.9_rt * rho_c_old));


            Real drho_he = (M_He - mass_he) /
                ((mass_he  - mass_he_old) / (rho_he - rho_he_old));

            rho_he_old = rho_he;
            mass_he_old = mass_he;

            rho_he = amrex::min(1.1_rt * rho_he_old,
                                amrex::max((rho_he + drho_he), 0.9_rt * rho_he_old));

            std::cout << "current mass = " << mass_wd / C::M_solar << " " << mass_he / C::M_solar << std::endl;
        }

    } // end mass constraint loop

    if (! mass_converged) {
        amrex::Error("ERROR: mass did not converge");
    }

    std::cout << "final masses: " << std::endl;
    std::cout << " mass WD: " << mass_wd / C::M_solar << std::endl;
    std::cout << " mass He: " << mass_he / C::M_solar << std::endl;
    std::cout << ihe_layer << std::endl;

    // store the model

    std::string dxstr = num_to_unitstring(dCoord);

    std::stringstream ss;
    ss << std::setprecision(2) << mass_wd / C::M_solar;
    std::string mass_wd_str = ss.str();

    ss.clear();
    ss << std::setprecision(3) << mass_he / C::M_solar;
    std::string mass_he_str = ss.str();

    std::string outfile = "sub_chandra.M_WD-" + mass_wd_str + ".M_He-" + mass_he_str + ".hse";

    if (mixed_co_wd) {
        outfile += ".CO";
    } else {
        outfile += ".C";
    }

    if (X_N14 > 0.0_rt) {
        outfile += ".N14";
    }

    outfile += "." + dxstr;

    std::ofstream of;
    of.open(outfile);

    of << "# npts = " << nx << std::endl;
    of << "# num of variables = " << nvar << std::endl;
    of << "# density" << std::endl;
    of << "# temperature" << std::endl;
    of << "# pressure" << std::endl;

    for (int n = 0; n < NumSpec; ++n) {
        of << "# " << spec_names_cxx[n] << std::endl;
    }

    for (int i = 0; i < nx; ++i) {
        of << std::setprecision(12) << std::setw(20) << xzn_hse(i) << " ";
        of << std::setprecision(12) << std::setw(20) << model_hse(i, idens) << " ";
        of << std::setprecision(12) << std::setw(20) << model_hse(i, itemp) << " ";
        of << std::setprecision(12) << std::setw(20) << model_hse(i, ipres) << " ";
        for (int n = 0; n < NumSpec; ++n) {
            of << std::setprecision(12) << std::setw(20) << model_hse(i, ispec+n) << " ";
        }
        of << std::endl;
    }

    of.close();

#if 0
  ! extra info
  if (mixed_co_wd) then
     outfile = "sub_chandra.M_WD-" // trim(adjustl(mass_wd_str)) // &
          ".M_He-" // trim(adjustl(mass_he_str)) // &
          ".extras.CO." // trim(adjustl(dxstr))
  else
     outfile = "sub_chandra.M_WD-" // trim(adjustl(mass_wd_str)) // &
          ".M_He-" // trim(adjustl(mass_he_str)) // &
          ".extras.C." // trim(adjustl(dxstr))
  endif

  open (unit=51, file=outfile, status="unknown")

  write (51,1001) "# npts = ", nx
  write (51,1002) "# cs"
  write (51,1002) "# entropy"

  do i = 1, nx
     write (51,1000) xzn_hse(i), cs_hse(i), s_hse(i)
  enddo

  close (51)

  ! compute the maximum HSE error
  max_hse_error = -1.e30

  do i = 2, nx-1
     g_zone = -Gconst*M_enclosed(i-1)/xznr(i-1)**2
     dpdr = (model_hse(i,ipres) - model_hse(i-1,ipres))/delx
     rhog = HALF*(model_hse(i,idens) + model_hse(i-1,idens))*g_zone

     print *, xzn_hse(i), g_zone*xzn_hse(i)

     if (dpdr /= ZERO .and. model_hse(i+1,idens) > low_density_cutoff) then
        max_hse_error = max(max_hse_error, abs(dpdr - rhog)/abs(dpdr))
     endif

  enddo

  print *, 'maximum HSE error = ', max_hse_error
  print *, ' '
#endif
}


#endif

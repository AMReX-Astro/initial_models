#ifndef INIT_1D_H
#define INIT_1D_H

#include <AMReX_Array.H>


#include <sstream>

#include <extern_parameters.H>
#include <fundamental_constants.H>

#include <coord_info.H>
#include <model_util.H>
#include <read_model.H>
#include <interpolate.H>

using namespace amrex;

// Take an initial model from a Lagrangian code and put it onto
// a uniform grid and make sure that it is happy with the EOS in
// our code.  The output is a .hse file that can be read directly
// by Maestro.
//
// The model is placed into HSE by the following differencing:
//
//   (1/dr) [ <P>_i - <P>_{i-1} ] = (1/2) [ <rho>_i + <rho>_{i-1} ] g
//
// We do a hybrid model: isentropic in the interior and Kepler's temperature
// structure outside of that.

using namespace amrex;

constexpr Real TOL_HSE = 1.e-10_rt;

constexpr int MAX_ITER = 250;


AMREX_INLINE void init_1d() {

    // read in the initial model

    initial_model_t kepler_model;
    read_file(problem_rp::model_file, kepler_model);

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

    Array2D<Real, 0, NPTS_MODEL-1, 0, model::nvar-1> model_kepler_hse;
    Array2D<Real, 0, NPTS_MODEL-1, 0, model::nvar-1> model_isentropic_hse;
    Array2D<Real, 0, NPTS_MODEL-1, 0, model::nvar-1> model_hybrid_hse;

    Array1D<Real, 0, NPTS_MODEL-1> M_enclosed;
    Array1D<Real, 0, NPTS_MODEL-1> entropy_want;

    // these are only needed for an irregular grid

    Array1D<Real, 0, NPTS_MODEL-1> delrl;
    Array1D<Real, 0, NPTS_MODEL-1> delrr;

    // compute the coordinates of the new gridded function

    fill_coord_arrays_irreg(nr, xzn_hse, xznl, xznr, delrl, delrr);

    // put the data on the initial grid

    eos_t eos_state;

    for (int i = 0; i < nr; ++i) {
        for (int n = 0; n < model::nvar; ++n) {
            model_kepler_hse(i,n) = interpolate(xzn_hse(i), n,  kepler_model, false);
        }

        // make sure the mass fractions sum to 1

        Real sum_x = 0.0_rt;
        for (int n = 0; n < NumSpec; ++n) {
            model_kepler_hse(i, model::ispec+n) = std::max(model_kepler_hse(i, model::ispec+n), problem_rp::smallx);
            sum_x += model_kepler_hse(i, model::ispec+n);
        }

        for (int n = 0; n < NumSpec; ++n) {
            model_kepler_hse(i, model::ispec+n) /= sum_x;
        }


        // fix the thermodynamics

        eos_state.rho = model_kepler_hse(i, model::idens);
        eos_state.T = model_kepler_hse(i, model::itemp);

        for (int n = 0; n < NumSpec; ++n) {
            eos_state.xn[n] = model_kepler_hse(i, model::ispec+n);
        }

        eos(eos_input_rt, eos_state);

        model_kepler_hse(i, model::ipres) = eos_state.p;
    }

    write_model("uniform", xzn_hse, model_kepler_hse);

    // iterate to find the central density

    Real dens_zone;
    Real temp_zone;
    Real pres_zone;
    Real xn[NumSpec];

    // because the Kepler model likely begins at a larger radius than
    // our first HSE model zone, simple interpolation will not do a
    // good job.  We want to integrate in from the zone that best
    // matches the first Kepler model zone, assuming HSE and constant
    // entropy.

    // find the zone in the uniformly gridded model that corresponds to the
    // first zone of the original model

    int ibegin{-1};

    for (int i = 0; i < nr; ++i) {
        if (xzn_hse(i) >= kepler_model.r(0)) {
            ibegin = i;
            break;;
        }
    }

    // store the central density.  We will iterate until the central density
    // converges

    Real central_density = model_kepler_hse(0, model::idens);

    std::cout << "interpolated central density = " << central_density << std::endl;

    bool converged_central_density{false};

    for (int iter_dens = 0; iter_dens < MAX_ITER; ++iter_dens) {

        // compute the enclosed mass

        Real dx = xzn_hse(1) - xzn_hse(0);
        M_enclosed(0) = (4.0_rt/3.0_rt) * M_PI * std::pow(dx, 3) * model_kepler_hse(0, model::idens);

        for (int i = 1; i <= ibegin; ++i) {
            M_enclosed(i) = M_enclosed(i-1) +
                (4.0_rt/3.0_rt) * M_PI * (xznr(i) - xznl(i)) *
                (std::pow(xznr(i), 2) + xznl(i) * xznr(i) + std::pow(xznl(i), 2)) *
                model_kepler_hse(i, model::idens);
        }

        // now start at ibegin and integrate inward

        eos_state.T = model_kepler_hse(ibegin, model::itemp);
        eos_state.rho = model_kepler_hse(ibegin, model::idens);
        for (int n = 0; n < NumSpec; ++n) {
            eos_state.xn[n] = model_kepler_hse(ibegin, model::ispec+n);
        }

       eos(eos_input_rt, eos_state);

       model_kepler_hse(ibegin, model::ipres) = eos_state.p;

       for (int i = 0; i < nr; ++i) {
           entropy_want(i) = eos_state.s;
       }

       for (int i = ibegin-1; i >= 0; --i) {

           // as the initial guess for the temperature and density, use
           // the previous zone

           dens_zone = model_kepler_hse(i+1, model::idens);
           temp_zone = model_kepler_hse(i+1, model::itemp);
           for (int n = 0; n < NumSpec; ++n) {
               xn[n] = model_kepler_hse(i, model::ispec+n);
           }

           Real delx{0};
           Real rfrac{0};
           if (problem_rp::use_irreg_grid) {
               delx = delrr(i+1) + delrl(i);
               rfrac = delrl(i+1) / delx;
           } else {
               delx = xzn_hse(i+1) - xzn_hse(i);
               rfrac = 0.5;
           }

          // compute the gravitational acceleration on the interface between zones
          // i and i+1

          Real g_zone = -C::Gconst * M_enclosed(i) / (xznr(i) * xznr(i));

          // iteration loop

          // start off the Newton loop by saying that the zone has not converged

          bool converged_hse = false;

          Real p_want;
          Real drho;
          Real dtemp;

          for (int iter = 0; iter < MAX_ITER; ++iter) {

             p_want = model_kepler_hse(i+1, model::ipres) -
                 delx * (rfrac * dens_zone + (1.0_rt - rfrac) * model_hybrid_hse(i+1, model::idens)) * g_zone;

             // (p, s) -> (T, rho)

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

             if (std::abs(drho) < TOL_HSE * dens_zone &&
                 std::abs(dtemp) < TOL_HSE * temp_zone) {
                 converged_hse = true;
                 break;
             }

          }

          if (! converged_hse) {
              std::cout << "Error zone " << i << " did not converge in init_1d" << std::endl;
              std::cout << "integrate down" << std::endl;
              std::cout << "dens_zone, temp_zone = " << dens_zone << " " << temp_zone << std::endl;
              std::cout << "p_want = " << p_want << std::endl;
              std::cout << "drho = " << drho << std::endl;
              amrex::Error("Error: HSE non-convergence");
          }

          // call the EOS one more time for this zone and then go on to the next
          // (t, rho) -> (p, s)

          eos_state.T = temp_zone;
          eos_state.rho = dens_zone;
          for (int n = 0; n < NumSpec; ++n) {
              eos_state.xn[n] = xn[n];
          }

          eos(eos_input_rt, eos_state);

          pres_zone = eos_state.p;

          // update the thermodynamics in this zone
          model_kepler_hse(i, model::idens) = dens_zone;
          model_kepler_hse(i, model::itemp) = temp_zone;
          model_kepler_hse(i, model::ipres) = pres_zone;
          model_kepler_hse(i, model::ientr) = eos_state.s;

       }

       if (std::abs(model_kepler_hse(0, model::idens) - central_density) < TOL_HSE * central_density) {
           converged_central_density = true;
           break;
       }

       central_density = model_kepler_hse(0, model::idens);

    }

    if (! converged_central_density) {
        amrex::Error("Error: non-convergence of central density");
    }

    std::cout << "converged central density = " << model_hybrid_hse(0, model::idens) << std::endl << std::endl;


    // now integrate the entire model

    std::cout << "creating HSE model..." << std::endl;

    bool fluff{false};
    int index_hse_fluff{-1};

    // keep track of the mass enclosed below the current zone

    M_enclosed(0) = (4.0_rt / 3.0_rt) * M_PI *
        (std::pow(xznr(0), 3) - std::pow(xznl(0), 3)) * model_kepler_hse(0, model::idens);


    // HSE solve for the full star
    // here we take the temperature from the initial model

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

        // as the initial guess for the density, use
        // the previous zone

        dens_zone = model_kepler_hse(i-1, model::idens);

        // get the composition for this zone

        for (int n = 0; n < NumSpec; ++n) {
            xn[n] = model_kepler_hse(i, model::ispec+n);
        }

        // compute the gravitational acceleration on the interface
        // between zones i-1 and i

        Real g_zone = -C::Gconst * M_enclosed(i-1) / (xznl(i) * xznl(i));

        // iteration loop

        // start off the Newton loop by saying that the zone has not converged
        bool converged_hse{false};

        if (! fluff) {

            Real p_want;
            Real drho;

            for (int iter = 0; iter < MAX_ITER; ++iter) {

                p_want = model_kepler_hse(i-1, model::ipres) +
                    delx * ((1.0_rt - rfrac) * dens_zone + rfrac * model_kepler_hse(i-1, model::idens)) * g_zone;

                // we take the temperature from the model

                temp_zone = model_kepler_hse(i, model::itemp);


                if (model_kepler_hse(i-1, model::idens) < problem_rp::temp_fluff_cutoff) {
                    temp_zone = problem_rp::temp_fluff;
                }

                // we need to find the density given this pressure and our model's temperature

                eos_state.T = temp_zone;
                eos_state.rho = dens_zone;  // initial guess
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

                if (dens_zone < problem_rp::low_density_cutoff) {
                    dens_zone = problem_rp::low_density_cutoff;
                    temp_zone = problem_rp::temp_fluff;
                    converged_hse = true;
                    fluff = true;
                    index_hse_fluff = i;
                    break;
                }

            }

            if (! converged_hse) {
                std::cout << "Error zone " << i <<  " did not converge in init_1d" << std::endl;
                std::cout << dens_zone << " " << temp_zone << std::endl;
                std::cout << p_want;
                std::cout << drho;
                amrex::Error("Error: HSE non-convergence");
            }

            if (temp_zone < problem_rp::temp_fluff) {
                temp_zone = problem_rp::temp_fluff;
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

        model_kepler_hse(i, model::idens) = dens_zone;
        model_kepler_hse(i, model::itemp) = temp_zone;
        model_kepler_hse(i, model::ipres) = pres_zone;
        model_kepler_hse(i, model::ientr) = eos_state.s;

        M_enclosed(i) = M_enclosed(i-1) +
            (4.0_rt / 3.0_rt) * M_PI * (xznr(i) - xznl(i)) *
            (std::pow(xznr(i), 2) + xznl(i) * xznr(i) + std::pow(xznl(i), 2)) * model_kepler_hse(i, model::idens);

    } // end loop over zones

    std::cout << "mass = " << M_enclosed(nr-1) / C::M_solar << std::endl;;

    write_model("hse", xzn_hse, model_kepler_hse);


    // compute the alternate model using the same central density and
    // temperature, but assuming that we are isentropic (and in HSE).

    std::cout << "creating isentropic model..." << std::endl;

    fluff = false;
    bool isentropic{true};

    // start by using the Kepler model as the initial guess

    for (int i = 0; i < nr; ++i) {
        for (int n = 0; n < model::nvar; ++n) {
            model_isentropic_hse(i, n) = model_kepler_hse(i, n);
        }
    }

    // we'll try to have every zone match the entropy of the initial zone

    for (int i = 0; i < nr; ++i) {
        entropy_want(i) = model_isentropic_hse(0, model::ientr);
    }

    // keep track of the mass enclosed below the current zone

    M_enclosed(0) = (4.0_rt / 3.0_rt) * M_PI *
        (std::pow(xznr(0), 3) - std::pow(xznl(0), 3)) * model_isentropic_hse(0, model::idens);


    // HSE solve for the full star
    // here we take the entropy as given

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

        // as the initial guess for the density and temperature, use
        // the previous zone

        dens_zone = model_isentropic_hse(i-1, model::idens);
        temp_zone = model_isentropic_hse(i-1, model::itemp);

        // get the composition for this zone

        for (int n = 0; n < NumSpec; ++n) {
            xn[n] = model_isentropic_hse(i, model::ispec+n);
        }

        // compute the gravitational acceleration on the interface
        // between zones i-1 and i

        Real g_zone = -C::Gconst * M_enclosed(i-1) / (xznl(i) * xznl(i));

        // iteration loop

        // start off the Newton loop by saying that the zone has not converged
        bool converged_hse{false};

        if (! fluff) {

            Real p_want;
            Real drho;
            Real dtemp;

            for (int iter = 0; iter < MAX_ITER; ++iter) {

                p_want = model_isentropic_hse(i-1, model::ipres) +
                    delx * ((1.0_rt - rfrac) * dens_zone + rfrac * model_isentropic_hse(i-1, model::idens)) * g_zone;

                if (isentropic) {

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

                    if (std::abs(drho) < TOL_HSE * dens_zone &&
                        std::abs(dtemp) < TOL_HSE * temp_zone) {
                        converged_hse = true;
                        break;
                    }

                } else {

                    eos_state.T = problem_rp::temp_fluff;
                    eos_state.rho = dens_zone;  // initial guess
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

            }

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

        model_isentropic_hse(i, model::idens) = dens_zone;
        model_isentropic_hse(i, model::itemp) = temp_zone;
        model_isentropic_hse(i, model::ipres) = pres_zone;
        model_isentropic_hse(i, model::ientr) = eos_state.s;

        M_enclosed(i) = M_enclosed(i-1) +
            (4.0_rt / 3.0_rt) * M_PI * (xznr(i) - xznl(i)) *
            (std::pow(xznr(i), 2) + xznl(i) * xznr(i) + std::pow(xznl(i), 2)) * model_isentropic_hse(i, model::idens);

    } // end loop over zones

    std::cout << "mass = " << M_enclosed(nr-1) / C::M_solar << std::endl;;

    write_model("isentropic", xzn_hse, model_isentropic_hse);


    // compute a hybrid model -- isentropic in the interior, Kepler's
    // temperature structure outside.

    std::cout << "creating hybrid model..." << std::endl;

    Real max_temp = std::numeric_limits<Real>::lowest();
    for (int i = 0; i < nr; ++i) {
        max_temp = std::max(max_temp, model_kepler_hse(i, model::itemp));
    }

    int i_isentropic{-1};

    // copy the isentropic model as our initial guess

    for (int i = 0; i < nr; ++i) {
        for (int n = 0; n < model::nvar; ++n) {
            model_hybrid_hse(i, n) = model_isentropic_hse(i, n);
        }

        // set the temperature -- basically once the Kepler temperature
        // really flattens out, we want to stop being isentropic and use
        // the original Kepler model temperature

        if (model_kepler_hse(i, model::itemp) > model_isentropic_hse(i, model::itemp)) {

            // there will be a few regions in the very, very center
            // where the kepler temperature may be slightly higher
            // than the isentropic, but we are still not done with
            // isentropic.  i_isentropic is an index that keeps track
            // of when we switch to the original kepler model "for
            // real".  This is used for diagnostics.  We require the
            // temperature to have dropped by 10% from the central
            // value at least...

            if (i_isentropic == -1 &&
                model_isentropic_hse(i, model::itemp) < 0.9_rt * max_temp) {
                i_isentropic = i;
                std::cout << "i_isentropic = " << i_isentropic << std::endl;
            }

            model_hybrid_hse(i, model::itemp) = model_kepler_hse(i, model::itemp);
        }

    }

    // the outer part of the star will be using the original kepler
    // temperature structure.  Because the hybrid model might hit the
    // fluff region earlier or later than the kepler model, reset the
    // temperatures in the fluff region to the last valid kepler zone.

    for (int i = index_hse_fluff; i < nr; ++i) {
        model_hybrid_hse(i, model::itemp) = model_kepler_hse(index_hse_fluff-1, model::itemp);
    }


    fluff = false;

    // keep track of the mass enclosed below the current zone

    M_enclosed(0) = (4.0_rt / 3.0_rt) * M_PI *
        (std::pow(xznr(0), 3) - std::pow(xznl(0), 3)) * model_hybrid_hse(0, model::idens);


    // HSE solve for the full star

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

        // as the initial guess for the density use the previous zone

        dens_zone = model_hybrid_hse(i-1, model::idens);
        temp_zone = model_hybrid_hse(i-1, model::itemp);

        // get the composition for this zone

        for (int n = 0; n < NumSpec; ++n) {
            xn[n] = model_hybrid_hse(i, model::ispec+n);
        }

        // compute the gravitational acceleration on the interface
        // between zones i-1 and i

        Real g_zone = -C::Gconst * M_enclosed(i-1) / (xznl(i) * xznl(i));

        // iteration loop

        // start off the Newton loop by saying that the zone has not converged
        bool converged_hse{false};

        if (! fluff) {

            Real p_want;
            Real drho;

            for (int iter = 0; iter < MAX_ITER; ++iter) {

                p_want = model_hybrid_hse(i-1, model::ipres) +
                    delx * ((1.0_rt - rfrac) * dens_zone + rfrac * model_hybrid_hse(i-1, model::idens)) * g_zone;

                temp_zone = model_hybrid_hse(i, model::itemp);

                eos_state.T = temp_zone;
                eos_state.rho = dens_zone;  // initial guess
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

                if (dens_zone < problem_rp::low_density_cutoff) {
                    dens_zone = problem_rp::low_density_cutoff;
                    temp_zone = problem_rp::temp_fluff;
                    converged_hse = true;
                    fluff = true;
                    break;
                }

            }

            if (! converged_hse) {
                std::cout << "Error zone " << i <<  " did not converge in init_1d" << std::endl;
                std::cout << dens_zone << " " << temp_zone << std::endl;
                std::cout << p_want;
                std::cout << drho;
                amrex::Error("Error: HSE non-convergence");
            }

            if (temp_zone < problem_rp::temp_fluff) {
                temp_zone = problem_rp::temp_fluff;
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

        model_hybrid_hse(i, model::idens) = dens_zone;
        model_hybrid_hse(i, model::itemp) = temp_zone;
        model_hybrid_hse(i, model::ipres) = pres_zone;
        model_hybrid_hse(i, model::ientr) = eos_state.s;

        M_enclosed(i) = M_enclosed(i-1) +
            (4.0_rt / 3.0_rt) * M_PI * (xznr(i) - xznl(i)) *
            (std::pow(xznr(i), 2) + xznl(i) * xznr(i) + std::pow(xznl(i), 2)) * model_hybrid_hse(i, model::idens);

    } // end loop over zones

    std::cout << "mass = " << M_enclosed(nr-1) / C::M_solar << std::endl;;

    write_model("hybrid", xzn_hse, model_hybrid_hse);


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

        Real dpdr = (model_hybrid_hse(i, model::ipres) - model_hybrid_hse(i-1, model::ipres)) / delx;
        Real rhog = ((1.0 - rfrac) * model_hybrid_hse(i, model::idens) + rfrac * model_hybrid_hse(i-1, model::idens)) * g_zone;

        if (dpdr != 0.0_rt && model_hybrid_hse(i+1, model::idens) > problem_rp::low_density_cutoff) {
            max_hse_error = amrex::max(max_hse_error, std::abs(dpdr - rhog) / std::abs(dpdr));
        }
    }

    std::cout << "maximum HSE error = " << max_hse_error << std::endl;
}
#endif


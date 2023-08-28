#ifndef INIT_1D_H
#define INIT_1D_H

#include <AMReX_Array.H>


#include <sstream>

#include <extern_parameters.H>
#include <fundamental_constants.H>

#include <coord_info.H>
#include <model_util.H>
#include <read_model.H>

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

    int nr{-1};

    if (problem_rp::use_irreg_grid) {
        nr = (3 * (problem_rp::nx/2 - 0.5_rt) * (problem_rp::nx/2 - 0.5_rt) - 0.75_rt) / 2.0_rt;
    } else {
        nr = problem_rp::nx;
    }

    if (nr > NPTS_MODEL) {
        amrex::Error("too many zones requested -- increase NPTS_MODEL");
    }

    Array1D<Real, 0, NPTS_MODEL-1> xzn_hse;
    Array1D<Real, 0, NPTS_MODEL-1> xznl;
    Array1D<Real, 0, NPTS_MODEL-1> xznr;

    Array2D<Real, 0, NPTS_MODEL-1, 0, nvar-1> model_kepler_hse;

    Array1D<Real, 0, NPTS_MODEL-1> M_enclosed;
    Array1D<Real, 0, NPTS_MODEL-1> entropy_want;

    // these are only needed for an irregular grid

    Array1D<Real, 0, NPTS_MODEL-1> delrl;
    Array1D<Real, 0, NPTS_MODEL-1> delrr;

    // compute the coordinates of the new gridded function

    // note: this uses the uniform grid dx, regardless of whether we are doing irregular
    // or uniform gridding

    Real dCoord = (problem_rp::xmax - problem_rp::xmin) / static_cast<Real>(problem_rp::nx);

    if (problem_rp::use_irreg_grid) {

        for (int i = 0; i < nr; ++i) {
            if (i == 0) {
                // set the first edge node to xmin
                xznl(i) = problem_rp::xmin;
            } else {
                xznl(i) = problem_rp::xmin + std::sqrt(0.75_rt + 2.0_rt * (static_cast<Real>(i) - 0.5_rt)) * dCoord;
            }

            xznr(i) = problem_rp::xmin + std::sqrt(0.75_rt + 2.0_rt * (static_cast<Real>(i) + 0.5_rt)) * dCoord;
            xzn_hse(i) = problem_rp::xmin + std::sqrt(0.75_rt + 2.0_rt * static_cast<Real>(i)) * dCoord;
            delrl(i) = xzn_hse(i) - xznl(i);
            delrr(i) = xznr(i) - xzn_hse(i);
        }
    } else {

        for (int i = 0; i < nr; ++i) {
            xznl(i) = problem_rp::xmin + static_cast<Real>(i) * dCoord;
            xznr(i) = problem_rp::xmin + (static_cast<Real>(i) + 1.0_rt) * dCoord;
            xzn_hse(i) = 0.5_rt * (xznl(i) + xznr(i));
        }
    }

    // put the data on the initial grid

    for (int i = 0; i < nr; ++i) {
        for (int n = 0; n < model::nvar; ++n) {

            if (n == model::itemp) {
                model_kepler_hse(i,n) = amrex::max(problem_rp::temp_cutoff,
                                                   interpolate(xzn_hse(i), n, kepler_model, true));
            } else {
                model_kepler_hse(i,n) = interpolate(xzn_hse(i), n,  kepler_model, true);
            }

        }

        // fix the thermodynamics

        eos_state.rho = model_kepler_hse(i, model::idens);
        eos_state.T = model_kepler_hse(i, model::itemp);

        for (int n = 0; n < NumSpec; ++n){
            eos_state.xn[n] = model_kepler_hse(i, model::ispec+n);
        }

        eos(eos_input_rt, eos_state);

        model_kepler_hse(i, model::ipres) = eos_state.p;
    }

    write_model("uniform", nr, xzn_hse, model_kepler_hse);

    // iterate to find the central density

    Real dens_zone;
    Real temp_zone;
    Real pres_zone;
    Real entropy;
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

    eos_t eos_state;

    for (int iter_dens = 0; iter_dens < MAX_ITER; ++iter_dens) {

        // compute the enclosed mass

        M_enclosed(0) = (4.0_rt/3.0_rt) * M_PI * std::pow(delx, 3) * model_kepler_hse(0, model::idens);

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
               delx = delrr(i-1) + delrl(i);
               rfrac = delrr(i-1) / delx;
           } else {
               delx = xzn_hse(i) - xzn_hse(i-1);
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
                 delx * (rfrac * dens_zone + (1.0_rt - rfrac) * model_hybid_mesa_hse(i+1, model::idens)) * g_zone;

             // (p, s) -> (T, rho)

             eos_state.T = temp_zone;  // initial guess
             eos_state.rho = dens_zone;  // initial guess
             eos_state.p = p_want;
             eos_state.s = entropy_want(i);

             for (int n = 0; n < NumSpec; ++n) {
                 eos_state.xn[n] = xn[n];
             }

             eos(eos_input_ps, eos_state);

             Real dens_new = eos_state.rho;
             Real temp_new = eos_state.T;

             drho = dens_new - dens_zone;
             dens_zone = dens_new;

             dtemp = temp_new - temp_zone;
             temp_zone = temp_new;

             if (std::abs(drho) < TOL * dens_zone && std::abs(dtemp) < TOL * temp_zone) {
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

       if (std::abs(model_kepler_hse(0, model::idens) - central_density) < TOL * central_density) {
           converged_central_density = true;
           break;
       }

       central_density = model_kepler_hse(0, model::idens);

    }

    if (! converged_central_density) {
        amrex::Error("Error: non-convergence of central density");
    }

    std::cout << "converged central density = " << model_hybdrid_hse(0, model::idens) << std::endl << std::endl;


    // now integrate the entire model

    bool fluff{false};
    int index_hse_fluff{-1};

    // keep track of the mass enclosed below the current zone

    M_enclosed(0) = (4.0_rt / 3.0_rt) * M_PI *
        (std::pow(xznr(0), 3) - std::pow(xznl(0), 3)) * model_kepler_hse(0,idens);


    // HSE + entropy solve for the full star

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

        dens_zone = model_hse(i-1, model::idens);

        // get the composition for this zone

        for (int n = 0; n < NumSpec; ++n) {
            xn[n] = model_hse(i, model::ispec+n);
        }

        // we take the temperature from the model

        temp_zone = model_hse(i, model::itemp);

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

                p_want = model_hse(i-1, ipres) +
                    delx * ((1.0_rt - rfrac) * dens_zone + rfrac * model_hse(i-1, idens)) * g_zone;

                if (model_hse(i-1, idens) < problem_rp::temp_fluff_cutoff) {
                    temp_zone = problem_rp::temp_fluff;
                }

                // we need to find the density given this pressure and our model's temperature

                eos_state.T = temp_zone;
                eos_state.rho = dens_zone;  // initial guess
                eos_state.p = p_want;
                for (int n = 0; n < NumSpec; ++n) {
                    eos_state.xn[n] = xn[n];
                }

                // (t, rho) -> (p, s)
                eos(eos_input_tp, eos_state);

                Real dens_new = eos_state.rho;
                drho = dens_new - dens_zone;
                dens_zone = dens_new;

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
            (std::pow(xznr(i), 2) + xznl(i) * xznr(i) + std::pow(xznl(i), 2)) * model_hse(i, idens);

    } // end loop over zones




    std::cout << "mass = " << M_enclosed(nr-1) / C::M_solar << std::endl;;

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
        Real rhog = ((1.0 - rfrac) * model_hse(i, idens) + rfrac * model_hse(i-1, idens)) * g_zone;

        if (dpdr != 0.0_rt && model_hse(i+1, model::idens) > problem_rp::low_density_cutoff) {
            max_hse_error = amrex::max(max_hse_error, std::abs(dpdr - rhog) / std::abs(dpdr));
        }
    }

    std::cout << "maximum HSE error = " << max_hse_error << std::endl;
}
#endif


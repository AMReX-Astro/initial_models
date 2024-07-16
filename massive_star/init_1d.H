#ifndef INIT_1D_H
#define INIT_1D_H

//  Take an initial model from a Lagrangian code and put it onto
//  a uniform grid and make sure that it is happy with the EOS in
//  our code.  The output is a .hse file that can be read directly
//  by Maestro.
//
//  The model is placed into HSE by the following differencing:
//
//   (1/dr) [ <P>_i - <P>_{i-1} ] = (1/2) [ <rho>_i + <rho>_{i-1} ] g_{i-1/2}
//
//  We take the temperature structure directly from the original
//  initial model.  For the composition, we interpolate Ye and X
//  from the initial model.  If we are in an NSE region, then we
//  take Ye and call the NSE table to get X.  If we are not in NSE,
//  then we use X to compute Ye.  Since we build with USE_NSE = TRUE,
//  the EOS always requires Ye and abar through the aux data.
//  We adjust the density and pressure according to HSE using the EOS.
//

#include <AMReX_Array.H>

#include <sstream>
#include <string>

#include <extern_parameters.H>
#include <fundamental_constants.H>

#include <network.H>
#include <nse_table.H>
#include <nse_table_check.H>

#include <eos.H>
#include <eos_composition.H>

#include <coord_info.H>
#include <read_model.H>
#include <interpolate.H>

using namespace amrex;


constexpr Real TOL = 1.e-10_rt;

constexpr int MAX_ITER = 250;

constexpr Real smallx = 1.e-10_rt;

constexpr Real low_density_cutoff = 1.e-7_rt;

// temp_fluff_cutoff is the density below which we hold the temperature
// constant for the MESA model

// MAESTRO
// constexpr Real temp_fluff_cutoff = 1.e-4_rt;

// CASTRO
constexpr Real temp_fluff_cutoff = 2.e-7_rt;
constexpr Real temp_fluff = 1.e5_rt;


AMREX_INLINE void
set_aux(eos_t& eos_state) {

    bool nse_check = in_nse(eos_state);

    if (nse_check) {

        // we are in NSE, so leave ye alone, but get abar and redo xn

        Real abar_pass;
        Real bea_pass;
        Real dyedt_pass;
        Real dabardt_pass;
        Real dbeadt_pass;
        Real e_nu_pass;

	nse_table_t nse_state;
	nse_state.T = eos_state.T;
	nse_state.rho = eos_state.rho;
	nse_state.Ye = eos_state.aux[AuxZero::iye];

        nse_interp(nse_state);

	for (int n = 0; n < NumSpec; ++n) {
	    eos_state.xn[n] = nse_state.X[n];
	}

        eos_state.aux[AuxZero::iabar] = nse_state.abar;

    } else {
        set_aux_comp_from_X(eos_state);
    }
}


AMREX_INLINE void init_1d() {

    // Create a 1-d uniform grid that is identical to the mesh that we
    // are mapping onto, and then we want to force it into HSE on that
    // mesh.

    std::cout << "here " << problem_rp::nx << " " << NPTS_MODEL << std::endl;

    if (problem_rp::nx > NPTS_MODEL) {
        amrex::Error("too many zones requested -- increase NPTS_MODEL");
    }

    Array1D<Real, 0, NPTS_MODEL-1> xzn_hse;
    Array1D<Real, 0, NPTS_MODEL-1> xznl;
    Array1D<Real, 0, NPTS_MODEL-1> xznr;

    Array2D<Real, 0, NPTS_MODEL-1, 0, model::nvar-1> model_mesa_hse;

    Array1D<Real, 0, NPTS_MODEL-1> M_enclosed;

    Array1D<Real, 0, NPTS_MODEL-1> entropy_want;

    // compute the coordinates of the new gridded function

    Real delx = (problem_rp::xmax - problem_rp::xmin) / static_cast<Real>(problem_rp::nx);

    fill_coord_arrays(xzn_hse, xznl, xznr);

    // read in the MESA model

    initial_model_t initial_model;

    read_file(problem_rp::model_file, initial_model);

    std::ofstream of;
    of.open("model.orig");

    of << "# initial model as read in" << std::endl;

    for (int i = 0; i < initial_model.npts; ++i) {
        of << std::setprecision(12) << std::setw(20) << initial_model.r(i);
        for (int j = 0; j < model::nvar; ++j) {
            of << std::setprecision(12) << std::setw(20) << initial_model.state(i,j);
        }
        of << std::endl;
    }

    of.close();

    // compute the mass of the initial model -- the data is
    // non-uniform, and appears to be node-centered finite difference.
    // We'll just do a very simple trapezoid rule.

    // the first node is not at r = 0, so just use a rectangle rule to get us started

    Real M_total = (4.0_rt / 3.0_rt) * M_PI * std::pow(initial_model.r(0), 3) *
        initial_model.state(0, model::idens);

    for (int i = 1; i < initial_model.npts-1; ++i) {

        // integrate 4pi r**2 rho using trapezoid

        Real rl = initial_model.r(i-1);
        Real rr = initial_model.r(i);

        M_total += 4.0_rt * M_PI * 0.5_rt * (rr - rl) *
            (rl * rl * initial_model.state(i-1, model::idens) +
             rr * rr * initial_model.state(i, model::idens));
    }

    std::cout << "total mass as read from original model = " << M_total << std::endl;

    // put the model onto our new uniform grid

    for (int i = 0; i < problem_rp::nx; ++i) {

       for (int n = 0; n < model::nvar; ++n) {

           if (xzn_hse(i) < initial_model.r(initial_model.npts-1)) {
               model_mesa_hse(i,n) = interpolate(xzn_hse(i), n, initial_model);
           } else {
               model_mesa_hse(i,n) = initial_model.state(initial_model.npts-1, n);
           }

       }

       // make sure that the species (mass fractions) sum to 1

       Real summ = 0.0_rt;
       for (int n = 0; n < NumSpec; ++n) {
           model_mesa_hse(i, model::ispec+n) =
               amrex::max(model_mesa_hse(i, model::ispec+n), smallx);
           summ += model_mesa_hse(i, model::ispec+n);
       }

       for (int n = 0; n < NumSpec; ++n) {
           model_mesa_hse(i, model::ispec+n) /= summ;
       }

    }

    write_model("uniform", xzn_hse, model_mesa_hse, true);


    // reset the composition if we are in NSE

    eos_t eos_state;

    for (int i = 0; i < problem_rp::nx; ++i) {

        eos_state.rho = model_mesa_hse(i, model::idens);
        eos_state.T = model_mesa_hse(i, model::itemp);

        for (int n = 0; n < NumAux; ++n) {
            eos_state.aux[n] = 0.0_rt;
        }

        eos_state.aux[AuxZero::iye] = model_mesa_hse(i, model::iyef);

        for (int n = 0; n < NumSpec; ++n) {
            eos_state.xn[n] = model_mesa_hse(i, model::ispec+n);
        }

        set_aux(eos_state);

        // copy the composition variables back

        for (int n = 0; n < NumSpec; ++n) {
            model_mesa_hse(i, model::ispec+n) = eos_state.xn[n];
        }

        model_mesa_hse(i,model::iyef) = eos_state.aux[AuxZero::iye];

    }

    write_model("composition", xzn_hse, model_mesa_hse, true);


    // iterate to find the central density

    Real dens_zone;
    Real temp_zone;
    Real pres_zone;
    Real entropy;
    Real ye;
    Real xn[NumSpec];

    // the MESA model may begin at a larger radius than our first HSE
    // model zone, so simple interpolation will not do a good job.  We
    // want to integrate in from the zone that best matches the first
    // MESA model zone, assuming HSE and constant entropy.

    // find the zone in the uniformly gridded model that corresponds to the
    // first zone of the original model

    int ibegin = -1;

    for (int i = 0; i < problem_rp::nx; ++i) {
        if (xzn_hse(i) >= initial_model.r(0)) {
            ibegin = i;
            break;;
        }
    }

    std::cout << "ibegin = " << ibegin << std::endl;

    if (ibegin > 0) {

	// store the central density.  We will iterate until the central density
	// converges

	Real central_density = model_mesa_hse(0, model::idens);

	std::cout << "interpolated central density = " << central_density << std::endl;

	bool converged_central_density = false;

	for (int iter_dens = 0; iter_dens < MAX_ITER; ++iter_dens) {

	    // compute the enclosed mass

	    M_enclosed(0) = (4.0_rt/3.0_rt) * M_PI * std::pow(delx, 3) * model_mesa_hse(0, model::idens);

	    for (int i = 1; i <= ibegin; ++i) {
		M_enclosed(i) = M_enclosed(i-1) +
		    (4.0_rt/3.0_rt) * M_PI * (xznr(i) - xznl(i)) *
		    (std::pow(xznr(i), 2) + xznl(i) * xznr(i) + std::pow(xznl(i), 2)) *
		    model_mesa_hse(i, model::idens);
	    }

	    // now start at ibegin and integrate inward

	    eos_state.T = model_mesa_hse(ibegin, model::itemp);
	    eos_state.rho = model_mesa_hse(ibegin, model::idens);
	    for (int n = 0; n < NumSpec; ++n) {
		eos_state.xn[n] = model_mesa_hse(ibegin, model::ispec+n);
	    }

	   eos_state.aux[AuxZero::iye] = model_mesa_hse(ibegin, model::iyef);

	   set_aux(eos_state);

	   eos(eos_input_rt, eos_state);

	   model_mesa_hse(ibegin, model::ipres) = eos_state.p;

	   for (int i = 0; i < problem_rp::nx; ++i) {
	       entropy_want(i) = eos_state.s;
	   }

	   for (int i = ibegin-1; i >= 0; --i) {

	       // as the initial guess for the temperature and density, use
	       // the previous zone

	       dens_zone = model_mesa_hse(i+1, model::idens);
	       temp_zone = model_mesa_hse(i+1, model::itemp);
	       for (int n = 0; n < NumSpec; ++n) {
		   xn[n] = model_mesa_hse(i, model::ispec+n);
	       }

	      ye = model_mesa_hse(i, model::iyef);

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

		 p_want = model_mesa_hse(i+1, model::ipres) -
		     delx * 0.5_rt * (dens_zone + model_mesa_hse(i+1, model::idens)) * g_zone;

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

		 eos_state.aux[AuxZero::iye] = ye;

		 set_aux(eos_state);

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
		 Real dAdrho = -0.5_rt * delx * g_zone - dpd;
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

	      eos_state.aux[AuxZero::iye] = ye;

	      set_aux(eos_state);

	      eos(eos_input_rt, eos_state);

	      pres_zone = eos_state.p;

	      // update the thermodynamics in this zone
	      model_mesa_hse(i, model::idens) = dens_zone;
	      model_mesa_hse(i, model::itemp) = temp_zone;
	      model_mesa_hse(i, model::ipres) = pres_zone;
	      model_mesa_hse(i, model::ientr) = eos_state.s;

	   }

	   if (std::abs(model_mesa_hse(0, model::idens) - central_density) < TOL*central_density) {
	       converged_central_density = true;
	       break;
	   }

	   central_density = model_mesa_hse(0, model::idens);

	}

	if (! converged_central_density) {
	    amrex::Error("Error: non-convergence of central density");
	}

	std::cout << "converged central density = " << model_mesa_hse(0, model::idens) << std::endl << std::endl;

    } else {

	// we will integrate from zone 0, so make sure that is
	// thermodynamically consistent

	eos_state.T = model_mesa_hse(0, model::itemp);
	eos_state.rho = model_mesa_hse(0, model::idens);
	for (int n = 0; n < NumSpec; ++n) {
	    eos_state.xn[n] = model_mesa_hse(0, model::ispec+n);
	}

	eos_state.aux[AuxZero::iye] = model_mesa_hse(0, model::iyef);

	set_aux(eos_state);
	eos(eos_input_rt, eos_state);

	model_mesa_hse(0, model::ipres) = eos_state.p;

    }

    // compute the full HSE model using our new central density and
    // temperature, and the temperature structure as dictated by the
    // MESA model.

    std::cout << "putting MESA model into HSE on our grid..." << std::endl;

    // compute the enclosed mass

    M_enclosed(0) = (4.0_rt/3.0_rt) * M_PI * std::pow(delx, 3) * model_mesa_hse(0, model::idens);

    bool fluff = false;
    int index_hse_fluff = -1;

    for (int i = 1; i < problem_rp::nx; ++i) {

        // use previous zone as initial guess for rho

        dens_zone = model_mesa_hse(i-1, model::idens);

        // we use the model value for temperature and compositon

        temp_zone = model_mesa_hse(i, model::itemp);

        for (int n = 0; n < NumSpec; ++n) {
            xn[n] = model_mesa_hse(i, model::ispec+n);
        }

       ye = model_mesa_hse(i, model::iyef);

       // compute the gravitational acceleration on the interface between zones
       // i-1 and i

       Real g_zone = -C::Gconst * M_enclosed(i-1) / (xznr(i-1) * xznr(i-1));

       // iteration loop

       // the goal here is to find the density that is consistent with HSE

       bool converged_hse = false;

       if (! fluff) {

           Real p_want;
           Real drho;

           for (int iter = 0; iter < MAX_ITER; ++iter) {

               // HSE differencing

               p_want = model_mesa_hse(i-1, model::ipres) +
                   delx * 0.5_rt * (dens_zone + model_mesa_hse(i-1, model::idens)) * g_zone;

               temp_zone = model_mesa_hse(i, model::itemp);

               if (model_mesa_hse(i-1, model::idens) < temp_fluff_cutoff) {
                   temp_zone = temp_fluff;
               }

               // (t, rho) -> (p)

               eos_state.T = temp_zone;
               eos_state.rho = dens_zone;
               for (int n = 0; n < NumSpec; ++n) {
                   eos_state.xn[n] = xn[n];
               }

               eos_state.aux[AuxZero::iye] = ye;

	       // note: if we are in NSE, then as we adjust density in
	       // this Newton loop, the composition will change, so
	       // convergence should mean that the density and
	       // composition have converged.

               set_aux(eos_state);

               eos(eos_input_rt, eos_state);

               pres_zone = eos_state.p;

               Real dpd = eos_state.dpdr;
               drho = (p_want - pres_zone) / (dpd - 0.5_rt * delx * g_zone);

               dens_zone =
                   amrex::max(0.9_rt * dens_zone,
                              amrex::min(dens_zone + drho, 1.1_rt * dens_zone));

               if (std::abs(drho) < TOL * dens_zone) {
                   converged_hse = true;
                   break;
               }

               if (dens_zone < low_density_cutoff) {
                   dens_zone = low_density_cutoff;
                   temp_zone = temp_fluff;
                   converged_hse = true;
                   fluff = true;
                   index_hse_fluff = i;
               }

           }

           if (! converged_hse) {
               std::cout << "Error zone " << i << " did not converge in init_1d" << std::endl;
               std::cout << "integrate up" << std::endl;
               std::cout << "dens_zone, temp_zone = " << dens_zone << " " << temp_zone << std::endl;
               std::cout << "p_want = " << p_want << std::endl;
               std::cout << "drho = " << drho << std::endl;
               amrex::Error("Error: HSE non-convergence");
           }


           if (temp_zone < temp_fluff) {
               temp_zone = temp_fluff;
           }

       } else {
           dens_zone = low_density_cutoff;
           temp_zone = temp_fluff;
       }

       // call the EOS one more time for this zone and then go on to the next
       // (t, rho) -> (p)

       eos_state.T = temp_zone;
       eos_state.rho = dens_zone;
       for (int n = 0; n < NumSpec; ++n) {
           eos_state.xn[n] = xn[n];
       }

       eos_state.aux[AuxZero::iye] = ye;

       set_aux(eos_state);

       // if we were in NSE, then this updated eos_state % xn(:), so copy that over

       for (int n = 0; n < NumSpec; ++n) {
           model_mesa_hse(i, model::ispec+n) = eos_state.xn[n];
       }

       // if we were not in NSE, then this updated ye

       model_mesa_hse(i, model::iyef) = eos_state.aux[AuxZero::iye];

       std::cout << "output: " << eos_state.T << " " << eos_state.rho << " " << eos_state.aux[AuxZero::iye] << std::endl;

       eos(eos_input_rt, eos_state);

       pres_zone = eos_state.p;

       // update the thermodynamics in this zone

       model_mesa_hse(i, model::idens) = dens_zone;
       model_mesa_hse(i, model::itemp) = temp_zone;
       model_mesa_hse(i, model::ipres) = pres_zone;
       model_mesa_hse(i, model::ientr) = eos_state.s;

       M_enclosed(i) = M_enclosed(i-1) +
           (4.0_rt/3.0_rt) * M_PI * (xznr(i) - xznl(i)) *
           (std::pow(xznr(i), 2) + xznl(i) * xznr(i) + std::pow(xznl(i), 2)) * model_mesa_hse(i, model::idens);

    }


    // output

    std::string model_name = "hse";

    write_model(model_name, xzn_hse, model_mesa_hse, true);

    std::cout << "total mass = " << M_enclosed(problem_rp::nx-1) << " g; " << M_enclosed(problem_rp::nx-1) / C::M_solar << " solar masses" << std::endl;

    // compute the maximum HSE error

    Real max_hse_error = -1.e30_rt;

    for (int i = 1; i < problem_rp::nx-1; ++i) {
        Real g_zone = -C::Gconst * M_enclosed(i-1) / std::pow(xznr(i-1), 2);
        Real dpdr = (model_mesa_hse(i, model::ipres) - model_mesa_hse(i-1, model::ipres)) / delx;
        Real rhog = 0.5_rt * (model_mesa_hse(i, model::idens) + model_mesa_hse(i-1, model::idens)) * g_zone;

        if (dpdr != 0.0_rt && model_mesa_hse(i+1, model::idens) > low_density_cutoff) {
            max_hse_error = amrex::max(max_hse_error, std::abs(dpdr - rhog) / std::abs(dpdr));
        }
    }

    std::cout << "maximum HSE error = " << max_hse_error << std::endl;
}
#endif


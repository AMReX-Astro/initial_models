#ifndef INIT_1D_H
#define INIT_1D_H

#include <AMReX_Array.H>


#include <cmath>
#include <sstream>
#include <network.H>
#include <eos.H>
#include <actual_rhs.H>
#include <burn_type.H>

#include <extern_parameters.H>
#include <fundamental_constants.H>

#include <model_util.H>

using namespace amrex;


constexpr Real TOL = 1.e-11_rt;

constexpr int MAX_ITER = 1000000;
constexpr int MAX_RETRY = 50;


AMREX_INLINE void fopt_urca_23(eos_t& eos_state, Real& fopt, Real& r_ecap, Real& r_beta){
    burn_t burn_state;
    rate_t rates;

    int ine23 = network_spec_index("neon-23");
    int ina23 = network_spec_index("sodium-23");

    composition(eos_state);
    eos_to_burn(eos_state, burn_state);
    constexpr int do_T_derivatives{0};
    evaluate_rates<do_T_derivatives, rate_t>(burn_state, rates);

    r_ecap = rates.screened_rates(k_na23__ne23);
    r_beta = rates.screened_rates(k_ne23__na23);

    // calculate e-capture and beta decay rates are in equilibrium
    Real xr_ecap = burn_state.xn[ina23] * r_ecap;
    Real xr_beta = burn_state.xn[ine23] * r_beta;

    fopt = xr_ecap - xr_beta;
}

AMREX_INLINE void floor_species(eos_t& eos_state) {
    // Species mass fractions below 1E-30 should just be 0.0
    // This prevents getting abundances of ~1E-100. Note that
    // the initial model writing omits the "E" in the exponent
    // if the exponent requires 3 digits.

    Real mass_fraction_floor = 1.0e-30;

    for (int i=0; i < NumSpec; ++i) {
        if (eos_state.xn[i] < mass_fraction_floor){
            eos_state.xn[i] = 0.0_rt;
        }
    }
}

AMREX_INLINE void renormalize_species(eos_t& eos_state) {
    // Renormalize the mass fractions so they sum to 1
    Real sumx=0.0_rt;
    for (int i=0; i < NumSpec; ++i){
        sumx += eos_state.xn[i];
    }

    for (int i=0; i < NumSpec; ++i){
        eos_state.xn[i] = eos_state.xn[i]/sumx;
    }

}

AMREX_INLINE void composition_equilibrium(eos_t& eos_state){
    burn_t burn_state;
    rate_t rates;
    Real fopt;

    // get the species indices
    int ic12  = network_spec_index("carbon-12");
    int io16  = network_spec_index("oxygen-16");
    int ine23 = network_spec_index("neon-23");
    int ina23 = network_spec_index("sodium-23");

    // Initialize mass fractions given "in" values
    for (int i=0; i < NumSpec; ++i){
        eos_state.xn[i] = 0.0_rt;
    }

    eos_state.xn[ic12] = problem_rp::c12_in;
    eos_state.xn[io16] = problem_rp::o16_in;
    eos_state.xn[ine23] = 0.5_rt * problem_rp::na_ne_23;
    eos_state.xn[ina23] = 0.5_rt * problem_rp::na_ne_23;

    // Estimate the mass fractions approximating the rates as
    // independent of ye.
    composition(eos_state);
    eos_to_burn(eos_state, burn_state);
    constexpr int do_T_derivatives{0};
    evaluate_rates<do_T_derivatives, rate_t>(burn_state, rates);

    Real r_ecap = rates.screened_rates(k_na23__ne23);
    Real r_beta = rates.screened_rates(k_ne23__na23);

    eos_state.xn[ine23] = problem_rp::na_ne_23/(1.0_rt + r_beta/r_ecap);
    eos_state.xn[ina23] = problem_rp::na_ne_23 - eos_state.xn[ine23];

    //Keep the mass fraction sum X(ne23) + X(na23) = na_ne_23
    // Find the A=23 mass fractions such that A=23 Urca rates are in equilibrium
    // Do Newton iterations approximating the rates as independent of mass fraction
    fopt_urca_23(eos_state, fopt, r_ecap, r_beta);

    Real rate_equilibrium_tol = 1.0e-10;
    int max_equilibrium_iters = 10000;
    int j=1;
    while (std::abs(fopt) > rate_equilibrium_tol && j < max_equilibrium_iters) {
        Real dx = -1.0_rt * fopt/(r_ecap + r_beta);
        if (fopt > 0.0_rt) {
            eos_state.xn[ina23] = eos_state.xn[ina23] + dx;
            eos_state.xn[ine23] = problem_rp::na_ne_23 - eos_state.xn[ina23];
            fopt_urca_23(eos_state, fopt, r_ecap, r_beta);
        } else {
            eos_state.xn[ine23] = eos_state.xn[ine23] - dx;
            eos_state.xn[ina23] = problem_rp::na_ne_23 - eos_state.xn[ine23];
            fopt_urca_23(eos_state, fopt, r_ecap, r_beta);
        }
        j++;
    }

    if (j == max_equilibrium_iters) {
        amrex::Error("species iteration did not converge!");
    }
}

AMREX_INLINE void set_urca_composition(eos_t& eos_state) {
    // Construct composition profiles
    composition_equilibrium(eos_state);

    // Floor species so abundances below 1E-30 are 0.0d0
    floor_species(eos_state);

    // Renormalize species
    renormalize_species(eos_state);

}


AMREX_INLINE void get_fluff_temperature(Real& temp, Real temp_previous){
    // There are 2 kinds of fluff temperature handling:
    // 1) fluff_type = "constant"   : fluff is at the temperature temp_fluff in the inputs
    // 2) fluff_type = "continuous" : fluff is at the same temperature as the material
    //                              immediately interior to the fluff.

    if (problem_rp::fluff_type == "constant") {
        temp = problem_rp::temp_fluff;
    } else if (problem_rp::fluff_type == "continuous") {
        temp = temp_previous;
    } else {
        amrex::Error("ERROR: invalid fluff_type");
    }
}


AMREX_INLINE void iter_dens_temp(Real& dens_zone, Real& temp_zone, eos_t& eos_state_inwards,
                          Real g_zone, Real delx, Real tol, Real dtol_fac,
                          bool& isentropic, bool test_hse_convergence, bool& converged_hse, bool fluff,
                          int& i_fluff, int& i_conv, int izone){

    Real dtemp, drho, dpd;
    eos_t eos_state;

    if (isentropic){
        Real p_want = eos_state_inwards.p + delx*0.5_rt*(dens_zone + eos_state_inwards.rho)*g_zone;

        // now we have two functions to zero:
        //   A = p_want - p(rho,T)
        //   B = entropy_want - s(rho,T)
        // We use a two dimensional Taylor expansion and find the deltas
        // for both density and temperature

        eos_state.T     = temp_zone;
        eos_state.rho   = dens_zone;
        set_urca_composition(eos_state);

        // (t, rho) -> (p, s)
        eos(eos_input_rt, eos_state);

        Real entropy = eos_state.s;
        Real pres_zone = eos_state.p;

        Real dpt = eos_state.dpdT;
        dpd = eos_state.dpdr;
        Real dst = eos_state.dsdT;
        Real dsd = eos_state.dsdr;

        Real A = p_want - pres_zone;
        Real B = eos_state_inwards.s - entropy;

        Real dAdT = -dpt;
        Real dAdrho = 0.5_rt*delx*g_zone - dpd;
        Real dBdT = -dst;
        Real dBdrho = -dsd;

        dtemp = (B - (dBdrho/dAdrho)*A)/
          ((dBdrho/dAdrho)*dAdT - dBdT);

        drho = -(A + dAdT*dtemp)/dAdrho;

        dens_zone = amrex::max((1.0_rt-dtol_fac)*dens_zone,
            amrex::min(dens_zone + drho, (1.0_rt+dtol_fac)*dens_zone));

        temp_zone = amrex::max((1.0_rt-dtol_fac)*temp_zone,
            amrex::min(temp_zone + dtemp, (1.0_rt+dtol_fac)*temp_zone));

        // check if the density falls below our minimum cut-off --
        // if so, floor it

        if (dens_zone < problem_rp::low_density_cutoff) {
            i_fluff = izone;

            dens_zone = problem_rp::low_density_cutoff;

            Real temp_before_fluff = eos_state_inwards.T;

            get_fluff_temperature(temp_zone, temp_before_fluff);

            converged_hse = true;
            fluff = true;
            return;
        }

        if (dens_zone < problem_rp::dens_conv_zone && isentropic){
            i_conv = izone;
            isentropic = false;
        }

        if (test_hse_convergence){
            // Check to see if HSE is satisfied to relative tolerance TOL

            // Call EOS again using the new values of dens_zone and temp_zone
            // to get pres_zone ...
            eos_state.T     = temp_zone;
            eos_state.rho   = dens_zone;
            set_urca_composition(eos_state);

            // (t, rho) -> (p, s)
            eos(eos_input_rt, eos_state);

            pres_zone = eos_state.p;

            // Compute dpdr
            Real dpdr = (pres_zone - eos_state_inwards.p)/delx;

            // Compute rho*g
            Real rhog = 0.5_rt*(dens_zone + eos_state_inwards.rho)*g_zone;

            if (std::abs(dpdr - rhog) < tol*std::abs(dpdr)){
                converged_hse = true;
                return;
            }
        } else {
            // Check to see if drho is under the relative tolerance TOL
            if (std::abs(drho) < tol*dens_zone){
                converged_hse = true;
                return;
            }
        }
    } else {
        // do isothermal
        Real p_want = eos_state_inwards.p +
            delx*0.5_rt*(dens_zone + eos_state_inwards.rho)*g_zone;

        temp_zone = eos_state_inwards.T;


        eos_state.T     = temp_zone;
        eos_state.rho   = dens_zone;
        set_urca_composition(eos_state);

        // (t, rho) -> (p, s)
        eos(eos_input_rt, eos_state);

        Real pres_zone = eos_state.p;

        dpd = eos_state.dpdr;

        drho = (p_want - pres_zone)/(dpd - 0.5_rt*delx*g_zone);

        dens_zone = amrex::max(0.9_rt*dens_zone,
            amrex::min(dens_zone + drho, 1.1_rt*dens_zone));

        if (dens_zone < problem_rp::low_density_cutoff) {
            i_fluff = izone;

            dens_zone = problem_rp::low_density_cutoff;

            Real temp_before_fluff = eos_state_inwards.T;

            get_fluff_temperature(temp_zone, temp_before_fluff);

            converged_hse = true;
            fluff = true;
            return;
        }

        if (test_hse_convergence){
            // Check to see if HSE is satisfied to relative tolerance TOL

            // Call EOS again using the new values of dens_zone and temp_zone
            // to get pres_zone ...
            eos_state.T     = temp_zone;
            eos_state.rho   = dens_zone;
            set_urca_composition(eos_state);

            // (t, rho) -> (p, s)
            eos(eos_input_rt, eos_state);

            pres_zone = eos_state.p;

            // Compute dpdr
            Real dpdr = (pres_zone - eos_state_inwards.p)/delx;

            // Compute rho*g
            Real rhog = 0.5_rt*(dens_zone + eos_state_inwards.rho)*g_zone;

            if (std::abs(dpdr - rhog) < tol*std::abs(dpdr)){
                converged_hse = true;
                return;
            }
        } else {
            // Check to see if drho is under the relative tolerance TOL
            if (std::abs(drho) < tol*dens_zone){
                converged_hse = true;
                return;
            }
        }
    }
}


// generate an initial model

AMREX_INLINE void init_1d() {


    // i think unnecessary. int i, j, n;

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

    //const Real TOL_HSE = 1.e-10_rt;
    const Real dtol_fac = 0.000001_rt;

    if (problem_rp::nx > NPTS_MODEL) {
        amrex::Error("too many zones requested -- increase NPTS_MODEL");
    }

    Array1D<Real, 0, NPTS_MODEL-1> xzn_hse;
    Array1D<Real, 0, NPTS_MODEL-1> xznl;
    Array1D<Real, 0, NPTS_MODEL-1> xznr;

    Array2D<Real, 0, NPTS_MODEL-1, 0, nvar-1> model_hse;

    Array1D<Real, 0, NPTS_MODEL-1> M_enclosed;
    Array1D<Real, 0, NPTS_MODEL-1> entropy_want;

    // we'll get the composition indices from the network module
    int ihe4 = network_spec_index("helium-4");
    int ic12 = network_spec_index("carbon-12");
    int io16 = network_spec_index("oxygen-16");

    int ine20 = network_spec_index("neon-20");
    int ine23 = network_spec_index("neon-23");
    int ina23 = network_spec_index("sodium-23");
    int img23 = network_spec_index("magnesium-23");

    if (ihe4 < 0 || ic12 < 0 || io16 < 0 || ine20 < 0 ||
        ine23 < 0 || ina23 < 0 || img23 < 0) {
        amrex::Error("ERROR: species not defined");
    }

    // check the species mass fractions
    if (problem_rp::c12_in < 0.0_rt || problem_rp::c12_in > 1.0_rt) {
       amrex::Error("ERROR: c12_in must be between 0 and 1");

    } else if (problem_rp::o16_in < 0.0_rt || problem_rp::o16_in > 1.0_rt) {
       amrex::Error("ERROR: o16_in must be between 0 and 1");

    } else if (problem_rp::na_ne_23 < 0.0_rt || problem_rp::na_ne_23 > 1.0_rt) {
       amrex::Error("ERROR: na_ne_23 must be between 0 and 1");
    }

    //Initialize the temperature interior to fluff
    //(used to calculate fluff temperature for various values of fluff_type)
    Real temp_before_fluff = problem_rp::temp_fluff;

    //use_eos_coulomb comes in from extern_probin_module -- override
    //here if desired
    use_eos_coulomb = true;


    //-----------------------------------------------------------------------------
    // Create a 1-d uniform grid that is identical to the mesh that we are
    // mapping onto, and then we want to force it into HSE on that mesh.
    //-----------------------------------------------------------------------------

    // compute the coordinates of the new gridded function

    Real delx = (problem_rp::xmax - problem_rp::xmin) / static_cast<Real>(problem_rp::nx);

    for (int i = 0; i < problem_rp::nx; ++i) {
        xznl(i) = problem_rp::xmin + static_cast<Real>(i) * delx;
        xznr(i) = problem_rp::xmin + (static_cast<Real>(i) + 1.0_rt) * delx;
        xzn_hse(i) = 0.5_rt * (xznl(i) + xznr(i));
    }

    // we start at the center of the WD and integrate outward.
    // Initialize the central conditions.

    eos_t eos_state;
    eos_state.T = problem_rp::temp_base;
    eos_state.rho = problem_rp::dens_base;
    set_urca_composition(eos_state);

    // (t, rho) -> (p, s)
    eos(eos_input_rt, eos_state);

    // make the initial guess be completely uniform

    for (int i = 0; i < problem_rp::nx; ++i) {
        model_hse(i,idens) = eos_state.rho;
        model_hse(i,itemp) = eos_state.T;
        model_hse(i,ipres) = eos_state.p;

        for (int n = 0; n < NumSpec; ++n) {
            model_hse(i,ispec+n) = eos_state.xn[n];
        }
    }

    // print out initial uniform model
    std::ofstream init_of;
    init_of.open("initial_model.out");

    init_of << "# npts = " << problem_rp::nx << std::endl;
    init_of << "# num init_of variables = " << nvar << std::endl;
    init_of << "# density" << std::endl;
    init_of << "# temperature" << std::endl;
    init_of << "# pressure" << std::endl;

    for (int n = 0; n < NumSpec; ++n) {
        init_of << "# " << spec_names_cxx[n] << std::endl;
    }

    for (int i = 0; i < problem_rp::nx; ++i) {
        init_of << std::setprecision(12) << std::setw(20) << xzn_hse(i) << " ";
        init_of << std::setprecision(12) << std::setw(20) << model_hse(i, idens) << " ";
        init_of << std::setprecision(12) << std::setw(20) << model_hse(i, itemp) << " ";
        init_of << std::setprecision(12) << std::setw(20) << model_hse(i, ipres) << " ";
        for (int n = 0; n < NumSpec; ++n) {
            init_of << std::setprecision(12) << std::setw(20) << model_hse(i, ispec+n) << " ";
        }
        init_of << std::endl;
    }

    init_of.close();

    for (int i = 0; i < problem_rp::nx; ++i) {
           entropy_want(i) = eos_state.s;
    }

    // keep track of the mass enclosed below the current zone

    M_enclosed(0) = (4.0_rt / 3.0_rt) * M_PI *
        (std::pow(xznr(0), 3) - std::pow(xznl(0), 3)) * model_hse(0,idens);

    //-----------------------------------------------------------------------------
    // HSE + entropy solve
    //-----------------------------------------------------------------------------

    bool isentropic = true;
    bool test_hse_convergence = true;
    bool fluff = false;
    int i_conv, i_fluff;
    Real dens_zone;
    Real temp_zone;
    Real pres_zone;


    for (int i = 1; i < problem_rp::nx; ++i) {

        delx = xzn_hse(i) - xzn_hse(i-1);

        //! as the initial guess for the temperature and density, use the previous
        //! zone
        dens_zone = model_hse(i-1, idens);
        temp_zone = model_hse(i-1, itemp);
        eos_state.rho = dens_zone;
        eos_state.T = temp_zone;
        set_urca_composition(eos_state);
        eos_state.p = model_hse(i-1, ipres);
        eos_state.s = entropy_want(i);

        Real g_zone = -C::Gconst*M_enclosed(i-1)/std::pow(xznl(i), 2);

        //-----------------------------------------------------------------------
        // iteration loop
        //-----------------------------------------------------------------------

        // start off the Newton loop by saying that the zone has not converged
        bool converged_hse = false;

        Real dpd;
        if (! fluff) {


            for (int iretry = 1; iretry < MAX_RETRY; ++iretry){
                Real try_tol = TOL*std::pow(1.1_rt,(iretry-1));

                for (int iter = 0; iter < MAX_ITER; ++iter) {
                    iter_dens_temp(dens_zone, temp_zone, eos_state,
                                  g_zone, delx, try_tol, dtol_fac,
                                  isentropic, test_hse_convergence, converged_hse, fluff,
                                  i_fluff, i_conv, i);
                    if (converged_hse) {
                        break;
                    }
                }
                if (converged_hse){
                    break;
                } else if (iretry == MAX_RETRY) {
                    amrex::Error("Error: HSE non-convergence with retries!");
                }
            }

            if (temp_zone < problem_rp::temp_fluff) {
                if (problem_rp::fluff_type == "constant") {
                    get_fluff_temperature(temp_zone, temp_before_fluff);
                }
                isentropic = false;
            }
        } else {
            dens_zone = problem_rp::low_density_cutoff;
            get_fluff_temperature(temp_zone, temp_before_fluff);

        }

        // call the EOS one more time for this zone and then go on
        // to the next

        eos_state.T = temp_zone;
        eos_state.rho = dens_zone;
        set_urca_composition(eos_state);

        // (t, rho) -> (p, s)
        eos(eos_input_rt, eos_state);

        pres_zone = eos_state.p;

        dpd = eos_state.dpdr;

        // update the thermodynamics in this zone

        model_hse(i, idens) = dens_zone;
        model_hse(i, itemp) = temp_zone;
        model_hse(i, ipres) = pres_zone;

        for (int n = 0; n < NumSpec; ++n) {
            model_hse(i, ispec+n) = eos_state.xn[n];
        }


        M_enclosed(i) = M_enclosed(i-1) +
                (4.0_rt / 3.0_rt) * M_PI * (xznr(i) - xznl(i)) *
                (std::pow(xznr(i), 2) + xznl(i) * xznr(i) + std::pow(xznl(i), 2)) * model_hse(i, idens);

        if ((M_enclosed(i) > problem_rp::M_conv_zone * C::M_solar) && isentropic) {
            i_conv = i;
            isentropic = false;

        }
    }

    std::cout << "mass = " << M_enclosed(problem_rp::nx-1) / C::M_solar << std::endl;

    //print/save model
    int ipos = problem_rp::prefix.find(".dat");
    if (ipos < 0) {
        ipos = problem_rp::prefix.find(".txt");
    }

    std::string outfile = problem_rp::prefix.substr(0, ipos) + ".hse." + std::to_string(problem_rp::nx);

    //std::cout << "writing " << model_name << " model to " << outfile << std::endl;

    std::ofstream of;
    of.open(outfile);

    of << "# npts = " << problem_rp::nx << std::endl;
    of << "# num of variables = " << nvar << std::endl;
    of << "# density" << std::endl;
    of << "# temperature" << std::endl;
    of << "# pressure" << std::endl;

    for (int n = 0; n < NumSpec; ++n) {
        of << "# " << spec_names_cxx[n] << std::endl;
    }

    for (int i = 0; i < problem_rp::nx; ++i) {
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


    // compute the maximum HSE error

    Real max_hse_error = -1.e30_rt;

    for (int i = 1; i < problem_rp::nx-1; ++i) {
        Real g_zone = -C::Gconst * M_enclosed(i-1) / (xznr(i-1) * xznr(i-1));
        Real dpdr = (model_hse(i, ipres) - model_hse(i-1, ipres)) / delx;
        Real rhog = 0.5_rt * (model_hse(i, idens) + model_hse(i-1, idens)) * g_zone;

        if (dpdr != 0.0_rt && model_hse(i+1, idens) > problem_rp::low_density_cutoff) {
            max_hse_error = amrex::max(max_hse_error, std::abs(dpdr - rhog) / std::abs(dpdr));
        }
    }

    // print out the stats
    std::cout << "maximum HSE error = " << max_hse_error << std::endl;
    std::cout << std::endl;

    std::cout << "total mass (Msun) = " << M_enclosed(i_fluff) / C::M_solar << std::endl;
    std::cout << "convective zone mass (Msun) = " << M_enclosed(i_conv) / C::M_solar << std::endl;
    std::cout << "convective zone boundary density (g/cm^3) = " << model_hse(i_conv, idens) << std::endl;
    std::cout << "convective zone boundary radius (cm) = " << xzn_hse(i_conv) << std::endl;
    std::cout << "convective zone boundary index = " << i_conv << std::endl;
}
#endif

#ifndef INIT_1D_H
#define INIT_1D_H

#include <AMReX_Array.H>


#include <sstream>

#include <extern_parameters.H>
#include <fundamental_constants.H>

#include <model_util.H>

using namespace amrex;


constexpr Real TOL = 1.e-11;

constexpr int MAX_ITER = 1000000;
constexpr int MAX_RETRY = 50;

// generate an initial model for an arbitrary-mass, isothermal C WD
// with an isentropic He envelope on the surface.

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

    const Real TOL_HSE = 1.e-10_rt;
    const Real dtol_fac = 0.000001_rt;

    Array1D<Real, 0, NPTS_MODEL-1> xzn_hse;
    Array1D<Real, 0, NPTS_MODEL-1> xznl;
    Array1D<Real, 0, NPTS_MODEL-1> xznr;
    
    Array2D<Real, 0, NPTS_MODEL-1, 0, model::nvar-1> model_hse;

    Array1D<Real, 0, NPTS_MODEL-1> M_enclosed;
    Array1D<Real, 0, NPTS_MODEL-1> entropy_want;    
    Array1D<Real, 0, NPTS_MODEL-1> entropy_store;

    cont Real M_sun = 1.9891e33 rt;

    //! we'll get the composition indices from the network module
    int ihe4 = network_spec_index("helium-4");
    int ic12 = network_spec_index("carbon-12");
    int io16 = network_spec_index("oxygen-16");

    int ine20 = network_spec_index("neon-20");
    int ine23 = network_spec_index("neon-23");
    int ina23 = network_spec_index("sodium-23");
    int img23 = network_spec_index("magnesium-23");

    if (ihe4 < 0 || ic12 < 0 || io16 < 0) {
        amrex::Error("ERROR: species not defined");
    }

    //Initialize the temperature interior to fluff
    //(used to calculate fluff temperature for various values of fluff_type)
    Real temp_before_fluff = temp_fluff;

    //initialize the EOS and network

    //use_eos_coulomb comes in from extern_probin_module -- override
    //here if desired
    bool use_eos_coulomb = true;

    // Initialize the composition module
    init_urca_composition();

    //-----------------------------------------------------------------------------
    // Create a 1-d uniform grid that is identical to the mesh that we are
    // mapping onto, and then we want to force it into HSE on that mesh.
    //-----------------------------------------------------------------------------

    // compute the coordinates of the new gridded function

    Real delx = (xmax - xmin) / static_cast<Real>(nx);

    for (int i = 0; i < nx; ++i) {
        xznl(i) = xmin + static_cast<Real>(i) * delx;
        xznr(i) = xmin + (static_cast<Real>(i) + 1.0_rt) * delx;
        xzn_hse(i) = 0.5_rt * (xznl(i) + xznr(i));
    }

    // we start at the center of the WD and integrate outward.
    // Initialize the central conditions.

    eos_t eos_state;
    eos_state.T = temp_core;
    eos_state.rho = rho_c;
    set_urca_composition(eos_state);

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

    for (int i = 0; i < nx; ++i) {
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


    for (int i = 1; i < nx; ++i) {

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

        if (! fluff) {
            //not sure about theses
            Real p_want;
            Real drho;
            Real dtemp; 
            Real dpd;

            for (int iretry = 1; iretry < MAX_RETRY; ++iretry){
                Real try_tol = TOL*(1.1_rt**(iretry-1));

                for (int iter = 0; iter < MAX_ITER; ++iter) {
                    iter_dens_temp(dens_zone, temp_zone, eos_state, 
                                  g_zone, delx, TOL, dtol_fac, dens_conv_zone, low_density_cutoff, temp_fluff, 
                                  isentropic, test_hse_convergence, converged_hse, fluff, 
                                  fluff_type, 
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

            if (temp_zone < temp_fluff){
                if ( fluff_type == "constant"){
                    get_fluff_temperature(temp_zone, temp_fluff, temp_before_fluff, fluff_type);
                }
                isentropic = false;
            }   
        } else {
        dens_zone = low_density_cutoff;
        get_fluff_temperature(temp_zone, temp_fluff, temp_before_fluff, fluff_type);

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
            model_hse(i, ispec+n) = xn[n];
        }


        M_enclosed(i) = M_enclosed(i-1) +
                (4.0_rt / 3.0_rt) * M_PI * (xznr(i) - xznl(i)) *
                (std::pow(xznr(i), 2) + xznl(i) * xznr(i) + std::pow(xznl(i), 2)) * model_hse(i, idens);

        if ((M_enclosed(i) > M_conv_zone*M_sun) && isentropic){
            i_conv = i;
            isentropic = false;

        }
    }    

    std:cout << 'mass = ' << M_enclosed(nx)/M_sun << endl;

    //print/save model
    int ipos = prefix.find(".dat");
    if (ipos < 0) {
        ipos = prefix.find(".txt");
    }

    std::string outfile = model_file.substr(0, ipos) + ".hse." + num;

    std::cout << "writing " << model_name << " model to " << outfile << std::endl;

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


    // compute the maximum HSE error

    Real max_hse_error = -1.e30_rt;

    for (int i = 1; i < nx-1; ++i) {
        Real g_zone = -C::Gconst * M_enclosed(i-1) / (xznr(i-1) * xznr(i-1));
        Real dpdr = (model_hse(i, ipres) - model_hse(i-1, ipres)) / delx;
        Real rhog = 0.5_rt * (model_hse(i, idens) + model_hse(i-1, idens)) * g_zone;

        if (dpdr != 0.0_rt && model_hse(i+1, idens) > low_density_cutoff) {
            max_hse_error = amrex::max(max_hse_error, std::abs(dpdr - rhog) / std::abs(dpdr));
        }
    }

    // print out the stats
    std::cout << "maximum HSE error = " << max_hse_error << std::endl;
    std:cout << std::endl;

    std:cout << 'total mass (Msun) = ' << M_enclosed(i_fluff)/M_sun
    std:cout << 'convective zone mass (Msun) = ' << M_enclosed(i_conv)/M_sun << std:endl;
    std:cout << 'convective zone boundary density (g/cm^3) = ' << model_hse(i_conv, idens) << std:endl;
    std:cout << 'convective zone boundary radius (cm) = ' << xzn_hse(i_conv) << std:endl;
    std:cout << 'convective zone boundary index = ' << i_conv << std:endl;
}


AMREX_INLINE void get_fluff_temperature(Real temp, Real temp_fluff, Real temp_previous, std::string fluff_type){
    // There are 2 kinds of fluff temperature handling:
    // 1) fluff_type = "constant"   : fluff is at the temperature temp_fluff in the inputs
    // 2) fluff_type = "continuous" : fluff is at the same temperature as the material
    //                              immediately interior to the fluff.

    if (fluff_type == "constant") {
        temp = temp_fluff;
    } else if (fluff_type == "continuous") {
        temp = temp_previous;
    } else {
        call amrex_error("ERROR: invalid fluff_type");
    }
}


AMREX_INLINE void iter_dens_temp(Real dens_zone, Real temp_zone, eos_t eos_state_inwards, 
                          Real g_zone, Real delx, Real tol, Real dtol_fac, Real dens_conv_zone, low_density_cutoff, temp_fluff, 
                          bool isentropic, bool test_hse_convergence, bool converged_hse, bool fluff, 
                          std::string fluff_type, 
                          int i_fluff, int i_conv, int izone){
    
    Real dtemp, drho;

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

        Real dpt = eos_state.dpdt;
        Real dpd = eos_state.dpdr;
        Real dst = eos_state.dsdt;
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

        if (dens_zone < low_density_cutoff) {
            i_fluff = izone;

            dens_zone = low_density_cutoff;

            temp_before_fluff = eos_state_inwards.T;

            get_fluff_temperature(temp_zone, temp_fluff, temp_before_fluff, fluff_type);

            converged_hse = true;
            fluff = true;
            break;
        }

        if (dens_zone < dens_conv_zone && isentropic){
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
            dpdr = (pres_zone - eos_state_inwards.p)/delx;

            // Compute rho*g
            rhog = 0.5_rt*(dens_zone + eos_state_inwards.rho)*g_zone;

            if (std::abs(dpdr - rhog) < TOL*std::abs(dpdr)){
                converged_hse = true;
                break;
            }
        } else {
            // Check to see if drho is under the relative tolerance TOL
            if (abs(drho) < TOL*dens_zone){
                converged_hse = true;
                break;
            }
        }
    } else {
        // do isothermal
        p_want = eos_state_inwards.p + 
            delx*0.5_rt*(dens_zone + eos_state_inwards.rho)*g_zone;

        temp_zone = eos_state_inwards.T;


        eos_state.T     = temp_zone;
        eos_state.rho   = dens_zone;
        set_urca_composition(eos_state);

        // (t, rho) -> (p, s)
        eos(eos_input_rt, eos_state);

        entropy = eos_state.s;
        pres_zone = eos_state.p;

        dpd = eos_state.dpdr;

        drho = (p_want - pres_zone)/(dpd - 0.5_rt*delx*g_zone);

        dens_zone = amrex::max(0.9_rt*dens_zone,
            amrex::min(dens_zone + drho, 1.1_rt*dens_zone));

        if (dens_zone < low_density_cutoff) {
            i_fluff = izone;

            dens_zone = low_density_cutoff;

            temp_before_fluff = eos_state_inwards.T;

            get_fluff_temperature(temp_zone, temp_fluff, 
                temp_before_fluff, fluff_type);

            converged_hse = true;
            fluff = true;
            break;
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
            dpdr = (pres_zone - eos_state_inwards.p)/delx;

            // Compute rho*g
            rhog = 0.5_rt*(dens_zone + eos_state_inwards.rho)*g_zone;

            if (std::abs(dpdr - rhog) < TOL*std::abs(dpdr)){
                converged_hse = true;
                break;
            }
        } else {
            // Check to see if drho is under the relative tolerance TOL
            if (abs(drho) < TOL*dens_zone){
                converged_hse = true;
                break;
            }
        }
    }
}
#include <AMReX_Array.H>

#include <sstream>

#include <extern_parameters.H>
#include <fundamental_constants.H>

#include <coord_info.H>
#include <read_model.H>
#include <model_util.H>

using namespace amrex;

AMREX_INLINE void init_1d()
{

    const Real TOL = 1.e-10_rt;

    const int MAX_ITER = 250;

    const Real smallx = 1.e-10_rt;

    // get the species indices

    bool species_defined{true};

    int ifuel1 = network_spec_index(problem_rp::fuel1_name);
    if (ifuel1 < 0)
        species_defined = false;

    int ifuel2{-1};
    if (problem_rp::fuel2_name != "")
    {
        ifuel2 = network_spec_index(problem_rp::fuel2_name);
        if (ifuel2 < 0)
            species_defined = false;
    }

    int ifuel3{-1};
    if (problem_rp::fuel3_name != "")
    {
        ifuel3 = network_spec_index(problem_rp::fuel3_name);
        if (ifuel3 < 0)
            species_defined = false;
    }

    int ifuel4{-1};
    if (problem_rp::fuel4_name != "")
    {
        ifuel4 = network_spec_index(problem_rp::fuel4_name);
        if (ifuel4 < 0)
            species_defined = false;
    }

    int ifuel5{-1};
    if (problem_rp::fuel5_name != "")
    {
        ifuel5 = network_spec_index(problem_rp::fuel5_name);
        if (ifuel5 < 0)
            species_defined = false;
    }

    int ifuel6{-1};
    if (problem_rp::fuel6_name != "")
    {
        ifuel6 = network_spec_index(problem_rp::fuel6_name);
        if (ifuel6 < 0)
            species_defined = false;
    }

    int ifuel7{-1};
    if (problem_rp::fuel7_name != "")
    {
        ifuel7 = network_spec_index(problem_rp::fuel7_name);
        if (ifuel7 < 0)
            species_defined = false;
    }

    int iash1 = network_spec_index(problem_rp::ash1_name);
    if (iash1 < 0)
        species_defined = false;

    int iash2{-1};
    if (problem_rp::ash2_name != "")
    {
        iash2 = network_spec_index(problem_rp::ash2_name);
        if (iash2 < 0)
            species_defined = false;
    }

    int iash3{-1};
    if (problem_rp::ash3_name != "")
    {
        iash3 = network_spec_index(problem_rp::ash3_name);
        if (iash3 < 0)
            species_defined = false;
    }

    int iash4{-1};
    if (problem_rp::ash4_name != "")
    {
        iash4 = network_spec_index(problem_rp::ash4_name);
        if (iash4 < 0)
            species_defined = false;
    }

    if (!species_defined)
    {
        std::cout << ifuel1 << " " << ifuel2 << " " << ifuel3 << " " << ifuel4 << " " << ifuel5 << " " << ifuel6 << " " << ifuel7 << std::endl;
        std::cout << iash1 << " " << iash2 << " " << iash3 << " " << iash4 << std::endl;
        amrex::Error("ERROR: species not defined");
    }

    Real xn[NumSpec] = {0.0_rt};

    // set the composition of the underlying star

    Real xn_star[NumSpec] = {smallx};

    xn_star[iash1] = problem_rp::ash1_frac;
    if (problem_rp::ash2_name != "")
        xn_star[iash2] = problem_rp::ash2_frac;
    if (problem_rp::ash3_name != "")
        xn_star[iash3] = problem_rp::ash3_frac;
    if (problem_rp::ash4_name != "")
        xn_star[iash4] = problem_rp::ash4_frac;

    // and the composition of the accreted layer

    Real xn_base[NumSpec] = {smallx};

    xn_base[ifuel1] = problem_rp::fuel1_frac;
    if (problem_rp::fuel2_name != "")
        xn_base[ifuel2] = problem_rp::fuel2_frac;
    if (problem_rp::fuel3_name != "")
        xn_base[ifuel3] = problem_rp::fuel3_frac;
    if (problem_rp::fuel4_name != "")
        xn_base[ifuel4] = problem_rp::fuel4_frac;
    if (problem_rp::fuel5_name != "")
        xn_base[ifuel5] = problem_rp::fuel5_frac;
    if (problem_rp::fuel6_name != "")
        xn_base[ifuel6] = problem_rp::fuel6_frac;
    if (problem_rp::fuel7_name != "")
        xn_base[ifuel7] = problem_rp::fuel7_frac;

    // check if star comp. adds to 1
    double sum_star{0.0};
    for (auto e : xn_star) {
        sum_star += e;
    }
    if (std::abs(sum_star - 1.0_rt) > NumSpec * smallx) {
        amrex::Error("ERROR: ash mass fractions don't sum to 1");
    }

    // check if base comp. adds to 1
    double sum_base{0.0};
    for (auto e : xn_base) {
        sum_base += e;
    }
    if (std::abs(sum_base - 1.0_rt) > NumSpec * smallx) {
        amrex::Error("ERROR: fuel mass fractions don't sum to 1");
    }

    // Create a 1-d uniform grid that is identical to the mesh that we are
    // mapping onto, and then we want to force it into HSE on that mesh.

    if (problem_rp::nx > NPTS_MODEL) {
        amrex::Error("too many zones requested -- increase NPTS_MODEL");
    }

    Array1D<Real, 0, NPTS_MODEL-1> xznl_hse;
    Array1D<Real, 0, NPTS_MODEL-1> xzn_hse;
    Array1D<Real, 0, NPTS_MODEL-1> xznr_hse;

    Array2D<Real, 0, NPTS_MODEL-1, 0, model::nvar-1> model_hse;

    // compute the coordinates of the new gridded function

    fill_coord_arrays(xzn_hse, xznl_hse, xznr_hse);

    // find the index of the base height
    int index_base = -1;

    for (int i = 0; i < problem_rp::nx; ++i) {
        if (xzn_hse(i) >= problem_rp::xmin + problem_rp::H_star + problem_rp::delta) {
            index_base = i+1;
            break;
        }
    }

    if (index_base == -1) {
        amrex::Error("ERROR: base_height not found on grid");
    }


    // put the model onto our new uniform grid

    bool fluff = false;

    // determine the conditions at the base

    eos_t eos_state;
    eos_state.T = problem_rp::T_base;
    eos_state.rho = problem_rp::dens_base;
    for (int n = 0; n < NumSpec; ++n) {
        eos_state.xn[n] = xn_base[n];
    }

    eos(eos_input_rt, eos_state);

    // store the conditions at the base -- we'll use the entropy later
    // to constrain the isentropic layer

    Real pres_base = eos_state.p;
    Real entropy_base = eos_state.s;

    // set an initial temperature profile and composition

    for (int i = 0; i < problem_rp::nx; ++i) {

        // hyperbolic tangent transition:

        for (int n = 0; n < NumSpec; ++n) {
            model_hse(i, model::ispec+n) = xn_star[n] +
                0.5_rt * (xn_base[n] - xn_star[n]) *
                (1.0_rt + std::tanh((xzn_hse(i) - (problem_rp::xmin + problem_rp::H_star - problem_rp::delta) + problem_rp::delta) / problem_rp::delta));
        }

        model_hse(i, model::itemp) = problem_rp::T_star + 0.5_rt * (problem_rp::T_base - problem_rp::T_star) *
            (1.0_rt + std::tanh((xzn_hse(i) - (problem_rp::xmin + problem_rp::H_star - problem_rp::delta) + problem_rp::delta) / problem_rp::delta));

        // the density and pressure will be determined via HSE,
        // for now, set them to the base conditions

        model_hse(i, model::idens) = problem_rp::dens_base;
        model_hse(i, model::ipres) = pres_base;

    }


    if (problem_rp::index_base_from_temp == 1) {
        // find the index of the base height -- look at the temperature for this
        index_base = -1;
        for (int i = 0; i < problem_rp::nx; ++i) {

            if (model_hse(i, model::itemp) > 0.9995 * problem_rp::T_base) {
                index_base = i+1;
                break;
            }
        }

        if (index_base == -1) {
            amrex::Error("ERROR: base_height not found on grid");
        }

    }

    std::cout << "index_base = " << index_base << std::endl;

    // make the base thermodynamics consistent for this base point -- that is
    // what we will integrate from!

    eos_state.rho = model_hse(index_base, model::idens);
    eos_state.T = model_hse(index_base, model::itemp);
    for (int n = 0; n < NumSpec; ++n) {
        eos_state.xn[n] = model_hse(index_base, model::ispec+n);
    }

    eos(eos_input_rt, eos_state);

    model_hse(index_base, model::ipres) = eos_state.p;


    // HSE + entropy solve

    // the HSE state will be done putting creating an isentropic state until
    // the temperature goes below T_lo -- then we will do isothermal.
    // also, once the density goes below low_density_cutoff, we stop HSE

    bool isentropic = true;


    Real dens_zone;
    Real temp_zone;
    Real pres_zone;

    Real p_want;
    Real drho;
    Real dtemp;
    Real entropy;

    // integrate up

    for (int i = index_base+1; i < problem_rp::nx; ++i) {

        Real delx = xzn_hse(i) - xzn_hse(i-1);

        // compute the gravitation acceleration at the lower edge

        Real g_zone;
        if (problem_rp::do_invsq_grav == 1) {
            g_zone = -C::Gconst * problem_rp::M_enclosed / std::pow(xznl_hse(i), 2);
        } else {
            g_zone = problem_rp::g_const;
        }

        // we've already set initial guesses for density, temperature, and
        // composition

        dens_zone = model_hse(i, model::idens);
        temp_zone = model_hse(i, model::itemp);
        for (int n = 0; n < NumSpec; ++n) {
            xn[n] = model_hse(i, model::ispec+n);
        }

        // iteration loop


        // start off the Newton loop by saying that the zone has not converged
        bool converged_hse = false;

        if (! fluff) {

            for (int iter = 0; iter < MAX_ITER; ++iter) {

                if (isentropic) {

                    // get the pressure we want from the HSE equation, just the
                    // zone below the current.  Note, we are using an average of
                    // the density of the two zones as an approximation of the
                    // interface value -- this means that we need to iterate for
                    // find the density and pressure that are consistent

                    // furthermore, we need to get the entropy that we need,
                    // which will come from adjusting the temperature in
                    // addition to the density.

                    // HSE differencing

                    p_want = model_hse(i-1, model::ipres) +
                        delx * 0.5_rt * (dens_zone + model_hse(i-1, model::idens)) * g_zone;

                    // now we have two functions to zero:
                    //   A = p_want - p(rho,T)
                    //   B = entropy_base - s(rho,T)
                    // We use a two dimensional Taylor expansion and find the deltas
                    // for both density and temperature

                    // now we know the pressure and the entropy that we want, so we
                    // need to find the temperature and density through a two
                    // dimensional root find

                    // (t, rho) -> (p, s)
                    eos_state.T = temp_zone;
                    eos_state.rho = dens_zone;
                    for (int n = 0; n < NumSpec; ++n) {
                        eos_state.xn[n] = xn[n];
                    }

                    eos(eos_input_rt, eos_state);

                    entropy = eos_state.s;
                    pres_zone = eos_state.p;

                    Real dpt = eos_state.dpdT;
                    Real dpd = eos_state.dpdr;
                    Real dst = eos_state.dsdT;
                    Real dsd = eos_state.dsdr;

                    Real A = p_want - pres_zone;
                    Real B = entropy_base - entropy;

                    dtemp = ((dsd / (dpd - 0.5_rt * delx * g_zone)) * A - B)/
                        (dsd * dpt / (dpd - 0.5_rt * delx * g_zone) - dst);

                    drho = (A - dpt * dtemp) / (dpd - 0.5_rt * delx * g_zone);

                    dens_zone = amrex::max(0.9_rt * dens_zone,
                                           amrex::min(dens_zone + drho, 1.1_rt * dens_zone));

                    temp_zone = amrex::max(0.9_rt * temp_zone,
                                           amrex::min(temp_zone + dtemp, 1.1_rt * temp_zone));

                    // check if the density falls below our minimum cut-off --
                    // if so, floor it
                    if (dens_zone < problem_rp::low_density_cutoff) {
                        dens_zone = problem_rp::low_density_cutoff;
                        temp_zone = problem_rp::T_lo;
                        converged_hse = true;
                        fluff = true;
                        break;
                    }

                    if (std::abs(drho) < TOL * dens_zone &&
                        std::abs(dtemp) < TOL*temp_zone) {
                        converged_hse = true;
                        break;
                    }

                } else {

                    // do isothermal
                    p_want = model_hse(i-1, model::ipres) +
                        delx * 0.5_rt * (dens_zone + model_hse(i-1, model::idens)) * g_zone;

                    temp_zone = problem_rp::T_lo;

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

                    drho = (p_want - pres_zone) / (dpd - 0.5_rt * delx * g_zone);

                    dens_zone = amrex::max(0.9_rt * dens_zone,
                                           amrex::min(dens_zone + drho, 1.1_rt * dens_zone));

                    if (std::abs(drho) < TOL * dens_zone) {
                        converged_hse = true;
                        break;
                    }

                    if (dens_zone < problem_rp::low_density_cutoff) {
                        dens_zone = problem_rp::low_density_cutoff;
                        temp_zone = problem_rp::T_lo;
                        converged_hse = true;
                        fluff = true;
                        break;
                    }

                }

                if (temp_zone < problem_rp::T_lo) {
                    temp_zone = problem_rp::T_lo;
                    isentropic = false;
                }

            }


            if (! converged_hse) {
                std::cout << "Error zone " << i << " did not converge in init_1d" << std::endl;
                std::cout << "integrate up" << std::endl;
                std::cout << dens_zone << " " << temp_zone << std::endl;
                std::cout << p_want << " " << entropy_base << " " << entropy << std::endl;
                std::cout << drho << " " << dtemp << std::endl;
                amrex::Error("Error: HSE non-convergence");
            }

        } else {
            dens_zone = problem_rp::low_density_cutoff;
            temp_zone = problem_rp::T_lo;
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

        // to make this process converge faster, set the density in the
        // next zone to the density in this zone
        // model_hse(i+1, model::idens) = dens_zone

    }


    // integrate down -- using the temperature profile defined above

    for (int i = index_base-1; i >= 0; --i) {

        Real delx = xzn_hse(i+1) - xzn_hse(i);

        // compute the gravitation acceleration at the upper edge

        Real g_zone;
        if (problem_rp::do_invsq_grav == 1) {
            g_zone = -C::Gconst * problem_rp::M_enclosed / std::pow(xznr_hse(i), 2);
        } else {
            g_zone = problem_rp::g_const;
        }

        // we already set the temperature and composition profiles
        temp_zone = model_hse(i, model::itemp);
        for (int n = 0; n < NumSpec; ++n) {
            xn[n] = model_hse(i, model::ispec+n);
        }

        // use our previous initial guess for density

        dens_zone = model_hse(i+1, model::idens);


        // iteration loop

        // start off the Newton loop by saying that the zone has not converged
        bool converged_hse = false;

        for (int iter = 0; iter < MAX_ITER; ++iter) {

            // get the pressure we want from the HSE equation, just the
            // zone below the current.  Note, we are using an average of
            // the density of the two zones as an approximation of the
            // interface value -- this means that we need to iterate for
            // find the density and pressure that are consistent

            // HSE differencing

            p_want = model_hse(i+1, model::ipres) -
                delx * 0.5_rt * (dens_zone + model_hse(i+1, model::idens)) * g_zone;

            // we will take the temperature already defined in model_hse
            // so we only need to zero:
            //   A = p_want - p(rho)

            // (t, rho) -> (p)
            eos_state.T = temp_zone;
            eos_state.rho = dens_zone;
            for (int n = 0; n < NumSpec; ++n) {
                eos_state.xn[n] = xn[n];
            }

            eos(eos_input_rt, eos_state);

            pres_zone = eos_state.p;

            Real dpd = eos_state.dpdr;

            Real A = p_want - pres_zone;

            drho = A / (dpd + 0.5_rt * delx * g_zone);

            dens_zone = amrex::max(0.9_rt * dens_zone,
                                   amrex::min(dens_zone + drho, 1.1_rt * dens_zone));

            if (std::abs(drho) < TOL * dens_zone) {
                converged_hse = true;
                break;
            }

        }

        if (! converged_hse) {
            std::cout << "Error zone " << i << " did not converge in init_1d" << std::endl;
            std::cout << "integrate down" << std::endl;
            std::cout << dens_zone << " " << temp_zone << std::endl;
            std::cout << p_want << std::endl;
            std::cout << drho << std::endl;
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

    }


    auto deltastr = num_to_unitstring(problem_rp::delta);
    Real dCoord = xzn_hse(1) - xzn_hse(0);
    auto dxstr = num_to_unitstring(dCoord);

    std::string outfile = problem_rp::model_prefix + ".hse.tanh.delta_" + deltastr + ".dx_" + dxstr;
    std::string outfile2 = outfile + ".extras";

    std::ofstream of(outfile);
    std::ofstream of2(outfile2);

    of << "# npts = " << problem_rp::nx << std::endl;
    of << "# num of variables = " << 3 + NumSpec << std::endl;
    of << "# density" << std::endl;
    of << "# temperature" << std::endl;
    of << "# pressure" << std::endl;

    for (int n = 0; n < NumSpec; ++n) {
        of << "# " << spec_names_cxx[n] << std::endl;
    }

    for (int i = 0; i < problem_rp::nx; ++i) {
        of << std::setprecision(12) << std::setw(20) << xzn_hse(i) << " ";
        of << std::setprecision(12) << std::setw(20) << model_hse(i, model::idens) << " ";
        of << std::setprecision(12) << std::setw(20) << model_hse(i, model::itemp) << " ";
        of << std::setprecision(12) << std::setw(20) << model_hse(i, model::ipres) << " ";
        for (int n = 0; n < NumSpec; ++n) {
            of << std::setprecision(12) << std::setw(20) << model_hse(i, model::ispec+n) << " ";
        }
        of << std::endl;
    }

    of.close();

    // some metadata
    of << "# generated by toy_atm" << std::endl;
    //of << "# inputs file: " << params_file << std::endl;

    // extras file

    of2 << "# npts = " << problem_rp::nx << std::endl;
    of2 << "# num of variables = " << 2 << std::endl;
    of2 << "# entropy" << std::endl;
    of2 << "# c_s" << std::endl;

    for (int i = 0; i < problem_rp::nx; ++i) {
        eos_state.rho = model_hse(i, model::idens);
        eos_state.T = model_hse(i, model::itemp);
        for (int n = 0; n < NumSpec; ++n) {
            eos_state.xn[n] = model_hse(i, model::ispec+n);
        }

        eos(eos_input_rt, eos_state);

        of2 << std::setprecision(12) << std::setw(20) << xzn_hse(i) << std::endl;
        of2 << std::setprecision(12) << std::setw(20) << eos_state.s << std::endl;
        of2 << std::setprecision(12) << std::setw(20) << eos_state.cs << std::endl;
    }

    // compute the maximum HSE error

    Real max_hse_error = -1.e30;

    for (int i = 1; i < problem_rp::nx-1; ++i) {

        Real delx = xzn_hse(i) - xzn_hse(i-1);

        // compute the gravitation acceleration at the lower edge

        Real g_zone;
        if (problem_rp::do_invsq_grav == 1) {
            g_zone = -C::Gconst * problem_rp::M_enclosed / std::pow(xznl_hse(i), 2);
        } else {
            g_zone = problem_rp::g_const;
        }

        Real dpdr = (model_hse(i, model::ipres) - model_hse(i-1, model::ipres)) / delx;
        Real rhog = 0.5_rt * (model_hse(i, model::idens) + model_hse(i-1, model::idens)) * g_zone;

        if (dpdr != 0.0_rt && model_hse(i+1, model::idens) > problem_rp::low_density_cutoff) {
            max_hse_error = amrex::max(max_hse_error, std::abs(dpdr - rhog) / std::abs(dpdr));
        }
    }

    std::cout << "maximum HSE error = " << max_hse_error << std::endl;

}

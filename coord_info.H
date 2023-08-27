#ifndef COORD_INFO_H
#define COORD_INFO_H

#include <AMReX_Array.H>
#include <extern_parameters.H>

using namespace amrex;

///
/// get the number of radial points corresponding to an exact mapping from
/// the 3D spherical domain to a 1D radial.  See the MAESTROeX paper
/// for details
///
inline int
get_irreg_nr() {

    int nr{-1};

    if (problem_rp::use_irreg_grid) {
        nr = (3 * (problem_rp::nx/2 - 0.5_rt) * (problem_rp::nx/2 - 0.5_rt) - 0.75_rt) / 2.0_rt;
    } else {
        nr = problem_rp::nx;
    }
    return nr;
}

///
/// compute the coordinates of the new gridded function with irregular spacing
///
inline
void
fill_coord_arrays_irreg(const int nr,
                        Array1D<Real, 0, NPTS_MODEL-1>& xzn_hse,
                        Array1D<Real, 0, NPTS_MODEL-1>& xznl,
                        Array1D<Real, 0, NPTS_MODEL-1>& xznr,
                        Array1D<Real, 0, NPTS_MODEL-1>& delrl,
                        Array1D<Real, 0, NPTS_MODEL-1>& delrr) {

    Real dCoord = (problem_rp::xmax - problem_rp::xmin) / static_cast<Real>(problem_rp::nx);

    // note: this uses the uniform grid dx, regardless of whether we are doing irregular
    // or uniform gridding

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

}

///
/// compute the coordinates of the new gridded function with regular spacing
///
inline
void
fill_coord_arrays(Array1D<Real, 0, NPTS_MODEL-1>& xzn_hse,
                  Array1D<Real, 0, NPTS_MODEL-1>& xznl,
                  Array1D<Real, 0, NPTS_MODEL-1>& xznr) {

    Real dCoord = (problem_rp::xmax - problem_rp::xmin) / static_cast<Real>(problem_rp::nx);

    for (int i = 0; i < problem_rp::nx; ++i) {
        xznl(i) = problem_rp::xmin + static_cast<Real>(i) * dCoord;
        xznr(i) = problem_rp::xmin + (static_cast<Real>(i) + 1.0_rt) * dCoord;
        xzn_hse(i) = 0.5_rt * (xznl(i) + xznr(i));
    }
}

#endif

#ifndef INTERPOLATE_H
#define INTERPOLATE_H

#include <read_model.H>

AMREX_INLINE AMREX_GPU_HOST_DEVICE
int
locate(const Real r, const initial_model_t& initial_model) {

    int loc;

    if (r <= initial_model.r(0)) {
       loc = 0;

    } else if (r > initial_model.r(initial_model.npts-2)) {
       loc = initial_model.npts-1;

    } else {

        int ilo = 0;
        int ihi = initial_model.npts-2;

        while (ilo+1 != ihi) {
            int imid = (ilo + ihi) / 2;

            if (r <= initial_model.r(imid)) {
                ihi = imid;
            } else {
                ilo = imid;
            }
        }

        loc = ihi;
    }

    return loc;
}


AMREX_INLINE AMREX_GPU_HOST_DEVICE
Real
interpolate(const Real r, const int var_index, const initial_model_t& initial_model) {

    // find the value of model_state component var_index at point r
    // using linear interpolation.  Eventually, we can do something
    // fancier here.

    int id = locate(r, initial_model);

    Real slope;
    Real interp;

    if (id == 0) {

       slope = (initial_model.state(id+1, var_index) - initial_model.state(id, var_index)) /
           (initial_model.r(id+1) - initial_model.r(id));
       interp = slope * (r - initial_model.r(id)) + initial_model.state(id, var_index);

       // safety check to make sure interp lies within the bounding points

       Real minvar = amrex::min(initial_model.state(id+1, var_index),
                                initial_model.state(id, var_index));
       Real maxvar = amrex::max(initial_model.state(id+1, var_index),
                                initial_model.state(id, var_index));
       interp = amrex::max(interp, minvar);
       interp = amrex::min(interp, maxvar);

    } else if (id == initial_model.npts-1) {

       slope = (initial_model.state(id, var_index) - initial_model.state(id-1, var_index)) /
           (initial_model.r(id) - initial_model.r(id-1));
       interp = slope * (r - initial_model.r(id)) + initial_model.state(id, var_index);

       // safety check to make sure interp lies within the bounding points

       Real minvar = amrex::min(initial_model.state(id-1, var_index),
                                initial_model.state(id, var_index));
       Real maxvar = amrex::max(initial_model.state(id-1, var_index),
                                initial_model.state(id, var_index));
       interp = amrex::max(interp, minvar);
       interp = amrex::min(interp, maxvar);

    } else {

        if (r >= initial_model.r(id)) {

            slope = (initial_model.state(id+1, var_index) -
                     initial_model.state(id, var_index)) /
                (initial_model.r(id+1) - initial_model.r(id));
            interp = slope * (r - initial_model.r(id)) + initial_model.state(id, var_index);

        } else {

            slope = (initial_model.state(id, var_index) -
                     initial_model.state(id-1, var_index)) /
                (initial_model.r(id) - initial_model.r(id-1));
            interp = slope * (r - initial_model.r(id)) + initial_model.state(id, var_index);

        }

    }

    return interp;

}

#endif

#ifndef CFMT_MODULE
#define CFMT_MODULE

AMREX_INLINE Real cfmt(Real x) {

    Real safe_x;
    
    safe_x = x;
    if (std::abs(x) < 1.0e-98) {
        safe_x = 1.0e-98;
    }

    return safe_x;
}
#endif

: dipole_pp.mod - creates point process mechanism Dipole
:
: v 1.9.1m0
: rev 2015-12-15 (SL: minor)
: last rev: (SL: added Qtotal back, used for par calc)
: v 2.0.0
: rev 2020-02-19 (BC)
: last rev: (BC: made thread-safe by removing POINTERs)

NEURON {
    THREADSAFE
    POINT_PROCESS Dipole
    RANGE ri, ia, Q, ztan
    POINTER pv

    : for POINT_PROCESS. Gets additions from dipole
    POINTER Qtotal
}

UNITS {
    (nA)   = (nanoamp)
    (mV)   = (millivolt)
    (Mohm) = (megaohm)
    (um)   = (micrometer)
    (Am)   = (amp meter)
    (fAm)  = (femto amp meter)
}

ASSIGNED {
    ia (nA)
    ri (Mohm)
    pv (mV)
    v (mV)
    ztan (um)
    Q (fAm)
    Qtotal (fAm)
}

: solve for v's first then use them
AFTER SOLVE {
    ia = (pv - v) / ri
    Q = ia * ztan
    Qtotal = Qtotal + Q
}

AFTER INITIAL {
    ia = (pv - v) / ri
    Q = ia * ztan
    Qtotal = Qtotal + Q
}

: following needed for POINT_PROCESS only but will work if also in SUFFIX
BEFORE INITIAL {
    Qtotal = 0
}

BEFORE BREAKPOINT {
    Qtotal = 0
}

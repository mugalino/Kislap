// This defines a list of commonly used constants in 
// astrophysics in both CGS and MKS unit systems
// 
// Usage: include "./constants.hpp"
//        
// To use a particular constant 
//        <unit_system>::<unit>
//        e.g. cgs::GRAV_CONSTANT 

#ifndef __CONSTANTS_HPP__
#define __CONSTANTS_HPP__

namespace constants {
    #include <cmath>

    // Physical constants
    double GRAV_CONSTANT = 6.67428e-8; // dyne cm2 g-2
    double SPEED_LIGHT = 2.99792e10; // cm s-1 
    double RADIATION = 7.56577e-15; // erg cm-3 K-4

    double ELECTRON_CHARGE = 4.80325e-10; // statcoulomb
    double ELECTRON_MASS = 9.10938e-28; // g 

    double ELECTRON_VOLT = 1.60218e-12; // erg
    double PLANCK_H = 6.62607e-27; // erg s
    double PLANCK_HBAR = PLANCK_H / 2 * M_PI; // erg s
    double BOLTZMANN_CONST = 1.38065e-16; // erg K-1
    double PROTON_MASS = 1.67262e-24; // g 
    
    double BOHR_MAGNETON = 9.27401e-21; // erg Gauss-1 
    double BOHR_RADIUS = 5.29177e-9; // cm

    double RYDBERG_CONST = 1.09737e5; // cm-1
    double RYDBERG_FREQ = 3.28984e15; // s-1

    double STEF_BOLTZ = 5.67040e-5; // erg s-1 cm-2

    double THOMSON_CSECTION = 6.65245e-25; // cm2
    double ATOM_MASS =  1.66054e-24; // g

    // Astronomical constants 
    double AU = 1.49598e13; // cm
    double PC = 3.0856e18; // cm
    double LY = 9.4606e17; // cm

    double SOLAR_MASS = 1.989e33; // g
    double SOLAR_RADS = 6.9558e10; // cm
    double SOLAR_LUMS = 3.826e33; // erg s-1

} // namespace constants

#endif

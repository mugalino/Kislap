// Defines prefix and unit system conversion
// 
// Usage: include "./conversion.hpp"
//        
// To perform a conversion from cgs to mks 
// (default for astronomy and astrophysics)
//        <conversion_type>::<>
//        e.g. conversion.length
//
// These are all forward conversions. To do 
// the inverse, simply invert, i.e. 1/x.

#ifndef __CONVERSION_HPP__
#define __CONVERSION_HPP__

namespace conversion {
    double length = 1.e-2; // cm to m
    double mass = 1.e-3; // g to kg
    double energy = 1.e-7; // erg to J
    double force =  1.e-5; // dyne to N
    double power = energy; // erg s-1 to W or J/s
    double charge = 1. / (3.e-9); // statcoulomb to C
    double current = charge; // statampere to A or C/s
    double efield = 3.e4; // statvolt cm-1 to volt m-1
    double bfield = 1.e-4; // gauss to tesla
    double resistance = 9.e11; // s cm-1 to ohm
    double voltage = 3.e2; // statvolt to volt
}

#endif

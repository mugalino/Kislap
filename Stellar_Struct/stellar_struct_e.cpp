// Assuming that a star is isothermal, there are two equations and one
// closure that uniquely solves the stellar structure equations. This
// C++ code provides an illustration of how to use RK4 to solve a
// system of coupled ODEs like the stellar structure equations
//
// To compile: g++ stellar_struct.cpp; make -j
// To run, ./<executable>

// CENTIMETER-GRAM-SECOND unit system

#include <cmath>
#include <iostream>
#include <fstream>
#include "./constants.hpp"
#include "./conversion.hpp"

using namespace std;

// Initialize constants related to the polytropic equation of state
// At high central densities, the cold Fermi EoS approaches a
// relativistic form : P =  1.2435 * 10 ** 15 * (Y_e * rho) ** (4/3) erg cm**-3
// (https://cococubed.com/code_pages/coldwd.shtml)
//
// We assume this closure for this illustration.
double POLY_CONSTANT = 1.2435 * std::pow(10., 15.);
double POLY_INDEX = 4. / 3.;

// Assume typical central density for near-Chandrasekhar mass white dwarf
double CENTRAL_DENSITY = 2.2e11;
double CENTRAL_PRESSURE = POLY_CONSTANT * std::pow(CENTRAL_DENSITY, POLY_INDEX);

// Initial dr, 1 kilometer in CGS
double INIT_STEP = 1.0e5; 

// NOTE:
// Formula for a 4th order Runge-Kutta integrator
// k_1 = f(x_n, y_n)
// k_2 = f(x_n + dx / 2, y_n + dx * k_1 / 2)
// k_3 = f(x_n + dx / 2, y_n + dx * k_2 / 2)
// k_4 = f(x_n + dx, y_n + dx * k_3)
// y_(n+1) = y_n + (dx / 6) * (k_1 + 2 * k_2 + 2 * k_3 + k_4)
// x_(n+1) = x_n + dx

double mass_function(double density, double radius) {
  double mass = 4 * M_PI * radius * radius * density;

  return mass;
}

double pressure_function(double density, double radius, double mass) {
  double pressure = -1.0 * density * constants::GRAV_CONSTANT * mass / radius / radius;

  return pressure;
}

double density_function(double pressure) {
  double density = std::pow(pressure / POLY_CONSTANT, 1. / POLY_INDEX);

  return density;
}

int main() {
  std::cout << "Writing profile data into file" << std::endl;
  ofstream profile("profile.txt");
  std::cout << "Start RK4 step with initial value" << std::endl;
  std::cout << "Step:" << 0 << "\t "
            << "Radius:" << 0.0 << "\t "
            << "Density:" << CENTRAL_DENSITY << "\t "
            << "Pressure:" << CENTRAL_PRESSURE << std::endl;

  double mass = 0.0;
  double dens = CENTRAL_DENSITY;
  double pres = CENTRAL_PRESSURE;
  double rads = 0.0;
  double drad = INIT_STEP;
  double k1, k2, k3, k4;

  for (int idx = 0; dens > 1.e5; idx ++) {
    idx++;
    rads += drad;

    // RK4 solve for MASS
    k1 = mass_function(dens, rads);
    k2 = mass_function(dens, rads + drad/2);
    k3 = mass_function(dens, rads + drad/2);
    k4 = mass_function(dens, rads + drad);
    
    mass += (drad/6.) * (k1 + 2*k2 + 2*k3 + k4);

    // RK4 solve for PRESSURE
    k1 = pressure_function(dens, rads, mass);
    k2 = pressure_function(dens, rads + drad/2, mass);
    k3 = pressure_function(dens, rads + drad/2, mass);
    k4 = pressure_function(dens, rads + drad, mass);

    pres += (drad/6) * (k1 + 2*k2 + 2*k3 + k4);

    dens = density_function(pres);

    std::cout << "Step:" << idx << "\t "
              << "Radius:" << rads << "\t "
              << "Density:" << dens << "\t "
              << "Pressure:" << pres << "\t "
              << "Total Mass Enclosed:" << mass << std::endl;
    profile << rads << "\t " << dens << "\t " << pres <<"\t " << mass << "\n ";
  } // while (dens > 1.e5) {}
  profile.close();
}

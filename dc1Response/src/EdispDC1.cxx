/** 
 * @file EdispDC1.cxx
 * @brief Implementation for LAT energy dispersion class.
 * @author J. Chiang
 *
 * $Header$
 */

#include <cmath>

#include <sstream>
#include <stdexcept>

#include "CLHEP/Random/RandGauss.h"

using CLHEP::RandGauss;

#include "astro/SkyDir.h"

#include "EdispDC1.h"

namespace {
   auto erfcc(double x) -> double {
/* (C) Copr. 1986-92 Numerical Recipes Software 0@.1Y.. */
      double t, z, ans;
      
      z=fabs(x);
      t=1.0/(1.0+0.5*z);
      ans = t*exp(-z*z-1.26551223+t*(1.00002368+t*(0.37409196+t*(0.09678418+
            t*(-0.18628806+t*(0.27886807+t*(-1.13520398+t*(1.48851587+
            t*(-0.82215223+t*0.17087277)))))))));
      return x >= 0.0 ? ans : 2.0-ans;
   }
} // unnamed namespace

namespace dc1Response {

EdispDC1::EdispDC1(const std::string &filename) : DC1(filename) {
   normalizeDists();
}
     
EdispDC1::EdispDC1(const std::string &filename, int hdu, int npars) 
   : DC1(filename, hdu, npars) {
   normalizeDists();
}

auto EdispDC1::value(double appEnergy,
                       double energy, 
                       const astro::SkyDir &srcDir,
                       const astro::SkyDir &scZAxis,
                       const astro::SkyDir &, 
                       double time) const -> double {
// Inclination wrt spacecraft z-axis in degrees.
   double theta = srcDir.difference(scZAxis)*180./M_PI;

// The azimuthal angle is not used by the DC1 irfs.
   double phi(0);

   return value(appEnergy, energy, theta, phi, time);
}

auto EdispDC1::value(double appEnergy, double energy,
                       double theta, double phi, double time) const -> double {
   (void)(phi);
   (void)(time);

   if (theta < 0) {
      std::ostringstream message;
      message << "dc1Response::EdispDC1"
              << "::value(double, double, double, double):\n"
              << "theta cannot be less than zero. "
              << "Value passed: " << theta;
      throw std::invalid_argument(message.str());
   }

   if (theta < m_theta[m_theta.size()-1]) {
      std::vector<double> pars = fitParams(energy, theta);
      double x = appEnergy/energy;
      double my_value = pars[0]*exp(-0.5*pow((x - pars[1])/pars[2], 2));
      return my_value;
   }
   return 0;
}

auto EdispDC1::appEnergy(double energy,
                           const astro::SkyDir &srcDir,
                           const astro::SkyDir &scZAxis,
                           const astro::SkyDir &, 
                           double time) const -> double {
   (void)(time);
// Inclination wrt spacecraft z-axis in degrees.
   double inc = srcDir.difference(scZAxis)*180./M_PI;

   std::vector<double> pars = fitParams(energy, inc);

   double my_energy = (RandGauss::shoot()*pars[2] + pars[1])*energy;

   if (my_energy < 0) {
//       return 0;
      return energy;
   }
   return my_energy;
}

auto EdispDC1::integral(double emin, double emax, double energy,
                          const astro::SkyDir &srcDir, 
                          const astro::SkyDir &scZAxis,
                          const astro::SkyDir &scXAxis, 
                          double time) const -> double {
   (void)(scXAxis);
   (void)(time);
   double phi(0);
   double theta = srcDir.difference(scZAxis)*180./M_PI;
   return integral(emin, emax, energy, theta, phi);
}
   
auto EdispDC1::integral(double emin, double emax, double energy, 
                          double theta, double phi, double time) const -> double {
   (void)(phi);
   (void)(time);
   if (theta < 0) {
      std::ostringstream message;
      message << "dc1Response::EdispDC1"
              << "::integral(double, double, double, double, double):\n"
              << "theta cannot be less than zero. "
              << "Value passed: " << theta;
      throw std::invalid_argument(message.str());
   }

   std::vector<double> pars = fitParams(energy, theta);
   double E0 = energy*pars[1];
   double sigma = energy*pars[2];
   double scaledEmin = (emin - E0)/sqrt(2.)/sigma;
   double scaledEmax = (emax - E0)/sqrt(2.)/sigma;
   double my_integral = (::erfcc(scaledEmin) - ::erfcc(scaledEmax))/2.;
   return my_integral;
}

void EdispDC1::normalizeDists() {
// Compute proper Gaussian function normalization for each set of
// parameters.
   for (unsigned int ipar = 0; ipar < m_pars.size(); ipar++) {
      m_pars[ipar][0] = 1./sqrt(2.*M_PI)/m_pars[ipar][2];
   }
}

} // namespace dc1Response

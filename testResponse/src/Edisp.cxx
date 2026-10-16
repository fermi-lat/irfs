/** 
 * @file Edisp.cxx
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

#include "Edisp.h"

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

namespace testResponse {

auto Edisp::value(double appEnergy,
                    double energy, 
                    const astro::SkyDir &srcDir,
                    const astro::SkyDir &scZAxis,
                    const astro::SkyDir &,
                    double time) const -> double {
// Inclination wrt spacecraft z-axis in degrees.
   double theta = srcDir.difference(scZAxis)*180./M_PI;

   double phi(0);

   return value(appEnergy, energy, theta, phi, time);
}

auto Edisp::value(double appEnergy, double energy,
                    double theta, double phi, double time) const -> double {
   (void)(phi);
   (void)(time);

   if (theta < 0) {
      std::ostringstream message;
      message << "testResponse::Edisp"
              << "::value(double, double, double, double):\n"
              << "theta cannot be less than zero. "
              << "Value passed: " << theta;
      throw std::invalid_argument(message.str());
   }

   double sigma = m_eSigma*energy;
   double arg = (appEnergy - energy)/sigma;
   double my_value = exp(-arg*arg/2.)/sqrt(2.*M_PI)/sigma;
   return my_value;
}

auto Edisp::appEnergy(double energy,
                        const astro::SkyDir &,
                        const astro::SkyDir &,
                        const astro::SkyDir &,
                        double time) const -> double {
   (void)(time);

   double xi = RandGauss::shoot();
   double my_energy = (xi*m_eSigma + 1.)*energy;

   if (my_energy < 0) {
//       return 0;
      return energy;
   }
   return my_energy;
}

auto Edisp::integral(double emin, double emax, double energy,
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
   
auto Edisp::integral(double emin, double emax, double energy, 
                       double theta, double phi, double time) const -> double {
   (void)(phi);
   (void)(time);
   if (theta < 0) {
      std::ostringstream message;
      message << "testResponse::Edisp"
              << "::integral(double, double, double, double, double):\n"
              << "theta cannot be less than zero. "
              << "Value passed: " << theta;
      throw std::invalid_argument(message.str());
   }

   double E0 = energy;
   double sigma = energy*m_eSigma;
   double scaledEmin = (emin - E0)/sqrt(2.)/sigma;
   double scaledEmax = (emax - E0)/sqrt(2.)/sigma;
   double my_integral = (::erfcc(scaledEmin) - ::erfcc(scaledEmax))/2.;
   return my_integral;
}

} // namespace testResponse

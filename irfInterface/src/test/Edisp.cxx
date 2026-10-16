/**
 * @file Edisp.cxx
 * @brief Gaussian energy dispersion function.
 * @author J. Chiang
 *
 * $Header$
 */

#include <cmath>

#include "astro/SkyDir.h"

#include "Edisp.h"

Edisp::Edisp(double resolution) 
   : irfInterface::IEdisp(), m_resolution(resolution) {}

auto Edisp::value(double appEnergy,
                    double energy, 
                    const astro::SkyDir & srcDir, 
                    const astro::SkyDir & scZAxis,
                    const astro::SkyDir & scXAxis, 
                    double time) const -> double {
   (void)(scXAxis);
   static double phi(0);

   double theta(srcDir.difference(scZAxis)*180./M_PI);

   return value(appEnergy, energy, theta, phi, time);
}

auto Edisp::value(double appEnergy, double energy, double theta, 
                    double phi, double time) const -> double {
   (void)(theta);
   (void)(phi);
   (void)(time);

   double sigma(m_resolution*energy);
   double de(appEnergy - energy);

   return std::exp(-de*de/2./sigma/sigma)/std::sqrt(2.*M_PI)/sigma;
}

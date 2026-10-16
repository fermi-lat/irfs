/**
 * @file Psf.cxx
 * @brief Simple Psf to perform simple checks of IPsf default member functions.
 * @author J. Chiang
 *
 * $Header$
 */

#include <cmath>

#include "Psf.h"

Psf::Psf(double maxSep) : irfInterface::IPsf(), m_maxSep(maxSep) {
   m_value = 1./2./M_PI/(1. - std::cos(m_maxSep*M_PI/180.));
}

auto Psf::value(const astro::SkyDir & appDir,
                  double energy,
                  const astro::SkyDir & srcDir,
                  const astro::SkyDir & scZAxis,
                  const astro::SkyDir & scXAxis,
                  double time) const -> double {
   (void)(energy);
   (void)(scZAxis);
   (void)(scXAxis);
   (void)(time);
   double sep(srcDir.difference(appDir)*180./M_PI);
   if (sep < m_maxSep) {
      return m_value;
   }
   return 0;
}

auto Psf::value(double sep, double energy, double theta, double phi, 
                  double time) const -> double {
   (void)(energy);
   (void)(theta);
   (void)(phi);
   (void)(time);
   if (sep < m_maxSep) {
      return m_value;
   }
   return 0;
}

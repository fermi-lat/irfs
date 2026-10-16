/**
 * @file AeffEpochDep.cxx
 * @brief Class to handle Aeffs that change from epoch to epoch.  This
 * should eventually be refactored to be the "standard"
 * implementation, rather than treated as a special case.
 *
 * @author J. Chiang <jchiang@slac.stanford.edu>
 * 
 * $Header$
 */

#include "AeffEpochDep.h"

namespace latResponse {

AeffEpochDep::AeffEpochDep() : m_upperLimit(0) {
}

AeffEpochDep::AeffEpochDep(const AeffEpochDep & other) 
   : irfInterface::IAeff(other),
     EpochDep(other),
     m_upperLimit(other.m_upperLimit) {
   for (size_t i(0); i < other.m_aeffs.size(); i++) {
      m_aeffs.push_back(other.m_aeffs[i]->clone());
   }
}

auto AeffEpochDep::operator=(const AeffEpochDep & rhs) -> AeffEpochDep & {
   if (this != &rhs) {
      irfInterface::IAeff::operator=(rhs);
      EpochDep::operator=(rhs);
      m_upperLimit = rhs.m_upperLimit;
      for (size_t i(0); i < m_aeffs.size(); i++) {
         delete m_aeffs[i];
      }
      m_aeffs.clear();
      for (size_t i(0); i < rhs.m_aeffs.size(); i++) {
         m_aeffs.push_back(rhs.m_aeffs[i]->clone());
      }
   }
   return *this;
}

AeffEpochDep::~AeffEpochDep() {
   for (size_t i(0); i < m_aeffs.size(); i++) {
      delete m_aeffs[i];
   }
}

auto AeffEpochDep::value(double energy,
                           const astro::SkyDir & srcDir,
                           const astro::SkyDir & scZDir,
                           const astro::SkyDir & scXDir,
                           double time) const -> double {
   size_t indx(index(time));
   return m_aeffs[indx]->value(energy, srcDir, scZDir, scXDir, time);
}

auto AeffEpochDep::value(double energy, double theta, double phi, 
                           double time) const -> double {
   size_t indx(index(time));
   return m_aeffs[indx]->value(energy, theta, phi, time);
}

auto AeffEpochDep::upperLimit() const -> double {
   return m_upperLimit;
}

void AeffEpochDep::addAeff(const irfInterface::IAeff & aeff,
                           double epoch_start) {
   appendEpoch(epoch_start);
   m_aeffs.push_back(const_cast<irfInterface::IAeff &>(aeff).clone());
   if (aeff.upperLimit() > m_upperLimit) {
      m_upperLimit = aeff.upperLimit();
   }
}

} // namespace latResponse

/**
 * @file PsfEpochDep.cxx
 * @brief Definition for class that handles Psfs that change from
 * epoch to epoch.
 *
 * @author J. Chiang <jchiang@slac.stanford.edu>
 *
 * $Header$
 */

#include "PsfEpochDep.h"

namespace latResponse {

PsfEpochDep::PsfEpochDep() {
}

PsfEpochDep::PsfEpochDep(const PsfEpochDep & other)
   : irfInterface::IPsf(other), EpochDep(other) {
   for (size_t i(0); i < other.m_psfs.size(); i++) {
      m_psfs.push_back(other.m_psfs[i]->clone());
   }
}

auto PsfEpochDep::operator=(const PsfEpochDep & rhs) -> PsfEpochDep & {
   if (this != &rhs) {
      irfInterface::IPsf::operator=(rhs);
      EpochDep::operator=(rhs);
      for (size_t i(0); i < m_psfs.size(); i++) {
         delete m_psfs[i];
      }
      m_psfs.clear();
      for (size_t i(0); i < rhs.m_psfs.size(); i++) {
         m_psfs.push_back(rhs.m_psfs[i]->clone());
      }
   }
   return *this;
}

PsfEpochDep::~PsfEpochDep() {
   for (size_t i(0); i < m_psfs.size(); i++) {
      delete m_psfs[i];
   }
}

auto PsfEpochDep::value(const astro::SkyDir & appDir, 
                          double energy, 
                          const astro::SkyDir & srcDir, 
                          const astro::SkyDir & scZAxis,
                          const astro::SkyDir & scXAxis, 
                          double time) const -> double {
   size_t indx(index(time));
   return m_psfs[indx]->value(appDir, energy, srcDir, scZAxis, scXAxis, time);
}

auto PsfEpochDep::value(double separation, double energy, double theta,
                          double phi, double time) const -> double {
   size_t indx(index(time));
   return m_psfs[indx]->value(separation, energy, theta, phi, time);
}

auto PsfEpochDep::
angularIntegral(double energy,
                const astro::SkyDir & srcDir,
                double theta, 
                double phi, 
                const AcceptanceConeVector_t & acceptanceCones, 
                double time) -> double {
   size_t indx(index(time));
   return m_psfs[indx]->angularIntegral(energy, srcDir, theta, phi, 
                                        acceptanceCones, time);
}

auto PsfEpochDep::
angularIntegral(double energy, double theta, double phi,
                double radius, double time) const -> double {
   size_t indx(index(time));
   return m_psfs[indx]->angularIntegral(energy, theta, phi, radius, time);
}

void PsfEpochDep::addPsf(const irfInterface::IPsf & psf,
                         double epoch_start) {
   appendEpoch(epoch_start);
   m_psfs.push_back(const_cast<irfInterface::IPsf &>(psf).clone());
}

} // namespace latResponse

/**
 * @file EfficiencyFactorEpochDep.cxx
 * @brief Epoch-dependent efficiency factor class.
 *
 * @author J. Chiang <jchiang@slac.stanford.edu>
 *
 * $Header$
 */

#include "EfficiencyFactorEpochDep.h"

namespace latResponse {

EfficiencyFactorEpochDep::EfficiencyFactorEpochDep() {
}

EfficiencyFactorEpochDep::
EfficiencyFactorEpochDep(const EfficiencyFactorEpochDep & other) 
   : irfInterface::IEfficiencyFactor(other),
     EpochDep(other) {
   for (size_t i(0); i < other.m_effs.size(); i++) {
      m_effs.push_back(other.m_effs[i]->clone());
   }
}

auto 
EfficiencyFactorEpochDep::operator=(const EfficiencyFactorEpochDep & rhs) -> EfficiencyFactorEpochDep & {
   if (this != &rhs) {
      irfInterface::IEfficiencyFactor::operator=(rhs);
      EpochDep::operator=(rhs);
      for (size_t i(0); i < m_effs.size(); i++) {
         delete m_effs[i];
      }
      m_effs.clear();
      for (size_t i(0); i < rhs.m_effs.size(); i++) {
         m_effs.push_back(rhs.m_effs[i]->clone());
      }
   }
   return *this;
}

EfficiencyFactorEpochDep::~EfficiencyFactorEpochDep() throw() {
   for (size_t i(0); i < m_effs.size(); i++) {
      delete m_effs[i];
   }
}

auto EfficiencyFactorEpochDep::operator()(double energy, double met) const -> double {
   size_t indx(index(met));
   return m_effs[indx]->operator()(energy, met);
}

auto EfficiencyFactorEpochDep::value(double energy, double livetimefrac,
                                       bool front, double met) const -> double {
   size_t indx(index(met));
   return m_effs[indx]->value(energy, livetimefrac, front, met);
}

void EfficiencyFactorEpochDep::
getLivetimeFactors(double energy, double & factor1, double & factor2,
                   double met) const {
   size_t indx(index(met));
   m_effs[indx]->value(energy, factor1, factor2, met);
}

void EfficiencyFactorEpochDep::add(irfInterface::IEfficiencyFactor & eff,
                                   double epoch_start) {
   appendEpoch(epoch_start);
   m_effs.push_back(eff.clone());
}

} // namespace latResponse

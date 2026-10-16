/**
 * @file IEfficiencyFactor.cxx
 * @brief Function object that returns the IRF-dependent efficiency
 * corrections as a function of livetime fraction.
 * 
 * @author J. Chiang
 *
 * $Header$
 */

#include <cmath>

#include <algorithm>
#include <iostream>
#include <stdexcept>

#include "facilities/Util.h"

#include "tip/IFileSvc.h"
#include "tip/Table.h"

#include "st_facilities/Util.h"

#include "irfInterface/IEfficiencyFactor.h"

namespace irfInterface {

auto IEfficiencyFactor::value(double energy, double livetimefrac,
                                double met) const -> double {
   // Since we do not have access to the front and back effective
   // areas separately, just return the average.
   if (::getenv("OMIT_EFFICIENCY_FACTOR")) {
      return 1;
   }
   return (value(energy, livetimefrac, true, met) + 
           value(energy, livetimefrac, false, met))/2.;
}

void IEfficiencyFactor::readFt2File(std::string ft2file) {
   facilities::Util::expandEnvVar(&ft2file);
   const tip::Table * scData = 
      tip::IFileSvc::instance().readTable(ft2file, "SC_DATA");
   tip::Table::ConstIterator it(scData->begin());
   tip::ConstTableRecord & row(*it);
   for ( ; it != scData->end(); ++it) {
      m_start.push_back(row["START"].get());
      m_stop.push_back(row["STOP"].get());
      double duration(m_stop.back() - m_start.back());
      m_livetimefrac.push_back(row["livetime"].get()/duration);
   }
   delete scData;
}

void IEfficiencyFactor::clearFt2Data() {
   m_start.clear();
   m_stop.clear();
   m_livetimefrac.clear();
}

} // namespace irfInterface

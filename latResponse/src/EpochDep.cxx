/**
 * @file EpochDep.cxx
 * @brief Mix-in class to provide functionality needed for
 * epoch-dependent IRFs
 *
 * @author J. Chiang <jchiang@slac.stanford.edu>
 *
 * $Header$
 */

#include <algorithm>
#include <stdexcept>

#include "facilities/Util.h"

#include "tip/IFileSvc.h"
#include "tip/Table.h"

#include "astro/JulianDate.h"

#include "EpochDep.h"

namespace latResponse {

EpochDep::EpochDep() : m_curr_index(0) {
}

auto EpochDep::epochStart(const std::string & fitsfile,
                            const std::string & extname) -> double {
   const tip::Table * table 
      = tip::IFileSvc::instance().readTable(fitsfile, extname);
   const tip::Header & header(table->getHeader());
   std::string validity_start_date;
   header["CVSD0001"].get(validity_start_date);
   std::string validity_start_time;
   header["CVST0001"].get(validity_start_time);
   delete table;

   std::vector<std::string> date_tokens;
   facilities::Util::stringTokenize(validity_start_date, "-", date_tokens);
   std::vector<std::string> time_tokens;
   facilities::Util::stringTokenize(validity_start_time, ":", time_tokens);
   double hours = (std::atof(time_tokens[0].c_str()) 
                   + 60.*(std::atof(time_tokens[1].c_str()) 
                          + 60.*std::atof(time_tokens[2].c_str())));
   astro::JulianDate jd(std::atoi(date_tokens[0].c_str()),
                        std::atoi(date_tokens[1].c_str()),
                        std::atoi(date_tokens[2].c_str()),
                        hours);
   double met = jd.seconds() - astro::JulianDate::missionStart().seconds();
   return met;
}

void EpochDep::appendEpoch(double epoch_start) {
   m_epochStart.push_back(epoch_start);
}

auto EpochDep::index(double met) const -> size_t {
   if ( (met >= m_epochStart[m_curr_index] && 
         m_curr_index == m_epochStart.size()-1) ||
        (met >= m_epochStart[m_curr_index] && 
         met <= m_epochStart[m_curr_index + 1]) ) {
      return m_curr_index;
   } else if (met < m_epochStart.front()) {
      throw std::runtime_error("Requested MET not covered by selected IRFs.");
   } else if (met >= m_epochStart.back()) {
      m_curr_index = m_epochStart.size() - 1;
   } else {
      std::vector<double>::const_iterator it 
         =  std::upper_bound(m_epochStart.begin(), m_epochStart.end(), met);
      m_curr_index = it - m_epochStart.end() - 1;
   }
   return m_curr_index;
}

} // namespace latResponse

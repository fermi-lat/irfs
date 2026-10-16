/** 
 * @file Aeff.cxx
 * @brief Implementation for the DC2 effective Area class.
 * @author J. Chiang
 * 
 * $Header$
 */

#include <cmath>

#include <stdexcept>

#include "tip/IFileSvc.h"
#include "tip/Table.h"

#include "st_facilities/Util.h"

#include "astro/SkyDir.h"

#include "Aeff.h"

namespace dc2Response {

Aeff::Aeff(const std::string & filename, const std::string & extname) 
   : DC2(filename, extname), m_aeffMax(0) {
   readData();
}

Aeff::Aeff(const Aeff & rhs) 
   : IAeff(rhs), DC2(rhs), m_effArea(rhs.m_effArea),
     m_logElo(rhs.m_logElo), m_logEhi(rhs.m_logEhi),
     m_logE(rhs.m_logE), m_cosinc(rhs.m_cosinc), m_aeffMax(rhs.m_aeffMax) {
}

void Aeff::readData() {
   tip::IFileSvc & fileSvc(tip::IFileSvc::instance());

   const tip::Table * effArea = fileSvc.readTable(m_filename, m_extname);

   tip::Table::ConstIterator it(effArea->begin());
   tip::ConstTableRecord & row(*it);

   std::vector<double> effarea;
   row["effarea"].get(effarea);

   std::vector<double> elo, ehi;
   row["energ_lo"].get(elo);
   row["energ_hi"].get(ehi);

   for (size_t k = 0; k < elo.size(); k++) {
      m_logElo.push_back(std::log(elo.at(k)));
      m_logEhi.push_back(std::log(ehi.at(k)));
      m_logE.push_back((m_logElo.back() + m_logEhi.back())/2.);
   }

   std::vector<double> theta;
   row["theta"].get(theta);

   size_t indx(0);
   for (size_t i = 0; i < theta.size(); i++) {
      m_cosinc.push_back(std::cos(theta.at(i)));
      std::vector<double> row;
      for (size_t k = 0; k < m_logElo.size(); k++, indx++) {
// convert to cm^2 and take log for inevitable extrapolations outside of grid
         if (effarea.at(indx) <= 0) {
            effarea.at(indx) = 1e-8;
         }
         if (effarea.at(indx) > m_aeffMax) {
            m_aeffMax = effarea.at(indx);
         }
         row.push_back(std::log(effarea.at(indx)*1e4));  
      }
      m_effArea.push_back(row);
   }
   m_aeffMax *= 1e4; // convert to cm^2
   delete effArea;
}

auto Aeff::value(double energy, 
                   const astro::SkyDir & srcDir, 
                   const astro::SkyDir & scZAxis,
                   const astro::SkyDir &,
                   double time) const -> double {
// Inclination wrt spacecraft z-axis in radians.
   double theta = srcDir.difference(scZAxis);
   theta *= 180./M_PI;
   return value(energy, theta, 0., time);
}

auto Aeff::value(double energy, double theta, double phi,
                   double time) const -> double {
   (void)(phi);
   (void)(time);

   if (theta < 0) {
      std::ostringstream message;
      message << "dc2Response::Aeff::value(double, double, double):\n"
              << "theta cannot be less than zero. "
              << "Value passed: " << theta;
      throw std::invalid_argument(message.str());
   }

   double mu(std::cos(theta*M_PI/180.));
   double logE(std::log(energy));

   if (mu < m_cosinc.front() || mu > m_cosinc.back() || 
       logE < m_logElo.front() || logE > m_logEhi.back()) {
      return 0;
   }

   double my_value = st_facilities::Util::bilinear(m_cosinc, mu,
                                                   m_logE, logE, m_effArea);
   return std::exp(my_value);
}

auto Aeff::upperLimit() const -> double {
   return m_aeffMax;
}

} // namespace dc2Response

/** 
 * @file AeffDC1.cxx
 * @brief Implementation for the DC1 effective Area class.
 * @author Johann Cohen-Tanugi
 * 
 * $Header$
 */

#include <cmath>

#include <sstream>
#include <stdexcept>

#include "astro/SkyDir.h"

#include "st_facilities/FitsUtil.h"

#include "AeffDC1.h"

namespace dc1Response {

AeffDC1::AeffDC1(const std::string &filename, int hdu) 
   : DC1(filename, hdu) {
   read_FITS_table();
}

AeffDC1::~AeffDC1() {}

AeffDC1::AeffDC1(const AeffDC1 &rhs) 
   : IAeff(rhs), DC1(rhs), m_aeffTable(rhs.m_aeffTable),
     m_aeffMax(rhs.m_aeffMax) {
}

void AeffDC1::read_FITS_table() {
   std::string extName;
   st_facilities::FitsUtil::getFitsHduName(m_filename, m_hdu, extName);
   st_facilities::FitsUtil::getRecordVector(m_filename, extName, "energy_lo", 
                                            m_energy);

   std::vector<double> upper_bounds;
   st_facilities::FitsUtil::getRecordVector(m_filename, extName, "energy_hi", 
                                            upper_bounds);
   m_energy.push_back( *(upper_bounds.end() - 1) );

   st_facilities::FitsUtil::getRecordVector(m_filename, extName, "theta_lo",
                                            m_theta);
   st_facilities::FitsUtil::getRecordVector(m_filename, extName, "theta_hi", 
                                            upper_bounds);
   m_theta.push_back( *(upper_bounds.end() - 1) );

   m_theta[0] = 0.;  // kludge until aeff_DC1.fits is fixed.
   st_facilities::FitsUtil::getRecordVector(m_filename, extName, "effarea",
                                            m_aeffTable);

   m_aeffMax = m_aeffTable.at(0);
   for (size_t i=1; i < m_aeffTable.size(); i++) {
      if (m_aeffMax < m_aeffTable.at(i)) {
         m_aeffMax = m_aeffTable.at(i);
      }
   }
}

auto AeffDC1::value(double energy, 
                      const astro::SkyDir &srcDir, 
                      const astro::SkyDir &scZAxis,
                      const astro::SkyDir &,
                      double time) const -> double {
// Inclination wrt spacecraft z-axis in radians.
   double theta = srcDir.difference(scZAxis);
   theta *= 180./M_PI;
   return value(energy, theta, 0., time);
}

auto AeffDC1::value(double energy, double theta, double phi, 
                      double time) const -> double {
   (void)(phi);
   (void)(time);

   if (theta < 0) {
      std::ostringstream message;
      message << "dc1Response::AeffDC1::value(double, double, double):\n"
              << "theta cannot be less than zero. "
              << "Value passed: " << theta;
      throw std::invalid_argument(message.str());
   }

   double my_value(0);

   int indx = getAeffIndex(energy, theta);
   if (indx >= 0) {
      my_value = m_aeffTable[indx];
   }
   return my_value;
}

auto AeffDC1::upperLimit() const -> double {
   return m_aeffMax;
}
   
auto AeffDC1::getAeffIndex(double energy, double theta) const -> int {
// Find the location of energy and inclination in grid boundary
// arrays.
   int ien, ith;
   if (energy < *m_energy.begin() || energy >= *(m_energy.end() - 1)
       || theta < *m_theta.begin() || theta >= *(m_theta.end() - 1)) {
// Desired grid point is outside the valid range.
      return -1;
   }
   ien = find_iterator(m_energy, energy) - m_energy.begin();
   ith = find_iterator(m_theta, theta) - m_theta.begin();
   int indx = ien*(m_theta.size()-1) + ith;
   return indx;
}

} // namespace dc1Response

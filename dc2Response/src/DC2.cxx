/** 
 * @file DC2.cxx
 * @brief Implementation for DC2 base class.
 * @author J. Chiang
 * 
 * $Header$
 */

#include <cmath>

#include <algorithm>
#include <stdexcept>
#include <vector>

#include "DC2.h"

namespace dc2Response {

DC2::DC2(const std::string & filename, const std::string & extname) 
   : m_filename(filename), m_extname(extname) {}

DC2::DC2(const DC2 & rhs) {
   m_filename = rhs.m_filename;
   m_extname = rhs.m_extname;
}

auto 
DC2::find_iterator(const std::vector<double> & gridValues, double target) -> std::vector<double>::const_iterator {
   std::vector<double>::const_iterator it;
   if (target < gridValues.front()) {
      it = gridValues.begin();
   } else if (target >= gridValues.back()) {
      it = gridValues.end() - 2;
   } else {
      it = std::upper_bound(gridValues.begin(), gridValues.end(), target) - 1;
      if (!(*it <= target && *(it+1) >= target)) {
         throw std::range_error("DC2::find_iterator: target value is " 
                                "out of range.");
      }
   }
   return it;
}

} // namespace dc2Response

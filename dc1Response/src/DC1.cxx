/** 
 * @file DC1.cxx
 * @brief Implementation for DC1 base class.
 * @author J. Chiang
 * 
 * $Header$
 */

#include <cassert>
#include <cmath>

#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <vector>

#include "fitsio.h"

#include "tip/IFileSvc.h"
#include "tip/Table.h"

#include "st_facilities/FitsUtil.h"

#include "DC1.h"

namespace dc1Response {

DC1::DC1(const std::string &filename, bool havePars) 
   : m_filename(filename), m_have_FITS_data(false) {
   if (havePars) {
      readGridBoundaries();
      readFitParams();
   }
}

DC1::DC1(const std::string &filename, int hdu, int npars) 
   : m_filename(filename), m_hdu(hdu+1), m_npars(npars), 
     m_have_FITS_data(true) {
   readGridBoundaries();
   readFitParams();
}

DC1::DC1(const std::string &filename, int hdu) 
   : m_filename(filename), m_hdu(hdu+1), m_have_FITS_data(true) {
}

DC1::DC1(const DC1 &rhs) {
   m_filename = rhs.m_filename;
   m_hdu = rhs.m_hdu;
   m_npars = rhs.m_npars;
   m_have_FITS_data = rhs.m_have_FITS_data;
   m_pars = rhs.m_pars;
   m_theta = rhs.m_theta;
   m_energy = rhs.m_energy;
}

void DC1::readFitParams() {

   std::string extensionName;
   std::vector<std::string> colNames;

   if (m_have_FITS_data) {
      std::vector< std::vector<double> > my_params(m_npars);
      st_facilities::FitsUtil::getFitsHduName(m_filename, m_hdu,
                                              extensionName);
      st_facilities::FitsUtil::getFitsColNames(m_filename, m_hdu, colNames);
      for (int icol = 0; icol < m_npars; icol++) {
         st_facilities::FitsUtil::getRecordVector(m_filename, extensionName, 
                         colNames[icol], my_params[icol]);
      }
      int nrows = my_params[0].size();
      m_pars.resize(nrows);
      for (int i = 0; i < nrows; i++) {
         m_pars[i].resize(m_npars);
         for (int j = 0; j < m_npars; j++) {
            m_pars[i][j] = my_params[j][i];
         }
      }
   } else {
      extensionName = "fitParams";
      tip::Table * fitParams =
         tip::IFileSvc::instance().editTable(m_filename, extensionName);
      colNames = fitParams->getValidFields();
      m_npars = colNames.size() - 3;

      int nrows = fitParams->getNumRecords();
      m_pars.resize(nrows);
   
      tip::Table::Iterator it = fitParams->begin();
      tip::Table::Record & row = *it;
      for (int i = 0; it != fitParams->end() && i < nrows; ++it, ++i) {
         m_pars[i].resize(m_npars);
         for (int j = 0; j < m_npars; j++) {
            row[colNames[j]].get(m_pars[i][j]);
         }
      }
      delete fitParams;
   }
}

auto 
DC1::find_iterator(const std::vector<double> &gridValues, 
                   double target) const -> std::vector<double>::const_iterator {
   std::vector<double>::const_iterator it;
   if (target < *(gridValues.begin())) {
      it = gridValues.begin();
   } else if (target >= *(gridValues.end()-1)) {
      it = gridValues.end() - 2;
   } else {
// Can't user std::lower_bound here since it will return the element
// before the target if the target is equal to the element desired.
      it = std::upper_bound(gridValues.begin(), gridValues.end(), 
                            target) - 1;
      assert(*it <= target && *(it+1) >= target);
   }
   return it;
}

auto DC1::fitParams(double energy, double inc) const -> const std::vector<double> & {
   unsigned int indx = getParamsIndex(energy, inc);
   assert(indx < m_pars.size());
   return m_pars[indx];
}

auto DC1::getParamsIndex(double energy, double inc) const -> int {
// Find the location of energy and inclination in grid boundary
// arrays.
   int ien = find_iterator(m_energy, energy) - m_energy.begin();
   int ith = find_iterator(m_theta, inc) - m_theta.begin();

// The number of grid boundary points in each dimension is one greater
// than the number of bins.  Energy is traversed in the inner loop in
// irfAnalysis/MakeDists, and angles are traversed in the outer loop.
   int indx = ith*(m_energy.size()-1) + ien;

   return indx;
}

void DC1::readGridBoundaries() {
   if (m_have_FITS_data) {
      std::string extName;
      st_facilities::FitsUtil::getFitsHduName(m_filename, m_hdu, extName);
      st_facilities::FitsUtil::getRecordVector(m_filename, extName, "theta",
                                               m_theta);
      st_facilities::FitsUtil::getRecordVector(m_filename, extName, "energy",
                                               m_energy);
   } else {
      st_facilities::FitsUtil::getTableVector(m_filename, "thetaGrid", "theta",
                                              m_theta);
      st_facilities::FitsUtil::getTableVector(m_filename, "energyGrid",
                                              "energy", m_energy);
   }
}

} // namespace dc1Response

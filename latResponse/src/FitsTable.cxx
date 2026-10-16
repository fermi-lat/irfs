/**
 * @file FitsTable.cxx
 * @brief Implementation for FitsTable for accessing IRF data.
 *
 * @author J. Chiang <jchiang@slac.stanford.edu>
 *
 * $Header$
 */

#include <cmath>

#include <algorithm>
#include <iostream>
#include <stdexcept>

#include "tip/IFileSvc.h"
#include "tip/Table.h"

#include "latResponse/Bilinear.h"
#include "latResponse/FitsTable.h"

namespace {
   auto binIndex(double x, const std::vector<double> & xx) -> size_t {
      std::vector<double>::const_iterator ix = 
         std::upper_bound(xx.begin(), xx.end(), x);
      return ix - xx.begin();
   }
}

namespace latResponse {

FitsTable::FitsTable(const std::string & filename,
                     const std::string & extname,
                     const std::string & tablename,
                     size_t nrow) : m_interpolator(0) {

   const tip::Table * table(tip::IFileSvc::instance().readTable(filename, 
                                                                extname));

   std::vector<double> elo, ehi;
   getVectorData(table, "ENERG_LO", elo, nrow);
   getVectorData(table, "ENERG_HI", ehi, nrow);
   for (size_t k(0); k < elo.size(); k++) {
      m_ebounds.push_back(std::log10(elo.at(k)));
      m_logEnergies.push_back(std::log10(std::sqrt(elo.at(k)*ehi.at(k))));
   }
   m_ebounds.push_back(std::log10(ehi.back()));

   std::vector<double> mulo, muhi;
   getVectorData(table, "CTHETA_LO", mulo, nrow);
   getVectorData(table, "CTHETA_HI", muhi, nrow);
   for (size_t i(0); i < muhi.size(); i++) {
      m_tbounds.push_back(mulo.at(i));
      m_mus.push_back((m_tbounds.at(i) + muhi.at(i))/2.);
   }
   m_tbounds.push_back(muhi.back());

   m_minCosTheta = mulo.front();

   getVectorData(table, tablename, m_values, nrow);
   m_maxValue = m_values.front();
   for (size_t i(1); i < m_values.size(); i++) {
      if (m_values.at(i) > m_maxValue) {
         m_maxValue = m_values.at(i);
      }
   }

// Replicate nasty THF2 and RootEval::Table behavior from handoff_response,
// by passing xlo, xhi, ylo, yhi values
   double xlo, xhi, ylo, yhi;
   m_interpolator = new Bilinear(m_logEnergies, m_mus, m_values,
                                 xlo=0., xhi=10., ylo=-1., yhi=1.);

   delete table;
}

FitsTable::FitsTable() : m_interpolator(0) {}

FitsTable::FitsTable(const FitsTable & rhs) 
   : m_interpolator(0), m_logEnergies(rhs.m_logEnergies), m_mus(rhs.m_mus),
     m_values(rhs.m_values), m_ebounds(rhs.m_ebounds),
     m_tbounds(rhs.m_tbounds), m_minCosTheta(rhs.m_minCosTheta), 
     m_maxValue(rhs.m_maxValue) {
   m_interpolator = new Bilinear(m_logEnergies, m_mus, m_values,
                                 0, 10, -1, 1);
}

auto FitsTable::operator=(const FitsTable & rhs) -> FitsTable & {
   if (this != &rhs) {
      delete m_interpolator;
      m_interpolator = new Bilinear(*rhs.m_interpolator);
      m_logEnergies = rhs.m_logEnergies;
      m_mus = rhs.m_mus;
      m_values = rhs.m_values;
      m_ebounds = rhs.m_ebounds;
      m_tbounds = rhs.m_tbounds;
      m_minCosTheta = rhs.m_minCosTheta;
      m_maxValue = rhs.m_maxValue;
   }
   return *this;
}

FitsTable::~FitsTable() { 
   delete m_interpolator;
}

auto FitsTable::
value(double logenergy, double costh, bool interpolate) const -> double {
   if (interpolate) {
      if (costh > m_mus.back()) {
         costh = m_mus.back();
      }
      return (*m_interpolator)(logenergy, costh);
   }

   // if (logenergy <= m_logEnergies.at(1)) { // use first bin
   //    logenergy = m_logEnergies.at(1);
   // }

   size_t ix = binIndex(logenergy, m_ebounds);
   if (ix == m_ebounds.size()) {
      ix -= 1;
   }
   if (ix == 0) {
      ix = 1;
   }
   size_t iy = binIndex(costh, m_tbounds);
   if (iy == 0) {
      iy = 1;
   }
   if (iy > m_mus.size()) {
      iy = m_mus.size();
   }
   size_t indx = (iy - 1)*m_logEnergies.size() + ix - 1;

   return m_values.at(indx);
}

void FitsTable::getValues(std::vector<double> & values) const {
   values.clear();
   for (size_t i(0); i < m_values.size(); i++) {
      values.push_back(m_values.at(i));
   }
}

void FitsTable::getCornerPars(double logE, double costh,
                              double & tt, double & uu,
                              std::vector<double> & cornerEnergies,
                              std::vector<double> & cornerPars) const {
   std::vector<double> corner_logE(4, 0);
   std::vector<double> corner_costh(4, 0);
   m_interpolator->getCorners(logE, costh, tt, uu, &corner_logE[0],
                              &corner_costh[0], &cornerPars[0]);
   for (size_t i(0); i < corner_logE.size(); i++) {
      cornerEnergies[i] = std::pow(10., corner_logE[i]);
   }
}

auto FitsTable::getPar(size_t ilogE, size_t icosth) const -> double {
   return m_interpolator->getPar(ilogE, icosth);
}

void FitsTable::setPar(size_t ilogE, size_t icosth, double value) {
   m_interpolator->setPar(ilogE, icosth, value);
   m_values.at(icosth*m_logEnergies.size() + ilogE) = value;
}

void FitsTable::getVectorData(const tip::Table * table,
                              const std::string & fieldName,
                              std::vector<double> & values,
                              size_t nrow) {
   values.clear();

   tip::Table::ConstIterator it(table->begin());
   tip::ConstTableRecord & row(*it);

   for (size_t i(0); i < nrow; i++) {
      ++it;
   }
   std::vector<float> my_values;
   row[fieldName].get(my_values);

   for (size_t i(0); i < my_values.size(); i++) {
      values.push_back(my_values.at(i));
   }
}

void FitsTable::setValues(const std::vector<double>& values) {
  if(values.size() != m_values.size())
    throw std::runtime_error("Wrong size for parameter array.");
  m_values = values;
  delete m_interpolator;
  m_interpolator = new Bilinear(m_logEnergies, m_mus, m_values,
				0, 10, -1, 1);
  m_maxValue = *std::max_element(m_values.begin(),m_values.end());
}

} // namespace latResponse

/**
 * @file EdispInterpolator.cxx
 * @brief Bilinear interoplator of energy dispersion distributions.
 * @author J. Chiang
 *
 * $Header$
 */

#include <cmath>

#include <algorithm>
#include <iostream>
#include <sstream>
#include <stdexcept>

#include "tip/IFileSvc.h"
#include "tip/Table.h"

#include "latResponse/FitsTable.h"

#include "latResponse/EdispInterpolator.h"

namespace {
   auto sqr(double x) -> double {
      return x*x;
   }
   class Array {
   public:
      Array(const std::vector<double> & values, size_t nx) 
         : m_values(values), m_nx(nx) {}
      auto operator()(size_t iy, size_t ix) const -> double {
         return m_values[iy*m_nx + ix];
      }
   private:
      const std::vector<double> & m_values;
      size_t m_nx;
   };
}

namespace latResponse {

EdispInterpolator::EdispInterpolator(const std::string & fitsfile,
                                     const std::string & extname,
                                     size_t nrow)
   : m_fitsfile(fitsfile), m_extname(extname), m_nrow(nrow),
     m_renormalized(false) {
   readFits();
}

EdispInterpolator::~EdispInterpolator() throw() {
}

void EdispInterpolator::readFits() {
   tip::IFileSvc & fileSvc(tip::IFileSvc::instance());
   const tip::Table * table(fileSvc.readTable(m_fitsfile, m_extname));
   const std::vector<std::string> & validFields(table->getValidFields());

   // The first four columns *must* be "ENERG_LO", "ENERG_HI", "CTHETA_LO",
   // "CTHETA_HI", in that order.
   const char * boundsName[] = {"energ_lo", "energ_hi", 
                                "ctheta_lo", "ctheta_hi"};
   size_t numBoundsCols(4);
   for (size_t i(0); i < numBoundsCols; i++) {
      if (validFields.at(i) != boundsName[i]) {
         std::ostringstream message;
         message << "latResponse::ParTables::ParTables: "
                 << "invalid header in " << m_fitsfile << "  "
                 << validFields.at(i) << "  " << i;
         throw std::runtime_error(message.str());
      }
   }

   // Push boundary values onto energy and theta arrays, replicating
   // parameter values along outer boundary.

   std::vector<double> elo, ehi;
   FitsTable::getVectorData(table, "ENERG_LO", elo, m_nrow);
   FitsTable::getVectorData(table, "ENERG_HI", ehi, m_nrow);
   std::vector<double> logEs;
   for (size_t k(0); k < elo.size(); k++) {
      logEs.push_back(std::log10(std::sqrt(elo[k]*ehi[k])));
   }

   std::vector<double> mulo, muhi;
   FitsTable::getVectorData(table, "CTHETA_LO", mulo, m_nrow);
   FitsTable::getVectorData(table, "CTHETA_HI", muhi, m_nrow);
   std::vector<double> cosths;
   for (size_t i(0); i < muhi.size(); i++) {
      cosths.push_back((mulo[i] + muhi[i])/2.);
   }

   size_t par_size(elo.size()*mulo.size());

   std::vector<double> values;
   for (size_t i(numBoundsCols); i < validFields.size(); i++) {
      const std::string & tablename(validFields[i]);
      FitsTable::getVectorData(table, tablename, values, m_nrow);
      if (values.size() != par_size) {
         std::ostringstream message;
         message << "Parameter array size does not match "
                 << "expected size based on energy and costheta arrays "
                 << "for table " << tablename
                 << " in  " << m_fitsfile;
         throw std::runtime_error(message.str());
      }

      std::vector<double> my_values;
      generateBoundaries(logEs, cosths, values, 
                         m_logEs, m_cosths, my_values);

      if (i == numBoundsCols) {
         m_parVectors.resize(m_logEs.size()*m_cosths.size(),
                             std::vector<double>());
      }
      for (size_t j(0); j < my_values.size(); j++) {
         m_parVectors[j].push_back(my_values[j]);
      }
   }
   for (size_t k(0); k < m_logEs.size(); k++) {
      m_energies.push_back(std::pow(10., m_logEs[k]));
   }
   for (size_t j(0); j < m_cosths.size(); j++) {
      m_thetas.push_back(std::acos(m_cosths[j])*180./M_PI);
   }
   if (m_parVectors[0].size() != (validFields.size()-numBoundsCols)) {
      std::ostringstream message;
      message << "Number of PSF parameters in "
              << m_fitsfile
              << " does no match the expected number of 6.";
      throw std::runtime_error(message.str());
   }
   delete table;
}

void EdispInterpolator::generateBoundaries(const std::vector<double> & x,
                                           const std::vector<double> & y,
                                           const std::vector<double> & values,
                                           std::vector<double> & xout,
                                           std::vector<double> & yout,
                                           std::vector<double> & values_out, 
                                           double xlo, double xhi,
                                           double ylo, double yhi) {
   xout.resize(x.size() + 2);
   std::copy(x.begin(), x.end(), xout.begin() + 1);
   xout.front() = xlo;
   xout.back() = xhi;

   yout.resize(y.size() + 2);
   std::copy(y.begin(), y.end(), yout.begin() + 1);
   yout.front() = ylo;
   yout.back() = yhi;

   Array array(values, x.size());
   values_out.push_back(array(0, 0));
   for (size_t i(0); i < x.size(); i++) {
      values_out.push_back(array(0, i));
   }
   values_out.push_back(array(0, x.size()-1));
   for (size_t j(0); j < y.size(); j++) {
      values_out.push_back(array(j, 0));
      for (size_t i(0); i < x.size(); i++) {
         values_out.push_back(array(j, i));
      }
      values_out.push_back(array(j, x.size()-1));
   }
   values_out.push_back(array(y.size()-1, 0));
   for (size_t i(0); i < x.size(); i++) {
      values_out.push_back(array(y.size()-1, i));
   }
   values_out.push_back(array(y.size()-1, x.size()-1));
}

void EdispInterpolator::getCornerPars(double energy, double theta,
                                      double phi, double time,
                                      double & tt, double & uu,
                                      std::vector<double> & cornerEnergies,
                                      std::vector<double> & cornerThetas,
                                      std::vector<size_t> & indx) const {
   (void)(phi);
   (void)(time);
   double logE(std::log10(energy));
   double costh(std::cos(theta*M_PI/180.));
   int i(findIndex(m_logEs, logE));
   int j(findIndex(m_cosths, costh));

   tt = (logE - m_logEs[i-1])/(m_logEs[i] - m_logEs[i-1]);
   uu = (costh - m_cosths[j-1])/(m_cosths[j] - m_cosths[j-1]);
   cornerEnergies[0] = m_energies[i-1];
   cornerEnergies[1] = m_energies[i];
   cornerEnergies[2] = m_energies[i];
   cornerEnergies[3] = m_energies[i-1];
   cornerThetas[0] = m_thetas[j-1];
   cornerThetas[1] = m_thetas[j-1];
   cornerThetas[2] = m_thetas[j];
   cornerThetas[3] = m_thetas[j];

   size_t xsize(m_energies.size());
   indx[0] = xsize*(j-1) + (i-1);
   indx[1] = xsize*(j-1) + (i);
   indx[2] = xsize*(j) + (i);
   indx[3] = xsize*(j) + (i-1);
}

auto EdispInterpolator::findIndex(const std::vector<double> & xx, double x) -> int {
   typedef std::vector<double>::const_iterator const_iterator_t;

   const_iterator_t ix(std::upper_bound(xx.begin(), xx.end(), x));
   if (ix == xx.end() && x != xx.back()) {
      std::cout << xx.front() << "  "
                << x << "  "
                << xx.back() << std::endl;
      throw std::invalid_argument("EdispInterpolator::findIndex: x out of range");
   }
   if (x == xx.back()) {
      ix = xx.end() - 1;
   } else if (x <= xx.front()) {
      ix = xx.begin() + 1;
   }
   int i(ix - xx.begin());
   return i;
}

auto EdispInterpolator::params(size_t indx) const -> std::vector<double> { 

  if(indx >= m_parVectors[0].size())
    throw std::runtime_error("Parameter index out of range.");

  std::vector<double> vals(m_parVectors.size(),0.0);
  std::vector< std::vector<double> >::const_iterator itr0 = m_parVectors.begin();
  std::vector<double>::iterator itr1 = vals.begin();
  for( ;itr0 != m_parVectors.end() && itr1 != vals.end();
       ++itr0, ++itr1)
    *itr1 = (*itr0)[indx];

  return vals;
}

void EdispInterpolator::setParams(size_t indx, const std::vector<double>& params) {
  
  if(indx >= m_parVectors[0].size())
    throw std::runtime_error("Parameter index out of range.");
  else if(params.size() != m_parVectors.size())
    throw std::runtime_error("Wrong size for parameter array.");

  std::vector<double>::const_iterator itr0 = params.begin();
  std::vector< std::vector<double> >::iterator itr1 = m_parVectors.begin();

  m_renormalized = false;

  for( ;itr0 != params.end() && itr1 != m_parVectors.end();
       ++itr0, ++itr1)
    (*itr1)[indx] = *itr0;
}

} // namespace latResponse

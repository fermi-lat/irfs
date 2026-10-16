/** 
 * @file Edisp.cxx
 * @brief Implementation for LAT energy dispersion class.
 * @author J. Chiang
 *
 * $Header$
 */

#include <cmath>

#include <algorithm>
#include <iostream>
#include <sstream>
#include <stdexcept>

#include "CLHEP/Random/RandFlat.h"

using CLHEP::RandFlat;

#include <tip/IFileSvc.h>
#include <tip/Table.h>

#include "st_facilities/dgaus8.h"

#include "astro/SkyDir.h"

#include "Edisp.h"

namespace {
   auto reverse_cmp(const double & x, const double & y) -> bool {
      return x > y;
   }
}

namespace dc2Response {

double Edisp::s_rwidth;
double Edisp::s_ltail;

std::vector<double> Edisp::s_xvals;

Edisp::Edisp(const std::string & fitsfile, const std::string & extname) 
   : DC2(fitsfile, extname) {
   readData();
}

auto Edisp::value(double appEnergy,
                    double energy, 
                    const astro::SkyDir &srcDir,
                    const astro::SkyDir &scZAxis,
                    const astro::SkyDir &,
                    double time) const -> double {
// Inclination wrt spacecraft z-axis in degrees.
   double theta = srcDir.difference(scZAxis)*180./M_PI;

// The azimuthal angle is not used by the DC2 irfs.
   double phi(0);

   return value(appEnergy, energy, theta, phi, time);
}

auto Edisp::value(double appEnergy, double energy,
                    double theta, double phi, double time) const -> double {
   (void)(phi);
   (void)(time);

   if (theta < 0) {
      std::ostringstream message;
      message << "dc2Response::Edisp"
              << "::value(double, double, double, double):\n"
              << "theta cannot be less than zero. "
              << "Value passed: " << theta;
      throw std::invalid_argument(message.str());
   }

   double mu(std::cos(theta*M_PI/180.));
   size_t indx = parIndex(energy, mu);

   double p1 = m_ltail.at(indx);
   double p2 = m_rwidth.at(indx);
   double norm = m_norms.at(indx);

   double x((appEnergy - energy)/energy);
   double arg(x/p2);
   if (arg > 40) {
      return norm*std::pow(1. + x, p1)*std::exp(-arg)/energy;
   }
   return norm*std::pow(1. + x, p1)/(1. + std::exp(arg))/energy;
}

auto Edisp::appEnergy(double energy,
                        const astro::SkyDir & srcDir,
                        const astro::SkyDir & scZAxis,
                        const astro::SkyDir &,
                        double time) const -> double {
   (void)(time);
   double mu(std::cos(srcDir.difference(scZAxis)));
   size_t indx = parIndex(energy, mu);

   float xi(RandFlat::shoot());

   size_t j = std::upper_bound(m_cumDists.at(indx).begin(), 
                                m_cumDists.at(indx).end(), xi)
      - m_cumDists.at(indx).begin() - 1;
   return energy*((xi - m_cumDists.at(indx).at(j))
                  /(m_cumDists.at(indx).at(j+1) - m_cumDists.at(indx).at(j))
                  *(s_xvals.at(j+1) - s_xvals.at(j)) + s_xvals.at(j))
      + energy;
}

auto Edisp::parIndex(double energy, double mu) const -> size_t {
   size_t k = std::upper_bound(m_eBounds.begin(), m_eBounds.end(),
                               energy) - m_eBounds.begin() - 1;
   size_t i = std::upper_bound(m_muBounds.begin(), m_muBounds.end(),
                               mu, ::reverse_cmp) - m_muBounds.begin() - 1;
   if (i > m_muBounds.size() - 2) {
      i = m_muBounds.size() - 2;
   }
   if (k > m_eBounds.size() - 2) {
      k = m_eBounds.size() - 2;
   }
   size_t indx(i*(m_eBounds.size()-1) + k);
   if (indx > m_cumDists.size() - 1) {
      throw std::runtime_error("dc2Response::Edisp::parIndex: "
                               "index out-of-range");
   }
   return indx;
}

// double Edisp::integral(double emin, double emax, double energy,
//                        const astro::SkyDir & srcDir, 
//                        const astro::SkyDir & scZAxis,
//                        const astro::SkyDir &,
//                        double time) const {
//    return integral(emin, emax, energy,
//                    srcDir.difference(scZAxis)*180./M_PI, 0, time);
// }
   
// double Edisp::integral(double emin, double emax, double energy, 
//                        double theta, double phi, double time) const {
//    (void)(phi);
//    (void)(time);
//    if (theta < 0) {
//       std::ostringstream message;
//       message << "dc2Response::Edisp"
//               << "::integral(double, double, double, double, double):\n"
//               << "theta cannot be less than zero. "
//               << "Value passed: " << theta;
//       throw std::invalid_argument(message.str());
//    }
//    double mu(std::cos(theta*M_PI/180.));
//    size_t indx = parIndex(energy, mu);
//    s_rwidth = m_rwidth.at(indx);
//    s_ltail = m_ltail.at(indx);
//    double lowerLim((emin - energy)/energy);
//    double upperLim((emax - energy)/energy);
//    double err(1e-5);
//    double my_integral;
//    long ierr;
//    dgaus8_(&edispIntegrand, &lowerLim, &upperLim, &err, &my_integral, &ierr);
//    return m_norms.at(indx)*my_integral;
// }

void Edisp::readData() {
   tip::IFileSvc & fileSvc(tip::IFileSvc::instance());

   const tip::Table * edisp = fileSvc.readTable(m_filename, m_extname);

   tip::Table::ConstIterator it(edisp->begin());
   tip::ConstTableRecord & row(*it);

   row["rwidth"].get(m_rwidth);
   row["ltail"].get(m_ltail);

   row["energ_lo"].get(m_eBounds);
   std::vector<double> ehi;
   row["energ_hi"].get(ehi);
   m_eBounds.push_back(ehi.back());

   row["costh_hi"].get(m_muBounds);
   std::vector<double> mulo;
   row["costh_lo"].get(mulo);
   m_muBounds.push_back(mulo.back());

   delete edisp;
   
   computeCumulativeDists();
}

void Edisp::computeCumulativeDists() {
   size_t npts(100);
   if (s_xvals.empty()) {
      double xstep(2./(npts-1));
      for (size_t j = 0; j < npts; j++) {
         s_xvals.push_back(xstep*j - 1.);
      }
   }
   m_norms.clear();
   m_norms.reserve(m_rwidth.size());
   m_cumDists.clear();
   m_cumDists.reserve(m_rwidth.size());
   for (size_t i = 0; i < m_rwidth.size(); i++) {
      s_rwidth = m_rwidth.at(i);
      s_ltail = m_ltail.at(i);
      double partialInt;
      double err(1e-5);
      long ierr;
      std::vector<float> row;
      row.push_back(0);
      for (size_t j = 1; j < npts; j++) {
         dgaus8_(&edispIntegrand, &s_xvals.at(j-1), &s_xvals.at(j), 
                 &err, &partialInt, &ierr);
         row.push_back(row.back() + partialInt);
      }
      m_norms.push_back(1./row.back());
      for (size_t j = 0; j < npts; j++) {
         row.at(j) *= m_norms.back();
      }
      m_cumDists.push_back(row);
   }
}

auto Edisp::edispIntegrand(double * xx) -> double {
   double & x(*xx);
   if (x/s_rwidth > 7.) { // This leads to a 0.1% error in the integrand value
      double logf = s_ltail*std::log(1. + x) - x/s_rwidth;
      return std::exp(logf);
   }
   return std::pow(1. + x, s_ltail)/(1. + std::exp(x/s_rwidth));   
}

} // namespace dc2Response

/**
 * @file EfficiencyFactor.cxx
 * @brief Function that returns the IRF-dependent efficiency
 * corrections as a function of livetime fraction.  This is based on 
 * ratios of the livetime fraction-dependent effective area (P6_v6_diff) to
 * the livetime averaged effective area (P6_V3_DIFFUSE).  See
 * http://confluence.slac.stanford.edu/display/DC2/Efficiency+Correction+Parametrization+as+a+function+of+livetime+fraction+using+MC+(P6_v6_diff)
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

#include "irfUtil/IrfHdus.h"

#include "latResponse/FitsTable.h"

#include "EfficiencyFactor.h"

namespace latResponse {

EfficiencyFactor::
EfficiencyFactor() : m_havePars(false) {
   char * parfile = ::getenv("EFFICIENCY_PAR_FILE");
   if (parfile != 0) {
      readPars(parfile);
      m_havePars = true;
   }
}

EfficiencyFactor::
EfficiencyFactor(const std::string & parfile) : m_havePars(true) {
   readPars(parfile);
}

EfficiencyFactor::
EfficiencyFactor(const irfUtil::IrfHdus & aeff_hdus, size_t iepoch) 
   : m_havePars(true) {
   const std::string & fitsfile(aeff_hdus("EFFICIENCY_PARS")[iepoch].first);
   const std::string & extname(aeff_hdus("EFFICIENCY_PARS")[iepoch].second);
   readFitsFile(fitsfile, extname);
}

void EfficiencyFactor::
readPars(std::string parfile) {
   facilities::Util::expandEnvVar(&parfile);
   st_facilities::Util::file_ok(parfile);

   if (st_facilities::Util::isFitsFile(parfile)) {
      readFitsFile(parfile);
      return;
   }

   std::vector<std::string> lines;
   st_facilities::Util::readLines(parfile, lines, "#", true);

   std::vector< std::vector<double> > parVectors;
   for (size_t i(0); i < lines.size(); i++) {
      std::vector<std::string> tokens;
      bool clear;
      facilities::Util::stringTokenize(lines.at(i), " ", tokens, clear=true);
      std::vector<double> pars;
      for (size_t j(0); j < tokens.size(); j++) {
         pars.push_back(std::atof(tokens.at(j).c_str()));
      }
      parVectors.push_back(pars);
   }

   m_p0 = EfficiencyParameter(parVectors.at(0));
   m_p1 = EfficiencyParameter(parVectors.at(1));
}

void EfficiencyFactor::readFitsFile(const std::string & fitsfile,
                                    const std::string & extname) {
   const tip::Table * table = 
      tip::IFileSvc::instance().readTable(fitsfile, extname);

   long nrows;
   table->getHeader()["NAXIS2"].get(nrows);

   bool all_zeros(true);
   std::vector< std::vector<double> > parVectors;
   for (size_t i(0); i < static_cast<unsigned long>(nrows); i++) {
      std::vector<double> fltValues;
      FitsTable::getVectorData(table, "EFFICIENCY_PARS", fltValues, i);
      for (size_t j(0); j < fltValues.size(); j++) {
         if (fltValues.at(j) != 0) {
            all_zeros = false;
         }
      }
      std::vector<double> values(fltValues.size(), 0);
      std::copy(fltValues.begin(), fltValues.end(), values.begin());
      parVectors.push_back(values);
   }
   if (all_zeros) {
      delete table;
      m_havePars = false;
      return;
   }
   m_p0 = EfficiencyParameter(parVectors.at(0));
   m_p1 = EfficiencyParameter(parVectors.at(1));
   delete table;
}

auto EfficiencyFactor::operator()(double energy, double met) const -> double {
   if (!m_havePars || m_start.empty()) {
      return 1;
   }
   double ltfrac;

// Find the interval corresponedng to the desired met.  Intervals may
// not be uniform, so must do a search.  
   double tmin(m_start.front());
   double tmax(m_stop.back());
   double tol(0);
   if (met < tmin - tol || met > tmax + tol) {
      std::ostringstream message;
      message << "Requested MET of " << met << " "
              << "lies outside the range of valid times in the "
              << "pointing/livetime history: " 
              << tmin << " to " << tmax << "MET s";
      throw std::runtime_error(message.str());
   }
   std::vector<double>::const_iterator it 
      = std::upper_bound(m_start.begin(), m_start.end(), met);
   size_t indx = it - m_start.begin() - 1;

   if (m_start.at(indx) <= met && met <= m_stop.at(indx)) {
      ltfrac = m_livetimefrac.at(indx);
   }

   double value = IEfficiencyFactor::value(energy, ltfrac);
   return value;
}

auto EfficiencyFactor::value(double energy, double livetimefrac,
                               bool front, double met) const -> double {
   /// front is no longer used with more general event_types in Pass 8
   /// but we retain it to keep the interface backwards-compatible.
   (void)(front);  
   (void)(met);
   if (!m_havePars) {
      return 1;
   }
   double logE(std::log10(energy));
   double value(m_p0(logE)*livetimefrac + m_p1(logE));
   return value;
}

void EfficiencyFactor::
getLivetimeFactors(double energy, double & factor1, double & factor2,
                   double met) const {
   (void)(met);
   if (!m_havePars) {
      factor1 = 1;
      factor2 = 0;
      return;
   }
   double logE(std::log10(energy));
   factor1 = m_p1(logE);
   factor2 = m_p0(logE);
}

} // namespace latResponse

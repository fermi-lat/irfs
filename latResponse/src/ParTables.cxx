/**
 * @file ParTables.cxx
 * @brief Implementation for class to manage tables of parameters for
 * IRF tables.
 *
 * @author J. Chiang
 *
 * $Header$
 */

#include <sstream>
#include <stdexcept>

#include "tip/IFileSvc.h"
#include "tip/Table.h"

#include "latResponse/ParTables.h"

namespace latResponse {

ParTables::ParTables(const std::string & fitsfile,
                     const std::string & extname,
                     size_t nrow) {
   tip::IFileSvc & fileSvc(tip::IFileSvc::instance());
   const tip::Table * table(fileSvc.readTable(fitsfile, extname));
   const std::vector<std::string> & validFields(table->getValidFields());

   // The first four columns *must* be "ENERG_LO", "ENERG_HI", "CTHETA_LO",
   // "CTHETA_HI", in that order.
   const char * boundsName[] = {"energ_lo", "energ_hi", 
                                "ctheta_lo", "ctheta_hi"};
   for (size_t i(0); i < 4; i++) {
      if (validFields.at(i) != boundsName[i]) {
         delete table;
         std::ostringstream message;
         message << "latResponse::ParTables::ParTables: "
                 << "invalid header in " << fitsfile << "  "
                 << validFields.at(i) << "  " << i;
         throw std::runtime_error(message.str());
      }
   }

   // Read in the table values for the remaining columns.
   for (size_t i(4); i < validFields.size(); i++) {
      const std::string & tablename(validFields.at(i));
      m_parNames.push_back(tablename);
      m_parIndices[i-4] = tablename;
      m_parTables.insert(
         std::map<std::string, FitsTable>::
         value_type(tablename, FitsTable(fitsfile, extname, tablename, nrow)));
   }
   delete table;
}

auto ParTables::operator[](const std::string & parName) const -> const FitsTable & {
   std::map<std::string, FitsTable>::const_iterator table =
      m_parTables.find(parName);
   if (table == m_parTables.end()) {
      throw std::runtime_error("latResponse::ParTables::operator[]: "
                               "table name not found.");
   }
   return table->second;
}

void ParTables::getPars(double loge, double costh, double * pars,
                        bool interpolate) const {
   for (size_t i(0); i < m_parNames.size(); i++) {
      pars[i] = operator[](m_parNames.at(i)).value(loge, costh, interpolate);
   }
}

void ParTables::getParVector(const std::string & parName,
                             std::vector<double> & pars) const {
   operator[](parName).getValues(pars);
}

void ParTables::
getCornerPars(double logE, double costh, double & tt, double & uu, 
              std::vector<double> & cornerEnergies,
              std::vector<std::vector<double> > & parVectors) const {
// Create parameter vectors for each of the 4 corners of the cell.
   parVectors.clear();
   for (size_t j(0); j < 4; j++) {
      parVectors.push_back(std::vector<double>());
   }
// Loop over parameter names, retrieve the values for each corner, and
// append them to their respective corner vectors.
   for (size_t i(0); i < m_parNames.size(); i++) {
      std::vector<double> pars(4);
      operator[](m_parNames.at(i)).getCornerPars(logE, costh, tt, uu, 
                                                 cornerEnergies, pars);
      for (size_t j(0); j < pars.size(); j++) {
         parVectors.at(j).push_back(pars[j]);
      }
   }
}

void ParTables::
getPars(size_t ilogE, size_t icosth, std::vector<double> & pars) const {
   pars.clear();
   for (size_t i(0); i < m_parNames.size(); i++) {
      const FitsTable & fitsTable(operator[](m_parNames.at(i)));
      pars.push_back(fitsTable.getPar(ilogE, icosth));
   }                                                       
}

void ParTables::
setPars(size_t ilogE, size_t icosth, const std::vector<double> & pars) {
   for (size_t i(0); i < m_parNames.size(); i++) {
      std::map<std::string, FitsTable>::iterator table 
         = m_parTables.find(m_parNames.at(i));
      table->second.setPar(ilogE, icosth, pars[i]);
   }
}

auto ParTables::params(size_t indx) -> std::vector<double> { 

  std::map<size_t, std::string>::const_iterator itr = m_parIndices.find(indx);
  if(itr == m_parIndices.end())
    throw std::runtime_error("Parameter index out of range.");

  std::string parName = itr->second;

  return m_parTables[parName].values();
}

void ParTables::setParams(size_t indx, const std::vector<double>& params) {
  
  std::map<size_t, std::string>::const_iterator itr = m_parIndices.find(indx);
  if(itr == m_parIndices.end())
    throw std::runtime_error("Parameter index out of range.");

  std::string parName = itr->second;
  m_parTables[parName].setValues(params);
}


} // namespace latResponse

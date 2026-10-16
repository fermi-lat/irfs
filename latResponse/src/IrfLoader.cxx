/**
 * @brief Implementation for post-handoff review IRFs
 * @author J. Chiang
 *
 * $Header$
 */

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <vector>

#include "facilities/commonUtilities.h"
#include "facilities/Util.h"

#include "tip/IFileSvc.h"
#include "tip/Table.h"

#include "st_stream/StreamFormatter.h"

#include "st_facilities/Env.h"
#include "st_facilities/FitsUtil.h"
#include "st_facilities/Util.h"

#include "irfInterface/IrfRegistry.h"
#include "irfInterface/Irfs.h"
#include "irfInterface/IrfsFactory.h"

#include "irfUtil/EventTypeMapper.h"
#include "irfUtil/HdCaldb.h"
#include "irfUtil/IrfHdus.h"
#include "irfUtil/Util.h"

#include "latResponse/IrfLoader.h"

#include "latResponse/Aeff.h"
#include "AeffEpochDep.h"
#include "CaldbDate.h"
#include "Edisp.h"
#include "Edisp2.h"
#include "latResponse/Edisp3.h"
#include "EdispEpochDep.h"
#include "EfficiencyFactor.h"
#include "EfficiencyFactorEpochDep.h"
#include "Irfs.h"
#include "Psf.h"
#include "Psf2.h"
#include "latResponse/Psf3.h"
#include "PsfEpochDep.h"

namespace {
   typedef std::map<std::string, std::pair<unsigned int, std::string> > 
   EventTypeMapping_t;
}

namespace latResponse {

bool IrfLoader::s_interpolate_edisp(true);

IrfLoader::IrfLoader() 
   : m_hdcaldb(new irfUtil::HdCaldb("GLAST", "LAT")) {
   read_caldb_indx();
   readCustomIrfNames();
}

IrfLoader::~IrfLoader() {
   delete m_hdcaldb;
}

void IrfLoader::registerEventClasses() const {
   for (size_t i(0); i < m_caldbNames.size(); i++) {
      try {
         registerEventClasses(m_caldbNames.at(i));
      } catch (std::runtime_error & eObj) {
         // do nothing
      }
   }
}

void IrfLoader::registerEventClasses(const std::string & irfName) const {
   irfInterface::IrfRegistry & registry(irfInterface::IrfRegistry::instance());

   irfUtil::EventTypeMapper & evMapper(irfUtil::EventTypeMapper::instance());

   const EventTypeMapping_t & event_type_mapping(evMapper.mapping(irfName));

   std::vector<std::string> partitions;
   evMapper.getPartitions(irfName, partitions);

   for (size_t ip(0); ip < partitions.size(); ip++) {
      std::vector<std::string> classNames;
      EventTypeMapping_t::const_iterator it(event_type_mapping.begin());
      for ( ; it != event_type_mapping.end(); ++it) {
         if (it->second.second == partitions[ip]) {
            classNames.push_back(irfName + "::" + it->first);
            registry.registerEventClass(classNames.back(),classNames.back());
         }
      }
      std::ostringstream partition_id;
      if (partitions[ip] != "none") {
         partition_id << " (" << partitions[ip] << ")";
      }
      if (!classNames.empty()) {
         registry.registerEventClasses(irfName + partition_id.str(), 
                                       classNames);
      }
   }
}

void IrfLoader::loadIrfs() const {
   for (size_t i(0); i < m_caldbNames.size(); i++) {
      try {
         loadIrfs(m_caldbNames[i]);
      } catch (std::runtime_error & eObj) {
         // do nothing.
      }
   }
}

void IrfLoader::loadIrfs(const std::string & irfName) const {
   if (std::find(m_caldbNames.begin(), m_caldbNames.end(), irfName) 
       == m_caldbNames.end()) {
      throw std::runtime_error("IRF " + irfName + "not found in CALDB");
   }
   const EventTypeMapping_t & event_type_mapping =
      irfUtil::EventTypeMapper::instance().mapping(irfName);
   EventTypeMapping_t::const_iterator it(event_type_mapping.begin());
   for ( ; it != event_type_mapping.end(); ++it) {
      addIrfs(irfName, it->first);
   }
}

void IrfLoader::addIrfs(const std::string & irf_name, 
                        const std::string & event_type) const {
   irfInterface::IrfsFactory * myFactory(irfInterface::IrfsFactory::instance());
   const std::vector<std::string> & irfNames(myFactory->irfNames());

// Check if this set of IRFs already exists.
   if (std::count(irfNames.begin(), irfNames.end(), irf_name+"::"+event_type)) {
      return;
   }

   Irfs * irfs = new Irfs(irf_name, event_type);
   myFactory->addIrfs(irf_name + "::" + event_type, irfs);
}

auto 
IrfLoader::aeff(const irfUtil::IrfHdus & aeff_hdus) -> irfInterface::IAeff * {
   if (aeff_hdus.numEpochs() == 1) {
      return new Aeff(aeff_hdus, 0);
   } else {
      AeffEpochDep * my_aeff(new AeffEpochDep());
      const std::vector<std::pair<std::string, std::string> > &
         filename_hdu_pairs(aeff_hdus("EFF_AREA"));
      for (size_t j(0); j < aeff_hdus.numEpochs(); j++) {
         double epoch_start(EpochDep::epochStart(filename_hdu_pairs[j].first,
                                                 filename_hdu_pairs[j].second));
         irfInterface::IAeff * aeff_component(new Aeff(aeff_hdus, j));
         my_aeff->addAeff(*aeff_component, epoch_start);
         delete aeff_component;
      }
      return my_aeff;
   }
   return 0;
}

auto 
IrfLoader::psf(const irfUtil::IrfHdus & psf_hdus) -> irfInterface::IPsf * {
   if (psf_hdus.numEpochs() == 1) {
      return psf(psf_hdus, 0);
   } else {
      PsfEpochDep * my_psf(new PsfEpochDep());
      const std::vector<std::pair<std::string, std::string> > &
         filename_hdu_pairs(psf_hdus("RPSF"));
      for (size_t j(0); j < psf_hdus.numEpochs(); j++) {
         double epoch_start(EpochDep::epochStart(filename_hdu_pairs[j].first,
                                                 filename_hdu_pairs[j].second));
         irfInterface::IPsf * psf_component(psf(psf_hdus, j));
         my_psf->addPsf(*psf_component, epoch_start);
         delete psf_component;
      }
      return my_psf;
   }
   return 0;
}

auto IrfLoader::psf(const irfUtil::IrfHdus & psf_hdus, 
                                    size_t iepoch) -> irfInterface::IPsf * {
   // Return the single epoch Psf.
   const std::vector<std::pair<std::string, std::string> > &
      filename_hdu_pairs(psf_hdus("RPSF"));
   const std::string & psf_file(filename_hdu_pairs[iepoch].first);
   const std::string & extname(filename_hdu_pairs[iepoch].second);
   switch (psfVersion(psf_file, extname)) {
   case 1:
      return new Psf(psf_file, psf_hdus.convType()==0, extname);
      break;
   case 2:
      return new Psf2(psf_file, psf_hdus.convType()==0, extname);
      break;
   case 3:
      return new Psf3(psf_hdus, iepoch);
      break;
   default:
      throw std::runtime_error("PSF version not found.");
   }
   return 0;
}

auto 
IrfLoader::edisp(const irfUtil::IrfHdus & edisp_hdus) -> irfInterface::IEdisp * {
   if (edisp_hdus.numEpochs() == 1) {
      return edisp(edisp_hdus, 0);
   } else {
      EdispEpochDep * my_edisp(new EdispEpochDep());
      const std::vector<std::pair<std::string, std::string> > &
         filename_hdu_pairs(edisp_hdus("EDISP"));
      for (size_t j(0); j < edisp_hdus.numEpochs(); j++) {
         double epoch_start(EpochDep::epochStart(filename_hdu_pairs[j].first,
                                                 filename_hdu_pairs[j].second));
         irfInterface::IEdisp * edisp_component(edisp(edisp_hdus, j));
         my_edisp->addEdisp(*edisp_component, epoch_start);
         delete edisp_component;
      }
      return my_edisp;
   }
   return 0;
}

auto 
IrfLoader::edisp(const irfUtil::IrfHdus & edisp_hdus, size_t iepoch) -> irfInterface::IEdisp * {
   const std::vector<std::pair<std::string, std::string> > &
      filename_hdu_pairs(edisp_hdus("EDISP"));
   const std::string & edisp_file(filename_hdu_pairs[iepoch].first);
   const std::string & extname(filename_hdu_pairs[iepoch].second);
   switch (edispVersion(edisp_file, extname)) {
   case 1:
      return new Edisp(edisp_file, extname);
      break;
   case 2:
      return new Edisp2(edisp_file, extname);
      break;
   case 3:
      return new Edisp3(edisp_hdus, iepoch);
      break;
   default:
      throw std::runtime_error("EDISP version not found.");
   }
   return 0;
}

auto
IrfLoader::efficiency_factor(const irfUtil::IrfHdus & aeff_hdus) -> irfInterface::IEfficiencyFactor * {
   irfInterface::IEfficiencyFactor * my_eff(0);
   if (aeff_hdus.numEpochs() == 1) {
      my_eff = new EfficiencyFactor(aeff_hdus, 0);
   } else {
      my_eff = new EfficiencyFactorEpochDep();
      const std::vector<std::pair<std::string, std::string> > &
         filename_hdu_pairs(aeff_hdus("EFFICIENCY_PARS"));
      for (size_t j(0); j < aeff_hdus.numEpochs(); j++) {
         double epoch_start(EpochDep::epochStart(filename_hdu_pairs[j].first,
                                                 filename_hdu_pairs[j].second));
         irfInterface::IEfficiencyFactor * 
            eff_component(new EfficiencyFactor(aeff_hdus, j));
         dynamic_cast<EfficiencyFactorEpochDep *>(my_eff)->add(*eff_component,
                                                               epoch_start);
         delete eff_component;
      }
   }
   return my_eff;
}

void IrfLoader::readCustomIrfNames() {
   char * custom_irf_dir(::getenv("CUSTOM_IRF_DIR"));
   if (!custom_irf_dir) {
      return;
   }
   m_customIrfDir = custom_irf_dir;

   char * custom_irf_names(::getenv("CUSTOM_IRF_NAMES"));

   if (!custom_irf_names) {
      throw std::runtime_error("CUSTOM_IRF_NAMES env var not set." );
   }

   facilities::Util::stringTokenize(custom_irf_names, ", ", m_customIrfNames);

//    std::cout << "Adding custom IRFs: " << std::endl;
//    for (size_t i(0); i < m_customIrfNames.size(); i++) {
//       std::cout << m_customIrfNames.at(i) << std::endl;
//    }
}

void IrfLoader::loadCustomIrfs() const {
   // const std::string & irfDir(m_customIrfDir);

   // irfInterface::IrfsFactory * myFactory(irfInterface::IrfsFactory::instance());
   // const std::vector<std::string> & irfNames(myFactory->irfNames());

   // for (size_t i(0); i < m_customIrfNames.size(); i++) {
   //    const std::string & irfName(m_customIrfNames.at(i));
   //    if (!std::count(irfNames.begin(), irfNames.end(), irfName + "::FRONT")) {
   //       std::string section("front");
   //       std::string aeff_file = st_facilities::Env
   //          ::appendFileName(irfDir, "aeff_"+irfName+"_"+section+".fits");
   //       std::string psf_file = st_facilities::Env
   //          ::appendFileName(irfDir, "psf_"+irfName+"_"+section+".fits");
   //       std::string edisp_file = st_facilities::Env
   //          ::appendFileName(irfDir,"edisp_"+irfName+"_"+section+".fits");
   //       addIrfs(aeff_file, psf_file, edisp_file, 0, irfName);

   //       section = "back";
   //       aeff_file = st_facilities::Env
   //          ::appendFileName(irfDir, "aeff_"+irfName+"_"+section+".fits");
   //       psf_file = st_facilities::Env
   //          ::appendFileName(irfDir, "psf_"+irfName+"_"+section+".fits");
   //       edisp_file = st_facilities::Env
   //          ::appendFileName(irfDir,"edisp_"+irfName+"_"+section+".fits");
   //       addIrfs(aeff_file, psf_file, edisp_file, 1, irfName);
   //    }
   // }
}

void IrfLoader::read_caldb_indx() {
   /// Read the caldb.indx file and extract the IRF names from the
   /// CAL_CBD column.
   m_caldbNames.clear();
   char * caldb_path = ::getenv("CALDB");
   if (!caldb_path) {
      throw std::runtime_error("CALDB env var not set");
   }
   std::string caldb_indx;
   find_cif(caldb_indx);
   caldb_indx = facilities::commonUtilities::joinPath(caldb_path, caldb_indx);

   tip::IFileSvc & fileSvc(tip::IFileSvc::instance());
   const tip::Table * table = fileSvc.readTable(caldb_indx, "CIF");
   
   CaldbDate cutoff_date("2007-01-01");

   tip::Table::ConstIterator it(table->begin());
   tip::ConstTableRecord & row(*it);
   for ( ; it != table->end(); ++it) {
      std::string cal_date;
      row["cal_date"].get(cal_date);
      CaldbDate caldbDate(cal_date);
      if (caldbDate > cutoff_date) {
         std::vector<std::string> cal_cbd;
         row["cal_cbd"].get(cal_cbd);
         std::vector<std::string> tokens;
         facilities::Util::stringTokenize(cal_cbd.front(), "()", tokens);
         const std::string & caldbName(tokens.at(1));
         if (!std::count(m_caldbNames.begin(), m_caldbNames.end(), caldbName)) {
            m_caldbNames.push_back(caldbName);
         }
      }
   }
   delete table;
}

void IrfLoader::find_cif(std::string & caldb_indx) const {
   /// This returns the path to the caldb.indx file starting from
   /// $CALDB.  For LAT data, this path is typically
   /// "data/glast/lat/caldb.indx".  This information is ascertained
   /// from the $CALDBCONFIG file.
   char * caldb_config = ::getenv("CALDBCONFIG");
   if (!caldb_config) {
      throw std::runtime_error("CALDBCONFIG env var not set");
   }
   std::vector<std::string> lines;
   st_facilities::Util::readLines(caldb_config, lines);
   for (size_t i(0); i < lines.size(); i++) {
      std::vector<std::string> tokens;
      facilities::Util::stringTokenize(lines.at(i), " \t", tokens);
      if (tokens.at(0) == "GLAST" && tokens.at(1) == "LAT") {
         caldb_indx = facilities::commonUtilities::joinPath(tokens.at(3),
                                                            tokens.at(4));
         return;
      }
   }
   throw std::runtime_error("GLAST LAT not found in caldb.config file");
}

auto IrfLoader::edispVersion(const std::string & fitsfile, 
                            const std::string & extname) -> int {
   const tip::Table * table = 
      tip::IFileSvc::instance().readTable(fitsfile, extname);
   int version(1);
   try {
      table->getHeader()["EDISPVER"].get(version);
   } catch (tip::TipException & eObj) {
      /// EDISPVER keyword is (probably) missing, so assume default version
   }
   delete table;
   return version;
}

auto IrfLoader::psfVersion(const std::string & fitsfile, 
                          const std::string & extname) -> int {
   const tip::Table * table = 
      tip::IFileSvc::instance().readTable(fitsfile, extname);
   int version(1);
   try {
      table->getHeader()["PSFVER"].get(version);
   } catch (tip::TipException & eObj) {
      /// PSFVER keyword is (probably) missing, so assume default version
   }
   delete table;
   return version;
}

} // namespace latResponse

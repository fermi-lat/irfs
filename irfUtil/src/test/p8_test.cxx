#include <iostream>
#include <string>
#include <vector>

#include "tip/IFileSvc.h"
#include "tip/Table.h"

#include "irfUtil/HdCaldb.h"
#include "irfUtil/Util.h"
#include "irfUtil/IrfHdus.h"

void print_hdus(const irfUtil::IrfHdus & my_hdus) {
   std::cout << "irfID: " << my_hdus.bitPos() << std::endl;
   const std::vector<std::string> & cnames(my_hdus.cnames());
   for (size_t i(0); i < cnames.size(); i++) {
      const std::vector<std::pair<std::string, std::string> > 
         & fh_pairs(my_hdus(cnames[i]));
      for (size_t j(0); j < fh_pairs.size(); j++) {
         std::cout << cnames[i] << "  " 
                   << fh_pairs[j].first << "  "
                   << fh_pairs[j].second << std::endl;
      }
   }
}

auto main() -> int {
//   std::string detname("PSF0");
   std::string detname("EDISP1");
   std::string irfName("P8_SOURCE_V5");
   std::vector<std::string> respnames;
   respnames.push_back("EFF_AREA");
   respnames.push_back("PHI_DEP");
   respnames.push_back("EFFICIENCY_PARS");
   respnames.push_back("EDISP");
   respnames.push_back("EDISP_SCALING");
   respnames.push_back("RPSF");
   respnames.push_back("PSF_SCALING");
   respnames.push_back("FISHEYE_CORR");
   irfUtil::HdCaldb hdcaldb;
   for (size_t j(0); j < respnames.size(); j++) {
      std::vector<std::string> files;
      std::vector<int> hdus;
      hdcaldb.getFiles(files, hdus, detname, respnames[j], irfName);
      std::cout << detname << "  "
                << respnames[j] << std::endl;
      for (size_t i(0); i < files.size(); i++) {
         std::cout << "  " << i << "  "
                   << files[i] << "  "
                   << hdus[i] << std::endl;
      }
   }

/// Exercise the operator() interface
   irfUtil::HdCaldb foo;
   std::pair<std::string, int> calfile;
   for (size_t j(0); j < respnames.size(); j++) {
      calfile = foo(detname, respnames[j], "VERSION.eq." + irfName,
                    "NONE", "2008-08-04", "00:00:00");
      std::cout << detname << "  "
                << respnames[j] << std::endl;
      std::cout << "  " << calfile.first << "  "
                << calfile.second << std::endl;
      std::cout << std::endl;
   }

/// Exercise the getCaldbFile function
   std::string irf_file;
   long hdu;
   for (size_t j(0); j < respnames.size(); j++) {
      irfUtil::Util::getCaldbFile(detname, respnames[j], irfName,
                                  irf_file, hdu, "GLAST", "LAT", "NONE", 
                                  "2009-08-04", "00:00:00");
      std::cout << irf_file << "  " << hdu << std::endl;
   }

// /// Exercise get_event_[class,type]_mapping
//    std::map<std::string, unsigned int> mapping;
//    irfUtil::Util::get_event_class_mapping(mapping);
//    for (std::map<std::string, unsigned int>::const_iterator it(mapping.begin());
//         it != mapping.end(); ++it) {
//       std::cout << it->first << "  " << it->second << std::endl;
//    }

//    std::map<std::string, std::pair<unsigned int, std::string> > evtype_mapping;

//    std::cout << "\nP7REP_SOURCE_V10 event types:\n";
//    irfUtil::Util::get_event_type_mapping("P7REP_SOURCE_V10", evtype_mapping);
//    for (std::map<std::string, std::pair<unsigned int, std::string> >::const_iterator 
//            it(evtype_mapping.begin()); it != evtype_mapping.end(); ++it) {
//       std::cout << it->first << "  " << it->second.first << std::endl;
//    }

//    std::cout << "\nP8_SOURCE_V5 event types:\n";
//    irfUtil::Util::get_event_type_mapping("P8_SOURCE_V5", evtype_mapping);
//    for (std::map<std::string, std::pair<unsigned int, std::string> >::const_iterator 
//            it(evtype_mapping.begin()); it != evtype_mapping.end(); ++it) {
//       std::cout << it->first << "  " << it->second.first << std::endl;
//    }

// /// Check that tip can read by extension number.
   
//    std::vector<std::string> files;
//    std::vector<int> hdus;

//    hdcaldb.getFiles(files, hdus, "PSF0", "EFFICIENCY_PARS", "P8_CLEAN_V5");
//    std::ostringstream extname;
//    extname << hdus[0];
//    const tip::Table * 
//       table(tip::IFileSvc::instance().readTable(files[0], extname.str()));
                            
//    tip::Table::ConstIterator it(table->begin());
//    tip::ConstTableRecord & row(*it);
//    for ( ; it != table->end(); ++it) {
//       std::vector<float> my_values;
//       row["EFFICIENCY_PARAMS"].get(my_values);
//       for (size_t i(0); i < my_values.size(); i++) {
//          std::cout << my_values[i] << "  ";
//       }
//       std::cout << std::endl;
//    }
//    delete table;

//    table = tip::IFileSvc::instance().readTable(files[0], 
//                                                "EFFICIENCY_PARAMS_PSF0");
//    it = table->begin();
//    for ( ; it != table->end(); ++it) {
//       std::vector<float> my_values;
//       row["EFFICIENCY_PARAMS"].get(my_values);
//       for (size_t i(0); i < my_values.size(); i++) {
//          std::cout << my_values[i] << "  ";
//       }
//       std::cout << std::endl;
//    }
//    delete table;

// /// Test IrfHdus class
//    // irfUtil::IrfHdus empty_hdus("Un_irf", "Un_event_type",
//    //                             irfUtil::IrfHdus::s_aeff_cnames);
//    // std::cout << "Contents of empty_hdus: " 
//    //           << empty_hdus.numEpochs() << std::endl;

//    std::cout << "\nP8_SOURCE_V5, PSF0, aeff:" << std::endl;
//    print_hdus(irfUtil::IrfHdus::aeff("P8_SOURCE_V5", "PSF0"));

//    std::cout << "\nP8_CLEAN_V5, EDISP1, psf:" << std::endl;
//    try {
//       print_hdus(irfUtil::IrfHdus::psf("P8_CLEAN_V5", "EDISP1"));
//    } catch (...) {
//    }

//    std::cout << "\nP8_ULTRACLEAN_V5, FRONT, edisp:" << std::endl;
//    try {
//       print_hdus(irfUtil::IrfHdus::edisp("P8_ULTRACLEAN_V5", "FRONT"));
//    } catch (...) {
//    }
   return 0;
}

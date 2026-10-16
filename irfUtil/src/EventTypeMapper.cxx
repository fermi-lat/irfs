/**
 * @file EventTypeMapper.cxx
 * @brief Implementation for a Singleton class to manage data on the
 * mapping between IRF name and event types.
 * @author J. Chiang
 * 
 * $Header$
 */

#include <stdexcept>

#include "facilities/commonUtilities.h"
#include "facilities/Util.h"

#include "st_facilities/Environment.h"

#include "tip/IFileSvc.h"
#include "tip/Table.h"

#include "irfUtil/EventTypeMapper.h"
#include "irfUtil/HdCaldb.h"
#include "irfUtil/Util.h"

namespace irfUtil {

EventTypeMapper * EventTypeMapper::s_instance(0);

EventTypeMapper::EventTypeMapper() {
   // Find the irf_index.fits file.
   std::string sub_path;
   Util::joinPaths("data glast lat bcf irf_index.fits", sub_path);
   std::string irf_index = facilities::commonUtilities::joinPath(
      st_facilities::Environment::getEnv("CALDB"), sub_path);
   
   // Open the bitmask_mapping extension.
   const tip::Table * evclass_map
      = tip::IFileSvc::instance().readTable(irf_index, "BITMASK_MAPPING");

   // Open event_type_mapping extension.
   const tip::Table * evtype_map
      = tip::IFileSvc::instance().readTable(irf_index, "EVENT_TYPE_MAPPING");

   // Loop over event_class names in BITMASK_MAPPING extension.
   tip::Table::ConstIterator evclass(evclass_map->begin());
   tip::ConstTableRecord & evclass_row(*evclass);
   for ( ; evclass != evclass_map->end(); ++evclass) {
      std::string event_class;
      unsigned int allowed_evtypes(0);
      evclass_row["event_class"].get(event_class);
      evclass_row["event_types"].get(allowed_evtypes);

      // Create mapping of event_type name to allowed bit positions
      // and mapping of partition name to full bit-masks for each
      // event_class.
      EventTypeMapper::EvTypeMapping_t mapping;
      std::map<std::string, unsigned int> full_bitmasks;

      // Loop over event_types.
      tip::Table::ConstIterator evtype;
      tip::ConstTableRecord & evtype_row(*evtype);
      for (evtype = evtype_map->begin(); evtype != evtype_map->end();
           ++evtype) {
         std::string event_type;
         int bitpos;
         std::string partition;
         evtype_row["event_type"].get(event_type);
         evtype_row["bitposition"].get(bitpos);
         evtype_row["event_type_partition"].get(partition);
         if ((allowed_evtypes >> bitpos) & 1) {
            mapping[event_type] = std::make_pair(bitpos, partition);
         }
         if (full_bitmasks.find(partition) == full_bitmasks.end()) {
            full_bitmasks[partition] = 0;
         }
      }
      m_mappings[event_class] = mapping;

      for (EventTypeMapper::EvTypeMapping_t::const_iterator 
              itor(mapping.begin()); itor != mapping.end(); ++itor) {
         full_bitmasks[itor->second.second] += (1 << itor->second.first);
      }
      m_full_bitmasks[event_class] = full_bitmasks;
   }
   delete evtype_map;
   delete evclass_map;
}

auto EventTypeMapper::instance() -> EventTypeMapper & {
   if (s_instance == 0) {
      s_instance = new EventTypeMapper();
   }
   return *s_instance;
}

void EventTypeMapper::deleteInstance() {
   delete s_instance;
   s_instance = 0;
}

auto
EventTypeMapper::bitPos(const std::string & irfName,
                        const std::string & event_type) const -> const std::pair<unsigned int, std::string> & {
   const EvTypeMapping_t & my_mapping = mapping(irfName);
   EvTypeMapping_t::const_iterator it = my_mapping.find(event_type);
   if (it == my_mapping.end()) {
      throw std::runtime_error(irfName + " and " + event_type + " not found.");
   }
   return it->second;
}

auto
EventTypeMapper::mapping(const std::string & irfName) const -> const EventTypeMapper::EvTypeMapping_t & {
   std::map<std::string, EvTypeMapping_t>::const_iterator it 
      = m_mappings.find(irfName);
   if (it == m_mappings.end()) {
      throw std::runtime_error("IRFs named " + irfName + 
                               " not found in CALDB.");
   }
   return it->second;
}

auto
EventTypeMapper::full_bitmasks_by_partition(const std::string & irfName) const -> const std::map<std::string, unsigned int> & {
   std::map<std::string, std::map<std::string, unsigned int> >::const_iterator
      it = m_full_bitmasks.find(irfName);
   if (it == m_full_bitmasks.end()) {
      throw std::runtime_error("IRFs named " + irfName + 
                               " not found in CALDB.");
   }
   return it->second;
}

void EventTypeMapper::
getPartitions(const std::string & irfName,
              std::vector<std::string> & partitions) const {
   const std::map<std::string, unsigned int> & full_bitmasks
      = full_bitmasks_by_partition(irfName);
   for (std::map<std::string, unsigned int>::const_iterator it 
           = full_bitmasks.begin(); it != full_bitmasks.end(); ++it) {
      partitions.push_back(it->first);
   }
}

} // namespace irfUtil

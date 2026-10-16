/**
 * $Header$
 */

#include <cstdlib>

#include <iostream>

#include "st_facilities/Environment.h"

#include "irfInterface/IrfsFactory.h"
#include "dc1Response/loadIrfs.h"

auto main() -> int {
   st_facilities::Environment::instance();
   if (!std::getenv("CALDB")) {
      std::cout << "CALDB not set, exiting." << std::endl;
      std::exit(0);
   }
   dc1Response::load_irfs();
   irfInterface::IrfsFactory * myFactory 
      = irfInterface::IrfsFactory::instance();
   irfInterface::Irfs * myIrfs = myFactory->create("DC1::Front");
   delete myIrfs;
   myIrfs = myFactory->create("DC1::Back");
   delete myIrfs;

// try loading again
   dc1Response::load_irfs();
   myIrfs = myFactory->create("DC1::Front");
   delete myIrfs;
   myIrfs = myFactory->create("DC1::Back");
   delete myIrfs;

// access upper limit on aeff
   std::cout << "Maximum effective areas:\n";
   dc1Response::load_irfs();
   myIrfs = myFactory->create("DC1::Front");
   std::cout << "Front: " << myIrfs->aeff()->upperLimit() << std::endl;
   delete myIrfs;
   myIrfs = myFactory->create("DC1::Back");
   std::cout << "Back: " << myIrfs->aeff()->upperLimit() << std::endl;
   delete myIrfs;
}


/**
 * @file test.cxx
 * @brief Test program for dc2Response
 * @author J. Chiang
 *
 * $Header$
 */

#include <cmath>
#include <cstdlib>
#include <iostream>
#include "st_facilities/Environment.h"
#include "irfInterface/IrfsFactory.h"
#include "dc2Response/loadIrfs.h"

void checkPsf(irfInterface::Irfs * irf, float energy, float theta, float phi) {

// Check limit when the separation goes to 0, for which a special
// formula is needed
   double value0=irf->psf()->value(0.,energy,theta,phi);
   double value1=irf->psf()->value(0.00001,energy,theta,phi);
   std::cout<<"values at 0 and close to it: "<<value0<<" "<<value1<<std::endl;

// Check that the integral is correctly normalized to 1
   float max = 50; //degrees
   float step=0.001;//degrees
   float Integral=0.;
   int run=0;
   while(run*step<max){
      float delta=run*step;
      float psfval = irf->psf()->value(delta,energy,theta,phi);
      Integral += psfval*std::sin(delta*M_PI/180.)*step;
      run++;
   }
   Integral*=M_PI/180.;//remormalize step in the integration above
   Integral*=2*M_PI;//No dependance on phi
   std::cout<<"Integral from 0 to 50 degrees : "<<Integral<<std::endl;
}

auto main(int argc, char* argv[]) -> int {
   st_facilities::Environment::instance();
   if (!std::getenv("CALDB")) {
      std::cout << "CALDB not set, exiting." << std::endl;
      std::exit(0);
   }

   float energy = 1000.;            //MeV
   float theta = 0.*M_PI/180.;
   float phi = 0.*M_PI/180.;  
   switch (argc){
   case 2 :
      energy = std::atof(argv[1]);
      break;
   case 3 :
      energy = std::atof(argv[1]);
      theta = std::atof(argv[2])*M_PI/180.;
      break;
   case 4 :
      energy = std::atof(argv[1]);
      theta = std::atof(argv[2])*M_PI/180.;
      phi = std::atof(argv[3])*M_PI/180.;
      break;
   default :
      break;
   }

   try {
      dc2Response::load_irfs();
      irfInterface::IrfsFactory * myFactory 
         = irfInterface::IrfsFactory::instance();
      irfInterface::Irfs * myIrfs(0);
      myIrfs = myFactory->create("DC2::FrontA");
      double totalAeff(0);
      std::cout << "FrontA: " << myIrfs->aeff()->upperLimit() << std::endl;
      totalAeff += myIrfs->aeff()->upperLimit();
      checkPsf(myIrfs,energy,theta,phi);
      delete myIrfs;
      myIrfs = myFactory->create("DC2::BackA");
      std::cout << "BackA: " << myIrfs->aeff()->upperLimit() << std::endl;
      totalAeff += myIrfs->aeff()->upperLimit();
      checkPsf(myIrfs,energy,theta,phi);
      delete myIrfs;
      myIrfs = myFactory->create("DC2::FrontB");
      std::cout << "FrontB: " << myIrfs->aeff()->upperLimit() << std::endl;
      totalAeff += myIrfs->aeff()->upperLimit();
      checkPsf(myIrfs,energy,theta,phi);
      delete myIrfs;
      myIrfs = myFactory->create("DC2::BackB");
      std::cout << "BackB: " << myIrfs->aeff()->upperLimit() << std::endl;
      totalAeff += myIrfs->aeff()->upperLimit();
      checkPsf(myIrfs,energy,theta,phi);
      delete myIrfs;
      std::cout << "Total: " << totalAeff << std::endl;
      std::exit(0);
   } catch (std::exception & eObj) {
      std::cout << eObj.what() << std::endl;
   } catch (std::string & what) {
      std::cout << what << std::endl;
   }
   std::exit(1);
}

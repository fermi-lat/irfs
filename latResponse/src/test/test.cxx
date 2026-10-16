/**
 * @file test.cxx
 * @brief Test program for latResponse package
 * @author J. Chiang
 *
* $Header$
*/

#ifdef TRAP_FPE
#include <fenv.h>
#endif

#include <cmath>
#include <cstdlib>

#include <iostream>
#include <iomanip>
#include <stdexcept>

#include <cppunit/ui/text/TextTestRunner.h>
#include <cppunit/extensions/HelperMacros.h>

#include "astro/SkyDir.h"

#include "facilities/commonUtilities.h"

#include "st_facilities/Environment.h"

#include "irfInterface/IrfsFactory.h"
#include "irfInterface/AcceptanceCone.h"

#include "latResponse/IrfLoader.h"

#include "latResponse/Aeff.h"
#include "latResponse/Psf3.h"
#include "Edisp2.h"
#include "EfficiencyFactor.h"
#include "AeffEpochDep.h"
#include "PsfEpochDep.h"
#include "EdispEpochDep.h"
#include "EfficiencyFactorEpochDep.h"

using facilities::commonUtilities;

namespace {
   auto getEnv(const std::string & envVarName) -> std::string {
      char * envvar(::getenv(envVarName.c_str()));
      if (envvar == 0) {
         throw std::runtime_error("Please set the " + envVarName 
                                  + " environment variable.");
      }
      return envvar;
   }
}

class LatResponseTests : public CppUnit::TestFixture {

   CPPUNIT_TEST_SUITE(LatResponseTests);

   CPPUNIT_TEST(irf_assignment);

   CPPUNIT_TEST(psf_zero_separation);
   CPPUNIT_TEST(psf_normalization);
   CPPUNIT_TEST(psf_roi_integral);

   CPPUNIT_TEST(edisp_normalization);
   CPPUNIT_TEST(edisp_sampling);

   CPPUNIT_TEST(epochDep_tests);

   CPPUNIT_TEST_SUITE_END();

public:

   void setUp();
   void tearDown();

   void irf_assignment();

   void psf_zero_separation();
   void psf_normalization();
   void psf_roi_integral();

   void edisp_normalization();
   void edisp_sampling();

   void epochDep_tests();

private:

   irfInterface::IrfsFactory * m_irfsFactory;
   std::vector<std::string> m_irfNames;

};

void LatResponseTests::setUp() {
   m_irfsFactory = irfInterface::IrfsFactory::instance();
   m_irfsFactory->getIrfsNames(m_irfNames);
//    m_irfNames.push_back("P6_V8_DIFFUSE::FRONT");
//    m_irfNames.push_back("P6_V8_DIFFUSE::BACK");
//    m_irfNames.push_back("P6_v8_diff::FRONT");
//    m_irfNames.push_back("P6_v8_diff::BACK");
}

void LatResponseTests::tearDown() {
   m_irfNames.clear();
}

void LatResponseTests::irf_assignment() {
// This assumes FRONT IRFs have even irfIDs and BACK IRFs odd irfIDs.
   for (std::vector<std::string>::const_iterator name(m_irfNames.begin());
        name != m_irfNames.end(); ++name) {
      if (name->find("FRONT") != std::string::npos) {
         irfInterface::Irfs * my_irfs(m_irfsFactory->create(*name));
         CPPUNIT_ASSERT(my_irfs->irfID() % 2 == 0);
         delete my_irfs;
      }
      if (name->find("BACK") != std::string::npos) {
         irfInterface::Irfs * my_irfs(m_irfsFactory->create(*name));
         CPPUNIT_ASSERT(my_irfs->irfID() % 2 == 1);
         delete my_irfs;
      }
   }
}

void LatResponseTests::psf_zero_separation() {
   double energy(1e3);
   double theta(0);
   double phi(0);

   double delta_sep(1e-4);
   double tol(1e-3);

   double time(239846401.);  // 08Aug2008 00:00:00
   
   for (std::vector<std::string>::const_iterator name = m_irfNames.begin();
        name != m_irfNames.end(); ++name) {
      irfInterface::Irfs * myIrfs(m_irfsFactory->create(*name));
      const irfInterface::IPsf & psf(*myIrfs->psf());
      for (theta = 0; theta < 70; theta += 5.) {
         double value0(psf.value(0, energy, theta, phi, time));
         double value1(psf.value(delta_sep, energy, theta, phi, time));
         if (std::fabs((value0 - value1)/value0) >= tol) {
            std::cout << *name << ": \n";
            std::cout << "theta = " << theta << ": "
                      << "psf value at 0 = " << value0 << "  "
                      << "psf value at 1e-4 = " << value1 << std::endl;
         }
         CPPUNIT_ASSERT(std::fabs((value0 - value1)/value0) < tol);
      }
      delete myIrfs;
   }
}

auto integrate_psf(const irfInterface::IPsf & psf,
		     const irfInterface::AcceptanceCone & cone,
		     double energy, double theta, double phi, double time, 
		     int nstep) -> double {

  const astro::SkyDir& roidir = cone.center();

  double x0 = roidir.ra();
  double y0 = roidir.dec();
  double r0 = ::sqrt(x0*x0+y0*y0);

  const double radius = cone.radius();
  const double dstep = 2*radius/double(nstep);

  double sum = 0;
  for (int i = 0; i < nstep; i++) {
    for (int j = 0; j < nstep; j++) {

      double x = x0 - radius + (i+0.5)*dstep;
      double y = y0 - radius + (j+0.5)*dstep;
      double r = ::sqrt(x*x+y*y);
      double deltar = ::sqrt(std::pow(x-x0,2) + std::pow(y-y0,2));

      if( deltar > cone.radius() )
	continue;

      double psf_val = psf.value(r, energy, theta, phi, time);
      sum += psf_val;
    }
  }

  sum *= std::pow(dstep,2)*std::pow(M_PI/180.,2);
  return sum;
}

void LatResponseTests::psf_roi_integral() {

  double phi(0);
  double time(239846401.);  // 08Aug2008 00:00:00
  astro::SkyDir srcdir(0, 0.0);
  astro::SkyDir roidir(0, 1.5);

  double tol(1e-1);
  
  std::vector<double> psi;

  std::vector<irfInterface::AcceptanceCone *> cones;

  irfInterface::AcceptanceCone * cone = new irfInterface::AcceptanceCone(roidir,1.3);
  cones.push_back(cone);

  int nstep = 100;

  std::vector<double> energies;
  double emin(5e1);
  double emax(5e5);
  size_t nee(13);
  double dee(std::log(emax/emin)/(nee-1));
  for (size_t i = 0; i < nee; i++) {
    energies.push_back(emin*std::exp(i*dee));
  }

  std::vector<double> thetas(3);
  thetas[0] = 10.0;
  thetas[1] = 30.0;
  thetas[2] = 60.0;

  std::vector<std::string> irfNames;
  irfNames.push_back("P8R2_SOURCE_V6::FRONT");
  irfNames.push_back("P8R2_SOURCE_V6::BACK");
  irfNames.push_back("P8R2_SOURCE_V6::PSF0");
  irfNames.push_back("P8R2_SOURCE_V6::PSF1");
  irfNames.push_back("P8R2_SOURCE_V6::PSF2");
  irfNames.push_back("P8R2_SOURCE_V6::PSF3");
  irfNames.push_back("P8R2_SOURCE_V6::EDISP0");
  irfNames.push_back("P8R2_SOURCE_V6::EDISP1");
  irfNames.push_back("P8R2_SOURCE_V6::EDISP2");
  irfNames.push_back("P8R2_SOURCE_V6::EDISP3");

  bool integralFailures(false);

  for (std::vector<std::string>::const_iterator name = irfNames.begin();
       name != irfNames.end(); ++name) {

    irfInterface::Irfs * myIrfs(m_irfsFactory->create(*name));
    irfInterface::IPsf & psf(*myIrfs->psf());
    
    for (std::vector<double>::const_iterator energy = energies.begin();
	 energy != energies.end(); ++energy) {

      for (std::vector<double>::const_iterator theta = thetas.begin();
	   theta != thetas.end(); ++theta) {

	double angIntRef = integrate_psf(psf,*cone,*energy,*theta,phi,time,nstep);
	double angInt(psf.angularIntegral(*energy, srcdir, *theta, phi, cones, time));
	double fdiff = (angInt-angIntRef)/angIntRef;
      
	if(fabs(fdiff) > tol) {
	  std::cout 
	    << std::setw(25) << *name
	    << std::setw(15) << *energy 
	    << std::setw(15) << *theta 
	    << std::setw(15) << angIntRef 
	    << std::setw(15) << angInt 
	    << std::setw(15) << fdiff 
	    << std::endl;
	  integralFailures = true;
	}
      }
    }
  }

  delete cone;

  CPPUNIT_ASSERT(!integralFailures);
}

void LatResponseTests::psf_normalization() {
   double phi(0);
   double time(239846401.);  // 08Aug2008 00:00:00

   double tol(2e-2);

   std::vector<double> psi;
   double psimin(1e-4);
   double psimax(90);
   size_t npsi(1000);
   double dpsi(std::log(psimax/psimin)/(npsi-1));
   psi.push_back(0);
   for (size_t i = 0; i < npsi; i++) {
      psi.push_back(psimin*std::exp(i*dpsi));
   }

   std::vector<double> energies;
   double emin(30);
   double emax(5e5);
   size_t nee(100);
   double dee(std::log(emax/emin)/(nee-1));
   for (size_t i = 0; i < nee; i++) {
      energies.push_back(emin*std::exp(i*dee));
   }

   std::vector<double> thetas;
   double thmin(0);
   double thmax(65);
   size_t nth(20);
   double dth((thmax - thmin)/(nth-1));
   for (size_t i(0); i < nth; i++) {
      thetas.push_back(i*dth + thmin);
   }

   std::cout << "PSF integral values that fail "
             << int(tol*100) << "% tolerance: \n"
             << "energy  inclination  integral est.  angularIntegral\n";

   bool integralFailures(false);

   for (std::vector<std::string>::const_iterator name = m_irfNames.begin();
        name != m_irfNames.end(); ++name) {
      std::cout << *name << ": \n";
      if (*name == "P6_V7_TRANSIENT::FRONT" || 
          *name == "P6_V7_TRANSIENT::BACK" || 
          *name == "P7TRANSIENT_V6::BACK" ||
          *name == "P7REP_TRANSIENT_V10::BACK" ||
          *name == "P7REP_TRANSIENT_V15::BACK") { 
         continue;
      }
      irfInterface::Irfs * myIrfs(m_irfsFactory->create(*name));
      const irfInterface::IPsf & psf(*myIrfs->psf());
      
      size_t jjj(0);
      for (std::vector<double>::const_iterator energy = energies.begin();
           energy != energies.end(); ++energy) {
         for (std::vector<double>::const_iterator theta = thetas.begin();
              theta != thetas.end(); ++theta, jjj++) {
           
            std::vector<double> psf_values;
            for (size_t i = 0; i < psi.size(); i++) {
               psf_values.push_back(psf.value(psi[i], *energy, 
                                              *theta, phi, time));
            }
            
            double integral(0);
            for (size_t i = 0; i < psi.size() - 1; i++) {
               integral += ((psf_values[i]*std::sin(psi[i]*M_PI/180.) + 
                        psf_values[i+1]*std::sin(psi[i+1]*M_PI/180.))/2.)
                  *(psi[i+1] - psi[i])*M_PI/180.;
            }
            integral *= 2.*M_PI;
            double angInt(psf.angularIntegral(*energy, *theta, phi, psimax,
                                              time));
            if (std::fabs(integral - 1.) >= tol) { 
               std::cout << *energy  << "      " 
                         << *theta   << "           "
                         << integral << "        "
                         << angInt << std::endl;
               integralFailures = true;
            }
//             CPPUNIT_ASSERT(std::fabs(integral - 1.) < tol);
//             CPPUNIT_ASSERT(std::fabs(angInt - 1.) < tol);
         }
      }
      delete myIrfs;
   }
   CPPUNIT_ASSERT(!integralFailures);
}

void LatResponseTests::edisp_normalization() {
   double time(239846401.);  // 08Aug2008 00:00:00

   std::vector<double> energies;
   double emin(30);
   double emax(3e5);
   size_t nee(10);
   double dee(std::log(emax/emin)/(nee-1));
   for (size_t i = 1; i < nee; i++) {
      energies.push_back(emin*std::exp(i*dee));
   }
   /// Add specific values where Pass 8 has problems.
   energies.push_back(133.4);
   energies.push_back(421.7);
   energies.push_back(749.9);
   energies.push_back(2371.4);

   std::vector<double> thetas;
   double thmin(0);
   double thmax(60);
   size_t nth(6);
   double dth((thmax - thmin)/(nth-1));
   for (size_t i = 0; i < nth; i++) {
      thetas.push_back(i*dth + thmin);
   }
   /// Add specific values where Pass 8 has problems.
   thetas.push_back(49.46);
   thetas.push_back(63.26);
//   thetas.push_back(75.52);

   double phi(0);

   double tol(3e-2);

   bool integralFailures(false);
   std::cout << "Energy dispersion integral values that fail "
             << int(tol*100) << "% tolerance: \n"
             << "energy  inclination  integral  (integral - 1) \n";

   for (std::vector<std::string>::const_iterator name(m_irfNames.begin());
        name != m_irfNames.end(); ++name) {
      // if ((*name).find("P8_SOURCE") == std::string::npos) {
      //    continue;
      // }
      std::cout << *name << ": \n";
      irfInterface::Irfs * myIrfs(m_irfsFactory->create(*name));
      const irfInterface::IEdisp & edisp(*myIrfs->edisp());
      for (std::vector<double>::const_iterator energy(energies.begin());
           energy != energies.end(); ++energy) {
         double efactor(10);
         double elower(*energy/efactor);
         double eupper(*energy*efactor);
         for (std::vector<double>::const_iterator theta(thetas.begin());
              theta != thetas.end(); ++theta) {
            double integral(edisp.integral(elower, eupper, *energy,
                                           *theta, phi, time));
            if (std::fabs(integral - 1.) >= tol) {
               std::cout << *energy << "     "
                         << *theta  << "          "
                         << integral << "     "
                         << (integral - 1) << std::endl;
               integralFailures = true;
            }
//            CPPUNIT_ASSERT(std::fabs(integral - 1.) < tol);
         }
      }
      delete myIrfs;
   }
   CPPUNIT_ASSERT(!integralFailures);
}

void LatResponseTests::edisp_sampling() {
// Just exercise the energy dispersion sampling to look for floating
// point exceptions.
   std::vector<double> energies;
   double emin(1e5); //(10);
   double emax(3e5);
   size_t nee(10);
   double dee(std::log(emax/emin)/(nee-1));
   for (size_t i = 0; i < nee; i++) {
      energies.push_back(emin*std::exp(i*dee));
   }

   std::vector<double> thetas;
   double thmin(0);
   double thmax(60);
   size_t nth(6);
   double dth((thmax - thmin)/(nth-1));
   for (size_t i = 0; i < nth; i++) {
      thetas.push_back(i*dth + thmin);
   }

   size_t nsamp(20);

   astro::SkyDir zAxis(0, 0);
   astro::SkyDir xAxis(90, 0);
   double time(239846401.);  // 08Aug2008 00:00:00

   for (std::vector<std::string>::const_iterator name(m_irfNames.begin());
        name != m_irfNames.end(); ++name) {
      irfInterface::Irfs * myIrfs(m_irfsFactory->create(*name));
      const irfInterface::IEdisp & edisp(*myIrfs->edisp());
      for (std::vector<double>::const_iterator theta(thetas.begin());
           theta != thetas.end(); ++theta) {
         astro::SkyDir srcDir(0, *theta);
         for (std::vector<double>::const_iterator energy(energies.begin());
              energy != energies.end(); ++energy) {
            for (size_t j = 0; j < nsamp; j++) {
               try {
                  edisp.appEnergy(*energy, srcDir, zAxis, xAxis, time);
               } catch(const std::exception & e) {
                  std::cerr << "caught sampling error, "
                            << (*energy)<<", "<<srcDir.dec() <<", " 
                            << zAxis.dec() << ", " <<xAxis.dec()  << std::endl;
                  throw;
               }
            }
         }
      }
   }
}

void LatResponseTests::epochDep_tests() {
   double energy(100);
   double theta(20);
   double phi(30);
   double met0(252460801.);  // 2009-01-01 00:00:00 (in epoch 0)
   double met1(283996802.);  // 2010-01-01 00:00:00 (in epoch 1)

   std::string dataPath(st_facilities::Environment::dataPath("latResponse"));

   // Effective area
   latResponse::AeffEpochDep aeff;

   std::string aeff0 = commonUtilities::joinPath(dataPath,
                                                 "aeff_epoch_0.fits");
   std::string aeff1 = commonUtilities::joinPath(dataPath,
                                                 "aeff_epoch_1.fits");

   latResponse::Aeff aeff_epoch0(aeff0);
   latResponse::Aeff aeff_epoch1(aeff1);

   aeff.addAeff(aeff_epoch0,
                latResponse::EpochDep::epochStart(aeff0, "EFFECTIVE AREA"));
   aeff.addAeff(aeff_epoch1,
                latResponse::EpochDep::epochStart(aeff1, "EFFECTIVE AREA"));

   CPPUNIT_ASSERT(aeff.value(energy, theta, phi, met0) == 
                  aeff_epoch0.value(energy, theta, phi, met0));
   CPPUNIT_ASSERT(aeff.value(energy, theta, phi, met0) != 
                  aeff_epoch1.value(energy, theta, phi, met0));
   CPPUNIT_ASSERT(aeff.value(energy, theta, phi, met1) == 
                  aeff_epoch1.value(energy, theta, phi, met1));
   CPPUNIT_ASSERT(aeff.value(energy, theta, phi, met1) != 
                  aeff_epoch0.value(energy, theta, phi, met1));

   // PSF
   latResponse::PsfEpochDep psf;

   std::string psf0(commonUtilities::joinPath(dataPath, "psf_epoch_0.fits"));
   std::string psf1(commonUtilities::joinPath(dataPath, "psf_epoch_1.fits"));
                                                
   latResponse::Psf3 psf_epoch0(psf0);
   latResponse::Psf3 psf_epoch1(psf1);
   psf.addPsf(psf_epoch0, latResponse::EpochDep::epochStart(psf0, "RPSF"));
   psf.addPsf(psf_epoch1, latResponse::EpochDep::epochStart(psf1, "RPSF"));

   double sep(15);

   CPPUNIT_ASSERT(psf.value(sep, energy, theta, phi, met0) == 
                  psf_epoch0.value(sep, energy, theta, phi, met0));
   CPPUNIT_ASSERT(psf.value(sep, energy, theta, phi, met0) != 
                  psf_epoch1.value(sep, energy, theta, phi, met0));
   CPPUNIT_ASSERT(psf.value(sep, energy, theta, phi, met1) == 
                  psf_epoch1.value(sep, energy, theta, phi, met1));
   CPPUNIT_ASSERT(psf.value(sep, energy, theta, phi, met1) != 
                  psf_epoch0.value(sep, energy, theta, phi, met1));

   // Energy dispersion
   std::string edisp0(commonUtilities::joinPath(dataPath,
                                                "edisp_epoch_0.fits"));
   std::string edisp1(commonUtilities::joinPath(dataPath,
                                                "edisp_epoch_1.fits"));
   latResponse::EdispEpochDep edisp;

   latResponse::Edisp2 edisp_epoch0(edisp0);
   latResponse::Edisp2 edisp_epoch1(edisp1);

   std::string extname("ENERGY DISPERSION");
   edisp.addEdisp(edisp_epoch0,
                  latResponse::EpochDep::epochStart(edisp0, extname));
   edisp.addEdisp(edisp_epoch1,
                  latResponse::EpochDep::epochStart(edisp1, extname));

   double measE(120.);

   CPPUNIT_ASSERT(edisp.value(measE, energy, theta, phi, met0) == 
                  edisp_epoch0.value(measE, energy, theta, phi, met0));
   CPPUNIT_ASSERT(edisp.value(measE, energy, theta, phi, met0) != 
                  edisp_epoch1.value(measE, energy, theta, phi, met0));
   CPPUNIT_ASSERT(edisp.value(measE, energy, theta, phi, met1) == 
                  edisp_epoch1.value(measE, energy, theta, phi, met1));
   CPPUNIT_ASSERT(edisp.value(measE, energy, theta, phi, met1) != 
                  edisp_epoch0.value(measE, energy, theta, phi, met1));

   // EfficiencyFactor
   latResponse::EfficiencyFactorEpochDep eff;
   latResponse::EfficiencyFactor eff_epoch0(aeff0);
   latResponse::EfficiencyFactor eff_epoch1(aeff1);
   eff.add(eff_epoch0, 
           latResponse::EpochDep::epochStart(aeff0, "EFFICIENCY_PARAMS"));
   eff.add(eff_epoch1, 
           latResponse::EpochDep::epochStart(aeff1, "EFFICIENCY_PARAMS"));
   CPPUNIT_ASSERT(eff(energy, met0) == eff_epoch0(energy, met0));
   CPPUNIT_ASSERT(eff(energy, met1) == eff_epoch1(energy, met1));
}

auto main(int iargc, char * argv[]) -> int {
#ifdef TRAP_FPE
// Add floating point exception traps.
   feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
#endif
   st_facilities::Environment::instance();
   latResponse::IrfLoader myLoader;
   myLoader.loadIrfs();
   
   if (iargc > 1 && std::string(argv[1]) == "-d") {
      LatResponseTests testObj;
      testObj.setUp();
      testObj.edisp_normalization();
      testObj.tearDown();
   } else {
      
      if (!std::getenv("CALDB")) {
         std::cout << "CALDB not set, exiting." << std::endl;
         std::exit(0);
      }
      
      CppUnit::TextTestRunner runner;
      runner.addTest(LatResponseTests::suite());
      bool result(runner.run());
      if (result) {
         return 0;
      } else {
         return 1;
      }
   }
}

/**
 * @file Psf3.h
 * @brief Psf3 class declaration.
 * @author J. Chiang
 *
 * $Header$
 */

#ifndef latResponse_Psf3_h
#define latResponse_Psf3_h

#include <array>
#include <string>
#include <vector>

#include "latResponse/PsfBase.h"

namespace irfUtil {
class IrfHdus;
}

namespace latResponse {

class PsfIntegralCache;

/**
 * @class Psf3
 *
 * @brief Revised PSF that is the sum of two King model functions.
 * See http://confluence.slac.stanford.edu/x/bADIAw.
 * In contrast to Psf2, this class interpolates the distributions
 * rather than the parameters.
 *
 */

class Psf3 : public PsfBase {

public:
  Psf3(const irfUtil::IrfHdus &psf_hdus, size_t iepoch, size_t nrow = 0);

  Psf3(const std::string &fitsfile, bool isFront = true,
       const std::string &extname = "RPSF", size_t nrow = 0);

  Psf3(const Psf3 &other);

  Psf3 &operator=(const Psf3 &rhs);

  virtual ~Psf3();

  /// A member function returning the point-spread function value.
  /// @param appDir Apparent (reconstructed) photon direction.
  /// @param energy True photon energy in MeV.
  /// @param srcDir True photon direction.
  /// @param scZAxis Spacecraft z-axis.
  /// @param scXAxis Spacecraft x-axis.
  /// @param time Photon arrival time (MET s)
  virtual double value(const astro::SkyDir &appDir, double energy,
                       const astro::SkyDir &srcDir,
                       const astro::SkyDir &scZAxis,
                       const astro::SkyDir &scXAxis, double time = 0) const;

  /// Return the psf as a function of instrument coordinates.
  /// @param separation Angle between apparent and true photon directions
  ///        (degrees).
  /// @param energy True photon energy (MeV).
  /// @param theta True photon inclination angle (degrees).
  /// @param phi True photon azimuthal angle measured wrt the instrument
  ///            X-axis (degrees).
  /// @param time Photon arrival time (MET s)
  virtual double value(double separation, double energy, double theta,
                       double phi, double time = 0) const;

  typedef std::vector<irfInterface::AcceptanceCone *> AcceptanceConeVector_t;

  /// Angular integral of the PSF over the intersection of acceptance
  /// cones.
  virtual double angularIntegral(double energy, const astro::SkyDir &srcDir,
                                 double theta, double phi,
                                 const AcceptanceConeVector_t &acceptanceCones,
                                 double time = 0);

  virtual double angularIntegral(double energy, double theta, double phi,
                                 double radius, double time = 0) const;

  virtual int nparams() const { return m_parVectors[0].size(); }
  virtual std::vector<double> params(size_t indx) const;

  /// Bin centers in energy
  virtual const std::vector<double> &energies() const { return m_energies; }

  /// Bin centers in log10(energy)
  virtual const std::vector<double> &logEnergies() const { return m_logEs; }

  /// Bin centers in cos(theta)
  virtual const std::vector<double> &costhetas() const { return m_cosths; }

  /// Bin centers in theta
  virtual const std::vector<double> &thetas() const { return m_thetas; }

  virtual irfInterface::IPsf *clone() { return new Psf3(*this); }

  void setParams(size_t indx, const std::vector<double> &params);

  static int findIndex(const std::vector<double> &xx, double x);

private:
  // PSF parameters, energy and cos(theta) bin defs.
  std::vector<double> m_logEs;
  std::vector<double> m_energies;
  std::vector<double> m_cosths;
  std::vector<double> m_thetas;
  std::vector<std::vector<double>> m_parVectors;

  PsfIntegralCache *m_integralCache;

  void readFits(const std::string &fitsfile,
                const std::string &extname = "RPSF", size_t nrow = 0);

  void normalize_pars(double radius = 90.);

  double evaluate(double energy, double sep, const double *pars) const;

  void getCornerPars(double energy, double theta, double &tt, double &uu,
                     std::array<double, 4> &cornerEnergies,
                     std::array<size_t, 4> &indx) const;

  double psf_base_integral(double energy, double radius,
                           const double *pars) const;

  double angularIntegral(double energy, double psi,
                         const std::vector<double> &pars);

  static void
  generateBoundaries(const std::vector<double> &x, const std::vector<double> &y,
                     const std::vector<double> &values,
                     std::vector<double> &xout, std::vector<double> &yout,
                     std::vector<double> &values_out, double xlo = 0,
                     double xhi = 10., double ylo = -1., double yhi = 1.);
};

} // namespace latResponse

#endif // latResponse_Psf3_h

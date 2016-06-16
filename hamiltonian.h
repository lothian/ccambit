#ifndef HAMILTONIAN_H
#define HAMILTONIAN_H

#include <boost/shared_ptr.hpp>
#include <libmints/mints.h>
#include <libtrans/integraltransform.h>

#include <ambit/tensor.h>

namespace psi { namespace ccambit {

using namespace ambit;

class Hamiltonian {
public:
  Hamiltonian(boost::shared_ptr<PSIO>, boost::shared_ptr<Wavefunction>, std::vector<boost::shared_ptr<MOSpace> >);
  virtual ~Hamiltonian();

protected:
  int nmo_;
  int nso_;
  ambit::Dimension nact_;
  int nfzc_;
  int nfzv_;
  double efzc_;

//  double **fock_;
  Tensor fock_;
  double ****ints_;
  double ****L_;

  friend class CCWfn;
  friend class HBAR;
  friend class CCLambda;
  friend class CCDensity;
  friend class CCPert;
}; // Hamiltonian

}} // psi::ccambit

#endif // HAMILTONIAN_H

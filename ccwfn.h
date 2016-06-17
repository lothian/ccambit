#ifndef CCWFN_H
#define CCWFN_H

#include "hamiltonian.h"
#include <libmints/mints.h>
#include <boost/shared_ptr.hpp>

#include <ambit/tensor.h>

namespace psi { namespace ccambit {

class CCWfn: public Wavefunction {
public:
  CCWfn(boost::shared_ptr<Wavefunction> reference, boost::shared_ptr<Hamiltonian> H, Options &options);
  virtual ~CCWfn();

protected:
  std::string wfn_;     // wfn type (CCSD, CCSD_T, etc.)
  double convergence_;  // conv. on RMS residual change between iterations
  int maxiter_;         // maximum number of iterations
  bool do_diis_;        // use DIIS algorithms?
  bool ooc_;            // Use out-of-core algorithms?
  std::string dertype_; // Gradient level -- needed only for (T) gradients

  int no_;  // Number of active occupied MOs
  int nv_;  // Number of active virtual MOs

  boost::shared_ptr<Hamiltonian> H_; // integrals and Fock matrix

  // Energy denominators
  Tensor D1_;
  Tensor D2_;

  // Ground-state T amplitudes
  Tensor t1_;      // Current T1
  Tensor t1old_;   // Previous iteration T1
  Tensor t2_;    // Current T2
  Tensor t2old_; // Previous iteration T2

  // DIIS-related vectors
  std::vector<double> t1diis_;
  std::vector<double> t2diis_;
  std::vector<double> t1err_;
  std::vector<double> t2err_;

  // Effective doubles 
  Tensor tau_;  // tau(ijab) = t2(ijab) + t1(ia) * t1(jb)
  Tensor ttau_; // ttau(ijab) = t2(ijab) + (1/2) t1(ia) * t1(jb)

  // Biorthogonal projection doubles
  Tensor t1s_;
  Tensor t2s_;
  
  // CCSD intermediates for amplitude equations (related to, but not the
  // same as corresponding HBAR quantities)
  Tensor Fvv_;     
  Tensor Foo_;     
  Tensor Fov_;     
  Tensor Woooo_; 
  Tensor Wovvo_; 
  Tensor Wovov_; 

  // Extra contributions for (T) gradients
  Tensor t3_; // only for in-core code
  Tensor l3_; // only for in-core code
  Tensor s1_;
  Tensor s2_;
  Tensor Doo_;
  Tensor Dvv_;
  Tensor Dov_;
  Tensor Gooov_;
  Tensor Gvvvo_;
  Tensor Goovv_;

//  double energy();
//  void build_tau();
//  void amp_save();
//  void build_F();
//  void build_W();
//  void build_t1();
//  void build_t2();
//  void build_diis_error();
//  void save_diis_vectors();
//  double t1norm();
//  double increment_amps();
//  void build_tstar();

//  double tcorr();
//  double tcorr_ooc();
//  double tcorr_ooc_TJL();
//  void tgrad();
//  void tgrad_ooc();

//  void t3_ijk(Tensor t3, int i, int j, int k, Tensor t2, Tensor fock, Tensor ints);
//  void W3_ijk(Tensor W3, int i, int j, int k, Tensor t2, Tensor ints);
//  void t3_abc(Tensor t3, int a, int b, int c, Tensor t2, Tensor fock, Tensor ints);
//  void M3_ijk(Tensor M3, int i, int j, int k, Tensor t2, Tensor fock, Tensor ints);
//  void M3_abc(Tensor M3, int a, int b, int c, Tensor t2, Tensor fock, Tensor ints);
//  void N3_ijk(Tensor N3, int i, int j, int k, Tensor t2, Tensor t1, Tensor fock, Tensor ints);
//  void N3_abc(Tensor N3, int a, int b, int c, Tensor t2, Tensor t1, Tensor fock, Tensor ints);

public:
//  double compute_energy();

  friend class HBAR;
  friend class CCLambda;
  friend class CCDensity;
  friend class CCPert;
}; // CCWfn

}} // psi::ugacc

#endif // CCWFN_H

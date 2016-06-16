#include "hamiltonian.h"
#include <libiwl/iwl.h>
#include <libmints/wavefunction.h>
#include <psifiles.h>
#include <libciomr/libciomr.h>
#include <libqt/qt.h>
#include <libtrans/integraltransform.h>
#include <libdpd/dpd.h>

#include <ambit/tensor.h>

#define ID(x) ints.DPD_ID(x)

namespace psi { namespace ccambit {

Hamiltonian::Hamiltonian(boost::shared_ptr<PSIO> psio, boost::shared_ptr<Wavefunction> ref, std::vector<boost::shared_ptr<MOSpace> > spaces)
{
  nmo_ = ref->nmo();
  nso_ = ref->nso();
  nfzc_ = ref->nfrzc();
  nfzv_ = 0;
  for(int i=0; i < ref->nirrep(); i++) 
    nfzv_ += ref->frzvpi()[i];
  nact_ = nmo_ - nfzc_ - nfzv_;
  size_t nact = nact_;

  int nact2 = nact_ * nact_;
  int nact3 = nact2 * nact_;

  int *mo_offset = init_int_array(ref->nirrep()); // Pitzer offsets
  for(int h=1; h < ref->nirrep(); h++) mo_offset[h] = mo_offset[h-1] + ref->nmopi()[h-1];

  int *map = init_int_array(nmo_); // Translates from Pitzer (including frozen docc) to QT
  reorder_qt((int *) ref->doccpi(), (int *) ref->soccpi(), (int *) ref->frzcpi(), (int *) ref->frzvpi(), 
             map, (int *) ref->nmopi(), ref->nirrep());

  Tensor fock_ = Tensor::build(CoreTensor, "Fock Matrix", {nact, nact});
  vector<double>& fockV = fock_.data();
  double *fockP = fockV.data();

  Tensor H = Tensor::build(CoreTensor, "Core Hamiltonian Matrix", {nact, nact});
  vector<double>& HV = H.data();
  double *HP = HV.data();

  // Prepare core Hamiltonian in MO basis in QT ordering

  // Prepare Fock matrix in MO basis in QT ordering
  SharedMatrix Fa = ref->Fa();
  SharedMatrix Hcore = ref->H();
  SharedMatrix Ca = ref->Ca();
  Fa->transform(Ca);
  Hcore->transform(Ca);
  for(int h=0; h < ref->nirrep(); h++) {
    int nmo = ref->nmopi()[h]; int nfv = ref->frzvpi()[h]; int nfc = ref->frzcpi()[h];
    for(int p=nfc; p < nmo-nfv; p++) {
      for(int q=nfc; q < nmo-nfv; q++) {
      int P = map[p+mo_offset[h]]; int Q = map[q+mo_offset[h]];
      fockP[(P-nfzc_)*nact_ + (Q-nfzc_)] = Fa->get(h,p,q);
      HP[(P-nfzc_)*nact_ + (Q-nfzc_)] = Hcore->get(h,p,q);
      }
    }
  }
  fock_.print(stdout);
  H.print(stdout);

  free(mo_offset);
  free(map);

  // Use reorder_qt() to generate a new mapping array w/o frozen core or virtual orbitals
  int *doccpi = init_int_array(ref->nirrep());
  int *nmopi = init_int_array(ref->nirrep());
  int *null = init_int_array(ref->nirrep());
  for(int h=0; h < ref->nirrep(); h++) {
    doccpi[h] = ref->doccpi()[h] - ref->frzcpi()[h];
    nmopi[h] = ref->nmopi()[h] - ref->frzcpi()[h] - ref->frzvpi()[h];
  }
  int *map2 = init_int_array(nact); // Translates from Pitzer (w/o frozen MOs) to QT
  reorder_qt(doccpi, (int *) ref->soccpi(), null, null, map2, nmopi, ref->nirrep());
  free(null); free(nmopi); free(doccpi);

  IntegralTransform ints(ref, spaces, IntegralTransform::Restricted, IntegralTransform::DPDOnly);
  ints.transform_tei(MOSpace::all, MOSpace::all, MOSpace::all, MOSpace::all);
  efzc_ = ints.get_frozen_core_energy();

  psio->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);
  dpd_set_default(ints.get_dpd_id());
  dpdbuf4 K;
  global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[A,A]"), ID("[A,A]"), ID("[A>=A]+"), ID("[A>=A]+"), 0, "MO Ints (AA|AA)");
  Tensor ints_ = Tensor::build(CoreTensor, "MO Two-Electron Integrals", {nact, nact, nact, nact});
  vector<double>& intsV = ints_.data();
  double *intsP = intsV.data();
  for(int h=0; h < ref->nirrep(); h++) {
    global_dpd_->buf4_mat_irrep_init(&K, h);
    global_dpd_->buf4_mat_irrep_rd(&K, h);
    for(int pq=0; pq < K.params->rowtot[h]; pq++) {
      int p = map2[ K.params->roworb[h][pq][0] ];
      int q = map2[ K.params->roworb[h][pq][1] ];
      for(int rs=0; rs < K.params->coltot[h]; rs++) {
        int r = map2[ K.params->colorb[h][rs][0] ];
        int s = map2[ K.params->colorb[h][rs][1] ];
        intsP[p*nact3 + r*nact2 + q*nact_ + s] = K.matrix[h][pq][rs];
      }
    }
    global_dpd_->buf4_mat_irrep_close(&K, h);
  }
  global_dpd_->buf4_close(&K);
  psio->close(PSIF_LIBTRANS_DPD, 1);

  ints_.print(stdout);

  // L(pqrs) = 2<pq|rs> - <pq|sr>  
  Tensor L_ = Tensor::build(CoreTensor, "2<pq|rs> - <pq|sr>", {nact, nact, nact, nact});
  L_("p,q,r,s") = 2.0 * ints_("p,q,r,s") - ints_("p,q,s,r");
  L_.print(stdout);
}

Hamiltonian::~Hamiltonian()
{
/*
  free_4d_array(ints_, nact_, nact_, nact_);
  free_4d_array(L_, nact_, nact_, nact_);
  free_block(fock_); 
*/
}

}} // namespace psi::ccambit

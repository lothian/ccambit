/*
 * @BEGIN LICENSE
 *
 * ccambit by Psi4 Developer, a plugin to:
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2016 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

#include <psi4-dec.h>
#include <libparallel/parallel.h>
#include <liboptions/liboptions.h>
#include <libmints/mints.h>
#include <libpsio/psio.hpp>
#include <ambit/tensor.h>

#include <libtrans/integraltransform.h>
#include <map>

#include "hamiltonian.h"
#include "ccwfn.h"

using namespace boost;
using namespace std;

namespace psi{ namespace ccambit {

extern "C"
int read_options(std::string name, Options& options)
{
  if (name == "CCAMBIT"|| options.read_globals()) {
    options.add_int("PRINT", 1);
    options.add_str("REFERENCE", "RHF");
    options.add_str("WFN", "CCSD");
    options.add_str("DERTYPE", "NONE");
    options.add_int("MAXITER", 100);
    options.add_bool("DIIS", true);
    options.add_double("R_CONVERGENCE", 1e-7);
    options.add_bool("OOC", false);
  }

  return true;
}

extern "C"
SharedWavefunction ccambit(SharedWavefunction ref, Options& options)
{
  outfile->Printf("\n");
  outfile->Printf("\t\t\t**************************\n");
  outfile->Printf("\t\t\t*                        *\n");
  outfile->Printf("\t\t\t*        CC-AMBIT        *\n");
  outfile->Printf("\t\t\t*                        *\n");
  outfile->Printf("\t\t\t**************************\n");
  outfile->Printf("\n");

  outfile->Printf("\tWave function  = %s\n", options.get_str("WFN").c_str());
  outfile->Printf("\tMaxiter        = %d\n", options.get_int("MAXITER"));
  outfile->Printf("\tConvergence    = %3.1e\n", options.get_double("R_CONVERGENCE"));
  outfile->Printf("\tDIIS           = %s\n", options.get_bool("DIIS") ?  "Yes" : "No");
  outfile->Printf("\tOut-of-core    = %s\n", options.get_bool("OOC") ?  "Yes" : "No");
  outfile->Printf("\tDertype        = %s\n", options.get_str("DERTYPE").c_str());

  // Error trapping – need closed-shell SCF in place
  if(!ref) throw PSIEXCEPTION("SCF has not been run yet!");
  if(options.get_str("REFERENCE") != "RHF")
    throw PSIEXCEPTION("Only for use with RHF references.");
  for(int h=0; h < ref->nirrep(); h++)
    if(ref->soccpi()[h]) throw PSIEXCEPTION("UGACC is for closed-shell systems only.");

  boost::shared_ptr<PSIO> psio(_default_psio_lib_);
  std::vector<boost::shared_ptr<MOSpace> > spaces;
  spaces.push_back(MOSpace::all);
  boost::shared_ptr<Hamiltonian> H(new Hamiltonian(psio, ref, spaces));

  boost::shared_ptr<CCWfn> cc(new CCWfn(ref, H, options));

  return ref;
}

}} // End namespaces


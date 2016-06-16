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

using namespace boost;

namespace psi{ namespace ccambit {

extern "C"
int read_options(std::string name, Options& options)
{
    if (name == "CCAMBIT"|| options.read_globals()) {
        /*- The amount of information printed to the output file -*/
        options.add_int("PRINT", 1);
    }

    return true;
}

extern "C"
SharedWavefunction ccambit(SharedWavefunction ref_wfn, Options& options)
{
    int print = options.get_int("PRINT");

    /* Your code goes here */
    using namespace ambit;
    Tensor A = Tensor::build(CoreTensor, "A", {5,5});
    A.print(stdout);

    // Typically you would build a new wavefunction and populate it with data
    return ref_wfn;
}

}} // End namespaces


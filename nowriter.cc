/*
 *@BEGIN LICENSE
 *
 * nowriter by Ward Poelmans, a plugin to:
 *
 * PSI4: an ab initio quantum chemistry software package
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
 *@END LICENSE
 */

#include "psi4-dec.h"
#include <libmints/mints.h>
#include "libciomr/libciomr.h"
#include <liboptions/liboptions.h>
#include "libchkpt/chkpt.h"
#include "libdpd/dpd.h"
#include "psifiles.h"
#include "libiwl/iwl.h"
#include "libmints/mints.h"
#include "libtrans/integraltransform.h"
#include "libiwl/iwl.hpp"
#include "libplugin/plugin.h"
#include <libmints/writer.h>
#include <algorithm>
#include "psifiles.h"

INIT_PLUGIN

using namespace boost;

namespace psi{ namespace nowriter {

extern "C" 
int read_options(std::string name, Options& options)
{
    if (name == "NOWRITER"|| options.read_globals()) {
        options.add_int("PRINT", 1);
        options.add_str_i("OUTPUT_FILENAME", "out.molden");
    }

    return true;
}

extern "C" 
PsiReturnType nowriter(Options& options)
{
    const int print = options.get_int("PRINT");
    std::string filename = options.get_str("OUTPUT_FILENAME");
    // PSI always converts to all upper cases for the moment...
    std::transform(filename.begin(), filename.end(), filename.begin(), ::toupper);

    shared_ptr<PSIO> psio(_default_psio_lib_);

    std::vector<shared_ptr<MOSpace> > spaces;

    std::string reference = options.get_str("REFERENCE");

    if(reference != "RHF")
    {
        outfile->Printf("No RHF reference, that's not gonna fly with this plugin...\n");
        return Failure;
    }

    boost::shared_ptr<Wavefunction> wfn = Process::environment.wavefunction();

    // all of this is for the ints.alpha_corr_to_pitzer()
    // No idea where else to get this?
    // Pitzer means: order by irrep, then orbital energy
    spaces.push_back(MOSpace::all);
    IntegralTransform ints(wfn,
            spaces,
            reference == "RHF" ? IntegralTransform::Restricted : IntegralTransform::Unrestricted,
            IntegralTransform::DPDOnly,
            IntegralTransform::QTOrder,
            IntegralTransform::None,
            false);
    dpd_set_default(ints.get_dpd_id());
    ints.set_print(print);
    ints.initialize();
    const int *aPitzer = ints.alpha_corr_to_pitzer();
    const int *bPitzer = ints.beta_corr_to_pitzer(); // This is the same as the alpha array for RHF references

    char **labels = Process::environment.molecule()->irrep_labels();
    int nso       = Process::environment.wavefunction()->nso();
    int nmo       = Process::environment.wavefunction()->nmo();
    int nIrreps   = Process::environment.wavefunction()->nirrep();
    int *orbspi   = Process::environment.wavefunction()->nmopi();
    int *frzcpi   = Process::environment.wavefunction()->frzcpi();
    int *frzvpi   = Process::environment.wavefunction()->frzvpi();
    int *clsdpi   = Process::environment.wavefunction()->doccpi();
    double eNuc   = Process::environment.molecule()->nuclear_repulsion_energy();
    SharedMatrix Ca = Process::environment.wavefunction()->Ca(); 

    outfile->Printf( "\n\n\t Irrep  frzcpi doccpi frzvpi  sopi\n");
    for(int h = 0; h < nIrreps; ++h)
        outfile->Printf( "\t  %3s    %3d    %3d    %3d    %3d\n",
                labels[h], frzcpi[h], clsdpi[h], frzvpi[h], orbspi[h]);

    int nTriMo = nmo * (nmo + 1) / 2;
    int nTriSo = nso * (nso + 1) / 2;

    if(print > 4)
    {
        outfile->Printf("Pitzer ordening:\n");
        for(int n = 0; n < nmo; ++n)
            outfile->Printf( "\tAlpha: %3d -> %3d   Beta: %3d -> %3d\n", n, aPitzer[n], n, bPitzer[n]);
    }

    Matrix moOpdm("MO OPDM", nIrreps, orbspi, orbspi);

    psio->open(PSIF_MO_OPDM, PSIO_OPEN_OLD);

    // I hope this is always the correct dimension...
    double **tempOPDM = block_matrix(nmo, nmo);
    psio->read_entry(PSIF_MO_OPDM, "MO-basis OPDM", (char *) tempOPDM[0], sizeof(double)*nmo*nmo);

    double *tempMo = new double[nTriMo];
    for(int p=0;p<nmo;++p)
        for(int q=0;q<=p;++q)
        {
            int P = aPitzer[p];
            int Q = aPitzer[q];
            size_t PQ = INDEX(P,Q);
            tempMo[PQ] = 0.5 * (tempOPDM[p][q] + tempOPDM[q][p]);
        }
    free_block(tempOPDM);
    psio->close(PSIF_MO_OPDM, 1);

    moOpdm.set(tempMo);
    delete[] tempMo;

    if(print > 4)
    {
        outfile->Printf( "The MO basis OPDM, in Pitzer order\n");
        moOpdm.print();
    }


    SharedMatrix eigvec = Matrix::create("Naturals", moOpdm.rowspi(), moOpdm.colspi());
    SharedVector eigv = Vector::create("Occ numbers", moOpdm.rowspi());
    moOpdm.diagonalize(eigvec, eigv, descending);

    if(print > 4)
    {
        outfile->Printf( "Eigs:\n");
        eigv->print("outfile");
        outfile->Printf( "Eig vec:\n");
        eigvec->print();
    }

    Ca->print();

    SharedMatrix soOpdm = moOpdm.clone();
    soOpdm->set_name("SO OPDM");
    soOpdm->zero();

    soOpdm->gemm(false, false, 1.0, Ca, eigvec, 0.0);

    soOpdm->print();

    boost::shared_ptr<psi::MoldenWriter> molden(new psi::MoldenWriter(wfn));

    SharedVector epsilon_zero(wfn->epsilon_a()->clone());
    epsilon_zero->zero();

    SharedVector eigv2(eigv->clone());
    eigv2->zero();

    molden->write(filename, soOpdm, soOpdm, epsilon_zero, epsilon_zero, eigv, eigv2);

    return Success;
}

}} // End namespaces


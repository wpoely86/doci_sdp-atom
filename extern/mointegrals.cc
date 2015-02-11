#include <libplugin/plugin.h>
#include <psi4-dec.h>
#include <libdpd/dpd.h>
#include <psifiles.h>
#include <libpsio/psio.hpp>
#include <libiwl/iwl.hpp>
#include <libtrans/integraltransform.h>
#include <libmints/wavefunction.h>
#include <libmints/mints.h>
#include <libciomr/libciomr.h>
#include <liboptions/liboptions.h>
#include <libchkpt/chkpt.h>
#include <hdf5.h>

#include <stdlib.h>
#include <iostream>

#include "Irreps.h"
#include "Hamiltonian.h"
#include "Options.h"

using namespace std;

// This allows us to be lazy in getting the spaces in DPD calls
#define ID(x) ints.DPD_ID(x)

// macro to help check return status of HDF5 functions
#define HDF5_STATUS_CHECK(status) { \
    if(status < 0) \
    fprintf(outfile, "%s:%d: Error with HDF5. status code=%d\n", __FILE__, __LINE__, status); \
} 

INIT_PLUGIN

namespace psi{ namespace mointegrals{

extern "C" int
read_options(std::string name, Options &options)
{
    if (name == "MOINTEGRALS"|| options.read_globals()) {

        /*- The amount of information printed
            to the output file -*/
        options.add_bool("PRINT_INTEGRALS", true);
        /*- Whether to compute two-electron integrals -*/
        options.add_bool("DO_TEI", true);
        // save to a HDF5 file
        options.add_bool("SAVEHDF5", false);
        options.add_str_i("HDF5_FILENAME", "integrals.h5");
    }

    return true;
}


extern "C" PsiReturnType
mointegrals(Options &options)
{
    /*
     * This plugin shows a simple way of obtaining MO basis integrals, directly from a DPD buffer.  It is also
     * possible to generate integrals with labels (IWL) formatted files, but that's not shown here.
     */
    bool print = options.get_bool("PRINT_INTEGRALS");
    bool doTei = options.get_bool("DO_TEI");
    bool savehdf5 = options.get_bool("SAVEHDF5");
    std::string filename = options.get_str("HDF5_FILENAME");
    boost::algorithm::to_lower(filename);
   
    // Grab the global (default) PSIO object, for file I/O
    boost::shared_ptr<PSIO> psio(_default_psio_lib_);

    // Now we want the reference (SCF) wavefunction
    boost::shared_ptr<Wavefunction> wfn = Process::environment.wavefunction();
    if(!wfn) throw PSIEXCEPTION("SCF has not been run yet!");
    
    /*MoldenWriter mollie(wfn);
    mollie.write("infoForMolden.txt");*/ //Please call the MoldenWriter from the psi4 input file from now on.
    
    // Quickly check that there are no open shell orbitals here...
    int nirrep  = wfn->nirrep();
    
    // For now, we'll just transform for closed shells and generate all integrals.  For more elaborate use of the
    // LibTrans object, check out the plugin_mp2 example.
    std::vector<boost::shared_ptr<MOSpace> > spaces;
    spaces.push_back(MOSpace::all);
    IntegralTransform ints(wfn, spaces, IntegralTransform::Restricted);
    ints.transform_tei(MOSpace::all, MOSpace::all, MOSpace::all, MOSpace::all);
    // Use the IntegralTransform object's DPD instance, for convenience
    dpd_set_default(ints.get_dpd_id());
    
    //Readin the MO OEI in moOei & print everything
    int nmo       = Process::environment.wavefunction()->nmo();
    int nIrreps   = Process::environment.wavefunction()->nirrep();
    int *orbspi   = Process::environment.wavefunction()->nmopi();
    int *docc     = Process::environment.wavefunction()->doccpi();
    int *socc     = Process::environment.wavefunction()->soccpi();

    shared_ptr<Molecule> molecule = Process::environment.molecule();
    int nelectrons = 0;
    for(int i=0;i<molecule->natom();i++)
        nelectrons += molecule->true_atomic_number(i);

    nelectrons -= molecule->molecular_charge();

    fprintf(outfile, "****  Molecular Integrals For CheMPS Start Here \n");
    
    std::string SymmLabel =  Process::environment.molecule()->sym_label();
    fprintf(outfile, "Symmetry Label = ");
    fprintf(outfile, SymmLabel.c_str());
    fprintf(outfile, " \n");
    fprintf(outfile, "Nirreps = %1d \n", nirrep);
    double NuclRepulsion =  Process::environment.molecule()->nuclear_repulsion_energy();
    fprintf(outfile, "Nuclear Repulsion Energy = %16.48f \n", NuclRepulsion);
    fprintf(outfile, "Number Of Molecular Orbitals = %2d \n", nmo);
    fprintf(outfile, "Irreps Of Molecular Orbitals = \n");
    for (int h=0; h<nirrep; ++h){
       for (int cnt=0; cnt<moOei.rowspi(h); ++cnt){
          fprintf(outfile, "%1d ",h);
       }
    }
    fprintf(outfile, "\n");
    fprintf(outfile, "DOCC = ");
    for (int h=0; h<nirrep; h++){
       fprintf(outfile, "%2d ",docc[h]);
    }
    fprintf(outfile, "\n");
    fprintf(outfile, "SOCC = ");
    for (int h=0; h<nirrep; h++){
       fprintf(outfile, "%2d ",socc[h]);
    }
    fprintf(outfile, "\n");

    
    int nTriMo = nmo * (nmo + 1) / 2;
    double *temp = new double[nTriMo];
    Matrix moOei("MO OEI", nIrreps, orbspi, orbspi);
    IWL::read_one(psio.get(), PSIF_OEI, PSIF_MO_OEI, temp, nTriMo, 0, 0, outfile);
    moOei.set(temp);
    
    int * orbitalIrreps = new int[nmo];

    int SyGroup = 0;
    bool stopFindGN = false;
    std::string SymmLabel = Process::environment.molecule()->sym_label();
    do {
        if (SymmLabel.compare(CheMPS2::Irreps::getGroupName(SyGroup))==0) stopFindGN = true;
        else SyGroup += 1;
    } while ((!stopFindGN) && (SyGroup<42));
    fprintf(outfile, "If anything went wrong: Is ");
    fprintf(outfile, SymmLabel.c_str());
    fprintf(outfile, " equal to ");
    fprintf(outfile, (CheMPS2::Irreps::getGroupName(SyGroup)).c_str());
    fprintf(outfile, " ?\n");

    int counterFillOrbitalIrreps = 0;
    for (int h=0; h<nirrep; ++h){
       for (int cnt=0; cnt<moOei.rowspi(h); ++cnt){
          orbitalIrreps[counterFillOrbitalIrreps] = h;
          counterFillOrbitalIrreps++;
       }
    }


    CheMPS2::Hamiltonian * Ham = new CheMPS2::Hamiltonian(nmo,SyGroup,orbitalIrreps);
    delete [] orbitalIrreps;
    
    double NuclRepulsion =  Process::environment.molecule()->nuclear_repulsion_energy();
    Ham->setEconst(NuclRepulsion);
    double EnergyHF = NuclRepulsion;
    
    fprintf(outfile, "****  MO OEI \n");

    int nTot = 0;
    for(int h = 0; h < nirrep; ++h){
       for (int cnt = 0; cnt < moOei.rowspi(h); cnt++){
          for (int cnt2 = cnt; cnt2 < moOei.colspi(h); cnt2++){
             fprintf(outfile, "%1d %1d %16.48f \n",
                               nTot+cnt, nTot+cnt2, moOei[h][cnt][cnt2]);
             Ham->setTmat(nTot+cnt, nTot+cnt2, moOei[h][cnt][cnt2]);

          }
          if (cnt <docc[h])                            EnergyHF += 2*moOei[h][cnt][cnt];
          if ((cnt>=docc[h]) && (cnt<docc[h]+socc[h])) EnergyHF +=   moOei[h][cnt][cnt];
       }
       nTot += moOei.rowspi(h);
    }

    fprintf(outfile, "****  MO TEI \n");
 
    /*
     * Now, loop over the DPD buffer, printing the integrals
     */
    dpdbuf4 K;
    psio->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);
    // To only process the permutationally unique integrals, change the ID("[A,A]") to ID("[A>=A]+")
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[A,A]"), ID("[A,A]"),
                  ID("[A>=A]+"), ID("[A>=A]+"), 0, "MO Ints (AA|AA)");
    for(int h = 0; h < nirrep; ++h){
        global_dpd_->buf4_mat_irrep_init(&K, h);
        global_dpd_->buf4_mat_irrep_rd(&K, h);
        for(int pq = 0; pq < K.params->rowtot[h]; ++pq){
            int p = K.params->roworb[h][pq][0];
            int q = K.params->roworb[h][pq][1];
            int psym = K.params->psym[p];
            int qsym = K.params->qsym[q];
            int prel = p - K.params->poff[psym];
            int qrel = q - K.params->qoff[qsym];
            for(int rs = 0; rs < K.params->coltot[h]; ++rs){
                int r = K.params->colorb[h][rs][0];
                int s = K.params->colorb[h][rs][1];
                int rsym = K.params->rsym[r];
                int ssym = K.params->ssym[s];
                int rrel = r - K.params->roff[rsym];
                int srel = s - K.params->soff[ssym];
                // Print out the absolute orbital numbers, the relative (within irrep)
                // numbers, the symmetries, and the integral itself
                fprintf(outfile, "%1d %1d %1d %1d %16.48f \n",
                                 p, q, r, s, K.matrix[h][pq][rs]);
                Ham->setVmat(p,r,q,s,K.matrix[h][pq][rs]);
                if ((p==q) && (r==s)){
                   if ((prel <docc[psym]) && (rrel < docc[rsym])) EnergyHF += 2*K.matrix[h][pq][rs];
                   if ((prel>=docc[psym]) && (prel < socc[psym] + docc[psym]) && (rrel < docc[rsym])) EnergyHF += K.matrix[h][pq][rs];
                   if ((prel <docc[psym]) && (rrel>= docc[rsym]) && (rrel < socc[rsym] + docc[rsym])) EnergyHF += K.matrix[h][pq][rs];
                   if ((prel>=docc[psym]) && (prel < socc[psym] + docc[psym]) && (rrel>= docc[rsym]) && (rrel < socc[rsym] + docc[rsym])) EnergyHF += 0.5*K.matrix[h][pq][rs];
                }
                if ((p==s) && (r==q)){
                   if ((prel <docc[psym]) && (rrel < docc[rsym])) EnergyHF -= K.matrix[h][pq][rs];
                   if ((prel>=docc[psym]) && (prel < socc[psym] + docc[psym]) && (rrel < docc[rsym])) EnergyHF -= 0.5*K.matrix[h][pq][rs];
                   if ((prel <docc[psym]) && (rrel>= docc[rsym]) && (rrel < socc[rsym] + docc[rsym])) EnergyHF -= 0.5*K.matrix[h][pq][rs];
                   if ((prel>=docc[psym]) && (prel < socc[psym] + docc[psym]) && (rrel>= docc[rsym]) && (rrel < socc[rsym] + docc[rsym])) EnergyHF -= 0.5*K.matrix[h][pq][rs];
                }
            }
        }
        global_dpd_->buf4_mat_irrep_close(&K, h);
    }
    global_dpd_->buf4_close(&K);
    psio->close(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);
    fprintf(outfile, "****  HF Energy = %16.48f \n", EnergyHF);

    fprintf(outfile, "****  The debug check of Hamiltonian ****\n");
    fprintf(outfile, "Econst = %16.24f \n", Ham->getEconst());
   
    double test = 0.0;
    double test2 = 0.0;
    for (int i=0; i<Ham->getL(); i++){
       for (int j=0; j<Ham->getL(); j++){
          test += Ham->getTmat(i,j);
          if (i<=j) test2 += Ham->getTmat(i,j);
       }
    }
    fprintf(outfile, "1-electron integrals: Sum over all elements  : %16.24f \n", test);
    fprintf(outfile, "1-electron integrals: Sum over Tij with i<=j : %16.24f \n", test2);
      
    test = 0.0;
    test2 = 0.0;
    for (int i=0; i<Ham->getL(); i++){
       for (int j=0; j<Ham->getL(); j++){
          for (int k=0; k<Ham->getL(); k++){
             for (int l=0; l<Ham->getL(); l++){
                test += Ham->getVmat(i,j,k,l);
                if ((i<=j) && (j<=k) && (k<=l)) test2 += Ham->getVmat(i,j,k,l);
             }
          }
       }
    }
    fprintf(outfile, "2-electron integrals: Sum over all elements          : %16.24f \n", test);
    fprintf(outfile, "2-electron integrals: Sum over Vijkl with i<=j<=k<=l : %16.24f \n", test2);
 
    cout.precision(15);
    if(savehdf5)
        Ham->save2(filename);

    delete Ham;
    delete [] temp;

    if(savehdf5)
    {
        hid_t file_id = H5Fopen(CheMPS2::HAMILTONIAN_ParentStorageName.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
        HDF5_STATUS_CHECK(file_id);

        hid_t group_id = H5Gopen(file_id, ""/Data"", H5P_DEFAULT);
        HDF5_STATUS_CHECK(group_id);

        hid_t dataspace_id = H5Screate(H5S_SCALAR);
        hid_t dataset_id = H5Dcreate(group_id, "nelectrons", H5T_STD_I32LE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        hid_t status = H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &nelectrons);
        HDF5_STATUS_CHECK(status);

        status = H5Dclose(dataset_id);
        HDF5_STATUS_CHECK(status);

        status = H5Gclose(group_id);
        HDF5_STATUS_CHECK(status);

        status = H5Fclose(file_id);
        HDF5_STATUS_CHECK(status);
    }

    return Success;
}

}} // End Namespaces

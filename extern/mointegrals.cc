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
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>


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
    outfile->Printf("%s:%d: Error with HDF5. status code=%d\n", __FILE__, __LINE__, status); \
} 

INIT_PLUGIN

namespace psi{ namespace mointegrals{

extern "C" int
read_options(std::string name, Options &options)
{
    if (name == "MOINTEGRALS"|| options.read_globals()) {

        /*- The amount of information printed
            to the output file -*/
        options.add_bool("PRINT_INTEGRALS", false);
        /*- Whether to compute two-electron integrals -*/
        options.add_bool("DO_TEI", true);
        // save the unitary transformation from SO to MO integrals
        options.add_bool("SAVE_MO", false);
        options.add_str_i("U_MO_FILENAME", "unitary-mo.h5");
        // save to a HDF5 file
        options.add_bool("SAVEHDF5", false);
        options.add_str_i("HDF5_FILENAME", "integrals.h5");
    }

    return true;
}


extern "C" PsiReturnType
mointegrals(Options &options)
{
    const bool print = options.get_bool("PRINT_INTEGRALS");
    const bool doTei = options.get_bool("DO_TEI");
    const bool save_mo = options.get_bool("SAVE_MO");
    const bool savehdf5 = options.get_bool("SAVEHDF5");
    std::string filename = options.get_str("HDF5_FILENAME");
    std::string mofilename = options.get_str("U_MO_FILENAME");
    boost::algorithm::to_lower(filename);
    boost::algorithm::to_lower(mofilename);
   
    // Grab the global (default) PSIO object, for file I/O
    boost::shared_ptr<PSIO> psio(_default_psio_lib_);

    // Now we want the reference (SCF) wavefunction
    boost::shared_ptr<Wavefunction> wfn = Process::environment.wavefunction();
    if(!wfn) throw PSIEXCEPTION("SCF has not been run yet!");
    
    const int nirrep  = wfn->nirrep();
    
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

    boost::shared_ptr<Molecule> molecule = Process::environment.molecule();
    int nelectrons = 0;
    for(int i=0;i<molecule->natom();i++)
        nelectrons += molecule->true_atomic_number(i);

    nelectrons -= molecule->molecular_charge();


    int nTriMo = nmo * (nmo + 1) / 2;
    double *temp = new double[nTriMo];
    Matrix moOei("MO OEI", nIrreps, orbspi, orbspi);
    IWL::read_one(psio.get(), PSIF_OEI, PSIF_MO_OEI, temp, nTriMo, 0, 0, "outfile");
    moOei.set(temp);

    outfile->Printf("****  Molecular Integrals For CheMPS Start Here \n");
    
    std::string SymmLabel =  Process::environment.molecule()->sym_label();
    outfile->Printf("Symmetry Label = ");
    outfile->Printf(SymmLabel.c_str());
    outfile->Printf(" \n");
    outfile->Printf("Nirreps = %1d \n", nirrep);
    double NuclRepulsion =  Process::environment.molecule()->nuclear_repulsion_energy();
    outfile->Printf("Nuclear Repulsion Energy = %16.48f \n", NuclRepulsion);
    outfile->Printf("Number Of Molecular Orbitals = %2d \n", nmo);
    outfile->Printf("Irreps Of Molecular Orbitals = \n");
    for (int h=0; h<nirrep; ++h){
       for (int cnt=0; cnt<moOei.rowspi(h); ++cnt){
          outfile->Printf("%1d ",h);
       }
    }
    outfile->Printf("\n");
    outfile->Printf("DOCC = ");
    for (int h=0; h<nirrep; h++){
       outfile->Printf("%2d ",docc[h]);
    }
    outfile->Printf("\n");
    outfile->Printf("SOCC = ");
    for (int h=0; h<nirrep; h++){
       outfile->Printf("%2d ",socc[h]);
    }
    outfile->Printf("\n");


    if(save_mo)
    {
        // contains the transformation from SO -> MO
        SharedMatrix ca = Process::environment.wavefunction()->Ca();
        ca->print();

        SharedMatrix S = Process::environment.wavefunction()->S();
        S->print();

        boost::shared_ptr<MatrixFactory> factory = Process::environment.wavefunction()->matrix_factory();

        SharedMatrix aotoso = Process::environment.wavefunction()->sobasisset()->petite_list()->aotoso();
        SharedMatrix sotoao = Process::environment.wavefunction()->sobasisset()->petite_list()->sotoao();

        boost::shared_ptr<IntegralFactory> so_integral = Process::environment.wavefunction()->integral();

        shared_ptr<OneBodySOInt> sOBI(so_integral->so_overlap());
        shared_ptr<OneBodySOInt> tOBI(so_integral->so_kinetic());
        shared_ptr<OneBodySOInt> vOBI(so_integral->so_potential());
        // Form the one-electron integral matrices from the matrix factory
        SharedMatrix sMat(factory->create_matrix("Overlap"));
        SharedMatrix tMat(factory->create_matrix("Kinetic"));
        SharedMatrix vMat(factory->create_matrix("Potential"));
        SharedMatrix hMat(factory->create_matrix("One Electron Ints"));
        SharedMatrix hMat2(factory->create_matrix("One Electron Ints 2"));

        sOBI->compute(sMat);
        tOBI->compute(tMat);
        vOBI->compute(vMat);

        // Form h = T + V by first cloning T and then adding V
        hMat->copy(tMat);
        hMat->add(vMat);

        if(print)
        {
            sMat->print();
            hMat->print();


            outfile->Printf("so2ao:\n");
            sotoao->print();
        }


        // Construct Shalf
        SharedMatrix eigvec= factory->create_shared_matrix("L");
        SharedMatrix temp= factory->create_shared_matrix("Temp");
        SharedMatrix temp2= factory->create_shared_matrix("Temp2");
        SharedVector eigval(factory->create_vector());

        S->diagonalize(eigvec, eigval);

        // Convert the eigenvales to 1/sqrt(eigenvalues)
        int *dimpi = eigval->dimpi();
        double min_S = fabs(eigval->get(0,0));
        for (int h=0; h<nIrreps; ++h)
            for (int i=0; i<dimpi[h]; ++i)
            {
                if (min_S > eigval->get(h,i))
                    min_S = eigval->get(h,i);
                //double scale = 1.0 / sqrt(eigval->get(h, i));
                double scale = sqrt(eigval->get(h, i));
                eigval->set(h, i, scale);
            }

        outfile->Printf("Lowest eigenvalue of overlap S = %14.10E\n", min_S);

        if(min_S < options.get_double("S_TOLERANCE") )
            outfile->Printf("WARNING: Min value of overlap below treshold!!!!\n");


        // Create a vector matrix from the converted eigenvalues
        temp2->set_diagonal(eigval);

        temp->gemm(false, true, 1.0, temp2, eigvec, 0.0);
        S->gemm(false, false, 1.0, eigvec, temp, 0.0);
        if(print)
        {
            outfile->Printf("S^1/2:\n");
            S->print();
        }

        // order of multi. is changed to have the output compatible to UnitaryMatrix class.
        // (columns/rows interchanged)
        temp->gemm(false, false, 1.0, S, ca, 0.0);
//        temp->gemm(false, false, 1.0, ca, S, 0.0);
//        temp->transpose_this();

//        S->invert();
        hMat2->transform(hMat, ca);

        if(print)
        {
            outfile->Printf("Mo test:\n");
            hMat2->print();

            outfile->Printf("Matrix to store:\n");
            temp->print();
        }

        boost::shared_ptr<PetiteList> pl = Process::environment.wavefunction()->sobasisset()->petite_list();

        // need dimensions
        const Dimension aos = pl->AO_basisdim();
        const Dimension sos = pl->SO_basisdim();
        const Dimension nmo = ca->colspi();

        SharedMatrix Ca_ao_mo(new Matrix("temp_ca AO x MO", aos, nmo));

        // do the half transform
        Ca_ao_mo->gemm(false, false, 1.0, aotoso, ca, 0.0);

        Ca_ao_mo->print();

        hid_t file_id = H5Fcreate(mofilename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
        hid_t group_id = H5Gcreate(file_id, "/Data", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        const Dimension &ca_dims = ca->rowspi();

        for (int irrep=0; irrep< nirrep; irrep++)
        {
            int norb = ca_dims[irrep];

            if(norb > 0)
            {
                std::stringstream irrepname;
                irrepname << "irrep_" << irrep;

                double *p_data = temp->pointer(irrep)[0];

                hsize_t dimarray      = norb * norb;
                hid_t dataspace_id    = H5Screate_simple(1, &dimarray, NULL);
                hid_t dataset_id      = H5Dcreate(group_id, irrepname.str().c_str(), H5T_IEEE_F64LE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
                H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, p_data);

                H5Dclose(dataset_id);
                H5Sclose(dataspace_id);
            }
        }

        H5Gclose(group_id);
        H5Fclose(file_id);

        SharedMatrix Ca_ao_ao(new Matrix("ca AO x AO", aos, aos));

        Ca_ao_ao->remove_symmetry(ca, sotoao);
        Ca_ao_ao->print();

        /*
        std::string mofilename2 = "S-";
        mofilename2 += mofilename;

        file_id = H5Fcreate(mofilename2.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
        group_id = H5Gcreate(file_id, "/Data", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        for (int irrep=0; irrep< nirrep; irrep++)
        {
            int norb = ca_dims[irrep];

            if(norb > 0)
            {
                std::stringstream irrepname;
                irrepname << "irrep_" << irrep;

                double *p_data = Ssqrt->pointer(irrep)[0];

                hsize_t dimarray      = norb * norb;
                hid_t dataspace_id    = H5Screate_simple(1, &dimarray, NULL);
                hid_t dataset_id      = H5Dcreate(group_id, irrepname.str().c_str(), H5T_IEEE_F64LE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
                H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, p_data);

                H5Dclose(dataset_id);
                H5Sclose(dataspace_id);
            }
        }

        H5Gclose(group_id);
        H5Fclose(file_id);
        */
    }

    
    
    int * orbitalIrreps = new int[nmo];

    int SyGroup = 0;
    bool stopFindGN = false;
    do {
        if (SymmLabel.compare(CheMPS2::Irreps::getGroupName(SyGroup))==0) stopFindGN = true;
        else SyGroup += 1;
    } while ((!stopFindGN) && (SyGroup<42));
    outfile->Printf("If anything went wrong: Is ");
    outfile->Printf(SymmLabel.c_str());
    outfile->Printf(" equal to ");
    outfile->Printf((CheMPS2::Irreps::getGroupName(SyGroup)).c_str());
    outfile->Printf(" ?\n");

    int counterFillOrbitalIrreps = 0;
    for (int h=0; h<nirrep; ++h){
       for (int cnt=0; cnt<moOei.rowspi(h); ++cnt){
          orbitalIrreps[counterFillOrbitalIrreps] = h;
          counterFillOrbitalIrreps++;
       }
    }


    CheMPS2::Hamiltonian * Ham = new CheMPS2::Hamiltonian(nmo,SyGroup,orbitalIrreps);
    delete [] orbitalIrreps;

    Ham->setNe(nelectrons);
    
    Ham->setEconst(NuclRepulsion);
    Ham->reset();

    double EnergyHF = NuclRepulsion;
    
    outfile->Printf("****  MO OEI \n");

    int nTot = 0;
    for(int h = 0; h < nirrep; ++h){
       for (int cnt = 0; cnt < moOei.rowspi(h); cnt++){
          for (int cnt2 = cnt; cnt2 < moOei.colspi(h); cnt2++){
              if(print)
                  outfile->Printf("%1d %1d %16.48f \n",
                          nTot+cnt, nTot+cnt2, moOei[h][cnt][cnt2]);

              Ham->setTmat(nTot+cnt, nTot+cnt2, moOei[h][cnt][cnt2]);

          }
          if (cnt <docc[h])                            EnergyHF += 2*moOei[h][cnt][cnt];
          if ((cnt>=docc[h]) && (cnt<docc[h]+socc[h])) EnergyHF +=   moOei[h][cnt][cnt];
       }
       nTot += moOei.rowspi(h);
    }

    outfile->Printf("****  MO TEI \n");
 
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
                if(print)
                    outfile->Printf("%1d %1d %1d %1d %16.48f \n",
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
    outfile->Printf("****  HF Energy = %16.48f \n", EnergyHF);

    outfile->Printf("****  The debug check of Hamiltonian ****\n");
    outfile->Printf("Econst = %16.24f \n", Ham->getEconst());
   
    double test = 0.0;
    double test2 = 0.0;
    for (int i=0; i<Ham->getL(); i++){
       for (int j=0; j<Ham->getL(); j++){
          test += Ham->getTmat(i,j);
          if (i<=j) test2 += Ham->getTmat(i,j);
       }
    }
    outfile->Printf("1-electron integrals: Sum over all elements  : %16.24f \n", test);
    outfile->Printf("1-electron integrals: Sum over Tij with i<=j : %16.24f \n", test2);
      
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
    outfile->Printf("2-electron integrals: Sum over all elements          : %16.24f \n", test);
    outfile->Printf("2-electron integrals: Sum over Vijkl with i<=j<=k<=l : %16.24f \n", test2);
 
    cout.precision(15);
    if(savehdf5)
        Ham->save2(filename);

    delete Ham;
    delete [] temp;

    return Success;
}

}} // End Namespaces

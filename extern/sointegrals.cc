#include <libplugin/plugin.h>
#include "psi4-dec.h"
#include <libparallel/parallel.h>
#include <liboptions/liboptions.h>
#include <libmints/mints.h>
#include <libmints/sointegral_twobody.h>
#include <libpsio/psio.h>

#include <hdf5.h>

#include "Irreps.h"
#include "Hamiltonian.h"

INIT_PLUGIN

// macro to help check return status of HDF5 functions
#define HDF5_STATUS_CHECK(status) { \
    if(status < 0) \
    fprintf(outfile, "%s:%d: Error with HDF5. status code=%d\n", __FILE__, __LINE__, status); \
} 

using namespace boost;

namespace psi{ namespace sointegrals{

extern "C" int
read_options(std::string name, Options &options)
{
    if (name == "SOINTEGRALS"|| options.read_globals()) {
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

class ERIPrinter
{
public:

    ERIPrinter() { count = 0; }

    shared_ptr<CheMPS2::Hamiltonian> Ham;

    unsigned int count;

    // Our functor...the nice thing about using a C++ functor is that the
    // code here is inlined by the compiler.
    void operator() (int pabs, int qabs, int rabs, int sabs,
                     int pirrep, int pso,
                     int qirrep, int qso,
                     int rirrep, int rso,
                     int sirrep, int sso,
                     double value)
    {
        fprintf(outfile, "%1d %1d %1d %1d %16.48f \n", 
			pabs, qabs, rabs, sabs, value);
        Ham->setVmat(pabs, rabs, qabs, sabs, value);

        count++;
    }
};

extern "C" PsiReturnType
sointegrals(Options &options)
{
    bool print = options.get_bool("PRINT_INTEGRALS");
    bool doTei = options.get_bool("DO_TEI");
    bool savehdf5 = options.get_bool("SAVEHDF5");
    std::string filename = options.get_str("HDF5_FILENAME");

    if(options.get_str("S_ORTHOGONALIZATION") != "SYMMETRIC")
    {
        fprintf(outfile, "This will only work with symmetric orthogonalisation: AO == OSO\n");
        return Failure;
    }

    shared_ptr<Molecule> molecule = Process::environment.molecule();

    // Form basis object:
    // Create a basis set parser object.
    shared_ptr<BasisSetParser> parser(new Gaussian94BasisSetParser());
    // Construct a new basis set.
    shared_ptr<BasisSet> aoBasis = BasisSet::construct(parser, molecule, "BASIS");

    // The integral factory oversees the creation of integral objects
    shared_ptr<IntegralFactory> integral(new IntegralFactory
            (aoBasis, aoBasis, aoBasis, aoBasis));
    
    // N.B. This should be called after the basis has been built, because
    // the geometry has not been
    // fully initialized until this time.
    molecule->print();
    // The basis set is also created from the information stored in the
    // checkpoint file
    // ...it needs to raid the checkpoint file to find those dimensions.
    // Create an SOBasis object using the basis set and integral factory.
    shared_ptr<SOBasisSet> soBasis(new SOBasisSet(aoBasis, integral));

    // Obtain block dimensions from the SO basis
    const Dimension dimension = soBasis->dimension();

    // The matrix factory can create matrices of the correct dimensions...
    shared_ptr<MatrixFactory> factory(new MatrixFactory);
    factory->init_with(dimension, dimension);

    // Form the one-electron integral objects from the integral factory
    shared_ptr<OneBodySOInt> sOBI(integral->so_overlap());
    shared_ptr<OneBodySOInt> tOBI(integral->so_kinetic());
    shared_ptr<OneBodySOInt> vOBI(integral->so_potential());
    // Form the one-electron integral matrices from the matrix factory
    SharedMatrix sMat(factory->create_matrix("Overlap"));
    SharedMatrix tMat(factory->create_matrix("Kinetic"));
    SharedMatrix vMat(factory->create_matrix("Potential"));
    SharedMatrix hMat(factory->create_matrix("One Electron Ints"));

    const int nmo = dimension.sum();
    const int nirreps   = dimension.n();
    const double NuclRepulsion =  molecule->nuclear_repulsion_energy();

    int nelectrons = 0;
    for(int i=0;i<molecule->natom();i++)
        nelectrons += molecule->true_atomic_number(i);

    nelectrons -= molecule->molecular_charge();

    fprintf(outfile, "****  Molecular Integrals Start Here \n");

    std::string SymmLabel =  molecule->sym_label();
    fprintf(outfile, "Symmetry Label = %s\n", SymmLabel.c_str());
    fprintf(outfile, "Nuclear Repulsion Energy = %16.48f \n", NuclRepulsion);
    fprintf(outfile, "Nelectrons = %1d \n", nelectrons);
    fprintf(outfile, "Nirreps = %1d \n", nirreps);
    fprintf(outfile, "Dimension of Irreps = ");
    for (int h=0; h<nirreps; ++h)
        fprintf(outfile, "%2d ", dimension[h]);
    fprintf(outfile, "\n");

    fprintf(outfile, "Number Of SO Orbitals = %2d \n", nmo);
    fprintf(outfile, "Irreps Of SO Orbitals = ");
    for (int h=0; h<nirreps; ++h)
        for (int i=0; i<dimension[h]; ++i)
            fprintf(outfile, "%2d ", h);
    fprintf(outfile, "\n");

    std::vector<int> OrbIrreps;
    OrbIrreps.reserve(nmo);
    for (int h=0; h<nirreps; ++h)
        for (int i=0; i<dimension[h]; ++i)
            OrbIrreps.push_back(h);

    // Compute the one electron integrals, telling each object where to
    // store the result
    sOBI->compute(sMat);
    tOBI->compute(tMat);
    vOBI->compute(vMat);

    if(print)
    {
        sMat->print();
        tMat->print();
        vMat->print();
    }
    // Form h = T + V by first cloning T and then adding V
    hMat->copy(tMat);
    hMat->add(vMat);

    // Construct Shalf
    /*
    SharedMatrix eigvec= factory->create_shared_matrix("L");
    SharedMatrix temp= factory->create_shared_matrix("Temp");
    SharedMatrix temp2= factory->create_shared_matrix("Temp2");
    SharedVector eigval(factory->create_vector());

    //Used to do this 3 times, now only once
    sMat->diagonalize(eigvec, eigval);

    // Convert the eigenvales to 1/sqrt(eigenvalues)
    int *dimpi = eigval->dimpi();
    double min_S = fabs(eigval->get(0,0));
    for (int h=0; h<nirreps; ++h)
        for (int i=0; i<dimpi[h]; ++i)
        {
            if (min_S > eigval->get(h,i))
                min_S = eigval->get(h,i);
            double scale = 1.0 / sqrt(eigval->get(h, i));
            eigval->set(h, i, scale);
        }

    fprintf(outfile, "Lowest eigenvalue of overlap S = %14.10E\n", min_S);

    if(min_S < options.get_double("S_TOLERANCE") )
        fprintf(outfile, "WARNING: Min value of overlap below treshold!!!!\n");


    // Create a vector matrix from the converted eigenvalues
    temp2->set_diagonal(eigval);

    SharedMatrix ca = Process::environment.wavefunction()->Ca();

//    temp->gemm(true, false, 1.0, ca, hMat, 0.0);
//    hMat->gemm(false, false, 1.0, temp, ca, 0.0);

    temp->gemm(false, true, 1.0, temp2, eigvec, 0.0);
    sMat->gemm(false, false, 1.0, eigvec, temp, 0.0);

    temp->gemm(false, true, 1.0, hMat, sMat, 0.0);
    hMat->gemm(false, false, 1.0, sMat, temp, 0.0);
    */

    hMat->print();

    int SyGroup = 0;
    bool stopFindGN = false;
    do {
        if (SymmLabel.compare(CheMPS2::Irreps::getGroupName(SyGroup))==0)
            stopFindGN = true;
        else
            SyGroup += 1;
    } 
    while ((!stopFindGN) && (SyGroup<42));

    fprintf(outfile, "If anything went wrong: Is %s equal to %s?\n", SymmLabel.c_str(), (CheMPS2::Irreps::getGroupName(SyGroup)).c_str());

    shared_ptr<CheMPS2::Hamiltonian> Ham(new CheMPS2::Hamiltonian(nmo, SyGroup, OrbIrreps.data()));
    Ham->setEconst(NuclRepulsion);
    Ham->setNe(nelectrons);
    Ham->reset();

    fprintf(outfile, "*** OEI\n");

    int count = 0;
    for (int h=0; h<nirreps; ++h)
    {
        for (int i=0; i<dimension[h]; ++i)
            for (int j=i; j<dimension[h]; ++j)
            {
                fprintf(outfile, "%1d %1d %16.48f \n", count+i, count+j, hMat->get(h,i,j));
                Ham->setTmat(count+i, count+j, hMat->get(h,i,j));
            }

        count += dimension[h];
    }


    if(doTei){
        // 1. Obtain an object that knows how to compute two-electron AO
        // integrals.
        shared_ptr<TwoBodyAOInt> tb(integral->eri());
        
        // 2. Create an object that knows how to convert any two-body AO
        // integral to SO.
        shared_ptr<TwoBodySOInt> eri(new TwoBodySOInt(tb, integral));
        
        // 3. Find out how many SO shells we have.
        int nsoshell = soBasis->nshell();

        // 4. We to create an instance of our ERIPrinter
        ERIPrinter printer;
        printer.Ham = Ham;

        fprintf(outfile, "*** TEI\n");
        
        // 5. Create an SOShellCombintationsIterator to step through the
        // necessary combinations
        SOShellCombinationsIterator shellIter(soBasis, soBasis, soBasis, soBasis);
        for (shellIter.first(); shellIter.is_done() == false; shellIter.next()) {
            // 6. Call the TwoBodySOInt object to compute integrals giving
            // it the
            // instance to our functor.
            eri->compute_shell(shellIter, printer);
        }

        fprintf(outfile, "number of TEI: %d\n", printer.count);
    }

    if(savehdf5)
        Ham->save2(filename);

    return Success;
}

}} // End namespaces

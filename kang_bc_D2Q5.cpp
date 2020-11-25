#include "palabos2D.h"
#include "palabos2D.hh"

using namespace std;
using namespace plb;

typedef double T;

#define ADESCRIPTOR descriptors::AdvectionDiffusionD2Q5Descriptor
#define ADYNAMICS AdvectionDiffusionBGKdynamics

template<typename T, template<typename U> class Descriptor>
class KangBCTop_BDY : public BoxProcessingFunctional2D_L<T,Descriptor> 
// {{{ 
{
public:
    KangBCTop_BDY(T kr_,T C_eq_)
        : kr(kr_),C_eq(C_eq_)
    { }
    virtual void process(   Box2D domain,
                            BlockLattice2D<T,Descriptor>& lattice)
    {
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            Cell <T,Descriptor>& cell = lattice.get(iX,iY);
            T temp_f4 = cell[4] + Descriptor<T>::t[4];
            cell[2] = (2.0*temp_f4+kr*C_eq)/(3.0*kr+1.0) - temp_f4-Descriptor<T>::t[2];
            }
        }
    }

    virtual KangBCTop_BDY<T,Descriptor>* clone() const {
        return new KangBCTop_BDY(*this);
    }

    void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        modified[0] = modif::staticVariables;
    }
    
    virtual BlockDomain::DomainT appliesTo() const {
        return BlockDomain::bulk;
    }
private:
    T kr;
    T C_eq;
};
//}}}

void assign_params( string f_name,
                    plint &maxIter,
                    T &epsilon,
                    plint &nx,
                    plint &ny,
                    T &Da,
                    T &D,
                    T &C_0,
                    T &C_eq)
// {{{
{
    try{
    XMLreader xmlFile(f_name);
    pcout << "CONFIGURATION" << endl;
    pcout << "=============" << endl; 
    xmlFile.print(0);
    xmlFile["inputs"]["maxIter"].read(maxIter);
    xmlFile["inputs"]["epsilon"].read(epsilon);
    xmlFile["inputs"]["nx"].read(nx);
    xmlFile["inputs"]["ny"].read(ny);
    xmlFile["inputs"]["Da"].read(Da);
    xmlFile["inputs"]["D"].read(D);
    xmlFile["inputs"]["C_0"].read(C_0);
    xmlFile["inputs"]["C_eq"].read(C_eq);
    pcout << "=============" << endl << endl; 
    } catch (PlbIOException& exception) { 
        pcout << exception.what() << endl;
    }
}
//}}}

void computeNormalizedConcentration(MultiBlockLattice2D<T,ADESCRIPTOR>& lattice,
                                    MultiScalarField2D<T>& normalizedConcentration,
                                    T C_0,
                                    T C_eq)
//{{{
{
    Box2D domain = lattice.getBoundingBox();
    computeDensity(lattice,normalizedConcentration,domain);
    // Page W12S14 from Kang
    subtractInPlace(normalizedConcentration,C_eq,domain);
    divideInPlace(normalizedConcentration,C_0-C_eq,domain);
}
//}}}

void writeVTK(  MultiBlockLattice2D<T,ADESCRIPTOR>& lattice,
                MultiScalarField2D<T> normalizedConcentration,
                plint iter,
                T C_0,
                T C_eq)
// {{{
    {
    VtkImageOutput2D<T> vtkOut(createFileName("vtk", iter, 6), 1);
    vtkOut.writeData<float>(*computeDensity(lattice), "concentration", 1);
    vtkOut.writeData<float>(normalizedConcentration, "normalizedSoluteConcentration", 1);
    }
// }}}

void rectangularDomainSetup(MultiBlockLattice2D<T,ADESCRIPTOR> &adLattice,
OnLatticeAdvectionDiffusionBoundaryCondition2D<T,ADESCRIPTOR>& boundaryCondition,
                            plint nx,
                            plint ny,
                            T C_0,
                            T C_eq,
                            T kr)
//{{{
    {
    // Defininig each boundary
    //                          The Corner will be reaction boundary?
    //                              |
    //                              v
    Box2D left_boundary(0, 0, 0, ny-1);
    Box2D top_boundary(0,nx-1,ny-1,ny-1);

    Array<T,2> u0((T)0,(T)0);                       // 2D velocity = 0 = (0,0)


    // Temperature == Concentration in this case -> is 0N correct?
    boundaryCondition.addTemperatureBoundary0N( left_boundary,
                                                adLattice);
                                                
    // Setting the density (concentration) to be a constant value on the left
    setBoundaryDensity( adLattice,
                        left_boundary,
                        C_0);


    plint processorLevel = 0;
    integrateProcessingFunctional(  new KangBCTop_BDY<T,ADESCRIPTOR>(kr,C_eq), 
                                    top_boundary,
                                    adLattice,
                                    processorLevel);

    initializeAtEquilibrium(adLattice,
                            adLattice.getBoundingBox(),
                            C_eq,
                            u0);

    adLattice.initialize();
    }
//}}}

int main(int argc, char* argv[])
//{{{
    {
    plbInit(&argc, &argv);
    global::directories().setOutputDir("./tmp/");

    // Instantiate variables
    plint maxIter;  // Maximum number of iterations
    T epsilon;      // Convergence Precision
    plint nx;       // x domain
    plint ny;       // y domain
    T D;            // Diffusivity
    T Da;           // Damkohler
    T C_0;          // Constant Concentration at left boundary
    T C_eq;         // Bulk Concentration

    // Read the parameters from params.xml file
    assign_params(  "params.xml",
                    maxIter,
                    epsilon,
                    nx,
                    ny,
                    Da,
                    D,
                    C_0,
                    C_eq);
     
    // Calculated values
    // Conversion to Tau from Diffusion Coefficient for D2Q5 (with rest
    // fraction == 0) -> C_q = 1/2
    T tau = 3.0*D;                  // Page 4 of Kang 2007: equation 18 for D2Q5 (not tau-1/2 !!!!!!!!!!)
    T adOmega = 1.0/tau;            // omega = 1/tau 
    T b = static_cast<T>(ny);
    T kr = (Da*D)/b;                // Rate constant (equation 50 from Kang)
    
    pcout << "b = " << b << endl;
    pcout << "tau = " << tau << endl;
    pcout << "omega = " << adOmega << endl;
    pcout << "rate constant = " << kr << endl;
    pcout << "============================================" << endl;

    // Instantiate the Lattice
    MultiBlockLattice2D<T,ADESCRIPTOR> adLattice(
                                    nx,
                                    ny,
                                    new ADYNAMICS<T, ADESCRIPTOR>(adOmega));

    // This scalar field will hold the values that are plotted in Kang's paper
    // on the 2D domain
    MultiScalarField2D<T> normalized_concentration(nx,ny);

    OnLatticeAdvectionDiffusionBoundaryCondition2D<T,ADESCRIPTOR>*
        boundaryCondition = createLocalAdvectionDiffusionBoundaryCondition2D<T,ADESCRIPTOR>();

    rectangularDomainSetup( adLattice,
                            *boundaryCondition,
                            nx,
                            ny,
                            C_0,
                            C_eq,
                            kr);

     
    
    plint convergenceIter = 250;
    plint iT = 0;
    //                           characteristic velocity of the system
    //                                |                  characteristic length of the system
    //                                |                   |  Precision of convergence
    //                                |                   |    |
    //                                v                   v    v
    util::ValueTracer<T> convergence(0.01*convergenceIter,b,epsilon);
    // ========================================================================
    // Time Stepping
    // ========================================================================
    while(!convergence.hasConverged() && iT < maxIter)
    //for (plint iT = 0; iT<=maxIter; iT++)
        {
        if (iT%convergenceIter==0) {
            computeNormalizedConcentration( adLattice,
                                            normalized_concentration,
                                            C_0,
                                            C_eq);
            convergence.takeValue(computeAverage(normalized_concentration),true);
            pcout << "Writing vtk_instance at iT = " << iT << std::endl;
            writeVTK(   adLattice,
                        normalized_concentration,
                        iT,
                        C_0,
                        C_eq);
        }
        adLattice.collideAndStream();
        ++iT;
    }

    computeNormalizedConcentration( adLattice,
                                    normalized_concentration,
                                    C_0,
                                    C_eq);
    plb_ofstream ofile("lbm_contour.dat");
    ofile << setprecision(10) << normalized_concentration << endl;
    }
// }}}

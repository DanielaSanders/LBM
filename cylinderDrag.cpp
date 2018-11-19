#include "palabos2D.h"
#include "palabos2D.hh"
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>

using namespace plb;
using namespace plb::descriptors;
using namespace std;

typedef double T;
#define DESCRIPTOR D2Q9Descriptor

//Velocity on parabolic Poiseuille profile
T poiseuilleVelocity(plint iY, IncomprFlowParam<T> const& parameters) {
    T y = (T)iY / parameters.getResolution();
    return 4.*parameters.getLatticeU() * (y-y*y); //velocity in lattice units
}

//Linearly decreasing pressure profile
T poiseuillePressure(plint iX, IncomprFlowParam<T> const& parameters){
    T Lx = parameters.getNx()-1; //number of Lattice cells in x-direction
    T Ly = parameters.getNy()-1; //number of Lattice cells in y-direction
    return 8.*parameters.getLatticeNu()*parameters.getLatticeU() 
            / (Ly*Ly) * (Lx/(T)2-(T)iX); //Nu: viscosity in lattice units, U: velocity in lattice units
}

//convert pressure to density
T poiseuilleDensity(plint iX, IncomprFlowParam<T> const& parameters){
    return poiseuillePressure(iX, parameters)*DESCRIPTOR<T>::invCs2 + (T)1;
}

//Initialize velocity for boundary conditions
template<typename T>
class PoiseuilleVelocity{
public:
    PoiseuilleVelocity(IncomprFlowParam<T> parameters_)
        : parameters(parameters_)
    { }
    void operator()(plint iX, plint iY, Array<T, 2>& u) const {
        u[0] = poiseuilleVelocity(iY, parameters);
        u[1] = T();
    }
private:
    IncomprFlowParam<T> parameters;
};

//Initialize pressure boundary to constant velocity
template<typename T>
class ConstantDensity {
public:
    ConstantDensity(T density_)
        : density(density_)
    { }
    T operator()(plint iX, plint iY) const{
        return density;
    }
private:
    T density;
};

// Create an initial condition for density and velocity
template<typename T>
class PoiseuilleVelocityAndDensity {
public:
    PoiseuilleVelocityAndDensity(IncomprFlowParam<T> parameters_)
        : parameters(parameters_)
    { }
    void operator()(plint iX, plint iY, T& rho, Array<T,2>& u) const {
        rho = poiseuilleDensity(iX,parameters);
        u[0] = poiseuilleVelocity(iY, parameters);
        u[1] = T();
    }
private:
    IncomprFlowParam<T> parameters;
};

// instantiate bounce-back nodes at the locations of the cylinder
template<typename T>
class CylinderShapeDomain2D : public plb::DomainFunctional2D {
public:
    CylinderShapeDomain2D(plb::plint cx_, plb::plint cy_, plb::plint radius)
        : cx(cx_),
          cy(cy_),
          radiusSqr(plb::util::sqr(radius))
    { }
    virtual bool operator() (plb::plint iX, plb::plint iY) const {
        return plb::util::sqr(iX-cx) + plb::util::sqr(iY-cy) <= radiusSqr;
    }
    virtual CylinderShapeDomain2D<T>* clone() const {
        return new CylinderShapeDomain2D<T>(*this);
    }
private: 
    plb::plint cx;
    plb::plint cy;
    plb::plint radiusSqr;
};


void cylinderSetup (MultiBlockLattice2D<T, DESCRIPTOR>& lattice,
                    IncomprFlowParam<T> const& parameters,
                    OnLatticeBoundaryCondition2D<T, DESCRIPTOR>& boundaryCondition,
                    Array<plint,2> forceIds )
{
    //Lattice dimensions
    const plint nx = parameters.getNx();
    const plint ny = parameters.getNy();

    //Boundary Conditions
    Box2D outlet(nx-1, nx-1, 1, ny-2);
    Box2D everythingButOutlet( 0, nx-2, 0, ny-1); 

    //Velocity Boundary conditions everywhere except right boundary
    boundaryCondition.setVelocityConditionOnBlockBoundaries(
            lattice, everythingButOutlet );
    boundaryCondition.setPressureConditionOnBlockBoundaries(
            lattice, outlet );
    
    setBoundaryVelocity (
            lattice, lattice.getBoundingBox(),
            PoiseuilleVelocity<T>(parameters) );
    setBoundaryDensity (
            lattice, outlet,
            ConstantDensity<T>(1.) );
    initializeAtEquilibrium(
        lattice, lattice.getBoundingBox(),
        PoiseuilleVelocityAndDensity<T>(parameters) );

//coordinates of cylinder
plint cx     = nx/4;   //position in x-direction
plint cy     = ny/2;   //position in y- direction
plint radius = cy/4;   //radius of cylinder

//Use of dynamic MomentumExchangeBounceBack for drag
defineDynamics(lattice, lattice.getBoundingBox(),
                new CylinderShapeDomain2D<T>(cx,cy,radius),
                new MomentumExchangeBounceBack<T,DESCRIPTOR>(forceIds));
initializeMomentumExchange(lattice, lattice.getBoundingBox(), 
                new CylinderShapeDomain2D<T>(cx, cy, radius) );

lattice.initialize();
}

void writeGif(MultiBlockLattice2D<T, DESCRIPTOR>& lattice, plint iter)
{
    const plint imSize = 600;

    ImageWriter<T> ImageWriter("leeloo");
    ImageWriter.writeScaledGif(createFileName("u", iter, 6), 
                            *computeVelocityNorm(lattice),
                            imSize, imSize );
}

void writeVTK(MultiBlockLattice2D<T,DESCRIPTOR>& lattice,
              IncomprFlowParam<T> const& parameters, plint iter)
{
    T dx = parameters.getDeltaX();
    T dt = parameters.getDeltaT();
    VtkImageOutput2D<T> vtkOut(createFileName("vtk", iter, 6), dx);
    vtkOut.writeData<float>(*computeVelocityNorm(lattice), "velocityNorm", dx/dt);
    vtkOut.writeData<2,float>(*computeVelocity(lattice), "velocity", dx/dt);
}

int main(int argc, char* argv[]) {
    plbInit(&argc, &argv);

    global::directories().setOutputDir("./tmp/");

    IncomprFlowParam<T> parameters(
        (T) 1e-2,   //uMax  maximum velocity
        (T) 1.,     //Re    Reynolds number
        200,        //N     number of discrete particle
        6.,         //lx    lattice size x
        1.          //ly    lattice size y
    );
    
    const T logT    = (T)0.05 ;
    const T imSave  = (T)0.1 ;
    const T vtkSave = (T)0.5 ;
    const T maxT    = (T)20.1 ;

    writeLogFile(parameters, "Poiseuille flow");

    MultiBlockLattice2D<T, DESCRIPTOR> lattice(
            parameters.getNx(), parameters.getNy(),
            new BGKdynamics<T, DESCRIPTOR>(parameters.getOmega()) );
    lattice.initialize();

    Array<plint,2> forceIds;
    forceIds[0] = lattice.internalStatSubscription().subscribeSum();
    forceIds[1] = lattice.internalStatSubscription().subscribeSum();

    
    OnLatticeBoundaryCondition2D<T,DESCRIPTOR>*
        boundaryCondition = createInterpBoundaryCondition2D<T, DESCRIPTOR>();

    cylinderSetup(lattice, parameters, *boundaryCondition, forceIds);

    // Main loop over time iterations
    for (plint iT=0; iT*parameters.getDeltaT()<maxT; ++iT) {
        
        if (iT%parameters.nStep(imSave)==0){
            pcout << "Saving Gif ..." << endl;
            writeGif(lattice, iT);
        }

        if (iT%parameters.nStep(vtkSave)==0 && iT>0) {
            pcout << "Saving VTK file ..." << endl;
            writeVTK(lattice, parameters, iT);
        }

        if (iT%parameters.nStep(logT)==0) {
            pcout << "step " << iT 
                  //<< "; lattice time =" << lattice.getTimeCounter().getTime()
                  << "; t =" << iT*parameters.getDeltaT();
        }

        // Lattice Boltzmann iteration step
        lattice.collideAndStream();

        if (iT%parameters.nStep(logT)==0) {
            pcout << "; av energy ="
                  << setprecision(10) << getStoredAverageEnergy<T>(lattice)
                  << "; av rho ="
                  << getStoredAverageDensity<T>(lattice) 
                  << "; drag force =" << lattice.getInternalStatistics().getSum(
                      forceIds[0])
                  << "; drag coefficient =" << (2* lattice.getInternalStatistics().getSum(forceIds[0])) 
                     / (getStoredAverageDensity<T>(lattice)*((1e-2)*(1e-2))*((parameters.getNx()+4)/8)*((parameters.getNx()+4)/8)*3.14159)
                  << endl;

        }
    }

    delete boundaryCondition;
}
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
T poiseuilleVelocity(plint iY, ImcomprFlowParam<T> const& parameters) {
    T y = (T)iY / parameters.getResolution();
    return 4.*parameters.getLatticeU() * (y-y*y); //velocity in lattice units
}

//Linearly decreasing pressure profile
T poiseuillePressure(plint iX, IncomprFlowParam<T> const& parameters){
    T Lx = parameters.getNx()-1; //number of Lattice cells in x-direction
    T Ly = parameters.getNy()-1; //number of Lattice cells in y-direction
    return 8.*parameters.getLatticeNu()*parameters.getLatticeU() 
            / (Ly*Ly) * (Lx/(T)2-(T)iX); //viscosity in lattice units
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
}

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
                    OnLatticeBoundaryCondition2D<T, DESCRIPTOR>& boundaryCondition)
{
    const plint nx = parameters.getNx();
    const plint ny = parameters.getNy();

    //Boundary Conditions
    Box2D inlet(0, 0, 1, ny-2 );
    Box2D outlet(nx-1, nx-1, 1, ny-2);
    Box2D bottomWall( 0, nx-1, 0, 0);
    Box2D topWall( 0, nx-1, ny-1, ny-1); 

    //Dirichlet condition
    boundaryCondition.setVelocityConditionOnBlockBoundaries(
            lattice, inlet );
    //outflow condition
    boundaryCondition.setVelocityConditionOnBlockBoundaries(
            lattice, outlet, boundary::outflow );
    //free-slip boundary conditions
    boundaryCondition.setVelocityConditionOnBlockBoundaries(
            lattice, bottomWall, boundary::freeslip );
    boundaryCondition.setVelocityConditionOnBlockBoundaries(
            lattice, topWall, boundary::freeslip );

    setBoundaryVelocity (
            lattice, lattice.getBoundingBox(),
            PoiseuilleVelocity<T>(parameters) );
    setBoundaryDensity (
            lattice, outlet,
            ConstantDensity<T>(1.) );
    initializeAtEquilibrium(
        lattice, lattice.getBoundingBox(),
        PoiseuilleVelocity<T>(parameters) );

plint cx     = nx/4;
plint cy     = ny/22; 
plint radius = cy/4;
defineDynamics(lattice, lattice.getBoundingBox(),
                new CylinderShapeDomain2D<T>(cx,cy,radius),
                new plb::BounceBack<T,DESCRIPTOR>);

lattice.initialize();
}

void writeGif(MultiBlockLattice2D<T, DESCRIPTOR>& lattice, plint iter)
{
    ImageWriter<T> ImageWriter("leeloo");
    imageWriter.writeScaleGif(createFileName("u", iter, 6), 
                            *computeVelocityNorm(lattice) );
}

void writeVTK(MultiBlockLattice2D<T, DESCRIPTOR>& lattice,
                IncomprFlowParam<T> const& parameters, plint iter)
{
    T dx = parameters.getDeltaT();
    T dy = parameters.getDeltaT();
    VTKImageOutput2D<T> vtkOut(createFileName("vtk", iter, 6), dx);
    vtkOut.writeData<float>(*computeVelocityNorm(lattice), "velocityNorm", dx/dt);
    vtkOut.writeData<2,float>(*computeVelocity(lattice), "velocity", dx/dt);
}

int main(int argc, char* argv[]) {
    plbInit(&argc, &argv);

    global::directories().setOutputDir("./tmp/");

    IncomprFlowParam<T> parameters(
        (T) ,    //uMax
        (T) 0.,     //Re
        100,        //N
        60.,        //lx
        10.,        //ly
    );
    const T logT    = (T) 0.02;
    const T imSave  = (T)0.2 ;
    const T vtkSav  = (T)1. ;
    const T maxT    = (T)10.1 ;

    wirteLogFile(paramters, "");

    MultiBlockLattice2D<T, DESCRIPTOR> lattice(
        parameters.getNx(), parameters.getNy(),
        new BGKdynamics<T, DESCRIPTOR>(parameters.getOmega()) );
    
    OnLatticeBoundaryCondition2D<T,DESCRIPTOR>*
        boundaryCondition = createLocalBoundaryCondition2D<T,DESCRIPTOR>();

    cylinderSetup(lattice, parameters, *boundaryCondition);

    // Main loop over time iterations.
    for (plint iT=0; iT*parameters.getDeltaT()<maxT; ++iT) {
        
        if (iT%parameters.nStep(imSave)==0){
            pcout << "Saving Gif ..." << endl;
            writeGif(lattice, iT);
        }

        if (iT%parameters.nStep(vtkSave)==0 && iT>0) {
            pcout << "Saving VTK file ..." << endl.;
            writeVTK(lattice, parameters, iT);
        }

        lattice.collideAndStream();

        if (iT%parameters.nStep(logT)==0) {
            pcout << "; av energy ="
                  << setprecision(10) << getStoredAverageEnergy<T>(lattice)
                  << "; av rho ="
                  << getStoredAverageDensity<T>(lattice) << endl;
        }
    }

    delete boundaryCondition;Box2D( , , , )
}
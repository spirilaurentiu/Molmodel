/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simbody(tm)                         *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2009 Stanford University and the Authors.           *
 * Authors: Christopher Bruns                                                 *
 * Contributors:                                                              *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */

#include "SimTKmolmodel.h"
#include "molmodel/internal/NoseHooverThermostat.h"
#include "molmodel/internal/AtomSubsystem.h"
#include "molmodel/internal/VanDerWaalsForce.h"

#include "SimTKcommon/Testing.h"

// define SHOW_VIZ for visualisation (for debugging)
// #define SHOW_VIZ 1

// Set which integrator you want to use
#define INTEGRATOR RungeKuttaMersonIntegrator
//#define INTEGRATOR VerletIntegrator
#define SETACCURACY(integ) // use default
//#define SETACCURACY(integ) integ.setAccuracy(1e-2)

using namespace SimTK;
using namespace std;

#define ASSERT(cond) {SimTK_ASSERT_ALWAYS(cond, "Assertion failed");}


class HarmonicOscillator
{
public:
    class OscillatorReporter : public PeriodicEventReporter 
    {
    public:
        mutable int eventCount;
        mutable Real sumEnergy;
        mutable Real sumEnergySquared;
        mutable Vec3 sumVelocity;
        mutable Real sumSpeedSq;
        mutable Real sumSpeed;
        mutable Vec3 sumPosition;
        mutable Real sumRMSVelPos;

        OscillatorReporter(HarmonicOscillator& oscillator, Real reportInterval) 
            : PeriodicEventReporter(reportInterval), oscillator(oscillator),
              eventCount(0), sumEnergy(0.0), sumEnergySquared(0.0), 
              sumVelocity(0), sumSpeedSq(0), sumSpeed(0),
              sumPosition(0.0), sumRMSVelPos(0.0)
        {}

        void handleEvent(const State& state) const 
        {
            // Equilibrate a bit before collecting data
            if (state.getTime() <= 10.0)
                return;

            eventCount ++;
            oscillator.updSystem().realize(state, Stage::Dynamics);

            Real energy = oscillator.getSystem().calcKineticEnergy(state);
            sumEnergy += energy;
            sumEnergySquared += energy*energy;

            Vec3 position = oscillator.getPosition(state);
            sumPosition += position;

            Vec3 velocity = oscillator.getVelocity(state);
            sumVelocity += velocity;
            const Real speedSq = ~velocity*velocity;
            sumSpeed += std::sqrt(speedSq);
            sumSpeedSq += speedSq;
        }
    private:
        HarmonicOscillator& oscillator;

    };

    HarmonicOscillator(int dim) :
        nDim(dim), system(), matter(system), forces(system), atoms(system), mass(1.0)
    {
        Vec3 station(0.0);

        // Use a slider to constrain the oscillator to one dimension of motion
        if (nDim == 1) {
            MobilizedBody::Slider atomBody( 
                                    matter.updGround(), Transform(),
                                    Body::Rigid(MassProperties(mass, station, Inertia(1))),
                                    Vec3(2,0,0));
            sliderIndex = atomBody.getMobilizedBodyIndex();
            Force::MobilityLinearSpring(forces, atomBody, 0, 10.0, 2.0);
        } else if (nDim == 2) {
            MobilizedBody::Slider dummy(matter.updGround(), Rotation(Pi/2,ZAxis),
                                    MassProperties(0,Vec3(0),Inertia(0)), 
                                    Rotation(Pi/2,ZAxis));
            MobilizedBody::Slider atomBody( 
                                    dummy, Transform(),
                                    Body::Rigid(MassProperties(mass, station, Inertia(1))),
                                    Vec3(2,2,0));
            sliderIndex = atomBody.getMobilizedBodyIndex();
            Force::MobilityLinearSpring(forces, dummy,    0, 10.0, 2.0);
            Force::MobilityLinearSpring(forces, atomBody, 0, 10.1, 2.0);
        } else {
            MobilizedBody::Cartesian atomBody( 
                                    matter.updGround(), Transform(),
                                    Body::Rigid(MassProperties(mass, station, Inertia(1))),
                                    Vec3(2,2,2));
            sliderIndex = atomBody.getMobilizedBodyIndex();
            Force::MobilityLinearSpring(forces, atomBody, 0, 10.0, 2.0);
            Force::MobilityLinearSpring(forces, atomBody, 1, 10.1, 2.0);
            Force::MobilityLinearSpring(forces, atomBody, 2, 10.2, 2.0);
        }

        AtomSubsystem::AtomIndex atomIx = 
                atoms.addAtom(mass);
        atoms.updAtom(atomIx).setStationInBodyFrame(station);
        atoms.setAtomMobilizedBodyIndex(atomIx, sliderIndex);
    }

    void simulate() {
        // Visualize - for debugging only
#ifdef SHOW_VIZ
        Visualizer viz(system);
        viz.setBackgroundType(Visualizer::SolidColor);
        system.addEventReporter( new Visualizer::Reporter(viz, 0.2) );
#endif

        reporter = new OscillatorReporter(*this, 0.25);
        system.addEventReporter(reporter);

        State state = system.realizeTopology();
        Random::Uniform randVel(-1e-4,1e-4);
        for (int i=0; i < state.getNU(); ++i)
            state.updU()[i] = randVel.getValue();

        // Simulate it.
        
        INTEGRATOR integ(system);
        SETACCURACY(integ);
        TimeStepper ts(system, integ);
        ts.initialize(state);
        ts.stepTo(1000.0);
    }

    void assertTemperature(Real temperature) const 
    {
        // ensure we collected some data
        SimTK_TEST(reporter->eventCount > 100);

        const int nAtoms = 1;
        const int nExcluded = 0;
        const int nThermalDofs = nDim - nExcluded; // 1, 2 or 3
        const int nAtomSamples = nAtoms * reporter->eventCount; 

        // Nicer math names.
        const Real T = temperature;
        const Real k = SimTK_BOLTZMANN_CONSTANT_MD;
        const int  N = nThermalDofs;
        const Real M = nAtoms * mass;
        const int  d = nDim;

        // Sanity checks

        // Mean position should be roughly zero
        const Vec3 expectedMeanPosition = Vec3(0);
        const Vec3 measuredMeanPosition = reporter->sumPosition / nAtomSamples;
        cout << "meanPos=" << measuredMeanPosition.norm() << endl;
        SimTK_TEST_EQ_TOL(expectedMeanPosition, measuredMeanPosition, 0.1);

        // Mean velocity should be roughly zero
        const Vec3 expectedMeanVelocity = Vec3(0);
        const Vec3 measuredMeanVelocity = reporter->sumVelocity / nAtomSamples;
        cout << "meanVel=" << measuredMeanVelocity.norm() << endl;
        SimTK_TEST_EQ_TOL(expectedMeanVelocity, measuredMeanVelocity, 0.1);

        // Check temperature
        const Real measuredMeanEnergy = reporter->sumEnergy/reporter->eventCount;
        const Real expectedMeanEnergy = N*k*T/2; // kT/2 per degree of freedom
        cout << "energy err=" << std::abs(1.0 - measuredMeanEnergy/expectedMeanEnergy) << endl;
        SimTK_TEST_EQ_TOL(expectedMeanEnergy, measuredMeanEnergy, 0.2);

        // Boltzmann distribution stuff

        // Mean squared speed should be v2bar=N*kT/M where N is number of thermal dofs and
        // M is the total mass. Usually you see this as 3*kT/m because N=3*natoms
        // and M=m*natoms. N is typically defined to have rigid body dofs removed
        // (up to 6) so really 3kT/m is the limit as natoms->infinity. You'll only
        // notice the different for small numbers of atoms, but doesn't everyone start
        // with those problems?

        const Real expectedMeanSpeedSq = N*k*T/M;
        const Real measuredMeanSpeedSq = reporter->sumSpeedSq / nAtomSamples;
        cout << "v^2 err=" << std::abs(1.0 - measuredMeanSpeedSq/expectedMeanSpeedSq) << endl;
        SimTK_TEST_EQ_TOL(expectedMeanSpeedSq, measuredMeanSpeedSq, 0.2);

        // OK, the mean speed (not speed squared) vbar is a LOT harder to come by! For one
        // thing a 1-chain Nose'-Hoover won't even get this right. 2-chain should work
        // and ours defaults to 3 so should work well. You might see this written as
        // vbar=sqrt(8*kT/(m*pi)) but that is really the natoms->infinity limit of
        // vbar=sqrt(8*NkT/(3*M*pi)) with N=3*natoms-6 and M=m*natoms.
        // When the system dimensionality (that is, # dofs per particle) is not 3,
        // the above formula is not correct. The formulas are
        //      vbar1d = sqrt(    2/pi * NkT/M )            (1)
        //      vbar2d = sqrt(   pi/4  * NkT/M )            (2)
        //      vbar3d = sqrt((8/3)/pi * NkT/M )            (3)
        //
        // DERIVATION
        // These can all be derived from the velocity vector probability density 
        // function which has this form for each coordinate (these are multiplied
        // together to give probability densities for multi-dof vectors).
        //      fv(x)=sqrt(m/(pi*2kT)) * exp(-m*x^2/(2kT)) (similarly for y and z;
        // and coordinates are signed quantities). That is
        //      fv1(x)     = fv(x)
        //      fv2(x,y)   = fv(x)*fv(y)
        //      fv3(x,y,z) = fv(x)*fv(y)*fv(z)
        // Then these are converted to speed v by integrating over constant speed 
        // sphere in 3d (area 4pi*v^2), constant speed circle in 2d (circumference 
        // 2pi*v), and just two points of opposite sign for 1d ("area" is just 2). 
        // So the *speed* probability density functions are
        //      f1(v) = 2       * sqrt(m/(pi*2kT))   * exp(-m*v^2/(2kT))
        //      f2(v) = 2pi*v   * sqrt(m/(pi*2kT))^2 * exp(-m*v^2/(2kT))
        //      f3(v) = 4pi*v^2 * sqrt(m/(pi*2kT))^3 * exp(-m*v^2/(2kT))
        // where v is the scalar speed in each case. The we get the average
        // speed by integrating v*f1(v), v*f2(v), v*f3(v) from v=0 to v=Inf.
        // That gives the following formulas for vbar:
        //      vbar1d = sqrt( 2/pi * kT/m )
        //      vbar2d = sqrt( pi/2 * kT/m )
        //      vbar2d = sqrt( 8/pi * kT/m )
        // Those are natoms->infinity limits; substituting the actual number of
        // thermal dofs N=d*natoms - excluded we get the formulas (1)-(3) above.
        //
        // (sherm 091211)

        const Real c = (d==1 ? 2/Pi : (d==2 ? Pi/4 : 8/(3*Pi)));
        const Real expectedMeanSpeed = std::sqrt( c*N*k*T/M );
        const Real measuredMeanSpeed = reporter->sumSpeed / nAtomSamples;
        cout << "v err=" << std::abs(1.0 - measuredMeanSpeed/expectedMeanSpeed) << endl;
        SimTK_TEST_EQ_TOL(expectedMeanSpeed, measuredMeanSpeed, 0.2);
    }

    MultibodySystem& updSystem() {return system;}
    const MultibodySystem& getSystem() const {return system;}

    SimbodyMatterSubsystem& updMatterSubsystem() {return matter;}
    const SimbodyMatterSubsystem& getMatterSubsystem() const {return matter;}

    GeneralForceSubsystem& updForceSubsystem() {return forces;}
    const GeneralForceSubsystem& getForceSubsystem() const {return forces;}

    AtomSubsystem& updAtomSubsystem() {return atoms;}
    const AtomSubsystem& getAtomSubsystem() const {return atoms;}

    // slider coordinate
    Vec3 getPosition(const State& state) const {
        return matter.getMobilizedBody(sliderIndex).getBodyOriginLocation(state);
    }

    Vec3 getVelocity(const State& state) const {
        return matter.getMobilizedBody(sliderIndex).getBodyOriginVelocity(state);
    }

    Real getTime(const State& state) const {
        return state.getTime();
    }

private:
    const int               nDim;
    MultibodySystem         system;
    SimbodyMatterSubsystem  matter;
    GeneralForceSubsystem   forces;
    AtomSubsystem           atoms;

    const units::md::mass_t mass;
    MobilizedBodyIndex      sliderIndex;

    OscillatorReporter*     reporter; // just a reference; don't delete
};

class ArgonGasSphere
{
public:

    class ArgonReporter : public PeriodicEventReporter 
    {
    public:
        mutable int eventCount;
        mutable Real sumEnergy;
        mutable Real sumEnergySquared;
        mutable Vec3 sumVelocity;
        mutable Real sumSpeedSq;
        mutable Real sumSpeed;
        mutable Vec3 sumPosition;
        mutable int  nDofs, nConstraints;

        ArgonReporter(ArgonGasSphere& argon, Real reportInterval) 
            : PeriodicEventReporter(reportInterval), argon(argon),
              eventCount(0), sumEnergy(0.0), sumEnergySquared(0.0),
              sumVelocity(0), sumSpeedSq(0), sumSpeed(0), sumPosition(0)
        {}

        void handleEvent(const State& state) const 
        {
            const MultibodySystem& system = argon.getSystem();
            const SimbodyMatterSubsystem& matter = system.getMatterSubsystem();

            if (state.getTime()==0) {
                nDofs = state.getNU();
                nConstraints = state.getNUDotErr();
            }

            // Equilibrate a bit before collecting data
            if (state.getTime() <= 10.0)
                return;

            eventCount ++;
            system.realize(state, Stage::Dynamics);

            Real energy = system.calcKineticEnergy(state);
            sumEnergy += energy;
            sumEnergySquared += energy*energy;

            // Sum all the atom speeds and speeds squared.
            for (MobilizedBodyIndex mbx(1); mbx < matter.getNumBodies(); ++mbx) {
                const MobilizedBody& atom = matter.getMobilizedBody(mbx);
                const Vec3& p = atom.getBodyOriginLocation(state);
                const Vec3& v = atom.getBodyOriginVelocity(state);
                sumPosition += p; sumVelocity += v;
                const Real v2 = ~v*v;
                sumSpeedSq += v2;
                sumSpeed += std::sqrt(v2);
            }
        }
    private:
        ArgonGasSphere& argon;

    };

    static Real getArgonMass() {return 39.948 * md::daltons;}

    ArgonGasSphere(int numberOfAtoms) :
        system(), matter(system), forces(system), atoms(system)
    {
        Vec3 station(0.0 * md::angstroms);
        Real atomMass = getArgonMass();

        VanDerWaalsForce* vdwForce = new VanDerWaalsForce(atoms);
        Force::Custom(forces, vdwForce);

        // Place atoms on a cubic lattice to start
        Real latSpacing = 5.0 * md::angstroms; // distance between atoms

        // How many atoms per edge of the lattice?
        Real l1 = std::pow(numberOfAtoms, 0.333333);
        int atomsPerLatticeEdge = (int)( std::ceil(std::pow(numberOfAtoms, 0.33333)) ); // max edge dimension

        int latX = 0;
        int latY = 0;
        int latZ = 0;
        Vec3 latCenter = Vec3(0.5, 0.5, 0.5) * latSpacing * (atomsPerLatticeEdge - 1);

        // Add a little randomness to the atom positions
        Random::Uniform jiggle( -0.5 * md::angstroms, 0.5 * md::angstroms );

        for (int a = 0; a < numberOfAtoms; ++a) 
        {
            MobilizedBody::Cartesian body(
                matter.updGround(),
                Body::Rigid( MassProperties(atomMass, station, Inertia(1)) )
            );

            Vec3 jiggleVec = Vec3( jiggle.getValue(), jiggle.getValue(), jiggle.getValue() );
            Vec3 location = Vec3(latX, latY, latZ) * latSpacing - latCenter + jiggleVec;
            body.setDefaultQ( location );

            AtomSubsystem::AtomIndex atomIx = atoms.addAtom(atomMass);
            atoms.updAtom(atomIx).setStationInBodyFrame(station);
            atoms.setAtomMobilizedBodyIndex(atomIx, body.getMobilizedBodyIndex());

            // forces
            vdwForce->addAtom(atomIx, 1.88 * md::angstroms, 0.23725 * md::kilocalories_per_mole);

            // Update lattice position for next atom
            ++latX;
            if (latX >= atomsPerLatticeEdge) {
                latX = 0;
                ++latY;
                if (latY >= atomsPerLatticeEdge) {
                    latY = 0;
                    ++latZ;
                }
            }
        }

        // This did seem to work.
        //if (numberOfAtoms==3)
        //    for (MobilizedBodyIndex mbx(1); mbx < numberOfAtoms; ++mbx) 
        //        Constraint::Rod(matter.updMobilizedBody(mbx),matter.updMobilizedBody(MobilizedBodyIndex(mbx+1)), 1.);

        matter.setShowDefaultGeometry(false);

        // (sherm 091210) I removed this because the system should be isolated
        // from ground in order to conserve linear and angular momentum.

        //vdwForce->addBoundarySphere(
        //    Vec3(0,0,0), 
        //    latSpacing * atomsPerLatticeEdge + 1.88 * md::angstroms,
        //    1.88 * md::angstroms,
        //    0.23725 * md::kilocalories_per_mole);

#ifdef SHOW_VIZ
        DecorationSubsystem decorations(system);
        vdwForce->decorateAtoms(decorations);
        vdwForce->decorateBoundingSpheres(decorations);
#endif

    }

    void simulate(const NoseHooverThermostat& thermo) {
        // Visualize - for debugging only
#ifdef SHOW_VIZ
        Visualizer viz(system);
        viz.setBackgroundType(Visualizer::SolidColor);
        system.addEventReporter( new Visualizer::Reporter(viz, 0.20) );
#endif

        reporter = new ArgonReporter(*this, 0.25);
        system.addEventReporter(reporter);

        State state = system.realizeTopology();

        Random::Uniform randVels(-1e-4,1e-4);
        for (int i=0; i < state.getNU(); ++i)
            state.updU()[i] = randVels.getValue();

        // Nuke COM pos and velocity (masses are all the same). There are expected
        // to be conserved in this free flying system.
        system.realize(state, Stage::Velocity);
        Vec3 com(0);
        Vec3 comv(0);
        int nAtoms = 0;
        for (MobilizedBodyIndex mbx(1); mbx < matter.getNumBodies(); ++mbx) {
            const MobilizedBody& body = matter.getMobilizedBody(mbx);
            const Vec3& p = body.getBodyOriginLocation(state);
            const Vec3& v = body.getBodyOriginVelocity(state);
            com += p;
            comv += v;
            ++nAtoms;
        }
        const Vec3 comErr = com / nAtoms;
        const Vec3 comvErr = comv / nAtoms;
        //cout << "BEFORE com=" << comErr.norm() << " comv=" << comvErr.norm() << endl;
        for (MobilizedBodyIndex mbx(1); mbx < matter.getNumBodies(); ++mbx) {
            const MobilizedBody& body = matter.getMobilizedBody(mbx);
            const MobilizedBody::Cartesian atom = MobilizedBody::Cartesian::downcast(body);
            atom.setU(state, atom.getU(state) - comvErr);
            atom.setQ(state, atom.getQ(state) - comErr);
            system.realize(state, Stage::Velocity);
        }

        com = 0;
        comv = 0;
        for (MobilizedBodyIndex mbx(1); mbx < matter.getNumBodies(); ++mbx) {
            const MobilizedBody& body = matter.getMobilizedBody(mbx);
            const Vec3& p = body.getBodyOriginLocation(state);
            const Vec3& v = body.getBodyOriginVelocity(state);
            com += p;
            comv += v;
        }
        //cout << "AFTER com=" << (com/nAtoms).norm() << " comv=" << (comv/nAtoms).norm() <<endl;

        // Simulate it.

        cout << "N-H says #dofs is " << thermo.getNumThermalDofs(state) << endl;
        
        INTEGRATOR integ(system);
        SETACCURACY(integ);
        TimeStepper ts(system, integ);
        ts.initialize(state);
        ts.stepTo(1000.0);
    }

    // SEE assertTemperature() in Oscillator above for comments to explain this.
    void assertTemperature(Real temperature) const 
    {
        // ensure we collected some data
        SimTK_TEST(reporter->eventCount > 100);

        const int nAtoms = getAtomSubsystem().getNumAtoms();
        const int nDofs = reporter->nDofs;
        const int nConstraints = reporter->nConstraints;
        SimTK_TEST(nDofs == 3*nAtoms);
        const int nExcluded = nAtoms==1 ? 3 : (nAtoms==2?5:6);
        const int nThermalDofs = (nDofs-nConstraints) - nExcluded;
        const int nDim = std::min(nThermalDofs, 3); // 1, 2, or 3
        const int nAtomSamples = nAtoms * reporter->eventCount; 

        cout << "Calculated nThermalDofs=" << nThermalDofs << endl;

        // Nicer math names.
        const Real T = temperature;
        const Real k = SimTK_BOLTZMANN_CONSTANT_MD;
        const int  N = nThermalDofs;
        const Real M = nAtoms * getArgonMass();
        const int  d = nDim;

        // Sanity checks

        // Mean position should be zero
        const Vec3 expectedMeanPosition = Vec3(0);
        const Vec3 measuredMeanPosition = reporter->sumPosition / nAtomSamples;
        cout << "meanPos=" << measuredMeanPosition.norm() << endl;
        SimTK_TEST_EQ_TOL(expectedMeanPosition, measuredMeanPosition, 0.01);

        // Mean velocity should be zero
        const Vec3 expectedMeanVelocity = Vec3(0);
        const Vec3 measuredMeanVelocity = reporter->sumVelocity / nAtomSamples;
        cout << "meanVel=" << measuredMeanVelocity.norm() << endl;
        SimTK_TEST_EQ_TOL(expectedMeanVelocity, measuredMeanVelocity, 0.01);

        // Check temperature
        const Real measuredMeanEnergy = reporter->sumEnergy/reporter->eventCount;
        const Real expectedMeanEnergy = N*k*T/2; // kT/2 per degree of freedom
        cout << "energy err=" << std::abs(1.0 - measuredMeanEnergy/expectedMeanEnergy) << endl;
        SimTK_TEST_EQ_TOL(expectedMeanEnergy, measuredMeanEnergy, 0.2);

        // Boltzmann distribution stuff

        // See comments in the other assertTemperature() routine.
        const Real expectedMeanSpeedSq = N*k*T/M;
        const Real measuredMeanSpeedSq = reporter->sumSpeedSq / nAtomSamples;
        cout << "v^2 err=" << std::abs(1.0 - measuredMeanSpeedSq/expectedMeanSpeedSq) << endl;
        SimTK_TEST_EQ_TOL(expectedMeanSpeedSq, measuredMeanSpeedSq, 0.2);

        // See comments in the other assertTemperature() routine.
        const Real c = (d==1 ? 2/Pi : (d==2 ? Pi/4 : 8/(3*Pi)));
        const Real expectedMeanSpeed = std::sqrt( c*N*k*T/M );
        const Real measuredMeanSpeed =  reporter->sumSpeed / nAtomSamples;
        cout << "v err=" << std::abs(1.0 - measuredMeanSpeed/expectedMeanSpeed) << endl;
        SimTK_TEST_EQ_TOL(expectedMeanSpeed, measuredMeanSpeed, 0.2);
    }

    MultibodySystem& updSystem() {return system;}
    const MultibodySystem& getSystem() const {return system;}

    SimbodyMatterSubsystem& updMatterSubsystem() {return matter;}
    const SimbodyMatterSubsystem& getMatterSubsystem() const {return matter;}

    GeneralForceSubsystem& updForceSubsystem() {return forces;}
    const GeneralForceSubsystem& getForceSubsystem() const {return forces;}

    AtomSubsystem& updAtomSubsystem() {return atoms;}
    const AtomSubsystem& getAtomSubsystem() const {return atoms;}

     Real getTime(const State& state) const {
        return state.getTime();
    }

private:
    MultibodySystem         system;
    SimbodyMatterSubsystem  matter;
    GeneralForceSubsystem   forces;
    AtomSubsystem           atoms;

    ArgonReporter*          reporter; // just a reference; don't delete
};

// The oscillator is Case study 12, page 155 in
// Understanding Molecular Simulation: From Algorithms to Applications
// Frenkel and Smit.
// (sherm 091210) Neither Chris Bruns nor I could find a source for the
// Boltzmann distribution equations we're using in the odd cases where
// the removed rigid body dofs matter and when the problem dimensionality
// is 1 or 2 (even though Frenkel and Smit do the oscillator as a
// 1-dimensional problem). See comments in the oscillator assertTemperature()
// method above.


// This is always 1 particle tethered to Ground, but it can be a 
// 1-dof, 2-dof, or 3-dof particle. When a system interacts with Ground
// there is no guarantee it will conserve momentum so the thermostat should
// be told not to drop the 6 rigid body dofs that it normally likes to do.
void testOscillatorTemperature(int dim, Real temperature)
{
    SimTK_ASSERT_ALWAYS(1<=dim && dim <=3, "bad universe");

    cout << dim << "-d OSCILLATOR T=" << temperature << endl;

    HarmonicOscillator oscillator(dim);
    GeneralForceSubsystem& forces = oscillator.updForceSubsystem();
    // 1ps relaxation time; no excluded dofs
    NoseHooverThermostat thermo(forces, oscillator.updMatterSubsystem(), temperature, 1, 0);
    //thermo.setDefaultNumChains(1);   // won't conserve momentum with just 1 chain
    oscillator.simulate();
    oscillator.assertTemperature(temperature);
}

// This is a free-flying system so does conserve momentum. However, for 2 atoms
// there are only 5 rigid body momenta conserved -- check McQuarrie's Statistical
// Mechanics for the case of a linear molecule if you don't believe me!
void testArgonTemperature(int nArgons, Real temperature)
{
    SimTK_ASSERT_ALWAYS(nArgons>0, "come on, don't waste my time");
    const int nExcluded = (nArgons==1?3:(nArgons==2?5:6));

    cout << nArgons << " ARGONS at T=" << temperature 
          << " (excluding " << nExcluded << " dofs)\n";

    ArgonGasSphere argon(nArgons);
    GeneralForceSubsystem& forces = argon.updForceSubsystem();
    NoseHooverThermostat thermo(forces, argon.updMatterSubsystem(), temperature, 1, nExcluded); 
    //thermo.setDefaultNumChains(1);   // won't conserve momentum with just 1 chain
    argon.simulate(thermo);
    argon.assertTemperature(temperature); 
}

int main() 
{
    SimTK_START_TEST("TestNoseHooverThermostat");

        SimTK_SUBTEST2(testOscillatorTemperature, 1, 50.0);
        SimTK_SUBTEST2(testOscillatorTemperature, 2, 50.0);
        SimTK_SUBTEST2(testOscillatorTemperature, 3, 50.0);
        // Skip these if in Debug mode
        #ifdef NDEBUG
            SimTK_SUBTEST2(testOscillatorTemperature, 1, 300.0);
            SimTK_SUBTEST2(testOscillatorTemperature, 3, 1000.0);
            SimTK_SUBTEST2(testArgonTemperature, 9, 1000.0);
        #endif

        // 1 dof should do nothing but doesn't work currently
        //SimTK_SUBTEST2(testArgonTemperature, 1, 300.0);
        SimTK_SUBTEST2(testArgonTemperature, 2, 300.0);
        SimTK_SUBTEST2(testArgonTemperature, 3, 300.0);
        SimTK_SUBTEST2(testArgonTemperature, 8, 300.0);
        // Skip this if in Debug mode
        #ifdef NDEBUG
            SimTK_SUBTEST2(testArgonTemperature, 17, 20.0);
        #endif
    SimTK_END_TEST();
}

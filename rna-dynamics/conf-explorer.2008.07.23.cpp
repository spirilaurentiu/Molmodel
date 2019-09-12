/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Molmodel                            *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2006-7 Stanford University and the Authors.         *
 * Authors: Samuel Flores, Christopher Bruns                                  *
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
#include "SimTKsimbody_aux.h"

//#include "/Users/samuelflores/rna-dynamics/Ligands.h"  
//#include "/Users/samuelflores/rna-dynamics/Water.h"    
//#include "/Users/samuelflores/rna-dynamics/WaterDroplet.h"

#include "Ligands.h"  
#include "Water.h"    
//#include "/home/scflores/rna-scratch/Water.h"    
#include "WaterDroplet.h"
//#include "/home/scflores/rna-scratch/WaterDroplet.h"

//#include "/Users/samuelflores/svn/molmodel/include/molmodel/internal/Ligands.h"
#include <iostream>
#include <vector>

using std::cout;
using std::endl;

using namespace SimTK;
using namespace std;

// define CREATE_VTK_WINDOW to see animated vtk window of simulation
// undefine for automated nightly builds
// #define CREATE_VTK_WINDOW
// Rubber


static const Real rubber_density = 1100.;  // kg/m^3
static const Real rubber_young   = 0.01e9; // pascals (N/m)
static const Real rubber_poisson = 0.5;    // ratio
static const Real rubber_planestrain =
                    rubber_young/(1.-rubber_poisson*rubber_poisson);
static const Real rubber_dissipation = 0.005;
// /

int main(int argc, char *argv[]) {
    int firstres=17 ;

    CompoundSystem system;
    SimbodyMatterSubsystem  matter(system);
    TinkerDuMMForceFieldSubsystem dumm(system);
    //DecorationSubsystem     artwork(system);
    //HuntCrossleyContact     contact(system);

    char* inFileName;
    char* myBondMobilityType;
    //Real timeInterval = 0.01;
    float myCoulombGlobalScaleFactor = 1.0;
    float myGbsaGlobalScaleFactor = 0;
    float myVdwGlobalScaleFactor = 1.0;
    float myReportingInterval    =0.01;
    float maxTime                =10000;//ps
    float thermostatTimeConstant = 0.0;
    int myUseMultithreadedComputation =1; 
    float myTemperature = 290;
    float waterRadius   = 3.50;
    int numWater         =0   ;
    int numDivalents     = 0  ;
    int makeWaterDroplet=1;	


 for  (int q =1; q<argc;q+=2)
    {   
        cout<<"argv["<<q<<"]  = "<<argv[q]<<endl;
        string key = argv[q];
        if      (key=="-T") {myTemperature = atof(argv[q+1]);}
        else if (key     == "-I") {myReportingInterval= atof(argv[q+1]);}
        else if (key     == "-S") {inFileName = argv[q+1];}
        else if (key     == "-V") {myVdwGlobalScaleFactor 	 = atof(argv[q+1]);}
        else if (key     == "-G") {myGbsaGlobalScaleFactor	 =atof(argv[q+1]);}
        else if (key     == "-C") {myCoulombGlobalScaleFactor    = atof(argv[q+1]);}
        else if (key     == "-M") {maxTime                       = atof(argv[q+1]);}
        else if (key     == "-TC") {thermostatTimeConstant       = atof(argv[q+1]);}
        else if (key     == "-MT") {myUseMultithreadedComputation=atoi(argv[q+1]);}
        else if (key     == "-ND") {numDivalents                 =atoi(argv[q+1]);}
        else if (key     == "-NW") {numWater                     =atoi(argv[q+1]);}
        else if (key     == "-WD") {makeWaterDroplet             =atoi(argv[q+1]);}
        else if (key     == "-WR") {waterRadius                  =atof(argv[q+1]);}
        else if (key     == "-BM") {myBondMobilityType           =(argv[q+1]);}
    }   

    dumm.setCoulombGlobalScaleFactor(myCoulombGlobalScaleFactor);	
    dumm.setGbsaGlobalScaleFactor(myGbsaGlobalScaleFactor);	
    dumm.setVdwGlobalScaleFactor(myVdwGlobalScaleFactor);	
    dumm.setUseMultithreadedComputation(myUseMultithreadedComputation);



    //contact.addHalfSpace(matter.Ground(), UnitVec3(0,1,0), -3, rubber_planestrain, rubber_dissipation);
    //artwork.addBodyFixedDecoration(GroundIndex, Transform(Vec3(0, -3, 0)), DecorativeBrick(Vec3(7,.02,7)).setColor(Yellow).setOpacity(0.25));

    //contact.addHalfSpace(matter.Ground(), UnitVec3(1,0,0), -3, rubber_planestrain, rubber_dissipation);
    //artwork.addBodyFixedDecoration(GroundIndex, Transform(Vec3(-3, 0, 0)), DecorativeBrick(Vec3(.02, 7 ,7)).setColor(Yellow).setOpacity(0.25));

    //contact.addHalfSpace(matter.Ground(), UnitVec3(1,0,0), 3, rubber_planestrain, rubber_dissipation);
    //artwork.addBodyFixedDecoration(GroundIndex, Transform(Vec3(3, 0, 0)), DecorativeBrick(Vec3(.02,7,7)).setColor(Yellow).setOpacity(0.25));

    //ifstream tinkerStream("/home/scflores/svn/molmodel/resources/tinker_amber99_sam.prm");
    ifstream tinkerStream("tinker_amber99_sam.prm");
    dumm.populateFromTinkerParameterFile(tinkerStream);
    tinkerStream.close();
    P12      methanol1(dumm);//, methanol2, methanol3;
    //Water    myWaterVec(dumm);//[numWater];
    Water * myWaterVec[numWater];
    for (int i=0; i<numWater; i++)  {
		myWaterVec[i]=new Water(dumm);
    }
    MagnesiumIon myMagnesiumIonVec[numDivalents]; 
    for (int i = 0; i < numDivalents; i++) {
	(myMagnesiumIonVec[i]).setAmberLikeParameters(dumm);
	myMagnesiumIonVec[i].setPdbResidueNumber(i);
	myMagnesiumIonVec[i].setPdbChainId('C');    
	myMagnesiumIonVec[i].setPdbResidueName("MG2");
	
    }	
    //myWaterVec.setPdbResidueName("H2O");
    //dumm.loadAmber99Parameters();
    //RNA myMolecule("G");//GCAGAUCUGAGCCUGGGAGCUCUCUGCC");
    RNA myMolecule("GGCAGAUCUGAGCCUGGGAGCUCUCUGCC");
    myMolecule.assignBiotypes();


    for (int q = 0; q <  (myMolecule.getNResidues()); q++)
         {   
          myMolecule.updResidue(Compound::Index(q)).setPdbChainId('B');
          cout<<"myMolecule.updResidue(Compound::Index(q)).setPdbResidueNumber(firstres+q);"<<std::endl;
          myMolecule.updResidue(Compound::Index(q)).setPdbResidueNumber(firstres+q);
        }   
    std::ifstream inFileStream("1ARJ.simple.pdb",ifstream::in);
    char * qvectorfilename = "qvector.dat";
   // ofstream qvectorstream;
    //qvectorstream.open(qvectorfilename);
    	

    ifstream qvectorstream;
    qvectorstream.open(qvectorfilename);
    int r=0;	
    if (qvectorstream.is_open()) {
	while(qvectorstream.good()) {
		float myq = qvectorstream.get();
		cout<<"q ("<<r<<") = "<<endl;
		r++;
		}
        qvectorstream.close();		
	} //else {	
    PdbStructure pdbStructure(inFileStream);
    Compound::AtomTargetLocations atomTargets = myMolecule.createAtomTargets(pdbStructure);
    std::cout<<"atomtargest.szie "<<atomTargets.size()<<"versus atoms in myMolecule = "<<myMolecule.getNAtoms()<<std::endl;    


    // Four steps to a perfect match
    myMolecule.matchDefaultAtomChirality(atomTargets);
    myMolecule.matchDefaultBondLengths(atomTargets);
    myMolecule.matchDefaultBondAngles(atomTargets);
    myMolecule.matchDefaultDihedralAngles(atomTargets);
    myMolecule.matchDefaultTopLevelTransform(atomTargets);
    Real residual = myMolecule.getTransformAndResidual(atomTargets).residual;
	//while (r < 
	//}
    for (int q = 0; q <  (myMolecule.getNResidues()); q++)  
		if (!((q<(26-17)) &&(q>(22-17)))) {
                  if ((q>0)&&(q<(myMolecule.getNResidues()))) {
                        stringstream ss1;
                        ss1<<q-1<<"/O3*";
                        stringstream ss2;
                        ss2<<q<<"/P";
                        myMolecule.setBondMobility(BondMobility::Rigid, ss1.str()    ,ss2.str()  );
                        //myMolecule.setBondMobility(BondMobility::Free , ss1.str()    ,ss2.str()  );
                        cout<<"rigidifying: "<< ss1.str()  <<" "<< ss2.str() <<endl;
                        }
                  for (int r =0 ; r<myMolecule.getResidue(Compound::Index(q)).getNBonds(); r++)    //)//(Compound::Index(q)).
                        {
                        //myMolecule.updResidue(Compound::Index(q)).setBondMobility(BondMobility::Free ,Compound::BondIndex(r));
                        myMolecule.updResidue(Compound::Index(q)).setBondMobility(BondMobility::Rigid,Compound::BondIndex(r));
                        }
                }

/*
    // hydroxyl hydrogen
    DuMM::AtomClassIndex amberHOAtomClassIndex(31);
    DuMM::ChargedAtomTypeIndex methanolHOAtomTypeIndex(4003);
    dumm.defineChargedAtomType(methanolHOAtomTypeIndex, "Methanol OH", amberHOAtomClassIndex, 0.4); // from serine OH
    dumm.setBiotypeChargedAtomType( methanolHOAtomTypeIndex, Biotype::get("Methanol", "HO").getIndex() );
*/

    //scf put this ligand back in later
    system.adoptCompound(methanol1, Vec3(-0.5, 2, 0));
    for (int i=0; i<numWater; i++)  {
    	(*myWaterVec[i]).setPdbResidueName("H2O");
    	(*myWaterVec[i]).setPdbResidueNumber(i);
    	(*myWaterVec[i]).setPdbChainId('D');    
	cout<< i<<","<< 10.0*sin(1.0*i/30*360*Deg2Rad )<<"," <<10.0*cos(1.0*i/30*360*Deg2Rad )<<endl;
    	system.adoptCompound(*myWaterVec[i],Vec3(   1.0*i/200.0-2.000, 3 *cos(1.0*i/50 *360*Deg2Rad ), 3 *sin(1.0*i/50 *360*Deg2Rad )));
    }	

    for (int i = 1; i<numDivalents; i++)
    {
   	system.adoptCompound(myMagnesiumIonVec[i],Vec3(i,0,0));
    }	 	
    //system.adoptCompound(methanol3, Vec3( 0.5, 0, 0));
    system.adoptCompound(myMolecule);

	
    //system.modelCompounds();        
    //State state = system.realizeTopology();
    GeneralForceSubsystem forces; //this is needed only for the VanderWallSphere in WaterDroplet
    //scf put back in later
    if (makeWaterDroplet) WaterDroplet myReturnInt(system,dumm,forces);	    
    //WaterDroplet myReturnInt(system,state,dumm);	    
    cout<<"check 1.0"<<endl;
    system.modelCompounds();        
    cout<<"check 1.5"<<endl;
    State state = system.realizeTopology();
    system.realize(state,Stage::Position);

    cout<<"check 2.0"<<endl;
    //const SimTK::Transform& mytransform2 = matter.getMobilizedBody(myMolecule.getAtomMobilizedBodyIndex(Compound::AtomIndex(0))).findBodyTransformInAnotherBody(state,matter.getMobilizedBody(myMolecule.getAtomMobilizedBodyIndex(Compound::AtomIndex(myMolecule.getNAtoms()-1))));
    cout<<"check 2.5"<<endl;
    //Constraint::Weld myConstraint(matter.updMobilizedBody(myMolecule.getAtomMobilizedBodyIndex(Compound::AtomIndex(0))),Transform(Vec3(0)) , matter.updMobilizedBody(myMolecule.getAtomMobilizedBodyIndex(Compound::AtomIndex(myMolecule.getNAtoms()-1))),mytransform2      );
    cout<<"check 3.0"<<endl;
/*
    ifstream qvectorstreamin ;
    qvectorstreamin.open(qvectorfilename);
    for (int s = 0 ; s < state.getNQ(); s++){//;//.	
	qvectorstreamin  >>      state.updQ()[s];//       <<endl;
//	cout <<      state.updQ()[s]       <<endl;
	}
    qvectorstreamout.close();//Vector myqVector = myMolecule.getQAsVector(state);

    ofstream qvectorstreamout;
    qvectorstreamout.open(qvectorfilename);
    for (int s = 0 ; s <      state.getNQ(); s++){//;//.	
	qvectorstreamout <<      state.updQ()[s]       <<endl;
	cout <<      state.updQ()[s]       <<endl;
	//cout<<myMolecule.getOneQ(state,s)<<endl;

	}
    qvectorstreamout.close();	
*/

#ifdef CREATE_VTK_WINDOW
    VTKVisualizer display(system, 0.1);
#endif

    VelocityRescalingThermostat * myVelocityRescalingThermostat = new VelocityRescalingThermostat(system,  myTemperature, myReportingInterval);
    ofstream output("mymovie.pdb"    );
    PeriodicPdbWriter * myPeriodicPdbWriter = new PeriodicPdbWriter(system,output,myReportingInterval);
    cout<<"check 4.0"<<endl;
    system.updDefaultSubsystem().addEventHandler(myVelocityRescalingThermostat);
    system.updDefaultSubsystem().addEventReporter(myPeriodicPdbWriter);
    cout<<"check 5.0"<<endl;





    RungeKuttaMersonIntegrator study(system);
    cout<<"check 6.0"<<endl;
    study.initialize(state);

    time_t rawtime;
    struct tm * timeinfo;


    //scf put back in later
    cout<<"about to start minimizing"<<endl;
    //LocalEnergyMinimizer::minimizeEnergy(system,state,15.0);

    system.realizeTopology();

    int p = 0;
    stringstream ss1;
    ss1<<"out."<<p<<".pdb";
    ofstream outputFrame(ss1.str().c_str());

    cout<<"p = "<<p<<endl;
    for (int i = 0; i< system.getNumCompounds(); i++)
		{
		cout<<"printing out Compound: "<<i<<endl;	
		system.updCompound(Compound::Index(i)).writePdb(study.getState(),outputFrame,Vec3(0));
		}
    cout<<"done printing out minimized structure"<<endl;	
    cout<<"multithreading set to : "<<dumm.getUseMultithreadedComputation()<<endl;
    cout<<"Gbsa set to : "<<myGbsaGlobalScaleFactor      <<endl;
#ifdef CREATE_VTK_WINDOW
    display.report(study.getState());
#endif

    TimeStepper ts(system,study);
    ts.initialize(system.getDefaultState());
    ts.stepTo(30);

   

    for (Real simTime=0.0; simTime < (10000 * myReportingInterval); simTime += myReportingInterval) // picoseconds
    {
	p++;
        study.stepTo(simTime);

	State currentState = study.getState();
        stringstream ss1;
        ss1<<"out."<<p<<".pdb";
        ofstream outputFrame(ss1.str().c_str());

        cout<<"p = "<<p<<"Last successful step size = "<<study.getPreviousStepSizeTaken()<<endl;

	cout << "current temperature = "<<(*myVelocityRescalingThermostat).getTemperature()<<endl;

        time ( &rawtime );
        timeinfo = localtime ( &rawtime );
        cout   <<"REMARK local time: "<<asctime (timeinfo) <<"REMARK elapsed time: "<<(clock()/CLOCKS_PER_SEC)<<endl;

	outputFrame<<"MODEL "<<p<<endl;
	for (int i = 0; i< system.getNumCompounds(); i++)
		{
		cout<<"printing out Compound: "<<i<<endl;	
		system.updCompound(Compound::Index(i)).writePdb(currentState,outputFrame,Vec3(0));
		}
		//system.updCompound(Compound::Index(i)).writePdb(currentState,outputFrame,Vec3(0));
	//methanol1.writePdb(currentState,outputFrame,Vec3(0));

        //for (int i=0; i<numWater; i++)  {
	// 	(*myWaterVec[i]).writePdb(  currentState,outputFrame,Vec3(0));
	//}
	
    //for (int i = 1; i<numDivalents; i++)
    {
	//trade study.getstate for staticlydefiend state
	//myMagnesiumIonVec[i].writePdb(  currentState  ,outputFrame,Vec3(0));	
    }	 	
        //myMolecule.writePdb(study.getState(),outputFrame,Vec3(0));
	outputFrame<<"ENDMDL"<<endl;
#ifdef CREATE_VTK_WINDOW
        display.report(study.getState());
#endif

    }
    for (int i=0; i<numWater; i++)  {
	delete myWaterVec[i];
    }

}




/* Portions copyright (c) 2006 Stanford University and Simbios.
 * Contributors: Pande Group
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject
 * to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE
 * LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 * OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 * WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

#include "cpuObcInterface.h"

#include "SimTKOpenMMLog.h"
#include "SimTKOpenMMUtilities.h"
#include "ObcParameters.h"
#include "CpuObc.h"

using namespace SimTK;

/**---------------------------------------------------------------------------------------

	Setup for Obc calculations like Gromacs

   @param numberOfAtoms            number of atoms

   @param obcScaleFactors          array of OBC scale factors (one entry each atom)

   @param atomicRadii              atomic radii in Angstrom (one entry each atom)

   @param includeAceApproximation  if true, then include nonpolar 
                                   ACE term in calculations

   @param soluteDielectric         solute dielectric

   @param solventDielectric        solvent dielectric

   @param log                      log reference -- if NULL, then errors/warnings
                                   output to stderr

   The method creates a CpuObc instance -- currently the OBC type II model is the
   default (see paper). If the OBC type I model is desired change

      ObcParameters* obcParameters  = new ObcParameters( numberOfAtoms, ObcParameters::ObcTypeII );
   to
      ObcParameters* obcParameters  = new ObcParameters( numberOfAtoms, ObcParameters::ObcTypeI  );

   The created object is a static member of the class CpuObc; 
   when the force routine, cpuCalculateObcForces(), is called, 
   the static object is used to compute the forces and energy 

   @return 0

   --------------------------------------------------------------------------------------- */

extern "C" int 
cpuSetObcParameters( int numberOfAtoms, RealOpenMM* atomicRadii, RealOpenMM* obcScaleFactors,
                     int includeAceApproximation,
                     RealOpenMM soluteDielectric, RealOpenMM solventDielectric, FILE* log ){

   // ---------------------------------------------------------------------------------------

   static const char* methodName   = "\ncpuSetObcParameters: ";

   // ---------------------------------------------------------------------------------------
   
   // set log file if not NULL

   if( log ){
      SimTKOpenMMLog::setSimTKOpenMMLog( log );
   }

   // set OBC parameters (Type II)

   ObcParameters* obcParameters  = new ObcParameters( numberOfAtoms, ObcParameters::ObcTypeII );
   obcParameters->setScaledRadiusFactors( obcScaleFactors );
   obcParameters->setAtomicRadii( atomicRadii, SimTKOpenMMCommon::KcalAngUnits );

   // dielectric constants

   obcParameters->setSolventDielectric( solventDielectric );
   obcParameters->setSoluteDielectric( soluteDielectric );

   // ---------------------------------------------------------------------------------------

   // create CpuObc instance that will calculate forces
  
   CpuObc* cpuObc = new CpuObc( obcParameters );

   // set static member for subsequent calls to calculate forces/energy 

   CpuImplicitSolvent::setCpuImplicitSolvent( cpuObc );

   // set base file name, ...

   //cpuObc->readInfoFile( "CpuImplicitSolventInfo" );

   // include/do not include ACE approximation (nonpolar solvation)

   cpuObc->setIncludeAceApproximation( includeAceApproximation );

   // ---------------------------------------------------------------------------------------

   // diagnostics
 
   if( log ){
      std::string state = cpuObc->getStateString( methodName );
      (void) fprintf( log, "\n%s\nDone w/ setup\n", state.c_str() );
      (void) fflush( log );
   }

   // ---------------------------------------------------------------------------------------

   return 0;
}

/**---------------------------------------------------------------------------------------

   Calculate implicit solvent forces and energy

   @param atomCoordinates   atom coordinates in Angstrom; format of array is
                            atomCoordinates[atom][3] in Angstrom

   @param partialCharges    partial charges

   @param forces            output forces in kcal/mol.A; format of array is 
                            forces[atom][3]

   @param energy            output energy in kcal/mol

   @param updateBornRadii   if set, then Born radii are updated for current configuration; 
                            otherwise radii correspond to configuration from previous iteration

   Function calls a static method in CpuImplicitSolvent class to calculate forces/energy

   @return result from CpuImplicitSolvent::computeImplicitSolventForces

   --------------------------------------------------------------------------------------- */

extern "C" int
cpuCalculateImplicitSolventForces( const RealOpenMM** atomCoordinates,
                                   const RealOpenMM* partialCharges,
                                   RealOpenMM** forces, RealOpenMM* energy, SimTK::Parallel2DExecutor* executor ){

   // ---------------------------------------------------------------------------------------

   //static const char* methodName   = "\ncpuCalculateImplicitSolventForces: ";

   // ---------------------------------------------------------------------------------------

   int status = CpuImplicitSolvent::getCpuImplicitSolvent()->computeImplicitSolventForces( atomCoordinates, partialCharges,
                                                                  forces, executor );

   *energy = CpuImplicitSolvent::getCpuImplicitSolvent()->getEnergy(); 
   // printf( "\ncpuCalculateImplicitSolventForcesE=%.5e", *energy );

   return status;

}

/**---------------------------------------------------------------------------------------

   Retrieve the calculated energy from the static class member
   The energy is calculated in cpuCalculateImplicitSolventForces()
   along w/ the forces

   @return the calculated energy from the static class member

   --------------------------------------------------------------------------------------- */

extern "C" RealOpenMM cpuGetImplicitSolventEnergy( void ){

   // ---------------------------------------------------------------------------------------

   //static const char* methodName   = "\ncpuGetImplicitSolventEnergy: ";

   // ---------------------------------------------------------------------------------------

   RealOpenMM energy =  CpuImplicitSolvent::getCpuImplicitSolvent()->getEnergy();
   // printf( "\ncpuGetImplicitSolventEnergy E=%.5e", energy );

   return energy;
}

/**---------------------------------------------------------------------------------------

   Delete the Obc associated object(s)

   @return 0 if static CpuObc object was set; else return -1

   --------------------------------------------------------------------------------------- */

extern "C" int cpuDeleteObcParameters( void ){
   return CpuImplicitSolvent::deleteCpuImplicitSolvent();
}

/**---------------------------------------------------------------------------------------

   Get OBC scale factors given masses

   @param numberOfAtoms number of atoms
   @param masses        input masses 
   @param scaleFactors  output atomic numbers

   @return SimTKOpenMMCommon::DefaultReturn

   --------------------------------------------------------------------------------------- */

extern "C" int getObcScaleFactorsGivenAtomMasses( int numberOfAtoms, const RealOpenMM* masses,
                                                  RealOpenMM* scaleFactors ){

   // ---------------------------------------------------------------------------------------

   static const std::string methodName = "\ngetObcScaleFactorsGivenAtomMasses";

   // ---------------------------------------------------------------------------------------

   for( int atomI = 0; atomI < numberOfAtoms; atomI++ ){

      double scaleFactor = 0.8;
      RealOpenMM mass    = masses[atomI];

      if ( mass < 1.2 && mass >= 1.0 ){        // hydrogen
         scaleFactor  = 0.85; 
      } else if( mass > 11.8 && mass < 12.2 ){ // carbon
         scaleFactor  = 0.72; 
      } else if( mass > 14.0 && mass < 15.0 ){ // nitrogen
         scaleFactor  = 0.79;
      } else if( mass > 15.5 && mass < 16.5 ){ // oxygen
         scaleFactor  = 0.85; 
      } else if( mass > 23.9 && mass < 24.5 ){ // magnesium
         scaleFactor  = 0.85; 
      } else if( mass > 31.5 && mass < 32.5 ){ // sulphur
         scaleFactor  = 0.96;
      } else if( mass > 29.5 && mass < 30.5 ){ // phosphorus
         scaleFactor  = 0.86;
      } else {
         std::stringstream message;
         message << methodName;
         message << " Warning: mass for atom " << atomI << " mass=" << mass << "> not recognized.";
         SimTKOpenMMLog::printMessage( message );
      }

      scaleFactors[atomI] = (RealOpenMM) scaleFactor;
   }

   return SimTKOpenMMCommon::DefaultReturn;
}

/**---------------------------------------------------------------------------------------

   Get OBC scale factors given atomic numbers

   @param numberOfAtoms number of atoms
   @param atomicNumber  input atomic number for each atom
   @param scaleFactors  output atomic numbers

   @return SimTKOpenMMCommon::DefaultReturn

   --------------------------------------------------------------------------------------- */

extern "C" int getObcScaleFactors( int numberOfAtoms, const int* atomicNumber, RealOpenMM* scaleFactors ){

   // ---------------------------------------------------------------------------------------

   static const std::string methodName = "getObcScaleFactors";

   // ---------------------------------------------------------------------------------------

   for( int atomI = 0; atomI < numberOfAtoms; atomI++ ){

      double scaleFactor;
      switch( atomicNumber[atomI] ){

         case 1: // hydrogen

            scaleFactor  = 0.85; 
            break;

         case 6: // carbon

            scaleFactor  = 0.72; 
            break;

         case 7: // nitrogen

            scaleFactor  = 0.79;
            break;

         case 8: // oxygen

            scaleFactor  = 0.85;
            break;

         case 9: // Fluorine
            
            scaleFactor = 0.88;
            break;

         case 12: // magnesium

            scaleFactor  = 0.85;
            break;

         case 15: // phosphorus

            scaleFactor  = 0.86;
            break;

         case 16: // sulphur

            scaleFactor  = 0.96;
            break;

         default:

            scaleFactor = 0.8;

            std::stringstream message;
            message << methodName << "(): ";
            message << "Warning: element #" << atomicNumber[atomI] << " for atom " << atomI 
                    << " not handled -- using default scale factor.\n";
            SimTKOpenMMLog::printMessage( message );
            break;
      }

      scaleFactors[atomI] = (RealOpenMM) scaleFactor;
   }

   return SimTKOpenMMCommon::DefaultReturn;
}

/**---------------------------------------------------------------------------------------

   Get GBSA radii

   @param numberOfAtoms             number of atoms

   @param atomicNumber              input atomic number for each atom

   @param numberOfCovalentPartners  input number of covalent partners for each atom
                                    1 for H, 2,3,4 for C, ...; 
                                    the values are only used for C, N & O

   @param atomicNumberOfHCovalentPartner    used only for H. E.g. if atom 22 is a H and 
                                    it is bonded to a nitrogen,
                                    then atomicNumberOfCovalentPartner[22] = 7

   @param gbsaRadii                 output GBSA radii

   @return SimTKOpenMMCommon::DefaultReturn

   --------------------------------------------------------------------------------------- */

extern "C" int getGbsaRadii( int numberOfAtoms, const int* atomicNumber, 
                             const int* numberOfCovalentPartners, 
                             const int* atomicNumberOfHCovalentPartner,
                             RealOpenMM* gbsaRadii ){

   // ---------------------------------------------------------------------------------------

   static const std::string methodName = "\ngetGbsaRadii";

   // ---------------------------------------------------------------------------------------

   // loop over atoms

   // Sulea T.A. has changed the values for most of these so they are similar to those used
   // in the Amber forcefield, which was used as input to Robosample at time of writing (Q1-2023) 
   // These are equivalent to the mbondi2 Generalized Born Parameter Set (iGBparm=6)

   for( int atomI = 0; atomI < numberOfAtoms; atomI++ ){

      // branch based on atomic number

      double radius;
      switch ( atomicNumber[atomI] ){

         case 1: // H

            // radius is modified if heavy atom is N or O
            assert(atomicNumberOfHCovalentPartner[atomI] > 0);
            if( atomicNumberOfHCovalentPartner[atomI] == 7 ){
               radius = 1.3;
            } else if( atomicNumberOfHCovalentPartner[atomI] == 8 ){
               radius = 1.2;
            }else if( atomicNumberOfHCovalentPartner[atomI] == 16 ){
               radius = 1.2;
            }
            else {
               radius = 1.2;
            }
            break;

         case 6: // C

            radius = 1.7;
            break;

         case 7: // N

            radius = 1.55;
            break;

         case 8: // O

            radius = 1.5;
            break;

         case 9: // F
            radius = 1.5;
            break;

         case 14: // Si
            radius = 2.1;
            break;

         case 15: // P
            radius = 1.85;
            break;
         case 16: // S
            radius = 1.8;
            break;
         case 17: // Cl
            radius = 1.7;
            break;

         default:
            radius = 1.5;
            break;
      }
         
      gbsaRadii[atomI] = (RealOpenMM) radius;
   }

   return SimTKOpenMMCommon::DefaultReturn;
}

#include <mpi.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include "dsmc_output.h"
#include "dsmc_parameter.h"

using namespace std ;


//==============================================================================================================
//==============================================================================================================

void OutputResult(	DSMC_DOMAIN		*h_pDomain ,
			DSMC_CELL		*h_Cell ,
			DSMC_WALLTYPE		*h_WallType ,
			DSMC_SURFACE		*h_Surface ,
			DSMC_DSMC		*h_pDSMC ,
			DSMC_RESULT		*h_pResult ,
			DSMC_SPECIES		*h_Species ,
			int 			TimestepNo ){
	
	MPI_Barrier( MPI_COMM_WORLD ) ;
	
	// Output simulation result for flow field.				
	OutputResultFlowfield(	h_pDomain , h_Cell , h_pResult , TimestepNo ) ;
	
	
	// Output simulation result of each species for flow field.
	if ( h_pDomain->SpeciesNum > 1 )
		OutputResultFlowfieldSpecies( h_pDomain , h_Cell , h_pResult , TimestepNo ) ;
	
	
	// Output siumation result for surface properties.				
	OutputResultSurface( h_pDomain , h_Cell , h_WallType , h_Surface , h_pDSMC , TimestepNo ) ;
	
	
	OutputCellInformation( h_pDomain , h_Cell , h_Species ) ;
}

//==============================================================================================================

void OutputResultFlowfield(	DSMC_DOMAIN		*h_pDomain ,
				DSMC_CELL		*h_Cell ,
				DSMC_RESULT		*h_pResult ,
				int 			TimestepNo ){
	
	int		MPISize , MPIMyID , CellNum ;
	ofstream	Output ;
	string		FileName , buffer1 ;
	char		buffer2[8] ;

	MPI_Comm_size( MPI_COMM_WORLD , &MPISize ) ;
	MPI_Comm_rank( MPI_COMM_WORLD , &MPIMyID ) ;
	

	CellNum		= h_pDomain->CellNum ;

	sprintf( buffer2 , "%d" , TimestepNo ) ;
	
	buffer1		= buffer2 ;
	FileName	= "Result-"+buffer1+".dat" ;


	for ( int i=0 ; i<MPISize ; i++ ){
		if ( MPIMyID == i ){
			//Output.open( FileName.c_str() , ios::out | ios::app ) ;
			if ( MPIMyID == 0 )
				Output.open( FileName.c_str() , ios::out | ios::trunc ) ;
			else
				Output.open( FileName.c_str() , ios::out | ios::app ) ;
				
	
			if ( !Output.fail() ){
				if ( MPIMyID == 0 ){
					Output << setw(20) << "CellNo" << setw(20) << "XCenter" << setw(20) << "YCenter" << setw(20) << "ZCenter" << setw(20) << "Density"
					<< setw(20) << "NumDensity" << setw(20) << "U-Vel" << setw(20) << "V-Vel" << setw(20) << "W-Vel" << setw(20) << "Temp" 
					<< setw(20) << "TransTemp" << setw(20) << "RotTemp" << setw(20) << "VibTemp" << setw(20) << "AveParticlesPerCell" 
					<< setw(20) << "mcs/mfp" << setw(20) << "ProcessorNo" << '\n' ;
				}

				for ( int j=0 ; j<CellNum ; j++ ){
					Output << setw(20) << h_Cell[j].Id << setw(20) << h_Cell[j].XCenter << setw(20) << h_Cell[j].YCenter << setw(20) << h_Cell[j].ZCenter << setw(20) << h_pResult->Density[j]
						<< setw(20) << h_pResult->NumDensity[j] << setw(20) << h_pResult->XVel[j] << setw(20) << h_pResult->YVel[j] << setw(20) << h_pResult->ZVel[j] 
						<< setw(20) << h_pResult->Temp[j] << setw(20) << h_pResult->TransTemp[j] << setw(20) << h_pResult->RotTemp[j] << setw(20) << h_pResult->VibTemp[j] 
						<< setw(20) << h_pResult->AveParticleNum[j] << setw(20) << h_pResult->MeanCollSpacingMeanFreePath[j] << setw(20) << i << '\n' ;
				}
			}else{
				cout << "Fail to open " << FileName << '\n' ;
			}

			Output.clear() ;
			Output.close() ;
		}
		
		MPI_Barrier( MPI_COMM_WORLD ) ;
	}
}

//==============================================================================================================

void OutputResultFlowfieldSpecies(	DSMC_DOMAIN		*h_pDomain ,
					DSMC_CELL		*h_Cell ,
					DSMC_RESULT		*h_pResult ,
					int 			TimestepNo ){
	
	int		MPISize , MPIMyID , CellNum , SpeciesNum ;
	ofstream	Output ;
	string		FileName , buffer1 , buffer2 ;
	char		buffer3[8] , buffer4[8] ;


	MPI_Comm_size( MPI_COMM_WORLD , &MPISize ) ;
	MPI_Comm_rank( MPI_COMM_WORLD , &MPIMyID ) ;


	CellNum		= h_pDomain->CellNum ;
	SpeciesNum	= h_pDomain->SpeciesNum ;


	for ( int i=0 ; i<SpeciesNum ; i++ ){
		sprintf( buffer3 , "%d" , i ) ;
		sprintf( buffer4 , "%d" , TimestepNo ) ;
	
		buffer1		= buffer3 ;
		buffer2		= buffer4 ;
		FileName	= "Result-Species"+buffer1+'-'+buffer2+".dat" ;
		
		for ( int p=0 ; p<MPISize ; p++ ){
			if ( MPIMyID == p ){
				//Output.open( FileName.c_str() , ios::out | ios::app ) ;
				if ( MPIMyID == 0 )
					Output.open( FileName.c_str() , ios::out | ios::trunc ) ;
				else
					Output.open( FileName.c_str() , ios::out | ios::app ) ;
		
		
				if ( !Output.fail() ){
					if ( MPIMyID == 0 ){
						Output << setw(20) << "CellNo" << setw(20) << "XCenter" << setw(20) << "YCenter" << setw(20) << "ZCenter" << setw(20) << "Density"
							<< setw(20) << "NumDensity" << setw(20) << "U-Vel" << setw(20) << "V-Vel" << setw(20) << "W-Vel" << setw(20) << "Temp" 
							<< setw(20) << "TransTemp" << setw(20) << "RotTemp" << setw(20) << "VibTemp" << setw(20) << "TransXTemp" 
							<< setw(20) << "TransYTemp" << setw(20) << "TransZTemp" << setw(20) << "AveParticlesPerCell" << '\n' ;
					}


					for ( int j=0 ; j<CellNum ; j++ ){
						Output << setw(20) << h_Cell[j].Id << setw(20) << h_Cell[j].XCenter << setw(20) << h_Cell[j].YCenter << setw(20) << h_Cell[j].ZCenter << setw(20) << h_pResult->DensitySpecies[i][j]
							<< setw(20) << h_pResult->NumDensitySpecies[i][j] << setw(20) << h_pResult->XVelSpecies[i][j] << setw(20) << h_pResult->YVelSpecies[i][j] << setw(20) << h_pResult->ZVelSpecies[i][j] 
							<< setw(20) << h_pResult->TempSpecies[i][j] << setw(20) << h_pResult->TransTempSpecies[i][j] << setw(20) << h_pResult->RotTempSpecies[i][j] << setw(20) << h_pResult->VibTempSpecies[i][j] 
							<< setw(20) << h_pResult->TransXTempSpecies[i][j] << setw(20) << h_pResult->TransYTempSpecies[i][j] << setw(20) << h_pResult->TransZTempSpecies[i][j] << setw(20) << h_pResult->AveParticleNumSpecies[i][j]	<< '\n' ;
					}
		
				}else{
					cout << "Fail to open " << FileName << '\n' ;
				}
		
				Output.clear() ;
				Output.close() ;
			}
			
			MPI_Barrier( MPI_COMM_WORLD ) ;
		}
	}
}

//==============================================================================================================

void OutputResultSurface(	DSMC_DOMAIN		*h_pDomain ,
				DSMC_CELL		*h_Cell ,
				DSMC_WALLTYPE		*h_WallType ,
				DSMC_SURFACE		*h_Surface ,
				DSMC_DSMC		*h_pDSMC ,
				int			TimestepNo ){
					
	int		MPISize , MPIMyID ;
	ofstream	Output ;
	string		FileName , buffer1 ;
	char		buffer2[8] ;
	int		CellNo , WallNo , FaceNo , WallFaceNum , SpeciesNum ;
	double		SurPro[17] , SamplingTime , Weight , Area , XCenter , YCenter , ZCenter ;
	
	
	MPI_Comm_size( MPI_COMM_WORLD , &MPISize ) ;
	MPI_Comm_rank( MPI_COMM_WORLD , &MPIMyID ) ;
	
	
	WallFaceNum	= h_pDomain->WallFaceNum ;
	SpeciesNum	= h_pDomain->SpeciesNum ;
	
	
	sprintf( buffer2 , "%d" , TimestepNo ) ;
	
	buffer1		= buffer2 ;
	FileName	= "SurfaceSample-"+buffer1+".dat" ;


	for ( int z=0 ; z<MPISize ; z++ ){
		if ( MPIMyID == z ){
			//Output.open( FileName.c_str() , ios::out | ios::app ) ;
			if ( MPIMyID == 0 )
				Output.open( FileName.c_str() , ios::out | ios::trunc ) ;
			else
				Output.open( FileName.c_str() , ios::out | ios::app ) ;


			if ( !Output.fail() ){
				if ( MPIMyID == 0 ){
					Output << setw(12) << "1-CellNo" << setw(16) << "2-XCenter" << setw(16) << "3-YCenter" << setw(16) << "4-ZCenter" << setw(20) << "5-SP0" 
						<< setw(20) << "6-SP1" << setw(20) << "7-SP2" << setw(20) << "8-SP3" << setw(20) << "9-SP4" << setw(20) << "10-SP5" 
						<< setw(20) << "11-SP6" << setw(20) << "12-SP7" << setw(20) << "13-SP8" << setw(20) << "14-SP9" << setw(20) << "15-SP10" 
						<< setw(20) << "16-SP11" << setw(20) << "17-SP12" << setw(20) << "18-SP13"<< setw(20) << "19-SP14" << setw(20) << "20-SP15" 
						<< setw(20) << "21-SP16" << setw(16) << "22-Area" << setw(16) << "23-Time" << setw(15) << "24-WallNo" << setw(15) << "25-FaceNo" << '\n' ;
				}
					

				for ( int i=0 ; i<WallFaceNum ; i++ ){
					CellNo		= h_Surface[i].CellNo ;
					WallNo		= h_WallType[h_Surface[i].WallNo].WallNo ;
					FaceNo		= h_Surface[i].FaceNo ;
					XCenter		= h_Surface[i].XCenter ;
					YCenter		= h_Surface[i].YCenter ;
					ZCenter		= h_Surface[i].ZCenter ;
					Area		= h_Surface[i].Area ;
			
					SamplingTime	= h_pDSMC->SamplingTimeEnd[CellNo] - h_pDSMC->SamplingTimeInit[CellNo] ;
					Weight		= h_Cell[CellNo].Weighting / SamplingTime ;
			
			
					for ( int j=0 ; j<17 ; j++ )
						SurPro[j]	= 0. ;
			
			
					for ( int j=0 ; j<SpeciesNum ; j++ ){
						SurPro[0]	+= h_pDSMC->SampleSurfaceParticleNum[i][j] ;
						SurPro[1]	+= h_pDSMC->SampleSurfaceParticleNum[i][j] ;
						SurPro[2]	+= h_pDSMC->SampleSurfaceStickingParticleNum[i][j] ;
						SurPro[3]	+= h_pDSMC->SampleSurfaceInNormMomentum[i][j] ;
						SurPro[4]	+= h_pDSMC->SampleSurfaceReNormMomentum[i][j] ;
						SurPro[5]	+= h_pDSMC->SampleSurfaceInXMomentum[i][j] ;
						SurPro[6]	+= h_pDSMC->SampleSurfaceReXMomentum[i][j] ;
						SurPro[7]	+= h_pDSMC->SampleSurfaceInYMomentum[i][j] ;
						SurPro[8]	+= h_pDSMC->SampleSurfaceReYMomentum[i][j] ;
						SurPro[9]	+= h_pDSMC->SampleSurfaceInZMomentum[i][j] ;
						SurPro[10]	+= h_pDSMC->SampleSurfaceReZMomentum[i][j] ;
						SurPro[11]	+= h_pDSMC->SampleSurfaceInTransEng[i][j] ;
						SurPro[12]	+= h_pDSMC->SampleSurfaceReTransEng[i][j] ;
						SurPro[13]	+= h_pDSMC->SampleSurfaceInRotEng[i][j] ;
						SurPro[14]	+= h_pDSMC->SampleSurfaceReRotEng[i][j] ;
						SurPro[15]	+= h_pDSMC->SampleSurfaceInVibEng[i][j] ;
						SurPro[16]	+= h_pDSMC->SampleSurfaceReVibEng[i][j] ;
					}
			
	
					for ( int j=1 ; j<17 ; j++ )
						SurPro[j]	= SurPro[j]*Weight/Area ;


					CellNo = h_Cell[CellNo].Id ;
				
					Output << setw(12) << CellNo << setw(16) << XCenter << setw(16) << YCenter << setw(16) << ZCenter << setw(20) << SurPro[0] 
						<< setw(20) << SurPro[1] << setw(20) << SurPro[2] << setw(20) << SurPro[3] << setw(20) << SurPro[4] << setw(20) << SurPro[5] 
						<< setw(20) << SurPro[6] << setw(20) << SurPro[7] << setw(20) << SurPro[8] << setw(20) << SurPro[9] << setw(20) << SurPro[10] 
						<< setw(20) << SurPro[11] << setw(20) << SurPro[12] << setw(20) << SurPro[13] << setw(20) << SurPro[14] << setw(20) << SurPro[15] 
						<< setw(20) << SurPro[16] << setw(16) << Area << setw(16) << SamplingTime << setw(15) << WallNo << setw(15) << FaceNo << '\n' ;
				}
			}else{
				cout << "Fail to open " << FileName << '\n' ;
			}
			
			Output.clear() ;
			Output.close() ;
		}
		
		MPI_Barrier( MPI_COMM_WORLD ) ;
	}
}

//==============================================================================================================

void OutputCellInformation( DSMC_DOMAIN *h_pDomain , DSMC_CELL *h_Cell , DSMC_SPECIES *h_Species ){
	
	int		MPISize , MPIMyID , CellNum ;
	double		*NewTimestep , *NewWeighting , Length ;
	ofstream	Output ;


	MPI_Comm_size( MPI_COMM_WORLD , &MPISize ) ;
	MPI_Comm_rank( MPI_COMM_WORLD , &MPIMyID ) ;
	

	CellNum		= h_pDomain->CellNum ;
	
	NewTimestep	= new double[CellNum] ;
	NewWeighting	= new double[CellNum] ;

	
	// Calculate new timestep and weighting based on local mean free path, old timestep, and old weighting.
	for ( int i=0 ; i<CellNum ; i++ ){
		if ( h_Cell[i].MeanFreePath < h_Cell[i].CharacteristicLength  )
			Length	= h_Cell[i].MeanFreePath ;
		else
			Length	= h_Cell[i].CharacteristicLength ;
		
		NewTimestep[i]	= Length/(3.*(h_Cell[i].Speed + sqrt(2.*BOLTZ*h_Cell[i].Temp/h_Species[0].Mass))) ;
		NewWeighting[i]	= h_Cell[i].Weighting * (NewTimestep[i]/h_Cell[i].Timestep) ;
	}


	for ( int i=0 ; i<MPISize ; i++ ){
		if ( MPIMyID == i ){
			if ( MPIMyID == 0 )
				Output.open( "Cell-Information.dat" , ios::out | ios::trunc ) ;
			else
				Output.open( "Cell-Information.dat" , ios::out | ios::app ) ;
					
	
			if ( !Output.fail() ){
				if ( MPIMyID == 0 )
					/*Output << setw(20) << "CellNo" << setw(20) << "MeanFreePath" << setw(20) << "MeanCollisionTime" << setw(20) << "Temp" << setw(20) << "Speed" 
						<< setw(20) << "Timestep" << setw(20) << "Weighting" << setw(20) << "T_Ratio" << setw(20) << "W_Ratio" << setw(20) << "SubcellNum" 
						<< setw(20) << "SelectNum" << setw(20) << "CollRatio" << '\n' ;*/
					Output << setw(20) << "CellNo" << setw(20) << "MeanFreePath" << setw(20) << "MeanCollisionTime" << setw(20) << "Timestep" << setw(20) << "Weighting" 
						<< setw(20) << "NewTimestep" << setw(20) << "NewWeighting" << setw(20) << "AveParticlesPerCell" << setw(20) << "mcs/mfp" << setw(20) << "SubcellNum" 
						<< setw(20) << "SelectNum" << setw(20) << "CollRatio" << '\n' ;


				for ( int j=0 ; j<CellNum ; j++ ){
					/*Output << setw(20) << h_Cell[j].Id << setw(20) << h_Cell[j].MeanFreePath << setw(20) << h_Cell[j].MeanCollisionTime << setw(20) << h_Cell[j].Temp 
						<< setw(20) << h_Cell[j].Speed << setw(20) << h_Cell[j].Timestep << setw(20) << h_Cell[j].Weighting 
						<< setw(20) << (h_Cell[j].Timestep/h_Cell[j].InitTimestep) << setw(20) << (h_Cell[j].Weighting/h_Cell[j].InitWeighting) 
						<< setw(20) << h_Cell[j].SubcellNum << setw(20) << h_Cell[j].SelectNum << setw(20) << h_Cell[j].CollRatio << '\n' ;*/
					Output << setw(20) << h_Cell[j].Id << setw(20) << h_Cell[j].MeanFreePath << setw(20) << h_Cell[j].MeanCollisionTime << setw(20) << h_Cell[j].Timestep 
						<< setw(20) << h_Cell[j].Weighting << setw(20) << NewTimestep[j] << setw(20) << NewWeighting[j] << setw(20) << h_Cell[j].AveParticleNum 
						<< setw(20) << h_Cell[j].MeanCollSpacingMeanFreePath << setw(20) << h_Cell[j].SubcellNum << setw(20) << h_Cell[j].SelectNum 
						<< setw(20) << h_Cell[j].CollRatio << '\n' ;
				}
			}else{
				cout << "Fail to open Cell-Information.dat" << '\n' ;
			}

			Output.clear() ;
			Output.close() ;
		}
		
		MPI_Barrier( MPI_COMM_WORLD ) ;
	}
	
	delete [] NewTimestep ;
	delete [] NewWeighting ;
}

//==============================================================================================================
//==============================================================================================================

void DumpParticlePerCell( DSMC_DOMAIN *h_pDomain , DSMC_DSMC *h_pDSMC , DSMC_CELL *h_Cell , int CellNo ){
	int		GroupNum , ParticleNum , M ;
	
	GroupNum	= h_pDomain->SpeciesGroupNum ;

	for ( int i=0 ; i<GroupNum ; i++ ){
		ParticleNum	= h_pDSMC->IndexCell2[i][CellNo] ;

		cout << "==================================================================\n" ;
		cout << "GroupNo: " << i << ", ParticleNum: " << ParticleNum << ", Timestep: " << h_Cell[CellNo].Timestep << '\n' ;
		cout << "==================================================================\n" ;		
		for ( int j=0 ; j<ParticleNum ; j++ ){
			M	= h_pDSMC->IndexParticle[h_pDSMC->IndexCell1[i][CellNo] + j] ;
			cout	<< setw(10) << M 
				<< setw(8)  << h_pDSMC->ParticleCellNo[M]
				<< setw(4)  << h_pDSMC->ParticleSpeciesNo[M]
				<< setw(16) << h_pDSMC->ParticleXCoord[M]
				<< setw(16) << h_pDSMC->ParticleYCoord[M]
				<< setw(16) << h_pDSMC->ParticleZCoord[M]
				<< setw(14) << h_pDSMC->ParticleXVel[M]
				<< setw(14) << h_pDSMC->ParticleYVel[M] 
				<< setw(14) << h_pDSMC->ParticleZVel[M]
				<< setw(14) << h_pDSMC->ParticleRotation[M]
				<< setw(14) << h_pDSMC->ParticleVibration[M]
				<< setw(8)  << h_pDSMC->ParticleVibLevel[M]
				<< '\n' ;
		}
	}
}

//==============================================================================================================

void DumpMultiSpecies( DSMC_MULTISPECIES *pMultiSpecies ){
	cout	<< setw(14) << pMultiSpecies->CrossSection
		<< setw(14) << pMultiSpecies->RefTemp
		<< setw(14) << pMultiSpecies->VisTempIndex
		<< setw(14) << pMultiSpecies->VSSParameter
		<< setw(14) << pMultiSpecies->ReduceMass
		<< setw(14) << pMultiSpecies->GammaValue
		<< '\n' ;
}

//==============================================================================================================

void DumpSimulationInformation( DSMC_DOMAIN *h_pDomain , DSMC_TIME *pTime ){
	cout << "\n\n\n====================================================================================================================\n\n" ;
	cout << setw(15) << "Timestep No:" << setw(10) << h_pDomain->TimestepNo << setw(20) << "Particle Number:" << setw(15) << h_pDomain->ParticleNum 
	     << setw(32) << "Ratio of Error Tracking:" << setw(15) << (double)h_pDomain->ErrorTrackingNum/h_pDomain->TimestepNo/h_pDomain->ParticleNum 
	     << setw(22) << "Error Tracking New:" << setw(10) << h_pDomain->ErrorTrackingNumNew << '\n' ;
	cout << "\n------------------------------------------------------------------------------------\n" ;
	cout << setw(30) << "The simulation time (second):" << setw(15) << pTime->Total << '\n' ;
	cout << setw(30) << "Initialization:" 		<< setw(15) << pTime->Init	<< '\n' ;
	cout << setw(30) << "Particle Movement:"	<< setw(15) << pTime->Move	<< '\n' ;
	cout << setw(30) << "Index:"			<< setw(15) << pTime->Index	<< '\n' ;
	cout << setw(30) << "Collision:"		<< setw(15) << pTime->Collision	<< '\n' ; 
	cout << setw(30) << "Sampling:"			<< setw(15) << pTime->Sample	<< '\n' ; 
	cout << setw(30) << "Calculate Result:"		<< setw(15) << pTime->CalResult	<< '\n' ;
}
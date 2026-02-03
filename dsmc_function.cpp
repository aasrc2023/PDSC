#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include "dsmc_function.h"
#include "dsmc_toolfunction.h"
#include "dsmc_parameter.h"
#include "dsmc_output.h"

using namespace std ;

//==============================================================================================================
//==============================================================================================================

void ParticleMovement(	DSMC_DOMAIN		*h_pDomain ,
			DSMC_PROCESSOR		*h_pProcessor ,
			DSMC_NODE		*h_Node ,
			DSMC_CELL		*h_Cell ,
			DSMC_INLET		*h_Inlet ,
			DSMC_WALLTYPE		*h_WallType ,
			DSMC_SPECIES		*h_Species ,
			DSMC_SURFACE		*h_Surface ,
			DSMC_DSMC		*h_pDSMC ,
			CELLMAPPING		*h_pMapping ,
			DSMC_MPI_PARTICLE	*h_MPIParticleOut , 
			DSMC_MPI_PARTICLE	*h_MPIParticleIn ,
			DSMC_MPI_DATATYPE	*pMPIDataType , 
			ofstream		&OutputDebug ){
	
	int		MPISize , MPIMyID ;
	int		Dimension , ParticleNum , *MPIParticleNumOut , *MPIParticleNumIn , MPIParticleNum ;
	
	
	MPI_Comm_size( MPI_COMM_WORLD , &MPISize ) ;
	MPI_Comm_rank( MPI_COMM_WORLD , &MPIMyID ) ;
	

	MPIParticleNumOut	= new int[MPISize] ;
	MPIParticleNumIn	= new int[MPISize] ;

	
	Dimension	= h_pDomain->Dimension ;


	MPIParticleNum	= 0 ;
	SetInitValue( MPIParticleNumOut , 0 , MPISize ) ;
	SetInitValue( MPIParticleNumIn , 0 , MPISize ) ;
	
	
	// Move all particles on the domain.
	if ( Dimension == 2 ){
		MoveAllParticle2D( h_pDomain , h_pProcessor , h_Node , h_Cell , h_WallType , h_Species , h_Surface , h_pDSMC , h_pMapping , 
				   h_MPIParticleIn , &MPIParticleNum , OutputDebug ) ;
	}else if ( Dimension == 3 ){
		MoveAllParticle3D( h_pDomain , h_pProcessor , h_Node , h_Cell , h_WallType , h_Species , h_Surface , h_pDSMC , h_pMapping , 
				   h_MPIParticleIn , &MPIParticleNum , OutputDebug ) ;
	}else if ( Dimension == 4 ){
		MoveAllParticleAxisymmetric( h_pDomain , h_pProcessor , h_Node , h_Cell , h_WallType , h_Species , h_Surface , h_pDSMC , h_pMapping , 
				   h_MPIParticleIn , &MPIParticleNum , OutputDebug ) ;
	}
	
	
	// Enter new particles from inlet boundary (Type: -3 and -4) and move them over a timestep*random[0-1].
	if ( h_pDomain->InletFaceNum >0 ){
		ParticleNum	= h_pDomain->ParticleNum ;
		
		if ( Dimension == 2 ){
			EnterNewPaticle2D( h_pDomain , h_Node , h_Inlet , h_Species , h_pDSMC ) ;
			MoveNewParticle2D( h_pDomain , h_pProcessor , h_Node , h_Cell , h_WallType , h_Species , h_Surface , h_pDSMC , 
					   h_pMapping , ParticleNum , h_MPIParticleIn , &MPIParticleNum , OutputDebug ) ;
		}else if ( Dimension == 3 ){
			EnterNewPaticle3D( h_pDomain , h_Node , h_Cell , h_Inlet , h_Species , h_pDSMC , h_pMapping ) ;
			MoveNewParticle3D( h_pDomain , h_pProcessor , h_Node , h_Cell , h_WallType , h_Species , h_Surface , h_pDSMC , 
					   h_pMapping , ParticleNum , h_MPIParticleIn , &MPIParticleNum , OutputDebug ) ;
		}else if ( Dimension == 4 ){
			EnterNewPaticleAxisymmetric( h_pDomain , h_Node , h_Inlet , h_Species , h_pDSMC ) ;
			MoveNewParticleAxisymmetric( h_pDomain , h_pProcessor , h_Node , h_Cell , h_WallType , h_Species , h_Surface , h_pDSMC , 
					   h_pMapping , ParticleNum , h_MPIParticleIn , &MPIParticleNum , OutputDebug ) ;
		}
	}


	// Move the particles which are from other processors.
	if ( Dimension == 2 ){
		MoveOtherParticle2D( h_pDomain , h_pProcessor , h_Node , h_Cell , h_WallType , h_Species , h_Surface , h_pDSMC , h_pMapping , 
				     h_MPIParticleOut , h_MPIParticleIn , MPIParticleNumOut , MPIParticleNumIn , &MPIParticleNum , pMPIDataType , OutputDebug ) ;
	}else if ( Dimension == 3 ){
		MoveOtherParticle3D( h_pDomain , h_pProcessor , h_Node , h_Cell , h_WallType , h_Species , h_Surface , h_pDSMC , h_pMapping , 
				     h_MPIParticleOut , h_MPIParticleIn , MPIParticleNumOut , MPIParticleNumIn , &MPIParticleNum , pMPIDataType , OutputDebug ) ;
	}else if ( Dimension == 4 ){
		MoveOtherParticleAxisymmetric( h_pDomain , h_pProcessor , h_Node , h_Cell , h_WallType , h_Species , h_Surface , h_pDSMC , h_pMapping , 
				     h_MPIParticleOut , h_MPIParticleIn , MPIParticleNumOut , MPIParticleNumIn , &MPIParticleNum , pMPIDataType , OutputDebug ) ;
	}
	
	
	delete [] MPIParticleNumOut ;
	delete [] MPIParticleNumIn ;
}

//==============================================================================================================

void MoveAllParticle2D(	DSMC_DOMAIN		*h_pDomain ,
			DSMC_PROCESSOR		*h_pProcessor ,
			DSMC_NODE		*h_Node ,
			DSMC_CELL		*h_Cell ,
			DSMC_WALLTYPE		*h_WallType ,
			DSMC_SPECIES		*h_Species ,
			DSMC_SURFACE		*h_Surface ,
			DSMC_DSMC		*h_pDSMC ,
			CELLMAPPING		*h_pMapping ,
			DSMC_MPI_PARTICLE	*h_MPIParticleIn ,
			int			*pMPIParticleNum , 
			ofstream		&OutputDebug ){
				
	int		MPISize , MPIMyID ;
	int		ParticleNum , SpeciesNo , FaceNo , ProcessorNo ;
	int		LocalCellNo , BeforeCellNo , NeighborCellNo ;
	int		TrackingNum , ErrorTrackingNum ;
	double		RemainderTime , Time , CollideTime ;
	double		XCoord , YCoord , XVel , YVel , ZVel , BeforeXCoord , BeforeYCoord ;
	bool		Tracking ;
	
	
	MPI_Comm_size( MPI_COMM_WORLD , &MPISize ) ;
	MPI_Comm_rank( MPI_COMM_WORLD , &MPIMyID ) ;
	
	
	ParticleNum		= h_pDomain->ParticleNum ;
	(*pMPIParticleNum)	= 0 ;
				
	Tracking	= true ;
	ErrorTrackingNum= 0 ;
				
	
	// Debug.			
	/*if ( MPIMyID == 9 ){
		for ( int i=0 ; i<h_pDomain->CellNum ; i++ ){
			OutputDebug << "LocalCellNo: " << i << ", GlobalCellNo: " << h_Cell[i].Id << ", LocalCellNo: " << h_pProcessor->LocalCellNo[h_Cell[i].Id] << endl ;
		}
		
		OutputDebug << "2395: " << h_Cell[2395].Id << ", " << h_pProcessor->LocalCellNo[h_Cell[2395].Id] << endl ;
	}*/
				
				
				
	for ( int i=0 ; i<ParticleNum ; i++ ){
		Tracking	= true ;
		BeforeCellNo	= -1 ;
		
		LocalCellNo	= h_pDSMC->ParticleCellNo[i] ;
		SpeciesNo	= h_pDSMC->ParticleSpeciesNo[i] ;
		XCoord		= h_pDSMC->ParticleXCoord[i] ;
		YCoord		= h_pDSMC->ParticleYCoord[i] ;
		XVel		= h_pDSMC->ParticleXVel[i] ;
		YVel		= h_pDSMC->ParticleYVel[i] ;
		ZVel		= h_pDSMC->ParticleZVel[i] ;
		
		
		RemainderTime	= h_Cell[LocalCellNo].Timestep ;
		TrackingNum	= 0 ;

		
		// Start to track particle position after one timestep.
		do{	
			Time		= RemainderTime ;
			BeforeXCoord	= XCoord ;
			BeforeYCoord	= YCoord ;
			
			
			// Update the particle position directly.
			MoveParticle2D( &XCoord , &YCoord , XVel , YVel , Time ) ;
			
			
			if ( InCell2D( h_Node , &h_Cell[LocalCellNo] , XCoord , YCoord , h_pMapping , 1 ) )
				break ;


			CalculateTimeCollideFace2D( &CollideTime , &FaceNo , Time , XCoord , YCoord , BeforeXCoord , BeforeYCoord , BeforeCellNo , &h_Cell[LocalCellNo] , h_pMapping , 43 ) ;
			
		
			if ( FaceNo == -1 ){
				XCoord	= h_Cell[LocalCellNo].XCenter ;
				YCoord	= h_Cell[LocalCellNo].YCenter ;
				ErrorTrackingNum++ ;
				break ;
			}
			
			
			XCoord	= BeforeXCoord ;
			YCoord	= BeforeYCoord ;
			
			// Move particle to the face of cell.
			MoveParticle2D( &XCoord , &YCoord , XVel , YVel , CollideTime ) ;
			
			
			// Calculate remainder time after moving particle onto face of cell.
			RemainderTime	= Time - CollideTime ;
			NeighborCellNo	= h_Cell[LocalCellNo].Neighbor[FaceNo] ;
			
			
			if ( NeighborCellNo < 0 ){
				Tracking	= ParticleCollideSurface2D( &XVel , &YVel , &ZVel , &h_Cell[LocalCellNo] , FaceNo , &i , &ParticleNum , &h_Species[SpeciesNo] , 
									    h_pDomain , h_Node , h_WallType , h_Surface , h_pDSMC , h_pMapping , OutputDebug ) ;
			}else{
				ProcessorNo	= h_pProcessor->CellProcessorNo[NeighborCellNo] ;
				
				
				// Particle still stay the same processor.
				if ( ProcessorNo == MPIMyID ){
					// Translate the global neighbor cell no. to local.
					NeighborCellNo	= h_pProcessor->LocalCellNo[NeighborCellNo] ;
					
					RemainderTime	= RemainderTime * (h_Cell[NeighborCellNo].Timestep/h_Cell[LocalCellNo].Timestep) ;
				
					// Update partice next cell no.
					BeforeCellNo	= h_Cell[LocalCellNo].Id ;
					LocalCellNo	= NeighborCellNo ;
					
				// Particle go to other processor.
				}else{
					BeforeCellNo			= h_Cell[LocalCellNo].Id ;
					h_pDSMC->ParticleXCoord[i]	= XCoord ;
					h_pDSMC->ParticleYCoord[i]	= YCoord ;
					h_pDSMC->ParticleXVel[i]	= XVel ;
					h_pDSMC->ParticleYVel[i]	= YVel ;
					h_pDSMC->ParticleZVel[i]	= ZVel ;
					
					TransferParticleToBuffer( h_pDSMC , &i , &ParticleNum , h_MPIParticleIn , pMPIParticleNum , ProcessorNo , 
								  NeighborCellNo , BeforeCellNo , RemainderTime , h_Cell[LocalCellNo].Timestep ) ;
					
					Tracking	= false ;
				}
			}
			
			
			TrackingNum++ ;
			
			
			// Debug.
			if ( TrackingNum == 6 ){
				XCoord	= h_Cell[LocalCellNo].XCenter ;
				YCoord	= h_Cell[LocalCellNo].YCenter ;
				break ;	
			}
		}while ( Tracking ) ;
		
		
		// Update particle information after tracking.
		if ( Tracking ){
			h_pDSMC->ParticleCellNo[i]	= LocalCellNo ;
			h_pDSMC->ParticleXCoord[i]	= XCoord ;
			h_pDSMC->ParticleYCoord[i]	= YCoord ;
			h_pDSMC->ParticleXVel[i]	= XVel ;
			h_pDSMC->ParticleYVel[i]	= YVel ;
			h_pDSMC->ParticleZVel[i]	= ZVel ;
		}
	}
	
	h_pDomain->ParticleNum		= ParticleNum ;	
	h_pDomain->TrackingNum		+=ParticleNum ;
	h_pDomain->ErrorTrackingNum	+=ErrorTrackingNum ;
	
	// Debug.
	//OutputDebug << "MoveAllParticle: " << ", Num: "<< setw(20) << ParticleNum << ", ERROR: " << setw(20) << ErrorTrackingNum << endl ;
}


void MoveAllParticle3D(	DSMC_DOMAIN		*h_pDomain ,
			DSMC_PROCESSOR		*h_pProcessor ,
			DSMC_NODE		*h_Node ,
			DSMC_CELL		*h_Cell ,
			DSMC_WALLTYPE		*h_WallType ,
			DSMC_SPECIES		*h_Species ,
			DSMC_SURFACE		*h_Surface ,
			DSMC_DSMC		*h_pDSMC ,
			CELLMAPPING		*h_pMapping ,
			DSMC_MPI_PARTICLE	*h_MPIParticleIn ,
			int			*pMPIParticleNum , 
			ofstream		&OutputDebug ){
				
	int		MPISize , MPIMyID ;
	int		ParticleNum , SpeciesNo , FaceNo , ProcessorNo ;
	int		LocalCellNo , BeforeCellNo , NeighborCellNo ;
	int		TrackingNum , ErrorTrackingNum ;
	double		RemainderTime , Time , CollideTime ;
	double		XCoord , YCoord , ZCoord , XVel , YVel , ZVel , BeforeXCoord , BeforeYCoord , BeforeZCoord ;
	bool		Tracking ;
	
	
	MPI_Comm_size( MPI_COMM_WORLD , &MPISize ) ;
	MPI_Comm_rank( MPI_COMM_WORLD , &MPIMyID ) ;
	
	
	ParticleNum		= h_pDomain->ParticleNum ;
	(*pMPIParticleNum)	= 0 ;
				
				
	Tracking	= true ;
	ErrorTrackingNum= 0 ;
				
				
	for ( int i=0 ; i<ParticleNum ; i++ ){
		Tracking	= true ;
		BeforeCellNo	= -1 ;
		
		LocalCellNo	= h_pDSMC->ParticleCellNo[i] ;
		SpeciesNo	= h_pDSMC->ParticleSpeciesNo[i] ;
		XCoord		= h_pDSMC->ParticleXCoord[i] ;
		YCoord		= h_pDSMC->ParticleYCoord[i] ;
		ZCoord		= h_pDSMC->ParticleZCoord[i] ;
		XVel		= h_pDSMC->ParticleXVel[i] ;
		YVel		= h_pDSMC->ParticleYVel[i] ;
		ZVel		= h_pDSMC->ParticleZVel[i] ;
		
		RemainderTime	= h_Cell[LocalCellNo].Timestep ;
		TrackingNum	= 0 ;


		// Start to track particle position after one timestep.
		do{	
			Time		= RemainderTime ;
			BeforeXCoord	= XCoord ;
			BeforeYCoord	= YCoord ;
			BeforeZCoord	= ZCoord ;
			
			
			// Update the particle position directly.
			MoveParticle3D( &XCoord , &YCoord , &ZCoord , XVel , YVel , ZVel , Time ) ;

			
			if ( InCell3D( h_Node , &h_Cell[LocalCellNo] , XCoord , YCoord , ZCoord , h_pMapping , 0 ) )
				break ;

			
			CalculateTimeCollideFace3D( &CollideTime , &FaceNo , Time , XCoord , YCoord , ZCoord , BeforeXCoord , BeforeYCoord , 
							BeforeZCoord , BeforeCellNo , &h_Cell[LocalCellNo] , h_pMapping , 0 ) ;
			
		
			if ( FaceNo == -1 ){
				XCoord	= h_Cell[LocalCellNo].XCenter ;
				YCoord	= h_Cell[LocalCellNo].YCenter ;
				ZCoord	= h_Cell[LocalCellNo].ZCenter ;
				ErrorTrackingNum++ ;
				break ;
			}			
			
			XCoord	= BeforeXCoord ;
			YCoord	= BeforeYCoord ;
			ZCoord	= BeforeZCoord ;
			
			
			// Move particle to the face of cell.
			MoveParticle3D( &XCoord , &YCoord , &ZCoord , XVel , YVel , ZVel , CollideTime ) ;
			
			
			// Calculate remainder time after moving particle onto face of cell.
			RemainderTime	= Time - CollideTime ;


			NeighborCellNo	= h_Cell[LocalCellNo].Neighbor[FaceNo] ;


			if ( NeighborCellNo < 0 ){
				Tracking	= ParticleCollideSurface3D( &XVel , &YVel , &ZVel , &h_Cell[LocalCellNo] , FaceNo , &i , &ParticleNum , &h_Species[SpeciesNo] , 
									    h_pDomain , h_Node , h_WallType , h_Surface , h_pDSMC , h_pMapping ) ;
									    
			}else{
				ProcessorNo	= h_pProcessor->CellProcessorNo[NeighborCellNo] ;
				
				
				// Particle still stay the same processor.
				if ( ProcessorNo == MPIMyID ){
					// Translate the global neighbor cell no. to local.
					NeighborCellNo	= h_pProcessor->LocalCellNo[NeighborCellNo] ;
				
					RemainderTime	= RemainderTime * (h_Cell[NeighborCellNo].Timestep/h_Cell[LocalCellNo].Timestep) ;
				
					// Update cell no.
					BeforeCellNo	= h_Cell[LocalCellNo].Id ; ;
					LocalCellNo	= NeighborCellNo ;
				
				// Particle go to other processor.
				}else{
					BeforeCellNo			= h_Cell[LocalCellNo].Id ;
					h_pDSMC->ParticleXCoord[i]	= XCoord ;
					h_pDSMC->ParticleYCoord[i]	= YCoord ;
					h_pDSMC->ParticleZCoord[i]	= ZCoord ;
					h_pDSMC->ParticleXVel[i]	= XVel ;
					h_pDSMC->ParticleYVel[i]	= YVel ;
					h_pDSMC->ParticleZVel[i]	= ZVel ;
					
					TransferParticleToBuffer( h_pDSMC , &i , &ParticleNum , h_MPIParticleIn , pMPIParticleNum , ProcessorNo , 
								  NeighborCellNo , BeforeCellNo , RemainderTime , h_Cell[LocalCellNo].Timestep ) ;
					
					Tracking	= false ;
				}
			}
			
			TrackingNum++ ;
			
			// Debug.
			if ( TrackingNum == 6 ){
				XCoord	= h_Cell[LocalCellNo].XCenter ;
				YCoord	= h_Cell[LocalCellNo].YCenter ;
				ZCoord	= h_Cell[LocalCellNo].ZCenter ;
				break ;	
			}
		}while ( Tracking ) ;
		
		// Update particle information after tracking.
		if ( Tracking ){
			h_pDSMC->ParticleCellNo[i]	= LocalCellNo ;
			h_pDSMC->ParticleXCoord[i]	= XCoord ;
			h_pDSMC->ParticleYCoord[i]	= YCoord ;
			h_pDSMC->ParticleZCoord[i]	= ZCoord ;
			h_pDSMC->ParticleXVel[i]	= XVel ;
			h_pDSMC->ParticleYVel[i]	= YVel ;
			h_pDSMC->ParticleZVel[i]	= ZVel ;
		}
	}
	
	h_pDomain->ParticleNum		= ParticleNum ;	
	h_pDomain->TrackingNum		+=ParticleNum ;
	h_pDomain->ErrorTrackingNum	+=ErrorTrackingNum ;
	
	// Debug.
//	OutputDebug << "1. MoveAllParticle: " << ", Num: "<< setw(20) << ParticleNum << ", ERROR: " << setw(20) << ErrorTrackingNum << endl ;
}


void MoveAllParticleAxisymmetric(	DSMC_DOMAIN		*h_pDomain ,
					DSMC_PROCESSOR		*h_pProcessor ,
					DSMC_NODE		*h_Node ,
					DSMC_CELL		*h_Cell ,
					DSMC_WALLTYPE		*h_WallType ,
					DSMC_SPECIES		*h_Species ,
					DSMC_SURFACE		*h_Surface ,
					DSMC_DSMC		*h_pDSMC ,
					CELLMAPPING		*h_pMapping ,
					DSMC_MPI_PARTICLE	*h_MPIParticleIn ,
					int			*pMPIParticleNum , 
					ofstream		&OutputDebug ){
				
	int		MPISize , MPIMyID ;
	int		ParticleNum , SpeciesNo , FaceNo , ProcessorNo , WallNo ;
	int		LocalCellNo , BeforeCellNo , NeighborCellNo ;
	int		TrackingNum , ErrorTrackingNum ;
	double		RemainderTime , Time , CollideTime ;
	double		XCoord , YCoord , yCoord , XVel , YVel , ZVel , BufferYVel ;
	double		BeforeXCoord , BeforeYCoord , BeforeXVel , BeforeYVel , BeforeZVel ;
	bool		Tracking ;
	
	
	MPI_Comm_size( MPI_COMM_WORLD , &MPISize ) ;
	MPI_Comm_rank( MPI_COMM_WORLD , &MPIMyID ) ;
	
	
	ParticleNum		= h_pDomain->ParticleNum ;
	(*pMPIParticleNum)	= 0 ;
				
	Tracking	= true ;
	ErrorTrackingNum= 0 ;		
				
				
	for ( int i=0 ; i<ParticleNum ; i++ ){
		Tracking	= true ;
		BeforeCellNo	= -1 ;
		
		LocalCellNo	= h_pDSMC->ParticleCellNo[i] ;
		SpeciesNo	= h_pDSMC->ParticleSpeciesNo[i] ;
		XCoord		= h_pDSMC->ParticleXCoord[i] ;
		YCoord		= h_pDSMC->ParticleYCoord[i] ;
		XVel		= h_pDSMC->ParticleXVel[i] ;
		YVel		= h_pDSMC->ParticleYVel[i] ;
		ZVel		= h_pDSMC->ParticleZVel[i] ;
		
		
		RemainderTime	= h_Cell[LocalCellNo].Timestep ;
		TrackingNum	= 0 ;

		
		// Start to track particle position after one timestep.
		do{	
			Time		= RemainderTime ;
			BeforeXCoord	= XCoord ;
			BeforeYCoord	= YCoord ;
			BeforeXVel	= XVel ;
			BeforeYVel	= YVel ;
			BeforeZVel	= ZVel ;
			
			
			// Update the particle position directly.
			yCoord = MoveParticleAxisymmetric( &XCoord , &YCoord , &XVel , &YVel , &ZVel , Time ) ;
			
			//OutputDebug << "yCoord: " << yCoord << ", YCoord: " << YCoord << endl ;
			
			
			if ( InCell2D( h_Node , &h_Cell[LocalCellNo] , XCoord , YCoord , h_pMapping , 1 ) )
				break ;


			CalculateTimeCollideFace2D( &CollideTime , &FaceNo , Time , XCoord , YCoord , BeforeXCoord , BeforeYCoord , BeforeCellNo , &h_Cell[LocalCellNo] , h_pMapping , 43 ) ;
			//CalculateTimeCollideFace2D( &CollideTime , &FaceNo , Time , XCoord , yCoord , BeforeXCoord , BeforeYCoord , BeforeCellNo , &h_Cell[LocalCellNo] , h_pMapping , 43 ) ;
			
		
			if ( FaceNo == -1 ){
				XCoord	= h_Cell[LocalCellNo].XCenter ;
				YCoord	= h_Cell[LocalCellNo].YCenter ;
				ErrorTrackingNum++ ;
				break ;
			}


			// Calculate remainder time after moving particle onto face of cell.
			RemainderTime	= Time - CollideTime ;
			NeighborCellNo	= h_Cell[LocalCellNo].Neighbor[FaceNo] ;

			/*if ( yCoord < 0. && NeighborCellNo >= 0 ){
				//OutputDebug << "yCoord: " << yCoord << endl ;
				
				ProcessorNo	= h_pProcessor->CellProcessorNo[NeighborCellNo] ;
				RemainderTime	= 0. ;
				
				// Particle still stay the same processor.
				if ( ProcessorNo == MPIMyID ){
					// Translate the global neighbor cell no. to local.
					NeighborCellNo	= h_pProcessor->LocalCellNo[NeighborCellNo] ;
					
					// Update partice next cell no.
					BeforeCellNo	= h_Cell[LocalCellNo].Id ;
					LocalCellNo	= NeighborCellNo ;
					
				// Particle go to other processor.
				}else{
					BeforeCellNo			= h_Cell[LocalCellNo].Id ;
					h_pDSMC->ParticleXCoord[i]	= XCoord ;
					h_pDSMC->ParticleYCoord[i]	= YCoord ;
					h_pDSMC->ParticleXVel[i]	= XVel ;
					h_pDSMC->ParticleYVel[i]	= YVel ;
					h_pDSMC->ParticleZVel[i]	= ZVel ;
					
					TransferParticleToBuffer( h_pDSMC , &i , &ParticleNum , h_MPIParticleIn , pMPIParticleNum , ProcessorNo , 
								  NeighborCellNo , BeforeCellNo , RemainderTime , h_Cell[LocalCellNo].Timestep ) ;
					
					Tracking	= false ;
				}
				break ;	
			}*/
			
			
			BufferYVel	= (YCoord - BeforeYCoord)/Time ;
			XCoord		= BeforeXCoord ;
			YCoord		= BeforeYCoord ;
			XVel		= BeforeXVel ;
			YVel		= BeforeYVel ;
			ZVel		= BeforeZVel ;
			
			
			// Move particle to the face of cell.
			MoveParticle2D( &XCoord , &YCoord , XVel , BufferYVel , CollideTime ) ;
			
			
			if ( NeighborCellNo < 0 ){
				// Debug.
				//if ( BeforeYVel*YVel < 0. ) OutputDebug << "Before YVel: " << BeforeYVel << ", YVel: " << YVel << endl ;
	
				AdjustVelocityAxisymmetric( BeforeYCoord , &YVel , &ZVel , CollideTime ) ;
				
				Tracking	= ParticleCollideSurface2D( &XVel , &YVel , &ZVel , &h_Cell[LocalCellNo] , FaceNo , &i , &ParticleNum , &h_Species[SpeciesNo] , 
									    h_pDomain , h_Node , h_WallType , h_Surface , h_pDSMC , h_pMapping , OutputDebug ) ;
			}else{
				AdjustVelocityAxisymmetric( BeforeYCoord , &YVel , &ZVel , CollideTime ) ;
				
				ProcessorNo	= h_pProcessor->CellProcessorNo[NeighborCellNo] ;
				
				
				// Particle still stay the same processor.
				if ( ProcessorNo == MPIMyID ){
					// Translate the global neighbor cell no. to local.
					NeighborCellNo	= h_pProcessor->LocalCellNo[NeighborCellNo] ;
					
					RemainderTime	= RemainderTime * (h_Cell[NeighborCellNo].Timestep/h_Cell[LocalCellNo].Timestep) ;
				
					// Update partice next cell no.
					BeforeCellNo	= h_Cell[LocalCellNo].Id ;
					LocalCellNo	= NeighborCellNo ;
					
				// Particle go to other processor.
				}else{
					BeforeCellNo			= h_Cell[LocalCellNo].Id ;
					h_pDSMC->ParticleXCoord[i]	= XCoord ;
					h_pDSMC->ParticleYCoord[i]	= YCoord ;
					h_pDSMC->ParticleXVel[i]	= XVel ;
					h_pDSMC->ParticleYVel[i]	= YVel ;
					h_pDSMC->ParticleZVel[i]	= ZVel ;
					
					TransferParticleToBuffer( h_pDSMC , &i , &ParticleNum , h_MPIParticleIn , pMPIParticleNum , ProcessorNo , 
								  NeighborCellNo , BeforeCellNo , RemainderTime , h_Cell[LocalCellNo].Timestep ) ;
					
					Tracking	= false ;
				}
			}
			
			
			TrackingNum++ ;
			
			
			// Debug.
			if ( TrackingNum == 6 ){
				XCoord	= h_Cell[LocalCellNo].XCenter ;
				YCoord	= h_Cell[LocalCellNo].YCenter ;
				ErrorTrackingNum++ ;
				break ;	
			}
		}while ( Tracking ) ;
		
		
		// Update particle information after tracking.
		if ( Tracking ){
			h_pDSMC->ParticleCellNo[i]	= LocalCellNo ;
			h_pDSMC->ParticleXCoord[i]	= XCoord ;
			h_pDSMC->ParticleYCoord[i]	= YCoord ;
			h_pDSMC->ParticleXVel[i]	= XVel ;
			h_pDSMC->ParticleYVel[i]	= YVel ;
			h_pDSMC->ParticleZVel[i]	= ZVel ;
		}
	}
	
	h_pDomain->ParticleNum		= ParticleNum ;	
	h_pDomain->TrackingNum		+=ParticleNum ;
	h_pDomain->ErrorTrackingNum	+=ErrorTrackingNum ;
	
	// Debug.
	//OutputDebug << "1. MoveAllParticle: " << ", Num: "<< setw(20) << ParticleNum << ", ERROR: " << setw(20) << ErrorTrackingNum << endl ;
}

//==============================================================================================================

void MoveNewParticle2D(	DSMC_DOMAIN		*h_pDomain ,
			DSMC_PROCESSOR		*h_pProcessor ,
			DSMC_NODE		*h_Node ,
			DSMC_CELL		*h_Cell ,
			DSMC_WALLTYPE		*h_WallType ,
			DSMC_SPECIES		*h_Species ,
			DSMC_SURFACE		*h_Surface ,
			DSMC_DSMC		*h_pDSMC ,
			CELLMAPPING		*h_pMapping , 
			int 			ParticleNo ,
			DSMC_MPI_PARTICLE	*h_MPIParticleIn ,
			int			*pMPIParticleNum , 
			ofstream		&OutputDebug ){
				
	int		MPISize , MPIMyID ;
	int		ParticleNum , SpeciesNo , FaceNo , ProcessorNo ;
	int		LocalCellNo , BeforeCellNo , NeighborCellNo ;
	int		TrackingNum , ErrorTrackingNum , DebugValue1 ;
	double		RemainderTime , Time , CollideTime ;
	double		XCoord , YCoord , XVel , YVel , ZVel , BeforeXCoord , BeforeYCoord ;
	bool		Tracking ;
	
	
	MPI_Comm_size( MPI_COMM_WORLD , &MPISize ) ;
	MPI_Comm_rank( MPI_COMM_WORLD , &MPIMyID ) ;
	
	ParticleNum	= h_pDomain->ParticleNum ;
				
	Tracking	= true ;
	ErrorTrackingNum= 0 ;
	DebugValue1	= ParticleNum - ParticleNo ;	
	
	
	for ( int i=ParticleNo ; i<ParticleNum ; i++ ){
		Tracking	= true ;
		BeforeCellNo	= -1 ;
		
		LocalCellNo	= h_pDSMC->ParticleCellNo[i] ;
		SpeciesNo	= h_pDSMC->ParticleSpeciesNo[i] ;
		XCoord		= h_pDSMC->ParticleXCoord[i] ;
		YCoord		= h_pDSMC->ParticleYCoord[i] ;
		XVel		= h_pDSMC->ParticleXVel[i] ;
		YVel		= h_pDSMC->ParticleYVel[i] ;
		ZVel		= h_pDSMC->ParticleZVel[i] ;
		
		RemainderTime	= h_Cell[LocalCellNo].Timestep*Randn() ;
		TrackingNum	= 0 ;

		// Start to track particle position after one timestep.
		do{	
			Time		= RemainderTime ;
			BeforeXCoord	= XCoord ;
			BeforeYCoord	= YCoord ;

			
			// Update the particle position directly.
			MoveParticle2D( &XCoord , &YCoord , XVel , YVel , Time ) ;
			
			
			if ( InCell2D( h_Node , &h_Cell[LocalCellNo] , XCoord , YCoord , h_pMapping , 1 ) )
				break ;


			CalculateTimeCollideFace2D( &CollideTime , &FaceNo , Time , XCoord , YCoord , BeforeXCoord , BeforeYCoord , BeforeCellNo , &h_Cell[LocalCellNo] , h_pMapping , 43 ) ;
		
		
			if ( FaceNo == -1 ){
				XCoord	= h_Cell[LocalCellNo].XCenter ;
				YCoord	= h_Cell[LocalCellNo].YCenter ;
				ErrorTrackingNum++ ;
				break ;
			}			
			
			XCoord	= BeforeXCoord ;
			YCoord	= BeforeYCoord ;
			
			// Move particle to the face of cell.
			MoveParticle2D( &XCoord , &YCoord , XVel , YVel , CollideTime ) ;
			
			// Calculate remainder time after moving particle onto face of cell.
			RemainderTime	= Time - CollideTime ;
			
			NeighborCellNo	= h_Cell[LocalCellNo].Neighbor[FaceNo] ;
			
			if ( NeighborCellNo < 0 ){
				Tracking	= ParticleCollideSurface2D( &XVel , &YVel , &ZVel , &h_Cell[LocalCellNo] , FaceNo , &i , &ParticleNum , &h_Species[SpeciesNo] , 
									    h_pDomain , h_Node , h_WallType , h_Surface , h_pDSMC , h_pMapping , OutputDebug ) ;
									    
			}else{
				ProcessorNo	= h_pProcessor->CellProcessorNo[NeighborCellNo] ;
				
				// Particle still stay the same processor.
				if ( ProcessorNo == MPIMyID ){
					// Translate the global neighbor cell no. to local.
					NeighborCellNo	= h_pProcessor->LocalCellNo[NeighborCellNo] ;
					
					RemainderTime	= RemainderTime * (h_Cell[NeighborCellNo].Timestep/h_Cell[LocalCellNo].Timestep) ;
				
					// Update partice next cell no.
					BeforeCellNo	= h_Cell[LocalCellNo].Id ;
					LocalCellNo	= NeighborCellNo ;
					
				// Particle go to other processor.
				}else{
					BeforeCellNo			= h_Cell[LocalCellNo].Id ;
					h_pDSMC->ParticleXCoord[i]	= XCoord ;
					h_pDSMC->ParticleYCoord[i]	= YCoord ;
					h_pDSMC->ParticleXVel[i]	= XVel ;
					h_pDSMC->ParticleYVel[i]	= YVel ;
					h_pDSMC->ParticleZVel[i]	= ZVel ;
					
					TransferParticleToBuffer( h_pDSMC , &i , &ParticleNum , h_MPIParticleIn , pMPIParticleNum , ProcessorNo , 
								  NeighborCellNo , BeforeCellNo , RemainderTime , h_Cell[LocalCellNo].Timestep ) ;
					
					Tracking	= false ;
				}
			}
			
			TrackingNum++ ;
			
			// Debug.
			if ( TrackingNum == 6 ){
				XCoord	= h_Cell[LocalCellNo].XCenter ;
				YCoord	= h_Cell[LocalCellNo].YCenter ;
				break ;	
			}
			
		}while ( Tracking ) ;
		
		// Update particle information after tracking.
		if ( Tracking ){
			h_pDSMC->ParticleCellNo[i]	= LocalCellNo ;
			h_pDSMC->ParticleXCoord[i]	= XCoord ;
			h_pDSMC->ParticleYCoord[i]	= YCoord ;
			h_pDSMC->ParticleXVel[i]	= XVel ;
			h_pDSMC->ParticleYVel[i]	= YVel ;
			h_pDSMC->ParticleZVel[i]	= ZVel ;
		}
	}
	
	h_pDomain->ParticleNum		= ParticleNum ;
	h_pDomain->TrackingNum		+=ParticleNum ;
	h_pDomain->ErrorTrackingNum	+=ErrorTrackingNum ;
	h_pDomain->ErrorTrackingNumNew	+=ErrorTrackingNum ;
	
	// Debug.
	//OutputDebug << "2. MoveNewParticle: " << " Num: " << DebugValue1 << ", ERROR: " << ErrorTrackingNum << endl ;
}


void MoveNewParticle3D(	DSMC_DOMAIN		*h_pDomain ,
			DSMC_PROCESSOR		*h_pProcessor ,
			DSMC_NODE		*h_Node ,
			DSMC_CELL		*h_Cell ,
			DSMC_WALLTYPE		*h_WallType ,
			DSMC_SPECIES		*h_Species ,
			DSMC_SURFACE		*h_Surface ,
			DSMC_DSMC		*h_pDSMC ,
			CELLMAPPING		*h_pMapping , 
			int 			ParticleNo ,
			DSMC_MPI_PARTICLE	*h_MPIParticleIn ,
			int			*pMPIParticleNum , 
			ofstream		&OutputDebug ){
				
	int		MPISize , MPIMyID ;
	int		ParticleNum , SpeciesNo , FaceNo , ProcessorNo ;
	int		LocalCellNo , BeforeCellNo , NeighborCellNo ;
	int		TrackingNum , ErrorTrackingNum , DebugValue1 ;
	double		RemainderTime , Time , CollideTime ;
	double		XCoord , YCoord , ZCoord , XVel , YVel , ZVel , BeforeXCoord , BeforeYCoord , BeforeZCoord ;
	bool		Tracking ;
	
	
	MPI_Comm_size( MPI_COMM_WORLD , &MPISize ) ;
	MPI_Comm_rank( MPI_COMM_WORLD , &MPIMyID ) ;
	
	
	ParticleNum	= h_pDomain->ParticleNum ;
				
	Tracking	= true ;
	ErrorTrackingNum= 0 ;
				
	
	// Debug.
	DebugValue1	= ParticleNum - ParticleNo ;			
	

	for ( int i=ParticleNo ; i<ParticleNum ; i++ ){
		Tracking	= true ;
		BeforeCellNo	= -1 ;
		
		LocalCellNo	= h_pDSMC->ParticleCellNo[i] ;
		SpeciesNo	= h_pDSMC->ParticleSpeciesNo[i] ;
		XCoord		= h_pDSMC->ParticleXCoord[i] ;
		YCoord		= h_pDSMC->ParticleYCoord[i] ;
		ZCoord		= h_pDSMC->ParticleZCoord[i] ;
		XVel		= h_pDSMC->ParticleXVel[i] ;
		YVel		= h_pDSMC->ParticleYVel[i] ;
		ZVel		= h_pDSMC->ParticleZVel[i] ;
		
		RemainderTime	= h_Cell[LocalCellNo].Timestep*Randn() ;
		TrackingNum	= 0 ;


		// Start to track particle position after one timestep.
		do{	
			Time		= RemainderTime ;
			BeforeXCoord	= XCoord ;
			BeforeYCoord	= YCoord ;
			BeforeZCoord	= ZCoord ;

			
			// Update the particle position directly.
			MoveParticle3D( &XCoord , &YCoord , &ZCoord , XVel , YVel , ZVel , Time ) ;

			
			if ( InCell3D( h_Node , &h_Cell[LocalCellNo] , XCoord , YCoord , ZCoord , h_pMapping , 0 ) )
				break ;
			

			CalculateTimeCollideFace3D( &CollideTime , &FaceNo , Time , XCoord , YCoord , ZCoord , BeforeXCoord , BeforeYCoord , 
							BeforeZCoord , BeforeCellNo , &h_Cell[LocalCellNo] , h_pMapping , 0 ) ;
			
		
			if ( FaceNo == -1 ){
				XCoord	= h_Cell[LocalCellNo].XCenter ;
				YCoord	= h_Cell[LocalCellNo].YCenter ;
				ZCoord	= h_Cell[LocalCellNo].ZCenter ;
				ErrorTrackingNum++ ;
				break ;
			}			
			
			XCoord	= BeforeXCoord ;
			YCoord	= BeforeYCoord ;
			ZCoord	= BeforeZCoord ;
			
			
			// Move particle to the face of cell.
			MoveParticle3D( &XCoord , &YCoord , &ZCoord , XVel , YVel , ZVel , CollideTime ) ;
			
			
			// Calculate remainder time after moving particle onto face of cell.
			RemainderTime	= Time - CollideTime ;
			
			
			NeighborCellNo	= h_Cell[LocalCellNo].Neighbor[FaceNo] ;
			
			
			if ( NeighborCellNo < 0 ){
				Tracking	= ParticleCollideSurface3D( &XVel , &YVel , &ZVel , &h_Cell[LocalCellNo] , FaceNo , &i , &ParticleNum , &h_Species[SpeciesNo] , 
									    h_pDomain , h_Node , h_WallType , h_Surface , h_pDSMC , h_pMapping ) ;
									    
			}else{
				ProcessorNo	= h_pProcessor->CellProcessorNo[NeighborCellNo] ;
				
				// Particle still stay the same processor.
				if ( ProcessorNo == MPIMyID ){
					// Translate the global neighbor cell no. to local.
					NeighborCellNo	= h_pProcessor->LocalCellNo[NeighborCellNo] ;
				
					RemainderTime	= RemainderTime * (h_Cell[NeighborCellNo].Timestep/h_Cell[LocalCellNo].Timestep) ;
				
					// Update partice next cell no.
					BeforeCellNo	= h_Cell[LocalCellNo].Id ;
					LocalCellNo	= NeighborCellNo ;
				
				// Particle go to other processor.
				}else{
					BeforeCellNo			= h_Cell[LocalCellNo].Id ;
					h_pDSMC->ParticleXCoord[i]	= XCoord ;
					h_pDSMC->ParticleYCoord[i]	= YCoord ;
					h_pDSMC->ParticleZCoord[i]	= ZCoord ;
					h_pDSMC->ParticleXVel[i]	= XVel ;
					h_pDSMC->ParticleYVel[i]	= YVel ;
					h_pDSMC->ParticleZVel[i]	= ZVel ;
					
					TransferParticleToBuffer( h_pDSMC , &i , &ParticleNum , h_MPIParticleIn , pMPIParticleNum , ProcessorNo , 
								  NeighborCellNo , BeforeCellNo , RemainderTime , h_Cell[LocalCellNo].Timestep ) ;
					
					Tracking	= false ;
				}
			}
			
			TrackingNum++ ;
			
			// Debug.
			if ( TrackingNum == 6 ){
				XCoord	= h_Cell[LocalCellNo].XCenter ;
				YCoord	= h_Cell[LocalCellNo].YCenter ;
				ZCoord	= h_Cell[LocalCellNo].ZCenter ;
				break ;	
			}
			
		}while ( Tracking ) ;
		
		// Update particle information after tracking.
		if ( Tracking ){
			h_pDSMC->ParticleCellNo[i]	= LocalCellNo ;
			h_pDSMC->ParticleXCoord[i]	= XCoord ;
			h_pDSMC->ParticleYCoord[i]	= YCoord ;
			h_pDSMC->ParticleZCoord[i]	= ZCoord ;
			h_pDSMC->ParticleXVel[i]	= XVel ;
			h_pDSMC->ParticleYVel[i]	= YVel ;
			h_pDSMC->ParticleZVel[i]	= ZVel ;
		}
	}
	
	h_pDomain->ParticleNum		= ParticleNum ;	
	h_pDomain->TrackingNum		+=ParticleNum ;
	h_pDomain->ErrorTrackingNum	+=ErrorTrackingNum ;
	h_pDomain->ErrorTrackingNumNew	+=ErrorTrackingNum ;
	
	//OutputDebug << "2. MoveNewParticle: " << ", Num: " << DebugValue1 << ", ERROR: " << ErrorTrackingNum << endl ;
}


void MoveNewParticleAxisymmetric(	DSMC_DOMAIN		*h_pDomain ,
					DSMC_PROCESSOR		*h_pProcessor ,
					DSMC_NODE		*h_Node ,
					DSMC_CELL		*h_Cell ,
					DSMC_WALLTYPE		*h_WallType ,
					DSMC_SPECIES		*h_Species ,
					DSMC_SURFACE		*h_Surface ,
					DSMC_DSMC		*h_pDSMC ,
					CELLMAPPING		*h_pMapping , 
					int 			ParticleNo ,
					DSMC_MPI_PARTICLE	*h_MPIParticleIn ,
					int			*pMPIParticleNum , 
					ofstream		&OutputDebug ){
				
	int		MPISize , MPIMyID ;
	int		ParticleNum , SpeciesNo , FaceNo , ProcessorNo , WallNo ;
	int		LocalCellNo , BeforeCellNo , NeighborCellNo ;
	int		TrackingNum , ErrorTrackingNum , DebugValue1 ;
	double		RemainderTime , Time , CollideTime ;
	double		XCoord , YCoord , yCoord , XVel , YVel , ZVel , BufferYVel ; 
	double		BeforeXCoord , BeforeYCoord , BeforeXVel , BeforeYVel , BeforeZVel ;
	bool		Tracking ;
	
	
	MPI_Comm_size( MPI_COMM_WORLD , &MPISize ) ;
	MPI_Comm_rank( MPI_COMM_WORLD , &MPIMyID ) ;
	
	ParticleNum	= h_pDomain->ParticleNum ;
				
	Tracking	= true ;
	ErrorTrackingNum= 0 ;
	DebugValue1	= ParticleNum - ParticleNo ;	
	
	
	for ( int i=ParticleNo ; i<ParticleNum ; i++ ){
		Tracking	= true ;
		BeforeCellNo	= -1 ;
		
		LocalCellNo	= h_pDSMC->ParticleCellNo[i] ;
		SpeciesNo	= h_pDSMC->ParticleSpeciesNo[i] ;
		XCoord		= h_pDSMC->ParticleXCoord[i] ;
		YCoord		= h_pDSMC->ParticleYCoord[i] ;
		XVel		= h_pDSMC->ParticleXVel[i] ;
		YVel		= h_pDSMC->ParticleYVel[i] ;
		ZVel		= h_pDSMC->ParticleZVel[i] ;
		
		RemainderTime	= h_Cell[LocalCellNo].Timestep*Randn() ;
		TrackingNum	= 0 ;

		// Start to track particle position after one timestep.
		do{	
			Time		= RemainderTime ;
			BeforeXCoord	= XCoord ;
			BeforeYCoord	= YCoord ;
			BeforeXVel	= XVel ;
			BeforeYVel	= YVel ;
			BeforeZVel	= ZVel ;
			
			
			// Update the particle position directly.
			yCoord = MoveParticleAxisymmetric( &XCoord , &YCoord , &XVel , &YVel , &ZVel , Time ) ;
			
			
			if ( InCell2D( h_Node , &h_Cell[LocalCellNo] , XCoord , YCoord , h_pMapping , 1 ) )
				break ;


			CalculateTimeCollideFace2D( &CollideTime , &FaceNo , Time , XCoord , YCoord , BeforeXCoord , BeforeYCoord , BeforeCellNo , &h_Cell[LocalCellNo] , h_pMapping , 43 ) ;
		
		
			if ( FaceNo == -1 ){
				XCoord	= h_Cell[LocalCellNo].XCenter ;
				YCoord	= h_Cell[LocalCellNo].YCenter ;
				ErrorTrackingNum++ ;
				break ;
			}
			
			
			// Calculate remainder time after moving particle onto face of cell.
			RemainderTime	= Time - CollideTime ;
			NeighborCellNo	= h_Cell[LocalCellNo].Neighbor[FaceNo] ;
			
			/*if ( yCoord < 0. && NeighborCellNo >= 0 ){
				ProcessorNo	= h_pProcessor->CellProcessorNo[NeighborCellNo] ;
				RemainderTime	= 0. ;
				
				// Particle still stay the same processor.
				if ( ProcessorNo == MPIMyID ){
					// Translate the global neighbor cell no. to local.
					NeighborCellNo	= h_pProcessor->LocalCellNo[NeighborCellNo] ;
					
					// Update partice next cell no.
					BeforeCellNo	= h_Cell[LocalCellNo].Id ;
					LocalCellNo	= NeighborCellNo ;
					
				// Particle go to other processor.
				}else{
					BeforeCellNo			= h_Cell[LocalCellNo].Id ;
					h_pDSMC->ParticleXCoord[i]	= XCoord ;
					h_pDSMC->ParticleYCoord[i]	= YCoord ;
					h_pDSMC->ParticleXVel[i]	= XVel ;
					h_pDSMC->ParticleYVel[i]	= YVel ;
					h_pDSMC->ParticleZVel[i]	= ZVel ;
					
					TransferParticleToBuffer( h_pDSMC , &i , &ParticleNum , h_MPIParticleIn , pMPIParticleNum , ProcessorNo , 
								  NeighborCellNo , BeforeCellNo , RemainderTime , h_Cell[LocalCellNo].Timestep ) ;
					
					Tracking	= false ;
				}
				break ;	
			}*/
			
			
			BufferYVel	= (YCoord - BeforeYCoord)/Time ;
			XCoord		= BeforeXCoord ;
			YCoord		= BeforeYCoord ;
			XVel		= BeforeXVel ;
			YVel		= BeforeYVel ;
			ZVel		= BeforeZVel ;
			
			
			// Move particle to the face of cell.
			MoveParticle2D( &XCoord , &YCoord , XVel , BufferYVel , CollideTime ) ;
			
			
			if ( NeighborCellNo < 0 ){
				
				AdjustVelocityAxisymmetric( BeforeYCoord , &YVel , &ZVel , CollideTime ) ;
				
				Tracking	= ParticleCollideSurface2D( &XVel , &YVel , &ZVel , &h_Cell[LocalCellNo] , FaceNo , &i , &ParticleNum , &h_Species[SpeciesNo] , 
									    h_pDomain , h_Node , h_WallType , h_Surface , h_pDSMC , h_pMapping , OutputDebug ) ;
									    
			}else{
				AdjustVelocityAxisymmetric( BeforeYCoord , &YVel , &ZVel , CollideTime ) ;
				
				ProcessorNo	= h_pProcessor->CellProcessorNo[NeighborCellNo] ;
				
				// Particle still stay the same processor.
				if ( ProcessorNo == MPIMyID ){
					// Translate the global neighbor cell no. to local.
					NeighborCellNo	= h_pProcessor->LocalCellNo[NeighborCellNo] ;
					
					RemainderTime	= RemainderTime * (h_Cell[NeighborCellNo].Timestep/h_Cell[LocalCellNo].Timestep) ;
				
					// Update partice next cell no.
					BeforeCellNo	= h_Cell[LocalCellNo].Id ;
					LocalCellNo	= NeighborCellNo ;
					
				// Particle go to other processor.
				}else{
					BeforeCellNo			= h_Cell[LocalCellNo].Id ;
					h_pDSMC->ParticleXCoord[i]	= XCoord ;
					h_pDSMC->ParticleYCoord[i]	= YCoord ;
					h_pDSMC->ParticleXVel[i]	= XVel ;
					h_pDSMC->ParticleYVel[i]	= YVel ;
					h_pDSMC->ParticleZVel[i]	= ZVel ;
					
					TransferParticleToBuffer( h_pDSMC , &i , &ParticleNum , h_MPIParticleIn , pMPIParticleNum , ProcessorNo , 
								  NeighborCellNo , BeforeCellNo , RemainderTime , h_Cell[LocalCellNo].Timestep ) ;
					
					Tracking	= false ;
				}
			}
			
			TrackingNum++ ;
			
			// Debug.
			if ( TrackingNum == 6 ){
				XCoord	= h_Cell[LocalCellNo].XCenter ;
				YCoord	= h_Cell[LocalCellNo].YCenter ;
				ErrorTrackingNum++ ;
				break ;	
			}
			
		}while ( Tracking ) ;
		
		// Update particle information after tracking.
		if ( Tracking ){
			h_pDSMC->ParticleCellNo[i]	= LocalCellNo ;
			h_pDSMC->ParticleXCoord[i]	= XCoord ;
			h_pDSMC->ParticleYCoord[i]	= YCoord ;
			h_pDSMC->ParticleXVel[i]	= XVel ;
			h_pDSMC->ParticleYVel[i]	= YVel ;
			h_pDSMC->ParticleZVel[i]	= ZVel ;
		}
	}
	
	h_pDomain->ParticleNum		= ParticleNum ;
	h_pDomain->TrackingNum		+=ParticleNum ;
	h_pDomain->ErrorTrackingNum	+=ErrorTrackingNum ;
	h_pDomain->ErrorTrackingNumNew	+=ErrorTrackingNum ;
	
	// Debug.
	//OutputDebug << "2. MoveNewParticle: " << ", Num: "<< setw(20) << DebugValue1 << ", ERROR: " << setw(20) << ErrorTrackingNum << endl ;
}

//==============================================================================================================

void MoveOtherParticle2D(	DSMC_DOMAIN		*h_pDomain ,
				DSMC_PROCESSOR		*h_pProcessor ,
				DSMC_NODE		*h_Node ,
				DSMC_CELL		*h_Cell ,
				DSMC_WALLTYPE		*h_WallType ,
				DSMC_SPECIES		*h_Species ,
				DSMC_SURFACE		*h_Surface ,
				DSMC_DSMC		*h_pDSMC ,
				CELLMAPPING		*h_pMapping , 
				DSMC_MPI_PARTICLE	*h_MPIParticleOut , 
				DSMC_MPI_PARTICLE	*h_MPIParticleIn ,
				int			*MPIParticleNumOut ,
				int			*MPIParticleNumIn ,
				int			*pMPIParticleNum ,
				DSMC_MPI_DATATYPE	*pMPIDataType ,
				ofstream		&OutputDebug ){
	
	int		MPISize , MPIMyID ;
	int		ParticleNo , ParticleNum , SpeciesNo , FaceNo , ProcessorNo , TotalTransferParticleNum , DebugValue1 , DebugValue2 ;
	int		LocalCellNo , BeforeCellNo , NeighborCellNo ;
	int		TrackingNum , ErrorTrackingNum , TransferNum ;
	double		RemainderTime , Time , CollideTime ;
	double		XCoord , YCoord , XVel , YVel , ZVel , BeforeXCoord , BeforeYCoord ;
	bool		Tracking , Moving ;
	
	
	MPI_Comm_size( MPI_COMM_WORLD , &MPISize ) ;
	MPI_Comm_rank( MPI_COMM_WORLD , &MPIMyID ) ;
	
	ParticleNum	= h_pDomain->ParticleNum ;
				
	Tracking			= true ;
	Moving				= true ;
	TransferNum			= 0 ;
	ErrorTrackingNum		= 0 ;
	TotalTransferParticleNum	= 0 ;
	DebugValue1			= 0 ;
	DebugValue2			= 0 ;
		

	do{
		TransferNum++ ;
		
		// Transfer particle data from buffer memory to other processor's memory.
		MPITransferParticle( h_pDomain , h_pProcessor , h_MPIParticleOut , h_MPIParticleIn , MPIParticleNumOut , MPIParticleNumIn , 
				     pMPIParticleNum , pMPIDataType , OutputDebug ) ;
				     
		
		// Add the particles from other processors.
		ParticleNo	= AddParticleFromOtherProcessor( h_pDomain , h_pProcessor , h_pDSMC , h_Cell , h_MPIParticleIn , pMPIParticleNum , OutputDebug ) ;
		ParticleNum	= h_pDomain->ParticleNum ;
		
		// Debug.
		//OutputDebug << "ParticleNo: " << ParticleNo << ", ParticleNum: " << ParticleNum << '\n' ;

		
		// Debug.
		if ( TransferNum == 1 ) DebugValue1 = ParticleNum - ParticleNo ;
		
		
		for ( int i=ParticleNo ; i<ParticleNum ; i++ ){
			Tracking	= true ;
			BeforeCellNo	= h_pDSMC->ParticleBeforeCellNo[i] ;
		
			LocalCellNo	= h_pDSMC->ParticleCellNo[i] ;
			SpeciesNo	= h_pDSMC->ParticleSpeciesNo[i] ;
			XCoord		= h_pDSMC->ParticleXCoord[i] ;
			YCoord		= h_pDSMC->ParticleYCoord[i] ;
			XVel		= h_pDSMC->ParticleXVel[i] ;
			YVel		= h_pDSMC->ParticleYVel[i] ;
			ZVel		= h_pDSMC->ParticleZVel[i] ;
		
			RemainderTime	= h_pDSMC->ParticleTimestep[i] ;
			TrackingNum	= 0 ;


			// Start to track particle position after one timestep.
			do{	
				Time		= RemainderTime ;
				BeforeXCoord	= XCoord ;
				BeforeYCoord	= YCoord ;

			
				// Update the particle position directly.
				MoveParticle2D( &XCoord , &YCoord , XVel , YVel , Time ) ;
			
			
				if ( InCell2D( h_Node , &h_Cell[LocalCellNo] , XCoord , YCoord , h_pMapping , 1 ) )
					break ;


				CalculateTimeCollideFace2D( &CollideTime , &FaceNo , Time , XCoord , YCoord , BeforeXCoord , BeforeYCoord , BeforeCellNo , &h_Cell[LocalCellNo] , h_pMapping , 43 ) ;
		
				// Debug.
				DebugValue2++ ;
		
		
				if ( FaceNo == -1 ){
					// Debug.
					//OutputDebug << i << ", TN: " << TrackingNum << ", Timestep: " << RemainderTime << ", t/tmax: " << RemainderTime/h_Cell[LocalCellNo].Timestep << endl ;
					//OutputDebug << "Before CellNo: " << BeforeCellNo << ", LocalCellNo: " << LocalCellNo << ", GlobalCellNo: " << h_Cell[LocalCellNo].Id << endl ;
					//OutputDebug << InCell2D( h_Node , &h_Cell[LocalCellNo] , XCoord , YCoord , h_pMapping , 1 ) << ", " << InCell2D( h_Node , &h_Cell[LocalCellNo] , BeforeXCoord , BeforeYCoord , h_pMapping , 1 ) << endl ;
					
					
					XCoord	= h_Cell[LocalCellNo].XCenter ;
					YCoord	= h_Cell[LocalCellNo].YCenter ;
					ErrorTrackingNum++ ;
					break ;
				}			
			
				XCoord	= BeforeXCoord ;
				YCoord	= BeforeYCoord ;
			
				// Move particle to the face of cell.
				MoveParticle2D( &XCoord , &YCoord , XVel , YVel , CollideTime ) ;
			
				// Calculate remainder time after moving particle onto face of cell.
				RemainderTime	= Time - CollideTime ;
			
				NeighborCellNo	= h_Cell[LocalCellNo].Neighbor[FaceNo] ;
			
				if ( NeighborCellNo < 0 ){
					Tracking	= ParticleCollideSurface2D( &XVel , &YVel , &ZVel , &h_Cell[LocalCellNo] , FaceNo , &i , &ParticleNum , &h_Species[SpeciesNo] , 
										    h_pDomain , h_Node , h_WallType , h_Surface , h_pDSMC , h_pMapping , OutputDebug ) ;
									    
				}else{
					ProcessorNo	= h_pProcessor->CellProcessorNo[NeighborCellNo] ;
				
					// Particle still stay the same processor.
					if ( ProcessorNo == MPIMyID ){
						// Translate the global neighbor cell no. to local.
						NeighborCellNo	= h_pProcessor->LocalCellNo[NeighborCellNo] ;
					
						RemainderTime	= RemainderTime * (h_Cell[NeighborCellNo].Timestep/h_Cell[LocalCellNo].Timestep) ;
				
						// Update partice next cell no.
						BeforeCellNo	= h_Cell[LocalCellNo].Id ;
						LocalCellNo	= NeighborCellNo ;
					
					// Particle go to other processor.
					}else{
						BeforeCellNo			= h_Cell[LocalCellNo].Id ;
						h_pDSMC->ParticleXCoord[i]	= XCoord ;
						h_pDSMC->ParticleYCoord[i]	= YCoord ;
						h_pDSMC->ParticleXVel[i]	= XVel ;
						h_pDSMC->ParticleYVel[i]	= YVel ;
						h_pDSMC->ParticleZVel[i]	= ZVel ;
					
						TransferParticleToBuffer( h_pDSMC , &i , &ParticleNum , h_MPIParticleIn , pMPIParticleNum , ProcessorNo , 
									  NeighborCellNo , BeforeCellNo , RemainderTime , h_Cell[LocalCellNo].Timestep ) ;
					
						Tracking	= false ;
					}
				}
			
			
				TrackingNum++ ;
			
				// Debug.
				if ( TrackingNum == 6 ){
					XCoord	= h_Cell[LocalCellNo].XCenter ;
					YCoord	= h_Cell[LocalCellNo].YCenter ;
					ErrorTrackingNum++ ;
					break ;	
				}
			
			}while ( Tracking ) ;
		
			// Update particle information after tracking.
			if ( Tracking ){
				h_pDSMC->ParticleCellNo[i]	= LocalCellNo ;
				h_pDSMC->ParticleXCoord[i]	= XCoord ;
				h_pDSMC->ParticleYCoord[i]	= YCoord ;
				h_pDSMC->ParticleXVel[i]	= XVel ;
				h_pDSMC->ParticleYVel[i]	= YVel ;
				h_pDSMC->ParticleZVel[i]	= ZVel ;
			}
		}
		
		h_pDomain->ParticleNum		= ParticleNum ;
		
		
		// Sum the transfered particle from each processor.
		MPI_Allreduce( pMPIParticleNum , &TotalTransferParticleNum , 1 , MPI_INT , MPI_SUM , MPI_COMM_WORLD ) ;
		
		
		//if ( TotalTransferParticleNum == 0 || TransferNum == 5 ){
		if ( TotalTransferParticleNum == 0 ){
			// Debug.
			//cout << "Fail to Tansfer Moveing! " << ParticleNo << ", " << ParticleNum << ", " << (*pMPIParticleNum) << '\n' ;
			
			Moving = false ;
		}
	}while ( Moving ) ;
	
	
	// Debug.
	//OutputDebug << "3. MoveOtherParticle: " << ", Num: "<< setw(20) << DebugValue1 << ", CheckNum: " << DebugValue2 << ", ERROR: " << setw(20) << ErrorTrackingNum << endl ;
}


void MoveOtherParticle3D(	DSMC_DOMAIN		*h_pDomain ,
				DSMC_PROCESSOR		*h_pProcessor ,
				DSMC_NODE		*h_Node ,
				DSMC_CELL		*h_Cell ,
				DSMC_WALLTYPE		*h_WallType ,
				DSMC_SPECIES		*h_Species ,
				DSMC_SURFACE		*h_Surface ,
				DSMC_DSMC		*h_pDSMC ,
				CELLMAPPING		*h_pMapping , 
				DSMC_MPI_PARTICLE	*h_MPIParticleOut , 
				DSMC_MPI_PARTICLE	*h_MPIParticleIn ,
				int			*MPIParticleNumOut ,
				int			*MPIParticleNumIn ,
				int			*pMPIParticleNum ,
				DSMC_MPI_DATATYPE	*pMPIDataType ,
				ofstream		&OutputDebug ){
	
	int		MPISize , MPIMyID ;
	int		ParticleNo , ParticleNum , SpeciesNo , FaceNo , ProcessorNo , TotalTransferParticleNum , DebugValue1 , DebugValue2 ;
	int		LocalCellNo , BeforeCellNo , NeighborCellNo ;
	int		TrackingNum , ErrorTrackingNum , TransferNum ;
	double		RemainderTime , Time , CollideTime ;
	double		XCoord , YCoord , ZCoord , XVel , YVel , ZVel , BeforeXCoord , BeforeYCoord , BeforeZCoord ;
	bool		Tracking , Moving ;
	
	
	MPI_Comm_size( MPI_COMM_WORLD , &MPISize ) ;
	MPI_Comm_rank( MPI_COMM_WORLD , &MPIMyID ) ;
	
	ParticleNum	= h_pDomain->ParticleNum ;
				
	Tracking			= true ;
	Moving				= true ;
	TransferNum			= 0 ;
	ErrorTrackingNum		= 0 ;
	TotalTransferParticleNum	= 0 ;
	DebugValue1			= 0 ;
	DebugValue2			= 0 ;
		

	do{
		TransferNum++ ;
		
		// Transfer particle data from buffer memory to other processor's memory.
		MPITransferParticle( h_pDomain , h_pProcessor , h_MPIParticleOut , h_MPIParticleIn , MPIParticleNumOut , MPIParticleNumIn , 
				     pMPIParticleNum , pMPIDataType , OutputDebug ) ;
				     
		
		// Add the particles from other processors.
		ParticleNo	= AddParticleFromOtherProcessor( h_pDomain , h_pProcessor , h_pDSMC , h_Cell , h_MPIParticleIn , pMPIParticleNum , OutputDebug ) ;
		ParticleNum	= h_pDomain->ParticleNum ;
		
		// Debug.
		//OutputDebug << "ParticleNo: " << ParticleNo << ", ParticleNum: " << ParticleNum << '\n' ;

		
		// Debug.
		if ( TransferNum == 1 ) DebugValue1 = ParticleNum - ParticleNo ;
		
		
		for ( int i=ParticleNo ; i<ParticleNum ; i++ ){
			Tracking	= true ;
			BeforeCellNo	= h_pDSMC->ParticleBeforeCellNo[i] ;
		
			LocalCellNo	= h_pDSMC->ParticleCellNo[i] ;
			SpeciesNo	= h_pDSMC->ParticleSpeciesNo[i] ;
			XCoord		= h_pDSMC->ParticleXCoord[i] ;
			YCoord		= h_pDSMC->ParticleYCoord[i] ;
			ZCoord		= h_pDSMC->ParticleZCoord[i] ;
			XVel		= h_pDSMC->ParticleXVel[i] ;
			YVel		= h_pDSMC->ParticleYVel[i] ;
			ZVel		= h_pDSMC->ParticleZVel[i] ;
		
			RemainderTime	= h_pDSMC->ParticleTimestep[i] ;
			TrackingNum	= 0 ;


			// Start to track particle position after one timestep.
			do{	
				Time		= RemainderTime ;
				BeforeXCoord	= XCoord ;
				BeforeYCoord	= YCoord ;
				BeforeZCoord	= ZCoord ;

			
				// Update the particle position directly.
				MoveParticle3D( &XCoord , &YCoord , &ZCoord , XVel , YVel , ZVel , Time ) ;
			
			
				if ( InCell3D( h_Node , &h_Cell[LocalCellNo] , XCoord , YCoord , ZCoord , h_pMapping , 0 ) )
					break ;


				CalculateTimeCollideFace3D( &CollideTime , &FaceNo , Time , XCoord , YCoord , ZCoord , BeforeXCoord , BeforeYCoord , 
							BeforeZCoord , BeforeCellNo , &h_Cell[LocalCellNo] , h_pMapping , 0 ) ;
		
		
				// Debug.
				DebugValue2++ ;
		
		
				if ( FaceNo == -1 ){
					XCoord	= h_Cell[LocalCellNo].XCenter ;
					YCoord	= h_Cell[LocalCellNo].YCenter ;
					ZCoord	= h_Cell[LocalCellNo].ZCenter ;
					ErrorTrackingNum++ ;
					break ;
				}			
			
				XCoord	= BeforeXCoord ;
				YCoord	= BeforeYCoord ;
				ZCoord	= BeforeZCoord ;
				
			
				// Move particle to the face of cell.
				MoveParticle3D( &XCoord , &YCoord , &ZCoord , XVel , YVel , ZVel , CollideTime ) ;
				
			
				// Calculate remainder time after moving particle onto face of cell.
				RemainderTime	= Time - CollideTime ;
			
				NeighborCellNo	= h_Cell[LocalCellNo].Neighbor[FaceNo] ;
			
				if ( NeighborCellNo < 0 ){
					Tracking	= ParticleCollideSurface3D( &XVel , &YVel , &ZVel , &h_Cell[LocalCellNo] , FaceNo , &i , &ParticleNum , &h_Species[SpeciesNo] , 
									  	    h_pDomain , h_Node , h_WallType , h_Surface , h_pDSMC , h_pMapping ) ;
									    
				}else{
					ProcessorNo	= h_pProcessor->CellProcessorNo[NeighborCellNo] ;
				
					// Particle still stay the same processor.
					if ( ProcessorNo == MPIMyID ){
						// Translate the global neighbor cell no. to local.
						NeighborCellNo	= h_pProcessor->LocalCellNo[NeighborCellNo] ;
					
						RemainderTime	= RemainderTime * (h_Cell[NeighborCellNo].Timestep/h_Cell[LocalCellNo].Timestep) ;
				
						// Update partice next cell no.
						BeforeCellNo	= h_Cell[LocalCellNo].Id ;
						LocalCellNo	= NeighborCellNo ;
					
					// Particle go to other processor.
					}else{
						BeforeCellNo			= h_Cell[LocalCellNo].Id ;
						h_pDSMC->ParticleXCoord[i]	= XCoord ;
						h_pDSMC->ParticleYCoord[i]	= YCoord ;
						h_pDSMC->ParticleZCoord[i]	= ZCoord ;
						h_pDSMC->ParticleXVel[i]	= XVel ;
						h_pDSMC->ParticleYVel[i]	= YVel ;
						h_pDSMC->ParticleZVel[i]	= ZVel ;
					
						TransferParticleToBuffer( h_pDSMC , &i , &ParticleNum , h_MPIParticleIn , pMPIParticleNum , ProcessorNo , 
									  NeighborCellNo , BeforeCellNo , RemainderTime , h_Cell[LocalCellNo].Timestep ) ;
					
						Tracking	= false ;
					}
				}
			
			
				TrackingNum++ ;
			
				// Debug.
				if ( TrackingNum == 6 ){
					XCoord	= h_Cell[LocalCellNo].XCenter ;
					YCoord	= h_Cell[LocalCellNo].YCenter ;
					ZCoord	= h_Cell[LocalCellNo].ZCenter ;
					ErrorTrackingNum++ ;
					break ;	
				}
			
			}while ( Tracking ) ;
		
			// Update particle information after tracking.
			if ( Tracking ){
				h_pDSMC->ParticleCellNo[i]	= LocalCellNo ;
				h_pDSMC->ParticleXCoord[i]	= XCoord ;
				h_pDSMC->ParticleYCoord[i]	= YCoord ;
				h_pDSMC->ParticleZCoord[i]	= ZCoord ;
				h_pDSMC->ParticleXVel[i]	= XVel ;
				h_pDSMC->ParticleYVel[i]	= YVel ;
				h_pDSMC->ParticleZVel[i]	= ZVel ;
			}
		}
		
		h_pDomain->ParticleNum		= ParticleNum ;
		
		
		// Sum the transfered particle from each processor.
		MPI_Allreduce( pMPIParticleNum , &TotalTransferParticleNum , 1 , MPI_INT , MPI_SUM , MPI_COMM_WORLD ) ;
		
		
		//if ( TotalTransferParticleNum == 0 || TransferNum == 5 ){
		if ( TotalTransferParticleNum == 0 ){
			// Debug.
			//cout << "Fail to Tansfer Moveing! " << ParticleNo << ", " << ParticleNum << ", " << (*pMPIParticleNum) << '\n' ;
			
			Moving = false ;
		}
	}while ( Moving ) ;
	
	
	// Debug.
	//OutputDebug << "3. MoveOtherParticle: " << ", Num: "<< setw(20) << DebugValue1 << ", CheckNum: " << DebugValue2 << ", ERROR: " << setw(20) << ErrorTrackingNum << endl ;
}


void MoveOtherParticleAxisymmetric(	DSMC_DOMAIN		*h_pDomain ,
					DSMC_PROCESSOR		*h_pProcessor ,
					DSMC_NODE		*h_Node ,
					DSMC_CELL		*h_Cell ,
					DSMC_WALLTYPE		*h_WallType ,
					DSMC_SPECIES		*h_Species ,
					DSMC_SURFACE		*h_Surface ,
					DSMC_DSMC		*h_pDSMC ,
					CELLMAPPING		*h_pMapping , 
					DSMC_MPI_PARTICLE	*h_MPIParticleOut , 
					DSMC_MPI_PARTICLE	*h_MPIParticleIn ,
					int			*MPIParticleNumOut ,
					int			*MPIParticleNumIn ,
					int			*pMPIParticleNum ,
					DSMC_MPI_DATATYPE	*pMPIDataType ,
					ofstream		&OutputDebug ){
	
	int		MPISize , MPIMyID ;
	int		ParticleNo , ParticleNum , SpeciesNo , FaceNo , ProcessorNo , TotalTransferParticleNum , DebugValue1 , DebugValue2 ;
	int		LocalCellNo , BeforeCellNo , NeighborCellNo , WallNo ;
	int		TrackingNum , ErrorTrackingNum , TransferNum ;
	double		RemainderTime , Time , CollideTime ;
	double		XCoord , YCoord , yCoord , XVel , YVel , ZVel , BufferYVel ;
	double		BeforeXCoord , BeforeYCoord , BeforeXVel , BeforeYVel , BeforeZVel ;
	bool		Tracking , Moving ;
	
	
	MPI_Comm_size( MPI_COMM_WORLD , &MPISize ) ;
	MPI_Comm_rank( MPI_COMM_WORLD , &MPIMyID ) ;
	
	ParticleNum	= h_pDomain->ParticleNum ;
				
	Tracking			= true ;
	Moving				= true ;
	TransferNum			= 0 ;
	ErrorTrackingNum		= 0 ;
	TotalTransferParticleNum	= 0 ;
	DebugValue1			= 0 ;
	DebugValue2			= 0 ;
		

	do{
		TransferNum++ ;
		
		// Transfer particle data from buffer memory to other processor's memory.
		MPITransferParticle( h_pDomain , h_pProcessor , h_MPIParticleOut , h_MPIParticleIn , MPIParticleNumOut , MPIParticleNumIn , 
				     pMPIParticleNum , pMPIDataType , OutputDebug ) ;
				     
		
		// Add the particles from other processors.
		ParticleNo	= AddParticleFromOtherProcessor( h_pDomain , h_pProcessor , h_pDSMC , h_Cell , h_MPIParticleIn , pMPIParticleNum , OutputDebug ) ;
		ParticleNum	= h_pDomain->ParticleNum ;
		
		// Debug.
		//OutputDebug << "ParticleNo: " << ParticleNo << ", ParticleNum: " << ParticleNum << '\n' ;

		
		// Debug.
		if ( TransferNum == 1 ) DebugValue1 = ParticleNum - ParticleNo ;
		
		
		for ( int i=ParticleNo ; i<ParticleNum ; i++ ){
			Tracking	= true ;
			BeforeCellNo	= h_pDSMC->ParticleBeforeCellNo[i] ;
		
			LocalCellNo	= h_pDSMC->ParticleCellNo[i] ;
			SpeciesNo	= h_pDSMC->ParticleSpeciesNo[i] ;
			XCoord		= h_pDSMC->ParticleXCoord[i] ;
			YCoord		= h_pDSMC->ParticleYCoord[i] ;
			XVel		= h_pDSMC->ParticleXVel[i] ;
			YVel		= h_pDSMC->ParticleYVel[i] ;
			ZVel		= h_pDSMC->ParticleZVel[i] ;
		
			RemainderTime	= h_pDSMC->ParticleTimestep[i] ;
			TrackingNum	= 0 ;


			// Start to track particle position after one timestep.
			do{	
				Time		= RemainderTime ;
				BeforeXCoord	= XCoord ;
				BeforeYCoord	= YCoord ;
				BeforeXVel	= XVel ;
				BeforeYVel	= YVel ;
				BeforeZVel	= ZVel ;

			
				// Update the particle position directly.
				yCoord = MoveParticleAxisymmetric( &XCoord , &YCoord , &XVel , &YVel , &ZVel , Time ) ;
			
			
				if ( InCell2D( h_Node , &h_Cell[LocalCellNo] , XCoord , YCoord , h_pMapping , 1 ) )
					break ;


				CalculateTimeCollideFace2D( &CollideTime , &FaceNo , Time , XCoord , YCoord , BeforeXCoord , BeforeYCoord , BeforeCellNo , &h_Cell[LocalCellNo] , h_pMapping , 43 ) ;
		
				// Debug.
				DebugValue2++ ;
		
		
				if ( FaceNo == -1 ){
					// Debug.
					//OutputDebug << i << ", TN: " << TrackingNum << ", Timestep: " << RemainderTime << ", t/tmax: " << RemainderTime/h_Cell[LocalCellNo].Timestep << endl ;
					//OutputDebug << "Before CellNo: " << BeforeCellNo << ", LocalCellNo: " << LocalCellNo << ", GlobalCellNo: " << h_Cell[LocalCellNo].Id << endl ;
					//OutputDebug << InCell2D( h_Node , &h_Cell[LocalCellNo] , XCoord , YCoord , h_pMapping , 1 ) << ", " << InCell2D( h_Node , &h_Cell[LocalCellNo] , BeforeXCoord , BeforeYCoord , h_pMapping , 1 ) << endl ;
					
					
					XCoord	= h_Cell[LocalCellNo].XCenter ;
					YCoord	= h_Cell[LocalCellNo].YCenter ;
					ErrorTrackingNum++ ;
					break ;
				}
				
				
				// Calculate remainder time after moving particle onto face of cell.
				RemainderTime	= Time - CollideTime ;
				NeighborCellNo	= h_Cell[LocalCellNo].Neighbor[FaceNo] ;		


				/*if ( yCoord < 0. && NeighborCellNo >= 0 ){
					ProcessorNo	= h_pProcessor->CellProcessorNo[NeighborCellNo] ;
					RemainderTime	= 0. ;
				
					// Particle still stay the same processor.
					if ( ProcessorNo == MPIMyID ){
						// Translate the global neighbor cell no. to local.
						NeighborCellNo	= h_pProcessor->LocalCellNo[NeighborCellNo] ;
					
						// Update partice next cell no.
						BeforeCellNo	= h_Cell[LocalCellNo].Id ;
						LocalCellNo	= NeighborCellNo ;
					
					// Particle go to other processor.
					}else{
						BeforeCellNo			= h_Cell[LocalCellNo].Id ;
						h_pDSMC->ParticleXCoord[i]	= XCoord ;
						h_pDSMC->ParticleYCoord[i]	= YCoord ;
						h_pDSMC->ParticleXVel[i]	= XVel ;
						h_pDSMC->ParticleYVel[i]	= YVel ;
						h_pDSMC->ParticleZVel[i]	= ZVel ;
						
						TransferParticleToBuffer( h_pDSMC , &i , &ParticleNum , h_MPIParticleIn , pMPIParticleNum , ProcessorNo , 
									  NeighborCellNo , BeforeCellNo , RemainderTime , h_Cell[LocalCellNo].Timestep ) ;
					
						Tracking	= false ;
					}
					break ;	
				}*/

			
				BufferYVel	= (YCoord - BeforeYCoord)/Time ;
				XCoord		= BeforeXCoord ;
				YCoord		= BeforeYCoord ;
				XVel		= BeforeXVel ;
				YVel		= BeforeYVel ;
				ZVel		= BeforeZVel ;
				
			
				// Move particle to the face of cell.
				MoveParticle2D( &XCoord , &YCoord , XVel , BufferYVel , CollideTime ) ;
				
			
				if ( NeighborCellNo < 0 ){
					
					AdjustVelocityAxisymmetric( BeforeYCoord , &YVel , &ZVel , CollideTime ) ;
					
					Tracking	= ParticleCollideSurface2D( &XVel , &YVel , &ZVel , &h_Cell[LocalCellNo] , FaceNo , &i , &ParticleNum , &h_Species[SpeciesNo] , 
										    h_pDomain , h_Node , h_WallType , h_Surface , h_pDSMC , h_pMapping , OutputDebug ) ;
		    
				}else{
					AdjustVelocityAxisymmetric( BeforeYCoord , &YVel , &ZVel , CollideTime ) ;
					
					ProcessorNo	= h_pProcessor->CellProcessorNo[NeighborCellNo] ;
				
					// Particle still stay the same processor.
					if ( ProcessorNo == MPIMyID ){
						// Translate the global neighbor cell no. to local.
						NeighborCellNo	= h_pProcessor->LocalCellNo[NeighborCellNo] ;
					
						RemainderTime	= RemainderTime * (h_Cell[NeighborCellNo].Timestep/h_Cell[LocalCellNo].Timestep) ;
				
						// Update partice next cell no.
						BeforeCellNo	= h_Cell[LocalCellNo].Id ;
						LocalCellNo	= NeighborCellNo ;
					
					// Particle go to other processor.
					}else{
						BeforeCellNo			= h_Cell[LocalCellNo].Id ;
						h_pDSMC->ParticleXCoord[i]	= XCoord ;
						h_pDSMC->ParticleYCoord[i]	= YCoord ;
						h_pDSMC->ParticleXVel[i]	= XVel ;
						h_pDSMC->ParticleYVel[i]	= YVel ;
						h_pDSMC->ParticleZVel[i]	= ZVel ;
					
						TransferParticleToBuffer( h_pDSMC , &i , &ParticleNum , h_MPIParticleIn , pMPIParticleNum , ProcessorNo , 
									  NeighborCellNo , BeforeCellNo , RemainderTime , h_Cell[LocalCellNo].Timestep ) ;
					
						Tracking	= false ;
					}
				}
			
			
				TrackingNum++ ;
			
				// Debug.
				if ( TrackingNum == 6 ){
					XCoord	= h_Cell[LocalCellNo].XCenter ;
					YCoord	= h_Cell[LocalCellNo].YCenter ;
					ErrorTrackingNum++ ;
					break ;	
				}
			
			}while ( Tracking ) ;
		
			// Update particle information after tracking.
			if ( Tracking ){
				h_pDSMC->ParticleCellNo[i]	= LocalCellNo ;
				h_pDSMC->ParticleXCoord[i]	= XCoord ;
				h_pDSMC->ParticleYCoord[i]	= YCoord ;
				h_pDSMC->ParticleXVel[i]	= XVel ;
				h_pDSMC->ParticleYVel[i]	= YVel ;
				h_pDSMC->ParticleZVel[i]	= ZVel ;
			}
		}
		
		h_pDomain->ParticleNum		= ParticleNum ;
		
		
		// Sum the transfered particle from each processor.
		MPI_Allreduce( pMPIParticleNum , &TotalTransferParticleNum , 1 , MPI_INT , MPI_SUM , MPI_COMM_WORLD ) ;
		
		
		//if ( TotalTransferParticleNum == 0 || TransferNum == 5 ){
		if ( TotalTransferParticleNum == 0 ){
			// Debug.
			//cout << "Fail to Tansfer Moveing! " << ParticleNo << ", " << ParticleNum << ", " << (*pMPIParticleNum) << '\n' ;
			
			Moving = false ;
		}
	}while ( Moving ) ;
}

//==============================================================================================================

void MoveParticle2D( double *pXCoord , double *pYCoord , double XVel , double YVel , double Time ){
	(*pXCoord)	+= XVel*Time ;
	(*pYCoord)	+= YVel*Time ;
}


void MoveParticle3D( double *pXCoord , double *pYCoord , double *pZCoord , double XVel , double YVel , double ZVel , double Time ){
	(*pXCoord)	+= XVel*Time ;
	(*pYCoord)	+= YVel*Time ;
	(*pZCoord)	+= ZVel*Time ;
}


double MoveParticleAxisymmetric( double *pXCoord , double *pYCoord , double *pXVel , double *pYVel , double *pZVel , double Time ){
	double		XCoord , YCoord , ZCoord , yCoord , XVel , YVel , ZVel ;
	
	XVel	= (*pXVel) ;
	YVel	= (*pYVel) ;
	ZVel	= (*pZVel) ;
	
	XCoord	= (*pXCoord) + XVel*Time ;
	YCoord	= (*pYCoord) + YVel*Time ;
	ZCoord	= ZVel*Time ;
	yCoord	= 0. ;
	
	(*pXCoord)	= XCoord ;
	(*pYCoord)	= sqrt((YCoord*YCoord) + (ZCoord*ZCoord)) ;
	
	(*pYVel)	= ((YVel*YCoord) + (ZVel*ZCoord))/(*pYCoord) ;
	(*pZVel)	= ((ZVel*YCoord) - (YVel*ZCoord))/(*pYCoord) ;
	
	//if ( YCoord < 0. ) (*pYCoord) *= -1. ;
	if ( YCoord < 0. ) yCoord = -1.*(*pYCoord) ;
		
	return	yCoord ;
}


void MoveParticleAxisymmetric( double *pXCoord , double *pYCoord , double XVel , double YVel , double ZVel , double Time ){
	double		YCoord , ZCoord ;
	
	YCoord	= (*pYCoord) + YVel*Time ;
	ZCoord	= ZVel*Time ;
	
	(*pXCoord)	= (*pXCoord) + XVel*Time ;
	(*pYCoord)	= sqrt((YCoord*YCoord) + (ZCoord*ZCoord)) ;
	
	if ( YCoord < 0. ) (*pYCoord) *= -1. ;
}

//==============================================================================================================

void AdjustVelocityAxisymmetric( double YCoord , double *pYVel , double *pZVel , double Time ){
	double		ZCoord , YVel , ZVel , DY ;
	
	YVel	= (*pYVel) ;
	ZVel	= (*pZVel) ;
	
	YCoord	= YCoord + YVel*Time ;
	ZCoord	= ZVel*Time ;
	
	DY	= sqrt((YCoord*YCoord) + (ZCoord*ZCoord)) ;
	
	(*pYVel)	= ((YVel*YCoord) + (ZVel*ZCoord))/DY ;
	(*pZVel)	= ((ZVel*YCoord) - (YVel*ZCoord))/DY ;
}

//==============================================================================================================

void CalculateTimeCollideFace2D( double		*pCollideTime , 
				 int		*pFaceNo ,
				 double		Time , 
				 double		XCoord , 
				 double		YCoord , 
				 double		BeforeXCoord , 
				 double		BeforeYCoord , 
				 int		BeforeCellNo ,
				 DSMC_CELL	*_Cell ,
				 CELLMAPPING	*h_pMapping ,
				 int		Debug ){
				 	
	int		Type , MinFaceNo ;
	double		Distance1 , Distance2 , CTime , MinCTime ;
	
	Type		= _Cell->Type ;	
					
	(*pCollideTime)	= 1.E+6 ;
	(*pFaceNo)	= -1 ;
	MinCTime	= 0. ;
	MinFaceNo	= -1 ;
	
	// Debug.
	/*if ( Debug == 43 ){
		cout << "BX: " << BeforeXCoord << ", BY: " << BeforeYCoord << '\n' ;
		cout << "X: " << XCoord << ", Y: " << YCoord << '\n' ;	
	}*/

	for ( int i=0 ; i<h_pMapping->SurfaceNum[Type] ; i++ ){
		Distance1	= BeforeXCoord*_Cell->FaceFA[i] + BeforeYCoord*_Cell->FaceFB[i] + _Cell->FaceFC[i] ;
		Distance2	= XCoord*_Cell->FaceFA[i] + YCoord*_Cell->FaceFB[i] + _Cell->FaceFC[i] ;
	
		if ( fabs(Distance1) == 0. || _Cell->Neighbor[i] == BeforeCellNo || _Cell->Neighbor[i] == -6 ){
		//	MinCTime	= Time*fabs( Distance1/(Distance1-Distance2) ) ;
		//	MinFaceNo	= i ;
			continue ;
		}
	
		// Debug.
		/*if ( Debug == 43 ){
			cout << "FaceNo: " << i << '\n' ;
			cout << "FA: " << _Cell->FaceFA[i] << ", FB: " << _Cell->FaceFB[i] << ", FC: " << _Cell->FaceFC[i] << '\n' ;
			cout << "D1: " << Distance1 << ", D2: " << Distance2 << '\n' ;	
		}*/
			
	
		if ( (Distance1*Distance2) < 0. ){
		//if ( (Distance1*Distance2) <= 0. ){
			CTime	= Time*fabs( Distance1/(Distance1-Distance2) ) ;
			
			if ( CTime < (*pCollideTime) ){
				(*pCollideTime)	= CTime ;
				(*pFaceNo)	= i ;
			}
		}
	}
}


void CalculateTimeCollideFace3D( double		*pCollideTime , 
				 int		*pFaceNo ,
				 double		Time , 
				 double		XCoord , 
				 double		YCoord , 
				 double		ZCoord ,
				 double		BeforeXCoord , 
				 double		BeforeYCoord , 
				 double		BeforeZCoord ,
				 int		BeforeCellNo ,
				 DSMC_CELL	*_Cell ,
				 CELLMAPPING	*h_pMapping ,
				 int		Debug ){
				 	
	int		Type , MinFaceNo ;
	double		Distance1 , Distance2 , CTime , MinCTime ;
	
	Type		= _Cell->Type ;	
					
	(*pCollideTime)	= 1.E+6 ;
	(*pFaceNo)	= -1 ;
	MinCTime	= 0. ;
	MinFaceNo	= -1 ;
	

	for ( int i=0 ; i<h_pMapping->SurfaceNum[Type] ; i++ ){
		Distance1	= BeforeXCoord*_Cell->FaceFA[i] + BeforeYCoord*_Cell->FaceFB[i] + BeforeZCoord*_Cell->FaceFC[i] + _Cell->FaceFD[i] ;
		Distance2	= XCoord*_Cell->FaceFA[i] + YCoord*_Cell->FaceFB[i] + ZCoord*_Cell->FaceFC[i] + _Cell->FaceFD[i] ;
	
		if ( fabs(Distance1) == 0. || _Cell->Neighbor[i] == BeforeCellNo )
			continue ;
	
		if ( (Distance1*Distance2) < 0. ){
			CTime	= Time*fabs( Distance1/(Distance1-Distance2) ) ;
			
			if ( CTime < (*pCollideTime) ){
				(*pCollideTime)	= CTime ;
				(*pFaceNo)	= i ;
			}
		}
	}
}

//==============================================================================================================

bool ParticleCollideSurface2D(	double		*pXVel ,
				double		*pYVel ,
				double		*pZVel ,
				DSMC_CELL	*_Cell ,
				int		FaceNo ,
				int		*pParticleNo ,
				int		*pParticleNum ,
				DSMC_SPECIES	*_pSpecies ,
				DSMC_DOMAIN	*h_pDomain ,
				DSMC_NODE	*h_Node ,
				DSMC_WALLTYPE	*h_WallType ,
				DSMC_SURFACE	*h_Surface ,
				DSMC_DSMC	*h_pDSMC ,
				CELLMAPPING	*h_pMapping , 
				ofstream	&OutputDebug ){
					
	bool		tracking ;
	int		SurfaceNo , WallNo , SpeciesNo , Type , Node0 , Node1 ;
	double		NormVel , ParaVel , MostProbableSpeed , WallTemp, A, B;
	
	
	tracking	= true ;
	
	
	if ( _Cell->Neighbor[FaceNo] == -3 || _Cell->Neighbor[FaceNo] == -4 ){
		RemoveParticle( pParticleNo , pParticleNum , h_pDSMC ) ;
		tracking	= false ;
		
	}else if ( _Cell->Neighbor[FaceNo] < -20 ){
		SpeciesNo	= _pSpecies->Id ;
		SurfaceNo	= _Cell->Surface[FaceNo] ;
		WallNo		= h_Surface[SurfaceNo].WallNo ;
		
		Type		= _Cell->Type ;
		Node0		= _Cell->Node[h_pMapping->Node[Type][FaceNo][0]] ;
		Node1		= _Cell->Node[h_pMapping->Node[Type][FaceNo][1]] ;
		
		
		// Calculate normal and parallel velocities on the face.
		InverseVelocityToNormPara2D( &NormVel , &ParaVel , (*pXVel) , (*pYVel) , h_Node[Node0] , h_Node[Node1] ) ;
		
		
		// Sample the properties of the incident particle on surface.
		SamplingSurface( 1 , *pXVel , *pYVel , *pZVel , h_pDSMC->ParticleRotation[*pParticleNo] , h_pDSMC->ParticleVibration[*pParticleNo] , SurfaceNo , _pSpecies , h_pDSMC ) ;
		SamplingSurface( 3 , NormVel , NormVel , NormVel , h_pDSMC->ParticleRotation[*pParticleNo] , h_pDSMC->ParticleVibration[*pParticleNo] , SurfaceNo , _pSpecies , h_pDSMC ) ;
	
	
		// Calculate velocity and internal energy of the particle after reflection.
		if ( h_WallType[WallNo].Type == 1 ){
			// Fully-specular reflection.
			NormVel		= -1.*NormVel ;
			
			// Debug.
			//OutputDebug << "Before-XVel: " << (*pXVel) << ", YVel: " << (*pYVel) << endl ;
			
			// Calculate velocities in Cartesian coordinate system.
			InverseVelocityToCartesian2D( pXVel , pYVel , NormVel , ParaVel , h_Node[Node0] , h_Node[Node1] ) ;
			
			h_pDSMC->ParticleLastCollide[*pParticleNo]	= -1 ;
			
			// Debug.
			//OutputDebug << "After-XVel: " << (*pXVel) << ", YVel: " << (*pYVel) << endl ;
			
		}else if ( h_WallType[WallNo].Type == 2 ){
			WallTemp		= h_WallType[WallNo].Temp ;
			MostProbableSpeed	= sqrt(2.*BOLTZ*WallTemp/_pSpecies->Mass) ;
			
			
			// Fully-diffusive reflection.
			NormVel		= -sqrt(-log(Randn()))*MostProbableSpeed ;
			RandVelocity( &ParaVel , pZVel , MostProbableSpeed ) ;

			
			// Calculate velocities in Cartesian coordinate system.
			InverseVelocityToCartesian2D( pXVel , pYVel , NormVel , ParaVel , h_Node[Node0] , h_Node[Node1] ) ;


			(*pXVel)	+= h_WallType[WallNo].XVel ;
			(*pYVel)	+= h_WallType[WallNo].YVel ;


			// Set rotational energy.
			if ( _pSpecies->RotDOF > 0. )
				h_pDSMC->ParticleRotation[*pParticleNo]	= RotationalEnergy( WallTemp , _pSpecies->RotDOF ) ; 
				

			// Set vibrational energy.
//			if ( h_pDSMC->EffVibDOF[_Cell->LocalId][SpeciesNo]  > 0.){
			if (_pSpecies->VibMode	> 0. ) {
				
				h_pDSMC->ParticleVibration[*pParticleNo]= VibrationalEnergy( WallTemp , &h_pDSMC->ParticleVibLevel[*pParticleNo] , &h_pDSMC->ParticleEffTemp[*pParticleNo] ,
											     h_pDSMC->EffVibDOF[_Cell->LocalId][SpeciesNo] , h_pDomain , _pSpecies->VibTemp ) ;

											     
			}
			
			
			h_pDSMC->ParticleLastCollide[*pParticleNo]	= -1 ;
			
		}else if ( h_WallType[WallNo].Type == 3 ){
			
			
			// Isotropic scattering 	
			MostProbableSpeed =sqrt((NormVel)*(NormVel)+(ParaVel)*(ParaVel));		
			A = 6.283185308*Randn()  ;
			NormVel		= - abs (MostProbableSpeed*cos(A)) ;
			ParaVel   =   MostProbableSpeed*sin(A);
	
			// Isotropic scattering 
/*	
			MostProbableSpeed =sqrt((NormVel)*(NormVel)+(ParaVel)*(ParaVel)+(*pZVel)*(*pZVel));		
			A = Randn() ;
			NormVel		= - MostProbableSpeed*sqrt(A) ;
			B = 6.283185308*Randn() ;
			ParaVel  = MostProbableSpeed*sqrt(1-A)*cos(B);
			*pZVel  = MostProbableSpeed*sqrt(1-A)*sin(B);
*/			

			// Calculate velocities in Cartesian coordinate system.
			InverseVelocityToCartesian2D( pXVel , pYVel , NormVel , ParaVel , h_Node[Node0] , h_Node[Node1] ) ;
			
				(*pXVel)	+= h_WallType[WallNo].XVel ;
				(*pYVel)	+= h_WallType[WallNo].YVel ;

			
			h_pDSMC->ParticleLastCollide[*pParticleNo]	= -1 ;
			
			
		}else if ( h_WallType[WallNo].Type == 4 ){
			
			if ( Randn() <= h_WallType[WallNo].StickingCoef ){
				SamplingSurface( 5 , *pXVel , *pYVel , *pZVel , h_pDSMC->ParticleRotation[*pParticleNo] , h_pDSMC->ParticleVibration[*pParticleNo] , SurfaceNo , _pSpecies , h_pDSMC ) ;
				
				
				RemoveParticle( pParticleNo , pParticleNum , h_pDSMC ) ;
				tracking	= false ;
				return	tracking ;
				
			}else{
				WallTemp		= h_WallType[WallNo].Temp ;
				MostProbableSpeed	= sqrt(2.*BOLTZ*WallTemp/_pSpecies->Mass) ;
			
			
			// Fully-diffusive reflection.
				NormVel		= -sqrt(-log(Randn()))*MostProbableSpeed ;
				RandVelocity( &ParaVel , pZVel , MostProbableSpeed ) ;

			
			// Calculate velocities in Cartesian coordinate system.
				InverseVelocityToCartesian2D( pXVel , pYVel , NormVel , ParaVel , h_Node[Node0] , h_Node[Node1] ) ;


				(*pXVel)	+= h_WallType[WallNo].XVel ;
				(*pYVel)	+= h_WallType[WallNo].YVel ;


				// Set rotational energy.
				if ( _pSpecies->RotDOF > 0. )
					h_pDSMC->ParticleRotation[*pParticleNo]	= RotationalEnergy( WallTemp , _pSpecies->RotDOF ) ; 

				
				// Set vibrational energy.
//				if ( h_pDSMC->EffVibDOF[_Cell->LocalId][SpeciesNo]  > 0. ){
				if (_pSpecies->VibMode	> 0. ) {
					h_pDSMC->ParticleVibration[*pParticleNo]= VibrationalEnergy( WallTemp , &h_pDSMC->ParticleVibLevel[*pParticleNo] , &h_pDSMC->ParticleEffTemp[*pParticleNo] ,
											     h_pDSMC->EffVibDOF[_Cell->LocalId][SpeciesNo] , h_pDomain , _pSpecies->VibTemp ) ;
				}


				h_pDSMC->ParticleLastCollide[*pParticleNo]	= -1 ;
				
			}
			
		}else if ( h_WallType[WallNo].Type == 5 ){
			
			double  AlphaN , SigmaT , EintCoef , AlphaT , R , Theta , Un , Ux , Angle,VPI ,VNI ,UPI,WPI,CTH,OM;
			
		
			WallTemp		= h_WallType[WallNo].Temp ;
			AlphaN      = h_WallType[WallNo].AlphaN ;
			SigmaT      = h_WallType[WallNo].SigmaT ;
			EintCoef        = h_WallType[WallNo].EintCoef ;
		
			MostProbableSpeed	= sqrt(2.*BOLTZ*WallTemp/_pSpecies->Mass) ;
				
			VNI = abs(NormVel / MostProbableSpeed) ;
			UPI = ParaVel / MostProbableSpeed  ;
			WPI = *pZVel  / MostProbableSpeed ;
			
			Angle = atan2(WPI,UPI);
				
			VPI = sqrt(UPI*UPI+WPI*WPI);
														
			R= sqrt(-AlphaN*log(Randn()));			
			Theta =6.283185308*Randn() ;
			Un = VNI*sqrt(1.-AlphaN);		
				
			NormVel = -MostProbableSpeed * sqrt ( R*R+Un*Un+2.*R*Un* cos(Theta)) ; 
			
			AlphaT=SigmaT*(2.-SigmaT);
			
			R= sqrt(-AlphaT*log(Randn()));
			Theta =6.283185308*Randn() ;
			
			Ux = VPI*sqrt(1-AlphaT);
			ParaVel =  MostProbableSpeed *((Ux+R*cos(Theta))*cos(Angle)-(R*sin(Theta)*sin(Angle))) ;
			*pZVel  =  MostProbableSpeed *((Ux+R*cos(Theta))*sin(Angle)+(R*sin(Theta)*cos(Angle))) ;


			// Calculate velocities in Cartesian coordinate system.
			InverseVelocityToCartesian2D( pXVel , pYVel , NormVel , ParaVel , h_Node[Node0] , h_Node[Node1] ) ;
			
				(*pXVel)	+= h_WallType[WallNo].XVel ;
				(*pYVel)	+= h_WallType[WallNo].YVel ;
				
			// Set rotational energy.
      OM = sqrt ( h_pDSMC->ParticleRotation[*pParticleNo]*(1.-EintCoef)/ ( BOLTZ * WallTemp ));
	    	
				if ( _pSpecies->RotDOF == 2.){
					
        	R = sqrt(- EintCoef * log( Randn()));
        
        	CTH = cos (6.283185308*Randn());

				}else{
									
				}
				
			h_pDSMC->ParticleRotation[*pParticleNo] = BOLTZ * WallTemp*(R*R+OM*OM+2*R*OM*CTH);
			
/*			
			// Set rotational energy.			
				if ( _pSpecies->RotDOF > 0.)
					
					h_pDSMC->ParticleRotation[*pParticleNo]	= RotationalEnergy( WallTemp , _pSpecies->RotDOF ) ; 

*/				
				
				// Set vibrational energy.
//				if ( h_pDSMC->EffVibDOF[_Cell->LocalId][SpeciesNo]  > 0. ){
				if (_pSpecies->VibMode	> 0. ) {
					h_pDSMC->ParticleVibration[*pParticleNo]= VibrationalEnergy( WallTemp , &h_pDSMC->ParticleVibLevel[*pParticleNo] , &h_pDSMC->ParticleEffTemp[*pParticleNo] ,
											     h_pDSMC->EffVibDOF[_Cell->LocalId][SpeciesNo] , h_pDomain , _pSpecies->VibTemp ) ;
				}

				h_pDSMC->ParticleLastCollide[*pParticleNo]	= -1 ;
		}
		
		
		// Sample the properties of the reflective particle on surface.
		SamplingSurface( 2 , *pXVel , *pYVel , *pZVel , h_pDSMC->ParticleRotation[*pParticleNo] , h_pDSMC->ParticleVibration[*pParticleNo] , SurfaceNo , _pSpecies , h_pDSMC ) ;
		SamplingSurface( 4 , NormVel , NormVel , NormVel , h_pDSMC->ParticleRotation[*pParticleNo] , h_pDSMC->ParticleVibration[*pParticleNo] , SurfaceNo , _pSpecies , h_pDSMC ) ;
	}
	
	return	tracking ;
}


bool ParticleCollideSurface3D(	double		*pXVel ,
				double		*pYVel ,
				double		*pZVel ,
				DSMC_CELL	*_Cell ,
				int		FaceNo ,
				int		*pParticleNo ,
				int		*pParticleNum ,
				DSMC_SPECIES	*_pSpecies ,
				DSMC_DOMAIN	*h_pDomain ,
				DSMC_NODE	*h_Node ,
				DSMC_WALLTYPE	*h_WallType ,
				DSMC_SURFACE	*h_Surface ,
				DSMC_DSMC	*h_pDSMC ,
				CELLMAPPING	*h_pMapping ){
					
	bool		tracking ;
	int		SurfaceNo , WallNo , SpeciesNo , Type , Node0 , Node1 ;
	double		NormVel , ParaVelX , ParaVelY , MostProbableSpeed , WallTemp ,A ,B ;
	
	
	tracking	= true ;
	
	
	if ( _Cell->Neighbor[FaceNo] == -3 || _Cell->Neighbor[FaceNo] == -4 ){
		RemoveParticle( pParticleNo , pParticleNum , h_pDSMC ) ;
		tracking	= false ;
		
	}else if ( _Cell->Neighbor[FaceNo] < -20 ){
		SpeciesNo	= _pSpecies->Id ;
		SurfaceNo	= _Cell->Surface[FaceNo] ;
		WallNo		= h_Surface[SurfaceNo].WallNo ;
		
		Type		= _Cell->Type ;
		
		// Calculate normal and parallel velocities on the face.
		InverseVelocityToNormPara3D( &NormVel , &ParaVelX , &ParaVelY , (*pXVel) , (*pYVel) , (*pZVel) , h_Node , 
						_Cell , FaceNo , h_pMapping ) ;
		
		
		// Sample the properties of the incident particle on surface.
		SamplingSurface( 1 , *pXVel , *pYVel , *pZVel , h_pDSMC->ParticleRotation[*pParticleNo] , h_pDSMC->ParticleVibration[*pParticleNo] , SurfaceNo , _pSpecies , h_pDSMC ) ;
		SamplingSurface( 3 , NormVel , NormVel , NormVel , h_pDSMC->ParticleRotation[*pParticleNo] , h_pDSMC->ParticleVibration[*pParticleNo] , SurfaceNo , _pSpecies , h_pDSMC ) ;
	
	
		// Calculate velocity and internal energy of the particle after reflection.
		if ( h_WallType[WallNo].Type == 1 ){
			// Fully-specular reflection.
			NormVel		= -1.*NormVel ;
			
			
			// Calculate velocities in Cartesian coordinate system.
			InverseVelocityToCartesian3D( pXVel , pYVel , pZVel , NormVel , ParaVelX , ParaVelY , h_Node , 
							_Cell , FaceNo , h_pMapping ) ;
			
			h_pDSMC->ParticleLastCollide[*pParticleNo]	= -1 ;
			
		}else if ( h_WallType[WallNo].Type == 2 ){
			WallTemp		= h_WallType[WallNo].Temp ;
			MostProbableSpeed	= sqrt(2.*BOLTZ*WallTemp/_pSpecies->Mass) ;
			
			// Fully-diffusive reflection.
			NormVel		= -sqrt(-log(Randn()))*MostProbableSpeed ;
			RandVelocity( &ParaVelX , &ParaVelY , MostProbableSpeed ) ;


			// Calculate velocities in Cartesian coordinate system.
			InverseVelocityToCartesian3D( pXVel , pYVel , pZVel , NormVel , ParaVelX , ParaVelY , h_Node , 
							_Cell , FaceNo , h_pMapping ) ;


			(*pXVel)	+= h_WallType[WallNo].XVel ;
			(*pYVel)	+= h_WallType[WallNo].YVel ;
			(*pZVel)	+= h_WallType[WallNo].ZVel ;


			// Set rotational energy.
			if ( _pSpecies->RotDOF > 0. )
				h_pDSMC->ParticleRotation[*pParticleNo]	= RotationalEnergy( WallTemp , _pSpecies->RotDOF ) ; 

				
			// Set vibrational energy.
//			if ( h_pDSMC->EffVibDOF[_Cell->LocalId][SpeciesNo]  > 0. ){
			if (_pSpecies->VibMode	> 0. ) {
				h_pDSMC->ParticleVibration[*pParticleNo]= VibrationalEnergy( WallTemp , &h_pDSMC->ParticleVibLevel[*pParticleNo] , &h_pDSMC->ParticleEffTemp[*pParticleNo] ,
											     h_pDSMC->EffVibDOF[_Cell->LocalId][SpeciesNo] , h_pDomain , _pSpecies->VibTemp ) ;
			}


			h_pDSMC->ParticleLastCollide[*pParticleNo]	= -1 ;
			
		}else if ( h_WallType[WallNo].Type == 3 ){
			
				MostProbableSpeed =sqrt((NormVel)*(NormVel)+(ParaVelX)*(ParaVelX)+(ParaVelY)*(ParaVelY));
		
			// Isotropic scattering 
			
			A = Randn() ;
			NormVel		= -MostProbableSpeed*sqrt(A) ;
			B = 6.283185308*Randn() ;
			ParaVelX  = MostProbableSpeed*sqrt(1-A)*cos(B);
			ParaVelY  = MostProbableSpeed*sqrt(1-A)*sin(B);

			// Calculate velocities in Cartesian coordinate system.
			InverseVelocityToCartesian3D( pXVel , pYVel , pZVel , NormVel , ParaVelX , ParaVelY , h_Node , 
							_Cell , FaceNo , h_pMapping ) ;
			
				(*pXVel)	+= h_WallType[WallNo].XVel ;
				(*pYVel)	+= h_WallType[WallNo].YVel ;
				(*pZVel)	+= h_WallType[WallNo].ZVel ;

			
			h_pDSMC->ParticleLastCollide[*pParticleNo]	= -1 ;
			
			
		}else if ( h_WallType[WallNo].Type == 4 ){
			if ( Randn() <= h_WallType[WallNo].StickingCoef ){
				SamplingSurface( 5 , *pXVel , *pYVel , *pZVel , h_pDSMC->ParticleRotation[*pParticleNo] , h_pDSMC->ParticleVibration[*pParticleNo] , SurfaceNo , _pSpecies , h_pDSMC ) ;
				
				
				RemoveParticle( pParticleNo , pParticleNum , h_pDSMC ) ;
				tracking	= false ;
				return	tracking ;
				
			}else{
				WallTemp		= h_WallType[WallNo].Temp ;
				MostProbableSpeed	= sqrt(2.*BOLTZ*WallTemp/_pSpecies->Mass) ;
			
				// Fully-diffusive reflection.
				NormVel		= -sqrt(-log(Randn()))*MostProbableSpeed ;
				RandVelocity( &ParaVelX , &ParaVelY , MostProbableSpeed ) ;


				// Calculate velocities in Cartesian coordinate system.
				InverseVelocityToCartesian3D( pXVel , pYVel , pZVel , NormVel , ParaVelX , ParaVelY , h_Node , _Cell , FaceNo , h_pMapping ) ;


				(*pXVel)	+= h_WallType[WallNo].XVel ;
				(*pYVel)	+= h_WallType[WallNo].YVel ;
				(*pZVel)	+= h_WallType[WallNo].ZVel ;


				// Set rotational energy.
				if ( _pSpecies->RotDOF > 0. )
					h_pDSMC->ParticleRotation[*pParticleNo]	= RotationalEnergy( WallTemp , _pSpecies->RotDOF ) ; 

				
				// Set vibrational energy.
//				if ( h_pDSMC->EffVibDOF[_Cell->LocalId][SpeciesNo]  > 0. ){
				if (_pSpecies->VibMode	> 0. ) {
					h_pDSMC->ParticleVibration[*pParticleNo]= VibrationalEnergy( WallTemp , &h_pDSMC->ParticleVibLevel[*pParticleNo] , &h_pDSMC->ParticleEffTemp[*pParticleNo] ,
											     h_pDSMC->EffVibDOF[_Cell->LocalId][SpeciesNo] , h_pDomain , _pSpecies->VibTemp ) ;
				}


				h_pDSMC->ParticleLastCollide[*pParticleNo]	= -1 ;
				
			}
		}else if ( h_WallType[WallNo].Type == 5 ){
			
			double  AlphaN , SigmaT , EintCoef , AlphaT , R , Theta , Un , Ux , OM ,CTH ;
			
			WallTemp		= h_WallType[WallNo].Temp ;
			AlphaN      = h_WallType[WallNo].AlphaN ;
			SigmaT      = h_WallType[WallNo].SigmaT ;
			EintCoef        = h_WallType[WallNo].EintCoef ;
			
			MostProbableSpeed	= sqrt(2.*BOLTZ*WallTemp/_pSpecies->Mass) ;
			
			AlphaT=SigmaT*(2-SigmaT);
			
			
			R= sqrt(-AlphaN*log(Randn()));			
			Theta =6.283185308*Randn() ;
			Un = abs( NormVel/ MostProbableSpeed )*sqrt(1-AlphaN);			
			NormVel = -MostProbableSpeed * sqrt ( R*R+Un*Un+2*R*Un* cos(Theta)) ; 
			
			R= sqrt(-AlphaT*log(Randn()));
			Theta =6.283185308*Randn() ;
			Ux =  ParaVelX/ MostProbableSpeed  *sqrt(1-AlphaT);
			ParaVelX =  MostProbableSpeed *(Ux+R*cos(Theta)) ;
			
			R= sqrt(-AlphaT*log(Randn()));
			Theta =6.283185308*Randn() ;	
			ParaVelY  =  MostProbableSpeed*R*cos(Theta);


			// Calculate velocities in Cartesian coordinate system.
			InverseVelocityToCartesian3D( pXVel , pYVel , pZVel , NormVel , ParaVelX , ParaVelY , h_Node , 
							_Cell , FaceNo , h_pMapping ) ;
			
				(*pXVel)	+= h_WallType[WallNo].XVel ;
				(*pYVel)	+= h_WallType[WallNo].YVel ;
				(*pZVel)	+= h_WallType[WallNo].ZVel ;
				
			// Set rotational energy.
      	OM = sqrt ( h_pDSMC->ParticleRotation[*pParticleNo]*(1.-EintCoef)/ ( BOLTZ * WallTemp ));
	    	
				if ( _pSpecies->RotDOF == 2.){
					
        	R = sqrt(- EintCoef * log( Randn()));
        
        	CTH = cos (6.283185308*Randn());

				}else{
									
				}
				
				h_pDSMC->ParticleRotation[*pParticleNo] = BOLTZ * WallTemp*(R*R+OM*OM+2*R*OM*CTH);
				// Set vibrational energy.
//				if ( h_pDSMC->EffVibDOF[_Cell->LocalId][SpeciesNo]  > 0. ){
				if (_pSpecies->VibMode	> 0. ) {
					h_pDSMC->ParticleVibration[*pParticleNo]= VibrationalEnergy( WallTemp , &h_pDSMC->ParticleVibLevel[*pParticleNo] , &h_pDSMC->ParticleEffTemp[*pParticleNo] ,
											     h_pDSMC->EffVibDOF[_Cell->LocalId][SpeciesNo] , h_pDomain , _pSpecies->VibTemp ) ;
				}
			
				h_pDSMC->ParticleLastCollide[*pParticleNo]	= -1 ;
		}
	
		// Sample the properties of the reflective particle on surface.
		SamplingSurface( 2 , *pXVel , *pYVel , *pZVel , h_pDSMC->ParticleRotation[*pParticleNo] , h_pDSMC->ParticleVibration[*pParticleNo] , SurfaceNo , _pSpecies , h_pDSMC ) ;
		SamplingSurface( 4 , NormVel , NormVel , NormVel , h_pDSMC->ParticleRotation[*pParticleNo] , h_pDSMC->ParticleVibration[*pParticleNo] , SurfaceNo , _pSpecies , h_pDSMC ) ;
	}
	
	return	tracking ;
}

//==============================================================================================================

void SamplingSurface(	int		Method ,
			double		XVel ,
			double		YVel ,
			double		ZVel ,
			double		RotEnergy ,
			double		VibEnergy ,
			int		SurfaceNo ,
			DSMC_SPECIES	*_pSpecies ,
			DSMC_DSMC	*h_pDSMC ){
				
	int		SpeciesNo ;
	double		Mass ;

	SpeciesNo	= _pSpecies->Id ;
	Mass		= _pSpecies->Mass ;

	if ( Method == 1 ){
		h_pDSMC->SampleSurfaceParticleNum[SurfaceNo][SpeciesNo]		+= 1. ;
		h_pDSMC->SampleSurfaceInXMomentum[SurfaceNo][SpeciesNo]		+= Mass*XVel ;
		h_pDSMC->SampleSurfaceInYMomentum[SurfaceNo][SpeciesNo]		+= Mass*YVel ;
		h_pDSMC->SampleSurfaceInZMomentum[SurfaceNo][SpeciesNo]		+= Mass*ZVel ;
		h_pDSMC->SampleSurfaceInTransEng[SurfaceNo][SpeciesNo]		+= 0.5*Mass*(XVel*XVel + YVel*YVel + ZVel*ZVel) ;
		h_pDSMC->SampleSurfaceInRotEng[SurfaceNo][SpeciesNo]		+= RotEnergy ;
		h_pDSMC->SampleSurfaceInVibEng[SurfaceNo][SpeciesNo]		+= VibEnergy ;
		
	}else if ( Method == 2 ){
		h_pDSMC->SampleSurfaceReXMomentum[SurfaceNo][SpeciesNo]		-= Mass*XVel ;
		h_pDSMC->SampleSurfaceReYMomentum[SurfaceNo][SpeciesNo]		-= Mass*YVel ;
		h_pDSMC->SampleSurfaceReZMomentum[SurfaceNo][SpeciesNo]		-= Mass*ZVel ; 
		h_pDSMC->SampleSurfaceReTransEng[SurfaceNo][SpeciesNo]		-= 0.5*Mass*(XVel*XVel + YVel*YVel + ZVel*ZVel) ;
		h_pDSMC->SampleSurfaceReRotEng[SurfaceNo][SpeciesNo]		-= RotEnergy ;
		h_pDSMC->SampleSurfaceReVibEng[SurfaceNo][SpeciesNo]		-= VibEnergy ;
		
	}else if ( Method == 3 ){
		h_pDSMC->SampleSurfaceInNormMomentum[SurfaceNo][SpeciesNo]	+= Mass*XVel ;
		
	}else if ( Method == 4 ){
		h_pDSMC->SampleSurfaceReNormMomentum[SurfaceNo][SpeciesNo]	-= Mass*XVel ;
		
	}else if ( Method == 5 ){
		h_pDSMC->SampleSurfaceStickingParticleNum[SurfaceNo][SpeciesNo]		+= 1. ;
	}
}

//==============================================================================================================

void EnterNewPaticle2D(	DSMC_DOMAIN		*h_pDomain ,
			DSMC_NODE		*h_Node ,
			DSMC_INLET		*h_Inlet ,
			DSMC_SPECIES		*h_Species ,
			DSMC_DSMC		*h_pDSMC  ){

	int		InletFaceNum , SpeciesNum , CellNo , ParticleNo , MaxParticleNum , EnterParticleNumINT ;
	double		EnterParticleNum , InXVel , InYVel , InZVel , InTemp , MostProbableSpeed ;
	double		InNormVel , InParaVel , InNormVelRatio , InParaVelRatio ;
	double		NormVel , ParaVel , XVel , YVel , ZVel ;
	double		FS1 , FS2 , U , UN , QA , A , randn ;
	DSMC_NODE	Node0 , Node1 ;
	
	InletFaceNum	= h_pDomain->InletFaceNum ;
	SpeciesNum	= h_pDomain->SpeciesNum ;
	ParticleNo	= h_pDomain->ParticleNum ;
	MaxParticleNum	= h_pDomain->MaxParticleNum ;
				
	
	for ( int i=0 ; i<InletFaceNum ; i++ ){
		InXVel		= h_Inlet[i].XVel ;
		InYVel		= h_Inlet[i].YVel ;
		InZVel		= h_Inlet[i].ZVel ;
		InTemp		= h_Inlet[i].Temp ;
		CellNo		= h_Inlet[i].CellNo ;
		Node0		= h_Node[h_Inlet[i].Node[0]] ;
		Node1		= h_Node[h_Inlet[i].Node[1]] ;
		
		
		for ( int j=0 ; j<SpeciesNum ; j++ ){
			MostProbableSpeed	= sqrt(2.*BOLTZ*InTemp/h_Species[j].Mass) ;
			
			EnterParticleNum			= h_pDSMC->InletEnterNum[i][j] + h_pDSMC->InletRemainderEnterNum[i][j] ;
			EnterParticleNumINT			= (int)EnterParticleNum ;
			h_pDSMC->InletRemainderEnterNum[i][j]	= EnterParticleNum - EnterParticleNumINT ;
			
			
			InverseVelocityToNormPara2D( &InNormVel , &InParaVel , InXVel , InYVel , Node0 , Node1 ) ;
				
			InNormVelRatio	= -InNormVel/MostProbableSpeed ;
			InParaVelRatio	= InParaVel/MostProbableSpeed ;

			FS1		= InNormVelRatio + sqrtf(InNormVelRatio*InNormVelRatio + 2.) ;
			FS2		= 0.5 * ( 1. + InNormVelRatio*(2.*InNormVelRatio-FS1) ) ;
			// The above constants are required for the entering distn. of eqn (12.5).
			
			
			h_pDomain->EnterParticleNum	+= EnterParticleNumINT ;
			

			// Set enter particles
			for ( int k=0 ; k<EnterParticleNumINT ; k++ ){
				
				if ( fabs( InNormVelRatio*MostProbableSpeed ) > 1.E-6 ){
					// U is a potential normalised thermal velocity component.
					// UN is a potential inward velocity component.
					// The inward normalised vel. component has been selected (eqn (12.5)).
					QA	= 3. ;
					if ( InNormVelRatio < -3. ) QA = fabs(InNormVelRatio) + 1. ;

					do{
						do{
							U	= -QA + 2.*QA*Randn() ;
							UN	= U + InNormVelRatio ;
						}while ( UN < 0. ) ;
						
						A	= (2.*UN/FS1) * exp(FS2-U*U) ;
					}while ( A < Randn() ) ;
					
					RandVelocity( &ParaVel , &ZVel , MostProbableSpeed ) ;
					
					NormVel	= -(UN*MostProbableSpeed) ;
					ParaVel	+= InParaVelRatio*MostProbableSpeed ;
				}else{
					RandVelocity( &ParaVel , &ZVel , MostProbableSpeed ) ;
					
					NormVel	= -(sqrt(-log(Randn()))*MostProbableSpeed + InNormVel*MostProbableSpeed) ;
					ParaVel += InParaVelRatio*MostProbableSpeed ;
				}
				
				InverseVelocityToCartesian2D( &XVel , &YVel , NormVel , ParaVel , Node0 , Node1 ) ;
				
				// Set velocity of particle.
				h_pDSMC->ParticleXVel[ParticleNo]	= XVel ;
				h_pDSMC->ParticleYVel[ParticleNo]	= YVel ;
				h_pDSMC->ParticleZVel[ParticleNo]	= ZVel ;
				
				
				// Set rotational energy of particle.
				if ( h_Species[j].RotDOF > 0. )
					h_pDSMC->ParticleRotation[ParticleNo]	= RotationalEnergy( InTemp , h_Species[j].RotDOF ) ; 

				// Set vibrational energy of particle.
//				if ( h_pDSMC->EffVibDOF[CellNo][j]  > 0. ){
				if (h_Species[j].VibMode	> 0. ) {
					h_pDSMC->ParticleVibration[ParticleNo]= VibrationalEnergy( InTemp , &h_pDSMC->ParticleVibLevel[ParticleNo] , &h_pDSMC->ParticleEffTemp[ParticleNo] ,
												   h_pDSMC->EffVibDOF[CellNo][j] , h_pDomain , h_Species[j].VibTemp ) ;
				}
				
				
				// Set position of particle.
				randn	= Randn() ;
				h_pDSMC->ParticleXCoord[ParticleNo]	= Node0.XCoord + (Node1.XCoord - Node0.XCoord)*randn ;
				h_pDSMC->ParticleYCoord[ParticleNo]	= Node0.YCoord + (Node1.YCoord - Node0.YCoord)*randn ;
				h_pDSMC->ParticleZCoord[ParticleNo]	= 0. ;
				
				
				// Set cell no. and species no. of particle.
				h_pDSMC->ParticleCellNo[ParticleNo]	= CellNo ;
				h_pDSMC->ParticleSpeciesNo[ParticleNo]	= j ;
				h_pDSMC->ParticleLastCollide[ParticleNo]= -1 ;
				
				ParticleNo++ ;
				
				if ( ParticleNo >= h_pDomain->MaxParticleNum ){
					cout << "Particle array size has reached the limit. Maximum_Particle_Number is " << MaxParticleNum << '.' << endl ;
					exit(1) ;
				}
			}
		}	
	}
	
	// Update particle number.
	h_pDomain->ParticleNum	= ParticleNo ;			
}


void EnterNewPaticle3D(	DSMC_DOMAIN		*h_pDomain ,
			DSMC_NODE		*h_Node ,
			DSMC_CELL		*h_Cell ,
			DSMC_INLET		*h_Inlet ,
			DSMC_SPECIES		*h_Species ,
			DSMC_DSMC		*h_pDSMC ,
			CELLMAPPING		*h_pMapping ){

	int		MPISize , MPIMyID ;
	int		InletFaceNum , SpeciesNum , CellNo , FaceNo , ParticleNo , MaxParticleNum , EnterParticleNumINT ;
	double		EnterParticleNum , InXVel , InYVel , InZVel , InTemp , MostProbableSpeed ;
	double		InNormVel , InParaVelX , InParaVelY , InNormVelRatio , InParaVelXRatio , InParaVelYRatio ;
	double		NormVel , ParaVelX , ParaVelY , XVel , YVel , ZVel , XCoord , YCoord , ZCoord ;
	double		FS1 , FS2 , U , UN , QA , A , randn ;
	
	
	MPI_Comm_size( MPI_COMM_WORLD , &MPISize ) ;
	MPI_Comm_rank( MPI_COMM_WORLD , &MPIMyID ) ;
	
	
	InletFaceNum	= h_pDomain->InletFaceNum ;
	SpeciesNum	= h_pDomain->SpeciesNum ;
	ParticleNo	= h_pDomain->ParticleNum ;
	MaxParticleNum	= h_pDomain->MaxParticleNum ;
				
	
	for ( int i=0 ; i<InletFaceNum ; i++ ){
		InXVel		= h_Inlet[i].XVel ;
		InYVel		= h_Inlet[i].YVel ;
		InZVel		= h_Inlet[i].ZVel ;
		InTemp		= h_Inlet[i].Temp ;
		CellNo		= h_Inlet[i].CellNo ;
		FaceNo		= h_Inlet[i].FaceNo ;
		
		
		for ( int j=0 ; j<SpeciesNum ; j++ ){
			if ( h_Inlet[i].CosineLawCoef < 0. )
				MostProbableSpeed	= sqrt(2.*BOLTZ*InTemp/h_Species[j].Mass) ;
			else if ( h_Inlet[i].CosineLawCoef >= 0. )
				MostProbableSpeed	= sqrt(8.*BOLTZ*InTemp/(PI*h_Species[j].Mass)) ;
				
			
			EnterParticleNum			= h_pDSMC->InletEnterNum[i][j] + h_pDSMC->InletRemainderEnterNum[i][j] ;
			EnterParticleNumINT			= (int)EnterParticleNum ;
			h_pDSMC->InletRemainderEnterNum[i][j]	= EnterParticleNum - EnterParticleNumINT ;
			
			
			InverseVelocityToNormPara3D( &InNormVel , &InParaVelX , &InParaVelY , InXVel , InYVel , InZVel , h_Node , 
							&h_Cell[CellNo] , FaceNo , h_pMapping ) ;
							
			// Debug.
			//cout << "InXVel: " << InXVel << ", InYVel: " << InYVel << ", InZVel: " << InZVel << '\n' ;
			//cout << "InNormVel: " << InNormVel << ", InParaVelX: " << InParaVelX << ", InParaVelY: " << InParaVelY << '\n' ;
			//getchar() ;
				
			InNormVelRatio	= -InNormVel/MostProbableSpeed ;
			InParaVelXRatio	= InParaVelX/MostProbableSpeed ;
			InParaVelYRatio	= InParaVelY/MostProbableSpeed ;

			FS1		= InNormVelRatio + sqrtf(InNormVelRatio*InNormVelRatio + 2.) ;
			FS2		= 0.5 * ( 1. + InNormVelRatio*(2.*InNormVelRatio-FS1) ) ;
			// The above constants are required for the entering distn. of eqn (12.5).
			
			
			h_pDomain->EnterParticleNum	+= EnterParticleNumINT ;
			

			// Set enter particles
			for ( int k=0 ; k<EnterParticleNumINT ; k++ ){
				if ( h_Inlet[i].CosineLawCoef < 0. ){
				
					if ( fabs( InNormVelRatio*MostProbableSpeed ) > 1.E-6 ){
						// U is a potential normalised thermal velocity component.
						// UN is a potential inward velocity component.
						// The inward normalised vel. component has been selected (eqn (12.5)).
						QA	= 3. ;
						if ( InNormVelRatio < -3. ) QA = fabs(InNormVelRatio) + 1. ;

						do{
							do{
								U	= -QA + 2.*QA*Randn() ;
								UN	= U + InNormVelRatio ;
							}while ( UN < 0. ) ;
						
							A	= (2.*UN/FS1) * exp(FS2-U*U) ;
						}while ( A < Randn() ) ;
					
						RandVelocity( &ParaVelX , &ParaVelY , MostProbableSpeed ) ;
					
						// Debug.
						//cout << "MPS: " << MostProbableSpeed << '\n' ;
						//cout << "FA: " << h_Cell[CellNo].FaceFA[FaceNo] << ", FB: " << h_Cell[CellNo].FaceFB[FaceNo] << ", FC: " << h_Cell[CellNo].FaceFC[FaceNo] << ", FD: " << h_Cell[CellNo].FaceFD[FaceNo] << '\n' ;
						//cout << "InNormVel: " << InNormVel << ", InParaVelX: " << InParaVelX << ", InParaVelY: " << InParaVelY << '\n' ;
						//cout << "ParaVelX: " << ParaVelX << ", ParaVelY: " << ParaVelY << '\n' ;
					
						NormVel		= -(UN*MostProbableSpeed) ;
						ParaVelX	+= InParaVelXRatio*MostProbableSpeed ;
						ParaVelY	+= InParaVelYRatio*MostProbableSpeed ;
					}else{
						RandVelocity( &ParaVelX , &ParaVelY , MostProbableSpeed ) ;
					
						//NormVel		= -(sqrt(-log(Randn()))*MostProbableSpeed + InNormVel*MostProbableSpeed) ;
						NormVel		= -(sqrt(-log(Randn()))*MostProbableSpeed) ;
						ParaVelX	+= InParaVelXRatio*MostProbableSpeed ;
						ParaVelY	+= InParaVelYRatio*MostProbableSpeed ;
					}
					
				}else if ( h_Inlet[i].CosineLawCoef >= 0. ){

					CosineLawVelocity( &NormVel , &ParaVelX , &ParaVelY , h_Inlet[i].CosineLawCoef , MostProbableSpeed ) ;
					NormVel	= -1.*NormVel ;
					
					// Debug.
					//cout << "Temp: " << InTemp << ", Mass: " << h_Species[j].Mass << ", : " << (8.*BOLTZ*InTemp/(PI*h_Species[j].Mass)) << '\n' ;
					//cout << "MPIMyID: " << MPIMyID << ", V: " << MostProbableSpeed << ", N: " << h_Inlet[i].CosineLawCoef << ", VN: " << NormVel << ", VPX: " << ParaVelX << ", VPY: " << ParaVelY << '\n' ;
					
				}


				InverseVelocityToCartesian3D( &XVel , &YVel , &ZVel , NormVel , ParaVelX , ParaVelY , h_Node , 
								&h_Cell[CellNo] , FaceNo , h_pMapping ) ;
				
				
				// Set velocity of particle.
				h_pDSMC->ParticleXVel[ParticleNo]	= XVel ;
				h_pDSMC->ParticleYVel[ParticleNo]	= YVel ;
				h_pDSMC->ParticleZVel[ParticleNo]	= ZVel ;
				
				
				// Set rotational energy of particle.
				if ( h_Species[j].RotDOF > 0. )
					h_pDSMC->ParticleRotation[ParticleNo]	= RotationalEnergy( InTemp , h_Species[j].RotDOF ) ; 

				// Set vibrational energy of particle.
//				if ( h_pDSMC->EffVibDOF[CellNo][j]  > 0. ){
				if (h_Species[j].VibMode > 0. ) {
					h_pDSMC->ParticleVibration[ParticleNo]= VibrationalEnergy( InTemp , &h_pDSMC->ParticleVibLevel[ParticleNo] , &h_pDSMC->ParticleEffTemp[ParticleNo] ,
												   h_pDSMC->EffVibDOF[CellNo][j] , h_pDomain , h_Species[j].VibTemp ) ;
				}
				
				
				// Set position of particle.
				CreateParticlePositionEnter3D( &XCoord , &YCoord , &ZCoord , &h_Inlet[i] , h_Node , h_Cell , h_pMapping ) ;
				
				h_pDSMC->ParticleXCoord[ParticleNo]	= XCoord ;
				h_pDSMC->ParticleYCoord[ParticleNo]	= YCoord ;
				h_pDSMC->ParticleZCoord[ParticleNo]	= ZCoord ;
				
				
				// Set cell no. and species no. of particle.
				h_pDSMC->ParticleCellNo[ParticleNo]	= CellNo ;
				h_pDSMC->ParticleSpeciesNo[ParticleNo]	= j ;
				h_pDSMC->ParticleLastCollide[ParticleNo]= -1 ;
				
				ParticleNo++ ;
				
				if ( ParticleNo >= h_pDomain->MaxParticleNum ){
					cout << "Particle array size has reached the limit. Maximum_Particle_Number is " << MaxParticleNum << '.' << endl ;
					exit(1) ;
				}
			}
		}	
	}
	
	// Update particle number.
	h_pDomain->ParticleNum	= ParticleNo ;			
}


void EnterNewPaticleAxisymmetric(	DSMC_DOMAIN		*h_pDomain ,
					DSMC_NODE		*h_Node ,
					DSMC_INLET		*h_Inlet ,
					DSMC_SPECIES		*h_Species ,
					DSMC_DSMC		*h_pDSMC  ){

	int		InletFaceNum , SpeciesNum , CellNo , ParticleNo , MaxParticleNum , EnterParticleNumINT ;
	double		EnterParticleNum , InXVel , InYVel , InZVel , InTemp , MostProbableSpeed ;
	double		InNormVel , InParaVel , InNormVelRatio , InParaVelRatio ;
	double		NormVel , ParaVel , XVel , YVel , ZVel ;
	double		FS1 , FS2 , U , UN , QA , A , randn ;
	double		FaceFA , FaceFB , FaceFC ;
	DSMC_NODE	Node0 , Node1 ;
	
	InletFaceNum	= h_pDomain->InletFaceNum ;
	SpeciesNum	= h_pDomain->SpeciesNum ;
	ParticleNo	= h_pDomain->ParticleNum ;
	MaxParticleNum	= h_pDomain->MaxParticleNum ;
				
	
	for ( int i=0 ; i<InletFaceNum ; i++ ){
		InXVel		= h_Inlet[i].XVel ;
		InYVel		= h_Inlet[i].YVel ;
		InZVel		= h_Inlet[i].ZVel ;
		InTemp		= h_Inlet[i].Temp ;
		CellNo		= h_Inlet[i].CellNo ;
		Node0		= h_Node[h_Inlet[i].Node[0]] ;
		Node1		= h_Node[h_Inlet[i].Node[1]] ;
		

		FaceFA	= Node1.YCoord - Node0.YCoord ;
		FaceFB	= Node0.XCoord - Node1.XCoord ;
		FaceFC	= -1. * (FaceFA*Node0.XCoord + FaceFB*Node0.YCoord) ;
		
		
		for ( int j=0 ; j<SpeciesNum ; j++ ){
			MostProbableSpeed	= sqrt(2.*BOLTZ*InTemp/h_Species[j].Mass) ;
			
			EnterParticleNum			= h_pDSMC->InletEnterNum[i][j] + h_pDSMC->InletRemainderEnterNum[i][j] ;
			EnterParticleNumINT			= (int)EnterParticleNum ;
			h_pDSMC->InletRemainderEnterNum[i][j]	= EnterParticleNum - EnterParticleNumINT ;
			
			
			InverseVelocityToNormPara2D( &InNormVel , &InParaVel , InXVel , InYVel , Node0 , Node1 ) ;
				
			InNormVelRatio	= -InNormVel/MostProbableSpeed ;
			InParaVelRatio	= InParaVel/MostProbableSpeed ;

			FS1		= InNormVelRatio + sqrtf(InNormVelRatio*InNormVelRatio + 2.) ;
			FS2		= 0.5 * ( 1. + InNormVelRatio*(2.*InNormVelRatio-FS1) ) ;
			// The above constants are required for the entering distn. of eqn (12.5).
			
			
			h_pDomain->EnterParticleNum	+= EnterParticleNumINT ;
			

			// Set enter particles
			for ( int k=0 ; k<EnterParticleNumINT ; k++ ){
				
				if ( fabs( InNormVelRatio*MostProbableSpeed ) > 1.E-6 ){
					// U is a potential normalised thermal velocity component.
					// UN is a potential inward velocity component.
					// The inward normalised vel. component has been selected (eqn (12.5)).
					QA	= 3. ;
					if ( InNormVelRatio < -3. ) QA = fabs(InNormVelRatio) + 1. ;

					do{
						do{
							U	= -QA + 2.*QA*Randn() ;
							UN	= U + InNormVelRatio ;
						}while ( UN < 0. ) ;
						
						A	= (2.*UN/FS1) * exp(FS2-U*U) ;
					}while ( A < Randn() ) ;
					
					RandVelocity( &ParaVel , &ZVel , MostProbableSpeed ) ;
					
					NormVel	= -(UN*MostProbableSpeed) ;
					ParaVel	+= InParaVelRatio*MostProbableSpeed ;
				}else{
					RandVelocity( &ParaVel , &ZVel , MostProbableSpeed ) ;
					
					NormVel	= -(sqrt(-log(Randn()))*MostProbableSpeed + InNormVel*MostProbableSpeed) ;
					ParaVel += InParaVelRatio*MostProbableSpeed ;
				}
				
				InverseVelocityToCartesian2D( &XVel , &YVel , NormVel , ParaVel , Node0 , Node1 ) ;
				
				// Set velocity of particle.
				h_pDSMC->ParticleXVel[ParticleNo]	= XVel ;
				h_pDSMC->ParticleYVel[ParticleNo]	= YVel ;
				h_pDSMC->ParticleZVel[ParticleNo]	= ZVel ;
				
				
				// Set rotational energy of particle.
				if ( h_Species[j].RotDOF > 0. )
					h_pDSMC->ParticleRotation[ParticleNo]	= RotationalEnergy( InTemp , h_Species[j].RotDOF ) ; 

				// Set vibrational energy of particle.
//				if ( h_pDSMC->EffVibDOF[CellNo][j]  > 0.){
				if (h_Species[j].VibMode	> 0. ) {
					h_pDSMC->ParticleVibration[ParticleNo]= VibrationalEnergy( InTemp , &h_pDSMC->ParticleVibLevel[ParticleNo] , &h_pDSMC->ParticleEffTemp[ParticleNo] ,
												   h_pDSMC->EffVibDOF[CellNo][j] , h_pDomain , h_Species[j].VibTemp ) ;
				}
				
				
				// Set position of particle.
				randn	= Randn() ;
				h_pDSMC->ParticleYCoord[ParticleNo]	= sqrt(Node0.YCoord*Node0.YCoord + (Node1.YCoord*Node1.YCoord - Node0.YCoord*Node0.YCoord)*randn) ;
					
				if ( fabs(FaceFA) <= 1.e-15 || fabs(FaceFB) <= 1.e-15 ){
					h_pDSMC->ParticleXCoord[ParticleNo]	= Node0.XCoord + (Node1.XCoord - Node0.XCoord)*randn ;
				}else{
					h_pDSMC->ParticleXCoord[ParticleNo]	= -1.*(FaceFB*h_pDSMC->ParticleYCoord[ParticleNo] + FaceFC)/FaceFA ;
				}
				
				h_pDSMC->ParticleZCoord[ParticleNo]	= 0. ;
				
				
				// Set cell no. and species no. of particle.
				h_pDSMC->ParticleCellNo[ParticleNo]	= CellNo ;
				h_pDSMC->ParticleSpeciesNo[ParticleNo]	= j ;
				h_pDSMC->ParticleLastCollide[ParticleNo]= -1 ;
				
				ParticleNo++ ;
				
				if ( ParticleNo >= h_pDomain->MaxParticleNum ){
					cout << "Particle array size has reached the limit. Maximum_Particle_Number is " << MaxParticleNum << '.' << endl ;
					exit(1) ;
				}
			}
		}	
	}
	
	// Update particle number.
	h_pDomain->ParticleNum	= ParticleNo ;			
}

//==============================================================================================================

void RemoveParticle( int *pParticleNo , int *pParticleNum , DSMC_DSMC *h_pDSMC ){
	h_pDSMC->ParticleCellNo[(*pParticleNo)]		= h_pDSMC->ParticleCellNo[(*pParticleNum)-1] ;
	h_pDSMC->ParticleSpeciesNo[(*pParticleNo)]	= h_pDSMC->ParticleSpeciesNo[(*pParticleNum)-1] ;
	h_pDSMC->ParticleLastCollide[(*pParticleNo)]	= h_pDSMC->ParticleLastCollide[(*pParticleNum)-1] ;
	h_pDSMC->ParticleXCoord[(*pParticleNo)]		= h_pDSMC->ParticleXCoord[(*pParticleNum)-1] ;
	h_pDSMC->ParticleYCoord[(*pParticleNo)]		= h_pDSMC->ParticleYCoord[(*pParticleNum)-1] ;
	h_pDSMC->ParticleZCoord[(*pParticleNo)]		= h_pDSMC->ParticleZCoord[(*pParticleNum)-1] ;
	h_pDSMC->ParticleXVel[(*pParticleNo)]		= h_pDSMC->ParticleXVel[(*pParticleNum)-1] ;
	h_pDSMC->ParticleYVel[(*pParticleNo)]		= h_pDSMC->ParticleYVel[(*pParticleNum)-1] ;
	h_pDSMC->ParticleZVel[(*pParticleNo)]		= h_pDSMC->ParticleZVel[(*pParticleNum)-1] ;
	h_pDSMC->ParticleRotation[(*pParticleNo)]	= h_pDSMC->ParticleRotation[(*pParticleNum)-1] ;
	h_pDSMC->ParticleVibration[(*pParticleNo)]	= h_pDSMC->ParticleVibration[(*pParticleNum)-1] ;
	h_pDSMC->ParticleVibLevel[(*pParticleNo)]	= h_pDSMC->ParticleVibLevel[(*pParticleNum)-1] ;
	
	(*pParticleNo)-- ;
	(*pParticleNum)-- ;
}

//==============================================================================================================

void TransferParticleToBuffer(	DSMC_DSMC		*h_pDSMC , 
				int			*pParticleNo , 
				int			*pParticleNum , 
				DSMC_MPI_PARTICLE	*h_MPIParticleIn , 
				int			*pMPIParticleNum , 
				int			ProcessorNo ,
				int			GlobalCellNo ,
				int			BeforeCellNo ,
				double			Timestep , 
				double			BeforeTimestep ){
	
	int		No ;
	
	No		= (*pMPIParticleNum) ;
			
	h_MPIParticleIn[No].ProcessorNo		= ProcessorNo ;
	h_MPIParticleIn[No].GlobalCellNo	= GlobalCellNo ;
	h_MPIParticleIn[No].BeforeCellNo	= BeforeCellNo ;
	h_MPIParticleIn[No].SpeciesNo		= h_pDSMC->ParticleSpeciesNo[(*pParticleNo)] ;
	h_MPIParticleIn[No].LastCollide		= h_pDSMC->ParticleLastCollide[(*pParticleNo)] ;
	h_MPIParticleIn[No].VibLevel		= h_pDSMC->ParticleVibLevel[(*pParticleNo)] ;
	h_MPIParticleIn[No].XCoord		= h_pDSMC->ParticleXCoord[(*pParticleNo)] ;
	h_MPIParticleIn[No].YCoord		= h_pDSMC->ParticleYCoord[(*pParticleNo)] ;
	h_MPIParticleIn[No].ZCoord		= h_pDSMC->ParticleZCoord[(*pParticleNo)] ;
	h_MPIParticleIn[No].XVel		= h_pDSMC->ParticleXVel[(*pParticleNo)] ;
	h_MPIParticleIn[No].YVel		= h_pDSMC->ParticleYVel[(*pParticleNo)] ;
	h_MPIParticleIn[No].ZVel		= h_pDSMC->ParticleZVel[(*pParticleNo)] ;
	h_MPIParticleIn[No].Rotation		= h_pDSMC->ParticleRotation[(*pParticleNo)] ;
	h_MPIParticleIn[No].Vibration		= h_pDSMC->ParticleVibration[(*pParticleNo)] ;
	h_MPIParticleIn[No].EffTemp		= h_pDSMC->ParticleEffTemp[(*pParticleNo)] ;
	h_MPIParticleIn[No].Timestep		= Timestep ;
	h_MPIParticleIn[No].BeforeTimestep	= BeforeTimestep ;
	h_MPIParticleIn[No].ReactionIndex		= h_pDSMC->ReactionIndex[(*pParticleNo)] ;		
			
	(*pMPIParticleNum)++ ;
			
	RemoveParticle( pParticleNo , pParticleNum , h_pDSMC ) ;
}

//==============================================================================================================
//==============================================================================================================

int AddParticleFromOtherProcessor(	DSMC_DOMAIN		*h_pDomain ,
					DSMC_PROCESSOR		*h_pProcessor , 
					DSMC_DSMC		*h_pDSMC ,
					DSMC_CELL		*h_Cell ,
					DSMC_MPI_PARTICLE	*h_MPIParticleIn ,
					int			*pMPIParticleNum , 
					ofstream		&OutputDebug ){

	int		MPISize , MPIMyID ;
	int		ParticleNo , ParticleNum , LocalCellNo ;


	MPI_Comm_size( MPI_COMM_WORLD , &MPISize ) ;
	MPI_Comm_rank( MPI_COMM_WORLD , &MPIMyID ) ;


	ParticleNo		= h_pDomain->ParticleNum ;
	ParticleNum		= ParticleNo ;


	for ( int i=0 ; i<(*pMPIParticleNum) ; i++ ){
		LocalCellNo	= h_pProcessor->LocalCellNo[h_MPIParticleIn[i].GlobalCellNo] ;
		
		
		// Debug.
		/*if ( MPIMyID == 9 ){
			if ( LocalCellNo >= h_pDomain->CellNum ){
				OutputDebug << "Hello, G: " << i << ", " << h_MPIParticleIn[i].GlobalCellNo << ", L: " << LocalCellNo << ", PN: " << h_pProcessor->CellProcessorNo[h_MPIParticleIn[i].GlobalCellNo] << endl ;
			}
		}*/
		
			
		h_pDSMC->ParticleCellNo[ParticleNum]		= LocalCellNo ;
		h_pDSMC->ParticleSpeciesNo[ParticleNum]		= h_MPIParticleIn[i].SpeciesNo ;
		h_pDSMC->ParticleLastCollide[ParticleNum]	= h_MPIParticleIn[i].LastCollide ;
		h_pDSMC->ParticleBeforeCellNo[ParticleNum]	= h_MPIParticleIn[i].BeforeCellNo ;
		h_pDSMC->ParticleXCoord[ParticleNum]		= h_MPIParticleIn[i].XCoord ;
		h_pDSMC->ParticleYCoord[ParticleNum]		= h_MPIParticleIn[i].YCoord ;
		h_pDSMC->ParticleZCoord[ParticleNum]		= h_MPIParticleIn[i].ZCoord ;
		h_pDSMC->ParticleXVel[ParticleNum]		= h_MPIParticleIn[i].XVel ;
		h_pDSMC->ParticleYVel[ParticleNum]		= h_MPIParticleIn[i].YVel ;
		h_pDSMC->ParticleZVel[ParticleNum]		= h_MPIParticleIn[i].ZVel ;
		h_pDSMC->ParticleRotation[ParticleNum]		= h_MPIParticleIn[i].Rotation ;
		h_pDSMC->ParticleVibration[ParticleNum]		= h_MPIParticleIn[i].Vibration ;
		h_pDSMC->ParticleEffTemp[ParticleNum]		= h_MPIParticleIn[i].EffTemp ;
		h_pDSMC->ParticleTimestep[ParticleNum]		= h_MPIParticleIn[i].Timestep * 
									  (h_Cell[LocalCellNo].Timestep/h_MPIParticleIn[i].BeforeTimestep) ;
		h_pDSMC->ParticleVibLevel[ParticleNum]		= h_MPIParticleIn[i].VibLevel ;
		h_pDSMC->ReactionIndex[ParticleNum]		= h_MPIParticleIn[i].ReactionIndex ;
			
		ParticleNum++ ;
	}


	// Update the total particle number.
	h_pDomain->ParticleNum	= ParticleNum ;
	(*pMPIParticleNum)	= 0 ;

	
	return	ParticleNo ;	
}

//==============================================================================================================

void MPITransferParticle(	DSMC_DOMAIN		*h_pDomain ,
				DSMC_PROCESSOR		*h_pProcessor , 
				DSMC_MPI_PARTICLE	*h_MPIParticleOut , 
				DSMC_MPI_PARTICLE	*h_MPIParticleIn ,
				int			*MPIParticleNumOut ,
				int			*MPIParticleNumIn ,
				int			*pMPIParticleNum ,
				DSMC_MPI_DATATYPE	*pMPIDataType ,
				ofstream		&OutputDebug ){
	
	int		MPISize , MPIMyID ;
	int		ProcessorNo , No ;
	int		*SumNumOut , *SumNumIn ;
	MPI_Status	istat[8] ;
	
	
	MPI_Comm_size( MPI_COMM_WORLD , &MPISize ) ;
	MPI_Comm_rank( MPI_COMM_WORLD , &MPIMyID ) ;


	SumNumOut	= new int[MPISize] ;
	SumNumIn	= new int[MPISize] ;
	

	SetInitValue( MPIParticleNumOut , 0 , MPISize ) ;
	SetInitValue( MPIParticleNumIn , 0 , MPISize ) ;

	
	for ( int i=0 ; i<(*pMPIParticleNum) ; i++ ){
		ProcessorNo	= h_MPIParticleIn[i].ProcessorNo ;
		
		MPIParticleNumIn[ProcessorNo]++ ;
	}
	
	
	SumNumOut[0]	= 0 ;
	for ( int i=1 ; i<MPISize ; i++ ){
		SumNumOut[i]	= SumNumOut[i-1] + MPIParticleNumIn[i-1] ;
	}
	
	
	for ( int i=0 ; i<(*pMPIParticleNum) ; i++ ){
		ProcessorNo	= h_MPIParticleIn[i].ProcessorNo ;
		No		= SumNumOut[ProcessorNo] + MPIParticleNumOut[ProcessorNo] ;
		
		h_MPIParticleOut[No]	= h_MPIParticleIn[i] ;
		
		MPIParticleNumOut[ProcessorNo]++ ;
	}
	
	
	// Debug.
	/*for ( int i=0 ; i<MPISize ; i++ ){
		OutputDebug << i << ", " << SumNumOut[i] << ", " << MPIParticleNumOut[i] << ", " << MPIParticleNumIn[i] << endl ;
	}*/
	
	
	SetInitValue( MPIParticleNumIn , 0 , MPISize ) ;
	
	
	// Send the number of transfered particles.
	for ( int i=0 ; i<MPISize ; i++ ){
		if ( MPIMyID != i ){
			MPI_Send( &MPIParticleNumOut[i] , 1 , MPI_INT , i , 0 , MPI_COMM_WORLD ) ;
		}else{
			for ( int j=0 ; j<MPISize ; j++ ){
				if ( MPIMyID != j ){
					MPI_Recv( &MPIParticleNumIn[j] , 1 , MPI_INT , j , 0 , MPI_COMM_WORLD , istat ) ;
				}
			}
		}
	}
	
	
	SumNumIn[0]	= 0 ;
	for ( int i=1 ; i<MPISize ; i++ ){
		SumNumIn[i]	= SumNumIn[i-1] + MPIParticleNumIn[i-1] ;
	}
	
	
	// Send and recive Particle data to buffer of others processor.
	for ( int i=0 ; i<MPISize ; i++ ){
		if ( MPIMyID != i ){
			MPI_Send( &h_MPIParticleOut[SumNumOut[i]] , MPIParticleNumOut[i] , pMPIDataType->MPI_PARTICLE , i , 1 , MPI_COMM_WORLD ) ;
		}else{
			for ( int j=0 ; j<MPISize ; j++ ){
				if ( MPIMyID != j ){
					MPI_Recv( &h_MPIParticleIn[SumNumIn[j]] , MPIParticleNumIn[j] , pMPIDataType->MPI_PARTICLE , j , 1 , MPI_COMM_WORLD , istat ) ;
				}
			}
		}	
	}
	
	
	(*pMPIParticleNum)	= SumNumIn[MPISize-1] + MPIParticleNumIn[MPISize-1] ;
	

	// Debug.
	/*for ( int i=0 ; i<MPISize ; i++ ){
		OutputDebug << i << ", " << SumNumIn[i] << ", " << MPIParticleNumIn[i] << endl ;
	}*/
	/*for ( int j=0 ; j<(*pMPIParticleNum) ; j++ ){
		//if ( h_MPIParticleIn[j].ProcessorNo != MPIMyID ){
			OutputDebug << setw(10) << j << setw(10) << h_MPIParticleIn[j].ProcessorNo << setw(10) << h_MPIParticleIn[j].GlobalCellNo << setw(10) << h_MPIParticleIn[j].XVel << endl ;
		//}
	}
	MPI_Barrier( MPI_COMM_WORLD ) ;*/
	//exit(1) ;

	
	delete [] SumNumOut ;
	delete [] SumNumIn ;
}

//==============================================================================================================
//==============================================================================================================

void Index(	DSMC_DOMAIN		*h_pDomain ,
		DSMC_SPECIES		*h_Species ,
		DSMC_DSMC		*h_pDSMC , 
		ofstream		&OutputDebug ){

	int	MPISize , MPIMyID ;
	int	*IndexGroup1 , *IndexGroup2 , GroupNum , CellNum , ParticleNum , CountParticleNum ;
	int	SpeciesNo , GroupNo , CellNo , ParticleNo ;


	MPI_Comm_size( MPI_COMM_WORLD , &MPISize ) ;
	MPI_Comm_rank( MPI_COMM_WORLD , &MPIMyID ) ;


	GroupNum	= h_pDomain->SpeciesGroupNum ;
	CellNum		= h_pDomain->CellNum ;
	ParticleNum	= h_pDomain->ParticleNum ;
	CountParticleNum= 0 ;


	IndexGroup1	= new int[GroupNum] ;
	IndexGroup2	= new int[GroupNum] ;


	for ( int i=0 ; i<GroupNum ; i++ ){
		IndexGroup1[i]	= 0 ;
		IndexGroup2[i]	= 0 ;
			
		for ( int j=0 ; j<CellNum ; j++ ){
			h_pDSMC->IndexCell2[i][j]	= 0 ;
		}
	}
	
	for ( int i=0 ; i<ParticleNum ; i++ ){
		
		if ( h_pDSMC->ReactionIndex[i] == 3 ){
			
	h_pDSMC->ParticleCellNo[i]		= h_pDSMC->ParticleCellNo[(ParticleNum)-1] ;
	h_pDSMC->ParticleSpeciesNo[i]= h_pDSMC->ParticleSpeciesNo[(ParticleNum)-1] ;
	h_pDSMC->ParticleLastCollide[i]	= h_pDSMC->ParticleLastCollide[(ParticleNum)-1] ;
	h_pDSMC->ParticleXCoord[i]	= h_pDSMC->ParticleXCoord[(ParticleNum)-1] ;
	h_pDSMC->ParticleYCoord[i]	= h_pDSMC->ParticleYCoord[(ParticleNum)-1] ;
	h_pDSMC->ParticleZCoord[i]	= h_pDSMC->ParticleZCoord[(ParticleNum)-1] ;
	h_pDSMC->ParticleXVel[i]	= h_pDSMC->ParticleXVel[(ParticleNum)-1] ;
	h_pDSMC->ParticleYVel[i]	= h_pDSMC->ParticleYVel[(ParticleNum)-1] ;
	h_pDSMC->ParticleZVel[i]		= h_pDSMC->ParticleZVel[(ParticleNum)-1] ;
	h_pDSMC->ParticleRotation[i]= h_pDSMC->ParticleRotation[(ParticleNum)-1] ;
	h_pDSMC->ParticleVibration[i]	= h_pDSMC->ParticleVibration[(ParticleNum)-1] ;
	h_pDSMC->ParticleVibLevel[i]= h_pDSMC->ParticleVibLevel[(ParticleNum)-1] ;

	ParticleNum-- ;
			
			}
	
	}

	// Count particle number in each cell and each group.
	for ( int i=0 ; i<ParticleNum ; i++ ){
		
		SpeciesNo	= h_pDSMC->ParticleSpeciesNo[i] ;
		
		GroupNo		= h_Species[SpeciesNo].GroupNo ;
		CellNo		= h_pDSMC->ParticleCellNo[i] ;
		
		
		IndexGroup2[GroupNo]			+= 1 ;
		h_pDSMC->IndexCell2[GroupNo][CellNo]	+= 1 ;
	}
	
	
	// Scan IndexGroup1
	CountParticleNum = 0 ;
	for ( int i=0 ; i<GroupNum ; i++ ){
		IndexGroup1[i]	 = CountParticleNum ;
		CountParticleNum+= IndexGroup2[i] ;
	}


	// Scan IndexCell1
	for ( int i=0 ; i<GroupNum ; i++ ){
		CountParticleNum = IndexGroup1[i] ;
		
		for ( int j=0 ; j<CellNum ; j++ ){
			h_pDSMC->IndexCell1[i][j]	= CountParticleNum ;
			CountParticleNum		+= h_pDSMC->IndexCell2[i][j] ;
			
			h_pDSMC->IndexCell2[i][j]	= 0 ;
		}

	}
	
	
	// Count particle number on each cell and group, and create mapping between particles and cells.
	for ( int i= 0 ; i<ParticleNum ; i++ ){
		
			
		SpeciesNo	= h_pDSMC->ParticleSpeciesNo[i] ;
		
		GroupNo		= h_Species[SpeciesNo].GroupNo ;
		CellNo		= h_pDSMC->ParticleCellNo[i] ;
		
		
		ParticleNo = h_pDSMC->IndexCell1[GroupNo][CellNo] + h_pDSMC->IndexCell2[GroupNo][CellNo] ;
		h_pDSMC->IndexParticle[ParticleNo]	= i ;
		h_pDSMC->IndexCell2[GroupNo][CellNo]	+= 1 ;
		h_pDSMC->ReactionIndex[ParticleNo]	= -1. ;
		
	}
	
		h_pDomain->ParticleNum = ParticleNum;

	delete [] IndexGroup1 ;
	delete [] IndexGroup2 ;			
}

//==============================================================================================================
//==============================================================================================================

void Collision(	DSMC_DOMAIN		*h_pDomain ,
		DSMC_CELL		*h_Cell ,
		DSMC_SPECIES		*h_Species ,
		DSMC_DSMC		*h_pDSMC ,
		DSMC_RESULT		*h_pResult ,
		DSMC_CHEMICALTCE *h_ChemicalTCE ,
		ofstream		&OutputDebug ){
	
	int			CellNum , GroupNum , SpeciesNum , CollisionPairs ,SelectPair, SelectParticleNo[2] , SpeciesNo[2],SubCellNo[2] ;
	int			**IndexCellSub1 , **IndexCellSub2 , *IndexParticleSub , *IndexParticleSubCellNo , SubCellNum , TotalSubCellNum , ParticleNum  ,*GroupParticleNum1,*GroupParticleNum2,**ReactionBufferSub ,*ReactionBuffer;
	double			SumParticleNum , SamplingNum , RCollisionPairs , MaxCrossSectionSpeed , CrossSectionSpeed ;
	double			RelativeVel[3] , RelativeVelSq , RelativeSpeed ;
	double      MeanCrossSectionSpeed , AverageCrossSectionSpeed1, AverageCrossSectionSpeed2;

	DSMC_MULTISPECIES	MultiSpecies ;
	
//	DSMC_CHEMICALTCE	PChemicalTCE ;
	
	CellNum		= h_pDomain->CellNum ;
	GroupNum	= h_pDomain->SpeciesGroupNum ;
	SpeciesNum	= h_pDomain->SpeciesNum ;
	SamplingNum	= (double)h_pDomain->SamplingNum ;
	
	RelativeVel[0]	= 0. ;
	RelativeVel[1]	= 0. ;
	RelativeVel[2]	= 0. ;
	RelativeVelSq	= 0. ;
	RelativeSpeed	= 0. ;
	SubCellNum	= 0 ;
	SelectPair  = 0 ;
	ReactionBuffer = new int[GroupNum] ;
	
	for ( int i=0 ; i<CellNum ; i++ ){
		// Debug.
//		OutputDebug << "Cell: " << i << endl ;
		
		
		for ( int ii=0 ; ii<GroupNum ; ii++ )		
		ReactionBuffer[ii] =  h_pDSMC->IndexCell2[ii][i] ;

		if ( h_pDomain->SubcellModel > 0 ){
			ParticleNum	= 0 ;
			for ( int ii=0 ; ii<GroupNum ; ii++ )
				ParticleNum	+= h_pDSMC->IndexCell2[ii][i] ;
			
			
			SubCellNum		= GetSubCellNum( &h_Cell[i] , ParticleNum , h_pDomain->SubcellModel , h_pDomain->Dimension ) ;
			h_Cell[i].SubcellNum	= SubCellNum ;


			// Debug.			
			//OutputDebug << "SubN: " << SubCellNum << ", PN: " << ParticleNum << endl ;
			
			if ( SubCellNum > 0 ){
				// Allocate memory for subcell
				IndexCellSub1	= new int*[GroupNum] ;
				IndexCellSub2	= new int*[GroupNum] ;
				ReactionBufferSub	= new int*[GroupNum] ;			
				
				for ( int ii=0 ; ii<GroupNum ; ii++ ){
					if ( h_pDomain->Dimension == 2 || h_pDomain->Dimension == 4 ){
						TotalSubCellNum		= SubCellNum*SubCellNum ;
						IndexCellSub1[ii]	= new int[TotalSubCellNum] ;
						IndexCellSub2[ii]	= new int[TotalSubCellNum] ;
						ReactionBufferSub[ii]	 = new int[TotalSubCellNum] ;
					}else if ( h_pDomain->Dimension == 3 ){
						TotalSubCellNum		= SubCellNum*SubCellNum*SubCellNum ;
						IndexCellSub1[ii]	= new int[TotalSubCellNum] ;
						IndexCellSub2[ii]	= new int[TotalSubCellNum] ;
						ReactionBufferSub[ii]	 = new int[TotalSubCellNum] ;
					}
				}
				IndexParticleSub	= new int[ParticleNum] ;
				IndexParticleSubCellNo	= new int[ParticleNum] ;
				GroupParticleNum1	= new int[GroupNum] ;
				GroupParticleNum2	= new int[GroupNum] ;
				// End of allocated memory for subcell
				
				
				if ( h_pDomain->Dimension == 2 || h_pDomain->Dimension == 4 ){
					// Debug
					
					IndexSubCell2D( h_pDomain , h_pDSMC , &h_Cell[i] , i , IndexCellSub1 , IndexCellSub2 , IndexParticleSub , IndexParticleSubCellNo , SubCellNum , 
							ParticleNum ,GroupParticleNum1,GroupParticleNum2, OutputDebug ) ;	
				}else if ( h_pDomain->Dimension == 3 ){
					IndexSubCell3D( h_pDomain , h_pDSMC , &h_Cell[i] , i , IndexCellSub1 , IndexCellSub2 , IndexParticleSub , IndexParticleSubCellNo , SubCellNum , 
							ParticleNum ,GroupParticleNum1,GroupParticleNum2, OutputDebug ) ;
				}
				
					for ( int ii=0 ; ii<GroupNum ; ii++ ){	
						for ( int j=0 ; j<TotalSubCellNum ; j++){
							ReactionBufferSub[ii][j]=IndexCellSub2[ii][j] ;
						}
					}
			}
		}
		
		
		for ( int j=0 ; j<GroupNum ; j++ ){
			for ( int k=0 ; k<GroupNum ; k++ ){
				SumParticleNum	= 0. ;
				
				for ( int m=0 ; m<SpeciesNum ; m++ )
					if ( h_Species[m].GroupNo == k ) SumParticleNum += h_pDSMC->SampleParticleNum[m][i] ;


				if ( SumParticleNum > 1. )
					SumParticleNum	/= SamplingNum ;
				else
					SumParticleNum	= h_pDSMC->IndexCell2[k][i] ;
/*				
				if (j==k)					
					SumParticleNum = h_pDSMC->IndexCell2[k][i]-1;				
				else				
					SumParticleNum = h_pDSMC->IndexCell2[k][i];	
*/				
			
					
				// RCollisionPairs is the number of pairs to be selected, see eqn (11.5)
				RCollisionPairs	= 0.5*h_pDSMC->IndexCell2[j][i]*SumParticleNum*h_Cell[i].Weighting*h_pDSMC->MaxCrossSectionSpeed[i][j][k]*
							h_Cell[i].Timestep/h_Cell[i].Volume + h_pDSMC->RemainderCollisionPair[i][j][k] ;
							
//				RCollisionPairs	= 0.5*h_pDSMC->IndexCell2[j][i]*(h_pDSMC->IndexCell2[k][i])*h_Cell[i].Weighting*h_pDSMC->MaxCrossSectionSpeed[i][j][k]*
//							h_Cell[i].Timestep/h_Cell[i].Volume + h_pDSMC->RemainderCollisionPair[i][j][k] ;
				
				CollisionPairs	= (int)RCollisionPairs ;
				h_pDSMC->RemainderCollisionPair[i][j][k]	= RCollisionPairs - CollisionPairs ;
				
				
				// Debug.
				//OutputDebug << "In Coll: CollPairs: " << RCollisionPairs << endl ;
				//OutputDebug << h_pDSMC->IndexCell2[j][i] << ", SM: " << SumParticleNum << ", Wei: " << h_Cell[i].Weighting << ", Cros: " 
				//	    << h_pDSMC->MaxCrossSectionSpeed[i][j][k] << ", time: " << h_Cell[i].Timestep << ", V: " << h_Cell[i].Volume << endl ;
				
				if ( CollisionPairs > 0 ){
					if ( ((j != k) && (ReactionBuffer[j]< 1 || ReactionBuffer[k] < 1 )) || 
					     ((j == k) && (ReactionBuffer[j]< 2)) ){
					     	
						h_pDSMC->RemainderCollisionPair[i][j][k]	+= CollisionPairs ;
						
					}else{
												
						MaxCrossSectionSpeed	= h_pDSMC->MaxCrossSectionSpeed[i][j][k] ;
					
						for ( int m=0 ; m<CollisionPairs ; m++ ){
							
							// Select two paritcles to collide and return (cross section*speed).
							if ( ((j != k) && (ReactionBuffer[j]< 1 || ReactionBuffer[k] < 1 )) || 
					     ((j == k) && (ReactionBuffer[j]< 2)) ){
					     
					     SelectPair = m ;	
						
								break;	
							}
							
							 SelectPair = m+1 ;					
							
							if ( SubCellNum == 0 ){
								CrossSectionSpeed	= SelectParticle( h_Species , h_pDSMC , RelativeVel , &RelativeVelSq , &RelativeSpeed , SelectParticleNo , 
												h_pDSMC->IndexCell1[j][i] , h_pDSMC->IndexCell2[j][i] , h_pDSMC->IndexCell1[k][i] , 
												h_pDSMC->IndexCell2[k][i] , &MultiSpecies ) ;
							}else if ( SubCellNum > 0 ){
								CrossSectionSpeed	= SelectParticleSubCell( h_Species , h_pDSMC , RelativeVel , &RelativeVelSq , &RelativeSpeed , SelectParticleNo , 
												h_pDSMC->IndexCell1[j][i] , h_pDSMC->IndexCell2[j][i] , h_pDSMC->IndexCell1[k][i] , 
												h_pDSMC->IndexCell2[k][i] , &MultiSpecies , IndexCellSub1 , IndexCellSub2 , IndexParticleSub , 
												IndexParticleSubCellNo , j , k , SubCellNum , GroupParticleNum1 ,SubCellNo,ReactionBufferSub,TotalSubCellNum  ) ;
							}
							
							
								SpeciesNo[0]	= h_pDSMC->ParticleSpeciesNo[SelectParticleNo[0]] ;
								SpeciesNo[1]	= h_pDSMC->ParticleSpeciesNo[SelectParticleNo[1]] ;

							
							// if necessary, the maximum product in MaxCrossSectionSpeed is upgraded
						  	if ( CrossSectionSpeed > MaxCrossSectionSpeed ) MaxCrossSectionSpeed = CrossSectionSpeed ;
//								if ( CrossSectionSpeed > h_pDSMC->MaxCrossSectionSpeed[i][j][k] ) h_pDSMC->MaxCrossSectionSpeed[i][j][k] = CrossSectionSpeed ;

							
							// the collision is accepted with the probability of eqn (11.6)
								if ( Randn() < (CrossSectionSpeed/h_pDSMC->MaxCrossSectionSpeed[i][j][k]) ){
									h_Cell[i].CollNum++ ;
									h_Cell[i].CollDistance	+= sqrt( (h_pDSMC->ParticleXCoord[SelectParticleNo[0]]-h_pDSMC->ParticleXCoord[SelectParticleNo[1]]) * (h_pDSMC->ParticleXCoord[SelectParticleNo[0]]-h_pDSMC->ParticleXCoord[SelectParticleNo[1]]) + 
												 (h_pDSMC->ParticleYCoord[SelectParticleNo[0]]-h_pDSMC->ParticleYCoord[SelectParticleNo[1]]) * (h_pDSMC->ParticleYCoord[SelectParticleNo[0]]-h_pDSMC->ParticleYCoord[SelectParticleNo[1]]) + 
												 (h_pDSMC->ParticleZCoord[SelectParticleNo[0]]-h_pDSMC->ParticleZCoord[SelectParticleNo[1]]) * (h_pDSMC->ParticleZCoord[SelectParticleNo[0]]-h_pDSMC->ParticleZCoord[SelectParticleNo[1]]) ) ;
												 
												 
								// Validation for vibration energy with Bird's DSMC0V.FOR code.
									h_pDSMC->CollisionNum[SpeciesNo[0]][SpeciesNo[1]]++ ;
									h_pDSMC->CollisionNum[SpeciesNo[1]][SpeciesNo[0]]++ ;
								
									h_pDSMC->TotalCrossSectionSpeed[i][j][k] += CrossSectionSpeed ;
								
									MeanCrossSectionSpeed = h_pDSMC->TotalCrossSectionSpeed[i][j][k]/(h_pDSMC->CollisionNum[j][k]/2);
									
//							  AverageCrossSectionSpeed1 = 2/sqrt(PI)*MultiSpecies.CrossSection*pow((2-MultiSpecies.VisTempIndex)*MultiSpecies.RefTemp,MultiSpecies.VisTempIndex)*GammaFunction(2-MultiSpecies.VisTempIndex)*(2*BOLTZ/MultiSpecies.ReduceMass)*pow((h_pDomain->Temp),(0.5-MultiSpecies.VisTempIndex)) ;
//								AverageCrossSectionSpeed2 = 2/sqrt(PI)*MultiSpecies.CrossSection*sqrt(2*BOLTZ*MultiSpecies.RefTemp/MultiSpecies.ReduceMass)*pow((h_pDomain->Temp/MultiSpecies.RefTemp),(1-MultiSpecies.VisTempIndex)) ;
								
								
									if (h_pDomain->ChemicalNum != 0)	{
										
																		
									ChemicalTCE(	h_pDomain , h_Cell, h_Species , h_pDSMC , SelectParticleNo , SpeciesNo , RelativeVel , &RelativeVelSq , 
										&RelativeSpeed , &MultiSpecies ,h_ChemicalTCE ,  i , h_pResult->Temp[i] , h_pResult->TransTemp[i] ,ReactionBuffer, ReactionBufferSub,SubCellNo,SubCellNum ,OutputDebug) ;
																		
									}else{	
								// Rotational and vibrational energy are redistributed.
										if ( (h_Species[SpeciesNo[0]].RotDOF + h_Species[SpeciesNo[1]].RotDOF + h_pDSMC->EffVibDOF[i][SpeciesNo[0]] + h_pDSMC->EffVibDOF[i][SpeciesNo[1]]) > 0.01  ){
										Inelrv(	h_pDomain , h_Species , h_pDSMC , SelectParticleNo , SpeciesNo , RelativeVel , &RelativeVelSq , 
										&RelativeSpeed , &MultiSpecies, i , h_pResult->Temp[i] , h_pResult->TransTemp[i] ,OutputDebug ) ;
										}
								
									Elastic( h_pDSMC , RelativeVel , &RelativeSpeed , SelectParticleNo , &MultiSpecies , h_Species[SpeciesNo[0]].Mass , h_Species[SpeciesNo[1]].Mass ) ;
									}
								}	
																				
						}// End collision number	
						
						h_Cell[i].SelectNum	+= SelectPair ;						
						h_pDSMC->RemainderCollisionPair[i][j][k] += (CollisionPairs-SelectPair) ;

						h_pDSMC->MaxCrossSectionSpeed[i][j][k]	= MaxCrossSectionSpeed ;							
					}
					
				}// End if (CollisionPairs > 0)
				
			}// End group number -- k
			
			
		}// End group number -- j
		
		
		// Deallocate memory for subcell.
		if ( SubCellNum > 0 ){
			for ( int ii=0 ; ii<GroupNum ; ii++ ){
				delete [] IndexCellSub1[ii] ;
				delete [] IndexCellSub2[ii] ;
				delete [] ReactionBufferSub[ii];
			}
			delete [] IndexCellSub1 ;
			delete [] IndexCellSub2 ;
			delete [] IndexParticleSub ;
			delete [] IndexParticleSubCellNo ;
			delete [] GroupParticleNum1;
			delete [] GroupParticleNum2;
			delete [] ReactionBufferSub;
		}
		
		
	}// End cell number -- i
	
	delete []ReactionBuffer;
/*	
			OutputDebug << "CollisionNumforward: " << h_pDSMC->CollisionNum[1][0]/2 << endl ;
			OutputDebug << "CollisionNumbackword: " << h_pDSMC->CollisionNum[1][1] << endl ;
			OutputDebug << "SuccessReactionbackword1: " << h_ChemicalTCE[0].SuccessReaction[0] << endl ;
			OutputDebug << "SuccessReactionbackword2: " << h_ChemicalTCE[0].SuccessReaction[1] << endl ;
			OutputDebug << "Probabilitybackword1: " << h_ChemicalTCE[0].SuccessReaction[0]/(h_pDSMC->CollisionNum[1][1]/2) << endl ;
		  OutputDebug << "Probabilitybackword2: " << h_ChemicalTCE[0].SuccessReaction[1]/(h_pDSMC->CollisionNum[1][1]) << endl ;
			OutputDebug << "MeanCrossSectionSpeed: " << MeanCrossSectionSpeed << endl ;
			OutputDebug << "AverageCrossSectionSpeed1: " << AverageCrossSectionSpeed1 << endl ;
			OutputDebug << "AverageCrossSectionSpeed2: " << AverageCrossSectionSpeed2 << endl ;
*/
}
//==============================================================================================================
// Start Ming-Chung Lo
void ChemicalTCE(	DSMC_DOMAIN		*h_pDomain ,
		DSMC_CELL		*h_Cell ,
		DSMC_SPECIES		*h_Species ,
		DSMC_DSMC		*h_pDSMC ,
		int			*SelectParticleNo ,
		int			*SelectSpeciesNo ,
		double			*RelativeVel , 
		double			*pRelativeVelSq ,
		double			*pRelativeSpeed ,
		DSMC_MULTISPECIES	*pMultiSpecies ,
		DSMC_CHEMICALTCE *h_ChemicalTCE ,
		int			CellNo ,
		double			Temp , 
		double			TransTemp,
		int   *ReactionBuffer,
		int   **ReactionBufferSub,
		int   *SubCellNo,
		int    SubCellNum, 
		ofstream	&OutputDebug ){
										
//		DSMC_MULTISPECIES MultiSpecies ;

	
	int k,kk,kkk,ThirdBody, Na, ReactionClassNum, SuccessReactionNum, BadReactionNum,SymmetryFactor,GroupNum,ColClassNo,ChemicalNum,iclass,ReactionClassNo ;
	double InitTransEnergy, InitRotEnergy, InitVibEnergy, TotalEnergy,SumDegreeofFreedom,AverageRotDOF,AverageEffVibDOF,ReactionProbability   ; 
	double pi,*Probability, ProbabilityCoff1,  ProbabilityCoff2,ProbabilityCoff3,RadomProbability ,RecombineProbability;
	int Ispec ,Jspec ,Kspec, Select ;
	double ThirdBodyXVel,ThirdBodyYVel,ThirdBodyZVel;
	double ThirdBodyTransEnergy,ThirdBodyRotEnergy ,ThirdBodyVibEnergy,ThirdBodyTotalEnergy;
	double Buffer1,Buffer2,Buffer3;
	double NetTranDOF,NetRotDOF,NetEffVibDOF,ETA;
	
	RecombineProbability = -1 ;
	pi = PI ;
	SymmetryFactor = 1 ;
	GroupNum	= h_pDomain->SpeciesGroupNum ;
			
	SuccessReactionNum = 0 ;
	BadReactionNum = 0 ; 
				
//	CalculateMultiSpecies( &MultiSpecies , &h_Species[SelectSpeciesNo[0]] , &h_Species[SelectSpeciesNo[1]] ) ;
		
  // First particle must be a diatomic  
  if ( h_Species[SelectSpeciesNo[0]].RotDOF < 0.01 && h_Species[SelectSpeciesNo[1]].RotDOF >0.01 ){ 
  	
  	 k = SelectSpeciesNo[1] ;
  	 SelectSpeciesNo[1] =	SelectSpeciesNo[0] ;
  	 SelectSpeciesNo[0] = k ; 
  	 k = SelectParticleNo[1] ; 	
  	 SelectParticleNo[1] = SelectParticleNo[0] ;
  	 SelectParticleNo[0] = k ;
  	 k = SubCellNo[1] ; 	
  	 SubCellNo[1] = SubCellNo[0] ;
  	 SubCellNo[0] = k ;
  	 
	}
	
	ColClassNo =  h_Species[SelectSpeciesNo[0]].ReactionColClass[SelectSpeciesNo[1]] ;
	
	
	if	(	ColClassNo >= 0 ){
		
		ReactionClassNum = h_ChemicalTCE[ColClassNo].ReactionClassNum	;

			
		if ( ReactionClassNum != 0 ){
			
			if ( ReactionClassNum > 0.) Probability = new double [ReactionClassNum] ;

			InitTransEnergy	= 0.5*pMultiSpecies->ReduceMass*(*pRelativeVelSq) ;
			InitRotEnergy	= h_pDSMC->ParticleRotation[SelectParticleNo[0]] + h_pDSMC->ParticleRotation[SelectParticleNo[1]] ;
			InitVibEnergy	= h_pDSMC->ParticleVibration[SelectParticleNo[0]] + h_pDSMC->ParticleVibration[SelectParticleNo[1]] ;
			TotalEnergy	= InitTransEnergy + InitRotEnergy + InitVibEnergy ;
			
			NetTranDOF = 2*( 2.5 - pMultiSpecies->VisTempIndex );
			NetRotDOF =  h_Species[SelectSpeciesNo[0]].RotDOF + h_Species[SelectSpeciesNo[1]].RotDOF ;
			NetEffVibDOF = h_pDSMC->EffVibDOF[CellNo][SelectSpeciesNo[0]] + h_pDSMC->EffVibDOF[CellNo][SelectSpeciesNo[1]];
			SumDegreeofFreedom= NetTranDOF + NetRotDOF + NetEffVibDOF ;
			Buffer1= 2.*pMultiSpecies->CrossSection/sqrt(PI) ;
			Buffer2= pow((2.*BOLTZ*pMultiSpecies->RefTemp/pMultiSpecies->ReduceMass),(pMultiSpecies->VisTempIndex-0.5)) ;
			Buffer3= pow((2.*BOLTZ/pMultiSpecies->ReduceMass),(1-pMultiSpecies->VisTempIndex)) ;
			ETA=Buffer1*Buffer2*Buffer3 ;
			
			if ( SelectSpeciesNo[0] == SelectSpeciesNo[1] ) ETA=ETA/2. ;
				
				
			if ( ReactionClassNum > 0 ){
				
				ReactionProbability = 0. ;
			
				for ( int i=0 ; i< ReactionClassNum ; i++ ){ 
									
					Probability[i] = 0 ;
				
					if ( ( TotalEnergy - h_ChemicalTCE[ColClassNo].ArrheniusActiveEnergy[i]) > 0. ){
						
						ProbabilityCoff2 = ( SumDegreeofFreedom / 2.)-1.;
						ProbabilityCoff1 =  h_ChemicalTCE[ColClassNo].ArrheniusTempExp[i]-0.5+pMultiSpecies->VisTempIndex-0.5+ProbabilityCoff2 ;
						Buffer1 = ProbabilityCoff2+1. ;
						Buffer2 = ProbabilityCoff1+1. ;
						Buffer3 = ProbabilityCoff1-ProbabilityCoff2 ;
						ProbabilityCoff3 =  GammaFunction(Buffer1)/GammaFunction(Buffer2)*h_ChemicalTCE[ColClassNo].ArrheniusConstant[i]/ETA*pow(BOLTZ,-Buffer3) ;
						Probability[i] = ProbabilityCoff3*pow(( TotalEnergy-h_ChemicalTCE[ColClassNo].ArrheniusActiveEnergy[i]),ProbabilityCoff1)*pow(TotalEnergy,-ProbabilityCoff2) ;
						ReactionProbability = ReactionProbability +	Probability[i] ;
					
					}			
				}
				
							 			
				if ( Randn() < ReactionProbability ){
					
					SuccessReactionNum = SuccessReactionNum + 1 ;
					
					if ( ReactionProbability > 1)			BadReactionNum = BadReactionNum + 1 ;
					
					if ( ReactionClassNum == 1 ){
						
						ReactionClassNo = 0 ;
						
					}else{
				
						RadomProbability = Randn() * ReactionProbability ;
						ReactionClassNo = 0 ;
						ReactionProbability = Probability[0] ;
						
						while( RadomProbability > ReactionProbability ){		
							
						ReactionClassNo = ReactionClassNo +1 ;
						ReactionProbability = ReactionProbability + Probability [ReactionClassNo] ;
					
						}
					}
								
					h_ChemicalTCE[ColClassNo].SuccessReaction[ReactionClassNo] ++;
					
					if ( h_ChemicalTCE[ColClassNo].ReactionType[ReactionClassNo] == 1){
						 					  
						
						Dissociation(	h_pDomain , h_Species , h_pDSMC , SelectParticleNo , SelectSpeciesNo , RelativeVel , pRelativeVelSq , 
										pRelativeSpeed , pMultiSpecies ,h_ChemicalTCE , CellNo , Temp , TransTemp, ColClassNo ,ReactionClassNo,TotalEnergy,ReactionBuffer,ReactionBufferSub,SubCellNo,SubCellNum,OutputDebug ) ;
					
					}else{
						
						if ( h_ChemicalTCE[ColClassNo].ReactionType[ReactionClassNo] == 2){
																				
						Exchange(	h_pDomain , h_Species , h_pDSMC , SelectParticleNo , SelectSpeciesNo , RelativeVel , pRelativeVelSq , 
										pRelativeSpeed , pMultiSpecies ,h_ChemicalTCE , CellNo , Temp , TransTemp, ColClassNo ,ReactionClassNo,TotalEnergy,ReactionBuffer, ReactionBufferSub,SubCellNo,SubCellNum,OutputDebug ) ;
						}
					}	
						
				}else{ 
					
					if ( (h_Species[SelectSpeciesNo[0]].RotDOF + h_Species[SelectSpeciesNo[1]].RotDOF + h_pDSMC->EffVibDOF[CellNo][SelectSpeciesNo[0]] + h_pDSMC->EffVibDOF[CellNo][SelectSpeciesNo[1]]) > 0.01  ){
									Inelrv(	h_pDomain , h_Species , h_pDSMC , SelectParticleNo , SelectSpeciesNo , RelativeVel , pRelativeVelSq , 
										pRelativeSpeed , pMultiSpecies , CellNo , Temp, TransTemp,OutputDebug ) ;
					}
								
					Elastic( h_pDSMC , RelativeVel , pRelativeSpeed , SelectParticleNo , pMultiSpecies , h_Species[SelectSpeciesNo[0]].Mass , h_Species[SelectSpeciesNo[1]].Mass ) ;
				}
			
			}else{	
				
				kk = 0 ;
				kkk= 0 ;
				
				for ( int i=0 ; i < GroupNum ; i++ ){ 
								
					kk = kk + h_pDSMC->IndexCell2[i][CellNo] ;
					kkk = kkk + ReactionBuffer[i] ;
				}
				
				if ( kkk > 2 ){
									
					do{
//						Select = Randn()*(kk-0.001) + 1. ;
//						Na = 0 ;					   
//						 for  ( int j=0 ; j< GroupNum ; j++ ){							
//							Na = Na + h_pDSMC->IndexCell2[j][CellNo];							
//							 if ( Select <= Na ){								
//								Select = Randn()*(h_pDSMC->IndexCell2[j][CellNo]-1E-3) ;						
//					      break ;					      					      
//							}															  								
//						}								
              Select = Randn()*(kk-0.001) ;	
					   	Select =  Select + h_pDSMC->IndexCell1[0][CellNo] ;	
							ThirdBody = h_pDSMC->IndexParticle[Select] ;
							Kspec = h_pDSMC->ParticleSpeciesNo[ThirdBody] ;
					
					}while( ThirdBody == SelectParticleNo[0]||ThirdBody == SelectParticleNo[1]||h_pDSMC->ReactionIndex[ThirdBody] ==3.);
				
					
																				  																				
					ReactionClassNo = Kspec ;
					
					ThirdBodyXVel		= h_pDSMC->ParticleXVel[ThirdBody] ;
					ThirdBodyYVel		= h_pDSMC->ParticleYVel[ThirdBody] ;
					ThirdBodyZVel		= h_pDSMC->ParticleZVel[ThirdBody] ;
					
					ThirdBodyTransEnergy = 0.5 * h_Species[Kspec].Mass * (ThirdBodyXVel*ThirdBodyXVel+ThirdBodyYVel*ThirdBodyYVel+ThirdBodyZVel*ThirdBodyZVel) ;	
					ThirdBodyRotEnergy = h_pDSMC->ParticleRotation[ThirdBody] ;
					ThirdBodyVibEnergy = h_pDSMC->ParticleVibration[ThirdBody] ;		
					ThirdBodyTotalEnergy = ThirdBodyTransEnergy + ThirdBodyRotEnergy + ThirdBodyVibEnergy ;
					
					ProbabilityCoff1 = h_ChemicalTCE[ColClassNo].ArrheniusTempExp[ReactionClassNo]-0.5+pMultiSpecies->VisTempIndex-0.5 ;
					
					Buffer1 = 7.0/2.0+h_Species[Kspec].RotDOF/2.-pMultiSpecies->VisTempIndex+0.5 ;
					Buffer2 = Buffer1+ProbabilityCoff1 ;
					
					ProbabilityCoff3 = h_Cell[CellNo].Weighting*GammaFunction(Buffer1)/GammaFunction(Buffer2)*h_ChemicalTCE[ColClassNo].ArrheniusConstant[ReactionClassNo]/ETA*pow(BOLTZ,-ProbabilityCoff1);
		      RecombineProbability =  kk/h_Cell[CellNo].Volume*ProbabilityCoff3*pow((TotalEnergy+ThirdBodyTotalEnergy),ProbabilityCoff1);					
	
//				RecombineProbability =  h_pDSMC->IndexCell2[Kspec][CellNo]/h_Cell[CellNo].Volume*ProbabilityCoff3*pow((TotalEnergy+ThirdBodyTotalEnergy),ProbabilityCoff1);
//			  RecombineProbability =  ReactionBuffer[Kspec]/h_Cell[CellNo].Volume*ProbabilityCoff3*pow((TotalEnergy+ThirdBodyTotalEnergy),ProbabilityCoff1);						
//			  RecombineProbability =  kk/h_Cell[CellNo].Volume*ProbabilityCoff3*pow((TotalEnergy+ThirdBodyTotalEnergy),ProbabilityCoff1)*0.5*pow(1.0/(2-pMultiSpecies->VisTempIndex+0.5),pMultiSpecies->VisTempIndex-0.5);
			 
					if ( Randn() < RecombineProbability ){ 	
											
//						SuccessReactionNum = SuccessReactionNum + 1 ;		
				
						h_ChemicalTCE[ColClassNo].SuccessReaction[ReactionClassNo] ++;
					
							if ( ReactionProbability > 1)			BadReactionNum = BadReactionNum + 1 ;
								
								Recombination (	h_pDomain , h_Species , h_pDSMC , SelectParticleNo , SelectSpeciesNo , RelativeVel , pRelativeVelSq , 
										pRelativeSpeed , pMultiSpecies ,h_ChemicalTCE , CellNo , Temp , TransTemp, ColClassNo ,ReactionClassNo,TotalEnergy,ThirdBody,&SelectParticleNo[1],&h_pDomain->ParticleNum ,ReactionBuffer,ReactionBufferSub,SubCellNo,SubCellNum,OutputDebug) ;
					
					}else{
						
						if ( (h_Species[SelectSpeciesNo[0]].RotDOF + h_Species[SelectSpeciesNo[1]].RotDOF + h_pDSMC->EffVibDOF[CellNo][SelectSpeciesNo[0]] + h_pDSMC->EffVibDOF[CellNo][SelectSpeciesNo[1]]) > 0.01  ){
									Inelrv(	h_pDomain , h_Species , h_pDSMC , SelectParticleNo , SelectSpeciesNo , RelativeVel , pRelativeVelSq , 
										pRelativeSpeed , pMultiSpecies , CellNo , Temp , TransTemp,OutputDebug ) ;
						}
								
						Elastic( h_pDSMC , RelativeVel , pRelativeSpeed , SelectParticleNo , pMultiSpecies , h_Species[SelectSpeciesNo[0]].Mass , h_Species[SelectSpeciesNo[1]].Mass ) ;	
																				
					}
					
				}else{
					
					if ( (h_Species[SelectSpeciesNo[0]].RotDOF + h_Species[SelectSpeciesNo[1]].RotDOF + h_pDSMC->EffVibDOF[CellNo][SelectSpeciesNo[0]] + h_pDSMC->EffVibDOF[CellNo][SelectSpeciesNo[1]]) > 0.01  ){
									Inelrv(	h_pDomain , h_Species , h_pDSMC , SelectParticleNo , SelectSpeciesNo , RelativeVel , pRelativeVelSq , 
										pRelativeSpeed , pMultiSpecies , CellNo , Temp , TransTemp,OutputDebug ) ;
					}
								
					Elastic( h_pDSMC , RelativeVel , pRelativeSpeed , SelectParticleNo , pMultiSpecies , h_Species[SelectSpeciesNo[0]].Mass , h_Species[SelectSpeciesNo[1]].Mass ) ;	
															
				}
			}	
			
			if ( ReactionClassNum > 0.) delete [] Probability ;	
		}
		
	}else{
					
		if ( (h_Species[SelectSpeciesNo[0]].RotDOF + h_Species[SelectSpeciesNo[1]].RotDOF + h_pDSMC->EffVibDOF[CellNo][SelectSpeciesNo[0]] + h_pDSMC->EffVibDOF[CellNo][SelectSpeciesNo[1]]) > 0.01  ){
			Inelrv(	h_pDomain , h_Species , h_pDSMC , SelectParticleNo , SelectSpeciesNo , RelativeVel , pRelativeVelSq , 
				pRelativeSpeed , pMultiSpecies , CellNo , Temp , TransTemp,OutputDebug ) ;
		}
								
		Elastic( h_pDSMC , RelativeVel , pRelativeSpeed , SelectParticleNo , pMultiSpecies , h_Species[SelectSpeciesNo[0]].Mass , h_Species[SelectSpeciesNo[1]].Mass ) ;	
														
	}
	
	
}

//==============================================================================================================

void Dissociation(	DSMC_DOMAIN		*h_pDomain ,
		DSMC_SPECIES		*h_Species ,
		DSMC_DSMC		*h_pDSMC ,
		int			*SelectParticleNo ,
		int			*SelectSpeciesNo ,
		double			*RelativeVel , 
		double			*pRelativeVelSq ,
		double			*pRelativeSpeed ,
		DSMC_MULTISPECIES	*pMultiSpecies ,
		DSMC_CHEMICALTCE	*h_ChemicalTCE ,
		int			CellNo ,
		double			Temp , 
		double			TransTemp,
		int 		ColClassNo,
		int 		ReactionClassNo,
		double  TotalEnergy,
		int *ReactionBuffer,
		int   **ReactionBufferSub,
		int   *SubCellNo,
		int    SubCellNum, 
		ofstream	&OutputDebug ){
		
		DSMC_MULTISPECIES	MultiSpecies ;
			
		int ParticleNum,VibLevel,MaxLevel,ispec,jspec,kspec,k ;
		double XIA,XIB,ProbabilityRatio,DividEnergy,ERM,BetweenRelativeEnergy,RotationEnergy,RelativeEnergy ;
		double B,Angle,RelativeVelPost[3],MassRatio[2],CenterMassVel[3],AzimuthAngle,RelativeSpeed;
				
		ispec = h_ChemicalTCE[ColClassNo].Post1Species[ReactionClassNo] ;
		jspec = h_ChemicalTCE[ColClassNo].Post2Species[ReactionClassNo] ;
		kspec = h_ChemicalTCE[ColClassNo].Post3Species[ReactionClassNo] ;
		
		
		if ( SelectSpeciesNo[1] != jspec ){ 
  	
  	 k = SelectSpeciesNo[1] ;
  	 SelectSpeciesNo[1] =	SelectSpeciesNo[0] ;
  	 SelectSpeciesNo[0] = k ; 
  	 
  	 k = SelectParticleNo[1] ; 	
  	 SelectParticleNo[1] = SelectParticleNo[0] ;
  	 SelectParticleNo[0] = k ;
  	 
  	 k = SubCellNo[1] ; 	
  	 SubCellNo[1] = SubCellNo[0] ;
  	 SubCellNo[0] = k ;
  	   	 
		}
				
		ReactionBuffer[h_Species[SelectSpeciesNo[0]].GroupNo] = ReactionBuffer[h_Species[SelectSpeciesNo[0]].GroupNo]-1 ;
		h_pDSMC->ReactionIndex[SelectParticleNo[0]] = 1 ;
		
		if (SubCellNum >0)
		ReactionBufferSub [h_Species[SelectSpeciesNo[0]].GroupNo][SubCellNo[0]] = ReactionBufferSub [h_Species[SelectSpeciesNo[0]].GroupNo][SubCellNo[0]]-1;
		
		ParticleNum	= h_pDomain->ParticleNum ;
		
		TotalEnergy	=  TotalEnergy + h_ChemicalTCE[ColClassNo].HeatofReaction[ReactionClassNo] ;
		
		
		XIA=0.5*h_Species[SelectSpeciesNo[1]].RotDOF ;
		XIB=( 2.5 - pMultiSpecies->VisTempIndex ) ;
    
//    if ( h_pDSMC->EffVibDOF[CellNo][SelectSpeciesNo[1]] > 0){
//    if ( h_pDSMC->EffVibDOF[CellNo][SelectSpeciesNo[1]] > 0.01){ 	
    if ( h_Species[SelectSpeciesNo[1]].VibMode > 0. ) {
    			
    	VibLevel = 0 ;
			MaxLevel = TotalEnergy / ( BOLTZ * h_Species[SelectSpeciesNo[1]].VibTemp );
			
			do {
				
				h_pDSMC->ParticleVibLevel[SelectParticleNo[1]] = Randn() * ( MaxLevel+0.99999 ) ;							
				h_pDSMC->ParticleVibration[SelectParticleNo[1]] = h_pDSMC->ParticleVibLevel[SelectParticleNo[1]] * (BOLTZ*h_Species[SelectSpeciesNo[1]].VibTemp);			
				ProbabilityRatio	= pow(( 1.0 - h_pDSMC->ParticleVibration[SelectParticleNo[1]]/TotalEnergy	),( 2*XIB+XIA-1 ));
			
			}while ( ProbabilityRatio < Randn() ) ;
			
			TotalEnergy	= TotalEnergy - h_pDSMC->ParticleVibration[SelectParticleNo[1]] ;						
	
		}else {
			
			h_pDSMC->ParticleVibLevel[SelectParticleNo[1]] =  0 ;
			h_pDSMC->ParticleVibration[SelectParticleNo[1]] = 0 ;
			
		}
			
		ERM	= LarsenBorgnakkeEnergyRatio( XIB-1 , XIA+XIB-1 ) ;	
		
		BetweenRelativeEnergy = TotalEnergy * ERM ;
			
		RotationEnergy = TotalEnergy - BetweenRelativeEnergy ;
			
		if (h_Species[SelectSpeciesNo[1]].RotDOF > 0.01 ){
			
			ERM	= LarsenBorgnakkeEnergyRatio( XIB-1 , XIA-1 ) ;				
			RelativeEnergy = RotationEnergy * ERM ;			
			h_pDSMC->ParticleRotation[SelectParticleNo[1]] = RotationEnergy - RelativeEnergy ;
			
		}else{
			
			RelativeEnergy = RotationEnergy ;
			h_pDSMC->ParticleRotation[SelectParticleNo[1]] = 0 ;
			
		}
			 
		MassRatio[0]	= pMultiSpecies->ReduceMass/h_Species[SelectSpeciesNo[1]].Mass ;
		MassRatio[1]	= pMultiSpecies->ReduceMass/h_Species[SelectSpeciesNo[0]].Mass ;	
			
		CenterMassVel[0]	= MassRatio[0]*(h_pDSMC->ParticleXVel[SelectParticleNo[0]]) + MassRatio[1]*(h_pDSMC->ParticleXVel[SelectParticleNo[1]]) ;
		CenterMassVel[1]	= MassRatio[0]*(h_pDSMC->ParticleYVel[SelectParticleNo[0]]) + MassRatio[1]*(h_pDSMC->ParticleYVel[SelectParticleNo[1]]) ;
		CenterMassVel[2]	= MassRatio[0]*(h_pDSMC->ParticleZVel[SelectParticleNo[0]]) + MassRatio[1]*(h_pDSMC->ParticleZVel[SelectParticleNo[1]]) ;
		
		RelativeSpeed = sqrt(2.*BetweenRelativeEnergy/pMultiSpecies->ReduceMass) ;

		B	= 2.*Randn() - 1. ;		
		Angle	= sqrt(1. - (B*B)) ;
			
		RelativeVelPost[0]	= B * RelativeSpeed ;

		// AzimuthAngle is a random azimuth angle
		AzimuthAngle		= 2.*PI*Randn() ;

		RelativeVelPost[1]	= Angle * cos(AzimuthAngle) * RelativeSpeed ;
		RelativeVelPost[2]	= Angle * sin(AzimuthAngle) * RelativeSpeed ;	
		
		h_pDSMC->ParticleXVel[SelectParticleNo[1]]	= CenterMassVel[0] - RelativeVelPost[0]*MassRatio[0] ;
		h_pDSMC->ParticleYVel[SelectParticleNo[1]]	= CenterMassVel[1] - RelativeVelPost[1]*MassRatio[0] ;
		h_pDSMC->ParticleZVel[SelectParticleNo[1]]	= CenterMassVel[2] - RelativeVelPost[2]*MassRatio[0] ; 

		CenterMassVel[0]	= CenterMassVel[0] + RelativeVelPost[0]*MassRatio[1] ;
		CenterMassVel[1]	= CenterMassVel[1] + RelativeVelPost[1]*MassRatio[1] ;
		CenterMassVel[2]	= CenterMassVel[2] + RelativeVelPost[2]*MassRatio[1] ;
				
		CalculateMultiSpecies( &MultiSpecies , &h_Species[ispec] , &h_Species[kspec] ) ;

		RelativeSpeed = sqrt(2.*RelativeEnergy/MultiSpecies.ReduceMass) ;
		
		B	= 2.*Randn() - 1. ;		
		Angle	= sqrt(1. - (B*B)) ;
		
		RelativeVelPost[0]	= B * RelativeSpeed ;

		// AzimuthAngle is a random azimuth angle
		AzimuthAngle		= 2.*PI*Randn() ;
		
		RelativeVelPost[1]	= Angle * cos(AzimuthAngle) * RelativeSpeed ;
		RelativeVelPost[2]	= Angle * sin(AzimuthAngle) * RelativeSpeed ;	
		
		MassRatio[0]	= MultiSpecies.ReduceMass/h_Species[kspec].Mass ;
		MassRatio[1]	= MultiSpecies.ReduceMass/h_Species[ispec].Mass ;	
			
		h_pDSMC->ParticleXVel[SelectParticleNo[0]]	= CenterMassVel[0] + RelativeVelPost[0]*MassRatio[1] ;
		h_pDSMC->ParticleYVel[SelectParticleNo[0]]	= CenterMassVel[1] + RelativeVelPost[1]*MassRatio[1] ;
		h_pDSMC->ParticleZVel[SelectParticleNo[0]]	= CenterMassVel[2] + RelativeVelPost[2]*MassRatio[1] ; 
		
		h_pDSMC->ParticleSpeciesNo[SelectParticleNo[0]] = ispec ;
		h_pDSMC->ParticleVibLevel [SelectParticleNo[0]] = 0 ;
		h_pDSMC->ParticleRotation [SelectParticleNo[0]] = 0 ;
		h_pDSMC->ParticleVibration[SelectParticleNo[0]] = 0 ;
				
		h_pDSMC->ParticleXVel[ParticleNum]	= CenterMassVel[0] - RelativeVelPost[0]*MassRatio[0] ;
		h_pDSMC->ParticleYVel[ParticleNum]	= CenterMassVel[1] - RelativeVelPost[1]*MassRatio[0] ;
		h_pDSMC->ParticleZVel[ParticleNum]	= CenterMassVel[2] - RelativeVelPost[2]*MassRatio[0] ; 
		
		h_pDSMC->ParticleSpeciesNo[ParticleNum] = kspec ;
		h_pDSMC->ParticleVibLevel [ParticleNum] = 0 ;
		h_pDSMC->ParticleRotation [ParticleNum] = 0 ;
		h_pDSMC->ParticleVibration[ParticleNum] = 0 ;
		
		h_pDSMC->ParticleXCoord[ParticleNum] = 	h_pDSMC->ParticleXCoord[SelectParticleNo[0]];
		h_pDSMC->ParticleYCoord[ParticleNum] =  h_pDSMC->ParticleYCoord[SelectParticleNo[0]];
		h_pDSMC->ParticleZCoord[ParticleNum] =	h_pDSMC->ParticleZCoord[SelectParticleNo[0]];
		h_pDSMC->ParticleCellNo[ParticleNum] = CellNo ;
				
		ParticleNum ++ ;
		
		h_pDomain->ParticleNum		= ParticleNum ;	
	
		
}
		
		
// End Ming-Chung Lo
//==============================================================================================================

void Exchange(	DSMC_DOMAIN		*h_pDomain ,
		DSMC_SPECIES		*h_Species ,
		DSMC_DSMC		*h_pDSMC ,
		int			*SelectParticleNo ,
		int			*SelectSpeciesNo ,
		double			*RelativeVel , 
		double			*pRelativeVelSq ,
		double			*pRelativeSpeed ,
		DSMC_MULTISPECIES	*pMultiSpecies ,
		DSMC_CHEMICALTCE	*h_ChemicalTCE ,
		int			CellNo ,
		double			Temp , 
		double			TransTemp,
		int 		ColClassNo,
		int 		ReactionClassNo,
		double  TotalEnergy,
		int *ReactionBuffer,
		int   **ReactionBufferSub,
		int   *SubCellNo,
		int    SubCellNum, 
		ofstream	&OutputDebug ){
			
		int VibLevel,MaxLevel,k ;
		double XIA,XIB,ProbabilityRatio,ERM,BetweenRelativeEnergy ;
		double B,Angle,RelativeVelPost[3],MassRatio[2],CenterMassVel[3],AzimuthAngle,RelativeSpeed;
		
		DSMC_MULTISPECIES	MultiSpecies ;
		

		 ReactionBuffer[h_Species[SelectSpeciesNo[0]].GroupNo] = ReactionBuffer[h_Species[SelectSpeciesNo[0]].GroupNo]-1 ;
		 ReactionBuffer[h_Species[SelectSpeciesNo[1]].GroupNo] = ReactionBuffer[h_Species[SelectSpeciesNo[1]].GroupNo]-1 ;
		  
		 h_pDSMC->ReactionIndex[SelectParticleNo[0]] = 2 ;
		 h_pDSMC->ReactionIndex[SelectParticleNo[1]] = 2 ;
		 
		if (SubCellNum >0) { 
		 ReactionBufferSub [h_Species[SelectSpeciesNo[0]].GroupNo][SubCellNo[0]] = ReactionBufferSub [h_Species[SelectSpeciesNo[0]].GroupNo][SubCellNo[0]]-1;
		 ReactionBufferSub [h_Species[SelectSpeciesNo[1]].GroupNo][SubCellNo[1]] = ReactionBufferSub [h_Species[SelectSpeciesNo[1]].GroupNo][SubCellNo[1]]-1;
		}


		TotalEnergy	= TotalEnergy + h_ChemicalTCE[ColClassNo].HeatofReaction[ReactionClassNo] ;
		
				 
		MassRatio[0]	= pMultiSpecies->ReduceMass/h_Species[SelectSpeciesNo[1]].Mass ;
		MassRatio[1]	= pMultiSpecies->ReduceMass/h_Species[SelectSpeciesNo[0]].Mass ;	
			
		CenterMassVel[0]	= MassRatio[0]*(h_pDSMC->ParticleXVel[SelectParticleNo[0]]) + MassRatio[1]*(h_pDSMC->ParticleXVel[SelectParticleNo[1]]) ;
		CenterMassVel[1]	= MassRatio[0]*(h_pDSMC->ParticleYVel[SelectParticleNo[0]]) + MassRatio[1]*(h_pDSMC->ParticleYVel[SelectParticleNo[1]]) ;
		CenterMassVel[2]	= MassRatio[0]*(h_pDSMC->ParticleZVel[SelectParticleNo[0]]) + MassRatio[1]*(h_pDSMC->ParticleZVel[SelectParticleNo[1]]) ;
		
	  SelectSpeciesNo[0] = h_ChemicalTCE[ColClassNo].Post1Species[ReactionClassNo] ;
	  SelectSpeciesNo[1] = h_ChemicalTCE[ColClassNo].Post2Species[ReactionClassNo] ;
	  
	  if ( h_Species[SelectSpeciesNo[0]].RotDOF > h_Species[SelectSpeciesNo[1]].RotDOF ){ 
  	
  		k = SelectSpeciesNo[1] ;
  	 	SelectSpeciesNo[1] =	SelectSpeciesNo[0] ;
  	 	SelectSpeciesNo[0] = k ; 
	 
		}
		
	  
	  h_pDSMC->ParticleSpeciesNo[SelectParticleNo[0]] =  SelectSpeciesNo[0] ;
		h_pDSMC->ParticleSpeciesNo[SelectParticleNo[1]] =  SelectSpeciesNo[1] ;
		
		CalculateMultiSpecies( &MultiSpecies , &h_Species[SelectSpeciesNo[0]] , &h_Species[SelectSpeciesNo[1]] ) ;
				
		XIA=0.5*h_Species[SelectSpeciesNo[1]].RotDOF ;
		XIB=( 2.5 - MultiSpecies.VisTempIndex ) ;
    
		VibLevel = 0 ;
		MaxLevel = TotalEnergy / ( BOLTZ * h_Species[SelectSpeciesNo[1]].VibTemp );
			
		do {
				
				h_pDSMC->ParticleVibLevel[SelectParticleNo[1]] = Randn() * ( MaxLevel+0.99999 ) ;							
				h_pDSMC->ParticleVibration[SelectParticleNo[1]] = h_pDSMC->ParticleVibLevel[SelectParticleNo[1]] * (BOLTZ*h_Species[SelectSpeciesNo[1]].VibTemp);			
				ProbabilityRatio	= pow(( 1.0 - h_pDSMC->ParticleVibration[SelectParticleNo[1]]/TotalEnergy	),( XIB+XIA-1 ));
			
		}while ( ProbabilityRatio < Randn() ) ;
			
		TotalEnergy	= TotalEnergy - h_pDSMC->ParticleVibration[SelectParticleNo[1]] ;			
		
		ERM	= LarsenBorgnakkeEnergyRatio( XIB-1 , XIA-1 ) ;	
		
		BetweenRelativeEnergy = TotalEnergy * ERM ;
			
		h_pDSMC->ParticleRotation[SelectParticleNo[1]] = TotalEnergy - BetweenRelativeEnergy ;
		
		RelativeSpeed = sqrt(2.*BetweenRelativeEnergy/MultiSpecies.ReduceMass) ;

		B	= 2.*Randn() - 1. ;		
		Angle	= sqrt(1. - (B*B)) ;
			
		RelativeVelPost[0]	= B * RelativeSpeed ;

		// AzimuthAngle is a random azimuth angle
		AzimuthAngle		= 2.*PI*Randn() ;

		RelativeVelPost[1]	= Angle * cos(AzimuthAngle) * RelativeSpeed ;
		RelativeVelPost[2]	= Angle * sin(AzimuthAngle) * RelativeSpeed ;	
		
		MassRatio[0]	= MultiSpecies.ReduceMass/h_Species[SelectSpeciesNo[1]].Mass ;
		MassRatio[1]	= MultiSpecies.ReduceMass/h_Species[SelectSpeciesNo[0]].Mass ;
		
			
		h_pDSMC->ParticleXVel[SelectParticleNo[0]]	= CenterMassVel[0] + RelativeVelPost[0]*MassRatio[1] ;
		h_pDSMC->ParticleYVel[SelectParticleNo[0]]	= CenterMassVel[1] + RelativeVelPost[1]*MassRatio[1] ;
		h_pDSMC->ParticleZVel[SelectParticleNo[0]]	= CenterMassVel[2] + RelativeVelPost[2]*MassRatio[1] ; 
		
		h_pDSMC->ParticleVibLevel [SelectParticleNo[0]] = 0 ;
		h_pDSMC->ParticleRotation [SelectParticleNo[0]] = 0 ;
		h_pDSMC->ParticleVibration[SelectParticleNo[0]] = 0 ;
				
		h_pDSMC->ParticleXVel[SelectParticleNo[1]]	= CenterMassVel[0] - RelativeVelPost[0]*MassRatio[0] ;
		h_pDSMC->ParticleYVel[SelectParticleNo[1]]	= CenterMassVel[1] - RelativeVelPost[1]*MassRatio[0] ;
		h_pDSMC->ParticleZVel[SelectParticleNo[1]]	= CenterMassVel[2] - RelativeVelPost[2]*MassRatio[0] ; 
			
}		
//==============================================================================================================

void Recombination (DSMC_DOMAIN		*h_pDomain ,
		DSMC_SPECIES		*h_Species ,
		DSMC_DSMC		*h_pDSMC ,
		int			*SelectParticleNo ,
		int			*SelectSpeciesNo ,
		double			*RelativeVel , 
		double			*pRelativeVelSq ,
		double			*pRelativeSpeed ,
		DSMC_MULTISPECIES	*pMultiSpecies ,
		DSMC_CHEMICALTCE	*h_ChemicalTCE ,
		int			CellNo ,
		double			Temp , 
		double			TransTemp,
		int 		ColClassNo,
		int 		ReactionClassNo,
		double  TotalEnergy,
		int ThirdBody,
		int		*pParticleNo ,
		int		*pParticleNum ,
		int *ReactionBuffer,
		int   **ReactionBufferSub,
		int   *SubCellNo,
		int    SubCellNum, 
		ofstream	&OutputDebug ){
			
		DSMC_MULTISPECIES	MultiSpecies ;

		double MassRatio[2],CenterMassVel[3],Velocity[3],InitTransEnergy,CollisionEnergy, ComplexRelativeEnergy;
		double  XIA, XIB,BetweenRelativeEnergy,ProbabilityRatio,ERM,RotationEnergy;
		double B,Angle,RelativeVelPost[3],AzimuthAngle,RelativeSpeed;
		int VibLevel,MaxLevel,Kspec,ParticleNum;
		
		
		 ReactionBuffer[h_Species[SelectSpeciesNo[0]].GroupNo] = ReactionBuffer[h_Species[SelectSpeciesNo[0]].GroupNo]-1 ;
		 ReactionBuffer[h_Species[SelectSpeciesNo[1]].GroupNo] = ReactionBuffer[h_Species[SelectSpeciesNo[1]].GroupNo]-1 ;
		 
		  
		 h_pDSMC->ReactionIndex[SelectParticleNo[0]] = 2 ;
		 h_pDSMC->ReactionIndex[SelectParticleNo[1]] = 3 ;
		 
		if (SubCellNum >0) { 
			
		 ReactionBufferSub [h_Species[SelectSpeciesNo[0]].GroupNo][SubCellNo[0]] = ReactionBufferSub [h_Species[SelectSpeciesNo[0]].GroupNo][SubCellNo[0]]-1;
		 ReactionBufferSub [h_Species[SelectSpeciesNo[1]].GroupNo][SubCellNo[1]] = ReactionBufferSub [h_Species[SelectSpeciesNo[1]].GroupNo][SubCellNo[1]]-1;
		}
		
		
		Kspec = h_pDSMC->ParticleSpeciesNo[ThirdBody];
		
		h_pDSMC->ParticleXCoord[SelectParticleNo[0]] = 	0.5*(h_pDSMC->ParticleXCoord[SelectParticleNo[0]]+h_pDSMC->ParticleXCoord[SelectParticleNo[1]]);
		h_pDSMC->ParticleYCoord[SelectParticleNo[0]] =  0.5*(h_pDSMC->ParticleYCoord[SelectParticleNo[0]]+h_pDSMC->ParticleYCoord[SelectParticleNo[1]]);
		h_pDSMC->ParticleZCoord[SelectParticleNo[0]] =	0.5*(h_pDSMC->ParticleZCoord[SelectParticleNo[0]]+h_pDSMC->ParticleZCoord[SelectParticleNo[1]]);

		MassRatio[0]	= pMultiSpecies->ReduceMass/h_Species[SelectSpeciesNo[1]].Mass ;
		MassRatio[1]	= pMultiSpecies->ReduceMass/h_Species[SelectSpeciesNo[0]].Mass ;	
		
		CenterMassVel[0]	= MassRatio[0]*(h_pDSMC->ParticleXVel[SelectParticleNo[0]]) + MassRatio[1]*(h_pDSMC->ParticleXVel[SelectParticleNo[1]]) ;
		CenterMassVel[1]	= MassRatio[0]*(h_pDSMC->ParticleYVel[SelectParticleNo[0]]) + MassRatio[1]*(h_pDSMC->ParticleYVel[SelectParticleNo[1]]) ;
		CenterMassVel[2]	= MassRatio[0]*(h_pDSMC->ParticleZVel[SelectParticleNo[0]]) + MassRatio[1]*(h_pDSMC->ParticleZVel[SelectParticleNo[1]]) ;
		
		InitTransEnergy	= 0.5*pMultiSpecies->ReduceMass*(*pRelativeVelSq) ;
		
//		OutputDebug << "SelectSpeciesNo1: " << SelectSpeciesNo[0] << endl ;

	  SelectSpeciesNo[0]= h_ChemicalTCE[ColClassNo].Post1Species[ReactionClassNo];
	  h_pDSMC->ParticleSpeciesNo[SelectParticleNo[0]] =  SelectSpeciesNo[0] ;
	  
//	  OutputDebug << "SelectSpeciesNo2: " << SelectSpeciesNo[0] << endl ;
//	  OutputDebug << "EffVibDOF: " << h_pDSMC->EffVibDOF[CellNo][SelectSpeciesNo[0]]<< endl ;
	    
	  TotalEnergy = h_ChemicalTCE[ColClassNo].HeatofReaction[ReactionClassNo];
	    
    Velocity[0] = CenterMassVel[0]-h_pDSMC->ParticleXVel[ThirdBody];
    Velocity[1] = CenterMassVel[1]-h_pDSMC->ParticleYVel[ThirdBody];
		Velocity[2] = CenterMassVel[2]-h_pDSMC->ParticleZVel[ThirdBody];
		
		CalculateMultiSpecies( &MultiSpecies , &h_Species[SelectSpeciesNo[0]] , &h_Species[Kspec] ) ;

    ComplexRelativeEnergy= 0.5* MultiSpecies.ReduceMass*(Velocity[0]*Velocity[0]+Velocity[1]*Velocity[1]+Velocity[2]*Velocity[2]);
    
    TotalEnergy	 = TotalEnergy + ComplexRelativeEnergy + InitTransEnergy + h_pDSMC->ParticleRotation[ThirdBody];
    
    XIA=0.5*(h_Species[SelectSpeciesNo[0]].RotDOF+h_Species[Kspec].RotDOF);
		XIB=( 2.5 - MultiSpecies.VisTempIndex ) ;
		
//			OutputDebug << "VibLevel1: " << h_pDSMC->ParticleVibLevel[SelectParticleNo[0]] << endl ;
//			OutputDebug << "Vibration1: " << h_pDSMC->ParticleVibration[SelectParticleNo[0]] << endl ;

//			if ( h_pDSMC->EffVibDOF[CellNo][SelectSpeciesNo[0]] > 0. ){
//			if ( h_pDSMC->EffVibDOF[CellNo][SelectSpeciesNo[0]] > 0.01 ){
  			if ( h_Species[SelectSpeciesNo[0]].VibMode > 0. ) {
		
    	VibLevel = 0 ;
			MaxLevel = TotalEnergy/( BOLTZ * h_Species[SelectSpeciesNo[0]].VibTemp );

			do {
				
				h_pDSMC->ParticleVibLevel[SelectParticleNo[0]] = Randn() * ( MaxLevel+0.99999 ) ;							
				h_pDSMC->ParticleVibration[SelectParticleNo[0]] = h_pDSMC->ParticleVibLevel[SelectParticleNo[0]] * (BOLTZ*h_Species[SelectSpeciesNo[0]].VibTemp);			
				ProbabilityRatio	= pow(( 1.0 - h_pDSMC->ParticleVibration[SelectParticleNo[0]]/TotalEnergy	),( XIB+XIA-1. ));
			
			}while ( ProbabilityRatio < Randn() ) ;
			
			TotalEnergy		= TotalEnergy - h_pDSMC->ParticleVibration[SelectParticleNo[0]] ;			
			
    }

//			OutputDebug << "VibLevel2: " << h_pDSMC->ParticleVibLevel[SelectParticleNo[0]] << endl ;
//			OutputDebug << "Vibration2: " << h_pDSMC->ParticleVibration[SelectParticleNo[0]] << endl ;

//			if ( h_pDSMC->EffVibDOF[CellNo][Kspec] > 0.){
//			if ( h_pDSMC->EffVibDOF[CellNo][Kspec] > 0.01 ){			
		  if ( h_Species[Kspec].VibMode > 0. ) {
		
		  TotalEnergy	= TotalEnergy + h_pDSMC->ParticleVibration[ThirdBody] ;	
    	VibLevel = 0 ;
			MaxLevel = TotalEnergy/( BOLTZ * h_Species[Kspec].VibTemp );

			do {
				
				h_pDSMC->ParticleVibLevel[ThirdBody] = Randn() * ( MaxLevel+0.99999 ) ;							
				h_pDSMC->ParticleVibration[ThirdBody] = h_pDSMC->ParticleVibLevel[ThirdBody] * (BOLTZ*h_Species[Kspec].VibTemp);			
				ProbabilityRatio	= pow(( 1.0 - h_pDSMC->ParticleVibration[ThirdBody]/TotalEnergy	),( XIB+XIA-1. ));
			
			}while ( ProbabilityRatio < Randn() ) ;
			
			TotalEnergy	 = TotalEnergy - h_pDSMC->ParticleVibration[ThirdBody] ;			
			
    }

		ERM	= LarsenBorgnakkeEnergyRatio( XIB-1. , XIA-1. ) ;
		
		BetweenRelativeEnergy = TotalEnergy * ERM  ;
			
		RotationEnergy = TotalEnergy - BetweenRelativeEnergy ;
		
//		OutputDebug << "Rotation1: " << h_pDSMC->ParticleRotation[SelectParticleNo[0]] << endl ;
//		OutputDebug << "RotDOF1: " << h_Species[SelectSpeciesNo[0]].RotDOF<< endl ;
		
		if (h_Species[Kspec].RotDOF > 0.01 ){
			
			ERM	= LarsenBorgnakkeEnergyRatio( 0.5*h_Species[SelectSpeciesNo[0]].RotDOF-1 , 0.5*h_Species[Kspec].RotDOF-1 ) ;				
			h_pDSMC->ParticleRotation[SelectParticleNo[0]] = RotationEnergy * ERM ;			
			h_pDSMC->ParticleRotation[ThirdBody] = RotationEnergy - h_pDSMC->ParticleRotation[SelectParticleNo[0]] ;
			
		}else{
			
			h_pDSMC->ParticleRotation[SelectParticleNo[0]] = RotationEnergy ;
			h_pDSMC->ParticleRotation[ThirdBody] = 0 ;
			
		}
		
//		OutputDebug << "Rotation2: " << h_pDSMC->ParticleRotation[SelectParticleNo[0]] << endl ;
//		OutputDebug << "RotDOF2: " << h_Species[SelectSpeciesNo[0]].RotDOF << endl ;

		MassRatio[0]	= MultiSpecies.ReduceMass/h_Species[Kspec].Mass ;
		MassRatio[1]	= MultiSpecies.ReduceMass/h_Species[SelectSpeciesNo[0]].Mass ;	

		CenterMassVel[0]	= MassRatio[0]*CenterMassVel[0] + MassRatio[1]*(h_pDSMC->ParticleXVel[ThirdBody]) ;
		CenterMassVel[1]	= MassRatio[0]*CenterMassVel[1] + MassRatio[1]*(h_pDSMC->ParticleYVel[ThirdBody]) ;
		CenterMassVel[2]	= MassRatio[0]*CenterMassVel[2] + MassRatio[1]*(h_pDSMC->ParticleZVel[ThirdBody]) ;
	
		RelativeSpeed = sqrt(2.*BetweenRelativeEnergy/MultiSpecies.ReduceMass) ;
		
		B	= 2.*Randn() - 1. ;		
		Angle	= sqrt(1. - (B*B)) ;
		
		RelativeVelPost[0]	= B * RelativeSpeed ;

		// AzimuthAngle is a random azimuth angle
		AzimuthAngle		= 2.*PI*Randn() ;
		
		RelativeVelPost[1]	= Angle * cos(AzimuthAngle) * RelativeSpeed ;
		RelativeVelPost[2]	= Angle * sin(AzimuthAngle) * RelativeSpeed ;	
		
		h_pDSMC->ParticleXVel[SelectParticleNo[0]]	= CenterMassVel[0] + RelativeVelPost[0]*MassRatio[1] ;
		h_pDSMC->ParticleYVel[SelectParticleNo[0]]	= CenterMassVel[1] + RelativeVelPost[1]*MassRatio[1] ;
		h_pDSMC->ParticleZVel[SelectParticleNo[0]]	= CenterMassVel[2] + RelativeVelPost[2]*MassRatio[1] ; 
		
		h_pDSMC->ParticleXVel[ThirdBody]	= CenterMassVel[0] - RelativeVelPost[0]*MassRatio[0] ;
		h_pDSMC->ParticleYVel[ThirdBody]	= CenterMassVel[1] - RelativeVelPost[1]*MassRatio[0] ;
		h_pDSMC->ParticleZVel[ThirdBody]	= CenterMassVel[2] - RelativeVelPost[2]*MassRatio[0] ;

		
		
}
//==============================================================================================================

double SelectParticle(	DSMC_SPECIES		*h_Species ,
			DSMC_DSMC		*h_pDSMC ,
			double			*RelativeVel , 
			double			*pRelativeVelSq ,
			double			*pRelativeSpeed ,
			int			*SelectParticleNo ,
			int			IndexCell1_Group1 ,
			int			IndexCell2_Group1 ,
			int			IndexCell1_Group2 ,
			int			IndexCell2_Group2 ,
			DSMC_MULTISPECIES	*pMultiSpecies ){
				
	int		Select , SpeciesNo1 , SpeciesNo2 ; 
	double		CrossSectionSpeed ;
	bool		buffer = true ;
	int  TrackingNum ;
	
	TrackingNum =0 ;
	do{
	// Select first particle.
			Select			= (IndexCell2_Group1 - 1.E-3) * Randn() ;
			SelectParticleNo[0]	= h_pDSMC->IndexParticle[IndexCell1_Group1 + Select] ;
	
	}while ( h_pDSMC->ReactionIndex[SelectParticleNo[0]] > 0.) ;
	
	/*do{
		// Select second particle.
		Select			= (IndexCell2_Group2 - 1.E-3) * Randn() ;
		SelectParticleNo[1]	= h_pDSMC->IndexParticle[IndexCell1_Group2 + Select] ;
		
		if ( SelectParticleNo[0] != SelectParticleNo[1] ){
			if ( h_pDSMC->ParticleLastCollide[SelectParticleNo[0]] != SelectParticleNo[1] || 
			     h_pDSMC->ParticleLastCollide[SelectParticleNo[1]] != SelectParticleNo[0] || 
			     h_pDSMC->ParticleLastCollide[SelectParticleNo[0]] < 0 || 
			     h_pDSMC->ParticleLastCollide[SelectParticleNo[1]] < 0 ){
				break ;
			}else if ( IndexCell2_Group2 <= 2 ){
				break ;
			}
		}
	}while ( buffer ) ;*/
		
	do{
		// Select second particle.
		Select			= (IndexCell2_Group2 - 1.E-3) * Randn() ;
		SelectParticleNo[1]	= h_pDSMC->IndexParticle[IndexCell1_Group2 + Select] ;
				
	}while ( SelectParticleNo[0] == SelectParticleNo[1] || h_pDSMC->ReactionIndex[SelectParticleNo[1]] > 0.) ;
	
	
	SpeciesNo1	= h_pDSMC->ParticleSpeciesNo[SelectParticleNo[0]] ;
	SpeciesNo2	= h_pDSMC->ParticleSpeciesNo[SelectParticleNo[1]] ;
	
	// Calculate properties and mean value between species 1 and species 2.
	CalculateMultiSpecies( pMultiSpecies , &h_Species[SpeciesNo1] , &h_Species[SpeciesNo2] ) ;
	
	
	// RelativeVel[3] are the components of the relative velocity
	RelativeVel[0]		= h_pDSMC->ParticleXVel[SelectParticleNo[0]] - h_pDSMC->ParticleXVel[SelectParticleNo[1]] ;
	RelativeVel[1]		= h_pDSMC->ParticleYVel[SelectParticleNo[0]] - h_pDSMC->ParticleYVel[SelectParticleNo[1]] ;
	RelativeVel[2]		= h_pDSMC->ParticleZVel[SelectParticleNo[0]] - h_pDSMC->ParticleZVel[SelectParticleNo[1]] ;
	
	(*pRelativeVelSq)	= RelativeVel[0]*RelativeVel[0] + RelativeVel[1]*RelativeVel[1] + RelativeVel[2]*RelativeVel[2] ;
	(*pRelativeSpeed)	= sqrt(*pRelativeVelSq) ;
	
	
	// the collision cross-section is based on eqn (4.63)
	CrossSectionSpeed	= (*pRelativeSpeed)*pMultiSpecies->CrossSection*
					(pow((2.*BOLTZ*pMultiSpecies->RefTemp/(pMultiSpecies->ReduceMass*(*pRelativeVelSq))) , (pMultiSpecies->VisTempIndex-0.5)))/
					pMultiSpecies->GammaValue ;

					
	return	CrossSectionSpeed ;
}

//==============================================================================================================

int GetSubCellNum( DSMC_CELL *_Cell , int ParticleNum , int SubcellModel , int Dimension ){
	int	SubCellNum , NewSubCellNum ;
	
	
	SubCellNum	= 0 ;
	NewSubCellNum	= 0 ;
	
	if ( Dimension == 2 || Dimension == 4 ){
		if ( ParticleNum < 9 ){
			SubCellNum	= 0 ;
		}else{
			SubCellNum	= sqrt((double)ParticleNum/2.) ;
		
			if ( SubcellModel == 1 ){
				return	SubCellNum ;	
			}else if ( SubcellModel == 2 ){
				NewSubCellNum	= (int)(_Cell->CharacteristicLength / (0.5*_Cell->InitMeanFreePath)) ;
			}
		
			if ( NewSubCellNum < 2 ){
				SubCellNum	= 0 ;	
			}else if ( SubCellNum > NewSubCellNum ){
				SubCellNum	= NewSubCellNum ;
			}
		}
	}else if ( Dimension ==3 ){
		if ( ParticleNum < 17 ){
			SubCellNum	= 0 ;
		}else{
			SubCellNum	= pow((double)ParticleNum/2. , 1./3.) ;
		
			if ( SubcellModel == 1 ){
				return	SubCellNum ;	
			}else if ( SubcellModel == 2 ){
				NewSubCellNum	= (int)(_Cell->CharacteristicLength / (0.5*_Cell->InitMeanFreePath)) ;
			}
		
			if ( NewSubCellNum < 2 ){
				SubCellNum	= 0 ;	
			}else if ( SubCellNum > NewSubCellNum ){
				SubCellNum	= NewSubCellNum ;
			}
		}
	}
	
	
	return	SubCellNum ;
}

//==============================================================================================================

void IndexSubCell2D( 	DSMC_DOMAIN	*h_pDomain , 
			DSMC_DSMC	*h_pDSMC ,
			DSMC_CELL	*_Cell ,
			int		CellNo , 
			int		**IndexCellSub1 , 
			int		**IndexCellSub2 , 
			int		*IndexParticleSub , 
			int		*IndexParticleSubCellNo , 
			int		SubCellNum , 
			int 		ParticleNum , 
			int		*GroupParticleNum1,
			int     *GroupParticleNum2,
			ofstream	&OutputDebug ){
				
	int	TotalSubCellNum , GroupNum , ParticleNo , XSubCellNo , YSubCellNo , SubCellNo ;
	int	*ParticleSub , *ParticleSubCellNo , No , CountParticleNum ; // *GroupParticleNum1 , *GroupParticleNum2 ,
	double	MaxXCoord , MinXCoord , MaxYCoord , MinYCoord , XSubCellSize , YSubCellSize ;
	
	
	ParticleSub		= new int[ParticleNum] ;
	ParticleSubCellNo	= new int[ParticleNum] ;

	
	
	GroupNum	= h_pDomain->SpeciesGroupNum ;
	TotalSubCellNum	= SubCellNum * SubCellNum ;
	MaxXCoord	= _Cell->MaxXCoord ;
	MinXCoord	= _Cell->MinXCoord ;
	MaxYCoord	= _Cell->MaxYCoord ;
	MinYCoord	= _Cell->MinYCoord ;
	XSubCellSize	= (MaxXCoord - MinXCoord)/SubCellNum ;
	YSubCellSize	= (MaxYCoord - MinYCoord)/SubCellNum ;
	
	
	// Debug.
	//OutputDebug << "MaxX: " << MaxXCoord << ", MinX: " << MinXCoord << ", MaxY: " << MaxYCoord << ", MinYCoord: " << MinYCoord << endl ;
	//OutputDebug << "TSC: " << TotalSubCellNum << ", CellNo: " << CellNo << endl ;
	
	
	for ( int i=0 ; i<GroupNum ; i++ ){
		GroupParticleNum1[i]	= 0 ;
		GroupParticleNum2[i]	= 0 ;
		
		for ( int j=0 ; j<TotalSubCellNum ; j++ ){
			IndexCellSub1[i][j]	= 0 ;
			IndexCellSub2[i][j]	= 0 ;
		}
	}
	for ( int i=0 ; i<ParticleNum ; i++ ){
		IndexParticleSub[i]		= 0 ;
		IndexParticleSubCellNo[i]	= 0 ;
		ParticleSub[i]			= 0 ;	
		ParticleSubCellNo[i]		= 0 ;
	}
	
	
	// Debug.
	//OutputDebug << "isc-1" << endl ;

	No	= 0 ;
	for ( int i=0 ; i<GroupNum ; i++ ){
		for ( int j=0 ; j<h_pDSMC->IndexCell2[i][CellNo] ; j++ ){
			// Debug.
			//OutputDebug << "isc-1.1, No:" << No << ", PN: " << ParticleNum << endl ;
			
			ParticleNo	= h_pDSMC->IndexCell1[i][CellNo] + j ;
			ParticleNo	= h_pDSMC->IndexParticle[ParticleNo] ;
			
			
			XSubCellNo	= (h_pDSMC->ParticleXCoord[ParticleNo] - MinXCoord) / XSubCellSize ;
			YSubCellNo	= (h_pDSMC->ParticleYCoord[ParticleNo] - MinYCoord) / YSubCellSize ;
			
			// Debug.
			//OutputDebug << "isc-1.2, XSub:" << XSubCellNo << ", YSub: " << YSubCellNo << endl ;
			
			
			if ( XSubCellNo < 0 )
				XSubCellNo = 0 ;
			if ( XSubCellNo >= SubCellNum )
				XSubCellNo = SubCellNum - 1 ;
			if ( YSubCellNo < 0 ) 
				YSubCellNo = 0 ;
			if ( YSubCellNo >= SubCellNum ) 
				YSubCellNo = SubCellNum - 1 ;
				
			SubCellNo	= YSubCellNo*SubCellNum + XSubCellNo ;
			
			// Debug.
			//OutputDebug << "isc-1.3, SubCellNo:" << SubCellNo << endl ;
			
			ParticleSub[No]		= ParticleNo ;
			ParticleSubCellNo[No]	= SubCellNo ;
			IndexCellSub2[i][SubCellNo]++ ;
			GroupParticleNum2[i]++ ;
			
			No++ ;
		}
	}
	
	
	// Debug.
	//OutputDebug << "isc-2" << endl ;
	
		CountParticleNum	= 0 ;
		for ( int i=0 ; i<GroupNum ; i++ ){
		GroupParticleNum1[i]	= CountParticleNum ;
		CountParticleNum	+= GroupParticleNum2[i] ;
	}
/*
	
	CountParticleNum	= 0 ;
	for ( int i=0 ; i<TotalSubCellNum  ; i++ ){
		for ( int j=0 ; j< GroupNum; j++ ){
		IndexCellSub1[j][i]	= CountParticleNum ;
		CountParticleNum	+= IndexCellSub2[j][i] ;
		IndexCellSub2[j][i]	= 0 ;
		}
	}
*/	
	// Debug.
	//OutputDebug << "isc-3" << endl ;
	
	
	for ( int i=0 ; i<GroupNum ; i++ ){
		CountParticleNum	= GroupParticleNum1[i] ;
		
		for ( int j=0 ; j<TotalSubCellNum ; j++ ){
			IndexCellSub1[i][j]	= CountParticleNum ;
			CountParticleNum	+= IndexCellSub2[i][j] ;
			
			IndexCellSub2[i][j]	= 0 ;
		}
	}
	
	
	// Debug.
	//OutputDebug << "isc-4" << endl ;
	
	No	= 0 ;
	for ( int i=0 ; i<GroupNum ; i++ ){
		for ( int j=0 ; j<h_pDSMC->IndexCell2[i][CellNo] ; j++ ){
			SubCellNo	= ParticleSubCellNo[No] ;
		
			ParticleNo	= IndexCellSub1[i][SubCellNo] + IndexCellSub2[i][SubCellNo] ;
		
			IndexParticleSub[ParticleNo]	= ParticleSub[No];
			//IndexParticleSubCellNo[ParticleNo]	= ParticleSubCellNo[No] ;
			IndexParticleSubCellNo[No]	= ParticleSubCellNo[No] ;
			IndexCellSub2[i][SubCellNo]++ ;
		
			No++ ;
		}
	}
	
	
	// Debug.
	/*No	= 0 ;
	OutputDebug << "=============" << endl ;
	for ( int i=0 ; i<GroupNum ; i++ ){
		for ( int j=0 ; j<h_pDSMC->IndexCell2[i][CellNo] ; j++ ){
			OutputDebug << "ParticleNo: " << setw(3) << No << setw(12) << h_pDSMC->IndexParticle[h_pDSMC->IndexCell1[i][CellNo]+j] << setw(4) << IndexParticleSubCellNo[No] 
				<< setw(12) << ParticleSub[No] << setw(12) << IndexParticleSub[No] << endl ;
		
		
		
			No++ ;
		}
	}
	OutputDebug << "==" << endl ;
	for ( int i=0 ; i<GroupNum ; i++ ){
		for ( int j=0 ; j<TotalSubCellNum ; j++ ){
			OutputDebug << setw(4) << j << setw(8) << IndexCellSub1[i][j] << setw(8) << IndexCellSub2[i][j] << endl ;
		}
	}
	OutputDebug << "=============" << endl ;*/
	
	// Debug.
	//OutputDebug << "isc-5" << endl ;
	
//	delete [] GroupParticleNum1 ;
//	delete [] GroupParticleNum2 ;
	delete [] ParticleSub ;
	delete [] ParticleSubCellNo ;
	
	
	// Debug.
	//OutputDebug << "isc-6" << endl ;
}


void IndexSubCell3D( 	DSMC_DOMAIN	*h_pDomain , 
			DSMC_DSMC	*h_pDSMC ,
			DSMC_CELL	*_Cell ,
			int		CellNo , 
			int		**IndexCellSub1 , 
			int		**IndexCellSub2 , 
			int		*IndexParticleSub , 
			int		*IndexParticleSubCellNo , 
			int		SubCellNum , 
			int 		ParticleNum , 
			int	*GroupParticleNum1 ,
			int *GroupParticleNum2 ,
			ofstream	&OutputDebug ){
				
	int	TotalSubCellNum , GroupNum , ParticleNo , XSubCellNo , YSubCellNo , ZSubCellNo , SubCellNo ; 
	int	*ParticleSub , *ParticleSubCellNo ,  No , CountParticleNum ;
	double	MaxXCoord , MinXCoord , MaxYCoord , MinYCoord , MaxZCoord , MinZCoord , XSubCellSize , YSubCellSize , ZSubCellSize ;
	
	
	ParticleSub		= new int[ParticleNum] ;
	ParticleSubCellNo	= new int[ParticleNum] ;
//	GroupParticleNum1	= new int[h_pDomain->SpeciesGroupNum] ;
//	GroupParticleNum2	= new int[h_pDomain->SpeciesGroupNum] ;
	
	
	GroupNum	= h_pDomain->SpeciesGroupNum ;
	TotalSubCellNum	= SubCellNum * SubCellNum * SubCellNum ;
	MaxXCoord	= _Cell->MaxXCoord ;
	MinXCoord	= _Cell->MinXCoord ;
	MaxYCoord	= _Cell->MaxYCoord ;
	MinYCoord	= _Cell->MinYCoord ;
	MaxZCoord	= _Cell->MaxZCoord ;
	MinZCoord	= _Cell->MinZCoord ;
	XSubCellSize	= (MaxXCoord - MinXCoord)/SubCellNum ;
	YSubCellSize	= (MaxYCoord - MinYCoord)/SubCellNum ;
	ZSubCellSize	= (MaxZCoord - MinZCoord)/SubCellNum ;
	XSubCellNo	= 0 ;
	YSubCellNo	= 0 ;
	ZSubCellNo	= 0 ;
	SubCellNo	= 0 ;
	
	
	for ( int i=0 ; i<GroupNum ; i++ ){
		GroupParticleNum1[i]	= 0 ;
		GroupParticleNum2[i]	= 0 ;
		
		for ( int j=0 ; j<TotalSubCellNum ; j++ ){
			IndexCellSub1[i][j]	= 0 ;
			IndexCellSub2[i][j]	= 0 ;
		}
	}
	for ( int i=0 ; i<ParticleNum ; i++ ){
		IndexParticleSub[i]		= 0 ;
		IndexParticleSubCellNo[i]	= 0 ;
		ParticleSub[i]			= 0 ;	
		ParticleSubCellNo[i]		= 0 ;
	}
	
	
	No	= 0 ;
	for ( int i=0 ; i<GroupNum ; i++ ){
		for ( int j=0 ; j<h_pDSMC->IndexCell2[i][CellNo] ; j++ ){
			ParticleNo	= h_pDSMC->IndexCell1[i][CellNo] + j ;
			ParticleNo	= h_pDSMC->IndexParticle[ParticleNo] ;
			
			
			XSubCellNo	= (h_pDSMC->ParticleXCoord[ParticleNo] - MinXCoord) / XSubCellSize ;
			YSubCellNo	= (h_pDSMC->ParticleYCoord[ParticleNo] - MinYCoord) / YSubCellSize ;
			ZSubCellNo	= (h_pDSMC->ParticleZCoord[ParticleNo] - MinZCoord) / ZSubCellSize ;
			
			
			if ( XSubCellNo < 0 )
				XSubCellNo = 0 ;
			if ( XSubCellNo >= SubCellNum )
				XSubCellNo = SubCellNum - 1 ;
			if ( YSubCellNo < 0 ) 
				YSubCellNo = 0 ;
			if ( YSubCellNo >= SubCellNum ) 
				YSubCellNo = SubCellNum - 1 ;
			if ( ZSubCellNo < 0 ) 
				ZSubCellNo = 0 ;
			if ( ZSubCellNo >= SubCellNum ) 
				ZSubCellNo = SubCellNum - 1 ;
				
				
			SubCellNo	= ZSubCellNo*SubCellNum*SubCellNum + YSubCellNo*SubCellNum + XSubCellNo ;
			
			
			ParticleSub[No]		= ParticleNo ;
			ParticleSubCellNo[No]	= SubCellNo ;
			IndexCellSub2[i][SubCellNo]++ ;
			GroupParticleNum2[i]++ ;
			
			No++ ;
		}
	}
	
	
	CountParticleNum	= 0 ;
	for ( int i=0 ; i<GroupNum ; i++ ){
		GroupParticleNum1[i]	= CountParticleNum ;
		CountParticleNum	+= GroupParticleNum2[i] ;
	}
	
	
	for ( int i=0 ; i<GroupNum ; i++ ){
		CountParticleNum	= GroupParticleNum1[i] ;
		
		for ( int j=0 ; j<TotalSubCellNum ; j++ ){
			IndexCellSub1[i][j]	= CountParticleNum ;
			CountParticleNum	+= IndexCellSub2[i][j] ;
			
			IndexCellSub2[i][j]	= 0 ;
		}
	}
	
	
	No	= 0 ;
	for ( int i=0 ; i<GroupNum ; i++ ){
		for ( int j=0 ; j<h_pDSMC->IndexCell2[i][CellNo] ; j++ ){
			SubCellNo	= ParticleSubCellNo[No] ;
		
			ParticleNo	= IndexCellSub1[i][SubCellNo] + IndexCellSub2[i][SubCellNo] ;
		
			IndexParticleSub[ParticleNo]	= ParticleSub[No] ;
			//IndexParticleSubCellNo[ParticleNo]	= ParticleSubCellNo[No] ;
			IndexParticleSubCellNo[No]	= ParticleSubCellNo[No] ;
			IndexCellSub2[i][SubCellNo]++ ;
		
			No++ ;
		}
	}
	
	
//	delete [] GroupParticleNum1 ;
//	delete [] GroupParticleNum2 ;
	delete [] ParticleSub ;
	delete [] ParticleSubCellNo ;
}

//==============================================================================================================

double SelectParticleSubCell(	DSMC_SPECIES		*h_Species ,
				DSMC_DSMC		*h_pDSMC ,
				double			*RelativeVel , 
				double			*pRelativeVelSq ,
				double			*pRelativeSpeed ,
				int			*SelectParticleNo ,
				int			IndexCell1_Group1 ,
				int			IndexCell2_Group1 ,
				int			IndexCell1_Group2 ,
				int			IndexCell2_Group2 ,
				DSMC_MULTISPECIES	*pMultiSpecies , 
				int			**IndexCellSub1 , 
				int			**IndexCellSub2 , 
				int			*IndexParticleSub ,
				int			*IndexParticleSubCellNo ,
				int			GroupNo1 ,
				int			GroupNo2 , 
				int			SubCellNum , 
				int         *GroupParticleNum1,
				int 		*SubCellNo,
				int			**ReactionBufferSub,
				int			TotalSubCellNum ){
				
	int		Select , SpeciesNo1 , SpeciesNo2 ;
	int		SubCellShift , InitSubCellNo ;
	double		CrossSectionSpeed ;
	bool		buffer = true ;
	int TrackingNum ;
	
	
	SubCellShift	= 0 ;
	TrackingNum  =0 ;
	
	do{
	// Select first particle.
	Select			= (IndexCell2_Group1 - 1.E-3) * Randn() ;
	SelectParticleNo[0]	= h_pDSMC->IndexParticle[IndexCell1_Group1 + Select] ;
	SubCellNo[0]		= IndexParticleSubCellNo[ GroupParticleNum1[GroupNo1] + Select] ;
			
	}while ( h_pDSMC->ReactionIndex[SelectParticleNo[0]] > 0.) ;

	
	/*if ( IndexCellSub2[GroupNo2][SubCellNo] >= 2 ){
		do{
			// Select second particle.
			//Select		= (IndexCell2_Group2 - 1.E-3) * Randn() ;
			//SelectParticleNo[1]	= h_pDSMC->IndexParticle[IndexCell1_Group2 + Select] ;
			Select			= (IndexCellSub2[GroupNo2][SubCellNo] - 1.E-3) * Randn() ;
			SelectParticleNo[1]	= IndexParticleSub[IndexCellSub1[GroupNo2][SubCellNo]+Select] ;
		
		}while ( SelectParticleNo[0] == SelectParticleNo[1] ) ;
	}else{
		InitSubCellNo	= SubCellNo ;
		
		do{
			do{
				if ( SubCellShift == 0 ){
					SubCellShift++ ;	
				}else if ( SubCellShift < 0 ){
					SubCellShift *= -1 ;
					SubCellShift++ ;
				}else if ( SubCellShift > 0 ){
					SubCellShift *= -1 ;	
				}
			
				SubCellNo	= InitSubCellNo + SubCellShift ;
				
			}while ( SubCellNo < 0 || SubCellNo >= TotalSubCellNum ) ;
		}while ( IndexCellSub2[GroupNo2][SubCellNo] == 0 ) ;
		
		Select			= (IndexCellSub2[GroupNo2][SubCellNo] - 1.E-3) * Randn() ;
		SelectParticleNo[1]	= IndexParticleSub[IndexCellSub1[GroupNo2][SubCellNo]+Select] ;
	}*/
	
	if ( ((GroupNo1 != GroupNo2) && (ReactionBufferSub[GroupNo2][SubCellNo[0]] >= 1)) || 
	     ((GroupNo1 == GroupNo2) && (ReactionBufferSub[GroupNo2][SubCellNo[0]] >= 2)) ){
		do{
			// Select second particle.
			Select			= (IndexCellSub2[GroupNo2][SubCellNo[0]] - 1.E-3) * Randn() ;
			SelectParticleNo[1]	= IndexParticleSub[IndexCellSub1[GroupNo2][SubCellNo[0]]+Select] ;
		
		}while ( SelectParticleNo[0] == SelectParticleNo[1] || h_pDSMC->ReactionIndex[SelectParticleNo[1]] > 0.) ;
		
		SubCellNo[1]= SubCellNo[0];
			 		
	}else{
		
		InitSubCellNo	= SubCellNo[0] ;

			do{
				do{
					if ( SubCellShift == 0 ){
						SubCellShift++ ;	
					}else if ( SubCellShift < 0 ){
						SubCellShift *= -1 ;
						SubCellShift++ ;
					}else if ( SubCellShift > 0 ){
						SubCellShift *= -1 ;	
					}			
					SubCellNo[1]	= InitSubCellNo + SubCellShift ;
				
				}while ( SubCellNo[1] < 0 || SubCellNo[1] >= TotalSubCellNum ) ;
			}while ( ReactionBufferSub[GroupNo2][SubCellNo[1]] == 0 ) ;
				
			do{
				
				Select			= (IndexCellSub2[GroupNo2][SubCellNo[1]] - 1.E-3) * Randn() ;
				SelectParticleNo[1]	= IndexParticleSub[IndexCellSub1[GroupNo2][SubCellNo[1]]+Select] ;
				
			}while ( h_pDSMC->ReactionIndex[SelectParticleNo[1]] > 0.) ;		
	}
	
	
	SpeciesNo1	= h_pDSMC->ParticleSpeciesNo[SelectParticleNo[0]] ;
	SpeciesNo2	= h_pDSMC->ParticleSpeciesNo[SelectParticleNo[1]] ;
	
	
	// Calculate properties and mean value between species 1 and species 2.
	CalculateMultiSpecies( pMultiSpecies , &h_Species[SpeciesNo1] , &h_Species[SpeciesNo2] ) ;
	
	
	// RelativeVel[3] are the components of the relative velocity
	RelativeVel[0]		= h_pDSMC->ParticleXVel[SelectParticleNo[0]] - h_pDSMC->ParticleXVel[SelectParticleNo[1]] ;
	RelativeVel[1]		= h_pDSMC->ParticleYVel[SelectParticleNo[0]] - h_pDSMC->ParticleYVel[SelectParticleNo[1]] ;
	RelativeVel[2]		= h_pDSMC->ParticleZVel[SelectParticleNo[0]] - h_pDSMC->ParticleZVel[SelectParticleNo[1]] ;
	
	(*pRelativeVelSq)	= RelativeVel[0]*RelativeVel[0] + RelativeVel[1]*RelativeVel[1] + RelativeVel[2]*RelativeVel[2] ;
	(*pRelativeSpeed)	= sqrt(*pRelativeVelSq) ;
	
	
	// the collision cross-section is based on eqn (4.63)
	CrossSectionSpeed	= (*pRelativeSpeed)*pMultiSpecies->CrossSection*
					(pow((2.*BOLTZ*pMultiSpecies->RefTemp/(pMultiSpecies->ReduceMass*(*pRelativeVelSq))) , (pMultiSpecies->VisTempIndex-0.5)))/
					pMultiSpecies->GammaValue ;
					
	return	CrossSectionSpeed ;
}

//==============================================================================================================

void CalculateMultiSpecies( DSMC_MULTISPECIES *pMultiSpecies , DSMC_SPECIES *_pSpecies1 , DSMC_SPECIES *_pSpecies2 ){
	
	pMultiSpecies->CrossSection	= 0.25*PI*(_pSpecies1->Diameter + _pSpecies2->Diameter)*(_pSpecies1->Diameter + _pSpecies2->Diameter) ;
	pMultiSpecies->RefTemp		= 0.5 * (_pSpecies1->RefTemp + _pSpecies2->RefTemp ) ;
	pMultiSpecies->VisTempIndex	= 0.5 * (_pSpecies1->VisTempIndex + _pSpecies2->VisTempIndex ) ;
	pMultiSpecies->VSSParameter	= 0.5 * (_pSpecies1->VSSParameter + _pSpecies2->VSSParameter ) ;
	pMultiSpecies->ReduceMass	= (_pSpecies1->Mass/(_pSpecies1->Mass+_pSpecies2->Mass)) * _pSpecies2->Mass ;
	pMultiSpecies->GammaValue	= GammaFunction( 2.5 - pMultiSpecies->VisTempIndex ) ;
}

//==============================================================================================================

void Inelrv(	DSMC_DOMAIN		*h_pDomain ,
		DSMC_SPECIES		*h_Species ,
		DSMC_DSMC		*h_pDSMC ,
		int			*SelectParticleNo ,
		int			*SelectSpeciesNo ,
		double			*RelativeVel , 
		double			*pRelativeVelSq ,
		double			*pRelativeSpeed ,
		DSMC_MULTISPECIES	*pMultiSpecies ,
		int			CellNo ,
		double			Temp , 
		double			TransTemp,
		ofstream	&OutputDebug ){

	int			ParticleNo , SpeciesNo , OtherSpeciesNo , MaxLevel ;
	double			InitTransEnergy , InitEnergy , FinalEnergy , DividEnergy , PostTransEnergy , CollTemp ;
	double			XIA , XIB , ERM , RedistributedProbability , ProbabilityRatio , Buffer , pi ,A,B;
	double      TotalEnergy , RefTemp ,SPVM ;
	DSMC_MULTISPECIES	MultiSpecies ;

	// InitTransEnergy is the initial translational energy
	// InitEnergy is the initial energy in the active rotational modes
	// FinalEnergy is the final energy in these modes
	// DividEnergy is the energy to be divided
	// XIB is th number of modes in the redistribution
	// IRTR is 0,1 if no,any rotational energy redistribution is made
	// IRTV is 0,1 if no,any vibrational energy redistribution is made

	InitTransEnergy	= 0.5*pMultiSpecies->ReduceMass*(*pRelativeVelSq) ;
	InitEnergy	= 0. ;
	FinalEnergy	= 0. ;
	TotalEnergy = InitTransEnergy ;
	DividEnergy	= InitTransEnergy ;
	XIB		= 2.5 - pMultiSpecies->VisTempIndex ;
	pi		= PI ;
	MaxLevel	= 0 ;
	RefTemp = 0.;
  SPVM    = 0.;

	// consider the molecules in turn
	for ( int i=0 ; i<2 ; i++ ){
		ParticleNo	= SelectParticleNo[i] ;
	
		if ( i == 0 ){
			SpeciesNo	= SelectSpeciesNo[0] ;
			OtherSpeciesNo	= SelectSpeciesNo[1] ;
		}else{
			SpeciesNo	= SelectSpeciesNo[1] ;
			OtherSpeciesNo	= SelectSpeciesNo[0] ;
		}

		CalculateMultiSpecies( &MultiSpecies , &h_Species[SpeciesNo] , &h_Species[OtherSpeciesNo] ) ;
		
		
		if ( h_Species[SpeciesNo].VibMode  > 0 ){
			// When h_Species[SpeciesNo].VibModel = 1, calculate the Zv by Parker's model (Eq. A5)
			// When h_Species[SpeciesNo].VibModel = 2, calculate the Zv by Millikan and White's model (Eq. 6.53 or A6)
			// When h_Species[SpeciesNo].VibModel = 3, assign the Zv as a constant
			
			if ( h_Species[SpeciesNo].VibModel == 1 ){
				// RedistributedProbability is the probability that vibration is redistributed to molecule "ParticleNo"
				Buffer				= 79.8/TransTemp ;
				RedistributedProbability	= (1.+pow(pi,1.5)/2.*sqrt(Buffer)+(pi*pi/4.+pi)*Buffer)/21. ;
				
			}else if ( h_Species[SpeciesNo].VibModel == 2 ){
				
/*				
				if ( h_pDSMC->EffVibDOF[CellNo][SpeciesNo] > 0.01 ){
					
					// the collision temperature CollTemp is calculated from eqn (11.34)
					//CollTemp	= (InitTransEnergy + h_pDSMC->ParticleVibration[ParticleNo])/
					//			((2.5 - MultiSpecies.VisTempIndex + 0.5*h_pDSMC->EffVibDOF[CellNo][SpeciesNo])*BOLTZ) ;
					CollTemp	= (DividEnergy + h_pDSMC->ParticleVibration[ParticleNo])/
								((2.5 - MultiSpecies.VisTempIndex + 0.5*h_pDSMC->EffVibDOF[CellNo][SpeciesNo])*BOLTZ) ;
				}else{
					CollTemp	= Temp ;
				}
*/
				if ( h_Species[SpeciesNo].VibMode > 0. ){
					
					MaxLevel = InitTransEnergy + h_pDSMC->ParticleVibration[ParticleNo]/(BOLTZ*h_Species[SpeciesNo].VibTemp);					
					CollTemp = MaxLevel*h_Species[SpeciesNo].VibTemp /(3.5-MultiSpecies.VisTempIndex);
				
				}else{
					CollTemp	= Temp ;
				}
				
				Buffer	= h_Species[SpeciesNo].VibConstant2[OtherSpeciesNo] * pow(CollTemp , -0.33333) ;
				
				// the vibrational collision number has been calculated from eqn (6.53)
				if ( Buffer < 50. ){
					RedistributedProbability	= 1./((h_Species[SpeciesNo].VibConstant1[OtherSpeciesNo]/pow(CollTemp , pMultiSpecies->VisTempIndex))*exp(Buffer)) ;
				}else{
					RedistributedProbability	= 1./1.E+7 ;
				}

			}else if ( h_Species[SpeciesNo].VibModel == 3 ){
				// the vibrational relaxation collision number is a constant
				RedistributedProbability	= 1./h_Species[SpeciesNo].VibConstant1[OtherSpeciesNo] ;
			
			}else if ( h_Species[SpeciesNo].VibModel == 4 ){
																		
				B = (h_Species[SpeciesNo].DisTemp) / (h_Species[SpeciesNo].VibTemp) ;	
				A = (h_Species[SpeciesNo].DisTemp) / (TransTemp);
					
				RedistributedProbability = 1. / ( pow( A , pMultiSpecies->VisTempIndex ) * pow(( h_Species[SpeciesNo].VibTemp * pow (B , -pMultiSpecies->VisTempIndex)), ((pow(A,0.3333333)-1)/(pow(B,0.33333)-1))));
			
			}else if ( h_Species[SpeciesNo].VibModel == 5 ){
				
				TotalEnergy	+= h_pDSMC->ParticleVibration[ParticleNo] ;
				MaxLevel	= TotalEnergy / (BOLTZ*h_Species[SpeciesNo].VibTemp) ;
				CollTemp = MaxLevel*h_Species[SpeciesNo].VibTemp/(3.5-MultiSpecies.VisTempIndex);
				RefTemp = 7*h_Species[SpeciesNo].VibTemp ;
				
				B = (h_Species[SpeciesNo].DisTemp) / (RefTemp) ;	
				A = (h_Species[SpeciesNo].DisTemp) / (CollTemp);
				
				SPVM=pow(h_Species[SpeciesNo].VibConstant1[OtherSpeciesNo]/RefTemp,MultiSpecies.VisTempIndex)*exp(h_Species[SpeciesNo].VibConstant2[OtherSpeciesNo] / pow(RefTemp , 0.33333333));
					
				RedistributedProbability = 1. / ( pow( A , pMultiSpecies->VisTempIndex ) * pow(( SPVM * pow (B , -pMultiSpecies->VisTempIndex)), ((pow(A,0.3333333)-1)/(pow(B,0.33333)-1))));
																		
			}// End Vibrational Model


			if ( RedistributedProbability > Randn() ){
				DividEnergy	+= h_pDSMC->ParticleVibration[ParticleNo] ;
				InitEnergy	+= h_pDSMC->ParticleVibration[ParticleNo] ;
				
				if ( h_pDomain->VibrationalModel == 1 ){
					// When h_Species[SpeciesNo].VibConstant2[OtherSpeciesNo]=-1, the Zv is a constant for classical version and calculated by Eq. 5.53
					if ( h_Species[SpeciesNo].VibConstant2[OtherSpeciesNo] == -1. ){
						for ( int j=0 ; j<h_Species[SpeciesNo].VibMode ; j++ )
							Buffer	= (h_Species[SpeciesNo].VibTemp/h_pDSMC->ParticleEffTemp[ParticleNo])/(exp(h_Species[SpeciesNo].VibTemp/h_pDSMC->ParticleEffTemp[ParticleNo])-1.) ;

						h_pDSMC->ParticleEffTemp[ParticleNo]	= (DividEnergy/BOLTZ)/((XIB+0.5*h_pDSMC->EffVibDOF[CellNo][SpeciesNo])+Buffer) ;
						h_pDSMC->EffVibDOF[CellNo][SpeciesNo]	= (2.*h_Species[SpeciesNo].VibTemp/h_pDSMC->ParticleEffTemp[ParticleNo])/(exp(h_Species[SpeciesNo].VibTemp/h_pDSMC->ParticleEffTemp[ParticleNo])-1.) ;
					}


					if (h_pDSMC->EffVibDOF[CellNo][SpeciesNo] == 2.){
						ERM	= 1. - pow(Randn() , (1./XIB)) ;
					}else{
						// apply the general Larsen-Borgnakke distribution function
						XIA	= 0.5 * h_pDSMC->EffVibDOF[CellNo][SpeciesNo] ;
						ERM	= LarsenBorgnakkeEnergyRatio( XIA-1. , XIB-1. ) ;
					}
					
					h_pDSMC->ParticleVibration[ParticleNo]	= ERM * DividEnergy ;
					
					// the available energy is reduced accordingly
					DividEnergy	-= h_pDSMC->ParticleVibration[ParticleNo] ;
					FinalEnergy	+= h_pDSMC->ParticleVibration[ParticleNo] ;
					
				}else if ( h_pDomain->VibrationalModel == 2 ){
					// MaxLevel is the maximum level within the available energy (eqn (5.62))
					MaxLevel	= DividEnergy / (BOLTZ*h_Species[SpeciesNo].VibTemp) ;
					
					do {
						h_pDSMC->ParticleVibLevel[ParticleNo]	= Randn()*(MaxLevel+0.99999) ;
						// the above statement chooses a level uniformly from 0 to "MaxLevel"
						
						h_pDSMC->ParticleVibration[ParticleNo]	= h_pDSMC->ParticleVibLevel[ParticleNo]*BOLTZ*h_Species[SpeciesNo].VibTemp ;
					
						// ProbabilityRatio is the probability ratio (eqn (5.61)) ;
						ProbabilityRatio	= pow((1.-h_pDSMC->ParticleVibration[ParticleNo]/DividEnergy) , (1.5-MultiSpecies.VisTempIndex)) ;

					}while ( ProbabilityRatio < Randn() ) ;
						
					DividEnergy	-= h_pDSMC->ParticleVibration[ParticleNo] ;
					FinalEnergy	+= h_pDSMC->ParticleVibration[ParticleNo] ;
										
				}
			}// End if ( RedistributedProbability > Randn() ) (Vibration)
		}// End if ( h_Species[SpeciesNo].VibMode  > 0 )

		
		if ( h_Species[SpeciesNo].RotDOF > 0.01 ){
			// When h_Species[SpeciesNo].RotModel = 1, calculate the Zv by Parker's model (Eq. A5)
			// When h_Species[SpeciesNo].RotModel = 2, assign the Zv as a constant)
			
			if ( h_Species[SpeciesNo].RotModel == 1 ) {
				Buffer				= 79.8/TransTemp ;
				RedistributedProbability	= (1.+pow(PI,1.5)/2.*sqrt(Buffer)+(PI*PI/4.+PI)*Buffer)/21. ;
				
			}else if ( h_Species[SpeciesNo].RotModel == 2 ){
				RedistributedProbability	= 1./h_Species[SpeciesNo].RotRelaxationNum[OtherSpeciesNo] ;
				
			}else if ( h_Species[SpeciesNo].RotModel == 3 ){
				Buffer				= 80.0/TransTemp ;
				RedistributedProbability	= (1.+(pow(PI,0.5)/2.)*sqrt(Buffer)+(PI*PI/4.+PI)*Buffer)/15.7 ;								
			}
			
			// RedistributedProbability is the probability that rotation is redistributed to molecule "ParticleNo"
			if ( RedistributedProbability > Randn() ){
				DividEnergy	+= h_pDSMC->ParticleRotation[ParticleNo] ;
				InitEnergy	+= h_pDSMC->ParticleRotation[ParticleNo] ;
				
				if (h_Species[SpeciesNo].RotDOF == 2.){
					ERM	= 1. - pow(Randn() , (1./XIB)) ;
				}else{
					// apply the general Larsen-Borgnakke distribution function
					XIA	= 0.5 * h_Species[SpeciesNo].RotDOF ;
					ERM	= LarsenBorgnakkeEnergyRatio( XIA-1. , XIB-1. ) ;
				}
				
				h_pDSMC->ParticleRotation[ParticleNo]		= ERM * DividEnergy ;
				
				// the available energy is reduced accordingly
				DividEnergy	-= h_pDSMC->ParticleRotation[ParticleNo] ;
				FinalEnergy	+= h_pDSMC->ParticleRotation[ParticleNo] ;
				
			}// End if ( RedistributedProbability > Randn() ) (Rotation)
		}// End if ( h_Species[SpeciesNo].RotDOF > 0.01 ) ;
	}// End of consider the molecules.

	// PostTransEnergy  is the post-collision translational energy
	PostTransEnergy		= InitTransEnergy + InitEnergy - FinalEnergy ;


	// adjust (*pRelativeSpeed) and, for the VSS model, RelativeSpeed[3] for the change in energy
	//(*pRelativeVelSq)	= 2. * PostTransEnergy / pMultiSpecies->ReduceMass ;
	Buffer			= sqrt(2.*PostTransEnergy/pMultiSpecies->ReduceMass) ;
	
	if ( fabs(pMultiSpecies->VSSParameter - 1.) < 1.E-3 ){
		(*pRelativeSpeed)	= Buffer ;
	}else{
		RelativeVel[0]		= RelativeVel[0] * Buffer / (*pRelativeSpeed) ;
		RelativeVel[1]		= RelativeVel[1] * Buffer / (*pRelativeSpeed) ;
		RelativeVel[2]		= RelativeVel[2] * Buffer / (*pRelativeSpeed) ;
		(*pRelativeSpeed)	= Buffer ;
	}
}

//==============================================================================================================

double LarsenBorgnakkeEnergyRatio( double XMA , double XMB ){
	double	ERM , P ;
	
	// selects a Larsen-Borgnakke energy ratio using eqn (11.9)
	do{
		ERM	= Randn() ;
		if ( XMA < 1.E-6 || XMB < 1.E-6 ){
			if ( XMA < 1.E-6 && XMB < 1.E-6 ) break ;
			if ( XMA < 1.E-6 ) P = pow((1.-ERM) , XMB) ;
			if ( XMB < 1.E-6 ) P = pow((1.-ERM) , XMA) ;
		}else{
			P	= pow(((XMA+XMB)*ERM/XMA) , XMA) * pow(((XMA+XMB)*(1.-ERM)/XMB) , XMB) ;
		}
	}while ( P < Randn() ) ;

	return	ERM ;
}

//==============================================================================================================

void Elastic(	DSMC_DSMC		*h_pDSMC ,
		double			*RelativeVel , 
		double			*pRelativeSpeed ,
		int			*SelectParticleNo ,
		DSMC_MULTISPECIES	*pMultiSpecies , 
		double			Mass1 ,
		double			Mass2 ){
	
	double		MassRatio[2] ,  RelativeVelPost[3] , CenterMassVel[3] ;
	double		B , D , AzimuthAngle , CosAzimuthAngle , SinAzimuthAngle , Angle ;

	// RelativeVelPost[3] are the post-collision components of the relative velocity
	RelativeVelPost[0]	= 0. ;
	RelativeVelPost[1]	= 0. ;
	RelativeVelPost[2]	= 0. ;

	AzimuthAngle		= 0. ;
	Angle			= 0. ;
	CosAzimuthAngle		= 0. ;
	SinAzimuthAngle		= 0. ;
	D			= 0. ;


	MassRatio[0]	= pMultiSpecies->ReduceMass/Mass2 ;
	MassRatio[1]	= pMultiSpecies->ReduceMass/Mass1 ;


	// CenterMassVel defines the components of the centre-of-mass velocity, eqn (2.1)
	CenterMassVel[0]	= MassRatio[0]*(h_pDSMC->ParticleXVel[SelectParticleNo[0]]) + MassRatio[1]*(h_pDSMC->ParticleXVel[SelectParticleNo[1]]) ;
	CenterMassVel[1]	= MassRatio[0]*(h_pDSMC->ParticleYVel[SelectParticleNo[0]]) + MassRatio[1]*(h_pDSMC->ParticleYVel[SelectParticleNo[1]]) ;
	CenterMassVel[2]	= MassRatio[0]*(h_pDSMC->ParticleZVel[SelectParticleNo[0]]) + MassRatio[1]*(h_pDSMC->ParticleZVel[SelectParticleNo[1]]) ;


	// use the VHS logic
	if ( fabs(pMultiSpecies->VSSParameter-1.) < 1.e-3 ){
		// B is the cosine of a random elevation angle
		B	= 2.*Randn() - 1. ;
		
		Angle	= sqrt(1. - (B*B)) ;

		RelativeVelPost[0]	= B * (*pRelativeSpeed) ;

		// AzimuthAngle is a random azimuth angle
		AzimuthAngle		= 2.*PI*Randn() ;

		RelativeVelPost[1]	= Angle * cos(AzimuthAngle) * (*pRelativeSpeed) ;
		RelativeVelPost[2]	= Angle * sin(AzimuthAngle) * (*pRelativeSpeed) ;
		
	// use the VSS logic
	}else{
		// B is the cosine of the deflection angle for the VSS model, eqn (11.8)
		B	= 2.*(pow(Randn() , pMultiSpecies->VSSParameter)) - 1. ;

		Angle	= sqrt(1. - (B*B)) ;

		AzimuthAngle	= 2.*PI*Randn() ;
		CosAzimuthAngle	= cos(AzimuthAngle) ;
		SinAzimuthAngle	= sin(AzimuthAngle) ;
		
		D	= sqrt(RelativeVel[1]*RelativeVel[1] + RelativeVel[2]*RelativeVel[2]) ;
		
		// the post-collision relative velocity components are based on eqn (2.22)
		// RelativeVelPost[3] are the components of the post-collision relative vel.
		if ( D > 1.E-6 ){
			RelativeVelPost[0]	= B*RelativeVel[0] + Angle*SinAzimuthAngle*D ;
			RelativeVelPost[1]	= B*RelativeVel[1] + Angle*((*pRelativeSpeed)*RelativeVel[2]*CosAzimuthAngle-RelativeVel[0]*RelativeVel[1]*SinAzimuthAngle)/D ;
			RelativeVelPost[2]	= B*RelativeVel[2] - Angle*((*pRelativeSpeed)*RelativeVel[1]*CosAzimuthAngle+RelativeVel[0]*RelativeVel[2]*SinAzimuthAngle)/D ;
		}else{
			RelativeVelPost[0]	= B * RelativeVel[0] ;
			RelativeVelPost[0]	= Angle * CosAzimuthAngle * RelativeVel[0] ;
			RelativeVelPost[0]	= Angle * SinAzimuthAngle * RelativeVel[0] ;
		}
	}// End if ( fabs(pMultiSpecies->VSSParameter-1.) < 1.e-3 )


	h_pDSMC->ParticleXVel[SelectParticleNo[0]]	= CenterMassVel[0] + RelativeVelPost[0]*MassRatio[1] ;
	h_pDSMC->ParticleYVel[SelectParticleNo[0]]	= CenterMassVel[1] + RelativeVelPost[1]*MassRatio[1] ;
	h_pDSMC->ParticleZVel[SelectParticleNo[0]]	= CenterMassVel[2] + RelativeVelPost[2]*MassRatio[1] ;

	h_pDSMC->ParticleXVel[SelectParticleNo[1]]	= CenterMassVel[0] - RelativeVelPost[0]*MassRatio[0] ;
	h_pDSMC->ParticleYVel[SelectParticleNo[1]]	= CenterMassVel[1] - RelativeVelPost[1]*MassRatio[0] ;
	h_pDSMC->ParticleZVel[SelectParticleNo[1]]	= CenterMassVel[2] - RelativeVelPost[2]*MassRatio[0] ;
}

//==============================================================================================================
//==============================================================================================================

void Sample(	DSMC_DOMAIN		*h_pDomain ,
		DSMC_SPECIES		*h_Species ,
		DSMC_DSMC		*h_pDSMC ,
		DSMC_CELL		*h_Cell ){
			
	int		GroupNum , CellNum , ParticleNum , SumParticleNum ;
	int		ParticleNo , SpeciesNo ;
	
	// Debug.
	//double		SampleVelSq[3] , SampleVel[3] ;
	
	
	GroupNum	= h_pDomain->SpeciesGroupNum ;
	CellNum		= h_pDomain->CellNum ;
	
	h_pDomain->SamplingNum++ ;
	for ( int i=0 ; i<GroupNum ; i++ ){
		for ( int j=0 ; j<CellNum ; j++ ){
			SumParticleNum	= h_pDSMC->IndexCell1[i][j] ;
			ParticleNum	= h_pDSMC->IndexCell2[i][j] ;
			
			// Debug.
			/*SampleVelSq[0]	= 0. ;
			SampleVelSq[1]	= 0. ;
			SampleVelSq[2]	= 0. ;
			SampleVel[0]	= 0. ;
			SampleVel[1]	= 0. ;
			SampleVel[2]	= 0. ;*/
			
			for ( int z=0 ; z<ParticleNum ; z++ ){
				ParticleNo	= h_pDSMC->IndexParticle[SumParticleNum + z] ;

				SpeciesNo	= h_pDSMC->ParticleSpeciesNo[ParticleNo] ;

				h_pDSMC->SampleParticleNum[SpeciesNo][j]++ ;
				h_pDSMC->SampleXVel[SpeciesNo][j]	+= h_pDSMC->ParticleXVel[ParticleNo] ;
				h_pDSMC->SampleYVel[SpeciesNo][j]	+= h_pDSMC->ParticleYVel[ParticleNo] ;
				h_pDSMC->SampleZVel[SpeciesNo][j]	+= h_pDSMC->ParticleZVel[ParticleNo] ;
				h_pDSMC->SampleXVelSq[SpeciesNo][j]	+= (h_pDSMC->ParticleXVel[ParticleNo] * h_pDSMC->ParticleXVel[ParticleNo]) ;
				h_pDSMC->SampleYVelSq[SpeciesNo][j]	+= (h_pDSMC->ParticleYVel[ParticleNo] * h_pDSMC->ParticleYVel[ParticleNo]) ;
				h_pDSMC->SampleZVelSq[SpeciesNo][j]	+= (h_pDSMC->ParticleZVel[ParticleNo] * h_pDSMC->ParticleZVel[ParticleNo]) ;
				
				// Debug.
				/*SampleVel[0]	+= h_pDSMC->ParticleXVel[ParticleNo] ;
				SampleVel[1]	+= h_pDSMC->ParticleYVel[ParticleNo] ;
				SampleVel[2]	+= h_pDSMC->ParticleZVel[ParticleNo] ;
				SampleVelSq[0]	+= (h_pDSMC->ParticleXVel[ParticleNo] * h_pDSMC->ParticleXVel[ParticleNo]) ;
				SampleVelSq[1]	+= (h_pDSMC->ParticleYVel[ParticleNo] * h_pDSMC->ParticleYVel[ParticleNo]) ;
				SampleVelSq[2]	+= (h_pDSMC->ParticleZVel[ParticleNo] * h_pDSMC->ParticleZVel[ParticleNo]) ;*/
				
				if ( h_Species[SpeciesNo].RotDOF > 0.01 ) 
					h_pDSMC->SampleRotation[SpeciesNo][j]	+= h_pDSMC->ParticleRotation[ParticleNo] ;
					
				if ( h_Species[SpeciesNo].VibMode > 0.01 ) 
					h_pDSMC->SampleVibration[SpeciesNo][j]	+= h_pDSMC->ParticleVibration[ParticleNo] ;
					
				if ( h_Species[SpeciesNo].VibMode == 0. ) 
					continue ;
										
				if ( h_pDomain->VibrationalModel == 2 ){
					
					if ( h_pDSMC->ParticleVibLevel[ParticleNo] == 0 ) h_pDSMC->SampleVibGroundNum[SpeciesNo][j]++ ;
					if ( h_pDSMC->ParticleVibLevel[ParticleNo] == 1 ) h_pDSMC->SampleVibLevelNum[SpeciesNo][j]++ ;
						
					h_pDSMC->SampleCellVibLevel[j]	+= h_pDSMC->ParticleVibLevel[ParticleNo] ;						
					h_pDSMC->SampleVibLevel[SpeciesNo][j]	+= h_pDSMC->ParticleVibLevel[ParticleNo] ;
				}
			}// End particle number
			
			// Debug.
			/*if ( h_Cell[j].Id == 1342 ){
				cout << "Sample-1: " << SampleVel[0] << ", " << SampleVel[1] << ", " << SampleVel[2] << ", " << ParticleNum << '\n' ;
				cout << "Sample-2: " << SampleVelSq[0] << ", " << SampleVelSq[1] << ", " << SampleVelSq[2] << ", " << ParticleNum << '\n' ;
			}*/
		}// End cell number
	}// End species group number
}

//==============================================================================================================

void SampleInit( DSMC_DOMAIN *h_pDomain , DSMC_DSMC *h_pDSMC , DSMC_CELL *h_Cell  ){
	int		CellNum , WallFaceNum , SpeciesNum ;
	
	CellNum		= h_pDomain->CellNum ;
	WallFaceNum	= h_pDomain->WallFaceNum ;
	SpeciesNum	= h_pDomain->SpeciesNum ;
	
	
	h_pDomain->SamplingNum		= 0 ;
	h_pDomain->EnterParticleNum	= 0 ;
	
	
	for ( int i=0 ; i<CellNum ; i++ ){
		
		h_pDSMC->SamplingTimeInit[i]	= h_pDSMC->SamplingTimeEnd[i] ;
		
		h_Cell[i].CollDistance = 0. ;
		h_Cell[i].CollNum = 0. ;	
		h_Cell[i].SelectNum = 0. ;
		
	}
	
	// Fluid field
	for ( int i=0 ; i<SpeciesNum ; i++ ){
		for ( int j=0 ; j<CellNum ; j++ ){
			h_pDSMC->SampleParticleNum[i][j]	= 1.E-6 ;
			h_pDSMC->SampleXVel[i][j]		= 0. ;
			h_pDSMC->SampleYVel[i][j]		= 0. ;
			h_pDSMC->SampleZVel[i][j]		= 0. ;
			h_pDSMC->SampleXVelSq[i][j]		= 0. ;
			h_pDSMC->SampleYVelSq[i][j]		= 0. ;
			h_pDSMC->SampleZVelSq[i][j]		= 0. ;
			h_pDSMC->SampleRotation[i][j]		= 0. ;
			h_pDSMC->SampleVibration[i][j]		= 0. ;
			h_pDSMC->SampleVibGroundNum[i][j]	= 0. ;
			h_pDSMC->SampleVibLevelNum[i][j]	= 0. ;
			h_pDSMC->SampleVibLevel [i][j]	= 0. ;
		}	
	}
	
	for ( int i=0 ; i<CellNum ; i++ )
		h_pDSMC->SampleCellVibLevel[i]	= 0. ;
	
	
	// Surface
	for ( int i=0 ; i<WallFaceNum ; i++ ){
		for ( int j=0 ; j<SpeciesNum ; j++ ){
			h_pDSMC->SampleSurfaceParticleNum[i][j]		= 1.E-6 ;
			h_pDSMC->SampleSurfaceStickingParticleNum[i][j]	= 0. ;
			h_pDSMC->SampleSurfaceInNormMomentum[i][j]	= 0. ;
			h_pDSMC->SampleSurfaceReNormMomentum[i][j]	= 0. ;
			h_pDSMC->SampleSurfaceInXMomentum[i][j]		= 0. ;
			h_pDSMC->SampleSurfaceReXMomentum[i][j]		= 0. ;
			h_pDSMC->SampleSurfaceInYMomentum[i][j]		= 0. ;
			h_pDSMC->SampleSurfaceReYMomentum[i][j]		= 0. ;
			h_pDSMC->SampleSurfaceInZMomentum[i][j]		= 0. ;
			h_pDSMC->SampleSurfaceReZMomentum[i][j]		= 0. ;
			h_pDSMC->SampleSurfaceInTransEng[i][j]		= 0. ;
			h_pDSMC->SampleSurfaceReTransEng[i][j]		= 0. ;
			h_pDSMC->SampleSurfaceInRotEng[i][j]		= 0. ;
			h_pDSMC->SampleSurfaceReRotEng[i][j]		= 0. ;
			h_pDSMC->SampleSurfaceInVibEng[i][j]		= 0. ;
			h_pDSMC->SampleSurfaceReVibEng[i][j]		= 0. ;
		}	
	}
}

//==============================================================================================================

void SamplePhysicalTime( DSMC_DOMAIN *h_pDomain , DSMC_CELL *h_Cell , DSMC_DSMC *h_pDSMC ){
	int		CellNum ;
	
	CellNum		= h_pDomain->CellNum ;
	
	for ( int i=0 ; i<CellNum ; i++ ){
		h_pDSMC->SamplingTimeEnd[i]	+= h_Cell[i].Timestep ;	
	}
}

//==============================================================================================================
//==============================================================================================================

void CalculateResult(	DSMC_DOMAIN		*h_pDomain ,
			DSMC_CELL		*h_Cell ,
			DSMC_SPECIES		*h_Species ,
			DSMC_DSMC		*h_pDSMC ,
			DSMC_RESULT		*h_pResult ){
	
	int		MPISize , MPIMyID ;
	int		CellNum , SamplingNum , SpeciesNum ;
	double		Weight , SumParticleNum , SumTotalParticleNum , SumMass , SumMassVel[3] , SumMassVelSq ;
	double		SumRotDOF , SumRotEnergy , SumVibDOF , SumVibEnergy , SumVibTemp ;
	double		NumDensity , Density , Vel[3] , VelSq , TransTemp , TransXTemp , TransYTemp , TransZTemp ;
	double		RotTemp , VibTemp , TotalTemp , AveParticleNum , MeanFreePath , MeanCollSpacingMeanFreePath ;
	double 		SumParticleRot , SumParticleVib, buffer1,buffer2, A ,B ,MeanCollisionTime ,SumVelSq	;
	
	
	MPI_Comm_size( MPI_COMM_WORLD , &MPISize ) ;
	MPI_Comm_rank( MPI_COMM_WORLD , &MPIMyID ) ;
	
/*	
	// Validation for vibration energy with Bird's DSMC0V.FOR code (Debug).
	ofstream	OutputTemp1 ;
	ofstream	OutputTemp2 ;
	ofstream	OutputTemp3 ;
	
	if ( MPIMyID == 0 )
		OutputTemp1.open( "Debug-Vibration-Species1.dat" , ios::out | ios::app ) ;
		OutputTemp2.open( "Debug-Vibration-Species2.dat" , ios::out | ios::app ) ;
		OutputTemp3.open( "Debug-Vibration.dat" , ios::out | ios::app ) ;	
*/	

	CellNum		= h_pDomain->CellNum ;
	SamplingNum	= h_pDomain->SamplingNum ;
	SpeciesNum	= h_pDomain->SpeciesNum ;
	
	
	for ( int i=0 ; i<CellNum ; i++ ){
		Weight	= h_Cell[i].Weighting/(h_Cell[i].Volume*SamplingNum) ;
		
		SumParticleNum		= 0. ;
		SumTotalParticleNum	= 0. ;
		SumMass			= 0. ;
		SumMassVel[0]		= 0. ;
		SumMassVel[1]		= 0. ;
		SumMassVel[2]		= 0. ;
		SumMassVelSq		= 0. ;
		SumRotDOF		= 0. ;
		SumRotEnergy		= 0. ;
		SumVibDOF		= 0. ;
		SumVibEnergy		= 0. ;
		SumVibTemp		= 0. ;
		SumParticleVib  = 0. ;
		SumVelSq   = 0. ;
		
		for ( int j=0 ; j<SpeciesNum ; j++ ){
			SumTotalParticleNum += h_pDSMC->SampleParticleNum[j][i] ;
			
			if ( h_Species[j].VibMode != 0. )	SumParticleVib+= h_pDSMC->SampleParticleNum[j][i] ;
		
		}
		
		for ( int j=0 ; j<SpeciesNum ; j++ ){
			SumParticleNum	+= h_pDSMC->SampleParticleNum[j][i] ;
			SumMass		+= (h_Species[j].Mass * h_pDSMC->SampleParticleNum[j][i]) ;
			SumMassVel[0]	+= (h_Species[j].Mass * h_pDSMC->SampleXVel[j][i]) ;
			SumMassVel[1]	+= (h_Species[j].Mass * h_pDSMC->SampleYVel[j][i]) ;
			SumMassVel[2]	+= (h_Species[j].Mass * h_pDSMC->SampleZVel[j][i]) ;
			SumMassVelSq	+= (h_Species[j].Mass * (h_pDSMC->SampleXVelSq[j][i]+h_pDSMC->SampleYVelSq[j][i]+h_pDSMC->SampleZVelSq[j][i])) ;
			SumVelSq	+= (h_pDSMC->SampleXVelSq[j][i]+h_pDSMC->SampleYVelSq[j][i]+h_pDSMC->SampleZVelSq[j][i]) ;
			
			
			SumRotEnergy	+= h_pDSMC->SampleRotation[j][i] ;
			SumRotDOF	+= (h_Species[j].RotDOF * h_pDSMC->SampleParticleNum[j][i]) ;
			
			
			
			if ( h_pDomain->VibrationalModel == 2 ){
				
				if ( h_pDSMC->SampleVibLevel[j][i] > 0.01 ){
					
//					if ( h_pDSMC->SampleVibGroundNum[j][i] > h_pDSMC->SampleVibLevelNum[j][i] ){						
//						VibTemp	= h_Species[j].VibTemp/log(h_pDSMC->SampleVibGroundNum[j][i]/h_pDSMC->SampleVibLevelNum[j][i]) ;					
//					}else{				
						// VibTemp the vib. temp. of species j is calculated from eqn (11.27 & 11.33)
//						VibTemp	= h_Species[j].VibTemp/log(1.+1./(h_pDSMC->SampleCellVibLevel[i]/SumTotalParticleNum)) ;	
						// VibTemp the vib. temp. of species j is calculated from new eqn (4.45)
						VibTemp	= h_Species[j].VibTemp/log(1.+1./(h_pDSMC->SampleVibLevel[j][i]/h_pDSMC->SampleParticleNum[j][i])) ;	
//					}
						//EffVibDOF is the effective number of vibrational degrees of freedom calculated from eqn (11.28)
//						h_pDSMC->EffVibDOF[i][j]	= 2.*h_pDSMC->SampleVibration[j][i]/(h_pDSMC->SampleParticleNum[j][i]*BOLTZ*VibTemp) ;	
						//EffVibDOF is the effective number of vibrational degrees of freedom calculated from new eqn (4.46)					
						h_pDSMC->EffVibDOF[i][j]	= 2.*(h_pDSMC->SampleVibLevel[j][i]/h_pDSMC->SampleParticleNum[j][i])*log(1.+1./(h_pDSMC->SampleVibLevel[j][i]/h_pDSMC->SampleParticleNum[j][i]));					
//						SumVibTemp			+= (VibTemp*h_pDSMC->EffVibDOF[i][j]*h_pDSMC->SampleParticleNum[j][i]) ;
						SumVibTemp			+= (VibTemp*h_pDSMC->SampleParticleNum[j][i]/SumParticleVib) ;	
								
				}else{
					
					VibTemp			= 1.E-6 ;
					h_pDSMC->EffVibDOF[i][j]= 0. ;
					
				}
				
				h_pResult->VibTempSpecies[j][i]	= VibTemp ;
			}
			
			SumVibEnergy	+= h_pDSMC->SampleVibration[j][i] ;
			SumVibDOF	+= (h_pDSMC->EffVibDOF[i][j] * h_pDSMC->SampleParticleNum[j][i]) ;
			
		}// End species number
		
		// NumDensity is the number density, see eqn (1.34)
		// Density is the density, see eqn (1.42)
		// Vel[3] are the stream velocity components, see eqn (1.43)
		NumDensity	= SumParticleNum * Weight ;
		Density		= NumDensity*SumMass/SumParticleNum ;
		Vel[0]		= SumMassVel[0]/SumMass ;
		Vel[1]		= SumMassVel[1]/SumMass ;
		Vel[2]		= SumMassVel[2]/SumMass ;
		
		// TransTemp is the translational temperature, see eqn (1.51)
		VelSq		= Vel[0]*Vel[0] + Vel[1]*Vel[1] + Vel[2]*Vel[2] ;
		TransTemp	= (SumMassVelSq - SumMass*VelSq)/(3.*BOLTZ*SumParticleNum) ;
		//TransTemp	= (SumMassVelSq - SumMass*VelSq)/(3.*1.3806504e-23*SumParticleNum) ;
		
		// RotTemp is the rotational temperature, see eqn (11.11)
		if ( SumRotDOF > 1.E-6 ){
			RotTemp	= (2./BOLTZ)*SumRotEnergy/SumRotDOF ;
		}else{
			RotTemp	= 0. ;
		}
		
		// VibTemp is the vibrational temperature, see eqn (11.11)	
		if ( h_pDomain->VibrationalModel == 1 ){
			if ( SumVibDOF > 1.E-6 )
				VibTemp	= (2./BOLTZ)*SumVibEnergy/SumVibDOF ;
			else
				VibTemp	= 0. ;
		}else if ( h_pDomain->VibrationalModel == 2 ){
			if ( SumVibDOF > 1.E-6 )

//				VibTemp	= SumVibTemp/SumVibDOF  ;
					VibTemp	= SumVibTemp  ;

			else
				VibTemp	= 0. ;
		}else{
			VibTemp	= 0. ;
		}
/*		
		// TotalTemp is the overall temperature, see eqn (11.12)
		TotalTemp	= (3.*TransTemp + (SumRotDOF/SumParticleNum)*RotTemp + (SumVibDOF/SumParticleNum)*VibTemp)/
				  (3.+(SumRotDOF+SumVibDOF)/SumParticleNum) ;
		 
		
		MeanFreePath	= pow((TotalTemp/h_Species[0].RefTemp) , (h_Species[0].VisTempIndex-0.5f)) / (4.4429*h_Species[0].Diameter*h_Species[0].Diameter*NumDensity) ;
		
		
		// Calculate the ratio of the mean collision spacing (the mean distance between two particles selected for collision) and the mean free path (mcs/mfp)
		if ( h_Cell[i].CollNum < 1. ){
			MeanCollSpacingMeanFreePath	= 0. ;
		}else{
			if ( MeanFreePath == 0. )
				MeanCollSpacingMeanFreePath	= 0. ;
			else
				MeanCollSpacingMeanFreePath	= (h_Cell[i].CollDistance/h_Cell[i].CollNum) / MeanFreePath ;
		}
		
	
		if ( h_Cell[i].SelectNum == 0. )
			h_Cell[i].CollRatio	= 0. ;
		else 
			h_Cell[i].CollRatio	= h_Cell[i].CollNum/h_Cell[i].SelectNum ;
*/		
		
		// Debug.
		/*if ( h_Cell[i].Id == 9209 ){
			cout << "Trans: " << TransTemp << setw(20) << h_pDSMC->SampleXVelSq[0][i] << setw(20) << h_pDSMC->SampleYVelSq[0][i] << setw(20) << h_pDSMC->SampleZVelSq[0][i] << '\n' ;
			cout << setw(26) << h_pDSMC->SampleXVelSq[0][i]/SamplingNum << setw(20) << h_pDSMC->SampleYVelSq[0][i]/SamplingNum << setw(20) << h_pDSMC->SampleZVelSq[0][i]/SamplingNum << '\n' ;
			cout << setw(30) << SumMass << setw(20) << SumParticleNum << setw(20) << h_pDSMC->SampleXVel[0][i] << setw(20) << h_pDSMC->SampleYVel[0][i] << setw(20) << h_pDSMC->SampleZVel[0][i] << '\n' ;
			cout << setw(30) << Vel[0] << setw(20) << Vel[1] << setw(20) << Vel[2] << '\n' ;
			cout << "mass: " << h_Species[0].Mass << ", " << SumMass << ", " << h_pDSMC->SampleParticleNum[0][i] << ", " << SumParticleNum << ", " << SamplingNum << '\n' ;
			cout << setw(30) << SumMassVelSq/(3.*1.3806504e-23*SumParticleNum) << setw(30) << SumMass*VelSq/(3.*1.3806504e-23*SumParticleNum) << setw(30) << ( SumMassVelSq/(3.*1.3806504e-23*SumParticleNum) - SumMass*VelSq/(3.*1.3806504e-23*SumParticleNum) ) << '\n' ;
		}*/
	
		if ( h_Cell[i].CollNum < 1. ){
			
			MeanCollSpacingMeanFreePath	= 0. ;
			MeanCollisionTime = 0. ;
			MeanFreePath = 0. ;
			
		}else{
			
			MeanCollisionTime = 0.5 * SumParticleNum * h_Cell[i].Timestep / (h_Cell[i].CollNum / h_pDomain->SamplingFrequency) ;
			MeanFreePath = 0.92132 * sqrt ( fabs( SumVelSq / SumParticleNum - VelSq ))* MeanCollisionTime;
//			MeanFreePath = 0.92132 * sqrt( Vel[0]*Vel[0] + Vel[1]*Vel[1] + Vel[2]*Vel[2] )* MeanCollisionTime;
			
			
			if ( MeanFreePath == 0. )
				MeanCollSpacingMeanFreePath	= 0. ;
			else
				MeanCollSpacingMeanFreePath	= (h_Cell[i].CollDistance/h_Cell[i].CollNum) / MeanFreePath ;
		}
		
		if ( h_Cell[i].SelectNum == 0. )
			h_Cell[i].CollRatio	= 0. ;
		else 
			h_Cell[i].CollRatio	= h_Cell[i].CollNum/h_Cell[i].SelectNum ;
		
		h_pResult->NumDensity[i]	= NumDensity ;
		h_pResult->Density[i]		= Density ;
		h_pResult->XVel[i]		= Vel[0] ;
		h_pResult->YVel[i]		= Vel[1] ;
		h_pResult->ZVel[i]		= Vel[2] ;
//		h_pResult->Temp[i]		= TotalTemp ;
		h_pResult->TransTemp[i]		= TransTemp ;
		h_pResult->RotTemp[i]		= RotTemp ;
		h_pResult->VibTemp[i]		= VibTemp ;
		h_pResult->AveParticleNum[i]	= SumParticleNum/SamplingNum ;
		h_pResult->MeanCollSpacingMeanFreePath[i] = MeanCollSpacingMeanFreePath ;
		
		h_Cell[i].MeanFreePath		= MeanFreePath ;
//		h_Cell[i].MeanCollisionTime	= MeanFreePath/sqrt(VelSq) ;
//		h_Cell[i].MeanCollisionTime	= MeanFreePath/sqrt(VelSq) ;
		h_Cell[i].MeanCollisionTime	= MeanCollisionTime ;
//		h_Cell[i].Temp			= TotalTemp ;
		h_Cell[i].Speed			= sqrt( Vel[0]*Vel[0] + Vel[1]*Vel[1] + Vel[2]*Vel[2] ) ;
		h_Cell[i].AveParticleNum	= h_pResult->AveParticleNum[i] ;
		h_Cell[i].MeanCollSpacingMeanFreePath = h_pResult->MeanCollSpacingMeanFreePath[i] ;

	}// End cell number
	


	// Calculate macroscopic properties for each species.
	for ( int i=0 ; i<SpeciesNum ; i++ ){
		for ( int j=0 ; j<CellNum ; j++ ){
			Weight		= h_Cell[j].Weighting/(h_Cell[j].Volume*SamplingNum) ;
			
			
			// NumDensity is the partial number density,only for one species
			// Density is the partial density, see eqn (1.13)
			// Vel[3] defines the average velocity of the species L molecules
			NumDensity	= h_pDSMC->SampleParticleNum[i][j]*Weight ;
			Density		= h_Species[i].Mass*NumDensity ;
			Vel[0]		= h_pDSMC->SampleXVel[i][j]/h_pDSMC->SampleParticleNum[i][j] ;
			Vel[1]		= h_pDSMC->SampleYVel[i][j]/h_pDSMC->SampleParticleNum[i][j] ;
			Vel[2]		= h_pDSMC->SampleZVel[i][j]/h_pDSMC->SampleParticleNum[i][j] ;
			
			// the component temperatures are based on eqn (1.30)
			VelSq		= Vel[0]*Vel[0] + Vel[1]*Vel[1] + Vel[2]*Vel[2] ;
			TransXTemp	= (h_Species[i].Mass/BOLTZ)*(h_pDSMC->SampleXVelSq[i][j]/h_pDSMC->SampleParticleNum[i][j]-Vel[0]*Vel[0]) ;
			TransYTemp	= (h_Species[i].Mass/BOLTZ)*(h_pDSMC->SampleYVelSq[i][j]/h_pDSMC->SampleParticleNum[i][j]-Vel[1]*Vel[1]) ;
			TransZTemp	= (h_Species[i].Mass/BOLTZ)*(h_pDSMC->SampleZVelSq[i][j]/h_pDSMC->SampleParticleNum[i][j]-Vel[2]*Vel[2]) ;
			
			// TransTemp is the translational temperature, see eqn (1.29)
			TransTemp	= (h_Species[i].Mass/(3.*BOLTZ))*
						((h_pDSMC->SampleXVelSq[i][j]+h_pDSMC->SampleYVelSq[i][j]+h_pDSMC->SampleZVelSq[i][j])/h_pDSMC->SampleParticleNum[i][j]-VelSq) ;
        		
        		// RotTemp is the rotational temperature, see eqn (11.10)
       			if ( h_Species[i].RotDOF > 0.01 ){
        			RotTemp	= 2.*h_pDSMC->SampleRotation[i][j]/(h_Species[i].RotDOF*BOLTZ*h_pDSMC->SampleParticleNum[i][j]) ;
        		}else{
        			RotTemp	= 0. ;	
        		}
        		
        		// VibTemp is the vibrational temperature, see eqn (11.10)
        		if ( h_pDomain->VibrationalModel == 1 ){
        			if ( h_pDSMC->EffVibDOF[j][i] > 0.01 ){
        				VibTemp	= 2.*h_pDSMC->SampleVibration[i][j]/(h_pDSMC->EffVibDOF[j][i]*BOLTZ*h_pDSMC->SampleParticleNum[i][j]) ;
        			}else{
        				VibTemp	= 0. ;
        			}
        			
        			h_pResult->VibTempSpecies[i][j]	= VibTemp ;
        		}
        		
        		TotalTemp	= (3.*TransTemp + h_Species[i].RotDOF*RotTemp + h_pDSMC->EffVibDOF[j][i]*h_pResult->VibTempSpecies[i][j])/
        					(3. + h_Species[i].RotDOF + h_pDSMC->EffVibDOF[j][i]) ;
        		
        		
      h_pResult->NumDensitySpecies[i][j]	= NumDensity ;
			h_pResult->DensitySpecies[i][j]		= Density ;
			h_pResult->XVelSpecies[i][j]		= Vel[0] ;
			h_pResult->YVelSpecies[i][j]		= Vel[1] ;
			h_pResult->ZVelSpecies[i][j]		= Vel[2] ;
			h_pResult->TempSpecies[i][j]		= TotalTemp ;
			h_pResult->TransTempSpecies[i][j]	= TransTemp ;
			h_pResult->TransXTempSpecies[i][j]	= TransXTemp ;
			h_pResult->TransYTempSpecies[i][j]	= TransYTemp ;
			h_pResult->TransZTempSpecies[i][j]	= TransZTemp ;
			h_pResult->RotTempSpecies[i][j]		= RotTemp ;
			h_pResult->AveParticleNumSpecies[i][j]	= h_pDSMC->SampleParticleNum[i][j]/SamplingNum ;
			
/*			
			// Validation for vibration energy with Bird's DSMC0V.FOR code (Debug).
//			if ( h_pDomain->TimestepNo <= h_pDomain->SamplingTime ){
			if ( h_pDomain->TimestepNo <= 50000 ){
				
				double		CollisionRate ;				
				CollisionRate	= 0. ;
				
				for ( int k=0 ; k<SpeciesNum ; k++ )
				
					CollisionRate	+= h_pDSMC->CollisionNum[i][k]*SamplingNum/h_pDSMC->SampleParticleNum[i][j] ; 

				if ( MPIMyID == 0 ){
					
					if (i==0 && ( h_pDomain->TimestepNo <=10 || ( h_pDomain->TimestepNo % 20) == 0))
					OutputTemp1 << setw(20) << h_pDomain->TimestepNo << setw(20) << h_pDSMC->CollisionNum[i][i]*SamplingNum/h_pDSMC->SampleParticleNum[i][j]<< setw(20) << h_pResult->TransTempSpecies[i][j] 
							<< setw(20) << h_pResult->RotTempSpecies[i][j] << setw(20) << h_pResult->VibTempSpecies[i][j]	 
							<< setw(20) << h_pResult->TempSpecies[i][j] <<  setw(20) <<  h_pResult->NumDensitySpecies[i][j]/h_pDomain->NumDen
							<<  setw(20) << h_pDSMC->EffVibDOF[j][i]<<'\n' ;
							
					if (i==1 && ( h_pDomain->TimestepNo <=10 || ( h_pDomain->TimestepNo % 20) == 0)	)	
					OutputTemp2 << setw(20) << h_pDomain->TimestepNo<< setw(20) << h_pDSMC->CollisionNum[i][i]*SamplingNum/h_pDSMC->SampleParticleNum[i][j] << setw(20) << h_pResult->TransTempSpecies[i][j] 
							<< setw(20) << h_pResult->RotTempSpecies[i][j] << setw(20) << h_pResult->VibTempSpecies[i][j]	 
							<< setw(20) << h_pResult->TempSpecies[i][j] <<  setw(20) <<  h_pResult->NumDensitySpecies[i][j]/h_pDomain->NumDen
							<<  setw(20) << h_pDSMC->EffVibDOF[j][i]<<'\n' ;

					if ( i==4 && ( h_pDomain->TimestepNo <=10 ||( h_pDomain->TimestepNo<=100 && ( h_pDomain->TimestepNo % 10) == 0) ||	
						
						( h_pDomain->TimestepNo<=1000 && ( h_pDomain->TimestepNo % 100) == 0)||( h_pDomain->TimestepNo<=10000 && ( h_pDomain->TimestepNo % 1000) == 0)))
					
					OutputTemp3 << setw(20) << h_pDomain->TimestepNo << setw(20) << h_pDSMC->EffVibDOF[j][0] << setw(20) << h_Cell[j].Temp 
							<<  setw(20) <<  h_pResult->NumDensitySpecies[0][j]/h_pDomain->NumDen<<  setw(20) <<  h_pResult->NumDensitySpecies[1][j]/h_pDomain->NumDen 
							<<  setw(20) <<  h_pResult->NumDensitySpecies[2][j]/h_pDomain->NumDen<<  setw(20) <<  h_pResult->NumDensitySpecies[3][j]/h_pDomain->NumDen 
							<<  setw(20) <<  h_pResult->NumDensitySpecies[4][j]/h_pDomain->NumDen
							<<'\n' ;	
											
				}
			}
*/
			
		}// End cell number
	}//End species number
  

	
	for ( int i=0 ; i< CellNum ; i++ ){
		
		  A = 0.;
  		B = 0.;
  		
		for ( int j=0 ; j< SpeciesNum; j++ ){
			
			A = A + ( 3. + h_Species[j].RotDOF + h_pDSMC->EffVibDOF[i][j] ) * h_pDSMC->SampleParticleNum[j][i] * h_pResult->TempSpecies[j][i];		
			B = B +( 3. + h_Species[j].RotDOF + h_pDSMC->EffVibDOF[i][j] ) * h_pDSMC->SampleParticleNum[j][i];
			 
		}
			 
			// TotalTemp is the overall temperature, see eqn (11.12)
		TotalTemp	= A / B ;
		 
//		MeanCollisionTime = 0.5 * h_Cell[i].AveParticleNum * h_Cell[i].Timestep / (h_Cell[i].CollNum/h_pDomain->TimestepNo) ;
//		MeanFreePath = 0.92132 * sqrt ( fabs((SumVelSq/SumTotalParticleNum - VelSq )* MeanCollisionTime)) ;
		
//		MeanFreePath	= pow((TotalTemp/h_Species[0].RefTemp) , (h_Species[0].VisTempIndex-0.5f)) / (4.4429*h_Species[0].Diameter*h_Species[0].Diameter*h_pResult->NumDensity[i]) ;
		
/*		
		// Calculate the ratio of the mean collision spacing (the mean distance between two particles selected for collision) and the mean free path (mcs/mfp)
		if ( h_Cell[i].CollNum < 1. ){
			
			MeanCollSpacingMeanFreePath	= 0. ;
			MeanCollisionTime = 0. ;
			MeanFreePath = 0. ;
			
		}else{
			
			MeanCollisionTime = 0.5 * SumTotalParticleNum * h_Cell[i].Timestep / (h_Cell[i].CollNum) ;
			MeanFreePath = 0.92132 * sqrt ( fabs( SumVelSq / SumTotalParticleNum - VelSq ))* MeanCollisionTime;
			
			if ( MeanFreePath == 0. )
				MeanCollSpacingMeanFreePath	= -1. ;
			else
				MeanCollSpacingMeanFreePath	= (h_Cell[i].CollDistance/h_Cell[i].CollNum) / MeanFreePath ;
		}
		
		
		if ( h_Cell[i].SelectNum == 0. )
			h_Cell[i].CollRatio	= 0. ;
		else 
			h_Cell[i].CollRatio	= h_Cell[i].CollNum/h_Cell[i].SelectNum ;
		
*/		
		// Debug.
		/*if ( h_Cell[i].Id == 9209 ){
			cout << "Trans: " << TransTemp << setw(20) << h_pDSMC->SampleXVelSq[0][i] << setw(20) << h_pDSMC->SampleYVelSq[0][i] << setw(20) << h_pDSMC->SampleZVelSq[0][i] << '\n' ;
			cout << setw(26) << h_pDSMC->SampleXVelSq[0][i]/SamplingNum << setw(20) << h_pDSMC->SampleYVelSq[0][i]/SamplingNum << setw(20) << h_pDSMC->SampleZVelSq[0][i]/SamplingNum << '\n' ;
			cout << setw(30) << SumMass << setw(20) << SumParticleNum << setw(20) << h_pDSMC->SampleXVel[0][i] << setw(20) << h_pDSMC->SampleYVel[0][i] << setw(20) << h_pDSMC->SampleZVel[0][i] << '\n' ;
			cout << setw(30) << Vel[0] << setw(20) << Vel[1] << setw(20) << Vel[2] << '\n' ;
			cout << "mass: " << h_Species[0].Mass << ", " << SumMass << ", " << h_pDSMC->SampleParticleNum[0][i] << ", " << SumParticleNum << ", " << SamplingNum << '\n' ;
			cout << setw(30) << SumMassVelSq/(3.*1.3806504e-23*SumParticleNum) << setw(30) << SumMass*VelSq/(3.*1.3806504e-23*SumParticleNum) << setw(30) << ( SumMassVelSq/(3.*1.3806504e-23*SumParticleNum) - SumMass*VelSq/(3.*1.3806504e-23*SumParticleNum) ) << '\n' ;
		}*/
	
		
//		h_pResult->NumDensity[i]	= NumDensity ;
//		h_pResult->Density[i]		= Density ;
//		h_pResult->XVel[i]		= Vel[0] ;
//		h_pResult->YVel[i]		= Vel[1] ;
//		h_pResult->ZVel[i]		= Vel[2] ;
			h_pResult->Temp[i]		= TotalTemp ;
//		h_pResult->TransTemp[i]		= TransTemp ;
//		h_pResult->RotTemp[i]		= RotTemp ;
//		h_pResult->VibTemp[i]		= VibTemp ;
//		h_pResult->AveParticleNum[i]	= SumParticleNum/SamplingNum ;
//			h_pResult->MeanCollSpacingMeanFreePath[i] = MeanCollSpacingMeanFreePath ;
		
//			h_Cell[i].MeanFreePath		= MeanFreePath ;
//			h_Cell[i].MeanCollisionTime	= MeanFreePath/sqrt(h_pResult->XVel[i]*h_pResult->XVel[i] + h_pResult->YVel[i]*h_pResult->YVel[i]+h_pResult->ZVel[i]*h_pResult->ZVel[i]) ;
//			h_Cell[i].MeanCollisionTime	= MeanCollisionTime ;
			h_Cell[i].Temp			= TotalTemp ;
//		h_Cell[i].Speed			= sqrt( Vel[0]*Vel[0] + Vel[1]*Vel[1] + Vel[2]*Vel[2] ) ;
//		h_Cell[i].AveParticleNum	= h_pResult->AveParticleNum[i] ;
//			h_Cell[i].MeanCollSpacingMeanFreePath = h_pResult->MeanCollSpacingMeanFreePath[i] ;
	}	
/*	
	buffer1=0 ;
	buffer2=0 ;
	
	for ( int i=0 ; i<CellNum ; i++ ){
		for ( int j=0 ; j<SpeciesNum ; j++ ){
			buffer1 += h_pResult->NumDensitySpecies[j][i]*(h_pResult->TempSpecies[j][i])*(3.+ h_Species[j].RotDOF + h_pDSMC->EffVibDOF[i][j]) ;
			buffer2 += h_pResult->NumDensitySpecies[j][i]*(3+h_Species[j].RotDOF+h_pDSMC->EffVibDOF[i][j]);
		}
		
		h_Cell[i].Temp = buffer1/buffer2;
		
		h_pDSMC->EffVibDOF[i][0]= 2*(2256/h_Cell[i].Temp)/(exp(2256/h_Cell[i].Temp)-1); 
		
	 
	}
*/	
	

	
/*	
	// Validation for vibration energy with Bird's DSMC0V.FOR code (Debug).
	if ( MPIMyID == 0 ){
		OutputTemp1.clear() ;
		OutputTemp1.close() ;
		OutputTemp2.clear() ;
		OutputTemp2.close() ;
		OutputTemp3.clear() ;
		OutputTemp3.close() ;
	}
*/
}

//==============================================================================================================
//==============================================================================================================

void CheckConvergence( DSMC_DOMAIN *h_pDomain , DSMC_RESULT *h_pResult , DSMC_CONVERGENCE *h_pConvergence , int TimestepNo , ofstream &OutputConvergence ){
	int		MPISize , MPIMyID , CellNum , TotalCellNum ;
	double		SumParticleNum , SumSpeed , SumTemp , AveParticleNum , AveSpeed , AveTemp ;
	ofstream	OutputStartSampling ;
	
	
	MPI_Comm_size( MPI_COMM_WORLD , &MPISize ) ;
	MPI_Comm_rank( MPI_COMM_WORLD , &MPIMyID ) ;
	
	
	SumParticleNum	= 0. ;
	SumSpeed	= 0. ;
	SumTemp		= 0. ;
	AveParticleNum	= 0. ;
	AveSpeed	= 0. ;
	AveTemp		= 0. ;
	CellNum		= h_pDomain->CellNum ;
	TotalCellNum	= h_pDomain->TotalCellNum ;
	
	
	//if ( TimestepNo <= h_pConvergence->DataNum ) h_pConvergence->DataNum = TimestepNo ;
	
	
	// Sum of particles, speed, and temperature.
	SumParticleNum	= h_pDomain->ParticleNum ;
	for ( int i=0 ; i<CellNum ; i++ ){
		SumSpeed	+= sqrt( h_pResult->XVel[i]*h_pResult->XVel[i] + h_pResult->YVel[i]*h_pResult->YVel[i] + h_pResult->ZVel[i]*h_pResult->ZVel[i] ) ;
		SumTemp		+= h_pResult->Temp[i] ;
	}
	MPI_Allreduce( &SumParticleNum , &AveParticleNum , 1 , MPI_DOUBLE , MPI_SUM , MPI_COMM_WORLD ) ;
	MPI_Allreduce( &SumSpeed , &AveSpeed , 1 , MPI_DOUBLE , MPI_SUM , MPI_COMM_WORLD ) ;
	MPI_Allreduce( &SumTemp , &AveTemp , 1 , MPI_DOUBLE , MPI_SUM , MPI_COMM_WORLD ) ;

	AveSpeed	/= TotalCellNum ;
	AveTemp		/= TotalCellNum ;

	if ( TimestepNo <= h_pConvergence->DataNum ){
		h_pConvergence->Sum( AveParticleNum , AveSpeed , AveTemp , TimestepNo ) ;
		
		if ( TimestepNo == h_pConvergence->DataNum ){
			h_pConvergence->Average() ;
			h_pConvergence->InitSum() ;
		}
	}else{
		h_pConvergence->Sum( AveParticleNum , AveSpeed , AveTemp , h_pConvergence->DataNum ) ;
		h_pConvergence->Average() ;
		
		/*if ( h_pConvergence->Convergence() && TimestepNo <= h_pDomain->SamplingTime ){
			h_pDomain->SamplingTime		= TimestepNo ;
			h_pDomain->TotalTimestep	= TimestepNo + 2000 ;
		}else{
			h_pConvergence->InitSum() ;
		}*/
		
		if ( h_pDomain->Convergence > 0 ){
			if ( h_pConvergence->Convergence() && TimestepNo <= h_pDomain->SamplingTime ){
				h_pDomain->SamplingTime		= TimestepNo ;
				h_pDomain->TotalTimestep	= TimestepNo + 5000 ;
				
				OutputStartSampling.open( "StartSampling.dat" , ios::out | ios::trunc ) ;
					
				OutputStartSampling << setw(30) << "TimestepNo:" << setw(25) << (h_pDomain->SamplingTime*20) << setw(25) << h_pDomain->SamplingTime << '\n' ;
				OutputStartSampling << setw(30) << "TotalTimestep:" << setw(25) << (h_pDomain->TotalTimestep*20) << setw(25) << h_pDomain->TotalTimestep << '\n' ;
					
				OutputStartSampling.clear() ;
				OutputStartSampling.close() ;
			}
		}
		
		h_pConvergence->InitSum() ;
		
		OutputConvergence << setw(15) << h_pDomain->TimestepNo << setw(15) << TimestepNo << setw(15) << h_pConvergence->OldParticleNum 
					<< setw(12) << h_pConvergence->OldSpeed << setw(12) << h_pConvergence->OldTemp << setw(15) << h_pConvergence->DiffParticleNum 
					<< setw(12) << h_pConvergence->DiffSpeed << setw(12) << h_pConvergence->DiffTemp << '\n' ;
	}
}

//==============================================================================================================
//==============================================================================================================


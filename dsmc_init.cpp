#include <mpi.h>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include "dsmc_parameter.h"
#include "dsmc_init.h"
#include "dsmc_toolfunction.h"

using namespace std ;


//==============================================================================================================

void Initialization( 	DSMC_DOMAIN		*h_pDomain ,
			DSMC_NODE		*h_Node ,
			DSMC_CELL		*h_Cell ,
			DSMC_INLET		*h_Inlet ,
			DSMC_WALLTYPE		*h_WallType ,
			DSMC_SPECIES		*h_Species ,
			DSMC_SURFACE		*h_Surface ,
			DSMC_DSMC		*h_pDSMC ,
			CELLMAPPING		*h_pMapping , 
			ofstream		&OutputDebug ){
				
	int		MPISize , MPIMyID ;
	
	MPI_Comm_size( MPI_COMM_WORLD , &MPISize ) ;
	MPI_Comm_rank( MPI_COMM_WORLD , &MPIMyID ) ;
	

	// Initialize the some species array value for vibrational model.
	InitialSpecieslValue( h_pDomain , h_Species , h_pDSMC ) ;


	// To create cell information, including volume, each face data...etc, and then return minimum cell volume.
	CreateCellInformation( h_pDomain , h_Node , h_Cell , h_pMapping ) ;


	// To calculate timestep, particle weighting, initial maximum (cross section)*(relative speed), and remainder collision pairs for each cell.
	CalculateTimestepWeighting( h_pDomain , h_Cell , h_Species , h_pDSMC , h_pMapping , OutputDebug ) ;


	// Debug.
	/*for ( int i=0 ; i<MPISize ; i++ ){
		if ( MPIMyID == i ){
			cout << "Timestep: " << h_pDomain->Timestep << ", Weighting: " << h_pDomain->ParticleWeighting << '\n' ;	
		}
		MPI_Barrier( MPI_COMM_WORLD ) ;
	}*/

	
	// To create inlet face information.
	if ( h_pDomain->InletFaceNum > 0 )
		CreateInletInformation( h_pDomain , h_Node , h_Cell , h_Inlet , h_Species , h_pDSMC , h_pMapping ) ;
	
	
	// Link surfaces and cells.
	LinkCellSurface( h_pDomain , h_Node , h_Cell , h_Surface , h_pMapping ) ;
	
	
	// To create initial particles based on M-B distribution and initial simulation conditions.
	CreateInitialParticle( h_pDomain , h_Node , h_Cell , h_Species , h_pDSMC , h_pMapping , OutputDebug ) ;
}

//==============================================================================================================
//==============================================================================================================

void InitialSpecieslValue( DSMC_DOMAIN *h_pDomain , DSMC_SPECIES *h_Species , DSMC_DSMC *h_pDSMC ){
	
	if ( h_pDomain->VibrationalModel == 0 ){
		for ( int i=0 ; i<h_pDomain->SpeciesNum ; i++ ){
			h_Species[i].VibMode		= 0 ;
			h_Species[i].VibModel		= 0 ;
			h_Species[i].VibDOF		= 0. ;
			
			for ( int j=0 ; j<h_pDomain->SpeciesNum ; j++ ){
				h_Species[i].VibConstant1[j]	= 0. ;
				h_Species[i].VibConstant2[j]	= 0. ;
			}
			
			h_Species[i].VibTemp		= 0. ;
			
		}
	}else{
		for ( int i=0 ; i<h_pDomain->CellNum ; i++ ){
			for ( int j=0 ; j<h_pDomain->SpeciesNum ; j++ )
				h_pDSMC->EffVibDOF[i][j]	= h_Species[j].VibDOF ;
		}
	}
}

//==============================================================================================================
//==============================================================================================================

void CreateCellInformation(	DSMC_DOMAIN		*h_pDomain ,
				DSMC_NODE		*h_Node ,
				DSMC_CELL		*h_Cell ,
				CELLMAPPING		*h_pMapping ){

	int		MPISize , MPIMyID , CellNum , Dimension ;
	double		MinCellVolume , MinCharacteristicLength , MinCenter ;
	
	MPI_Comm_size( MPI_COMM_WORLD , &MPISize ) ;
	MPI_Comm_rank( MPI_COMM_WORLD , &MPIMyID ) ;
	
	CellNum			= h_pDomain->CellNum ;
	Dimension		= h_pDomain->Dimension ;
	MinCellVolume		= 1.E+6 ;
	MinCharacteristicLength	= 1.E+6 ;
	MinCenter		= 1.E+6 ;


	for ( int i=0 ; i<CellNum ; i++ ){
		// Calculate coordinate of cell center .
		CalculateCellCenter( h_Node , &h_Cell[i] , h_pMapping ) ;
		
		
		// Calculate cell volume, each face function (2D and 3D).
		if ( Dimension == 2 ){
			CalculateCellVolume2D( h_Node , &h_Cell[i] , h_pMapping ) ;
			CreateCellFaceFunction2D( h_Node , &h_Cell[i] , h_pMapping ) ;
		}else if ( Dimension ==3 ){
			CalculateCellVolume3D( h_Node , &h_Cell[i] , h_pMapping ) ;
			CreateCellFaceFunction3D( h_Node , &h_Cell[i] , h_pMapping ) ;
		}else if ( Dimension == 4 ){
			CalculateCellVolumeAxisymmetric( h_Node , &h_Cell[i] , h_pMapping ) ;
			CreateCellFaceFunction2D( h_Node , &h_Cell[i] , h_pMapping ) ;
		}
		
		
		// Find out the minimum volume of cell and characteristic length.
		if ( h_Cell[i].Volume < MinCellVolume ) MinCellVolume = h_Cell[i].Volume ;
		if ( h_Cell[i].CharacteristicLength < MinCharacteristicLength ) MinCharacteristicLength = h_Cell[i].CharacteristicLength ;
		if ( h_Cell[i].YCenter < MinCenter ) MinCenter = h_Cell[i].YCenter ;
		// Debug.
		//h_Cell[i].Dump(i) ;
		//getchar() ;
	}
	//h_Cell[0].Dump(0) ;
	//getchar() ;


	// Update the minimum cell volume of cell.
	MPI_Allreduce( &MinCellVolume , &h_pDomain->MinVolume , 1 , MPI_DOUBLE , MPI_MIN , MPI_COMM_WORLD ) ;
	MPI_Allreduce( &MinCharacteristicLength , &h_pDomain->MinLength , 1 , MPI_DOUBLE , MPI_MIN , MPI_COMM_WORLD ) ;
	MPI_Allreduce( &MinCenter , &h_pDomain->MinCenter , 1 , MPI_DOUBLE , MPI_MIN , MPI_COMM_WORLD ) ;
	
	// Debug.
	/*for ( int i=0 ; i<MPISize ; i++ ){
		if ( MPIMyID == i )
			cout << "Processor: " << MPIMyID << ", MinVolume: " << h_pDomain->MinVolume << ", MinLength: " << h_pDomain->MinLength << '\n' ;
		MPI_Barrier( MPI_COMM_WORLD ) ;
	}*/
}

//==============================================================================================================

void CalculateCellCenter( DSMC_NODE *h_Node , DSMC_CELL *_Cell , CELLMAPPING *h_pMapping ){
	int		Type , NodeNo ;
	double		XCenter , YCenter , ZCenter , MaxXCoord , MinXCoord , MaxYCoord , MinYCoord , MaxZCoord , MinZCoord ;
	
	Type		= _Cell->Type ;
	XCenter		= 0. ;
	YCenter		= 0. ;
	ZCenter		= 0. ;
	MaxXCoord	= -1.E+9 ;
	MinXCoord	= 1.E+9 ;
	MaxYCoord	= -1.E+9 ;
	MinYCoord	= 1.E+9 ;
	MaxZCoord	= -1.E+9 ;
	MinZCoord	= 1.E+9 ;
	

	for ( int i=0 ; i<h_pMapping->NodeNum[Type] ; i++ ){
		XCenter	+= h_Node[_Cell->Node[i]].XCoord ;
		YCenter	+= h_Node[_Cell->Node[i]].YCoord ;
		ZCenter	+= h_Node[_Cell->Node[i]].ZCoord ;
		
		if ( h_Node[_Cell->Node[i]].XCoord > MaxXCoord ) MaxXCoord = h_Node[_Cell->Node[i]].XCoord ;
		if ( h_Node[_Cell->Node[i]].XCoord < MinXCoord ) MinXCoord = h_Node[_Cell->Node[i]].XCoord ;
			
		if ( h_Node[_Cell->Node[i]].YCoord > MaxYCoord ) MaxYCoord = h_Node[_Cell->Node[i]].YCoord ;
		if ( h_Node[_Cell->Node[i]].YCoord < MinYCoord ) MinYCoord = h_Node[_Cell->Node[i]].YCoord ;
			
		if ( h_Node[_Cell->Node[i]].ZCoord > MaxZCoord ) MaxZCoord = h_Node[_Cell->Node[i]].ZCoord ;
		if ( h_Node[_Cell->Node[i]].ZCoord < MinZCoord ) MinZCoord = h_Node[_Cell->Node[i]].ZCoord ;
		
		// Debug.
		//cout << "Node-Coord: " << setw(12) << h_Node[_Cell->Node[i]].XCoord << setw(12) << h_Node[_Cell->Node[i]].YCoord << setw(12) << h_Node[_Cell->Node[i]].ZCoord << '\n' ;
	}

	_Cell->XCenter		= XCenter/h_pMapping->NodeNum[Type] ;
	_Cell->YCenter		= YCenter/h_pMapping->NodeNum[Type] ;
	_Cell->ZCenter		= ZCenter/h_pMapping->NodeNum[Type] ;
	_Cell->MaxXCoord	= MaxXCoord ;
	_Cell->MinXCoord	= MinXCoord ;
	_Cell->MaxYCoord	= MaxYCoord ;
	_Cell->MinYCoord	= MinYCoord ;
	_Cell->MaxZCoord	= MaxZCoord ;
	_Cell->MinZCoord	= MinZCoord ;
	
	
	// Debug.
	//cout << "Center: " << setw(12) << _Cell->XCenter << setw(12) << _Cell->YCenter << setw(12) << _Cell->ZCenter << '\n' ;
	//cout << "MinX: " << _Cell->MinXCoord << ", MaxX: " << _Cell->MaxXCoord << '\n' ;
	//cout << "MinY: " << _Cell->MinYCoord << ", MaxY: " << _Cell->MaxYCoord << '\n' ;
	//cout << "MinZ: " << _Cell->MinZCoord << ", MaxZ: " << _Cell->MaxZCoord << '\n' ;
	//getchar() ;
}

//==============================================================================================================

void CalculateCellVolume2D( DSMC_NODE *h_Node , DSMC_CELL *_Cell , CELLMAPPING *h_pMapping ){
	int		Type ;
	double		Volume , x0 , x1 , y0 , y1 ;
	
	Type	= _Cell->Type ;
	Volume	= 0. ;
	
	for ( int i=0 ; i<h_pMapping->NodeNum[Type] ; i++ ){
		x0	= h_Node[_Cell->Node[i]].XCoord ;
		y0	= h_Node[_Cell->Node[i]].YCoord ;
		
		if ( i == (h_pMapping->NodeNum[Type]-1) ){
			x1	= h_Node[_Cell->Node[0]].XCoord ;
			y1	= h_Node[_Cell->Node[0]].YCoord ;
		}else{
			x1	= h_Node[_Cell->Node[i+1]].XCoord ;
			y1	= h_Node[_Cell->Node[i+1]].YCoord ;
		}
		
		Volume	+= ((x0*y1)-(y0*x1)) ;
	}
	if ( Volume <= 0. ){
		cout << "Cell volume is error. Id: " << _Cell->Id << ", Volume: " << Volume << endl ;
		exit(1) ;	
	}
	
	// calculate cell volume and characteristic length.
	_Cell->Volume			= Volume/2. ;
	_Cell->CharacteristicLength	= sqrt(_Cell->Volume) ;
}


void CalculateCellVolume3D( DSMC_NODE *h_Node , DSMC_CELL *_Cell , CELLMAPPING *h_pMapping ){
	int		Type , Node[2] ;
	double		Volume , VectorX[3] , VectorY[3] , VectorZ[3] ;
	
	Type	=	_Cell->Type ;
	Volume	=	0. ;
	
	for ( int i=0 ; i<h_pMapping->TetraCellNum[Type] ; i++ ){
		Node[0]		= h_pMapping->TetraCell[Type][i][0] ;
		Node[0]		= _Cell->Node[Node[0]] ;
		
		for ( int j=0 ; j<3 ; j++ ){
			Node[1]		= h_pMapping->TetraCell[Type][i][j+1] ;
			Node[1]		= _Cell->Node[Node[1]] ;
			
			VectorX[j]	= h_Node[Node[1]].XCoord - h_Node[Node[0]].XCoord ;
			VectorY[j]	= h_Node[Node[1]].YCoord - h_Node[Node[0]].YCoord ;
			VectorZ[j]	= h_Node[Node[1]].ZCoord - h_Node[Node[0]].ZCoord ;
		}
		
		Volume	+= fabs(( VectorX[0]*(VectorY[1]*VectorZ[2] - VectorZ[1]*VectorY[2]) 
				+ VectorY[0]*(VectorZ[1]*VectorX[2] - VectorX[1]*VectorZ[2]) 
				+ VectorZ[0]*(VectorX[1]*VectorY[2] - VectorY[1]*VectorX[2])) / 6.) ;
	}
	if ( Volume <= 0. ){
		cout << "Cell volume is error. Id: " << _Cell->Id << ", Volume: " << Volume << endl ;
		exit(1) ;	
	}
	
	// calculate cell volume and characteristic length.
	_Cell->Volume			= Volume ;
	_Cell->CharacteristicLength	= pow(_Cell->Volume , 1./3.) ;
}


void CalculateCellVolumeAxisymmetric( DSMC_NODE *h_Node , DSMC_CELL *_Cell , CELLMAPPING *h_pMapping ){
	int		Type ;
	double		Volume , x0 , x1 , y0 , y1 ;
	
	Type	= _Cell->Type ;
	Volume	= 0. ;
	
	for ( int i=0 ; i<h_pMapping->NodeNum[Type] ; i++ ){
		x0	= h_Node[_Cell->Node[i]].XCoord ;
		y0	= h_Node[_Cell->Node[i]].YCoord ;
		
		if ( i == (h_pMapping->NodeNum[Type]-1) ){
			x1	= h_Node[_Cell->Node[0]].XCoord ;
			y1	= h_Node[_Cell->Node[0]].YCoord ;
		}else{
			x1	= h_Node[_Cell->Node[i+1]].XCoord ;
			y1	= h_Node[_Cell->Node[i+1]].YCoord ;
		}
		
		Volume	+= ((x0*y1)-(y0*x1)) ;
	}
	if ( Volume <= 0. ){
		cout << "Cell volume is error. Id: " << _Cell->Id << ", Volume: " << Volume << endl ;
		exit(1) ;	
	}
	
	
	Volume		= Volume/2. ;
	// calculate cell volume and characteristic length.
	_Cell->Volume			= 2.*PI*_Cell->YCenter*Volume ;
	//_Cell->CharacteristicLength	= pow((_Cell->Volume/(2.*PI)) , 1./3.) ;
	_Cell->CharacteristicLength	= sqrt(Volume) ;
}

//==============================================================================================================

void CreateCellFaceFunction2D( DSMC_NODE *h_Node , DSMC_CELL *_Cell , CELLMAPPING *h_pMapping ){
	int		Type , Node0 , Node1 ;
	double		x0 , y0 , x1 , y1 ;

	Type	= _Cell->Type ;

	for ( int i=0 ; i<h_pMapping->SurfaceNum[Type] ; i++ ){
		Node0	= h_pMapping->Node[Type][i][0] ;
		Node1	= h_pMapping->Node[Type][i][1] ;
	
	
		Node0	= _Cell->Node[Node0] ;
		Node1	= _Cell->Node[Node1] ;
	
		x0	= h_Node[Node0].XCoord ;
		y0	= h_Node[Node0].YCoord ;
		x1	= h_Node[Node1].XCoord ;
		y1	= h_Node[Node1].YCoord ;

		
		_Cell->FaceFA[i]	= y1-y0 ;
		_Cell->FaceFB[i]	= x0-x1 ;
		_Cell->FaceFC[i]	= -1. * (_Cell->FaceFA[i]*x0 + _Cell->FaceFB[i]*y0) ;
	}
}


void CreateCellFaceFunction3D( DSMC_NODE *h_Node , DSMC_CELL *_Cell , CELLMAPPING *h_pMapping ){
	int		Type , Node[2] ;
	double		VectorX[2] , VectorY[2] , VectorZ[2] , FaceFA , FaceFB , FaceFC , Norm ;
	
	Type	= _Cell->Type ;
	
	// Create each face function of a cell.
	for ( int i=0 ; i<h_pMapping->SurfaceNum[Type] ; i++ ){
		Node[0]		= h_pMapping->Node[Type][i][0] ;
		Node[0]		= _Cell->Node[Node[0]] ;
		
		for ( int j=0 ; j<2 ; j++ ){
			Node[1]	= h_pMapping->Node[Type][i][j+1] ;
			Node[1]	= _Cell->Node[Node[1]] ;
			
			VectorX[j]	= h_Node[Node[1]].XCoord - h_Node[Node[0]].XCoord ;
			VectorY[j]	= h_Node[Node[1]].YCoord - h_Node[Node[0]].YCoord ;
			VectorZ[j]	= h_Node[Node[1]].ZCoord - h_Node[Node[0]].ZCoord ;
		}
	
	
		FaceFA	= VectorY[0]*VectorZ[1] - VectorZ[0]*VectorY[1] ;
		FaceFB	= VectorZ[0]*VectorX[1] - VectorX[0]*VectorZ[1] ;
		FaceFC	= VectorX[0]*VectorY[1] - VectorY[0]*VectorX[1] ;
		Norm	= sqrt( FaceFA*FaceFA + FaceFB*FaceFB + FaceFC*FaceFC ) ;
		
		_Cell->FaceFA[i]	=  FaceFA/Norm ;
		_Cell->FaceFB[i]	=  FaceFB/Norm ;
		_Cell->FaceFC[i]	=  FaceFC/Norm ;
		_Cell->FaceFD[i]	= -1. * ( _Cell->FaceFA[i]*h_Node[Node[0]].XCoord + 
						  _Cell->FaceFB[i]*h_Node[Node[0]].YCoord + 
						  _Cell->FaceFC[i]*h_Node[Node[0]].ZCoord ) ;
	}
}

//==============================================================================================================
//==============================================================================================================

void CalculateTimestepWeighting( DSMC_DOMAIN		*h_pDomain ,
				 DSMC_CELL		*h_Cell ,
				 DSMC_SPECIES		*h_Species ,
				 DSMC_DSMC		*h_pDSMC ,
				 CELLMAPPING		*h_pMapping , 
				 ofstream		&OutputDebug ){

	int		CellNum , Dimension , SubCellNum ;
	double		Timestep , ParticleWeighting , XVel , YVel , ZVel , Temp , Mass , CrossSection , MinLength , MinVolume , MinCenter ; 
	double		Coeff , MinX , ParticleNum ;
	double Tref ,Wab ,Mr, *CollisionTime;
	
	Dimension	= h_pDomain->Dimension ;
	CellNum		= h_pDomain->CellNum ;
	XVel		= h_pDomain->XVel ;
	YVel		= h_pDomain->YVel ;
	ZVel		= h_pDomain->ZVel ;
	Temp		= h_pDomain->Temp ;
	Mass		= h_Species[0].Mass ;
	MinLength	= h_pDomain->MinLength ;
	MinVolume	= h_pDomain->MinVolume ;
	MinCenter	= h_pDomain->MinCenter ;
	
	CrossSection	= 0. ;
	Coeff		= 0. ;
  CollisionTime = new double [h_pDomain->SpeciesGroupNum]; 
	
	// To calculate the timestep and particle weighting based on the minimum cell volume.
	if ( h_pDomain->TimestepRatio == 0 ){
		Timestep		= h_pDomain->Timestep ;
	}else{
		Timestep		= h_pDomain->MinLength/(3.*(sqrt(XVel*XVel+YVel*YVel+ZVel*ZVel)+sqrt(2.*BOLTZ*Temp/Mass))) * h_pDomain->TimestepRatio ;
	}
	
	if ( h_pDomain->WeightingRatio == 0 ){
		ParticleWeighting	= h_pDomain->ParticleWeighting ;
	}else{
		ParticleWeighting	= h_pDomain->NumDen*MinVolume/h_pDomain->WeightingRatio ;
	}
	

	// Debug.
	//OutputDebug << "In Init-CalTimeWeight, Timestep: " << Timestep << ", Weighting: " << ParticleWeighting << endl ;
	//OutputDebug << "NumDen: " << h_pDomain->NumDen << ", MinVolume: " << h_pDomain->MinVolume << ", WeightRatio: " << h_pDomain->WeightingRatio << endl ;


	for ( int i=0 ; i<CellNum ; i++ ){
		for ( int j=0 ; j<h_pDomain->SpeciesGroupNum ; j++ ){
			
			CollisionTime[j]= 0.;
			
			for ( int k=0 ; k<h_pDomain->SpeciesGroupNum ; k++ ){
				
				CrossSection	= 0.25*PI*pow((h_Species[j].Diameter+h_Species[k].Diameter),2.) ;
				
				Tref = 0.5*(h_Species[j].RefTemp + h_Species[k].RefTemp) ;
				Wab = 0.5*(h_Species[j].VisTempIndex + h_Species[k].VisTempIndex) ;
				Mr = h_Species[j].Mass* h_Species[k].Mass/ (h_Species[j].Mass + h_Species[k].Mass) ;
								
				CollisionTime[j] += 2.* CrossSection * h_pDomain->NumDen*h_Species[k].Fraction*pow((Temp/Tref),(1-Wab))*sqrt(2*BOLTZ*Tref/(Mr*PI));			
				
				h_pDSMC->RemainderCollisionPair[i][j][k]	= Randn() ;
				h_pDSMC->MaxCrossSectionSpeed[i][j][k]		= CrossSection*300.*sqrt(Temp/300.) ;
			}
		}

//		h_Cell[i].Timestep	= 0.2/CollisionTime[0];

		h_Cell[i].Timestep	= Timestep ;
		h_Cell[i].Weighting	= ParticleWeighting ;
	
		// Variable Timestep scheme.
		if ( h_pDomain->Dimension != 4  ){
			if ( h_pDomain->VariableTimestepScheme == 1 ){
				h_Cell[i].Timestep	= Timestep * (h_Cell[i].CharacteristicLength/MinLength) ;
				h_Cell[i].Weighting	= ParticleWeighting * (h_Cell[i].CharacteristicLength/MinLength) ;
			}else if ( h_pDomain->VariableTimestepScheme == 2 ){
				
	
				
				
			}else if ( h_pDomain->VariableTimestepScheme == 3 ){


				
			}
		}else if ( h_pDomain->Dimension == 4 ){
			if ( h_pDomain->VariableTimestepScheme == 4 ){
				h_Cell[i].Timestep	= Timestep * (h_Cell[i].YCenter/MinCenter) * (h_Cell[i].CharacteristicLength/MinLength) * h_pDomain->AxisAdjustFactor ;
				h_Cell[i].Weighting	= ParticleWeighting * (h_Cell[i].YCenter/MinCenter) * (h_Cell[i].CharacteristicLength/MinLength) * h_pDomain->AxisAdjustFactor ;
			}else if ( h_pDomain->VariableTimestepScheme == 5 ){
				//h_Cell[i].Timestep	= Timestep * sqrt(h_Cell[i].YCenter/MinCenter) * (h_Cell[i].CharacteristicLength/MinLength) * h_pDomain->AxisAdjustFactor ;
				//h_Cell[i].Weighting	= ParticleWeighting * sqrt(h_Cell[i].YCenter/MinCenter) * (h_Cell[i].CharacteristicLength/MinLength) * h_pDomain->AxisAdjustFactor ;
				
				h_Cell[i].Timestep	= h_Cell[i].InitTimestep ;
				h_Cell[i].Weighting	= h_Cell[i].InitWeighting/h_pDomain->WeightingRatio ;
				
			}else if ( h_pDomain->VariableTimestepScheme == 6 ){
				h_Cell[i].Timestep	= h_Cell[i].InitTimestep ;
				h_Cell[i].Weighting	= h_Cell[i].InitWeighting/h_pDomain->WeightingRatio ;
				
				if ( h_Cell[i].MeanCollSpacingMeanFreePath >= 1. ){
					SubCellNum	= h_Cell[i].CharacteristicLength/(0.5*h_Cell[i].InitMeanFreePath) ;
					ParticleNum	= SubCellNum*SubCellNum*2 + 2 ;
				}else{
					ParticleNum	= 10. ;
				}
					
				if ( h_Cell[i].AveParticleNum < ParticleNum ){
					h_Cell[i].Timestep	= h_Cell[i].Timestep/(ParticleNum/h_Cell[i].AveParticleNum) ;
					h_Cell[i].Weighting	= h_Cell[i].Weighting/(ParticleNum/h_Cell[i].AveParticleNum) ;
				}
			}
		}
		
		
		/*if ( h_pDomain->VariableTimestepScheme == 1 || h_pDomain->VariableTimestepScheme == 2 ){
			
			if ( Dimension == 4 ){
				h_Cell[i].Timestep	= Timestep * (h_Cell[i].YCenter/MinCenter) * (h_Cell[i].CharacteristicLength/MinLength) * 0.005 ;
				h_Cell[i].Weighting	= ParticleWeighting * (h_Cell[i].YCenter/MinCenter) * (h_Cell[i].CharacteristicLength/MinLength) * 0.005 ;
				//h_Cell[i].Timestep	= Timestep * (h_Cell[i].Volume/MinVolume) ;
				//h_Cell[i].Weighting	= ParticleWeighting * (h_Cell[i].Volume/MinVolume) ;
			}else{
				h_Cell[i].Timestep	= Timestep * (h_Cell[i].CharacteristicLength/MinLength) ;
				h_Cell[i].Weighting	= ParticleWeighting * (h_Cell[i].CharacteristicLength/MinLength) ;
			}
			
		}else if ( h_pDomain->VariableTimestepScheme == 3 ){
			
			if ( h_Cell[i].CharacteristicLength < h_Cell[i].InitMeanFreePath )
				MinX	= h_Cell[i].CharacteristicLength ;
			else
				MinX	= h_Cell[i].InitMeanFreePath ;
			
			//h_Cell[i].Timestep	= h_Cell[i].InitMeanFreePath/(3.*(h_Cell[i].InitSpeed+sqrt(2.*BOLTZ*h_Cell[i].InitTemp/Mass))) ;
			h_Cell[i].Timestep	= MinX/(3.*(h_Cell[i].InitSpeed+sqrt(2.*BOLTZ*h_Cell[i].InitTemp/Mass))) ;
			Coeff			= h_Cell[i].Timestep / Timestep ;
			h_Cell[i].Weighting	= ParticleWeighting * Coeff ;
			
		}else{
			h_Cell[i].Timestep	= Timestep ;
			h_Cell[i].Weighting	= ParticleWeighting ;
		}*/
		
		// Debug.
		//cout << "CellNo: " << i  << ", Timestep: " << h_Cell[i].Timestep << ", Weighting: " << h_Cell[i].Weighting << '\n' ;
		//if ( i%250 == 0 ) getchar() ;
	}
	
	OutputDebug << "CollisionTime: " << CollisionTime[0] << endl ;
	
	h_pDomain->Timestep		= Timestep ;
	h_pDomain->ParticleWeighting	= ParticleWeighting ;
	delete [] CollisionTime ;	
	// Debug.
	//cout << "General Timestep: " << h_pDomain->Timestep << ", Weighting: " << h_pDomain->ParticleWeighting << '\n' ;
}

//==============================================================================================================
//==============================================================================================================

void CreateInletInformation( 	DSMC_DOMAIN		*h_pDomain ,
				DSMC_NODE		*h_Node ,
				DSMC_CELL		*h_Cell ,
				DSMC_INLET		*h_Inlet ,
				DSMC_SPECIES		*h_Species ,
				DSMC_DSMC		*h_pDSMC ,
				CELLMAPPING		*h_pMapping ){
					
	int		MPISize , MPIMyID , Type , FaceNum , CellNo , NodeNo[4] , Dimension , DebugFaceNum ;
	int		*FaceCellNo , *FaceCellCheck , *FaceLocalCellCheck ;
	double		XVel , YVel , ZVel , NumDen , Temp , Timestep , Weighting ;
	double		MostProbableSpeed , NormalVel , SpeedRatio , spi , Flux ;
	bool		CheckNode ;
	
	MPI_Comm_size( MPI_COMM_WORLD , &MPISize ) ;
	MPI_Comm_rank( MPI_COMM_WORLD , &MPIMyID ) ;
	
	FaceNum		= 0 ;
	DebugFaceNum	= 0 ;
	Dimension	= h_pDomain->Dimension ;
	spi		= sqrt(PI) ;
	
	FaceCellNo		= new int[h_pDomain->InletFaceNum] ;
	FaceCellCheck		= new int[h_pDomain->InletFaceNum] ;
	FaceLocalCellCheck	= new int[h_pDomain->InletFaceNum] ;


	for ( int i=0 ; i<h_pDomain->InletFaceNum ; i++ ){
		FaceCellNo[i]		= 0 ;
		FaceCellCheck[i]	= 0 ;
		FaceLocalCellCheck[i]	= 0 ;
	}


	for ( int i=0 ; i<h_pDomain->CellNum ; i++ ){
		Type	= h_Cell[i].Type ;
		for ( int j=0 ; j<h_pMapping->SurfaceNum[Type] ; j++ ){
			if ( h_Cell[i].Neighbor[j] == -3 || h_Cell[i].Neighbor[j] == -4 ){
				FaceCellNo[FaceNum]	= i ;
				//FaceLocalCellCheck[FaceNum]	= 0 ;
				FaceNum++ ;
			}
		}	
	}
	

	// Debug.
	//cout << "MyID: " << MPIMyID << "FaceNum: " << FaceNum << '\n' ;
	//MPI_Barrier( MPI_COMM_WORLD ) ;
	

	// To create the relation between inlet faces and cells.
	for ( int f1=0 ; f1<h_pDomain->InletFaceNum ; f1++ ){
		//for ( int f2=0 ; f2<h_pDomain->InletFaceNum ; f2++ ){
		for ( int f2=0 ; f2<FaceNum ; f2++ ){
			//if ( FaceCellCheck[f2] == 1 ) continue ;
			if ( FaceLocalCellCheck[f2] == 1 ) continue ;
				
			CellNo	= FaceCellNo[f2] ;
			Type	= h_Cell[CellNo].Type ;
		
			for ( int i=0 ; i<h_pMapping->SurfaceNum[Type] ; i++ ){
				if ( h_Cell[CellNo].Neighbor[i] == -3 || h_Cell[CellNo].Neighbor[i] == -4 ){
					CheckNode	= true ;
					
					for ( int j=0 ; j<h_pMapping->SurfaceNode[Type][i] ; j++ ){
						NodeNo[j]	= h_Cell[CellNo].Node[h_pMapping->Node[Type][i][j]] ;
						if ( h_Inlet[f1].Node[j] != NodeNo[j]) CheckNode = false ;
					}
					
					if ( CheckNode ){
						h_Inlet[f1].CellNo	= CellNo ;
						h_Inlet[f1].FaceNo	= i ;
						FaceCellCheck[f1]	= 1 ;
						FaceLocalCellCheck[f2]	= 1 ;
						DebugFaceNum++ ;
						break ;
					}
				}
			}
			
			if ( FaceLocalCellCheck[f2] == 1 ) break ;
		}
	}


	// Debug.
	//cout << "MyID: " << MPIMyID << ", DebugFaceNum: " << DebugFaceNum << '\n' ;
	
	// Debug.
	/*for ( int i=0 ; i<MPISize ; i++ ){
		if ( MPIMyID == i ){
			cout << "Processor: " << MPIMyID << '\n' ;
			//for ( int j=0 ; j<h_pDomain->InletFaceNum ; j++ )
			//	cout << setw(12) << h_Inlet[j].CellNo << setw(12) << h_Inlet[j].FaceNo << setw(12) << FaceCellCheck[j] << '\n' ;	
			for ( int j=0 ; j<FaceNum ; j++ )
				cout << setw(12) << FaceLocalCellCheck[j] << '\n' ;	
		}
		MPI_Barrier( MPI_COMM_WORLD ) ;	
	}*/
	
	FaceNum	= 0 ;
	for ( int i=0 ; i<h_pDomain->InletFaceNum ; i++ ){
		if ( FaceCellCheck[i] == 1 ){
			h_Inlet[FaceNum]	= h_Inlet[i] ;
			FaceNum++ ;
		}
	}
	//Update the inlet face number of each processor.
	h_pDomain->InletFaceNum	= FaceNum ;
	
	// Debug.
	/*for ( int i=0 ; i<MPISize ; i++ ){
		if ( MPIMyID == i ){
			cout << "Processor: " << MPIMyID << '\n' ;
			for ( int j=0 ; j<h_pDomain->InletFaceNum ; j++ )
				cout << setw(12) << h_Inlet[j].CellNo << setw(12) << h_Inlet[j].FaceNo << setw(12) << '\n' ;	
		}
		MPI_Barrier( MPI_COMM_WORLD ) ;	
	}*/
	// Debug.
	//MPI_Barrier( MPI_COMM_WORLD ) ;	
	//cout << "MyID: " << MPIMyID << ", FaceNum: " << h_pDomain->InletFaceNum << '\n' ;
	

	
	for ( int i=0 ; i<h_pDomain->InletFaceNum ; i++ ){
		// Calculate the area of inlet face.
		if ( Dimension == 2 ){
			CalculateInletFaceArea2D( h_Node , &h_Inlet[i] ) ;
		}else if ( Dimension == 3 ){
			CalculateInletFaceArea3D( h_Node , &h_Inlet[i] ) ;
		}else if ( Dimension == 4 ){
			CalculateInletFaceAreaAxisymmetric( h_Node , &h_Inlet[i] ) ;
		}
		
		
		XVel		= h_Inlet[i].XVel ;
		YVel		= h_Inlet[i].YVel ;
		ZVel		= h_Inlet[i].ZVel ;
		//NumDen		= h_Inlet[i].NumDen ;
		Temp		= h_Inlet[i].Temp ;
		CellNo		= h_Inlet[i].CellNo ;
		Timestep	= h_Cell[CellNo].Timestep ;
		Weighting	= h_Cell[CellNo].Weighting;
		
		
		// Calculate the flux from inlet face for each species
		for ( int j=0 ; j<h_pDomain->SpeciesNum ; j++ ){
			if ( h_pDomain->InletSpecifiedNumDenNum == 1 ){
				NumDen	= h_Inlet[i].NumDen[0] ;
			}else if ( h_pDomain->InletSpecifiedNumDenNum > 1){
				NumDen	= h_Inlet[i].NumDen[j] ;
			}
			
			
			MostProbableSpeed	= sqrt(2.*BOLTZ*Temp/h_Species[j].Mass) ;
			
			if ( Dimension == 2 || Dimension == 4 ){
				// Calculate a normal velocity from inlet face into the cell.
				NormalVel	= -1.*CalculateNormalVelocity2D( XVel ,YVel , &h_Node[h_Inlet[i].Node[0]] , &h_Node[h_Inlet[i].Node[1]] ) ;
			}else if ( Dimension == 3 ){
				NormalVel	= -1.*CalculateNormalVelocity3D( XVel , YVel , ZVel , &h_Inlet[i] , &h_Cell[CellNo] ) ;
			}
			
			SpeedRatio	= NormalVel/MostProbableSpeed ;
		
			// the non-dimensional flux of eqn (4.22)
			if ( fabs(SpeedRatio) < 10.1 )	Flux = (exp(-SpeedRatio*SpeedRatio) + spi*SpeedRatio*(1.+ErrorFunction(SpeedRatio))) / (2.*spi) ;
			if ( SpeedRatio > 10. )		Flux = SpeedRatio ;
			if ( SpeedRatio < -10. )	Flux = 0. ;
			
			
			if ( h_Inlet[i].CosineLawCoef < 0. ){
				h_pDSMC->InletEnterNum[i][j]	= NumDen*h_Species[j].Fraction*Flux*MostProbableSpeed*Timestep*h_Inlet[i].Area/Weighting ;
			}else if ( h_Inlet[i].CosineLawCoef >= 0. ){
				h_pDSMC->InletEnterNum[i][j]	= NumDen*h_Species[j].Fraction*Timestep*h_Inlet[i].Area/Weighting ;
			}
		}
	}

	delete [] FaceCellNo ;
	delete [] FaceCellCheck ;
	delete [] FaceLocalCellCheck ;
}


//==============================================================================================================

void CalculateInletFaceArea2D( DSMC_NODE *h_Node , DSMC_INLET *_Inlet ){
	double		x0 , y0 , x1 , y1 , a , b ;
	
	x0	= h_Node[_Inlet->Node[0]].XCoord ;
	y0	= h_Node[_Inlet->Node[0]].YCoord ;
	x1	= h_Node[_Inlet->Node[1]].XCoord ;
	y1	= h_Node[_Inlet->Node[1]].YCoord ;
	a	= x1 - x0 ;
	b	= y1 - y0 ;
	
	_Inlet->Area	= sqrt((a*a)+(b*b)) ;
}


void CalculateInletFaceArea3D( DSMC_NODE *h_Node , DSMC_INLET *_Inlet ){
	int		NodeNum ;
	double		Area , VectorX[2] , VectorY[2] , VectorZ[2] , Buffer[3] ;
	
	NodeNum	= _Inlet->NodeNum ;
	Area	= 0. ;
	
	VectorX[0]	= h_Node[_Inlet->Node[1]].XCoord - h_Node[_Inlet->Node[0]].XCoord ;
	VectorY[0]	= h_Node[_Inlet->Node[1]].YCoord - h_Node[_Inlet->Node[0]].YCoord ;
	VectorZ[0]	= h_Node[_Inlet->Node[1]].ZCoord - h_Node[_Inlet->Node[0]].ZCoord ;
	
	VectorX[1]	= h_Node[_Inlet->Node[NodeNum-1]].XCoord - h_Node[_Inlet->Node[0]].XCoord ;
	VectorY[1]	= h_Node[_Inlet->Node[NodeNum-1]].YCoord - h_Node[_Inlet->Node[0]].YCoord ;
	VectorZ[1]	= h_Node[_Inlet->Node[NodeNum-1]].ZCoord - h_Node[_Inlet->Node[0]].ZCoord ;
	
	Buffer[0]	= VectorY[0]*VectorZ[1] - VectorZ[0]*VectorY[1] ;
	Buffer[1]	= VectorZ[0]*VectorX[1] - VectorX[0]*VectorZ[1] ;
	Buffer[2]	= VectorX[0]*VectorY[1] - VectorY[0]*VectorX[1] ;
	
	Area	+= (sqrt( Buffer[0]*Buffer[0] + Buffer[1]*Buffer[1] + Buffer[2]*Buffer[2] )/2.) ;
	
	if ( NodeNum == 4 ){
		VectorX[0]	= h_Node[_Inlet->Node[1]].XCoord - h_Node[_Inlet->Node[2]].XCoord ;
		VectorY[0]	= h_Node[_Inlet->Node[1]].YCoord - h_Node[_Inlet->Node[2]].YCoord ;
		VectorZ[0]	= h_Node[_Inlet->Node[1]].ZCoord - h_Node[_Inlet->Node[2]].ZCoord ;
	
		VectorX[1]	= h_Node[_Inlet->Node[NodeNum-1]].XCoord - h_Node[_Inlet->Node[2]].XCoord ;
		VectorY[1]	= h_Node[_Inlet->Node[NodeNum-1]].YCoord - h_Node[_Inlet->Node[2]].YCoord ;
		VectorZ[1]	= h_Node[_Inlet->Node[NodeNum-1]].ZCoord - h_Node[_Inlet->Node[2]].ZCoord ;
	
		Buffer[0]	= VectorY[0]*VectorZ[1] - VectorZ[0]*VectorY[1] ;
		Buffer[1]	= VectorZ[0]*VectorX[1] - VectorX[0]*VectorZ[1] ;
		Buffer[2]	= VectorX[0]*VectorY[1] - VectorY[0]*VectorX[1] ;
	
		Area	+= (sqrt( Buffer[0]*Buffer[0] + Buffer[1]*Buffer[1] + Buffer[2]*Buffer[2] )/2.) ;
	}
	
	
	// Debug.
	/*cout << _Inlet->CellNo << ", " << _Inlet->FaceNo << '\n' ;
	cout << _Inlet->Node[0] << ", " << _Inlet->Node[1] << ", " << _Inlet->Node[2] << ", " << _Inlet->Node[3] << '\n' ;
	cout << "Area: " << Area << '\n' ;
	getchar() ;*/
	
	_Inlet->Area	= Area ;
}


void CalculateInletFaceAreaAxisymmetric( DSMC_NODE *h_Node , DSMC_INLET *_Inlet ){
	double		x0 , y0 , x1 , y1 , a , b , Area ;
	
	x0	= h_Node[_Inlet->Node[0]].XCoord ;
	y0	= h_Node[_Inlet->Node[0]].YCoord ;
	x1	= h_Node[_Inlet->Node[1]].XCoord ;
	y1	= h_Node[_Inlet->Node[1]].YCoord ;
	a	= x1 - x0 ;
	b	= y1 - y0 ;
	
	Area		= sqrt((a*a)+(b*b)) ;
	_Inlet->Area	= (y0+y1)*Area*PI ;
}

//==============================================================================================================
//==============================================================================================================

void LinkCellSurface(	DSMC_DOMAIN		*h_pDomain ,
			DSMC_NODE		*h_Node ,
			DSMC_CELL		*h_Cell ,
			DSMC_SURFACE		*h_Surface ,
			CELLMAPPING		*h_pMapping ){
				
	int		Type , SurfaceNo , Dimension ;
	
	Dimension	= h_pDomain->Dimension ;
	SurfaceNo	= 0 ;
	
	for ( int i=0 ; i<h_pDomain->CellNum ; i++ ){
		Type	= h_Cell[i].Type ;
		
		for ( int j=0 ; j<h_pMapping->SurfaceNum[Type] ; j++ ){
			if ( h_Cell[i].Neighbor[j] < -20 ){
				h_Cell[i].Surface[j]	= SurfaceNo ;
				
				h_Surface[SurfaceNo].WallNo	= abs(h_Cell[i].Neighbor[j] + 21) ;
				h_Surface[SurfaceNo].CellNo	= i ;
				h_Surface[SurfaceNo].FaceNo	= j ;
				
				if ( Dimension == 2 ){
					// Calculate area and center of the surface. 
					CalculateSurfaceAreaCenter2D( h_Node , &h_Cell[i] , &h_Surface[SurfaceNo] , h_pMapping , j ) ;
				}else if ( Dimension == 3 ){
					CalculateSurfaceAreaCenter3D( h_Node , &h_Cell[i] , &h_Surface[SurfaceNo] , h_pMapping , j ) ;
				}else if ( Dimension ==4 ){
					CalculateSurfaceAreaCenterAxisymmetric( h_Node , &h_Cell[i] , &h_Surface[SurfaceNo] , h_pMapping , j ) ;
				}
				
				SurfaceNo++ ;
			}	
		}
	}
	
	if ( SurfaceNo > h_pDomain->WallFaceNum ){
		cout << "Wall_Face_Number is error. It is " << h_pDomain->WallFaceNum << ". But the wall face number is actually " << SurfaceNo << '.' << endl ;
		exit(1) ;
	}
	
	h_pDomain->WallFaceNum	= SurfaceNo ;
	
	// Debug.
	//cout << "WallFaceNum: " << h_pDomain->WallFaceNum << '\n' ;

	
	// Debug.
	//for ( int i=0 ; i<h_pDomain->WallFaceNum ; i++ ) h_Surface[i].Dump(i) ;
}

//==============================================================================================================

void CalculateSurfaceAreaCenter2D(	DSMC_NODE	*h_Node , 
					DSMC_CELL	*_Cell , 
					DSMC_SURFACE	*_Surface ,
					CELLMAPPING	*h_pMapping , 
					int		FaceNo ){
	int		Type , Node0 , Node1 ;
	double		x0 , y0 , x1 , y1 , a , b ;
	
	Type	= _Cell->Type ;
	Node0	= h_pMapping->Node[Type][FaceNo][0] ;
	Node1	= h_pMapping->Node[Type][FaceNo][1] ;
	
	x0	= h_Node[_Cell->Node[Node0]].XCoord ;
	y0	= h_Node[_Cell->Node[Node0]].YCoord ;
	x1	= h_Node[_Cell->Node[Node1]].XCoord ;
	y1	= h_Node[_Cell->Node[Node1]].YCoord ;
	a	= x1 - x0 ;
	b	= y1 - y0 ;
	
	_Surface->XCenter	= (x0+x1)/2. ;
	_Surface->YCenter	= (y0+y1)/2. ;
	_Surface->Area		= sqrt((a*a)+(b*b)) ;
}


void CalculateSurfaceAreaCenter3D(	DSMC_NODE	*h_Node , 
					DSMC_CELL	*_Cell , 
					DSMC_SURFACE	*_Surface ,
					CELLMAPPING	*h_pMapping , 
					int		FaceNo ){
	
	int		Type , Node[2] ;
	double		Area , XCenter , YCenter , ZCenter , VectorX[2] , VectorY[2] , VectorZ[2] , Buffer[3] ;
	
	Area	= 0. ;
	XCenter	= 0. ;
	YCenter	= 0. ;
	ZCenter	= 0. ;
	Type	= _Cell->Type ;
	
	Node[0]		= h_pMapping->Node[Type][FaceNo][0] ;
	Node[0]		= _Cell->Node[Node[0]] ;
	for ( int i=0 ; i<h_pMapping->SurfaceNode[Type][FaceNo] ; i++ ){
		Node[1]		= h_pMapping->Node[Type][FaceNo][i] ;
		Node[1]		= _Cell->Node[Node[1]] ;
		
		XCenter		+= h_Node[Node[1]].XCoord ;
		YCenter		+= h_Node[Node[1]].YCoord ;
		ZCenter		+= h_Node[Node[1]].ZCoord ;
		
		// Calculate Vector "0".
		VectorX[0]	= h_Node[Node[1]].XCoord - h_Node[Node[0]].XCoord ;
		VectorY[0]	= h_Node[Node[1]].YCoord - h_Node[Node[0]].YCoord ;
		VectorZ[0]	= h_Node[Node[1]].ZCoord - h_Node[Node[0]].ZCoord ;
		
		if ( i>1 ){
			Buffer[0]	= VectorY[0]*VectorZ[1] - VectorZ[0]*VectorY[1] ;
			Buffer[1]	= VectorZ[0]*VectorX[1] - VectorX[0]*VectorZ[1] ;
			Buffer[2]	= VectorX[0]*VectorY[1] - VectorY[0]*VectorX[1] ;
	
			Area	+= (sqrt( Buffer[0]*Buffer[0] + Buffer[1]*Buffer[1] + Buffer[2]*Buffer[2] )/2.) ;
		}
		
		VectorX[1]	= VectorX[0] ;
		VectorY[1]	= VectorY[0] ;
		VectorZ[1]	= VectorZ[0] ;
	}
	
	
	_Surface->XCenter	= XCenter/h_pMapping->SurfaceNode[Type][FaceNo] ;
	_Surface->YCenter	= YCenter/h_pMapping->SurfaceNode[Type][FaceNo] ;
	_Surface->ZCenter	= ZCenter/h_pMapping->SurfaceNode[Type][FaceNo] ;
	_Surface->Area		= Area ;
	
	// Debug.
	//cout << "XC: " << _Surface->XCenter << ", YC: " << _Surface->YCenter << ", ZC: " << _Surface->ZCenter << '\n' ;
	//cout << "Area: " << _Surface->Area << '\n' ;
	//getchar() ;
}


void CalculateSurfaceAreaCenterAxisymmetric(	DSMC_NODE	*h_Node , 
						DSMC_CELL	*_Cell , 
						DSMC_SURFACE	*_Surface ,
						CELLMAPPING	*h_pMapping , 
						int		FaceNo ){
							
	int		Type , Node0 , Node1 ;
	double		x0 , y0 , x1 , y1 , a , b , Area ;
	
	Type	= _Cell->Type ;
	Node0	= h_pMapping->Node[Type][FaceNo][0] ;
	Node1	= h_pMapping->Node[Type][FaceNo][1] ;
	
	x0	= h_Node[_Cell->Node[Node0]].XCoord ;
	y0	= h_Node[_Cell->Node[Node0]].YCoord ;
	x1	= h_Node[_Cell->Node[Node1]].XCoord ;
	y1	= h_Node[_Cell->Node[Node1]].YCoord ;
	a	= x1 - x0 ;
	b	= y1 - y0 ;
	
	Area	= sqrt((a*a)+(b*b)) ;
	
	_Surface->XCenter	= (x0+x1)/2. ;
	_Surface->YCenter	= (y0+y1)/2. ;
	_Surface->Area		= (y0+y1)*Area*PI ;
}

//==============================================================================================================
//==============================================================================================================

void CreateInitialParticle(	DSMC_DOMAIN		*h_pDomain ,
				DSMC_NODE		*h_Node ,
				DSMC_CELL		*h_Cell ,
				DSMC_SPECIES		*h_Species ,
				DSMC_DSMC		*h_pDSMC ,
				CELLMAPPING		*h_pMapping , 
				ofstream		&OutputDebug ){
	
	int		ParticleNo , MaxParticleNum  , Dimension , ParticleNumINT ;
	double		ParticleNum , RemainderParticleNum , NumDen , Temp ;
	double		XCoord , YCoord , ZCoord , MostProbableSpeed , buffer ;
	
	Dimension		= h_pDomain->Dimension ;
	NumDen			= h_pDomain->NumDen ;
	Temp			= h_pDomain->Temp ;
	MaxParticleNum		= h_pDomain->MaxParticleNum ;
	
	ParticleNo		= 0 ;
	ParticleNum		= 0. ;
	RemainderParticleNum	= 0. ;
	XCoord			= 0. ;
	YCoord			= 0. ;
	ZCoord			= 0. ;
	buffer			= 0. ;
	
	for ( int i=0 ; i<h_pDomain->CellNum ; i++ ){
		for ( int j=0 ; j<h_pDomain->SpeciesNum ; j++ ){
			ParticleNum		= (NumDen*h_Cell[i].Volume*h_Species[j].Fraction/h_Cell[i].Weighting) + RemainderParticleNum ;
			ParticleNumINT		= (int)ParticleNum ;
			RemainderParticleNum	= ParticleNum - ParticleNumINT ;

			
			MostProbableSpeed	= sqrt(2.*BOLTZ*Temp/h_Species[j].Mass) ;
			

			for ( int k=0 ; k<ParticleNumINT ; k++ ){
				h_pDSMC->ParticleCellNo[ParticleNo]	= i ;
				h_pDSMC->ParticleSpeciesNo[ParticleNo]	= j ;
				h_pDSMC->ParticleLastCollide[ParticleNo]= -1 ;
				
				if ( Dimension == 2 || Dimension == 4 ){
					CreateParticlePosition2D( &XCoord  , &YCoord , &ZCoord , h_Node , &h_Cell[i] , h_pMapping ) ;
				}else if ( Dimension == 3 ){
					CreateParticlePosition3D( &XCoord  , &YCoord , &ZCoord , h_Node , &h_Cell[i] , h_pMapping ) ;
				}
				
				h_pDSMC->ParticleXCoord[ParticleNo]	= XCoord ;
				h_pDSMC->ParticleYCoord[ParticleNo]	= YCoord ;
				h_pDSMC->ParticleZCoord[ParticleNo]	= ZCoord ;
				
				RandVelocity( &h_pDSMC->ParticleXVel[ParticleNo] , &buffer , MostProbableSpeed ) ;
				RandVelocity( &h_pDSMC->ParticleYVel[ParticleNo] , &buffer , MostProbableSpeed ) ;
				RandVelocity( &h_pDSMC->ParticleZVel[ParticleNo] , &buffer , MostProbableSpeed ) ;
				
				h_pDSMC->ParticleXVel[ParticleNo]	+= h_pDomain->XVel ;
				h_pDSMC->ParticleYVel[ParticleNo]	+= h_pDomain->YVel ;
				h_pDSMC->ParticleZVel[ParticleNo]	+= h_pDomain->ZVel ;
				
				// Set rotational energy of particle
				if ( h_Species[j].RotDOF > 0. )
					h_pDSMC->ParticleRotation[ParticleNo] 	= RotationalEnergy( Temp , h_Species[j].RotDOF ) ; 
					
				// Set vibrational energy of particle
				if ( h_pDSMC->EffVibDOF[i][j]  > 0. ) {
					 
//  				h_pDSMC->EffVibDOF[i][j] = 2.*(h_Species[j].VibTemp/Temp)/(exp(h_Species[j].VibTemp/Temp)-1.) ;
					h_pDSMC->ParticleVibration[ParticleNo]	= VibrationalEnergy( Temp , &h_pDSMC->ParticleVibLevel[ParticleNo] , &h_pDSMC->ParticleEffTemp[ParticleNo] , h_pDSMC->EffVibDOF[i][j] , h_pDomain , h_Species[j].VibTemp ) ;
					
				}
				// Validation for vibration energy with Bird's DSMC0V.FOR code (Debug).
				//h_pDSMC->ParticleRotation[ParticleNo] 	= RotationalEnergy( 0. , h_Species[j].RotDOF ) ; 
    		//h_pDSMC->ParticleVibration[ParticleNo]	= VibrationalEnergy( 0. , &h_pDSMC->ParticleVibLevel[ParticleNo] , &h_pDSMC->ParticleEffTemp[ParticleNo] , h_pDSMC->EffVibDOF[i][j] , h_pDomain , h_Species[j].VibTemp ) ;
			
			
				ParticleNo++ ;
				
				if ( ParticleNo >= MaxParticleNum ){
					cout << "Particle array size has reached the limit. Maximum_Particle_Number is " << MaxParticleNum << '.' << endl ;
					exit(1) ;
				}
			}
		}	
	}
	
	h_pDomain->ParticleNum	= ParticleNo ;
	
	// Debug.
	//cout << "Initial ParticleNum: " << h_pDomain->ParticleNum << '\n' ;
}

//==============================================================================================================
//==============================================================================================================

void InitializationTwoStage(	DSMC_DOMAIN		*h_pDomain ,
				DSMC_NODE		*h_Node ,
				DSMC_CELL		*h_Cell ,
				DSMC_INLET		*h_Inlet ,
				DSMC_WALLTYPE		*h_WallType ,
				DSMC_SPECIES		*h_Species ,
				DSMC_SURFACE		*h_Surface ,
				DSMC_DSMC		*h_pDSMC ,
				CELLMAPPING		*h_pMapping , 
				ofstream		&OutputDebug ){
				
	int		MPISize , MPIMyID ;
	
	MPI_Comm_size( MPI_COMM_WORLD , &MPISize ) ;
	MPI_Comm_rank( MPI_COMM_WORLD , &MPIMyID ) ;
	

	AdjustTimestepWeighting( h_pDomain , h_Cell , h_Species , h_pDSMC , h_pMapping , OutputDebug ) ;

	
	CreateInitialParticle( h_pDomain , h_Node , h_Cell , h_Species , h_pDSMC , h_pMapping , OutputDebug ) ;
}

//==============================================================================================================

void AdjustTimestepWeighting(	DSMC_DOMAIN		*h_pDomain ,
				DSMC_CELL		*h_Cell ,
				DSMC_SPECIES		*h_Species ,
				DSMC_DSMC		*h_pDSMC ,
				CELLMAPPING		*h_pMapping , 
				ofstream		&OutputDebug ){

	int		CellNum , Dimension , SubCellNum ;
	double		Timestep , ParticleWeighting , XVel , YVel , ZVel , Temp , Mass , CrossSection ; 
	double		Coeff , MinX , ParticleNum ;
	
	Dimension	= h_pDomain->Dimension ;
	CellNum		= h_pDomain->CellNum ;
	XVel		= h_pDomain->XVel ;
	YVel		= h_pDomain->YVel ;
	ZVel		= h_pDomain->ZVel ;
	Temp		= h_pDomain->Temp ;
	Mass		= h_Species[0].Mass ;
	
	CrossSection	= 0. ;
	Coeff		= 0. ;

	
	for ( int i=0 ; i<CellNum ; i++ ){
		/*for ( int j=0 ; j<h_pDomain->SpeciesGroupNum ; j++ ){
			for ( int k=0 ; k<h_pDomain->SpeciesGroupNum ; k++ ){
				CrossSection	= 0.25*PI*pow((h_Species[j].Diameter+h_Species[k].Diameter),2.) ;
				
				h_pDSMC->RemainderCollisionPair[i][j][k]	= Randn() ;
				h_pDSMC->MaxCrossSectionSpeed[i][j][k]		= CrossSection*300.*sqrt(Temp/300.) ;
			}
		}*/
	
		ParticleNum	= 10. ;
	
	
		if ( h_Cell[i].MeanCollSpacingMeanFreePath >= 1. ){
			//SubCellNum	= h_Cell[i].CharacteristicLength/(0.5*h_Cell[i].InitMeanFreePath) ;
			//SubCellNum	= h_Cell[i].CharacteristicLength/h_Cell[i].InitMeanFreePath ;
			SubCellNum	= h_Cell[i].CharacteristicLength/h_Cell[i].MeanFreePath ;
			ParticleNum	= SubCellNum*SubCellNum*2 + 2 ;
		}
		
		if ( ParticleNum < 10. ) ParticleNum = 10. ;
			
		if ( h_Cell[i].AveParticleNum < ParticleNum ){
			h_Cell[i].Timestep	= h_Cell[i].Timestep/(ParticleNum/h_Cell[i].AveParticleNum) ;
			h_Cell[i].Weighting	= h_Cell[i].Weighting/(ParticleNum/h_Cell[i].AveParticleNum) ;
		}
	}
}

//==============================================================================================================
//==============================================================================================================
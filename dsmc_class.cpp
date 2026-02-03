#include <mpi.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include "dsmc_class.h"

using namespace std ;

//==============================================================================================================

DSMC_NODE::DSMC_NODE(){
	Id		= 0 ;
	XCoord		= 0. ;
	YCoord		= 0. ;
	ZCoord		= 0. ;
}


void DSMC_NODE::Dump( int No ){
	cout << setw(10) << No << setw(10) << Id << setw(16) << XCoord << setw(16) << YCoord << setw(16) << ZCoord << '\n' ;
}

//==============================================================================================================

DSMC_CELL::DSMC_CELL(){
	Id		= 0 ;
	LocalId		= 0 ;
	Type		= -1 ;
	for ( int i=0 ; i<8 ; i++ ) Node[i]	= 0 ;
	for ( int i=0 ; i<6 ; i++ ) Neighbor[i]	= 0 ;
	for ( int i=0 ; i<6 ; i++ ) Surface[i]	= -1 ;
	ProcessorNo	= -1 ;
	XCenter		= 0. ;
	YCenter		= 0. ;
	ZCenter		= 0. ;
	Volume		= 0. ;
	Timestep	= 0. ;
	Weighting	= 0. ;
	CharacteristicLength	= 0. ;
	for ( int i=0 ; i<6 ; i++ ){
		FaceFA[i]= 0. ;
		FaceFB[i]= 0. ;
		FaceFC[i]= 0. ;
		FaceFD[i]= 0. ;
	}
	MaxXCoord	= 0. ;
	MinXCoord	= 0. ;
	MaxYCoord	= 0. ;
	MinYCoord	= 0. ;
	MaxZCoord	= 0. ;
	MinZCoord	= 0. ;
	
	SelectNum	= 0. ;
	CollNum		= 0. ;
	CollDistance	= 0. ;
	CollRatio	= 0. ;
	MeanFreePath	= 0. ;
	MeanCollisionTime = 0. ;
	Temp		= 0. ;
	Speed		= 0. ;
	InitMeanFreePath= 0. ;
	InitTemp	= 0. ;
	InitSpeed	= 0. ;
	InitTimestep	= 1. ;
	InitWeighting	= 1. ;
	AveParticleNum	= 0. ;
	MeanCollSpacingMeanFreePath = 0. ;
	SubcellNum	= 0 ;
}


void DSMC_CELL::Dump( int No ){
	cout << setw(10) << No << setw(10) << Id << setw(4) << Type << '\n' ;
	cout << "---------------------------------\n" ;
	cout << setw(16) << "Nodes: " << setw(10) << Node[0] << setw(10) << Node[1] << setw(10) << Node[2] << setw(10) << Node[3] << setw(10) << Node[4] << setw(10) << Node[5] << setw(10) << Node[6] << setw(10) << Node[7] << '\n' ;
	cout << setw(16) << "Neighbors: " << setw(10) << Neighbor[0] << setw(10) << Neighbor[1] << setw(10) << Neighbor[2] << setw(10) << Neighbor[3] << setw(10) << Neighbor[4] << setw(10) << Neighbor[5] << '\n' ;
	cout << "-------------\n" ;
	cout << setw(16) << "Center, V: " << setw(14) << XCenter << setw(14) << YCenter << setw(14) << ZCenter << setw(14) << Volume << setw(14) << CharacteristicLength << '\n' ;
	cout << setw(16) << "T, W: " << setw(14) << Timestep << setw(14) << Weighting << '\n' ;
	cout << setw(16) << "Face-1: " << setw(14) << FaceFA[0] << setw(14) << FaceFB[0] << setw(14) << FaceFC[0] << setw(14) << FaceFD[0] << '\n' ;
	cout << setw(16) << "Face-2: " << setw(14) << FaceFA[1] << setw(14) << FaceFB[1] << setw(14) << FaceFC[1] << setw(14) << FaceFD[1] << '\n' ;
	cout << setw(16) << "Face-3: " << setw(14) << FaceFA[2] << setw(14) << FaceFB[2] << setw(14) << FaceFC[2] << setw(14) << FaceFD[2] << '\n' ;
	cout << setw(16) << "Face-4: " << setw(14) << FaceFA[3] << setw(14) << FaceFB[3] << setw(14) << FaceFC[3] << setw(14) << FaceFD[3] << '\n' ;
	cout << setw(16) << "Face-5: " << setw(14) << FaceFA[4] << setw(14) << FaceFB[4] << setw(14) << FaceFC[4] << setw(14) << FaceFD[4] << '\n' ;
	cout << setw(16) << "Face-6: " << setw(14) << FaceFA[5] << setw(14) << FaceFB[5] << setw(14) << FaceFC[5] << setw(14) << FaceFD[5] << '\n' ;	
	cout << setw(12) << "MaxXCoord" << setw(12) << MaxXCoord << setw(12) << "MinXCoord" << setw(12) << MinXCoord << setw(12) << "MaxYCoord" << setw(12) << MaxYCoord << setw(12) << "MinYCoord" << setw(12) << MinYCoord << setw(12) << "MaxZCoord" << setw(12) << MaxZCoord << setw(12) << "MinZCoord" << setw(12) << MinZCoord << '\n' ;
}

//==============================================================================================================

DSMC_INLET::DSMC_INLET(){
	Id		= 0 ;
	CellNo		= 0 ;
	FaceNo		= -1 ;
	NodeNum		= 0 ;
	for ( int i=0 ; i<4 ; i++ ) Node[i] = 0 ;
	SpeciesNo	= -1 ;
	Area		= 0. ;
	XVel		= 0. ;
	YVel		= 0. ;
	ZVel		= 0. ;
	Temp		= 0. ;
	CosineLawCoef	= -1. ;
	NumDen		= NULL ;
}


DSMC_INLET::~DSMC_INLET(){
	//delete [] NumDen ;	
}


void DSMC_INLET::AllocateMemory( int SpeciesNum ){
	NumDen	= new double[SpeciesNum] ;
	
	for ( int i=0 ; i<SpeciesNum ; i++ )
		NumDen[i]	= 0. ;
}


void DSMC_INLET::DeleteMemory(){
	delete [] NumDen ;	
}


void DSMC_INLET::Dump( int No ){
	cout << "No: " << No << ", Id: " << Id << ", CellNo: " << CellNo << ", FaceNo: " << FaceNo << ", NodeNum: " << NodeNum << '\n' ;
	for ( int i=0 ; i<NodeNum ; i++ ) cout << setw(10) << Node[i] ;
	cout << '\n' ;
	cout << "Species No: " << SpeciesNo << '\n' ;
	cout << setw(12) << Area << setw(12) << XVel << setw(12) << YVel << setw(12) << ZVel << setw(12) << NumDen << setw(12) << Temp << setw(12) << CosineLawCoef << '\n' ;
}

//==============================================================================================================

DSMC_WALLTYPE::DSMC_WALLTYPE(){
	Id		= 0 ;
	WallNo		= 0 ;
	Type		= 0 ;
	Temp		= 0. ;
	XVel		= 0. ;
	YVel		= 0. ;
	ZVel		= 0. ;
	DiffRatio	= 0. ;
	StickingCoef	= 0. ;
	AlphaN	= 0. ;
	SigmaT  = 0. ;
	EintCoef = 0. ;
	
}


void DSMC_WALLTYPE::Dump( int No ){
	cout << "No: " << No << ", Id: " << Id << '\n' ;
	cout << setw(5) << WallNo << setw(4) << Type << setw(10) << Temp << setw(12) << XVel << setw(12) << YVel << setw(12) << ZVel << setw(12) << DiffRatio << setw(12) << StickingCoef 
			<< setw(12) << AlphaN << setw(12) << SigmaT << setw(12) << EintCoef << '\n' ;
} 

//==============================================================================================================

DSMC_SPECIES::DSMC_SPECIES(){
	Id		= 0 ;
	GroupNo		= 0 ;
	RotModel	= 0 ;
	VibMode		= 0 ;
	VibModel	= 0 ; 
	
	Fraction	= 0. ;
	Diameter	= 0. ;
	RefTemp		= 0. ;
	VisTempIndex	= 0. ;
	VSSParameter	= 0. ;
	Mass		= 0. ;
	
	RotDOF		= 0. ;
	RotRelaxationNum= NULL ;
	
	VibDOF		= 0. ;
	VibConstant1	= NULL ;
	VibConstant2	= NULL ;
	VibTemp		= 0. ;
	DisTemp   = 0. ;//M.C.
	
	ReactionColClass = NULL ;//M.C.
	
}


DSMC_SPECIES::~DSMC_SPECIES(){
	delete [] RotRelaxationNum ;
	delete [] VibConstant1 ;
	delete [] VibConstant2 ;
	delete [] ReactionColClass ;	
}


void DSMC_SPECIES::AllocateMemory( int SpeciesNum ){
	RotRelaxationNum	= new double[SpeciesNum] ;
	VibConstant1		= new double[SpeciesNum] ;
	VibConstant2		= new double[SpeciesNum] ;
	ReactionColClass		= new int[SpeciesNum] ;
	
	for ( int i=0 ; i<SpeciesNum ; i++ ){
		RotRelaxationNum[i]	= 0. ;
		VibConstant1[i]		= 0. ;
		VibConstant2[i]		= 0. ;
		ReactionColClass[i]		= -1 ;
	}
}


void DSMC_SPECIES::Dump( int No ){
	cout << "No: " << No << ", Id: " << Id << '\n' ; 
	cout << setw(24) << "GroupNo" 		<< setw(14) << GroupNo << '\n' ; 
	cout << setw(24) << "RotModel" 		<< setw(14) << RotModel << '\n' ;  
	cout << setw(24) << "VibMode"		<< setw(14) << VibMode << '\n' ; 
	cout << setw(24) << "VibModel"		<< setw(14) << VibModel << '\n' ; 
	
	cout << setw(24) << "Fraction" 		<< setw(14) << Fraction << '\n' ; 
	cout << setw(24) << "Diameter" 		<< setw(14) << Diameter << '\n' ; 
	cout << setw(24) << "RefTemp" 		<< setw(14) << RefTemp << '\n' ; 
	cout << setw(24) << "VisTempIndex"	<< setw(14) << VisTempIndex << '\n' ; 
	cout << setw(24) << "VSSParameter"	<< setw(14) << VSSParameter << '\n' ; 
	cout << setw(24) << "Mass"		<< setw(14) << Mass << '\n' ; 
	
	cout << setw(24) << "RotDOF" 		<< setw(14) << RotDOF << '\n' ; 
	cout << setw(24) << "RotRelaxationNum" 	<< setw(14) << RotRelaxationNum << '\n' ; 
	
	cout << setw(24) << "VibDOF" 		<< setw(14) << VibDOF << '\n' ; 
	cout << setw(24) << "VibConstant1" 	<< setw(14) << VibConstant1 << '\n' ; 
	cout << setw(24) << "VibConstant2" 	<< setw(14) << VibConstant2 << '\n' ; 
	cout << setw(24) << "VibTemp" 		<< setw(14) << VibTemp << '\n' ;
	cout << setw(24) << "DisTemp" 		<< setw(14) << DisTemp << '\n' ;
	cout << setw(24) << "ReactionColClass" 		<< setw(14) << ReactionColClass << '\n' ; 
}

//==============================================================================================================

DSMC_MULTISPECIES::DSMC_MULTISPECIES(){
	CrossSection	= 0. ;
	RefTemp		= 0. ;
	VisTempIndex	= 0. ;
	VSSParameter	= 0. ;
	ReduceMass	= 0. ;
	GammaValue	= 0. ;
}		

//==============================================================================================================
// Start Ming-Chung Lo

DSMC_CHEMICALTCE::DSMC_CHEMICALTCE(){
		 	
		ReactionClassNum	= 0 ;		 		
		ReactionType	= NULL ;
		
//	ReactionColClass	= NULL ;

		Pre1Species	= NULL ; 
		Pre2Species = NULL ;
		Pre3Species = NULL ; 
		Post1Species	= NULL ;
		Post2Species	= NULL ;
		Post3Species	= NULL ;
		
		ArrheniusConstant	= NULL ;
		ArrheniusTempExp = NULL ;
		ArrheniusActiveEnergy = NULL ;
		HeatofReaction = NULL ;
		
		SuccessReaction = NULL ;	
}

DSMC_CHEMICALTCE::~DSMC_CHEMICALTCE(){
				
		delete [] ReactionType ;
		
		delete [] Pre1Species ;
		delete [] Pre2Species ;
		delete [] Pre3Species ;
		delete [] Post1Species ;
		delete [] Post2Species ;
		delete [] Post3Species ;
	
		delete [] ArrheniusConstant ;
		delete [] ArrheniusTempExp ;
		delete [] ArrheniusActiveEnergy ;
		delete [] HeatofReaction ;	
		
		delete [] SuccessReaction ;
}


void DSMC_CHEMICALTCE::AllocateMemory ( int SpeciesNum ){
		
	  ReactionType	= new int[SpeciesNum] ;
		
		Pre1Species	= new int[SpeciesNum] ;
		Pre2Species	= new int[SpeciesNum] ;
		Pre3Species	= new int[SpeciesNum] ;
		Post1Species	= new int[SpeciesNum] ;
		Post2Species	= new int[SpeciesNum] ;
		Post3Species	= new int[SpeciesNum] ;

		ArrheniusConstant	= new double[SpeciesNum] ;
		ArrheniusTempExp	= new double[SpeciesNum] ;
		ArrheniusActiveEnergy	= new double[SpeciesNum] ;
		HeatofReaction	= new double[SpeciesNum] ;	
	
		SuccessReaction = new int[SpeciesNum] ;
		
		for ( int i=0 ; i<SpeciesNum ; i++ ){
			
			ReactionType[i]= -1 ;
							
			Pre1Species[i]= -1 ;
			Pre2Species[i]= -1 ;
			Pre3Species[i]= -1 ;
			Post1Species[i]= -1 ;
			Post2Species[i]= -1 ;
			Post3Species[i]= -1 ;

			ArrheniusConstant[i]= -1. ;
			ArrheniusTempExp[i]= -1. ;
			ArrheniusActiveEnergy[i]= -1. ;
			HeatofReaction[i]= -1. ;
			
			SuccessReaction[i]= 0. ;
			
		}
}

void DSMC_CHEMICALTCE::Dump( int No ){
	
	cout << "No: " << No  << '\n' ; 
	cout << setw(24) << "ReactionClassNum" 		<< setw(14) << ReactionClassNum << '\n' ; 
	cout << setw(24) << "ReactionType" 		<< setw(14) << ReactionType << '\n' ; 
	 
	cout << setw(24) << "Pre1Species"		<< setw(14) << Pre1Species << '\n' ; 
	cout << setw(24) << "Pre2Species"		<< setw(14) << Pre2Species << '\n' ; 	
	cout << setw(24) << "Pre3Species" 		<< setw(14) << Pre3Species << '\n' ; 
	cout << setw(24) << "Post1Species" 		<< setw(14) << Post1Species << '\n' ; 
	cout << setw(24) << "Post2Species" 		<< setw(14) << Post2Species << '\n' ; 
	cout << setw(24) << "Post3Species"	<< setw(14) << Post3Species << '\n' ;
	 
	cout << setw(24) << "ArrheniusConstant"	<< setw(14) << ArrheniusConstant << '\n' ; 
	cout << setw(24) << "ArrheniusTempExp"		<< setw(14) << ArrheniusTempExp << '\n' ; 	
	cout << setw(24) << "ArrheniusActiveEnergy" 		<< setw(14) << ArrheniusActiveEnergy << '\n' ; 	
	cout << setw(24) << "HeatofReaction" 		<< setw(14) << HeatofReaction << '\n' ; 
	
	cout << setw(24) << "SuccessReaction" 	<< setw(14) << SuccessReaction << '\n' ; 

}

// End Ming-Chung Lo
//==============================================================================================================

DSMC_SURFACE::DSMC_SURFACE(){
	WallNo	= 0 ;
	CellNo	= 0 ;
	FaceNo	= -1 ;
	XCenter	= 0. ;
	YCenter	= 0. ;
	ZCenter	= 0. ;
	Area	= 0. ;
}


void DSMC_SURFACE::Dump( int No ){
	cout << "No: " << No << ", Wall No: " << WallNo << ", CellNo: " << CellNo << ", FaceNo: " << FaceNo << '\n' ;
	cout << setw(12) << "Area" << setw(12) << Area << setw(12) << "XCenter" << setw(12) << XCenter << setw(12) << "YCenter" << setw(12) << YCenter << setw(12) << "ZCenter" << setw(12) << ZCenter << '\n' ;
}

//==============================================================================================================

DSMC_DSMC::DSMC_DSMC(){
	ParticleCellNo		= NULL ;
	ParticleSpeciesNo	= NULL ;
	ParticleLastCollide	= NULL ;
	ParticleBeforeCellNo	= NULL ;
	ParticleXCoord		= NULL ;
	ParticleYCoord		= NULL ;
	ParticleZCoord		= NULL ;
	ParticleXVel		= NULL ;
	ParticleYVel		= NULL ;
	ParticleZVel		= NULL ;
	ParticleRotation	= NULL ;
	ParticleVibration	= NULL ;
	ParticleEffTemp		= NULL ;
	ParticleTimestep	= NULL ;
	ParticleVibLevel	= NULL ;

	IndexCell1		= NULL ;
	IndexCell2		= NULL ;
	IndexParticle		= NULL ;
	ReactionIndex   = NULL ;
	
	SampleParticleNum	= NULL ;
	SampleXVel		= NULL ;
	SampleYVel		= NULL ;
	SampleZVel		= NULL ;
	SampleXVelSq		= NULL ;
	SampleYVelSq		= NULL ;
	SampleZVelSq		= NULL ;
	SampleRotation		= NULL ;
	SampleVibration		= NULL ;
	SampleVibLevel		= NULL ;
	SampleCellVibLevel		= NULL ;
	SampleVibGroundNum	= NULL ;
	SampleVibLevelNum	= NULL ;
	
	SampleSurfaceParticleNum	= NULL ;
	SampleSurfaceStickingParticleNum= NULL ;
	SampleSurfaceInNormMomentum	= NULL ;
	SampleSurfaceReNormMomentum	= NULL ;
	SampleSurfaceInXMomentum	= NULL ;
	SampleSurfaceReXMomentum	= NULL ;
	SampleSurfaceInYMomentum	= NULL ;
	SampleSurfaceReYMomentum	= NULL ;
	SampleSurfaceInZMomentum	= NULL ;
	SampleSurfaceReZMomentum	= NULL ; 
	SampleSurfaceInTransEng		= NULL ;
	SampleSurfaceReTransEng		= NULL ;
	SampleSurfaceInRotEng		= NULL ;
	SampleSurfaceReRotEng		= NULL ;
	SampleSurfaceInVibEng		= NULL ;
	SampleSurfaceReVibEng		= NULL ;
	
	MaxCrossSectionSpeed	= NULL ;
	TotalCrossSectionSpeed	= NULL ;
	RemainderCollisionPair	= NULL ;
	
	InletEnterNum		= NULL ;
	InletRemainderEnterNum	= NULL ;
	
	EffVibDOF 		= NULL ;
	
	SamplingTimeEnd		= NULL ;
	SamplingTimeInit	= NULL ;
	
	CollisionNum		= NULL ;
}


DSMC_DSMC::~DSMC_DSMC(){
	delete [] ParticleCellNo ;
	delete [] ParticleSpeciesNo ;
	delete [] ParticleLastCollide ;
	delete [] ParticleBeforeCellNo ;
	delete [] ParticleXCoord ;
	delete [] ParticleYCoord ;
	delete [] ParticleZCoord ;
	delete [] ParticleXVel ;
	delete [] ParticleYVel ;
	delete [] ParticleZVel ;
	delete [] ParticleRotation ;
	delete [] ParticleVibration ;
	delete [] ParticleEffTemp ;
	delete [] ParticleTimestep ;
	delete [] ParticleVibLevel ;

	delete [] IndexCell1 ;
	delete [] IndexCell2 ;
	delete [] IndexParticle ;
	delete [] ReactionIndex	;
	
	
	delete [] SampleParticleNum ;
	delete [] SampleXVel ;
	delete [] SampleYVel ;
	delete [] SampleZVel ;
	delete [] SampleXVelSq ;
	delete [] SampleYVelSq ;
	delete [] SampleZVelSq ;
	delete [] SampleRotation ;
	delete [] SampleVibration ;
	delete [] SampleVibLevel ;
	delete [] SampleCellVibLevel ;
	delete [] SampleVibGroundNum ;
	delete [] SampleVibLevelNum ;
	
	delete [] SampleSurfaceParticleNum ;
	delete [] SampleSurfaceStickingParticleNum ;
	delete [] SampleSurfaceInNormMomentum ;
	delete [] SampleSurfaceReNormMomentum ;
	delete [] SampleSurfaceInXMomentum ;
	delete [] SampleSurfaceReXMomentum ;
	delete [] SampleSurfaceInYMomentum ;
	delete [] SampleSurfaceReYMomentum ;
	delete [] SampleSurfaceInZMomentum ;
	delete [] SampleSurfaceReZMomentum ; 
	delete [] SampleSurfaceInTransEng ;
	delete [] SampleSurfaceReTransEng ;
	delete [] SampleSurfaceInRotEng ;
	delete [] SampleSurfaceReRotEng ;
	delete [] SampleSurfaceInVibEng ;
	delete [] SampleSurfaceReVibEng ;
	
	delete [] MaxCrossSectionSpeed ;
	delete [] TotalCrossSectionSpeed ;
	
	delete [] RemainderCollisionPair ;
	
	delete [] InletEnterNum ;
	delete [] InletRemainderEnterNum ;
	
	delete [] EffVibDOF ;
	
	delete [] SamplingTimeEnd ;
	delete [] SamplingTimeInit ;
	
	delete [] CollisionNum ;
}


void DSMC_DSMC::AllocateMemory( int MaxParticleNum , int CellNum , int WallFaceNum , int InletFaceNum , int SpeciesNum , int SpeciesGroupNum ){
	ParticleCellNo		= new int[MaxParticleNum] ;
	ParticleSpeciesNo	= new int[MaxParticleNum] ;
	ParticleLastCollide	= new int[MaxParticleNum] ;
	ParticleBeforeCellNo	= new int[MaxParticleNum] ;
	ParticleXCoord		= new double[MaxParticleNum] ;
	ParticleYCoord		= new double[MaxParticleNum] ;
	ParticleZCoord		= new double[MaxParticleNum] ;
	ParticleXVel		= new double[MaxParticleNum] ;
	ParticleYVel		= new double[MaxParticleNum] ;
	ParticleZVel		= new double[MaxParticleNum] ;
	ParticleRotation	= new double[MaxParticleNum] ;
	ParticleVibration	= new double[MaxParticleNum] ;
	ParticleEffTemp		= new double[MaxParticleNum] ;
	ParticleTimestep	= new double[MaxParticleNum] ;
	ParticleVibLevel	= new int[MaxParticleNum] ;

	IndexCell1		= new int*[SpeciesGroupNum] ;
	IndexCell2		= new int*[SpeciesGroupNum] ;
	IndexParticle		= new int[MaxParticleNum] ;
	ReactionIndex	= new int[MaxParticleNum] ;
		
	for ( int i=0 ; i<SpeciesGroupNum ; i++ ){
		IndexCell1[i]		= new int[CellNum] ;
		IndexCell2[i]		= new int[CellNum] ;
	}
	
	
	SampleParticleNum	= new double*[SpeciesNum] ;
	SampleXVel		= new double*[SpeciesNum] ;
	SampleYVel		= new double*[SpeciesNum] ;
	SampleZVel		= new double*[SpeciesNum] ;
	SampleXVelSq		= new double*[SpeciesNum] ;
	SampleYVelSq		= new double*[SpeciesNum] ;
	SampleZVelSq		= new double*[SpeciesNum] ;
	SampleRotation		= new double*[SpeciesNum] ;
	SampleVibration		= new double*[SpeciesNum] ;
	SampleVibGroundNum	= new double*[SpeciesNum] ;
	SampleVibLevelNum	= new double*[SpeciesNum] ;
	SampleVibLevel    = new double*[SpeciesNum] ;
	
	for ( int i=0 ; i<SpeciesNum ; i++ ){
		SampleParticleNum[i]	= new double[CellNum] ;
		SampleXVel[i]		= new double[CellNum] ;
		SampleYVel[i]		= new double[CellNum] ;
		SampleZVel[i]		= new double[CellNum] ;
		SampleXVelSq[i]		= new double[CellNum] ;
		SampleYVelSq[i]		= new double[CellNum] ;
		SampleZVelSq[i]		= new double[CellNum] ;
		SampleRotation[i]	= new double[CellNum] ;
		SampleVibration[i]	= new double[CellNum] ;
		SampleVibGroundNum[i]	= new double[CellNum] ;
		SampleVibLevelNum[i]	= new double[CellNum] ;
		SampleVibLevel[i]		= new double[CellNum] ;
	}
	
	
	
	SampleSurfaceParticleNum	= new double*[WallFaceNum] ;
	SampleSurfaceStickingParticleNum= new double*[WallFaceNum] ;
	SampleSurfaceInNormMomentum	= new double*[WallFaceNum] ;
	SampleSurfaceReNormMomentum	= new double*[WallFaceNum] ;
	SampleSurfaceInXMomentum	= new double*[WallFaceNum] ;
	SampleSurfaceReXMomentum	= new double*[WallFaceNum] ;
	SampleSurfaceInYMomentum	= new double*[WallFaceNum] ;
	SampleSurfaceReYMomentum	= new double*[WallFaceNum] ;
	SampleSurfaceInZMomentum	= new double*[WallFaceNum] ;
	SampleSurfaceReZMomentum	= new double*[WallFaceNum] ; 
	SampleSurfaceInTransEng		= new double*[WallFaceNum] ;
	SampleSurfaceReTransEng		= new double*[WallFaceNum] ;
	SampleSurfaceInRotEng		= new double*[WallFaceNum] ;
	SampleSurfaceReRotEng		= new double*[WallFaceNum] ;
	SampleSurfaceInVibEng		= new double*[WallFaceNum] ;
	SampleSurfaceReVibEng		= new double*[WallFaceNum] ;
	for ( int i=0 ; i<WallFaceNum ; i++ ){
		SampleSurfaceParticleNum[i]		= new double[SpeciesNum] ;
		SampleSurfaceStickingParticleNum[i]	= new double[SpeciesNum] ;
		SampleSurfaceInNormMomentum[i]		= new double[SpeciesNum] ;
		SampleSurfaceReNormMomentum[i]		= new double[SpeciesNum] ;
		SampleSurfaceInXMomentum[i]		= new double[SpeciesNum] ;
		SampleSurfaceReXMomentum[i]		= new double[SpeciesNum] ;
		SampleSurfaceInYMomentum[i]		= new double[SpeciesNum] ;
		SampleSurfaceReYMomentum[i]		= new double[SpeciesNum] ;
		SampleSurfaceInZMomentum[i]		= new double[SpeciesNum] ;
		SampleSurfaceReZMomentum[i]		= new double[SpeciesNum] ; 
		SampleSurfaceInTransEng[i]		= new double[SpeciesNum] ;
		SampleSurfaceReTransEng[i]		= new double[SpeciesNum] ;
		SampleSurfaceInRotEng[i]		= new double[SpeciesNum] ;
		SampleSurfaceReRotEng[i]		= new double[SpeciesNum] ;
		SampleSurfaceInVibEng[i]		= new double[SpeciesNum] ;
		SampleSurfaceReVibEng[i]		= new double[SpeciesNum] ;
	}
	
	
	MaxCrossSectionSpeed	= new double**[CellNum] ;
	TotalCrossSectionSpeed	= new double**[CellNum] ;
	RemainderCollisionPair	= new double**[CellNum] ;
	for ( int i=0 ; i<CellNum ; i++ ){
		MaxCrossSectionSpeed[i]		= new double*[SpeciesGroupNum] ;
		TotalCrossSectionSpeed[i]		= new double*[SpeciesGroupNum] ;
		RemainderCollisionPair[i]	= new double*[SpeciesGroupNum] ;
		
		for ( int j=0 ; j<SpeciesGroupNum ; j++ ){
			MaxCrossSectionSpeed[i][j]	= new double[SpeciesGroupNum] ;
			TotalCrossSectionSpeed[i][j]	= new double[SpeciesGroupNum] ;
			RemainderCollisionPair[i][j]	= new double[SpeciesGroupNum] ;	
		}
	}
	
	InletEnterNum		= new double*[InletFaceNum] ;
	InletRemainderEnterNum	= new double*[InletFaceNum] ;
	for ( int i=0 ; i<InletFaceNum ; i++ ){
		InletEnterNum[i]		= new double[SpeciesNum] ;
		InletRemainderEnterNum[i]	= new double[SpeciesNum] ;
	}
	
	EffVibDOF		= new double*[CellNum] ;
	for ( int i=0 ; i<CellNum ; i++ )
		EffVibDOF[i]	= new double[SpeciesNum] ;
		
	SamplingTimeEnd		= new double[CellNum] ;
	SamplingTimeInit	= new double[CellNum] ;
	
	SampleCellVibLevel	= new double[CellNum] ;
	
	CollisionNum		= new double*[SpeciesNum] ;
	for ( int i=0 ; i<SpeciesNum ; i++ )
		CollisionNum[i]	= new double[SpeciesNum] ;
}


void DSMC_DSMC::DeleteMemory( int CellNum , int WallFaceNum , int InletFaceNum , int SpeciesNum , int SpeciesGroupNum ){
	for ( int i=0 ; i<SpeciesGroupNum ; i++ ){
		delete [] IndexCell1[i] ;
		delete [] IndexCell2[i] ;
	}

	for ( int i=0 ; i<SpeciesNum ; i++ ){
		delete [] SampleParticleNum[i] ;
		delete [] SampleXVel[i] ;
		delete [] SampleYVel[i] ;
		delete [] SampleZVel[i] ;
		delete [] SampleXVelSq[i] ;
		delete [] SampleYVelSq[i] ;
		delete [] SampleZVelSq[i] ;
		delete [] SampleRotation[i] ;
		delete [] SampleVibration[i] ;
		delete [] SampleVibGroundNum[i] ;
		delete [] SampleVibLevelNum[i] ;
		delete [] SampleVibLevel[i] ;
	}
	
	for ( int i=0 ; i<WallFaceNum ; i++ ){
		delete [] SampleSurfaceParticleNum[i] ;
		delete [] SampleSurfaceStickingParticleNum[i] ;
		delete [] SampleSurfaceInNormMomentum[i] ;
		delete [] SampleSurfaceReNormMomentum[i] ;
		delete [] SampleSurfaceInXMomentum[i] ;
		delete [] SampleSurfaceReXMomentum[i] ;
		delete [] SampleSurfaceInYMomentum[i] ;
		delete [] SampleSurfaceReYMomentum[i] ;
		delete [] SampleSurfaceInZMomentum[i] ;
		delete [] SampleSurfaceReZMomentum[i] ; 
		delete [] SampleSurfaceInTransEng[i] ;
		delete [] SampleSurfaceReTransEng[i] ;
		delete [] SampleSurfaceInRotEng[i] ;
		delete [] SampleSurfaceReRotEng[i] ;
		delete [] SampleSurfaceInVibEng[i] ;
		delete [] SampleSurfaceReVibEng[i] ;
	}
	
	for ( int i=0 ; i<CellNum ; i++ ){
		for ( int j=0 ; j<SpeciesGroupNum ; j++ ){
			delete [] MaxCrossSectionSpeed[i][j] ;
			delete [] TotalCrossSectionSpeed[i][j] ;
			delete [] RemainderCollisionPair[i][j] ;
		}
	}
	for ( int i=0 ; i<CellNum ; i++ ){
		delete [] MaxCrossSectionSpeed[i] ;
		delete [] TotalCrossSectionSpeed[i] ;
		delete [] RemainderCollisionPair[i] ;
	}
	
	for ( int i=0 ; i<InletFaceNum ; i++ ){
		delete [] InletEnterNum[i] ;
		delete [] InletRemainderEnterNum[i] ;
	}
	
	for ( int i=0 ; i<CellNum ; i++ )
		delete [] EffVibDOF[i] ;
		
	for ( int i=0 ; i<SpeciesNum ; i++ )	
		delete [] CollisionNum[i] ;
}


void DSMC_DSMC::InitValue( int MaxParticleNum , int CellNum , int WallFaceNum , int InletFaceNum , int SpeciesNum , int SpeciesGroupNum ){
	for ( int i=0 ; i<MaxParticleNum ; i++ ){
		ParticleCellNo[i]	= 0 ;
		ParticleSpeciesNo[i]	= 0 ;
		ParticleLastCollide[i]	= -1 ;
		ParticleBeforeCellNo[i]	= 0 ;
		ParticleXCoord[i]	= 0. ;
		ParticleYCoord[i]	= 0. ;
		ParticleZCoord[i]	= 0. ;
		ParticleXVel[i]		= 0. ;
		ParticleYVel[i]		= 0. ;
		ParticleZVel[i]		= 0. ;
		ParticleRotation[i]	= 0. ;
		ParticleVibration[i]	= 0. ;
		ParticleEffTemp[i]	= 0. ;
		ParticleTimestep[i]	= 0. ;
		ParticleVibLevel[i]	= 0 ;
		
		IndexParticle[i]	= 0 ;
		ReactionIndex[i]	= 0 ;
	}


	for ( int i=0 ; i<SpeciesGroupNum ; i++ ){
		for ( int j=0 ; j<CellNum ; j++ ){
			IndexCell1[i][j]	= 0 ;
			IndexCell2[i][j]	= 0 ;
		}
	}
	
	for ( int i=0 ; i<SpeciesNum ; i++ ){
		for ( int j=0 ; j<CellNum ; j++ ){
			SampleParticleNum[i][j]	= 0. ;
			SampleXVel[i][j]	= 0. ;
			SampleYVel[i][j]	= 0. ;
			SampleZVel[i][j]	= 0. ;
			SampleXVelSq[i][j]	= 0. ;
			SampleYVelSq[i][j]	= 0. ;
			SampleZVelSq[i][j]	= 0. ;
			SampleRotation[i][j]	= 0. ;
			SampleVibration[i][j]	= 0. ;
			SampleVibGroundNum[i][j]= 0. ;
			SampleVibLevelNum[i][j]	= 0. ;
			SampleVibLevel[i][j]	= 0. ;
		}
	}
	
	for ( int i=0 ; i<CellNum ; i++ ){
		SampleCellVibLevel[i]	= 0. ;	
	}

	for ( int i=0 ; i<WallFaceNum ; i++ ){
		for ( int j=0 ; j<SpeciesNum ; j++ ){
			SampleSurfaceParticleNum[i][j]	= 0. ;
			SampleSurfaceStickingParticleNum[i][j]	= 0. ;
			SampleSurfaceInNormMomentum[i][j]= 0. ;
			SampleSurfaceReNormMomentum[i][j]= 0. ;
			SampleSurfaceInXMomentum[i][j]	= 0. ;
			SampleSurfaceReXMomentum[i][j]	= 0. ;
			SampleSurfaceInYMomentum[i][j]	= 0. ;
			SampleSurfaceReYMomentum[i][j]	= 0. ;
			SampleSurfaceInZMomentum[i][j]	= 0. ;
			SampleSurfaceReZMomentum[i][j]	= 0. ; 
			SampleSurfaceInTransEng[i][j]	= 0. ;
			SampleSurfaceReTransEng[i][j]	= 0. ;
			SampleSurfaceInRotEng[i][j]	= 0. ;
			SampleSurfaceReRotEng[i][j]	= 0. ;
			SampleSurfaceInVibEng[i][j]	= 0. ;
			SampleSurfaceReVibEng[i][j]	= 0. ;
		}
	}
	
	for ( int i=0 ; i<CellNum ; i++ ){
		for ( int j=0 ; j<SpeciesGroupNum ; j++ ){
			for ( int k=0 ; k<SpeciesGroupNum ; k++ ){
				MaxCrossSectionSpeed[i][j][k]	= 0. ;
				TotalCrossSectionSpeed[i][j][k]	= 0. ;
				RemainderCollisionPair[i][j][k]	= 0. ;
			}
		}
	}
	
	for ( int i=0 ; i<InletFaceNum ; i++ ){
		for ( int j=0 ; j<SpeciesNum ; j++ ){
			InletEnterNum[i][j]		= 0. ;
			InletRemainderEnterNum[i][j]	= 0. ;
		}
	}
	
	for ( int i=0 ; i<CellNum ; i++ ){
		for ( int j=0 ; j<SpeciesNum ; j++ )
			EffVibDOF[i][j]	= 0. ;
	}
	
	for ( int i=0 ; i<CellNum ; i++ ){
		SamplingTimeEnd[i]	= 0. ;
		SamplingTimeInit[i]	= 0. ;	
	}
	
	for ( int i=0 ; i<SpeciesNum ; i++ ){
		for ( int j=0 ; j<SpeciesNum ; j++ )
			CollisionNum[i][j]	= 0. ;	
	}
}


void DSMC_DSMC::ParticleDump( int No ){
	cout << setw(20) << "Particle No"	<< setw(15) << No ;
	cout << setw(20) << "Cell No"		<< setw(15) << *(ParticleCellNo+No) << '\n' ;
	cout <<	setw(20) << "Species No"	<< setw(15) << *(ParticleSpeciesNo+No) << '\n' ;
	cout <<	setw(20) << "Last Collide"	<< setw(15) << *(ParticleLastCollide+No) << '\n' ;
	cout << setw(20) << "Before Cell No"	<< setw(15) << *(ParticleBeforeCellNo+No) << '\n' ;
	cout <<	setw(20) << "XCoord"		<< setw(15) << *(ParticleXCoord+No) << '\n' ;
	cout <<	setw(20) << "YCoord"		<< setw(15) << *(ParticleYCoord+No) << '\n' ;
	cout <<	setw(20) << "ZCoord"		<< setw(15) << *(ParticleZCoord+No) << '\n' ;
	cout <<	setw(20) << "XVel"		<< setw(15) << *(ParticleXVel+No) << '\n' ;
	cout <<	setw(20) << "YVel"		<< setw(15) << *(ParticleYVel+No) << '\n' ;
	cout <<	setw(20) << "ZVel"		<< setw(15) << *(ParticleZVel+No) << '\n' ;
	cout <<	setw(20) << "Rotation"		<< setw(15) << *(ParticleRotation+No) << '\n' ;
	cout <<	setw(20) << "Vibration"		<< setw(15) << *(ParticleVibration+No) << '\n' ;
	cout <<	setw(20) << "Vibration"		<< setw(15) << *(ParticleEffTemp+No) << '\n' ;
	cout << setw(20) << "Timestep"		<< setw(15) << *(ParticleTimestep+No) << '\n' ;
	cout <<	setw(20) << "VibLevel"		<< setw(15) << *(ParticleVibLevel+No) << '\n' ;
}

//==============================================================================================================

DSMC_RESULT::DSMC_RESULT(){
	NumDensity		= NULL ;
	Density			= NULL ;
	XVel			= NULL ;
	YVel			= NULL ;
	ZVel			= NULL ;
	Temp			= NULL ;
	TransTemp		= NULL ;
	RotTemp			= NULL ;
	VibTemp			= NULL ;
	TransXTemp		= NULL ;
	TransYTemp		= NULL ;
	TransZTemp		= NULL ;
	AveParticleNum		= NULL ;
	MeanCollSpacingMeanFreePath = NULL ;

	NumDensitySpecies	= NULL ;
	DensitySpecies		= NULL ;
	XVelSpecies		= NULL ;
	YVelSpecies		= NULL ;
	ZVelSpecies		= NULL ;
	TempSpecies		= NULL ;
	TransTempSpecies	= NULL ;
	TransXTempSpecies	= NULL ;
	TransYTempSpecies	= NULL ;
	TransZTempSpecies	= NULL ;
	RotTempSpecies		= NULL ;
	VibTempSpecies		= NULL ;
	AveParticleNumSpecies = NULL ;
	
}


DSMC_RESULT::~DSMC_RESULT(){
	delete [] NumDensity ;
	delete [] Density ;
	delete [] XVel ;
	delete [] YVel ;
	delete [] ZVel ;
	delete [] Temp ;
	delete [] TransTemp ;
	delete [] RotTemp ;
	delete [] VibTemp ;
	delete [] TransXTemp ;
	delete [] TransYTemp ;
	delete [] TransZTemp ;
	delete [] AveParticleNum ;
	delete [] MeanCollSpacingMeanFreePath ;
	
	delete [] NumDensitySpecies ;        
	delete [] DensitySpecies ;
	delete [] XVelSpecies ;        
	delete [] YVelSpecies ;        
	delete [] ZVelSpecies ;        
	delete [] TempSpecies ;        
	delete [] TransTempSpecies ;        
	delete [] TransXTempSpecies ;        
	delete [] TransYTempSpecies ;        
	delete [] TransZTempSpecies ;        
	delete [] RotTempSpecies ;
	delete [] VibTempSpecies ;
	delete [] AveParticleNumSpecies  ;
}


void DSMC_RESULT::AllocateMemory( int CellNum , int SpeciesNum ){
	NumDensity		= new double[CellNum] ;
	Density			= new double[CellNum] ;
	XVel			= new double[CellNum] ;
	YVel			= new double[CellNum] ;
	ZVel			= new double[CellNum] ;
	Temp			= new double[CellNum] ;
	TransTemp		= new double[CellNum] ;
	RotTemp			= new double[CellNum] ;
	VibTemp			= new double[CellNum] ;
	TransXTemp		= new double[CellNum] ;
	TransYTemp		= new double[CellNum] ;
	TransZTemp		= new double[CellNum] ;
	AveParticleNum		= new double[CellNum] ;
	MeanCollSpacingMeanFreePath = new double[CellNum] ;
	
	NumDensitySpecies	= new double*[SpeciesNum] ;
	DensitySpecies		= new double*[SpeciesNum] ;
	XVelSpecies		= new double*[SpeciesNum] ;
	YVelSpecies		= new double*[SpeciesNum] ;
	ZVelSpecies		= new double*[SpeciesNum] ;
	TempSpecies		= new double*[SpeciesNum] ;
	TransTempSpecies	= new double*[SpeciesNum] ;
	TransXTempSpecies	= new double*[SpeciesNum] ;
	TransYTempSpecies	= new double*[SpeciesNum] ;
	TransZTempSpecies	= new double*[SpeciesNum] ;
	RotTempSpecies		= new double*[SpeciesNum] ;
	VibTempSpecies		= new double*[SpeciesNum] ;
	AveParticleNumSpecies = new double*[SpeciesNum] ;
	
	for ( int i=0 ; i<SpeciesNum ; i++ ){
		NumDensitySpecies[i]	= new double[CellNum] ;
		DensitySpecies[i]	= new double[CellNum] ;
		XVelSpecies[i]		= new double[CellNum] ;
		YVelSpecies[i]		= new double[CellNum] ;
		ZVelSpecies[i]		= new double[CellNum] ;
		TempSpecies[i]		= new double[CellNum] ;
		TransTempSpecies[i]	= new double[CellNum] ;
		TransXTempSpecies[i]	= new double[CellNum] ;
		TransYTempSpecies[i]	= new double[CellNum] ;
		TransZTempSpecies[i]	= new double[CellNum] ;
		RotTempSpecies[i]	= new double[CellNum] ;
		VibTempSpecies[i]	= new double[CellNum] ;
		AveParticleNumSpecies [i]	= new double[CellNum] ;
	}
}


void DSMC_RESULT::DeleteMemory( int SpeciesNum ){
	for ( int i=0 ; i<SpeciesNum ; i++ ){
		delete [] NumDensitySpecies[i] ;        
		delete [] DensitySpecies[i] ;
		delete [] XVelSpecies[i] ;        
		delete [] YVelSpecies[i] ;        
		delete [] ZVelSpecies[i] ;        
		delete [] TempSpecies[i] ;        
		delete [] TransTempSpecies[i] ;        
		delete [] TransXTempSpecies[i] ;        
		delete [] TransYTempSpecies[i] ;        
		delete [] TransZTempSpecies[i] ;        
		delete [] RotTempSpecies[i] ;
		delete [] VibTempSpecies[i] ;	
		delete [] AveParticleNumSpecies [i] ;		
	}
}


void DSMC_RESULT::InitValue( int CellNum , int SpeciesNum , double InitTemp ){
		for ( int i=0 ; i<CellNum ; i++ ){
			NumDensity[i]		= 0. ;
			Density[i]		= 0. ;
			XVel[i]			= 0. ;
			YVel[i]			= 0. ;
			ZVel[i]			= 0. ;
			Temp[i]			= InitTemp ;
			TransTemp[i]		= InitTemp ;
			RotTemp[i]		= 0. ;
			VibTemp[i]		= 0. ;
			TransXTemp[i]		= 0. ;
			TransYTemp[i]		= 0. ;
			TransZTemp[i]		= 0. ;
			AveParticleNum[i]	= 0. ;
			MeanCollSpacingMeanFreePath[i] = 0. ;
		}
	
		for ( int i=0 ; i<SpeciesNum ; i++ ){
			for ( int j=0 ; j<CellNum ; j++ ){
				NumDensitySpecies[i][j]	= 0. ;
				DensitySpecies[i][j]	= 0. ;
				XVelSpecies[i][j]	= 0. ;
				YVelSpecies[i][j]	= 0. ;
				ZVelSpecies[i][j]	= 0. ;
				TempSpecies[i][j]	= 0. ;
				TransTempSpecies[i][j]	= 0. ;
				TransXTempSpecies[i][j]	= 0. ;
				TransYTempSpecies[i][j]	= 0. ;
				TransZTempSpecies[i][j]	= 0. ;
				RotTempSpecies[i][j]	= 0. ;
				VibTempSpecies[i][j]	= 0. ;	
				AveParticleNumSpecies [i][j]	= 0. ;
			}	
		}
}

//==============================================================================================================

DSMC_DOMAIN::DSMC_DOMAIN(){
	Dimension		= 0 ;
	NodeNum			= 0 ;
	CellNum			= 0 ;
	LocalCellNum		= 0 ;
	TotalCellNum		= 0 ;
	InletFaceNum		= 0 ;
	WallFaceNum		= 0 ;
	WallTypeNum		= 0 ;
	MaxParticleNum		= 0 ;
	TransferParticleNum	= 0 ;
	InletSpecifiedNumDenNum	= 0 ;

	SamplingFrequency 	= 0 ;
	OutputFrequency		= 0 ;
	SamplingTime		= 0 ;
	TotalTimestep		= 0 ;

	TimestepRatio		= 0. ;
	WeightingRatio		= 0. ;
	XVel			= 0. ;
	YVel			= 0. ;
	ZVel			= 0. ;
	NumDen			= 0. ;
	Temp			= 0. ;

	SpeciesNum		= 0 ;
	SpeciesGroupNum		= 0 ;

	VariableTimestepScheme	= 0 ;
	SubcellModel		= 0 ;
	DynamicDomainDecompisition= 0 ;
	VibrationalModel	= 0 ;
	SimulationStage		= 0 ;
	
	Convergence		= 0 ;
	ParticleNumCheck	= 100. ;
	SpeedCheck		= 100. ;
	TempCheck		= 100. ;


	TimestepNo		= 0 ; 
	ParticleNum		= 0 ;
	SamplingNum		= 0 ;
	EnterParticleNum	= 0 ;
	
	MinVolume		= 1.E+6 ; 
	MinLength		= 1.E+6 ;
	MinCenter		= 1.E+6 ;
	Timestep		= 0. ;
	ParticleWeighting	= 0. ;
	AxisAdjustFactor	= 1. ;
	
	TrackingNum		= 0 ;
	ErrorTrackingNum	= 0 ;
	ErrorTrackingNumNew	= 0 ;
	
	ChemicalNum = 0 ;
}


void DSMC_DOMAIN::Dump(){
	cout << setw(28) << "Dimension"				<< setw(14) << Dimension << '\n' ;
	cout << setw(28) << "NodeNum"				<< setw(14) << NodeNum << '\n' ;
	cout << setw(28) << "CellNum" 				<< setw(14) << CellNum << '\n' ;
	cout << setw(28) << "LocalCellNum" 			<< setw(14) << LocalCellNum << '\n' ;
	cout << setw(28) << "TotalCellNum" 			<< setw(14) << TotalCellNum << '\n' ;
	cout << setw(28) << "InletFaceNum" 			<< setw(14) << InletFaceNum << '\n' ;
	cout << setw(28) << "WallFaceNum" 			<< setw(14) << WallFaceNum << '\n' ;
	cout << setw(28) << "WallTypeNum" 			<< setw(14) << WallTypeNum << '\n' ;
	cout << setw(28) << "MaxParticleNum" 			<< setw(14) << MaxParticleNum << '\n' ;
	cout << setw(28) << "TransferParticleNum"		<< setw(14) << TransferParticleNum << '\n' ;
	cout << setw(28) << "InletSpecifiedNumDenNum"		<< setw(14) << InletSpecifiedNumDenNum << '\n' ;
	cout << setw(28) << "SamplingFrequency" 		<< setw(14) << SamplingFrequency << '\n' ;
	cout << setw(28) << "OutputFrequency" 			<< setw(14) << OutputFrequency << '\n' ;
	cout << setw(28) << "SamplingTime" 			<< setw(14) << SamplingTime << '\n' ;
	cout << setw(28) << "TotalTimestep" 			<< setw(14) << TotalTimestep << '\n' ;
	cout << setw(28) << "TimestepRatio" 			<< setw(14) << TimestepRatio << '\n' ;
	cout << setw(28) << "WeightingRatio" 			<< setw(14) << WeightingRatio << '\n' ;
	cout << setw(28) << "XVel" 				<< setw(14) << XVel << '\n' ;
	cout << setw(28) << "YVel" 				<< setw(14) << YVel << '\n' ;
	cout << setw(28) << "ZVel" 				<< setw(14) << ZVel << '\n' ;
	cout << setw(28) << "NumDen" 				<< setw(14) << NumDen << '\n' ;
	cout << setw(28) << "Temp" 				<< setw(14) << Temp << '\n' ;
	cout << setw(28) << "SpeciesNum" 			<< setw(14) << SpeciesNum << '\n' ;
	cout << setw(28) << "SpeciesGroupNum" 			<< setw(14) << SpeciesGroupNum << '\n' ;
	cout << setw(28) << "VariableTimestepScheme" 		<< setw(14) << VariableTimestepScheme << '\n' ;
	cout << setw(28) << "SubcellModel" 			<< setw(14) << SubcellModel << '\n' ;
	cout << setw(28) << "DynamicDomainDecompisition"	<< setw(14) << DynamicDomainDecompisition << '\n' ;
	cout << setw(28) << "VibrationalModel" 			<< setw(14) << VibrationalModel << '\n' ;
	cout << setw(28) << "SimulationStage" 			<< setw(14) << SimulationStage << '\n' ;
	cout << setw(28) << "Convergence" 			<< setw(14) << Convergence << '\n' ;
	cout << setw(28) << "ParticleNumCheck" 			<< setw(14) << ParticleNumCheck << '\n' ;
	cout << setw(28) << "SpeedCheck" 			<< setw(14) << SpeedCheck << '\n' ;
	cout << setw(28) << "TempCheck" 			<< setw(14) << TempCheck << '\n' ;
	cout << setw(28) << "TimestepNo" 			<< setw(14) << TimestepNo << '\n' ;
	cout << setw(28) << "ParticleNum" 			<< setw(14) << ParticleNum << '\n' ;
	cout << setw(28) << "SamplingNum" 			<< setw(14) << SamplingNum << '\n' ;
	cout << setw(28) << "EnterParticleNum" 			<< setw(14) << EnterParticleNum << '\n' ;
	cout << setw(28) << "MinVolume" 			<< setw(14) << MinVolume << '\n' ;
	cout << setw(28) << "MinLength" 			<< setw(14) << MinLength << '\n' ;
	cout << setw(28) << "MinCenter" 			<< setw(14) << MinCenter << '\n' ;
	cout << setw(28) << "Timestep" 				<< setw(14) << Timestep << '\n' ;
	cout << setw(28) << "ParticleWeighting" 		<< setw(14) << ParticleWeighting << '\n' ;
	cout << setw(28) << "AxisAdjustFactor" 			<< setw(14) << AxisAdjustFactor << '\n' ;
	cout << setw(28) << "TrackingNum"	 		<< setw(14) << TrackingNum << '\n' ;
	cout << setw(28) << "ErrorTrackingNum"	 		<< setw(14) << ErrorTrackingNum << '\n' ;
	cout << setw(28) << "ErrorTrackingNumNew" 		<< setw(14) << ErrorTrackingNumNew << '\n' ;
	cout << setw(28) << "ChemicalNum" 		<< setw(14) << ChemicalNum << '\n' ;
}

//==============================================================================================================

void DSMC_DOMAIN::DumpFile( ofstream &OutputInfo ){
	OutputInfo << "========================\n" ;
	OutputInfo << setw(28) << "Dimension"			<< setw(14) << Dimension << '\n' ;
	OutputInfo << setw(28) << "NodeNum"			<< setw(14) << NodeNum << '\n' ;
	OutputInfo << setw(28) << "CellNum" 			<< setw(14) << CellNum << '\n' ;
	OutputInfo << setw(28) << "LocalCellNum" 		<< setw(14) << LocalCellNum << '\n' ;
	OutputInfo << setw(28) << "TotalCellNum" 		<< setw(14) << TotalCellNum << '\n' ;
	OutputInfo << setw(28) << "InletFaceNum" 		<< setw(14) << InletFaceNum << '\n' ;
	OutputInfo << setw(28) << "WallFaceNum" 		<< setw(14) << WallFaceNum << '\n' ;
	OutputInfo << setw(28) << "WallTypeNum" 		<< setw(14) << WallTypeNum << '\n' ;
	OutputInfo << setw(28) << "MaxParticleNum" 		<< setw(14) << MaxParticleNum << '\n' ;
	OutputInfo << setw(28) << "TransferParticleNum"		<< setw(14) << TransferParticleNum << '\n' ;
	OutputInfo << setw(28) << "InletSpecifiedNumDenNum"	<< setw(14) << InletSpecifiedNumDenNum << '\n' ;
	OutputInfo << setw(28) << "SamplingFrequency" 		<< setw(14) << SamplingFrequency << '\n' ;
	OutputInfo << setw(28) << "OutputFrequency" 		<< setw(14) << OutputFrequency << '\n' ;
	OutputInfo << setw(28) << "SamplingTime" 		<< setw(14) << SamplingTime << '\n' ;
	OutputInfo << setw(28) << "TotalTimestep" 		<< setw(14) << TotalTimestep << '\n' ;
	OutputInfo << setw(28) << "TimestepRatio"		<< setw(14) << TimestepRatio << '\n' ;
	OutputInfo << setw(28) << "WeightingRatio" 		<< setw(14) << WeightingRatio << '\n' ;
	OutputInfo << setw(28) << "XVel" 			<< setw(14) << XVel << '\n' ;
	OutputInfo << setw(28) << "YVel" 			<< setw(14) << YVel << '\n' ;
	OutputInfo << setw(28) << "ZVel" 			<< setw(14) << ZVel << '\n' ;
	OutputInfo << setw(28) << "NumDen" 			<< setw(14) << NumDen << '\n' ;
	OutputInfo << setw(28) << "Temp" 			<< setw(14) << Temp << '\n' ;
	OutputInfo << setw(28) << "SpeciesNum" 			<< setw(14) << SpeciesNum << '\n' ;
	OutputInfo << setw(28) << "SpeciesGroupNum" 		<< setw(14) << SpeciesGroupNum << '\n' ;
	OutputInfo << setw(28) << "VariableTimestepScheme" 	<< setw(14) << VariableTimestepScheme << '\n' ;
	OutputInfo << setw(28) << "SubcellModel" 		<< setw(14) << SubcellModel << '\n' ;
	OutputInfo << setw(28) << "DynamicDomainDecompisition"	<< setw(14) << DynamicDomainDecompisition << '\n' ;
	OutputInfo << setw(28) << "VibrationalModel" 		<< setw(14) << VibrationalModel << '\n' ;
	OutputInfo << setw(28) << "SimulationStage" 		<< setw(14) << SimulationStage << '\n' ;
	OutputInfo << setw(28) << "Convergence" 			<< setw(14) << Convergence << '\n' ;
	OutputInfo << setw(28) << "ParticleNumCheck" 			<< setw(14) << ParticleNumCheck << '\n' ;
	OutputInfo << setw(28) << "SpeedCheck" 			<< setw(14) << SpeedCheck << '\n' ;
	OutputInfo << setw(28) << "TempCheck" 			<< setw(14) << TempCheck << '\n' ;
	OutputInfo << setw(28) << "TimestepNo" 			<< setw(14) << TimestepNo << '\n' ;
	OutputInfo << setw(28) << "ParticleNum" 		<< setw(14) << ParticleNum << '\n' ;
	OutputInfo << setw(28) << "SamplingNum" 		<< setw(14) << SamplingNum << '\n' ;
	OutputInfo << setw(28) << "MinVolume" 			<< setw(14) << MinVolume << '\n' ;
	OutputInfo << setw(28) << "MinLength" 			<< setw(14) << MinLength << '\n' ;
	OutputInfo << setw(28) << "MinCenter" 			<< setw(14) << MinCenter << '\n' ;
	OutputInfo << setw(28) << "Timestep" 			<< setw(14) << Timestep << '\n' ;
	OutputInfo << setw(28) << "ParticleWeighting" 		<< setw(14) << ParticleWeighting << '\n' ;
	OutputInfo << setw(28) << "AxisAdjustFactor" 		<< setw(14) << AxisAdjustFactor << '\n' ;
	//OutputInfo << setw(28) << "TrackingNum"	 	<< setw(14) << TrackingNum << '\n' ;
	//OutputInfo << setw(28) << "ErrorTrackingNum"	 	<< setw(14) << ErrorTrackingNum << '\n' ;
	//OutputInfo << setw(28) << "ErrorTrackingNumNew" 	<< setw(14) << ErrorTrackingNumNew << '\n' ;
	OutputInfo << "========================\n" ;
}

//==============================================================================================================

DSMC_CONVERGENCE::DSMC_CONVERGENCE(){
	DataNum		= 200 ;
	CurrentNo	= 0 ;
	ConditionParticleNum = 100. ;
	ConditionSpeed	= 100. ;
	ConditionTemp	= 100. ;
	ParticleNum	= NULL ;
	Speed		= NULL ;
	Temp		= NULL ;
	OldParticleNum	= 0. ;
	OldSpeed	= 0. ;
	OldTemp		= 0. ;
	DiffParticleNum	= 100. ;
	DiffSpeed	= 100. ;
	DiffTemp	= 100. ;
}


DSMC_CONVERGENCE::~DSMC_CONVERGENCE(){
	delete [] ParticleNum ;
	delete [] Speed ;
	delete [] Temp ;
}


void DSMC_CONVERGENCE::AllocateMemory(){
	ParticleNum	= new double[DataNum] ;
	Speed		= new double[DataNum] ;
	Temp		= new double[DataNum] ;	
}


void DSMC_CONVERGENCE::InitValue( double ParticleNumCheck , double SpeedCheck , double TempCheck ){
	ConditionParticleNum	= ParticleNumCheck ;
	ConditionSpeed		= SpeedCheck ;
	ConditionTemp		= TempCheck ;
	
	for ( int i=0 ; i<DataNum ; i++ ){
		ParticleNum[i]	= 0. ;
		Speed[i]	= 0. ;
		Temp[i]		= 0. ;
	}
}


bool DSMC_CONVERGENCE::Convergence(){
	bool		a = true ;
	
	DiffParticleNum	= fabs( ParticleNum[CurrentNo] - OldParticleNum )/OldParticleNum ;
	//DiffSpeed	= fabs( Speed[CurrentNo] - OldSpeed ) ;
	//DiffTemp	= fabs( Temp[CurrentNo] - OldTemp ) ;
	DiffSpeed	= fabs( Speed[CurrentNo] - OldSpeed )/OldSpeed ;
	DiffTemp	= fabs( Temp[CurrentNo] - OldTemp )/OldTemp ;
	
	if ( DiffParticleNum > ConditionParticleNum ) a = false ;
	if ( DiffSpeed > ConditionSpeed ) a = false ;
	if ( DiffTemp > ConditionTemp ) a = false ;
	
	return	a ;
}


void DSMC_CONVERGENCE::Sum( double AveParticleNum , double AveSpeed , double AveTemp , int Num ){
	for ( int i=0 ; i<Num ; i++ ){
		ParticleNum[i]	+= AveParticleNum ;
		Speed[i]	+= AveSpeed ;
		Temp[i]		+= AveTemp ;
	}
}


void DSMC_CONVERGENCE::InitSum(){
	OldParticleNum		= ParticleNum[CurrentNo] ;
	OldSpeed		= Speed[CurrentNo] ;
	OldTemp			= Temp[CurrentNo] ;
	ParticleNum[CurrentNo]	= 0. ;
	Speed[CurrentNo]	= 0. ;
	Temp[CurrentNo]		= 0. ;
	CurrentNo		= (CurrentNo+1)%DataNum ;
}


void DSMC_CONVERGENCE::Average(){
	ParticleNum[CurrentNo]	/= DataNum ;
	Speed[CurrentNo]	/= DataNum ;
	Temp[CurrentNo]		/= DataNum ;
}

//==============================================================================================================

DSMC_PROCESSOR::DSMC_PROCESSOR(){
	NeighborNum	= 0 ;
	
	for ( int i=0 ; i<64 ; i++ )
		Neighbor[i]	= -999 ;
		
	MaxNeighborNum	= 64 ;
	LocalCellNo	= NULL ;
	CellProcessorNo	= NULL ;
}


void DSMC_PROCESSOR::AllocateMemory( int TotalCellNum ){
	LocalCellNo	= new int[TotalCellNum]	;
	CellProcessorNo	= new int[TotalCellNum]	;
}


void DSMC_PROCESSOR::InitValue( int TotalCellNum ){
	for ( int i=0 ; i<TotalCellNum ; i++ ){
		LocalCellNo[i]		= -1 ;
		CellProcessorNo[i]	= -1 ;
	}
}


void DSMC_PROCESSOR::DeleteMemory(){
	delete [] LocalCellNo ;
	delete [] CellProcessorNo ;
}


void DSMC_PROCESSOR::Dump( int ID , int CellNum ){
	cout << "Processor ID: " << ID << ", CellNum: " << CellNum << ", NeighborNum: " << NeighborNum ;
	for ( int i=0 ; i<NeighborNum ; i++ )
		cout << setw(8) << Neighbor[i] ;
	cout << '\n' ;
}

//==============================================================================================================

DSMC_MPI_PARTICLE::DSMC_MPI_PARTICLE(){
	ProcessorNo	= 0 ;
	GlobalCellNo	= 0 ;
	BeforeCellNo	= 0 ;
	SpeciesNo	= 0 ;
	LastCollide	= 0 ;
	VibLevel	= 0 ;
	
	XCoord		= 0. ;
	YCoord		= 0. ;
	ZCoord		= 0. ;
	XVel		= 0. ;
	YVel		= 0. ;
	ZVel		= 0. ;
	Rotation	= 0. ;
	Vibration	= 0. ;
	EffTemp		= 0. ;
	Timestep	= 0. ;
	BeforeTimestep	= 0. ;
}

//==============================================================================================================

void DSMC_MPI_DATATYPE::InitMPIDataType( DSMC_DOMAIN *pDomain ){
	DomainVariableNum	= 47 ;
	
	
	for ( int i=0 ; i<DomainVariableNum ; i++ )
		DomainVariableLen[i]	= 1 ;
	
	
	for ( int i=0 ; i<DomainVariableNum ; i++ )
		DomainOldType[i]	= MPI_INT ;
	
	DomainOldType[15]	= MPI_DOUBLE ;
	DomainOldType[16]	= MPI_DOUBLE ;
	DomainOldType[17]	= MPI_DOUBLE ;
	DomainOldType[18]	= MPI_DOUBLE ;
	DomainOldType[19]	= MPI_DOUBLE ;
	DomainOldType[20]	= MPI_DOUBLE ;
	DomainOldType[21]	= MPI_DOUBLE ;
	
	DomainOldType[30]	= MPI_DOUBLE ;
	DomainOldType[31]	= MPI_DOUBLE ;
	DomainOldType[32]	= MPI_DOUBLE ;
	
	DomainOldType[37]	= MPI_DOUBLE ;
	DomainOldType[38]	= MPI_DOUBLE ;
	DomainOldType[39]	= MPI_DOUBLE ;
	DomainOldType[40]	= MPI_DOUBLE ;
	DomainOldType[41]	= MPI_DOUBLE ;
	DomainOldType[42]	= MPI_DOUBLE ;
	
	
	MPI_Address( (void*)&pDomain->Dimension		, &DomainDisp[0] ) ;
	MPI_Address( (void*)&pDomain->NodeNum		, &DomainDisp[1] ) ;
	MPI_Address( (void*)&pDomain->CellNum		, &DomainDisp[2] ) ;
	MPI_Address( (void*)&pDomain->LocalCellNum	, &DomainDisp[3] ) ;
	MPI_Address( (void*)&pDomain->TotalCellNum	, &DomainDisp[4] ) ;
	
	MPI_Address( (void*)&pDomain->InletFaceNum	, &DomainDisp[5] ) ;
	MPI_Address( (void*)&pDomain->WallFaceNum	, &DomainDisp[6] ) ;
	MPI_Address( (void*)&pDomain->WallTypeNum	, &DomainDisp[7] ) ;
	MPI_Address( (void*)&pDomain->MaxParticleNum	, &DomainDisp[8] ) ;
	MPI_Address( (void*)&pDomain->TransferParticleNum, &DomainDisp[9] ) ;
	MPI_Address( (void*)&pDomain->InletSpecifiedNumDenNum, &DomainDisp[10] ) ;
	
	MPI_Address( (void*)&pDomain->SamplingFrequency	, &DomainDisp[11] ) ;
	MPI_Address( (void*)&pDomain->OutputFrequency	, &DomainDisp[12] ) ;
	MPI_Address( (void*)&pDomain->SamplingTime	, &DomainDisp[13] ) ;
	MPI_Address( (void*)&pDomain->TotalTimestep	, &DomainDisp[14] ) ;
	
	MPI_Address( (void*)&pDomain->TimestepRatio	, &DomainDisp[15] ) ;
	MPI_Address( (void*)&pDomain->WeightingRatio	, &DomainDisp[16] ) ;
	MPI_Address( (void*)&pDomain->XVel		, &DomainDisp[17] ) ;
	MPI_Address( (void*)&pDomain->YVel		, &DomainDisp[18] ) ;
	MPI_Address( (void*)&pDomain->ZVel		, &DomainDisp[19] ) ;
	MPI_Address( (void*)&pDomain->NumDen		, &DomainDisp[20] ) ;
	MPI_Address( (void*)&pDomain->Temp		, &DomainDisp[21] ) ;
	
	MPI_Address( (void*)&pDomain->SpeciesNum	, &DomainDisp[22] ) ;
	MPI_Address( (void*)&pDomain->SpeciesGroupNum	, &DomainDisp[23] ) ;
	
	MPI_Address( (void*)&pDomain->VariableTimestepScheme	, &DomainDisp[24] ) ;
	MPI_Address( (void*)&pDomain->SubcellModel		, &DomainDisp[25] ) ;
	MPI_Address( (void*)&pDomain->DynamicDomainDecompisition, &DomainDisp[26] ) ;
	MPI_Address( (void*)&pDomain->VibrationalModel		, &DomainDisp[27] ) ;
	MPI_Address( (void*)&pDomain->SimulationStage		, &DomainDisp[28] ) ;
	
	MPI_Address( (void*)&pDomain->Convergence		, &DomainDisp[29] ) ;
	MPI_Address( (void*)&pDomain->ParticleNumCheck		, &DomainDisp[30] ) ;
	MPI_Address( (void*)&pDomain->SpeedCheck		, &DomainDisp[31] ) ;
	MPI_Address( (void*)&pDomain->TempCheck			, &DomainDisp[32] ) ;

	MPI_Address( (void*)&pDomain->TimestepNo	, &DomainDisp[33] ) ;
	MPI_Address( (void*)&pDomain->ParticleNum	, &DomainDisp[34] ) ;
	MPI_Address( (void*)&pDomain->SamplingNum	, &DomainDisp[35] ) ;
	MPI_Address( (void*)&pDomain->EnterParticleNum	, &DomainDisp[36] ) ;
	
	MPI_Address( (void*)&pDomain->MinVolume		, &DomainDisp[37] ) ;
	MPI_Address( (void*)&pDomain->MinLength		, &DomainDisp[38] ) ;
	MPI_Address( (void*)&pDomain->MinCenter		, &DomainDisp[39] ) ;
	MPI_Address( (void*)&pDomain->Timestep		, &DomainDisp[40] ) ;
	MPI_Address( (void*)&pDomain->ParticleWeighting	, &DomainDisp[41] ) ;
	MPI_Address( (void*)&pDomain->AxisAdjustFactor	, &DomainDisp[42] ) ;
	
	MPI_Address( (void*)&pDomain->TrackingNum		, &DomainDisp[43] ) ;
	MPI_Address( (void*)&pDomain->ErrorTrackingNum		, &DomainDisp[44] ) ;
	MPI_Address( (void*)&pDomain->ErrorTrackingNumNew	, &DomainDisp[45] ) ;
	
	MPI_Address( (void*)&pDomain->ChemicalNum	, &DomainDisp[46] ) ;
	
	
	for ( int i=(DomainVariableNum-1) ; i>=0 ; i-- )
		DomainDisp[i] -= DomainDisp[0] ;
	
	
	MPI_Type_struct( DomainVariableNum , DomainVariableLen, DomainDisp , DomainOldType , &MPI_DOMAIN ) ;
	MPI_Type_commit( &MPI_DOMAIN ) ;
}


void DSMC_MPI_DATATYPE::InitMPIDataType( DSMC_MPI_PARTICLE *pMPIParticle ){
	ParticleVariableNum	= 18 ;
	
	for ( int i=0 ; i<ParticleVariableNum ; i++ )
		ParticleVariableLen[i]	= 1 ;
		
	
	for ( int i=0 ; i<ParticleVariableNum ; i++ )
		ParticleOldType[i]	= MPI_INT ;
		
	
	ParticleOldType[6]	= MPI_DOUBLE ;
	ParticleOldType[7]	= MPI_DOUBLE ;
	ParticleOldType[8]	= MPI_DOUBLE ;
	ParticleOldType[9]	= MPI_DOUBLE ;
	ParticleOldType[10]	= MPI_DOUBLE ;
	ParticleOldType[11]	= MPI_DOUBLE ;
	ParticleOldType[12]	= MPI_DOUBLE ;
	ParticleOldType[13]	= MPI_DOUBLE ;
	ParticleOldType[14]	= MPI_DOUBLE ;
	ParticleOldType[15]	= MPI_DOUBLE ;
	ParticleOldType[16]	= MPI_DOUBLE ;
	
	
	MPI_Address( (void*)&pMPIParticle->ProcessorNo	, &ParticleDisp[0] ) ;
	MPI_Address( (void*)&pMPIParticle->GlobalCellNo	, &ParticleDisp[1] ) ;
	MPI_Address( (void*)&pMPIParticle->BeforeCellNo	, &ParticleDisp[2] ) ;
	MPI_Address( (void*)&pMPIParticle->SpeciesNo	, &ParticleDisp[3] ) ;
	MPI_Address( (void*)&pMPIParticle->LastCollide	, &ParticleDisp[4] ) ;
	MPI_Address( (void*)&pMPIParticle->VibLevel	, &ParticleDisp[5] ) ;
	MPI_Address( (void*)&pMPIParticle->XCoord	, &ParticleDisp[6] ) ;
	MPI_Address( (void*)&pMPIParticle->YCoord	, &ParticleDisp[7] ) ;
	MPI_Address( (void*)&pMPIParticle->ZCoord	, &ParticleDisp[8] ) ;
	MPI_Address( (void*)&pMPIParticle->XVel		, &ParticleDisp[9] ) ;
	MPI_Address( (void*)&pMPIParticle->YVel		, &ParticleDisp[10] ) ;
	MPI_Address( (void*)&pMPIParticle->ZVel		, &ParticleDisp[11] ) ;
	MPI_Address( (void*)&pMPIParticle->Rotation	, &ParticleDisp[12] ) ;
	MPI_Address( (void*)&pMPIParticle->Vibration	, &ParticleDisp[13] ) ;
	MPI_Address( (void*)&pMPIParticle->EffTemp	, &ParticleDisp[14] ) ;
	MPI_Address( (void*)&pMPIParticle->Timestep	, &ParticleDisp[15] ) ;
	MPI_Address( (void*)&pMPIParticle->BeforeTimestep, &ParticleDisp[16] ) ;
	MPI_Address( (void*)&pMPIParticle->ReactionIndex, &ParticleDisp[17] ) ;
	
	
	for ( int i=(ParticleVariableNum-1) ; i>=0 ; i-- )
		ParticleDisp[i] -= ParticleDisp[0] ;
	
	
	MPI_Type_struct( ParticleVariableNum , ParticleVariableLen, ParticleDisp , ParticleOldType , &MPI_PARTICLE ) ;
	MPI_Type_commit( &MPI_PARTICLE ) ;
}

//==============================================================================================================

CELLMAPPING::CELLMAPPING(){
	for ( int i=0 ; i<6 ; i++ ){
		for ( int j=0 ; j<6 ; j++ ){
			SurfaceNode[i][j]	= 0 ;
			
			for ( int k=0 ; k<4 ; k++ )
				Node[i][j][k]	= 0 ;
		}	
		
		for ( int j=0 ; j<5 ; j++ ){
			for ( int k=0 ; k<4 ; k++ )
				TetraCell[i][j][k]	= 0 ;
		}	
	}
	
	
	// For Tetra Cell Type
	NodeNum[0]		= 4 ;
	SurfaceNum[0]		= 4 ;
	TetraCellNum[0]		= 1 ; 

	SurfaceNode[0][0]	= 3 ;
	SurfaceNode[0][1]	= 3 ;
	SurfaceNode[0][2]	= 3 ;
	SurfaceNode[0][3]	= 3 ;

	Node[0][0][0]		= 1 ;
	Node[0][0][1]		= 2 ;
	Node[0][0][2]		= 4 ;

	Node[0][1][0]		= 2 ;
	Node[0][1][1]		= 3 ;
	Node[0][1][2]		= 4 ;

	Node[0][2][0]		= 3 ;
	Node[0][2][1]		= 1 ;
	Node[0][2][2]		= 4 ;

	Node[0][3][0]		= 3 ;
	Node[0][3][1]		= 2 ;
	Node[0][3][2]		= 1 ;
	
	TetraCell[0][0][0]	= 1 ;
	TetraCell[0][0][1]	= 2 ;
	TetraCell[0][0][2]	= 3 ;
	TetraCell[0][0][3]	= 4 ;
	// End (tetra cell)
	


	// For Hex. Cell Type
	NodeNum[1]		= 8 ;
	SurfaceNum[1]		= 6 ;
	TetraCellNum[1]		= 5 ; 
                                  
	SurfaceNode[1][0]	= 4 ;
	SurfaceNode[1][1]	= 4 ;
	SurfaceNode[1][2]	= 4 ;
	SurfaceNode[1][3]	= 4 ;
	SurfaceNode[1][4]	= 4 ;
	SurfaceNode[1][5]	= 4 ;
                                  
	Node[1][0][0]		= 1 ;
	Node[1][0][1]		= 2 ;
	Node[1][0][2]		= 6 ;
	Node[1][0][3]		= 5 ;
                                  
	Node[1][1][0]		= 2 ;
	Node[1][1][1]		= 3 ;
	Node[1][1][2]		= 7 ;
	Node[1][1][3]		= 6 ;
                                  
	Node[1][2][0]		= 3 ;
	Node[1][2][1]		= 4 ;
	Node[1][2][2]		= 8 ;
	Node[1][2][3]		= 7 ;
                                  
	Node[1][3][0]		= 4 ;
	Node[1][3][1]		= 1 ;
	Node[1][3][2]		= 5 ;
	Node[1][3][3]		= 8 ;
                                  
	Node[1][4][0]		= 3 ;
	Node[1][4][1]		= 2 ;
	Node[1][4][2]		= 1 ;
	Node[1][4][3]		= 4 ;
                                  
	Node[1][5][0]		= 6 ;
	Node[1][5][1]		= 7 ;
	Node[1][5][2]		= 8 ;
	Node[1][5][3]		= 5 ;
	                          
	TetraCell[1][0][0]	= 1 ;
	TetraCell[1][0][1]	= 2 ;
	TetraCell[1][0][2]	= 3 ;
	TetraCell[1][0][3]	= 6 ;
	                          
	TetraCell[1][1][0]	= 1 ;
	TetraCell[1][1][1]	= 3 ;
	TetraCell[1][1][2]	= 4 ;
	TetraCell[1][1][3]	= 8 ;
	                          
	TetraCell[1][2][0]	= 1 ;
	TetraCell[1][2][1]	= 3 ;
	TetraCell[1][2][2]	= 6 ;
	TetraCell[1][2][3]	= 8 ;
	                          
	TetraCell[1][3][0]	= 3 ;
	TetraCell[1][3][1]	= 6 ;
	TetraCell[1][3][2]	= 7 ;
	TetraCell[1][3][3]	= 8 ;
	                          
	TetraCell[1][4][0]	= 1 ;
	TetraCell[1][4][1]	= 5 ;
	TetraCell[1][4][2]	= 6 ;
	TetraCell[1][4][3]	= 8 ;
	// End (hex. cell)        
                                  
                                  
                                  
	// For Pyramid Cell Type  
	NodeNum[2]		= 5 ;
	SurfaceNum[2]		= 5 ;
	TetraCellNum[2]		= 2 ;
                                  
	SurfaceNode[2][0]	= 3 ;
	SurfaceNode[2][1]	= 3 ;
	SurfaceNode[2][2]	= 3 ;
	SurfaceNode[2][3]	= 3 ;
	SurfaceNode[2][4]	= 4 ;
                                  
	Node[2][0][0]		= 1 ;
	Node[2][0][1]		= 2 ;
	Node[2][0][2]		= 5 ;
                                  
	Node[2][1][0]		= 2 ;
	Node[2][1][1]		= 3 ;
	Node[2][1][2]		= 5 ;
                                  
	Node[2][2][0]		= 3 ;
	Node[2][2][1]		= 4 ;
	Node[2][2][2]		= 5 ;
                                  
	Node[2][3][0]		= 4 ;
	Node[2][3][1]		= 1 ;
	Node[2][3][2]		= 5 ;
                                  
	Node[2][4][0]		= 3 ;
	Node[2][4][1]		= 2 ;
	Node[2][4][2]		= 1 ;
	Node[2][4][3]		= 4 ;
	
	TetraCell[2][0][0]	= 1 ;
	TetraCell[2][0][1]	= 2 ;
	TetraCell[2][0][2]	= 4 ;
	TetraCell[2][0][3]	= 5 ;
	
	TetraCell[2][1][0]	= 2 ;
	TetraCell[2][1][1]	= 3 ;
	TetraCell[2][1][2]	= 4 ;
	TetraCell[2][1][3]	= 5 ;
	// End (pyramid cell)     
                                  
                                  
                                  
	// For Prism Cell Type    
	NodeNum[3]		= 6 ;
	SurfaceNum[3]		= 5 ;
	TetraCellNum[3]		= 3 ;
                                  
	SurfaceNode[3][0]	= 4 ;
	SurfaceNode[3][1]	= 3 ;
	SurfaceNode[3][2]	= 4 ;
	SurfaceNode[3][3]	= 3 ;
	SurfaceNode[3][4]	= 4 ;
                                  
	Node[3][0][0]		= 1 ;
	Node[3][0][1]		= 2 ;
	Node[3][0][2]		= 6 ;
	Node[3][0][3]		= 5 ;
                                  
	Node[3][1][0]		= 2 ;
	Node[3][1][1]		= 3 ;
	Node[3][1][2]		= 6 ;
                                  
	Node[3][2][0]		= 3 ;
	Node[3][2][1]		= 4 ;
	Node[3][2][2]		= 5 ;
	Node[3][2][3]		= 6 ;
                                  
	Node[3][3][0]		= 4 ;
	Node[3][3][1]		= 1 ;
	Node[3][3][2]		= 5 ;
                                  
	Node[3][4][0]		= 3 ;
	Node[3][4][1]		= 2 ;
	Node[3][4][2]		= 1 ;
	Node[3][4][3]		= 4 ;
	
	TetraCell[3][0][0]	= 1 ;
	TetraCell[3][0][1]	= 4 ;
	TetraCell[3][0][2]	= 5 ;
	TetraCell[3][0][3]	= 2 ;
	
	TetraCell[3][1][0]	= 4 ;
	TetraCell[3][1][1]	= 3 ;
	TetraCell[3][1][2]	= 5 ;
	TetraCell[3][1][3]	= 2 ;
	
	TetraCell[3][2][0]	= 2 ;
	TetraCell[3][2][1]	= 3 ;
	TetraCell[3][2][2]	= 5 ;
	TetraCell[3][2][3]	= 6 ;
	// End (prism cell)



	// For 2-D Triangular Cell Type
	NodeNum[4]		= 3 ;
	SurfaceNum[4]		= 3 ;
	TetraCellNum[4]		= 0 ;
                                  
	SurfaceNode[4][0]	= 2 ;
	SurfaceNode[4][1]	= 2 ;
	SurfaceNode[4][2]	= 2 ;
                                  
	Node[4][0][0]		= 1 ;
	Node[4][0][1]		= 2 ;
                                  
	Node[4][1][0]		= 2 ;
	Node[4][1][1]		= 3 ;
                                  
	Node[4][2][0]		= 3 ;
	Node[4][2][1]		= 1 ;
	// End (2-d triangular cell)



	// For 2-D Quadrilateral Cell Type
	NodeNum[5]		= 4 ;
	SurfaceNum[5]		= 4 ;
	TetraCellNum[5]		= 0 ;
                                  
	SurfaceNode[5][0]	= 2 ;
	SurfaceNode[5][1]	= 2 ;
	SurfaceNode[5][2]	= 2 ;
	SurfaceNode[5][3]	= 2 ;
                                  
	Node[5][0][0]		= 1 ;
	Node[5][0][1]		= 2 ;
                                  
	Node[5][1][0]		= 2 ;
	Node[5][1][1]		= 3 ;
                                  
	Node[5][2][0]		= 3 ;
	Node[5][2][1]		= 4 ;
                                  
	Node[5][3][0]		= 4 ;
	Node[5][3][1]		= 1 ;
	// End (2-d quadrilateral cell)



	for ( int i=0 ; i<6 ; i++ ){
		for ( int j=0 ; j<6 ; j++ ){
			for ( int k=0 ; k<4 ; k++ )
				Node[i][j][k]-- ;
		}
	}
	
	for ( int i=0 ; i<6 ; i++ ){
		for ( int j=0 ; j<5 ; j++ ){
			for ( int k=0 ; k<4 ; k++ )
				TetraCell[i][j][k]-- ;
		}
	}
}

//==============================================================================================================

DSMC_POST_SURFACE::DSMC_POST_SURFACE(){
	Id		= 0 ;
	CellNo		= 0 ;
	WallNo		= 0 ;
	FaceNo		= 0 ;
	
	for ( int i=0 ; i<17 ; i++ )
		SurfacePro[i]	= 0. ;
		
	for ( int i=0 ; i<3 ; i++ )
		Center[i]	= 0. ;
		
	Area		= 0. ;
	SimTime		= 0. ;
	DepositionRate	= 0. ;
	Cp		= 0. ;
	Cf		= 0. ;
	Ch		= 0. ;
	Pressure	= 0. ;
	HeatFlux	= 0. ;
}


void DSMC_POST_SURFACE::Dump( int No ){
	cout << "==============================================\n" ;
	cout << No << '\n' ;
	cout << Id << '\n' ;
	cout << CellNo << '\n' ;
	cout << WallNo << '\n' ;
	
	for ( int i=0 ; i<15 ; i++ )
		cout << SurfacePro[i] << '\n' ;
		
	for ( int i=0 ; i<4 ; i++ )
		cout << Center[i] << '\n' ;
		
	cout << Area << '\n' ;
	cout << SimTime << '\n' ;
	cout << DepositionRate << '\n' ;
	cout << Cp << '\n' ;
	cout << Cf << '\n' ;
	cout << Ch << '\n' ;
	cout << Pressure << '\n' ;
}

//==============================================================================================================

DSMC_POST_SURFACE_CONDITION::DSMC_POST_SURFACE_CONDITION( string Filename ){
	ifstream		Input ;
	string			get_line , word , number ;
	int			wordstart , wordend , wallno ;

	wallno	= 0 ;

	Input.open( Filename.c_str() , ios::in ) ;

	if ( !Input ){
		cout << "Failed to open file " << Filename << endl ;
		exit(1) ;
	}
	
	while ( getline(Input, get_line) ){
		if ( get_line[0] != '#' ){
			wordend = get_line.find(' ') ;
			word = get_line.substr(0,wordend) ;


			if ( word == "Density_Free_Stream"){
				wordstart	= get_line.find('=') ;
				number		= get_line.substr(wordstart+2) ;

				Density	= atof( number.c_str() ) ;
			
			}else if ( word == "Velocity_Free_Stream"){
				wordstart	= get_line.find('=') ;
				number		= get_line.substr(wordstart+2) ;

				Velocity = atof( number.c_str() ) ;
				
			}else if ( word == "Mass_Species"){
				wordstart	= get_line.find('=') ;
				number		= get_line.substr(wordstart+2) ;

				Mass = atof( number.c_str() ) ;

			}else if ( word == "Output_Wall_Number"){
				wordstart	= get_line.find('=') ;
				number		= get_line.substr(wordstart+2) ;

				WallNum = atoi( number.c_str() ) ;
				
				WallNo	= new int[WallNum] ;

			}else if ( word == "Output_Wall_No"){
				wordstart	= get_line.find('=') ;
				number		= get_line.substr(wordstart+2) ;

				WallNo[wallno]	= atoi( number.c_str() ) ;
				
				wallno++ ;

			}else if ( word == "Sort_Name"){
				wordstart	= get_line.find('=') ;
				number		= get_line.substr(wordstart+2) ;

				SortName = number ;
			}
		}
	}
	
	// To close the file.
	Input.clear() ;
	Input.close() ;

	Drag		= 0. ;
	HeatFlux	= 0. ;
	CoeffAxial = 0. ;
	CoeffNomal = 0. ;
}


DSMC_POST_SURFACE_CONDITION::~DSMC_POST_SURFACE_CONDITION(){
	delete [] WallNo ;	
}


void DSMC_POST_SURFACE_CONDITION::Dump( int WallFaceNum ){
	cout << "\n==============================================================================\n" ;
	cout << setw(30) << "Density in free-stream:"	<< setw(16) << Density << '\n' ;
	cout << setw(30) << "Velocity in free-stream:"	<< setw(16) << Velocity << '\n' ;
	cout << setw(30) << "Mass of Species:"		<< setw(16) << Mass << '\n' ;
	cout << setw(30) << "Drag:"			<< setw(16) << Drag << '\n' ;
	cout << setw(30) << "HeatFlux:"			<< setw(16) << HeatFlux << '\n' ;
	cout << setw(30) << "Post-processed wall number:" << setw(16) << WallNum << '\n' ;
	for ( int i=0 ; i<WallNum ; i++ )
		cout << setw(26) << "WallNo[" << setw(2) << i << "]:" << setw(16) << WallNo[i] << '\n' ;
	cout << setw(30) << "Sort Name:" 		<< setw(16) << SortName << '\n' ;
	cout << setw(30) << "Wall Face number:"		<< setw(16) << WallFaceNum << '\n' ;
	cout << "==============================================================================\n\n" ;
}

//==============================================================================================================

DSMC_POST_NODECELL::~DSMC_POST_NODECELL(){
	CellRelation.clear() ;
}

//==============================================================================================================

DSMC_TIME::DSMC_TIME(){
	Total		= 0. ;
	Init		= 0. ;
	Move		= 0. ;
	Index		= 0. ;
	Collision	= 0. ;
	Sample		= 0. ;
	CalResult	= 0. ;

	time		= 0. ;
	tstart		= 0. ;
}


void DSMC_TIME::Start(){
	gettimeofday(&start,0) ;
	tstart = (long double)(start.tv_sec+start.tv_usec*1e-6) ;
}


void DSMC_TIME::End(){
	gettimeofday(&end,0) ;
	Total = (long double)(end.tv_sec+end.tv_usec*1e-6) - tstart ; 
}


void DSMC_TIME::Time(){
	struct timeval	tv ;
	gettimeofday(&tv,0) ;
	time = (long double)(tv.tv_sec+tv.tv_usec*1e-6) ;
}


void DSMC_TIME::Time( long double *t ){
	struct timeval	tv ;
	gettimeofday(&tv,0) ;
	(*t) += (long double)(tv.tv_sec+tv.tv_usec*1e-6) - time ;
}


void DSMC_TIME::PrintFile(){
	ofstream	Output ;
	
	Output.open( "SimulationTime.dat" , ios::out | ios::trunc ) ;
	
	if ( !Output.fail() ){
		Output << setw(30) << "The simulation time (second):"	<< setw(15) << Total	<< '\n' ;
		Output << setw(30) << "Initialization:" 		<< setw(15) << Init	<< '\n' ;
		Output << setw(30) << "Particle Movement:"		<< setw(15) << Move	<< '\n' ;
		Output << setw(30) << "Index:"				<< setw(15) << Index	<< '\n' ;
		Output << setw(30) << "Collision:"			<< setw(15) << Collision	<< '\n' ; 
		Output << setw(30) << "Sampling:"			<< setw(15) << Sample	<< '\n' ; 
		Output << setw(30) << "Calculate Result:"		<< setw(15) << CalResult	<< '\n' ;
	}else{
		cout << "Fail to open SimTime.dat\n" ;
	}

	
	Output.clear() ;
	Output.close() ;
}

//==============================================================================================================
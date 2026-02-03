#include <iostream>
#include <iomanip>
#include <cstdlib>
#include "dsmc_class_pre.h"

using namespace std ;

//==============================================================================================================

DSMC_PRE_NODE::DSMC_PRE_NODE(){
	Id		= 0 ;
	XCoord		= 0. ;
	YCoord		= 0. ;
	ZCoord		= 0. ;
}

//==============================================================================================================

DSMC_PRE_CELL::DSMC_PRE_CELL(){
	Id		= 0 ;
	Type		= -1 ;
	ProcessorNo	= -1 ;

	for ( int i=0 ; i<6 ; i++ ) Neighbor[i]	= -999 ;
	for ( int i=0 ; i<8 ; i++ ) pNode[i]	= NULL ;

	Weight		= 0 ;
	MeanFreePath	= 0. ;
	Temp		= 0. ;
	Speed		= 0. ;
	Timestep	= 0. ;
	Weighting	= 0. ;
	TimestepRatio	= 0. ;
	WeightingRatio	= 0. ;
	SubcellNum	= 0 ;
}


void DSMC_PRE_CELL::DumpNode( int CellNo , int NodeNum ){
	cout << "CellNo: " << CellNo << ", ID: " << Id << '\n' ;
	for ( int i=0 ; i<NodeNum ; i++ ){
		cout << setw(12) << pNode[i]->Id << setw(16) << pNode[i]->XCoord << setw(16) << pNode[i]->YCoord << setw(16) << pNode[i]->ZCoord << '\n' ;
	}
}

//==============================================================================================================

DSMC_PRE_BCFACE::DSMC_PRE_BCFACE(){
	Id		= 0 ;
	NodeNum		= 0 ;
	CountNum	= 0 ;

	for ( int i=0 ; i<4 ; i++ ) Node[i]	= 0 ;
	for ( int i=0 ; i<2 ; i++ ) NewNode[i]	= 0 ;
	
	XCenter		= 0. ;
	YCenter		= 0. ;
	ZCenter		= 0. ;
	XVel		= 0. ;
	YVel		= 0. ;
	ZVel		= 0. ;
	Temp		= 0. ;
}


DSMC_PRE_BCFACE::~DSMC_PRE_BCFACE(){
	NumDen.clear() ;
}

//==============================================================================================================

DSMC_PRE_BCTYPE::DSMC_PRE_BCTYPE(){
	Type		= 0 ;
	FaceNum		= 0 ;
	TypeName	= "NoType" ;
	XVel		= 0. ;
	YVel		= 0. ;
	ZVel		= 0. ;
	NumDen		= 0. ;
	Temp		= 0. ;
	CosineLawCoef	= -1. ;
}


void DSMC_PRE_BCTYPE::Dump( int BCTypeNo ){
	cout << "Face No.: " << setw(8) << BCTypeNo << '\n' ;
	cout << "===============================\n" ;
	cout << setw(25) << "TypeName" << setw(8) << "Type" << setw(10) << "XVel" << setw(10) << "YVel" << setw(10) << "ZVel" << setw(16) << "NumDen" << setw(8) << "Temp" << setw(17) << "CosineLawCoef\n" ; 
	cout << setw(25) << TypeName << setw(8) << Type << setw(10) << XVel << setw(10) << YVel << setw(10) << ZVel << setw(16) << NumDen << setw(8) << Temp << setw(16) << CosineLawCoef << '\n' ;
	cout << "=====================================================================================================\n" ;
}


void DSMC_PRE_BCTYPE::Dump(){
	cout << setw(20) << TypeName << setw(8) << Type << setw(15) << FaceNum << '\n' ;
}

//==============================================================================================================

DSMC_PRE_NODECELL::~DSMC_PRE_NODECELL(){
	CellRelation.clear() ;
	FaceRelation.clear() ;
}


void DSMC_PRE_NODECELL::Init(){
	CellRelation.clear() ;
	FaceRelation.clear() ;
}

//==============================================================================================================

DSMC_PRE_GRAPHCELL::DSMC_PRE_GRAPHCELL(){
	Num		= 0 ;
	for ( int i=0 ; i<6 ; i++ ){
		Check[i]	= -999 ;
		Neighbor[i]	= -999 ;
	}
}

//==============================================================================================================

DSMC_PRE_PROCESSOR::~DSMC_PRE_PROCESSOR(){
	Cell.clear() ;
	Neighbor.clear() ;
}

//==============================================================================================================

DSMC_PRE_INLET::DSMC_PRE_INLET(){
	XCoord	= 0. ;
	YCoord	= 0. ;
	ZCoord	= 0. ;
	XVel	= 0. ;
	YVel	= 0. ;
	ZVel	= 0. ;
	NumDen	= NULL ;
	Temp	= 0. ;
}


DSMC_PRE_INLET::~DSMC_PRE_INLET(){
	delete [] NumDen ;
}


void DSMC_PRE_INLET::AllocateMemory( int Num ){
	NumDen	= new double[Num] ;	
	for ( int i=0 ; i<Num ; i++ ) NumDen[i] = 0. ;
}

//==============================================================================================================

DSMC_PRE_RESULT::DSMC_PRE_RESULT(){
	XCenter			= 0. ;
	YCenter			= 0. ;
	ZCenter			= 0. ;
	NumDensity		= 0. ;
	Density			= 0. ;
	XVel			= 0. ;
	YVel			= 0. ;
	ZVel			= 0. ;
	Temp			= 0. ;
	TransTemp		= 0. ;
	RotTemp			= 0. ;
	VibTemp			= 0. ;
	AveParticleNum		= 0. ;
	MeanCollSpacingMeanFreePath = 0. ;
	ProcessorNo		= 0 ;
}

//==============================================================================================================

DSMC_PRE_OPTION::DSMC_PRE_OPTION(){
	OpenScale			= 0 ;
	OpenInputInlet			= 0 ;
	OpenDomainRedecomposition	= 0 ;
	OpenCell			= 0 ;
	
	
	ProcessorNum		= 0 ;	
	Scale			= 1. ; 
	InletFileName		= "NoFile" ;
	InitResultFileName	= "NoFile" ;
	InitMeshFileName	= "NoFile" ;
	CellFileName		= "NoFile" ;
}


void DSMC_PRE_OPTION::Dump(){
	cout << setw(40) << "Processor Number:" << "  " << ProcessorNum << '\n' ;
	cout << setw(40) << "Scale:" << "  " << Scale << '\n' ;
	cout << setw(40) << "Inlet Bounday Condition File:" << "  " << InletFileName << "\n" ;
	cout << setw(40) << "The Last Result File:" << "  " << InitResultFileName << "\n" ;
	cout << setw(40) << "The Last Mesh File:" << "  " << InitMeshFileName << "\n" ;
	cout << setw(40) << "Cell-Information File:" << "  " << CellFileName << "\n\n" ;
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


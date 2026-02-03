#include <fstream>
#include "dsmc_class.h"


using namespace std ;

#if !defined(__DSMC_INIT_H)
#define __DSMC_INIT_H

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
			ofstream		&OutputDebug ) ;

//==============================================================================================================
//==============================================================================================================		

void InitialSpecieslValue( DSMC_DOMAIN *h_pDomain , DSMC_SPECIES *h_Species , DSMC_DSMC *h_pDSMC ) ;

//==============================================================================================================
//==============================================================================================================

void CreateCellInformation(	DSMC_DOMAIN		*h_pDomain ,
				DSMC_NODE		*h_Node ,
				DSMC_CELL		*h_Cell ,
				CELLMAPPING		*h_pMapping ) ;

//==============================================================================================================

void CalculateCellCenter( DSMC_NODE *h_Node , DSMC_CELL *_Cell , CELLMAPPING *h_pMapping ) ;

//==============================================================================================================

void CalculateCellVolume2D( DSMC_NODE *h_Node , DSMC_CELL *_Cell , CELLMAPPING *h_pMapping ) ;


void CalculateCellVolume3D( DSMC_NODE *h_Node , DSMC_CELL *_Cell , CELLMAPPING *h_pMapping ) ;


void CalculateCellVolumeAxisymmetric( DSMC_NODE *h_Node , DSMC_CELL *_Cell , CELLMAPPING *h_pMapping ) ;

//==============================================================================================================

void CreateCellFaceFunction2D( DSMC_NODE *h_Node , DSMC_CELL *_Cell , CELLMAPPING *h_pMapping ) ;


void CreateCellFaceFunction3D( DSMC_NODE *h_Node , DSMC_CELL *_Cell , CELLMAPPING *h_pMapping ) ;

//==============================================================================================================
//==============================================================================================================

void CalculateTimestepWeighting( DSMC_DOMAIN		*h_pDomain ,
				 DSMC_CELL		*h_Cell ,
				 DSMC_SPECIES		*h_Species ,
				 DSMC_DSMC		*h_pDSMC ,
				 CELLMAPPING		*h_pMapping , 
				 ofstream		&OutputDebug ) ;

//==============================================================================================================
//==============================================================================================================

void CreateInletInformation( 	DSMC_DOMAIN		*h_pDomain ,
				DSMC_NODE		*h_Node ,
				DSMC_CELL		*h_Cell ,
				DSMC_INLET		*h_Inlet ,
				DSMC_SPECIES		*h_Species ,
				DSMC_DSMC		*h_pDSMC ,
				CELLMAPPING		*h_pMapping ) ;

//==============================================================================================================			

void CalculateInletFaceArea2D( DSMC_NODE *h_Node , DSMC_INLET *_Inlet ) ;


void CalculateInletFaceArea3D( DSMC_NODE *h_Node , DSMC_INLET *_Inlet ) ;


void CalculateInletFaceAreaAxisymmetric( DSMC_NODE *h_Node , DSMC_INLET *_Inlet ) ;

//==============================================================================================================
//==============================================================================================================

void LinkCellSurface(	DSMC_DOMAIN		*h_pDomain ,
			DSMC_NODE		*h_Node ,
			DSMC_CELL		*h_Cell ,
			DSMC_SURFACE		*h_Surface ,
			CELLMAPPING		*h_pMapping ) ;

//==============================================================================================================

void CalculateSurfaceAreaCenter2D(	DSMC_NODE	*h_Node , 
					DSMC_CELL	*_Cell , 
					DSMC_SURFACE	*_Surface ,
					CELLMAPPING	*h_pMapping , 
					int		FaceNo ) ;
					
					
void CalculateSurfaceAreaCenter3D(	DSMC_NODE	*h_Node , 
					DSMC_CELL	*_Cell , 
					DSMC_SURFACE	*_Surface ,
					CELLMAPPING	*h_pMapping , 
					int		FaceNo ) ;
					
					
void CalculateSurfaceAreaCenterAxisymmetric(	DSMC_NODE	*h_Node , 
						DSMC_CELL	*_Cell , 
						DSMC_SURFACE	*_Surface ,
						CELLMAPPING	*h_pMapping , 
						int		FaceNo ) ;

//==============================================================================================================					
//==============================================================================================================
				
void CreateInitialParticle(	DSMC_DOMAIN		*h_pDomain ,
				DSMC_NODE		*h_Node ,
				DSMC_CELL		*h_Cell ,
				DSMC_SPECIES		*h_Species ,
				DSMC_DSMC		*h_pDSMC ,
				CELLMAPPING		*h_pMapping , 
				ofstream		&OutputDebug ) ;
				
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
				ofstream		&OutputDebug ) ;
				
//==============================================================================================================
				
void AdjustTimestepWeighting(	DSMC_DOMAIN		*h_pDomain ,
				DSMC_CELL		*h_Cell ,
				DSMC_SPECIES		*h_Species ,
				DSMC_DSMC		*h_pDSMC ,
				CELLMAPPING		*h_pMapping , 
				ofstream		&OutputDebug ) ;
				
//==============================================================================================================					
//==============================================================================================================

#endif
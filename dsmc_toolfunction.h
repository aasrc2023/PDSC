#include <string>
#include "dsmc_class.h"

using namespace std ;

#if !defined(__DSMC_TOOLFUNCTION_H)
#define __DSMC_TOOLFUNCTION_H


//==============================================================================================================

double CalculateNormalVelocity2D( double XVel , double YVel , DSMC_NODE *_Node0 , DSMC_NODE *_Node1 ) ;


double CalculateNormalVelocity3D( double XVel , double YVel , double ZVel , DSMC_INLET *_Inlet , DSMC_CELL *_Cell ) ;

//==============================================================================================================

void InverseVelocityToNormPara2D( double *pNormVel , double *pParaVel , double XVel , double YVel , DSMC_NODE Node0 , DSMC_NODE Node1 ) ;


void InverseVelocityToNormPara3D(	double		*pNormVel , 
					double		*pParaVelX , 
					double		*pParaVelY , 
					double		XVel , 
					double		YVel , 
					double		ZVel , 
					DSMC_NODE	*h_Node , 
					DSMC_CELL	*_Cell , 
					int		FaceNo , 
					CELLMAPPING	*h_pMapping ) ;

//==============================================================================================================

void InverseVelocityToCartesian2D( double *pXVel , double *pYVel , double NormVel , double ParaVel , DSMC_NODE Node0 , DSMC_NODE Node1 ) ;


void InverseVelocityToCartesian3D(	double		*pXVel , 
					double		*pYVel , 
					double		*pZVel , 
					double		NormVel , 
					double		ParaVelX , 
					double		ParaVelY , 
					DSMC_NODE	*h_Node , 
					DSMC_CELL	*_Cell , 
					int		FaceNo , 
					CELLMAPPING	*h_pMapping ) ;

//==============================================================================================================

bool InCell2D( DSMC_NODE *h_Node , DSMC_CELL *_Cell , double XCoord , double YCoord , CELLMAPPING *h_pMapping , int Debug ) ;


bool InCell3D( DSMC_NODE *h_Node , DSMC_CELL *_Cell , double XCoord , double YCoord , double ZCoord , CELLMAPPING *h_pMapping , int Debug ) ;
	
//==============================================================================================================

void CreateParticlePosition2D(	double		*XCoord  , 
				double		*YCoord , 
				double		*ZCoord , 
				DSMC_NODE	*h_Node , 
				DSMC_CELL	*_Cell , 
				CELLMAPPING	*pMapping ) ;
				
				
void CreateParticlePosition3D(	double		*XCoord  , 
				double		*YCoord , 
				double		*ZCoord , 
				DSMC_NODE	*h_Node , 
				DSMC_CELL	*_Cell , 
				CELLMAPPING	*h_pMapping ) ;
				
//==============================================================================================================

void CreateParticlePositionEnter3D(	double		*XCoord , 
					double		*YCoord ,
					double		*ZCoord ,
					DSMC_INLET	*_Inlet ,
					DSMC_NODE	*h_Node ,  
					DSMC_CELL	*h_Cell ,
					CELLMAPPING	*h_pMapping ) ;

//==============================================================================================================

double CalculateArea3D( DSMC_NODE *_Node1 , DSMC_NODE *_Node2 , DSMC_NODE *_Node3 ) ;

//==============================================================================================================

void CalculateCellFaceArea3D( double *Area , double *Area1 , double *Area2 , DSMC_NODE *h_Node , DSMC_CELL *_Cell , int FaceNo , CELLMAPPING *h_pMapping ) ;

//==============================================================================================================

void RandVelocity( double *u , double *v , double MostProbableSpeed ) ;


//==============================================================================================================

void CosineLawVelocity( double *NormVel , double *ParaVelX , double *ParaVelY , double CosineLawCoef , double Velocity ) ;

//==============================================================================================================

double RotationalEnergy( double Temp , int RotDOF ) ;

//==============================================================================================================

double VibrationalEnergy( double Temp , int *_ParticleVibLevel , double *_ParticleEffTemp , double VibDOF , DSMC_DOMAIN *h_pDomain , double VibTemp ) ;

//==============================================================================================================

int FindMaxMin( int *Min , int *Max , int *Value , int Num ) ;

//==============================================================================================================

double Randn() ;

//==============================================================================================================

double ErrorFunction( double S ) ;

//==============================================================================================================

double GammaFunction( double x ) ;

//==============================================================================================================

string IntToString( int N ) ;

//==============================================================================================================

void SetInitValue( int *Array , int Value , int Num ) ;

//==============================================================================================================

#endif
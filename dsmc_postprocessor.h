#include <string>
#include "dsmc_class.h"

using namespace std ;

#if !defined(__DSMC_POSTPROCESSOR_H)
#define __DSMC_POSTPROCESSOR_H

//==============================================================================================================
//==============================================================================================================

void OutptuResultTec( DSMC_DOMAIN *pDomain , DSMC_NODE *Node , DSMC_CELL *Cell , DSMC_RESULT *pResult , CELLMAPPING *pMapping ) ;

//==============================================================================================================

void OutptuResultSpeciesTec( DSMC_DOMAIN *pDomain , DSMC_NODE *Node , DSMC_CELL *Cell , DSMC_RESULT *pResult , CELLMAPPING *pMapping , string SpeciesName ) ;
	
//==============================================================================================================

void OutptuResultTec2D( DSMC_DOMAIN *pDomain , DSMC_NODE *Node , DSMC_CELL *Cell , DSMC_RESULT *pResult , CELLMAPPING *pMapping ) ;


void OutptuResultTec3D( DSMC_DOMAIN *pDomain , DSMC_NODE *Node , DSMC_CELL *Cell , DSMC_RESULT *pResult , CELLMAPPING *pMapping ) ;

//==============================================================================================================

void OutptuResultSpeciesTec2D( DSMC_DOMAIN *pDomain , DSMC_NODE *Node , DSMC_CELL *Cell , DSMC_RESULT *pResult , CELLMAPPING *pMapping , string SpeciesName ) ;


void OutptuResultSpeciesTec3D( DSMC_DOMAIN *pDomain , DSMC_NODE *Node , DSMC_CELL *Cell , DSMC_RESULT *pResult , CELLMAPPING *pMapping , string SpeciesName ) ;

//==============================================================================================================
//==============================================================================================================

void RemoveSurface( DSMC_DOMAIN *pDomain , DSMC_POST_SURFACE *Surface , DSMC_POST_SURFACE_CONDITION *pSurfaceCondition ) ;

//==============================================================================================================

void CalculateAngle( DSMC_DOMAIN *pDomain , DSMC_POST_SURFACE *Surface ) ;

//==============================================================================================================

double CalculateAngleDeg( double a , double b , double c ) ;

//==============================================================================================================

void SortSurface( DSMC_DOMAIN *pDomain , DSMC_POST_SURFACE *Surface , DSMC_POST_SURFACE_CONDITION *pSurfaceCondition ) ;

//==============================================================================================================

void CalculateSurfaceProperty( DSMC_DOMAIN *pDomain , DSMC_POST_SURFACE *Surface , DSMC_POST_SURFACE_CONDITION *pSurfaceCondition ) ;

//==============================================================================================================

void CalculateSurfaceProperty2D( DSMC_DOMAIN *pDomain , DSMC_POST_SURFACE *Surface , DSMC_POST_SURFACE_CONDITION *pSurfaceCondition ) ;


void CalculateSurfaceProperty3D( DSMC_DOMAIN *pDomain , DSMC_POST_SURFACE *Surface , DSMC_POST_SURFACE_CONDITION *pSurfaceCondition ) ;

//==============================================================================================================

void OutputSurface( DSMC_NODE	*Node , DSMC_CELL *Cell , DSMC_DOMAIN *pDomain , DSMC_POST_SURFACE *Surface , CELLMAPPING *pMapping ) ;

//==============================================================================================================

void OutputSurface2D( DSMC_DOMAIN *pDomain , DSMC_POST_SURFACE *Surface ) ;


void OutputSurface3D( DSMC_NODE	*Node , DSMC_CELL *Cell , DSMC_POST_SURFACE *Surface , DSMC_DOMAIN *pDomain , CELLMAPPING *pMapping ) ;

//==============================================================================================================
//==============================================================================================================

void OutputCellTec( DSMC_DOMAIN *pDomain , DSMC_NODE *Node , DSMC_CELL *Cell , CELLMAPPING *pMapping ) ;

//==============================================================================================================

void OutputCellTec2D( DSMC_DOMAIN *pDomain , DSMC_NODE *Node , DSMC_CELL *Cell , CELLMAPPING *pMapping ) ;


void OutputCellTec3D( DSMC_DOMAIN *pDomain , DSMC_NODE *Node , DSMC_CELL *Cell , CELLMAPPING *pMapping ) ;

//==============================================================================================================
//==============================================================================================================

void ReadResult( DSMC_DOMAIN *pDomain , DSMC_CELL *Cell , DSMC_RESULT *pResult , string Filename ) ;

//==============================================================================================================

void ReadResultSpecies( DSMC_DOMAIN *pDomain , DSMC_CELL *Cell , DSMC_RESULT *pResult , string Filename ) ;

//==============================================================================================================

void ReadResultSurface( DSMC_DOMAIN *pDomain , DSMC_POST_SURFACE *Surface , string Filename ) ;

//==============================================================================================================

void ReadCellInformation( DSMC_DOMAIN *pDomain , DSMC_CELL *Cell , string Filename ) ;

//==============================================================================================================
//==============================================================================================================

#endif
#include "dsmc_class.h"

#if !defined(__DSMC_OUTPUT_H)
#define __DSMC_OUTPUT_H

//==============================================================================================================
//==============================================================================================================

void OutputResult(	DSMC_DOMAIN		*h_pDomain ,
			DSMC_CELL		*h_Cell ,
			DSMC_WALLTYPE		*h_WallType ,
			DSMC_SURFACE		*h_Surface ,
			DSMC_DSMC		*h_pDSMC ,
			DSMC_RESULT		*h_pResult ,
			DSMC_SPECIES		*h_Species ,
			int 			TimestepNo ) ;
			
//==============================================================================================================

void OutputResultFlowfield(	DSMC_DOMAIN		*h_pDomain ,
				DSMC_CELL		*h_Cell ,
				DSMC_RESULT		*h_pResult ,
				int 			TimestepNo ) ;

//==============================================================================================================
				
void OutputResultFlowfieldSpecies(	DSMC_DOMAIN		*h_pDomain ,
					DSMC_CELL		*h_Cell ,
					DSMC_RESULT		*h_pResult ,
					int 			TimestepNo ) ;				

//==============================================================================================================

void OutputResultSurface(	DSMC_DOMAIN		*h_pDomain ,
				DSMC_CELL		*h_Cell ,
				DSMC_WALLTYPE		*h_WallType ,
				DSMC_SURFACE		*h_Surface ,
				DSMC_DSMC		*h_pDSMC ,
				int			TimestepNo ) ;

//==============================================================================================================

void OutputCellInformation( DSMC_DOMAIN *h_pDomain , DSMC_CELL *h_Cell , DSMC_SPECIES *h_Species ) ;
		
//==============================================================================================================
//==============================================================================================================

void DumpParticlePerCell( DSMC_DOMAIN *h_pDomain , DSMC_DSMC *h_pDSMC , DSMC_CELL *h_Cell , int CellNo ) ;

//==============================================================================================================

void DumpMultiSpecies( DSMC_MULTISPECIES *pMultiSpecies ) ;

//==============================================================================================================

void DumpSimulationInformation( DSMC_DOMAIN *h_pDomain , DSMC_TIME *pTime ) ;

//==============================================================================================================

#endif
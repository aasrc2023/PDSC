#include <fstream>
#include <string>

using namespace std ;

#if !defined(__DSMC_READFILE_H)
#define __DSMC_READFILE_H


//==============================================================================================================
//==============================================================================================================

void ReadInput( DSMC_DOMAIN *pDomain , string Filename ) ;

//==============================================================================================================
//==============================================================================================================

void ReadSimulationCondition(	DSMC_DOMAIN		*h_pDomain ,
				DSMC_PROCESSOR		*h_pProcessor ,
				DSMC_NODE		*h_Node ,
				DSMC_CELL		*h_Cell ,
				DSMC_INLET		*h_Inlet ,
				DSMC_WALLTYPE		*h_WallType ,
				DSMC_SPECIES		*h_Species ,
				CELLMAPPING		*h_pMapping , 
				ofstream		&OutputDebug ) ;

//==============================================================================================================

void ReadPartitionMPI( DSMC_DOMAIN *h_pDomain , DSMC_PROCESSOR *h_pProcessor , string Filename ) ;

//==============================================================================================================

void ReadMeshMPI( DSMC_NODE *h_Node , DSMC_CELL *h_Cell , DSMC_DOMAIN *h_pDomain , DSMC_PROCESSOR *h_pProcessor , CELLMAPPING *h_pMapping , string Filename , ofstream &OutputDebug ) ;


void ReadMesh( DSMC_NODE *h_Node , DSMC_CELL *h_Cell , DSMC_DOMAIN *h_pDomain , CELLMAPPING *h_pMapping , string Filename ) ;

//==============================================================================================================

void ReadInlet( DSMC_INLET *h_Inlet , DSMC_DOMAIN *h_pDomain , string Filename ) ;

//==============================================================================================================

void ReadWallType( DSMC_WALLTYPE *h_WallType , DSMC_DOMAIN *h_pDomain , string Filename ) ;

//==============================================================================================================

void ReadSpecies( DSMC_SPECIES *h_Species , DSMC_DOMAIN *h_pDomain , string Filename , string FileSpeciesName ) ;

//==============================================================================================================

void ReadCellInformation( DSMC_CELL *h_Cell , DSMC_DOMAIN *h_pDomain , DSMC_PROCESSOR *h_pProcessor , string Filename ) ;

//==============================================================================================================
//==============================================================================================================

// Tool Functions.
bool FindString( string Word1 , string Word2 ) ;

void GetSpeciesName( string *SpeciesName , DSMC_DOMAIN *h_pDomain , string Filename ) ;

#endif
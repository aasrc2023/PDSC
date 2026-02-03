#include <string>
#include <metis.h>
#include "dsmc_class_pre.h"

using namespace std ;

#if !defined(__DSMC_PREPROCESSOR_H)
#define __DSMC_PREPROCESSOR_H


void InitOption( DSMC_PRE_OPTION *Option , int argc , char* argv[] ) ;

void GetArraySizeFieldview( int *pNodeNum , int *pCellNum , int *pBCTypeNum , int *pBCFaceNum , string Filename ) ;

bool FindString( string Word1 , string Word2 ) ;

void ReadGridFieldview( int		NodeNum , 
			int		CellNum , 
			int		BCTypeNum , 
			int		BCFaceNum , 
			DSMC_PRE_NODE	*Node , 
			DSMC_PRE_CELL	*Cell , 
			DSMC_PRE_BCTYPE	*BCType , 
			DSMC_PRE_BCFACE	*BCFace ,
			CELLMAPPING	*pMapping , 
			string		Filename ) ;

void ReadBoundaryCondition( int BCTypeNum , DSMC_PRE_BCTYPE *BCType ) ;

bool Check2DMesh( int BCTypeNum , DSMC_PRE_BCTYPE *BCType ) ;

int Transfer3DTo2D(	int		NodeNum , 
			int		CellNum , 
			int		BCFaceNum ,
			DSMC_PRE_NODE	*Node , 
			DSMC_PRE_CELL	*Cell , 
			DSMC_PRE_BCFACE	*BCFace , 
			DSMC_PRE_BCTYPE	*BCType ,  
			CELLMAPPING	*pMapping ) ;

void NodeCCW( DSMC_PRE_NODE *Node , DSMC_PRE_CELL *_Cell , CELLMAPPING *pMapping ) ;

double CalculateAngleDeg( double a , double b , double c ) ;

template <class TYPE>
void BubbleSort( TYPE *a , int *Index , int Num ) ;

template <class TYPE>
void Swap( TYPE *a , TYPE *b ) ;

void CreateNodeCellRelation( int CellNum , DSMC_PRE_CELL *Cell , DSMC_PRE_NODECELL *NodeCell , CELLMAPPING *pMapping ) ;

void CreateNodeFaceRelation( int BCFaceNum , DSMC_PRE_BCFACE *BCFace , DSMC_PRE_BCTYPE *BCType , DSMC_PRE_NODECELL *NodeCell , CELLMAPPING *pMapping ) ;

void DomainPartition(		int			CellNum , 
				int			NodeNum , 
				DSMC_PRE_NODE		*Node ,
				DSMC_PRE_CELL		*Cell , 
				DSMC_PRE_PROCESSOR	*Processor ,
				DSMC_PRE_OPTION		*Option ,
				CELLMAPPING		*pMapping ) ;

int CreateGraph( int			CellNum , 
		 int			NodeNum , 
		 DSMC_PRE_CELL		*Cell ,
		 DSMC_PRE_GRAPHCELL	*GraphCell ,
		 CELLMAPPING		*pMapping ) ;
		 
void ResortCellNo(	DSMC_PRE_CELL		*Cell , 
			DSMC_PRE_PROCESSOR	*Processor , 
			int			CellNum , 
			int			ProcessorNum , 
			idx_t			*MetisCellProcessorNo ) ;

void CreateCellNeighbor( 	int		CellNum , 
				int		NodeNum , 
				int		BCFaceNum ,
				DSMC_PRE_NODE	*Node ,
				DSMC_PRE_CELL	*Cell , 
				DSMC_PRE_BCFACE	*BCFace ,
				DSMC_PRE_BCTYPE	*BCType ,
				CELLMAPPING	*pMapping ) ;

bool CheckNeighbor( DSMC_PRE_CELL *_CellA , DSMC_PRE_CELL *_CellB , int SurfaceA , int SurfaceB , CELLMAPPING *pMapping ) ;

bool CheckNeighbor( DSMC_PRE_CELL *_CellA , DSMC_PRE_BCFACE *_BCFace, int SurfaceA , CELLMAPPING *pMapping ) ;

void CreateProcessorNeighbor( int ProcessorNum , DSMC_PRE_CELL *Cell , DSMC_PRE_PROCESSOR *Processor , CELLMAPPING *pMapping ) ;

void CalculateBoundaryNum( int BCTypeNum , DSMC_PRE_BCTYPE *BCType , int *pBCIONum , int *pBCWallNum ) ;

void CreateNonUniformInletBC( int BCFaceNum , int *pInletSpecifiedNumDenNum , DSMC_PRE_NODE *Node , DSMC_PRE_BCFACE *BCFace , DSMC_PRE_BCTYPE *BCType , string FileName ) ;

int GetInletNo( DSMC_PRE_BCFACE *_BCFace , DSMC_PRE_INLET *Inlet , int InletNum ) ;

void OutputDSMCInputFile(	int 			NodeNum , 
				int			CellNum ,  
				int			BCFaceNum , 
				int			BCIONum ,
				int			ProcessorNum , 
				int			InletSpecifiedNumDenNum ,
				DSMC_PRE_NODE		*Node , 
				DSMC_PRE_CELL		*Cell , 
				DSMC_PRE_BCFACE		*BCFace , 
				DSMC_PRE_BCTYPE		*BCType ,
				DSMC_PRE_PROCESSOR	*Processor ,
				CELLMAPPING 		*pMapping ,
				DSMC_PRE_OPTION		*pOption ) ;

void OutputInletTec(	int			NodeNum ,
			int			BCFaceNum , 
			int			InletSpecifiedNumDenNum ,
			DSMC_PRE_NODE		*Node , 
			DSMC_PRE_BCFACE		*BCFace , 
			DSMC_PRE_BCTYPE		*BCType ,
			DSMC_PRE_OPTION		*pOption , 
			CELLMAPPING 		*pMapping ) ;			

void OutputMeshTec( int NodeNum , int CellNum , DSMC_PRE_NODE *Node , DSMC_PRE_CELL *Cell , CELLMAPPING *pMapping , int Dimension , double Scale ) ;

void ReadResult( int ResultNum , DSMC_PRE_CELL *Cell , DSMC_PRE_RESULT *Result , string Filename ) ;

void ReadMesh( int CellNum , DSMC_PRE_NODE *Node , DSMC_PRE_CELL *Cell , CELLMAPPING *pMapping , string Filename ) ;

void ReadCellInformation( int CellNum , DSMC_PRE_CELL *Cell , string Filename ) ;

int GetFileNum( string Filename ) ;

#endif

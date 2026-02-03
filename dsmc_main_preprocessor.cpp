#include <cstdlib>
#include <iostream>
#include <iomanip>
#include "dsmc_class_pre.h"
#include "dsmc_preprocessor.h"

using namespace std ;

int main( int argc , char* argv[] ){
	int			NodeNum , CellNum , BCTypeNum , BCFaceNum , BCIONum , BCWallNum , Dimension , ResultNum , InletSpecifiedNumDenNum ;
	DSMC_PRE_NODE		*Node ;
	DSMC_PRE_CELL		*Cell ;	
	DSMC_PRE_BCTYPE		*BCType ;
	DSMC_PRE_BCFACE		*BCFace ;
	DSMC_PRE_PROCESSOR	*Processor ;
	DSMC_PRE_RESULT		*Result ;
	CELLMAPPING		Mapping ;
	DSMC_PRE_OPTION		Option ;


	Dimension		= 3 ;
	ResultNum		= 0 ;
	InletSpecifiedNumDenNum	= 1 ;

	
	InitOption( &Option , argc , argv ) ;
	Option.Dump() ;

	if ( Option.ProcessorNum == 0 ){
		cout << "Please input processor (cpu) number: np=<processor number>\n\n" ;
		exit(1) ;	
	}


	// To get the array size from grids information (fieldview format).
	GetArraySizeFieldview( &NodeNum , &CellNum , &BCTypeNum , &BCFaceNum , argv[1] ) ;
	
	
	if ( Option.OpenDomainRedecomposition == 1 )
		ResultNum	= GetFileNum( Option.InitResultFileName ) - 1 ;
	

	cout << setw(20) << "Node Num:"		<< setw(12) << NodeNum << '\n' ;
	cout << setw(20) << "Cell Num:" 	<< setw(12) << CellNum << '\n' ;
	cout << setw(20) << "Boundary Num:" 	<< setw(12) << BCTypeNum << '\n' ;
	cout << setw(20) << "Face Num:" 	<< setw(12) << BCFaceNum << '\n' ;
	cout << setw(20) << "Result Num:"	<< setw(12) << ResultNum << '\n' ;
	cout << "===============================================================\n\n" ;
	//exit(1) ;


	// To allocate the memory.
	Node		= new DSMC_PRE_NODE[NodeNum] ;
	Cell		= new DSMC_PRE_CELL[CellNum] ;
	BCType		= new DSMC_PRE_BCTYPE[BCTypeNum] ;
	BCFace		= new DSMC_PRE_BCFACE[BCFaceNum] ;
	Processor	= new DSMC_PRE_PROCESSOR[Option.ProcessorNum] ;
	Result		= new DSMC_PRE_RESULT[ResultNum] ;


	// Read Girds Information, including node, cell and face.
	ReadGridFieldview( NodeNum , CellNum , BCTypeNum , BCFaceNum , Node , Cell , BCType , BCFace , &Mapping , argv[1] ) ;


	// Read Boundary Conditions (Boundary.txt).
	ReadBoundaryCondition( BCTypeNum , BCType ) ;


	// Debug.
	//for ( int i=0 ; i<BCTypeNum ; i++ ) BCType[i].Dump(i) ;
	//exit(1) ;
	
	
	if ( Check2DMesh( BCTypeNum , BCType ) ){
		cout << "Transfer 3D Mesh to 2D........\n" ;
		NodeNum = Transfer3DTo2D( NodeNum , CellNum , BCFaceNum , Node , Cell , BCFace , BCType , &Mapping ) ;
		Dimension = 2 ;
	}
	
	
	if ( Option.OpenDomainRedecomposition == 1 ){
		ReadMesh( CellNum , Node , Cell , &Mapping , Option.InitMeshFileName ) ;
		ReadResult( ResultNum , Cell , Result , Option.InitResultFileName ) ;
	}
	
	
	if ( Option.OpenCell == 1 ){
		ReadMesh( CellNum , Node , Cell , &Mapping , Option.InitMeshFileName ) ;
		ReadCellInformation( CellNum , Cell , Option.CellFileName ) ;
	}


	cout << "Domain Partition.................\n" ;
	DomainPartition( CellNum , NodeNum , Node , Cell , Processor , &Option , &Mapping ) ;
	

	cout << "Creating the neighbor cells of each cell.................\n" ;
	CreateCellNeighbor( CellNum , NodeNum , BCFaceNum , Node , Cell , BCFace , BCType , &Mapping ) ;
	

	cout << "Creating the neighbor processors of each processor.......\n" ;
	CreateProcessorNeighbor( Option.ProcessorNum , Cell , Processor , &Mapping ) ;


	cout << "===============================================================\n" ;
	cout << setw(20) << "Boundary Name" << setw(8) << "Type" << setw(16) << "Face Number\n" ;
	for ( int i=0 ; i<BCTypeNum ; i++ ) BCType[i].Dump() ;

	
	CalculateBoundaryNum( BCTypeNum , BCType , &BCIONum , &BCWallNum ) ;

	
	// Create non-uniform inlet boundary condition.
	if ( Option.OpenInputInlet == 1 ){
		cout << "Creating a non-uniform inlet boundary condition................\n" ;
		CreateNonUniformInletBC( BCFaceNum , &InletSpecifiedNumDenNum , Node , BCFace , BCType , Option.InletFileName ) ;
	}


	// Output DSMC inptufile (Mesh.inp and Inlet.inp).
	cout << "Outputting DSMC Inputfiles (Mesh.inp, Inlet.inp, and Partition.inp)..........\n" ;
	OutputDSMCInputFile( NodeNum , CellNum , BCFaceNum , BCIONum , Option.ProcessorNum , InletSpecifiedNumDenNum , Node , Cell , BCFace , BCType , Processor , &Mapping , &Option ) ;


	// Debug.
	OutputMeshTec( NodeNum , CellNum , Node , Cell , &Mapping , Dimension , Option.Scale ) ;
	if ( Dimension == 3 ){
		cout << "Outputting Specified Inlet Boundary Conditions (Inlet-tec.dat).......\n" ;
		OutputInletTec( NodeNum , BCFaceNum , InletSpecifiedNumDenNum , Node , BCFace , BCType , &Option , &Mapping ) ;
	}



	cout << "===============================================================\n" ;
	cout << "********************************************************\n" ;
	cout << "********************************************************\n" ;
	cout << setw(18) << "Partition Num:"	<< setw(13) << Option.ProcessorNum << '\n' ;
	cout << setw(18) << "Node Num:"		<< setw(13) << NodeNum << '\n' ;
	cout << setw(18) << "Cell Num:"		<< setw(13) << CellNum << '\n' ;
	cout << setw(18) << "BCType Num:"	<< setw(13) << BCTypeNum << '\n' ;
	cout << setw(18) << "BCFace Num:"	<< setw(13) << (BCIONum+BCWallNum) << '\n' ;
	cout << setw(18) << "I/O Face Num:"	<< setw(13) << BCIONum << '\n' ;
	cout << setw(18) << "Wall Face Num:"	<< setw(13) << BCWallNum << '\n' ;
	cout << "===============================================================\n" ;


	delete [] Node ;
	delete [] Cell ;
	delete [] BCType ;
	delete [] BCFace ;
	delete [] Processor ;
	delete [] Result ;

	return	0 ;
}

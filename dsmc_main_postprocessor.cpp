#include <cstdlib>
#include <iostream>
#include <iomanip>
#include "dsmc_class.h"
#include "dsmc_readfile.h"
#include "dsmc_postprocessor.h"


using namespace std ;

int main( int argc , char** argv ){
	int				ResultType ;
	DSMC_DOMAIN			Domain ;
	DSMC_NODE			*Node ;
	DSMC_CELL			*Cell ;
	DSMC_RESULT			Result ;
	DSMC_POST_SURFACE		*Surface ;
	CELLMAPPING			Mapping ;
	DSMC_POST_SURFACE_CONDITION	SurfaceCondition( "Surface-post.txt" ) ;
	string				SpeciesName ;
	
	
	if ( argc != 2 ){
		cout << '\n' ;
		cout << "Post-process DSMC results. Usage: \n"
				"     ./dsmc_post <ResultFileName>\n\n" ;
		exit(1) ;
	}
	cout << "The reslut is \n" ;
	cout << "    1: Flow field result.\n" ;
	cout << "    2: Flow field result for one spaces.\n" ;
	cout << "    3: Surface result.\n" ;
	cout << "    4: Output cell information.\n" ;
	cin >> ResultType ;
	if ( ResultType != 1 && ResultType != 2 && ResultType != 3 && ResultType != 4 ){
		cout << "Error. Please input \"1, 2, or 3\n" ;
		exit(1) ;
	}
	
	// Read input file.
	ReadInput( &Domain , "Input.txt" ) ;
	if ( Domain.Dimension == 4 ) Domain.Dimension = 2 ;
	
	
	Domain.CellNum	= Domain.TotalCellNum ;

	
	// Allocate memory for nodes and cells.
	Node		= new DSMC_NODE[Domain.NodeNum] ;
	Cell		= new DSMC_CELL[Domain.CellNum] ;
	Result.AllocateMemory( Domain.CellNum , Domain.SpeciesNum ) ;
	Result.InitValue( Domain.CellNum , Domain.SpeciesNum , Domain.Temp ) ;
	
	Surface		= new DSMC_POST_SURFACE[Domain.WallFaceNum] ;

	
	// Read mesh file.
	cout << "Reading mesh information\n\n" ;
	ReadMesh( Node , Cell , &Domain , &Mapping , "Mesh.inp" ) ;
	
	
	if ( ResultType == 1 ){
		cout << "Reading simulation result information\n\n" ;
		ReadResult( &Domain , Cell , &Result , argv[1] ) ;
		
		cout << "Output the simulation result (Tecplot format: Result-tec.dat)\n\n" ;
		OutptuResultTec( &Domain , Node , Cell , &Result , &Mapping ) ;
		
	}else if ( ResultType == 2 ){
		cout << "Please input species name\n" ;
		cin >> SpeciesName ;
		cout << '\n' ;
		
		
		cout << "Reading simulation result information\n\n" ;
		ReadResultSpecies( &Domain , Cell , &Result , argv[1] ) ;
		
		
		cout << "Output the simulation result (Tecplot format: Result-SpeciesName-tec.dat)\n\n" ;
		OutptuResultSpeciesTec( &Domain , Node , Cell , &Result , &Mapping , SpeciesName ) ;
		
	}else if ( ResultType == 3 ){
		cout << "Reading sampling surface properties\n\n" ;
		ReadResultSurface( &Domain , Surface , argv[1] ) ;
		
		RemoveSurface( &Domain , Surface , &SurfaceCondition ) ;
		
		SurfaceCondition.Dump( Domain.WallFaceNum ) ;
		
		
		if ( SurfaceCondition.SortName == "Angle" ) CalculateAngle( &Domain , Surface ) ;
		
		if ( SurfaceCondition.SortName != "3D" ) SortSurface( &Domain , Surface , &SurfaceCondition ) ;
	
	
		cout << "Calculating coefficient (Cp, Cf, and Ch) and drag\n\n" ;
		CalculateSurfaceProperty( &Domain , Surface , &SurfaceCondition ) ;
	
	
		cout << "Output surface propertices (Result-Surface.dat)\n\n" ;
		OutputSurface( Node , Cell , &Domain , Surface , &Mapping ) ;
		
		
		cout << "===========================\n" ;
		cout << "Results: \n" ;
		cout << setw(20) << "Drag (N):"		<< setw(20) << SurfaceCondition.Drag << '\n' ;
		cout << setw(20) << "Heat Flux (W/m^2):"<< setw(20) << SurfaceCondition.HeatFlux << '\n' ;
		cout << setw(20) << "CA:"<< setw(20) << SurfaceCondition.CoeffAxial << '\n' ;
		cout << setw(20) << "CN:"<< setw(20) << SurfaceCondition.CoeffNomal << '\n' ;
		
	}else if ( ResultType == 4 ) {
		cout << "Reading cell information\n\n" ;
		ReadCellInformation( &Domain , Cell , argv[1] ) ;
		
		cout << "Output the cell information (Tecplot format: Cell-tec.dat)\n\n" ;
		OutputCellTec( &Domain , Node , Cell , &Mapping ) ;
		
	}

	
	delete [] Node ;
	delete [] Cell ;
	Result.DeleteMemory( Domain.SpeciesNum ) ;
	
	delete [] Surface ;
	
	cout << "\n=================== End of Post-processing ===================\n" ;
}


#include <mpi.h>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include "dsmc_class.h"
#include "dsmc_readfile.h"
#include "dsmc_init.h"
#include "dsmc_function.h"
#include "dsmc_output.h"
#include "dsmc_toolfunction.h"


using namespace std ;

int main( int argc , char** argv ){
	MPI_Init( &argc , &argv ) ;
	
	int			MPISize , MPIMyID , Debug = 0 ;
	int			TotalParticleNum , *ParticleNum , EnterParticleNum , MaxParticleNum , MinParticleNum , AveParticleNum ;
	int			TotalTimestepNo , SimulationStage , BufferTotalTimestep ;
	DSMC_DOMAIN		h_Domain ;
	DSMC_PROCESSOR		h_Processor ;
	DSMC_NODE		*h_Node ;
	DSMC_CELL		*h_Cell ;
	DSMC_INLET		*h_Inlet ;
	DSMC_WALLTYPE		*h_WallType ;
	DSMC_SPECIES		*h_Species ;
	DSMC_SURFACE		*h_Surface ;
	DSMC_DSMC		h_DSMC ;
	DSMC_RESULT		h_Result ;
	CELLMAPPING		h_Mapping ;
	DSMC_MPI_PARTICLE	*h_MPIParticleOut , *h_MPIParticleIn ;
	DSMC_TIME		Time ;
	DSMC_MPI_DATATYPE	MPIDataType ;
	DSMC_CONVERGENCE	h_Convergence ;
	ofstream		OutputParticleNum , OutputInfo , OutputConvergence , OutputDebug ;
	string			ProcessorID , StageNo , FileName ;
	DSMC_CHEMICALTCE *h_ChemicalTCE ;
	
  int MCLASS ;

	MPI_Barrier( MPI_COMM_WORLD ) ;
	Time.Start() ;
	
	// Obtain the processor number and processor ID.
	MPI_Comm_size( MPI_COMM_WORLD , &MPISize ) ;
	MPI_Comm_rank( MPI_COMM_WORLD , &MPIMyID ) ;
	
	Debug = 0 ;
	
	ProcessorID	= IntToString( MPIMyID ) ;
	FileName	= "Debug-" + ProcessorID + ".dat" ;
	if ( Debug == 1 )
		OutputDebug.open( FileName.c_str() , ios::out | ios::trunc ) ;
	
	
	// Create a MPI Type for DSMC_DOMAIN.
	MPIDataType.InitMPIDataType( &h_Domain ) ;
	

	// Read simulation conditions from the file (Input.txt) for Processor-0.
	if ( MPIMyID == 0 )
		ReadInput( &h_Domain , "Input.txt" ) ;
		
	h_Domain.ChemicalNum = 0 ;

	MPI_Bcast( &h_Domain , 1 , MPIDataType.MPI_DOMAIN , 0 , MPI_COMM_WORLD ) ;


	/*for ( int i=0 ; i<MPISize ; i++ ){
		if ( MPIMyID == i ){
			cout << "MyID: " << i << '\n' ;
			cout << "***************************\n" ;
			h_Domain.Dump() ;
		}
		MPI_Barrier( MPI_COMM_WORLD ) ;
	}
	MPI_Barrier( MPI_COMM_WORLD ) ;
	exit(1) ;*/
	
	MCLASS = 15 ;// M.C.modifify
	
	// Allocate Memory for Node, Cell, Boundary, Particle, Species,and Sampling...etc.
	h_Node			= new DSMC_NODE[h_Domain.NodeNum] ;
	h_Cell			= new DSMC_CELL[h_Domain.LocalCellNum] ;
	h_Inlet			= new DSMC_INLET[h_Domain.InletFaceNum] ;
	h_WallType		= new DSMC_WALLTYPE[h_Domain.WallTypeNum] ;
	h_Species		= new DSMC_SPECIES[h_Domain.SpeciesNum] ;
	h_Surface		= new DSMC_SURFACE[h_Domain.WallFaceNum] ;
	ParticleNum		= new int[MPISize] ;
	
	h_ChemicalTCE		= new DSMC_CHEMICALTCE[MCLASS] ;
	
	
	for ( int i=0 ; i<h_Domain.SpeciesNum ; i++ )
		h_Species[i].AllocateMemory( h_Domain.SpeciesNum ) ;
		
	for ( int i=0 ; i<h_Domain.InletFaceNum ; i++ )
		h_Inlet[i].AllocateMemory( h_Domain.SpeciesNum ) ;
		
	for ( int i=0 ; i<MCLASS; i++ )
		h_ChemicalTCE[i].AllocateMemory(h_Domain.SpeciesNum ) ;
		
	
	h_DSMC.AllocateMemory( h_Domain.MaxParticleNum , h_Domain.LocalCellNum , h_Domain.WallFaceNum , h_Domain.InletFaceNum , h_Domain.SpeciesNum , h_Domain.SpeciesGroupNum ) ;
	h_DSMC.InitValue( h_Domain.MaxParticleNum , h_Domain.LocalCellNum , h_Domain.WallFaceNum , h_Domain.InletFaceNum , h_Domain.SpeciesNum , h_Domain.SpeciesGroupNum ) ;
	
	h_Result.AllocateMemory( h_Domain.LocalCellNum , h_Domain.SpeciesNum ) ;
	h_Result.InitValue( h_Domain.LocalCellNum , h_Domain.SpeciesNum , h_Domain.Temp ) ;
	
	h_Processor.AllocateMemory( h_Domain.TotalCellNum ) ;
	h_Processor.InitValue( h_Domain.TotalCellNum ) ;
	
	h_Convergence.AllocateMemory() ;
	h_Convergence.InitValue( h_Domain.ParticleNumCheck , h_Domain.SpeedCheck , h_Domain.TempCheck ) ;
	// End of allocated memory.
	
	
	// Read simulation conditions, including partition, mesh, inlet/oulet face, solid face, and species information 
	// (Partition.inp, Mesh.inp, Inlet.inp, WallType.txt, and Species.txt).
	ReadSimulationCondition( &h_Domain , &h_Processor , h_Node , h_Cell , h_Inlet , h_WallType , h_Species , &h_Mapping , OutputDebug ) ;




	// Allocate memory for communication of all processors.
	h_MPIParticleOut	= new DSMC_MPI_PARTICLE[h_Domain.TransferParticleNum] ;
	h_MPIParticleIn		= new DSMC_MPI_PARTICLE[h_Domain.TransferParticleNum] ;
	
	// Create a MPI Type for DSMC_DOMAIN.
	MPIDataType.InitMPIDataType( h_MPIParticleOut ) ;
	// End of allocated memory and MPI type.

	
	// Initialized the flow field base on the simulation conditions.
	//Time.Time() ;
	Initialization( &h_Domain , h_Node , h_Cell , h_Inlet , h_WallType , h_Species , h_Surface , &h_DSMC , &h_Mapping , OutputDebug ) ;
	//Time.Time( &Time.Init ) ;


	TotalTimestepNo	= 0 ;
	SimulationStage	= 0 ;
	
//--------------------------------------------------------	
/*
	h_Species[0].ReactionColClass[0] = 0 ;	
	
	h_ChemicalTCE[0].ReactionClassNum = 1 ;	
	
	h_ChemicalTCE[0].ReactionType[0] = 1 ;	
	
	h_ChemicalTCE[0].Pre1Species[0] = 0 ;
	h_ChemicalTCE[0].Pre2Species[0] = 0 ;
	h_ChemicalTCE[0].Pre3Species[0] = -1 ;
	
	h_ChemicalTCE[0].Post1Species[0] = 3 ;
	h_ChemicalTCE[0].Post2Species[0] = 0 ;
	h_ChemicalTCE[0].Post3Species[0] = 3 ;	
	
	h_ChemicalTCE[0].ArrheniusConstant[0] =3.188E-13 ;
	h_ChemicalTCE[0].ArrheniusTempExp[0] = -0.5 ;	
	h_ChemicalTCE[0].ArrheniusActiveEnergy[0] = 1.561E-18 ;
	h_ChemicalTCE[0].HeatofReaction[0] = -1.561E-18 ;
	
//--------------------------------------------------------		

	h_Species[0].ReactionColClass[1] = 1 ;	
	h_Species[1].ReactionColClass[0] = 1 ;
	
	h_ChemicalTCE[1].ReactionClassNum = 2 ;	
	
	h_ChemicalTCE[1].ReactionType[0] = 1 ;
		
	h_ChemicalTCE[1].Pre1Species[0] = 0 ;
	h_ChemicalTCE[1].Pre2Species[0] = 1 ;
	h_ChemicalTCE[1].Pre3Species[0] = -1 ;
	
	h_ChemicalTCE[1].Post1Species[0] = 3 ;
	h_ChemicalTCE[1].Post2Species[0] = 1 ;
	h_ChemicalTCE[1].Post3Species[0] = 3 ;
	
	h_ChemicalTCE[1].ArrheniusConstant[0] = 3.188E-13;
	h_ChemicalTCE[1].ArrheniusTempExp[0] = -0.5 ;	
	h_ChemicalTCE[1].ArrheniusActiveEnergy[0] = 1.561E-18 ;
	h_ChemicalTCE[1].HeatofReaction[0] = -1.561E-18 ;
	
	h_ChemicalTCE[1].ReactionType[1] = 1 ;
		
	h_ChemicalTCE[1].Pre1Species[1] = 1 ;
	h_ChemicalTCE[1].Pre2Species[1] = 0 ;
	h_ChemicalTCE[1].Pre3Species[1] = -1 ;
	
	h_ChemicalTCE[1].Post1Species[1] = 4 ;
	h_ChemicalTCE[1].Post2Species[1] = 0 ;
	h_ChemicalTCE[1].Post3Species[1] = 4 ;
	
	h_ChemicalTCE[1].ArrheniusConstant[1] = 5.995E-12 ;
	h_ChemicalTCE[1].ArrheniusTempExp[1] = -1 ;	
	h_ChemicalTCE[1].ArrheniusActiveEnergy[1] = 8.201E-19 ;
	h_ChemicalTCE[1].HeatofReaction[1] = -8.201E-19 ;

//--------------------------------------------------------	
	
	h_Species[0].ReactionColClass[2] = 2 ;	
	h_Species[2].ReactionColClass[0] = 2 ;
	
	h_ChemicalTCE[2].ReactionClassNum = 2 ;	
	
	h_ChemicalTCE[2].ReactionType[0] = 1 ;
		
	h_ChemicalTCE[2].Pre1Species[0] = 0 ;
	h_ChemicalTCE[2].Pre2Species[0] = 2 ;
	h_ChemicalTCE[2].Pre3Species[0] = -1 ;
	
	h_ChemicalTCE[2].Post1Species[0] = 3 ;
	h_ChemicalTCE[2].Post2Species[0] = 2 ;
	h_ChemicalTCE[2].Post3Species[0] = 3 ;
	
	h_ChemicalTCE[2].ArrheniusConstant[0] = 3.188E-13 ;
	h_ChemicalTCE[2].ArrheniusTempExp[0] = -0.5 ;	
	h_ChemicalTCE[2].ArrheniusActiveEnergy[0] = 1.561E-18 ;
	h_ChemicalTCE[2].HeatofReaction[0] = -1.561E-18 ;
	
	h_ChemicalTCE[2].ReactionType[1] = 1 ;
		
	h_ChemicalTCE[2].Pre1Species[1] = 2 ;
	h_ChemicalTCE[2].Pre2Species[1] = 0 ;
	h_ChemicalTCE[2].Pre3Species[1] = -1 ;
	
	h_ChemicalTCE[2].Post1Species[1] = 3 ;
	h_ChemicalTCE[2].Post2Species[1] = 0 ;
	h_ChemicalTCE[2].Post3Species[1] = 4 ;
	
	h_ChemicalTCE[2].ArrheniusConstant[1] = 6.592E-10 ;
	h_ChemicalTCE[2].ArrheniusTempExp[1] = -1.5 ;	
	h_ChemicalTCE[2].ArrheniusActiveEnergy[1] = 1.044E-18 ;
	h_ChemicalTCE[2].HeatofReaction[1] = -1.044E-18 ;
	
//-----------------------------------------------------------

	h_Species[0].ReactionColClass[3] = 3 ;
	
	h_ChemicalTCE[3].ReactionClassNum = 1 ;	
	
	h_ChemicalTCE[3].ReactionType[0] = 1 ;
	
	h_ChemicalTCE[3].Pre1Species[0] = 0 ;
	h_ChemicalTCE[3].Pre2Species[0] = 3 ;
	h_ChemicalTCE[3].Pre3Species[0] = -1 ;
	
	h_ChemicalTCE[3].Post1Species[0] = 3 ;
	h_ChemicalTCE[3].Post2Species[0] = 3 ;
	h_ChemicalTCE[3].Post3Species[0] = 3 ;
	
	h_ChemicalTCE[3].ArrheniusConstant[0] = 6.891E-8 ;
	h_ChemicalTCE[3].ArrheniusTempExp[0] = -1.5 ;	
	h_ChemicalTCE[3].ArrheniusActiveEnergy[0] = 1.561E-18 ;
	h_ChemicalTCE[3].HeatofReaction[0] = -1.561E-18 ;

//--------------------------------------------------------	
	
	h_Species[0].ReactionColClass[4] = 4 ;	
	
	h_ChemicalTCE[4].ReactionClassNum = 2 ;	
	
	h_ChemicalTCE[4].ReactionType[0] = 1 ;
		
	h_ChemicalTCE[4].Pre1Species[0] = 0 ;
	h_ChemicalTCE[4].Pre2Species[0] = 4 ;
	h_ChemicalTCE[4].Pre3Species[0] = -1 ;
	
	h_ChemicalTCE[4].Post1Species[0] = 3 ;
	h_ChemicalTCE[4].Post2Species[0] = 4 ;
	h_ChemicalTCE[4].Post3Species[0] = 3 ;
	
	h_ChemicalTCE[4].ArrheniusConstant[0] = 3.188E-13 ;
	h_ChemicalTCE[4].ArrheniusTempExp[0] = -0.5 ;	
	h_ChemicalTCE[4].ArrheniusActiveEnergy[0] = 1.561E-18 ;
	h_ChemicalTCE[4].HeatofReaction[0] = -1.561E-18 ;
	
	h_ChemicalTCE[4].ReactionType[1] = 2 ;
		
	h_ChemicalTCE[4].Pre1Species[1] = 0 ;
	h_ChemicalTCE[4].Pre2Species[1] = 4 ;
	h_ChemicalTCE[4].Pre3Species[1] = -1 ;
	
	h_ChemicalTCE[4].Post1Species[1] = 3 ;
	h_ChemicalTCE[4].Post2Species[1] = 2 ;
	h_ChemicalTCE[4].Post3Species[1] = -1 ;
	
	h_ChemicalTCE[4].ArrheniusConstant[1] = 1.121E-16 ;
	h_ChemicalTCE[4].ArrheniusTempExp[1] = 0 ;	
	h_ChemicalTCE[4].ArrheniusActiveEnergy[1] = 5.177E-19 ;
	h_ChemicalTCE[4].HeatofReaction[1] = -5.177E-19 ;
	
//-----------------------------------------------------------

	h_Species[1].ReactionColClass[1] = 5 ;
	
	h_ChemicalTCE[5].ReactionClassNum = 1 ;	
	
	h_ChemicalTCE[5].ReactionType[0] = 1 ;
	
	h_ChemicalTCE[5].Pre1Species[0] = 1 ;
	h_ChemicalTCE[5].Pre2Species[0] = 1 ;
	h_ChemicalTCE[5].Pre3Species[0] = -1 ;
	
	h_ChemicalTCE[5].Post1Species[0] = 4 ;
	h_ChemicalTCE[5].Post2Species[0] = 1 ;
	h_ChemicalTCE[5].Post3Species[0] = 4 ;
	
	h_ChemicalTCE[5].ArrheniusConstant[0] = 5.995E-12 ;
	h_ChemicalTCE[5].ArrheniusTempExp[0] = -1 ;	
	h_ChemicalTCE[5].ArrheniusActiveEnergy[0] = 8.201E-19 ;
	h_ChemicalTCE[5].HeatofReaction[0] = -8.201E-19 ;

//--------------------------------------------------------	

	h_Species[1].ReactionColClass[2] = 6 ;	
	h_Species[2].ReactionColClass[1] = 6 ;	
	
	h_ChemicalTCE[6].ReactionClassNum = 2 ;	
	
	h_ChemicalTCE[6].ReactionType[0] = 1 ;
		
	h_ChemicalTCE[6].Pre1Species[0] = 1 ;
	h_ChemicalTCE[6].Pre2Species[0] = 2 ;
	h_ChemicalTCE[6].Pre3Species[0] = -1 ;
	
	h_ChemicalTCE[6].Post1Species[0] = 4 ;
	h_ChemicalTCE[6].Post2Species[0] = 2 ;
	h_ChemicalTCE[6].Post3Species[0] = 4 ;
	
	h_ChemicalTCE[6].ArrheniusConstant[0] = 5.995E-12 ;
	h_ChemicalTCE[6].ArrheniusTempExp[0] = -1 ;	
	h_ChemicalTCE[6].ArrheniusActiveEnergy[0] = 8.201E-19 ;
	h_ChemicalTCE[6].HeatofReaction[0] = -8.201E-19 ;
	
	h_ChemicalTCE[6].ReactionType[1] = 1 ;
		
	h_ChemicalTCE[6].Pre1Species[1] = 2 ;
	h_ChemicalTCE[6].Pre2Species[1] = 1 ;
	h_ChemicalTCE[6].Pre3Species[1] = -1 ;
	
	h_ChemicalTCE[6].Post1Species[1] = 3 ;
	h_ChemicalTCE[6].Post2Species[1] = 1 ;
	h_ChemicalTCE[6].Post3Species[1] = 4 ;
	
	h_ChemicalTCE[6].ArrheniusConstant[1] = 6.592E-10 ;
	h_ChemicalTCE[6].ArrheniusTempExp[1] = -1.5 ;	
	h_ChemicalTCE[6].ArrheniusActiveEnergy[1] = 1.044E-18 ;
	h_ChemicalTCE[6].HeatofReaction[1] = -1.044E-18 ;
	
//-----------------------------------------------------------

	h_Species[1].ReactionColClass[3] = 7 ;	
	
	h_ChemicalTCE[7].ReactionClassNum = 2 ;	
	
	h_ChemicalTCE[7].ReactionType[0] = 1 ;
		
	h_ChemicalTCE[7].Pre1Species[0] = 1 ;
	h_ChemicalTCE[7].Pre2Species[0] = 3 ;
	h_ChemicalTCE[7].Pre3Species[0] = -1 ;
	
	h_ChemicalTCE[7].Post1Species[0] = 4 ;
	h_ChemicalTCE[7].Post2Species[0] = 3 ;
	h_ChemicalTCE[7].Post3Species[0] = 4 ;
	
	h_ChemicalTCE[7].ArrheniusConstant[0] = 5.995E-12 ;
	h_ChemicalTCE[7].ArrheniusTempExp[0] = -1 ;	
	h_ChemicalTCE[7].ArrheniusActiveEnergy[0] = 8.201E-19 ;
	h_ChemicalTCE[7].HeatofReaction[0] = -8.201E-19 ;
	
	h_ChemicalTCE[7].ReactionType[1] = 2 ;
		
	h_ChemicalTCE[7].Pre1Species[1] = 1 ;
	h_ChemicalTCE[7].Pre2Species[1] = 3 ;
	h_ChemicalTCE[7].Pre3Species[1] = -1 ;
	
	h_ChemicalTCE[7].Post1Species[1] = 4 ;
	h_ChemicalTCE[7].Post2Species[1] = 2 ;
	h_ChemicalTCE[7].Post3Species[1] = -1 ;
	
	h_ChemicalTCE[7].ArrheniusConstant[1] = 1.599E-18 ;
	h_ChemicalTCE[7].ArrheniusTempExp[1] = 0.5 ;	
	h_ChemicalTCE[7].ArrheniusActiveEnergy[1] = 4.97E-20 ;
	h_ChemicalTCE[7].HeatofReaction[1] = 2.72E-19 ;

//-----------------------------------------------------------

	h_Species[1].ReactionColClass[4] = 8 ;
	
	h_ChemicalTCE[8].ReactionClassNum = 1 ;	
	
	h_ChemicalTCE[8].ReactionType[0] = 1 ;
	
	h_ChemicalTCE[8].Pre1Species[0] = 1 ;
	h_ChemicalTCE[8].Pre2Species[0] = 4 ;
	h_ChemicalTCE[8].Pre3Species[0] = -1 ;
	
	h_ChemicalTCE[8].Post1Species[0] = 4 ;
	h_ChemicalTCE[8].Post2Species[0] = 4 ;
	h_ChemicalTCE[8].Post3Species[0] = 4 ;
	
	h_ChemicalTCE[8].ArrheniusConstant[0] = 5.995E-12 ;
	h_ChemicalTCE[8].ArrheniusTempExp[0] = -1 ;	
	h_ChemicalTCE[8].ArrheniusActiveEnergy[0] = 8.201E-19 ;
	h_ChemicalTCE[8].HeatofReaction[0] = -8.201E-19 ;

//--------------------------------------------------------	

	h_Species[2].ReactionColClass[2] = 9 ;
	
	h_ChemicalTCE[9].ReactionClassNum = 1 ;	
	
	h_ChemicalTCE[9].ReactionType[0] = 1 ;
	
	h_ChemicalTCE[9].Pre1Species[0] = 2 ;
	h_ChemicalTCE[9].Pre2Species[0] = 2 ;
	h_ChemicalTCE[9].Pre3Species[0] = -1 ;
	
	h_ChemicalTCE[9].Post1Species[0] = 3 ;
	h_ChemicalTCE[9].Post2Species[0] = 2 ;
	h_ChemicalTCE[9].Post3Species[0] = 4 ;
	
	h_ChemicalTCE[9].ArrheniusConstant[0] = 6.592E-10 ;
	h_ChemicalTCE[9].ArrheniusTempExp[0] = -1.5 ;	
	h_ChemicalTCE[9].ArrheniusActiveEnergy[0] = 1.044E-18 ;
	h_ChemicalTCE[9].HeatofReaction[0] = -1.044E-18 ;

//--------------------------------------------------------	

	h_Species[2].ReactionColClass[3] = 10 ;	
	
	h_ChemicalTCE[10].ReactionClassNum = 2 ;	
	
	h_ChemicalTCE[10].ReactionType[0] = 1 ;
		
	h_ChemicalTCE[10].Pre1Species[0] = 2 ;
	h_ChemicalTCE[10].Pre2Species[0] = 3 ;
	h_ChemicalTCE[10].Pre3Species[0] = -1 ;
	
	h_ChemicalTCE[10].Post1Species[0] = 3 ;
	h_ChemicalTCE[10].Post2Species[0] = 3 ;
	h_ChemicalTCE[10].Post3Species[0] = 4 ;
	
	h_ChemicalTCE[10].ArrheniusConstant[0] = 6.592E-10 ;
	h_ChemicalTCE[10].ArrheniusTempExp[0] = -1.5 ;	
	h_ChemicalTCE[10].ArrheniusActiveEnergy[0] = 1.044E-18 ;
	h_ChemicalTCE[10].HeatofReaction[0] = -1.044E-18 ;
	
	h_ChemicalTCE[10].ReactionType[1] = 2 ;
		
	h_ChemicalTCE[10].Pre1Species[1] = 2 ;
	h_ChemicalTCE[10].Pre2Species[1] = 3 ;
	h_ChemicalTCE[10].Pre3Species[1] = -1 ;
	
	h_ChemicalTCE[10].Post1Species[1] = 4 ;
	h_ChemicalTCE[10].Post2Species[1] = 0 ;
	h_ChemicalTCE[10].Post3Species[1] = -1 ;
	
	h_ChemicalTCE[10].ArrheniusConstant[1] = 2.491E-17 ;
	h_ChemicalTCE[10].ArrheniusTempExp[1] = 0 ;	
	h_ChemicalTCE[10].ArrheniusActiveEnergy[1] = 0 ;
	h_ChemicalTCE[10].HeatofReaction[1] = 5.177E-19 ;

//-----------------------------------------------------------
	
	h_Species[2].ReactionColClass[4] = 11 ;	
	
	h_ChemicalTCE[11].ReactionClassNum = 2 ;	
	
	h_ChemicalTCE[11].ReactionType[0] = 1 ;
		
	h_ChemicalTCE[11].Pre1Species[0] = 2 ;
	h_ChemicalTCE[11].Pre2Species[0] = 4 ;
	h_ChemicalTCE[11].Pre3Species[0] = -1 ;
	
	h_ChemicalTCE[11].Post1Species[0] = 3 ;
	h_ChemicalTCE[11].Post2Species[0] = 4 ;
	h_ChemicalTCE[11].Post3Species[0] = 4 ;
	
	h_ChemicalTCE[11].ArrheniusConstant[0] = 6.592E-10 ;
	h_ChemicalTCE[11].ArrheniusTempExp[0] = -1.5 ;	
	h_ChemicalTCE[11].ArrheniusActiveEnergy[0] = 1.044E-18 ;
	h_ChemicalTCE[11].HeatofReaction[0] = -1.044E-18 ;
	
	h_ChemicalTCE[11].ReactionType[1] = 2 ;
		
	h_ChemicalTCE[11].Pre1Species[1] = 2 ;
	h_ChemicalTCE[11].Pre2Species[1] = 4 ;
	h_ChemicalTCE[11].Pre3Species[1] = -1 ;
	
	h_ChemicalTCE[11].Post1Species[1] = 3 ;
	h_ChemicalTCE[11].Post2Species[1] = 1 ;
	h_ChemicalTCE[11].Post3Species[1] = -1 ;
	
	h_ChemicalTCE[11].ArrheniusConstant[1] = 5.281E-21 ;
	h_ChemicalTCE[11].ArrheniusTempExp[1] = 1 ;	
	h_ChemicalTCE[11].ArrheniusActiveEnergy[1] = 2.72E-19 ;
	h_ChemicalTCE[11].HeatofReaction[1] = -2.72E-19 ;

//-----------------------------------------------------------
	
	h_Species[3].ReactionColClass[3] = 12 ;	
	
	h_ChemicalTCE[12].ReactionClassNum = -5 ;	
	
	h_ChemicalTCE[12].ReactionType[0] = 3 ;
		
	h_ChemicalTCE[12].Pre1Species[0] = 3 ;
	h_ChemicalTCE[12].Pre2Species[0] = 3 ;
	h_ChemicalTCE[12].Pre3Species[0] = 0 ;
	
	h_ChemicalTCE[12].Post1Species[0] = 0 ;
	h_ChemicalTCE[12].Post2Species[0] = -1 ;
	h_ChemicalTCE[12].Post3Species[0] = 0 ;
	
	h_ChemicalTCE[12].ArrheniusConstant[0] = 3.006E-44 ;
	h_ChemicalTCE[12].ArrheniusTempExp[0] = -0.5 ;	
	h_ChemicalTCE[12].ArrheniusActiveEnergy[0] = 0 ;
	h_ChemicalTCE[12].HeatofReaction[0] = 1.561E-18 ;
	
	h_ChemicalTCE[12].ReactionType[1] = 3 ;
		
	h_ChemicalTCE[12].Pre1Species[1] = 3 ;
	h_ChemicalTCE[12].Pre2Species[1] = 3 ;
	h_ChemicalTCE[12].Pre3Species[1] = 1 ;
	
	h_ChemicalTCE[12].Post1Species[1] = 0 ;
	h_ChemicalTCE[12].Post2Species[1] = -1 ;
	h_ChemicalTCE[12].Post3Species[1] = 1 ;
	
	h_ChemicalTCE[12].ArrheniusConstant[1] = 3.006E-44 ;
	h_ChemicalTCE[12].ArrheniusTempExp[1] = -0.5 ;	
	h_ChemicalTCE[12].ArrheniusActiveEnergy[1] = 0 ;
	h_ChemicalTCE[12].HeatofReaction[1] = 1.561E-18 ;

	h_ChemicalTCE[12].ReactionType[2] = 3 ;
		
	h_ChemicalTCE[12].Pre1Species[2] = 3 ;
	h_ChemicalTCE[12].Pre2Species[2] = 3 ;
	h_ChemicalTCE[12].Pre3Species[2] = 2 ;
	
	h_ChemicalTCE[12].Post1Species[2] = 0 ;
	h_ChemicalTCE[12].Post2Species[2] = -1 ;
	h_ChemicalTCE[12].Post3Species[2] = 2 ;
	
	h_ChemicalTCE[12].ArrheniusConstant[2] = 3.0066E-44 ;
	h_ChemicalTCE[12].ArrheniusTempExp[2] = -0.5 ;	
	h_ChemicalTCE[12].ArrheniusActiveEnergy[2] = 0 ;
	h_ChemicalTCE[12].HeatofReaction[2] = 1.561E-18 ;

	h_ChemicalTCE[12].ReactionType[3] = 3 ;
		
	h_ChemicalTCE[12].Pre1Species[3] = 3 ;
	h_ChemicalTCE[12].Pre2Species[3] = 3 ;
	h_ChemicalTCE[12].Pre3Species[3] = 3 ;
	
	h_ChemicalTCE[12].Post1Species[3] = 0 ;
	h_ChemicalTCE[12].Post2Species[3] = -1 ;
	h_ChemicalTCE[12].Post3Species[3] = 3 ;
	
	h_ChemicalTCE[12].ArrheniusConstant[3] = 6.397E-39 ;
	h_ChemicalTCE[12].ArrheniusTempExp[3] = -1.5 ;	
	h_ChemicalTCE[12].ArrheniusActiveEnergy[3] = 0 ;
	h_ChemicalTCE[12].HeatofReaction[3] = 1.561E-18 ;
	
	h_ChemicalTCE[12].ReactionType[4] = 3 ;
		
	h_ChemicalTCE[12].Pre1Species[4] = 3 ;
	h_ChemicalTCE[12].Pre2Species[4] = 3 ;
	h_ChemicalTCE[12].Pre3Species[4] = 4 ;
	
	h_ChemicalTCE[12].Post1Species[4] = 0 ;
	h_ChemicalTCE[12].Post2Species[4] = -1 ;
	h_ChemicalTCE[12].Post3Species[4] = 4 ;
	
	h_ChemicalTCE[12].ArrheniusConstant[4] = 3.006E-44 ;
	h_ChemicalTCE[12].ArrheniusTempExp[4] = -0.5 ;	
	h_ChemicalTCE[12].ArrheniusActiveEnergy[4] = 0 ;
	h_ChemicalTCE[12].HeatofReaction[4] = 1.561E-18 ;

//-----------------------------------------------------------
	
	h_Species[4].ReactionColClass[4] = 13 ;	
	
	h_ChemicalTCE[13].ReactionClassNum = -5 ;	
	
	h_ChemicalTCE[13].ReactionType[0] = 3 ;
		
	h_ChemicalTCE[13].Pre1Species[0] = 4 ;
	h_ChemicalTCE[13].Pre2Species[0] = 4 ;
	h_ChemicalTCE[13].Pre3Species[0] = 0 ;
	
	h_ChemicalTCE[13].Post1Species[0] = 1 ;
	h_ChemicalTCE[13].Post2Species[0] = -1 ;
	h_ChemicalTCE[13].Post3Species[0] = 0 ;
	
	h_ChemicalTCE[13].ArrheniusConstant[0] = 8.3E-45 ;
	h_ChemicalTCE[13].ArrheniusTempExp[0] = -0.5 ;	
	h_ChemicalTCE[13].ArrheniusActiveEnergy[0] = 0 ;
	h_ChemicalTCE[13].HeatofReaction[0] = 8.201E-19 ;
	
	h_ChemicalTCE[13].ReactionType[1] = 3 ;
		
	h_ChemicalTCE[13].Pre1Species[1] = 4 ;
	h_ChemicalTCE[13].Pre2Species[1] = 4 ;
	h_ChemicalTCE[13].Pre3Species[1] = 1 ;
	
	h_ChemicalTCE[13].Post1Species[1] = 1 ;
	h_ChemicalTCE[13].Post2Species[1] = -1 ;
	h_ChemicalTCE[13].Post3Species[1] = 1 ;
	
	h_ChemicalTCE[13].ArrheniusConstant[1] =  8.3E-45 ;
	h_ChemicalTCE[13].ArrheniusTempExp[1] = -0.5 ;	
	h_ChemicalTCE[13].ArrheniusActiveEnergy[1] = 0 ;
	h_ChemicalTCE[13].HeatofReaction[1] = 8.201E-19 ;

	h_ChemicalTCE[13].ReactionType[2] = 3 ;
		
	h_ChemicalTCE[13].Pre1Species[2] = 4 ;
	h_ChemicalTCE[13].Pre2Species[2] = 4 ;
	h_ChemicalTCE[13].Pre3Species[2] = 2 ;
	
	h_ChemicalTCE[13].Post1Species[2] = 1 ;
	h_ChemicalTCE[13].Post2Species[2] = -1 ;
	h_ChemicalTCE[13].Post3Species[2] = 2 ;
	
	h_ChemicalTCE[13].ArrheniusConstant[2] =  8.3E-45  ;
	h_ChemicalTCE[13].ArrheniusTempExp[2] = -0.5 ;	
	h_ChemicalTCE[13].ArrheniusActiveEnergy[2] = 0 ;
	h_ChemicalTCE[13].HeatofReaction[2] = 8.201E-19 ;

	h_ChemicalTCE[13].ReactionType[3] = 3 ;
		
	h_ChemicalTCE[13].Pre1Species[3] = 4 ;
	h_ChemicalTCE[13].Pre2Species[3] = 4 ;
	h_ChemicalTCE[13].Pre3Species[3] = 3 ;
	
	h_ChemicalTCE[13].Post1Species[3] = 1 ;
	h_ChemicalTCE[13].Post2Species[3] = -1 ;
	h_ChemicalTCE[13].Post3Species[3] = 3 ;
	
	h_ChemicalTCE[13].ArrheniusConstant[3] =  8.3E-45  ;
	h_ChemicalTCE[13].ArrheniusTempExp[3] = -0.5 ;	
	h_ChemicalTCE[13].ArrheniusActiveEnergy[3] = 0 ;
	h_ChemicalTCE[13].HeatofReaction[3] = 8.201E-19 ;
	
	h_ChemicalTCE[13].ReactionType[4] = 3 ;
		
	h_ChemicalTCE[13].Pre1Species[4] = 4 ;
	h_ChemicalTCE[13].Pre2Species[4] = 4 ;
	h_ChemicalTCE[13].Pre3Species[4] = 4 ;
	
	h_ChemicalTCE[13].Post1Species[4] = 1 ;
	h_ChemicalTCE[13].Post2Species[4] = -1 ;
	h_ChemicalTCE[13].Post3Species[4] = 4 ;
	
	h_ChemicalTCE[13].ArrheniusConstant[4] =  8.3E-45  ;
	h_ChemicalTCE[13].ArrheniusTempExp[4] = -0.5 ;	
	h_ChemicalTCE[13].ArrheniusActiveEnergy[4] = 0 ;
	h_ChemicalTCE[13].HeatofReaction[4] = 8.201E-19 ;

//-----------------------------------------------------------

	h_Species[3].ReactionColClass[4] = 14 ;	
	h_Species[4].ReactionColClass[3] = 14 ;	
	
	h_ChemicalTCE[14].ReactionClassNum = -5 ;	
	
	h_ChemicalTCE[14].ReactionType[0] = 3 ;
		
	h_ChemicalTCE[14].Pre1Species[0] = 3 ;
	h_ChemicalTCE[14].Pre2Species[0] = 4 ;
	h_ChemicalTCE[14].Pre3Species[0] = 0 ;
	
	h_ChemicalTCE[14].Post1Species[0] = 2 ;
	h_ChemicalTCE[14].Post2Species[0] = -1 ;
	h_ChemicalTCE[14].Post3Species[0] = 0 ;
	
	h_ChemicalTCE[14].ArrheniusConstant[0] = 2.785E-40 ;
	h_ChemicalTCE[14].ArrheniusTempExp[0] = -1.5 ;	
	h_ChemicalTCE[14].ArrheniusActiveEnergy[0] = 0 ;
	h_ChemicalTCE[14].HeatofReaction[0] = 1.044E-18 ;
	
	h_ChemicalTCE[14].ReactionType[1] = 3 ;
		
	h_ChemicalTCE[14].Pre1Species[1] = 3 ;
	h_ChemicalTCE[14].Pre2Species[1] = 4 ;
	h_ChemicalTCE[14].Pre3Species[1] = 1 ;
	
	h_ChemicalTCE[14].Post1Species[1] = 2 ;
	h_ChemicalTCE[14].Post2Species[1] = -1 ;
	h_ChemicalTCE[14].Post3Species[1] = 1 ;
	
	h_ChemicalTCE[14].ArrheniusConstant[1] = 2.785E-40 ;
	h_ChemicalTCE[14].ArrheniusTempExp[1] = -1.5 ;	
	h_ChemicalTCE[14].ArrheniusActiveEnergy[1] = 0 ;
	h_ChemicalTCE[14].HeatofReaction[1] = 1.044E-18 ;

	h_ChemicalTCE[14].ReactionType[2] = 3 ;
		
	h_ChemicalTCE[14].Pre1Species[2] = 3 ;
	h_ChemicalTCE[14].Pre2Species[2] = 4 ;
	h_ChemicalTCE[14].Pre3Species[2] = 2 ;
	
	h_ChemicalTCE[14].Post1Species[2] = 2 ;
	h_ChemicalTCE[14].Post2Species[2] = -1 ;
	h_ChemicalTCE[14].Post3Species[2] = 2 ;
	
	h_ChemicalTCE[14].ArrheniusConstant[2] = 2.785E-40 ;
	h_ChemicalTCE[14].ArrheniusTempExp[2] = -1.5 ;	
	h_ChemicalTCE[14].ArrheniusActiveEnergy[2] = 0 ;
	h_ChemicalTCE[14].HeatofReaction[2] = 1.044E-18 ;

	h_ChemicalTCE[14].ReactionType[3] = 3 ;
		
	h_ChemicalTCE[14].Pre1Species[3] = 3 ;
	h_ChemicalTCE[14].Pre2Species[3] = 4 ;
	h_ChemicalTCE[14].Pre3Species[3] = 3 ;
	
	h_ChemicalTCE[14].Post1Species[3] = 2 ;
	h_ChemicalTCE[14].Post2Species[3] = -1 ;
	h_ChemicalTCE[14].Post3Species[3] = 3 ;
	
	h_ChemicalTCE[14].ArrheniusConstant[3] = 2.785E-40 ;
	h_ChemicalTCE[14].ArrheniusTempExp[3] = -1.5 ;	
	h_ChemicalTCE[14].ArrheniusActiveEnergy[3] = 0 ;
	h_ChemicalTCE[14].HeatofReaction[3] = 1.044E-18 ;
	
	h_ChemicalTCE[14].ReactionType[4] = 3 ;
		
	h_ChemicalTCE[14].Pre1Species[4] = 3 ;
	h_ChemicalTCE[14].Pre2Species[4] = 4 ;
	h_ChemicalTCE[14].Pre3Species[4] = 4 ;
	
	h_ChemicalTCE[14].Post1Species[4] = 2 ;
	h_ChemicalTCE[14].Post2Species[4] = -1 ;
	h_ChemicalTCE[14].Post3Species[4] = 4 ;
	
	h_ChemicalTCE[14].ArrheniusConstant[4] = 2.785E-40 ;
	h_ChemicalTCE[14].ArrheniusTempExp[4] = -1.5 ;	
	h_ChemicalTCE[14].ArrheniusActiveEnergy[4] = 0 ;
	h_ChemicalTCE[14].HeatofReaction[4] = 1.044E-18 ;
*/
//-----------------------------------------------------------
//-----------------------------------------------------------	
	do{
		SimulationStage++ ;
		
		StageNo	= IntToString( SimulationStage ) ;
		
		if ( h_Domain.SimulationStage == 2 ){
			if ( SimulationStage == 1 ){
				BufferTotalTimestep	= h_Domain.TotalTimestep ;
				h_Domain.TotalTimestep	= h_Domain.SamplingTime + 1000 ;
		
			}else if ( SimulationStage == h_Domain.SimulationStage ){
				h_Domain.TotalTimestep	= BufferTotalTimestep ;
				h_Domain.SamplingTime	= 3000 ;
				
				
				//InitializationTwoStage( &h_Domain , h_Node , h_Cell , h_Inlet , h_WallType , h_Species , h_Surface , &h_DSMC , 
				//			&h_Mapping , OutputDebug ) ;
				AdjustTimestepWeighting( &h_Domain , h_Cell , h_Species , &h_DSMC , &h_Mapping , OutputDebug ) ;
			}
		}
		

		if ( MPIMyID == 0 ){
			FileName	= "ParticlesNum-" + StageNo + ".dat" ;
			OutputParticleNum.open( FileName.c_str() , ios::out | ios::trunc ) ;
				
			FileName	= "Information-" + StageNo + ".dat" ;
			OutputInfo.open( FileName.c_str() , ios::out | ios::trunc ) ;
				
			FileName	= "Convergence-" + StageNo + ".dat" ;	
			OutputConvergence.open( FileName.c_str() , ios::out | ios::trunc ) ;
			
			OutputParticleNum << setw(25) << "TotalTimestepNo" << setw(25) << "TimestepNo" << setw(25) << "ParticleNum" << '\n' ;
			OutputConvergence << setw(15) << "TimestepNo-1" << setw(15) << "TimestepNo-2" << setw(15) << "ParticleNum" << setw(12) << "Speed" 
						<< setw(12) << "Temp" << setw(15) << "DiffParticle" << setw(12) << "DiffSpeed" << setw(12) << "DiffTemp" << '\n' ; 
		
			// Debug.
			h_Domain.DumpFile( OutputInfo ) ;
		}
		//MPI_Barrier( MPI_COMM_WORLD ) ;
		//exit(1) ;
	
	
		h_Domain.TimestepNo	= 0 ;
		
	
		for ( int i=1 ; i<=h_Domain.TotalTimestep ; i++ ){
			// Reset sampling data.
			if ( i<=h_Domain.SamplingTime ) SampleInit( &h_Domain , &h_DSMC , h_Cell) ;
			
		
			for ( int j=1 ; j<=h_Domain.OutputFrequency ; j++){
			
				for ( int k=1 ; k<=h_Domain.SamplingFrequency ; k++ ){
					h_Domain.TimestepNo++ ;
					TotalTimestepNo++ ;
					
				
					MPI_Reduce( &h_Domain.ParticleNum , &TotalParticleNum , 1 , MPI_INT , MPI_SUM , 0 , MPI_COMM_WORLD ) ;
				
				
					if ( MPIMyID == 0 ){
						OutputParticleNum << setw(25) << TotalTimestepNo << setw(25) << h_Domain.TimestepNo << setw(25) << TotalParticleNum << '\n' ;
					
						// Debug.
						cout << setw(20) << "StageNo: " << StageNo << setw(20) << "TimestepNo:" << setw(10) << h_Domain.TimestepNo << setw(20) << "ParticleNum:" << setw(15) << TotalParticleNum << '\n' ; 
						//cout << setw(40) << h_Domain.TimestepNo << setw(40) << TotalParticleNum << '\n' ;
					}
				
				
					// Sample physical computing time.
					SamplePhysicalTime( &h_Domain , h_Cell , &h_DSMC ) ;
				
				
					// Debug.
					OutputDebug << setw(12) << h_Domain.TimestepNo << ", Before Move" << endl ;
					//MPI_Barrier( MPI_COMM_WORLD ) ;
				
				
					// Particle movement.
					//Time.Time() ;
					ParticleMovement( &h_Domain , &h_Processor , h_Node , h_Cell , h_Inlet , h_WallType , h_Species , h_Surface , &h_DSMC , 
							  &h_Mapping , h_MPIParticleOut , h_MPIParticleIn , &MPIDataType , OutputDebug ) ;
					//Time.Time( &Time.Move ) ;
				
					// Debug.
					OutputDebug << setw(12) << h_Domain.TimestepNo << ", After Move" << endl ;
					//MPI_Barrier( MPI_COMM_WORLD ) ;
				
				
					// Debug.
					OutputDebug << setw(12) << h_Domain.TimestepNo << ", Before Index" << endl ;
					//MPI_Barrier( MPI_COMM_WORLD ) ;
				
				
					// Index.
					//Time.Time() ;
					Index( &h_Domain , h_Species , &h_DSMC , OutputDebug ) ;
					//Time.Time( &Time.Index ) ;
				
				
					// Debug.
					OutputDebug << setw(12) << h_Domain.TimestepNo << ", After Index" << endl ;
					//MPI_Barrier( MPI_COMM_WORLD ) ;
				
				
					// Debug.
					OutputDebug << setw(12) << h_Domain.TimestepNo << ", Before Collision" << endl ;
					//MPI_Barrier( MPI_COMM_WORLD ) ;
				
					
					// Collision.
					//Time.Time() ;
					//if (  )
					Collision( &h_Domain , h_Cell , h_Species , &h_DSMC , &h_Result, h_ChemicalTCE , OutputDebug ) ;
					//Time.Time( &Time.Collision ) ;
				
				
					// Debug.
					OutputDebug << setw(12) << h_Domain.TimestepNo << ", After Collision" << endl ;
					//MPI_Barrier( MPI_COMM_WORLD ) ;
				}	
			
				// Debug.
				OutputDebug << setw(12) << h_Domain.TimestepNo << ", Before Sample" << endl ;
				//MPI_Barrier( MPI_COMM_WORLD ) ;
				
			
				Index( &h_Domain , h_Species , &h_DSMC , OutputDebug ) ;

				// Sampling.
				//Time.Time() ;
				Sample(	&h_Domain , h_Species , &h_DSMC , h_Cell ) ;
				//Time.Time( &Time.Sample ) ;
			
			
				// Debug.
				OutputDebug << setw(12) << h_Domain.TimestepNo << ", After Sample" << endl ;
				//MPI_Barrier( MPI_COMM_WORLD ) ;
			}
		
			// Debug.
			OutputDebug << setw(12) << h_Domain.TimestepNo << ", Before Calculate Result" << endl ;
			//MPI_Barrier( MPI_COMM_WORLD ) ;
		
		
			// Calculate macroscopic properties.
			//Time.Time() ;
			CalculateResult( &h_Domain , h_Cell , h_Species , &h_DSMC , &h_Result ) ;
			//Time.Time( &Time.CalResult ) ;
			
			
			//if ( i <= h_Domain.SamplingTime )
				CheckConvergence( &h_Domain , &h_Result , &h_Convergence , i , OutputConvergence ) ;


			// Ouput simulation results.
			if ( i > (h_Domain.SamplingTime+50) ){
				if ( (i%500) == 0 || i == h_Domain.TotalTimestep ){
					OutputResult( &h_Domain , h_Cell , h_WallType , h_Surface , &h_DSMC , &h_Result , h_Species , i ) ;
				}
			}
		
		
			// Output information.
			MPI_Allgather( &h_Domain.ParticleNum , 1 , MPI_INT , ParticleNum , 1 , MPI_INT , MPI_COMM_WORLD ) ;
			MPI_Reduce( &h_Domain.EnterParticleNum , &EnterParticleNum , 1 , MPI_INT , MPI_SUM , 0 , MPI_COMM_WORLD ) ;
		
		
			if ( MPIMyID == 0 ){
				// Find the maximum and minimum number of particles across processors.
				AveParticleNum	= FindMaxMin( &MinParticleNum , &MaxParticleNum , ParticleNum , MPISize ) ;
			
				OutputInfo << "Timestep No: " << h_Domain.TimestepNo << '\n' ;
				OutputInfo << "Enter Particle Number:         " << EnterParticleNum << '\n' ;
				OutputInfo << "Average Enter Particle Number: " << double(EnterParticleNum)/double(h_Domain.SamplingNum*2) << '\n' ;
				OutputInfo << "Maximum Particle Number:       " << MaxParticleNum << '\n' ;
				OutputInfo << "Minimum Particle Number:       " << MinParticleNum << '\n' ;
				OutputInfo << "Average Particle Number:       " << AveParticleNum << '\n' ;
				OutputInfo << "Non-Uniformity (%):            " << double(MaxParticleNum - MinParticleNum)/double(AveParticleNum) * 100. << '\n' ;
			
			
				for ( int i=0 ; i<MPISize ; i++ )
					OutputInfo << "Processor No:" << setw(8) << i << ", Particle Number: " << setw(16) << ParticleNum[i] << '\n' ;
			}
		
			// Output particle number and simulation time.
			//Time.End() ;
			//DumpSimulationInformation( &h_Domain , &Time ) ;
		}
	
		
		if ( MPIMyID == 0 ){
			OutputInfo << setw(25) << "SamplingTime:" << setw(25) << h_Domain.SamplingTime << setw(25) << "TotalTimestep:" << setw(25) << h_Domain.TotalTimestep << '\n' ;
			
			
			OutputParticleNum.clear() ;
			OutputParticleNum.close() ;
		
			OutputInfo.clear() ;
			OutputInfo.close() ;
			
			OutputConvergence.clear() ;
			OutputConvergence.close() ;
		}
	}while ( SimulationStage < h_Domain.SimulationStage ) ;
	
	
	if ( Debug == 1 ){
		OutputDebug.clear() ;
		OutputDebug.close() ;
	}
	
	
	// Delete Memory.
	delete [] h_Node ;
	delete [] h_Cell ;
	for ( int i=0 ; i<h_Domain.InletFaceNum ; i++ )
		h_Inlet[i].DeleteMemory() ;
	delete [] h_Inlet ;
	delete [] h_WallType ;
	delete [] h_Species ;
	delete [] h_Surface ;
	delete [] h_ChemicalTCE ; // M.C.
	
	h_DSMC.DeleteMemory( h_Domain.LocalCellNum , h_Domain.WallFaceNum , h_Domain.InletFaceNum , h_Domain.SpeciesNum , h_Domain.SpeciesGroupNum ) ;
	h_Result.DeleteMemory( h_Domain.SpeciesNum ) ;
	h_Processor.DeleteMemory() ;
	delete [] h_MPIParticleOut ;
	delete [] h_MPIParticleIn ;
	delete [] ParticleNum ;
	// End of deleted memory.


	MPI_Barrier( MPI_COMM_WORLD ) ;
	if ( MPIMyID == 0 )
		cout << "\n===================== END OF DSMC SIMULATION =====================\n\n" ;
	
	
	MPI_Barrier( MPI_COMM_WORLD ) ;
	Time.End() ;
	
	if ( MPIMyID == 0 )
		Time.PrintFile() ;
	
	MPI_Finalize() ;

	return	0 ;
}

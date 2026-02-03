#include <mpi.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include "dsmc_class.h"
#include "dsmc_readfile.h"


using namespace std ;


//==============================================================================================================
//==============================================================================================================

// Read Input File (Input.txt)
void ReadInput( DSMC_DOMAIN *pDomain , string Filename ){
	ifstream		Input ;
	string			get_line , word , number ;
	int			wordstart , wordend ;

	Input.open( Filename.c_str() , ios::in ) ;

	if ( Input.fail() ){
		cout << "Failed to open file " << Filename << endl ;
		exit(1) ;
	}


	while ( getline(Input, get_line) ){
		if ( get_line[0] != '#' ){
			wordend = get_line.find(' ') ;
			word = get_line.substr(0,wordend) ;


			if ( word == "Dimension"){
				wordstart = get_line.find('=') ;
				number = get_line.substr(wordstart+2) ;

				pDomain->Dimension = atoi( number.c_str() ) ;
			
			}else if ( word == "Node_Number"){
				wordstart = get_line.find('=') ;
				number = get_line.substr(wordstart+2) ;

				pDomain->NodeNum = atoi( number.c_str() ) ;

			}else if ( word == "Cell_Number"){
				wordstart = get_line.find('=') ;
				number = get_line.substr(wordstart+2) ;

				pDomain->CellNum = atoi( number.c_str() ) ;
			
			}else if ( word == "Local_Cell_Number"){
				wordstart = get_line.find('=') ;
				number = get_line.substr(wordstart+2) ;

				pDomain->LocalCellNum = atoi( number.c_str() ) ;
			
			}else if ( word == "Total_Cell_Number"){
				wordstart = get_line.find('=') ;
				number = get_line.substr(wordstart+2) ;

				pDomain->TotalCellNum = atoi( number.c_str() ) ;

			}else if ( word == "Inlet_Face_Number"){
				wordstart = get_line.find('=') ;
				number = get_line.substr(wordstart+2) ;

				pDomain->InletFaceNum = atoi( number.c_str() ) ;

			}else if ( word == "Wall_Face_Number"){
				wordstart = get_line.find('=') ;
				number = get_line.substr(wordstart+2) ;

				pDomain->WallFaceNum = atoi( number.c_str() ) ;

			}else if ( word == "Wall_Type_Number"){
				wordstart = get_line.find('=') ;
				number = get_line.substr(wordstart+2) ;

				pDomain->WallTypeNum = atoi( number.c_str() ) ;

			}else if ( word == "Sampling_Frequency"){
				wordstart = get_line.find('=') ;
				number = get_line.substr(wordstart+2) ;

				pDomain->SamplingFrequency = atoi( number.c_str() ) ;

			}else if ( word == "Output_Frequency"){
				wordstart = get_line.find('=') ;
				number = get_line.substr(wordstart+2) ;

				pDomain->OutputFrequency = atoi( number.c_str() ) ;

			}else if ( word == "Sampling_Time"){
				wordstart = get_line.find('=') ;
				number = get_line.substr(wordstart+2) ;

				pDomain->SamplingTime = atoi( number.c_str() ) ;

			}else if ( word == "Total_Timestep" ){
				wordstart = get_line.find('=') ;
				number = get_line.substr(wordstart+2) ;

				pDomain->TotalTimestep = atoi( number.c_str() ) ;

			}else if ( word == "Timestep_Ratio" ){
				wordstart = get_line.find('=') ;
				number = get_line.substr(wordstart+2) ;

				pDomain->TimestepRatio = atof( number.c_str() ) ;
				
			}else if ( word == "Timestep" ){
				wordstart = get_line.find('=') ;
				number = get_line.substr(wordstart+2) ;

				pDomain->Timestep = atof( number.c_str() ) ;

			}else if ( word == "Weighting_Ratio" ){
				wordstart = get_line.find('=') ;
				number = get_line.substr(wordstart+2) ;

				pDomain->WeightingRatio = atof( number.c_str() ) ;
				
			}else if ( word == "Weighting" ){
				wordstart = get_line.find('=') ;
				number = get_line.substr(wordstart+2) ;

				pDomain->ParticleWeighting = atof( number.c_str() ) ;

			}else if ( word == "Maximum_Particle_Number" ){
				wordstart = get_line.find('=') ;
				number = get_line.substr(wordstart+2) ;

				pDomain->MaxParticleNum = atoi( number.c_str() ) ;
				
			}else if ( word == "Transfer_Particle_Number" ){
				wordstart = get_line.find('=') ;
				number = get_line.substr(wordstart+2) ;

				pDomain->TransferParticleNum = atoi( number.c_str() ) ;

			}else if ( word == "X_Velocity" ){
				wordstart = get_line.find('=') ;
				number = get_line.substr(wordstart+2) ;

				pDomain->XVel = atof( number.c_str() ) ;

			}else if ( word == "Y_Velocity" ){
				wordstart = get_line.find('=') ;
				number = get_line.substr(wordstart+2) ;

				pDomain->YVel = atof( number.c_str() ) ;

			}else if ( word == "Z_Velocity" ){
				wordstart = get_line.find('=') ;
				number = get_line.substr(wordstart+2) ;

				pDomain->ZVel = atof( number.c_str() ) ;

			}else if ( word == "Number_Density" ){
				wordstart = get_line.find('=') ;
				number = get_line.substr(wordstart+2) ;

				pDomain->NumDen = atof( number.c_str() ) ;

			}else if ( word == "Temperature" ){
				wordstart = get_line.find('=') ;
				number = get_line.substr(wordstart+2) ;

				pDomain->Temp = atof( number.c_str() ) ;

			}else if ( word == "Species_Number" ){
				wordstart = get_line.find('=') ;
				number = get_line.substr(wordstart+2) ;

				pDomain->SpeciesNum = atoi( number.c_str() ) ;

			}else if ( word == "Species_Group_Number" ){
				wordstart = get_line.find('=') ;
				number = get_line.substr(wordstart+2) ;

				pDomain->SpeciesGroupNum = atoi( number.c_str() ) ;

			}else if ( word == "Variable_Timestep_Scheme" ){
				wordstart = get_line.find('=') ;
				number = get_line.substr(wordstart+2) ;

				pDomain->VariableTimestepScheme = atoi( number.c_str() ) ;

			}else if ( word == "Subcell_Model" ){
				wordstart = get_line.find('=') ;
				number = get_line.substr(wordstart+2) ;

				pDomain->SubcellModel = atoi( number.c_str() ) ;
			
			}else if ( word == "Dynamic_Domain_Decompisition" ){
				wordstart = get_line.find('=') ;
				number = get_line.substr(wordstart+2) ;

				pDomain->DynamicDomainDecompisition = atoi( number.c_str() ) ;

			}else if ( word == "Vibrational_Model" ){
				wordstart = get_line.find('=') ;
				number = get_line.substr(wordstart+2) ;

				pDomain->VibrationalModel = atoi( number.c_str() ) ;
				
			}else if ( word == "Simulation_Stage" ){
				wordstart = get_line.find('=') ;
				number = get_line.substr(wordstart+2) ;

				pDomain->SimulationStage = atoi( number.c_str() ) ;
			
			}else if ( word == "Check_Convergence" ){
				wordstart = get_line.find('=') ;
				number = get_line.substr(wordstart+2) ;

				pDomain->Convergence = atoi( number.c_str() ) ;
			
			}else if ( word == "Particle_Number_Check" ){
				wordstart = get_line.find('=') ;
				number = get_line.substr(wordstart+2) ;

				pDomain->ParticleNumCheck = atof( number.c_str() ) ;
			
			}else if ( word == "Speed_Check" ){
				wordstart = get_line.find('=') ;
				number = get_line.substr(wordstart+2) ;

				pDomain->SpeedCheck = atof( number.c_str() ) ;
			
			}else if ( word == "Temperature_Check" ){
				wordstart = get_line.find('=') ;
				number = get_line.substr(wordstart+2) ;

				pDomain->TempCheck = atof( number.c_str() ) ;
			
			}else if ( word == "Axisymmetric_Adjusting_Factor" ){
				wordstart = get_line.find('=') ;
				number = get_line.substr(wordstart+2) ;

				pDomain->AxisAdjustFactor = atof( number.c_str() ) ;
			}
		}
	}
	// To close the file.
	Input.clear() ;
	Input.close() ;
	
	
	if ( pDomain->SubcellModel == 2 ){
		if ( pDomain->VariableTimestepScheme != 2 || pDomain->VariableTimestepScheme != 3 )
			pDomain->SubcellModel	= 1 ;
	}
	
	
	if ( pDomain->Convergence == 1 ){
		pDomain->SamplingTime	= 50000 ;
		pDomain->TotalTimestep	= 52000 ;
	}
}

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
				ofstream 		&OutputDebug ){

	int		MPISize , MPIMyID ;
	
	MPI_Comm_size( MPI_COMM_WORLD , &MPISize ) ;
	MPI_Comm_rank( MPI_COMM_WORLD , &MPIMyID ) ;
					
	
	// Synchronization for all processors.
	MPI_Barrier( MPI_COMM_WORLD ) ;
	
	
	// Read processor information from the file (Partition.inp)
	ReadPartitionMPI( h_pDomain , h_pProcessor , "Partition.inp" ) ;
	
	
	// Read mesh information from the file (Mesh.inp).
	ReadMeshMPI( h_Node , h_Cell , h_pDomain , h_pProcessor , h_pMapping , "Mesh.inp" , OutputDebug ) ;
	
	
	// Debug.
	/*for ( int i=0 ; i<MPISize ; i++ ){
		if ( MPIMyID == i ){
			cout << "Processor - " << MPIMyID << '\n' ;
			for ( int j=0 ; j<h_pDomain->CellNum ; j++ ){
				cout << setw(12) << h_Cell[j].Id << setw(12) << h_Cell[j].LocalId ;
				//for ( int k=0 ; k<h_pMapping->NodeNum[h_Cell[j].Type] ; k++ )
				//	cout << setw(12) << h_Cell[j].Node[k] ;
				for ( int k=0 ; k<h_pMapping->SurfaceNum[h_Cell[j].Type] ; k++ )
					cout << setw(12) << h_Cell[j].Neighbor[k] ;
				cout << '\n' ;
			}
		}
		MPI_Barrier( MPI_COMM_WORLD ) ;
	}*/

	
	// Read inlet/outlet face information from the file (Inlet.inp).
	ReadInlet( h_Inlet , h_pDomain , "Inlet.inp" ) ;
	
	
	// Debug.
	/*for ( int i=0 ; i<h_pDomain->InletFaceNum ; i++ ) h_Inlet[i].Dump(i) ;
	MPI_Barrier( MPI_COMM_WORLD ) ;
	exit(1) ;*/
	
	
	// Read solid face information from the file (WallType.txt).
	ReadWallType( h_WallType , h_pDomain , "WallType.txt" ) ;
	
	// Debug.
	//for ( int i=0 ; i<h_pDomain.WallTypeNum ; i++ ) h_WallType[i].Dump(i) ;
	
	
	// Read species information from the database file (Species.txt).
	ReadSpecies( h_Species , h_pDomain , "Species.txt" , "Input.txt" ) ;

	// Debug.
	//for ( int i=0 ; i<h_pDomain.SpeciesNum ; i++ ) h_Species[i].Dump(i) ;
	
	
	// Read cell information, including local mean free path, temp, speed...etc.
	//if ( h_pDomain->VariableTimestepScheme == ( 2 || 3 || 5 || 6 ) ){
	//	ReadCellInformation( h_Cell , h_pDomain ,  h_pProcessor , "Cell-Information.inp" ) ;
	//}
}

//==============================================================================================================

void ReadPartitionMPI( DSMC_DOMAIN *h_pDomain , DSMC_PROCESSOR *h_pProcessor , string Filename ){
	int			MPISize , MPIMyID , ProcessorNo ;
	ifstream		Input ;
	string			Buffer ;

	// Obtain number and ID of processors.
	MPI_Comm_size( MPI_COMM_WORLD , &MPISize ) ;
	MPI_Comm_rank( MPI_COMM_WORLD , &MPIMyID ) ;


	Input.open( Filename.c_str() , ios::in ) ;

	if ( Input.fail() ){
		cout << "Failed to open file " << Filename << endl ;
		exit(1) ;
	}

	
	for ( int i=0 ; i<MPISize ; i++ ){
		if ( MPIMyID == i ){
			Input >> ProcessorNo >> h_pDomain->CellNum >> h_pProcessor->NeighborNum ;
			
			if ( h_pDomain->CellNum > h_pDomain->LocalCellNum ){
				cout << "The setup of \"Local_Cell_Number\" is very small. Actually the cell number on this processor is " << h_pDomain->CellNum << '.' << endl ;
				exit(1) ;
			}else if ( h_pProcessor->NeighborNum > h_pProcessor->MaxNeighborNum ){
				cout << "Neighbor number of processor is more than maximum number. Processor ID: " << MPIMyID << ", Neighbor Number: " << h_pProcessor->NeighborNum << endl ;
				exit(1) ;
			}else{
				for ( int j=0 ; j<h_pProcessor->NeighborNum ; j++ ){
					Input >> h_pProcessor->Neighbor[j] ;
				}
			}
			break ;
		}else{
			getline( Input , Buffer ) ;
		}	
	}
	
	
	// Debug.
	/*for ( int i=0 ; i<MPISize ; i++ ){
		if ( MPIMyID == i ){
			h_pProcessor->Dump( MPIMyID , h_pDomain->CellNum ) ;
		}
		MPI_Barrier( MPI_COMM_WORLD ) ;
	}
	MPI_Barrier( MPI_COMM_WORLD ) ;
	getchar() ;*/
	
	
	
	// To close the file.
	Input.clear() ;
	Input.close() ;
}

//==============================================================================================================

// Read mesh file (Mesh.inp)
void ReadMeshMPI( DSMC_NODE *h_Node , DSMC_CELL *h_Cell , DSMC_DOMAIN *h_pDomain , DSMC_PROCESSOR *h_pProcessor , CELLMAPPING *h_pMapping , string Filename , ofstream &OutputDebug ){
	int			MPISize , MPIMyID , NodeNum , CellNum , TotalCellNum , CellNo , Type ;
	int			*GlobalCellNum , *SumGlobalCellNum ;
	ifstream		Input ;
	string			buffer ;
	
	MPI_Comm_size( MPI_COMM_WORLD , &MPISize ) ;
	MPI_Comm_rank( MPI_COMM_WORLD , &MPIMyID ) ;
	

	NodeNum		= 0 ;
	CellNum		= h_pDomain->CellNum ;
	TotalCellNum	= 0 ;
	
	GlobalCellNum		= new int[MPISize] ;
	SumGlobalCellNum	= new int[MPISize] ;
	
	
	for ( int i=0 ; i<MPISize ; i++ ){
		GlobalCellNum[i]	= 0 ;
		SumGlobalCellNum[i]	= 0 ;
	}
	
	
	//MPI_Allgather( &CellNum , 1 , MPI_INT , &GlobalCellNum[0] , 1 , MPI_INT , MPI_COMM_WORLD ) ;
	MPI_Gather( &CellNum , 1 , MPI_INT , GlobalCellNum , 1 , MPI_INT , 0 , MPI_COMM_WORLD ) ;
	MPI_Bcast( GlobalCellNum , MPISize , MPI_INT , 0 , MPI_COMM_WORLD ) ;
	
	
	SumGlobalCellNum[0]	= 0 ;
	for ( int i=1 ; i<MPISize ; i++ )
		SumGlobalCellNum[i]	= SumGlobalCellNum[i-1] + GlobalCellNum[i-1] ;
		
	// Debug.
	//for ( int i=0 ; i<MPISize ; i++ )
	//	OutputDebug << setw(10) << i << setw(10) << SumGlobalCellNum[i] << setw(10) << GlobalCellNum[i] << endl ;
		
	// Debug.
	//OutputDebug << MPISize << ", TotalCellNum: " << (SumGlobalCellNum[MPISize-1] + GlobalCellNum[MPISize-1]) << endl ;
	//MPI_Barrier( MPI_COMM_WORLD ) ;
	//exit(1) ;
	
	// Debug.
	/*for ( int i=0 ; i<MPISize ; i++ ){
		if ( MPIMyID == i ){
			for ( int j=0 ; j<MPISize ; j++ )	
				cout << setw(20) << SumGlobalCellNum[j] ;
				
			cout << '\n' ;
		}	
		MPI_Barrier( MPI_COMM_WORLD ) ;
	}*/
	
	
	for ( int i=0 ; i<MPISize ; i++ ){
		if ( MPIMyID == i ){
			Input.open( Filename.c_str() , ios::in ) ;

			if ( Input.fail() ){
				cout << "Failed to open file " << Filename << endl ;
				exit(1) ;
			}

			while ( getline( Input , buffer ) ){
				if ( FindString(buffer , "Nodes") ){
					Input >> NodeNum ;
			
					// Check node number is right.
					if ( NodeNum != h_pDomain->NodeNum ){
						cout << "Node number is error. Domain->NodeNum: " << h_pDomain->NodeNum << ", NodeNum: " << NodeNum << endl ;
						exit(1) ;
					}
			
					for ( int i=0 ; i<NodeNum ; i++ )
						Input >> h_Node[i].Id >> h_Node[i].XCoord >>  h_Node[i].YCoord >> h_Node[i].ZCoord ;
				}


				if ( FindString(buffer , "Cells") ){
					Input >> TotalCellNum ;
		
					// Check node number is right.
					if ( TotalCellNum != h_pDomain->TotalCellNum ){
						cout << "Cell number is error. Domain->TotalCellNum: " << h_pDomain->TotalCellNum << ", TotalCellNum: " << TotalCellNum << endl ;
						exit(1) ;
					}
					getline( Input , buffer ) ;
		
					for ( int i=0 ; i<TotalCellNum ; i++ ){
						if ( i >= SumGlobalCellNum[MPIMyID] ){
							for ( int j=0 ; j<CellNum ; j++ ){
								Input >> h_Cell[j].Id >> h_Cell[j].Type ;
				
								for ( int k=0 ; k<h_pMapping->NodeNum[h_Cell[j].Type] ; k++)
									Input >> h_Cell[j].Node[k] ;

								h_Cell[j].LocalId			= j ;
								h_pProcessor->LocalCellNo[h_Cell[j].Id]	= j ;
								h_pProcessor->CellProcessorNo[h_Cell[j].Id] = MPIMyID ;
							}

							break ;
						}else{
							getline( Input , buffer ) ;
						}
					}
				}


				if ( FindString(buffer , "Neighbors") ){
					for ( int i=0 ; i<TotalCellNum ; i++ ){
						if ( i >= SumGlobalCellNum[MPIMyID] ){
							for ( int j=0 ; j<CellNum ; j++ ){
								Input >> CellNo ;
								Type	= h_Cell[j].Type ;
								
								// Debug.
								//cout << "CellNo: " << CellNo << ", ID: " << h_Cell[j].Id << ", Type: " << h_Cell[j].Type << '\n' ;
				
								for ( int k=0 ; k<h_pMapping->SurfaceNum[Type] ; k++ )
									Input >> h_Cell[j].Neighbor[k] ;
							}
							
							break ;
						}else{
							getline( Input , buffer ) ;
						}
					}
				}
			}
			
			// To close the file.
			Input.clear() ;
			Input.close() ;
		}
		
		MPI_Barrier( MPI_COMM_WORLD ) ;
	} // End of each Processor
	
	
	// Debug.
	/*for ( int i=0 ; i<MPISize ; i++ ){
		if ( MPIMyID == i ){
			cout << "Processor: " << MPIMyID << '\n' ;
			cout << "==============================================\n" ;
			for ( int j=0 ; j<h_pDomain->TotalCellNum ; j++ ){
				cout << setw(15) << h_pProcessor->LocalCellNo[j] << setw(15) << h_pProcessor->CellProcessorNo[j] << '\n' ;
			}
		}
		MPI_Barrier( MPI_COMM_WORLD ) ;
	}*/
	
	
	for ( int i=0 ; i<MPISize ; i++ ){
		MPI_Bcast( &h_pProcessor->LocalCellNo[SumGlobalCellNum[i]] , GlobalCellNum[i] , MPI_INT , i , MPI_COMM_WORLD ) ;
		MPI_Bcast( &h_pProcessor->CellProcessorNo[SumGlobalCellNum[i]] , GlobalCellNum[i] , MPI_INT , i , MPI_COMM_WORLD ) ;
	}
	
	
	// Debug.
	/*for ( int i=0 ; i<MPISize ; i++ ){
		if ( MPIMyID == i ){
			cout << "Processor: " << MPIMyID << '\n' ;
			cout << "==============================================\n" ;
			for ( int j=0 ; j<h_pDomain->TotalCellNum ; j++ ){
				cout << setw(15) << h_pProcessor->LocalCellNo[j] << setw(15) << h_pProcessor->CellProcessorNo[j] << '\n' ;
			}
		}
		MPI_Barrier( MPI_COMM_WORLD ) ;
	}*/
	
	
	delete [] GlobalCellNum ;
	delete [] SumGlobalCellNum ;
}


// Read mesh file (Mesh.inp)
void ReadMesh( DSMC_NODE *h_Node , DSMC_CELL *h_Cell , DSMC_DOMAIN *h_pDomain , CELLMAPPING *h_pMapping , string Filename ){
	int			NodeNum , CellNum , CellNo , Type ;
	ifstream		Input ;
	string			buffer ;

	NodeNum	= 0 ;
	CellNum	= 0 ;


	Input.open( Filename.c_str() , ios::in ) ;

	if ( !Input ){
		cout << "Failed to open file " << Filename << endl ;
		exit(1) ;
	}


	while ( getline( Input , buffer ) ){
		if ( FindString(buffer , "Nodes") ){
			Input >> NodeNum ;
			
			// Check node number is right.
			if ( NodeNum != h_pDomain->NodeNum ){
				cout << "Node number is error. Domain->NodeNum: " << h_pDomain->NodeNum << ", NodeNum: " << NodeNum << endl ;
				exit(1) ;
			}
			
			for ( int i=0 ; i<NodeNum ; i++ )
				Input >> h_Node[i].Id >> h_Node[i].XCoord >>  h_Node[i].YCoord >> h_Node[i].ZCoord ;
		}


		if ( FindString(buffer , "Cells") ){
			Input >> CellNum ;
		
			// Check node number is right.
			if ( CellNum != h_pDomain->CellNum ){
				cout << "Cell number is error. Domain->CellNum: " << h_pDomain->CellNum << ", CellNum: " << CellNum << endl ;
				exit(1) ;
			}
		
			for ( int i=0 ; i<CellNum ; i++ ){
				Input >> h_Cell[i].Id >> h_Cell[i].Type ;
				
				for ( int j=0 ; j<h_pMapping->NodeNum[h_Cell[i].Type] ; j++)
					Input >> h_Cell[i].Node[j] ;
			}
		}


		if ( FindString(buffer , "Neighbors") ){
			for ( int i=0 ; i<CellNum ; i++ ){
				Input >> CellNo ;
				Type	= h_Cell[i].Type ;
				
				for ( int j=0 ; j<h_pMapping->SurfaceNum[Type] ; j++ )
					Input >> h_Cell[i].Neighbor[j] ;
			}
		}
	}
	// To close the file.
	Input.clear() ;
	Input.close() ;
}

//==============================================================================================================

void ReadInlet( DSMC_INLET *h_Inlet , DSMC_DOMAIN *h_pDomain , string Filename ){
	int			FaceNum , Type ;
	ifstream		Input ;
	string			buffer ;

	FaceNum	= 0 ;


	Input.open( Filename.c_str() , ios::in ) ;

	if ( Input.fail() ){
		cout << "Failed to open file " << Filename << endl ;
		exit(1) ;
	}

	
	Input >> FaceNum >> h_pDomain->InletSpecifiedNumDenNum ;
	if ( FaceNum != h_pDomain->InletFaceNum ){
		cout << "Inlet face number is error. Domain->InletFaceNum: " << h_pDomain->InletFaceNum << ", FaceNum: " << FaceNum << endl ;
		exit(1) ;
	}else if ( h_pDomain->SpeciesNum < h_pDomain->InletSpecifiedNumDenNum ){
		cout << "Species number is error. Domain->SpeciesNum: " << h_pDomain->SpeciesNum << ", h_Domain->InletSpecifiedNumDenNum: " << h_pDomain->InletSpecifiedNumDenNum << '\n' ;
		exit(1) ;
	}
	
	
	for ( int i=0 ; i<FaceNum ; i++ ){
		h_Inlet[i].Id	= i ;
		
		Input >> h_Inlet[i].XVel >> h_Inlet[i].YVel >> h_Inlet[i].ZVel ;
		
		for ( int j=0 ; j<h_pDomain->InletSpecifiedNumDenNum ; j++ )
			Input >> h_Inlet[i].NumDen[j] ;
			
		Input >> h_Inlet[i].Temp >> h_Inlet[i].CosineLawCoef >> h_Inlet[i].NodeNum ;
		
		for ( int j=0 ; j<h_Inlet[i].NodeNum ; j++ )
			Input >> h_Inlet[i].Node[j] ;
	}
	
	// To close the file.
	Input.clear() ;
	Input.close() ;
}

//==============================================================================================================

void ReadWallType( DSMC_WALLTYPE *h_WallType , DSMC_DOMAIN *h_pDomain , string Filename ){
	ifstream		Input ;
	string			buffer ;

	Input.open( Filename.c_str() , ios::in ) ;

	if ( Input.fail() ){
		cout << "Failed to open file " << Filename << endl ;
		exit(1) ;
	}

	getline( Input , buffer ) ;
	getline( Input , buffer ) ;
	for ( int i=0 ; i<h_pDomain->WallTypeNum ; i++ ){
		h_WallType[i].Id = i ;
		Input >> h_WallType[i].WallNo >> h_WallType[i].Type >> h_WallType[i].Temp >> h_WallType[i].XVel >> h_WallType[i].YVel >> h_WallType[i].ZVel >> h_WallType[i].DiffRatio >> h_WallType[i].StickingCoef
					>> h_WallType[i].AlphaN >> h_WallType[i].SigmaT >> h_WallType[i].EintCoef ;
	}
	// To close the file.
	Input.clear() ;
	Input.close() ;
}

//==============================================================================================================

void ReadSpecies( DSMC_SPECIES *h_Species , DSMC_DOMAIN *h_pDomain , string Filename , string FileSpeciesName ){
	int			SpeciesNum , wordstart , wordend , SpeciesNoRot , SpeciesNoVibC1 , SpeciesNoVibC2 ;
	ifstream		InputSpecies ;
	string			get_line , word , number , *SpeciesName , Name ;

	SpeciesNum	= h_pDomain->SpeciesNum ;
	SpeciesName	= new string[SpeciesNum] ;

	// Read species's name.
	GetSpeciesName( SpeciesName , h_pDomain , FileSpeciesName ) ;
	
	// Debug.
	/*for ( int i=0 ; i<SpeciesNum ; i++ ){
		cout << "Species Name[" << i << "] = " << SpeciesName[i] << '\n' ;
		cout << "        , size: " << SpeciesName[i].size() << '\n' ;	
	}
	getchar() ;*/


	InputSpecies.open( Filename.c_str() , ios::in ) ;


	if ( InputSpecies.fail() ){
		cout << "Failed to open file " << Filename << endl ;
		exit(1) ;
	}

	// Read species information from database.
	while ( getline(InputSpecies , get_line) ){
		if ( get_line[0] != '#' ){
			wordend = get_line.find(' ') ;
			word = get_line.substr(0,wordend) ;

			if ( word == "Species_Name" ){
				wordstart	= get_line.find('=') ;
				Name		= get_line.substr(wordstart+2) ;

				for ( int i=0 ; i<SpeciesNum ; i++ ){
					if ( FindString( Name , SpeciesName[i] ) ){
						h_Species[i].Id	= i ;
						
						SpeciesNoRot	= 0 ;
						SpeciesNoVibC1	= 0 ;
						SpeciesNoVibC2	= 0 ;
					
						while ( getline(InputSpecies , get_line) ){
							wordend = get_line.find(' ') ;
							word = get_line.substr(0,wordend) ;
						
							if ( word == "Fraction" ){
								wordstart = get_line.find('=') ;
								number = get_line.substr(wordstart+2) ;

								h_Species[i].Fraction = atof( number.c_str() ) ;
							
							}else if ( word == "Group_No" ){
								wordstart = get_line.find('=') ;
								number = get_line.substr(wordstart+2) ;

								h_Species[i].GroupNo = atoi( number.c_str() ) ;
			
							}else if ( word == "Diameter" ){
								wordstart = get_line.find('=') ;
								number = get_line.substr(wordstart+2) ;

								h_Species[i].Diameter = atof( number.c_str() ) ;
		
							}else if ( word == "Reference_Temperature" ){
								wordstart = get_line.find('=') ;
								number = get_line.substr(wordstart+2) ;

								h_Species[i].RefTemp = atof( number.c_str() ) ;
								
							}else if ( word == "Viscosity_Index" ){
								wordstart = get_line.find('=') ;
								number = get_line.substr(wordstart+2) ;

								h_Species[i].VisTempIndex = atof( number.c_str() ) ;
								
							}else if ( word == "VSS_Parameter" ){
								wordstart = get_line.find('=') ;
								number = get_line.substr(wordstart+2) ;

								h_Species[i].VSSParameter = atof( number.c_str() ) ;
								
							}else if ( word == "Mass" ){
								wordstart = get_line.find('=') ;
								number = get_line.substr(wordstart+2) ;

								h_Species[i].Mass = atof( number.c_str() ) ;
								
							}else if ( word == "Rotational_Degrees_of_Freedom" ){
								wordstart = get_line.find('=') ;
								number = get_line.substr(wordstart+2) ;

								h_Species[i].RotDOF = atof( number.c_str() ) ;
								
							}else if ( word == "Rotational_Model" ){
								wordstart = get_line.find('=') ;
								number = get_line.substr(wordstart+2) ;

								h_Species[i].RotModel = atoi( number.c_str() ) ;
								
							}else if ( word == "Rotational_Relaxation_Number" ){
								wordstart = get_line.find('=') ;
								number = get_line.substr(wordstart+2) ;

								h_Species[i].RotRelaxationNum[SpeciesNoRot] = atof( number.c_str() ) ;
								
								SpeciesNoRot++ ;
								
							}else if ( word == "Vibrational_Mode" ){
								wordstart = get_line.find('=') ;
								number = get_line.substr(wordstart+2) ;

								h_Species[i].VibMode = atoi( number.c_str() ) ;

							}else if ( word == "Vibrational_Model" ){
								wordstart = get_line.find('=') ;
								number = get_line.substr(wordstart+2) ;

								h_Species[i].VibModel = atoi( number.c_str() ) ;
								
							}else if ( word == "Vibrational_Degrees_of_Freedom" ){
								wordstart = get_line.find('=') ;
								number = get_line.substr(wordstart+2) ;

								h_Species[i].VibDOF = atof( number.c_str() ) ;
								
							}else if ( word == "Vibrational_Constant_1" ){
								wordstart = get_line.find('=') ;
								number = get_line.substr(wordstart+2) ;

								h_Species[i].VibConstant1[SpeciesNoVibC1] = atof( number.c_str() ) ;
								
								SpeciesNoVibC1++ ;
								
							}else if ( word == "Vibrational_Constant_2" ){
								wordstart = get_line.find('=') ;
								number = get_line.substr(wordstart+2) ;

								h_Species[i].VibConstant2[SpeciesNoVibC2] = atof( number.c_str() ) ;
								
								SpeciesNoVibC2++ ;
								
							}else if ( word == "Vibrational_Characteristic_Temperature" ){
								wordstart = get_line.find('=') ;
								number = get_line.substr(wordstart+2) ;

								h_Species[i].VibTemp = atof( number.c_str() ) ;
								
							}else if ( word == "Vibrational_Dissociation_Temperature" ){
								wordstart = get_line.find('=') ;
								number = get_line.substr(wordstart+2) ;

								h_Species[i].DisTemp = atof( number.c_str() ) ;
								
								break ;
							}
						}
					}	
				}
			}
		}
	}
	// To close the file.
	InputSpecies.clear() ;
	InputSpecies.close() ;
	

	delete [] SpeciesName ;
}

//==============================================================================================================

void ReadCellInformation( DSMC_CELL *h_Cell , DSMC_DOMAIN *h_pDomain , DSMC_PROCESSOR *h_pProcessor , string Filename ){
	int			MPISize , MPIMyID ;
	int			CellNo , ProcessorNo , TotalCellNum ;
	double			a ;
	ifstream		Input ;
	string			buffer ;


	MPI_Comm_size( MPI_COMM_WORLD , &MPISize ) ;
	MPI_Comm_rank( MPI_COMM_WORLD , &MPIMyID ) ;


	TotalCellNum	= h_pDomain->TotalCellNum ;
	

	Input.open( Filename.c_str() , ios::in ) ;

	if ( !Input ){
		cout << "Failed to open file " << Filename << endl ;
		exit(1) ;
	}


	getline( Input , buffer ) ;
	for ( int i=0 ; i<TotalCellNum ; i++ ){
		Input >> CellNo ;
		
		// Global cell number to local cell number.
		ProcessorNo	= h_pProcessor->CellProcessorNo[CellNo] ;
		if ( ProcessorNo == MPIMyID ){
			CellNo	= h_pProcessor->LocalCellNo[CellNo] ;
			Input >> h_Cell[CellNo].InitMeanFreePath >> a >> a >> a >> h_Cell[CellNo].InitTimestep >> h_Cell[CellNo].InitWeighting 
				>> h_Cell[CellNo].AveParticleNum >> h_Cell[CellNo].MeanCollSpacingMeanFreePath >> a >> a >> a ;
		}else{
			getline( Input , buffer ) ;
		}
	}

	
	// Debug.
	/*for ( int i=0 ; i<MPISize ; i++ ){
		if ( MPIMyID == i ){
			cout << "CellNum: " << h_pDomain->CellNum << '\n' ;
			for ( int j=0 ; j<h_pDomain->CellNum ; j++ ){
				cout << setw(20) << h_Cell[j].Id << setw(20) << h_Cell[j].InitMeanFreePath << setw(20) << h_Cell[j].InitTemp << setw(20) << h_Cell[j].InitSpeed 
					<< setw(20) << h_Cell[j].InitTimestep << setw(20) << h_Cell[j].InitWeighting << '\n' ;
			}
		}
		
		MPI_Barrier( MPI_COMM_WORLD ) ;
	}*/
	
	
	// To close the file.
	Input.clear() ;
	Input.close() ;
}

//==============================================================================================================

bool FindString( string Word1 , string Word2 ){
	bool	a = true ;
	//for ( int i=0 ; i<Word2.size() ; i++ ){
	for ( int i=0 ; i<(Word2.size()-1) ; i++ ){
		if ( Word1[i] != Word2[i] ){
			a = false ;
			break ;
		}
	}

	return	a ;
}

//==============================================================================================================

void GetSpeciesName( string *SpeciesName , DSMC_DOMAIN *h_pDomain , string Filename ){
	int			SpeciesNum , SpeciesNo , wordstart , wordend ;
	ifstream		Input ;
	string			get_line ;
	string			word ;

	SpeciesNum	= 0 ;
	SpeciesNo	= 0 ;

	Input.open( Filename.c_str() , ios::in ) ;

	if ( Input.fail() ){
		cout << "Failed to open file " << Filename << endl ;
		exit(1) ;
	}

	while ( getline(Input, get_line) ){
		if ( get_line[0] != '#' ){
			wordend = get_line.find(' ') ;
			word = get_line.substr(0,wordend) ;

			if ( word == "Species_Name" ) SpeciesNum++ ;
		}
	}
	// To close the file.
	Input.clear() ;
	Input.close() ;
	
	if ( SpeciesNum != h_pDomain->SpeciesNum ){
		cout << "Species number is error. Species_Number: " << h_pDomain->SpeciesNum << ", SpeciesNum: " << SpeciesNum << endl ;
		exit(1) ;
	}
	
	
	Input.open( Filename.c_str() , ios::in ) ;

	if ( Input.fail() ){
		cout << "Failed to open file " << Filename << endl ;
		exit(1) ;
	}

	while ( getline(Input, get_line) ){
		if ( get_line[0] != '#' ){
			wordend = get_line.find(' ') ;
			word = get_line.substr(0,wordend) ;

			if ( word == "Species_Name" ){
				wordstart		= get_line.find('=') ;
				SpeciesName[SpeciesNo]	= get_line.substr(wordstart+2) ;
				SpeciesNo++ ;
			}
		}
	}
	// To close the file.
	Input.clear() ;
	Input.close() ;
}

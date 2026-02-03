#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <cstdlib>
#include "dsmc_postprocessor.h"

using namespace std ;

//==============================================================================================================
//==============================================================================================================

void OutptuResultTec( DSMC_DOMAIN *pDomain , DSMC_NODE *Node , DSMC_CELL *Cell , DSMC_RESULT *pResult , CELLMAPPING *pMapping ){
	if ( pDomain->Dimension == 2 ){
		OutptuResultTec2D( pDomain , Node , Cell , pResult , pMapping ) ;
	}else if ( pDomain->Dimension == 3 ){
		OutptuResultTec3D( pDomain , Node , Cell , pResult , pMapping ) ;
	}
}

//==============================================================================================================

void OutptuResultSpeciesTec( DSMC_DOMAIN *pDomain , DSMC_NODE *Node , DSMC_CELL *Cell , DSMC_RESULT *pResult , CELLMAPPING *pMapping , string SpeciesName ){
	if ( pDomain->Dimension == 2 ){
		OutptuResultSpeciesTec2D( pDomain , Node , Cell , pResult , pMapping , SpeciesName ) ;
	}else if ( pDomain->Dimension == 3 ){
		OutptuResultSpeciesTec3D( pDomain , Node , Cell , pResult , pMapping , SpeciesName ) ;
	}
}

//==============================================================================================================

void OutptuResultTec2D( DSMC_DOMAIN *pDomain , DSMC_NODE *Node , DSMC_CELL *Cell , DSMC_RESULT *pResult , CELLMAPPING *pMapping ){
	ofstream	Output ;
	int		NodeNum , CellNum , OutNum ;
	
	NodeNum	= pDomain->NodeNum ;
	CellNum	= pDomain->CellNum ;
	OutNum	= 0 ;

	Output.open( "Result-tec.dat" , ios::out | ios::trunc ) ;

	if ( !Output.fail() ){
		Output << "TITLE = \"Finite volume dataset\"\n" ;
		Output << "VARIABLES = \"X\", \"Y\", \"Density\", \"NumDensity\", \"U-Vel\", \"V-Vel\", \"W-Vel\", \"Temp\", \"TransTemp\", \"RotTemp\", \"VibTemp\", \"#P\", \"mcs/mfp\", \"ProcessorNo\"\n" ;
		Output << "ZONE N=" << NodeNum << ", E=" << CellNum << ", DATAPACKING=BLOCK, ZONETYPE=FEQUADRILATERAL\n" ;
		Output << "VARLOCATION=([1-2]=NODAL, [3-14]=CELLCENTERED)\n" ;

		
		for ( int i=0 ; i<NodeNum ; i++ ){
			OutNum++ ;
			Output << setw(18) << Node[i].XCoord ;
			if ( OutNum%6 == 0 ) Output << '\n' ;
		}
		
		for ( int i=0 ; i<NodeNum ; i++ ){
			OutNum++ ;
			Output << setw(18) << Node[i].YCoord ;
			if ( OutNum%6 == 0 ) Output << '\n' ;
		}

		for ( int i=0 ; i<CellNum ; i++ ){
			OutNum++ ;
			Output << setw(18) << pResult->Density[i] ;
			if ( OutNum%6 == 0 ) Output << '\n' ;
		}

		for ( int i=0 ; i<CellNum ; i++ ){
			OutNum++ ;
			Output << setw(18) << pResult->NumDensity[i] ;
			if ( OutNum%6 == 0 ) Output << '\n' ;
		}

		for ( int i=0 ; i<CellNum ; i++ ){
			OutNum++ ;
			Output << setw(18) << pResult->XVel[i] ;
			if ( OutNum%6 == 0 ) Output << '\n' ;
		}
		
		for ( int i=0 ; i<CellNum ; i++ ){
			OutNum++ ;
			Output << setw(18) << pResult->YVel[i] ;
			if ( OutNum%6 == 0 ) Output << '\n' ;
		}
		
		for ( int i=0 ; i<CellNum ; i++ ){
			OutNum++ ;
			Output << setw(18) << pResult->ZVel[i] ;
			if ( OutNum%6 == 0 ) Output << '\n' ;
		}

		for ( int i=0 ; i<CellNum ; i++ ){
			OutNum++ ;
			Output << setw(18) << pResult->Temp[i] ;
			if ( OutNum%6 == 0 ) Output << '\n' ;
		}

		for ( int i=0 ; i<CellNum ; i++ ){
			OutNum++ ;
			Output << setw(18) << pResult->TransTemp[i] ;
			if ( OutNum%6 == 0 ) Output << '\n' ;
		}

		for ( int i=0 ; i<CellNum ; i++ ){
			OutNum++ ;
			Output << setw(18) << pResult->RotTemp[i] ;
			if ( OutNum%6 == 0 ) Output << '\n' ;
		}

		for ( int i=0 ; i<CellNum ; i++ ){
			OutNum++ ;
			Output << setw(18) << pResult->VibTemp[i] ;
			if ( OutNum%6 == 0 ) Output << '\n' ;
		}

		for ( int i=0 ; i<CellNum ; i++ ){
			OutNum++ ;
			Output << setw(18) << pResult->AveParticleNum[i] ;
			if ( OutNum%6 == 0 ) Output << '\n' ;
		}
		
		for ( int i=0 ; i<CellNum ; i++ ){
			OutNum++ ;
			Output << setw(18) << pResult->MeanCollSpacingMeanFreePath[i] ;
			if ( OutNum%6 == 0 ) Output << '\n' ;
		}
		
		for ( int i=0 ; i<CellNum ; i++ ){
			OutNum++ ;
			Output << setw(18) << Cell[i].ProcessorNo ;
			if ( OutNum%6 == 0 ) Output << '\n' ;
		}


		Output << "\n\n" ;
		for ( int i=0 ; i<CellNum ; i++ ){
			for ( int j=0 ; j<pMapping->NodeNum[Cell[i].Type] ; j++ )
				Output << setw(18) << (Cell[i].Node[j]+1) ;
			if ( Cell[i].Type == 4 ) Output << setw(18) << (Cell[i].Node[2]+1) ;
			Output << '\n' ;
		}
	}else{
		cout << "Fail to open Result-tec.dat \n" ;
		exit(1) ;
	}

	Output.clear() ;
	Output.close() ;
}


void OutptuResultTec3D( DSMC_DOMAIN *pDomain , DSMC_NODE *Node , DSMC_CELL *Cell , DSMC_RESULT *pResult , CELLMAPPING *pMapping ){
	ofstream	Output ;
	int		NodeNum , CellNum , OutNum ;
	
	NodeNum	= pDomain->NodeNum ;
	CellNum	= pDomain->CellNum ;
	OutNum	= 0 ;

	Output.open( "Result-tec.dat" , ios::out | ios::trunc ) ;

	if ( !Output.fail() ){
		Output << "TITLE = \"Finite volume dataset\"\n" ;
		Output << "VARIABLES = \"X\", \"Y\", \"Z\", \"Density\", \"NumDensity\", \"U-Vel\", \"V-Vel\", \"W-Vel\", \"Temp\", \"TransTemp\", \"RotTemp\", \"VibTemp\", \"#P\", \"mcs/mfp\", \"ProcessorNo\"\n" ;
		Output << "ZONE N=" << NodeNum << ", E=" << CellNum << ", DATAPACKING=BLOCK, ZONETYPE=FEBRICK\n" ;
		Output << "VARLOCATION=([1-3]=NODAL, [4-15]=CELLCENTERED)\n" ;

		
		for ( int i=0 ; i<NodeNum ; i++ ){
			OutNum++ ;
			Output << setw(18) << Node[i].XCoord ;
			if ( OutNum%6 == 0 ) Output << '\n' ;
		}
		
		for ( int i=0 ; i<NodeNum ; i++ ){
			OutNum++ ;
			Output << setw(18) << Node[i].YCoord ;
			if ( OutNum%6 == 0 ) Output << '\n' ;
		}
		
		for ( int i=0 ; i<NodeNum ; i++ ){
			OutNum++ ;
			Output << setw(18) << Node[i].ZCoord ;
			if ( OutNum%6 == 0 ) Output << '\n' ;
		}

		for ( int i=0 ; i<CellNum ; i++ ){
			OutNum++ ;
			Output << setw(18) << pResult->Density[i] ;
			if ( OutNum%6 == 0 ) Output << '\n' ;
		}

		for ( int i=0 ; i<CellNum ; i++ ){
			OutNum++ ;
			Output << setw(18) << pResult->NumDensity[i] ;
			if ( OutNum%6 == 0 ) Output << '\n' ;
		}

		for ( int i=0 ; i<CellNum ; i++ ){
			OutNum++ ;
			Output << setw(18) << pResult->XVel[i] ;
			if ( OutNum%6 == 0 ) Output << '\n' ;
		}
		
		for ( int i=0 ; i<CellNum ; i++ ){
			OutNum++ ;
			Output << setw(18) << pResult->YVel[i] ;
			if ( OutNum%6 == 0 ) Output << '\n' ;
		}
		
		for ( int i=0 ; i<CellNum ; i++ ){
			OutNum++ ;
			Output << setw(18) << pResult->ZVel[i] ;
			if ( OutNum%6 == 0 ) Output << '\n' ;
		}

		for ( int i=0 ; i<CellNum ; i++ ){
			OutNum++ ;
			Output << setw(18) << pResult->Temp[i] ;
			if ( OutNum%6 == 0 ) Output << '\n' ;
		}

		for ( int i=0 ; i<CellNum ; i++ ){
			OutNum++ ;
			Output << setw(18) << pResult->TransTemp[i] ;
			if ( OutNum%6 == 0 ) Output << '\n' ;
		}

		for ( int i=0 ; i<CellNum ; i++ ){
			OutNum++ ;
			Output << setw(18) << pResult->RotTemp[i] ;
			if ( OutNum%6 == 0 ) Output << '\n' ;
		}

		for ( int i=0 ; i<CellNum ; i++ ){
			OutNum++ ;
			Output << setw(18) << pResult->VibTemp[i] ;
			if ( OutNum%6 == 0 ) Output << '\n' ;
		}

		for ( int i=0 ; i<CellNum ; i++ ){
			OutNum++ ;
			Output << setw(18) << pResult->AveParticleNum[i] ;
			if ( OutNum%6 == 0 ) Output << '\n' ;
		}
		
		for ( int i=0 ; i<CellNum ; i++ ){
			OutNum++ ;
			Output << setw(18) << pResult->MeanCollSpacingMeanFreePath[i] ;
			if ( OutNum%6 == 0 ) Output << '\n' ;
		}
		
		for ( int i=0 ; i<CellNum ; i++ ){
			OutNum++ ;
			Output << setw(18) << Cell[i].ProcessorNo ;
			if ( OutNum%6 == 0 ) Output << '\n' ;
		}


		Output << "\n\n" ;
		for ( int i=0 ; i<CellNum ; i++ ){
			if ( Cell[i].Type == 0 ){
				Output << setw(18) << Cell[i].Node[0]+1 ;
				Output << setw(18) << Cell[i].Node[1]+1 ;
				Output << setw(18) << Cell[i].Node[2]+1 ;
				Output << setw(18) << Cell[i].Node[2]+1 ;
				Output << setw(18) << Cell[i].Node[3]+1 ;
				Output << setw(18) << Cell[i].Node[3]+1 ;
				Output << setw(18) << Cell[i].Node[3]+1 ;
				Output << setw(18) << Cell[i].Node[3]+1 ;

			// Hex. Cell
			}else if ( Cell[i].Type == 1 ){
				for ( int j=0 ; j<pMapping->NodeNum[1] ; j++ )
					Output << setw(18) << Cell[i].Node[j]+1 ;

			// Pyramid Cell
			}else if ( Cell[i].Type == 2 ){
				Output << setw(18) << Cell[i].Node[0]+1 ;
				Output << setw(18) << Cell[i].Node[1]+1 ;
				Output << setw(18) << Cell[i].Node[2]+1 ;
				Output << setw(18) << Cell[i].Node[3]+1 ;
				Output << setw(18) << Cell[i].Node[4]+1 ;
				Output << setw(18) << Cell[i].Node[4]+1 ;
				Output << setw(18) << Cell[i].Node[4]+1 ;
				Output << setw(18) << Cell[i].Node[4]+1 ;


			// Prism Cell
			}else if ( Cell[i].Type == 3 ){
				Output << setw(18) << Cell[i].Node[0]+1 ;
				Output << setw(18) << Cell[i].Node[3]+1 ;
				Output << setw(18) << Cell[i].Node[4]+1 ;
				Output << setw(18) << Cell[i].Node[4]+1 ;
				Output << setw(18) << Cell[i].Node[1]+1 ;
				Output << setw(18) << Cell[i].Node[2]+1 ;
				Output << setw(18) << Cell[i].Node[5]+1 ;
				Output << setw(18) << Cell[i].Node[5]+1 ;


			// 2D Triangular Cell
			/*}else if ( Cell[i].Type == 4 ){
				Output << setw(18) << Cell[i].Node[0]+1 ;
				Output << setw(18) << Cell[i].Node[1]+1 ;
				Output << setw(18) << Cell[i].Node[2]+1 ;
				Output << setw(18) << Cell[i].Node[2]+1 ;


			// 2D Quadrilateral Cell
			}else if ( Cell[i].Type == 5 ){
				Output << setw(18) << Cell[i].Node[0]+1 ;
				Output << setw(18) << Cell[i].Node[1]+1 ;
				Output << setw(18) << Cell[i].Node[2]+1 ;
				Output << setw(18) << Cell[i].Node[3]+1 ;*/
			}
			
			Output << '\n' ;
		}
	}else{
		cout << "Fail to open Result-tec.dat \n" ;
		exit(1) ;
	}

	Output.clear() ;
	Output.close() ;
}

//==============================================================================================================

void OutptuResultSpeciesTec2D( DSMC_DOMAIN *pDomain , DSMC_NODE *Node , DSMC_CELL *Cell , DSMC_RESULT *pResult , CELLMAPPING *pMapping , string SpeciesName ){
	ofstream	Output ;
	int		NodeNum , CellNum , OutNum ;
	string		FileName ;
	
	NodeNum		= pDomain->NodeNum ;
	CellNum		= pDomain->CellNum ;
	OutNum		= 0 ;
	
	FileName	= "Result-" + SpeciesName + "-tec.dat" ;

	Output.open( FileName.c_str() , ios::out | ios::trunc ) ;

	if ( !Output.fail() ){
		Output << "TITLE = \"Finite volume dataset\"\n" ;
		Output << "VARIABLES = \"X\", \"Y\", \"Density\", \"NumDensity\", \"U-Vel\", \"V-Vel\", \"W-Vel\", \"Temp\", \"TransTemp\", \"RotTemp\", \"VibTemp\", \"TransXTemp\", \"TransYTemp\", \"TransZTemp\", \"#P\"\n" ;
		Output << "ZONE N=" << NodeNum << ", E=" << CellNum << ", DATAPACKING=BLOCK, ZONETYPE=FEQUADRILATERAL\n" ;
		Output << "VARLOCATION=([1-2]=NODAL, [3-15]=CELLCENTERED)\n" ;

		
		for ( int i=0 ; i<NodeNum ; i++ ){
			OutNum++ ;
			Output << setw(18) << Node[i].XCoord ;
			if ( OutNum%6 == 0 ) Output << '\n' ;
		}
		
		for ( int i=0 ; i<NodeNum ; i++ ){
			OutNum++ ;
			Output << setw(18) << Node[i].YCoord ;
			if ( OutNum%6 == 0 ) Output << '\n' ;
		}

		for ( int i=0 ; i<CellNum ; i++ ){
			OutNum++ ;
			Output << setw(18) << pResult->Density[i] ;
			if ( OutNum%6 == 0 ) Output << '\n' ;
		}

		for ( int i=0 ; i<CellNum ; i++ ){
			OutNum++ ;
			Output << setw(18) << pResult->NumDensity[i] ;
			if ( OutNum%6 == 0 ) Output << '\n' ;
		}

		for ( int i=0 ; i<CellNum ; i++ ){
			OutNum++ ;
			Output << setw(18) << pResult->XVel[i] ;
			if ( OutNum%6 == 0 ) Output << '\n' ;
		}
		
		for ( int i=0 ; i<CellNum ; i++ ){
			OutNum++ ;
			Output << setw(18) << pResult->YVel[i] ;
			if ( OutNum%6 == 0 ) Output << '\n' ;
		}
		
		for ( int i=0 ; i<CellNum ; i++ ){
			OutNum++ ;
			Output << setw(18) << pResult->ZVel[i] ;
			if ( OutNum%6 == 0 ) Output << '\n' ;
		}

		for ( int i=0 ; i<CellNum ; i++ ){
			OutNum++ ;
			Output << setw(18) << pResult->Temp[i] ;
			if ( OutNum%6 == 0 ) Output << '\n' ;
		}

		for ( int i=0 ; i<CellNum ; i++ ){
			OutNum++ ;
			Output << setw(18) << pResult->TransTemp[i] ;
			if ( OutNum%6 == 0 ) Output << '\n' ;
		}

		for ( int i=0 ; i<CellNum ; i++ ){
			OutNum++ ;
			Output << setw(18) << pResult->RotTemp[i] ;
			if ( OutNum%6 == 0 ) Output << '\n' ;
		}

		for ( int i=0 ; i<CellNum ; i++ ){
			OutNum++ ;
			Output << setw(18) << pResult->VibTemp[i] ;
			if ( OutNum%6 == 0 ) Output << '\n' ;
		}

		for ( int i=0 ; i<CellNum ; i++ ){
			OutNum++ ;
			Output << setw(18) << pResult->TransXTemp[i] ;
			if ( OutNum%6 == 0 ) Output << '\n' ;
		}
		
		for ( int i=0 ; i<CellNum ; i++ ){
			OutNum++ ;
			Output << setw(18) << pResult->TransYTemp[i] ;
			if ( OutNum%6 == 0 ) Output << '\n' ;
		}
		
		for ( int i=0 ; i<CellNum ; i++ ){
			OutNum++ ;
			Output << setw(18) << pResult->TransZTemp[i] ;
			if ( OutNum%6 == 0 ) Output << '\n' ;
		}
		
		for ( int i=0 ; i<CellNum ; i++ ){
			OutNum++ ;
			Output << setw(18) << pResult->AveParticleNum[i] ;
			if ( OutNum%6 == 0 ) Output << '\n' ;
		}



		Output << "\n\n" ;
		for ( int i=0 ; i<CellNum ; i++ ){
			for ( int j=0 ; j<pMapping->NodeNum[Cell[i].Type] ; j++ )
				Output << setw(18) << (Cell[i].Node[j]+1) ;
			if ( Cell[i].Type == 4 ) Output << setw(18) << (Cell[i].Node[2]+1) ;
			Output << '\n' ;
		}
	}else{
		cout << "Fail to open " << FileName << '\n' ;
		exit(1) ;
	}

	Output.clear() ;
	Output.close() ;
	
}


void OutptuResultSpeciesTec3D( DSMC_DOMAIN *pDomain , DSMC_NODE *Node , DSMC_CELL *Cell , DSMC_RESULT *pResult , CELLMAPPING *pMapping , string SpeciesName ){
	ofstream	Output ;
	int		NodeNum , CellNum , OutNum ;
	string		FileName ;
	
	NodeNum		= pDomain->NodeNum ;
	CellNum		= pDomain->CellNum ;
	OutNum		= 0 ;
	
	FileName	= "Result-" + SpeciesName + "-tec.dat" ;

	Output.open( FileName.c_str() , ios::out | ios::trunc ) ;

	if ( !Output.fail() ){
		Output << "TITLE = \"Finite volume dataset\"\n" ;
		Output << "VARIABLES = \"X\", \"Y\", \"Z\", \"Density\", \"NumDensity\", \"U-Vel\", \"V-Vel\", \"W-Vel\", \"Temp\", \"TransTemp\", \"RotTemp\", \"VibTemp\", \"TransXTemp\", \"TransYTemp\", \"TransZTemp\", \"#P\"\n" ;
		Output << "ZONE N=" << NodeNum << ", E=" << CellNum << ", DATAPACKING=BLOCK, ZONETYPE=FEBRICK\n" ;
		Output << "VARLOCATION=([1-3]=NODAL, [4-16]=CELLCENTERED)\n" ;

		
		for ( int i=0 ; i<NodeNum ; i++ ){
			OutNum++ ;
			Output << setw(18) << Node[i].XCoord ;
			if ( OutNum%6 == 0 ) Output << '\n' ;
		}
		
		for ( int i=0 ; i<NodeNum ; i++ ){
			OutNum++ ;
			Output << setw(18) << Node[i].YCoord ;
			if ( OutNum%6 == 0 ) Output << '\n' ;
		}
		
		for ( int i=0 ; i<NodeNum ; i++ ){
			OutNum++ ;
			Output << setw(18) << Node[i].ZCoord ;
			if ( OutNum%6 == 0 ) Output << '\n' ;
		}

		for ( int i=0 ; i<CellNum ; i++ ){
			OutNum++ ;
			Output << setw(18) << pResult->Density[i] ;
			if ( OutNum%6 == 0 ) Output << '\n' ;
		}

		for ( int i=0 ; i<CellNum ; i++ ){
			OutNum++ ;
			Output << setw(18) << pResult->NumDensity[i] ;
			if ( OutNum%6 == 0 ) Output << '\n' ;
		}

		for ( int i=0 ; i<CellNum ; i++ ){
			OutNum++ ;
			Output << setw(18) << pResult->XVel[i] ;
			if ( OutNum%6 == 0 ) Output << '\n' ;
		}
		
		for ( int i=0 ; i<CellNum ; i++ ){
			OutNum++ ;
			Output << setw(18) << pResult->YVel[i] ;
			if ( OutNum%6 == 0 ) Output << '\n' ;
		}
		
		for ( int i=0 ; i<CellNum ; i++ ){
			OutNum++ ;
			Output << setw(18) << pResult->ZVel[i] ;
			if ( OutNum%6 == 0 ) Output << '\n' ;
		}

		for ( int i=0 ; i<CellNum ; i++ ){
			OutNum++ ;
			Output << setw(18) << pResult->Temp[i] ;
			if ( OutNum%6 == 0 ) Output << '\n' ;
		}

		for ( int i=0 ; i<CellNum ; i++ ){
			OutNum++ ;
			Output << setw(18) << pResult->TransTemp[i] ;
			if ( OutNum%6 == 0 ) Output << '\n' ;
		}

		for ( int i=0 ; i<CellNum ; i++ ){
			OutNum++ ;
			Output << setw(18) << pResult->RotTemp[i] ;
			if ( OutNum%6 == 0 ) Output << '\n' ;
		}

		for ( int i=0 ; i<CellNum ; i++ ){
			OutNum++ ;
			Output << setw(18) << pResult->VibTemp[i] ;
			if ( OutNum%6 == 0 ) Output << '\n' ;
		}

		for ( int i=0 ; i<CellNum ; i++ ){
			OutNum++ ;
			Output << setw(18) << pResult->TransXTemp[i] ;
			if ( OutNum%6 == 0 ) Output << '\n' ;
		}
		
		for ( int i=0 ; i<CellNum ; i++ ){
			OutNum++ ;
			Output << setw(18) << pResult->TransYTemp[i] ;
			if ( OutNum%6 == 0 ) Output << '\n' ;
		}
		
		for ( int i=0 ; i<CellNum ; i++ ){
			OutNum++ ;
			Output << setw(18) << pResult->TransZTemp[i] ;
			if ( OutNum%6 == 0 ) Output << '\n' ;
		}
		
		for ( int i=0 ; i<CellNum ; i++ ){
			OutNum++ ;
			Output << setw(18) << pResult->AveParticleNum[i]  ;
			if ( OutNum%6 == 0 ) Output << '\n' ;
		}


		Output << "\n\n" ;
		for ( int i=0 ; i<CellNum ; i++ ){
			if ( Cell[i].Type == 0 ){
				Output << setw(18) << Cell[i].Node[0]+1 ;
				Output << setw(18) << Cell[i].Node[1]+1 ;
				Output << setw(18) << Cell[i].Node[2]+1 ;
				Output << setw(18) << Cell[i].Node[2]+1 ;
				Output << setw(18) << Cell[i].Node[3]+1 ;
				Output << setw(18) << Cell[i].Node[3]+1 ;
				Output << setw(18) << Cell[i].Node[3]+1 ;
				Output << setw(18) << Cell[i].Node[3]+1 ;

			// Hex. Cell
			}else if ( Cell[i].Type == 1 ){
				for ( int j=0 ; j<pMapping->NodeNum[1] ; j++ )
					Output << setw(18) << Cell[i].Node[j]+1 ;

			// Pyramid Cell
			}else if ( Cell[i].Type == 2 ){
				Output << setw(18) << Cell[i].Node[0]+1 ;
				Output << setw(18) << Cell[i].Node[1]+1 ;
				Output << setw(18) << Cell[i].Node[2]+1 ;
				Output << setw(18) << Cell[i].Node[3]+1 ;
				Output << setw(18) << Cell[i].Node[4]+1 ;
				Output << setw(18) << Cell[i].Node[4]+1 ;
				Output << setw(18) << Cell[i].Node[4]+1 ;
				Output << setw(18) << Cell[i].Node[4]+1 ;


			// Prism Cell
			}else if ( Cell[i].Type == 3 ){
				Output << setw(18) << Cell[i].Node[0]+1 ;
				Output << setw(18) << Cell[i].Node[3]+1 ;
				Output << setw(18) << Cell[i].Node[4]+1 ;
				Output << setw(18) << Cell[i].Node[4]+1 ;
				Output << setw(18) << Cell[i].Node[1]+1 ;
				Output << setw(18) << Cell[i].Node[2]+1 ;
				Output << setw(18) << Cell[i].Node[5]+1 ;
				Output << setw(18) << Cell[i].Node[5]+1 ;
			}
			
			Output << '\n' ;
		}
	}else{
		cout << "Fail to open " << FileName << '\n' ;
		exit(1) ;
	}

	Output.clear() ;
	Output.close() ;
}

//==============================================================================================================
//==============================================================================================================

void RemoveSurface( DSMC_DOMAIN *pDomain , DSMC_POST_SURFACE *Surface , DSMC_POST_SURFACE_CONDITION *pSurfaceCondition ){
	int		WallFaceNum , WallNo ;
	
	WallFaceNum	= 0 ;
	
	for ( int i=0 ; i<pDomain->WallFaceNum ; i++ ){
		WallNo	= Surface[i].WallNo ;
		
		for ( int j=0 ; j<pSurfaceCondition->WallNum ; j++ ){
			if ( WallNo == pSurfaceCondition->WallNo[j] ){
				Surface[WallFaceNum]	= Surface[i] ;
				WallFaceNum++ ;
				
				break ;
			}
		}
	}
	
	pDomain->WallFaceNum	= WallFaceNum ;
}

//==============================================================================================================

void CalculateAngle( DSMC_DOMAIN *pDomain , DSMC_POST_SURFACE *Surface ){
	int		WallFaceNum ;
	double		a , b , c , x0 , y0 , x1 , y1 ;
	
	x0	= 9.03E-4 ;
	y0	= 0. ;
	c	= 9.03E-4 ;
	
	WallFaceNum	= pDomain->WallFaceNum ;
	
	for ( int i=0 ; i<WallFaceNum ; i++ ){
		x1	= Surface[i].Center[0] ;
		y1	= Surface[i].Center[1] ;
		
		a	= x1 - x0 ;
		b	= y1 - y0 ;
		
		Surface[i].Center[3]	= CalculateAngleDeg( a , b , c ) ;
		Surface[i].Center[3]	= 180. - Surface[i].Center[3] ;
	}
}

//==============================================================================================================

double CalculateAngleDeg( double a , double b , double c ){
	double	pi	=	3.141592654 ;
	double	angle ;
	
	angle	=	asin(b/c) ;
	angle	=	angle/pi*180. ;
	if ( a < 0. )
		angle = 180. - angle ;
	else if ( angle < 0. )
		angle = 360. + angle ;
	
	return	angle ;
}

//==============================================================================================================

void SortSurface( DSMC_DOMAIN *pDomain , DSMC_POST_SURFACE *Surface , DSMC_POST_SURFACE_CONDITION *pSurfaceCondition ){
	DSMC_POST_SURFACE	Buffer ;
	int			WallFaceNum , Direction ;
	double			Coord[2] ;
	
	WallFaceNum	= pDomain->WallFaceNum ;
	
	if ( pSurfaceCondition->SortName == "X-Direction" ){
		Direction	= 0 ;
	}else if ( pSurfaceCondition->SortName == "Y-Direction" ){
		Direction	= 1 ;
	}else if ( pSurfaceCondition->SortName == "Z-Direction" ){
		Direction	= 2 ;
	}else if ( pSurfaceCondition->SortName == "Angle" ){
		Direction	= 3 ;
	}else{
		cout << "SortName is ERROR. SortName is " << pSurfaceCondition->SortName << ". \n" ;
		exit(1) ;
	}
	
	for ( int i=0 ; i<WallFaceNum ; i++ ){
		for ( int j=0 ; j<(WallFaceNum-i-1) ; j++ ){
			Coord[0]	= Surface[j].Center[Direction] ;
			Coord[1]	= Surface[j+1].Center[Direction] ;
			
			if ( Coord[0] > Coord[1] ){
				Buffer		= Surface[j] ;
				Surface[j]	= Surface[j+1] ;
				Surface[j+1]	= Buffer ;
			}
		}
	}
}

//==============================================================================================================

void CalculateSurfaceProperty( DSMC_DOMAIN *pDomain , DSMC_POST_SURFACE *Surface , DSMC_POST_SURFACE_CONDITION *pSurfaceCondition ){
	if ( pDomain->Dimension == 2 ){
		CalculateSurfaceProperty2D( pDomain , Surface , pSurfaceCondition ) ;
	}else if ( pDomain->Dimension == 3 ){
		CalculateSurfaceProperty3D( pDomain , Surface , pSurfaceCondition ) ;
	}
}

//==============================================================================================================

void CalculateSurfaceProperty2D( DSMC_DOMAIN *pDomain , DSMC_POST_SURFACE *Surface , DSMC_POST_SURFACE_CONDITION *pSurfaceCondition ){
	int		WallFaceNum ;
	double		Density , Velocity ;
	double		PressureX , PressureY , PressureNorm , PressurePara , HeatFlux , Drag , PeakHeatFlux ;
	
	WallFaceNum	= pDomain->WallFaceNum ;
	Density		= pSurfaceCondition->Density ;
	Velocity	= pSurfaceCondition->Velocity ;
	Drag		= 0. ;
	PeakHeatFlux	= 0. ;
	
	for ( int i=0 ; i<WallFaceNum ; i++ ){
		PressureNorm	= Surface[i].SurfacePro[3] + Surface[i].SurfacePro[4] ;
		PressureX	= Surface[i].SurfacePro[5] + Surface[i].SurfacePro[6] ;
		PressureY	= Surface[i].SurfacePro[7] + Surface[i].SurfacePro[8] ;
		HeatFlux	= Surface[i].SurfacePro[11] + Surface[i].SurfacePro[12] + 
				  Surface[i].SurfacePro[13] + Surface[i].SurfacePro[14] + 
				  Surface[i].SurfacePro[15] + Surface[i].SurfacePro[16] ;
		
		Drag		+= PressureX * Surface[i].Area ;
		if ( HeatFlux > PeakHeatFlux ) PeakHeatFlux = HeatFlux ;
		
		PressurePara	= sqrt( fabs( pow(PressureX,2) + pow(PressureY,2) - pow(PressureNorm,2) ) ) ;
		
		Surface[i].Cp	= PressureNorm / (0.5*Density*Velocity*Velocity) ;
		Surface[i].Cf	= PressurePara / (0.5*Density*Velocity*Velocity) ;
		Surface[i].Ch	= HeatFlux / (0.5*Density*pow(Velocity,3)) ;
		Surface[i].Pressure	= PressureNorm ;
		Surface[i].HeatFlux	= HeatFlux ;
	}
	
	pSurfaceCondition->Drag		= Drag ;
	pSurfaceCondition->HeatFlux	= PeakHeatFlux ;
}

//==============================================================================================================

void CalculateSurfaceProperty3D( DSMC_DOMAIN *pDomain , DSMC_POST_SURFACE *Surface , DSMC_POST_SURFACE_CONDITION *pSurfaceCondition ){
	int		WallFaceNum ;
	double		Density , Velocity ,Area ;
	double		DepositionRate , PressureX , PressureY , PressureZ , PressureNorm , PressurePara , HeatFlux , Drag , PeakHeatFlux, Lift ;
	
	WallFaceNum	= pDomain->WallFaceNum ;
	Density		= pSurfaceCondition->Density ;
	Velocity	= pSurfaceCondition->Velocity ;
	DepositionRate	= 0. ;
	Drag		= 0. ;
	Lift    = 0. ;
	PeakHeatFlux	= 0. ;
	PressurePara	= 0. ;
	Area    = 0.59385 ;

	
	for ( int i=0 ; i<WallFaceNum ; i++ ){
		DepositionRate	= Surface[i].SurfacePro[2]*pSurfaceCondition->Mass ;
		PressureNorm	= Surface[i].SurfacePro[3] + Surface[i].SurfacePro[4] ;
		PressureX	= Surface[i].SurfacePro[5] + Surface[i].SurfacePro[6] ;
		PressureY	= Surface[i].SurfacePro[7] + Surface[i].SurfacePro[8] ;
		PressureZ	= Surface[i].SurfacePro[9] + Surface[i].SurfacePro[10] ;
		HeatFlux	= Surface[i].SurfacePro[11] + Surface[i].SurfacePro[12] + 
				  Surface[i].SurfacePro[13] + Surface[i].SurfacePro[14] + 
				  Surface[i].SurfacePro[15] + Surface[i].SurfacePro[16] ;
		Lift    += PressureY * Surface[i].Area ;
		Drag		+= PressureX * Surface[i].Area ;
//		Area    += Surface[i].Area ;
		if ( HeatFlux > PeakHeatFlux ) PeakHeatFlux = HeatFlux ;
		
		//PressurePara	= sqrt( fabs( pow(PressureX,2) + pow(PressureY,2) - pow(PressureNorm,2) ) ) ;
		
		Surface[i].Cp			= PressureNorm / (0.5*Density*Velocity*Velocity) ;
		//Surface[i].Cf			= PressurePara / (0.5*Density*Velocity*Velocity) ;
		Surface[i].Ch			= HeatFlux / (0.5*Density*pow(Velocity,3)) ;
		Surface[i].DepositionRate	= DepositionRate ;
		Surface[i].Pressure		= PressureNorm ;
		Surface[i].HeatFlux		= HeatFlux ;
		
		
		// Debug.
		//cout << "DS: " << Surface[i].DepositionRate << ", " << Surface[i].SurfacePro[2] << ", " << pSurfaceCondition->Mass << '\n' ;
	}
	
//	Drag = Drag /(0.5*Density*Velocity*Velocity*Area);
	
	pSurfaceCondition->Drag		= Drag ;
	pSurfaceCondition->CoeffAxial = Drag /(0.5*Density*Velocity*Velocity*Area);
	pSurfaceCondition->CoeffNomal = Lift /(0.5*Density*Velocity*Velocity*Area);
	pSurfaceCondition->HeatFlux	= PeakHeatFlux ;
}

//==============================================================================================================

void OutputSurface( DSMC_NODE	*Node , DSMC_CELL *Cell , DSMC_DOMAIN *pDomain , DSMC_POST_SURFACE *Surface , CELLMAPPING *pMapping ){
	if ( pDomain->Dimension == 2 ){
		OutputSurface2D( pDomain , Surface ) ;
	}else if ( pDomain->Dimension == 3 ){
		OutputSurface3D( Node , Cell , Surface , pDomain , pMapping ) ;
	}
}

//==============================================================================================================

void OutputSurface2D( DSMC_DOMAIN *pDomain , DSMC_POST_SURFACE *Surface ){
	ofstream	Output ;
	int		WallFaceNum ;
	
	WallFaceNum	= pDomain->WallFaceNum ;
	
	
	Output.open( "Result-Surface2D.dat" , ios::out | ios::trunc ) ;
		
		
	if ( !Output.fail() ){
		Output	<< setw(16) << "XCoord" << setw(16) << "YCoord" << setw(16) << "Angle" << setw(16) << "Pressure" << setw(16) << "HeatingRate" 
			<< setw(16) << "Cp" << setw(16) << "Cf" << setw(16) << "Ch" << setw(16) << "Area" << setw(16) << "WallNo" << '\n' ;
			
		for ( int i=0 ; i<WallFaceNum ; i++ ){
			Output	<< setw(16) << Surface[i].Center[0] << setw(16) << Surface[i].Center[1] << setw(16) << Surface[i].Center[3] 
				<< setw(16) << Surface[i].Pressure << setw(16) << Surface[i].HeatFlux << setw(16) << Surface[i].Cp 
				<< setw(16) << Surface[i].Cf << setw(16) << Surface[i].Ch << setw(16) << Surface[i].Area 
				<< setw(16) << Surface[i].WallNo << '\n' ;
		}
	}else{
		cout << "Fail to open Result-Surface2D.dat \n" ;
		exit(1) ;
	}
	Output.clear() ;
	Output.close() ;
	
	
	Output.open( "Result-Surface2D-tec.dat" , ios::out | ios::trunc ) ;
		
	if ( !Output.fail() ){
		Output << "VARIABLES = \"XCoord\" \"YCoord\" \"Angle\" \"Pressure\" \"HeatingRate\" \"Cp\" \"Cf\" \"Ch\" \"Area\" \"WallNo\"\n" ;
		Output << "ZONE\n" ;
		Output << "DT=(DOUBLE,DOUBLE,DOUBLE,DOUBLE,DOUBLE,DOUBLE,DOUBLE,DOUBLE,DOUBLE,DOUBLE)\n" ;
			
			
		for ( int i=0 ; i<WallFaceNum ; i++ ){
			Output	<< setw(16) << Surface[i].Center[0] << setw(16) << Surface[i].Center[1] << setw(16) << Surface[i].Center[3] 
				<< setw(16) << Surface[i].Pressure << setw(16) << Surface[i].HeatFlux << setw(16) << Surface[i].Cp 
				<< setw(16) << Surface[i].Cf << setw(16) << Surface[i].Ch << setw(16) << Surface[i].Area 
				<< setw(16) << Surface[i].WallNo << '\n' ;
		}
	}else{
		cout << "Fail to open Result-Surface2D-tec.dat \n" ;
		exit(1) ;
	}
	Output.clear() ;
	Output.close() ;
}


void OutputSurface3D( DSMC_NODE	*Node , DSMC_CELL *Cell , DSMC_POST_SURFACE *Surface , DSMC_DOMAIN *pDomain , CELLMAPPING *pMapping ){
	int			CellNum , NodeNum , CellNo , NodeNo , FaceNo , Type , OutNum ;
	DSMC_CELL		*sCell ;
	DSMC_POST_NODECELL	*NodeCell ;
	ofstream		Output ;
	
	
	sCell		= new DSMC_CELL[pDomain->WallFaceNum] ;
	NodeCell	= new DSMC_POST_NODECELL[pDomain->NodeNum] ;
	
	
	CellNum		= 0 ;
	NodeNum		= 0 ;
	FaceNo		= 0 ;
	Type		= 0 ;
	OutNum		= 0 ;
	
	
	for ( int i=0 ; i<pDomain->WallFaceNum ; i++ ){
		CellNo		= Surface[i].CellNo ;
		FaceNo		= Surface[i].FaceNo ;
		sCell[i]	= Cell[CellNo] ;
		Type		= Cell[CellNo].Type ;
		
		
		if ( pMapping->SurfaceNode[Type][FaceNo] == 3 )
			sCell[i].Type	= 4 ;
		else if ( pMapping->SurfaceNode[Type][FaceNo] == 4 )
			sCell[i].Type	= 5 ;


		for ( int j=0 ; j<pMapping->SurfaceNode[Type][FaceNo] ; j++ )
			sCell[i].Node[j]	= Cell[CellNo].Node[pMapping->Node[Type][FaceNo][j]] ;
			
		CellNum++ ;
	}
	
	
	
	for ( int i=0 ; i<CellNum ; i++ ){
		Type	= sCell[i].Type ;
		
		for ( int j=0 ; j<pMapping->NodeNum[Type] ; j++ ){
			NodeNo	= sCell[i].Node[j] ;

			NodeCell[NodeNo].CellRelation.push_back(i) ;
		}
	}
	
	
	for ( int i=0 ; i<pDomain->NodeNum ; i++ ){
		if ( NodeCell[i].CellRelation.size() > 0 ){
			for ( int j=0 ; j<NodeCell[i].CellRelation.size() ; j++ ){
				CellNo	= NodeCell[i].CellRelation[j] ;
				Type	= sCell[CellNo].Type ;
				
				
				for ( int k=0 ; k<pMapping->NodeNum[Type] ; k++ ){
					if ( sCell[CellNo].Node[k] == i ){
						sCell[CellNo].Node[k]	= NodeNum ;
					}
				}
			}
			
			Node[NodeNum]	= Node[i] ;
			
			NodeNum++ ;
		}
	}
	
	
	Output.open( "Result-Surface3D-tec.dat" , ios::out | ios::trunc ) ;
		
		
	if ( !Output.fail() ){
		Output << "TITLE = \"Finite volume dataset\"\n" ;
		Output << "VARIABLES = \"X\", \"Y\", \"Z\", \"DepositionRate\", \"Pressure\", \"HeatingRate\", \"Cp\", \"Cf\", \"Ch\"\n" ;
		Output << "ZONE N=" << NodeNum << ", E=" << CellNum << ", DATAPACKING=BLOCK, ZONETYPE=FEQUADRILATERAL\n" ;
		Output << "VARLOCATION=([1-3]=NODAL, [4-9]=CELLCENTERED)\n" ;
		
		
		
		for ( int i=0 ; i<NodeNum ; i++ ){
			OutNum++ ;
			Output << setw(18) << Node[i].XCoord ;
			if ( OutNum%6 == 0 ) Output << '\n' ;
		}
		
		for ( int i=0 ; i<NodeNum ; i++ ){
			OutNum++ ;
			Output << setw(18) << Node[i].YCoord ;
			if ( OutNum%6 == 0 ) Output << '\n' ;
		}
		
		for ( int i=0 ; i<NodeNum ; i++ ){
			OutNum++ ;
			Output << setw(18) << Node[i].ZCoord ;
			if ( OutNum%6 == 0 ) Output << '\n' ;
		}
		
		for ( int i=0 ; i<CellNum ; i++ ){
			OutNum++ ;
			Output << setw(18) << Surface[i].DepositionRate ;
			if ( OutNum%6 == 0 ) Output << '\n' ;
		}
		
		for ( int i=0 ; i<CellNum ; i++ ){
			OutNum++ ;
			Output << setw(18) << Surface[i].Pressure ;
			if ( OutNum%6 == 0 ) Output << '\n' ;
		}
		
		for ( int i=0 ; i<CellNum ; i++ ){
			OutNum++ ;
			Output << setw(18) << Surface[i].HeatFlux ;
			if ( OutNum%6 == 0 ) Output << '\n' ;
		}
		
		for ( int i=0 ; i<CellNum ; i++ ){
			OutNum++ ;
			Output << setw(18) << Surface[i].Cp ;
			if ( OutNum%6 == 0 ) Output << '\n' ;
		}
		
		for ( int i=0 ; i<CellNum ; i++ ){
			OutNum++ ;
			Output << setw(18) << Surface[i].Cf ;
			if ( OutNum%6 == 0 ) Output << '\n' ;
		}
		
		for ( int i=0 ; i<CellNum ; i++ ){
			OutNum++ ;
			Output << setw(18) << Surface[i].Ch ;
			if ( OutNum%6 == 0 ) Output << '\n' ;
		}
		
		
		Output << "\n\n" ;
		for ( int i=0 ; i<CellNum ; i++ ){
			for ( int j=0 ; j<pMapping->NodeNum[sCell[i].Type] ; j++ )
				Output << setw(18) << (sCell[i].Node[j]+1) ;
			if ( sCell[i].Type == 4 ) Output << setw(18) << (sCell[i].Node[2]+1) ;
			Output << '\n' ;
		}
	}
	
	Output.clear() ;
	Output.close() ;
	
	
	delete [] sCell ;
	delete [] NodeCell ;
}

//==============================================================================================================

void OutputCellTec( DSMC_DOMAIN *pDomain , DSMC_NODE *Node , DSMC_CELL *Cell , CELLMAPPING *pMapping ){

	if ( pDomain->Dimension == 2 ){
		OutputCellTec2D( pDomain , Node , Cell , pMapping ) ;
	}else if ( pDomain->Dimension == 3 ){
		OutputCellTec3D( pDomain , Node , Cell , pMapping ) ;
	}
	
}

//==============================================================================================================

void OutputCellTec2D( DSMC_DOMAIN *pDomain , DSMC_NODE *Node , DSMC_CELL *Cell , CELLMAPPING *pMapping ){
	ofstream	Output ;
	int		NodeNum , CellNum , OutNum ;
	
	NodeNum	= pDomain->NodeNum ;
	CellNum	= pDomain->CellNum ;
	OutNum	= 0 ;

	Output.open( "Cell-tec.dat" , ios::out | ios::trunc ) ;

	if ( !Output.fail() ){
		Output << "TITLE = \"Finite volume dataset\"\n" ;
		Output << "VARIABLES = \"X\", \"Y\", \"MFP\", \"MCT\", \"Timestep\", \"Weighting\", \"NewTimestep\", \"NewWeighting\", \"SubcellNum\", \"SelectNum\", \"CollRatio\"\n" ;
		Output << "ZONE N=" << NodeNum << ", E=" << CellNum << ", DATAPACKING=BLOCK, ZONETYPE=FEQUADRILATERAL\n" ;
		Output << "VARLOCATION=([1-2]=NODAL, [3-11]=CELLCENTERED)\n" ;

		
		for ( int i=0 ; i<NodeNum ; i++ ){
			OutNum++ ;
			Output << setw(18) << Node[i].XCoord ;
			if ( OutNum%6 == 0 ) Output << '\n' ;
		}
		
		for ( int i=0 ; i<NodeNum ; i++ ){
			OutNum++ ;
			Output << setw(18) << Node[i].YCoord ;
			if ( OutNum%6 == 0 ) Output << '\n' ;
		}

		for ( int i=0 ; i<CellNum ; i++ ){
			OutNum++ ;
			Output << setw(18) << Cell[i].MeanFreePath ;
			if ( OutNum%6 == 0 ) Output << '\n' ;
		}
		
		for ( int i=0 ; i<CellNum ; i++ ){
			OutNum++ ;
			Output << setw(18) << Cell[i].MeanCollisionTime ;
			if ( OutNum%6 == 0 ) Output << '\n' ;
		}
		
		for ( int i=0 ; i<CellNum ; i++ ){
			OutNum++ ;
			Output << setw(18) << Cell[i].Timestep ;
			if ( OutNum%6 == 0 ) Output << '\n' ;
		}
		
		for ( int i=0 ; i<CellNum ; i++ ){
			OutNum++ ;
			Output << setw(18) << Cell[i].Weighting ;
			if ( OutNum%6 == 0 ) Output << '\n' ;
		}
		
		for ( int i=0 ; i<CellNum ; i++ ){
			OutNum++ ;
			Output << setw(18) << Cell[i].InitTimestep ;
			if ( OutNum%6 == 0 ) Output << '\n' ;
		}
		
		for ( int i=0 ; i<CellNum ; i++ ){
			OutNum++ ;
			Output << setw(18) << Cell[i].InitWeighting ;
			if ( OutNum%6 == 0 ) Output << '\n' ;
		}
		
		for ( int i=0 ; i<CellNum ; i++ ){
			OutNum++ ;
			Output << setw(18) << Cell[i].SubcellNum ;
			if ( OutNum%6 == 0 ) Output << '\n' ;
		}
		
		for ( int i=0 ; i<CellNum ; i++ ){
			OutNum++ ;
			Output << setw(18) << Cell[i].SelectNum ;
			if ( OutNum%6 == 0 ) Output << '\n' ;
		}
		
		for ( int i=0 ; i<CellNum ; i++ ){
			OutNum++ ;
			Output << setw(18) << Cell[i].CollRatio ;
			if ( OutNum%6 == 0 ) Output << '\n' ;
		}


		Output << "\n\n" ;
		for ( int i=0 ; i<CellNum ; i++ ){
			for ( int j=0 ; j<pMapping->NodeNum[Cell[i].Type] ; j++ )
				Output << setw(18) << (Cell[i].Node[j]+1) ;
			if ( Cell[i].Type == 4 ) Output << setw(18) << (Cell[i].Node[2]+1) ;
			Output << '\n' ;
		}
	}else{
		cout << "Fail to open Cell-tec.dat \n" ;
		exit(1) ;
	}

	Output.clear() ;
	Output.close() ;	
}


void OutputCellTec3D( DSMC_DOMAIN *pDomain , DSMC_NODE *Node , DSMC_CELL *Cell , CELLMAPPING *pMapping ){
	ofstream	Output ;
	int		NodeNum , CellNum , OutNum ;
	
	NodeNum	= pDomain->NodeNum ;
	CellNum	= pDomain->CellNum ;
	OutNum	= 0 ;

	Output.open( "Cell-tec.dat" , ios::out | ios::trunc ) ;

	if ( !Output.fail() ){
		Output << "TITLE = \"Finite volume dataset\"\n" ;
		Output << "VARIABLES = \"X\", \"Y\", \"Z\", \"MFP\", \"MCT\", \"Temp\", \"Speed\", \"Timestep\", \"Weighting\", \"T_Ratio\", \"W_Ratio\", \"SubcellNum\", \"SelectNum\", \"CollRatio\"\n" ;
		Output << "ZONE N=" << NodeNum << ", E=" << CellNum << ", DATAPACKING=BLOCK, ZONETYPE=FEBRICK\n" ;
		Output << "VARLOCATION=([1-3]=NODAL, [4-14]=CELLCENTERED)\n" ;

		
		for ( int i=0 ; i<NodeNum ; i++ ){
			OutNum++ ;
			Output << setw(18) << Node[i].XCoord ;
			if ( OutNum%6 == 0 ) Output << '\n' ;
		}
		
		for ( int i=0 ; i<NodeNum ; i++ ){
			OutNum++ ;
			Output << setw(18) << Node[i].YCoord ;
			if ( OutNum%6 == 0 ) Output << '\n' ;
		}
		
		for ( int i=0 ; i<NodeNum ; i++ ){
			OutNum++ ;
			Output << setw(18) << Node[i].ZCoord ;
			if ( OutNum%6 == 0 ) Output << '\n' ;
		}
		
		for ( int i=0 ; i<CellNum ; i++ ){
			OutNum++ ;
			Output << setw(18) << Cell[i].MeanFreePath ;
			if ( OutNum%6 == 0 ) Output << '\n' ;
		}
		
		for ( int i=0 ; i<CellNum ; i++ ){
			OutNum++ ;
			Output << setw(18) << Cell[i].MeanCollisionTime ;
			if ( OutNum%6 == 0 ) Output << '\n' ;
		}
		
		for ( int i=0 ; i<CellNum ; i++ ){
			OutNum++ ;
			Output << setw(18) << Cell[i].Temp ;
			if ( OutNum%6 == 0 ) Output << '\n' ;
		}
		
		for ( int i=0 ; i<CellNum ; i++ ){
			OutNum++ ;
			Output << setw(18) << Cell[i].Speed ;
			if ( OutNum%6 == 0 ) Output << '\n' ;
		}
		
		for ( int i=0 ; i<CellNum ; i++ ){
			OutNum++ ;
			Output << setw(18) << Cell[i].Timestep ;
			if ( OutNum%6 == 0 ) Output << '\n' ;
		}
		
		for ( int i=0 ; i<CellNum ; i++ ){
			OutNum++ ;
			Output << setw(18) << Cell[i].Weighting ;
			if ( OutNum%6 == 0 ) Output << '\n' ;
		}
		
		for ( int i=0 ; i<CellNum ; i++ ){
			OutNum++ ;
			Output << setw(18) << Cell[i].InitTimestep ;
			if ( OutNum%6 == 0 ) Output << '\n' ;
		}
		
		for ( int i=0 ; i<CellNum ; i++ ){
			OutNum++ ;
			Output << setw(18) << Cell[i].InitWeighting ;
			if ( OutNum%6 == 0 ) Output << '\n' ;
		}
		
		for ( int i=0 ; i<CellNum ; i++ ){
			OutNum++ ;
			Output << setw(18) << Cell[i].SubcellNum ;
			if ( OutNum%6 == 0 ) Output << '\n' ;
		}
		
		for ( int i=0 ; i<CellNum ; i++ ){
			OutNum++ ;
			Output << setw(18) << Cell[i].SelectNum ;
			if ( OutNum%6 == 0 ) Output << '\n' ;
		}
		
		for ( int i=0 ; i<CellNum ; i++ ){
			OutNum++ ;
			Output << setw(18) << Cell[i].CollRatio ;
			if ( OutNum%6 == 0 ) Output << '\n' ;
		}


		Output << "\n\n" ;
		for ( int i=0 ; i<CellNum ; i++ ){
			if ( Cell[i].Type == 0 ){
				Output << setw(18) << Cell[i].Node[0]+1 ;
				Output << setw(18) << Cell[i].Node[1]+1 ;
				Output << setw(18) << Cell[i].Node[2]+1 ;
				Output << setw(18) << Cell[i].Node[2]+1 ;
				Output << setw(18) << Cell[i].Node[3]+1 ;
				Output << setw(18) << Cell[i].Node[3]+1 ;
				Output << setw(18) << Cell[i].Node[3]+1 ;
				Output << setw(18) << Cell[i].Node[3]+1 ;

			// Hex. Cell
			}else if ( Cell[i].Type == 1 ){
				for ( int j=0 ; j<pMapping->NodeNum[1] ; j++ )
					Output << setw(18) << Cell[i].Node[j]+1 ;

			// Pyramid Cell
			}else if ( Cell[i].Type == 2 ){
				Output << setw(18) << Cell[i].Node[0]+1 ;
				Output << setw(18) << Cell[i].Node[1]+1 ;
				Output << setw(18) << Cell[i].Node[2]+1 ;
				Output << setw(18) << Cell[i].Node[3]+1 ;
				Output << setw(18) << Cell[i].Node[4]+1 ;
				Output << setw(18) << Cell[i].Node[4]+1 ;
				Output << setw(18) << Cell[i].Node[4]+1 ;
				Output << setw(18) << Cell[i].Node[4]+1 ;


			// Prism Cell
			}else if ( Cell[i].Type == 3 ){
				Output << setw(18) << Cell[i].Node[0]+1 ;
				Output << setw(18) << Cell[i].Node[3]+1 ;
				Output << setw(18) << Cell[i].Node[4]+1 ;
				Output << setw(18) << Cell[i].Node[4]+1 ;
				Output << setw(18) << Cell[i].Node[1]+1 ;
				Output << setw(18) << Cell[i].Node[2]+1 ;
				Output << setw(18) << Cell[i].Node[5]+1 ;
				Output << setw(18) << Cell[i].Node[5]+1 ;

			}
			
			Output << '\n' ;
		}
	}else{
		cout << "Fail to open Cell-tec.dat \n" ;
		exit(1) ;
	}

	Output.clear() ;
	Output.close() ;
}

//==============================================================================================================
//==============================================================================================================

void ReadResult( DSMC_DOMAIN *pDomain , DSMC_CELL *Cell , DSMC_RESULT *pResult , string Filename ){
	ifstream	Input ;
	string		buffer ;
	int		CellNum , CellNo ;
	
	
	CellNum	= pDomain->CellNum ;


	Input.open( Filename.c_str() , ios::in ) ;
		
	if ( !Input ){
		cout << "Fail to open " << Filename << '\n' ;	
		exit(1) ;
	}
	
	getline( Input , buffer ) ;
	for ( int i=0 ; i<CellNum ; i++ ){
		Input >> CellNo ;
		Input >> Cell[CellNo].XCenter >> Cell[CellNo].YCenter >> Cell[CellNo].ZCenter >> pResult->Density[CellNo]
			>> pResult->NumDensity[CellNo] >> pResult->XVel[CellNo] >> pResult->YVel[CellNo] >> pResult->ZVel[CellNo] 
			>> pResult->Temp[CellNo] >> pResult->TransTemp[CellNo] >> pResult->RotTemp[CellNo] >> pResult->VibTemp[CellNo] 
			>> pResult->AveParticleNum[CellNo] >> pResult->MeanCollSpacingMeanFreePath[CellNo] >> Cell[CellNo].ProcessorNo ;
	}
	
	Input.clear() ;
	Input.close() ;
}

//==============================================================================================================

void ReadResultSpecies( DSMC_DOMAIN *pDomain , DSMC_CELL *Cell , DSMC_RESULT *pResult , string Filename ){
	ifstream	Input ;
	string		buffer ;
	int		CellNum , CellNo ;
	
	
	CellNum	= pDomain->CellNum ;


	Input.open( Filename.c_str() , ios::in ) ;
		
	if ( !Input ){
		cout << "Fail to open " << Filename << '\n' ;	
		exit(1) ;
	}
	
	getline( Input , buffer ) ;
	for ( int i=0 ; i<CellNum ; i++ ){
		Input >> CellNo ;
		Input >> Cell[CellNo].XCenter >> Cell[CellNo].YCenter >> Cell[CellNo].ZCenter >> pResult->Density[CellNo]
			>> pResult->NumDensity[CellNo] >> pResult->XVel[CellNo] >> pResult->YVel[CellNo] >> pResult->ZVel[CellNo] 
			>> pResult->Temp[CellNo] >> pResult->TransTemp[CellNo] >> pResult->RotTemp[CellNo] >> pResult->VibTemp[CellNo] 
			>> pResult->TransXTemp[CellNo] >> pResult->TransYTemp[CellNo] >> pResult->TransZTemp[CellNo]>> pResult->AveParticleNum[CellNo] ;
	}
	
	Input.clear() ;
	Input.close() ;
}

//==============================================================================================================

void ReadResultSurface( DSMC_DOMAIN *pDomain , DSMC_POST_SURFACE *Surface , string Filename ){
	ifstream	Input ;
	string		buffer ;
	int		WallFaceNum ;
	
	WallFaceNum	= pDomain->WallFaceNum ;
	
	Input.open( Filename.c_str() , ios::in ) ;
		
	if ( !Input ){
		cout << "Fail to open " << Filename << '\n' ;	
		exit(1) ;
	}
	
	getline( Input , buffer ) ;
	for ( int i=0 ; i<WallFaceNum ; i++ ){
		Surface[i].Id	= i ;
		Input >> Surface[i].CellNo >> Surface[i].Center[0] >> Surface[i].Center[1] >> Surface[i].Center[2] 
			>> Surface[i].SurfacePro[0] >> Surface[i].SurfacePro[1] >> Surface[i].SurfacePro[2] >> Surface[i].SurfacePro[3] 
			>> Surface[i].SurfacePro[4] >> Surface[i].SurfacePro[5] >> Surface[i].SurfacePro[6] >> Surface[i].SurfacePro[7] 
			>> Surface[i].SurfacePro[8] >> Surface[i].SurfacePro[9] >> Surface[i].SurfacePro[10]>> Surface[i].SurfacePro[11] 
			>> Surface[i].SurfacePro[12] >> Surface[i].SurfacePro[13] >> Surface[i].SurfacePro[14] >> Surface[i].SurfacePro[15]
			>> Surface[i].SurfacePro[16] >> Surface[i].Area >> Surface[i].SimTime >> Surface[i].WallNo >> Surface[i].FaceNo ;
			
		// Debug.
		//Surface[i].Dump(i) ;
		//getchar() ;
	}
	
	Input.clear() ;
	Input.close() ;
}

//==============================================================================================================

void ReadCellInformation( DSMC_DOMAIN *pDomain , DSMC_CELL *Cell , string Filename ){
	int			CellNo , CellNum ;
	double			Buffer ;
	ifstream		Input ;
	string			buffer ;


	CellNum	= pDomain->CellNum ;


	Input.open( Filename.c_str() , ios::in ) ;

	if ( !Input ){
		cout << "Failed to open file " << Filename << endl ;
		exit(1) ;
	}


	getline( Input , buffer ) ;
	for ( int i=0 ; i<CellNum ; i++ ){
		Input >> CellNo ;
		Input >> Cell[CellNo].MeanFreePath >> Cell[CellNo].MeanCollisionTime >> Cell[CellNo].Timestep >> Cell[CellNo].Weighting 
			>> Cell[CellNo].InitTimestep >> Cell[CellNo].InitWeighting >> Buffer >> Buffer >> Cell[CellNo].SubcellNum 
			>> Cell[CellNo].SelectNum >> Cell[CellNo].CollRatio ;
	}
	
	
	// To close the file.
	Input.clear() ;
	Input.close() ;
}

//==============================================================================================================
//==============================================================================================================

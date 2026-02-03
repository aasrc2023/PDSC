#include <cstdlib>
#include <fstream>
#include <string>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <metis.h>
#include "dsmc_class_pre.h"
#include "dsmc_preprocessor.h"

using namespace std ;


//==============================================================================================================

void InitOption( DSMC_PRE_OPTION *Option , int argc , char* argv[] ){
	int		OptionNum , FindNo ;
	string		*option , word ;
	
	
	OptionNum	= argc - 2 ;
	
	if ( argc > 2 ){
		option	= new string[argc-2] ;
		
		for ( int i=0 ; i<OptionNum ; i++ )
			option[i]	= argv[i+2] ;
		
		
		for ( int i=0 ; i<OptionNum ; i++ ){
			FindNo	= option[i].find('=') ;
			word	= option[i].substr( 0 , FindNo ) ;
			
			
			if ( word == "np" ){
				word	= option[i].substr( FindNo+1 ) ;
				
				Option->ProcessorNum	= atoi( word.c_str() ) ;
				
			}else if ( word == "scale" ){
				word	= option[i].substr( FindNo+1 ) ;

				Option->OpenScale	= 1 ;
				Option->Scale		= atof( word.c_str() ) ;
			
				//continue ;
			}else if ( word == "inlet" ){
				word	= option[i].substr( FindNo+1 ) ;

				Option->OpenInputInlet	= 1 ;
				Option->InletFileName	= word ;
			
				//continue ;
			}else if ( word == "drd" ){
				word	= option[i].substr( FindNo+1 ) ;
				
				Option->OpenDomainRedecomposition	= 1 ;
				Option->InitResultFileName		= word ;
				Option->InitMeshFileName		= "Mesh_Before.txt" ;
				
			}else if ( word == "cell" ){
				word	= option[i].substr( FindNo+1 ) ;

				Option->OpenCell	= 1 ;
				Option->CellFileName	= word ;
				Option->InitMeshFileName= "Mesh_Before.txt" ;
			}
		}		

		delete [] option ;
	}
}

//==============================================================================================================

void GetArraySizeFieldview( int *pNodeNum , int *pCellNum , int *pBCTypeNum , int *pBCFaceNum , string Filename ){
	string		buffer ;
	ifstream	Input ;

	(*pNodeNum)	= 0 ;
	(*pCellNum)	= 0 ;
	(*pBCTypeNum)	= 0 ;
	(*pBCFaceNum)	= 0 ;


	Input.open( Filename.c_str() , ios::in ) ;
	
	if ( Input.fail() ){
 		cout << "Fail to open file " << Filename << '\n';
		exit(1) ;
	}
	
	while ( getline( Input , buffer ) ){
		if ( FindString(buffer , "Boundary Table") ){
			getline( Input , buffer ) ;
			(*pBCTypeNum) = atoi(buffer.c_str()) ;
		}

		if ( FindString(buffer , "Nodes") ){
			getline( Input , buffer ) ;
			(*pNodeNum) = atoi(buffer.c_str()) ;
		}

		if ( FindString(buffer , "Boundary Faces") ){
			getline( Input , buffer ) ;
			(*pBCFaceNum) = atoi(buffer.c_str()) ;
		}

		if ( FindString(buffer , "Elements") ){
			while ( getline( Input , buffer ) ){
				if ( !FindString(buffer , "Variables") )
					(*pCellNum)++ ;
				else
					break ;
			}
		}
	}

	Input.clear() ;
	Input.close() ;	
}

//==============================================================================================================

bool FindString( string Word1 , string Word2 ){
	bool	a = true ;
	for ( int i=0 ; i<Word2.size() ; i++ ){
		if ( Word1[i] != Word2[i] ){
			a = false ;
			break ;
		}
	}

	return	a ;
}

//==============================================================================================================

void ReadGridFieldview( int		NodeNum , 
			int		CellNum , 
			int		BCTypeNum , 
			int		BCFaceNum , 
			DSMC_PRE_NODE	*Node , 
			DSMC_PRE_CELL	*Cell , 
			DSMC_PRE_BCTYPE	*BCType , 
			DSMC_PRE_BCFACE	*BCFace ,
			CELLMAPPING	*pMapping , 
			string		Filename ){

	int		NodeNo ;
	string		buffer ;
	ifstream	Input ;

	Input.open( Filename.c_str() , ios::in ) ;
	
	if ( Input.fail() ){
 		cout << "Fail to open file " << Filename << '\n';
		exit(1) ;
	}
	
	while ( getline( Input , buffer ) ){
		if ( FindString(buffer , "Boundary Table") ){
			getline( Input , buffer ) ;
			for ( int i=0 ; i<BCTypeNum ; i++ )
				Input >> buffer >> BCType[i].TypeName ;
		}


		if ( FindString(buffer , "Nodes") ){
			getline( Input , buffer ) ;
			for ( int i=0 ; i<NodeNum ; i++ ){
				Node[i].Id	=	i ;
				Input >> Node[i].XCoord >> Node[i].YCoord >> Node[i].ZCoord ;
			}
		}


		if ( FindString(buffer , "Boundary Faces") ){
			getline( Input , buffer ) ;
			for ( int i=0 ; i<BCFaceNum ; i++ ){
				Input >> BCFace[i].Id >> BCFace[i].NodeNum ;
				for ( int j=0 ; j<BCFace[i].NodeNum ; j++ ){
					Input >> BCFace[i].Node[j] ;
					BCFace[i].Node[j]-- ;
				}
				BCFace[i].Id-- ;
			}
		}


		if ( FindString(buffer , "Elements") ){
			for ( int i=0 ; i<CellNum ; i++ ){
				Input >> Cell[i].Type >> NodeNo ;
				
				Cell[i].Id	=	i ;
				
				if ( Cell[i].Type == 1 )
					Cell[i].Type = 0 ;
				else if ( Cell[i].Type == 2 )
					Cell[i].Type = 1 ;
				else if ( Cell[i].Type == 3 )
					Cell[i].Type = 3 ;
				else if ( Cell[i].Type == 4 )
					Cell[i].Type = 2 ;
				else if ( Cell[i].Type == 5 )
					Cell[i].Type = 0 ;
				else if ( Cell[i].Type == 6 )
					Cell[i].Type = 0 ;


				for ( int j=0 ; j<pMapping->NodeNum[Cell[i].Type] ; j++ ){
					Input >> NodeNo ;
					Cell[i].pNode[j] = &Node[NodeNo-1] ;
				}

				// For Tetra Cell Type
				if ( Cell[i].Type == 0 ){

					NodeNo = Cell[i].pNode[1]->Id ;
					Cell[i].pNode[1] = &Node[Cell[i].pNode[2]->Id] ;
					Cell[i].pNode[2] = &Node[NodeNo] ;

				// For Hex. Cell Type
				}else if ( Cell[i].Type == 1 ){

					NodeNo = Cell[i].pNode[2]->Id ;
					Cell[i].pNode[2] = &Node[Cell[i].pNode[3]->Id] ;
					Cell[i].pNode[3] = &Node[NodeNo] ;

					NodeNo = Cell[i].pNode[6]->Id ;
					Cell[i].pNode[6] = &Node[Cell[i].pNode[7]->Id] ;
					Cell[i].pNode[7] = &Node[NodeNo] ;
					
				// For Pyramid Cell Type
				//}else if ( Cell[i].Type == 2 ){


				// For Prism Cell Type
				}else if ( Cell[i].Type == 3 ){

					NodeNo = Cell[i].pNode[4]->Id ;
					Cell[i].pNode[4] = &Node[Cell[i].pNode[5]->Id] ;
					Cell[i].pNode[5] = &Node[NodeNo] ;
				}
			}
			break ;
		}
	}

	Input.clear() ;
	Input.close() ;
}

//==============================================================================================================

void ReadBoundaryCondition( int BCTypeNum , DSMC_PRE_BCTYPE *BCType ){
	ifstream	Input ;
	string		Buffer ;

	Input.open( "Boundary.txt" , ios::in ) ;
	
	if ( Input.fail() ){
 		cout << "Fail to open file Boundary.txt\n" ;
		exit(1) ;
	}


	getline( Input , Buffer ) ;
	getline( Input , Buffer ) ;
	for ( int i=0 ; i<BCTypeNum ; i++ ){
		Input >> Buffer ;

		for ( int j=0 ; j<BCTypeNum ; j++ ){
			if ( FindString( BCType[j].TypeName , Buffer ) ){
				Input >> BCType[j].Type ;

				if ( BCType[j].Type == -3 ){
					Input >> BCType[j].XVel 
					      >> BCType[j].YVel
					      >> BCType[j].ZVel
					      >> BCType[j].NumDen
					      >> BCType[j].Temp 
					      >> BCType[j].CosineLawCoef ;
				}
				
				break ;
			}
		}
	}

	Input.clear() ;
	Input.close() ;
}

//==============================================================================================================

bool Check2DMesh( int BCTypeNum , DSMC_PRE_BCTYPE *BCType ){
	bool		IS2D = false ;

	for ( int i=0 ; i<BCTypeNum ; i++ ){
		if ( BCType[i].Type == 1 ){
			IS2D	=	true ;
			break ;
		}
	}

	return	IS2D ;
}

//==============================================================================================================

int Transfer3DTo2D(	int		NodeNum , 
			int		CellNum , 
			int		BCFaceNum ,
			DSMC_PRE_NODE	*Node , 
			DSMC_PRE_CELL	*Cell , 
			DSMC_PRE_BCFACE	*BCFace , 
			DSMC_PRE_BCTYPE	*BCType ,  
			CELLMAPPING	*pMapping ){

	int			NodeNo , CellNo , NewNodeNum , FaceNo , Type ;
	DSMC_PRE_NODECELL	*NodeCell ;
	
	NodeCell	=	new DSMC_PRE_NODECELL[NodeNum] ;
	

	for ( int i=0 ; i<BCFaceNum ; i++ ){
		// Find out the normal face.
		Type = BCType[BCFace[i].Id].Type ;
		if ( Type == 1 ){
			for ( int j=0 ; j<CellNum ; j++ ){
				// Set type of j-th cell
				if ( BCFace[i+j].NodeNum == 3 ){
					// 2-d triangular cell.
					Cell[j].Type = 4 ;
				}else if ( BCFace[i+j].NodeNum == 4 ){
					// 2-d quadrilateral cell.
					Cell[j].Type = 5 ;
				}
				

				// Set Node of j-th Cell
				for ( int k=0 ; k<BCFace[i+j].NodeNum ; k++ ){
					NodeNo		= BCFace[i+j].Node[k] ;
					Cell[j].pNode[k]= &Node[NodeNo] ;
				}
			}

			break ;
		}
	}
	
	
	// Node counterclockwise.
	for ( int i=0 ; i<CellNum ; i++ ) NodeCCW( Node , &Cell[i] , pMapping ) ;
	

	CreateNodeCellRelation( CellNum , Cell , NodeCell , pMapping ) ;
	CreateNodeFaceRelation( BCFaceNum , BCFace , BCType , NodeCell , pMapping ) ;


	NewNodeNum = 0 ;
	for ( int i=0 ; i<NodeNum ; i++ ){
		if ( NodeCell[i].CellRelation.size() > 0 ){
			Node[i].Id	=	NewNodeNum ;
			Node[NewNodeNum]=	Node[i] ;

			// Update the Node for cell.
			for ( int j=0 ; j<NodeCell[i].CellRelation.size() ; j++ ){
				CellNo	=	NodeCell[i].CellRelation[j] ;

				for ( int k=0 ; k<pMapping->NodeNum[Cell[CellNo].Type] ; k++ ){		
					if ( Cell[CellNo].pNode[k] == &Node[i] )
						Cell[CellNo].pNode[k] = &Node[NewNodeNum] ;
				}
			}
			
			
			if ( NodeCell[i].FaceRelation.size() > 0 ){
				for ( int j=0 ; j<NodeCell[i].FaceRelation.size() ; j++ ){
					FaceNo	=	NodeCell[i].FaceRelation[j] ;
					for ( int k=0 ; k<BCFace[FaceNo].NodeNum ; k++ ){
						if ( BCFace[FaceNo].Node[k] == i ){
							BCFace[FaceNo].NewNode[BCFace[FaceNo].CountNum] = NewNodeNum ;
							BCFace[FaceNo].CountNum++ ;
							break ;
						}
					}
				}
			}
			

			NewNodeNum++ ;
		}
	}
	
	
	// Update node of each face.
	for ( int i=0 ; i<BCFaceNum ; i++ ){
		if ( BCType[BCFace[i].Id].Type < 0 ){
			BCFace[i].NodeNum	= 2 ;
			for ( int j=0 ; j<BCFace[i].NodeNum ; j++ )
				BCFace[i].Node[j]	=	BCFace[i].NewNode[j] ;
		}	
	}
	
	
	delete [] NodeCell ;
	
	return	NewNodeNum ;
}

//==============================================================================================================

void NodeCCW( DSMC_PRE_NODE *Node , DSMC_PRE_CELL *_Cell , CELLMAPPING *pMapping ){
	int		Type , *NodeNo , NodeNoMaxXCoord ;
	double		x0 , y0	, x1 , y1 , *Angle , a , b , c ;
	DSMC_PRE_NODE	*pNodeBuffer ;
	
	
	Type	=	_Cell->Type ;
	Angle	=	new double[pMapping->NodeNum[Type]] ;
	NodeNo	=	new int[pMapping->NodeNum[Type]] ;
	
	Angle[0]=	0. ;
	NodeNo[0]=	0 ;

	NodeNoMaxXCoord	= 0 ;
	
	for ( int i=1 ; i<pMapping->NodeNum[Type] ; i++ ){
		if ( _Cell->pNode[i]->XCoord > _Cell->pNode[NodeNoMaxXCoord]->XCoord ){
			NodeNoMaxXCoord	= i ;
		}
	}
	pNodeBuffer			= _Cell->pNode[0] ;
	_Cell->pNode[0]			= _Cell->pNode[NodeNoMaxXCoord] ;
	_Cell->pNode[NodeNoMaxXCoord]	= pNodeBuffer ;
	
		
	x0	=	_Cell->pNode[0]->XCoord ;
	y0	=	_Cell->pNode[0]->YCoord ;
	for ( int i=1 ; i<pMapping->NodeNum[Type] ; i++ ){
		x1	=	_Cell->pNode[i]->XCoord ;
		y1	=	_Cell->pNode[i]->YCoord ;
		a	=	x1 - x0 ;
		b	=	y1 - y0 ;
		c	=	sqrt(a*a + b*b) ;
		
		Angle[i]=	CalculateAngleDeg( a , b , c ) ;
		NodeNo[i]=	_Cell->pNode[i]->Id ;
	}
	
	
	// Sort the angle from the smallest to the biggest.
	BubbleSort( &Angle[1] , &NodeNo[1] , (pMapping->NodeNum[Type]-1) ) ;

	for ( int i=1 ; i<pMapping->NodeNum[Type] ; i++ )
		_Cell->pNode[i]	=	&Node[NodeNo[i]] ;

	
	delete [] Angle ;
	delete [] NodeNo ;
}

//==============================================================================================================

template <class TYPE>
void BubbleSort( TYPE *a , int *Index , int Num ){
	for ( int i=0 ; i<(Num-1) ; i++ ){
		for ( int j=0 ; j<(Num-(1+i)) ; j++ ){
			if ( a[j] > a[j+1] ){
				Swap( &a[j] , &a[j+1] ) ;
				Swap( &Index[j] , &Index[j+1] ) ;
			}
		}
	}
}

//==============================================================================================================

template <class TYPE>
void Swap( TYPE *a , TYPE *b ){
	TYPE	buffer ;
	buffer	=	(*a) ;
	(*a)	=	(*b) ;
	(*b)	=	buffer ;
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

void CreateNodeCellRelation( int CellNum , DSMC_PRE_CELL *Cell , DSMC_PRE_NODECELL *NodeCell , CELLMAPPING *pMapping ){
	int		NodeNo , Type ;

	for ( int i=0 ; i<CellNum ; i++ ){
		Type	=	Cell[i].Type ;
		for ( int j=0 ; j<pMapping->NodeNum[Type] ; j++ ){
			NodeNo	=	Cell[i].pNode[j]->Id ;
			
			NodeCell[NodeNo].CellRelation.push_back(i) ;
		}
	}
}

//==============================================================================================================

void CreateNodeFaceRelation( int BCFaceNum , DSMC_PRE_BCFACE *BCFace , DSMC_PRE_BCTYPE *BCType , DSMC_PRE_NODECELL *NodeCell , CELLMAPPING *pMapping ){
	int		NodeNo , Type ;

	for ( int i=0 ; i<BCFaceNum ; i++ ){
		Type	=	BCType[BCFace[i].Id].Type ;
		
		// To check the BCFace is the boundary face.
		if ( Type < 0 ){
			for ( int j=0 ; j<BCFace[i].NodeNum ; j++ ){
				NodeNo	=	BCFace[i].Node[j] ;
			
				NodeCell[NodeNo].FaceRelation.push_back(i) ;
			}
		}
	}
}

//==============================================================================================================

void DomainPartition(		int			CellNum , 
				int			NodeNum , 
				DSMC_PRE_NODE		*Node ,
				DSMC_PRE_CELL		*Cell , 
				DSMC_PRE_PROCESSOR	*Processor ,
				DSMC_PRE_OPTION		*Option ,
				CELLMAPPING		*pMapping ){
	
	DSMC_PRE_GRAPHCELL	*GraphCell ;
	int			Num , MetisState ;
	idx_t			*MetisCellProcessorNo , *MetisXadj , *MetisAdjncy , *MetisVwgt ; 
	idx_t			MetisCellNum , MetisProcessorNum , MetisObjval , MetisNcon ;

	
	GraphCell	= new DSMC_PRE_GRAPHCELL[CellNum] ;
	
	// Create Graph and return 2x edge number.
	Num	= CreateGraph( CellNum , NodeNum ,  Cell , GraphCell , pMapping ) ;


	MetisState		= 0 ;
	MetisCellNum		= CellNum ;
	MetisProcessorNum	= Option->ProcessorNum ;
	MetisObjval		= 0 ;
	MetisNcon		= 1 ;


	MetisCellProcessorNo	= new idx_t[CellNum] ;
	MetisXadj		= new idx_t[CellNum+1] ;
	MetisAdjncy		= new idx_t[Num] ;
	MetisVwgt		= new idx_t[CellNum] ; 
	
	
	// Create CSR format for graph data.
	MetisXadj[0]	= 0 ;
	for ( int i=0 ; i<CellNum ; i++ ){
		MetisXadj[i+1]	=	MetisXadj[i] ;

		for ( int j=0 ; j<GraphCell[i].Num ; j++ ){
			MetisAdjncy[MetisXadj[i+1]]	= GraphCell[i].Neighbor[j] ;
			MetisXadj[i+1]++ ;
		}
	}
	
	
	// Meitis API for graph partition.
	if ( Option->OpenDomainRedecomposition == 0 ){
		MetisState	= METIS_PartGraphRecursive( &MetisCellNum , &MetisNcon , MetisXadj , MetisAdjncy , NULL , NULL , 
							NULL , &MetisProcessorNum , NULL , NULL , NULL , &MetisObjval , MetisCellProcessorNo ) ;
	}else if ( Option->OpenDomainRedecomposition == 1 ){
		for ( int i=0 ; i<CellNum ; i++ ){
			MetisVwgt[i]	= Cell[i].Weight ;
			//cout << MetisVwgt[i] << '\n' ;
		}
			
		MetisState	= METIS_PartGraphRecursive( &MetisCellNum , &MetisNcon , MetisXadj , MetisAdjncy , MetisVwgt , NULL , 
							NULL , &MetisProcessorNum , NULL , NULL , NULL , &MetisObjval , MetisCellProcessorNo ) ;
	}
							
	
	ResortCellNo( Cell , Processor , CellNum , Option->ProcessorNum , MetisCellProcessorNo ) ;
	
	
	// Debug.
	//for ( int i=0 ; i<Option->ProcessorNum ; i++ )
	//	cout << "CellNum: " << Processor[i].Cell.size() << '\n' ;
	
	
	delete [] GraphCell ;
	delete [] MetisCellProcessorNo ;
	delete [] MetisXadj ;
	delete [] MetisAdjncy ;
	delete [] MetisVwgt ;
}

//==============================================================================================================

int CreateGraph( int			CellNum , 
		 int			NodeNum , 
		 DSMC_PRE_CELL		*Cell ,
		 DSMC_PRE_GRAPHCELL	*GraphCell ,
		 CELLMAPPING		*pMapping ){

	int			c1 , c2 , j2 , Type1 , Type2 , Num ;
	DSMC_PRE_NODECELL	*NodeCell ;

	Num	= 0 ;

	NodeCell	= new DSMC_PRE_NODECELL[NodeNum] ;

	// To create the relation between nodes and cells.
	CreateNodeCellRelation( CellNum , Cell , NodeCell , pMapping ) ;
	

	for ( int i=0 ; i<NodeNum ; i++ ){
		for ( int j1=0 ; j1<NodeCell[i].CellRelation.size() ; j1++ ){
			c1	= NodeCell[i].CellRelation[j1] ;

			Type1	= Cell[c1].Type ;
			for ( int k1=0 ; k1<pMapping->SurfaceNum[Type1] ; k1++ ){
				if ( GraphCell[c1].Check[k1] == -999 ){
					j2	= 0 ;
					
					while ( j2 < NodeCell[i].CellRelation.size() ){
						if ( j1 != j2  ){
							c2	= NodeCell[i].CellRelation[j2] ;

							Type2	= Cell[c2].Type ;
							for ( int k2=0 ; k2<pMapping->SurfaceNum[Type2] ; k2++ ){
								if ( CheckNeighbor( &Cell[c1] , &Cell[c2] , k1 , k2 , pMapping ) ){
									GraphCell[c1].Check[k1]	= 0 ;

									GraphCell[c1].Neighbor[GraphCell[c1].Num]	=	c2 ;
									GraphCell[c1].Num++ ;
									Num++ ;

									j2	=	NodeCell[i].CellRelation.size() ;
									break ;
								}
							}
						}

						j2++ ;
					}
				}
			}
		}
	}
	
	delete [] NodeCell ;

	return	Num ;
}

//==============================================================================================================

void ResortCellNo(	DSMC_PRE_CELL		*Cell , 
			DSMC_PRE_PROCESSOR	*Processor , 
			int			CellNum , 
			int			ProcessorNum , 
			idx_t			*MetisCellProcessorNo ){
				
	int		ProcessorNo , CellNo ;


	for ( int i=0 ; i<CellNum ; i++ ){
		ProcessorNo		= MetisCellProcessorNo[i] ;
		Cell[i].ProcessorNo	= ProcessorNo ;

		Processor[ProcessorNo].Cell.push_back( Cell[i] ) ;
	}
	
	
	CellNo	=	0 ;
	for ( int i=0 ; i<ProcessorNum ; i++ ){
		for ( int j=0 ; j<Processor[i].Cell.size() ; j++ ){
			Cell[CellNo]	=	Processor[i].Cell[j] ;
			CellNo++ ;
		}
	}
}

//==============================================================================================================

void CreateCellNeighbor( 	int		CellNum , 
				int		NodeNum , 
				int		BCFaceNum ,
				DSMC_PRE_NODE	*Node ,
				DSMC_PRE_CELL	*Cell , 
				DSMC_PRE_BCFACE	*BCFace ,
				DSMC_PRE_BCTYPE	*BCType ,
				CELLMAPPING	*pMapping ){

	int			c1 , c2 , j2 , Type1 , Type2 , NodeNo , FaceNo ;
	DSMC_PRE_NODECELL	*NodeCell ;

	NodeCell	= new DSMC_PRE_NODECELL[NodeNum] ;

	// To create the relation between nodes and cells.
	CreateNodeCellRelation( CellNum , Cell , NodeCell , pMapping ) ;


	for ( int i=0 ; i<NodeNum ; i++ ){
		for ( int j1=0 ; j1<NodeCell[i].CellRelation.size() ; j1++ ){
			c1	= NodeCell[i].CellRelation[j1] ;

			Type1	= Cell[c1].Type ;
			for ( int k1=0 ; k1<pMapping->SurfaceNum[Type1] ; k1++ ){
				if ( Cell[c1].Neighbor[k1] == -999 ){
					j2	= 0 ;
					while ( j2 < NodeCell[i].CellRelation.size() ){
						if ( j1 != j2  ){
							c2	= NodeCell[i].CellRelation[j2] ;

							Type2	= Cell[c2].Type ;
							for ( int k2=0 ; k2<pMapping->SurfaceNum[Type2] ; k2++ ){
								if ( CheckNeighbor( &Cell[c1] , &Cell[c2] , k1 , k2 , pMapping ) ){
									
									Cell[c1].Neighbor[k1]	= c2 ;
									Cell[c2].Neighbor[k2]	= c1 ;

									j2	= NodeCell[i].CellRelation.size() ;
									break ;
								}
							}
						}

						j2++ ;
					}
				}
			}
		}
	}

	// Debug.
	//cout << "before creating the relation between nodes and faces\n" ;

	// To create the relation between node and face.
	CreateNodeFaceRelation( BCFaceNum , BCFace , BCType , NodeCell , pMapping ) ;

	// Debug.
	//cout << "before creating the boundary face\n" ;

	for ( int i=0 ; i<CellNum ; i++ ){
		Type1	= Cell[i].Type ;
		for ( int j=0 ; j<pMapping->SurfaceNum[Type1] ; j++ ){
			if ( Cell[i].Neighbor[j] == -999 ){
				NodeNo	= pMapping->Node[Type1][j][0] ;
				NodeNo	= Cell[i].pNode[NodeNo]->Id ;

				for ( int k=0 ; k<NodeCell[NodeNo].FaceRelation.size() ; k++ ){
					FaceNo	= NodeCell[NodeNo].FaceRelation[k] ;

					if ( CheckNeighbor( &Cell[i] , &BCFace[FaceNo], j , pMapping ) ){
						Cell[i].Neighbor[j]	= BCType[BCFace[FaceNo].Id].Type ;

						BCType[BCFace[FaceNo].Id].FaceNum++ ;

						// Output the node of face as cells
						for ( int m=0 ; m<BCFace[FaceNo].NodeNum ; m++ ){
							int kk			= pMapping->Node[Type1][j][m] ;
							BCFace[FaceNo].Node[m]	= Cell[i].pNode[kk]->Id ;
						}
						break ;
					}
				}
			}
		}
	}

	delete [] NodeCell ;
}

//==============================================================================================================

bool CheckNeighbor( DSMC_PRE_CELL *_CellA , DSMC_PRE_CELL *_CellB , int SurfaceA , int SurfaceB , CELLMAPPING *pMapping ){
	bool		Check	=	false ;
	int		TypeA , TypeB , Num , NodeA , NodeB ;

	TypeA	=	_CellA->Type ;
	TypeB	=	_CellB->Type ;
	Num	=	0 ;

	for ( int i=0 ; i<pMapping->SurfaceNode[TypeA][SurfaceA] ; i++ ){
		NodeA	= pMapping->Node[TypeA][SurfaceA][i] ;
		NodeA	= _CellA->pNode[NodeA]->Id ;

		for ( int j=0 ; j<pMapping->SurfaceNode[TypeB][SurfaceB] ; j++ ){
			NodeB	= pMapping->Node[TypeB][SurfaceB][j] ;
			NodeB	= _CellB->pNode[NodeB]->Id ;

			if ( NodeA == NodeB ){
				Num++ ;
				break ;
			}
		}
	}


	if ( Num == pMapping->SurfaceNode[TypeA][SurfaceA] ) Check = true ;

	return	Check ;
}

//==============================================================================================================

bool CheckNeighbor( DSMC_PRE_CELL *_CellA , DSMC_PRE_BCFACE *_BCFace, int SurfaceA , CELLMAPPING *pMapping ){
	bool		Check	=	false ;
	int		TypeA , Num , NodeA , NodeB ;

	TypeA	=	_CellA->Type ;
	Num	=	0 ;

	for ( int i=0 ; i<pMapping->SurfaceNode[TypeA][SurfaceA] ; i++ ){
		NodeA	= pMapping->Node[TypeA][SurfaceA][i] ;
		NodeA	= _CellA->pNode[NodeA]->Id ;

		for ( int j=0 ; j<_BCFace->NodeNum ; j++ ){
			NodeB	= _BCFace->Node[j] ;

			if ( NodeA == NodeB ){
				Num++ ;
				break ;
			}
		}
	}


	if ( Num == pMapping->SurfaceNode[TypeA][SurfaceA] ) Check = true ;

	return	Check ;
}

//==============================================================================================================

void CreateProcessorNeighbor( int ProcessorNum , DSMC_PRE_CELL *Cell , DSMC_PRE_PROCESSOR *Processor , CELLMAPPING *pMapping ){
	int		CellNo , ProcessorNo ;
	int		**CheckProcessorNeighbor ;

	CellNo	= 0 ;

	// Allocate Memory.
	CheckProcessorNeighbor	= new int*[ProcessorNum] ;
	for ( int i=0 ; i<ProcessorNum ; i++ )
		CheckProcessorNeighbor[i]	= new int[ProcessorNum] ;


	for ( int i=0 ; i<ProcessorNum ; i++ ){
		for ( int j=0 ; j<ProcessorNum ; j++ )
			CheckProcessorNeighbor[i][j]	=	-999 ;
	}


	for ( int i=0 ; i<ProcessorNum ; i++ ){
		for ( int j=0 ; j<Processor[i].Cell.size() ; j++ ){
			for ( int k=0 ; k<pMapping->SurfaceNum[Cell[CellNo].Type] ; k++ ){
				if ( Cell[CellNo].Neighbor[k] < 0 ) continue ;

				ProcessorNo	= Cell[Cell[CellNo].Neighbor[k]].ProcessorNo ;

				if ( ProcessorNo != i && CheckProcessorNeighbor[i][ProcessorNo] == -999 ){

					Processor[i].Neighbor.push_back( ProcessorNo ) ;
					CheckProcessorNeighbor[i][ProcessorNo]	= 0 ;
				}
			}
			CellNo++ ;
		}
	}


	// Deallocate Memory.
	for ( int i=0 ; i<ProcessorNum ; i++ )
		delete [] CheckProcessorNeighbor[i] ;
	delete [] CheckProcessorNeighbor ;	
	
}

//==============================================================================================================

void CalculateBoundaryNum( int BCTypeNum , DSMC_PRE_BCTYPE *BCType , int *pBCIONum , int *pBCWallNum ){
	(*pBCIONum)	= 0 ;
	(*pBCWallNum)	= 0 ;

	for ( int i=0 ; i<BCTypeNum ; i++ ){
		if ( BCType[i].Type> 0) continue ;
		if ( BCType[i].Type == -3 || BCType[i].Type == -4 )
			(*pBCIONum)	+= BCType[i].FaceNum ;
		else if ( BCType[i].Type < -20 )
			(*pBCWallNum)	+= BCType[i].FaceNum ;
	}
}

//==============================================================================================================

void CreateNonUniformInletBC( int BCFaceNum , int *pInletSpecifiedNumDenNum , DSMC_PRE_NODE *Node , DSMC_PRE_BCFACE *BCFace , DSMC_PRE_BCTYPE *BCType , string FileName ){
	ifstream		Input ;
	string			GetLine ;
	int			InletNum , InletNo , NodeNo ;
	DSMC_PRE_INLET		*Inlet ;
	
	
	InletNum	= 0 ;
	InletNo		= 0 ;
	
	
	Input.open( FileName.c_str() , ios::in ) ;
	
	if ( Input.fail() ){
 		cout << "Fail to open file " << FileName << '\n';
		exit(1) ;
	}
	
	
	while ( getline( Input , GetLine ) ){
		InletNum++ ;
	}
	Input.clear() ;
	Input.close() ;	
	
	
	InletNum -= 3 ;
	cout << "Non-Uniform Inlet Boundary Number: " << InletNum << '\n' ;
	
	
	Inlet	= new DSMC_PRE_INLET[InletNum] ;
	
	
	Input.open( FileName.c_str() , ios::in ) ;
	
	if ( Input.fail() ){
 		cout << "Fail to open file " << FileName << '\n';
		exit(1) ;
	}
	
	
	// Read inlet information.
	getline( Input , GetLine ) ;
	Input >> (*pInletSpecifiedNumDenNum) ;
	getline( Input , GetLine ) ;
	cout << "The Number of Specified Number Density of Species: " << (*pInletSpecifiedNumDenNum) << '\n' ;
	
	
	for ( int i=0 ; i<InletNum ; i++ ){
		Inlet[i].AllocateMemory( (*pInletSpecifiedNumDenNum) ) ;	
	}
	
	
	getline( Input , GetLine ) ;
	for ( int i=0 ; i<InletNum ; i++ ){
		Input >> Inlet[i].XCoord >> Inlet[i].YCoord >> Inlet[i].ZCoord >> Inlet[i].XVel >> Inlet[i].YVel >> Inlet[i].ZVel ;
		
		for ( int j=0 ; j<(*pInletSpecifiedNumDenNum) ; j++ )
			Input >> Inlet[i].NumDen[j] ;
			
		Input >> Inlet[i].Temp ;
	}
	Input.clear() ; 
	Input.close() ;	


	for ( int i=0 ; i<BCFaceNum ; i++ ){
		if ( BCType[BCFace[i].Id].Type == -3 ){
			// Calculate the coordinate of the center of inlet face.
			for ( int j=0 ; j<BCFace[i].NodeNum ; j++ ){
				NodeNo	= BCFace[i].Node[j] ;
				
				BCFace[i].XCenter	+= Node[NodeNo].XCoord ;
				BCFace[i].YCenter	+= Node[NodeNo].YCoord ;
				BCFace[i].ZCenter	+= Node[NodeNo].ZCoord ;
			}
			BCFace[i].XCenter	/= BCFace[i].NodeNum ;
			BCFace[i].YCenter	/= BCFace[i].NodeNum ;
			BCFace[i].ZCenter	/= BCFace[i].NodeNum ;
		
		
		
			// Get a closest inlet no.
			InletNo	= GetInletNo( &BCFace[i] , Inlet , InletNum ) ;
		
		
			// Debug.
			//cout << "InletNo: " << InletNo << '\n' ;
		
		
			BCFace[i].XVel		= Inlet[InletNo].XVel ;
			BCFace[i].YVel		= Inlet[InletNo].YVel ;
			BCFace[i].ZVel		= Inlet[InletNo].ZVel ;
			BCFace[i].Temp		= Inlet[InletNo].Temp ;
			for ( int k=0 ; k<(*pInletSpecifiedNumDenNum) ; k++ )
				BCFace[i].NumDen.push_back( Inlet[InletNo].NumDen[k] ) ;
				
			// Debug.
			//cout << BCFace[i].Temp << '\n' ;
		}
	}
	
	
	delete [] Inlet ;
}

//==============================================================================================================

int GetInletNo( DSMC_PRE_BCFACE *_BCFace , DSMC_PRE_INLET *Inlet , int InletNum ){
	int		InletNo ;
	double		MinDistance , dx , dy , dz , Distance ;
	
	InletNo		= 0 ;
	MinDistance	= 1.E+10 ;
	
	for ( int i=0 ; i<InletNum ; i++ ){
		dx	= (_BCFace->XCenter - Inlet[i].XCoord)*(_BCFace->XCenter - Inlet[i].XCoord) ;
		dy	= (_BCFace->YCenter - Inlet[i].YCoord)*(_BCFace->YCenter - Inlet[i].YCoord) ;
		dz	= (_BCFace->ZCenter - Inlet[i].ZCoord)*(_BCFace->ZCenter - Inlet[i].ZCoord) ;
		
		Distance	= dx + dy + dz ;
		
		if ( Distance < MinDistance ){
			InletNo		= i ;
			MinDistance	= Distance ;
		}
	}
	
	return	InletNo ;
}

//==============================================================================================================

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
				DSMC_PRE_OPTION		*pOption ){

	ofstream	Output ;
	double		Scale ;
	
	Scale		= pOption->Scale ;


	// Output mesh (Mesh.inp).
	Output.open( "Mesh.inp" , ios::out | ios::trunc ) ;
	if ( !Output.fail() ){
		Output << "Mesh file for DSMC input\n" ;

		Output << "Nodes\n" ;
		Output << NodeNum << '\n' ;
		//Output << setw(12) << "ID" << setw(22) << "XCoord" << setw(22) << "YCoord" << setw(23) << "ZCoord" << '\n' ;
		for ( int i=0 ; i<NodeNum ; i++ ){
			Output	<< setw(12) << Node[i].Id
				<< setw(22) << Node[i].XCoord/Scale 
				<< setw(22) << Node[i].YCoord/Scale 
				<< setw(22) << Node[i].ZCoord/Scale << '\n' ;
		}
		Output << "\n\n" ;

		Output << "Cells\n" ;
		Output << CellNum << '\n' ;
		/*Output << setw(12) << "ID" << setw(8) << "Type" 
		       << setw(12) << "Node-0" << setw(12) << "Node-1" 
		       << setw(12) << "Node-2" << setw(12) << "Node-3" 
		       << setw(12) << "Node-4" << setw(12) << "Node-5"
		       << setw(12) << "Node-6" << setw(13) << "Node-7\n" ;*/
		for ( int i=0 ; i<CellNum ; i++ ){
			//Output << setw(12) << Cell[i].Id << setw(8) << Cell[i].Type ;
			Output << setw(12) << i << setw(8) << Cell[i].Type ;
			
			for ( int j=0 ; j<pMapping->NodeNum[Cell[i].Type] ; j++ )
				Output << setw(12) << Cell[i].pNode[j]->Id ;
			Output << '\n' ;
		}
		Output << "\n\n" ;

		Output << "Neighbors\n" ;
		/*Output << setw(12) << "ID" 
		       << setw(12) << "NBR-0" << setw(12) << "NBR-1" 
		       << setw(12) << "NBR-2" << setw(12) << "NBR-3" 
		       << setw(12) << "NBR-4" << setw(13) << "NBR-5\n" ;*/
		for ( int i=0 ; i<CellNum ; i++ ){
			//Output << setw(12) << Cell[i].Id ;
			Output << setw(12) << i ;
			
			for ( int j=0 ; j<pMapping->SurfaceNum[Cell[i].Type] ; j++ )
				Output << setw(12) << Cell[i].Neighbor[j] ;
			Output << '\n' ;
		}
	}
	Output.clear() ;
	Output.close() ;
	
	
	// Output I/O (Inlet.inp).
	Output.open( "Inlet.inp" , ios::out | ios::trunc ) ;
	if ( !Output.fail() ){
		Output << BCIONum << "          " << InletSpecifiedNumDenNum << '\n' ;


		for ( int i=0 ; i<BCFaceNum ; i++ ){
			int	j = BCFace[i].Id ;

			if ( BCType[j].Type == -3 || BCType[j].Type == -4 ){
				if ( pOption->OpenInputInlet == 0 ){
					Output << setw(11) << BCType[j].XVel 
						<< setw(11) << BCType[j].YVel 
						<< setw(11) << BCType[j].ZVel 
						<< setw(16) << BCType[j].NumDen
						<< setw(14) << BCType[j].Temp
						<< setw(8)  << BCType[j].CosineLawCoef ;
					
				}else if ( pOption->OpenInputInlet == 1 ){
					if ( BCType[j].Type == -3 ){
						Output << setw(11) << BCFace[i].XVel 
							<< setw(11) << BCFace[i].YVel 
							<< setw(11) << BCFace[i].ZVel ;
							
						for ( int k=0 ; k<InletSpecifiedNumDenNum ; k++ )	
							Output << setw(16) << BCFace[i].NumDen[k] ;
							
						Output << setw(14) << BCFace[i].Temp ;
					
					}else if ( BCType[j].Type == -4 ){
						Output << setw(11) << BCType[j].XVel 
							<< setw(11) << BCType[j].YVel 
							<< setw(11) << BCType[j].ZVel ;
							
						for ( int k=0 ; k<InletSpecifiedNumDenNum ; k++ )
							Output << setw(16) << "0." ;
							
						Output << setw(14) << BCType[j].Temp ;
					}
					
					Output << setw(8)  << BCType[j].CosineLawCoef ;
				}
				

				Output << setw(10) << BCFace[i].NodeNum ;
				for ( int k=0 ; k<BCFace[i].NodeNum ; k++ )
					Output << setw(10) << BCFace[i].Node[k] ;
				Output << '\n' ;
			}
		}
	}
	Output.clear() ;
	Output.close() ;
	
	
	
	// Output partition information (Partition.inp).
	Output.open( "Partition.inp" , ios::out | ios::trunc ) ;
	if ( !Output.fail() ){
		for ( int i=0 ; i<ProcessorNum ; i++ ){
			Output << setw(8) << i << setw(16) << Processor[i].Cell.size() << setw(6) << Processor[i].Neighbor.size() ;
			
			for ( int j=0 ; j<Processor[i].Neighbor.size() ; j++ )
				Output << setw(6) << Processor[i].Neighbor[j] ;
			Output << '\n' ;
		}
	}
	Output.clear() ;
	Output.close() ;
	
	
	// Output cell information (Cell-Information.inp).
	if ( pOption->OpenCell ){
		Output.open( "Cell-Information.inp" , ios::out | ios::trunc ) ;
		if ( !Output.fail() ){
			Output << setw(20) << "CellNo" << setw(20) << "MeanFreePath" << setw(20) << "Temp" << setw(20) << "Speed" << setw(20) << "Timestep" 
				<< setw(20) << "Weighting" << setw(20) << "T_Ratio" << setw(20) << "W_Ratio" << setw(20) << "SubcellNum" << '\n' ; 
			
			
			for ( int i=0 ; i<CellNum ; i++ )
				Output << setw(20) << i << setw(20) << Cell[i].MeanFreePath << setw(20) << Cell[i].Temp << setw(20) << Cell[i].Speed << setw(20) << Cell[i].Timestep 
					<< setw(20) << Cell[i].Weighting << setw(20) << Cell[i].TimestepRatio << setw(20) << Cell[i].WeightingRatio << setw(20) << Cell[i].SubcellNum << '\n' ; 
		}
		Output.clear() ;
		Output.close() ;
	}
}

//==============================================================================================================

void OutputMeshTec( int NodeNum , int CellNum , DSMC_PRE_NODE *Node , DSMC_PRE_CELL *Cell , CELLMAPPING *pMapping , int Dimension , double Scale ){
	int		OutNum = 0 ;
	ofstream	Output ;

	Output.open( "Mesh-tec.dat" , ios::out | ios::trunc ) ;

	if ( !Output.fail() ){
		Output << "TITLE = \"Finite volume dataset\"\n" ;
		if ( Dimension == 3 ){
			Output << "VARIABLES = \"X\", \"Y\", \"Z\", \"CPU_No\"\n" ;
			Output << "ZONE N=" << NodeNum << ", E=" << CellNum << ", DATAPACKING=BLOCK, ZONETYPE=FEBRICK\n" ;
			Output << "VARLOCATION=([1-3]=NODAL, [4]=CELLCENTERED)\n" ;
		}else if ( Dimension == 2 ){
			Output << "VARIABLES = \"X\", \"Y\", \"CPU_No\"\n" ;
			Output << "ZONE N=" << NodeNum << ", E=" << CellNum << ", DATAPACKING=BLOCK, ZONETYPE=FEQUADRILATERAL\n" ;
			Output << "VARLOCATION=([1-2]=NODAL, [3]=CELLCENTERED)\n" ;
		}

		for ( int i=0 ; i<NodeNum ; i++ ){
			OutNum++ ;
			Output << setw(15) << Node[i].XCoord/Scale ;
			if ( OutNum%6==0 ) Output << '\n' ;
		}

		for ( int i=0 ; i<NodeNum ; i++ ){
			OutNum++ ;
			Output << setw(15) << Node[i].YCoord/Scale ;
			if ( OutNum%6==0 ) Output << '\n' ;
		}


		// For 3D mesh
		if ( Dimension == 3 ){
			for ( int i=0 ; i<NodeNum ; i++ ){
				OutNum++ ;
				Output << setw(15) << Node[i].ZCoord/Scale ;
				if ( OutNum%6==0 ) Output << '\n' ;
			}
		}
		
		
		for ( int i=0 ; i<CellNum ; i++ ){
			OutNum++ ;
			Output << setw(15) << Cell[i].ProcessorNo ;
			if ( OutNum%6==0 ) Output << '\n' ;
		}

		Output << '\n' ;
		for ( int i=0 ; i<CellNum ; i++ ){
			// Tetra Cell.
			if ( Cell[i].Type == 0 ){
				Output << setw(12) << Cell[i].pNode[0]->Id+1 ;
				Output << setw(12) << Cell[i].pNode[1]->Id+1 ;
				Output << setw(12) << Cell[i].pNode[2]->Id+1 ;
				Output << setw(12) << Cell[i].pNode[2]->Id+1 ;
				Output << setw(12) << Cell[i].pNode[3]->Id+1 ;
				Output << setw(12) << Cell[i].pNode[3]->Id+1 ;
				Output << setw(12) << Cell[i].pNode[3]->Id+1 ;
				Output << setw(12) << Cell[i].pNode[3]->Id+1 ;

			// Hex. Cell
			}else if ( Cell[i].Type == 1 ){
				for ( int j=0 ; j<pMapping->NodeNum[1] ; j++ )
					Output << setw(12) << Cell[i].pNode[j]->Id+1 ;

			// Pyramid Cell
			}else if ( Cell[i].Type == 2 ){
				Output << setw(12) << Cell[i].pNode[0]->Id+1 ;
				Output << setw(12) << Cell[i].pNode[1]->Id+1 ;
				Output << setw(12) << Cell[i].pNode[2]->Id+1 ;
				Output << setw(12) << Cell[i].pNode[3]->Id+1 ;
				Output << setw(12) << Cell[i].pNode[4]->Id+1 ;
				Output << setw(12) << Cell[i].pNode[4]->Id+1 ;
				Output << setw(12) << Cell[i].pNode[4]->Id+1 ;
				Output << setw(12) << Cell[i].pNode[4]->Id+1 ;


			// Prism Cell
			}else if ( Cell[i].Type == 3 ){
				Output << setw(12) << Cell[i].pNode[0]->Id+1 ;
				Output << setw(12) << Cell[i].pNode[3]->Id+1 ;
				Output << setw(12) << Cell[i].pNode[4]->Id+1 ;
				Output << setw(12) << Cell[i].pNode[4]->Id+1 ;
				Output << setw(12) << Cell[i].pNode[1]->Id+1 ;
				Output << setw(12) << Cell[i].pNode[2]->Id+1 ;
				Output << setw(12) << Cell[i].pNode[5]->Id+1 ;
				Output << setw(12) << Cell[i].pNode[5]->Id+1 ;


			// 2D Triangular Cell
			}else if ( Cell[i].Type == 4 ){
				Output << setw(12) << Cell[i].pNode[0]->Id+1 ;
				Output << setw(12) << Cell[i].pNode[1]->Id+1 ;
				Output << setw(12) << Cell[i].pNode[2]->Id+1 ;
				Output << setw(12) << Cell[i].pNode[2]->Id+1 ;


			// 2D Quadrilateral Cell
			}else if ( Cell[i].Type == 5 ){
				Output << setw(12) << Cell[i].pNode[0]->Id+1 ;
				Output << setw(12) << Cell[i].pNode[1]->Id+1 ;
				Output << setw(12) << Cell[i].pNode[2]->Id+1 ;
				Output << setw(12) << Cell[i].pNode[3]->Id+1 ;
			}


			Output << '\n' ;
		}
	}

	
	Output.clear() ;
	Output.close() ;
}

//==============================================================================================================

void OutputInletTec(	int			NodeNum ,
			int			BCFaceNum , 
			int			InletSpecifiedNumDenNum ,
			DSMC_PRE_NODE		*Node , 
			DSMC_PRE_BCFACE		*BCFace , 
			DSMC_PRE_BCTYPE		*BCType ,
			DSMC_PRE_OPTION		*pOption , 
			CELLMAPPING 		*pMapping ){
	
	int				NodeNo , CellNo , OutNum ;
	double				Scale ;
	vector<DSMC_PRE_BCFACE>		iBCFace ;
	vector<DSMC_PRE_NODE>		iNode ;
	DSMC_PRE_NODECELL		*NodeCell ;
	ofstream			Output ;
	string				*WordNumDen , *WordCountNum ;
	char				buffer[8] ;
	
	
	NodeCell	= new DSMC_PRE_NODECELL[NodeNum] ;
	WordNumDen	= new string[InletSpecifiedNumDenNum] ;
	WordCountNum	= new string[InletSpecifiedNumDenNum] ;
	
	OutNum		= 0 ;
	Scale		= pOption->Scale ;
	
	
	for ( int i=0 ; i<BCFaceNum ; i++ ){
		int	j = BCFace[i].Id ;

		if ( BCType[j].Type == -3 ){
			iBCFace.push_back( BCFace[i] ) ;
			
			if ( pOption->OpenInputInlet == 0 ){
				iBCFace[iBCFace.size()-1].XVel		= BCType[j].XVel ;
				iBCFace[iBCFace.size()-1].YVel		= BCType[j].YVel ;
				iBCFace[iBCFace.size()-1].ZVel		= BCType[j].ZVel ;
				//iBCFace[iBCFace.size()-1].NumDen	= BCType[j].NumDen ;
				iBCFace[iBCFace.size()-1].Temp		= BCType[j].Temp ;
				
				iBCFace[iBCFace.size()-1].NumDen.push_back( BCType[j].NumDen ) ;
			}
		
		}
	}
	
	
	// Debug.
	//cout << "inlet cell: " << iBCFace.size() << '\n' ;
	
	
	
	for ( int i=0 ; i<iBCFace.size() ; i++ ){
		for ( int j=0 ; j<iBCFace[i].NodeNum ; j++ ){
			NodeNo	= iBCFace[i].Node[j] ;
			
			NodeCell[NodeNo].CellRelation.push_back(i) ;
		}	
	}
	
	
	NodeNo	= 0 ;
	for ( int i=0 ; i<NodeNum ; i++ ){
		if ( NodeCell[i].CellRelation.size() > 0 ){
			iNode.push_back( Node[i] ) ;
			
			for ( int j=0 ; j<NodeCell[i].CellRelation.size() ; j++ ){
				CellNo	= NodeCell[i].CellRelation[j] ;
				
				for ( int z=0 ; z<iBCFace[CellNo].NodeNum ; z++ ){
					if ( iBCFace[CellNo].Node[z] == i )
						iBCFace[CellNo].Node[z]	= NodeNo ;
				}
			}
			
			NodeNo++ ;
		}
	}
		
	//cout << "iNode Num: " << iNode.size() << '\n' ;
	
	
	for ( int i=0 ; i<InletSpecifiedNumDenNum ; i++ ){
		sprintf( buffer , "%d" , i+1 ) ;
		
		WordCountNum[i]	= buffer ;
		WordNumDen[i]	= "NumDen-" + WordCountNum[i] ;
	}
	
	
	Output.open( "Inlet-tec.dat" , ios::out | ios::trunc ) ;
		

	if ( !Output.fail() ){
		Output << "TITLE = \"Finite volume dataset\"\n" ;
		Output << "VARIABLES = \"X\", \"Y\", \"Z\", \"U-Vel\", \"V-Vel\", \"W-Vel\", \"Temp\"" ;
		for ( int i=0 ; i<InletSpecifiedNumDenNum ; i++ )
			Output << ", \"" << WordNumDen[i] << "\"" ;
		Output << '\n' ;
		Output << "ZONE N=" << iNode.size() << ", E=" << iBCFace.size() << ", DATAPACKING=BLOCK, ZONETYPE=FEQUADRILATERAL\n" ;
		Output << "VARLOCATION=([1-3]=NODAL, [4-" << (7+InletSpecifiedNumDenNum) << "]=CELLCENTERED)\n" ;
		
		
		for ( int i=0 ; i<iNode.size() ; i++ ){
			OutNum++ ;
			Output << setw(18) << iNode[i].XCoord/Scale ;
			if ( OutNum%6 == 0 ) Output << '\n' ;
		}
		
		for ( int i=0 ; i<iNode.size() ; i++ ){
			OutNum++ ;
			Output << setw(18) << iNode[i].YCoord/Scale ;
			if ( OutNum%6 == 0 ) Output << '\n' ;
		}
		
		for ( int i=0 ; i<iNode.size() ; i++ ){
			OutNum++ ;
			Output << setw(18) << iNode[i].ZCoord/Scale ;
			if ( OutNum%6 == 0 ) Output << '\n' ;
		}
		
		for ( int i=0 ; i<iBCFace.size() ; i++ ){
			OutNum++ ;
			Output << setw(18) << iBCFace[i].XVel ;
			if ( OutNum%6 == 0 ) Output << '\n' ;	
		}
		
		for ( int i=0 ; i<iBCFace.size() ; i++ ){
			OutNum++ ;
			Output << setw(18) << iBCFace[i].YVel ;
			if ( OutNum%6 == 0 ) Output << '\n' ;	
		}
		
		for ( int i=0 ; i<iBCFace.size() ; i++ ){
			OutNum++ ;
			Output << setw(18) << iBCFace[i].ZVel ;
			if ( OutNum%6 == 0 ) Output << '\n' ;	
		}
		
		for ( int i=0 ; i<iBCFace.size() ; i++ ){
			OutNum++ ;
			Output << setw(18) << iBCFace[i].Temp ;
			if ( OutNum%6 == 0 ) Output << '\n' ;	
		}
		
		for ( int j=0 ; j<InletSpecifiedNumDenNum ; j++ ){
			for ( int i=0 ; i<iBCFace.size() ; i++ ){
				OutNum++ ;
				Output << setw(18) << iBCFace[i].NumDen[j] ;
				if ( OutNum%6 == 0 ) Output << '\n' ;	
			}
		}
		
		
		Output << "\n\n" ;
		for ( int i=0 ; i<iBCFace.size() ; i++ ){
			for ( int j=0 ; j<iBCFace[i].NodeNum ; j++ )
				Output << setw(18) << ( iBCFace[i].Node[j]+1 ) ;
				
			if ( iBCFace[i].NodeNum == 3 )
				Output << setw(18) << ( iBCFace[i].Node[2]+1 ) ;
				
			Output << '\n' ;
		}
	}
	
	Output.clear() ;
	Output.close() ;
	
	
	iBCFace.clear() ;
	iNode.clear() ;
	
	delete [] NodeCell ;
	delete [] WordNumDen ;
	delete [] WordCountNum ;
}

//==============================================================================================================

void ReadResult( int ResultNum , DSMC_PRE_CELL *Cell , DSMC_PRE_RESULT *Result , string Filename ){
	ifstream	Input ;
	string		buffer ;
	int		CellNo ;
	
	
	Input.open( Filename.c_str() , ios::in ) ;
		
	if ( !Input ){
		cout << "Fail to open " << Filename << '\n' ;	
		exit(1) ;
	}
	
	getline( Input , buffer ) ;
	for ( int i=0 ; i<ResultNum ; i++ ){
		Input >> CellNo ;
		Input >> Result[CellNo].XCenter >> Result[CellNo].YCenter >> Result[CellNo].ZCenter >> Result[CellNo].Density
			>> Result[CellNo].NumDensity >> Result[CellNo].XVel >> Result[CellNo].YVel >> Result[CellNo].ZVel 
			>> Result[CellNo].Temp >> Result[CellNo].TransTemp >> Result[CellNo].RotTemp >> Result[CellNo].VibTemp 
			>> Result[CellNo].AveParticleNum >> Result[CellNo].MeanCollSpacingMeanFreePath >> Result[CellNo].ProcessorNo ;
			
		Cell[CellNo].Weight	= (int)(Result[CellNo].AveParticleNum*100.) ;
	}
	
	Input.clear() ;
	Input.close() ;
}

//==============================================================================================================

void ReadMesh( int CellNum , DSMC_PRE_NODE *Node , DSMC_PRE_CELL *Cell , CELLMAPPING *pMapping , string Filename ){
	int			cellnum , NodeNo , Type ;
	ifstream		Input ;
	string			buffer ;


	Input.open( Filename.c_str() , ios::in ) ;

	if ( !Input ){
		cout << "Failed to open file " << Filename << endl ;
		exit(1) ;
	}


	while ( getline( Input , buffer ) ){
		if ( FindString(buffer , "Cells") ){
			Input >> cellnum ;
		
			// Check cell number is right.
			if ( CellNum != cellnum ){
				cout << "Cell number is error. CellNum: " << CellNum << ", cellnum: " << cellnum << endl ;
				exit(1) ;
			}
		
			for ( int i=0 ; i<CellNum ; i++ ){
				Input >> Cell[i].Id >> Cell[i].Type ;
				
				for ( int j=0 ; j<pMapping->NodeNum[Cell[i].Type] ; j++){
					Input >> NodeNo ;
					Cell[i].pNode[j]	= &Node[NodeNo] ;
				}
			}
		}
	}
	// To close the file.
	Input.clear() ;
	Input.close() ;
}

//==============================================================================================================

void ReadCellInformation( int CellNum , DSMC_PRE_CELL *Cell , string Filename ){
	int			CellNo ;
	ifstream		Input ;
	string			buffer ;


	Input.open( Filename.c_str() , ios::in ) ;

	if ( !Input ){
		cout << "Failed to open file " << Filename << endl ;
		exit(1) ;
	}


	getline( Input , buffer ) ;
	for ( int i=0 ; i<CellNum ; i++ ){
		Input >> CellNo ;
		Input >> Cell[CellNo].MeanFreePath >> Cell[CellNo].Temp >> Cell[CellNo].Speed >> Cell[CellNo].Timestep >> Cell[CellNo].Weighting >> Cell[CellNo].TimestepRatio 
			>> Cell[CellNo].WeightingRatio >> Cell[CellNo].SubcellNum ;
	}
	
	
	// To close the file.
	Input.clear() ;
	Input.close() ;
}

//==============================================================================================================

int GetFileNum( string Filename ){
	ifstream	Input ;
	int		FileNum ;
	string		GetLine ;
	
	
	FileNum	= 0 ;

	
	Input.open( Filename.c_str() , ios::in ) ;

		
	if ( !Input ){
		cout << "Fail to open " << Filename << '\n' ;	
		exit(1) ;
	}


	while ( getline( Input , GetLine ) ){
		FileNum++ ;
	}
	
	Input.clear() ;
	Input.close() ;

	return	FileNum ;
}

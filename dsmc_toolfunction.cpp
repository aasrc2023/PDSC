#include <cmath>
#include <cstdlib>
#include <iostream>
#include "dsmc_toolfunction.h"
#include "dsmc_parameter.h"


using namespace std ;

//==============================================================================================================

double CalculateNormalVelocity2D( double XVel , double YVel , DSMC_NODE *_Node0 , DSMC_NODE *_Node1 ){
	double		NormalVel , Vector_x , Vector_y , NormalVector_x , NormalVector_y ;
	
	Vector_x	= _Node1->XCoord - _Node0->XCoord ;
	Vector_y	= _Node1->YCoord - _Node0->YCoord ;
	
	NormalVector_x	= Vector_y ;
	NormalVector_y	= -Vector_x ;
	
	NormalVel	= (XVel*NormalVector_x+YVel*NormalVector_y)/sqrt((NormalVector_x*NormalVector_x)+(NormalVector_y*NormalVector_y)) ;
	
	return	NormalVel ;
}


double CalculateNormalVelocity3D( double XVel , double YVel , double ZVel , DSMC_INLET *_Inlet , DSMC_CELL *_Cell ){
	int		FaceNo ;
	double		FaceFA , FaceFB , FaceFC , NormalVel ;

	FaceNo		= _Inlet->FaceNo ;
	
	FaceFA		= _Cell->FaceFA[FaceNo] ;
	FaceFB		= _Cell->FaceFB[FaceNo] ;
	FaceFC		= _Cell->FaceFC[FaceNo] ;
	
	NormalVel	= (XVel*FaceFA + YVel*FaceFB + ZVel*FaceFC) ;
	
	// Debug.
	//cout << "FA: " << FaceFA << ", FB: " << FaceFB << ", FC: " << FaceFC << '\n' ;
	//cout << "NormalVel: " << NormalVel << '\n' ;
	
	return	NormalVel ;
}

//==============================================================================================================

void InverseVelocityToNormPara2D( double *pNormVel , double *pParaVel , double XVel , double YVel , DSMC_NODE Node0 , DSMC_NODE Node1 ){
	double		VecPara_x ,VecPara_y , VecNorm_x , VecNorm_y , Norm ;
	
	VecPara_x	= Node1.XCoord - Node0.XCoord ;
	VecPara_y	= Node1.YCoord - Node0.YCoord ;
	VecNorm_x	= VecPara_y ;
	VecNorm_y	= -VecPara_x ;
	Norm		= sqrt( VecPara_x*VecPara_x + VecPara_y*VecPara_y ) ;
	
	(*pNormVel)	= (XVel*VecNorm_x + YVel*VecNorm_y)/Norm ;
	(*pParaVel)	= (XVel*VecPara_x + YVel*VecPara_y)/Norm ;
}


void InverseVelocityToNormPara3D(	double		*pNormVel , 
					double		*pParaVelX , 
					double		*pParaVelY , 
					double		XVel , 
					double		YVel , 
					double		ZVel , 
					DSMC_NODE	*h_Node , 
					DSMC_CELL	*_Cell , 
					int		FaceNo , 
					CELLMAPPING	*h_pMapping ){
	
	int		Type , Node[2] ;			
	double		VecParaX[2] ,VecParaY[2] , VecParaZ[2] , VecNormX , VecNormY , VecNormZ , Norm ;
	
	Type		= _Cell->Type ;
	
	Node[0]		= h_pMapping->Node[Type][FaceNo][0] ;
	Node[1]		= h_pMapping->Node[Type][FaceNo][1] ;
	Node[0]		= _Cell->Node[Node[0]] ;
	Node[1]		= _Cell->Node[Node[1]] ;
	
	VecParaX[0]	= h_Node[Node[1]].XCoord - h_Node[Node[0]].XCoord ;
	VecParaY[0]	= h_Node[Node[1]].YCoord - h_Node[Node[0]].YCoord ;
	VecParaZ[0]	= h_Node[Node[1]].ZCoord - h_Node[Node[0]].ZCoord ;
	Norm		= sqrt( VecParaX[0]*VecParaX[0] + VecParaY[0]*VecParaY[0] + VecParaZ[0]*VecParaZ[0]) ;
	
	VecParaX[0]	/= Norm ;
	VecParaY[0]	/= Norm ;
	VecParaZ[0]	/= Norm ;
	
	VecNormX	= _Cell->FaceFA[FaceNo] ;
	VecNormY	= _Cell->FaceFB[FaceNo] ;
	VecNormZ	= _Cell->FaceFC[FaceNo] ;
	
	VecParaX[1]	= VecNormY*VecParaZ[0] - VecNormZ*VecParaY[0] ;
	VecParaY[1]	= VecNormZ*VecParaX[0] - VecNormX*VecParaZ[0] ;
	VecParaZ[1]	= VecNormX*VecParaY[0] - VecNormY*VecParaX[0] ;
	
	(*pNormVel)	= XVel*VecNormX + YVel*VecNormY + ZVel*VecNormZ ;
	(*pParaVelX)	= XVel*VecParaX[0] + YVel*VecParaY[0] + ZVel*VecParaZ[0] ;
	(*pParaVelY)	= XVel*VecParaX[1] + YVel*VecParaY[1] + ZVel*VecParaZ[1] ;
	
	// Debug.
	//cout << "NormVel: " << (*pNormVel) << ", ParaVelX: " << (*pParaVelX) << ", ParaVelY: " << (*pParaVelY) << '\n' ;
	//cout << "XVel: " << XVel << ", YVel: " << YVel << ", ZVel: " << ZVel << '\n' ;
	//getchar() ;
}

//==============================================================================================================

void InverseVelocityToCartesian2D( double *pXVel , double *pYVel , double NormVel , double ParaVel , DSMC_NODE Node0 , DSMC_NODE Node1 ){
	double		VecPara_x ,VecPara_y , VecNorm_x , VecNorm_y , Norm ;
	
	VecPara_x	= Node1.XCoord - Node0.XCoord ;
	VecPara_y	= Node1.YCoord - Node0.YCoord ;
	VecNorm_x	= VecPara_y ;
	VecNorm_y	= -VecPara_x ;
	Norm		= sqrt( VecPara_x*VecPara_x + VecPara_y*VecPara_y ) ;
	
	(*pXVel)	= (NormVel*VecNorm_x + ParaVel*VecPara_x)/Norm ;
	(*pYVel)	= (NormVel*VecNorm_y + ParaVel*VecPara_y)/Norm ;
}


void InverseVelocityToCartesian3D(	double		*pXVel , 
					double		*pYVel , 
					double		*pZVel , 
					double		NormVel , 
					double		ParaVelX , 
					double		ParaVelY , 
					DSMC_NODE	*h_Node , 
					DSMC_CELL	*_Cell , 
					int		FaceNo , 
					CELLMAPPING	*h_pMapping ){
	
	int		Type , Node[2] ;			
	double		VecParaX[2] ,VecParaY[2] , VecParaZ[2] , VecNormX , VecNormY , VecNormZ , Norm ;
	
	Type		= _Cell->Type ;
	
	Node[0]		= h_pMapping->Node[Type][FaceNo][0] ;
	Node[1]		= h_pMapping->Node[Type][FaceNo][1] ;
	Node[0]		= _Cell->Node[Node[0]] ;
	Node[1]		= _Cell->Node[Node[1]] ;
	
	VecParaX[0]	= h_Node[Node[1]].XCoord - h_Node[Node[0]].XCoord ;
	VecParaY[0]	= h_Node[Node[1]].YCoord - h_Node[Node[0]].YCoord ;
	VecParaZ[0]	= h_Node[Node[1]].ZCoord - h_Node[Node[0]].ZCoord ;
	Norm		= sqrt( VecParaX[0]*VecParaX[0] + VecParaY[0]*VecParaY[0] + VecParaZ[0]*VecParaZ[0]) ;
	
	VecParaX[0]	/= Norm ;
	VecParaY[0]	/= Norm ;
	VecParaZ[0]	/= Norm ;
	
	VecNormX	= _Cell->FaceFA[FaceNo] ;
	VecNormY	= _Cell->FaceFB[FaceNo] ;
	VecNormZ	= _Cell->FaceFC[FaceNo] ;
	
	VecParaX[1]	= VecNormY*VecParaZ[0] - VecNormZ*VecParaY[0] ;
	VecParaY[1]	= VecNormZ*VecParaX[0] - VecNormX*VecParaZ[0] ;
	VecParaZ[1]	= VecNormX*VecParaY[0] - VecNormY*VecParaX[0] ;
	
	
	(*pXVel)	= NormVel*VecNormX + ParaVelX*VecParaX[0] + ParaVelY*VecParaX[1] ;
	(*pYVel)	= NormVel*VecNormY + ParaVelX*VecParaY[0] + ParaVelY*VecParaY[1] ;
	(*pZVel)	= NormVel*VecNormZ + ParaVelX*VecParaZ[0] + ParaVelY*VecParaZ[1] ;
	
	// Debug.
	//cout << "NormVel: " << NormVel << ", ParaVelX: " << ParaVelX << ", ParaVelY: " << ParaVelY << '\n' ;
	//cout << "XVel: " << (*pXVel) << ", YVel: " << (*pYVel) << ", ZVel: " << (*pZVel) << '\n' ;
	//getchar() ;
}

//==============================================================================================================

/*bool InCell2D( DSMC_NODE *h_Node , DSMC_CELL *_Cell , double XCoord , double YCoord , CELLMAPPING *pMapping ){
	bool	InCell ;
	int	Type , Num , OnPoint ;
	double	PX ; 
	
	InCell	= false ;
	Type	= _Cell->Type ;
	Num	= 0 ;
	OnPoint	= 0 ;

	if ( XCoord < _Cell->MinXCoord || XCoord > _Cell->MaxXCoord || YCoord < _Cell->MinYCoord || YCoord > _Cell->MaxYCoord ){
		return	InCell ;	
	}
	
	
	for ( int i=0 ; i<pMapping->SurfaceNum[Type] ; i++ ){
		if ( _Cell->FaceFA[i] ) <= 1.E-9 ) continue ;
		
		
		PX	= (_Cell->FaceFB[i]*YCoord + _Cell->FaceFC[i])/(-1.*_Cell->FaceFA[i]) ;

		if ( PX >= XCoord && PX <= _Cell->MaxXCoord ){
			Num++ ;
		}
	}
	if ( Num >= 3 ) return	InCell ;
	if ( Num%2 == 1 ) InCell = true ;
	
	return	InCell ;
}*/


bool InCell2D( DSMC_NODE *h_Node , DSMC_CELL *_Cell , double XCoord , double YCoord , CELLMAPPING *h_pMapping , int Debug ){
	bool	InCell ;
	int	Type , NodeA , NodeB ;
	double	VectorA_X , VectorA_Y , VectorB_X , VectorB_Y , Cross ; 
	
	InCell	= true ;
	Type	= _Cell->Type ;

	
	for ( int i=0 ; i<h_pMapping->SurfaceNum[Type] ; i++ ){
		NodeA		= _Cell->Node[h_pMapping->Node[Type][i][0]] ;
		NodeB		= _Cell->Node[h_pMapping->Node[Type][i][1]] ;
		
		VectorA_X	= h_Node[NodeA].XCoord - XCoord ;
		VectorA_Y	= h_Node[NodeA].YCoord - YCoord ;
		VectorB_X	= h_Node[NodeB].XCoord - XCoord ;
		VectorB_Y	= h_Node[NodeB].YCoord - YCoord ;

		Cross	= (VectorA_X*VectorB_Y) - (VectorA_Y*VectorB_X) ;
		
		
		//if ( Debug == 1 ) cout << "Cross: " << Cross << '\n' ;
		
		if ( Cross < 0. ){
			InCell	= false ;
			break ;	
		}
	}
	
	//if ( Debug == 1 ){
	//	cout << "Cell No: " << _Cell->Id << ", XCoord: " << XCoord << ", YCoord: " << YCoord << '\n' ;
	//}

	return	InCell ;
}


bool InCell3D( DSMC_NODE *h_Node , DSMC_CELL *_Cell , double XCoord , double YCoord , double ZCoord , CELLMAPPING *h_pMapping , int Debug ){
	bool	InCell ;
	int	Type ;
	double	Cross ; 
	
	InCell	= true ;
	Type	= _Cell->Type ;

	
	for ( int i=0 ; i<h_pMapping->SurfaceNum[Type] ; i++ ){
		Cross	= -1.*( _Cell->FaceFA[i]*XCoord + _Cell->FaceFB[i]*YCoord + _Cell->FaceFC[i]*ZCoord + _Cell->FaceFD[i] ) ;
		
		if ( Cross < 0. ){
			InCell	= false ;
			break ;	
		}
	}
	
	return	InCell ;
}

//==============================================================================================================

void CreateParticlePosition2D(	double		*XCoord  , 
				double		*YCoord , 
				double		*ZCoord , 
				DSMC_NODE	*h_Node , 
				DSMC_CELL	*_Cell , 
				CELLMAPPING	*h_pMapping ){
	(*ZCoord)	= 0. ;
				
	do{
		(*XCoord)	= _Cell->MinXCoord + (_Cell->MaxXCoord - _Cell->MinXCoord)*Randn() ;
		(*YCoord)	= _Cell->MinYCoord + (_Cell->MaxYCoord - _Cell->MinYCoord)*Randn() ;
		
	}while( !InCell2D( h_Node , _Cell , (*XCoord) , (*YCoord) , h_pMapping , 0 ) ) ;
}


void CreateParticlePosition3D(	double		*XCoord  , 
				double		*YCoord , 
				double		*ZCoord , 
				DSMC_NODE	*h_Node , 
				DSMC_CELL	*_Cell , 
				CELLMAPPING	*h_pMapping ){
				
	do{
		(*XCoord)	= _Cell->MinXCoord + (_Cell->MaxXCoord - _Cell->MinXCoord)*Randn() ;
		(*YCoord)	= _Cell->MinYCoord + (_Cell->MaxYCoord - _Cell->MinYCoord)*Randn() ;
		(*ZCoord)	= _Cell->MinZCoord + (_Cell->MaxZCoord - _Cell->MinZCoord)*Randn() ;
		
	}while( !InCell3D( h_Node , _Cell , (*XCoord) , (*YCoord) , (*ZCoord) , h_pMapping , 0 ) ) ;
}

//==============================================================================================================

void CreateParticlePositionEnter3D(	double		*XCoord , 
					double		*YCoord ,
					double		*ZCoord ,
					DSMC_INLET	*_Inlet ,
					DSMC_NODE	*h_Node ,  
					DSMC_CELL	*h_Cell ,
					CELLMAPPING	*h_pMapping ){
						
	int		CellNo , FaceNo , Type , Node[3] ;
	double		Rand1 , Rand2 , Area , Area1 , Area2 ;
	
	Rand1	= sqrt( Randn() ) ;
	Rand2	= Randn() ;
	
	CellNo	= _Inlet->CellNo ;
	FaceNo	= _Inlet->FaceNo ;
	Area	= _Inlet->Area ;
	Area1	= 0. ;
	
	Type	= h_Cell[CellNo].Type ;
	
	if ( _Inlet->NodeNum == 3 ){
		Node[0]		= h_pMapping->Node[Type][FaceNo][0] ;
		Node[1]		= h_pMapping->Node[Type][FaceNo][1] ;
		Node[2]		= h_pMapping->Node[Type][FaceNo][2] ;
		
		Node[0]		= h_Cell[CellNo].Node[Node[0]] ;
		Node[1]		= h_Cell[CellNo].Node[Node[1]] ;
		Node[2]		= h_Cell[CellNo].Node[Node[2]] ;
		
		
		(*XCoord)	= h_Node[Node[0]].XCoord + 
				 (h_Node[Node[1]].XCoord - h_Node[Node[0]].XCoord)*Rand1 + 
				 (h_Node[Node[2]].XCoord - h_Node[Node[1]].XCoord)*Rand1*Rand2 ;
		(*YCoord)	= h_Node[Node[0]].YCoord + 
				 (h_Node[Node[1]].YCoord - h_Node[Node[0]].YCoord)*Rand1 + 
				 (h_Node[Node[2]].YCoord - h_Node[Node[1]].YCoord)*Rand1*Rand2 ;
		(*ZCoord)	= h_Node[Node[0]].ZCoord + 
				 (h_Node[Node[1]].ZCoord - h_Node[Node[0]].ZCoord)*Rand1 + 
				 (h_Node[Node[2]].ZCoord - h_Node[Node[1]].ZCoord)*Rand1*Rand2 ;

	}else if ( _Inlet->NodeNum == 4 ){
		
		Node[0]		= h_pMapping->Node[Type][FaceNo][0] ;
		Node[1]		= h_pMapping->Node[Type][FaceNo][1] ;
		Node[2]		= h_pMapping->Node[Type][FaceNo][2] ;

		Node[0]		= h_Cell[CellNo].Node[Node[0]] ;
		Node[1]		= h_Cell[CellNo].Node[Node[1]] ;
		Node[2]		= h_Cell[CellNo].Node[Node[2]] ;
		
		Area1	= CalculateArea3D( &h_Node[Node[0]] , &h_Node[Node[1]] , &h_Node[Node[2]] ) ;
	
		if ( Randn() <= (Area1/Area) ){
			(*XCoord)	= h_Node[Node[0]].XCoord + 
					 (h_Node[Node[1]].XCoord - h_Node[Node[0]].XCoord)*Rand1 + 
					 (h_Node[Node[2]].XCoord - h_Node[Node[1]].XCoord)*Rand1*Rand2 ;
			(*YCoord)	= h_Node[Node[0]].YCoord + 
					 (h_Node[Node[1]].YCoord - h_Node[Node[0]].YCoord)*Rand1 + 
					 (h_Node[Node[2]].YCoord - h_Node[Node[1]].YCoord)*Rand1*Rand2 ;
			(*ZCoord)	= h_Node[Node[0]].ZCoord + 
					 (h_Node[Node[1]].ZCoord - h_Node[Node[0]].ZCoord)*Rand1 + 
					 (h_Node[Node[2]].ZCoord - h_Node[Node[1]].ZCoord)*Rand1*Rand2 ;
		}else{
			Node[0]		= h_pMapping->Node[Type][FaceNo][2] ;
			Node[1]		= h_pMapping->Node[Type][FaceNo][3] ;
			Node[2]		= h_pMapping->Node[Type][FaceNo][0] ;

			Node[0]		= h_Cell[CellNo].Node[Node[0]] ;
			Node[1]		= h_Cell[CellNo].Node[Node[1]] ;
			Node[2]		= h_Cell[CellNo].Node[Node[2]] ;
			
			
			(*XCoord)	= h_Node[Node[0]].XCoord + 
					 (h_Node[Node[1]].XCoord - h_Node[Node[0]].XCoord)*Rand1 + 
					 (h_Node[Node[2]].XCoord - h_Node[Node[1]].XCoord)*Rand1*Rand2 ;
			(*YCoord)	= h_Node[Node[0]].YCoord + 
					 (h_Node[Node[1]].YCoord - h_Node[Node[0]].YCoord)*Rand1 + 
					 (h_Node[Node[2]].YCoord - h_Node[Node[1]].YCoord)*Rand1*Rand2 ;
			(*ZCoord)	= h_Node[Node[0]].ZCoord + 
					 (h_Node[Node[1]].ZCoord - h_Node[Node[0]].ZCoord)*Rand1 + 
					 (h_Node[Node[2]].ZCoord - h_Node[Node[1]].ZCoord)*Rand1*Rand2 ;
		}
	}
}

//==============================================================================================================

double CalculateArea3D( DSMC_NODE *_Node1 , DSMC_NODE *_Node2 , DSMC_NODE *_Node3 ){
	double 	Area , VectorX[2] , VectorY[2] , VectorZ[2] , Buffer[3] ;;
	
	Area	= 0. ;
	
	VectorX[0]	= _Node1->XCoord - _Node2->XCoord ;
	VectorY[0]	= _Node1->YCoord - _Node2->YCoord ;
	VectorZ[0]	= _Node1->ZCoord - _Node2->ZCoord ;
	
	VectorX[1]	= _Node3->XCoord - _Node2->XCoord ;
	VectorY[1]	= _Node3->YCoord - _Node2->YCoord ;
	VectorZ[1]	= _Node3->ZCoord - _Node2->ZCoord ;
	
	Buffer[0]	= VectorY[0]*VectorZ[1] - VectorZ[0]*VectorY[1] ;
	Buffer[1]	= VectorZ[0]*VectorX[1] - VectorX[0]*VectorZ[1] ;
	Buffer[2]	= VectorX[0]*VectorY[1] - VectorY[0]*VectorX[1] ;
	
	Area	= (sqrt( Buffer[0]*Buffer[0] + Buffer[1]*Buffer[1] + Buffer[2]*Buffer[2] )/2.) ;
	
	return	Area ;
}

//==============================================================================================================

void CalculateCellFaceArea3D( double *Area , double *Area1 , double *Area2 , DSMC_NODE *h_Node , DSMC_CELL *_Cell , int FaceNo , CELLMAPPING *h_pMapping ){
	

	int		Type , Node[3] ;
	double		VectorX[2] , VectorY[2] , VectorZ[2] , Buffer[3] ;
	
	Type		= _Cell->Type ;
	(*Area)		= 0. ;
	(*Area1)	= 0. ;
	(*Area2)	= 0. ;
	
	
	// Calculate area of first triangular cell.
	Node[0]		= h_pMapping->Node[Type][FaceNo][0] ;
	Node[1]		= h_pMapping->Node[Type][FaceNo][1] ;
	Node[2]		= h_pMapping->Node[Type][FaceNo][3] ;

	Node[0]		= _Cell->Node[Node[0]] ;
	Node[1]		= _Cell->Node[Node[1]] ;
	Node[2]		= _Cell->Node[Node[2]] ;
	
	
	VectorX[0]	= h_Node[Node[1]].XCoord - h_Node[Node[0]].XCoord ;
	VectorY[0]	= h_Node[Node[1]].YCoord - h_Node[Node[0]].YCoord ;
	VectorZ[0]	= h_Node[Node[1]].ZCoord - h_Node[Node[0]].ZCoord ;
	
	VectorX[1]	= h_Node[Node[2]].XCoord - h_Node[Node[0]].XCoord ;
	VectorY[1]	= h_Node[Node[2]].YCoord - h_Node[Node[0]].YCoord ;
	VectorZ[1]	= h_Node[Node[2]].ZCoord - h_Node[Node[0]].ZCoord ;
	
	Buffer[0]	= VectorY[0]*VectorZ[1] - VectorZ[0]*VectorY[1] ;
	Buffer[1]	= VectorZ[0]*VectorX[1] - VectorX[0]*VectorZ[1] ;
	Buffer[2]	= VectorX[0]*VectorY[1] - VectorY[0]*VectorX[1] ;
	
	(*Area1)	= (sqrt( Buffer[0]*Buffer[0] + Buffer[1]*Buffer[1] + Buffer[2]*Buffer[2] )/2.) ;
	
	
	// Calculate area of second triangular cell.
	Node[0]		= h_pMapping->Node[Type][FaceNo][2] ;
	Node[1]		= h_pMapping->Node[Type][FaceNo][1] ;
	Node[2]		= h_pMapping->Node[Type][FaceNo][3] ;

	Node[0]		= _Cell->Node[Node[0]] ;
	Node[1]		= _Cell->Node[Node[1]] ;
	Node[2]		= _Cell->Node[Node[2]] ;
	
	
	VectorX[0]	= h_Node[Node[1]].XCoord - h_Node[Node[0]].XCoord ;
	VectorY[0]	= h_Node[Node[1]].YCoord - h_Node[Node[0]].YCoord ;
	VectorZ[0]	= h_Node[Node[1]].ZCoord - h_Node[Node[0]].ZCoord ;
	
	VectorX[1]	= h_Node[Node[2]].XCoord - h_Node[Node[0]].XCoord ;
	VectorY[1]	= h_Node[Node[2]].YCoord - h_Node[Node[0]].YCoord ;
	VectorZ[1]	= h_Node[Node[2]].ZCoord - h_Node[Node[0]].ZCoord ;
	
	Buffer[0]	= VectorY[0]*VectorZ[1] - VectorZ[0]*VectorY[1] ;
	Buffer[1]	= VectorZ[0]*VectorX[1] - VectorX[0]*VectorZ[1] ;
	Buffer[2]	= VectorX[0]*VectorY[1] - VectorY[0]*VectorX[1] ;
	
	(*Area2)	= (sqrt( Buffer[0]*Buffer[0] + Buffer[1]*Buffer[1] + Buffer[2]*Buffer[2] )/2.) ;

	(*Area)		= (*Area1) + (*Area2) ;
}

//==============================================================================================================

void RandVelocity( double *u , double *v , double MostProbableSpeed ){
	double		a , b ;

	a	= sqrt(-log(Randn())) ;
	b	= 6.283185308*Randn() ;
	
	(*u)	= a*sin(b)*MostProbableSpeed ;
	(*v)	= a*cos(b)*MostProbableSpeed ;
}

//==============================================================================================================

void CosineLawVelocity( double *NormVel , double *ParaVelX , double *ParaVelY , double CosineLawCoef , double Velocity ){
	double		AzimuthalAng , PolarAng ;
	
	AzimuthalAng	= 6.283185308*Randn() ;
	PolarAng	= acos(pow( 1.-Randn() , 1./(CosineLawCoef+1.) )) ;
	
	(*NormVel)	= Velocity * cos( PolarAng ) ;
	(*ParaVelX)	= Velocity * sin( PolarAng ) * cos( AzimuthalAng ) ;
	(*ParaVelY)	= Velocity * sin( PolarAng ) * sin( AzimuthalAng ) ;
}

//==============================================================================================================

double RotationalEnergy( double Temp , int RotDOF ){
	double		RotEnergy = 0. ;
	double		a , b , ERM ;

	if ( RotDOF == 2 ){
		// For 2 degrees of freedom, the sampling is directly from eqn (11.22).
		RotEnergy	=	-log(Randn())*BOLTZ*Temp ;
	}else{
		// Otherwise apply the acceptance-rejection method to eqn (11.23).
		a = 0.5*RotDOF - 1. ;
		do{
			// The cut-off internal energy is 10 kT.
			ERM = Randn()*10. ;
			b = (pow((ERM/a),a)) * exp(a-ERM) ;
		}while ( b < Randn() ) ;

		RotEnergy = ERM*BOLTZ*Temp ;
	}
	
	return	RotEnergy ;
}

//==============================================================================================================

double VibrationalEnergy( double Temp , int *_ParticleVibLevel , double *_ParticleEffTemp , double VibDOF , DSMC_DOMAIN *h_pDomain , double VibTemp ){
	double		VibEnergy ;
	double		a , b , ERM ;
	
	VibEnergy	= 0. ;
	
	if ( h_pDomain->VibrationalModel == 1 ){
		if ( VibDOF == 2 ){
			// For 2 degrees of freedom, the sampling is directly from eqn (11.22).
			VibEnergy	= -log(Randn())*BOLTZ*Temp ;
		}else{
			// Otherwise apply the acceptance-rejection method to eqn (11.23).
			a	= 0.5*VibDOF-1. ;
			
			do{
				// The cut-off internal energy is 10 kT.
				ERM	= Randn()*10. ;
				b	= (pow((ERM/a),a)) * exp(a-ERM) ;
			}while ( b < Randn() ) ;
			VibEnergy	= ERM*BOLTZ*Temp ;	
		}
		(*_ParticleVibLevel)	= 0 ;
		(*_ParticleEffTemp)	= Temp ;
		
	}else if ( h_pDomain->VibrationalModel == 2 ){
		// Eqn. (11.24) and (5.57) are used to set the vibrational level and energy.
		(*_ParticleVibLevel)	= -log(Randn())*Temp/VibTemp ;
		VibEnergy		= (*_ParticleVibLevel)*BOLTZ*VibTemp ;
	}
	
	return	VibEnergy ;
}

//==============================================================================================================

int FindMaxMin( int *Min , int *Max , int *Value , int Num ){
	int		Average ;
	
	
	(*Min)	= Value[0] ;
	(*Max)	= Value[0] ;
	Average	= Value[0] ;
	
	
	for ( int i=1 ; i<Num ; i++ ){
		if ( (*Min) > Value[i] ) (*Min) = Value[i] ;
		if ( (*Max) < Value[i] ) (*Max) = Value[i] ;
		
		Average	+= Value[i] ;	
	}
	Average	/= Num ;
	
	return	Average ;
}

//==============================================================================================================

double Randn(){
	double	random ;

	random	= (double)rand()/(double)RAND_MAX ;

	return random*(1. - 1.e-10) + 1.e-10 ;
}

//==============================================================================================================

double ErrorFunction( double S ){
	double	error_x ;
	double	B , C , D , T ;

	B	=	fabs(S) ;
	if ( B > 4. ){
		D	=	1. ;
	}else{
		C	=	exp(-B*B) ;
		T	=	1./(1.+0.3275911*B) ;
		D	=	1. - (0.254829592*T - 0.284496736*T*T + 1.421413741*T*T*T - 1.453152027*T*T*T*T + 1.061405429*T*T*T*T*T )*C ;
	}
	if ( S < 0. ) D = -D ;
	error_x = D ;
	
	return	error_x ;	
}

//==============================================================================================================

double GammaFunction( double x ){
	double		aa , y , gamma ;

	aa	= 1. ;
	y	= x ;

	if( y < 1. ){
		aa	= aa/y ;
	}else{
		y	= y - 1. ;
		while( y >= 1. ){
			aa	= aa * y ;
			y	= y - 1. ;
		}
	}

	gamma	= aa * (1.-0.5748646*y + 0.9512363*y*y - 0.6998588*y*y*y + 0.4245549*y*y*y*y - 0.1010678*y*y*y*y*y);

	return	gamma ;
}

//==============================================================================================================

string IntToString( int N ){
	char	Char_N[32] ;
	string	A ;	
	
	sprintf( Char_N , "%d" , N ) ;
	
	A	= Char_N ;
	
	return	A ;
}

//==============================================================================================================

void SetInitValue( int *Array , int Value , int Num ){
	for ( int i=0 ; i<Num ; i++ )
		Array[i]	= Value ;
}

//==============================================================================================================


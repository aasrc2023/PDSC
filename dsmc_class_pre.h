#include <string>
#include <vector>

using namespace std ;

#if !defined(__DSMC_CLASS_PRE_H)
#define __DSMC_CLASS_PRE_H

//==============================================================================================================

class DSMC_PRE_NODE{
	public:
		int		Id ;
		double		XCoord , YCoord , ZCoord ;

		DSMC_PRE_NODE() ;
} ;

//==============================================================================================================

class DSMC_PRE_CELL{
	public:
		int		Id , Type , ProcessorNo , Neighbor[6] , Weight ;
		DSMC_PRE_NODE	*pNode[8] ;
		double		MeanFreePath , Temp , Speed , Timestep , Weighting ,TimestepRatio , WeightingRatio ;
		int		SubcellNum ;


		DSMC_PRE_CELL() ;
		void		DumpNode( int CellNo , int NodeNum ) ;
} ;

//==============================================================================================================

class DSMC_PRE_BCFACE{
	public:
		int		Id , NodeNum , Node[4] , CountNum , NewNode[2] ;
		double		XCenter , YCenter , ZCenter ;
		double		XVel , YVel , ZVel , Temp ;
		vector<double>	NumDen ;
		

		DSMC_PRE_BCFACE() ;
		~DSMC_PRE_BCFACE() ;
} ;

//==============================================================================================================

class DSMC_PRE_BCTYPE{
	public:
		int		Type , FaceNum ;
		string		TypeName ;
		double		XVel , YVel , ZVel , NumDen , Temp , CosineLawCoef ;

		DSMC_PRE_BCTYPE() ;
		void		Dump( int BCTypeNo ) ;
		void		Dump() ;
} ;

//==============================================================================================================

class DSMC_PRE_NODECELL{
	public:
		vector<int>	CellRelation , FaceRelation ;

		~DSMC_PRE_NODECELL() ;
		void		Init() ;
} ;

//==============================================================================================================

class DSMC_PRE_GRAPHCELL{
	public:	
		int		Check[6] , Neighbor[6] , Num ;
		
		DSMC_PRE_GRAPHCELL() ;
} ;

//==============================================================================================================

class DSMC_PRE_PROCESSOR{
	public:	
		vector<DSMC_PRE_CELL>	Cell ;
		vector<int>		Neighbor ;

		~DSMC_PRE_PROCESSOR() ;
} ;

//==============================================================================================================

class DSMC_PRE_INLET{
	public:
		double		XCoord , YCoord , ZCoord , XVel , YVel , ZVel , Temp ;
		double		*NumDen ;
		
		DSMC_PRE_INLET() ;
		~DSMC_PRE_INLET() ;
		void		AllocateMemory( int Num ) ;
} ;

//==============================================================================================================

class DSMC_PRE_RESULT{
	public:
		double		XCenter , YCenter , ZCenter ;
		double		NumDensity , Density , XVel , YVel , ZVel , Temp ;
		double		TransTemp , RotTemp , VibTemp ;
		double		AveParticleNum , MeanCollSpacingMeanFreePath ;
		int		ProcessorNo ;

		DSMC_PRE_RESULT() ;
} ;

//==============================================================================================================

class DSMC_PRE_OPTION{
	public:
		int		OpenScale , OpenInputInlet , OpenDomainRedecomposition , OpenCell ;
		
		int		ProcessorNum ;
		double		Scale ;
		string		InletFileName ;
		string		InitResultFileName ;
		string		InitMeshFileName ;
		string		CellFileName ;
		
		
		DSMC_PRE_OPTION() ;
		void		Dump() ;
} ;

//==============================================================================================================

class CELLMAPPING{
	public:
		int		NodeNum[6] , SurfaceNum[6] , TetraCellNum[6] ; 
		int		SurfaceNode[6][6] , Node[6][6][4] , TetraCell[6][5][4] ;

		CELLMAPPING() ;
} ;


#endif

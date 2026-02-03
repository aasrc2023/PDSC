#include <mpi.h>
#include <string>
#include <vector>
#include <fstream>
#include <sys/time.h>

using namespace std ;

#if !defined(__DSMC_CLASS_H)
#define __DSMC_CLASS_H

//==============================================================================================================

class DSMC_NODE{
	public:
		int		Id ;
		double		XCoord , YCoord , ZCoord ;

		DSMC_NODE() ;
		void		Dump( int No ) ;
} ;

//==============================================================================================================

class DSMC_CELL{
	public:
		int		Id , LocalId , Type , Node[8] , Neighbor[6] , Surface[6] , ProcessorNo ;
		double		XCenter , YCenter , ZCenter , Volume , Timestep , Weighting , CharacteristicLength ;
		double		FaceFA[6] , FaceFB[6] , FaceFC[6] , FaceFD[6] ;
		double		MaxXCoord , MinXCoord , MaxYCoord , MinYCoord , MaxZCoord , MinZCoord ;
		
		double		SelectNum , CollNum , CollDistance , CollRatio , MeanFreePath , MeanCollisionTime ;
		double		Temp , Speed , InitMeanFreePath , InitTemp , InitSpeed , InitTimestep , InitWeighting ;
		double		AveParticleNum , MeanCollSpacingMeanFreePath ;
		int		SubcellNum ;

		DSMC_CELL() ;
		void		Dump( int No ) ;
		//void		DumpNodeCoord( DSMC_NODE *h_Node ) ;
} ;

//==============================================================================================================

class DSMC_INLET{
	public:
		int		Id , CellNo , FaceNo , NodeNum , Node[4] , SpeciesNo ;
		double		Area , XVel , YVel , ZVel , Temp , CosineLawCoef ;
		double		*NumDen ;

		DSMC_INLET() ;
		~DSMC_INLET() ;
		void		AllocateMemory( int SpeciesNum ) ;
		void		DeleteMemory() ;
		void		Dump( int No ) ;
} ;

//==============================================================================================================

class DSMC_WALLTYPE{
	public:
		int		Id , WallNo , Type ;
		double		Temp , XVel , YVel , ZVel , DiffRatio , StickingCoef, AlphaN , SigmaT , EintCoef ;
		
		DSMC_WALLTYPE() ;
		void		Dump( int No ) ;
} ;

//==============================================================================================================

class DSMC_SPECIES{
	public:
		int		Id , GroupNo , RotModel , VibMode , VibModel ;
		int   *ReactionColClass ;
		double		Fraction , Diameter , RefTemp , VisTempIndex , VSSParameter , Mass ;
		double		RotDOF , *RotRelaxationNum ;
		double		VibDOF , *VibConstant1 , *VibConstant2 , VibTemp , DisTemp ;
		
		
		DSMC_SPECIES() ;
		~DSMC_SPECIES() ;
		void		AllocateMemory( int SpeciesNum ) ;
		void		Dump( int No ) ;
} ;

//==============================================================================================================

class DSMC_MULTISPECIES{
	public:
		double		CrossSection , RefTemp , VisTempIndex , VSSParameter , ReduceMass , GammaValue ;
		
		DSMC_MULTISPECIES() ;
} ;

//==============================================================================================================
// Start Ming-Chung Lo

class DSMC_CHEMICALTCE{
	public:		
		int		*ReactionType , ReactionClassNum ,*SuccessReaction;
		
		int		*Pre1Species , *Pre2Species , *Pre3Species , *Post1Species , *Post2Species , *Post3Species ;
		double		*ArrheniusConstant , *ArrheniusTempExp , *ArrheniusActiveEnergy , *HeatofReaction ;
		
		DSMC_CHEMICALTCE() ;
		~DSMC_CHEMICALTCE() ;
		
		void		AllocateMemory( int SpeciesNum ) ;
		void		DeleteMemory( int SpeciesNum ) ;
		void		InitValue( int SpeciesNum  ) ;
				
		// Print information.
		void		Dump( int No ) ;
} ;

// End Ming-Chung Lo
//==============================================================================================================

class DSMC_SURFACE{
	public:
		int		WallNo , CellNo , FaceNo ;
		double		XCenter , YCenter , ZCenter , Area ;
		
		DSMC_SURFACE() ;
		void		Dump( int No ) ;
} ;

//==============================================================================================================

class DSMC_DSMC{
	public:
		// Particle information.
		int		*ParticleCellNo , *ParticleSpeciesNo , *ParticleLastCollide , *ParticleBeforeCellNo ;
		double		*ParticleXCoord , *ParticleYCoord , *ParticleZCoord ;
		double		*ParticleXVel , *ParticleYVel , *ParticleZVel ;
		double		*ParticleRotation , *ParticleVibration , *ParticleEffTemp , *ParticleTimestep ;
		int		*ParticleVibLevel ;


		// Index information [SpeciesGroupNum][CellNum].
		int		**IndexCell1 , **IndexCell2 , *IndexParticle , *ReactionIndex	;


		// Sample information for flow field [Species-i][CellNum-j].
		double		**SampleParticleNum ;
		double		**SampleXVel , **SampleYVel , **SampleZVel ;
		double		**SampleXVelSq , **SampleYVelSq , **SampleZVelSq ;
		double		**SampleRotation , **SampleVibration ;
		double		**SampleVibLevel , **SampleVibGroundNum , **SampleVibLevelNum ;
		
		
		// Sample information for surface [WallFaceNum-i][SpecieNum-j].
		double		**SampleSurfaceParticleNum , **SampleSurfaceStickingParticleNum ;
		double		**SampleSurfaceInNormMomentum , **SampleSurfaceInXMomentum , **SampleSurfaceInYMomentum , **SampleSurfaceInZMomentum ;
		double		**SampleSurfaceReNormMomentum , **SampleSurfaceReXMomentum , **SampleSurfaceReYMomentum , **SampleSurfaceReZMomentum ; 
		double		**SampleSurfaceInTransEng , **SampleSurfaceInRotEng , **SampleSurfaceInVibEng ;
		double		**SampleSurfaceReTransEng , **SampleSurfaceReRotEng , **SampleSurfaceReVibEng ;
		


		// Collision information [CellNum-i][Group-j][Group-k].
		double		***MaxCrossSectionSpeed , ***TotalCrossSectionSpeed ,***RemainderCollisionPair ;


		// Inlet information for enter particles. [InletFaceNum-i][SpeciesNum-j]
		double		**InletEnterNum , **InletRemainderEnterNum ;
		
		
		// Species information for vibrational model (effective number of vibrational degree of freedom). [CellNum-i][SpeciesNum-j].
		double		**EffVibDOF ;
		
		
		// Sample physical computing time [CellNum].
		double		*SamplingTimeEnd , *SamplingTimeInit , *SampleCellVibLevel ;
		
		
		// Sample collision number for Species.
		double		**CollisionNum ;


		DSMC_DSMC() ;
		~DSMC_DSMC() ;
		void		AllocateMemory( int MaxParticleNum , int CellNum , int WallFaceNum , int InletFaceNum , int SpeciesNum , int SpeciesGroupNum ) ;
		void		DeleteMemory( int CellNum , int WallFaceNum , int InletFaceNum , int SpeciesNum , int SpeciesGroupNum ) ;
		void		InitValue( int MaxParticleNum , int CellNum , int WallFaceNum , int InletFaceNum , int SpeciesNum , int SpeciesGroupNum ) ;
		
		// Print information.
		void		ParticleDump( int No ) ;
} ;

//==============================================================================================================

class DSMC_RESULT{
	public:
		double		*NumDensity , *Density , *XVel , *YVel , *ZVel , *Temp ;
		double		*TransTemp , *RotTemp , *VibTemp , *TransXTemp , *TransYTemp , *TransZTemp ;
		double		*AveParticleNum , *MeanCollSpacingMeanFreePath ;

		// [SpeciesNum][CellNum]
		double		**NumDensitySpecies , **DensitySpecies , **XVelSpecies , **YVelSpecies , **ZVelSpecies , **TempSpecies ;
		double		**TransTempSpecies , **TransXTempSpecies , **TransYTempSpecies , **TransZTempSpecies , **RotTempSpecies , **VibTempSpecies ,**AveParticleNumSpecies;
		
		
		DSMC_RESULT() ;
		~DSMC_RESULT() ;
		void		AllocateMemory( int CellNum , int SpeciesNum ) ;
		void		DeleteMemory( int SpeciesNum ) ;
		void		InitValue( int CellNum , int SpeciesNum , double Temp ) ;
} ;

//==============================================================================================================

class DSMC_DOMAIN{
	public:
		// These Values are from Input.txt
		int		Dimension , NodeNum , CellNum , LocalCellNum , TotalCellNum ; 
		int		InletFaceNum , WallFaceNum , WallTypeNum , MaxParticleNum , TransferParticleNum , InletSpecifiedNumDenNum ;
		int		SamplingFrequency , OutputFrequency , SamplingTime , TotalTimestep ;
		double		TimestepRatio , WeightingRatio , XVel , YVel , ZVel , NumDen , Temp ;
		int		SpeciesNum , SpeciesGroupNum ;
		int		VariableTimestepScheme , SubcellModel , DynamicDomainDecompisition , VibrationalModel , SimulationStage ;
		
		int		Convergence ;
		double		ParticleNumCheck , SpeedCheck , TempCheck ;

		// These values are depended on variable of flowfield.
		int		TimestepNo , ParticleNum , SamplingNum , EnterParticleNum ;
		double		MinVolume , MinLength , MinCenter , Timestep , ParticleWeighting , AxisAdjustFactor ;

		// Debug information.
		int		TrackingNum , ErrorTrackingNum , ErrorTrackingNumNew ;
		
		// These values are depended on variable of Chemical.		
		int   ChemicalNum ;


		DSMC_DOMAIN() ;
		void		Dump() ;
		void		DumpFile( ofstream &OutputInfo ) ;
} ;

//==============================================================================================================

class DSMC_CONVERGENCE{
	public:
		int		DataNum , CurrentNo ;
		double		ConditionParticleNum , ConditionSpeed , ConditionTemp ;
		double		*ParticleNum , *Speed , *Temp ;
		double		OldParticleNum , OldSpeed , OldTemp , DiffParticleNum , DiffSpeed , DiffTemp ;
	
	
		DSMC_CONVERGENCE() ;
		~DSMC_CONVERGENCE() ;
		void		AllocateMemory() ;
		void		InitValue( double ParticleNumCheck , double SpeedCheck , double TempCheck ) ;
		bool		Convergence() ;
		void		Sum( double AveParticleNum , double AveSpeed , double AveTemp , int Num ) ;
		void		InitSum() ;
		void 		Average() ;
} ;

//==============================================================================================================

class DSMC_PROCESSOR{
	public:	
		int		NeighborNum , Neighbor[64] , MaxNeighborNum ;
		int		*LocalCellNo , *CellProcessorNo ;
		
		DSMC_PROCESSOR() ;
		void		AllocateMemory( int TotalCellNum ) ;
		void		InitValue( int TotalCellNum ) ;
		void		DeleteMemory() ;
		void		Dump( int ID , int CellNum ) ;
} ;

//==============================================================================================================

class DSMC_MPI_PARTICLE{
	public:
		int		ProcessorNo , GlobalCellNo , BeforeCellNo , SpeciesNo , LastCollide , VibLevel ,ReactionIndex ;
		double		XCoord , YCoord , ZCoord ;
		double		XVel , YVel , ZVel ;
		double		Rotation , Vibration , EffTemp ;
		double		Timestep , BeforeTimestep ;
		
		DSMC_MPI_PARTICLE() ;
} ;

//==============================================================================================================

class DSMC_MPI_DATATYPE{
	public:
		// DSMC_DOMAIN
		int		DomainVariableNum , DomainVariableLen[47] ;
		MPI_Datatype	DomainOldType[47] , MPI_DOMAIN ;
		MPI_Aint	DomainDisp[47] ;

		void InitMPIDataType( DSMC_DOMAIN *pDomain ) ;
		
		
		
		// DSMC_MPI_PARTICLE
		int		ParticleVariableNum , ParticleVariableLen[18] ;
		MPI_Datatype	ParticleOldType[18] , MPI_PARTICLE ;
		MPI_Aint	ParticleDisp[18] ;
		
		void InitMPIDataType( DSMC_MPI_PARTICLE *pMPIParticle ) ;
} ;

//==============================================================================================================

class CELLMAPPING{
	public:
		int		NodeNum[6] , SurfaceNum[6] , TetraCellNum[6] ; 
		int		SurfaceNode[6][6] , Node[6][6][4] , TetraCell[6][5][4] ;

		CELLMAPPING() ;
} ;

//==============================================================================================================
//==============================================================================================================
// For DSMC Processor
//==============================================================================================================
//==============================================================================================================

class DSMC_POST_SURFACE{
	public:
		int		Id , CellNo , WallNo , FaceNo ;
		double		SurfacePro[17] , Center[4] , Area , SimTime ;
		double		DepositionRate , Cp , Cf , Ch , Pressure , HeatFlux ;
		
		DSMC_POST_SURFACE() ;
		void		 Dump( int No ) ;
} ;

//==============================================================================================================

class DSMC_POST_SURFACE_CONDITION{
	public:
		string		SortName ;
		int		WallNum , *WallNo ;		
		double		Density , Velocity , Mass , Drag , HeatFlux , CoeffAxial , CoeffNomal ;
		
		DSMC_POST_SURFACE_CONDITION( string Filename ) ;
		~DSMC_POST_SURFACE_CONDITION() ;
		void		Dump( int WallFaceNum ) ;
} ;

//==============================================================================================================

class DSMC_POST_NODECELL{
	public:
		vector<int>	CellRelation ;

		~DSMC_POST_NODECELL() ;
} ;

//==============================================================================================================

class DSMC_TIME{
	private:
		struct timeval	start , end ;
		long double	time , tstart ;

	public:
		long double	Total , Init , Move , Index , Collision , Sample , CalResult ;
		
		DSMC_TIME() ;
		void		Start() ;
		void		End() ;
		void		Time() ;
		void		Time( long double *t ) ;
		void		PrintFile() ;
} ;

//==============================================================================================================

#endif

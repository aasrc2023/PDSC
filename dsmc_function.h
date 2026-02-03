#include <fstream>
#include "dsmc_class.h"



#if !defined(__DSMC_FUNCTION_H)
#define __DSMC_FUNCTION_H


//==============================================================================================================
//==============================================================================================================

void ParticleMovement(	DSMC_DOMAIN		*h_pDomain ,
			DSMC_PROCESSOR		*h_pProcessor ,
			DSMC_NODE		*h_Node ,
			DSMC_CELL		*h_Cell ,
			DSMC_INLET		*h_Inlet ,
			DSMC_WALLTYPE		*h_WallType ,
			DSMC_SPECIES		*h_Species ,
			DSMC_SURFACE		*h_Surface ,
			DSMC_DSMC		*h_pDSMC ,
			CELLMAPPING		*h_pMapping , 
			DSMC_MPI_PARTICLE	*h_MPIParticleOut , 
			DSMC_MPI_PARTICLE	*h_MPIParticleIn ,
			DSMC_MPI_DATATYPE	*pMPIDataType , 
			ofstream		&OutputDebug ) ;

//==============================================================================================================

void MoveAllParticle2D(	DSMC_DOMAIN		*h_pDomain ,
			DSMC_PROCESSOR		*h_pProcessor ,
			DSMC_NODE		*h_Node ,
			DSMC_CELL		*h_Cell ,
			DSMC_WALLTYPE		*h_WallType ,
			DSMC_SPECIES		*h_Species ,
			DSMC_SURFACE		*h_Surface ,
			DSMC_DSMC		*h_pDSMC ,
			CELLMAPPING		*h_pMapping ,
			DSMC_MPI_PARTICLE	*h_MPIParticleIn ,
			int			*pMPIParticleNum , 
			ofstream		&OutputDebug ) ;
			

void MoveAllParticle3D(	DSMC_DOMAIN		*h_pDomain ,
			DSMC_PROCESSOR		*h_pProcessor ,
			DSMC_NODE		*h_Node ,
			DSMC_CELL		*h_Cell ,
			DSMC_WALLTYPE		*h_WallType ,
			DSMC_SPECIES		*h_Species ,
			DSMC_SURFACE		*h_Surface ,
			DSMC_DSMC		*h_pDSMC ,
			CELLMAPPING		*h_pMapping ,
			DSMC_MPI_PARTICLE	*h_MPIParticleIn ,
			int			*pMPIParticleNum , 
			ofstream		&OutputDebug ) ;
			
			
void MoveAllParticleAxisymmetric(	DSMC_DOMAIN		*h_pDomain ,
					DSMC_PROCESSOR		*h_pProcessor ,
					DSMC_NODE		*h_Node ,
					DSMC_CELL		*h_Cell ,
					DSMC_WALLTYPE		*h_WallType ,
					DSMC_SPECIES		*h_Species ,
					DSMC_SURFACE		*h_Surface ,
					DSMC_DSMC		*h_pDSMC ,
					CELLMAPPING		*h_pMapping ,
					DSMC_MPI_PARTICLE	*h_MPIParticleIn ,
					int			*pMPIParticleNum , 
					ofstream		&OutputDebug ) ;			

//==============================================================================================================
		
void MoveNewParticle2D(	DSMC_DOMAIN		*h_pDomain ,
			DSMC_PROCESSOR		*h_pProcessor ,
			DSMC_NODE		*h_Node ,
			DSMC_CELL		*h_Cell ,
			DSMC_WALLTYPE		*h_WallType ,
			DSMC_SPECIES		*h_Species ,
			DSMC_SURFACE		*h_Surface ,
			DSMC_DSMC		*h_pDSMC ,
			CELLMAPPING		*h_pMapping , 
			int 			ParticleNo ,
			DSMC_MPI_PARTICLE	*h_MPIParticleIn ,
			int			*pMPIParticleNum , 
			ofstream		&OutputDebug ) ;
			
			
void MoveNewParticle3D(	DSMC_DOMAIN		*h_pDomain ,
			DSMC_PROCESSOR		*h_pProcessor ,
			DSMC_NODE		*h_Node ,
			DSMC_CELL		*h_Cell ,
			DSMC_WALLTYPE		*h_WallType ,
			DSMC_SPECIES		*h_Species ,
			DSMC_SURFACE		*h_Surface ,
			DSMC_DSMC		*h_pDSMC ,
			CELLMAPPING		*h_pMapping , 
			int 			ParticleNo ,
			DSMC_MPI_PARTICLE	*h_MPIParticleIn ,
			int			*pMPIParticleNum , 
			ofstream		&OutputDebug ) ;
			

void MoveNewParticleAxisymmetric(	DSMC_DOMAIN		*h_pDomain ,
					DSMC_PROCESSOR		*h_pProcessor ,
					DSMC_NODE		*h_Node ,
					DSMC_CELL		*h_Cell ,
					DSMC_WALLTYPE		*h_WallType ,
					DSMC_SPECIES		*h_Species ,
					DSMC_SURFACE		*h_Surface ,
					DSMC_DSMC		*h_pDSMC ,
					CELLMAPPING		*h_pMapping , 
					int 			ParticleNo ,
					DSMC_MPI_PARTICLE	*h_MPIParticleIn ,
					int			*pMPIParticleNum , 
					ofstream		&OutputDebug ) ;

//==============================================================================================================

void MoveOtherParticle2D(	DSMC_DOMAIN		*h_pDomain ,
				DSMC_PROCESSOR		*h_pProcessor ,
				DSMC_NODE		*h_Node ,
				DSMC_CELL		*h_Cell ,
				DSMC_WALLTYPE		*h_WallType ,
				DSMC_SPECIES		*h_Species ,
				DSMC_SURFACE		*h_Surface ,
				DSMC_DSMC		*h_pDSMC ,
				CELLMAPPING		*h_pMapping , 
				DSMC_MPI_PARTICLE	*h_MPIParticleOut , 
				DSMC_MPI_PARTICLE	*h_MPIParticleIn ,
				int			*MPIParticleNumOut ,
				int			*MPIParticleNumIn ,
				int			*pMPIParticleNun ,
				DSMC_MPI_DATATYPE	*pMPIDataType ,
				ofstream		&OutputDebug ) ;
				
				
void MoveOtherParticle3D(	DSMC_DOMAIN		*h_pDomain ,
				DSMC_PROCESSOR		*h_pProcessor ,
				DSMC_NODE		*h_Node ,
				DSMC_CELL		*h_Cell ,
				DSMC_WALLTYPE		*h_WallType ,
				DSMC_SPECIES		*h_Species ,
				DSMC_SURFACE		*h_Surface ,
				DSMC_DSMC		*h_pDSMC ,
				CELLMAPPING		*h_pMapping , 
				DSMC_MPI_PARTICLE	*h_MPIParticleOut , 
				DSMC_MPI_PARTICLE	*h_MPIParticleIn ,
				int			*MPIParticleNumOut ,
				int			*MPIParticleNumIn ,
				int			*pMPIParticleNum ,
				DSMC_MPI_DATATYPE	*pMPIDataType ,
				ofstream		&OutputDebug ) ;
				

void MoveOtherParticleAxisymmetric(	DSMC_DOMAIN		*h_pDomain ,
					DSMC_PROCESSOR		*h_pProcessor ,
					DSMC_NODE		*h_Node ,
					DSMC_CELL		*h_Cell ,
					DSMC_WALLTYPE		*h_WallType ,
					DSMC_SPECIES		*h_Species ,
					DSMC_SURFACE		*h_Surface ,
					DSMC_DSMC		*h_pDSMC ,
					CELLMAPPING		*h_pMapping , 
					DSMC_MPI_PARTICLE	*h_MPIParticleOut , 
					DSMC_MPI_PARTICLE	*h_MPIParticleIn ,
					int			*MPIParticleNumOut ,
					int			*MPIParticleNumIn ,
					int			*pMPIParticleNum ,
					DSMC_MPI_DATATYPE	*pMPIDataType ,
					ofstream		&OutputDebug ) ;
				
//==============================================================================================================

void MoveParticle2D( double *pXCoord , double *pYCoord , double XVel , double YVel , double Time ) ;

void MoveParticle3D( double *pXCoord , double *pYCoord , double *pZCoord , double XVel , double YVel , double ZVel , double Time ) ;

double MoveParticleAxisymmetric( double *pXCoord , double *pYCoord , double *pXVel , double *pYVel , double *pZVel , double Time ) ;

void MoveParticleAxisymmetric( double *pXCoord , double *pYCoord , double XVel , double YVel , double ZVel , double Time ) ;

//==============================================================================================================

void AdjustVelocityAxisymmetric( double YCoord , double *pYVel , double *pZVel , double Time ) ;

//==============================================================================================================

void CalculateTimeCollideFace2D( double		*pCollideTime , 
				 int		*pFaceNo ,
				 double		Time , 
				 double		XCoord , 
				 double		YCoord , 
				 double		BeforeXCoord , 
				 double		BeforeYCoord ,
				 int		BeforeCellNo , 
				 DSMC_CELL	*_Cell ,
				 CELLMAPPING	*h_pMapping ,
				 int		Debug ) ;
				 
				 
void CalculateTimeCollideFace3D( double		*pCollideTime , 
				 int		*pFaceNo ,
				 double		Time , 
				 double		XCoord , 
				 double		YCoord , 
				 double		ZCoord ,
				 double		BeforeXCoord , 
				 double		BeforeYCoord , 
				 double		BeforeZCoord ,
				 int		BeforeCellNo ,
				 DSMC_CELL	*_Cell ,
				 CELLMAPPING	*h_pMapping ,
				 int		Debug ) ;

//==============================================================================================================

bool ParticleCollideSurface2D(	double		*pXVel ,
				double		*pYVel ,
				double		*pZVel ,
				DSMC_CELL	*_Cell ,
				int		FaceNo ,
				int		*pParticleNo ,
				int		*pParticleNum ,
				DSMC_SPECIES	*_pSpecies ,
				DSMC_DOMAIN	*h_pDomain ,
				DSMC_NODE	*h_Node ,
				DSMC_WALLTYPE	*h_WallType ,
				DSMC_SURFACE	*h_Surface ,
				DSMC_DSMC	*h_pDSMC ,
				CELLMAPPING	*h_pMapping , 
				ofstream	&OutputDebug ) ;
				
				
bool ParticleCollideSurface3D(	double		*pXVel ,
				double		*pYVel ,
				double		*pZVel ,
				DSMC_CELL	*_Cell ,
				int		FaceNo ,
				int		*pParticleNo ,
				int		*pParticleNum ,
				DSMC_SPECIES	*_pSpecies ,
				DSMC_DOMAIN	*h_pDomain ,
				DSMC_NODE	*h_Node ,
				DSMC_WALLTYPE	*h_WallType ,
				DSMC_SURFACE	*h_Surface ,
				DSMC_DSMC	*h_pDSMC ,
				CELLMAPPING	*h_pMapping ) ;

//==============================================================================================================

void SamplingSurface(	int		Method ,
			double		XVel ,
			double		YVel ,
			double		ZVel ,
			double		RotEnergy ,
			double		VibEnergy ,
			int		SurfaceNo ,
			DSMC_SPECIES	*_pSpecies ,
			DSMC_DSMC	*h_pDSMC ) ;

//==============================================================================================================

void EnterNewPaticle2D(	DSMC_DOMAIN		*h_pDomain ,
			DSMC_NODE		*h_Node ,
			DSMC_INLET		*h_Inlet ,
			DSMC_SPECIES		*h_Species ,
			DSMC_DSMC		*h_pDSMC  ) ;
			
			
void EnterNewPaticle3D(	DSMC_DOMAIN		*h_pDomain ,
			DSMC_NODE		*h_Node ,
			DSMC_CELL		*h_Cell ,
			DSMC_INLET		*h_Inlet ,
			DSMC_SPECIES		*h_Species ,
			DSMC_DSMC		*h_pDSMC ,
			CELLMAPPING		*h_pMapping ) ;
			

void EnterNewPaticleAxisymmetric(	DSMC_DOMAIN		*h_pDomain ,
					DSMC_NODE		*h_Node ,
					DSMC_INLET		*h_Inlet ,
					DSMC_SPECIES		*h_Species ,
					DSMC_DSMC		*h_pDSMC  ) ;

//==============================================================================================================

void RemoveParticle( int *pParticleNo , int *pParticleNum , DSMC_DSMC *h_pDSMC ) ;

//==============================================================================================================

/*void TransferParticleToBuffer(	DSMC_DOMAIN		*h_pDomain ,
				DSMC_PROCESSOR		*h_pProcessor , 
				DSMC_DSMC		*h_pDSMC , 
				int			*pParticleNo , 
				int			*pParticleNum , 
				DSMC_MPI_PARTICLE	*h_MPIParticleIn , 
				int			*pMPIParticleNum , 
				int			ProcessorNo ,
				int			GlobalCellNo ,
				int			BeforeCellNo ,
				double			Timestep , 
				double			BeforeTimestep ) ;*/
				
void TransferParticleToBuffer(	DSMC_DSMC		*h_pDSMC , 
				int			*pParticleNo , 
				int			*pParticleNum , 
				DSMC_MPI_PARTICLE	*h_MPIParticleIn , 
				int			*pMPIParticleNum , 
				int			ProcessorNo ,
				int			GlobalCellNo ,
				int			BeforeCellNo ,
				double			Timestep , 
				double			BeforeTimestep ) ;
				
//==============================================================================================================
//==============================================================================================================

int AddParticleFromOtherProcessor(	DSMC_DOMAIN		*h_pDomain ,
					DSMC_PROCESSOR		*h_pProcessor , 
					DSMC_DSMC		*h_pDSMC ,
					DSMC_CELL		*h_Cell ,
					DSMC_MPI_PARTICLE	*h_MPIParticleIn ,
					int			*pMPIParticleNum , 
					ofstream		&OutputDebug ) ;
					
//==============================================================================================================

void MPITransferParticle(	DSMC_DOMAIN		*h_pDomain , 
				DSMC_PROCESSOR		*h_pProcessor , 
				DSMC_MPI_PARTICLE	*h_MPIParticleOut , 
				DSMC_MPI_PARTICLE	*h_MPIParticleIn ,
				int			*MPIParticleNumOut ,
				int			*MPIParticleNumIn ,
				int			*pMPIParticleNum ,
				DSMC_MPI_DATATYPE	*pMPIDataType , 
				ofstream		&OutputDebug ) ;

//==============================================================================================================
//==============================================================================================================

void Index(	DSMC_DOMAIN		*h_pDomain ,
		DSMC_SPECIES		*h_Species ,
		DSMC_DSMC		*h_pDSMC , 
		ofstream		&OutputDebug ) ;

//==============================================================================================================
//==============================================================================================================

void Collision(	DSMC_DOMAIN		*h_pDomain ,
		DSMC_CELL		*h_Cell ,
		DSMC_SPECIES		*h_Species ,
		DSMC_DSMC		*h_pDSMC ,
		DSMC_RESULT		*h_pResult ,
		DSMC_CHEMICALTCE *h_ChemicalTCE , 
		ofstream		&OutputDebug ) ;
	
//==============================================================================================================
		
double SelectParticle(	DSMC_SPECIES		*h_Species ,
			DSMC_DSMC		*h_pDSMC ,
			double			*RelativeVel , 
			double			*pRelativeVelSq ,
			double			*pRelativeSpeed ,
			int			*SelectParticleNo ,
			int			IndexCell1_Group1 ,
			int			IndexCell2_Group1 ,
			int			IndexCell1_Group2 ,
			int			IndexCell2_Group2 ,
			DSMC_MULTISPECIES	*pMultiSpecies ) ;
			
//==============================================================================================================

int GetSubCellNum( DSMC_CELL *_Cell , int ParticleNum , int SubcellModel , int Dimension ) ;

//==============================================================================================================

void IndexSubCell2D( 	DSMC_DOMAIN	*h_pDomain , 
			DSMC_DSMC	*h_pDSMC ,
			DSMC_CELL	*_Cell ,
			int		CellNo , 
			int		**IndexCellSub1 , 
			int		**IndexCellSub2 , 
			int		*IndexParticleSub , 
			int		*IndexParticleSubCellNo , 
			int		SubCellNum , 
			int 		ParticleNum ,
			int		*GroupParticleNum1,
			int     *GroupParticleNum2,			
			ofstream	&OutputDebug ) ;
			
			
void IndexSubCell3D( 	DSMC_DOMAIN	*h_pDomain , 
			DSMC_DSMC	*h_pDSMC ,
			DSMC_CELL	*_Cell ,
			int		CellNo , 
			int		**IndexCellSub1 , 
			int		**IndexCellSub2 , 
			int		*IndexParticleSub , 
			int		*IndexParticleSubCellNo , 
			int		SubCellNum , 
			int 		ParticleNum , 
			int		*GroupParticleNum1,
			int     *GroupParticleNum2,		
			ofstream	&OutputDebug ) ;

//==============================================================================================================

double SelectParticleSubCell(	DSMC_SPECIES		*h_Species ,
				DSMC_DSMC		*h_pDSMC ,
				double			*RelativeVel , 
				double			*pRelativeVelSq ,
				double			*pRelativeSpeed ,
				int			*SelectParticleNo ,
				int			IndexCell1_Group1 ,
				int			IndexCell2_Group1 ,
				int			IndexCell1_Group2 ,
				int			IndexCell2_Group2 ,
				DSMC_MULTISPECIES	*pMultiSpecies , 
				int			**IndexCellSub1 , 
				int			**IndexCellSub2 , 
				int			*IndexParticleSub ,
				int			*IndexParticleSubCellNo ,
				int			GroupNo1 ,
				int			GroupNo2 , 
				int			SubCellNum , 
				int         *GroupParticleNum1,
				int 		*SubCellNo,
				int			**ReactionBufferSub,
				int			TotalSubCellNum ) ;
				
//==============================================================================================================

void CalculateMultiSpecies( DSMC_MULTISPECIES *pMultiSpecies , DSMC_SPECIES *_pSpecies1 , DSMC_SPECIES *_pSpecies2 ) ;

//==============================================================================================================

void ChemicalTCE(	DSMC_DOMAIN		*h_pDomain ,
		DSMC_CELL		*h_Cell ,		
		DSMC_SPECIES		*h_Species ,
		DSMC_DSMC		*h_pDSMC ,
		int			*SelectParticleNo ,
		int			*SelectSpeciesNo ,
		double			*RelativeVel , 
		double			*pRelativeVelSq ,
		double			*pRelativeSpeed ,
		DSMC_MULTISPECIES	*pMultiSpecies ,
		DSMC_CHEMICALTCE *h_ChemicalTCE ,
		int			CellNo ,
		double			Temp , 
		double			TransTemp,
		int   *ReactionBuffer,
		int   **ReactionBufferSub,
		int   *SubCellNo,
		int   SubCellNum,
		ofstream	&OutputDebug ) ;
		
//==============================================================================================================

void Dissociation(	DSMC_DOMAIN		*h_pDomain ,
		DSMC_SPECIES		*h_Species ,
		DSMC_DSMC		*h_pDSMC ,
		int			*SelectParticleNo ,
		int			*SelectSpeciesNo ,
		double			*RelativeVel , 
		double			*pRelativeVelSq ,
		double			*pRelativeSpeed ,
		DSMC_MULTISPECIES	*pMultiSpecies ,
		DSMC_CHEMICALTCE	*pChemicalTCE ,
		int			CellNo ,
		double			Temp , 
		double			TransTemp,
		int 		ColClassNo,
		int 		ReactionClassNo,
		double      TotalEnergy,
		int   *ReactionBuffer,
		int   **ReactionBufferSub,
		int   *SubCellNo,
		int   SubCellNum,
		ofstream	&OutputDebug) ;

//==============================================================================================================

void Exchange(	DSMC_DOMAIN		*h_pDomain ,
		DSMC_SPECIES		*h_Species ,
		DSMC_DSMC		*h_pDSMC ,
		int			*SelectParticleNo ,
		int			*SelectSpeciesNo ,
		double			*RelativeVel , 
		double			*pRelativeVelSq ,
		double			*pRelativeSpeed ,
		DSMC_MULTISPECIES	*pMultiSpecies ,
		DSMC_CHEMICALTCE	*pChemicalTCE ,
		int			CellNo ,
		double			Temp , 
		double			TransTemp,
		int 		ColClassNo,
		int 		ReactionClassNo,
		double      TotalEnergy,
		int   *ReactionBuffer,
		int   **ReactionBufferSub,
		int   *SubCellNo,
		int   SubCellNum,
		ofstream	&OutputDebug) ;

//==============================================================================================================

void Recombination(	DSMC_DOMAIN		*h_pDomain ,
		DSMC_SPECIES		*h_Species ,
		DSMC_DSMC		*h_pDSMC ,
		int			*SelectParticleNo ,
		int			*SelectSpeciesNo ,
		double			*RelativeVel , 
		double			*pRelativeVelSq ,
		double			*pRelativeSpeed ,
		DSMC_MULTISPECIES	*pMultiSpecies ,
		DSMC_CHEMICALTCE	*pChemicalTCE ,
		int			CellNo ,
		double			Temp , 
		double			TransTemp,
		int 		ColClassNo,
		int 		ReactionClassNo,
		double      TotalEnergy,
		int ThirdBody,
		int		*pParticleNo ,
		int		*pParticleNum ,
		int   *ReactionBuffer,
		int   **ReactionBufferSub,
		int   *SubCellNo,
		int   SubCellNum,
		ofstream	&OutputDebug) ;

//==============================================================================================================
void Inelrv(	DSMC_DOMAIN		*h_pDomain ,
		DSMC_SPECIES		*h_Species ,
		DSMC_DSMC		*h_pDSMC ,
		int			*SelectParticleNo ,
		int			*SelectSpeciesNo ,
		double			*RelativeVel , 
		double			*pRelativeVelSq ,
		double			*pRelativeSpeed ,
		DSMC_MULTISPECIES	*pMultiSpecies ,
		int			CellNo ,
		double			Temp , 
		double			TransTemp,
		ofstream		&OutputDebug ) ;

//==============================================================================================================

double LarsenBorgnakkeEnergyRatio( double XMA , double XMB ) ;

//==============================================================================================================

void Elastic(	DSMC_DSMC		*h_pDSMC ,
		double			*RelativeVel , 
		double			*pRelativeSpeed ,
		int			*SelectParticleNo ,
		DSMC_MULTISPECIES	*pMultiSpecies , 
		double			Mass1 ,
		double			Mass2 ) ;

//==============================================================================================================
//==============================================================================================================

void Sample(	DSMC_DOMAIN		*h_pDomain ,
		DSMC_SPECIES		*h_Species ,
		DSMC_DSMC		*h_pDSMC , 
		DSMC_CELL		*h_Cell ) ;

//==============================================================================================================

void SampleInit( DSMC_DOMAIN *h_pDomain , DSMC_DSMC *h_pDSMC, DSMC_CELL *h_Cell ) ;

//==============================================================================================================

void SamplePhysicalTime( DSMC_DOMAIN *h_pDomain , DSMC_CELL *h_Cell , DSMC_DSMC *h_pDSMC ) ;

//==============================================================================================================
//==============================================================================================================

void CalculateResult(	DSMC_DOMAIN		*h_pDomain ,
			DSMC_CELL		*h_Cell ,
			DSMC_SPECIES		*h_Species ,
			DSMC_DSMC		*h_pDSMC ,
			DSMC_RESULT		*h_pResult ) ;

//==============================================================================================================
//==============================================================================================================

void CheckConvergence( DSMC_DOMAIN *h_pDomain , DSMC_RESULT *h_pResult , DSMC_CONVERGENCE *h_pConvergence , int TimestepNo , ofstream &OutputConvergence ) ;

//==============================================================================================================
//==============================================================================================================

#endif

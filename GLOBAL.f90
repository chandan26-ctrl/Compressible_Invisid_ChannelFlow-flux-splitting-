!****************************************************************************                  
! *   FILE         = MULTIPLE FILES                                         *
! *   AUTHOR       = CHANDAN KUMAR (CHANDANKR@IITKGP.AC.IN)                 *
! *   INSTITUTE    = INDIAN INSTITUTE OF TECHNOLOGY (IIT), KHARAGPUR        *
!****************************************************************************
! *  THIS PROGRAM IS DISTRIBUTED IN A HOPE THAT IT WILL BE USEFUL.          * 
!**************************************************************************** 
MODULE GLOBAL
  IMPLICIT NONE
           DOUBLE PRECISION, PARAMETER :: UINF=300.0D0, PINF=101325.0D0, TINF=288.0D0, PI=4.0D0*DATAN(1.0D0) 
           DOUBLE PRECISION, PARAMETER :: XL=0.50D0, THETA=14.0632D0, RS=287.0D0, ROEEP=1.0D-3,GAMMAN=1.4D0,CFL=0.6D0
           DOUBLE PRECISION, ALLOCATABLE :: P(:,:), U(:,:), V(:,:), RHO(:,:), RHOU(:,:), RHOV(:,:), E(:,:),SOS(:,:)
           DOUBLE PRECISION, ALLOCATABLE :: DFI(:,:,:), DFJ(:,:,:), DFF(:,:,:), DFE(:,:,:), DU(:,:,:), MACH(:,:) , T(:,:)
!           DOUBLE PRECISION, ALLOCATABLE :: RHOO(:,:),RHOUO(:,:),RHOVO(:,:),EO(:,:),RHON(:,:),RHOUN(:,:),RHOVN(:,:),EN(:,:)
           DOUBLE PRECISION :: TTOT, PTOT, RHOTOT, TOUT, POUT, RHOOUT, COUT, UOUT, RHOINF,DT,DTB,VOL
           DOUBLE PRECISION :: YL1, YL2, RTHETA, DX, DY, MINF, CINF, LH, M2, M2N, FX, FXX           
           DOUBLE PRECISION :: UUJ, CC, XX, PPF, UUJP, UUJN, DDT, RATIO,DT1
           INTEGER :: I, J, K, L, N, KM, NY2, EXNY, ITER, REMAIN, II,JJ,ITER2
           INTEGER, PARAMETER :: NX=101, NY1=41, NTIMES=5000
           REAL, PARAMETER :: ERROR=0.00001
!           DOUBLE PRECISION :: MAXERR1,MAXERR2,MAXERR3,MAXERR4,ERROR1,ERROR2,ERROR3,ERROR4,FERR
END MODULE GLOBAL


!****************************************************************************                  
! *   FILE         = MULTIPLE FILES                                         *
! *   AUTHOR       = CHANDAN KUMAR (CHANDANKR@IITKGP.AC.IN)                 *
! *   INSTITUTE    = INDIAN INSTITUTE OF TECHNOLOGY (IIT), KHARAGPUR        *
!****************************************************************************
! *  THIS PROGRAM IS DISTRIBUTED IN A HOPE THAT IT WILL BE USEFUL.          * 
!****************************************************************************  
 SUBROUTINE FXI
!==========================================================================
! INITIALISE INVISID FLUX IN X DIRECTION
!==========================================================================
        USE GLOBAL
        IMPLICIT NONE

          DO I= 1,NX
          DO J= 1,NY2
              UUJ=U(I,J)
              CC = DSQRT(GAMMAN*P(I,J)/RHO(I,J))
              XX = (DSQRT((ROEEP*CC)**2+UUJ**2)+CC)
           
              PPF  = 0.5D0*P(I,J)
              UUJP = 0.5D0*(UUJ + XX)
              UUJN = 0.5D0*(UUJ - XX)

              DFI(I,J,1)=(RHO(I,J))*UUJP
              DFI(I,J,2)=(RHOU(I,J))*UUJP+PPF
              DFI(I,J,3)=(RHOV(I,J))*UUJP
              DFI(I,J,4)=(E(I,J))*UUJP+PPF*UUJ

              DFJ(I,J,1)=(RHO(I,J))*UUJN
              DFJ(I,J,2)=(RHOU(I,J))*UUJN+PPF
              DFJ(I,J,3)=(RHOV(I,J))*UUJN
              DFJ(I,J,4)=(E(I,J))*UUJN+PPF*UUJ
           END DO
           END DO

        END SUBROUTINE FXI

        SUBROUTINE DFXI
!==========================================================================
! COMPUTE DERIVATIVE FLUX IN X-DECTION
!==========================================================================
        USE GLOBAL
        IMPLICIT NONE

          DO N=1,4
          DO I= 2,NX-1
          DO J= 2,NY2-1
              DFF(I,J,N) = -(DFI(I,J,N)-DFI(I-1,J,N))/DX-(DFJ(I+1,J,N)-DFJ(I,J,N))/DX               
          END DO
          END DO
          END DO
        END SUBROUTINE DFXI


        SUBROUTINE FYI
!==========================================================================
! INITIALISE INVISID FLUX IN Y DIRECTION
!==========================================================================
        USE GLOBAL
        IMPLICIT NONE

              !________________________FLUX IN Y-DIRECTION______________________

          DO I= 1,NX
          DO J= 1,NY2
              UUJ=V(I,J)
              CC = DSQRT(GAMMAN*P(I,J)/RHO(I,J))
              XX = (DSQRT((ROEEP*CC)**2+UUJ**2)+CC)
           
              PPF  = 0.5D0*P(I,J)
              UUJP = 0.5D0*(UUJ + XX)
              UUJN = 0.5D0*(UUJ - XX)

              DFI(I,J,1)=(RHO(I,J))*UUJP
              DFI(I,J,2)=RHOU(I,J)*UUJP
              DFI(I,J,3)=RHOV(I,J)*UUJP+PPF
              DFI(I,J,4)=E(I,J)*UUJP+PPF*UUJ

              DFJ(I,J,1)=RHO(I,J)*UUJN
              DFJ(I,J,2)=RHOU(I,J)*UUJN
              DFJ(I,J,3)=RHOV(I,J)*UUJN+PPF
              DFJ(I,J,4)=E(I,J)*UUJN+PPF*UUJ
           END DO
           END DO


        END SUBROUTINE FYI

        SUBROUTINE DFYI
!==========================================================================
! COMPUTE DERIVATIVE FLUX IN Y-DECTION 
!==========================================================================
        USE GLOBAL
        IMPLICIT NONE


          DO N=1,4
          DO I= 2,NX-1
          DO J= 2,NY2-1
              DFE(I,J,N) = (-(DFI(I,J,N)-DFI(I,J-1,N))/DY-(DFJ(I,J+1,N)-DFJ(I,J,N))/DY)
          END DO
          END DO
          END DO
          
        END SUBROUTINE DFYI


        SUBROUTINE DUT
!==========================================================================
! CONTRIBUTION OF FLUXES TO DU
!==========================================================================
        USE GLOBAL
        IMPLICIT NONE

           
          DO N=1,4
          DO I= 2,NX-1
          DO J= 2,NY2-1
              DU(I,J,N)=DFF(I,J,N)+DFE(I,J,N)
           END DO
           END DO
           END DO

        END SUBROUTINE DUT

        SUBROUTINE NEW_U
!==========================================================================
! COMPUTE VARIABLE AT N+1 EXPLICITLY
!==========================================================================
        USE GLOBAL
        IMPLICIT NONE

          DO I= 2,NX-1
          DO J= 2,NY2-1
              RHO(I,J) = RHO(I,J)+DU(I,J,1)*DT
              RHOU(I,J)= RHOU(I,J)+DU(I,J,2)*DT
              RHOV(I,J)= RHOV(I,J)+DU(I,J,3)*DT
              E(I,J)   = E(I,J)+DU(I,J,4)*DT
           END DO
           END DO


        END SUBROUTINE NEW_U


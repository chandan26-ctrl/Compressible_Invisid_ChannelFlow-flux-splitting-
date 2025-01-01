!****************************************************************************                  
! *   FILE         = MULTIPLE FILES                                         *
! *   AUTHOR       = CHANDAN KUMAR (CHANDANKR@IITKGP.AC.IN)                 *
! *   INSTITUTE    = INDIAN INSTITUTE OF TECHNOLOGY (IIT), KHARAGPUR        *
!****************************************************************************
! *  THIS PROGRAM IS DISTRIBUTED IN A HOPE THAT IT WILL BE USEFUL.          * 
!****************************************************************************  

        SUBROUTINE INITIAL
!==========================================================================
! INTIALISE INITIAL AND INLET CONDITION
!==========================================================================
        USE GLOBAL
        IMPLICIT NONE



!____________________________INITIAL CONDITION_______________________________

       DO I= 1,NX
        DO J= 1,NY2
           U(I,J)=0.0D0
           V(I,J)=0.0D0
           P(I,J)=PINF
           RHO(I,J)=RHOINF
           RHOU(I,J)=RHO(I,J)*U(I,J)
           RHOV(I,J)=RHO(I,J)*V(I,J)
           E(I,J)=P(I,J)/(GAMMAN-1.0D0)+0.5D0*RHO(I,J)*(U(I,J)**2+V(I,J)**2)
        END DO
        END DO
     

!________________________INLET BOUNDARY CONDITION____________________________

        DO J=1,NY2         
           U(1,J)=UINF
           V(1,J)=0.0D0
           P(1,J)=PINF
           RHO(1,J)=RHOINF
           RHOU(1,J)=RHO(1,J)*U(1,J)
           RHOV(1,J)=RHO(1,J)*V(1,J)
           E(1,J)=P(1,J)/(GAMMAN-1.0D0)+0.5D0*RHO(1,J)*(U(1,J)**2+V(1,J)**2)         
        END DO
       
           
        END SUBROUTINE INITIAL


        SUBROUTINE BC
!==========================================================================
! BOUNDARY CONDITION
!==========================================================================
        USE GLOBAL
        IMPLICIT NONE

!____________________LOWER WALL (INVISCID WALL B C)_________________________

           J=1
        DO I=2,NX-1
          
           U(I,J)=U(I,J+1)
           V(I,J)=V(I,J+1)
           P(I,J)=P(I,J+1) 
           RHO(I,J)=RHO(I,J+1)
           RHOU(I,J)=RHO(I,J)*U(I,J)
           RHOV(I,J)=RHO(I,J)*V(I,J)
           E(I,J)= E(I,J+1)  
        END DO
!_____________________UPPER WALL (INVISCID WALL B C)_________________________

           J=NY2
        DO I=2,NX-1
           U(I,J)=U(I,J-1)
           V(I,J)=V(I,J-1)
           P(I,J)=P(I,J-1)
           RHO(I,J)=RHO(I,J-1)
           RHOU(I,J)=RHO(I,J)*U(I,J)
           RHOV(I,J)=RHO(I,J)*V(I,J)
           E(I,J)= E(I,J-1)  
        END DO

       
!________________________OUTLET BOUNDARY CONDITION___________________________

        DO J=1,NY2
           U(NX,J)=U(NX-1,J)
           V(NX,J)=V(NX-1,J)
           P(NX,J)=PINF
!           P(NX,J)=P(NX-1,J)
           RHO(NX,J)=RHO(NX-1,J)
           RHOU(NX,J)=RHOU(NX-1,J)
           RHOV(NX,J)=RHOV(NX-1,J)
           E(NX,J)=P(NX,J)/(GAMMAN-1.0D0)+0.5D0*RHO(NX,J)*(U(NX,J)**2+V(NX,J)**2) 
        END DO 

          
          


       DO I= 1,NX
       DO J= 1,NY2        
           T(I,J)= P(I,J)/(RS*RHO(I,J))
        END DO
        END DO

        END SUBROUTINE BC


        SUBROUTINE GEOMETRY
!==========================================================================
! GEOMETRY AND GRID POINT
!==========================================================================
        USE GLOBAL
        IMPLICIT NONE

           DX=XL/(NX-1)
           DY = DX

           YL1=DY*(NY1-1)
!   
           NY2=NY1
        END SUBROUTINE GEOMETRY


        SUBROUTINE ARRAY_ALLOCATION
!==========================================================================
! ALLOCATE THE ARRAY
!==========================================================================
        USE GLOBAL
        IMPLICIT NONE

           ALLOCATE(U(NX,NY2),V(NX,NY2),P(NX,NY2),E(NX,NY2))
           ALLOCATE(RHO(NX,NY2),RHOU(NX,NY2),RHOV(NX,NY2))
           ALLOCATE(DFI(NX,NY2,4),DFJ(NX,NY2,4),DFF(NX,NY2,4))
           ALLOCATE(DFE(NX,NY2,4), DU(NX,NY2,4))
           ALLOCATE(SOS(NX,NY2),MACH(NX,NY2))
           ALLOCATE(T(NX,NY2))

          

        END SUBROUTINE ARRAY_ALLOCATION



        SUBROUTINE NEW_VARIABLES
!==========================================================================
! COMPUTE REMAING VARIABLES N+1
!==========================================================================
        USE GLOBAL
        IMPLICIT NONE

           DO I= 2,NX-1
           DO J= 2,NY2-1  
              U(I,J)=RHOU(I,J)/RHO(I,J)
              V(I,J)=RHOV(I,J)/RHO(I,J)
              P(I,J)=(E(I,J)-0.5*RHO(I,J)*(U(I,J)**2+V(I,J)**2))*(GAMMAN-1.0D0)  
           END DO
           END DO
        

        END SUBROUTINE NEW_VARIABLES


        SUBROUTINE TIMESTEP
!==========================================================================
! COMPUTE DT USING CFL CONDITION
!==========================================================================
        USE GLOBAL
        IMPLICIT NONE

             DT1=1.0D0
           DO I= 2,NX-1
           DO J= 2,NY2-1 
              SOS(I,J)=DSQRT(GAMMAN*P(I,J)/RHO(I,J))
              DDT = CFL*DX/(DABS(U(I,J))+SOS(I,J))
              DT1 = DMIN1(DDT,DT1)
           END DO
           END DO

 
              DT=DT1
              
        END SUBROUTINE TIMESTEP


        SUBROUTINE MACH_NO
!==========================================================================
! COMPUTE MACH NO. AT EACH NODE
!==========================================================================
        USE GLOBAL
        IMPLICIT NONE

        DO I= 1,NX
        DO J= 1,NY2
           SOS(I,J)=DSQRT(GAMMAN*P(I,J)/RHO(I,J))
           MACH(I,J)=DSQRT(U(I,J)**2+V(I,J)**2)/SOS(I,J)
        END DO
        END DO

        END SUBROUTINE MACH_NO

        SUBROUTINE WRITE_FILE
!==========================================================================
! WRITE VARIABLES IN FILE
!==========================================================================
        USE GLOBAL
        IMPLICIT NONE
        CHARACTER(LEN=20) :: FILE_ID
        CHARACTER(LEN=50) :: FILE_NAME
         

        WRITE(FILE_ID, '(I6)') ITER
        FILE_NAME = 'PRINT' // TRIM(ADJUSTL(FILE_ID)) // '.dat'



        OPEN(FILE=TRIM(FILE_NAME), UNIT = 4)
        WRITE(4,*) ' VARIABLES = "X", "Y", "U", "V", "P", "RHO", "MACH"'
        WRITE(4,*) 'ZONE T="FRAME 0",', 'I=',NX, 'J=',NY2

        DO J=1,NY2
        DO I=1,NX
           
           WRITE(4,*) I,J,U(I,J), V(I,J), P(I,J), RHO(I,J),MACH(I,J)
        END DO
        END DO

        CLOSE(4)

        END SUBROUTINE WRITE_FILE


        SUBROUTINE READ_BY_FILE
!==========================================================================
! READ FROM FILE
!==========================================================================
        USE GLOBAL
        IMPLICIT NONE

         DO J=1,NY2
         DO I=1,NX
            OPEN(2,FILE='READ.dat')
            READ(2,*) II,JJ,U(I,J), V(I,J), P(I,J), RHO(I,J),MACH(I,J)
         END DO
         END DO
         CLOSE(2)

        DO J=1,NY2
        DO I=1,NX
            RHOU(I,J)=RHO(I,J)*U(I,J)
            RHOV(I,J)=RHO(I,J)*V(I,J)
            T(I,J)=P(I,J)/(RS*RHO(I,J))
            E(I,J)=P(I,J)/(GAMMAN-1.0D0)+0.5D0*RHO(I,J)*(U(I,J)**2+V(I,J)**2)
         END DO
         END DO


        END SUBROUTINE READ_BY_FILE






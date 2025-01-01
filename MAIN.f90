!****************************************************************************                  
! *   FILE         = MULTIPLE FILES                                         *
! *   AUTHOR       = CHANDAN KUMAR (CHANDANKR@IITKGP.AC.IN)                 *
! *   INSTITUTE    = INDIAN INSTITUTE OF TECHNOLOGY (IIT), KHARAGPUR        *
!****************************************************************************
! *  THIS PROGRAM IS DISTRIBUTED IN A HOPE THAT IT WILL BE USEFUL.          * 
!****************************************************************************     
	PROGRAM CHANNEL
        USE GLOBAL
        IMPLICIT NONE
           !__________INLET VARIABLES_______________________________
           CINF=SQRT(GAMMAN*RS*TINF)
           MINF=UINF/CINF
           RHOINF=PINF/(RS*TINF)
           PRINT *," PINF:",PINF, "RHOINF:",RHOINF, "TINF:",TINF
           PRINT *," CINF:",CINF, "MINF:",MINF, "UINF:",UINF
           !__________________________STAGNATION CONDITION______________________________
           TTOT= TINF*(1.0D0+0.5D0*(GAMMAN-1.0D0)*MINF**2)
           PTOT= PINF/(TINF/TTOT)**(GAMMAN/(GAMMAN-1.0D0)) 
           RHOTOT= RHOINF*(1.0D0+0.5D0*(GAMMAN-1.0D0)*MINF**2)**(1.0D0/(GAMMAN-1.0D0))  
           !_______________________GEOMETRY AND GRID POINT______________________________
           CALL GEOMETRY
           !_____________________________COMPUTATION OF DT______________________________
           VOL = DX*DY
           DT = CFL*DX/(UINF+CINF)
           PRINT *,"DT_CFL:",DT," CFL:",CFL
           !____________________________ARRAY ALLOCATION________________________________
           CALL ARRAY_ALLOCATION
           CALL INITIAL
           ITER=0
          CALL WRITE_FILE
           CALL BC
!           CALL READ_BY_FILE
        !___________________________MAIN PROCESS_____________________________________
        DO ITER= 1,NTIMES

          CALL FXI
          CALL DFXI
          CALL FYI
          CALL DFYI
          CALL DUT
          CALL NEW_U
          !__________________REMAINING VARIABLE AT N+1______________________
          CALL NEW_VARIABLES
          CALL BC
          !____________________DT USING CFL CONDITION_______________________

          CALL TIMESTEP
              
          !______________________PRINT ITERATION____________________________
          REMAIN= MOD(ITER,50)
          IF (REMAIN.EQ.0) THEN
             PRINT*, ITER, DT
             CALL WRITE_FILE
          END IF

        END DO

        !__________________COMPUTATION OF MACH NO. AT EACH NODE_______________
        CALL MACH_NO

        END


!******************************************************************************|
! diablo.f -> DNS In A Box, Laptop Optimized                       VERSION 2.0
!
! This Fortran 77 code computes incompressible flow in a box.
!
! Primative variables (u,v,w,p) are used, and continuity is enforced with a
! fractional step algorithm.
!
! SPATIAL DERIVATIVES:
!   0, 1, 2, or 3 directions are taken to be periodic and handled spectrally
!   (these cases are referred to as the "periodic", "channel", "duct", and
!    "cavity" cases respectively).
!   The remaining directions are taken to be bounded by walls and handled with
!   momentum- and energy-conserving second-order central finite differences.
!
! TIME ADVANCEMENT
!   Two main approaches are implemented:
!     1. RKW3 on nonlinear terms and CN on viscous terms over each RK substep.
!     2. RKW3 on y-derivative terms and CN on other terms over each RK substep.
!
! The emphasis in this introductory code is on code simplicity:
!   -> All variables are in core.
!   -> The code is not explicitly designed for use with either MPI or SMP.
!   -> Overindexing is not used.
! A few simple high-performance programming constructs are used:
!   -> The inner 2 loops are broken out in such a way as to enable out-of-order
!      execution of these loops as much as possible, thereby leveraging
!      vector and superscalar CPU architectures.
!   -> The outer loops are fairly long (including as many operations as
!      possible inside on a single J plane of data) in order to make effective
!      use of cache.
! Multiple time advancement algorithms are implemented for the periodic,
! channel, duct, and cavity cases in order to compare their efficiency for
! various flows on different computational architectures.  In a few of the
! algorithms, the code is broken out fully into a PASS1/PASS2 architecture
! to maximize the efficient use of cache.
!
!******************************************************************************|
!
! This code is free software; you can redistribute it and/or modify it
! under the terms of the GNU General Public License as published by the
! Free Software Foundation; either version 2 of the License, or (at your
! option) any later version. This code is distributed in the hope that it
! will be useful, but WITHOUT ANY WARRANTY; without even the implied
! warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details. You should have received a
! copy of the GNU General Public License along with this code; if not,
! write to the Free Software Foundation, Inc., 59 Temple Place - Suite
! 330, Boston, MA 02111-1307, USA.
!
!******************************************************************************|

!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      PROGRAM DIABLO
        INCLUDE 'header'
        INTEGER N
        LOGICAL FLAG

        CALL INITIALIZE
! Initialize START_TIME for run timing
        call WALL_TIME(START_TIME)

! A flag to determine if we are considering the first time-step
        FIRST_TIME=.TRUE.
        CALL SAVE_STATS(.FALSE.)
        CALL SAVE_FLOW(.FALSE.)
        DO TIME_STEP = TIME_STEP+1, TIME_STEP+N_TIME_STEPS
          IF (RANK.EQ.0)
     &          WRITE(6,*) 'Now beginning TIME_STEP = ',TIME_STEP

          DO RK_STEP=1,3
            IF (NUM_PER_DIR.EQ.3) THEN
              IF (TIME_AD_METH.EQ.1) CALL RK_PER_1
              IF (TIME_AD_METH.EQ.2) CALL RK_PER_2
            ELSEIF (NUM_PER_DIR.EQ.2) THEN
              IF (TIME_AD_METH.EQ.1) CALL RK_CHAN_1
              IF (TIME_AD_METH.EQ.2) CALL RK_CHAN_2
            ELSEIF (NUM_PER_DIR.EQ.1) THEN
              IF (TIME_AD_METH.EQ.1) CALL RK_DUCT_1
              IF (TIME_AD_METH.EQ.2) CALL RK_DUCT_2
            ELSEIF (NUM_PER_DIR.EQ.0) THEN
              IF (TIME_AD_METH.EQ.1) CALL RK_CAV_1
              IF (TIME_AD_METH.EQ.2) CALL RK_CAV_2
            END IF
          END DO
          TIME=TIME+DELTA_T
          FIRST_TIME=.FALSE.

! Optionally apply a filter to the scalar field
          DO N=1,N_TH
            IF (FILTER_TH(N)
     &         .AND.(MOD(TIME_STEP,FILTER_INT(N)).EQ.0)) THEN
              CALL FILTER_CHAN(N)
            END IF
          END DO

! Save statistics to an output file
          IF (MOD(TIME_STEP,SAVE_STATS_INT).EQ.0) THEN
            CALL SAVE_STATS(.FALSE.)
          END IF
! Save the flow to a restart file
          IF (MOD(TIME_STEP,SAVE_FLOW_INT).EQ.0) THEN
            CALL SAVE_FLOW(.FALSE.)
          END IF

          IF (USE_MPI) THEN
            CALL END_RUN_MPI(FLAG)
          ELSE
            CALL END_RUN(FLAG)
          END IF
          IF ( FLAG ) THEN
            EXIT
          END IF
        END DO

! Calculate and display the runtime for the simulation
        call WALL_TIME(END_TIME)
        IF (RANK.EQ.0) THEN
          WRITE(*,*) 'Elapsed Time (sec): ',end_time-start_time
          WRITE(*,*) 'Seconds per Iteration: '
     &         ,(end_time-start_time)/N_TIME_STEPS
        END IF

        TIME_STEP=TIME_STEP-1
        CALL SAVE_FLOW(.TRUE.)
        CALL SAVE_STATS(.TRUE.)
        IF (RANK.EQ.0) THEN
          WRITE(6,*)
          WRITE(6,*) '    ****** Hello world!  Have a nice day! ******'
          WRITE(6,*)
        END IF
      END



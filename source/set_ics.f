!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE CREATE_FLOW_CHAN
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
        INCLUDE 'header'


        INTEGER I,J,K,N
        REAL*8 RNUM1,RNUM2,RNUM3
        REAL(kind=8) :: X0, X1, Y0, Y1, xcomp, vert_comp, V0
        INTEGER, DIMENSION(:), ALLOCATABLE :: seed

! Initialize the random number generator
        CALL RANDOM_SEED(SIZE = K)
        Allocate (seed(1:K))
        do k=1,K
          seed(k)=RANK+k+999
        end do
        CALL RANDOM_SEED(PUT = seed)

! UBULK0 and KICK should be set in input.dat

! IC_TYPE is set in input_chan.dat and can be used to easily
! control which initial condition is used.  A few examples
! are given here. These can be modified, or new types can be
! added

        IF (IC_TYPE.eq.0) then
! Parabolic profile for laminar closed channel flow
          DO J=0,NY
            DO K=0,NZP-1
              DO I=0,NXM
                U1(I,K,J)=(3./2.)*UBULK0*(1.d0-GYF(J)**2.)
                U2(I,K,J)=0.
                U3(I,K,J)=0.
              END DO
            END DO
          END DO
        else if (IC_TYPE.eq.1) then
! Laminar profile for open channel flow :
          DO K=0,NZP-1
            DO I=0,NXM
              DO J=1,NY
                U1(I,K,J)=-(3./2.)*UBULK0*GYF(J)**2.+3.*UBULK0*GYF(J)
                U2(I,K,J)=0.
                U3(I,K,J)=0.
              END DO
              U1(I,K,0)=0.
              U3(I,K,0)=0.
              U1(I,K,NY+1)=0.
              U3(I,K,NY+1)=0.
            END DO
          END DO
        else if (IC_TYPE.eq.2) then
! Linear profile for laminar Couette flow:
          DO J=0,NY
            DO K=0,NZP-1
              DO I=0,NXM
                U1(I,K,J)=gyf(j)
                U2(I,K,J)=0.
                U3(I,K,J)=0.
              END DO
            END DO
          END DO
        else if (IC_TYPE.eq.3) then
! Tanh shear layer
          DO J=0,NY
            DO K=0,NZP-1
              DO I=0,NXM
                U1(I,K,J)=TANH(GYF(J))
                U2(I,K,J)=0.d0
                U3(I,K,J)=0.d0
              END DO
            END DO
          END DO
        else if (IC_TYPE.eq.4) then
! For Front
! Initialize in thermal wind balance:
          DO J=0,NY
            DO K=0,NZP-1
              DO I=0,NXM
                U1(I,K,J)=0.d0
                U2(I,K,J)=0.d0
                U3(I,K,J)=0.d0
              END DO
            END DO
          END DO


        else if (IC_TYPE.eq.5) then
          X0 = LX/2
          X1 = 1600
          Y0 = 0
          Y1 = 80
          V0 = 0.85d0
          DO J=0,NY
            vert_comp = EXP(-(GYF(J)-Y0)**2/Y1**2)
            DO K=0,NZP-1
              !print*,"CHOR K", K, GZF(K), GZ(K)
              DO I=0,NXM
                xcomp = ERF(-(GX(I) - X0)/X1) *
     &          exp(-(GX(I) - X0)**2/X1**2)
                U1(I,K,J) = 0.d0
                U2(I,K,J) = 0.d0
                U3(I,K,J) = V0 * vert_comp * xcomp
              END DO
            END DO
          END DO


        else if (IC_TYPE.eq.6) then
          X0 = LX/2
          X1 = 1600
          Y0 = 0
          Y1 = 80
          V0 = 0.4d0
          DO J=0,NY
            vert_comp = EXP(-(GYF(J)-Y0)**2/Y1**2)
            DO K=0,NZP-1
              DO I=0,NXM
        xcomp=exp(-(GZ(K)-5000)**2/X1**2) - exp(-(GZ(K)-11000)**2/X1**2)
                U1(I,K,J) = V0 * vert_comp * xcomp
                if(vert_comp>0.8.and.xcomp>.8)print*,"CHOR U",U1(I,K,J)
                U2(I,K,J) = 0.d0
                U3(I,K,J) = 0.d0
              END DO
            END DO
          END DO

        else
          WRITE(*,*) 'WARNING, unsupported IC_TYPE in CREATE_FLOW'
        end if
! Add random noise in physical space
        CALL RANDOM_NUMBER(RNUM1)
        CALL RANDOM_NUMBER(RNUM1)
        CALL RANDOM_NUMBER(RNUM1)

        DO J=0,NY+1
          DO K=0,NZP-1
            DO I=0,NXM
              CALL RANDOM_NUMBER(RNUM1)
              U1(I,K,J)=U1(I,K,J)+KICK*(RNUM1-0.5d0)
              CALL RANDOM_NUMBER(RNUM1)
              U2(I,K,J)=U2(I,K,J)+KICK*(RNUM1-0.5d0)
              CALL RANDOM_NUMBER(RNUM1)
              U3(I,K,J)=U3(I,K,J)+KICK*(RNUM1-0.5d0)
            END DO
          END DO
        END DO

! Zero the ghost cells
        IF (.NOT.USE_MPI) THEN
          DO K=0,NZM
            DO I=0,NXM
              U1(I,K,0)=0.
              U2(I,K,0)=0.
              U3(I,K,0)=0.
              U1(I,K,NY+1)=0.
              U2(I,K,NY+1)=0.
              U3(I,K,NY+1)=0.
            END DO
          END DO
        END IF

! Convert to Fourier space
        CALL FFT_XZ_TO_FOURIER(U1,CU1,0,NY+1)
        CALL FFT_XZ_TO_FOURIER(U2,CU2,0,NY+1)
        CALL FFT_XZ_TO_FOURIER(U3,CU3,0,NY+1)
        CALL FFT_XZ_TO_FOURIER(P,CP,0,NY+1)

        IF (USE_MPI) THEN
          CALL GHOST_CHAN_MPI
        END IF
! Apply Boundary conditions to velocity field
        IF (USE_MPI) THEN
          CALL APPLY_BC_VEL_MPI
        ELSE
          CALL APPLY_BC_VEL_LOWER
          CALL APPLY_BC_VEL_UPPER
        END IF

!! Remove the divergence of the velocity field
!        CALL REM_DIV_CHAN

!        IF (USE_MPI) THEN
!          CALL GHOST_CHAN_MPI
!        END IF

!! Save various statistics to keep track of the initial condition
!        CALL SAVE_STATS_CHAN(.FALSE.)

        RETURN
      END


!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE CREATE_TH_CHAN
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
! Initialize the scalar fields
! In this subroutine, you should initialize each scalar field for the
! particular problem of interest

        INCLUDE 'header'
        INTEGER I,J,K,N
! A variable for Front case...
        REAL*8 RI_B(0:NY+1)
        REAL*8 DRDX
        REAL*8 X0, X1, Y0, Y1, dBdy, B_inf, AUX0, AUX1, vert_comp, V0

        DO N=1,N_TH
          IF (CREATE_NEW_TH(N)) THEN

            IF (IC_TYPE.eq.0) THEN
! As an example, initialize TH1 with a sine in x
              DO K=0,NZP-1
                DO I=0,NXM
                  DO J=1,NY
                    TH(I,K,J,N)=sin(2.d0*PI*GX(I)/LX)/(4.d0*PI**2.d0)
                  END DO
                END DO
              END DO
            ELSE IF ((IC_TYPE.eq.1).or.(IC_TYPE.eq.2)) THEN
! Initialize with a linear profile using the bcs
              DO K=0,NZP-1
                DO I=0,NXM
                IF ((TH_BC_YMIN(N).EQ.0).AND.(TH_BC_YMAX(N).EQ.0)) THEN
                    DO J=1,NY
                      IF (GYF(J).LE.2.0) THEN
                        TH(I,K,J,N)=(TH_BC_YMAX_C1(N)-TH_BC_YMIN_C1(N))
     &                       *(GYF(J)+1.)/2.0+TH_BC_YMIN_C1(N)
                      ELSE
                        TH(I,K,J,N)=TH_BC_YMAX_C1(N)
                      END IF
                    END DO
                  ELSE IF ((TH_BC_YMIN(N).EQ.1)
     &                   .AND.(TH_BC_YMAX(N).EQ.1)) THEN
                    DO J=1,NY
! Linear profile with slope corresponding to lower value
                      TH(I,K,J,N)=TH_BC_YMIN_C1(N)*GYF(J)
                    END DO
                  ELSE
                    IF (RANK.EQ.0) then
                    WRITE(*,*) 'WARNING, THETA INITIALIZED TO ZERO ...'
               WRITE(*,*) 'CREATE AN INITIAL VALUE IN CREATE_FLOW_CHAN'
                    end if
                  END IF
                END DO
              END DO
            ELSE IF (IC_TYPE.eq.3) then
! Tanh profile
              DO K=0,NZP-1
                DO I=0,NXM
                  DO J=1,NY
                    TH(I,K,J,N)=TANH(GYF(J))
                  END DO
                END DO
              END DO
            ELSE IF (IC_TYPE.eq.4) THEN
              IF (N.eq.1) THEN
! For Front case, specify given RI_B profile
                DO K=0,NZP-1
                  DO I=0,NXM
                    TH(I,K,0,N)=0.d0
                    DO J=1,NY
!            RI_B(J)=4.d3*(0.5*tanh((-80.d0-GYF(J))/10.d0)+0.5d0)
!         TH(I,K,J,N)=4.d3*(0.5d0*LOG(COSH((-80.d0-GYF(J))/10.d0))*-10.d0
!     &                      +0.5d0*GYF(J))
!     &                *(RI(N)*0.00000003d0)**2.d0
!     &                    /I_RO**2.d0/RI(N)
                      DRDX=DRHODX(N)
                      if (DRDX.eq.0) then
                        DRDX=3.d-8
                      end if
                      if (GYF(J).lt.-60.d0) then
                        RI_B(J)=20.d0
                        TH(I,K,J,N)=(GYF(J)-GYF(1))*
     &                             RI_B(J)*(RI(N)*DRDX)**2.d0
     &                             /I_RO**2.d0/RI(N)
                      else
                        RI_B(J)=0.0d0
                        TH(I,K,J,N)=(GYF(J)+60.d0)*
     &                         RI_B(J)*(RI(N)*DRDX)**2.d0
     &                         /I_RO**2.d0/RI(N)
     &                        +(-60+140.d0)*20.d0*(RI(N)*DRDX)**2.d0
     &                         /I_RO**2.d0/RI(N)
                      end if
                    END DO
                  END DO
                END DO
              ELSE IF (N.eq.2) then
                DO K=0,NZP-1
                  DO I=0,NXM
                    TH(I,K,0,N)=0.d0
                    DO J=1,NY
                      if (GYF(J).lt.-60.d0) then
                        TH(I,K,J,N)=1
                      else
                        TH(I,K,J,N)=0
                      end if
                    END DO
                  END DO
                END DO
              ELSE
! Passive tracers
                DO K=0,NZP-1
                  DO I=0,NXM
                    DO J=1,NY
                      TH(I,K,J,N)=0*exp(GYF(J)/10.d0)
                    END DO
                  END DO
                END DO
              END IF



            ELSE IF (IC_TYPE.eq.5) THEN
              X0 = LX/2
              X1 = 1600
              Y0 = 0
              Y1 = 80
              V0 = 0.85d0
              DO K=0,NZP-1
                DO I=0,NXM
                  DO J=1,NY
                    vert_comp = EXP(-(GYF(J)-Y0)**2/Y1**2)
                    AUX1 = 2*I_RO*V0*(Y0 - GYF(J)) * vert_comp / Y1**2
                    AUX0 = sqrt(pi)*X1*(1 - ERF(-(GX(I) - X0)/X1)**2)/4
                    dBdy = 2e-4 * 25/1000 ! alpha * 25C / km
                    B_inf = 9.81*(1 + dBdy * GYF(j))
                    TH(I,K,J,N)=B_inf + AUX1 * AUX0
                    !TH(I,K,J,N) = B_inf + AUX1
                  END DO
                END DO
              END DO



            ELSE IF (IC_TYPE.eq.6) THEN
              X0 = LX/2
              X1 = 1600
              Y0 = 0
              Y1 = 80
              V0 = 0.4d0
              DO K=0,NZP-1
                DO I=0,NXM
                  DO J=1,NY
                    vert_comp = EXP(-(GYF(J)-Y0)**2/Y1**2)
                    AUX1 = (sqrt(pi)*X1*erf(-5000/X1 + GZ(K)/X1) -
     & sqrt(pi)*X1*erf(-11000/X1 + GZ(K)/X1))
               AUX0=(Y0-GYF(J))*I_RO*v0*exp(-(Y0-GYF(J))**2/Y1**2)/Y1**2
                    dBdy = 2e-4 * 25/1000 ! alpha * 25C / km
                    B_inf = 9.81*(1 + dBdy * GYF(j))
                    TH(I,K,J,N)=B_inf - AUX1 * AUX0
                    !TH(I,K,J,N) = B_inf + AUX1
                  END DO
                END DO
              END DO



            ELSE
              WRITE(*,*) 'WARNING, unsupported IC_TYPE in CREATE_FLOW'
            END IF


            S1(:,:,:)=TH(:,:,:,N)
            CALL FFT_XZ_TO_FOURIER(S1,CS1,0,NY+1)
            CTH(:,:,:,N)=CS1(:,:,:)

          END IF
        END DO

        RETURN
      END


!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE CREATE_FLOW_PER
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
        INCLUDE 'header'
        INTEGER I, J, K
        REAL*8 RNUM1,RNUM2,RNUM3
        REAL*8 K0

! For an initial vortex, define the location of the centerline
        REAL*8 XC(0:NY+1),ZC(0:NY+1)

        WRITE(6,*) 'Creating new flow from scratch.'

! Initialize random number generator
        CALL RANDOM_SEED

        IF (IC_TYPE.eq.0) THEN
! Initizlize the flow using a Taylor-Green vortex
! Nondimensionalize with U0 and 1/kappa
          DO J=0,NYM
            DO K=0,NZM
              DO I=0,NXM
                U1(I,K,J)=cos(2*pi*(GY(J))/LY)
     &                   *cos(2*pi*(GX(I))/LX)
     &                   *SIN(2*pi*(GZ(K))/LZ)
                U2(I,K,J)=0.d0
                U3(I,K,J)=-cos(2*pi*(GY(J))/LY)
     &                   *sin(2*pi*(GX(I))/LX)
     &                   *COS(2*pi*(GZ(K))/LZ)
              END DO
            END DO
          END DO
        ELSE IF (IC_TYPE.eq.1) THEN
! Start with an ideal vortex centered in the domain
          DO J=0,NYM
            XC(J)=LX/2.
            ZC(J)=LZ/2.
            DO K=0,NZM
              DO I=0,NXM
                IF ((GX(I)-XC(j))**2.+(GZ(K)-ZC(j))**2..gt.0.1) then
! If we aren't too close to the vortex center
                  U1(I,K,J)=-1.d0*(GZ(K)-ZC(j))
     &                    /((GX(I)-XC(j))**2.+(GZ(K)-ZC(j))**2.)
                  U3(I,K,J)=1.d0*(GX(I)-XC(j))
     &                    /((GX(I)-XC(j))**2.+(GZ(K)-ZC(j))**2.)
                  U2(I,K,J)=0.d0
                ELSE
! Otherwise:
                  U1(I,K,J)=-1.d0*(GZ(K)-ZC(j))
     &                    /0.1
                  U3(I,K,J)=1.d0*(GX(I)-XC(j))
     &                    /0.1
                  U2(I,K,J)=0.d0
                END IF
              END DO
            END DO
          END DO
        END IF
! Add random noise in Fourier space

        CALL FFT_XZY_TO_FOURIER(U1,CU1)
        CALL FFT_XZY_TO_FOURIER(U2,CU2)
        CALL FFT_XZY_TO_FOURIER(U3,CU3)
        DO J=1,TNKY
          DO K=1,TNKZ
            DO I=1,NKX
              CALL RANDOM_NUMBER(RNUM1)
              CALL RANDOM_NUMBER(RNUM2)
              CALL RANDOM_NUMBER(RNUM3)
              K0=sqrt(KX(I)**2.d0+KY(J)**2.d0+KZ(K)**2.d0)
     &          /sqrt(KX(1)**2.d0+KY(1)**2.d0+KZ(1)**2.d0)
              CU1(I,K,J)=CU1(I,K,J)+(RNUM1-0.5)*KICK/K0
              CU2(I,K,J)=CU2(I,K,J)+(RNUM1-0.5)*KICK/K0
              CU3(I,K,J)=CU3(I,K,J)+(RNUM1-0.5)*KICK/K0
            end do
          end do
        end do
! get the initial energy in low wavenumbers
        CALL FFT_XZY_TO_PHYSICAL(CU1,U1)
        CALL FFT_XZY_TO_PHYSICAL(CU2,U2)
        CALL FFT_XZY_TO_PHYSICAL(CU3,U3)
        EK0=0.d0
        DO J=0,NYM
          DO K=0,NZM
            DO I=0,NXM
              EK0=EK0+U1(I,K,J)**2.d0+U2(I,K,J)**2.d0+U3(I,K,J)**2.d0
            END DO
          END DO
        END DO
        write(*,*) 'EK0: ',EK0
        IF (N_TH.gt.0) THEN
!      EPSILON_TARGET=((1.d0/DX(1))**4.d0)*(NU**3.d0)*(PR(1))**(-2.d0)
        EPSILON_TARGET=((1.d0/DX(1))**4.d0)*(NU**3.d0)*(100.d0)**(-2.d0)
          write(*,*) 'EPSILON_TARGET: ',EPSILON_TARGET
          write(*,*) 'TARGET KOLMOGOROV SCALE: ',
     &             (NU**3.d0/epsilon_target)**(0.25d0)
        END IF
        CALL FFT_XZY_TO_FOURIER(U1,CU1)
        CALL FFT_XZY_TO_FOURIER(U2,CU2)
        CALL FFT_XZY_TO_FOURIER(U3,CU3)


        CALL REM_DIV_PER
        CALL POISSON_P_PER

        CALL SAVE_STATS_PER(.FALSE.)

        RETURN
      END
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE CREATE_TH_PER
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
        INCLUDE 'header'
        INTEGER I,J,K,N

! Note, Since stratification is not permitted in the periodic flow field
! Any background stratification must be added to the governing equations

        DO N=1,N_TH
          IF (CREATE_NEW_TH(N)) THEN
            DO J=0,NYM
              DO K=0,NZM
                DO I=0,NXM
! Example: Gaussian patch centered in the domain
                  TH(I,K,J,N)=EXP(-((GX(I)-LX/2)*10.d0)**2.d0
     &                            -((GY(J)-LY/2)*10.d0)**2.d0
     &                            -((GZ(K)-LZ/2)*10.d0)**2.d0)
                END DO
              END DO
            END DO
            CALL FFT_XZY_TO_FOURIER(TH(0,0,0,N),CTH(0,0,0,N))

          END IF
        END DO

        RETURN
      END

!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE CREATE_FLOW_DUCT
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|

        RETURN
      END


!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE CREATE_FLOW_CAV
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|

        RETURN
      END






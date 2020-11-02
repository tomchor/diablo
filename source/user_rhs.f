      SUBROUTINE USER_RHS_CHAN_PHYSICAL
        include 'header'
! Here, you can add terms to the right hand side
! of the momentum and scalar equations.
! The right hand side forcing arrays, CF2, CF2, CF3, CFTH
! are in Fourier space.  The velocity and scalars are available
! in physical space.
! S1 is available as a working variable

        integer i,j,k,n

        real*8 alpha

        RETURN
      END


      SUBROUTINE USER_RHS_CHAN_FOURIER
        include 'header'
! Here, you can add terms to the right hand side
! of the momentum and scalar equations.
! The right hand side forcing arrays, CF2, CF2, CF3, CFTH
! are in Fourier space.  The velocity and scalars are available
! in Fourier space.
! S1 is available as a working variable

        integer i,j,k,n

! Advection owing to thermal wind
        IF ((FLAVOR.eq.'Front').and.(I_RO.ne.0.d0)) THEN
          DO N=1,N_TH
! Loop over all scalars

! Add thermal wind advection to the momentum equations
            do j=JSTART,JEND
              do k=0,TNKZ
                do i=0,NXP-1
                  CF1(I,K,J)=CF1(I,K,J)
     &                     -(DRHODX(N)*RI(N)*(GYF(J)+LY/2.d0)/I_RO)
     &                                      *CIKZ(K)*CU1(I,K,J)
     &                  -(-1.d0*DRHODZ(N)*RI(N)*(GYF(J)+LY/2.d0)/I_RO)
     &                                      *CIKX(I)*CU1(I,K,J)
     &                             -(-1.d0*DRHODZ(N)*RI(N)/I_RO)
     &                                *0.5d0*(CU2(I,K,J)+CU2(I,K,J+1))
                  CF3(I,K,J)=CF3(I,K,J)
     &                     -(DRHODX(N)*RI(N)*(GYF(J)+LY/2.d0)/I_RO)
     &                                     *CIKZ(K)*CU3(I,K,J)
     &                 -(-1.d0*DRHODZ(N)*RI(N)*(GYF(J)+LY/2.d0)/I_RO)
     &                                     *CIKX(I)*CU3(I,K,J)
     &                             -(DRHODX(N)*RI(N)/I_RO)
     &                               *0.5d0*(CU2(I,K,J)+CU2(I,K,J+1))
                end do
              end do
            end do

            do j=2,NY
              do k=0,TNKZ
                do i=0,NXP-1
                  CF2(I,K,J)=CF2(I,K,J)
     &                      -(DRHODX(N)*RI(N)*(GY(J)+LY/2.d0)/I_RO)
     &                               *CIKZ(K)*CU2(I,K,J)
     &                  -(-1.d0*DRHODZ(N)*RI(N)*(GY(J)+LY/2.d0)/I_RO)
     &                               *CIKX(I)*CU2(I,K,J)
                end do
              end do
            end do

! Add advection by thermal wind to the scalar equations
            DO J=JSTART_TH(N),JEND_TH(N)
              DO K=0,TNKZ
                DO I=0,NXP-1
                  CFTH(I,K,J,N)=CFTH(I,K,J,N)
     &            -(RI(N)/I_RO)*DRHODX(N)*(GYF(J)+LY/2.d0)
     &                       *CIKZ(K)*CTH(I,K,J,N)
     &            -(RI(N)/I_RO)*-1.d0*DRHODZ(N)*(GYF(J)+LY/2.d0)
     &                       *CIKX(I)*CTH(I,K,J,N)
                END DO
              END DO
            END DO

! End do N_TH
          END DO

        END IF

! Add sponge layer forcing

! For dense tracers, only apply damping to the buoyancy (N=1)
!      N=1
        CALL SPONGE_TH(N)

        CALL SPONGE_VEL

        RETURN
      END


      SUBROUTINE USER_RHS_PER_PHYSICAL
        include 'header'
C Optionally, add forcing terms to the right hand side
C of the momentum and scalar evolution equations
C Here, CF1, CF2, CF3, and CFTH are the forcing arrays in Fourier space
C U1, U2, U3, and TH are the velocity and scalar variables in Physical space

        integer i,j,k,n

C For forced homogeneous turbulence:
C Add some forcing to the system to keep the Batchelor scale fixed
        EK=0.d0
        DO J=0,NYM
          DO K=0,NZM
            DO I=0,NXM
              EK=EK+U1(I,K,J)**2.d0+U2(I,K,J)**2.d0+U3(I,K,J)**2.d0
            END DO
          END DO
        END DO
C Note, that each cell has the same volume, so we can just average over all points
        EK=EK/dble(NX*NY*NZ)

! Scale EK by an amount to compensate for dissipation from 2/3 de-aliasing:
        EK=0.8d0*EK
        DO J=0,NYM
          DO K=0,NZM
            DO I=0,NXM
              S1(I,K,J)=(EPSILON_TARGET/EK)*U1(I,K,J)
            END DO
          END DO
        END DO
        CALL FFT_XZY_TO_FOURIER(S1,CS1)
        DO J=0,TNKY
          DO K=0,TNKZ
            DO I=0,NKX
              CF1(I,K,J)=CF1(I,K,J)+CS1(I,K,J)
            END DO
          END DO
        END DO
        DO J=0,NYM
          DO K=0,NZM
            DO I=0,NXM
              S1(I,K,J)=(EPSILON_TARGET/EK)*U2(I,K,J)
            END DO
          END DO
        END DO
        CALL FFT_XZY_TO_FOURIER(S1,CS1)
        DO J=0,TNKY
          DO K=0,TNKZ
            DO I=0,NKX
              CF2(I,K,J)=CF2(I,K,J)+CS1(I,K,J)
            END DO
          END DO
        END DO
        DO J=0,NYM
          DO K=0,NZM
            DO I=0,NXM
              S1(I,K,J)=(EPSILON_TARGET/EK)*U3(I,K,J)
            END DO
          END DO
        END DO
        CALL FFT_XZY_TO_FOURIER(S1,CS1)
        DO J=0,TNKY
          DO K=0,TNKZ
            DO I=0,NKX
              CF3(I,K,J)=CF3(I,K,J)+CS1(I,K,J)
            END DO
          END DO
        END DO

        RETURN
      END

      SUBROUTINE USER_RHS_DUCT_PHYSICAL
        include 'header'
        RETURN
      END

      SUBROUTINE USER_RHS_CAVITY_PHYSICAL
        include 'header'
        RETURN
      END


      SUBROUTINE SPONGE_TH(N)
! This subroutine applies a sponge relaxation (Rayleigh damping) towards a
! specified background state for the temperature field
! The intention is to allow an open boundary
        include 'header'
        integer i,j,k,n
        real*8 L_sponge,L_bottom
        real*8 SPONGE_AMP

! The following variables will store the background state
        real*8 TH_0(-1:NY+1)

        real*8 RI_B(0:NY+1)

! This variable will hold the forcing rate
        real*8 SPONGE_SIGMA(0:NY+1)

! Set the amplitude of the sponge
        SPONGE_AMP=0.01d0
! Set the top of the sponge layer in physical units
        L_sponge=-500 + LY/6
! Set the bottom of the computational domain in physical units
        L_bottom=-500
        DO J=0,NY+1
! Quadratic damping at lower wall
          if (GYF(J).lt.L_sponge) then
            SPONGE_SIGMA(j)=SPONGE_AMP*((L_sponge-GYF(J))
     &        /(L_sponge-L_bottom))**2.d0
          else
            SPONGE_SIGMA(j)=0.d0
          end if
        END DO

! Set the profile for relaxing the mean TH
        DO J=0,NY+1
          TH_0(J)=TH_BC_YMIN_C1(N)*GYF(J)
        END DO

! For front, specify RI_B profile
        if (n.eq.1) then
          TH_0(0)=0.d0
          do j=1,NY
            if (GYF(j).lt.-60.d0) then
              RI_B(j)=20.d0
              TH_0(j)=(GYF(j)-GYF(1))*
     &                      RI_B(j)*(RI(n)*DRHODX(n))**2.d0
     &                      /I_RO**2.d0/RI(n)
            else
              RI_B(j)=1.0d0
              TH_0(j)=(GYF(j)+60.d0)*
     &                 RI_B(j)*(RI(n)*DRHODX(n))**2.d0
     &                 /I_RO**2.d0/RI(n)
     &                 +(-60+140.d0)*20.d0*(RI(n)*DRHODX(n))**2.d0
     &                 /I_RO**2.d0/RI(N)
            end if
          end do
        else
          do j=0,NY+1
            TH_0(j)=0.d0
          end do
        end if

! Add damping to R-K terms
! Damp the perturbations towards 0
        do k=0,TNKZ
          do i=0,NXP-1
            if ((RANKZ.ne.0).or.(i.ne.0).or.(k.ne.0)) then
              do j=JSTART_TH(N),JEND_TH(N)
                CFTH(i,k,j,n)=CFTH(i,k,j,n)
     &                   -SPONGE_SIGMA(j)*(CTH(i,k,j,n)-0.)
              end do
            end if
          end do
        end do
! Damp the mean gradient towards TH_0
        IF (RANKZ.eq.0) THEN
          do j=JSTART_TH(N),JEND_TH(N)
            CFTH(0,0,j,n)=CFTH(0,0,j,n)-SPONGE_SIGMA(j)
     &            *(CTH(0,0,j,n)-TH_0(J))
          end do
        end if

        return
      end


      SUBROUTINE SPONGE_VEL
! This subroutine applies a sponge relaxation (Rayleigh damping) towards a
! specified background state
! The intention is to allow an open boundary
        include 'header'
        integer i,j,k

        real*8 L_sponge,L_bottom
        real*8 SPONGE_AMP

! The following variables will store the background state
        real*8 U1_0(-1:NY+1), U2_0(0:NY+1), U3_0(-1:NY+1)

! This variable will hold the forcing rate
        real*8 SPONGE_SIGMA(0:NY+1)

! Set the amplitude of the sponge
        SPONGE_AMP=0.01d0
        !SPONGE_AMP=0.0001d0
! Set the top of the sponge layer in physical units
        L_sponge=-500 + LY/6
! Set the bottom of the computational domain in physical units
        L_bottom=-500
        DO J=0,NY+1
! Quadratic damping at lower wall
          if (GYF(J).lt.L_sponge) then
            SPONGE_SIGMA(j)=SPONGE_AMP*((L_sponge-GYF(J))
     &        /(L_sponge-L_bottom))**2.d0
          else
            SPONGE_SIGMA(j)=0.d0
          end if
        END DO

! Set the background state
! Here, set the background to be geostrophic, with a linear temperature profile
        do j=0,NY+1
          U1_0(j)=0.d0
          U3_0(j)=0.d0
        end do
        do j=0,NY+1
          U2_0(j)=0.d0
        end do

! Add damping function to explicit R-K
        do k=0,TNKZ
          do i=0,NKX
            if ((i.ne.0).or.(k.ne.0)) then
              do j=jstart,jend
                CF1(I,K,J)=CF1(I,K,J)-SPONGE_SIGMA(j)*(CU1(i,k,j)-0.d0)
                CF3(I,K,J)=CF3(I,K,J)-SPONGE_SIGMA(j)*(CU3(i,k,j)-0.d0)
              end do
              do j=1,NY
                CF2(I,K,J)=CF2(I,K,J)-
     &         0.5*(SPONGE_SIGMA(j)+SPONGE_SIGMA(j+1))*(CU2(i,k,j)-0.d0)
              end do
            end if
          end do
        end do
! Damp mean flow
        do j=jstart,jend
          CF1(0,0,j)=CF1(0,0,j)-SPONGE_SIGMA(j)*(CU1(0,0,j)-U1_0(j))
          CF3(0,0,j)=CF3(0,0,j)-SPONGE_SIGMA(j)*(CU3(0,0,j)-U3_0(j))
        end do
        do j=1,NY
          CF2(0,0,j)=CF2(0,0,j)-SPONGE_SIGMA(j)*(CU2(0,0,j)-U2_0(j))
        end do

        return
      end




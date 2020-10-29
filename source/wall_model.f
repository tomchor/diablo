      SUBROUTINE WALL_MODEL_LOWER
! This subroutine creates an artificial slip condition for the velocity
! field using a modification of the wall model by Schumann (1975) and
! Grotzbach (1987) to account for a constant angle between the outer
! flow and the surface stress owing to system rotation (for an Ekman layer)
! The velocity field should be in **PHYSICAL** space upon calling
! Returns the wall stress in U_BC_LOWER_NWM(0:NXM,0:NZM), etc.

      include 'header'

      integer i,j,k

! alpha0 is the angle of surface stress with respect to the velocity
! at the first gridpoint (in degrees)
      real*8 alpha_0
      parameter (alpha_0=0.d0)

! Constants for the law of the wall
      real*8 ikappa,B,z0
! Smooth wall:
      parameter (ikappa=2.44d0,B=5.2d0,z0=0.d0)     
! Typical options for a rough wall
!      parameter (ikappa=2.44d0,B=0.d0,z0=2.3422d-6)     

! Location at which to apply log law (in wall units)
      real*8 zplus

      real*8 U1_mean,U3_mean,U_mean
      real*8 tauwall_mean_x,tauwall_mean_z
      real*8 tauwall_x(0:NX+1,0:NZP+1)
      real*8 tauwall_z(0:NX+1,0:NZP+1)
      real*8 tauwall_z_test
      real*8 coef_x1,coef_x2,coef_z1,coef_z2
      integer index_x1,index_x2,index_z1,index_z2
      real*8 shift_amp,shift_x,shift_z

! Variables for Newton-Raphson solver
      real*8 f,fp
      real*8 res,error
      integer iterations

! Here gy(2) is the wall location 
! GYF(2) is the first horizontal velocity point away from the wall
! For smooth wall
      zplus=(gyf(2)-gy(2))/NU
! For rough wall
!      zplus=(gyf(2)-gy(2))/z0

! Get the plane average of the horizontal velocity at zplus     
      U1_mean=SUM(U1(0:NXM,0:NZP-1,2))
      U3_mean=SUM(U3(0:NXM,0:NZP-1,2))

      call mpi_allreduce(mpi_in_place,U1_mean,1,MPI_DOUBLE_PRECISION,
     &     MPI_SUM,MPI_COMM_Z,ierror)
      call mpi_allreduce(mpi_in_place,U3_mean,1,MPI_DOUBLE_PRECISION,
     &     MPI_SUM,MPI_COMM_Z,ierror)

      U1_mean=U1_mean/dble(NX*NZ)
      U3_mean=U3_mean/dble(NX*NZ)

      U_mean=sqrt(U1_mean**2.d0+U3_mean**2.d0)      

! Smooth wall
! Make a first guess at the mean utau 
      utau_mean_lower=U_mean         
      iterations=0
      error=1.d-9
      res=1. 
      do while (res.ge.error.and.iterations.lt.50)
        iterations=iterations+1
        utau_mean_lower=abs(utau_mean_lower)
        f=ikappa*utau_mean_lower*log(zplus*utau_mean_lower)
     &          +B*utau_mean_lower-U_mean
        fp=ikappa*log(zplus*utau_mean_lower)+ikappa+B
        res=abs(f/sqrt(utau_mean_lower))
        utau_mean_lower=utau_mean_lower-f/fp
      end do
      if (res.gt.error) then
        write(*,*) 'WARNING: Failed convergence in wall_model'
      end if
! Schumann model
!      utau_mean_lower=1.0d0
! Rough wall
!      utau_mean_lower=U_mean/(ikappa*log(zplus))

! magnitude of the wall stress
      tauwall_mean=utau_mean_lower**2.d0

      do k=0,NZP-1
        do i=0,NXM
! MKP Model, no shift
           tauwall_x(i,k)=tauwall_mean*U1_mean
     &             /sqrt(U1_mean**2.d0+U3_mean**2.d0)
           tauwall_z(i,k)=tauwall_mean*U3_mean
     &             /sqrt(U1_mean**2.d0+U3_mean**2.d0)

! Set the velocity gradient at the wall based on the local wall stress
          U_BC_LOWER_NWM(i,k)=tauwall_x(i,k)/NU
          W_BC_LOWER_NWM(i,k)=tauwall_z(i,k)/NU

! Constant wall stress
!          U_BC_LOWER_NWM(i,k)=tauwall_mean_x/NU
!          W_BC_LOWER_NWM(i,k)=tauwall_mean_z/NU

        end do
      end do

C Calculate the Fourier transform of the NWM boundary condition
C Store the boundary condition in one plane of CS1 which is sized
C correctly for the FFT calls
        DO K=0,NZP-1
          DO I=0,NXM
            S1(I,K,0)=U_BC_LOWER_NWM(I,K)
          END DO
        END DO
        call FFT_XZ_TO_FOURIER(S1,CS1,0,0)
        DO K=0,TNKZ
          DO I=0,NXP-1
            CU_BC_LOWER_NWM(I,K)=CS1(I,K,0)
          END DO
        END DO
C Store the boundary condition in one plane of CS1 which is sized
C correctly for the FFT calls
        DO K=0,NZP-1
          DO I=0,NXM
            S1(I,K,0)=W_BC_LOWER_NWM(I,K)
          END DO
        END DO
        call FFT_XZ_TO_FOURIER(S1,CS1,0,0)
        DO K=0,TNKZ
          DO I=0,NXP-1
            CW_BC_LOWER_NWM(I,K)=CS1(I,K,0)
          END DO
        END DO


      return
      end 


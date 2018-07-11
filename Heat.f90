PROGRAM Heat

  USE ParametersModule, ONLY: &
       DP, nE, nX, nY, nZ, Pi, &
       X_L, X_R, Y_L, Y_R, Z_L, Z_R, t_end, C
  USE DiscretizationModule, ONLY: &
       ComputeIncrement_Heat

  USE OMP_LIB
  
  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER :: iE, iX, iY, iZ, mpierr, iCycle, threads
  REAL(DP) :: dX, dY, dZ, t, dt, wTime, err
  REAL(DP), ALLOCATABLE, DIMENSION(:,:,:,:) :: U, dU, analytic, diff, Ureal

  CALL MPI_INIT( mpierr )
  PRINT*, "mpierr = ", mpierr

  dX = ( X_R - X_L ) / DBLE( nX )
  dY = ( Y_R - Y_L ) / DBLE( nY )
  dZ = ( Z_R - Z_L ) / DBLE( nZ )
  
  ALLOCATE( U (1:nE,0:nX+1,0:nY+1,0:nZ+1) )
  ALLOCATE( dU(1:nE,1:nX+0,1:nY+0,1:nZ+0) )
  ALLOCATE( analytic(1:nE,1:nX,1:nY,1:nZ) )
  ALLOCATE( Ureal(1:nE,1:nX,1:nY,1:nZ) )
  ALLOCATE( diff(1:nE,1:nX,1:nY,1:nZ) )
  
  wTime = MPI_WTIME( )

 CALL OMP_SET_NUM_THREADS(2)

  !$OMP PARALLEL DO COLLAPSE(4)
  DO iZ = 1, nZ
    DO iY = 1, nY
      DO iX = 1, nX
        DO iE = 1, nE
            
           U(iE,iX,iY,iZ) = SIN( 2.0_DP * Pi * ( X_L + (DBLE(iX)-0.5_DP)*dX ) )
           analytic(iE,iX,iY,iZ) = SIN( 2.0_DP * Pi * ( X_L + (DBLE(iX)-0.5_DP)*dX ) )*exp(-4*pi*pi*t_end)
           threads = OMP_GET_NUM_THREADS()
            
        ENDDO
      ENDDO
    ENDDO
 ENDDO
!$OMP END PARALLEL DO

  wTime = MPI_WTIME( ) - wTime


  PRINT*, "wTime (allocation loop) = ", wTime
 
  t = 0.0_DP
  dt= C*dx*dx
 
  wTime = MPI_WTIME( )

  iCycle = 0

  
  DO WHILE( t < t_end )

    iCycle = iCycle + 1
  
    CALL ComputeIncrement_Heat &
         ( nE, nX, nY, nZ, 1, nE, 1, nX, 1, nY, 1, nZ, dX, dY, dZ, U, dU )

    U(1:nE,1:nX,1:nY,1:nZ) = U(1:nE,1:nX,1:nY,1:nZ) + dt * dU(1:nE,1:nX,1:nY,1:nZ)
    
    t = t + dt

  END DO

  wTime = MPI_WTIME( ) - wTime

  PRINT*, "wTime (computational loop) = ", wTime
  
  Ureal = U(1:nE,1:nX,1:nY,1:nz)
  diff = abs(Ureal - analytic)
  err = maxval(diff)
 
  OPEN(unit=10, file = 'HeatFinal.dat')
  WRITE(10,*) U(1,:,1,1)
  CLOSE(10)
  DEALLOCATE( U, dU, analytic, diff, Ureal )
  
  
  CALL MPI_FINALIZE( mpierr )

  WRITE(*,*) 'Advanced to time',t,'in',iCycle, 'steps'
  WRITE(*,*) 'nE=',nE
  WRITE(*,*) 'nX=',nX
  WRITE(*,*) 'nY=',nY 
  WRITE(*,*) 'nZ=',nZ
  WRITE(*,*) 'dt=',dt
  WRITE(*,*) 'error=',err

END PROGRAM Heat


MODULE DiscretizationModule

  USE ParametersModule, ONLY: &
       DP
  USE OMP_LIB
  
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ComputeIncrement_Heat

CONTAINS


  SUBROUTINE ComputeIncrement_Heat &
       ( nE, nX, nY, nZ, iE_L, iE_R, iX_L, iX_R, iY_L, iY_R, iZ_L, iZ_R, dX, dY, dZ, U, dU )

    INTEGER,  INTENT(in)     :: nE, nX, nY, nZ
    INTEGER,  INTENT(in)     :: iE_L, iX_L, iY_L, iZ_L
    INTEGER,  INTENT(in)     :: iE_R, iX_R, iY_R, iZ_R
    REAL(DP), INTENT(in)    :: dX, dY, dZ
    REAL(DP), INTENT(inout) :: U (1:nE,0:nX+1,0:nY+1,0:nZ+1)
    REAL(DP), INTENT(out)   :: dU(1:nE,1:nX+0,1:nY+0,1:nZ+0)

    INTEGER :: iE, iX, iY, iZ


    CALL ApplyBoundaryConditions( nE, nX, nY, nZ, iE_L, iE_R, iX_L, iX_R, iY_L, iY_R, iZ_L, iZ_R, U )
! $OMP TARGET DATA MAP(to:U) MAP(from:dU)  !GPU only
! $OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(4) !GPU only
!$OMP PARALLEL DO COLLAPSE(4) !CPU only
    DO iZ = iZ_L, iZ_R
      DO iY = iY_L, iY_R
        DO iX = iX_L, iX_R
          DO iE = iE_L, iE_R
              
             dU(iE,iX,iY,iZ) = (1/(dx*dx))*(U(iE,iX-1,iY,iZ) - 2*U(iE,iX,iY,iZ) + U(iE,iX+1,iY,iZ)) &
                + (1/(dy*dy))*(U(iE,iX,iY-1,iZ) - 2*U(iE,iX,iY,iZ) + U(iE,iX,iY+1,iZ)) &
                + (1/(dz*dz))*(U(iE,iX,iY,iZ-1) - 2*U(iE,iX,iY,iZ) + U(iE,iX,iY,iZ+1))
             
          ENDDO
        ENDDO
      ENDDO
    ENDDO
!$OMP END PARALLEL DO
! $OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO
! $OMP END TARGET DATA
    
  END SUBROUTINE ComputeIncrement_Heat


  SUBROUTINE ApplyBoundaryConditions &
    ( nE, nX, nY, nZ, iE_L, iE_R, iX_L, iX_R, iY_L, iY_R, iZ_L, iZ_R, U )

    INTEGER,  INTENT(in)    :: nE, nX, nY, nZ
    INTEGER,  INTENT(in)    :: iE_L, iX_L, iY_L, iZ_L
    INTEGER,  INTENT(in)    :: iE_R, iX_R, iY_R, iZ_R
    REAL(DP), INTENT(inout) :: U(1:nE,0:nX+1,0:nY+1,0:nZ+1)

    ! --- X-Dimension
    
      IF (iX_L == 1) THEN
         U(iE_L:iE_R, 0, iY_L:iY_R, iZ_L:iZ_R) = U(iE_L:iE_R, nX, iY_L:iY_R, iZ_L:iZ_R)
      ENDIF

      IF (iX_R == nX) THEN
         U(iE_L:iE_R, nX+1, iY_L:iY_R, iZ_L:iZ_R) = U(iE_L:iE_R, 1, iY_L:iY_R, iZ_L:iZ_R)
      ENDIF
    
    ! --- Y-Dimesnion

      IF (iY_L == 1) THEN
         U(iE_L:iE_R, iX_L:iX_R, 0, iZ_L:iZ_R) = U(iE_L:iE_R, iX_L:iX_R, nY, iZ_L:iZ_R)
      ENDIF

      IF (iY_R == nY) THEN
         U(iE_L:iE_R, iX_L:iX_R, nY+1, iZ_L:iZ_R) = U(iE_L:iE_R, iX_L:iX_R, 1, iZ_L:iZ_R)
      ENDIF

    ! --- Z-Dimension

      IF (iZ_L == 1) THEN
         U(iE_L:iE_R, iX_L:iX_R, iY_L:iY_R, 0) = U(iE_L:iE_R, iX_L:iX_R, iY_L:iY_R, nZ)
      ENDIF

      IF (iZ_R == nZ) THEN
         U(iE_L:iE_R, iX_L:iX_R, iY_L:iY_R, nZ+1) = U(iE_L:iE_R, iX_L:iX_R, iY_L:iY_R, 1)
      ENDIF
    
    
  END SUBROUTINE ApplyBoundaryConditions

  
END MODULE DiscretizationModule

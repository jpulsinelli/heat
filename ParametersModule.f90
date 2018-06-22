MODULE ParametersModule

  IMPLICIT NONE
  PRIVATE

  INTEGER, PARAMETER, PUBLIC :: DP = KIND( 1.d0 )

  REAL(DP), PARAMETER, PUBLIC :: Pi = ACOS( - 1.0_DP )

  INTEGER, PARAMETER, PUBLIC :: nE = 40
  INTEGER, PARAMETER, PUBLIC :: nX = 64
  INTEGER, PARAMETER, PUBLIC :: nY = 64
  INTEGER, PARAMETER, PUBLIC :: nZ = 64

  REAL(DP), PARAMETER, PUBLIC :: E_L = 0.0_DP
  REAL(DP), PARAMETER, PUBLIC :: E_R = 1.0_DP
  REAL(DP), PARAMETER, PUBLIC :: X_L = 0.0_DP
  REAL(DP), PARAMETER, PUBLIC :: X_R = 1.0_DP
  REAL(DP), PARAMETER, PUBLIC :: Y_L = 0.0_DP
  REAL(DP), PARAMETER, PUBLIC :: Y_R = 1.0_DP
  REAL(DP), PARAMETER, PUBLIC :: Z_L = 0.0_DP
  REAL(DP), PARAMETER, PUBLIC :: Z_R = 1.0_DP
  REAL(DP), PARAMETER, PUBLIC :: t_end = 0.1_DP
  REAL(DP), PARAMETER, PUBLIC :: C = 0.49_DP
  
END MODULE ParametersModule
!> \file GridTransform.F90
!> \author Sayop Kim

MODULE GridTransform_m
   USE Parameters_m, ONLY: wp
   USE io_m, ONLY: iControl
   USE SimulationVars_m, ONLY: nmax, error, ErrMax
   USE GridTransformSetup_m, ONLY: InitializeArrays, CalculateA123, &
                                   CalculatePiPsi, ThomasLoop, &
                                   CopyFrontTOBack, CalculateGridJacobian
   USE GridSetup_m, ONLY: GenerateInteriorPoints
   IMPLICIT NONE

CONTAINS

!-----------------------------------------------------------------------------!
   SUBROUTINE GridTransform()
!-----------------------------------------------------------------------------!
   IMPLICIT NONE
   INTEGER :: n

   CALL InitializeArrays
   IF ( iControl == 1) CALL CalculatePiPsi
   DO n = 1, nmax
      CALL CalculateA123  
      CALL ThomasLoop
      IF (error <= ErrMax) return
   ENDDO
   CALL CopyFrontTOBack
   CALL GenerateInteriorPoints
   CALL CalculateGridJacobian
END SUBROUTINE GridTransform

END MODULE

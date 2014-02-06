!> \file GridTransform.F90
!> \author Sayop Kim

MODULE GridTransform_m
   USE Parameters_m, ONLY: wp
   USE io_m, ONLY: iControl, WriteRMSlog, filenameLength
   USE SimulationVars_m, ONLY: nmax
   USE GridTransformSetup_m, ONLY: InitializeArrays, CalculateA123, &
                                   CalculatePiPsi, ThomasLoop, &
                                   CopyFrontTOBack, CalculateGridJacobian, &
                                   RMSres, RMScrit
   USE GridSetup_m, ONLY: GenerateInteriorPoints
   IMPLICIT NONE
   CHARACTER(LEN=filenameLength) :: RMSlogfile = 'RMSlog.dat'
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
      CALL WriteRMSlog(n,RMSlogfile)
      IF (RMSres <= RMScrit) EXIT
   ENDDO
   CALL CopyFrontTOBack
   CALL GenerateInteriorPoints
   CALL CalculateGridJacobian
   END SUBROUTINE GridTransform

END MODULE

!> \file: io.F90
!> \author: Sayop Kim
!> \brief: Provides routines to read input and write output

MODULE io_m
   USE Parameters_m, ONLY: wp
   USE SimulationVars_m, ONLY: nmax
   USE GridTransformSetup_m, ONLY: RMScrit
   IMPLICIT NONE

   PUBLIC :: filenameLength, Gpnts, FEsize, GeoSize, DCsize, &
             ReadGridInput, WriteTecPlot, WriteRMSlog, width, &
             iControl

   REAL(KIND=wp), DIMENSION(3,2) :: Gpnts    ! Geometry points(start,end)
   REAL(KIND=wp) :: width   ! width: domain width
   INTEGER :: FEsize, GeoSize, DCsize
   INTEGER :: iControl
   INTEGER, PARAMETER :: IOunit = 10, filenameLength = 64
   CHARACTER(LEN=50) :: prjTitle
   
CONTAINS

!-----------------------------------------------------------------------------!
   SUBROUTINE ReadGridInput()
!-----------------------------------------------------------------------------!
! Read input files for transformation 1:
!-----------------------------------------------------------------------------!
   USE SimulationVars_m, ONLY: imax, jmax, kmax,&
                               xblkV, cy1, cy2, cy3, cy4, cy5, cy6
   IMPLICIT NONE
   INTEGER :: ios, i, j
   CHARACTER(LEN=8) :: inputVar
   
   OPEN(IOunit, FILE = 'input.dat', FORM = 'FORMATTED', ACTION = 'READ', &
         STATUS = 'OLD', IOSTAT = ios)
   IF(ios /= 0) THEN
      WRITE(*,'(a)') ""
      WRITE(*,'(a)') "Fatal error: Could not open the input data file."
      RETURN
   ELSE
      WRITE(*,'(a)') ""
      WRITE(*,'(a)') "Reading input file for transformation 1"
   ENDIF

   READ(IOunit,*)
   READ(IOunit,'(a)') prjTitle
   WRITE(*,'(4a)') 'Project Title:', '"',TRIM(prjTitle),'"'
   READ(IOunit,*) inputVar, imax
   WRITE(*,'(a,i6)') inputVar,imax
   READ(IOunit,*) inputVar, jmax
   WRITE(*,'(a,i6)') inputVar,jmax
   READ(IOunit,*) inputVar, kmax
   WRITE(*,'(a,i6)') inputVar,kmax
   READ(IOunit,*)
   READ(IOunit,*) inputVar, xblkV(1,1), xblkV(2,1), xblkV(3,1)
   WRITE(*,'(a,3f6.3)') inputVar, xblkV(1,1), xblkV(2,1), xblkV(3,1)
   READ(IOunit,*) inputVar, xblkV(1,2), xblkV(2,2), xblkV(3,2)
   WRITE(*,'(a,3f6.3)') inputVar, xblkV(1,2), xblkV(2,2), xblkV(3,2)
   READ(IOunit,*) inputVar, xblkV(1,3), xblkV(2,3), xblkV(3,3)
   WRITE(*,'(a,3f6.3)') inputVar, xblkV(1,3), xblkV(2,3), xblkV(3,3)
   READ(IOunit,*) inputVar, xblkV(1,4), xblkV(2,4), xblkV(3,4)
   WRITE(*,'(a,3f6.3)') inputVar, xblkV(1,4), xblkV(2,4), xblkV(3,4)
   ! Set remaining corner points for 3D
   DO i = 1, 4
      DO j = 1, 3, 2
         xblkV(j,i+4) = xblkV(j,i)
      ENDDO
      xblkV(2,i+4) = width
   ENDDO
   ! Read Airfoil corner point data
   READ(IOunit,*) inputVar, Gpnts(1,1), Gpnts(2,1), Gpnts(3,1)
   WRITE(*,'(a,3f6.3)') inputVar, Gpnts(1,1), Gpnts(2,1), Gpnts(3,1)
   READ(IOunit,*) inputVar, Gpnts(1,2), Gpnts(2,2), Gpnts(3,2)
   WRITE(*,'(a,3f6.3)') inputVar, Gpnts(1,2), Gpnts(2,2), Gpnts(3,2)
   ! Read Airfoil grid resolution data (at bottom edge)
   READ(IOunit,*) inputVar, FEsize
   WRITE(*,'(a,i6)') inputVar, FEsize
   READ(IOunit,*) inputVar, GeoSize
   WRITE(*,'(a,i6)') inputVar, GeoSize
   READ(IOunit,*) inputVar, DCsize
   WRITE(*,'(a,i6)') inputVar, DCsize
   READ(IOunit,*) inputVar, width
   WRITE(*,'(a,f6.3)') inputVar,width
   READ(IOunit,*)
   READ(IOunit,*)
   READ(IOunit,*)
   READ(IOunit,*)
   READ(IOunit,*)
   READ(IOunit,*)
   READ(IOunit,*)
   READ(IOunit,*) inputVar, cy1
   WRITE(*,'(a,f6.3)') inputVar, cy1
   READ(IOunit,*) inputVar, cy2
   WRITE(*,'(a,f6.3)') inputVar, cy2
   READ(IOunit,*) inputVar, cy3
   WRITE(*,'(a,f6.3)') inputVar, cy3
   READ(IOunit,*) inputVar, cy4
   WRITE(*,'(a,f6.3)') inputVar, cy4
   READ(IOunit,*) inputVar, cy5
   WRITE(*,'(a,f6.3)') inputVar, cy5
   READ(IOunit,*) inputVar, cy6
   WRITE(*,'(a,f6.3)') inputVar, cy6
   READ(IOunit,*)
   READ(IOunit,*) inputVar, nmax
   WRITE(*,'(a,i6)') inputVar, nmax
   READ(IOunit,*)
   READ(IOunit,*) inputVar, RMScrit
   WRITE(*,'(a,f6.3)') inputVar, RMScrit
   READ(IOunit,*)
   READ(IOunit,*) inputVar, iControl
   WRITE(*,'(a,i6)') inputVar, iControl
   ! Set remaining corner points for 3D
   DO i = 1, 4
      DO j = 1, 3, 2
         xblkV(j,i+4) = xblkV(j,i)
      ENDDO
      xblkV(2,i+4) = width
   ENDDO

   CLOSE(IOunit)
   END SUBROUTINE ReadGridInput

   
!-----------------------------------------------------------------------------!
   SUBROUTINE WriteTecPlot(fileName,varList)
!-----------------------------------------------------------------------------!
! Write Tecplot file
!-----------------------------------------------------------------------------!
   USE SimulationVars_m, ONLY: imax, jmax, kmax,&
                               xp, inverseJacobian
   USE GridTransformSetup_m, ONLY: Pi, Psi
   IMPLICIT NONE
   CHARACTER(LEN=filenameLength), INTENT(IN) :: fileName
   CHARACTER(LEN=*), INTENT(IN) :: varList
   INTEGER :: i, j, k

   OPEN(IOunit, File = fileName, FORM = 'FORMATTED', ACTION = 'WRITE')
   ! writes the two line TECPLOT header
   WRITE(IOunit,'(a)') 'Title="' // TRIM(prjTitle) // '"'
   WRITE(IOunit,'(a)') 'Variables=' // TRIM(varList) 

   WRITE(IOunit,'(a)') ""
   WRITE(IOunit,'(a,i6,a,i6,a,i6,a)') 'Zone I=', imax, ', J=', jmax, ', K=', kmax, ', F=POINT'

   DO k = 1, kmax
      DO j = 1, jmax
         DO i = 1, imax
            WRITE(IOunit,'(6g15.6)') xp(1,i,j,k), xp(2,i,j,k), xp(3,i,j,k), &
                                     inverseJacobian(i,j,k), Pi(i,j,k), Psi(i,j,k)
         ENDDO
      ENDDO
   ENDDO
   CLOSE(IOunit)

   END SUBROUTINE WriteTecPlot

!-----------------------------------------------------------------------------!
   SUBROUTINE WriteRMSlog(nIter,fileName)
!-----------------------------------------------------------------------------!
! Write Tecplot file
!-----------------------------------------------------------------------------!
   USE GridTransformSetup_m, ONLY: RMSres
   IMPLICIT NONE
   CHARACTER(LEN=filenameLength), INTENT(IN) :: fileName
   INTEGER :: nIter

   IF ( nIter == 1 ) THEN
      OPEN(IOunit, File = fileName, FORM = 'FORMATTED', ACTION = 'WRITE')
   ELSE
      OPEN(IOunit, File = fileName, FORM = 'FORMATTED', ACTION = 'WRITE', &
           POSITION = 'APPEND')
   ENDIF
   write(IOunit,'(i6,g15.6)') nIter, RMSres
   CLOSE(IOunit)
   END SUBROUTINE WriteRMSlog
END MODULE io_m

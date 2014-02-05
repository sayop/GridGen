!> \file: GridSetup.F90
!> \author: Sayop Kim

MODULE GridSetup_m
   USE Parameters_m, ONLY: wp
   IMPLICIT NONE

   PUBLIC :: InitializeGrid, GenerateInteriorPoints

CONTAINS
!-----------------------------------------------------------------------------!
   SUBROUTINE InitializeGrid()
!-----------------------------------------------------------------------------!
   USE io_m, ONLY: ReadGridInput
   USE SimulationVars_m, ONLY: imax, jmax, kmax,&
                               xblkV, cy
   IMPLICIT NONE

   ! Create Bottom Edge coordinate values
   CALL ReadGridInput
   CALL InitializeGridArrays
   CALL CreateBottomEdge
   CALL SetEdgePnts
   CALL GridPntsAlgbra
   CALL GenerateInteriorPoints

   END SUBROUTINE

!-----------------------------------------------------------------------------!
   SUBROUTINE InitializeGridArrays()
!-----------------------------------------------------------------------------!
      ! imax: number of grid points in i-drection
      ! jmax: number of grid points in j-direction
      ! kmax: number of grid points in k-direction
      ! xp(3,imax,jmax,kmax): curvilinear coordinates in physical space
      USE SimulationVars_m, ONLY: imax, jmax, kmax, &
                                  xp, inverseJacobian
      IMPLICIT NONE

      WRITE(*,'(a)') ""
      WRITE(*,'(a)') "Initializing data arrays..."
      ALLOCATE(xp(3,imax,jmax,kmax))
      ALLOCATE(inverseJacobian(imax,jmax,kmax))
      xp = 0.0_wp
      inverseJacobian = 0.0_wp
   END SUBROUTINE


!-----------------------------------------------------------------------------!
   SUBROUTINE CreateBottomEdge()
!-----------------------------------------------------------------------------!
   USE io_m, ONLY: width, FEsize, GeoSize, DCsize, &
                   Gpnts
   USE SimulationVars_m, ONLY: imax, jmax, kmax,&
                               xblkV
   USE SimulationVars_m, ONLY: BOTedge
   USE SimulationSetup_m, ONLY: UniformSpacing
   IMPLICIT NONE
   INTEGER :: i

   ALLOCATE(BOTedge(3,imax))
   WRITE(*,*) ""
   WRITE(*,*) "Creating Bottome edge point values with Airfoil geometry"
   DO i = 2, FEsize
      BOTedge(1,i) = UniformSpacing(xblkV(1,1), Gpnts(1,1), i, FEsize)
      !BOTedge(2,i) = UniformSpacing(xblkV(2,1), Gpnts(2,1), i, FEsize)
      BOTedge(3,i) = UniformSpacing(xblkV(3,1), Gpnts(3,1), i, FEsize)
   ENDDO
   DO i = FEsize + 1, FEsize + GeoSize - 1
      BOTedge(1,i) = UniformSpacing(Gpnts(1,1), Gpnts(1,2), i-FEsize+1, GeoSize)
      !BOTedge(2,i) = UniformSpacing(Gpnts(2,1), Gpnts(2,2), i-FEsize+1, GeoSize)
      !BOTedge(3,i) = UniformSpacing(Gpnts(3,1), Gpnts(3,2), i-FEsize+1, GeoSize)
      BOTedge(3,i) = Airfoil(BOTedge(1,i))
   ENDDO
   DO i = FEsize + GeoSize, imax - 1
      BOTedge(1,i) = UniformSpacing(Gpnts(1,2), xblkV(1,2), i-FEsize-GeoSize+2, DCsize)
      !BOTedge(2,i) = UniformSpacing(Gpnts(2,2), xblkV(2,2), i-FEsize-GeoSize+2, DCsize)
      BOTedge(3,i) = UniformSpacing(Gpnts(3,2), xblkV(3,2), i-FEsize-GeoSize+2, DCsize)
   ENDDO

   END SUBROUTINE


!-----------------------------------------------------------------------------!
   FUNCTION Airfoil(xx) RESULT(yx)
!-----------------------------------------------------------------------------!
   IMPLICIT NONE
   REAL(KIND=wp) xint, thick, xx, yx
   xint = 1.008930411365_wp
   thick = 0.15_wp
   yx = 0.2969_wp * sqrt(xint * xx) - 0.126_wp * xint * xx - 0.3516_wp * &
        (xint * xx)**2 + 0.2843_wp * (xint * xx)**3 - 0.1015_wp * (xint * xx)**4
   yx = 5.0_wp * thick * yx

   END FUNCTION Airfoil


!-----------------------------------------------------------------------------!
   SUBROUTINE SetEdgePnts()
!-----------------------------------------------------------------------------!
      USE SimulationVars_m, ONLY: imax, jmax, kmax, &
                                  xp, xblkV, BOTedge
      USE SimulationSetup_m, ONLY: UniformSpacing
      IMPLICIT NONE
      INTEGER :: i

      WRITE(*,'(a)') ""
      WRITE(*,'(a)') "Setting Boundary Conditions..."
      !+++++++++++++++++++++++++++++++++++++++++++++++++++!
      ! Assign coordinates value in xblkV(8,3)            !
      ! Below shows 8 vertices defined in one single block!
      !                                                   !
      !         7--------------8                          !
      !        /|             /|                          !
      !       / |            / |                          !
      !      3--------------4  |    z  y                  !
      !      |  |           |  |    | /                   !
      !      |  5-----------|--6    |/                    !
      !      | /            | /     --- x                 !
      !      |/             |/                            !
      !      1--------------2                             !
      !                                                   !
      !+++++++++++++++++++++++++++++++++++++++++++++++++++!
      ! Vertex (1)
      !xblkV(1,1) = 0.0
      !xblkV(2,1) = 0.0
      !xblkV(3,1) = 0.0
      DO i = 1, 3
         xp(i,1,1,1) = xblkV(i,1)
      ENDDO
      ! Vertex (2)
      !xblkV(1,2) = 0.0
      !xblkV(2,2) = 0.0
      !xblkV(3,2) = 0.0
      DO i = 1, 3
         xp(i,imax,1,1) = xblkV(i,2)
      ENDDO
      ! Vertex (3)
      !xblkV(1,3) = 0.0
      !xblkV(2,3) = 0.0
      !xblkV(3,3) = 0.0
      DO i = 1, 3
         xp(i,1,1,kmax) = xblkV(i,3)
      ENDDO
      ! Vertex (4)
      !xblkV(1,4) = 0.0
      !xblkV(2,4) = 0.0
      !xblkV(3,4) = 0.0
      DO i = 1, 3
         xp(i,imax,1,kmax) = xblkV(i,4)
      ENDDO
      ! Vertex (5)
      !xblkV(1,5) = 0.0
      !xblkV(2,5) = 0.0
      !xblkV(3,5) = 0.0
      DO i = 1, 3
         xp(i,1,jmax,1) = xblkV(i,5)
      ENDDO
      ! Vertex (6)
      !xblkV(1,6) = 0.0
      !xblkV(2,6) = 0.0
      !xblkV(3,6) = 0.0
      DO i = 1, 3
         xp(i,imax,jmax,1) = xblkV(i,6)
      ENDDO
      ! Vertex (7)
      !xblkV(1,7) = 0.0
      !xblkV(2,7) = 0.0
      !xblkV(3,7) = 0.0
      DO i = 1, 3
         xp(i,1,jmax,kmax) = xblkV(i,7)
      ENDDO
      ! Vertex (8)
      !xblkV(1,8) = 0.0
      !xblkV(2,8) = 0.0
      !xblkV(3,8) = 0.0
      DO i = 1, 3
         xp(i,imax,jmax,kmax) = xblkV(i,8)
      ENDDO
      !+++++++++++++++++++++++++++++++++++++++++++++++++++
      ! Set up boundary point coordinates at every edge
      !
      !                 
      !         +--------(8)---------+
      !        /|                   /|
      !     (11)|                (12)|
      !      /  |                 /  |
      !     +---------(4)--------+  (6)
      !     |  (5)               |   |
      !     |   |                |   |     z  y
      !     |   |                |   |     | /
      !    (1)  +-------(7)------|---+     |/
      !     |  /                (2) /      ---x
      !     |(9)                 |(10)    
      !     |/                   |/      
      !     +--------(3)---------+
      !
      !+++++++++++++++++++++++++++++++++++++++++++++++++++
      ! edge (1)
      DO i = 2, kmax - 1
         xp(1,1,1,i) = UniformSpacing(xblkV(1,1), xblkV(1,3), i, kmax)
         xp(2,1,1,i) = UniformSpacing(xblkV(2,1), xblkV(2,3), i, kmax)
         xp(3,1,1,i) = GridStretching(xblkV(3,1), xblkV(3,3), i, kmax)
      ENDDO
      ! edge (2)
      DO i = 2, kmax - 1
         xp(1,imax,1,i) = UniformSpacing(xblkV(1,2), xblkV(1,4), i, kmax)
         xp(2,imax,1,i) = UniformSpacing(xblkV(2,2), xblkV(2,4), i, kmax)
         xp(3,imax,1,i) = GridStretching(xblkV(3,2), xblkV(3,4), i, kmax)
      ENDDO
      ! edige (3)
      DO i = 2, imax - 1
         !xp(1,i,1,1) = UniformSpacing(xblkV(1,1), xblkV(1,2), i, imax)
         xp(2,i,1,1) = UniformSpacing(xblkV(2,1), xblkV(2,2), i, imax)
         !xp(3,i,1,1) = UniformSpacing(xblkV(3,1), xblkV(3,2), i, imax)
         xp(1,i,1,1) = BOTedge(1,i)
         xp(3,i,1,1) = BOTedge(3,i)
      ENDDO
      ! edge (4)
      DO i = 2, imax - 1
         xp(1,i,1,kmax) = UniformSpacing(xblkV(1,3), xblkV(1,4), i, imax)
         xp(2,i,1,kmax) = UniformSpacing(xblkV(2,3), xblkV(2,4), i, imax)
         xp(3,i,1,kmax) = UniformSpacing(xblkV(3,3), xblkV(3,4), i, imax)
      ENDDO
      ! edge (5)
      DO i = 2, kmax - 1
         xp(1,1,jmax,i) = xp(1,1,1,i)
         xp(2,1,jmax,i) = UniformSpacing(xblkV(2,5), xblkV(2,7), i, kmax)
         xp(3,1,jmax,i) = xp(3,1,1,i)
      ENDDO
      ! edge (6)
      DO i = 2, kmax - 1
         xp(1,imax,jmax,i) = xp(1,imax,1,i)
         xp(2,imax,jmax,i) = UniformSpacing(xblkV(2,6), xblkV(2,8), i, kmax)
         xp(3,imax,jmax,i) = xp(3,imax,1,i)
      ENDDO
      ! edge (7)
      DO i = 2, imax - 1
         !xp(1,i,jmax,1) = UniformSpacing(xblkV(1,5), xblkV(1,6), i, imax)
         xp(2,i,jmax,1) = UniformSpacing(xblkV(2,5), xblkV(2,6), i, imax)
         !xp(3,i,jmax,1) = UniformSpacing(xblkV(3,5), xblkV(3,6), i, imax)
         xp(1,i,jmax,1) = BOTedge(1,i)
         xp(3,i,jmax,1) = BOTedge(3,i)
      ENDDO
      ! edge (8) 
      DO i = 2, imax - 1
         xp(1,i,jmax,kmax) = xp(1,i,1,kmax)
         xp(2,i,jmax,kmax) = UniformSpacing(xblkV(2,7), xblkV(2,8), i, imax)
         xp(3,i,jmax,kmax) = xp(3,i,1,kmax)
      ENDDO
      ! edge (9)
      DO i = 2, jmax - 1
         xp(1,1,i,1) = UniformSpacing(xblkV(1,1), xblkV(1,5), i, jmax)
         xp(2,1,i,1) = UniformSpacing(xblkV(2,1), xblkV(2,5), i, jmax)
         xp(3,1,i,1) = UniformSpacing(xblkV(3,1), xblkV(3,5), i, jmax)
      ENDDO
      ! edge (10)
      DO i = 2, jmax - 1
         xp(1,imax,i,1) = UniformSpacing(xblkV(1,2), xblkV(1,6), i, jmax)
         xp(2,imax,i,1) = UniformSpacing(xblkV(2,2), xblkV(2,6), i, jmax)
         xp(3,imax,i,1) = UniformSpacing(xblkV(3,2), xblkV(3,6), i, jmax)
      ENDDO
      ! edge (11)
      DO i = 2, jmax - 1
         xp(1,1,i,kmax) = UniformSpacing(xblkV(1,3), xblkV(1,7), i, jmax)
         xp(2,1,i,kmax) = UniformSpacing(xblkV(2,3), xblkV(2,7), i, jmax)
         xp(3,1,i,kmax) = UniformSpacing(xblkV(3,3), xblkV(3,7), i, jmax)
      ENDDO
      ! edge (12)
      DO i = 2, jmax - 1
         xp(1,imax,i,kmax) = UniformSpacing(xblkV(1,4), xblkV(1,8), i, jmax)
         xp(2,imax,i,kmax) = UniformSpacing(xblkV(2,4), xblkV(2,8), i, jmax)
         xp(3,imax,i,kmax) = UniformSpacing(xblkV(3,4), xblkV(3,8), i, jmax)
      ENDDO
   END SUBROUTINE

!-----------------------------------------------------------------------------!
   SUBROUTINE GridPntsAlgbra()
!-----------------------------------------------------------------------------!
      USE SimulationVars_m, ONLY: imax, jmax, kmax, &
                                  xp, xblkV
      USE SimulationSetup_m, ONLY: UniformSpacing
      IMPLICIT NONE
      INTEGER :: i, j, k

      WRITE(*,'(a)') ""
      WRITE(*,'(a)') "Writing grid points on block surface..."

      !+++++++++++++++++++++++++++++++++++++++++
      ! "front plane"
      !     +------------+
      !     |            |   z(k)
      !     | i-k plane  |    |
      !     |  (j = 1)   |    |
      !     1------------+    ---- x(i)
      !
      !k=kmax @---@---@---@---@
      !       |   |   |   |   |  @: edge points (known)
      !       @---o---o---o---@  o: interior points (unknown)
      !       |   |   |   |   |
      !       @---o---o---o---@
      !       |   |   |   |   |
      !   k=1 @---@---@---@---@
      !      i=1             i=imax
      ! x-coordinate is determined along the i=const lines
      ! y-coordinate is same as y of corner (1)
      ! z-coordinate is determined along the k=const lines
      !+++++++++++++++++++++++++++++++++++++++++
      DO i = 2, imax - 1
         DO k = 2, kmax - 1
            xp(1,i,1,k) = UniformSpacing(xp(1,i,1,1), xp(1,i,1,kmax), k, kmax)
            xp(2,i,1,k) = UniformSpacing(xp(2,i,1,1), xp(2,i,1,kmax), k, kmax)
            xp(3,i,1,k) = GridStretching(xp(3,i,1,1), xp(3,i,1,kmax), k, kmax)
         ENDDO
      ENDDO
      !+++++++++++++++++++++++++++++++++++++++++
      ! "back plane"
      !     +------------+
      !     |            |   z(k)
      !     | i-k plane  |    |
      !     | (j = jmax) |    |
      !     5------------+    ---- x(i)
      ! x-coordinate is determined along the i=const lines
      ! y-coordinate is same as y of corner (5)
      ! z-coordinate is determined along the k=const lines
      !+++++++++++++++++++++++++++++++++++++++++
      DO i = 2, imax - 1
         DO k = 2, kmax - 1
            xp(1,i,jmax,k) = xp(1,i,1,k)
            xp(2,i,jmax,k) = UniformSpacing(xp(2,i,jmax,1), xp(2,i,jmax,kmax), k, kmax)
            xp(3,i,jmax,k) = xp(3,i,1,k)
         ENDDO
      ENDDO
      !+++++++++++++++++++++++++++++++++++++++++
      ! "left plane"
      !                 +
      !                /|
      !               / |   j-k plane (i = 1)
      !              /  |
      !             +   +
      !             |  /  z(k) y(j)
      !             | /     |  /
      !             |/      | /
      !             1       |/ 
      ! x-coordinate is same as x of corner (1)
      ! y-coordinate is determined along the j=const lines
      ! z-coordinate is determined along the k=const lines
      !+++++++++++++++++++++++++++++++++++++++++
      DO j = 2, jmax - 1
         DO k = 2, kmax - 1
            xp(1,1,j,k) = UniformSpacing(xp(1,1,j,1), xp(1,1,j,kmax), k, kmax)
            xp(2,1,j,k) = UniformSpacing(xp(2,1,j,1), xp(2,1,j,kmax), k, kmax)
            xp(3,1,j,k) = GridStretching(xp(3,1,j,1), xp(3,1,j,kmax), k, kmax)
         ENDDO
      ENDDO
      !+++++++++++++++++++++++++++++++++++++++++
      ! "right plane"
      !                 +
      !                /|
      !               / |   j-k plane (i = imax)
      !              /  |
      !             +   +
      !             |  /  z(k) y(j)
      !             | /     |  /
      !             |/      | /
      !             2       |/ 
      ! x-coordinate is same as x of corner (2)
      ! y-coordinate is determined along the j=const lines
      ! z-coordinate is determined along the k=const lines
      !+++++++++++++++++++++++++++++++++++++++++
      DO j = 2, jmax - 1
         DO k = 2, kmax - 1
            xp(1,imax,j,k) = UniformSpacing(xp(1,imax,j,1), xp(1,imax,j,kmax), k, kmax)
            xp(2,imax,j,k) = xp(2,1,j,k)
            xp(3,imax,j,k) = xp(3,1,j,k)
         ENDDO
      ENDDO
      !+++++++++++++++++++++++++++++++++++++++++
      ! "bottom plane"
      !           +-------------+
      !          /             /   y(j)
      !         /  i-j plane  /   /
      !        /  (k = 1)    /   /
      !       1-------------+    ---->x(i)
      ! x-coordinate is determined along the i=const lines
      ! y-coordinate is determined along the j=const lines
      ! z-coordinate is same as z of corner (1)
      !+++++++++++++++++++++++++++++++++++++++++
      DO i = 2, imax - 1
         DO j = 2, jmax - 1
            xp(1,i,j,1) = UniformSpacing(xp(1,i,1,1), xp(1,i,jmax,1), j, jmax)
            xp(2,i,j,1) = UniformSpacing(xp(2,i,1,1), xp(2,i,jmax,1), j, jmax)
            xp(3,i,j,1) = xp(3,i,1,1)
         ENDDO
      ENDDO
      !+++++++++++++++++++++++++++++++++++++++++
      ! "top plane"
      !           +-------------+
      !          /             /   y(j)
      !         /  i-j plane  /   /
      !        /  (k = kmax) /   /
      !       3-------------+    ---->x(i)
      ! x-coordinate is determined along the i=const lines
      ! y-coordinate is determined along the j=const lines
      ! z-coordinate is same as z of corner (3)
      !+++++++++++++++++++++++++++++++++++++++++
      DO i = 2, imax - 1
         DO j = 2, jmax - 1
            xp(1,i,j,kmax) = xp(1,i,1,kmax)
            xp(2,i,j,kmax) = xp(2,i,j,1)
            xp(3,i,j,kmax) = UniformSpacing(xp(3,i,1,kmax), xp(3,i,jmax,kmax), j, jmax)
         ENDDO
      ENDDO
   END SUBROUTINE


!-----------------------------------------------------------------------------!
   SUBROUTINE GenerateInteriorPoints()
!-----------------------------------------------------------------------------!
      USE SimulationVars_m, ONLY: imax, jmax, kmax, &
                                  xp, xblkV
      USE SimulationSetup_m, ONLY: UniformSpacing
      IMPLICIT NONE
      INTEGER :: i, j, k

      WRITE(*,'(a)') ""
      WRITE(*,'(a)') "Writing interior grid points..."
      DO i = 2, imax -1
         DO k = 2, kmax - 1
            DO j = 2, jmax - 1
               xp(1,i,j,k) = UniformSpacing(xp(1,i,1,k), xp(1,i,jmax,k), j, jmax)
               xp(2,i,j,k) = UniformSpacing(xp(2,i,1,k), xp(2,i,jmax,k), j, jmax)
               xp(3,i,j,k) = GridStretching(xp(3,i,1,k), xp(3,i,jmax,k), j, jmax)
            ENDDO
         ENDDO
      ENDDO

   END SUBROUTINE

!-----------------------------------------------------------------------------!
   FUNCTION GridStretching(xmin,xmax,indx,indxMax) RESULT(outcome)
!-----------------------------------------------------------------------------!
      !Distribute interior grid points based on stretching coefficient
      !Interpolateion is made by referring to (i,j,k) indices
      USE SimulationVars_m, ONLY: cy
 
      IMPLICIT NONE
      REAL(KIND=wp), INTENT(IN) :: xmin, xmax
      INTEGER, INTENT(IN) :: indx, indxMax
      REAL(KIND=wp) :: outcome, coef
      coef = log(1.0_wp + (exp(-cy) - 1.0_wp) * REAL(indx - 1) / REAL(indxMax - 1))
      outcome = xmin - coef * (xmax - xmin) / cy
   END FUNCTION GridStretching

END MODULE GridSetup_m

!> \file GridTransform.F90
!> \author Sayop Kim

MODULE GridTransformSetup_m
   USE Parameters_m, ONLY: wp
   USE SimulationVars_m, ONLY: imax, jmax, kmax, &
                               xp, cy2, cy3
   IMPLICIT NONE

   PUBLIC CalculateA123, CalculatePiPsi, ThomasLoop, &
          RMScrit, RMSres

   REAL(KIND=wp), ALLOCATABLE, DIMENSION(:,:,:,:,:) :: InverseGridMetrics
   REAL(KIND=wp), ALLOCATABLE, DIMENSION(:,:,:) :: A1, A2, A3, Pi, Psi
   REAL(KIND=wp) :: RMScrit, RMSres
CONTAINS

!-----------------------------------------------------------------------------!
   SUBROUTINE InitializeArrays()
!-----------------------------------------------------------------------------!
   IMPLICIT NONE

   ALLOCATE(A1(imax,1,kmax))
   ALLOCATE(A2(imax,1,kmax))
   ALLOCATE(A3(imax,1,kmax))
   A1 = 0.0_wp
   A2 = 0.0_wp
   A3 = 0.0_wp

   ALLOCATE(Pi(imax,1,kmax)) 
   ALLOCATE(Psi(imax,1,kmax))
   Pi = 0.0_wp
   Psi = 0.0_wp


   END SUBROUTINE InitializeArrays

!-----------------------------------------------------------------------------!
   SUBROUTINE CalculateGridMetrics()
!-----------------------------------------------------------------------------!
   IMPLICIT NONE
   INTEGER :: i, j, k, l, m
  
   ALLOCATE(InverseGridMetrics(imax,jmax,kmax,3,3))
   ! Inverse Grid Metrics
   !    _                      _
   !   |  x_ksi  x_eta  x_zeta  |
   !   |  y_ksi  y_eta  y_zeta  |
   !   |_ z_ksi  z_eta  z_zeta _|
   !
   
   InverseGridMetrics = 0.0

   DO i = 1, imax
      DO j = 1, jmax
         DO k = 1, kmax
            !x,y,z loops
            DO l = 1, 3
               DO m = 1, 3
               ! l indicates physical space coordinates
               ! m indicates computational domain
               
               ENDDO
            ENDDO
         ENDDO
      ENDDO
   ENDDO      
      
   END SUBROUTINE CalculateGridMetrics

!-----------------------------------------------------------------------------!
   SUBROUTINE CalculateA123()
!-----------------------------------------------------------------------------!
! Evaluate A1, A2, A3 coefficients before looping Thomas method
! A1 = (x_k)^2 + (z_k)^2
! A2 = (x_i)*(x_k) + (z_i)*(z_k)
! A3 = (x_i)^2 + (z_i)^2

   IMPLICIT NONE
   INTEGER :: i, j, k
   ! _i: derivative w.r.t ksi
   ! _k: derivative w.r.t zeta
   REAL(KIND=wp) :: x_i, z_i, x_k, z_k

   ! Evaluate only on j=1 surface (2D i-k front plane)
   j = 1

   !WRITE(*,'(a)') ""
   !WRITE(*,'(a)') "Calculating A1, A2, A3 coefficients for the governing equation..."
   DO i = 2, imax - 1
      DO k = 2, kmax - 1
         x_i = 0.5_wp * (xp(1,i+1,j,k) - xp(1,i-1,j,k))
         z_i = 0.5_wp * (xp(3,i+1,j,k) - xp(3,i-1,j,k))
         x_k = 0.5_wp * (xp(1,i,j,k+1) - xp(1,i,j,k-1))
         z_k = 0.5_wp * (xp(3,i,j,k+1) - xp(3,i,j,k-1))
         A1(i,j,k) = x_k**2 + z_k**2
         A2(i,j,k) = x_i*x_k + z_i*z_k
         A3(i,j,k) = x_i**2 + z_i**2
      ENDDO
   ENDDO 

   END SUBROUTINE CalculateA123 

!-----------------------------------------------------------------------------!
   SUBROUTINE CalculatePiPsi()
!-----------------------------------------------------------------------------!
! Initialize Pi and Psy value before moving into pseudo time loop.
   USE SimulationSetup_m, ONLY: UniformSpacing, GridStretching

   IMPLICIT NONE
   INTEGER :: i, j, k
   ! _i: derivative w.r.t ksi
   ! _k: derivative w.r.t zeta
   REAL(KIND=wp) :: x_i, z_i, x_ii, z_ii, &
                    x_k, z_k, x_kk, z_kk

   ! Evaluate only on j=1 surface (2D i-k front plane)
   j = 1

   WRITE(*,'(a)') ""
   WRITE(*,'(a)') "Calculating Pi and Psi variables for controling elliptic grid..."
   ! Evaluate Psi on the boundaries (i=1, i=imax)
   DO i = 1, imax, imax - 1
      DO k = 2, kmax - 1
         x_k = 0.5_wp * (xp(1,i,j,k+1) - xp(1,i,j,k-1))
         z_k = 0.5_wp * (xp(3,i,j,k+1) - xp(3,i,j,k-1))
         x_kk = xp(1,i,j,k+1) - 2.0_wp * xp(1,i,j,k) + xp(1,i,j,k-1)
         z_kk = xp(3,i,j,k+1) - 2.0_wp * xp(3,i,j,k) + xp(3,i,j,k-1)
         IF(abs(x_k) > abs(z_k)) THEN
            Psi(i,j,k) = -x_kk / x_k
         ELSE
            Psi(i,j,k) = -z_kk / z_k
         ENDIF
      ENDDO
   ENDDO
   ! Evaluate Pi on the boundaries (k=1, k=kmax)
   DO k = 1, kmax, kmax - 1
      DO i = 2, imax - 1
         x_i = 0.5_wp * (xp(1,i+1,j,k) - xp(1,i-1,j,k))
         z_i = 0.5_wp * (xp(3,i+1,j,k) - xp(3,i-1,j,k))
         x_ii = xp(1,i+1,j,k) - 2.0_wp * xp(1,i,j,k) + xp(1,i-1,j,k)
         z_ii = xp(3,i+1,j,k) - 2.0_wp * xp(3,i,j,k) + xp(3,i-1,j,k)
         IF(abs(x_i) > abs(z_i)) THEN
            Pi(i,j,k) = -x_ii / x_i
         ELSE
            Pi(i,j,k) = -z_ii / z_i
         ENDIF
      ENDDO
   ENDDO

   ! Evaluate Pi and Psi at interior points
   DO i = 2, imax - 1
      DO k = 2, kmax - 1
         !Psi(i,j,k) = UniformSpacing(Psi(1,j,k), Psi(imax,j,k), i, imax)
         Psi(i,j,k) = GridStretching(Psi(1,j,k), Psi(imax,j,k), i, imax, cy3)
         !Pi(i,j,k) = UniformSpacing(Pi(i,j,1), Pi(i,j,kmax), k, kmax)
         Pi(i,j,k) = GridStretching(Pi(i,j,1), Pi(i,j,kmax), k, kmax, cy2)
      ENDDO
   ENDDO
   END SUBROUTINE CalculatePiPsi


!-----------------------------------------------------------------------------!
   SUBROUTINE ThomasLoop()
!-----------------------------------------------------------------------------!
! Thomas method for solving tridiagonal matrix system
! This subroutine should be run in a pseudo time loop
   IMPLICIT NONE
   INTEGER :: i, j, k

   REAL(KIND=wp), DIMENSION(imax) :: a, b, c, d
   REAL(KIND=wp) :: x_ik, x_k, z_ik, z_k
   RMSres = 0.0_wp
   j = 1

   DO k = 2, kmax - 1
      ! Calculate governing equation w.r.t x-coordinate
      DO i = 1, imax
         IF( i == 1 .or. i == imax ) THEN
             a(i) = 0.0_wp
             b(i) = 1.0_wp
             c(i) = 0.0_wp
             d(i) = xp(1,i,j,k)
         ELSE
             a(i) = A1(i,j,k) * (1.0_wp - 0.5_wp * Pi(i,j,k))
             b(i) = -2.0_wp * (A1(i,j,k) + A3(i,j,k))
             c(i) = A1(i,j,k) * (1.0_wp + 0.5_wp * Pi(i,j,k))
             x_k = 0.5_wp * (xp(1,i,j,k+1) - xp(1,i,j,k-1))
             x_ik = 0.25_wp * ( xp(1,i+1,j,k+1) - xp(1,i+1,j,k-1) &
                               -xp(1,i-1,j,k+1) + xp(1,i-1,j,k-1) )
             d(i) = 2.0_wp * A2(i,j,k) * x_ik - A3(i,j,k) * ( xp(1,i,j,k+1) + &
                                                              xp(1,i,j,k-1) + &
                                                              Psi(i,j,k) * x_k )
         ENDIF
      ENDDO
      ! Call Thomas method solver
      CALL SY(1, imax, a, b, c, d)
      ! Update values at n+1 pseudo time
      DO i = 1, imax
         RMSres = RMSres + (d(i) - xp(1,i,j,k)) ** 2
         xp(1,i,j,k) = d(i)
      ENDDO

      ! Calculate governing equation w.r.t x-coordinate
      DO i =1, imax
         IF( i == 1 .or. i == imax ) THEN
             a(i) = 0.0_wp
             b(i) = 1.0_wp
             c(i) = 0.0_wp
             d(i) = xp(3,i,j,k)
         ELSE
             a(i) = A1(i,j,k) * (1.0_wp - 0.5_wp * Pi(i,j,k))
             b(i) = -2.0_wp * (A1(i,j,k) + A3(i,j,k))
             c(i) = A1(i,j,k) * (1.0_wp + 0.5_wp * Pi(i,j,k))
             z_k = 0.5_wp * (xp(3,i,j,k+1) - xp(3,i,j,k-1))
             z_ik = 0.25_wp * ( xp(3,i+1,j,k+1) - xp(3,i+1,j,k-1) &
                               -xp(3,i-1,j,k+1) + xp(3,i-1,j,k-1) )
             d(i) = 2.0_wp * A2(i,j,k) * z_ik - A3(i,j,k) * ( xp(3,i,j,k+1) + &
                                                              xp(3,i,j,k-1) + &
                                                              Psi(i,j,k) * z_k )
         ENDIF
      ENDDO
      ! Call Thomas method solver
      CALL SY(1, imax, a, b, c, d)
      ! Update values at n+1 pseudo time
      DO i = 1, imax
         RMSres = RMSres + (d(i) - xp(3,i,j,k)) ** 2
         xp(3,i,j,k) = d(i)
      ENDDO
   ENDDO
   END SUBROUTINE ThomasLoop


!-----------------------------------------------------------------------------!
   SUBROUTINE SY(IL,IU,BB,DD,AA,CC)
!-----------------------------------------------------------------------------!
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: IL, IU
   REAL(KIND=wp), DIMENSION(IL:IU), INTENT(IN) :: AA, BB
   REAL(KIND=wp), DIMENSION(IL:IU), INTENT(INOUT) :: CC, DD

   INTEGER :: LP, I, J
   REAL(KIND=wp) :: R

   LP = IL + 1

   DO I = LP, IU
      R = BB(I) / DD(I-1)
      DD(I) = DD(I) - R*AA(I-1)
      CC(I) = CC(I) - R*CC(I-1)
   ENDDO

   CC(IU) = CC(IU)/DD(IU)
   DO I = LP, IU
      J = IU - I + IL
      CC(J) = (CC(J) - AA(J)*CC(J+1))/DD(J)
   ENDDO
   END SUBROUTINE SY


!-----------------------------------------------------------------------------!
   SUBROUTINE CopyFrontTOBack()
!-----------------------------------------------------------------------------!
   USE SimulationVars_m, ONLY: imax, jmax, kmax, xp
   USE SimulationSetup_m, ONLY: UniformSpacing

   IMPLICIT NONE
   INTEGER :: i, k

   DO i = 2, imax - 1
      DO k = 2, kmax - 1
         xp(1,i,jmax,k) = xp(1,i,1,k)
         xp(2,i,jmax,k) = UniformSpacing(xp(2,i,jmax,1), xp(2,i,jmax,kmax), k, kmax)
         xp(3,i,jmax,k) = xp(3,i,1,k)
      ENDDO
   ENDDO

   END SUBROUTINE CopyFrontTOBack


!-----------------------------------------------------------------------------!
   SUBROUTINE CalculateGridJacobian()
!-----------------------------------------------------------------------------!
   USE SimulationVars_m, ONLY: imax, jmax, kmax, xp, inverseJacobian

   IMPLICIT NONE
   INTEGER :: i, j, k
   ! xgst, ygst, zgst: arbitrary ghost cell points
   REAL(KIND=wp) :: x_i, y_i, z_i, x_j, y_j, z_j, x_k, y_k, z_k, &
                    xgst, ygst, zgst

   DO i = 1, imax 
      DO j = 1, jmax
         DO k = 1, kmax
            ! calculate x_i, y_i, z_i
            IF ( i == 1 ) THEN
               xgst = xp(1,i,j,k) - (xp(1,i+1,j,k) - xp(1,i,j,k))
               ygst = xp(2,i,j,k) - (xp(2,i+1,j,k) - xp(2,i,j,k))
               zgst = xp(3,i,j,k) - (xp(3,i+1,j,k) - xp(3,i,j,k))
               x_i = 0.5_wp * (xp(1,i+1,j,k) - xgst)
               y_i = 0.5_wp * (xp(2,i+1,j,k) - ygst)
               z_i = 0.5_wp * (xp(3,i+1,j,k) - zgst)
            ELSEIF ( i == imax ) THEN
               xgst = xp(1,i,j,k) + (xp(1,i,j,k) - xp(1,i-1,j,k))
               ygst = xp(2,i,j,k) + (xp(2,i,j,k) - xp(2,i-1,j,k))
               zgst = xp(3,i,j,k) + (xp(3,i,j,k) - xp(3,i-1,j,k))
               x_i = 0.5_wp * (xgst - xp(1,i-1,j,k))
               y_i = 0.5_wp * (ygst - xp(2,i-1,j,k))
               z_i = 0.5_wp * (zgst - xp(3,i-1,j,k))
            ELSE
               x_i = 0.5_wp * (xp(1,i+1,j,k) - xp(1,i-1,j,k))
               y_i = 0.5_wp * (xp(2,i+1,j,k) - xp(2,i-1,j,k))
               z_i = 0.5_wp * (xp(3,i+1,j,k) - xp(3,i-1,j,k))
            ENDIF
            ! calculate x_j, y_j, z_j
            IF ( j == 1 ) THEN
               xgst = xp(1,i,j,k) - (xp(1,i,j+1,k) - xp(1,i,j,k))
               ygst = xp(2,i,j,k) - (xp(2,i,j+1,k) - xp(2,i,j,k))
               zgst = xp(3,i,j,k) - (xp(3,i,j+1,k) - xp(3,i,j,k))
               x_j = 0.5_wp * (xp(1,i,j+1,k) - xgst)
               y_j = 0.5_wp * (xp(2,i,j+1,k) - ygst)
               z_j = 0.5_wp * (xp(3,i,j+1,k) - zgst)
            ELSEIF ( j == jmax ) THEN
               xgst = xp(1,i,j,k) + (xp(1,i,j,k) - xp(1,i,j-1,k))
               ygst = xp(2,i,j,k) + (xp(2,i,j,k) - xp(2,i,j-1,k))
               zgst = xp(3,i,j,k) + (xp(3,i,j,k) - xp(3,i,j-1,k))
               x_j = 0.5_wp * (xgst - xp(1,i,j-1,k))
               y_j = 0.5_wp * (ygst - xp(2,i,j-1,k))
               z_j = 0.5_wp * (zgst - xp(3,i,j-1,k))
            ELSE
               x_j = 0.5_wp * (xp(1,i,j+1,k) - xp(1,i,j-1,k))
               y_j = 0.5_wp * (xp(2,i,j+1,k) - xp(2,i,j-1,k))
               z_j = 0.5_wp * (xp(3,i,j+1,k) - xp(3,i,j-1,k))
            ENDIF
            ! calculate x_k, y_k, z_k
            IF ( k == 1 ) THEN
               xgst = xp(1,i,j,k) - (xp(1,i,j,k+1) - xp(1,i,j,k))
               ygst = xp(2,i,j,k) - (xp(2,i,j,k+1) - xp(2,i,j,k))
               zgst = xp(3,i,j,k) - (xp(3,i,j,k+1) - xp(3,i,j,k))
               x_k = 0.5_wp * (xp(1,i,j,k+1) - xgst)
               y_k = 0.5_wp * (xp(2,i,j,k+1) - ygst)
               z_k = 0.5_wp * (xp(3,i,j,k+1) - zgst)
            ELSEIF ( k == kmax ) THEN
               xgst = xp(1,i,j,k) + (xp(1,i,j,k) - xp(1,i,j,k-1))
               ygst = xp(2,i,j,k) + (xp(2,i,j,k) - xp(2,i,j,k-1))
               zgst = xp(3,i,j,k) + (xp(3,i,j,k) - xp(3,i,j,k-1))
               x_k = 0.5_wp * (xgst - xp(1,i,j,k-1))
               y_k = 0.5_wp * (ygst - xp(2,i,j,k-1))
               z_k = 0.5_wp * (zgst - xp(3,i,j,k-1))
            ELSE
               x_k = 0.5_wp * (xp(1,i,j,k+1) - xp(1,i,j,k-1))
               y_k = 0.5_wp * (xp(2,i,j,k+1) - xp(2,i,j,k-1))
               z_k = 0.5_wp * (xp(3,i,j,k+1) - xp(3,i,j,k-1))
            ENDIF
            ! Calculate 1/J: Inverse of grid Jacobian           
            inverseJacobian(i,j,k) = x_i * (y_j * z_k - y_k * z_j) - &
                                     x_j * (y_i * z_k - y_k * z_i) + &
                                     x_k * (y_i * z_j - y_j * z_i)
         ENDDO
      ENDDO
   ENDDO

   END SUBROUTINE CalculateGridJacobian


END MODULE

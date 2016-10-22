!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This module contains various driver subroutines which use LAPACK standard subroutines
! to diagonalize/invert matrices of various type.
!
! Written by Giuseppe Santoro and Angelo Russomanno
! Version 02: With a new subroutine Rinvert written by Tommaso Zanca 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE 	mod_matrices

IMPLICIT	NONE

!.......................... Useful private parameters here .............................

INTEGER,        PRIVATE, PARAMETER 	:: idp = KIND(1.0D0)
REAL(KIND=idp), PRIVATE, PARAMETER  	:: ZERO = 0.0D0, ONE = 1.0D0
REAL(KIND=idp), PRIVATE, PARAMETER  	:: Pi = 3.141592653589793D0

!.......................... Useful common variables ....................................

!.......... VECTORS WHICH WE WANT TO ALLOCATE/DEALLOCATE FREELY ........................

INTEGER,           PRIVATE, ALLOCATABLE 	:: IWORK_RH(:), IWORK_RTH(:)
INTEGER,           PRIVATE, ALLOCATABLE 	:: IFAIL_RH(:), IFAIL_RTH(:)
INTEGER,           PRIVATE, ALLOCATABLE 	:: IPIV_CI(:)
INTEGER,           PRIVATE, ALLOCATABLE         :: IPIV_RI(:)
REAL(KIND=idp),    PRIVATE, ALLOCATABLE  	:: AP_RH(:)
REAL(KIND=idp),    PRIVATE, ALLOCATABLE  	:: RWORK_RH(:), RWORK_RTH(:), RWORK_CH(:)
REAL(KIND=idp),    PRIVATE, ALLOCATABLE         :: RWORK_RI(:)
COMPLEX(KIND=idp), PRIVATE, ALLOCATABLE  	:: CWORK_CH(:), CWORK_CI(:)

CONTAINS

!-------------> Useful general subroutines below  <--------------------------------------
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subroutine to diagonalize an NxN Real Hermitean (symmetric) matrix H 
! using lapack DSPEVX
!
! Finds its M lowest eigenvalues   (OUTPUT in E)
! and  (if JOBZ='V') eigenvectors  (OUTPUT in W, by columns)
!
! INPUT: JOBZ = 'V' (Eigenvectors wanted) or JOBZ = 'N' (Only eigenvalues)
!        N (integer) dimension of matrix H(N,N)
!        M (integer) number of lowest eigenvalues wanted
!        H(N,N) (real*8) the symmetric matrix to be diagonalized
! OUTPUT: E(N) (real*8) the eigenvalues
!         W(N,M) (real*8) the M eigenvectors (by columns) if JOBVZ = 'V'
!
! If UPLO = 'L', AP_RH(i + (j-1)*(2*n-j)/2) = H(i,j) for j<=i<=n.
!
! The needed workspace is allocated/deallocated using RHdiagonalize_ALLOC/DEALLOC 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
SUBROUTINE RHdiagonalize( JOBZ, N, M, H, E, W )

IMPLICIT   NONE

INTEGER,          PARAMETER   	:: IL = 1   	! GS always included
CHARACTER(LEN=1), PARAMETER 	:: UPLO  = 'L' 	! Lower packed mode
CHARACTER(LEN=1), PARAMETER 	:: RANGE = 'I' 	! Only selected (lowest) 
                                               	! eigenpairs from IL to IU 
                                       	   	! are calculated
CHARACTER, INTENT(IN)       	:: JOBZ*1
INTEGER,   INTENT(IN)          	:: N, M
REAL(KIND=idp), INTENT(IN)     	:: H(N,N)
REAL(KIND=idp), INTENT(OUT)    	:: E(N), W(N,M)

INTEGER            	  	:: M_OUT, INFO, IU, i, j
REAL(KIND=idp)                  :: ABSTOL, DLAMCH
REAL(KIND=idp)         	  	:: VL, VU  	! Not referenced for RANGE=I

ABSTOL = 2.D0*DLAMCH('S') 	! Asier Laguna said it seems more accurate then ABSTOL=0.D0
IU = M         			! Maximum eigenvalue calculated

DO i=1,N			! Constructs vector AP_RH, the matrix in Packed mode
   DO j=1,i
      AP_RH(i + (j-1)*(2*N-j)/2) = H(i,j)
   ENDDO
ENDDO
 
CALL DSPEVX( JOBZ, RANGE, UPLO, N, AP_RH, VL, VU, IL, IU, &
             ABSTOL, M_OUT, E, W, N, RWORK_RH, IWORK_RH, IFAIL_RH, INFO )

IF (INFO.ne.0) THEN
    WRITE(6,*) ' WARNING: DSPEVX INFO = ', INFO
    DO i=1,M
       WRITE(6,*) '          i, IFAIL_RH(i) = ', i, IFAIL_RH(i)
    ENDDO
ENDIF

RETURN
END SUBROUTINE RHdiagonalize

!.................................................................................
! To Allocate Workspace for RHdiagonalize
!.................................................................................
 
SUBROUTINE RHdiagonalize_ALLOC( N )

IMPLICIT NONE
INTEGER     	:: N

ALLOCATE( IFAIL_RH(N)      )
ALLOCATE( IWORK_RH(5*N)    )
ALLOCATE( RWORK_RH(8*N)    )
ALLOCATE( AP_RH(N*(N+1)/2) )

RETURN
END SUBROUTINE RHdiagonalize_ALLOC

!.................................................................................
! To Deallocate Workspace for RHdiagonalize
!.................................................................................

SUBROUTINE RHdiagonalize_DEALLOC()

IMPLICIT NONE

DEALLOCATE( IFAIL_RH )
DEALLOCATE( IWORK_RH )
DEALLOCATE( RWORK_RH )
DEALLOCATE( AP_RH    )

RETURN
END SUBROUTINE RHdiagonalize_DEALLOC

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subroutine to diagonalize an NxN Real Tridiagonal Hermitean (symmetric) matrix
! using lapack DSTEVX
!
! NB: The case of a COMPLEX Tridiagonal Hermitean matrix can be easily reduced
!     to the REAL case by a diagonal unitary transformation: simply rescale-out
!     the phase-factors appearing in the subdiagonal and then re-insert them
!     when constructing the final eigenvectors
!
! Finds its M lowest eigenvalues (which are output in E)
! and optionally (JOBZ='V') eigenvectors (which are output in W by columns)
!
! DIAG(1..N)    contains the diagonal
! SDIAG(1..N-1) contains the subdiagonal
!
! JOBZ = 'V'  IF eigenvectors wanted, 'N' IF NOT
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE RTHdiagonalize( JOBZ, N, M, DIAG, SDIAG, E, W )

IMPLICIT   NONE

INTEGER,          PARAMETER   	:: IL = 1   		! GS always included
REAL(KIND=idp)                	:: ABSTOL, DLAMCH 
CHARACTER(LEN=1), PARAMETER 	:: RANGE = 'I' 		! Only selected (lowest) 
                                         	 	! eigenpairs from IL to IU 
                                          	   	! are calculated
CHARACTER, INTENT(IN)       	:: JOBZ*1
INTEGER,   INTENT(IN)          	:: N, M
REAL(KIND=idp), INTENT(IN)  	:: DIAG(N), SDIAG(N)
REAL(KIND=idp), INTENT(OUT)  	:: E(N), W(N,M)

INTEGER          	   	:: M_OUT, INFO, IU, i
REAL(KIND=idp)         	  	:: VL, VU  		! Not referenced for RANGE=I

ABSTOL = 2.D0*DLAMCH('S') 	! Asier Laguna said it seems more accurate then ABSTOL=0.D0
IU = M         			! Maximum eigenvalue calculated
    
CALL DSTEVX( JOBZ, RANGE, N, DIAG, SDIAG, VL, VU, IL, IU, & 
 	     ABSTOL, M_OUT, E, W, N, RWORK_RTH, IWORK_RTH, IFAIL_RTH, INFO )

IF (INFO.ne.0) THEN
   WRITE(6,*) ' WARNING: DSTEVX INFO = ', INFO
   DO i=1,M
      WRITE(6,*) '          i, IFAIL_RTH(i) = ', i, IFAIL_RTH(i)
   ENDDO
ENDIF

RETURN
END SUBROUTINE RTHdiagonalize

!.................................................................................
! To Allocate Workspace for RTH_Diagonalize
!.................................................................................
 
SUBROUTINE RTHdiagonalize_ALLOC( N )

IMPLICIT NONE
INTEGER     	:: N

ALLOCATE( IFAIL_RTH(N)   )
ALLOCATE( IWORK_RTH(5*N) )
ALLOCATE( RWORK_RTH(5*N) )

RETURN
END SUBROUTINE RTHdiagonalize_ALLOC

!.................................................................................
! To Deallocate Workspace for RTH_Diagonalize
!.................................................................................

SUBROUTINE RTHdiagonalize_DEALLOC()

IMPLICIT NONE

DEALLOCATE( IFAIL_RTH )
DEALLOCATE( IWORK_RTH )
DEALLOCATE( RWORK_RTH )

RETURN
END SUBROUTINE RTHdiagonalize_DEALLOC
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subroutine to diagonalize an NxN Complex Hermitean matrix
! using lapack ZHEEV
! Computes ALL eigenvalues and (optionally, if JOBZ='V') eigenvectors of H
! INPUT: JOBZ = 'V' (Eigenvectors wanted) or JOBZ = 'N' (Only eigenvalues)
!
! We directly set UPLO='L' for Lower triangle of matrix H. 
! Warning: On exit H is DESTROYED and substituted (if JOBZ='V') by eigenvectors (by colums)
!
! Workspace Allocated/Deallocated by CHdiagonalize_ALLOC/DEALLOC
! RWORK_CH(real*8) dimension should be at least 3*N-2
! CWORK_CH(complex*16) dimension should be LWORK >= 2*N-1
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
SUBROUTINE CHdiagonalize(JOBZ, N, H, E)
 
IMPLICIT NONE

 CHARACTER*1,  INTENT(IN)	:: JOBZ
INTEGER,      INTENT(IN) 	:: N
COMPLEX(idp), INTENT(INOUT) 	:: H(N,N) 
REAL(idp),    INTENT(OUT) 	:: E(N) 

INTEGER          		:: INFO, LWORK
!..................................................................
! To get optimal LWORK call with LWORK=-1 and CWORK_CH(1) will tell 
 LWORK = -1
 CALL   ZHEEV( JOBZ, 'L', N, H, N, E, CWORK_CH, LWORK, RWORK_CH, INFO )
!..................................................................
LWORK = 2*N-1 			! Optimal(see above) INT(CWORK_CH(1)) 
 CALL ZHEEV( JOBZ, 'L', N, H, N, E, CWORK_CH, LWORK, RWORK_CH, INFO )
IF (INFO.ne.0) WRITE(6,*) ' WARNING: ZHHEV INFO = ', INFO
 
RETURN
END SUBROUTINE CHdiagonalize

!.................................................................................
! To Allocate Workspace for CHdiagonalize
!.................................................................................

SUBROUTINE CHdiagonalize_ALLOC( N )

IMPLICIT NONE
INTEGER     	:: N

ALLOCATE( RWORK_CH(3*N-2) )
ALLOCATE( CWORK_CH(2*N-1) )

RETURN
END SUBROUTINE CHdiagonalize_ALLOC

!.................................................................................
! To Deallocate Workspace for CHdiagonalize
!.................................................................................

SUBROUTINE CHdiagonalize_DEALLOC()

IMPLICIT NONE

DEALLOCATE( RWORK_CH )
DEALLOCATE( CWORK_CH )

RETURN
END SUBROUTINE CHdiagonalize_DEALLOC

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Translates a complex Hermitean matrix H into a Real symmetric one HR = |H|
! which is equivalent if you only want the eigenvalues 
! Can be made more efficient by directly calling DSPEVX
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
SUBROUTINE C2Rdiagonalize(N, H, E)
 
IMPLICIT NONE
INTEGER,      INTENT(IN)	:: N
COMPLEX(idp), INTENT(IN)  	:: H(N,N)
REAL(idp),    INTENT(OUT) 	:: E(N)
REAL(idp) 			:: HR(N,N), W(N,N)

HR(:,:) = CDABS( H(:,:) )

CALL RHdiagonalize_ALLOC( N )
CALL RHdiagonalize('N', N, N, HR, E, W )
CALL RHdiagonalize_DEALLOC()
 
RETURN
END SUBROUTINE C2Rdiagonalize
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subroutine to calculate the inverse of a complex matrix using Lapack
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
SUBROUTINE Cinvert(N, H)
 
IMPLICIT NONE
INTEGER,      INTENT(IN) 	:: N
COMPLEX(idp), INTENT(INOUT) 	:: H(N,N)

INTEGER          	:: INFO, LWORK
 
LWORK = N 		! On exit (INFO=0) optimal LWORK=INT( CWORK_CI(1) )
CALL ZGETRF( N, N, H, N, IPIV_CI, INFO )
CALL ZGETRI( N, H, N, IPIV_CI, CWORK_CI, LWORK, INFO )
 
RETURN
END SUBROUTINE Cinvert

!.................................................................................
! To Allocate Workspace for Cinvert
!.................................................................................

SUBROUTINE Cinvert_ALLOC( N )

IMPLICIT NONE
INTEGER     	:: N

ALLOCATE( IPIV_CI(N)  )
ALLOCATE( CWORK_CI(N) )

RETURN
END SUBROUTINE Cinvert_ALLOC

!.................................................................................
! To Deallocate Workspace for Cinvert
!.................................................................................

SUBROUTINE Cinvert_DEALLOC()

IMPLICIT NONE

DEALLOCATE( IPIV_CI  )
DEALLOCATE( CWORK_CI )

RETURN
END SUBROUTINE Cinvert_DEALLOC

! INIZIO TOMMASO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subroutine to calculate the inverse of a real matrix using Lapack
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
SUBROUTINE Rinvert(N, H)
 
IMPLICIT NONE
INTEGER,      INTENT(IN) 	:: N
REAL(KIND=idp), INTENT(INOUT) 	:: H(N,N)

INTEGER          	:: INFO, LWORK
 
LWORK = N 		! On exit (INFO=0) optimal LWORK=INT( CWORK_CI(1) )
CALL DGETRF( N, N, H, N, IPIV_RI, INFO )
CALL DGETRI( N, H, N, IPIV_RI, RWORK_RI, LWORK, INFO )
 
RETURN
END SUBROUTINE Rinvert

!.................................................................................
! To Allocate Workspace for Rinvert
!.................................................................................

SUBROUTINE Rinvert_ALLOC( N )

IMPLICIT NONE
INTEGER     	:: N

ALLOCATE( IPIV_RI(N)  )
ALLOCATE( RWORK_RI(N) )

RETURN
END SUBROUTINE Rinvert_ALLOC

!.................................................................................
! To Deallocate Workspace for Rinvert
!.................................................................................

SUBROUTINE Rinvert_DEALLOC()

IMPLICIT NONE

DEALLOCATE( IPIV_RI  )
DEALLOCATE( RWORK_RI )

RETURN
END SUBROUTINE Rinvert_DEALLOC
! FINE TOMMASO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subroutine to calculate the inverse of a complex matrix 
! Slower than the Lapack one above!
! Inputs:
!   A       Matrix A to be inverted
!   N       Elements used in matrix A (N by N)
!   MAXN    Matrix dimensions as A(MAXN,MAXN)
! Outputs:
!  Ainv     Inverse of matrix A
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
SUBROUTINE Cinvert2( A, N, MAXN, Ainv )

IMPLICIT NONE
integer*4, PARAMETER	:: MAXMAXN = 2000
integer*4 		:: N, MAXN
complex*16 		:: A(MAXN,MAXN), Ainv(MAXN,MAXN)
integer*4 		:: i, j, k, index(MAXMAXN), jPivot, indexJ
real*8 			:: scale(MAXMAXN), scaleMax, ratio, ratioMax
complex*16 		:: AA(MAXMAXN,MAXMAXN), B(MAXMAXN,MAXMAXN), coeff, sum

if( MAXN .gt. MAXMAXN ) then
        write(*,*) 'ERROR in cinv: Matrix too large'
        stop
endif

      !* Matrix B is initialized to the identity matrix
      do i=1,N
       do j=1,N
         AA(i,j) = A(i,j)  ! Copy matrix so as not to overwrite
         B(i,j) = 0.0
       enddo
       B(i,i) = 1.0
      enddo

      !* Set scale factor, scale(i) = max( |a(i,j)| ), for each row
      do i=1,N
        index(i) = i     ! Initialize row index list
        scaleMax = 0.0
        do j=1,N
          if( abs(AA(i,j)) .gt. scaleMax ) then
            scaleMax = abs(AA(i,j))
          endif
        enddo
        scale(i) = scaleMax
      enddo

      !* Loop over rows k = 1, ..., (N-1)
      do k=1,(N-1)
        !* Select pivot row from max( |a(j,k)/s(j)| )
        ratiomax = 0.0
        jPivot = k
        do i=k,N
          ratio = abs(AA(index(i),k))/scale(index(i))
          if( ratio .gt. ratiomax ) then
            jPivot=i
            ratiomax = ratio
          endif
        enddo
        !* Perform pivoting using row index list
        indexJ = index(k)
        if( jPivot .ne. k ) then     ! Pivot
          indexJ = index(jPivot)
          index(jPivot) = index(k)   ! Swap index jPivot and k
          index(k) = indexJ
        endif
        !* Perform forward elimination
        do i=k+1,N
          coeff = AA(index(i),k)/AA(indexJ,k)
          do j=k+1,N
            AA(index(i),j) = AA(index(i),j) - coeff*AA(indexJ,j)
          enddo
          AA(index(i),k) = coeff
          do j=1,N
            B(index(i),j) = B(index(i),j) - AA(index(i),k)*B(indexJ,j)
          enddo
        enddo
      enddo

      !* Perform backsubstitution
      do k=1,N
        Ainv(N,k) = B(index(N),k)/AA(index(N),N)
        do i=N-1,1,-1
          sum = B(index(i),k)
          do j=i+1,N
            sum = sum - AA(index(i),j)*Ainv(j,k)
          enddo
          Ainv(i,k) = sum/AA(index(i),i)
        enddo
      enddo

RETURN
END SUBROUTINE Cinvert2
 
!.................. END OF SUBROUTINES INCLUDED IN THE MODULE ...............................
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE mod_matrices

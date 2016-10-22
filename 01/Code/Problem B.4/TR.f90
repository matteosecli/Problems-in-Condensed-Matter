! the functions f and g

complex(8) function f(kx,ky) result(res)
real(8),intent(inout)::kx,ky
real(8)::t,a
t=1.0
a=1.
res=complex(2*cos(kx*a*0.5)*cos(sqrt(3.0)*ky*a/6.0),2*cos(kx*a*0.5)*sin(sqrt(3.0)*ky*a/6.0))
end function f


complex(8) function g(kx,ky) result(res)
real(8),intent(inout)::kx,ky
real(8)::t,a
t=1.0
a=1.
res=complex(cos(sqrt(3.0)*ky*a/3.0),sin(sqrt(3.0)*ky*a/3.0))
end function g

! ************************************************************
!*************************************************************
  PROGRAM  ARMCHAIR
USE mod_matrices
  IMPLICIT NONE
!********************************************************************
	  INTEGER  ::i,N  ! Matrix Size
	  complex(8),ALLOCATABLE, DIMENSION(:,:)::H,HMAT  !our matrix
          REAL(8)::Eatom,a,kx,ky,t,kxmin,kxmax,hsteps,pi
	  complex(8),external::f,g
	  REAL(8),ALLOCATABLE, DIMENSION(:,:) :: ORTHO
	  INTEGER ::m,j,k,Nsteps,q
          integer::Ncells
	  CHARACTER(LEN=*),PARAMETER :: formatstring='(30(F10.5))' !string to format output
!*****************************************************************
! ****************************************************************
! ************       variables for LAPACK   **********************
          INTEGER :: INFO       !variable to get info from lapack 
	  INTEGER :: LWORK                   !will be the size of workspace
	  INTEGER:: LWMAX
          REAL,allocatable,dimension(:)::RWORK_CH,RWORK  !maximum size of workspace
!	  REAL(8), DIMENSION(LWMAX) :: WORK  !workspace of lapack subroutine
	  CHARACTER(LEN=1) :: JOBZ, UPLO       !characters to interact with lapack

         real(8),ALLOCATABLE,dimension(:)::E!eigenvalues

! ****************************************************************
pi=acos(-1.0)

!The BZ

kxmin=-pi
kxmax=pi   

!Parameters
a=1.0
Eatom=0.0
t=1.0
	 

!The number of cells we call N


	 Ncells=50
	 
	 
	  N=4*Ncells
    ALLOCATE(HMAT(N,N),H(N,N),E(N))
	t=1.0
         open(8,file='eigenval.out',status='unknown',position='append')



! The BZ discretization
Nsteps=600

hsteps=(kxmax-kxmin)/real(Nsteps)

!do q=1,Ncells+1
!m=-Ncells+(q-1)
m=0
ky=2*acos(-1.0)*m/sqrt(3.0)/real(Ncells)

	  do k=1,Nsteps

 kx=kxmin+(k-1)*hsteps




        !****************************************************
        !  THE HAMILTONIAN MATRIX  
        !****************************************************
        !


DO I=1,N
    DO J=1,N

     
       if(I==J) then
      Hmat(I,J)=Eatom

      elseif(I==N .and.j==1) then
      Hmat(I,J)=-t*g(kx,ky)
       elseif(I==1 .and.j==N) then
      Hmat(I,J)=-t*CONJG(g(kx,ky))
     else if (j==i+1.and.mod(i,2).ne.0) then
     Hmat(I,J)=-t*f(kx,ky)
     elseif (j==i+1.and.mod(i,2)==0) then
     Hmat(I,J)=-t*g(kx,ky)
      elseif( i==j+1.and.mod(i,2).ne.0) then
     Hmat(I,J)=-t*CONJG(g(kx,ky))
     elseif (i==j+1.and.mod(i,2)==0) then
      Hmat(I,J)=-t*CONJG(f(kx,ky))
     else
     Hmat(I,J)=0.0
      end if
     END DO
  END DO


!*************************************************
!DEFINE THE LOWER TRIANGLE MATRIX
!*************************************************

DO I=1,N

    DO J=1,N
     if(I>=J) then
      H(I,J)=Hmat(I,J)
      else
     H(I,J)=0.0
      end if
     END DO
  END DO


!!!!

JOBZ='N'
 Call CHdiagonalize_ALLOC( N )
 Call CHdiagonalize(JOBZ, N, H, E)
 Call CHdiagonalize_DEALLOC()


   Do i=1,N
   WRITE(8,*)kx,E(I)
   End do
    
end do
!end do

!DEALLOCATE(H,E) 

END PROGRAM









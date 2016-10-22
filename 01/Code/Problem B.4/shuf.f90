program shuff
implicit none


integer::NSTEPS,u,Ndata,Ncells,N


  INTEGER ::j,k,i,rowindex,colindex
  REAL(8), ALLOCATABLE,DIMENSION(:) ::X,W
  CHARACTER(LEN=*),PARAMETER :: formatstring='(30(F10.5))' !string to format output
!*****************************************************************


  write(*,*) 'ENTER THE NUMBER OF UNIT CELLS'
  read*,Ncells
  N=Ncells*4
Nsteps=400
Ndata=N*Nsteps
   ALLOCATE(X(Ndata),W(Ndata))

open(unit=4000,file='eigenval.out',status='unknown',action='read')


do i=1,Ndata

read(4000,*)x(i),W(i)

end do




u=1
do j=1,N

do i=j,Ndata,N

!read(1002,*)x(i),W(i)

write(u,'(2(f12.6))')x(i),W(i)

end do

u=u+1
if (u==5.or.u==6)then
u=u+2
end if

end do


End program

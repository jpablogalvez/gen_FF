!======================================================================!
!
! This subroutine performs the exact diagonalization of a real symmetryc         
!  matrix using the DSYEV subroutine from the LAPACK library
!
       subroutine diagonalize(n,matrix,eigenval,eigenvec)
!
       use utils, only: print_end
!
       implicit none
!
! Input/output variables
! 
       real(kind=8),dimension(n,n),intent(in)   ::  matrix 
       real(kind=8),dimension(n,n),intent(out)  ::  eigenvec
       real(kind=8),dimension(n),intent(out)    ::  eigenval
       integer,intent(in)                       ::  n
!
! Local variables
!
       real(kind=8),dimension(:,:),allocatable  ::  tmp
       real(kind=8),dimension(:),allocatable    ::  work
       integer                                  ::  lwork
       integer                                  ::  info
!
! Diagonalizng the input real symmetryc matrix
!
       allocate(tmp(n,n))
       allocate(work(1))
!
       tmp(:,:) = matrix(:,:)
! We set LWORK to -1 in order to get back the optimal size of Work
       call DSYEV('V','U',n,tmp,n,eigenval,work,-1,info) 
!
       if ( info .ne. 0 ) then
         write(*,*)
         write(*,'(2X,68("="))')
         write(*,'(3X,A)') 'ERROR:  Diagonalization failure (DSYEV)'
         write(*,*) 
         write(*,'(3X,A)') 'Error in the calculation of the optima'//  &
                                              'l size of the WORK array'
         write(*,'(2X,68("="))')
         write(*,*) 
!
         deallocate(work,tmp)
!
         call print_end()
       end if
! As lwork was set to -1 WORK(1) returns the optimal LWORK
       lwork = int(work(1)) 
!                      
       deallocate(work)
!
       allocate(work(lwork))
! The matrix_tmp could have changed and it would not be the entering
!  input matrix anymore
       tmp(:,:) = matrix(:,:) 
!
       call DSYEV('V','U',n,tmp,n,eigenval,work,lwork,info)
!
       if ( info .lt. 0 ) then
         write(*,*)
         write(*,'(2X,68("="))')
         write(*,'(3X,A)') 'ERROR:  Diagonalization failure (DSYEV)'
         write(*,*) 
         write(*,'(3X,A)') 'Wrong argument introduced in the DSYEV'//  &
                                                           ' subroutine'
         write(*,'(3X,A,1X,I3)') 'The argument with illegal value '//  &
                                               'is the number',abs(info)
         write(*,'(2X,68("="))')
         write(*,*) 
         call print_end()
       else if ( info .gt. 0 ) then
         write(*,*)
         write(*,'(2X,68("="))')
         write(*,'(3X,A)') 'ERROR:  Diagonalization failure (DSYEV)'
         write(*,*) 
         write(*,'(3X,A)') 'Convergence not reached'
         write(*,'(2X,68("="))')
         write(*,*) 
         call print_end()
       end if
!
       eigenvec(:,:) = tmp(:,:)
!
       deallocate(work,tmp)
!
       return
       end subroutine diagonalize
!
!======================================================================!
!
! This subroutine computes the singular value decomposition of a
!  real general rectangular matrix A(M,N) using the DGESVD subroutine
!  from the LAPACK library
!
! The SVD is written as
!
!   A = U*S*Vt
!
       subroutine SVD(flag,N,M,A,U,S,Vt)
!
       use utils, only: lowercase,                                     &
                        print_end
!
       implicit none
!
! Input/output variables
!
       real(kind=8),dimension(N,M),intent(in)   ::  A
       real(kind=8),dimension(M,M),intent(out)  ::  U
       real(kind=8),dimension(N,N),intent(out)  ::  Vt
       real(kind=8),dimension(N),intent(out)    ::  S
       integer,intent(in)                       ::  N
       integer,intent(in)                       ::  M
       character(len=1)                         ::  flag
!
! Local variables
!
       real(kind=8),dimension(:,:),allocatable  ::  tmp
       real(kind=8),dimension(:),allocatable    ::  work
       integer                                  ::  lwork
       integer                                  ::  info
!
! Computing the singular value decomposition
!
       allocate(tmp(N,N))
       tmp(:,:) = A(:,:)
!
! Finding optimal size for temporary arrays
!
       allocate(work(1))
       lwork = -1
!
       call DGESVD('A','A',M,N,tmp,M,S,U,M,Vt,N,work,lwork,info)
!
       if ( info .ne. 0 ) then
         write(*,*)
         write(*,'(2X,68("="))')
         write(*,'(3X,A)') 'ERROR:  Singular value decomposition f'//  &
                                                       'ailure (DGESVD)'
         write(*,*) 
         write(*,'(3X,A)') 'Error in the calculation of the optima'//  &
                                              'l size of the WORK array'
         write(*,'(2X,68("="))')
         write(*,*) 
! 
         deallocate(work,tmp)
!
         call print_end()
       end if
!
       lwork = int(work(1))
!
       deallocate(work)
!
! Computing SVD
!
       allocate(work(lwork))
! The matrix_tmp could have changed and it would not be the entering  ! FLAG: check if this happens
       tmp(:,:) = A(:,:) 
!
       call DGESVD('A','A',M,N,tmp,M,S,U,M,Vt,N,work,lwork,info)
!
       deallocate(work,tmp)
!
       if ( info .lt. 0 ) then
         write(*,*)
         write(*,'(2X,68("="))')
         write(*,'(3X,A)') 'ERROR:  Singular value decomposition f'//  &
                                                       'ailure (DGESVD)'
         write(*,*) 
         write(*,'(3X,A)') 'Wrong argument introduced in the DSYEV'//  &
                                                           ' subroutine'
         write(*,'(3X,A,1X,I3)') 'The argument with illegal value '//  &
                                               'is the number',abs(info)
         write(*,'(2X,68("="))')
         write(*,*) 
         call print_end()
       else if ( info .gt. 0 ) then
         write(*,*)
         write(*,'(2X,68("="))')
         write(*,'(3X,A)') 'ERROR:  Singular value decomposition f'//  &
                                                       'ailure (DGESVD)'
         write(*,*) 
         write(*,'(3X,A)') 'Convergence not reached'
         write(*,'(2X,68("="))')
         write(*,*) 
         call print_end()
       end if
! 
       if ( lowercase(flag) .eq. 'n' ) then       ! returns V
         Vt = transpose(Vt)
       else if ( lowercase(flag) .eq. 't' ) then  ! returns Vt
         return
       else
         write(*,*)
         write(*,'(2X,68("="))')
         write(*,'(3X,A)') 'ERROR:  Subroutine SVD called incorrectly'
         write(*,*) 
         write(*,'(3X,A)') 'Option '//flag//' not known'
         write(*,'(3X,A)') 'Please, select between "n" and "t"'
         write(*,'(2X,68("="))')
         write(*,*) 
         call print_end() 
       end if
!
       return
       end subroutine SVD
!
!======================================================================!
!
! Solve the linear system A.x = b where A is a NxN matrix
! and x and x are vectors of size N
!
       subroutine linear_solve(N,A,b,x)
!
       implicit none
!
! Input/output variables
! 
       real(kind=8),dimension(N,N),intent(in)  ::  A
       real(kind=8),dimension(N),intent(in)    ::  b
       real(kind=8),dimension(N),intent(out)   ::  x
       integer,intent(in)                      ::  N
!
! Local variables
!
       integer,dimension(:),allocatable        ::  ipiv
       real(kind=8),dimension(:),allocatable   ::  work
       integer                                 ::  lwork
       integer                                 ::  info
!
! Solving the linear system
!
       allocate(ipiv(N),work(N*N))
       lwork = size(work)
!
       x = b
!
       call dsysv('U',N,1,A,N,ipiv,x,N,work,lwork,info)
!
       if (info /= 0) then
         print *,  info
         stop 'error in linear_solve (dsysv)!!'
       endif
!
       return
       end subroutine linear_solve
!
!======================================================================!
!
! This function returns the inverse of the square matrix A(N,N) 
!
       function inverse(N,A) result(B)
!
       implicit none
!
! Input/output variables
!
       real(kind=8),dimension(N,N),intent(in)   ::  A
       real(kind=8),dimension(N,N)              ::  B
       integer,intent(in)                       ::  N
!
! Local variables
!
       integer,dimension(:),allocatable         ::  ipiv
       real(kind=8),dimension(:),allocatable    ::  work
       integer                                  ::  lwork
       integer                                  ::  info
!
! Computing the inverse of the input matrix
!
       allocate(ipiv(N),work(N*N))
       lwork = size(work)
!
       B(:,:) = A(:,:)
!
       call DGETRF(N,N,B,N,ipiv,info)
!
       if (info /= 0) then
         print*,info
         stop 'error in inverse (dgetrf)!!'
       endif
!
       call DGETRI(N,B,N,ipiv,work,lwork,info)
!
       if (info /= 0) then
         print *,  info
         stop 'error in inverse (dgetri)!!'
       endif
!
       deallocate(ipiv,work)
!
       return
       end function inverse
!
!======================================================================!
!
! This function returns the determinant of the square matrix A(N,N) 
!
       subroutine determinant(N,A,det)
!
       use utils, only: print_end
!
       implicit none
!
! Input/output variables
!
       real(kind=8),dimension(N,N),intent(in)   ::  A
       real(kind=8),intent(out)                 ::  det
       integer,intent(in)                       ::  N
!
! Local variables
!
       real(kind=8),dimension(:,:),allocatable  ::  tmp
       integer,dimension(:),allocatable         ::  ipiv
       integer                                  ::  info
       integer                                  ::  i
!
! Computing the inverse of the input matrix
!
       allocate(tmp(N,N),ipiv(N))
!
       tmp(:,:) = A(:,:)
       ipiv(:)  = 0
!
       call DGETRF(N,N,tmp,N,ipiv,info)
!
       if ( info .lt. 0 ) then
         write(*,*)
         write(*,'(2X,68("="))')
         write(*,'(3X,A)') 'ERROR:  LU decomposition failure (DGETRF)'
         write(*,*) 
         write(*,'(3X,A)') 'Wrong argument introduced in the DGETR'//  &
                                                          'F subroutine'
         write(*,'(3X,A,1X,I3)') 'The argument with illegal value '//  &
                                               'is the number',abs(info)
         write(*,'(2X,68("="))')
         write(*,*) 
         call print_end()
       else if ( info .gt. 0 ) then
         write(*,*)
         write(*,'(2X,68("="))')
         write(*,'(3X,A)') 'ERROR:  LU decomposition failure (DGETRF)'
         write(*,*) 
         write(*,'(3X,A)') 'The factor U is exactly singular'
         write(*,'(2X,68("="))')
         write(*,*) 
         call print_end()
       end if
!
       det = 1.0d0
!
       do i = 1, N
         if ( ipiv(i) .eq. i ) then
           det =  det*tmp(i,i)
         else 
           det = -det*tmp(i,i)
         end if
       end do
!
       deallocate(tmp,ipiv)
!
       return
       end subroutine determinant
!
!======================================================================!

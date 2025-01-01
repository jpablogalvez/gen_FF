!======================================================================!
!
       module sorting
       implicit none
!
       contains
!
!======================================================================!
!
! LMRCIVQSORT - Logical Matrix Rows and Columns Integer Vector Quicksort
!
! This subroutine performs the permutations in the interval [LOW,HIGH]
!  of the rows and columns of a logical matrix ADJ(N,N) of dimension N
!  and the elements of an integer vector IVEC(N) based on the order gi-
!  ven in the array REF(N)
!
! Seen on: https://cs.stackexchange.com/questions/104816/implementation-of-quicksort-to-handle-duplicates
!
       recursive subroutine lmrcivqsort(n,ref,adj,ivec,low,high)
!
! Input/output variables
!
       integer,dimension(n),intent(inout)    ::  ref
       integer,dimension(n),intent(inout)    ::  ivec
       logical,dimension(n,n),intent(inout)  ::  adj
       integer,intent(in)                    ::  n
       integer,intent(in)                    ::  low
       integer,intent(in)                    ::  high
!
! Local variables
!
       logical,dimension(n)                  ::  lvaux
       real(kind=4)                          ::  rndm
       integer                               ::  pivot
       integer                               ::  left
       integer                               ::  right
       integer                               ::  upper
       integer                               ::  iaux
!
! Sorting elements with indexes in the interval [low,high]
! --------------------------------------------------------
!
       if ( low .lt. high ) then
! Choice a random pivot (not best performance, but avoids worst-case)
         call random_number(rndm)
         pivot = ref(low + FLOOR((high+1-low)*rndm))
!
! Partitioning the reference array
!
         left  = low
         right = low
         upper = high
! Sorting block of elements with respect to the pivot
         do while ( right .le. upper ) 
           if ( ref(right) .lt. pivot ) then
! Sorting aggregates as they were found
             iaux       = ref(left)
             ref(left)  = ref(right)
             ref(right) = iaux
! Permuting columns of the adjacency matrix
             lvaux(:)     = adj(:,left)
             adj(:,left)  = adj(:,right)
             adj(:,right) = lvaux(:)
! Permuting rows of the adjacency matrix
             lvaux(:)     = adj(right,:)
             adj(right,:) = adj(left,:)
             adj(left,:)  = lvaux(:)
! Permuting elements of the integer vector
             iaux        = ivec(left)
             ivec(left)  = ivec(right)
             ivec(right) = iaux
! Updating partitions
             left  = left  + 1
             right = right + 1
           else if ( ref(right) .gt. pivot ) then
! Sorting aggregates as they were found
             iaux       = ref(upper)
             ref(upper) = ref(right)
             ref(right) = iaux
! Permuting columns of the adjacency matrix
             lvaux(:)     = adj(:,upper)
             adj(:,upper) = adj(:,right)
             adj(:,right) = lvaux(:)
! Permuting rows of the adjacency matrix
             lvaux(:)     = adj(right,:)
             adj(right,:) = adj(upper,:)
             adj(upper,:) = lvaux(:)
! Permuting elements of the integer vector
             iaux        = ivec(upper)
             ivec(upper) = ivec(right)
             ivec(right) = iaux
! Updating partitions
             upper = upper - 1
           else
! Updating partitions
             right = right + 1
           end if
         end do
!  
! Sorting the partitions not containing duplicates
!
         call lmrcivqsort(n,ref,adj,ivec,low,left-1)
         call lmrcivqsort(n,ref,adj,ivec,right,high)
       end if
!
       end subroutine lmrcivqsort
!
!======================================================================!
!
! IQSORT - Integer Quicksort
!
! Seen on: https://cs.stackexchange.com/questions/104816/implementation-of-quicksort-to-handle-duplicates
!
       recursive subroutine iqsort(n,ref,low,high)
!
! Input/output variables
!
       integer,dimension(n),intent(inout)  ::  ref
       integer,intent(in)                  ::  n
       integer,intent(in)                  ::  low
       integer,intent(in)                  ::  high
!
! Local variables
!
       real(kind=4)                        ::  rndm
       integer                             ::  pivot
       integer                             ::  left
       integer                             ::  right
       integer                             ::  upper
       integer                             ::  iaux
!
! Sorting elements with indexes in the interval [low,high]
! --------------------------------------------------------
!
       if ( low .lt. high ) then
! Choice a random pivot (not best performance, but avoids worst-case)
         call random_number(rndm)
         pivot = ref(low + FLOOR((high+1-low)*rndm))
!
! Partitioning the reference array
!
         left  = low
         right = low
         upper = high
! Sorting block of elements with respect to the pivot
         do while ( right .le. upper ) 
           if ( ref(right) .lt. pivot ) then
! Sorting reference array
             iaux       = ref(left)
             ref(left)  = ref(right)
             ref(right) = iaux
! Updating partitions
             left  = left  + 1
             right = right + 1
           else if ( ref(right) .gt. pivot ) then
! Sorting reference array
             iaux       = ref(upper)
             ref(upper) = ref(right)
             ref(right) = iaux
! Updating partitions
             upper = upper - 1
           else
! Updating partitions
             right = right + 1
           end if
         end do
!  
! Sorting the partitions not containing duplicates
!
         call iqsort(n,ref,low,left-1)
         call iqsort(n,ref,right,high)

       end if
!
       end subroutine iqsort
!
!======================================================================!
!
! IINSERTASC - Integer INSERTion sort ASCendent order
!
       subroutine iinsertasc(n,ref)
!
! Input/output variables
!
       integer,dimension(n),intent(inout)  ::  ref
       integer,intent(in)                  ::  n
!
! Local variables
!
       integer                             ::  iaux
       integer                             ::  i,j
!
! Sorting elements with indexes in the reference array
! ----------------------------------------------------
!
       do i = 2, n
         iaux = ref(i)
         j = i - 1
         do while ( j .ge. 1 )
           if ( ref(j) .le. iaux ) exit
             ref(j+1) = ref(j)
             j = j - 1
         end do
         ref(j+1) = iaux
       end do
!
       end subroutine iinsertasc
!
!======================================================================!
!
! IINSERTDESC - Integer INSERTion sort DESCendent order
!
       subroutine iinsertdesc(n,ref)
!
! Input/output variables
!
       integer,dimension(n),intent(inout)  ::  ref
       integer,intent(in)                  ::  n
!
! Local variables
!
       integer                             ::  iaux
       integer                             ::  i,j
!
! Sorting elements with indexes in the reference array
! ----------------------------------------------------
!
       do i = 2, n
         iaux = ref(i)
         j = i - 1
         do while ( j .ge. 1 )
           if ( ref(j) .le. iaux ) exit
             ref(j+1) = ref(j)
             j = j - 1
         end do
         ref(j+1) = iaux
       end do
!
       end subroutine iinsertdesc
!
!======================================================================!
!
! DQSORT - Double precision Quicksort
!
! Seen on: https://cs.stackexchange.com/questions/104816/implementation-of-quicksort-to-handle-duplicates
!
       recursive subroutine dqsort(n,ref,low,high)
!
! Input/output variables
!
       real(kind=8),dimension(n),intent(inout)  ::  ref
       integer,intent(in)                       ::  n
       integer,intent(in)                       ::  low
       integer,intent(in)                       ::  high
!
! Local variables
!
       real(kind=8)                             ::  pivot
       real(kind=8)                             ::  daux
       real(kind=4)                             ::  rndm
       integer                                  ::  left
       integer                                  ::  right
       integer                                  ::  upper
!
! Sorting elements with indexes in the interval [low,high]
! --------------------------------------------------------
!
       if ( low .lt. high ) then
! Choice a random pivot (not best performance, but avoids worst-case)
         call random_number(rndm)
         pivot = ref(low + FLOOR((high+1-low)*rndm))
!
! Partitioning the reference array
!
         left  = low
         right = low
         upper = high
! Sorting block of elements with respect to the pivot
         do while ( right .le. upper ) 
           if ( ref(right) .lt. pivot ) then
! Sorting reference array
             daux       = ref(left)
             ref(left)  = ref(right)
             ref(right) = daux
! Updating partitions
             left  = left  + 1
             right = right + 1
           else if ( ref(right) .gt. pivot ) then
! Sorting reference array
             daux       = ref(upper)
             ref(upper) = ref(right)
             ref(right) = daux
! Updating partitions
             upper = upper - 1
           else 
! Updating partitions
             right = right + 1
           end if
         end do
!  
! Sorting the partitions not containing duplicates
!
         call dqsort(n,ref,low,left-1)
         call dqsort(n,ref,right,high)

       end if
!
       end subroutine dqsort
!
!======================================================================!
!
! IVQSORT - Integer Vector Quicksort
!
! Seen on: https://cs.stackexchange.com/questions/104816/implementation-of-quicksort-to-handle-duplicates
!
       recursive subroutine ivqsort(n,ref,targ,low,high)
!
! Input/output variables
!
       integer,dimension(n),intent(inout)  ::  ref
       integer,dimension(n),intent(inout)  ::  targ
       integer,intent(in)                  ::  n
       integer,intent(in)                  ::  low
       integer,intent(in)                  ::  high
!
! Local variables
!
       real(kind=4)                        ::  rndm
       integer                             ::  pivot
       integer                             ::  left
       integer                             ::  right
       integer                             ::  upper
       integer                             ::  iaux
!
! Sorting elements with indexes in the interval [low,high]
! --------------------------------------------------------
!
       if ( low .lt. high ) then
! Choice a random pivot (not best performance, but avoids worst-case)
         call random_number(rndm)
         pivot = ref(low + FLOOR((high+1-low)*rndm))
!
! Partitioning the reference array
!
         left  = low
         right = low
         upper = high
! Sorting block of elements with respect to the pivot
         do while ( right .le. upper ) 
           if ( ref(right) .lt. pivot ) then
! Sorting reference array
             iaux       = ref(left)
             ref(left)  = ref(right)
             ref(right) = iaux
! Sorting target array according to the reference array 
             iaux        = targ(left)
             targ(left)  = targ(right)
             targ(right) = iaux
! Updating partitions
             left  = left  + 1
             right = right + 1
           else if ( ref(right) .gt. pivot ) then
! Sorting reference array
             iaux       = ref(upper)
             ref(upper) = ref(right)
             ref(right) = iaux
! Sorting target array according to the reference array 
             iaux        = targ(upper)
             targ(upper) = targ(right)
             targ(right) = iaux
! Updating partitions
             upper = upper - 1
           else
! Updating partitions
             right = right + 1
           end if
         end do
!  
! Sorting the partitions not containing duplicates
!
         call ivqsort(n,ref,targ,low,left-1)
         call ivqsort(n,ref,targ,right,high)

       end if
!
       end subroutine ivqsort
!
!======================================================================!
!
! IVVQSORT - Integer Two Vectors Quicksort
!
! Seen on: https://cs.stackexchange.com/questions/104816/implementation-of-quicksort-to-handle-duplicates
!
       recursive subroutine ivvqsort(n,ref,targ1,targ2,low,high)
!
! Input/output variables
!
       integer,dimension(n),intent(inout)  ::  ref
       integer,dimension(n),intent(inout)  ::  targ1
       integer,dimension(n),intent(inout)  ::  targ2
       integer,intent(in)                  ::  n
       integer,intent(in)                  ::  low
       integer,intent(in)                  ::  high
!
!  Local variables
!
       real(kind=4)                        ::  rndm
       integer                             ::  pivot
       integer                             ::  left
       integer                             ::  right
       integer                             ::  upper
       integer                             ::  iaux
!
! Sorting elements with indexes in the interval [low,high]
! --------------------------------------------------------
!
       if ( low .lt. high ) then
! Choice a random pivot (not best performance, but avoids worst-case)
         call random_number(rndm)
         pivot = ref(low + FLOOR((high+1-low)*rndm))
!
! Partitioning the reference array
!
         left  = low
         right = low
         upper = high
! Sorting block of elements with respect to the pivot
         do while ( right .le. upper ) 
           if ( ref(right) .lt. pivot ) then
! Sorting reference array
             iaux       = ref(left)
             ref(left)  = ref(right)
             ref(right) = iaux
! Sorting first target array according to the reference array 
             iaux         = targ1(left)
             targ1(left)  = targ1(right)
             targ1(right) = iaux
! Sorting second target array according to the reference array 
             iaux         = targ2(left)
             targ2(left)  = targ2(right)
             targ2(right) = iaux
! Updating partitions
             left  = left  + 1
             right = right + 1
           else if ( ref(right) .gt. pivot ) then
! Sorting reference array
             iaux       = ref(upper)
             ref(upper) = ref(right)
             ref(right) = iaux
! Sorting first target array according to the reference array 
             iaux         = targ1(upper)
             targ1(upper) = targ1(right)
             targ1(right) = iaux
! Sorting second target array according to the reference array 
             iaux         = targ2(upper)
             targ2(upper) = targ2(right)
             targ2(right) = iaux
! Updating partitions
             upper = upper - 1
           else
! Updating partitions
             right = right + 1
           end if
         end do
!  
! Sorting the partitions not containing duplicates
!
         call ivvqsort(n,ref,targ1,targ2,low,left-1)
         call ivvqsort(n,ref,targ1,targ2,right,high)

       end if
!
       end subroutine ivvqsort
!
!======================================================================!
!
! IdVQSORT - Integer Vector Double precision Vector Quicksort
!
! Seen on: https://cs.stackexchange.com/questions/104816/implementation-of-quicksort-to-handle-duplicates
!
       recursive subroutine ivdvqsort(n,ref,targ,dtarg,low,high)
!
! Input/output variables
!
       real(kind=8),dimension(n),intent(inout)  ::  dtarg
       integer,dimension(n),intent(inout)       ::  ref
       integer,dimension(n),intent(inout)       ::  targ
       integer,intent(in)                       ::  n
       integer,intent(in)                       ::  low
       integer,intent(in)                       ::  high
!
! Local variables
!
       real(kind=8)                             ::  daux
       real(kind=4)                             ::  rndm
       integer                                  ::  pivot
       integer                                  ::  left
       integer                                  ::  right
       integer                                  ::  upper
       integer                                  ::  iaux
!
! Sorting elements with indexes in the interval [low,high]
! --------------------------------------------------------
!
       if ( low .lt. high ) then
! Choice a random pivot (not best performance, but avoids worst-case)
         call random_number(rndm)
         pivot = ref(low + FLOOR((high+1-low)*rndm))
!
! Partitioning the reference array
!
         left  = low
         right = low
         upper = high
! Sorting block of elements with respect to the pivot
         do while ( right .le. upper ) 
           if ( ref(right) .lt. pivot ) then
! Sorting reference array
             iaux       = ref(left)
             ref(left)  = ref(right)
             ref(right) = iaux
! Sorting target array according to the reference array 
             iaux        = targ(left)
             targ(left)  = targ(right)
             targ(right) = iaux
! Sorting double precision target array according to the reference array 
             daux         = dtarg(left)
             dtarg(left)  = dtarg(right)
             dtarg(right) = daux
! Updating partitions
             left  = left  + 1
             right = right + 1
           else if ( ref(right) .gt. pivot ) then
! Sorting reference array
             iaux       = ref(upper)
             ref(upper) = ref(right)
             ref(right) = iaux
! Sorting target array according to the reference array 
             iaux        = targ(upper)
             targ(upper) = targ(right)
             targ(right) = iaux
! Sorting target array according to the reference array 
             daux         = dtarg(upper)
             dtarg(upper) = dtarg(right)
             dtarg(right) = daux
! Updating partitions
             upper = upper - 1
           else
! Updating partitions
             right = right + 1
           end if
         end do
!  
! Sorting the partitions not containing duplicates
!
         call ivdvqsort(n,ref,targ,dtarg,low,left-1)
         call ivdvqsort(n,ref,targ,dtarg,right,high)

       end if
!
       end subroutine ivdvqsort
!
!======================================================================!
!
       end module sorting
!
!======================================================================!

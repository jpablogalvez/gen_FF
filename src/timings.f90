!======================================================================!
!
       module timings
!
       implicit none
!
       save
!
! Declaration of time control variables
!
       real(kind=8)  ::  tcpu      !  Total CPU time
!
! Declaration of system_clock variables 
!
       integer       ::  count_rate
       integer       ::  count_max 
!
       end module timings
!
!======================================================================!

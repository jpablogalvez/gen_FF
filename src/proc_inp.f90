!======================================================================!
!
       module proc_inp
!
       use lengths,  only:  leninp,lenline
!
       implicit none
!
       contains
!
!======================================================================!
!
       subroutine find_key(uni,key,line,iost)
!
       implicit none
!
! Input/output variables
!
       character(len=*),intent(in)         ::  key     !  
       character(len=lenline),intent(out)  ::  line    !  Line read
       integer,intent(out)                 ::  iost    !  Reading status
       integer,intent(in)                  ::  uni     !  Input unit  
!
! Local variables
!
       integer                             ::  keylen  !  Keyword length
       integer                             ::  io      !  Input/Output status
!
! Finding the first line starting with the specified keyword 
!
       keylen = len(key)
       iost   = 0
!
       do
         read(uni,'(A)',iostat=io) line
         if ( io /= 0 ) exit
         line = adjustl(line)
         if ( line(:keylen) .eq. trim(key) ) return
       end do
!
       iost = 1
!
       return
       end subroutine find_key
!
!======================================================================!
!
       subroutine find_last(uni,key,line,iost)
!
       implicit none
!
! Input/output variables
!
       character(len=*),intent(in)         ::  key     !  
       character(len=lenline),intent(out)  ::  line    !  Line read
       integer,intent(in)                  ::  uni     !
       integer,intent(out)                 ::  iost    !  Reading status
!
! Local variables
!
       character(len=lenline)              ::  straux  !  Auxliary string
       integer                             ::  keylen  !  Keyword length
       integer                             ::  io      !  Input/Output status
!
! Finding the last line starting with the specified keyword 
!
       keylen = len_trim(key)
       iost   = 1
!
       do
         read(uni,'(A)',iostat=io) straux
         if ( io /= 0 ) exit
         straux = adjustl(straux)
         if ( straux(:keylen) .eq. trim(key) ) then
           iost = 0 
           line = straux
         end if
       end do
!
       return
       end subroutine find_last
!
!======================================================================!
!
       subroutine split_line(n,inline,str)
!
       use printings,  only:  print_end
!
       implicit none
!
! Input/output variables
!
       character(len=10),dimension(n),intent(out)  ::  str     !  
       character(len=lenline),intent(in)           ::  inline  !  Line read
       integer,intent(in)                          ::  n       !
!
! Local variables
!
       character(len=lenline)                      ::  line    !  Line read
       integer                                     ::  io      !  Scan status
       integer                                     ::  i       !  Index
!
! Splitting input line into strings
!
       line = inline
!
       line = adjustl(line)
       io = scan(line,' ')
!
       i = 0
       do while ( io .ne. 0 )
!
         i = i + 1  
!
         if ( i .le. n ) then
           str(i)  = line(:io-1)
         end if
!
         line = line(io+1:)  
         line = adjustl(line)
         if ( len_trim(line) .eq. 0 ) exit
!
         io = scan(line,' ')
!
       end do
!
       if ( i .gt. n ) then
         write(*,*)
         write(*,'(2X,68("="))')
         write(*,'(3X,A)')    'ERROR:  Line size larger than expected'
         write(*,*)
         write(*,'(3X,A)')    'Error while splitting line'
         write(*,'(3X,A,I3)') 'Maximum number of strings :',n
         write(*,'(3X,A,I3)') 'Number of actual strings  :',i
         write(*,*)
         write(*,'(2X,68("="))')
         write(*,*)  
         call print_end()
       end if
!
       return
       end subroutine split_line
!
!======================================================================!
!
       subroutine count_line(inline,n)
!
       use printings,  only:  print_end
!
       implicit none
!
! Input/output variables
!
       character(len=lenline),intent(in)  ::  inline  !  Line read
       integer,intent(out)                ::  n       !
!
! Local variables
!
       character(len=lenline)             ::  line    !  Line read
       integer                            ::  io      !  Scan status
!
! Splitting input line into strings
! 
       line = inline
!
       line = adjustl(line)
       io = scan(line,' ')
!
       n = 0
       do while ( io .ne. 0 )
!
         n = n + 1  
!         
         line = line(io+1:)  
         line = adjustl(line)
         if ( len_trim(line) .eq. 0 ) return
!
         io = scan(line,' ')
!
       end do
!
       return
       end subroutine count_line
!
!======================================================================!
!
       subroutine print_badread(inp,key)
!
       implicit none
!
! Input/output variables
!
       character(len=*),intent(in)       ::  key     !  Input file name
       character(len=leninp),intent(in)  ::  inp     !  Input file name
!
       write(*,*)
       write(*,'(2X,68("="))')
       write(*,'(3X,A)')      'ERROR:  Missing information in the '//  &
                              'input file'
       write(*,*)
       write(*,'(3X,2(A))')   'Keyword    : ',trim(adjustl(key))
       write(*,'(3X,2(A))')   'Input file : ',trim(inp)
       write(*,'(2X,68("="))')
       write(*,*)
       call exit(0)
!
       return
       end subroutine print_badread
!
!======================================================================!
!
       end module proc_inp
!
!======================================================================!

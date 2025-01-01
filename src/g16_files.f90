!======================================================================!
!
       module g16_files
!
       use proc_inp
       use lengths,  only:  leninp,lenarg,lenline,lentag,lenlab
       use units,    only:  uniinp
!
       implicit none
!
       private
       public  ::  chk_log,                                            &
                   read_log,                                           &
                   readwiberg
!
       contains
!
!======================================================================!
!
       subroutine read_log(inp,nat,coord,lab,znum,mass)
!
       use printings,  only:  print_missinp
!
       implicit none
!
! Input/output variables
!
       character(len=leninp),intent(in)                  ::  inp    ! 
       character(len=lenlab),dimension(nat),intent(out)  ::  lab    !  Atomic labels
       real(kind=8),dimension(3,nat),intent(out)         ::  coord  !
       real(kind=8),dimension(nat),intent(out)           ::  mass   !
       integer,dimension(nat),intent(out)                ::  znum   !  
       integer,intent(in)                                ::  nat    !  
!
! Local variables
!
       integer                                           ::  io     !  Input/Output status
!
! Reading information from Gaussian16 output file
!
       open(unit=uniinp,file=trim(inp),action='read',                  &
            status='old',iostat=io)
       if ( io .ne. 0 ) call print_missinp(inp)
!
! Reading nuclei information
!
       call read_nucinf(inp,nat,lab,znum,mass)
!
! Reading optimized geometry
!
       call read_opt(inp,nat,coord)
!
       close(uniinp)
!
       return
       end subroutine read_log
!
!======================================================================!
!
       subroutine chk_log(inp,nat)
!
       use printings,  only:  print_missinp
!
       implicit none
!
! Input/output variables
!
       character(len=leninp),intent(in)  ::  inp     !  Input file name
       integer,intent(out)               ::  nat     !  Number of atoms
!
! Local variables
!
       character(len=lenline)            ::  line    !  Line read
       character(len=lenarg)             ::  straux  !  Auxiliary string
       integer                           ::  io      !  Input/Output status
!
! Reading information from Gaussian16 output file
!
       open(unit=uniinp,file=trim(inp),action='read',status='old',     &
            iostat=io)
!
       if ( io .ne. 0 ) call print_missinp(inp)
!
! Fatal errors check
!
       call chk_errter(inp)
!
       rewind(uniinp)
!
       call chk_opt(inp)
       call chk_imagi(inp)
!
       rewind(uniinp)
! 
! Reading number of atoms
!
       call find_key(uniinp,'NAtoms=',line,io)
       if ( io .ne. 0 ) call print_badread(inp,'NAtoms=')
!
       read(line,*) straux,nat,line
!
       close(uniinp)
!
       return
       end subroutine chk_log
!
!======================================================================!
!
       subroutine chk_errter(inp)
!
       use printings,  only:  print_end
!
       implicit none
!
! Input/output variables
!
       character(len=leninp),intent(in)  ::  inp     ! 
!
! Local variables
!
       character(len=lenline)            ::  line    !  Line read
       integer                           ::  io      !  Input/Output status
!
! Checking if the calculation therminated abnormally
!
       call find_key(uniinp,'Error termination',line,io)
!
       if ( io .eq. 0 ) then
         write(*,*)
         write(*,'(2X,68("="))')
         write(*,'(3X,A)')      'ERROR:  Gaussian job terminated a'//  &
                                                             'bnormally'
         write(*,*)
         write(*,'(3X,A)')      'Please, check the following input'//  &
                                                                ' file:'
         write(*,*)
         write(*,'(4X,A)')       trim(inp)
         write(*,'(2X,68("="))')
         write(*,*)  
         call print_end()
       end if
!
       return
       end subroutine chk_errter
!
!======================================================================!
!
       subroutine chk_opt(inp)
!
       use printings,  only:  print_end
!
       implicit none
!
! Input/output variables
!
       character(len=*),intent(in)  ::  inp     ! 
!
! Local variables
!
       character(len=lenline)       ::  line    !  Line read
       integer                      ::  io      !  Input/Output status
!
! Checking if the optimization therminated abnormally
! 
       call find_key(uniinp,'-- Stationary point found',line,io)
!
       if ( io .ne. 0 ) then
         write(*,*)
         write(*,'(2X,68("="))')
         write(*,'(3X,A)')      'ERROR:  Geometry optimization ter'//  &
                                                    'minated abnormally'
         write(*,*)
         write(*,'(3X,A)')      'Please, check the following input'//  &
                                                                ' file:'
         write(*,*)
         write(*,'(4X,A)')       trim(inp)
         write(*,'(2X,68("="))')
         write(*,*)  
         call print_end()
       end if
!
       return
       end subroutine chk_opt
!
!======================================================================!
!
       subroutine chk_imagi(inp)
!
       use printings,  only:  print_end
!
       implicit none
!
! Input/output variables
!
       character(len=leninp),intent(in)  ::  inp     !  Input file name
!
! Local variables
!
       character(len=lenline)            ::  line    !  Line read
       integer                           ::  io      !  Input/Output status
!
! Checking if there are imaginary frequencies present
!
       call find_key(uniinp,'****** ',line,io)
!
       if ( io .eq. 0 ) then
         write(*,*)
         write(*,'(2X,68("="))')
         write(*,'(3X,A)')      'ERROR:  Geometry optimization did'//  &
                                       ' not finished in a true minimum'
         write(*,*)
         write(*,'(3X,A)')      'Please, check the following input'//  &
                                                                ' file:'
         write(*,*)
         write(*,'(4X,A)')       trim(inp)
         write(*,'(2X,68("="))')
         write(*,*)  
         call print_end()
       end if
!
       return
       end subroutine chk_imagi
!
!======================================================================!
!
       subroutine read_nucinf(inp,nat,lab,znum,mass)
!
       implicit none
!
! Input/output variables
!
       character(len=leninp),intent(in)                  ::  inp     ! 
       character(len=lenlab),dimension(nat),intent(out)  ::  lab     !  Atom names
       real(kind=8),dimension(nat),intent(out)           ::  mass    !  Atomic masses
       integer,dimension(nat),intent(out)                ::  znum    !  Atomic number
       integer,intent(in)                                ::  nat     !  Number of atoms
!
! Local variables
!
       real(kind=8),dimension(nat)                       ::  dvec    !
       character(len=lenline)                            ::  line    !  Line read
       character(len=lenarg)                             ::  straux  !  Auxiliary string
       integer                                           ::  io      !  Input/Output status
       integer                                           ::  i       !  Index
!
       integer                                           ::  ilower    !
       integer                                           ::  iupper    !
       integer,parameter                                 ::  num = 10  ! 
!
! Reading atomic numbers
!
       do ilower = 1, nat, num
!
         iupper = min(ilower+num-1,nat)
!
         call find_key(uniinp,'AtmWgt=',line,io)
         if ( io .ne. 0 ) call print_badread(inp,'AtmWgt=')
         read(line,*) straux,(mass(i),i=ilower,iupper)
!
         call find_key(uniinp,'AtZNuc=',line,io)
         if ( io .ne. 0 ) call print_badread(inp,'AtZNuc=')
         read(line,*) straux,(dvec(i),i=ilower,iupper)
!
       end do
!
       znum = int(dvec)
!
       call find_key(uniinp,'Distance matrix (angstroms):',line,io)
       if ( io .ne. 0 ) call print_badread(inp,'Distance matrix (a'//  &
                                                           'ngstroms):')
!
       read(uniinp,*) line
!
       do i = 1, nat
         read(uniinp,'(A)') line
         read(line,*) straux,lab(i),line
       end do
!
       return
       end subroutine read_nucinf
!
!======================================================================!
!
       subroutine read_opt(inp,nat,coord)
!
       implicit none
!
! Input/output variables
!
       character(len=*),intent(in)                ::  inp     ! 
       real(kind=8),dimension(3,nat),intent(out)  ::  coord   !
       integer,intent(in)                         ::  nat     !  Number of atoms
!
! Local variables
!
       integer                                    ::  io      !  Input/Output status
!
! Reading optimized geometry
!
       call find_coord(nat,coord,io)
       if ( io .ne. 0 ) call print_badread(inp,'Standard orientation:')
!
       return
       end subroutine read_opt
!======================================================================!
!
       subroutine find_coord(nat,coord,iost)
!
       implicit none
!
! Input/output variables
!
       real(kind=8),dimension(3,nat),intent(out)  ::  coord   !
       integer,intent(in)                         ::  nat     ! 
       integer,intent(out)                        ::  iost    !  Reading status
!
! Local variables
!
       character(len=leninp)                      ::  key     !  
       character(len=lenline)                     ::  line    !  Line read
       integer                                    ::  keylen  !  Keyword length
       integer                                    ::  iat     !
       integer                                    ::  io      !  Input/Output status
!
! Reading coordinates
!
       key    = 'Standard orientation:'
       keylen = len_trim(key)
       iost   = 1
!
       do
         read(uniinp,'(A)',iostat=io) line
         if ( io /= 0 ) exit
         line = adjustl(line)
         if ( line(:keylen) .eq. trim(key) ) then
           iost = 0 
!
           read(uniinp,*) line
           read(uniinp,*) line
           read(uniinp,*) line
           read(uniinp,*) line
!
           do iat = 1, nat
!
             read(uniinp,'(A)') line
! Skipping Center Number column
             line = adjustl(line)
             io   = scan(line,' ')
             line = line(io+1:)
! Skipping Atomic Number column
             line = adjustl(line)
             io   = scan(line,' ')
             line = line(io+1:)
! Skipping Atomic Type column
             line = adjustl(line)
             io   = scan(line,' ')
             line = line(io+1:)
! Reading coordinates
             read(line,*) coord(:,iat)
!
           end do
         end if
       end do
!
       return
       end subroutine find_coord
!
!======================================================================!
!
! READWIBERG - READ WIBERG bond index matrix
!
! This subroutine 
!
       subroutine readwiberg(opt,nat,wiberg)    
!
       implicit none
!
! Input/output variables
!
       character(len=leninp),intent(in)             ::  opt     ! 
       real(kind=8),dimension(nat,nat),intent(out)  ::  wiberg  !  Wiberg bond index matrix
       integer,intent(in)                           ::  nat     !  Number of nodes
!
! Local variables
! 
       character(len=lenline)                       ::  line    !  Line read
       character(len=lenline)                       ::  key     !  
       character(len=lentag)                        ::  str1    !  Auxiliary string
       character(len=lentag)                        ::  str2    !  Auxiliary string
       integer                                      ::  keylen  !  Keyword length
       integer                                      ::  io      !  Input/Output status
       integer                                      ::  ikey    !
       integer                                      ::  i,j     !
!
       integer                                      ::  ilower  !
       integer                                      ::  iupper  !
       integer                                      ::  num     ! 
       integer,parameter                            ::  ncol = 9  ! 
!
! Reading Wiberg bond index matrix
! --------------------------------
!
! Finding the last line starting with the specified keyword 
!
       open(unit=20,file=trim(opt),action='read',status='old',iostat=io)
!
       key    = 'Wiberg bond index matrix in the NAO basis'
       keylen = len_trim(key)
       ikey   = 1
!
       do
         read(20,'(A)',iostat=io) line
         if ( io /= 0 ) exit
         line = adjustl(line)
         if ( line(:keylen) .eq. trim(key) ) then
!
           ikey = 0
!
! Reading the matrix column by column
!
           do ilower = 1, nat, ncol
!
             iupper = min(ilower + ncol - 1,nat)
             num    = iupper - ilower + 1
!
             read(20,'(A)') line
             read(20,'(A)') line
             read(20,'(A)') line
!~              write(*,'(3X,10(9X,I6))') (j,j=ilower,iupper)
!
             do i = 1, nat
!~                do j=ilower,iupper
               read(20,'(A)') line
               read(line,*) str1,str2,(wiberg(i,j),j=ilower,iupper)
!~                  B(j-ilower+1) = Wib(i,j)

!~                end do
!
!~                write(*,'(I7,10F15.8)') i,(B(j),j=1,num)
             end do
!
           end do
!
         end if
       end do
!
       close(20)
!
!~        do i = 1, nat
!~          write(*,'(100F8.4)') (wiberg(i,j),j=1,nat)
!~        end do
!
       if ( ikey .ne. 0 ) call print_badread(opt,'Wiberg bond inde'//  &
                                            'x matrix in the NAO basis')
!
       return
       end subroutine readwiberg
!
!======================================================================!
!
       end module g16_files
!
!======================================================================!

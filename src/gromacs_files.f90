!======================================================================!
!
       module gromacs_files
!
       use proc_inp
       use lengths,  only:  leninp,lenarg,lenline,lenlab
!
       implicit none
!
       private
       public  ::  read_top,                                           &
                   read_bond,                                          &
                   read_angle,                                         &
                   read_dihe,                                          &
                   print_top,                                          &
                   print_bond,                                         &
                   print_angle,                                        &
                   print_dihe,                                         &
                   print_head,                                         &
                   print_tail,                                         &
                   print_pairs,                                        &
                   print_exclusions
!
       contains
!
!======================================================================!
!
       subroutine read_top(top,itype,intop,unitop)
!
       use datatypes,  only:  grotop
!
       use utils,      only:  findcv
       use printings,  only:  print_missinp,print_end
!
       implicit none
!
! Input/output variables
!
       type(grotop),intent(out)                      ::  top      !
       character(len=leninp),intent(in)              ::  intop    !
       integer,intent(in)                            ::  unitop   !
       integer,dimension(:),allocatable,intent(out)  ::  itype    !
!
! Local variables
!
       character(len=lenline)                        ::  line     !  Line read
       character(len=lenlab)                         ::  resname  !
       integer                                       ::  nstr     !
       integer                                       ::  io       !
       integer                                       ::  i        !
!
! Reading Gromacs topology information
! ------------------------------------
!
       open(unit=unitop,file=trim(intop),action='read',                &
            status='old',iostat=io)
       if ( io .ne. 0 ) call print_missinp(intop)
!
! Counting the number of atoms in the molecule
!
       call find_key(unitop,'[ atoms ]',line,io)
       if ( io .ne. 0 ) call print_badread(intop,'[ atoms ]')
!
       top%nat = 0
       do
         read(unitop,'(A)',iostat=io) line
         if ( io /= 0 ) exit
         if ( len_trim(line) .eq. 0 ) cycle
         line = adjustl(line)
         if ( line(1:1) .eq. ';' ) cycle
         if ( line(1:1) .eq. '#' ) cycle
         if ( line(1:1) .eq. '[' ) exit
         top%nat = top%nat + 1
       end do
!
       rewind(unitop)
!
! Allocating atom types information
!
       allocate(top%attype%atname(top%nat),top%attype%bond(top%nat),   &
                top%attype%ptype(top%nat),top%attype%atnum(top%nat),   &
                top%attype%mass(top%nat),top%attype%charge(top%nat),   &
                top%attype%sig(top%nat),top%attype%eps(top%nat))
!
! Allocating atoms information
!
       allocate(top%atom%attype(top%nat),top%atom%residue(top%nat),    &
                top%atom%atom(top%nat),top%atom%mass(top%nat),         &
                top%atom%charge(top%nat),top%atom%cgnr(top%nat),       &
                top%atom%atnr(top%nat),top%atom%resnr(top%nat),        &
                top%atom%itype(top%nat))
!
       allocate(itype(top%nat))     
!
       top%atom%attype(:)   = ''
!
       top%attype%atname(:) = ''
       top%attype%bond(:)   = ''
       top%attype%ptype(:)  = '' 
!
! Reading defaults section
!
       call find_key(unitop,'[ defaults ]',line,io)
       if ( io .ne. 0 ) call print_badread(intop,'[ defaults ]')
!
       do
         read(unitop,'(A)',iostat=io) line
         if ( io /= 0 ) exit
         if ( len_trim(line) .eq. 0 ) cycle
         line = adjustl(line)
         if ( (line(1:1).ne.';') .and. (line(1:1).ne.'#') ) exit ! Find first no blank line with no comments
       end do
!
       read(line,*) top%def%nbfunc,top%def%comrule,top%def%genpair,    &
                    top%def%fudgelj,top%def%fudgeqq
!
       rewind(unitop)
!
! Reading atomtypes section
! 
       call find_key(unitop,'[ atomtypes ]',line,io)
       if ( io .ne. 0 ) call print_badread(intop,'[ atomtypes ]')
!
       top%attype%ntype = 0
       do
!
         read(unitop,'(A)',iostat=io) line
         if ( io /= 0 ) exit
!
         io = scan(line,';')
         if ( io .ne. 0 ) then
           line = line(:io-1)
           line = adjustl(line)
         end if
!
         io = scan(line,'#')
         if ( io .ne. 0 ) then
           line = line(:io-1)
           line = adjustl(line)
         end if
!
         if ( len_trim(line) .eq. 0 ) cycle
         line = adjustl(line)
         if ( line(1:1) .eq. ';' ) cycle
         if ( line(1:1) .eq. '#' ) cycle
         if ( line(1:1) .eq. '[' ) exit
!
         top%attype%ntype = top%attype%ntype + 1
!
         call count_line(line,nstr)
!
         if ( nstr .eq. 6 ) then
           read(line,*) top%attype%atname(top%attype%ntype),             &
                        top%attype%mass(top%attype%ntype),               &
                        top%attype%charge(top%attype%ntype),             &
                        top%attype%ptype(top%attype%ntype),              &
                        top%attype%sig(top%attype%ntype),                &
                        top%attype%eps(top%attype%ntype) 
         else if ( nstr .eq. 7 ) then
           read(line,*) top%attype%atname(top%attype%ntype),             &
                        top%attype%bond(top%attype%ntype),               & 
                        top%attype%mass(top%attype%ntype),               &
                        top%attype%charge(top%attype%ntype),             &
                        top%attype%ptype(top%attype%ntype),              &
                        top%attype%sig(top%attype%ntype),                &
                        top%attype%eps(top%attype%ntype) 
         else
           write(*,'(2X,68("="))')
           write(*,'(3X,A)') 'ERROR:  Invalid format in [ atomtype'//  &
                                          's ] section of topology file'
           write(*,*)
           write(*,'(3X,A)') 'Topology file  : '//trim(intop)
           write(*,'(3X,A)') 'Invalid line   : '//trim(line)
           write(*,*)
           write(*,'(2X,68("="))')
           write(*,*)
           call print_end()
         end if
!
       end do
!
       rewind(unitop)
!
! Reading moleculetype section
! 
       call find_key(unitop,'[ moleculetype ]',line,io)
       if ( io .ne. 0 ) call print_badread(intop,'[ moleculetype ]')
!
       do
         read(unitop,'(A)',iostat=io) line
         if ( io /= 0 ) exit
         if ( len_trim(line) .eq. 0 ) cycle
         line = adjustl(line)
         if ( (line(1:1).ne.';') .and. (line(1:1).ne.'#') ) exit ! Find first no blank line with no comments
       end do
!
       read(line,*) top%mol%resname,top%mol%nrexcl
!
       rewind(unitop)
!
! Reading atoms section
! 
       call find_key(unitop,'[ atoms ]',line,io)
       if ( io .ne. 0 ) call print_badread(intop,'[ atoms ]')
!
       itype(:) = 0
       top%atom%itype(:) = 0
!
       top%atom%nat = top%nat
!
       i = 0
       do while ( i .lt. top%nat )
!
         read(unitop,'(A)',iostat=io) line
         if ( io /= 0 ) exit
         if ( len_trim(line) .eq. 0 ) cycle
         line = adjustl(line)
         if ( line(1:1) .eq. ';' ) cycle
         if ( line(1:1) .eq. '#' ) cycle
!
         i = i + 1
!
         read(line,*) top%atom%atnr(i),top%atom%attype(i),             &
                      top%atom%resnr(i),top%atom%residue(i),           &
                      top%atom%atom(i),top%atom%cgnr(i),               &
                      top%atom%charge(i),top%atom%mass(i) ! TODO: read top mass or use qc output mass
!
         io = findcv(top%attype%ntype,                                 &
                     top%attype%atname(:top%attype%ntype),             &
                     top%atom%attype(i))
!
         if ( io .le. 0 ) then
           write(*,*)
           write(*,'(2X,68("="))')
           write(*,'(3X,A)') 'ERROR:  Atomtype not specified'
           write(*,*)
           write(*,'(3X,A)') 'Please, check the following atom: '//    &
                                                      top%atom%attype(i)
           write(*,*)
           write(*,'(2X,68("="))')
           write(*,*)  
           call print_end()
         end if
!
         itype(i) = io
         top%atom%itype(i) = io
!
       end do
!
       rewind(unitop)
!
! Reading system section
! 
       call find_key(unitop,'[ system ]',line,io)
       if ( io .ne. 0 ) call print_badread(intop,'[ system ]')
!
       do
         read(unitop,'(A)',iostat=io) line
         if ( io /= 0 ) exit
         if ( len_trim(line) .eq. 0 ) cycle
         line = adjustl(line)
         if ( (line(1:1).ne.';') .and. (line(1:1).ne.'#') ) exit ! Find first no blank line with no comments
       end do
!
       read(line,*) top%mol%sysname
!
! Reading molecules section
! 
       call find_key(unitop,'[ molecules ]',line,io)
       if ( io .ne. 0 ) call print_badread(intop,'[ molecules ]')
!
       do
         read(unitop,'(A)',iostat=io) line
         if ( io /= 0 ) exit
         if ( len_trim(line) .eq. 0 ) cycle
         line = adjustl(line)
         if ( (line(1:1).ne.';') .and. (line(1:1).ne.'#') ) exit ! Find first no blank line with no comments
       end do
!
       read(line,*) resname,top%mol%nmol
!
       if ( trim(top%mol%resname) .ne. trim(resname) ) then
         write(*,*)
         write(*,'(2X,68("="))')
         write(*,'(3X,A)') 'ERROR:  Residue name mismatch in input topology'
         write(*,*)
         write(*,'(3X,A)') 'Residue name specified after [ moleculetype ] :',trim(top%mol%resname)
         write(*,'(3X,A)') 'Residue name specified after [ molecules ]    :',trim(resname)
         write(*,*)
         write(*,'(2X,68("="))')
         write(*,*)  
         call print_end()
       end if
!
       close(unitop)
!
       return
       end subroutine read_top
!
!======================================================================!
!
       subroutine read_bond(bonded,intop,unitop)
!
       use datatypes,  only:  grobonded
!
       use printings,  only:  print_missinp
!
       implicit none
!
! Input/output variables
!
       type(grobonded),intent(inout)     ::  bonded  !
       character(len=leninp),intent(in)  ::  intop   !
       integer,intent(in)                ::  unitop  !
!
! Local variables
!
       character(len=lenline)            ::  line    !  Line read
       integer                           ::  io      !
       integer                           ::  i       !
!
! Reading Gromacs bonds information
! ---------------------------------
!
       open(unit=unitop,file=trim(intop),action='read',                &
            status='old',iostat=io)
       if ( io .ne. 0 ) call print_missinp(intop)
!
! Counting the number of bonds in the molecule
!
       call find_key(unitop,'[ bonds ]',line,io)
       if ( io .ne. 0 ) call print_badread(intop,'[ bonds ]')
!
       bonded%nbond = 0
       do
         read(unitop,'(A)',iostat=io) line
         if ( io /= 0 ) exit
         if ( len_trim(line) .eq. 0 ) cycle
         line = adjustl(line)
         if ( line(1:1) .eq. ';' ) cycle
         if ( line(1:1) .eq. '#' ) cycle
         if ( line(1:1) .eq. '[' ) exit
         bonded%nbond = bonded%nbond + 1
       end do
!
       rewind(unitop)
!
! Allocating bonds information
!
       allocate(bonded%bond(bonded%nbond),bonded%kbond(bonded%nbond),  &
                bonded%fbond(bonded%nbond),bonded%ibond(2,bonded%nbond))
!
! Reading bonds section  ! TODO: only valid for harmonic bond and G96 bond interactions
! 
       call find_key(unitop,'[ bonds ]',line,io)
       if ( io .ne. 0 ) call print_badread(intop,'[ bonds ]')
!
       i = 0
       do while ( i .lt. bonded%nbond )
!
         read(unitop,'(A)',iostat=io) line
         if ( io /= 0 ) exit
         if ( len_trim(line) .eq. 0 ) cycle
         line = adjustl(line)
         if ( line(1:1) .eq. ';' ) cycle
         if ( line(1:1) .eq. '#' ) cycle
!
         i = i + 1
!
         read(line,*) bonded%ibond(:,i),bonded%fbond(i),               &
                      bonded%bond(i),bonded%kbond(i)
!
       end do
!
       close(unitop)
!
       return
       end subroutine read_bond
!
!======================================================================!
!
       subroutine read_angle(bonded,intop,unitop)
!
       use datatypes,  only:  grobonded
!
       use printings,  only:  print_missinp
!
       implicit none
!
! Input/output variables
!
       type(grobonded),intent(inout)     ::  bonded  !
       character(len=leninp),intent(in)  ::  intop   !
       integer,intent(in)                ::  unitop  !
!
! Local variables
!
       character(len=lenline)            ::  line    !  Line read
       integer                           ::  io      !
       integer                           ::  i       !
!
! Reading Gromacs angles information
! ----------------------------------
!
       open(unit=unitop,file=trim(intop),action='read',                &
            status='old',iostat=io)
       if ( io .ne. 0 ) call print_missinp(intop)
!
! Counting the number of angles in the molecule
!
       call find_key(unitop,'[ angles ]',line,io)
       if ( io .ne. 0 ) then
         write(*,*)
         write(*,'(2X,68("*"))')
         write(*,'(3X,A)') 'NOTE:  [ angles ] section not found in'//  & ! TODO: check if nat = 2
                                                        ' topology file'
         write(*,*)
         write(*,'(2X,68("*"))')
         write(*,*)  
         close(unitop)
         return
       end if
!
       bonded%nang = 0
       do
         read(unitop,'(A)',iostat=io) line
         if ( io /= 0 ) exit
         if ( len_trim(line) .eq. 0 ) cycle
         line = adjustl(line)
         if ( line(1:1) .eq. ';' ) cycle
         if ( line(1:1) .eq. '#' ) cycle
         if ( line(1:1) .eq. '[' ) exit
         bonded%nang = bonded%nang + 1
       end do
!
       rewind(unitop)
!
! Allocating angles information
!
       allocate(bonded%ang(bonded%nang),bonded%kang(bonded%nang),      &
                bonded%fang(bonded%nang),bonded%iang(3,bonded%nang))
!
! Reading angles section  ! TODO: only valid for harmonic bond angle interactions
! 
       call find_key(unitop,'[ angles ]',line,io)
       if ( io .ne. 0 ) call print_badread(intop,'[ angles ]')
!
       i = 0
       do while ( i .lt. bonded%nang )
!
         read(unitop,'(A)',iostat=io) line
         if ( io /= 0 ) exit
         if ( len_trim(line) .eq. 0 ) cycle
         line = adjustl(line)
         if ( line(1:1) .eq. ';' ) cycle
         if ( line(1:1) .eq. '#' ) cycle
!
         i = i + 1
!
         read(line,*) bonded%iang(:,i),bonded%fang(i),                 & 
                      bonded%ang(i),bonded%kang(i)
!
       end do
!
       close(unitop)
!
       return
       end subroutine read_angle
!
!======================================================================!
!
       subroutine read_dihe(bonded,intop,unitop)
!
       use datatypes, only: grobonded
!
       use printings, only: print_missinp
!
       implicit none
!
! Input/output variables
!
       type(grobonded),intent(inout)                      ::  bonded  !
       character(len=leninp),intent(in)                   ::  intop   !
       integer,intent(in)                                 ::  unitop  !
!
! Local variables
!
       character(len=10),dimension(11)                   ::  str      !  
       character(len=lenline)                            ::  line     !  Line read
       integer                                           ::  io       !
       integer                                           ::  i        !
!
! Reading Gromacs dihedrals information
! -------------------------------------
!
       open(unit=unitop,file=trim(intop),action='read',                &
            status='old',iostat=io)
       if ( io .ne. 0 ) call print_missinp(intop)
!
! Counting the number of dihedrals in the molecule
!
       bonded%ndihe = 0
       do
!
         call find_key(unitop,'[ dihedrals ]',line,io)
         if ( io .ne. 0 ) exit
!
         do
           read(unitop,'(A)',iostat=io) line
           if ( io /= 0 ) exit
           if ( len_trim(line) .eq. 0 ) cycle
           line = adjustl(line)
           if ( line(1:1) .eq. ';' ) cycle
           if ( line(1:1) .eq. '#' ) cycle
           if ( line(1:1) .eq. '[' ) exit
           bonded%ndihe = bonded%ndihe + 1
         end do
!
         if ( io /= 0 ) exit
!
       end do
!
       if ( bonded%ndihe .eq. 0 ) then
         write(*,*)
         write(*,'(2X,68("*"))')
         write(*,'(3X,A)') 'NOTE:  [ dihedrals ] section not foun'//  &
                                                    'd in topology file'
         write(*,*)
         write(*,'(2X,68("*"))')
         write(*,*)  
         close(unitop)
         return
       end if
!
       rewind(unitop)
!
! Allocating dihedrals information
!
       allocate(bonded%dihe(bonded%ndihe),bonded%fdihe(bonded%ndihe),  &
                bonded%kdihe(bonded%ndihe),bonded%multi(bonded%ndihe), &
                bonded%idihe(4,bonded%ndihe))
!
       allocate(bonded%c0(bonded%ndihe),bonded%c1(bonded%ndihe),       &
                bonded%c2(bonded%ndihe),bonded%c3(bonded%ndihe),       &
                bonded%c4(bonded%ndihe),bonded%c5(bonded%ndihe))
!
       bonded%idihe(:,:) = 0
       bonded%dihe(:)    = 0.0d0
       bonded%kdihe(:)   = 0.0d0
       bonded%multi(:)   = 0.0d0
!
       bonded%c0(:) = 0.0d0
       bonded%c1(:) = 0.0d0
       bonded%c2(:) = 0.0d0
       bonded%c3(:) = 0.0d0
       bonded%c4(:) = 0.0d0
       bonded%c5(:) = 0.0d0
!
! Reading dihedrals section
! 
       call find_key(unitop,'[ dihedrals ]',line,io)
       if ( io .ne. 0 ) call print_badread(intop,'[ dihedrals ]')
!
       i = 0
       do while ( i .lt. bonded%ndihe )
!
         read(unitop,'(A)',iostat=io) line
         if ( io /= 0 ) exit
         if ( len_trim(line) .eq. 0 ) cycle
         line = adjustl(line)
         if ( line(1:1) .eq. ';' ) cycle
         if ( line(1:1) .eq. '#' ) cycle
!
         i = i + 1
!
         io = scan(line,';')
         if ( io .ne. 0 ) then
           line = line(:io-1)
           line = adjustl(line)
         end if
!
         io = scan(line,'#')
         if ( io .ne. 0 ) then
           line = line(:io-1)
           line = adjustl(line)
         end if
!
         call split_line(11,line,str)
!
         read(str(1),*) bonded%idihe(1,i)
         read(str(2),*) bonded%idihe(2,i)
         read(str(3),*) bonded%idihe(3,i)
         read(str(4),*) bonded%idihe(4,i)
         read(str(5),*) bonded%fdihe(i)
!
         if ( bonded%fdihe(i) .eq. 1 ) then        !  Proper dihedral
           read(str(6),*) bonded%dihe(i)
           read(str(7),*) bonded%kdihe(i)
           read(str(8),*) bonded%multi(i)
         else if ( bonded%fdihe(i) .eq. 2 ) then   !  Improper dihedral
           read(str(6),*) bonded%dihe(i)
           read(str(7),*) bonded%kdihe(i)
         else if ( bonded%fdihe(i) .eq. 3 ) then   !  Ryckaert-Bellemans dihedral
           read(str(6),*)  bonded%c0(i)
           read(str(7),*)  bonded%c1(i)
           read(str(8),*)  bonded%c2(i)
           read(str(9),*)  bonded%c3(i)
           read(str(10),*) bonded%c4(i)
           read(str(11),*) bonded%c5(i)
         else if ( bonded%fdihe(i) .eq. 4 ) then   !  Periodic improper dihedral
           read(str(6),*) bonded%dihe(i)
           read(str(7),*) bonded%kdihe(i)
         else if ( bonded%fdihe(i) .eq. 5 ) then   !  Fourier dihedral
           read(str(6),*)  bonded%c1(i)
           read(str(7),*)  bonded%c2(i)
           read(str(8),*)  bonded%c3(i)
           read(str(9),*)  bonded%c4(i)
           read(str(10),*) bonded%c5(i)
         else if ( bonded%fdihe(i) .eq. 9 ) then   !  Multiple proper dihedral
           read(str(6),*) bonded%dihe(i)
           read(str(7),*) bonded%kdihe(i)
           read(str(8),*) bonded%multi(i)
         else if ( bonded%fdihe(i) .eq. 10 ) then  !  Restricted dihedral
           read(str(6),*) bonded%dihe(i)
           read(str(7),*) bonded%kdihe(i)
         end if
!
       end do
!
       close(unitop)
!
       return
       end subroutine read_dihe
!
!======================================================================!
!
! This subroutine 
!
       subroutine print_top(uni,nat,itype,mindis,top,dihe,bas,geo,     &
                            intop,topout,fpairs,fexcl,debug)
!
       use datatypes,   only:  grotop,                                 &
                               dihedrals
!
       implicit none
!
! Input/output variables
!
       type(grotop),intent(in)                ::  top      !
       type(dihedrals),intent(in)             ::  dihe     !
       character(len=leninp),intent(in)       ::  bas      !
       character(len=leninp),intent(in)       ::  intop    !
       character(len=leninp),intent(in)       ::  topout   !
       character(len=leninp),intent(in)       ::  geo      !
       integer,dimension(nat,nat),intent(in)  ::  mindis   !
       integer,dimension(nat),intent(in)      ::  itype    !
       integer,intent(in)                     ::  nat      !
       integer,intent(in)                     ::  uni      !
       logical,intent(in)                     ::  fpairs   !  
       logical,intent(in)                     ::  fexcl    !  
!
       logical,intent(in)                     ::  debug    !  Debug mode
!
! Printing Gromacs topology
! -------------------------
!
! Opening output files
!
       open(unit=uni,file=trim(topout),action='write')
!
! Printing tolopogy file header
!
       call print_head(uni,topout,intop,geo,top%def,top%attype,        &
                       top%mol,top%atom)
!
! Printing bonded terms
!
       call print_bond(uni,topout,top%bonded)
       if ( top%bonded%nang .gt. 1 ) call print_angle(uni,topout,top%bonded)
       if ( dihe%ndihe .gt. 1 ) call print_dihe(uni,topout,dihe)
!
! Printing pairs section
!
       if ( fpairs ) call print_pairs(uni,nat,itype,mindis,top%def,    &
                                     top%attype,top%mol,top%atom)
!
! Printing exclusions section
!
       if ( fexcl) call print_exclusions(uni,nat)
!
! Printing topology file tail
!
       call print_tail(uni,top%mol%sysname,top%mol%resname,top%mol%nmol)
!
       close(uni)
!
       return
       end subroutine print_top
!
!======================================================================!
!
       subroutine print_bond(uni,top,bonded)
!
       use datatypes, only: grobonded
!
       implicit none
!
! Input/output variables
!
       type(grobonded)                   ::  bonded  !  
       integer,intent(in)                ::  uni     !
       character(len=leninp),intent(in)  ::  top     !
!
! Local variables
!
       integer                           ::  i       !
!
! Printing Gromacs bonds section
! ------------------------------
!
       write(uni,'(A)') '; Stretchings'     
       write(uni,'(A)') '[ bonds ]'     
!
       do i = 1, bonded%nbond
         if ( bonded%fbond(i) .eq. 1 ) then
           write(uni,'(2X,3(1X,I3),2(3X,F12.4),12X,A,2X,I4)')           &
                      bonded%ibond(:,i),bonded%fbond(i),               & 
                      bonded%bond(i),bonded%kbond(i)!,';',bonded%sbond(i)
         else if ( bonded%fbond(i) .eq. 5 ) then
           write(uni,'(2X,3(1X,I3))') bonded%ibond(:,i),bonded%fbond(i)
         end if
       end do
       write(uni,*)
!
       return
       end subroutine print_bond
!
!======================================================================!
!
       subroutine print_angle(uni,top,bonded)
!
       use datatypes, only: grobonded
!
       implicit none
!
! Input/output variables
!
       type(grobonded)                   ::  bonded  !  
       integer,intent(in)                ::  uni     !
       character(len=leninp),intent(in)  ::  top     !
!
! Local variables
!
       integer                           ::  i       !
!
! Printing Gromacs angles section
! -------------------------------
!
       write(uni,'(A)') '; Bendings'
       write(uni,'(A)') '[ angles ]'
!
       do i = 1, bonded%nang    
         write(uni,'(2X,4(1X,I3),2(3X,F9.4),5X,A,2X,I4)')              &
                         bonded%iang(:,i),bonded%fang(i),              &
                         bonded%ang(i),bonded%kang(i)!,';',bonded%sang(i)
       end do
       write(uni,*)
!
       return
       end subroutine print_angle
!
!======================================================================!
!
       subroutine print_dihe(uni,top,dihe)
!
       use datatypes, only: dihedrals
!
       implicit none
!
! Input/output variables
!
       type(dihedrals)                   ::  dihe  !  
       integer,intent(in)                ::  uni   !
       character(len=leninp),intent(in)  ::  top   !
!
! Local variables
!
       integer                           ::  i,j   !
!
! Printing Gromacs dihedrals section
! ----------------------------------
!
       write(uni,'(A)') '; Torsions'
       write(uni,'(A)') '[ dihedrals ]'
!
       if ( dihe%nimpro .gt. 0 ) then
         do i = 1, dihe%nimpro    
            write(uni,'(2X,5(1X,I3),2(1X,F9.4),5X,A,2X,I4)')           &
                dihe%iimpro(:,i),dihe%fimpro(i),dihe%dimpro(i),        &
                                       dihe%kimpro(i)!,';',dihe%simpro(i)
         end do
       end if
!
       if ( dihe%ninv .gt. 0 ) then
         do i = 1, dihe%ninv    
            write(uni,'(2X,5(1X,I3),2(1X,F9.4),5X,A,2X,I4)')           &
                      dihe%iinv(:,i),dihe%finv(i),dihe%dinv(i),        &
                                           dihe%kinv(i)!,';',dihe%sinv(i)
         end do
       end if
!
       if ( dihe%nrigid .gt. 0 ) then
         do i = 1, dihe%nrigid    
            write(uni,'(2X,5(1X,I3),2(1X,F9.4),5X,A,2X,I4)')           &
                      dihe%irigid(:,i),dihe%frigid(i),dihe%drigid(i),  &
                                       dihe%krigid(i)!,';',dihe%srigid(i)
         end do
       end if
!
       if ( dihe%nflexi .gt. 0 ) then
         do i = 1, dihe%nflexi
           do j = 1, dihe%flexi(i)%ntor     
              write(uni,'(2X,5(1X,I3),2(1X,F9.4),1X,I4,5X,A,2X,I4)')   &
                dihe%flexi(i)%itor(:),dihe%fflexi(i),                  &
                dihe%flexi(i)%tor(j)%phase,dihe%flexi(i)%tor(j)%vtor,  &
                dihe%flexi(i)%tor(j)%multi!,';',dihe%flexi(i)%tor(j)%stor
           end do 
         end do
       end if
       write(uni,*)
!
       return
       end subroutine print_dihe
!
!======================================================================!
!
       subroutine print_head(uni,outp,top,geo,def,attype,mol,atom)
!
       use datatypes, only: grodefaults,                               &
                            groattype,                                 &
                            gromolecule,                               &
                            groatoms
!
       use utils,     only: print_host
       use printings, only: print_missinp
!
       implicit none
!
! Input/output variables
!
       type(grodefaults)                 ::  def      !  Default parameters
       type(groattype)                   ::  attype   !  Atom types information
       type(groatoms)                    ::  atom     !  Definition of the molecule
       type(gromolecule)                 ::  mol      !  Definition of the molecule
!
       integer,intent(in)                ::  uni      !
       character(len=leninp),intent(in)  ::  outp     !
       character(len=leninp),intent(in)  ::  top      !
       character(len=leninp),intent(in)  ::  geo      !
!
! Local variables
!
       integer                           ::  i        !
!
! Printing Gromacs topology header
! --------------------------------
!
! Printing title
!
       if ( len_trim(top) .eq. 0 ) then
         write(uni,'(A)') '; Topology generated from scratch'
       else
         write(uni,'(A)') '; Topology generated from '//trim(top)
       end if
       write(uni,'(A)') '; Equilibrium coordinates taken from '//      &
                                                               trim(geo)
       write(uni,'(A)') '; '//trim(outp)//' created by atomtypes at '  &
                                         //fdate()//' on '//print_host()
       write(uni,*)
!
! Printing defaults section
!
       write(uni,'(A)') '[ defaults ]'
       write(uni,'(A)') '; nbfunc        comb-rule       gen-pairs'//  &
                                                '       fudgeLJ fudgeQQ'
       write(uni,'(I1,13X,I1,13X,A,10X,F6.1,5X,F6.4)') def%nbfunc,     &
                         def%comrule,def%genpair,def%fudgelj,def%fudgeqq
       write(uni,*)
!
! Printing atomtypes section  ! TODO: handle virtual sites
! 
       write(uni,'(A)') '[ atomtypes ]'
       write(uni,'(A)') ';name   bond_type     mass     charge   p'//  &
                                          'type   sigma         epsilon'
       do i = 1, attype%ntype
         write(uni,'(A,5X,A,5X,F7.3,2X,F7.3,3X,A,5X,E12.6,3X,E12.6)')  &
                       attype%atname(i),attype%bond(i),attype%mass(i), &
                       attype%charge(i),attype%ptype(i),               &
                       attype%sig(i),attype%eps(i)
       end do 
       write(uni,*)
!
! Printing moleculetype section
! 
       write(uni,'(A)') '[ moleculetype ]'
       write(uni,'(A)') ';name            nrexcl'
       write(uni,'(A,13X,I1)') mol%resname,mol%nrexcl
       write(uni,*)
!
! Printing atoms section
! 
       write(uni,'(A)') '[ atoms ]'
       write(uni,'(A)') ';   nr  type  resi  res  atom  cgnr     c'//  &
                                                       'harge      mass'
!
       do i = 1, atom%nat
         write(uni,'(1X,I5,1X,A,1X,I4,2(1X,A),1X,I4,4X,F9.6,4X,F9.5)') &
                         atom%atnr(i),atom%attype(i),atom%resnr(i),    &
                         atom%residue(i),atom%atom(i),atom%cgnr(i),    &
                         atom%charge(i),atom%mass(i) ! TODO: read top mass or use qc output mass
       end do                                                            
       write(uni,*)
!
       return
       end subroutine print_head
!
!======================================================================!
!
       subroutine print_tail(uni,sysname,resname,nmol)
!
       implicit none
!
! Input/output variables
!
       integer,intent(in)                ::  uni      !
       integer,intent(in)                ::  nmol     !
       character(len=lenarg),intent(in)  ::  sysname  !
       character(len=lenlab),intent(in)  ::  resname  !
!
! Printing Gromacs topology tail
! ------------------------------
!
       write(uni,'(A)') '[ system ]'
       write(uni,'(1X,A)') trim(sysname)
       write(uni,*)
!
       write(uni,'(A)') '[ molecules ]'
       write(uni,'(A)') '; Compound        nmols'
       write(uni,'(1X,A,7X,I6)') resname,nmol
       write(uni,*)
!
       return
       end subroutine print_tail
!
!======================================================================!
!
       subroutine print_pairs(uni,nat,itype,mindis,def,attype,mol,atom)
!
       use datatypes, only: grodefaults,                               &
                            groattype,                                 &
                            gromolecule,                               &
                            groatoms
!
       implicit none
!
! Input/output variables
!
       type(grodefaults),intent(in)            ::  def      !
       type(groattype),intent(in)              ::  attype   !
       type(gromolecule),intent(in)            ::  mol      !
       type(groatoms),intent(in)               ::  atom     !
       integer,dimension(nat,nat),intent(in)   ::  mindis   !
       integer,dimension(nat),intent(in)       ::  itype    !
       integer,intent(in)                      ::  uni      !
       integer,intent(in)                      ::  nat      !
!
! Input/output variables
!
       real(kind=8)                            ::  V,W      !
       integer                                 ::  i,j      !
!
! Parameters
!
       real(kind=8)                            ::  one = 1.0d0        
!
! Printing pairs section in Gromacs topology
! ------------------------------------------
!
       write(uni,'(A)') '; Nonbonded terms'          
       write(uni,'(A)') '[ pairs ]'
       write(uni,'(A)') '; 1-4 interactions' ! TODO: think if print only 1-4 belonging to flexible dihedrals 
!
       do i = 1, nat-1
         do j = i+1,nat
           if ( mindis(i,j) .eq. 3 ) then
!
             if ( def%comrule .eq. 1 ) then
               stop 'combrule 1 not yet implemented'
             else if ( def%comrule .eq. 2 ) then
               V = 0.5d0*(attype%sig(itype(i)) + attype%sig(itype(j))) 
               W = dsqrt(attype%eps(itype(i))*attype%eps(itype(j)))    &
                                                            *def%fudgelj  ! TODO: joyce uses scaled 1-4 (without fudgelj)
             else if ( def%comrule .eq. 3 ) then
               V = dsqrt(attype%sig(itype(i))*attype%sig(itype(j)))     
               W = dsqrt(attype%eps(itype(i))*attype%eps(itype(j)))    &
                                                            *def%fudgelj  ! TODO: joyce uses scaled 1-4 (without fudgelj)
             end if
!
             write(uni,'(3X,3(1X,I3),3(1X,F8.5),1X,F9.4,1X,F11.5)')    &
                     i,j,2,def%fudgeqq,atom%charge(i),atom%charge(j),V,W
!
           end if
         end do
       end do
!
       write(uni,'(A)') '; other nb interactions'
!
       do i = 1, nat-1
         do j = i+1,nat
           if ( mindis(i,j) .gt. mol%nrexcl ) then
!
             if ( def%comrule .eq. 1 ) then
               stop 'combrule 1 not yet implemented'
             else if ( def%comrule .eq. 2 ) then
               V = 0.5d0*(attype%sig(itype(i)) + attype%sig(itype(j))) 
               W = dsqrt(attype%eps(itype(i))*attype%eps(itype(j))) 
             else if ( def%comrule .eq. 3 ) then
               V = dsqrt(attype%sig(itype(i))*attype%sig(itype(j))) 
               W = dsqrt(attype%eps(itype(i))*attype%eps(itype(j))) 
             end if
!
             write(uni,'(3X,3(1X,I3),3(1X,F8.5),1X,F9.4,1X,F11.5)')    &
                             i,j,2,one,atom%charge(i),atom%charge(j),V,W
!
           end if
         end do
       end do
       write(uni,*)
!
       return
       end subroutine print_pairs
!
!======================================================================!
!
       subroutine print_exclusions(uni,nat)
!
       implicit none
!
! Input/output variables
!
       integer,intent(in)                     ::  uni     !
       integer,intent(in)                     ::  nat     !
!
! Input/output variables
!
       integer                                ::  i,j     !
!
! Printing exclusions section in Gromacs topology
! -----------------------------------------------
!
       write(uni,'(A)') '; Exclusions from default nonbonded'
       write(uni,'(A)') '[ exclusions ]'
       write(uni,'(A)') '; ai        aj'
!
       do i = 1, nat-1
         do j = i+1,nat
           write(uni,'(I6,1X,I6)') i,j
         end do
       end do
       write(uni,*)
!
       return
       end subroutine print_exclusions
!
!======================================================================!
!
       end module gromacs_files
!
!======================================================================!

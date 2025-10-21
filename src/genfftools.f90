!======================================================================!
!
       module genfftools
!
       use lengths,       only:  leninp,lenlab
       use units,         only:  uniic,unideps,unitmp,unidihe,         &
                                 unijoyce,uninb,unitop,unicorr
!
       implicit none
!
       private
       public  ::  genrigidlist,                                       &
                   gencyclelist,                                       &
                   genheavylist,                                       &
                   findarcycles,                                       &
                   bonded2dihe
!
       contains
!
!======================================================================!
!
! GENCYCLELIST - GENerate CYCLE LIST
!
! This subroutine 
!
       subroutine gencyclelist(nat,r,adj,mcycle,ncycle,cycles,lcycle)
!
       implicit none
!
! Input/output variables
!
       logical,dimension(nat,nat),intent(in)   ::  adj     !  Adjacency matrix
       logical,dimension(nat,nat),intent(out)  ::  lcycle  !  Bonds belonging to rings
       integer,dimension(r,nat),intent(in)     ::  cycles  !  Cycles information
       integer,dimension(r),intent(in)         ::  ncycle  !  Number of atoms in each cycle
       integer,intent(in)                      ::  mcycle  !  Number of cycles
       integer,intent(in)                      ::  nat     !  Number of nodes
       integer,intent(in)                      ::  r       !  Cyclic rank
!
! Local variables
! 
       integer                                 ::  icycle  !
       integer                                 ::  i,j     !
       integer                                 ::  ii,jj   !
!
! Generating blacklist with bonds belonging to rings
! --------------------------------------------------
!
       lcycle(:,:) = .FALSE.
!
       do icycle = 1, mcycle           !  Loop over each cycle
         if ( ncycle(icycle) .lt. 8 ) then
           do i = 1, ncycle(icycle)-1  !  Loop over each pair of atoms in cycle icycle
             do j = i+1, ncycle(icycle)
!
               ii = cycles(icycle,i)
               jj = cycles(icycle,j)
!
               if ( adj(ii,jj) .and. (.NOT.lcycle(jj,ii)) ) then
!
                 lcycle(ii,jj) = .TRUE.
                 lcycle(jj,ii) = .TRUE.
!
               end if
!
             end do
           end do
         end if
       end do
!
       return
       end subroutine gencyclelist
!
!======================================================================!
!
! GENRIGIDLIST - GENerate RIGID LIST
!
! This subroutine 
!
       subroutine genrigidlist(nat,wiberg,lcycle,lrigid,laroma,latar)   
!
       implicit none
!
! Input/output variables
!
       real(kind=8),dimension(nat,nat),intent(in)  ::  wiberg  !
       logical,dimension(nat,nat),intent(in)       ::  lcycle  !  Bonds belonging to rings
       logical,dimension(nat,nat),intent(out)      ::  lrigid  !  Rigid bonds
       logical,dimension(nat,nat),intent(out)      ::  laroma  !  Rigid bonds
       logical,dimension(nat),intent(out)          ::  latar   !  Aromatic atom
       integer,intent(in)                          ::  nat     !  Number of nodes
!
! Local variables
! 
       integer                                 ::  i,j     !
       integer                                 ::  ii,jj   !
!
! Generating adjacency matrix with rigid bonds information
! --------------------------------------------------------
!
       lrigid(:,:) = .FALSE.  ! TODO: check if all dihedrals within atoms forming a cycle are planar
       laroma(:,:) = .FALSE.
       latar(:)    = .FALSE.
!
! We consider rigid bond 
!  if Wiberg index is greater than 1.5 it is a double bond
!  if Wiberg index is greater than 1.1 and the bond belongs to a cycle
!   then it belongs to an aromatic cycle
!
       do i = 1, nat-1
         do j = i+1, nat
!
           if ( wiberg(j,i) .ge. 1.5d0 ) then
             lrigid(i,j) = .TRUE.
             lrigid(j,i) = .TRUE.
           else if ( (wiberg(j,i).ge.1.1d0).and.lcycle(j,i) ) then 
             lrigid(i,j) = .TRUE.
             lrigid(j,i) = .TRUE.
!
             laroma(i,j) = .TRUE.
             laroma(j,i) = .TRUE.
!
             latar(i) = .TRUE.
             latar(j) = .TRUE.
           end if
!
         end do
       end do
!
       return
       end subroutine genrigidlist
!
!======================================================================!
!
! GENHEAVYLIST - GENerate HEAVY LIST
!
! This subroutine 
!
       subroutine genheavylist(nat,mass,lheavy)
!
       implicit none
!
! Input/output variables
!
       logical,dimension(nat),intent(out)              ::  lheavy   !  Adjacency matrix
       real(kind=8),dimension(nat),intent(in)          ::  mass     !  Atomic masses
       integer,intent(in)                              ::  nat      !  Number of nodes
!
! Local variables
!
       integer                                         ::  i,j      !  Index
!
! Building H-deplected adjacency matrix
!
       lheavy(:) = .FALSE.
!
       do i = 1, nat
         if ( mass(i) .gt. 3.5 ) lheavy(i) = .TRUE.
       end do
!
       return
       end subroutine genheavylist
!
!======================================================================!
!
! FINDARCYCLES - FIND ARomatic CYCLES
!
! This subroutine 
!
       subroutine findarcycles(nat,r,latar,mcycle,ncycle,cycles,       &
                               maroma,naroma,aroma,marunit,narunit,    &
                               arunit)
!
       use sorting, only: iqsort
!
       implicit none
!
! Input/output variables
!
       logical,dimension(nat),intent(in)       ::  latar    !  Atoms belonging to aromatic rings
       integer,dimension(r,nat),intent(in)     ::  cycles   !  Cycles information
       integer,dimension(r,nat),intent(out)    ::  aroma    !  Aromatic cycles information
       integer,dimension(r,nat),intent(out)    ::  arunit   !  Aromatic units information
       integer,dimension(r),intent(in)         ::  ncycle   !  Number of atoms in each cycle
       integer,dimension(r),intent(out)        ::  naroma   !  Number of atoms in each aromatic cycle
       integer,dimension(r),intent(out)        ::  narunit  !  Number of atoms in each aromatic unit
       integer,intent(in)                      ::  mcycle   !  Number of cycles
       integer,intent(out)                     ::  maroma   !  Number of aromatic cycles
       integer,intent(out)                     ::  marunit  !  Number of aromatic units
       integer,intent(in)                      ::  nat      !  Number of nodes
       integer,intent(in)                      ::  r        !  Cyclic rank
!
! Local variables
! 
       logical,dimension(nat,nat)              ::  lblist   !
       logical,dimension(mcycle)               ::  notvis   !
       logical                                 ::  lcheck   !
       logical                                 ::  lnew     !
       integer,dimension(mcycle)               ::  queue    !
       integer                                 ::  iqueue   !
       integer                                 ::  nqueue   !
       integer                                 ::  i,j,k    !
       integer                                 ::  ii,jj    !
!
! Generating array representation of aromatic cycles
! --------------------------------------------------
!
       maroma = 0
!
       naroma(:)  = -1
       aroma(:,:) = -1
!
       do i = 1, mcycle
!
         lcheck = .TRUE.    
!
         do j = 1, ncycle(i)
           if ( .NOT. latar(cycles(i,j)) ) then
             lcheck = .FALSE.
             exit
           end if
         end do
!
         if ( lcheck ) then
!
           maroma = maroma + 1
!
           naroma(maroma)  = ncycle(i)
           aroma(maroma,:) = cycles(i,:)    
!
         end if
!
       end do
!
! Generating array representation of aromatic units
! -------------------------------------------------
!
       marunit = 0
!
       narunit(:)  = -1
       arunit(:,:) = -1
!
       notvis(:) = .TRUE.
!
! Outer loop over every node
!
       do i = 1, maroma
         if ( notvis(i) ) then
!
           notvis(i) = .FALSE.
! Initializing aromatic units information
           marunit = marunit + 1
!
           narunit(marunit)  = naroma(i)
           arunit(marunit,:) = cycles(i,:)
! Initializing queue
           queue(:) = 0
           queue(1) = i
           iqueue   = 1  ! actual position in the queue
           nqueue   = 2  ! next position in the queue
!
           lblist(:,:) = .FALSE.
!
! Inner loop over queue elements
!
           do while ( iqueue .lt. nqueue )
!
             k = queue(iqueue)
!
! Adding edges in queue cycle to blacklist  
!
             do ii = 1, naroma(k)
               if ( ii .lt. naroma(k) ) then
                 jj = ii + 1
               else
                 jj = 1
               end if
               lblist(aroma(k,ii),aroma(k,jj)) = .TRUE.
               lblist(aroma(k,jj),aroma(k,ii)) = .TRUE.
             end do
!
! Checking if rest of cycles share a common edge
!
             do j = i+1, maroma
               if ( notvis(j) ) then
!
                 lcheck = .FALSE.
!
                 do ii = 1, naroma(j)
                   if ( ii .lt. naroma(j) ) then
                     jj = ii + 1
                   else
                     jj = 1
                   end if
                   if ( lblist(aroma(j,ii),aroma(j,jj)) ) then
                     lcheck = .TRUE.
                     exit
                   end if
                 end do
!
                 if ( lcheck ) then
!
! Updating queue
!
                   notvis(j)     = .FALSE.
                   queue(nqueue) = j
                   nqueue        = nqueue + 1
!
! Updating aromatic units information
!
                   do jj = 1, naroma(j)
                     lnew = .TRUE.
                     do ii = 1, narunit(marunit)
                       if ( arunit(marunit,ii) .eq. aroma(j,jj) ) then                       
                         lnew = .FALSE.
                         exit
                       end if
                     end do
                     if ( lnew ) then
                       narunit(marunit) = narunit(marunit) + 1
                       arunit(marunit,narunit(marunit)) = aroma(j,jj)
                     end if
                   end do
!
                 end if
!
               end if
             end do
!
             iqueue = iqueue + 1
!
           end do
!
         end if
       end do
!
       do i = 1, marunit
         call iqsort(narunit(i),arunit(i,:),1,narunit(i))
       end do
!
       return
       end subroutine findarcycles
!
!======================================================================!
!
! BONDED2DIHE - BONDED TO DIHEdral datatype
!
! This subroutine 
!
       subroutine bonded2dihe(ndihe,bonded,dihe,nat,adj)
!
       use datatypes,  only:  grobonded,                               &
                              dihedrals
!
       use printings,  only:  print_end
!
       implicit none
!
! Input/output variables
!
       type(grobonded),intent(in)             ::  bonded   !
       type(dihedrals),intent(out)            ::  dihe     !
       logical,dimension(nat,nat),intent(in)  ::  adj      !
       integer,intent(in)                     ::  ndihe    !
       integer,intent(in)                     ::  nat      !
!
! Local variables
!
       logical,dimension(ndihe)               ::  ldihe    !
       logical,dimension(ndihe)               ::  visited  !
       integer                                ::  i,j,k    !
!
!  Setting DIHE datatype from BONDED datatype
! -------------------------------------------
!
! Allocating information
!
       allocate(dihe%iimpro(4,ndihe),dihe%dimpro(ndihe),               &
                dihe%kimpro(ndihe),dihe%fimpro(ndihe))
!
       allocate(dihe%iinv(4,ndihe),dihe%dinv(ndihe),                   &
                dihe%kinv(ndihe),dihe%finv(ndihe))
!
       allocate(dihe%irigid(4,ndihe),dihe%drigid(ndihe),               &
                dihe%krigid(ndihe),dihe%frigid(ndihe))
!
       allocate(dihe%flexi(ndihe))
       allocate(dihe%iflexi(4,ndihe),dihe%dflexi(ndihe),               &
                dihe%kflexi(ndihe),dihe%fflexi(ndihe))
!
       dihe%ndihe  = ndihe 
       dihe%nflexi = 0 
       dihe%nrigid = 0 
       dihe%nimpro = 0 
       dihe%ninv   = 0 
!
       visited(:) = .FALSE.
!
! Setting up improper, inversion, and flexible dihedral
!
       do i = 1, ndihe
!
         if ( visited(i) ) cycle
!         
         if ( bonded%fdihe(i) .eq. 1 ) then        !  Proper dihedral
!
           visited(i) = .TRUE.
!
           dihe%nflexi = dihe%nflexi + 1
!
           dihe%flexi(dihe%nflexi)%ntor = 1
           allocate(dihe%flexi(dihe%nflexi)%tor(1)) 
!
           dihe%fflexi(dihe%nflexi)     = 1
!
           dihe%flexi(dihe%nflexi)%itor(:) = bonded%idihe(:,i)
           dihe%iflexi(:,dihe%nflexi)      = bonded%idihe(:,i)
!
           dihe%flexi(dihe%nflexi)%tor(1)%phase = bonded%dihe(i)
           dihe%flexi(dihe%nflexi)%tor(1)%vtor  = bonded%kdihe(i)
           dihe%flexi(dihe%nflexi)%tor(1)%multi = bonded%multi(i)
!
         else if ( bonded%fdihe(i) .eq. 2 ) then   !  Improper dihedral
!
! Checking if the improper dihedral corresponds to an 
!  out of plane bending or a double (rigid) bond
!
           visited(i) = .TRUE.
!
           if ( (adj(bonded%idihe(1,i),bonded%idihe(2,i))) .and.       &
                (adj(bonded%idihe(2,i),bonded%idihe(3,i))) .and.       &
                (adj(bonded%idihe(3,i),bonded%idihe(4,i))) ) then
!
             dihe%nrigid = dihe%nrigid + 1
!
             dihe%frigid(dihe%nrigid)   = 2
             dihe%irigid(:,dihe%nrigid) = bonded%idihe(:,i)
             dihe%drigid(dihe%nrigid)   = bonded%dihe(i)
             dihe%krigid(dihe%nrigid)   = bonded%kdihe(i)
!
           else
!
             dihe%nimpro = dihe%nimpro + 1
!
             dihe%fimpro(dihe%nimpro)   = 2
             dihe%iimpro(:,dihe%nimpro) = bonded%idihe(:,i)
             dihe%dimpro(dihe%nimpro)   = bonded%dihe(i)
             dihe%kimpro(dihe%nimpro)   = bonded%kdihe(i)
!
           end if
!
         else if ( bonded%fdihe(i) .eq. 3 ) then   !  Ryckaert-Bellemans dihedral
!
! TODO: convert to linear combination of proper dihedrals
!
           stop 'Ryckaert-Bellemans dihedral not supported'
!
         else if ( bonded%fdihe(i) .eq. 4 ) then   !  Periodic improper dihedral
!
           stop 'Periodic improper dihedral not supported'
!
         else if ( bonded%fdihe(i) .eq. 5 ) then   !  Fourier dihedral
!
! TODO: convert to linear combination of proper dihedrals
!
           stop 'Fourier dihedral not supported'
!
         else if ( bonded%fdihe(i) .eq. 9 ) then   !  Multiple proper dihedral
!
! Finding number of terms in the Fouier expansion and unique dihedrals
!
           visited(i) = .TRUE.
!
           dihe%nflexi = dihe%nflexi + 1
!
           dihe%fflexi(dihe%nflexi)   = 9
           dihe%iflexi(:,dihe%nflexi) = bonded%idihe(:,i)
!
           dihe%flexi(dihe%nflexi)%itor(:) = bonded%idihe(:,i)
!
           ldihe(:) = .FALSE.
           ldihe(i) = .TRUE.
!
! Counting the proper dihedrals functions sharing the same quadruplet
!
           dihe%flexi(dihe%nflexi)%ntor = 1
           do j = i+1, ndihe 
!
             if ( visited(j) ) cycle
             if ( (bonded%idihe(1,i).eq.bonded%idihe(1,j)) .and.       &
                  (bonded%idihe(2,i).eq.bonded%idihe(2,j)) .and.       &
                  (bonded%idihe(3,i).eq.bonded%idihe(3,j)) .and.       &
                  (bonded%idihe(4,i).eq.bonded%idihe(4,j)) ) then
!
               visited(j) = .TRUE.
               ldihe(j)   = .TRUE.
!
               dihe%flexi(dihe%nflexi)%ntor = dihe%flexi(dihe%nflexi)%ntor + 1
!                 
             end if
!
           end do
!
           allocate(dihe%flexi(dihe%nflexi)%tor(dihe%flexi(dihe%nflexi)%ntor)) 
!
! Setting multiple dihedrals sharing the same quadruplet 
!
           k = 0
           do j = 1, ndihe
             if ( ldihe(j) ) then
               k = k + 1
!
               dihe%flexi(dihe%nflexi)%tor(k)%phase = bonded%dihe(j)
               dihe%flexi(dihe%nflexi)%tor(k)%vtor  = bonded%kdihe(j)
               dihe%flexi(dihe%nflexi)%tor(k)%multi = bonded%multi(j)
             end if
           end do
!
           if ( k .ne. dihe%flexi(dihe%nflexi)%ntor ) then
             write(*,*)
             write(*,'(2X,68("="))')
             write(*,'(3X,A)') 'ERROR:  Datatype transformation is wrong'
             write(*,*)
             write(*,'(3X,A,4I4)') 'Target quadruplet :',bonded%idihe(:,i)
             write(*,'(3X,A,I3)')  'Number of expected terms :',dihe%flexi(dihe%nflexi)%ntor
             write(*,'(3X,A,I3)')  'Number of terms found    :',k
             write(*,'(2X,68("="))')
             write(*,*)
             call print_end()  
           end if
!
         else if ( bonded%fdihe(i) .eq. 10 ) then  !  Restricted dihedral
!
           stop 'Restricted dihedral not supported'
!
         end if
!
       end do
!
       return
       end subroutine bonded2dihe
!
!======================================================================!
!
       end module genfftools
!
!======================================================================!

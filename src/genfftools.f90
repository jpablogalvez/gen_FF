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
                   genmethlist,                                        &
                   genquad,                                            &
                   findarcycles,                                       &
                   bonded2dihe,                                        &
                   selectquad,                                         &
                   screenquad,                                         &
                   calc_angle,                                         &
                   diedro
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
       subroutine genrigidlist(nat,r,wiberg,coord,cycles,ncycle,       &
                               mcycle,lcycle,lrigid,laroma,latar)   
!
       implicit none
!
! Input/output variables
!
       real(kind=8),dimension(nat,nat),intent(in)  ::  wiberg  !
       real(kind=8),dimension(3,nat),intent(in)    ::  coord   !
       logical,dimension(nat,nat),intent(in)       ::  lcycle  !  Bonds belonging to rings
       logical,dimension(nat,nat),intent(out)      ::  lrigid  !  Rigid bonds
       logical,dimension(nat,nat),intent(out)      ::  laroma  !  Rigid bonds
       logical,dimension(nat),intent(out)          ::  latar   !  Aromatic atom
       integer,dimension(r,nat),intent(in)         ::  cycles  !
       integer,dimension(r),intent(in)             ::  ncycle  !
       integer,intent(in)                          ::  nat     !  Number of nodes
       integer,intent(in)                          ::  r       !  Graph rank
       integer,intent(in)                          ::  mcycle  !  Number of cycles
!
! Local variables
! 
       real(kind=8)                                ::  daux    !
       logical                                     ::  lcheck  !
       integer,dimension(4)                        ::  ivaux   !
       integer                                     ::  i,j,k   !
       integer                                     ::  ii      !
!
       real(kind=8),parameter                      ::  pi =  4*atan(1.0_8) 

!
! Generating adjacency matrix with rigid bonds information
! --------------------------------------------------------
!
       lrigid(:,:) = .FALSE.  ! TODO: check if all dihedrals within atoms forming a cycle are planar
       laroma(:,:) = .FALSE.
       latar(:)    = .FALSE.
!
! We consider rigid bond 
!  if Wiberg index is greater than 1.1 and the bond belongs to a cycle
!   then it belongs to an aromatic cycle
!  if Wiberg index is greater than 1.5 then it is a double bond
!
       do i = 1, nat-1
         do j = i+1, nat
!
           if ( (wiberg(j,i).ge.1.1d0).and.lcycle(j,i) ) then 
             lrigid(i,j) = .TRUE.
             lrigid(j,i) = .TRUE.
!
             laroma(i,j) = .TRUE.
             laroma(j,i) = .TRUE.
!
             latar(i) = .TRUE.
             latar(j) = .TRUE.
           else if ( wiberg(j,i) .ge. 1.5d0 ) then
             lrigid(i,j) = .TRUE.
             lrigid(j,i) = .TRUE.
           end if
!
         end do
       end do
!
! Bonds in fully planar rings are also considered rigid/aromatic
!
       do i = 1, mcycle
         lcheck = .TRUE.
         do j = 1, ncycle(i)
!
           do k = -1, 2
             ii = modulo(j + k - 1, ncycle(i)) + 1
             ivaux(k+2) = cycles(i,ii)
           end do
!
           call Diedro(coord(:,ivaux(1)),coord(:,ivaux(2)),            &
                       coord(:,ivaux(3)),coord(:,ivaux(4)),daux)
           daux = daux*180.0d0/pi
!
           if ( abs(daux) .gt. 25.0d0 ) then
             lcheck = .FALSE.
             exit
           end if
!
         end do
! 
         if ( lcheck ) then
           do j = 1, ncycle(i)
             if ( j .lt. ncycle(i) ) then
               k = j + 1
             else
               k = 1
             end if
!
             lrigid(cycles(i,k),cycles(i,j)) = .TRUE.
             lrigid(cycles(i,j),cycles(i,k)) = .TRUE.
!
             laroma(cycles(i,k),cycles(i,j)) = .TRUE.
             laroma(cycles(i,j),cycles(i,k)) = .TRUE.
!
             latar(cycles(i,k)) = .TRUE.
             latar(cycles(i,j)) = .TRUE.
!
           end do
         end if
!
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
       integer                                         ::  i        !  Index
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
! GENMETHLIST - GENerate METHyl LIST
!
! This subroutine 
!
       subroutine genmethlist(nat,adj,znum,dihe,nflexi,lch3,ich3,debug)
!
       use datatypes,  only:  dihedrals
!
       implicit none
!
! Input/output variables
!
       type(dihedrals),intent(in)               ::  dihe    !
       logical,dimension(nat,nat),intent(in)    ::  adj     !  Adjacency matrix
       logical,dimension(nflexi),intent(out)    ::  lch3    !
       integer,dimension(3,nflexi),intent(out)  ::  ich3    !
       integer,dimension(nat),intent(in)        ::  znum    !
       integer,intent(in)                       ::  nat     !  
       integer,intent(in)                       ::  nflexi  !  Number of nodes
       logical,intent(in)                       ::  debug   !
!
! Local variables
!
       logical                                  ::  match   !
       integer                                  ::  idx1    !
       integer                                  ::  idx2    !
       integer                                  ::  idx3    !
       integer                                  ::  i,j     !  Index
!
! Finding torsions associated to CH3 rotations
!
       lch3(:) = .FALSE.
!
       idx3 = -1
!
       do i = 1, dihe%nflexi
!
         idx1 = dihe%iflexi(2,i)
         idx2 = dihe%iflexi(3,i)
!
         match = .TRUE.
!
         if ( znum(idx1) .eq. 6 ) then
           do j = 1, nat
             if ( j .eq. idx2 ) cycle
             if ( adj(idx1,j) ) then
               if ( znum(j) .ne. 1 ) then
                 match = .FALSE.
                 exit
               end if                 
             end if
           end do
           idx3 = dihe%iflexi(4,i)
           if ( match ) then
             lch3(i) = .TRUE.
             ich3(1,i) = idx1
             ich3(2,i) = idx2
             ich3(3,i) = idx3  
             cycle
           end if
         end if
!
         if ( znum(idx2) .eq. 6 ) then
           do j = 1, nat
             if ( j .eq. idx1 ) cycle
             if ( adj(idx2,j) ) then
               if ( znum(j) .ne. 1 ) then
                 match = .FALSE.
                 exit
               end if                 
             end if
           end do
           idx3 = dihe%iflexi(1,i) 
           if ( match ) then
             lch3(i) = .TRUE.
             ich3(1,i) = idx2
             ich3(2,i) = idx1
             ich3(3,i) = idx3  
           end if
         end if
! 
       end do
!
       if ( debug ) then
         write(*,'(1X,A)') 'List of CH3 rotations'
         write(*,'(1X,A)') '---------------------'
         do i = 1, dihe%nflexi
           if ( lch3(i) ) write(*,'(1X,I4,1X,A,4I4,A,4I4)')            &
                                  i,'CH3 rotation',dihe%iflexi(:,i),   &
                                             ' with backbone ',ich3(:,i)
         end do
         write(*,*)
       end if
!
       return
       end subroutine genmethlist
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
       dihe%ntor = ndihe
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
! SELECTQUAD - SELECT QUADruplet
!
! This subroutine 
!
       subroutine selectquad(nat,ndihe,dihe,idat,nidat,debug)
!
       use datatypes,  only:  dihedrals
!
       use printings,  only:  print_end
!
       implicit none
!
! Input/output variables
!
       type(dihedrals),intent(inout)      ::  dihe    !
       integer,dimension(nat),intent(in)  ::  idat    !
       integer,intent(in)                 ::  nat     !
       integer,intent(in)                 ::  nidat   !
       integer,intent(in)                 ::  ndihe   !
       logical,intent(in)                 ::  debug   !
!
! Local variables
!
       real(kind=8),dimension(ndihe)      ::  dtmp     !
       logical,dimension(ndihe)           ::  visited  !
       integer,dimension(4,ndihe)         ::  itmp     !
       integer,dimension(ndihe)           ::  ftmp     !
       integer,dimension(ndihe)           ::  iddihe   !
       integer,dimension(4)               ::  vaux     !
       integer                            ::  ntmp     !
       integer                            ::  i,j,k    !
!
!  Selecting quadruplets with unique identifiers
! ----------------------------------------------
!
       do i = 1, dihe%nquad
         vaux(:) = dihe%iquad(:,i)
         iddihe(i) = min(idat(vaux(2)),idat(vaux(3)))*nidat            &
                     + max(idat(vaux(2)),idat(vaux(3)))
       end do
!
       visited(:) = .FALSE.
!
       itmp(:,:) = -1
       dtmp(:)   = 0.0d0
       ftmp(:)   = -1
       ntmp      = 0
!
       k = 0
       do i = 1, dihe%nquad
!
         if ( visited(i) ) cycle
         k = k + 1
         visited(i) = .TRUE.
!
         itmp(:,k) = dihe%iquad(:,i)
         dtmp(k)   = dihe%dquad(i)
         ftmp(k)   = dihe%fquad(i)
         ntmp = k      
!
         do j = i+1, dihe%nquad
           if ( visited(j) ) cycle
           if ( iddihe(i) .eq. iddihe(j) ) then
             visited(j) = .TRUE.
             exit
           end if
         end do
       end do
!
       dihe%iquad(:,:) = itmp(:,:)
       dihe%dquad(:)   = dtmp(:)
       dihe%fquad(:)   = ftmp(:)
       dihe%nquad = ntmp  
!
       if ( debug ) then
         write(*,'(1X,A)') 'Unique quadruplets'
         write(*,'(1X,A)') '------------------'
         do i = 1, dihe%nquad
           write(*,'(1X,A,4I4)') 'Quadruplet',dihe%iquad(:,i)
         end do
         write(*,*)
       end if
!
       return
       end subroutine selectquad
!
!======================================================================!
!
! SCREENQUAD - SCREEN QUADruplet
!
! This subroutine 
!
       subroutine screenquad(dihe,nflexi,lch3,ich3,ndihe,nat,idat)
!
       use datatypes,  only:  dihedrals
!
       use printings,  only:  print_end
!
       implicit none
!
! Input/output variables
!
       type(dihedrals),intent(inout)           ::  dihe    !
       logical,dimension(nflexi),intent(in)    ::  lch3    !
       integer,dimension(3,nflexi),intent(in)  ::  ich3    !
       integer,dimension(nat),intent(in)       ::  idat    !
       integer,intent(in)                      ::  nflexi  !
       integer,intent(in)                      ::  nat     !
       integer,intent(in)                      ::  ndihe   !
!
! Local variables
!
       type(dihedrals)                         ::  tmpdihe  !
       integer,dimension(4)                    ::  ivaux1   !
       integer,dimension(4)                    ::  ivaux2   !
       integer                                 ::  i,j,k    !
! 
!  Including only principal quadruplets 
! -------------------------------------
!
! Generating new representation with a reduced set of quadruplets
!
       allocate(tmpdihe%iflexi(4,ndihe),tmpdihe%dflexi(ndihe),         &
                tmpdihe%kflexi(ndihe),tmpdihe%fflexi(ndihe),           &
                tmpdihe%flexi(ndihe))
!
       tmpdihe%iflexi(:,:) = -1
       tmpdihe%dflexi(:)   = 999.0d0
       tmpdihe%kflexi(:)   = 0.0d0
       tmpdihe%fflexi(:)   = -1
!
       k = 0
       do i = 1, dihe%nflexi
         ivaux1(:) = dihe%iflexi(:,i)
         if ( lch3(i) ) then
!
write(*,*) 'Dihedral',i,'is a CH3 rotation',dihe%iflexi(:,i)
           do j = 1, dihe%nquad
             ivaux2(:) = dihe%iquad(:,j)
             if ( ((idat(ich3(1,i)).eq.idat(ivaux2(2))).and.             &
                   (idat(ich3(2,i)).eq.idat(ivaux2(3))).and.             &
                   (idat(ich3(3,i)).eq.idat(ivaux2(4)))                  &
             .or. ((idat(ich3(1,i)).eq.idat(ivaux2(3))).and.             &
                   (idat(ich3(2,i)).eq.idat(ivaux2(2))).and.             &
                   (idat(ich3(3,i)).eq.idat(ivaux2(1))))) ) then
write(*,*) '  Overlap with princial quadruplet',j,':',dihe%iquad(:,j) 
write(*,*)
!
               k = k + 1
!
               tmpdihe%iflexi(:,k) = dihe%iflexi(:,i)
               tmpdihe%dflexi(k)   = dihe%dflexi(i)
               tmpdihe%kflexi(k)   = dihe%kflexi(i)
               tmpdihe%fflexi(k)   = dihe%fflexi(i)
!
               tmpdihe%flexi(k)%ntor = dihe%flexi(i)%ntor
               allocate(tmpdihe%flexi(k)%tor(tmpdihe%flexi(k)%ntor))
!
               tmpdihe%flexi(k)%itor(:) = dihe%flexi(i)%itor(:)
               tmpdihe%flexi(k)%tor     = dihe%flexi(i)%tor
! 
               exit
!
             end if
           end do
!
         else
!
           do j = 1, dihe%nquad
             ivaux2(:) = dihe%iquad(:,j)
             if ( ((ivaux1(1).eq.ivaux2(1)).and.                       &
                   (ivaux1(2).eq.ivaux2(2)).and.                       &
                   (ivaux1(3).eq.ivaux2(3)).and.                       &
                   (ivaux1(4).eq.ivaux2(4))) ) then
!
               k = k + 1
!
               tmpdihe%iflexi(:,k) = dihe%iflexi(:,i)
               tmpdihe%dflexi(k)   = dihe%dflexi(i)
               tmpdihe%kflexi(k)   = dihe%kflexi(i)
               tmpdihe%fflexi(k)   = dihe%fflexi(i)
!
               tmpdihe%flexi(k)%ntor = dihe%flexi(i)%ntor
               allocate(tmpdihe%flexi(k)%tor(tmpdihe%flexi(k)%ntor))
!
               tmpdihe%flexi(k)%itor(:) = dihe%flexi(i)%itor(:)
               tmpdihe%flexi(k)%tor     = dihe%flexi(i)%tor
!
             end if
           end do
!         
         end if
       end do
!
       tmpdihe%nflexi = k
!
! Storing new representation in the original one
!
       do i = 1, dihe%nflexi
         deallocate(dihe%flexi(i)%tor)
       end do
!
       tmpdihe%ntor = 0
       do i = 1, tmpdihe%nflexi
         allocate(dihe%flexi(i)%tor(tmpdihe%flexi(i)%ntor))
         dihe%flexi(i) = tmpdihe%flexi(i)
       end do
!
       dihe%nflexi = tmpdihe%nflexi
! 
       dihe%iflexi(:,:) = tmpdihe%iflexi(:,:)
       dihe%dflexi(:)   = tmpdihe%dflexi(:)
       dihe%kflexi(:)   = tmpdihe%kflexi(:)
       dihe%fflexi(:)   = tmpdihe%fflexi(:)
!
       return
       end subroutine screenquad
!
!======================================================================!
!
! GENQUAD - GENerate QUADruplet
!
! This subroutine 
!
       subroutine genquad(nat,znum,dihed,ndihe,debug)
!
       use datatypes,  only:  dihedrals
!
       implicit none
!
! Input/output variables
!
       type(dihedrals),intent(inout)      ::  dihed    !
       integer,dimension(nat),intent(in)  ::  znum     !
       integer,intent(in)                 ::  nat      !
       integer,intent(in)                 ::  ndihe    !
       logical,intent(in)                 ::  debug    !
!
! Local variables
!
       type(dihedrals)                    ::  tmpdi    !            
       real(kind=8)                       ::  daux
       logical,dimension(ndihe)           ::  visited  !
       logical                            ::  flag     !
       integer,dimension(ndihe)           ::  imap     !  
       integer,dimension(ndihe)           ::  iflexi   !  
       integer,dimension(ndihe)           ::  neqquad  !  
       integer,dimension(4)               ::  vaux     !  
       integer,dimension(4)               ::  vaux1    !  
       integer,dimension(4)               ::  vaux2    !  
       integer                            ::  meqquad  !  
       integer                            ::  nmap     !
       integer                            ::  i,j,k    !
!
!  Generating principal quadruplets 
! ---------------------------------
!
! Sorting flexible dihedrals by central bond
!
       allocate(tmpdi%iflexi(4,ndihe),tmpdi%dflexi(ndihe),             &
                tmpdi%fflexi(ndihe))
!
       visited(:) = .FALSE.
! 
       neqquad(:) = 0
       meqquad = 0
!
       k = 0
       do i = 1, dihed%nflexi
         if ( visited(i) ) cycle 
!
         meqquad = meqquad + 1
         neqquad(meqquad) = neqquad(meqquad) + 1
!
         k = k + 1
         tmpdi%iflexi(:,k) = dihed%iflexi(:,i)
         tmpdi%fflexi(k)   = dihed%fflexi(i)
         tmpdi%dflexi(k)   = dihed%dflexi(i)
!
         vaux1(:) = dihed%iflexi(:,i)
         visited(i) = .TRUE.
!
         do j = 1, dihed%nflexi
           if ( visited(j) ) cycle
!
           vaux2(:) = dihed%iflexi(:,j) 
!
           if ( ((vaux1(2).eq.vaux2(2)).and.(vaux1(3).eq.vaux2(3))) .or. &
                ((vaux1(2).eq.vaux2(3)).and.(vaux1(3).eq.vaux2(2))) ) then
             neqquad(meqquad) = neqquad(meqquad) + 1
             visited(j) = .TRUE.
             k = k + 1
             tmpdi%iflexi(:,k) = dihed%iflexi(:,j)
             tmpdi%fflexi(k)   = dihed%fflexi(j)
             tmpdi%dflexi(k)   = dihed%dflexi(j)             
           end if
! 
         end do
       end do
!
       imap(:) = -1
!
       k = 0
       do i = 1, meqquad
!
! Removing H-related dihedrals if possible
!
         nmap = 0
!
         flag = .TRUE.
         do j = 1, neqquad(i)
           vaux(:) = tmpdi%iflexi(:,k+j)
           if ( (znum(vaux(1)).ne.1).and.(znum(vaux(4)).ne.1) ) then 
             imap(nmap+1) = k + j
             nmap = nmap + 1
           end if
         end do        
!
         if ( nmap .eq. 0 ) then
           flag = .TRUE.
           do j = 1, neqquad(i)
             vaux(:) = tmpdi%iflexi(:,k+j)
             if ( (znum(vaux(1)).ne.1).or.(znum(vaux(4)).ne.1) ) then 
               imap(nmap+1) = k + j
               nmap = nmap + 1
             end if
           end do  
         end if
!
         if ( nmap .eq. 0 ) then
           flag = .TRUE.
           do j = 1, neqquad(i)
             vaux(:) = tmpdi%iflexi(:,k+j)
             if ( (znum(vaux(1)).eq.1).and.(znum(vaux(4)).eq.1) ) then 
               imap(nmap+1) = k + j
               nmap = nmap + 1
             end if
           end do  
         end if
!
! Find leading quadruplet
!
         flag = .FALSE.
!
         do j = 1, nmap
           daux = abs(tmpdi%dflexi(imap(j)))
           if ( (daux.ge.0.0d0) .and. (daux.le.35.0d0) ) then
             flag = .TRUE.
             iflexi(i) = imap(j)
             k = k + neqquad(i)
             exit
           end if
         end do
         if ( flag ) cycle
!
         do j = 1, nmap
           daux = abs(tmpdi%dflexi(imap(j)))
           if ( (daux.gt.150.0d0) .and. (daux.lt.181.0d0) ) then
             flag = .TRUE.
             iflexi(i) = imap(j)
             k = k + neqquad(i)
             exit
           end if

         end do
         if ( flag ) cycle
!
         do j = 1, nmap
           daux = abs(tmpdi%dflexi(imap(j)))
           if ( (daux.ge.70.0d0) .and. (daux.le.105.0d0) ) then
             flag = .TRUE.
             iflexi(i) = imap(j)
             k = k + neqquad(i)
             exit
           end if
         end do
         if ( flag ) cycle
!
         do j = 1, nmap
           daux = abs(tmpdi%dflexi(imap(j)))
           if ( (daux.ge.35.0d0) .and. (daux.le.70.0d0) ) then
             flag = .TRUE.
             iflexi(i) = imap(j)
             k = k + neqquad(i)
             exit
           end if
         end do
         if ( flag ) cycle
!
         do j = 1, nmap
           daux = abs(tmpdi%dflexi(imap(j)))
           if ( (daux.ge.105.0d0) .and. (daux.le.150.0d0) ) then
             flag = .TRUE.
             iflexi(i) = imap(j)
             k = k + neqquad(i)
             exit
           end if
         end do
         if ( flag ) cycle
!
         k = k + neqquad(i)
!
         neqquad(i) = nmap
!
       end do
!
! Storing information of selected quadruplets
!
       dihed%nquad = meqquad
       allocate(dihed%iquad(4,meqquad),dihed%dquad(meqquad),           &
                dihed%fquad(meqquad),dihed%mapquad(meqquad))
!
       do i = 1, meqquad
         dihed%iquad(:,i) = tmpdi%iflexi(:,iflexi(i))
         dihed%dquad(i)   = tmpdi%dflexi(iflexi(i))
         dihed%fquad(i)   = tmpdi%fflexi(iflexi(i))
         dihed%mapquad(i) = iflexi(i)
       end do
!
       if ( debug ) then
         write(*,'(1X,A)') 'Principal quadruplets'
         write(*,'(1X,A)') '---------------------'
         do i = 1, dihed%nquad
           write(*,'(1X,A,4I4)') 'Quadruplet',dihed%iquad(:,i)
         end do
         write(*,*)
       end if
!
       return
       end subroutine genquad
!
!======================================================================!
!
    function calc_angle(r1,r2,r3) result(a)

        double precision,dimension(1:3),intent(in)::r1,r2,r3

        real(8) :: a

        !Local
        double precision b1,b2
        double precision,dimension(1:3)::raux1,raux2
        integer :: i

        do i=1,3
            raux1(i)=r1(i)-r2(i)
            raux2(i)=r3(i)-r2(i)
        enddo
        raux1 = r1-r2
        raux2 = r3-r2

        b1=dsqrt(dot_product(raux1,raux1))
        b2=dsqrt(dot_product(raux2,raux2))
        a=dacos(dot_product(raux1,raux2)/b1/b2)

        return

    end function calc_angle
!
!======================================================================!
!
! Subroutine taken from the Joyce code
!
      Subroutine Diedro (ci,cj,ck,cl,val)
!----------------------------------------------------------------------!

!     computes the dihedral betewwn four (in degrees)

      implicit double precision (a-h,o-z)
      dimension ci(3),cj(3),ck(3),cl(3)
      tol=1.d-06
!
      xjk=cj(1)-ck(1)
      yjk=cj(2)-ck(2)
      zjk=cj(3)-ck(3)
      a=dsqrt(xjk**2 + yjk**2 + zjk**2)

      xjk=xjk/a
      yjk=yjk/a
      zjk=zjk/a
!
      xji=cj(1)-ci(1)
      yji=cj(2)-ci(2)
      zji=cj(3)-ci(3)
      scal=xji*xjk+yji*yjk+zji*zjk
      xji=xji-xjk*scal
      yji=yji-yjk*scal
      zji=zji-zjk*scal
      a=dsqrt(xji**2 + yji**2 + zji**2)

      xji=xji/a
      yji=yji/a
      zji=zji/a
!
      xkl=ck(1)-cl(1)
      ykl=ck(2)-cl(2)
      zkl=ck(3)-cl(3)
!~       a=dsqrt(xkl**2 + ykl**2 + zkl**2)
!
      scal=xkl*xjk+ykl*yjk+zkl*zjk
      xkl=xkl-xjk*scal
      ykl=ykl-yjk*scal
      zkl=zkl-zjk*scal
      a=dsqrt(xkl**2 + ykl**2 + zkl**2)

      xkl=xkl/a
      ykl=ykl/a
      zkl=zkl/a
!
      coseno=xji*xkl+yji*ykl+zji*zkl
      seno  =xjk*(yji*zkl-ykl*zji)+yjk*(xkl*zji-xji*zkl)+              &
            zjk*(xji*ykl-yji*xkl)
      val=-atan2(seno,coseno)

      end
!
!======================================================================!
!
       end module genfftools
!
!======================================================================!

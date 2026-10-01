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
       subroutine genrigidlist(nat,r,wiberg,coord,znum,cycles,ncycle,  &
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
       integer,dimension(nat),intent(in)            ::  znum    !
       integer,intent(in)                          ::  nat     !  Number of nodes
       integer,intent(in)                          ::  r       !  Graph rank
       integer,intent(in)                          ::  mcycle  !  Number of cycles
!
! Local variables
! 
       real(kind=8)                                ::  daux    !
       real(kind=8)                                ::  thr     !
       logical                                     ::  lcheck  !
       integer,dimension(4)                        ::  ivaux   !
       integer                                     ::  i,j,k   !
       integer                                     ::  ii      !
!
       real(kind=8),parameter                      ::  pi =  4*atan(1.0_8)
       real(kind=8),parameter                      ::  thr_cc = 1.25d0
       real(kind=8),parameter                      ::  thr_cn = 1.20d0
       real(kind=8),parameter                      ::  thr_nn = 1.20d0
       real(kind=8),parameter                      ::  thr_xx = 1.50d0

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
!  if Wiberg index is greater than an element-pair threshold then it is
!   treated as a double/conjugated rigid bond. Lower thresholds for C-C,
!   C-N and N-N catch partial double bonds such as amides.
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
           else
             if ( (znum(i).eq.6).and.(znum(j).eq.6) ) then
               thr = thr_cc
             else if ( ((znum(i).eq.6).and.(znum(j).eq.7)) .or.       &
                       ((znum(i).eq.7).and.(znum(j).eq.6)) ) then
               thr = thr_cn
             else if ( (znum(i).eq.7).and.(znum(j).eq.7) ) then
               thr = thr_nn
             else
               thr = thr_xx
             end if
!
             if ( wiberg(j,i) .ge. thr ) then
             lrigid(i,j) = .TRUE.
             lrigid(j,i) = .TRUE.
             end if
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
       integer,dimension(ndihe)           ::  qmap     !
       integer,dimension(ndihe)           ::  stmp     !
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
       qmap(:)   = 0
       stmp(:)   = 0
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
         stmp(k)   = dihe%squad(i)
         qmap(i)   = k
         ntmp = k      
!
         do j = i+1, dihe%nquad
           if ( visited(j) ) cycle
           if ( iddihe(i) .eq. iddihe(j) ) then
             if ( (dihe%squad(i).eq.1).and.(dihe%squad(j).eq.1) ) cycle
             visited(j) = .TRUE.
             qmap(j) = k
             exit
           end if
         end do
       end do
!
       dihe%iquad(:,:) = itmp(:,:)
       dihe%dquad(:)   = dtmp(:)
       dihe%fquad(:)   = ftmp(:)
       dihe%squad(:)   = stmp(:)
       if ( allocated(dihe%depquad) ) then
         do i = 1, dihe%nflexi
           if ( dihe%depquad(i) .gt. 0 ) then
             if ( qmap(dihe%depquad(i)) .gt. 0 ) then
               dihe%depquad(i) = qmap(dihe%depquad(i))
             else
               dihe%depquad(i) = 0
             end if
           end if
         end do
       end if
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
       logical                                 ::  ldep     !
       integer                                 ::  i,j,k    !
! 
!  Including only principal quadruplets 
! -------------------------------------
!
! Generating new representation with a reduced set of quadruplets
!
       allocate(tmpdihe%iflexi(4,ndihe),tmpdihe%dflexi(ndihe),         &
                tmpdihe%kflexi(ndihe),tmpdihe%fflexi(ndihe),           &
                tmpdihe%flexi(ndihe),tmpdihe%depquad(ndihe))
!
       tmpdihe%iflexi(:,:) = -1
       tmpdihe%dflexi(:)   = 999.0d0
       tmpdihe%kflexi(:)   = 0.0d0
       tmpdihe%fflexi(:)   = -1
       tmpdihe%depquad(:)  = 0
!
       k = 0
       do i = 1, dihe%nflexi
         ivaux1(:) = dihe%iflexi(:,i)
         ldep = allocated(dihe%depquad)
         if ( ldep ) ldep = dihe%depquad(i) .gt. 0
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
               if ( allocated(dihe%depquad) )                          &
                 tmpdihe%depquad(k) = dihe%depquad(i)
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
           if ( ldep ) then
!
             k = k + 1
!
             tmpdihe%iflexi(:,k) = dihe%iflexi(:,i)
             tmpdihe%dflexi(k)   = dihe%dflexi(i)
             tmpdihe%kflexi(k)   = dihe%kflexi(i)
             tmpdihe%fflexi(k)   = dihe%fflexi(i)
             tmpdihe%depquad(k)  = dihe%depquad(i)
!
             tmpdihe%flexi(k)%ntor = dihe%flexi(i)%ntor
             allocate(tmpdihe%flexi(k)%tor(tmpdihe%flexi(k)%ntor))
!
             tmpdihe%flexi(k)%itor(:) = dihe%flexi(i)%itor(:)
             tmpdihe%flexi(k)%tor     = dihe%flexi(i)%tor
!
             cycle
!
           end if
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
               if ( allocated(dihe%depquad) )                          &
                 tmpdihe%depquad(k) = dihe%depquad(i)
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
       if ( allocated(dihe%depquad) ) dihe%depquad(:) = tmpdihe%depquad(:)
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
       subroutine genquad(nat,coord,adj,ideg,lcycle,lrigid,znum,dihed, &
                          ndihe,debug)
!
       use datatypes,  only:  dihedrals
!
       implicit none
!
! Input/output variables
!
       type(dihedrals),intent(inout)             ::  dihed    !
       real(kind=8),dimension(3,nat),intent(in)  ::  coord    !
       logical,dimension(nat,nat),intent(in)     ::  adj      !
       logical,dimension(nat,nat),intent(in)     ::  lcycle   !
       logical,dimension(nat,nat),intent(in)     ::  lrigid   !
       integer,dimension(nat),intent(in)         ::  ideg     !
       integer,dimension(nat),intent(in)         ::  znum     !
       integer,intent(in)                        ::  nat      !
       integer,intent(in)                        ::  ndihe    !
       logical,intent(in)                        ::  debug    !
!
! Local variables
!
       type(dihedrals)                    ::  tmpdi     !
       logical,dimension(ndihe)           ::  visited   !
       logical                            ::  flag      !
       integer,dimension(ndihe)           ::  imap      !
       integer,dimension(ndihe)           ::  iadd      !
       integer,dimension(ndihe)           ::  iflexi    !
       integer,dimension(ndihe)           ::  iselect   !
       integer,dimension(ndihe)           ::  sselect   !
       integer,dimension(ndihe)           ::  tmpmap    !
       integer,dimension(ndihe)           ::  neqquad   !
       integer,dimension(4)               ::  vaux      !
       integer,dimension(4)               ::  vaux1     !
       integer,dimension(4)               ::  vaux2     !
       integer                            ::  meqquad   !
       integer                            ::  nmap      !
       integer                            ::  nadd      !
       integer                            ::  nselect   !
       integer                            ::  i,j,k     !
!
       real(kind=8)                       ::  daux      !
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
         tmpmap(k)         = i
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
             tmpmap(k)         = j
           end if
! 
         end do
       end do
!
       imap(:) = -1
       iflexi(:)  = -1
       iselect(:) = -1
       sselect(:) = 0
       nselect = 0
!
       if ( allocated(dihed%depquad) ) deallocate(dihed%depquad)
       allocate(dihed%depquad(ndihe))
       dihed%depquad(:) = 0
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
! Select multiple ring-exocyclic quadruplets when the exocyclic atom is
! geometrically compatible with the requested cis/cis-trans cases.
!
         call select_ring_exocyclic(nat,coord,adj,ideg,lcycle,lrigid, &
                                    tmpdi,k,neqquad(i),ndihe,imap,    &
                                    nmap,iadd,nadd)
!
         if ( nadd .gt. 0 ) then
           nselect = nselect + 1
!
           iselect(nselect) = imap(1)
!
           do j = 1, nmap
             daux = abs(tmpdi%dflexi(imap(j)))
             if ( (daux.ge.0.0d0) .and. (daux.le.35.0d0) ) then
               iselect(nselect) = imap(j)
               GOTO 1000
             end if
           end do
!
           do j = 1, nmap
             daux = abs(tmpdi%dflexi(imap(j)))
             if ( (daux.gt.150.0d0) .and. (daux.lt.181.0d0) ) then
               iselect(nselect) = imap(j)
               GOTO 1000
             end if
           end do
!
           do j = 1, nmap
             daux = abs(tmpdi%dflexi(imap(j)))
             if ( (daux.ge.70.0d0) .and. (daux.le.105.0d0) ) then
               iselect(nselect) = imap(j)
               GOTO 1000
             end if
           end do
!
           do j = 1, nmap
             daux = abs(tmpdi%dflexi(imap(j)))
             if ( (daux.ge.35.0d0) .and. (daux.le.70.0d0) ) then
               iselect(nselect) = imap(j)
               GOTO 1000
             end if
           end do
!
           do j = 1, nmap
             daux = abs(tmpdi%dflexi(imap(j)))
             if ( (daux.ge.105.0d0) .and. (daux.le.150.0d0) ) then
               iselect(nselect) = imap(j)
               GOTO 1000
             end if
           end do
!
1000       continue
!
           sselect(nselect) = 0
           do j = 1, nadd
             if ( iadd(j) .ne. iselect(nselect) ) then
               dihed%depquad(tmpmap(iadd(j))) = nselect
             end if
           end do
           k = k + neqquad(i)
           cycle
         end if
!
! Find leading quadruplet
!
         iflexi(i) = leading_quad(tmpdi,imap,nmap,ndihe)
         nselect = nselect + 1
         iselect(nselect) = iflexi(i)
         sselect(nselect) = 0
!
         k = k + neqquad(i)
!
       end do
!
! Storing information of selected quadruplets
!
       dihed%nquad = nselect
       allocate(dihed%iquad(4,nselect),dihed%dquad(nselect),           &
                dihed%fquad(nselect),dihed%mapquad(nselect),           &
                dihed%squad(nselect))
!
       do i = 1, nselect
         dihed%iquad(:,i) = tmpdi%iflexi(:,iselect(i))
         dihed%dquad(i)   = tmpdi%dflexi(iselect(i))
         dihed%fquad(i)   = tmpdi%fflexi(iselect(i))
         dihed%mapquad(i) = iselect(i)
         dihed%squad(i)   = sselect(i)
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
!
       end subroutine genquad
!
!======================================================================!
!
       subroutine select_ring_exocyclic(nat,coord,adj,ideg,lcycle,    &
                                        lrigid,tmpdi,first,nitem,ndihe,&
                                        imap,nmap,iout,nout)

       use datatypes, only: dihedrals
!
       implicit none
!
       type(dihedrals),intent(in)             ::  tmpdi   !
       real(kind=8),dimension(3,nat),intent(in) ::  coord !
       logical,dimension(nat,nat),intent(in)  ::  adj     !
       logical,dimension(nat,nat),intent(in)  ::  lcycle  !
       logical,dimension(nat,nat),intent(in)  ::  lrigid  !
       integer,dimension(nat),intent(in)      ::  ideg    !
       integer,intent(in)                    ::  nat    !
       integer,intent(in)                    ::  ndihe  !
       integer,dimension(ndihe),intent(in)   ::  imap   !
       integer,dimension(ndihe),intent(out)  ::  iout   !
       integer,intent(in)                    ::  first  !
       integer,intent(in)                    ::  nitem  !
       integer,intent(in)                    ::  nmap   !
       integer,intent(out)                   ::  nout   !
!
       real(kind=8)                          ::  daux   !
       logical                               ::  ring1  !
       logical                               ::  ring2  !
       integer                               ::  c1     !
       integer                               ::  c2     !
       integer                               ::  exo    !
       integer                               ::  ring   !
       integer                               ::  ncis   !
       integer                               ::  ntrans !
       integer                               ::  j      !
!
       iout(:) = -1
       nout = 0
!
       if ( nmap .le. 1 ) return
!
       c1 = tmpdi%iflexi(2,first+1)
       c2 = tmpdi%iflexi(3,first+1)
!
       ring1 = has_ring_terminal(tmpdi,lcycle,c1,c2,first,nitem)
       ring2 = has_ring_terminal(tmpdi,lcycle,c2,c1,first,nitem)
!
       if ( ring1 .eqv. ring2 ) return
!
       if ( ring1 ) then
         ring = c1
         exo  = c2
       else
         ring = c2
         exo  = c1
       end if
!
       if ( is_sp2_like_exocyclic(nat,coord,adj,ideg,lrigid,exo,ring) ) then
         do j = 1, nmap
           daux = abs(tmpdi%dflexi(imap(j)))
           if ( (daux.ge.0.0d0) .and. (daux.le.35.0d0) ) then
             nout = nout + 1
             iout(nout) = imap(j)
           end if
         end do
         if ( nout .eq. 2 ) return
         nout = 0
       end if
!
       if ( is_valid_degree2_exocyclic(nat,coord,adj,ideg,lrigid,exo,ring) ) then
         ncis   = 0
         ntrans = 0
         do j = 1, nmap
           daux = abs(tmpdi%dflexi(imap(j)))
           if ( (daux.ge.0.0d0) .and. (daux.le.35.0d0) ) then
             ncis = ncis + 1
           else if ( (daux.gt.150.0d0) .and. (daux.lt.181.0d0) ) then
             ntrans = ntrans + 1
           end if
           nout = nout + 1
           iout(nout) = imap(j)
         end do
         if ( (nout.eq.2).and.(ncis.eq.1).and.(ntrans.eq.1) ) return
         nout = 0
       end if
!
       return
       end subroutine select_ring_exocyclic
!
!======================================================================!
!
       logical function has_ring_terminal(tmpdi,lcycle,center,other,  &
                                          first,nitem)

       use datatypes, only: dihedrals
!
       implicit none
!
       type(dihedrals),intent(in)          ::  tmpdi   !
       logical,dimension(:,:),intent(in)   ::  lcycle  !
       integer,intent(in)          ::  center  !
       integer,intent(in)          ::  other   !
       integer,intent(in)          ::  first   !
       integer,intent(in)          ::  nitem   !
!
       logical                     ::  match   !
       integer,dimension(4)        ::  iquad   !
       integer                     ::  tcenter !
       integer                     ::  tother  !
       integer                     ::  j       !
!
       has_ring_terminal = .FALSE.
!
       do j = 1, nitem
         iquad(:) = tmpdi%iflexi(:,first+j)
         call get_terminals(iquad,center,other,tcenter,tother,match)
         if ( match .and. lcycle(center,tcenter) ) then
           has_ring_terminal = .TRUE.
           return
         end if
       end do
!
       return
       end function has_ring_terminal
!
!======================================================================!
!
       subroutine get_terminals(iquad,center,other,tcenter,tother,match)
!
       implicit none
!
       logical,intent(out)               ::  match   !
       integer,dimension(4),intent(in)   ::  iquad   !
       integer,intent(in)                ::  center  !
       integer,intent(in)                ::  other   !
       integer,intent(out)               ::  tcenter !
       integer,intent(out)               ::  tother  !
!
       match = .FALSE.
       tcenter = -1
       tother  = -1
!
       if ( (iquad(2).eq.center).and.(iquad(3).eq.other) ) then
         tcenter = iquad(1)
         tother  = iquad(4)
         match = .TRUE.
       else if ( (iquad(3).eq.center).and.(iquad(2).eq.other) ) then
         tcenter = iquad(4)
         tother  = iquad(1)
         match = .TRUE.
       end if
!
       return
       end subroutine get_terminals
!
!======================================================================!
!
       logical function is_sp2_like_exocyclic(nat,coord,adj,ideg,     &
                                               lrigid,exo,ring)
!
       implicit none
!
       real(kind=8),dimension(3,nat),intent(in)  ::  coord   !
       logical,dimension(nat,nat),intent(in)     ::  adj     !
       logical,dimension(nat,nat),intent(in)     ::  lrigid  !
       integer,dimension(nat),intent(in)         ::  ideg    !
       integer,intent(in)                        ::  nat     !
       integer,intent(in)               ::  exo   !
       integer,intent(in)               ::  ring  !
!
       is_sp2_like_exocyclic = .FALSE.
!
       if ( ideg(exo) .ne. 3 ) return
!
       is_sp2_like_exocyclic = is_planar_degree3(nat,coord,adj,exo)   &
                               .or. has_rigid_exo_bond(nat,adj,lrigid, &
                                                       exo,ring)
!
       return
       end function is_sp2_like_exocyclic
!
!======================================================================!
!
       logical function is_planar_degree3(nat,coord,adj,idx)
!
       implicit none
!
       real(kind=8),dimension(3,nat),intent(in) ::  coord  !
       logical,dimension(nat,nat),intent(in)    ::  adj    !
       integer,intent(in)                       ::  nat    !
       integer,intent(in)                 ::  idx    !
!
       real(kind=8),dimension(3)          ::  angle  !
       real(kind=8)                       ::  asum   !
       integer,dimension(3)               ::  nei    !
       integer                            ::  n      !
       integer                            ::  j      !
!
       real(kind=8),parameter             ::  pi = 4*atan(1.0_8)
!
       is_planar_degree3 = .FALSE.
!
       n = 0
       do j = 1, nat
         if ( adj(idx,j) ) then
           n = n + 1
           if ( n .le. 3 ) nei(n) = j
         end if
       end do
!
       if ( n .ne. 3 ) return
!
       angle(1) = calc_angle(coord(:,nei(1)),coord(:,idx),coord(:,nei(2))) &
                  * 180.0d0 / pi
       angle(2) = calc_angle(coord(:,nei(1)),coord(:,idx),coord(:,nei(3))) &
                  * 180.0d0 / pi
       angle(3) = calc_angle(coord(:,nei(2)),coord(:,idx),coord(:,nei(3))) &
                  * 180.0d0 / pi
!
       asum = angle(1) + angle(2) + angle(3)
!
       is_planar_degree3 = (abs(asum-360.0d0).le.25.0d0)               &
                          .and. (maxval(angle).lt.170.0d0)
!
       return
       end function is_planar_degree3
!
!======================================================================!
!
       logical function has_rigid_exo_bond(nat,adj,lrigid,exo,ring)
!
       implicit none
!
       logical,dimension(nat,nat),intent(in) ::  adj     !
       logical,dimension(nat,nat),intent(in) ::  lrigid  !
       integer,intent(in)                    ::  nat     !
       integer,intent(in)        ::  exo   !
       integer,intent(in)        ::  ring  !
       integer                   ::  j     !
!
       has_rigid_exo_bond = .FALSE.
!
       do j = 1, nat
         if ( j .eq. ring ) cycle
         if ( adj(exo,j) .and. lrigid(exo,j) ) then
           has_rigid_exo_bond = .TRUE.
           return
         end if
       end do
!
       return
       end function has_rigid_exo_bond
!
!======================================================================!
!
       logical function is_valid_degree2_exocyclic(nat,coord,adj,ideg,&
                                                    lrigid,exo,ring)
!
       implicit none
!
       real(kind=8),dimension(3,nat),intent(in) ::  coord   !
       logical,dimension(nat,nat),intent(in)    ::  adj     !
       logical,dimension(nat,nat),intent(in)    ::  lrigid  !
       integer,dimension(nat),intent(in)        ::  ideg    !
       integer,intent(in)                       ::  nat     !
       integer,intent(in)             ::  exo    !
       integer,intent(in)             ::  ring   !
!
       real(kind=8)                   ::  angle  !
       integer                        ::  other  !
!
       real(kind=8),parameter         ::  pi = 4*atan(1.0_8)
!
       is_valid_degree2_exocyclic = .FALSE.
!
       if ( ideg(exo) .ne. 2 ) return
!
       other = other_exocyclic_neighbor(nat,adj,exo,ring)
       if ( other .lt. 1 ) return
!
       angle = calc_angle(coord(:,ring),coord(:,exo),coord(:,other))    &
               * 180.0d0 / pi
!
       is_valid_degree2_exocyclic = lrigid(exo,other)                  &
                                    .or. (angle.ge.160.0d0)
!
       return
       end function is_valid_degree2_exocyclic
!
!======================================================================!
!
       integer function other_exocyclic_neighbor(nat,adj,exo,ring)
!
       implicit none
!
       logical,dimension(nat,nat),intent(in) ::  adj    !
       integer,intent(in)                    ::  nat    !
       integer,intent(in)        ::  exo   !
       integer,intent(in)        ::  ring  !
       integer                   ::  j     !
!
       other_exocyclic_neighbor = -1
!
       do j = 1, nat
         if ( j .eq. ring ) cycle
         if ( adj(exo,j) ) then
           other_exocyclic_neighbor = j
           return
         end if
       end do
!
       return
       end function other_exocyclic_neighbor
!
!======================================================================!
!
       integer function leading_quad(tmpdi,imap,nmap,ndihe)

       use datatypes, only: dihedrals

       implicit none

       type(dihedrals),intent(in)           ::  tmpdi  !
       integer,intent(in)                   ::  nmap   !
       integer,intent(in)                   ::  ndihe  !
       integer,dimension(ndihe),intent(in)  ::  imap   !

       real(kind=8)                         ::  daux   !
       integer                              ::  j      !

       leading_quad = imap(1)

       do j = 1, nmap
         daux = abs(tmpdi%dflexi(imap(j)))
         if ( (daux.ge.0.0d0) .and. (daux.le.35.0d0) ) then
           leading_quad = imap(j)
           return
         end if
       end do

       do j = 1, nmap
         daux = abs(tmpdi%dflexi(imap(j)))
         if ( (daux.gt.150.0d0) .and. (daux.lt.181.0d0) ) then
           leading_quad = imap(j)
           return
         end if
       end do

       do j = 1, nmap
         daux = abs(tmpdi%dflexi(imap(j)))
         if ( (daux.ge.70.0d0) .and. (daux.le.105.0d0) ) then
           leading_quad = imap(j)
           return
         end if
       end do

       do j = 1, nmap
         daux = abs(tmpdi%dflexi(imap(j)))
         if ( (daux.ge.35.0d0) .and. (daux.le.70.0d0) ) then
           leading_quad = imap(j)
           return
         end if
       end do

       do j = 1, nmap
         daux = abs(tmpdi%dflexi(imap(j)))
         if ( (daux.ge.105.0d0) .and. (daux.le.150.0d0) ) then
           leading_quad = imap(j)
           return
         end if
       end do

       return
       end function leading_quad
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

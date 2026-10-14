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
                   findarcycles,                                       &
                   bonded2dihe,                                        &
                   apply_ref_flexible_terms,                           &
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
       dihe%iimpro(:,:) = 0
       dihe%dimpro(:)   = 0.0d0
       dihe%kimpro(:)   = 0.0d0
       dihe%fimpro(:)   = 0
       dihe%iinv(:,:)   = 0
       dihe%dinv(:)     = 0.0d0
       dihe%kinv(:)     = 0.0d0
       dihe%finv(:)     = 0
       dihe%irigid(:,:) = 0
       dihe%drigid(:)   = 0.0d0
       dihe%krigid(:)   = 0.0d0
       dihe%frigid(:)   = 0
       dihe%iflexi(:,:) = 0
       dihe%dflexi(:)   = 0.0d0
       dihe%kflexi(:)   = 0.0d0
       dihe%fflexi(:)   = 0
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
           dihe%fflexi(dihe%nflexi) = 1
!
           dihe%flexi(dihe%nflexi)%itor(:) = bonded%idihe(:,i)
           dihe%iflexi(:,dihe%nflexi)      = bonded%idihe(:,i)
           dihe%dflexi(dihe%nflexi)        = bonded%dihe(i)
           dihe%kflexi(dihe%nflexi)        = bonded%kdihe(i)
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
           dihe%dflexi(dihe%nflexi)   = bonded%dihe(i)
           dihe%kflexi(dihe%nflexi)   = bonded%kdihe(i)
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
! APPLY_REF_FLEXIBLE_TERMS - use reference topology flexible terms
!
! This subroutine preserves the Fourier multiplicities explicitly present in a
! reference topology.  It is useful in Step 2, where bonded terms are generated
! again from QM data but the user may have removed some soft torsional terms
! after Step 1.
!
       subroutine apply_ref_flexible_terms(dihe,refbonded)
!
       use datatypes,  only:  grobonded,                               &
                              dihedrals,                               &
                              torsion
!
       implicit none
!
! Input/output variables
!
       type(dihedrals),intent(inout)       ::  dihe      !
       type(grobonded),intent(in)          ::  refbonded !
!
! Local variables
!
       logical,dimension(:),allocatable    ::  matched  !
       integer                             ::  i,j,k     !
       integer                             ::  nmatch    !
       integer                             ::  source    !
!
! Use only explicitly matching proper torsional terms from the reference
! topology.  The matching is done by atom quadruplet and function type; atom
! order can be reversed.
!
       if ( dihe%nflexi .le. 0 ) return
!
       allocate(matched(dihe%nflexi))
       matched(:) = .FALSE.
!
       do i = 1, dihe%nflexi
!
         nmatch = 0
         do j = 1, refbonded%ndihe
           if ( refbonded%fdihe(j) .ne. dihe%fflexi(i) ) cycle
           if ( same_flexible_quad(refbonded%idihe(:,j),               &
                                   dihe%iflexi(:,i)) ) then
             nmatch = nmatch + 1
           end if
         end do
!
         if ( nmatch .le. 0 ) cycle
!
         matched(i) = .TRUE.
!
         if ( allocated(dihe%flexi(i)%tor) ) deallocate(dihe%flexi(i)%tor)
         dihe%flexi(i)%ntor = nmatch
         allocate(dihe%flexi(i)%tor(nmatch))
         dihe%flexi(i)%itor(:) = dihe%iflexi(:,i)
!
         k = 0
         do j = 1, refbonded%ndihe
           if ( refbonded%fdihe(j) .ne. dihe%fflexi(i) ) cycle
           if ( same_flexible_quad(refbonded%idihe(:,j),               &
                                   dihe%iflexi(:,i)) ) then
             k = k + 1
             dihe%flexi(i)%tor(k)%phase = refbonded%dihe(j)
             dihe%flexi(i)%tor(k)%vtor  = refbonded%kdihe(j)
             dihe%flexi(i)%tor(k)%multi = refbonded%multi(j)
           end if
         end do
!
       end do
!
! If a topology contains only the independent member of a dependent flexible
! torsion group, keep the dependent member synchronized with the same
! multiplicities.  This preserves user-edited Fourier truncations between
! Step 1 and Step 2 while keeping the dependent term's own phase.
!
       if ( allocated(dihe%depquad) ) then
         do i = 1, dihe%nflexi
           if ( matched(i) ) cycle
!
           source = 0
           if ( dihe%depquad(i) .gt. 0 ) then
             source = find_quad_flex_term(dihe,                       &
                                          dihe%iquad(:,dihe%depquad(i)))
             if ( source .gt. 0 ) then
               if ( .NOT. matched(source) ) source = 0
             end if
           else
             do j = 1, dihe%nflexi
               if ( .NOT. matched(j) ) cycle
               if ( dihe%depquad(j) .le. 0 ) cycle
               if ( same_flexible_quad(dihe%iflexi(:,i),              &
                    dihe%iquad(:,dihe%depquad(j))) ) then
                 source = j
                 exit
               end if
             end do
           end if
!
           if ( source .gt. 0 ) call mirror_reference_terms(dihe,i,source)
!
         end do
       end if
!
       deallocate(matched)
!
       return
!
       contains
!
       logical function same_flexible_quad(ia,ib)
!
       integer,dimension(4),intent(in)     ::  ia !
       integer,dimension(4),intent(in)     ::  ib !
!
       same_flexible_quad =                                          &
         ( (ia(1).eq.ib(1)).and.(ia(2).eq.ib(2)).and.                 &
           (ia(3).eq.ib(3)).and.(ia(4).eq.ib(4)) ) .or.               &
         ( (ia(1).eq.ib(4)).and.(ia(2).eq.ib(3)).and.                 &
           (ia(3).eq.ib(2)).and.(ia(4).eq.ib(1)) )
!
       return
       end function same_flexible_quad
!
!----------------------------------------------------------------------!
!
       integer function find_quad_flex_term(dihe,iquad)
!
       type(dihedrals),intent(in)          ::  dihe  !
       integer,dimension(4),intent(in)     ::  iquad !
       integer                             ::  k     !
!
       find_quad_flex_term = 0
       do k = 1, dihe%nflexi
         if ( same_flexible_quad(dihe%iflexi(:,k),iquad) ) then
           find_quad_flex_term = k
           return
         end if
       end do
!
       return
       end function find_quad_flex_term
!
!----------------------------------------------------------------------!
!
       subroutine mirror_reference_terms(dihe,target,source)
!
       type(dihedrals),intent(inout)       ::  dihe    !
       integer,intent(in)                  ::  target  !
       integer,intent(in)                  ::  source  !
!
       type(torsion),dimension(:),allocatable :: oldtor !
       integer                             ::  old_n   !
       integer                             ::  k       !
!
       if ( source .eq. target ) return
       if ( dihe%flexi(source)%ntor .le. 0 ) return
!
       old_n = dihe%flexi(target)%ntor
       if ( old_n .gt. 0 ) then
         allocate(oldtor(old_n))
         oldtor = dihe%flexi(target)%tor
       end if
!
       if ( allocated(dihe%flexi(target)%tor) )                        &
         deallocate(dihe%flexi(target)%tor)
!
       dihe%flexi(target)%ntor = dihe%flexi(source)%ntor
       allocate(dihe%flexi(target)%tor(dihe%flexi(target)%ntor))
       dihe%flexi(target)%itor(:) = dihe%iflexi(:,target)
!
       do k = 1, dihe%flexi(target)%ntor
         dihe%flexi(target)%tor(k)%multi =                            &
                                  dihe%flexi(source)%tor(k)%multi
         dihe%flexi(target)%tor(k)%phase = term_phase(oldtor,old_n,   &
                                  dihe%flexi(source)%tor(k)%multi,    &
                                  dihe%dflexi(target))
         dihe%flexi(target)%tor(k)%vtor = term_vtor(oldtor,old_n,     &
                                  dihe%flexi(source)%tor(k)%multi)
       end do
!
       if ( allocated(oldtor) ) deallocate(oldtor)
!
       return
       end subroutine mirror_reference_terms
!
!----------------------------------------------------------------------!
!
       real(kind=8) function term_phase(oldtor,old_n,multi,default)
!
       type(torsion),dimension(:),allocatable,intent(in) :: oldtor !
       integer,intent(in)                  ::  old_n   !
       integer,intent(in)                  ::  multi   !
       real(kind=8),intent(in)             ::  default !
       integer                             ::  k       !
!
       term_phase = default
       if ( old_n .le. 0 ) return
       do k = 1, old_n
         if ( oldtor(k)%multi .eq. multi ) then
           term_phase = oldtor(k)%phase
           return
         end if
       end do
!
       return
       end function term_phase
!
!----------------------------------------------------------------------!
!
       real(kind=8) function term_vtor(oldtor,old_n,multi)
!
       type(torsion),dimension(:),allocatable,intent(in) :: oldtor !
       integer,intent(in)                  ::  old_n   !
       integer,intent(in)                  ::  multi   !
       integer                             ::  k       !
!
       term_vtor = 0.0d0
       if ( old_n .le. 0 ) return
       do k = 1, old_n
         if ( oldtor(k)%multi .eq. multi ) then
           term_vtor = oldtor(k)%vtor
           return
         end if
       end do
!
       return
       end function term_vtor
!
       end subroutine apply_ref_flexible_terms
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

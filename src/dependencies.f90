!======================================================================!
!
       module dependencies
!
       use units,  only:  unitmp,unideps,uniscr
!
       implicit none
!
       private
       public  ::  symffbonded
!
       contains
!
!======================================================================!
!
! SYMFFBONDED - SYMmetrize Force Field BONDED interactions
!
! This subroutine
!
       subroutine symffbonded(nat,nidat,idat,newlab,bonded,dihed,      &
                              r,marunit,narunit,arunit,lheavy,adj,     &
                              fsymm,fadhoc,debug)
!
       use datatypes,   only: grobonded,                               &
                              dihedrals
       use lengths,     only: lenlab
       use graphtools,  only: inormlabels
!
       implicit none
!
! Input/output variables
!
       type(grobonded),intent(inout)                    ::  bonded    !
       type(dihedrals),intent(inout)                    ::  dihed     !
       logical,dimension(nat,nat),intent(in)            ::  adj       !
       logical,dimension(nat),intent(in)                ::  lheavy    !
       integer,dimension(r,nat),intent(in)              ::  arunit    !
       integer,dimension(nat),intent(in)                ::  idat      !
       integer,dimension(r),intent(in)                  ::  narunit   !
       character(len=lenlab),dimension(nat),intent(in)  ::  newlab    !  New atom type name
!
       integer,intent(in)                               ::  nat       !
       integer,intent(in)                               ::  marunit   !
       integer,intent(in)                               ::  r         !
       integer,intent(in)                               ::  nidat     !
!
       logical,intent(in)                               ::  fsymm     !
       logical,intent(in)                               ::  fadhoc !
       logical,intent(in)                               ::  debug     !  Debug mode
!
! Local variables
!
       character(len=20)                                ::  cnum      !
       integer,dimension(:),allocatable                 ::  ivaux     !
       integer,dimension(4)                             ::  vaux      !
       integer                                          ::  iaux      !
       integer                                          ::  iterm     !
       integer                                          ::  i,j       !
       integer                                          ::  orderidx !
!
! Symmetrizing bonded terms ! TODO: option to reorder dihedral atoms
! -------------------------
!
!  Finding equivalent bond-stretching terms
!
       allocate(bonded%sbond(bonded%nbond),                            &
                bonded%idbond(bonded%nbond),                           &
                bonded%labbond(bonded%nbond))
!
       call labbond(nat,idat,nidat,bonded%nbond,bonded%ibond,          &
                    bonded%idbond,bonded%nidbond,debug)
!
       iterm = 0
       do i = 1, bonded%nbond
         iterm = iterm + 1
         bonded%sbond(i) = iterm
         bonded%labbond(i) = trim(newlab(bonded%ibond(1,i)))//'-'//    &
                             trim(newlab(bonded%ibond(2,i)))
       end do
!
       call symterm(bonded%nbond,bonded%idbond,bonded%bond,            &
                    bonded%sbond,bonded%labbond,fsymm)
!
!  Finding equivalent angle-bending terms
!
       if ( bonded%nang .eq. 0 ) return
!
       allocate(bonded%sang(bonded%nang),                              &
                bonded%idang(bonded%nang),                             &
                bonded%labang(bonded%nang))
!
       call labangle(nat,idat,nidat,bonded%nang,bonded%iang,           &
                     bonded%idang,bonded%nidang,debug)
!
       do i = 1, bonded%nang
         iterm = iterm + 1
         bonded%sang(i) = iterm
         bonded%labang(i) = trim(newlab(bonded%iang(1,i)))//'-'//      &
                            trim(newlab(bonded%iang(2,i)))//'-'//      &
                            trim(newlab(bonded%iang(3,i)))
       end do
!
       call symterm(bonded%nang,bonded%idang,bonded%ang,bonded%sang,   &
                    bonded%labang,fsymm)
!
!  Finding equivalent dihedrals terms
!
       if ( dihed%ndihe .eq. 0 ) return
!
       allocate(ivaux(dihed%ndihe))
!
       allocate(dihed%srigid(dihed%ndihe),                             &
                dihed%idrigid(dihed%ndihe),                            &
                dihed%labrigid(dihed%ndihe))
!
       do i = 1, dihed%nrigid
         iterm = iterm + 1
         vaux(:) = dihed%irigid(:,i)
!
         iaux = idat(vaux(1))*nidat**3 + idat(vaux(2))*nidat**2        &
                                   + idat(vaux(3))*nidat + idat(vaux(4))
!
!~ write(*,*) 'RIGID',vaux(:),';',idat(vaux(1)),idat(vaux(2)),idat(vaux(3)),idat(vaux(4))
!
         dihed%srigid(i)  = iterm
         dihed%idrigid(i) = iaux
!
         dihed%labrigid(i) = trim(newlab(dihed%irigid(1,i)))//'-'//    &
                             trim(newlab(dihed%irigid(2,i)))//'-'//    &
                             trim(newlab(dihed%irigid(3,i)))//'-'//    &
                             trim(newlab(dihed%irigid(4,i)))
       end do
!
       allocate(dihed%simpro(dihed%ndihe),                             &
                dihed%idimpro(dihed%ndihe),                            &
                dihed%labimpro(dihed%ndihe))
!
! TODO: (option) impropers with same central atom are equivalent
!

!
       do i = 1, dihed%nimpro
         iterm = iterm + 1
         vaux(:) = dihed%iimpro(:,i)
!
         iaux = idat(vaux(1))*nidat**3 + idat(vaux(2))*nidat**2        &
                                   + idat(vaux(3))*nidat + idat(vaux(4))
!
         dihed%simpro(i)  = iterm
         dihed%idimpro(i) = iaux
!
         dihed%labimpro(i) = trim(newlab(dihed%iimpro(1,i)))//'-'//    &
                             trim(newlab(dihed%iimpro(2,i)))//'···'//  &
                             trim(newlab(dihed%iimpro(3,i)))//'···'//  &
                             trim(newlab(dihed%iimpro(4,i)))
       end do
!
! TODO: (option) inversion dihedrals with same central atom are equivalent
!

!
       allocate(dihed%sinv(dihed%ndihe),                               &
                dihed%idinv(dihed%ndihe),                              &
                dihed%labinv(dihed%ndihe))
!
       do i = 1, dihed%ninv
         iterm = iterm + 1
         vaux(:) = dihed%iinv(:,i)
!
         iaux = idat(vaux(1))*nidat**3 + idat(vaux(2))*nidat**2        &
                                   + idat(vaux(3))*nidat + idat(vaux(4))
!
         dihed%sinv(i)  = iterm
         dihed%idinv(i) = iaux
!
         dihed%labinv(i) = trim(newlab(dihed%iinv(1,i)))//'-'//        &
                           trim(newlab(dihed%iinv(2,i)))//'···'//      &
                           trim(newlab(dihed%iinv(3,i)))//'···'//      &
                           trim(newlab(dihed%iinv(4,i)))
       end do
!
       allocate(dihed%sflexi(dihed%ndihe),                             &
                dihed%idflexi(dihed%ndihe),                            &
                dihed%labflexi(dihed%ndihe))
!
       do i = 1, dihed%nflexi
         vaux(:) = dihed%iflexi(:,i)
         iaux = min(idat(vaux(2)),idat(vaux(3)))*nidat                 &
              + max(idat(vaux(2)),idat(vaux(3)))
!
         dihed%idflexi(i) = iaux
!
!
         dihed%labflexi(i) = trim(newlab(dihed%iflexi(1,i)))//'-'//    &
                             trim(newlab(dihed%iflexi(2,i)))//'-'//    &
                             trim(newlab(dihed%iflexi(3,i)))//'-'//    &
                             trim(newlab(dihed%iflexi(4,i)))
       end do
!
       dihed%ntor = 0
       do orderidx = 1, dihed%nflexi
         i=orderidx
         if (allocated(dihed%flex_order)) i=dihed%flex_order(orderidx)
         do j = 1, dihed%flexi(i)%ntor
           iterm = iterm + 1
           dihed%flexi(i)%tor(j)%stor = iterm
!
           dihed%ntor = dihed%ntor + 1
!
           write(cnum,*) dihed%flexi(i)%tor(j)%multi
           cnum = adjustl(cnum)
           dihed%flexi(i)%tor(j)%labtor = trim(dihed%labflexi(i))//    &
                                                       '_n='//trim(cnum)
         end do
       end do
!
       ivaux(:) = dihed%idrigid(:)
       call inormlabels(dihed%nrigid,ivaux(:dihed%nrigid),             &
                        dihed%idrigid(:dihed%nrigid),iaux)
!
       ivaux(:) = dihed%idimpro(:)
       call inormlabels(dihed%nimpro,ivaux(:dihed%nimpro),             &
                        dihed%idimpro(:dihed%nimpro),iaux)
!
       ivaux(:) = dihed%idinv(:)
       call inormlabels(dihed%ninv,ivaux(:dihed%ninv),                 &
                        dihed%idinv(:dihed%ninv),iaux)
!
       ivaux(:) = dihed%idflexi(:)
       call inormlabels(dihed%nflexi,ivaux(:dihed%nflexi),             &
                        dihed%idflexi(:dihed%nflexi),iaux)
!
! Symmetrizing dihedral terms
!
!~        if ( dihed%nrigid .gt. 0 ) then ! TODO: problem symmetrizing -179 and 179 dihedrals
!~          call symdihe(dihed%nrigid,dihed%idrigid,dihed%drigid,dihed%srigid,dihed%labrigid,.TRUE.,fsymm,debug)
!~        end if
!~ !
!~        if ( dihed%nimpro .gt. 0 ) then
!~          call symdihe(dihed%nimpro,dihed%idimpro,dihed%dimpro,dihed%simpro,dihed%labimpro,.TRUE.,fsymm,debug)
!~        end if
!
       if ( (dihed%nrigid+dihed%nimpro) .gt. 0 ) then
         call setdeps(nat,r,idat,dihed%nrigid+dihed%nimpro,            &
                      dihed%nrigid,dihed%irigid,dihed%idrigid,         &
                      dihed%drigid,dihed%srigid,dihed%labrigid,        &
                      dihed%nimpro,dihed%iimpro,dihed%idimpro,         &
                      dihed%dimpro,dihed%simpro,dihed%labimpro,        &
                      marunit,narunit,arunit,adj,lheavy,fadhoc)
       end if
!
       if ( dihed%ninv .gt. 0 ) then
         call symoop(dihed%ninv,dihed%idinv,dihed%sinv,dihed%labinv)
       end if
!
       if ( dihed%nflexi .gt. 0 ) then
         call symtor(dihed%ntor,dihed)
       end if
!
       deallocate(ivaux)
!
       return
       end subroutine symffbonded
!
!======================================================================!
!
! LABBOND - LABeling BONDs
!
! This subroutine
!
       subroutine labbond(nat,idat,nidat,nbond,ibond,idbond,nidbond,   &
                          debug)
!
       use graphtools, only:  inormlabels
!
       implicit none
!
! Input/output variables
!
       integer,dimension(2,nbond),intent(inout)  ::  ibond    !
       integer,dimension(nat),intent(in)         ::  idat     !
       integer,dimension(nbond),intent(out)      ::  idbond   !
       integer,intent(in)                        ::  nat      !
       integer,intent(in)                        ::  nbond    !
       integer,intent(in)                        ::  nidat    !
       integer,intent(out)                       ::  nidbond  !
       logical,intent(in)                        ::  debug    !
!
! Local variables
!
       integer,dimension(nbond)                  ::  ivaux    !
       integer                                   ::  itmp     !
       integer                                   ::  i,j      !
!
!  Combining labels of adjacent nodes
! -----------------------------------
!
       do i = 1, nbond
         if ( idat(ibond(1,i)) .gt. idat(ibond(2,i)) ) then
           itmp       = ibond(1,i)
           ibond(1,i) = ibond(2,i)
           ibond(2,i) = itmp
         end if
!
         ivaux(i) = idat(ibond(1,i))*nidat + idat(ibond(2,i))
       end do
!
       if ( debug ) then
         write(*,*) 'Non-normalized bond labels'
         write(*,*) '--------------------------'
         do i = 1, nbond
           write(*,'(200(1X,I3))') i,(idat(ibond(j,i)),j=1,2),ivaux(i)
         end do
         write(*,*)
       end if
!
       call inormlabels(nbond,ivaux,idbond,nidbond)
!
       if ( debug ) then
         write(*,*) 'Normalized bond labels'
         write(*,*) '----------------------'
         do i = 1, nbond
           write(*,'(200(1X,I3))') i,ivaux(i),idbond(i)
         end do
         write(*,*)
       end if
!
       return
       end subroutine labbond
!
!======================================================================!
!
! LABANGLE - LABeling ANGLEs
!
! This subroutine
!
       subroutine labangle(nat,idat,nidat,nang,iang,idang,nidang,debug)
!
       use graphtools, only:  inormlabels
!
       implicit none
!
! Input/output variables
!
       integer,dimension(3,nang),intent(inout)  ::  iang    !
       integer,dimension(nat),intent(in)        ::  idat    !
       integer,dimension(nang),intent(out)      ::  idang   !
       integer,intent(in)                       ::  nat     !
       integer,intent(in)                       ::  nang    !
       integer,intent(in)                       ::  nidat   !
       integer,intent(out)                      ::  nidang  !
       logical,intent(in)                       ::  debug   !
!
! Local variables
!
       integer,dimension(nang)                  ::  ivaux   !
       integer                                  ::  itmp    !
       integer                                  ::  i,j     !
!
!  Combining labels of adjacent nodes
! -----------------------------------
!
       do i = 1, nang
         if ( idat(iang(1,i)) .gt. idat(iang(3,i)) ) then
           itmp      = iang(1,i)
           iang(1,i) = iang(3,i)
           iang(3,i) = itmp
         end if
!
         ivaux(i) = idat(iang(2,i))*nidat**2                           &
                               + idat(iang(1,i))*nidat + idat(iang(3,i))
       end do
!
       if ( debug ) then
         write(*,*) 'Non-normalized angle labels'
         write(*,*) '---------------------------'
         do i = 1, nang
           write(*,'(200(1X,I3))') i,(idat(iang(j,i)),j=1,3),ivaux(i)
         end do
         write(*,*)
       end if
!
       call inormlabels(nang,ivaux,idang,nidang)
!
       if ( debug ) then
         write(*,*) 'Normalized angle labels'
         write(*,*) '-----------------------'
         do i = 1, nang
           write(*,'(200(1X,I3))') i,ivaux(i),idang(i)
         end do
         write(*,*)
       end if
!
       return
       end subroutine labangle
!
!======================================================================!
!
! SYMTERM - SYMmetrizing force field TERMs
!
! This subroutine
!
       subroutine symterm(nbond,idbond,dbond,sbond,labbond,fsymm)
!
       use printings
       use graphtools, only:  blockdiag,inormlabels
!
       implicit none
!
! Input/output variables
!
       character(len=50),dimension(nbond),intent(in)  ::  labbond !
       real(kind=8),dimension(nbond),intent(inout)    ::  dbond   !
       integer,dimension(nbond),intent(out)           ::  idbond  !
       integer,dimension(nbond),intent(in)            ::  sbond   !
       integer,intent(in)                             ::  nbond   !
       logical,intent(in)                             ::  fsymm   !
!
! Equivalent atoms information
!
       integer,dimension(nbond)                       ::  mol      !
       integer,dimension(nbond)                       ::  agg      !
       integer,dimension(nbond)                       ::  tag      !
       integer,dimension(nbond)                       ::  imol     !
       integer,dimension(nbond)                       ::  iagg     !
       integer,dimension(nbond)                       ::  nmol     !
       integer,dimension(nbond)                       ::  nagg     !
       integer                                        ::  magg     !  Number of aggregates
       integer                                        ::  nsize    !  Maximum aggregate size
!
! Local variables
!
       logical,dimension(nbond,nbond)                 ::  adj      !
       real(kind=8)                                   ::  daux     !
       integer                                        ::  i,j      !
       integer                                        ::  k,l      !
!
!  Symmetrizing bond distances
! ----------------------------
!
! Storing information of equivalent bond terms in an adjacency matrix
!
       adj(:,:) = .FALSE.
!
       do i = 1, nbond
         do j = 1, i
           if ( idbond(i) .eq. idbond(j) ) then
             adj(i,j) = .TRUE.
             adj(j,i) = .TRUE.
           end if
         end do
       end do
!
! Block-diagonalization of the adjacency matrix yields an array repre-
!  sentation with the equivalent bond terms information
!
       call blockdiag(nbond,adj,mol,tag,agg,nsize,nagg,iagg,           &
                      nmol,imol,magg)
!
       if ( fsymm ) then
         do i = 2, nsize
           k = imol(i)
           if ( nagg(i) .ne. 0 ) then
             do j = 1, nagg(i)
!
               daux = 0
               do l = k+1, k+i
                  daux = daux + dbond(mol(l))
               end do
!
               daux = daux/real(i)
               do l = k+1, k+i
                 dbond(mol(l)) = daux
               end do
!
               k = k + i
!
             end do
           end if
         end do
       end if
!
       do i = 2, nsize
         k = imol(i)
         if ( nagg(i) .ne. 0 ) then
           do j = 1, nagg(i)
!
             do l = k+2, k+i
               write(unideps,'(3X,I4,1X,A,1X,I4,A)')                   &
                      sbond(mol(l)),'=',sbond(mol(k+1)),'*1.d0 ; '//   &
                   trim(labbond(mol(l)))//' = '//trim(labbond(mol(k+1)))
               write(unitmp,'(3X,I4,1X,A,1X,I4,A)')                    &
                      sbond(mol(l)),'=',sbond(mol(k+1)),'*1.d0 ; '//   &
                   trim(labbond(mol(l)))//' = '//trim(labbond(mol(k+1)))
             end do
!
             k = k + i
!
           end do
         end if
       end do
!
       return
       end subroutine symterm
!
!======================================================================!
!
! SYMDIHE - SYMmetrizing force field DIHEdrals
!
! This subroutine
!
       subroutine symoop(ndihe,iddihe,sdihe,labdihe)
!
       use printings
       use graphtools, only:  blockdiag,inormlabels
!
       implicit none
!
! Input/output variables
!
       character(len=50),dimension(ndihe),intent(in)  ::  labdihe !
       integer,dimension(ndihe),intent(out)           ::  iddihe  !
       integer,dimension(ndihe),intent(in)            ::  sdihe   !
       integer,intent(in)                             ::  ndihe   !
!
! Equivalent atoms information
!
       integer,dimension(ndihe)                       ::  mol      !
       integer,dimension(ndihe)                       ::  agg      !
       integer,dimension(ndihe)                       ::  tag      !
       integer,dimension(ndihe)                       ::  imol     !
       integer,dimension(ndihe)                       ::  iagg     !
       integer,dimension(ndihe)                       ::  nmol     !
       integer,dimension(ndihe)                       ::  nagg     !
       integer                                        ::  magg     !  Number of aggregates
       integer                                        ::  nsize    !  Maximum aggregate size
!
! Local variables
!
       logical,dimension(ndihe,ndihe)                 ::  adj      !
       integer                                        ::  i,j      !
       integer                                        ::  k,l      !
!
!  Symmetrizing bond distances
! ----------------------------
!
! Storing information of equivalent bond terms in an adjacency matrix
!
       adj(:,:) = .FALSE.
!
       do i = 1, ndihe
         do j = 1, i
           if ( iddihe(i) .eq. iddihe(j)  ) then
             adj(i,j) = .TRUE.
             adj(j,i) = .TRUE.
           end if
         end do
       end do
!
! Block-diagonalization of the adjacency matrix yields an array repre-
!  sentation with the equivalent bond terms information
!
       call blockdiag(ndihe,adj,mol,tag,agg,nsize,nagg,iagg,           &
                      nmol,imol,magg)
!
       do i = 2, nsize
         k = imol(i)
         if ( nagg(i) .ne. 0 ) then
           do j = 1, nagg(i)
!
!             if ( fsymm ) then ! TODO: symmetrize absolute value but keep sign
!             daux = 0
!             do l = k+1, k+i
!                daux = daux + ddihe(mol(l))
!             end do
!
!             daux = daux/real(i)
!             do l = k+1, k+i
!               ddihe(mol(l)) = daux
!             end do
!             end if
!
             do l = k+2, k+i
               write(unideps,'(3X,I4,1X,A,1X,I4,A)')                   &
                      sdihe(mol(l)),'=',sdihe(mol(k+1)),'*1.d0 ; '//   &
                   trim(labdihe(mol(l)))//' = '//trim(labdihe(mol(k+1)))
               write(unitmp,'(3X,I4,1X,A,1X,I4,A)')                    &
                      sdihe(mol(l)),'=',sdihe(mol(k+1)),'*1.d0 ; '//   &
                   trim(labdihe(mol(l)))//' = '//trim(labdihe(mol(k+1)))
             end do
!
             k = k + i
!
           end do
         end if
       end do
!
       return
       end subroutine symoop
!
!======================================================================!
!
! SYMTOR - SYMmetrizing force field TORsions
!
! This subroutine
!
       subroutine symtor(ndihe,dihe)
!
       use datatypes,  only:  dihedrals
!
       use printings
       use graphtools, only:  inormlabels
!
       implicit none
!
! Input/output variables
!
       type(dihedrals),intent(in)                     ::  dihe    !
       integer,intent(in)                             ::  ndihe   !
!
! Local variables
!
       character(len=50),dimension(ndihe)             ::  labdihe  !
       logical,dimension(ndihe)                       ::  lquad    !
       logical,dimension(ndihe)                       ::  ldep     !
       integer,dimension(4,ndihe)                     ::  idihe    !
       integer,dimension(ndihe)                       ::  sdihe    !
       integer,dimension(ndihe)                       ::  mdihe    !
       integer,dimension(ndihe)                       ::  iddihe   !
       integer,dimension(ndihe)                       ::  iflexterm!
       integer,dimension(ndihe,0:ndihe)               ::  qterm    !
       integer,dimension(ndihe)                       ::  ivaux    !
       integer                                        ::  idep     !
       integer                                        ::  niddihe  !
       integer                                        ::  nterm    !
       integer                                        ::  i,j,k    !
!
!  Symmetrizing flexible dihedrals
! --------------------------------
!
! Generating unique representation with all torsions
!
       k = 0
       do i = 1, dihe%nflexi
         do j = 1, dihe%flexi(i)%ntor
           k = k + 1
!
           sdihe(k)   = dihe%flexi(i)%tor(j)%stor
           mdihe(k)   = dihe%flexi(i)%tor(j)%multi
           idihe(:,k) = dihe%flexi(i)%itor(:)
           labdihe(k) = dihe%flexi(i)%tor(j)%labtor
           iflexterm(k) = i
!
           ivaux(k) = dihe%idflexi(i)*10 + mdihe(k)
!~ write(*,*) k,trim(labdihe(k)),dihe%idflexi(i),mdihe(k),ivaux(k)
!
         end do
       end do
       nterm = k
!
       call inormlabels(nterm,ivaux(:nterm),iddihe(:nterm),niddihe)
!~        k = 0
!~        do i = 1, dihe%nflexi
!~          do j = 1, dihe%flexi(i)%ntor
!~            k = k + 1
!~ !
!~ write(*,*)'normalized', k,trim(labdihe(k)),dihe%idflexi(i),mdihe(k),iddihe(k)
!~ !
!~          end do
!~        end do
!
! Finding position of principal quadruplets in the new list
!
       lquad(:) = .FALSE.
       ldep(:)  = .FALSE.
       qterm(:,:) = 0
       do i = 1, dihe%nquad
         do j = 1, nterm
           if ( (idihe(1,j).eq.dihe%iquad(1,i))                        &
                 .and. (idihe(2,j).eq.dihe%iquad(2,i))                 &
                 .and. (idihe(3,j).eq.dihe%iquad(3,i))                 &
                 .and. (idihe(4,j).eq.dihe%iquad(4,i)) ) then
             lquad(j) = .TRUE.
             qterm(i,mdihe(j)) = j
           end if
         end do
       end do
!
! Printing explicit dependencies for extra ring-exocyclic quadruplets
! selected by genquad but not used as independent scan coordinates.
!
       if ( allocated(dihe%depquad) ) then
         do i = 1, nterm
           idep = dihe%depquad(iflexterm(i))
           if ( (idep.gt.0).and.(idep.le.dihe%nquad) ) then
             if ( qterm(idep,mdihe(i)) .gt. 0 ) then
               write(unideps,'(3X,I4,1X,A,1X,I4,A)')                   &
                         sdihe(i),'=',sdihe(qterm(idep,mdihe(i))),     &
                                  '*1.d0 ; '//trim(labdihe(i))//' = '//&
                                  trim(labdihe(qterm(idep,mdihe(i))))
               write(uniscr,'(3X,I4,1X,A,1X,I4,A)')                    &
                         sdihe(i),'=',sdihe(qterm(idep,mdihe(i))),     &
                                  '*1.d0 ; '//trim(labdihe(i))//' = '//&
                                  trim(labdihe(qterm(idep,mdihe(i))))
               write(unitmp,'(3X,I4,1X,A,1X,I4,A)')                    &
                         sdihe(i),'=',sdihe(qterm(idep,mdihe(i))),     &
                                  '*1.d0 ; '//trim(labdihe(i))//' = '//&
                                  trim(labdihe(qterm(idep,mdihe(i))))
               ldep(i) = .TRUE.
             end if
           end if
         end do
       end if
!
! Printing equivalent dihedral terms with same identifier based on
!  atomtypes and multiplicities
!
       do i = 1, nterm
         if ( lquad(i) ) then
           do j = 1, nterm
             if ( lquad(j) ) cycle
             if ( ldep(j) ) cycle
             if ( iddihe(i) .eq. iddihe(j) ) then
!~ write(*,*) 'matching',i,iddihe(i),'with',j,iddihe(j)
               write(unideps,'(3X,I4,1X,A,1X,I4,A)')                   &
                                  sdihe(j),'=',sdihe(i),'*1.d0 ; '//   &
                               trim(labdihe(j))//' = '//trim(labdihe(i))
               write(uniscr,'(3X,I4,1X,A,1X,I4,A)')                    &
                                  sdihe(j),'=',sdihe(i),'*1.d0 ; '//   &
                               trim(labdihe(j))//' = '//trim(labdihe(i))
               write(unitmp,'(3X,I4,1X,A,1X,I4,A)')                    &
                                  sdihe(j),'=',sdihe(i),'*1.d0 ; '//   &
                               trim(labdihe(j))//' = '//trim(labdihe(i))
! A later principal of the same type must not assign this target again.
               ldep(j) = .TRUE.
             end if
           end do
         end if
       end do
!
       return
       end subroutine symtor
!
!======================================================================!
!
! SETDEPS - SET DEPendencieS
!
! This subroutine
!
       subroutine setdeps(nat,r,idat,ndihe,nrigid,irigid,idrigid,      &
                           drigid,srigid,labrigid,nimpro,iimpro,       &
                           idimpro,dimpro,simpro,labimpro,marunit,     &
                           narunit,arunit,adj,lheavy,fadhoc)
!
       use printings
       use graphtools, only:  blockdiag,inormlabels
!
       implicit none
!
! Input/output variables
!
       real(kind=8),dimension(nrigid),intent(in)       ::  drigid    !
       real(kind=8),dimension(nimpro),intent(in)       ::  dimpro    !
       character(len=50),dimension(nrigid),intent(in)  ::  labrigid  !
       character(len=50),dimension(nimpro),intent(in)  ::  labimpro  !
       integer,dimension(4,nrigid),intent(in)          ::  irigid    !
       integer,dimension(4,nimpro),intent(in)          ::  iimpro    !
       integer,dimension(nat),intent(in)               ::  idat      !
       integer,dimension(nrigid),intent(in)            ::  idrigid   !
       integer,dimension(nimpro),intent(in)            ::  idimpro   !
       integer,dimension(nrigid),intent(in)            ::  srigid    !
       integer,dimension(nimpro),intent(in)            ::  simpro    !
       integer,dimension(r,nat),intent(in)             ::  arunit    !
       integer,dimension(r),intent(in)                 ::  narunit   !
       integer,intent(in)                              ::  nat       !
       integer,intent(in)                              ::  r         !
       integer,intent(in)                              ::  ndihe     !
       integer,intent(in)                              ::  nrigid    !
       integer,intent(in)                              ::  nimpro    !
       integer,intent(in)                              ::  marunit   !
       logical,dimension(nat,nat),intent(in)           ::  adj       !
       logical,dimension(nat),intent(in)               ::  lheavy    !
       logical,intent(in)                              ::  fadhoc !
!
! Equivalent atoms information
!
       real(kind=8),dimension(ndihe)                   ::  ddihe     !
       character(len=50),dimension(ndihe)              ::  labdihe   !
       integer,dimension(4,ndihe)                      ::  idihe     !
       integer,dimension(ndihe)                        ::  iddihe    !
       integer,dimension(ndihe)                        ::  sdihe     !
       integer,dimension(ndihe)                        ::  ivaux     !
       integer,dimension(ndihe)                        ::  mol       !
       integer,dimension(ndihe)                        ::  agg       !
       integer,dimension(ndihe)                        ::  tag       !
       integer,dimension(ndihe)                        ::  imol      !
       integer,dimension(ndihe)                        ::  iagg      !
       integer,dimension(ndihe)                        ::  nmol      !
       integer,dimension(ndihe)                        ::  nagg      !
       integer                                         ::  magg      !  Number of aggregates
       integer                                         ::  nsize     !  Maximum aggregate size
       integer                                         ::  niddihe   !
!
! Local variables
!
       logical,dimension(ndihe,ndihe)                     ::  eqvadj    !
       logical,dimension(ndihe)                           ::  visitdi   !
       logical,dimension(nat)                             ::  visited   !
       logical,dimension(nat)                             ::  lunit     !
       logical,dimension(nimpro)                          ::  chkimpro  !
       integer,dimension(2,nat*nat)                       ::  ibond     !
       integer,dimension(4,ndihe)                         ::  iquad     !
       integer,dimension(ndihe)                           ::  idxdihe   !
       integer,dimension(ndihe)                           ::  itype     !
       integer,dimension(nat)                             ::  nei1      !
       integer,dimension(nat)                             ::  nei2      !
       integer,dimension(4)                               ::  iidx      !
       integer                                            ::  nquad     !
       integer                                            ::  nbond     !
       integer                                            ::  inei1     !
       integer                                            ::  inei2     !
       integer                                            ::  itmp      !
       integer                                            ::  i,j       !
       integer                                            ::  k,l       !
!
! Finding equivalencies between improper and rigid dihedrals
! ..........................................................
!
! iimpro(1:2) is the deterministic center--reference bond selected in
! setdihe.  The aromatic matching below intentionally uses that pair.
! Generating unique representation with rigid and improper dihedrals
!
! Merge rigid dihedrals and OOP impropers into a single list so dependency
! groups can connect them in the same Joyce $dependence block.
!
       do i = 1, nrigid
         ivaux(i)   = idrigid(i)*2 + 1
         sdihe(i)   = srigid(i)
         ddihe(i)   = drigid(i)
         idihe(:,i) = irigid(:,i)
         labdihe(i) = labrigid(i)
       end do
!
       do i = 1, nimpro
         ivaux(nrigid+i)   = idimpro(i)*2 + 2
         sdihe(nrigid+i)   = simpro(i)
         ddihe(nrigid+i)   = dimpro(i)
         idihe(:,nrigid+i) = iimpro(:,i)
         labdihe(nrigid+i) = labimpro(i)
       end do
!
       call inormlabels(ndihe,ivaux,iddihe,niddihe)
!
       do i = 1, ndihe
         iddihe(i) = iddihe(i) + 4
       end do
       niddihe = niddihe + 4
!
! Storing information of equivalent dihedral terms based on atomtypes
!  in an adjacency matrix
!
! Terms with the same normalized atom-type identifier and compatible cis/trans
! geometry are initially considered equivalent.
!
       eqvadj(:,:) = .FALSE.
!
       do i = 1, ndihe
         do j = 1, i
           if ( (iddihe(i).eq.iddihe(j)) .and.                         &
                      (abs(abs(ddihe(i))-abs(ddihe(j))).lt.5.0d0) ) then ! cis/trans dihedrals have different type
             eqvadj(i,j) = .TRUE.
             eqvadj(j,i) = .TRUE.
           end if
         end do
       end do
!
! Setting ad hoc dependencies
! ---------------------------
!
       if ( fadhoc ) then

! The ad hoc aromatic pass adds chemical dependencies that are not captured by
! simple atom-type equivalence, especially rigid aromatic terms coupled to OOP
! impropers around the same center--reference bond.

       chkimpro(:) = .TRUE.
!
       do i = 1, marunit
!
! Generating list of atoms belonging to the aromatic unit
!
! lunit marks the current fused/aromatic unit so generated quadruplets can be
! classified by how many atoms belong to that unit.
!
         lunit(:) = .FALSE.
!
         do j = 1, narunit(i)
           lunit(arunit(i,j)) = .TRUE.
         end do
!
! Generating list of bonds belonging to the aromatic unit
!
         nbond = 0
         do j = 1, narunit(i)-1
           do k = j+1, narunit(i)
             if ( adj(arunit(i,j),arunit(i,k)) ) then
!
               nbond = nbond + 1
               ibond(1,nbond) = arunit(i,j)
               ibond(2,nbond) = arunit(i,k)
!
               if ( idat(ibond(1,nbond)) .gt. idat(ibond(2,nbond)) ) then
                 itmp           = ibond(1,nbond)
                 ibond(1,nbond) = ibond(2,nbond)
                 ibond(2,nbond) = itmp
               end if
!
             end if
           end do
         end do
!
! Generating all possible quadruplets from the bonds in the aromatic unit
!
! For every aromatic-unit bond, combine every neighbour of the first central
! atom with every neighbour of the second central atom.  The resulting list is
! a superset; invalid or absent terms are discarded later.
!
         nquad = 0
         do j = 1, nbond
! Finding atoms bonded to first central atom
           visited(:) = .FALSE.
!
           visited(ibond(1,j)) = .TRUE.
           visited(ibond(2,j)) = .TRUE.
!
           inei1 = 0
           do k = 1, nat
             if ( .NOT. visited(k) ) then
               if ( adj(k,ibond(1,j)) ) then
                 inei1 = inei1 + 1
                 nei1(inei1) = k
               end if
             end if
           end do
! Finding atoms bonded to second central atom
           visited(:) = .FALSE.
!
           visited(ibond(1,j)) = .TRUE.
           visited(ibond(2,j)) = .TRUE.
!
           inei2 = 0
           do k = 1, nat
             if ( .NOT. visited(k) ) then
               if ( adj(k,ibond(2,j)) ) then
                 inei2 = inei2 + 1
                 nei2(inei2) = k
               end if
             end if
           end do
!
           do k = 1, inei1
             do l = 1, inei2
               nquad = nquad + 1
               iquad(1,nquad) = nei1(k)
               iquad(2,nquad) = ibond(1,j)
               iquad(3,nquad) = ibond(2,j)
               iquad(4,nquad) = nei2(l)
             end do
           end do
!
         end do
!
! Finding equivalencies between selected ICs
!
! Normalize each generated quadruplet to the same orientation used by the real
! rigid/improper list, then try to find the actual term index.
!
         visitdi(:) = .FALSE.
         idxdihe(:) = -1
!
         itype(:) = -1
!
         do j = 1, nquad
!
           iidx(:) = iquad(:,j)
!
           if ( idat(iidx(2)) .gt. idat(iidx(3)) ) then
!
             itmp    = iidx(2)
             iidx(2) = iidx(3)
             iidx(3) = itmp
!
             itmp    = iidx(1)
             iidx(1) = iidx(4)
             iidx(4) = itmp
!
           else if ( idat(iidx(2)) .eq. idat(iidx(3)) ) then
             if ( idat(iidx(1)) .gt. idat(iidx(4)) ) then
!
               itmp    = iidx(2)
               iidx(2) = iidx(3)
               iidx(3) = itmp
!
               itmp    = iidx(1)
               iidx(1) = iidx(4)
               iidx(4) = itmp
!
             end if
           end if
!
           iquad(:,j) = iidx(:)
!
!!!write(*,*) 'TARGET',iidx(:),':',idat(iidx(1)),idat(iidx(2)),idat(iidx(3)),idat(iidx(4))
! Find index of actual quadruplet in original list
           do k = 1, nrigid
!!!write(*,*) 'checking',idihe(:,k)
             if ( .NOT. visitdi(k) ) then
               if ( (iidx(1).eq.idihe(1,k))                            &
                    .and. (iidx(2).eq.idihe(2,k))                      &
                    .and. (iidx(3).eq.idihe(3,k))                      &
                    .and. (iidx(4).eq.idihe(4,k)) ) then
                 idxdihe(j) = k
                 visitdi(k) = .TRUE.
!!!write(*,*) 'ORIGINAL DIHEDRAL FOUND'
                 exit
               end if
             end if
           end do ! FIXME: if idxdihe(j) == -1 there is a problem ; target dihedral is not found in original list
!!!write(*,*)
!
! Assign a type to the actual quadruplet
! --------------------------------------
!
           ! A generated aromatic quadruplet can be absent from the
           ! rigid list after geometric filtering. It has no ad hoc
           ! rigid/improper counterpart to classify in that case.
           if ( idxdihe(j) .lt. 1 ) cycle

!   i) backbone (equivalent to C-C-C-X and impropers)
!  ii) interunit (C-C-C-C [trans], C-C-C-H [cis])
! iii) C-C-C-H [trans] (equivalent to impropers)
!  iv) H-C-C-H
!   v) X-C-C-Y (based on atomtypes)
!
           if ( lunit(iidx(1)) .and. lunit(iidx(2))                    &
                        .and. lunit(iidx(3)) .and. lunit(iidx(4)) ) then
!
! If all the atoms belong to the same aromatic unit the dihedral can be
!  i) backbone   (cis)
! ii) inter-unit (trans)
!
             if ( abs(abs(ddihe(idxdihe(j)))-180.0d0) .lt. 5.0d0 ) then
! Dihedral associated to inter-unit flexibility (trans)
               itype(idxdihe(j)) = 2
!~ write(*,*) 'Dihedral ',trim(labdihe(idxdihe(j))),' is type',2
             else
! Dihedral is an ad hoc backbone (cis)
               itype(idxdihe(j)) = 1
!~ write(*,*) 'Dihedral ',trim(labdihe(idxdihe(j))),' is type',1
             end if
!
! If three atoms belong to the same aromatic unit the dihedral can be
!  i) C-C-C-H equivalent to oop
! ii) C-C-C-X equivalent to backbone
!
           else if ( (lunit(iidx(1)) .and. lunit(iidx(2))              &
                                             .and. lunit(iidx(3)))     &
                 .or. (lunit(iidx(2)) .and. lunit(iidx(3))             &
                                            .and. lunit(iidx(4))) ) then
!~ write(*,*) 'Dihedral ',trim(labdihe(idxdihe(j))),' has 3 atoms in the aromatic unit',(lheavy(iidx(l)),l=1,4)
             if ( abs(abs(ddihe(idxdihe(j)))-180.0d0) .lt. 5.0d0 ) then
!~ write(*,*) 'Dihedral ',trim(labdihe(idxdihe(j))),' is trans'
! Dihedrals C-C-C-X associated with impropers (trans)
               if ( (.NOT.lheavy(iidx(1))) .or. (.NOT.lheavy(iidx(4))) ) then
!~ write(*,*) 'Dihedral ',trim(labdihe(idxdihe(j))),' is C-C-C-H'
! C-C-C-H dihedrals are equivalent to the corresponding oop improper
                 itype(idxdihe(j)) = 3
!~ write(*,*) 'Dihedral ',trim(labdihe(idxdihe(j))),' is type',3
!
! Set up equivalencies with corresponding improper dihedrals
!
                 do k = 1, nimpro
                   if ( chkimpro(k) ) then ! TODO: alternatively check if iimpro(2) is H or nonH
                     if ( ((iidx(2).eq.iimpro(1,k)).and.(iidx(1).eq.iimpro(2,k))) &
                         .or. ((iidx(3).eq.iimpro(1,k)).and.(iidx(4).eq.iimpro(2,k))) ) then
                       itype(nrigid+k) = 3
!~ write(*,*) '  Improper ',trim(labdihe(nrigid+k)),' is type',3
                       chkimpro(k) = .FALSE.
                     end if
                   end if
                 end do
!
               else
! C-C-C-X trans dihedrals are ad hoc backbone
!~ write(*,*) 'Dihedral ',trim(labdihe(idxdihe(j))),' is C-C-C-X'
                 itype(idxdihe(j)) = 1
!~ write(*,*) 'Dihedral ',trim(labdihe(idxdihe(j))),' is type',1
!
! Set up equivalencies with corresponding improper dihedrals
!
                 do k = 1, nimpro ! TODO: alternatively check if iimpro(2) is H or nonH
                   if ( chkimpro(k) ) then
                     if ( ((iidx(2).eq.iimpro(1,k)).and.(iidx(1).eq.iimpro(2,k))) &
                         .or. ((iidx(3).eq.iimpro(1,k)).and.(iidx(4).eq.iimpro(2,k))) ) then
                       itype(nrigid+k) = 1
!~ write(*,*) '  Improper ',trim(labdihe(nrigid+k)),' is type',1
                       chkimpro(k) = .FALSE.
                     end if
                   end if
                 end do
!
               end if
             else
! Dihedral C1-C2-C2-X associated to inter-unit flexibility (cis)
!~ write(*,*) 'Dihedral ',trim(labdihe(idxdihe(j))),' is cis'
               itype(idxdihe(j)) = 2
!~ write(*,*) 'Dihedral ',trim(labdihe(idxdihe(j))),' is type',2
             end if
!
           else
!
! X-C-C-Y type dihedrals
!  i) H-C-C-H dihedrals are equivalent
! ii) H-C-C-X dihedrals are based on atomtypes
!
!~ write(*,*) 'Dihedral ',trim(labdihe(idxdihe(j))),' only central atoms in the aromatic unit',(lheavy(iidx(l)),l=1,4)
             if ( (.NOT.lheavy(iidx(1))) .and. (.NOT.lheavy(iidx(4))) ) then
               itype(idxdihe(j)) = 4
!~ write(*,*) 'Dihedral ',trim(labdihe(idxdihe(j))),' is H-C-C-H'
!~ write(*,*) 'Dihedral ',trim(labdihe(idxdihe(j))),' is type',4
             end if
!
           end if
!
!~ write(*,*)
!write(*,*) 'Dihedral ',idxdihe(j),':',iidx(:),trim(labdihe(idxdihe(j))),' is type',itype(idxdihe(j))
!
         end do
!
! Setting dependencies between dihedrals within the same aromatic unit
!
! Terms assigned the same ad hoc type are connected in eqvadj.  The later block
! diagonalization converts each connected component into Joyce dependencies.
!
         do j = 1, ndihe
           if ( itype(j) .eq. -1 ) cycle
           do k = 1, j
             if ( itype(k) .eq. -1 ) cycle
             if ( itype(j) .eq. itype(k)  ) then
               eqvadj(j,k) = .TRUE.
               eqvadj(k,j) = .TRUE.
             end if
           end do
         end do
!
       end do
       end if
!
! Block-diagonalization of the adjacency matrix yields an array repre-
!  sentation with the equivalent bond terms information
!
! Each connected component of eqvadj becomes one dependency family.  The first
! term in each family is kept independent and all others are assigned to it.
!
       call blockdiag(ndihe,eqvadj,mol,tag,agg,nsize,nagg,iagg,        &
                      nmol,imol,magg)
!
       do i = 2, nsize
         k = imol(i)
         if ( nagg(i) .ne. 0 ) then
           do j = 1, nagg(i)
!
             do l = k+2, k+i
               write(unideps,'(3X,I4,1X,A,1X,I4,A)')                   &
                      sdihe(mol(l)),'=',sdihe(mol(k+1)),'*1.d0 ; '//   &
                   trim(labdihe(mol(l)))//' = '//trim(labdihe(mol(k+1)))
               write(unitmp,'(3X,I4,1X,A,1X,I4,A)')                    &
                      sdihe(mol(l)),'=',sdihe(mol(k+1)),'*1.d0 ; '//   &
                   trim(labdihe(mol(l)))//' = '//trim(labdihe(mol(k+1)))
             end do
!
             k = k + i
!
           end do
         end if
       end do
!
!
       return
       end subroutine setdeps
!
!======================================================================!
!
       end module dependencies
!
!======================================================================!

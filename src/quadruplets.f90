!======================================================================!
!
       module quadruplets
!
       implicit none
!
       private
       public  ::  genquad,                                            &
                   selectquad,                                         &
                   screenquad
!
       contains
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
!~ write(*,*) 'Dihedral',i,'is a CH3 rotation',dihe%iflexi(:,i)
           do j = 1, dihe%nquad
             ivaux2(:) = dihe%iquad(:,j)
             if ( ((idat(ich3(1,i)).eq.idat(ivaux2(2))).and.             &
                   (idat(ich3(2,i)).eq.idat(ivaux2(3))).and.             &
                   (idat(ich3(3,i)).eq.idat(ivaux2(4)))                  &
             .or. ((idat(ich3(1,i)).eq.idat(ivaux2(3))).and.             &
                   (idat(ich3(2,i)).eq.idat(ivaux2(2))).and.             &
                   (idat(ich3(3,i)).eq.idat(ivaux2(1))))) ) then
!~ write(*,*) '  Overlap with princial quadruplet',j,':',dihe%iquad(:,j)
!~ write(*,*)
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
       imap(:)    = -1
       iflexi(:)  = -1
       iselect(:) = -1
       sselect(:) = 0
       nselect    = 0
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
! Each group contains all flexible torsions around the same central bond.
! Prefer a scan coordinate with heavy atoms at both ends, because such a
! torsion usually gives the most chemically informative relaxed scan.
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
!
! If no heavy-heavy terminal choice exists, allow one terminal hydrogen.  This
! still avoids the least informative H--X--Y--H choice when possible.
!
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
!
! Last fallback: all available candidates contain hydrogens at both terminal
! positions, so keep them rather than losing the rotatable bond entirely.
!
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
! geometrically compatible with the requested cis/cis-trans cases
!
! The special ring-exocyclic logic may return two chemically required torsions:
! one independent scan coordinate and one dependent torsion preserved in the
! topology through depquad
!
         call selringexo(nat,coord,adj,ideg,lcycle,lrigid,tmpdi,k,     &
                         neqquad(i),ndihe,imap,nmap,iadd,nadd)
!
         if ( nadd .gt. 0 ) then
           nselect = nselect + 1
!
           iselect(nselect) = imap(1)
!
! Choose the independent member of the ring-exocyclic pair using the same
! angular priority as leading_quad.  The remaining returned members become
! explicit dependencies of this selected scan.
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
!
! depquad is indexed by the original flexible-torsion list
! tmpmap translates the compact grouped index back to that original index
!
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
       subroutine selringexo(nat,coord,adj,ideg,lcycle,lrigid,tmpdi,  &
                             first,nitem,ndihe,imap,nmap,iout,nout)

       use datatypes, only: dihedrals
!
       implicit none
!
       type(dihedrals),intent(in)               ::  tmpdi   !
       real(kind=8),dimension(3,nat),intent(in) ::  coord   !
       logical,dimension(nat,nat),intent(in)    ::  adj     !
       logical,dimension(nat,nat),intent(in)    ::  lcycle  !
       logical,dimension(nat,nat),intent(in)    ::  lrigid  !
       integer,dimension(nat),intent(in)        ::  ideg    !
       integer,intent(in)                       ::  nat     !
       integer,intent(in)                       ::  ndihe   !
       integer,dimension(ndihe),intent(in)      ::  imap    !
       integer,dimension(ndihe),intent(out)     ::  iout    !
       integer,intent(in)                       ::  first   !
       integer,intent(in)                       ::  nitem   !
       integer,intent(in)                       ::  nmap    !
       integer,intent(out)                      ::  nout    !
!
       real(kind=8)                             ::  daux    !
       logical                                  ::  ring1   !
       logical                                  ::  ring2   !
       integer                                  ::  c1      !
       integer                                  ::  c2      !
       integer                                  ::  exo     !
       integer                                  ::  ring    !
       integer                                  ::  ncis    !
       integer                                  ::  ntrans  !
       integer                                  ::  j       !
!
       iout(:) = -1
       nout = 0
!
! A single candidate cannot create a dependent ring-exocyclic pair.
!
       if ( nmap .le. 1 ) return
!
! The central bond is common to every item in this group.  Identify which side
! belongs to a ring and which side points to the exocyclic substituent.
!
       c1 = tmpdi%iflexi(2,first+1)
       c2 = tmpdi%iflexi(3,first+1)
!
       ring1 = lringtermi(tmpdi,lcycle,c1,c2,first,nitem)
       ring2 = lringtermi(tmpdi,lcycle,c2,c1,first,nitem)
!
! The special rule applies only to true ring--exocyclic bonds.  If both sides
! are ring-like or both are non-ring, fall back to ordinary leading_quad.
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
       if ( lsp2exo(nat,coord,adj,ideg,lrigid,exo,ring) ) then
!
! Degree-3 sp2-like exocyclic centers are handled by selecting the two almost
! cis candidates.  They represent the pair that must be kept together.
!
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
! Checking if candidate is a valid exocyclic atom with degree 2
!
       if ( ldeg2exo(nat,adj,ideg,exo,ring) ) then
!
! Degree-2 exocyclic centers are accepted only when the complete group is an
! unambiguous cis/trans pair.
!
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
       end subroutine selringexo
!
!======================================================================!
!
       logical function lringtermi(tmpdi,lcycle,center,other,  &
                                   first,nitem)

       use datatypes, only: dihedrals
!
       implicit none
!
! Input/output variables
!
       type(dihedrals),intent(in)          ::  tmpdi   !
       logical,dimension(:,:),intent(in)   ::  lcycle  !
       integer,intent(in)                  ::  center  !
       integer,intent(in)                  ::  other   !
       integer,intent(in)                  ::  first   !
       integer,intent(in)                  ::  nitem   !
!
! Local variables
!
       logical                             ::  match   !
       integer,dimension(4)                ::  iquad   !
       integer                             ::  tcenter !
       integer                             ::  tother  !
       integer                             ::  j       !
!
       lringtermi = .FALSE.
!
       do j = 1, nitem
         iquad(:) = tmpdi%iflexi(:,first+j)
!
! Work out the terminal atom attached to "center" regardless of whether the
! stored quadruplet is oriented center--other or other--center.
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
         if ( match .and. lcycle(center,tcenter) ) then
!
! A terminal connected through a ring bond marks this side as the ring side of
! the rotatable bond.
!
           lringtermi = .TRUE.
           return
         end if
       end do
!
       return
       end function lringtermi
!
!======================================================================!
!
       logical function lsp2exo(nat,coord,adj,ideg,lrigid,exo,ring)
!
       use genfftools,  only:  calc_angle
!
       implicit none
!
! Input/output variables
!
       real(kind=8),dimension(3,nat),intent(in)  ::  coord        !
       logical,dimension(nat,nat),intent(in)     ::  adj          !
       logical,dimension(nat,nat),intent(in)     ::  lrigid       !
       integer,dimension(nat),intent(in)         ::  ideg         !
       integer,intent(in)                        ::  nat          !
       integer,intent(in)                        ::  exo          !
       integer,intent(in)                        ::  ring         !
!
! Local variables
!
       real(kind=8),dimension(3)                 ::  angle        !
       real(kind=8)                              ::  asum         !
       integer,dimension(3)                      ::  nei          !
       integer                                   ::  n            !
       integer                                   ::  j            !
       logical                                   ::  lplanardeg3  !
       logical                                   ::  exobond      !
!
       real(kind=8),parameter                    ::  pi = 4*atan(1.0_8)
!
! Checking if exocyclic atom has sp2 hibridization
!
       lsp2exo = .FALSE.
!
! This route is only meant for trigonal exocyclic atoms.  Degree-2 centers are
! handled separately by ldeg2exo.
!
       if ( ideg(exo) .ne. 3 ) return
!
! Checking if exocyclic atom is planar and has degree 3
!
       lplanardeg3 = .FALSE.
!
       n = 0
       do j = 1, nat
         if ( adj(exo,j) ) then
!
! Store the three neighbours used to test planarity around the exocyclic atom.
!
           n = n + 1
           if ( n .le. 3 ) nei(n) = j
         end if
       end do
!
       if ( n .ne. 3 ) return
!
       angle(1) = calc_angle(coord(:,nei(1)),coord(:,exo),             &
                             coord(:,nei(2)))*180.0d0/pi
       angle(2) = calc_angle(coord(:,nei(1)),coord(:,exo),             &
                             coord(:,nei(3)))*180.0d0/pi
       angle(3) = calc_angle(coord(:,nei(2)),coord(:,exo),             &
                             coord(:,nei(3)))*180.0d0/pi
!
       asum = angle(1) + angle(2) + angle(3)
!
! A nearly planar degree-3 center is treated as sp2-like.  The max-angle guard
! avoids accepting almost-linear geometries as trigonal planar.
!
       lplanardeg3 = (abs(asum-360.0d0).le.25.0d0)                     &
                                        .and. (maxval(angle).lt.170.0d0)
!
! Checking if cycle has rigid exocyclic bond
!
       exobond = .FALSE.
!
       do j = 1, nat
         if ( j .eq. ring ) cycle
         if ( adj(exo,j) .and. lrigid(exo,j) ) then
!
! A rigid bond to the non-ring substituent is a second indicator of sp2-like
! conjugation even if the angular planarity test is marginal.
!
           exobond = .TRUE.
           exit
         end if
       end do
!
! Checking if exocyclic atom has sp2 hybridization
!
       lsp2exo = lplanardeg3 .or. exobond
!
       return
       end function lsp2exo
!
!======================================================================!
!
       logical function ldeg2exo(nat,adj,ideg,exo,ring)
!
       implicit none
!
! Input/Output variables
!
       logical,dimension(nat,nat),intent(in)  ::  adj    !
       integer,dimension(nat),intent(in)      ::  ideg   !
       integer,intent(in)                     ::  nat    !
       integer,intent(in)                     ::  exo    !
       integer,intent(in)                     ::  ring   !
!
! Local variables
!
       integer                                ::  other  !
       integer                                ::  j      !
!
!
!
       ldeg2exo = .FALSE.
!
       if ( ideg(exo) .ne. 2 ) return
!
       other = -1
!
! For a degree-2 exocyclic center, require exactly one additional neighbour
! besides the ring atom.  The chemical identity of that neighbour is not fixed.
!
       do j = 1, nat
         if ( j .eq. ring ) cycle
         if ( adj(exo,j) ) then
           other = j
           exit
         end if
       end do
!
       if ( other .lt. 1 ) return
!
! For degree-2 ring-exocyclic bonds the cis/trans selection is decided
!  by selringexo itself.  Reaching this point only states that
!  the exocyclic atom has exactly one additional neighbour.
!
       ldeg2exo = .TRUE.
!
       return
       end function ldeg2exo
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
       end module quadruplets
!
!======================================================================!

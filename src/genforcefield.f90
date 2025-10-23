!======================================================================!
!
       module genforcefield
!
       use lengths,       only:  leninp,lenlab
       use units,         only:  uniic,unideps,unitmp,uniscr,unidihe,  &
                                 unijoyce,uninb,unitop,unicorr
!
       implicit none
!
       private
       public  ::  gentypes,                                           &
                   genffbonded,                                        &
                   symffbonded,                                        &
                   eqvatms,                                            &
                   nearnei
!
       contains
!
!======================================================================!
!
! EQVATOMS - chemically EQuiValent AToMS algorithm
!
! This subroutine 
!
       subroutine eqvatms(nat,adj,ilab,nlab,eqv)
!
       use isomorphism
       use graphtools,  only : calcdegundir,dnormlabels
!
       implicit none
!
! Input/output variables
!
       logical,dimension(nat,nat),intent(in)   ::  adj     !  Adjacency matrix
       logical,dimension(nat,nat),intent(out)  ::  eqv     !  Equivalent atoms relationships
       integer,dimension(nat),intent(in)       ::  ilab    !  Integers associated to atomic labels
       integer,intent(in)                      ::  nlab    !  Number of different nuclei
       integer,intent(in)                      ::  nat     !  Number of nodes
!
! Local variables
!
       real(kind=8),dimension(nat,nat)         ::  dadj    !  Double precision adjacency matrix
       real(kind=8),dimension(nat,nat)         ::  eigvec  !  Eigenvectors of the adjacency matrix
       real(kind=8),dimension(nat)             ::  eigval  !  Eigenvalues of the adjacency matrix
       real(kind=8),dimension(nat)             ::  vec     !  Eigenvector of greatest eigenvalue
       logical,dimension(nat)                  ::  check   !
       logical                                 ::  match   !
       integer,dimension(nat)                  ::  Order   !
       integer,dimension(nat)                  ::  ideg    !  Nodes degrees
       integer,dimension(nat)                  ::  ival    !  Integers associated eigenvalue centrality
       integer                                 ::  nval    !  Number of different eigenvalue centralities
       integer                                 ::  i,j     !  Indexes
!
       real(kind=8),parameter                  ::  thr = 1.0E-6
!
!   Chemically equivalent atoms strategy
!   ....................................
!
! Finding eigenvector centrality of the nodes
!
       dadj(:,:) = 0.0d0
       do j = 1, nat
         do i = 1, nat
           if ( adj(i,j) ) dadj(i,j) = 1.0d0
         end do
       end do
!
       call diagonalize(nat,dadj,eigval,eigvec)
!
       vec(:) = eigvec(:,nat)
!
! Creating new labels accounting for the eigenvector centrality
!
       call dnormlabels(nat,vec,ival,nval) 
!
!   Computing degrees
!
       ideg(:) = calcdegundir(nat,adj)
!
!   Combining eigenvector centrality, atom labels and degrees in a 
!    single label
!
       do i = 1, nat
         ideg(i) = ideg(i)*nlab*nval + ilab(i)*nval + ival(i)
       end do
!
! Finding equivalent nodes in the molecular graph
!
       check(:) = .TRUE.
       eqv(:,:) = .FALSE.
!
       do i = 1, nat-1
!
         check(i) = .FALSE.
!
!   Generate the order of the nodes of graph 1
!
!   Then, using a depth-first search, it produces an edge-visitation
!    list for Graph 1. 
!   A depth-first search on Graph 1 assigns labels to vertices based on
!    their order of visitation. This labeling, represented by the array 
!    ORDER, replaces the original vertex labeling
!
         call DFS(i,nat,adj,Order)  
!
         do j = i+1, nat
           if ( check(j) ) then
!
!   Checking (degree,label,ECV) sequence
!
             if ( ideg(i) .eq. ideg(j) ) then
!
!   Isomorphism testing
!
               match = VF2_Equivalent(nat,nat,adj,adj,ideg,ideg,order,i,j,check,eqv)
!~ write(*,*) 'checking',i,j,match
!
               if ( match ) then
                 check(j) = .FALSE.
                 eqv(i,j) = .TRUE.
                 eqv(j,i) = .TRUE.
               end if
!
             end if
           end if
         end do
!
       end do
!
       return
       end subroutine eqvatms
!
!======================================================================!
!
! NEARNEI - k-NEARest NEIghbors algorithm
!
! This subroutine 
!
       subroutine nearnei(nat,adj,eqv,nlab,ilab,knei,mindis,           &
                          r,mcycle,ncycle,cycles,debug)
!
       use isomorphism
       use sorting
!
       implicit none
!
! Input/output variables
!
       logical,dimension(nat,nat),intent(in)   ::  adj     !  Adjacency matrix
       logical,dimension(nat,nat),intent(out)  ::  eqv     !  Equivalent atoms relationships
       integer,dimension(nat,nat),intent(in)   ::  mindis  !  Minimum distance matrix
       integer,dimension(nat,r),intent(in)     ::  cycles  !  Cycles information
       integer,dimension(r),intent(in)         ::  ncycle  !  Number of atoms in each cycle
       integer,dimension(nat),intent(in)       ::  ilab    !  Integers associated to atomic labels
       integer,intent(in)                      ::  nlab    !  Number of different nuclei
       integer,intent(in)                      ::  mcycle  !  Number of cycles
       integer,intent(in)                      ::  knei    !  Maximum distance from source node
       integer,intent(in)                      ::  r       !  Cyle rank
       integer,intent(in)                      ::  nat     !  Number of nodes
       logical,intent(in)                      ::  debug   !  Debug mode
!
! Local variables
!
       real(kind=8),dimension(nat,nat)         ::  dval    !
       real(kind=8),dimension(nat,nat)         ::  dadj    !  Double precision adjacency matrix
       real(kind=8),dimension(nat,nat)         ::  eigvec  !  Eigenvectors of the adjacency matrix
       real(kind=8),dimension(nat)             ::  eigval  !  Eigenvalues of the adjacency matrix
       logical,dimension(nat,nat)              ::  adj1    !  Adjacency matrix of reference subgraph
       logical,dimension(nat,nat)              ::  adj2    !  Adjacency matrix of target subgraph
       logical,dimension(nat)                  ::  check   !  Nodes already processed
       logical,dimension(nat)                  ::  addcycle
       logical                                 ::  match   !
       integer,dimension(nat,nat)              ::  ideg    !  Subgraphs degrees
       integer,dimension(nat,nat)              ::  inei    !  Indexes of the atoms in each subgraph
       integer,dimension(nat)                  ::  Order   !
       integer,dimension(nat)                  ::  nnei    !  Number of atoms in each subgraph
       integer                                 ::  isubg   !  Subgraph index
       integer                                 ::  jsubg   !  Subgraph index
       integer                                 ::  icycle  !  Cycle index
       integer                                 ::  isize   !  Atom-in-cycle index
       integer                                 ::  n       !  Subgraph size
       integer                                 ::  i,j,k   !  Indexes
!
       real(kind=8),parameter                  ::  thr = 1.0E-6
!
!   K-Nearest neighbors algorithm
!   .............................
!
! Generating BFS-based subgraphs
!
       inei(:,:) = 0
       nnei(:)   = 1
!
       do j = 1, nat
         inei(1,j) = j
       end do
!
       do j = 1, nat-1
         do i = j+1, nat
           if ( mindis(i,j) .le. knei ) then
!
             nnei(j) = nnei(j) + 1
             nnei(i) = nnei(i) + 1
!
             inei(nnei(j),j) = i
             inei(nnei(i),i) = j
!
           end if
         end do
       end do
!
       if ( debug ) then
         write(*,*) 'Subgraphs information'
         write(*,*) '---------------------'
         write(*,*)
         write(*,'(1X,A)') 'Number of vertices:'
         write(*,'(200(1X,I3))') nnei(:)
         write(*,'(1X,20("-"))') 
         do i = 1, nat
           write(*,'(200(1X,I3))') (inei(i,j),j=1,nat)
         end do
         write(*,*)
       end if
!
! If any atom in the subgraph is present in a cycle then add all atoms
!  belonging to such cycle
!
       if ( mcycle .ne. 0 ) then
!
         if ( debug ) then
           write(*,*) 'Cycle information'
           write(*,*) '-----------------'
           do i = 1, mcycle
             write(*,*) 'Cycle Number', i, ':',(cycles(i,j),j=1,ncycle(i))
             write(*,*)
           end do
         end if
!
         do i = 1, nat  ! Loop over subgraphs
!
! Initializing the state of the nodes belonging to the subgraph
!
           check(:) = .FALSE.
           do j = 1, nnei(i)
             check(inei(j,i)) = .TRUE.
           end do
!
! Check if any atom in the subgraph belongs to a cycle
!
           addcycle(:) = .FALSE.
           do j = 1, nnei(i)  ! Loop over nodes in the subgraph i
             do icycle = 1, mcycle
               if ( .NOT. addcycle(icycle) ) then
                 do isize = 1, ncycle(icycle) 
                   if ( cycles(icycle,isize) .eq. inei(j,i) ) then
!~ write(*,*) inei(j,i),'is vertex',j,'in subgraph',i,'is in cycle',icycle
                     addcycle(icycle) = .TRUE.
!~ write(*,*) 'addcycle',addcycle(:)
                     exit
                   end if
                 end do     
               end if 
             end do      
           end do
!
           do icycle = 1, mcycle
             if ( addcycle(icycle) ) then
               do isize = 1, ncycle(icycle)
                 if ( .NOT. check(cycles(icycle,isize)) ) then
                   nnei(i) = nnei(i) + 1
                   inei(nnei(i),i) = cycles(icycle,isize)
                   check(cycles(icycle,isize)) = .TRUE.
                 end if
               end do
             end if
           end do
!
!~        write(*,*) 'Updated subgraphs'
!~        write(*,*) '-----------------'
!~        write(*,*)
!~        write(*,'(1X,A)') 'Number of vertices:'
!~        write(*,'(200(1X,I3))') nnei(:)
!~        write(*,'(1X,20("-"))') 
!~        do k = 1, nat
!~          write(*,'(200(1X,I3))') (inei(k,j),j=1,nat)
!~        end do
!~        write(*,*)
         end do
!
       end if
!
       if ( debug ) then
         write(*,*) 'Subgraphs information including cycles'
         write(*,*) '--------------------------------------'
         write(*,*)
         write(*,'(1X,A)') 'Number of vertices:'
         write(*,'(200(1X,I3))') nnei(:)
         write(*,'(1X,20("-"))') 
         do i = 1, nat
           write(*,'(200(1X,I3))') (inei(i,j),j=1,nat)
         end do
         write(*,*)
       end if
!
! Computing graph invariants of the generated subgraphs
!
!  Computing degrees
!
       do i = 1, nat  ! Loop over subgraphs
         do j = 1, nnei(i)
           ideg(j,i) = 0
           do k = 1, nnei(i)
             if ( adj(inei(j,i),inei(k,i)) ) ideg(j,i) = ideg(j,i) + 1
           end do
         end do            
       end do
!
!~        write(*,*) 'Degrees information'
!~        write(*,*) '-------------------'
!~        do i = 1, nat
!~          write(*,'(20(1X,I10))') (ideg(i,j),j=1,nat)
!~        end do
!~        write(*,*)
!
!  Computing eigenvector centrality of the nodes
!
!~        ival(:,:) = 0
       dval(:,:) = 0.0d0
!
       do isubg = 1, nat
!
         n = nnei(isubg)
!
         dadj(:,:) = 0.0d0
         do i = 1, n
           do j = 1, n
             if ( adj(inei(i,isubg),inei(j,isubg)) ) dadj(j,i) = 1.0d0
           end do 
         end do
!
         call diagonalize(n,dadj(:n,:n),eigval(:n),eigvec(:n,:n))
!~          vec(:n) = eigvec(:n,n)
         dval(:n,isubg) = abs(eigvec(:n,n))
!
!~          nval(isubg) = 1
!
!~          do i = 1, n
!~            match = .FALSE.
!~            do j = 1, i-1
!~              if ( abs(vec(i)-vec(j)) .lt. thr ) then
!~                match = .TRUE.
!~                ival(i,isubg) = ival(j,isubg)
!~                exit
!~              end if
!~            end do
!~            if ( .not. match ) then
!~              ival(i,isubg) = nval(isubg)
!~              nval(isubg)   = nval(isubg) + 1
!~            end if
!~          end do

!
       end do
!~ !
!~ ! Combining atom labels and degrees in a single variable
!~ !
!~        do i = 1, nat
!~          do j = 1, nnei(i)
!~            ideg(j,i) = ideg(j,i)*nlab*nval(i) + ilab(inei(j,i))*nval(i) + ival(j,i)
!~          end do
!~        end do
!
! Combining atom labels and degrees in a single variable
!
       do i = 1, nat
         do j = 1, nnei(i)
           ideg(j,i) = ideg(j,i)*nlab + ilab(inei(j,i))
         end do
       end do
!
!       write(*,*) 'Unique integer identifier'
!       write(*,*) '-------------------------'
!       do i = 1, nat
!         write(*,'(20(1X,I10))') (ideg(i,j),j=1,nat)
!       end do
!       write(*,*)
!
! Sorting subgraphs nodes according to the unique integer identifiers
!
       do i = 1, nat
!~          call ivqsort(nnei(i),ideg(:,i),inei(:,i),2,nnei(i))
         call ivdvqsort(nnei(i),ideg(:,i),inei(:,i),dval(:,i),2,nnei(i)) ! TODO: sorting must be based on dval
       end do
!
!~        write(*,*) 'Sorted subgraphs information'
!~        write(*,*) '--------.-------------------'
!~        write(*,*)
!~        write(*,'(1X,A,20(1X,I3))') 'Number of vertices:',nnei
!~        do i = 1, nat
!~          write(*,'(20(1X,I10))') (inei(i,j),j=1,nat)
!~        end do
!~        write(*,*)
!~ !
!~        write(*,*) 'Sorted unique integer identifier'
!~        write(*,*) '--------------------------------'
!~        do i = 1, nat
!~          write(*,'(20(1X,I10))') (ideg(i,j),j=1,nat)
!~        end do
!~        write(*,*)
!
! Solving graph isomorphism problem directly among subgraphs
!
       check(:) = .TRUE.
       eqv(:,:) = .FALSE.
!
       do isubg = 1, nat-1
!
         n = nnei(isubg)
         check(isubg) = .FALSE.
!
!   Building adjacency matrix of the reference subgraph
!
         adj1(:,:) = .FALSE.
         do i = 1, n-1
           do j = i+1, n
             if ( adj(inei(i,isubg),inei(j,isubg)) ) then
               adj1(i,j) =  .TRUE.
               adj1(j,i) =  .TRUE.
             end if
           end do 
         end do
!
! Generar el orden de los nodos del grafo 1
!
!   Then, using a depth-first search, it produces an edge-visitation
!    list for Graph 1. 
!   A depth-first search on Graph 1 assigns labels to vertices based on
!    their order of visitation. This labeling, represented by the array 
!    NUM, replaces the original vertex labeling
!
         call DFS(1,n,adj1(:n,:n),Order(:n))
!
         do jsubg = isubg+1, nat
           if ( check(jsubg) ) then
!~              if ( n .eq. nnei(jsubg) ) then
             if ( (n.eq.nnei(jsubg)) .and.                             &
                       (abs(dval(1,isubg)-dval(1,jsubg)).lt.thr)  ) then
!~ write(*,*) 'checking',isubg,jsubg,dval(1,isubg),dval(1,jsubg), &
!~ dval(1,isubg)-dval(1,jsubg),1.0E-6,((dval(1,isubg)-dval(1,jsubg)).lt.thr)
!
!   Prechecking (degree,label) sequence
!
               match = .TRUE.
               do i = 1, n
!~ write(*,*) 'checking',isubg,jsubg,i,dval(i,isubg),dval(i,jsubg),(dval(i,isubg)-dval(i,jsubg))
                 if ( ideg(i,isubg) .ne. ideg(i,jsubg) ) then  ! TODO: after sorting dval check dval
!~                  if ( (ideg(i,isubg).ne.ideg(i,jsubg)) .or.            &
!~                         (.NOT.(dval(i,isubg)-dval(i,jsubg)).lt.1.0E-6 ) ) then
!
                   match = .FALSE.
                   exit
                 end if
               end do
!
               if ( .NOT. match ) cycle
!
!   Building adjacency matrix of the target subgraph
!
               adj2(:,:) = .FALSE.
               do i = 1, n-1
                 do j = i+1, n
                   if ( adj(inei(i,jsubg),inei(j,jsubg)) ) then
                     adj2(i,j) =  .TRUE.
                     adj2(j,i) =  .TRUE.
                   end if
                 end do 
               end do
!
!   Isomorphism testing
!
               match = VF2_EquivalentEVC(n,n,adj1(:n,:n),adj2(:n,:n),  &
                                      ideg(:n,isubg),ideg(:n,jsubg),   &
                                      dval(:n,isubg),dval(:n,jsubg),   &
                                      order(:n),1,1)
!
               if ( match ) then
                 check(jsubg)     = .FALSE.
                 eqv(isubg,jsubg) = .TRUE.
                 eqv(jsubg,isubg) = .TRUE.
               end if
!
             end if
           end if
         end do
!
       end do
!
       return
       end subroutine nearnei
!
!======================================================================!
!
! GENTYPES - GENerate atom TYPES
!
! This subroutine 
!
       subroutine gentypes(nat,eqv,lab,newlab,ilab,mlab,eqlab,idat,    &
                           nidat,debug)
!
       use printings
       use graphtools
!
       implicit none
!
! Input/output variables
!
       character(len=lenlab),dimension(nat),intent(in)   ::  lab      !
       character(len=lenlab),dimension(nat),intent(out)  ::  newlab   !
       integer,dimension(nat),intent(in)                 ::  ilab     !  Integers associated to atomic labels
       integer,dimension(nat),intent(out)                ::  eqlab    !  Integers associated to each type of nuclei
       integer,dimension(nat),intent(out)                ::  idat     !  Integers associated to atomtypes
       integer,intent(in)                                ::  nat      !  Number of atoms
       integer,intent(in)                                ::  mlab     !  Number of different nuclei
       integer,intent(out)                               ::  nidat    !  Number of different atomtypes
       logical,dimension(nat,nat),intent(in)             ::  eqv      !  Equivalent atoms relationships
       logical,intent(in)                                ::  debug    !  Debug mode
!
! Equivalent atoms information
!
       integer,dimension(nat)                            ::  mol      !
       integer,dimension(nat)                            ::  agg      !
       integer,dimension(nat)                            ::  tag      !
       integer,dimension(nat)                            ::  imol     !
       integer,dimension(nat)                            ::  iagg     !
       integer,dimension(nat)                            ::  itag     !
       integer,dimension(nat)                            ::  nmol     !
       integer,dimension(nat)                            ::  nagg     !
       integer,dimension(nat)                            ::  ntag     !
       integer                                           ::  magg     !  Number of aggregates
       integer                                           ::  nsize    !  Maximum aggregate size
!
! Local variables
! 
       character(len=20)                                 ::  str      !
       integer,dimension(mlab)                           ::  nlab     !  Number of types of each different nuclei
       integer,dimension(nat)                            ::  ilabaux  !
       integer                                           ::  i,j,k    !  Indexes
       integer                                           ::  ii,kk    !  Indexes
! 
! Block diagonalizing equivalence matrix
!
       call blockdiag(nat,eqv,mol,tag,agg,nsize,nagg,iagg,             &
                      nmol,imol,magg)
!
! Generating new atomtypes
!
       nlab(:)  = 0
       eqlab(:) = 0
!
       do i = 1, nsize  !  TODO: check
         k = imol(i)
         do j = 1, nagg(i)
           nlab(ilab(mol(k+1))) = nlab(ilab(mol(k+1))) + 1
           kk = nlab(ilab(mol(k+1)))
           do ii = 1, i
             eqlab(mol(k+ii)) = kk
           end do
           k = k + i
         end do
       end do
!
! Combining atom-name label with equivalent-atom label 
!
       do i = 1, nat
         ilabaux(i) = mlab*eqlab(i) + ilab(i)
       end do
!
       call inormlabels(nat,ilabaux,idat,nidat)
!
! Generating new atom labels
!
       do i = 1, nat
         write(str,*) eqlab(i)
         str       = adjustl(str)
         newlab(i) = trim(lab(i))//trim(str)
       end do
!
! Printing atom typing information
!
       if ( debug ) then
         write(*,*) 'Equivalence matrix information'
         write(*,*) '------------------------------'
         write(*,'(2X,A,X,I6)')     'Total number of entities : ',magg
         write(*,*)
         write(*,'(2X,A,20(X,I6))') 'Aggregates of each type  : ',     &
                                                            nagg(:nsize)
         write(*,'(2X,A,20(X,I6))') '  Accumulation           : ',     &
                                                            iagg(:nsize)
         write(*,*)
         write(*,'(2X,A,20(X,I6))') 'Molecules of each type   : ',     &
                                                            nmol(:nsize)
         write(*,'(2X,A,20(X,I6))') '  Accumulation           : ',     &
                                                            imol(:nsize)
         write(*,*)
!
         call print_info(0,nat,agg,tag,mol,'agg','tag','mol')
       end if
!
       return
       end subroutine gentypes
!
!======================================================================!
!
! GENFFBONDED - GENerate Force Field BONDED interactions
!
! This subroutine 
!
       subroutine genffbonded(nat,idat,nidat,coord,adj,ideg,lcycle,    &
                              lrigid,znum,bonded,dihed,iroute,debug)
!
       use datatypes, only: grobonded,                                 &
                            dihedrals
!
       implicit none
!
! Input/output variables
!
       type(grobonded),intent(inout)                      ::  bonded    !
       type(dihedrals),intent(inout)                      ::  dihed     !
       real(kind=8),dimension(3,nat),intent(in)           ::  coord     !  Atomic coordinates
       integer,dimension(nat),intent(in)                  ::  idat      ! 
       integer,dimension(nat),intent(in)                  ::  ideg      ! 
       logical,dimension(nat,nat),intent(in)              ::  adj       !  Boolean adjacency 
       logical,dimension(nat,nat),intent(in)              ::  lcycle    !  Bonds belonging to rings
       logical,dimension(nat,nat),intent(in)              ::  lrigid    ! 
!
       integer,dimension(nat),intent(in)                  ::  znum      !  
       integer,intent(in)                                 ::  nidat     !
       integer,intent(in)                                 ::  nat       !  Number of atoms
!
       integer,intent(in)                                 ::  iroute    !
       logical,intent(in)                                 ::  debug     !  Debug mode
!
! Local variables
!  
       logical,dimension(:,:),allocatable                 ::  adjbond   !
       logical,dimension(:,:),allocatable                 ::  adjang    !
       integer,dimension(:,:),allocatable                 ::  edgeang   !
!
       integer                                            ::  iterm     !
!
! Generating bonded terms
! -----------------------
!
      call genbonds(nat,coord,adj,bonded%nbond,bonded%ibond,           &
                    bonded%bond,bonded%kbond,bonded%fbond,             &
                    adjbond,debug)
!
      if ( nat .gt. 2 ) then
        call genangles(nat,coord,bonded%nbond,bonded%ibond,adjbond,    &
                       bonded%nang,bonded%iang,bonded%ang,bonded%kang, &
                       bonded%fang,adjang,edgeang,debug)
       end if
!
       if ( nat .gt. 3 ) then
         call gendihe(nat,idat,nidat,coord,znum,ideg,lcycle,lrigid,    &
                      bonded%nbond,bonded%ibond,bonded%nang,           &
                      bonded%iang,adjang,edgeang,dihed,iroute,debug)
       end if
!
       deallocate(adjbond,adjang,edgeang)
!
       return
       end subroutine genffbonded
!
!======================================================================!
!
! GENBONDS - GENerate force field BONDS
!
! This subroutine 
!
       subroutine genbonds(nat,coord,adj,nbond,ibond,bond,kbond,fbond, &
                           adjbond,debug)
!
       use graphtools
!
       implicit none
!
! Input/output variables
!
       real(kind=8),dimension(:),allocatable,intent(out)  ::  bond     !  Equilibrium bond distance
       real(kind=8),dimension(:),allocatable,intent(out)  ::  kbond    !  Stretchings force constant
       integer,dimension(:,:),allocatable,intent(out)     ::  ibond    !  Bond index
       integer,dimension(:),allocatable,intent(out)       ::  fbond    !  Bond type
       integer,intent(out)                                ::  nbond    !  Number of bonds
!
       logical,dimension(nat,nat),intent(in)              ::  adj       !  Boolean adjacency 
       logical,dimension(:,:),allocatable,intent(out)     ::  adjbond   !  Boolean adjacency 
       real(kind=8),dimension(3,nat),intent(in)           ::  coord     !  Atomic coordinates
       integer,intent(in)                                 ::  nat       !  Number of atoms
!
       logical,intent(in)                                 ::  debug     !  Debug mode
!
! Local variables
!
       integer                                            ::  i
!
! Generating bond-stretching terms
! --------------------------------
!
!  Obtaining bonds information from the molecule (atoms representation)
!
       if ( debug ) then
         write(*,*) 'GENERATING BOND-STRETCHING TERMS'
         write(*,*) '................................'
         write(*,*) 
       end if
!
       call countedges(nat,adj,nbond)
!
       allocate(ibond(2,nbond),bond(nbond),kbond(nbond),fbond(nbond))
       allocate(adjbond(nbond,nbond))
!
       call adj2edge(nat,adj,nbond,ibond)
!
!  Changing from atoms to bonds representation
!
       call changerep(nat,adj,nbond,adjbond,ibond,debug)
!
! Computing bond-stretching equilibrium terms
!
       do i = 1, nbond
         bond(i) = sqrt((coord(1,ibond(1,i))-coord(1,ibond(2,i)))**2   &
                        +(coord(2,ibond(1,i))-coord(2,ibond(2,i)))**2  & 
                   +(coord(3,ibond(1,i))-coord(3,ibond(2,i)))**2)/10.0d0
       end do
!
       kbond(:) = 0.0d0
       fbond(:) = 1
!
       return
       end subroutine genbonds
!
!======================================================================!
!
! GENANGLES - GENerate force field ANGLES
!
! This subroutine 
!
       subroutine genangles(nat,coord,nbond,ibond,adjbond,nang,iang,   &
                            ang,kang,fang,adjang,edgeang,debug)
!
       use graphtools
       use genfftools
!
       implicit none
!
! Input/output variables
!
       real(kind=8),dimension(:),allocatable,intent(out)  ::  ang      !  Equilibrium angle  
       real(kind=8),dimension(:),allocatable,intent(out)  ::  kang     !  Bending force constant
       integer,dimension(:,:),allocatable,intent(out)     ::  iang     !  Angle index
       integer,dimension(:),allocatable,intent(out)       ::  fang     !  Angle type
       integer,intent(out)                                ::  nang     !  Angles number
!
       logical,dimension(nbond,nbond),intent(in)          ::  adjbond  !  Boolean adjacency 
       logical,dimension(:,:),allocatable,intent(out)     ::  adjang   !  Boolean adjacency 
       real(kind=8),dimension(3,nat),intent(in)           ::  coord    !  Atomic coordinates
       integer,dimension(:,:),allocatable,intent(out)     ::  edgeang  !
       integer,dimension(2,nbond),intent(in)              ::  ibond    !
       integer,intent(in)                                 ::  nbond    !  Number of bonds
       integer,intent(in)                                 ::  nat      !  Number of atoms
!
       logical,intent(in)                                 ::  debug    !  Debug mode
!
! Local variables
!
       integer                                            ::  i
!
       real(kind=8),parameter                             ::  pi =  4*atan(1.0_8) 
!
! Generating angle-bending terms
! --------------------------------
!
!  Obtaining angles information from the molecule (bonds representation)
!
       if ( debug ) then
         write(*,*) 'GENERATING ANGLE-BENDING TERMS'
         write(*,*) '..............................'
         write(*,*) 
       end if
!
       call countedges(nbond,adjbond,nang)
!
       allocate(iang(3,nang),ang(nang),kang(nang),fang(nang))
       allocate(edgeang(2,nang),adjang(nang,nang))
!
       call adj2edge(nbond,adjbond,nang,edgeang)
!
!  Changing from bonds to angles representation
!
       call changerep(nbond,adjbond,nang,adjang,edgeang,debug)
!
       call setang(nbond,ibond,nang,edgeang,iang)
!
! Computing angle-bending equilibrium terms
!
       do i = 1, nang
         ang(i) = calc_angle(coord(:,iang(1,i)),coord(:,iang(2,i)),    &
                                                     coord(:,iang(3,i)))
         ang(i) = ang(i)*180.0d0/pi
       end do
!
       kang(:) = 0.0d0
       fang(:) = 1
!
       return
       end subroutine genangles
!
!======================================================================!
!
! GENDIED - GENerate force field DIHEDrals
!
! This subroutine 
!
       subroutine gendihe(nat,idat,nidat,coord,znum,ideg,lcycle,       &
                          lrigid,nbond,ibond,nang,iang,adjang,         &
                          edgeang,dihed,iroute,debug)
!
       use datatypes, only: dihedrals
       use graphtools
!
       implicit none
!
! Input/output variables
!
       type(dihedrals),intent(out)                        ::  dihed    !  Equilibrium dihedral angle
       real(kind=8),dimension(3,nat),intent(in)           ::  coord    !  Atomic coordinates
       logical,dimension(nat,nat),intent(in)              ::  lcycle   !  Bonds belonging to rings
       logical,dimension(nat,nat),intent(in)              ::  lrigid   !  Rigid bonds
       integer,dimension(nat),intent(in)                  ::  idat     ! 
       integer,dimension(nat),intent(in)                  ::  znum     ! 
       integer,dimension(nat),intent(in)                  ::  ideg     ! 
       integer,intent(in)                                 ::  nat      !  Number of atoms
       integer,intent(in)                                 ::  nidat    !
!
       logical,dimension(nang,nang),intent(in)            ::  adjang   !  Boolean adjacency 
       integer,dimension(2,nang),intent(in)               ::  edgeang  !
       integer,dimension(2,nbond),intent(in)              ::  ibond    !
       integer,dimension(3,nang),intent(in)               ::  iang     !
       integer,intent(in)                                 ::  nang     !  Number of angles
       integer,intent(in)                                 ::  nbond    !  Number of bonds
!
       integer,intent(in)                                 ::  iroute   !
       logical,intent(in)                                 ::  debug    !  Debug mode
!
! Local variables
!
       integer,dimension(:,:),allocatable                 ::  edges
       integer                                            ::  ndihe    !  Number of dihedrals
       integer                                            ::  i
!
       real(kind=8),parameter                             ::  pi =  4*atan(1.0_8) 
!
! Generating proper dihedral or torsional terms and improper dihedrals
! --------------------------------------------------------------------
!
!  Obtaining dihedral information from the molecule
!
       if ( debug ) then
         write(*,*) 'GENERATING DIHEDRAL TERMS'
         write(*,*) '.........................'
         write(*,*) 
       end if
!
       call countedges(nang,adjang,ndihe)
!
       dihed%ndihe = ndihe
!
       allocate(dihed%iimpro(4,ndihe),dihed%dimpro(ndihe),             &
                dihed%kimpro(ndihe),dihed%fimpro(ndihe))
!
       dihed%iimpro(:,:) = 0
       dihed%fimpro(:)   = 0
       dihed%dimpro(:)   = 0.0d0
       dihed%kimpro(:)   = 0.0d0
       dihed%nimpro      = 0
!
       allocate(dihed%iinv(4,ndihe),dihed%dinv(ndihe),                 &
                dihed%kinv(ndihe),dihed%finv(ndihe))
!
       dihed%iinv(:,:) = 0
       dihed%finv(:)   = 0
       dihed%dinv(:)   = 0.0d0
       dihed%kinv(:)   = 0.0d0
       dihed%ninv      = 0
!
       allocate(dihed%irigid(4,ndihe),dihed%drigid(ndihe),             &
                dihed%krigid(ndihe),dihed%frigid(ndihe))
       dihed%irigid(:,:) = 0
       dihed%frigid(:)   = 0
       dihed%drigid(:)   = 0.0d0
       dihed%krigid(:)   = 0.0d0
       dihed%nrigid      = 0
!
       allocate(dihed%iflexi(4,ndihe),dihed%dflexi(ndihe),             &
                dihed%kflexi(ndihe),dihed%flexi(ndihe),dihed%fflexi(ndihe))
       dihed%iflexi(:,:) = 0
       dihed%fflexi(:)   = 0
       dihed%dflexi(:)   = 0.0d0
       dihed%kflexi(:)   = 0.0d0
       dihed%nflexi      = 0
!
       allocate(edges(2,ndihe))
!
       call adj2edge(nang,adjang,ndihe,edges)
!
! Classifying adjacent angles as proper or improper dihedrals
!
       call setdihe(nat,idat,nidat,coord,nbond,ibond,nang,edgeang,     &
                    iang,ndihe,dihed,edges,lcycle,lrigid,znum,ideg,    &
                    iroute)
!
! Removing 3-member cycle dihedrals  ! TODO
!

!
       deallocate(edges)
!
       return
       end subroutine gendihe
!
!======================================================================!
!
! SYMFFBONDED - SYMmetrize Force Field BONDED interactions
!
! This subroutine 
!
       subroutine symffbonded(nat,nidat,idat,newlab,bonded,dihed,      &
                              r,marunit,narunit,arunit,coord,          &
                              lheavy,adj,fsymm,debug)
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
       real(kind=8),dimension(3,nat),intent(in)         ::  coord     !
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
       logical,intent(in)                               ::  debug     !  Debug mode
!
! Local variables
!  
       character(len=20)                                ::  cnum      !
       real(kind=8)                                     ::  daux      !
       integer,dimension(:),allocatable                 ::  ivaux     !
       integer,dimension(4)                             ::  vaux      !
       integer                                          ::  iaux      !
       integer                                          ::  iterm     !
       integer                                          ::  itmp      !
       integer                                          ::  i,j       !
!
       real(kind=8),parameter                           ::  pi =  4*atan(1.0_8) 
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
                    bonded%sbond,bonded%labbond,fsymm,debug)
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
                    bonded%labang,fsymm,debug)
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
       do i = 1, dihed%nflexi
         do j = 1, dihed%flexi(i)%ntor
           iterm = iterm + 1
           dihed%flexi(i)%tor(j)%stor = iterm
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
                      marunit,narunit,arunit,adj,lheavy,fsymm,debug)
       end if
!
       if ( dihed%ninv .gt. 0 ) then
         call symoop(dihed%ninv,dihed%idinv,dihed%dinv,dihed%sinv,     &
                     dihed%labinv,fsymm,debug)
       end if 
!
       if ( dihed%nflexi .gt. 0 ) then
         call symtor(dihed%nflexi*10,dihed)
       end if
!
       deallocate(ivaux)
!
       return
       end subroutine symffbonded
!
!======================================================================!
!
! CHANGEREP - CHANGE REPresentation
!
! This subroutine 
!
       subroutine changerep(nat,adj,nbond,adjbond,edge,debug)
!
       use sorting,   only : iinsertdesc
       use graphtools
!
       implicit none
!
! Input/output variables
!
       logical,dimension(nat,nat),intent(in)       ::  adj       !  Boolean adjacency matrix
       logical,dimension(nbond,nbond),intent(out)  ::  adjbond   !
       integer,dimension(2,nbond),intent(in)       ::  edge      !
       integer,intent(in)                          ::  nat       !  Number of atoms
       integer,intent(in)                          ::  nbond     !
       logical,intent(in)                          ::  debug     !  Debug mode
!
! Local variables
!  
       logical,dimension(nbond,nat)                ::  incbond   !
       integer                                     ::  i,j       !  Indexes
!
! Generating bond-stretching terms
! --------------------------------
!
!  Obtaining incidence information from the adjacency matrix
!
       call adj2inc(nat,adj,nbond,incbond)
!
       if ( debug ) then
         write(*,*) 'Incidence matrix'
         write(*,*) '----------------'
         write(*,'(4X,200(1X,I3))') (j,j=1,nat)
         do i = 1, nbond
           write(*,'(1X,I3,200(1X,L3))') i,(incbond(i,j),j=1,nat)
         end do
         write(*,*) 
       end if     
!
!  Building adjacency matrix in the new representation
!
       call inc2adj(nat,nbond,incbond,edge,adjbond)
!
       if ( debug ) then
         write(*,*) 'Adjacency matrix in the new representation'
         write(*,*) '------------------------------------------'
         write(*,'(4X,200(1X,I3))') (j,j=1,nbond)
         do i = 1, nbond
           write(*,'(1X,I3,200(1X,L3))') i,(adjbond(i,j),j=1,nbond)
         end do
         write(*,*)
       end if
!
       return
       end subroutine changerep
!
!======================================================================!
!
! SETANG - SET ANGles
!
! This subroutine 
!
       subroutine setang(nbond,ibond,nang,edgeang,iang)
!
       implicit none
!
! Input/output variables
!
       integer,dimension(2,nbond),intent(in)  ::  ibond    !  
       integer,dimension(2,nang),intent(in)   ::  edgeang  !  
       integer,dimension(3,nang),intent(out)  ::  iang     !  
       integer,intent(in)                     ::  nbond    !
       integer,intent(in)                     ::  nang     !  
!
! Local variables
!  
       integer,dimension(2)                   ::  bond1    !
       integer,dimension(2)                   ::  bond2    !
       integer                                ::  i,j,k    !  Indexes
!
! Finding atoms belonging to each angle
! -------------------------------------
!
! Finding atoms belonging to each angle
!
       do i = 1, nang
!
         bond1(:) = ibond(:,edgeang(1,i))
         bond2(:) = ibond(:,edgeang(2,i))
!
         do j = 1, 2
           do k = 1, 2
             if ( bond1(j) .eq. bond2(k) ) then
!
               iang(2,i) = bond1(j)
!
               if ( bond1(modulo(j,2)+1) .lt. bond2(modulo(k,2)+1) )   & ! TODO: idat(bond1(modulo(j,2)+1))
                                                                    then
                 iang(1,i) = bond1(modulo(j,2)+1)
                 iang(3,i) = bond2(modulo(k,2)+1)
               else
                 iang(1,i) = bond2(modulo(k,2)+1)
                 iang(3,i) = bond1(modulo(j,2)+1)
               end if
!
               GO TO 1000
!
             end if
           end do
         end do
1000     continue
!
       end do
!
       return
       end subroutine setang
!
!======================================================================!
!
! SETDIHE - SET DIHEdrals
!
! This subroutine 
!
       subroutine setdihe(nat,idat,nidat,coord,nbond,ibond,nang,       &
                          edgeang,iang,ndihe,dihed,edgedihe,lcycle,    &
                          lrigid,znum,ideg,iroute)
!
       use datatypes,   only: dihedrals
       use genfftools
!
       implicit none
!
! Input/output variables
!
       type(dihedrals),intent(inout)              ::  dihed     !                            
       real(kind=8),dimension(3,nat),intent(in)   ::  coord     !
       logical,dimension(nat,nat),intent(in)      ::  lcycle    !  Bonds belonging to rings
       logical,dimension(nat,nat),intent(in)      ::  lrigid    !  Bonds belonging to rings
       integer,dimension(nat),intent(in)          ::  idat      ! 
       integer,dimension(nat),intent(in)          ::  znum      ! 
       integer,dimension(nat),intent(in)          ::  ideg      ! 
       integer,intent(in)                         ::  nidat     !
       integer,intent(in)                         ::  nat       !
!
       integer,dimension(2,nbond),intent(in)      ::  ibond     !  
       integer,intent(in)                         ::  nbond     !
!
       integer,dimension(2,nang),intent(in)       ::  edgeang   !  
       integer,dimension(3,nang),intent(in)       ::  iang      !  
       integer,intent(in)                         ::  nang      !
!
       integer,dimension(2,ndihe),intent(in)      ::  edgedihe  !  
       integer,intent(in)                         ::  ndihe     ! 
!
       integer,intent(in)                         ::  iroute    ! 
!
! Local variables
!  
       type(dihedrals)                            ::  tmpdi     !                            
       real(kind=8),dimension(ndihe)              ::  dimpro    !
       integer,dimension(4,ndihe)                 ::  iimpro    !
       integer,dimension(ndihe)                   ::  fimpro    !
       integer                                    ::  nimpro    !
       real(kind=8),dimension(ndihe)              ::  dinv      !
       integer,dimension(4,ndihe)                 ::  iinv      !
       integer,dimension(ndihe)                   ::  finv      !
       integer                                    ::  ninv      !
       logical,dimension(ndihe)                   ::  torsion   !
       logical,dimension(ndihe)                   ::  visited   !
       logical,dimension(4)                       ::  lcheck    !
       logical                                    ::  flag      !
       logical                                    ::  lfound    !
       real(kind=8)                               ::  daux      !
       real(kind=8)                               ::  daux1     !
       real(kind=8)                               ::  daux2     !
       integer,dimension(4,ndihe)                 ::  auxdihe   !  
       integer,dimension(ndihe)                   ::  neqquad   !  
       integer                                    ::  meqquad   !  
       integer,dimension(ndihe)                   ::  imap      !  
       integer,dimension(ndihe)                   ::  iflexi    !  
       integer,dimension(ndihe)                   ::  ivaux     !  
       integer,dimension(4)                       ::  vaux      !  
       integer,dimension(4)                       ::  vaux1     !  
       integer,dimension(4)                       ::  vaux2     !  
       integer,dimension(2,4)                     ::  bonds     !
       integer,dimension(2)                       ::  rbond     !
       integer,dimension(2)                       ::  bond1     !
       integer,dimension(2)                       ::  bond2     !
       integer,dimension(2)                       ::  bond3     !
       integer,dimension(2)                       ::  bond4     !
       integer,dimension(3)                       ::  ang1      !
       integer,dimension(3)                       ::  ang2      !
       integer                                    ::  nmap      !
       integer                                    ::  id1       !
       integer                                    ::  id2       !
       integer                                    ::  id11      !
       integer                                    ::  id12      !
       integer                                    ::  id21      !
       integer                                    ::  id22      !
       integer                                    ::  ntor      !
       integer                                    ::  iaux      !
       integer                                    ::  iaux1     !
       integer                                    ::  iaux2     !
       integer                                    ::  iaux3     !
       integer                                    ::  iaux4     !
       integer                                    ::  itmp      !
       integer                                    ::  i,j,k     !  Indexes
       integer                                    ::  ii,jj,kk  !  Indexes
!
       real(kind=8),parameter                    ::  pi =  4*atan(1.0_8) 
!
! Finding atoms belonging to each dihedral
! ----------------------------------------
!
       torsion(:) = .TRUE.
!
!  Finding improper out of plane dihedrals 
!  .......................................  
!
! If adjacent angles have the same central atom (geminal angles) then 
!  the dihedral is an improper out of plane tetrahedral angle (out of
!  plane angles)
!
       do i = 1, ndihe
!
         ang1(:) = iang(:,edgedihe(1,i))
         ang2(:) = iang(:,edgedihe(2,i))
!
         if ( ang1(2) .eq. ang2(2) ) then
!
           torsion(i) = .FALSE.
!
           vaux(1) = ang1(2)
!
           if ( ang1(1) .eq. ang2(1) ) then !> sorting k,l atoms by canonical order
             vaux(2) = ang1(1)
             call checklab(vaux,ang1(3),ang2(3))
           else if ( ang1(1) .eq. ang2(3) ) then
             vaux(2) = ang1(1)
             call checklab(vaux,ang1(3),ang2(1))
           else if ( ang1(3) .eq. ang2(1) ) then
             vaux(2) = ang1(3)
             call checklab(vaux,ang1(1),ang2(3))
           else if ( ang1(3) .eq. ang2(3) ) then
             vaux(2) = ang1(3)
             call checklab(vaux,ang1(1),ang2(1))
           end if
!
! Computing equilibrium value
!
           if ( idat(vaux(3)) .gt. idat(vaux(4)) ) then ! TODO: move to position 1 central atom
             itmp    = vaux(3)
             vaux(3) = vaux(4)
             vaux(4) = itmp           
           end if
!
           call Diedro(coord(:,vaux(1)),coord(:,vaux(2)),              & 
                       coord(:,vaux(3)),coord(:,vaux(4)),daux)
           daux = daux*180.0d0/pi
!
! If dihedral angle is planar classify it as improper
!
           if ( dabs(daux) .lt. 25.0d0 ) then
!            
             flag = .FALSE.
             if ( iroute .eq. 1 ) then      ! planar rings, ethene, ketone, etc.
               flag = .NOT. lcycle(vaux(1),vaux(2)) 
             else if ( iroute .eq. 2 ) then ! planar rings and double bonds
               flag = ( (.NOT.lcycle(vaux(1),vaux(2)))                 & ! double bonds
                               .and. (lrigid(vaux(1),vaux(2))) ) .or.  &
                        ( ((.NOT.lcycle(vaux(1),vaux(2))))             & ! planar rings
                                      .and. (lcycle(vaux(1),vaux(3)))  &
                                       .and. (lcycle(vaux(1),vaux(4))) )
             else if ( iroute .eq. 3 ) then ! only in planar rings
               flag = ((.NOT.lcycle(vaux(1),vaux(2))))                 & ! planar rings
                                      .and. (lcycle(vaux(1),vaux(3)))  &
                                         .and. (lcycle(vaux(1),vaux(4)))
             else if ( iroute .eq. 4 ) then ! only in double bonds
               flag = (.NOT.lcycle(vaux(1),vaux(2)))                   & ! double bonds
                                         .and. (lrigid(vaux(1),vaux(2)))
             else if ( iroute .eq. 4 ) then ! only in ketone, imine, etc.
               flag = (.NOT.lcycle(vaux(1),vaux(2)))                   &
                                        .and. lrigid(vaux(1),vaux(2))  &
                                              .and. (znum(vaux(2)).gt.6) 
             end if
!
             if ( flag ) then
!
               dihed%nimpro = dihed%nimpro + 1 
!
               dihed%iimpro(:,dihed%nimpro) = vaux(:)
               dihed%fimpro(dihed%nimpro)   = 2
               dihed%dimpro(dihed%nimpro)   = 0.0d0    ! TODO: option to pick equilibrium value or set to 0
!~                dihed%dimpro(dihed%nimpro)   = daux  ! TODO: option to pick equilibrium value or set to 0
!
             end if
!
! If dihedral angle deviates from planarity classify it as inversion  
!
           else if ( znum(vaux(1)) .eq. 7 ) then ! TODO: only inversion in amines
!
             dihed%ninv = dihed%ninv + 1
!
             dihed%iinv(:,dihed%ninv) = vaux(:)
!~              dihed%finv(dihed%ninv)   = 3  !  TODO: check best function for amines
             dihed%finv(dihed%ninv)   = 2     !  TODO: check best function for amines
             dihed%dinv(dihed%ninv)   = daux  !  TODO: how to handle equilibrium value
!
           end if
!
         end if
!
       end do
!
! Keep only one improper dihedral per central atom
!
       nimpro = 0
       do i = 1, dihed%nimpro
         flag = .FALSE.
         do j = 1, nimpro
           if ( dihed%iimpro(1,i) .eq. iimpro(1,j) ) then
             flag = .TRUE.
!
             vaux1(:) = dihed%iimpro(:,i)
             vaux2(:) = iimpro(:,j)
!
             id1 = idat(vaux1(1))*nidat**3 + idat(vaux1(2))*nidat**2   &
                   + idat(vaux1(3))*nidat + idat(vaux1(4))
!
             id2 = idat(vaux2(1))*nidat**3 + idat(vaux2(2))*nidat**2   &
                   + idat(vaux2(3))*nidat + idat(vaux2(4))     
!
             if ( id1 .gt. id2 ) then
               dihed%iimpro(:,i) = iimpro(:,nimpro)
               dihed%dimpro(i)   = dimpro(nimpro)
               dihed%fimpro(i)   = fimpro(nimpro) 
             end if   
!
             exit
           end if
         end do
         if ( .not. flag ) then
!
           nimpro = nimpro + 1
!
           iimpro(:,nimpro) = dihed%iimpro(:,i) 
           dimpro(nimpro)   = dihed%dimpro(i) 
           fimpro(nimpro)   = dihed%fimpro(i) 
!
         end if
       end do
!
       dihed%nimpro      = nimpro
       dihed%iimpro(:,:) = iimpro(:,:)
       dihed%dimpro(:)   = dimpro(:)
       dihed%fimpro(:)   = fimpro(:)
!
! Keep only one inversion dihedral per central atom
!
       ninv = 0
       do i = 1, dihed%ninv
         flag = .FALSE.
         do j = 1, ninv
           if ( dihed%iinv(1,i) .eq. iinv(1,j) ) then
             flag = .TRUE.
!
             vaux1(:) = dihed%iimpro(:,i)
             vaux2(:) = iimpro(:,i)
!
             id1 = idat(vaux1(1))*nidat**3 + idat(vaux1(2))*nidat**2   &
                   + idat(vaux1(3))*nidat + idat(vaux1(4))
!
             id2 = idat(vaux2(1))*nidat**3 + idat(vaux2(2))*nidat**2   &
                   + idat(vaux2(3))*nidat + idat(vaux2(4))     
!
             if ( id1 .gt. id2 ) then
               dihed%iimpro(:,i) = iimpro(:,nimpro)
               dihed%dimpro(i)   = dimpro(nimpro)
               dihed%fimpro(i)   = fimpro(nimpro) 
             end if   
!
             exit
           end if
         end do
         if ( .not. flag ) then
!
           ninv = ninv + 1
!
           iinv(:,nimpro) = dihed%iinv(:,i) 
           dinv(nimpro)   = dihed%dinv(i) 
           finv(nimpro)   = dihed%finv(i) 
!
         end if
       end do
!
       dihed%ninv      = ninv
       dihed%iinv(:,:) = iinv(:,:)
       dihed%dinv(:)   = dinv(:)
       dihed%finv(:)   = finv(:)
!
!  Finding torsional dihedrals
!  ...........................
!
       do i = 1, ndihe
         if ( torsion(i) ) then
!
           lcheck(:) = .TRUE.
!
           bonds(:,1) = ibond(:,edgeang(1,edgedihe(1,i)))
           bonds(:,2) = ibond(:,edgeang(2,edgedihe(1,i)))
           bonds(:,3) = ibond(:,edgeang(1,edgedihe(2,i)))
           bonds(:,4) = ibond(:,edgeang(2,edgedihe(2,i)))
!
           lfound = .FALSE.
! Finding two identical bonds
           do j = 1, 3
             do k = j+1, 4
               if ( (bonds(1,j).eq.bonds(1,k)) .and. ((bonds(2,j).eq.bonds(2,k))) ) then
                 lcheck(j) = .FALSE.
                 lcheck(k) = .FALSE.
!
                 rbond(1)  = bonds(1,j)
                 rbond(2)  = bonds(2,j)
!
                 vaux(2) = rbond(1)
                 vaux(3) = rbond(2)
!                 
                 lfound = .TRUE.
                 exit
               end if
             end do
             if ( lfound ) exit
           end do
!
! FIXME: if (.NOT.lfound) there was a problem
!
! Finding matching indexes in the other 2 bonds
           do j = 1, 4
             if ( lcheck(j) ) then
               if ( bonds(1,j) .eq. rbond(1) ) then
                 vaux(1) = bonds(2,j)
                 lcheck(j) = .FALSE.
               else if ( bonds(2,j) .eq. rbond(1) ) then
                 vaux(1) = bonds(1,j)
                 lcheck(j) = .FALSE.
               else if ( bonds(1,j) .eq. rbond(2) ) then
                 vaux(4) = bonds(2,j)
                 lcheck(j) = .FALSE.
               else if ( bonds(2,j) .eq. rbond(2) ) then
                 vaux(4) = bonds(1,j)
                 lcheck(j) = .FALSE.
               end if
             end if
           end do
!
! Computing equilibrium value
!
           if ( idat(vaux(2)) .gt. idat(vaux(3)) ) then
!
             itmp    = vaux(2)
             vaux(2) = vaux(3)
             vaux(3) = itmp
!
             itmp    = vaux(1)
             vaux(1) = vaux(4)
             vaux(4) = itmp
!
           else if ( idat(vaux(2)) .eq. idat(vaux(3)) ) then
             if ( idat(vaux(1)) .gt. idat(vaux(4)) ) then
!
               itmp    = vaux(2)
               vaux(2) = vaux(3)
               vaux(3) = itmp
!
               itmp    = vaux(1)
               vaux(1) = vaux(4)
               vaux(4) = itmp
!
             end if
           end if
!
           call Diedro(coord(:,vaux(1)),coord(:,vaux(2)),              & 
                       coord(:,vaux(3)),coord(:,vaux(4)),daux)
           daux = daux*180.0d0/pi
!
! If dihedral angle is rigid classify it as improper 
!
           if ( lrigid(vaux(2),vaux(3)) ) then
!
             dihed%nrigid = dihed%nrigid + 1
!
             dihed%irigid(:,dihed%nrigid) = vaux(:)
             dihed%frigid(dihed%nrigid)   = 2
!~              drigid(nrigid)   = daux     ! TODO: option to pick equilibrium value or set to 0/180
             if ( abs(daux) .le. 90 ) then  ! TODO: option to pick equilibrium value or set to 0/180
               dihed%drigid(dihed%nrigid) = 0.0d0
             else
               dihed%drigid(dihed%nrigid) = 180.0d0
             end if
!
! If dihedral angle belongs to a cycle classify it as inversion 
!
           else if ( lcycle(vaux(2),vaux(3)) ) then
!
             dihed%ninv = dihed%ninv + 1
!
             dihed%iinv(:,dihed%ninv) = vaux(:)
             dihed%dinv(dihed%ninv)   = daux  !  TODO: how to handle equilibrium value
             dihed%finv(dihed%ninv)   = 3
!
           else
!
             dihed%nflexi = dihed%nflexi + 1
!
             dihed%iflexi(:,dihed%nflexi) = vaux(:)
             dihed%dflexi(dihed%nflexi)   = daux ! TODO: option to pick equilibrium value or set to ideal value
!
           end if
!
         end if
       end do
!
! Generating Fourier series for each flexible dihedral
! ----------------------------------------------------
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
! TODO: option to keep only one quadruplet per torsion
!
! Generating Fourier series
!
       do i = 1, dihed%nflexi
!         
         dihed%flexi(i)%ntor = 7 ! TODO: choose number of functions in Fourier expansion
         allocate(dihed%flexi(i)%tor(dihed%flexi(i)%ntor)) 
!
         dihed%flexi(i)%itor(:) = dihed%iflexi(:,i)
         dihed%fflexi(i)        = 9
!
         do j = 1, dihed%flexi(i)%ntor
           dihed%flexi(i)%tor(j)%vtor  = 0.0d0
           dihed%flexi(i)%tor(j)%phase = dihed%dflexi(i) ! TODO: shift phase for even multiplicities
           dihed%flexi(i)%tor(j)%multi = j-1
         end do
!
       end do
!
       dihed%ndihe = dihed%nrigid + dihed%ninv + dihed%nimpro + dihed%nflexi
!
       return
       end subroutine setdihe
!
!======================================================================!
!
! CHECKLAB - CHECK LABel
!
! This subroutine 
!
       subroutine checklab(idihe,dp1,dp2)
!
       implicit none
!
! Input/output variables
!
       integer,dimension(4),intent(inout)  ::  idihe  !  
       integer,intent(in)                  ::  dp1    !  
       integer,intent(in)                  ::  dp2    !  
!
!  
! ----------------------------------------
!
       if ( dp1 .lt. dp2 ) then
         idihe(3) = dp1
         idihe(4) = dp2
       else
         idihe(3) = dp2
         idihe(4) = dp1
       end if
!
       return
       end subroutine checklab
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
! LABDIHE - LABeling DIHEdrals
!
! This subroutine 
!
       subroutine labdihe(nat,idat,nidat,ndihe,idihe,iddihe,niddihe,   &
                          nimpro,nrigid,nflexi,debug)
!
       use graphtools, only:  inormlabels
!
       implicit none
!
! Input/output variables
!
       integer,dimension(4,ndihe),intent(in)  ::  idihe    !
       integer,dimension(nat),intent(in)      ::  idat     ! 
       integer,dimension(ndihe),intent(out)   ::  iddihe   !
       integer,intent(in)                     ::  nat      !  
       integer,intent(in)                     ::  ndihe    !  
       integer,intent(in)                     ::  nidat    !  
       integer,intent(out)                    ::  niddihe  !  
       logical,intent(in)                     ::  debug    !
       integer,intent(in)                     ::  nimpro   !
       integer,intent(in)                     ::  nrigid   !
       integer,intent(in)                     ::  nflexi   !
!
! Local variables
!
       integer,dimension(ndihe)               ::  ivaux    !
       integer                                ::  i,j      !
!
!  Combining labels of adjacent nodes
! -----------------------------------
!
!~        do i = 1, nimpro
!~          ivaux(i) = min(idat(idihe(1,i)),idat(idihe(2,i)))*nidat**3    &
!~                   + max(idat(idihe(1,i)),idat(idihe(2,i)))*nidat**2    &
!~                   + min(idat(idihe(3,i)),idat(idihe(4,i)))*nidat       &
!~                   + max(idat(idihe(3,i)),idat(idihe(4,i)))
!~        end do
!~ !
!~        do i = nimpro+1, ndihe
!~          ivaux(i) = min(idat(idihe(2,i)),idat(idihe(3,i)))*nidat**3    &
!~                   + max(idat(idihe(2,i)),idat(idihe(3,i)))*nidat**2    &
!~                   + min(idat(idihe(1,i)),idat(idihe(4,i)))*nidat       &
!~                   + max(idat(idihe(1,i)),idat(idihe(4,i)))
!~        end do
!
!~        if ( debug ) then
!~          write(*,*) 'Non-normalized combined labels'
!~          write(*,*) '------------------------------'
!~          do i = 1, ndihe
!~            write(*,'(200(1X,I3))') i,(idat(idihe(j,i)),j=1,4),ivaux(i)
!~          end do
!~          write(*,*)
!~        end if
!
       call inormlabels(ndihe,ivaux,iddihe,niddihe)
!
!~        if ( debug ) then
!~          write(*,*) 'Normalized combined labels'
!~          write(*,*) '--------------------------'
!~          do i = 1, ndihe
!~            write(*,'(200(1X,I3))') i,ivaux(i),iddihe(i)
!~          end do
!~          write(*,*)
!~        end if
!
! Including dihedral type information in labels
!
       do i = 1, nimpro
         ivaux(i) = 3*iddihe(i) + 1
       end do

       do i = nimpro+1, nimpro+nrigid
         ivaux(i) = 3*iddihe(i) + 2
       end do
!
       do i = nimpro+nrigid+1, ndihe
         ivaux(i) = 3*iddihe(i) + 3
       end do
!
       call inormlabels(ndihe,ivaux,iddihe,niddihe)
!
       if ( debug ) then
         write(*,*) 'Normalized dihedral labels'
         write(*,*) '--------------------------'
         do i = 1, ndihe
           write(*,'(200(1X,I3))') i,ivaux(i),iddihe(i)
         end do
         write(*,*)
       end if
!
       return
       end subroutine labdihe
!
!======================================================================!
!
! SYMTERM - SYMmetrizing force field TERMs
!
! This subroutine 
!
       subroutine symterm(nbond,idbond,dbond,sbond,labbond,fsymm,debug)
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
       logical,intent(in)                             ::  debug   !
       logical,intent(in)                             ::  fsymm   !
!
! Equivalent atoms information
!
       integer,dimension(nbond)                       ::  mol      !
       integer,dimension(nbond)                       ::  agg      !
       integer,dimension(nbond)                       ::  tag      !
       integer,dimension(nbond)                       ::  imol     !
       integer,dimension(nbond)                       ::  iagg     !
       integer,dimension(nbond)                       ::  itag     !
       integer,dimension(nbond)                       ::  nmol     !
       integer,dimension(nbond)                       ::  nagg     !
       integer,dimension(nbond)                       ::  ntag     !
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
       subroutine symoop(ndihe,iddihe,ddihe,sdihe,labdihe,fsymm,debug)
!
       use printings
       use graphtools, only:  blockdiag,inormlabels
!
       implicit none
!
! Input/output variables
!
       character(len=50),dimension(ndihe),intent(in)  ::  labdihe !
       real(kind=8),dimension(ndihe),intent(inout)    ::  ddihe   !
       integer,dimension(ndihe),intent(out)           ::  iddihe  !
       integer,dimension(ndihe),intent(in)            ::  sdihe   !
       integer,intent(in)                             ::  ndihe   !  
       logical,intent(in)                             ::  fsymm   ! 
       logical,intent(in)                             ::  debug   !
!
! Equivalent atoms information
!
       integer,dimension(ndihe)                       ::  mol      !
       integer,dimension(ndihe)                       ::  agg      !
       integer,dimension(ndihe)                       ::  tag      !
       integer,dimension(ndihe)                       ::  imol     !
       integer,dimension(ndihe)                       ::  iagg     !
       integer,dimension(ndihe)                       ::  itag     !
       integer,dimension(ndihe)                       ::  nmol     !
       integer,dimension(ndihe)                       ::  nagg     !
       integer,dimension(ndihe)                       ::  ntag     !
       integer                                        ::  magg     !  Number of aggregates
       integer                                        ::  nsize    !  Maximum aggregate size
!
! Local variables
!
       logical,dimension(ndihe,ndihe)                 ::  adj      !
       real(kind=8)                                   ::  daux     !
       integer                                        ::  i,j      !
       integer                                        ::  k,l      !
!
       real(kind=8)                                   ::  thr = 5.0d0
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
!             daux = 0     ! TODO: symmetrize absolute value but keep sign
!             do l = k+1, k+i
!                daux = daux + ddihe(mol(l))
!             end do
!
!             daux = daux/real(i)
!             do l = k+1, k+i
!               ddihe(mol(l)) = daux
!             end do
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
       integer,dimension(4,ndihe)                     ::  idihe    !
       integer,dimension(ndihe)                       ::  sdihe    !
       integer,dimension(ndihe)                       ::  mdihe    !
       integer,dimension(ndihe)                       ::  iddihe   !
       integer,dimension(ndihe)                       ::  ivaux    !
       integer                                        ::  niddihe  !
       integer                                        ::  nterm    !
       integer                                        ::  i,j,k    !
!
       real(kind=8)                                   ::  thr = 5.0d0
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
!
           ivaux(k) = dihe%idflexi(i)*10 + mdihe(k)
!~ write(*,*) k,trim(labdihe(k)),dihe%idflexi(i),mdihe(k),ivaux(k)
!
         end do
       end do
       nterm = k
!
       call inormlabels(ndihe,ivaux,iddihe,niddihe)
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
       do i = 1, dihe%nquad
         do j = 1, nterm
           if ( lquad(j) ) cycle
           if ( (idihe(1,j).eq.dihe%iquad(1,i))                        &
                 .and. (idihe(2,j).eq.dihe%iquad(2,i))                 &
                 .and. (idihe(3,j).eq.dihe%iquad(3,i))                 &
                 .and. (idihe(4,j).eq.dihe%iquad(4,i)) ) then                 
             lquad(j) = .TRUE.
           end if
         end do
       end do
!
! Printing equivalent dihedral terms with same identifier based on 
!  atomtypes and multiplicities 
!
       do i = 1, nterm
         if ( lquad(i) ) then
           do j = 1, nterm
             if ( lquad(j) ) cycle
             if ( iddihe(i) .eq. iddihe(j) ) then
!~ write(*,*) 'matching',i,iddihe(i),'with',j,iddihe(j)
               write(unideps,'(3X,I4,1X,A,1X,I4,A)')                   &
                                  sdihe(j),'=',sdihe(i),'*1.d0 ; '//   &
                               trim(labdihe(j))//' = '//trim(labdihe(i))
               write(uniscr,'(3X,I4,1X,A,1X,I4,A)')                    &
                                  sdihe(j),'=',sdihe(i),'*1.d0 ; '//   &
                               trim(labdihe(j))//' = '//trim(labdihe(i))
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
                           narunit,arunit,adj,lheavy,fsymm,debug)
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
       logical,intent(in)                              ::  fsymm     !   
       logical,intent(in)                              ::  debug     !
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
       integer,dimension(ndihe)                        ::  itag      !
       integer,dimension(ndihe)                        ::  nmol      !
       integer,dimension(ndihe)                        ::  nagg      !
       integer,dimension(ndihe)                        ::  ntag      !
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
       real(kind=8)                                       ::  daux      !
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
! Generating unique representation with rigid and improper dihedrals
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
       chkimpro(:) = .TRUE.
!
       do i = 1, marunit
!
! Generating list of atoms belonging to the aromatic unit
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
!
! Block-diagonalization of the adjacency matrix yields an array repre-
!  sentation with the equivalent bond terms information
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
       end module genforcefield
!
!======================================================================!

!======================================================================!
!
       module genforcefield
!
       use lengths,       only:  leninp,lenlab
       use units,         only:  uniic,unideps,unitmp,unidihe,         &
                                 unijoyce,uninb,unitop,unicorr
!
       implicit none
!
       private
       public  ::  gentypes,                                           &
                   genffbonded,                                        &
                   symffbonded,                                        &
                   eqvatms,                                            &
                   nearnei,                                            &
                   genrigidlist,                                       &
                   gencyclelist,                                       &
                   bonded2dihe
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
       subroutine genffbonded(nat,coord,adj,ideg,lcycle,lrigid,        &
                              znum,bonded,dihed,iroute,debug)
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
       integer,dimension(nat),intent(in)                  ::  ideg      ! 
       logical,dimension(nat,nat),intent(in)              ::  adj       !  Boolean adjacency 
       logical,dimension(nat,nat),intent(in)              ::  lcycle    !  Bonds belonging to rings
       logical,dimension(nat,nat),intent(in)              ::  lrigid    ! 
!
       integer,dimension(nat),intent(in)                  ::  znum      !  
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
         call gendihe(nat,coord,znum,ideg,lcycle,lrigid,bonded%nbond,  &
                      bonded%ibond,bonded%nang,bonded%iang,adjang,     &
                      edgeang,dihed,iroute,debug)
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
       subroutine gendihe(nat,coord,znum,ideg,lcycle,lrigid,nbond,     &
                          ibond,nang,iang,adjang,edgeang,dihed,        &
                          iroute,debug)
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
       integer,dimension(nat),intent(in)                  ::  znum     ! 
       integer,dimension(nat),intent(in)                  ::  ideg     ! 
       integer,intent(in)                                 ::  nat      !  Number of atoms
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
       dihed%iimpro(:,:) = 0
       dihed%fimpro(:)   = 0
       dihed%dimpro(:)   = 0.0d0
       dihed%kimpro(:)   = 0.0d0
       dihed%nimpro      = 0
!
       allocate(dihed%iinv(4,ndihe),dihed%dinv(ndihe),                 &
                dihed%kinv(ndihe),dihed%finv(ndihe))
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
       call setdihe(nat,coord,nbond,ibond,nang,edgeang,iang,           &
                    ndihe,dihed,edges,lcycle,lrigid,znum,ideg,iroute)
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
       subroutine symffbonded(nat,nidat,idat,bonded,dihed,debug)
!
       use datatypes,   only: grobonded,                               &
                              dihedrals
       use graphtools,  only: inormlabels
!
       implicit none
!
! Input/output variables
!
       type(grobonded),intent(inout)                ::  bonded    !
       type(dihedrals),intent(inout)                ::  dihed     !
       integer,dimension(nat),intent(in)            ::  idat      !
!
       integer,intent(in)                           ::  nat       ! 
       integer,intent(in)                           ::  nidat     ! 
!
       logical,intent(in)                           ::  debug     !  Debug mode
!
! Local variables
!  
       integer,dimension(:),allocatable             ::  ivaux     !
       integer,dimension(4)                         ::  vaux      !
       integer                                      ::  iaux      !
       integer                                      ::  iterm     !
       integer                                      ::  i,j       !
!
! Symmetrizing bonded terms
! -------------------------
!
!  Finding equivalent bond-stretching terms
!
       allocate(bonded%sbond(bonded%nbond),bonded%idbond(bonded%nbond))
!
       call labbond(nat,idat,nidat,bonded%nbond,bonded%ibond,          &
                    bonded%idbond,bonded%nidbond,debug)
!
       iterm = 0
       do i = 1, bonded%nbond
         iterm = iterm + 1
         bonded%sbond(i) = iterm
       end do
!
       call symterm(bonded%nbond,bonded%idbond,bonded%bond,            &
                    bonded%sbond,debug)
!
!  Finding equivalent angle-bending terms
!
       if ( bonded%nang .eq. 0 ) return
!
       allocate(bonded%sang(bonded%nang),bonded%idang(bonded%nang))
!
       call labangle(nat,idat,nidat,bonded%nang,bonded%iang,           &
                     bonded%idang,bonded%nidang,debug)
!
       do i = 1, bonded%nang
         iterm = iterm + 1
         bonded%sang(i) = iterm
       end do
!
       call symterm(bonded%nang,bonded%idang,bonded%ang,bonded%sang,   &
                    debug)
!
!  Finding equivalent dihedrals terms
!
       if ( bonded%ndihe .eq. 0 ) return
!
       allocate(ivaux(dihed%ndihe))
!
       allocate(dihed%simpro(dihed%ndihe),dihed%idimpro(dihed%ndihe))
!
       do i = 1, dihed%nimpro
         iterm = iterm + 1
         vaux(:) = dihed%iimpro(:,i)
         iaux = min(idat(vaux(1)),idat(vaux(2)))*nidat**3              &
                + max(idat(vaux(1)),idat(vaux(2)))*nidat**2            &
                + min(idat(vaux(3)),idat(vaux(4)))*nidat               &
                + max(idat(vaux(3)),idat(vaux(4)))
!
         dihed%simpro(i)  = iterm  
         dihed%idimpro(i) = iaux  
       end do
!
       allocate(dihed%sinv(dihed%ndihe),dihed%idinv(dihed%ndihe))
!
       do i = 1, dihed%ninv
         iterm = iterm + 1
         vaux(:) = dihed%iinv(:,i)
         iaux = min(idat(vaux(1)),idat(vaux(2)))*nidat**3              &
                + max(idat(vaux(1)),idat(vaux(2)))*nidat**2            &
                + min(idat(vaux(3)),idat(vaux(4)))*nidat               &
                + max(idat(vaux(3)),idat(vaux(4)))
!
         dihed%sinv(i)  = iterm  
         dihed%idinv(i) = iaux  
       end do
!
       allocate(dihed%srigid(dihed%ndihe),dihed%idrigid(dihed%ndihe))
!
       do i = 1, dihed%nrigid
         iterm = iterm + 1
         vaux(:) = dihed%irigid(:,i)
         iaux = min(idat(vaux(2)),idat(vaux(3)))*nidat**3              &
              + max(idat(vaux(2)),idat(vaux(3)))*nidat**2              &
              + min(idat(vaux(1)),idat(vaux(4)))*nidat                 &
              + max(idat(vaux(1)),idat(vaux(4)))
!
         dihed%srigid(i)  = iterm  
         dihed%idrigid(i) = iaux  
       end do
!
       allocate(dihed%sflexi(dihed%ndihe),dihed%idflexi(dihed%ndihe))
!
       do i = 1, dihed%nflexi
!~          iterm = iterm + 1
         vaux(:) = dihed%iflexi(:,i)
         iaux = min(idat(vaux(2)),idat(vaux(3)))*nidat**3              &
              + max(idat(vaux(2)),idat(vaux(3)))*nidat**2              &
              + min(idat(vaux(1)),idat(vaux(4)))*nidat                 &
              + max(idat(vaux(1)),idat(vaux(4)))
!
!~          dihed%sflexi(i)  = iterm  
         dihed%idflexi(i) = iaux  
       end do
!
!~        do i = 1, dihed%nimpro
!~          iterm = iterm + 1
!~          dihed%simpro(i) = iterm
!~        end do
!~ !
!~        do i = 1, dihed%ninv
!~          iterm = iterm + 1
!~          dihed%sinv(i) = iterm
!~        end do
!~ !
!~        do i = 1, dihed%nrigid
!~          iterm = iterm + 1
!~          dihed%srigid(i) = iterm
!~        end do
!
       do i = 1, dihed%nflexi
         do j = 1, dihed%flexi(i)%ntor
           iterm = iterm + 1
           dihed%flexi(i)%tor(j)%stor = iterm
         end do
       end do
!
       ivaux(:) = dihed%idimpro(:)
       call inormlabels(dihed%nimpro,ivaux(:dihed%nimpro),dihed%idimpro(:dihed%nimpro),iaux)
!
       ivaux(:) = dihed%idinv(:)
       call inormlabels(dihed%ninv,ivaux(:dihed%ninv),dihed%idinv(:dihed%ninv),iaux)
!
       ivaux(:) = dihed%idrigid(:)
       call inormlabels(dihed%nrigid,ivaux(:dihed%nrigid),dihed%idrigid(:dihed%nrigid),iaux)
!
       ivaux(:) = dihed%idflexi(:)
       call inormlabels(dihed%nflexi,ivaux(:dihed%nflexi),dihed%idflexi(:dihed%nflexi),iaux)
!
! Symmetrizing dihedral terms
!
! FIXME: problem symmetrizing -179 and 179 dihedrals
       if ( dihed%nimpro .gt. 0 ) then
         call symdihe(dihed%nimpro,dihed%idimpro,dihed%dimpro,dihed%simpro,.TRUE.,debug)
       end if 
!
       if ( dihed%ninv .gt. 0 ) then
         call symdihe(dihed%ninv,dihed%idinv,dihed%dinv,dihed%sinv,.TRUE.,debug)
       end if 
!
       if ( dihed%nrigid .gt. 0 ) then
         call symdihe(dihed%nrigid,dihed%idrigid,dihed%drigid,dihed%srigid,.TRUE.,debug)
       end if 
!
       if ( dihed%nflexi .gt. 0 ) then
         call symdihe(dihed%nflexi,dihed%idflexi,dihed%dflexi,dihed%sflexi,.FALSE.,debug)
       end if
!
       deallocate(ivaux)
!
       return
       end subroutine symffbonded
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
       subroutine genrigidlist(nat,wiberg,lcycle,lrigid)   
!
       implicit none
!
! Input/output variables
!
       real(kind=8),dimension(nat,nat),intent(in)  ::  wiberg  !
       logical,dimension(nat,nat),intent(in)       ::  lcycle  !  Bonds belonging to rings
       logical,dimension(nat,nat),intent(out)      ::  lrigid  !  Rigid bonds
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
!
! We consider rigid bond 
!  if Wiberg index is greater than 1.5 it is a double bond
!  if Wiberg index is greater than 1.1 and the bond belongs to a cycle
!   then it belongs to an aromatic cycle
!
       do i = 1, nat-1
         do j = i+1, nat
!
           if ( (wiberg(j,i).ge.1.5d0) .or.                            &
                         ((wiberg(j,i).ge.1.1d0).and.lcycle(j,i)) ) then 
             lrigid(i,j) = .TRUE.
             lrigid(j,i) = .TRUE.
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
       subroutine setdihe(nat,coord,nbond,ibond,nang,edgeang,iang,     &
                          ndihe,dihed,edgedihe,lcycle,lrigid,znum,     &
                          ideg,iroute)
!
       use datatypes,   only: dihedrals
!
       implicit none
!
! Input/output variables
!
       type(dihedrals),intent(inout)              ::  dihed     !                            
       real(kind=8),dimension(3,nat),intent(in)   ::  coord     !
       logical,dimension(nat,nat),intent(in)      ::  lcycle    !  Bonds belonging to rings
       logical,dimension(nat,nat),intent(in)      ::  lrigid    !  Bonds belonging to rings
       integer,dimension(nat),intent(in)          ::  znum      ! 
       integer,dimension(nat),intent(in)          ::  ideg      ! 
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
       logical,dimension(ndihe)                   ::  torsion   !
       logical                                    ::  flag      !
       real(kind=8)                               ::  daux      !
       integer,dimension(4,ndihe)                 ::  auxdihe   !  
       integer,dimension(ndihe)                   ::  ivaux     !  
       integer,dimension(4)                       ::  vaux      !  
       integer,dimension(2)                       ::  bond1     !
       integer,dimension(2)                       ::  bond2     !
       integer,dimension(2)                       ::  bond3     !
       integer,dimension(2)                       ::  bond4     !
       integer,dimension(3)                       ::  ang1      !
       integer,dimension(3)                       ::  ang2      !
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
               dihed%dimpro(dihed%nimpro)   = 0.0d0  ! TODO: option to pick equilibrium value or set to 0
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
             dihed%finv(dihed%ninv)   = 3
             dihed%dinv(dihed%ninv)   = daux  !  TODO: how to handle equilibrium value
!
           end if
!
         end if
!
       end do
!
!  Finding torsional dihedrals
!  ...........................
!
       do i = 1, ndihe
         if ( torsion(i) ) then
!
           bond1(:) = ibond(:,edgeang(1,edgedihe(1,i)))
           bond2(:) = ibond(:,edgeang(2,edgedihe(1,i)))
           bond3(:) = ibond(:,edgeang(1,edgedihe(2,i)))
           bond4(:) = ibond(:,edgeang(2,edgedihe(2,i)))
!
           do j = 1, 2
             do k = 1, 2
!
               if ( bond1(j) .eq. bond2(k) ) then
                 iaux1 = j
                 iaux2 = k
               end if
!
               if ( bond3(j) .eq. bond4(k) ) then
                 iaux3 = j
                 iaux4 = k
               end if
!
             end do 
           end do
!
           if ( iaux1 .eq. 1 ) then
             iaux1 = bond1(2)
           else
             iaux1 = bond1(1)
           end if
!
           if ( iaux2 .eq. 1 ) then
             iaux2 = bond2(2)
           else
             iaux2 = bond2(1)
           end if
!
           if ( iaux3 .eq. 1 ) then
             iaux3 = bond3(2)
           else
             iaux3 = bond3(1)
           end if
!
           if ( iaux4 .eq. 1 ) then
             iaux4 = bond4(2)
           else
             iaux4 = bond4(1)
           end if
!
           if ( (bond1(1).eq.bond3(1)) .and. (bond1(2).eq.bond3(2)) ) then
!
             vaux(1) = iaux2
             vaux(2) = bond1(1)
             vaux(3) = bond1(2)
             vaux(4) = iaux4
!
           else if ( (bond1(1).eq.bond4(1)) .and. (bond1(2).eq.bond4(2)) ) then
!
             vaux(1) = iaux2
             vaux(2) = bond1(1)
             vaux(3) = bond1(2)
             vaux(4) = iaux3
!
           else if ( (bond2(1).eq.bond3(1)) .and. (bond2(2).eq.bond3(2)) ) then
!
             vaux(1) = iaux1
             vaux(2) = bond2(1)
             vaux(3) = bond2(2)
             vaux(4) = iaux4
!
           else if ( (bond2(1).eq.bond4(1)) .and. (bond2(2).eq.bond4(2)) ) then
!
             vaux(1) = iaux1
             vaux(2) = bond2(1)
             vaux(3) = bond2(2)
             vaux(4) = iaux3
!
           end if
!
! Computing equilibrium value
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
!~              drigid(nrigid)   = daux ! TODO: option to pick equilibrium value or set to 0/180
             if ( abs(daux) .le. 90 ) then
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
       do i = 1, dihed%nflexi
!         
         dihed%flexi(i)%ntor = 6 ! TODO: choose number of functions in Fourier expansion
         allocate(dihed%flexi(i)%tor(dihed%flexi(i)%ntor)) 
!
         dihed%flexi(i)%itor(:) = dihed%iflexi(:,i)
!
         do j = 1, dihed%flexi(i)%ntor
           dihed%flexi(i)%tor(j)%vtor  = 0.0d0
           dihed%flexi(i)%tor(j)%phase = dihed%dflexi(i) ! TODO: shift phase for even multiplicities
           dihed%flexi(i)%tor(j)%multi = j
         end do
!
       end do
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
       integer,dimension(2,nbond),intent(in)   ::  ibond    !
       integer,dimension(nat),intent(in)       ::  idat     ! 
       integer,dimension(nbond),intent(out)    ::  idbond   !
       integer,intent(in)                      ::  nat      !  
       integer,intent(in)                      ::  nbond    !  
       integer,intent(in)                      ::  nidat    !  
       integer,intent(out)                     ::  nidbond  !  
       logical,intent(in)                      ::  debug    !
!
! Local variables
!
       integer,dimension(nbond)                ::  ivaux    !
       integer                                 ::  i,j      !
!
!  Combining labels of adjacent nodes
! -----------------------------------
!
       do i = 1, nbond
         ivaux(i) = min(idat(ibond(1,i)),idat(ibond(2,i)))*nidat       &
                                + max(idat(ibond(1,i)),idat(ibond(2,i)))
       end do
!
       if ( debug ) then
         write(*,*) 'Non-normalized combined labels'
         write(*,*) '------------------------------'
         do i = 1, nbond
           write(*,'(200(1X,I3))') i,(idat(ibond(j,i)),j=1,2),ivaux(i)
         end do
         write(*,*)
       end if
!
       call inormlabels(nbond,ivaux,idbond,nidbond)
!
       if ( debug ) then
         write(*,*) 'Normalized combined labels'
         write(*,*) '--------------------------'
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
       integer,dimension(3,nang),intent(in)  ::  iang    !
       integer,dimension(nat),intent(in)     ::  idat    ! 
       integer,dimension(nang),intent(out)   ::  idang   !
       integer,intent(in)                    ::  nat     !  
       integer,intent(in)                    ::  nang    !  
       integer,intent(in)                    ::  nidat   !  
       integer,intent(out)                   ::  nidang  !  
       logical,intent(in)                    ::  debug   !
!
! Local variables
!
       integer,dimension(nang)               ::  ivaux   !
       integer                               ::  i,j     !
!
!  Combining labels of adjacent nodes
! -----------------------------------
!
       do i = 1, nang
         ivaux(i) = idat(iang(2,i))*nidat**2                           &
                        + min(idat(iang(1,i)),idat(iang(3,i)))*nidat   &
                                  + max(idat(iang(1,i)),idat(iang(3,i)))
       end do
!
       if ( debug ) then
         write(*,*) 'Non-normalized combined labels'
         write(*,*) '------------------------------'
         do i = 1, nang
           write(*,'(200(1X,I3))') i,(idat(iang(j,i)),j=1,3),ivaux(i)
         end do
         write(*,*)
       end if
!
       call inormlabels(nang,ivaux,idang,nidang)
!
       if ( debug ) then
         write(*,*) 'Normalized combined labels'
         write(*,*) '--------------------------'
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
         write(*,*) 'Normalized specific labels'
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
       subroutine symterm(nbond,idbond,dbond,sbond,debug)
!
       use printings
       use graphtools, only:  blockdiag,inormlabels
!
       implicit none
!
! Input/output variables
!
       real(kind=8),dimension(nbond),intent(inout)  ::  dbond   !
       integer,dimension(nbond),intent(out)         ::  idbond  !
       integer,dimension(nbond),intent(in)          ::  sbond   !
       integer,intent(in)                           ::  nbond   !  
       logical,intent(in)                           ::  debug   !
!
! Equivalent atoms information
!
       integer,dimension(nbond)                     ::  mol      !
       integer,dimension(nbond)                     ::  agg      !
       integer,dimension(nbond)                     ::  tag      !
       integer,dimension(nbond)                     ::  imol     !
       integer,dimension(nbond)                     ::  iagg     !
       integer,dimension(nbond)                     ::  itag     !
       integer,dimension(nbond)                     ::  nmol     !
       integer,dimension(nbond)                     ::  nagg     !
       integer,dimension(nbond)                     ::  ntag     !
       integer                                      ::  magg     !  Number of aggregates
       integer                                      ::  nsize    !  Maximum aggregate size
!
! Local variables
!
       logical,dimension(nbond,nbond)               ::  adj      !
       real(kind=8)                                 ::  daux     !
       integer                                      ::  i,j      !
       integer                                      ::  k,l      !
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
!~              do l = k+2, k+i
!~                write(unideps,'(3X,I4,1X,A,1X,I4,A)') sbond(mol(l)),'=',sbond(mol(k+1)),'*1.d0'
!~                write(unitmp,'(3X,I4,1X,A,1X,I4,A)') sbond(mol(l)),'=',sbond(mol(k+1)),'*1.d0'
!~              end do
!
             k = k + i
!
           end do
         end if
       end do
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
         call print_info(0,nbond,agg,tag,mol,'agg','tag','mol')
       end if
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
       subroutine symdihe(ndihe,iddihe,ddihe,sdihe,ftmp,debug)
!
       use printings
       use graphtools, only:  blockdiag,inormlabels
!
       implicit none
!
! Input/output variables
!
       real(kind=8),dimension(ndihe),intent(inout)  ::  ddihe   !
       integer,dimension(ndihe),intent(out)         ::  iddihe  !
       integer,dimension(ndihe),intent(in)          ::  sdihe   !
       integer,intent(in)                           ::  ndihe   !  
       logical,intent(in)                           ::  ftmp    !  
       logical,intent(in)                           ::  debug   !
!
! Equivalent atoms information
!
       integer,dimension(ndihe)                     ::  mol      !
       integer,dimension(ndihe)                     ::  agg      !
       integer,dimension(ndihe)                     ::  tag      !
       integer,dimension(ndihe)                     ::  imol     !
       integer,dimension(ndihe)                     ::  iagg     !
       integer,dimension(ndihe)                     ::  itag     !
       integer,dimension(ndihe)                     ::  nmol     !
       integer,dimension(ndihe)                     ::  nagg     !
       integer,dimension(ndihe)                     ::  ntag     !
       integer                                      ::  magg     !  Number of aggregates
       integer                                      ::  nsize    !  Maximum aggregate size
!
! Local variables
!
       logical,dimension(ndihe,ndihe)               ::  adj      !
       real(kind=8)                                 ::  daux     !
       integer                                      ::  i,j      !
       integer                                      ::  k,l      !
!
       real(kind=8)                                 ::  thr = 5.0d0
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
!           if ( (iddihe(i).eq.iddihe(j)) .and.                         &
!                                  (abs(ddihe(i)-ddihe(j)).lt.thr) ) then
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
!             daux = 0
!             do l = k+1, k+i
!                daux = daux + ddihe(mol(l))
!             end do
!
!             daux = daux/real(i)
!             do l = k+1, k+i
!               ddihe(mol(l)) = daux
!             end do
!
!~              do l = k+2, k+i
!~                write(unideps,'(3X,I4,1X,A,1X,I4,A)') sdihe(mol(l)),'=',sdihe(mol(k+1)),'*1.d0'
!~                if ( ftmp )write(unitmp,'(3X,I4,1X,A,1X,I4,A)') sdihe(mol(l)),'=',sdihe(mol(k+1)),'*1.d0'
!~              end do
!
             k = k + i
!
           end do
         end if
       end do
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
         call print_info(0,ndihe,agg,tag,mol,'agg','tag','mol')
       end if
!
       return
       end subroutine symdihe
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
    function calc_atm_dihed(atom1,atom2,atom3,atom4) result(dh)

        !----------------------------------------------------------------
        ! Esta función coloca el dihedro en una orientación estandard
        ! y mide el dihdro. Se obtienen valores con signo de acuerod
        ! con el criterio establecido.
        ! NOTE:
        ! The SR works with individual x(), y(), z() vectors, not atomic
        ! 3D vectors. Thus the SR is first interfaced to inpunt 3D atomic
        ! vectors as default in this module.
        !----------------------------------------------------------------

        real(kind=8),dimension(3),intent(in) :: atom1, atom2, atom3, atom4
        real(kind=8)::dh
        real(kind=8),dimension(1:4)::x,y,z
        real(kind=8)::dx,dy,dz,z_aux,y_aux,theta
        integer::i
        real(kind=8),parameter                   ::  pi =  4*atan(1.0_8) 

        !Interface: atom to vectors (in this case x,y,z)
        x(1:4)=(/atom1(1),atom2(1),atom3(1),atom4(1)/)
        y(1:4)=(/atom1(2),atom2(2),atom3(2),atom4(2)/)
        z(1:4)=(/atom1(3),atom2(3),atom3(3),atom4(3)/)


        !1) Center the system at Atom3
        dx=x(3)
        dy=y(3)
        dz=z(3)
        do i=1,4
            x(i)=x(i)-dx
            y(i)=y(i)-dy
            z(i)=z(i)-dz
        enddo


        !2) Lie Atom2--Atom3 on the z-axis.
        !   a. Rotation around X
        theta=atan(y(2)/z(2))
        z_aux=y(2)*sin(theta)+z(2)*cos(theta)
        if (z_aux>0) theta=theta+pi ! Note: 2 is behind 3
        do i=1,4
            z_aux=z(i)
            z(i)=y(i)*sin(theta)+z_aux*cos(theta)
            y(i)=y(i)*cos(theta)-z_aux*sin(theta)
        enddo
        !   b. Rotation around Y
        theta=atan(x(2)/z(2))
        z_aux=x(2)*sin(theta)+z(2)*cos(theta)
        if (z_aux>0) theta=theta+pi
        do i=1,4
            z_aux=z(i)
            z(i)=x(i)*sin(theta)+z_aux*cos(theta)
            x(i)=x(i)*cos(theta)-z_aux*sin(theta)
        enddo

        !3) Put Atom1 on the YZ plane (y>0) 
        !   a. Rotation around Z
        theta=atan(x(1)/y(1))
        y_aux=x(1)*sin(theta)+y(1)*cos(theta)
        if (y_aux<0) theta=theta+pi ! Note: 1 has positive y
        do i=1,4
            y_aux=y(i)
            y(i)=x(i)*sin(theta)+y_aux*cos(theta)
            x(i)=x(i)*cos(theta)-y_aux*sin(theta)
        enddo
        theta=atan(x(1)/y(1))
        y_aux=x(1)*sin(theta)+y(1)*cos(theta)
        if (y_aux<0) theta=theta+pi ! Note: 1 has positive y
        do i=1,4
            y_aux=y(i)
            y(i)=x(i)*sin(theta)+y_aux*cos(theta)
            x(i)=x(i)*cos(theta)-y_aux*sin(theta)
        enddo

        dh=atan(x(4)/y(4))
        ! Comprobamos que de los dos ángulos posibles (en +/- 180)
        ! obtenemos el que nos interesa comparando con el criterio
        ! según el cual signo(dh)=signo(x(4))
        if (dh>0 .and. x(4)<0) dh=dh-pi
        if (dh<0 .and. x(4)>0) dh=dh+pi

        ! Sale con el signo cambiado!
        dh=-dh
        
        return

     end function calc_atm_dihed
!
!======================================================================!
!
    function calc_improper(r1,r2,r3,r4) result(im)

        !----------------------------------------------------------------
        ! Ańgulo formado por el vector 4-1 con el plano 2-3-4
        !        2
        !       /
        !   1--4
        !       \
        !       3
        !----------------------------------------------------------------

        implicit none

        double precision,dimension(1:3),intent(in)::r1,r2,r3,r4

        real(8) :: im

        !Local
        real(kind=8),parameter          ::  pi =  4*atan(1.0_8) 
        real(8)            :: X1,Y1,Z1,&
                              X2,Y2,Z2,&
                              X3,Y3,Z3,&
                              X4,Y4,Z4
        real(8) :: v1,v1x,v1y,v1z,&
                      v2x,v2y,v2z,&
                   vn,vnx,vny,vnz
        integer::i

       !Interface: vectors to atom (in this case x,y,z)
        X1 = r1(1)
        X2 = r2(1)
        X3 = r3(1)
        X4 = r4(1)
!
        Y1 = r1(2)
        Y2 = r2(2)
        Y3 = r3(2)
        Y4 = r4(2)
!
        Z1 = r1(3)
        Z2 = r2(3)
        Z3 = r3(3)
        Z4 = r4(3)
!

       !Vectors 4-2 y 4-3
        v1x = X2 - X4
        v1y = Y2 - Y4
        v1z = Z2 - Z4
        v2x = X3 - X4
        v2y = Y3 - Y4
        v2z = Z3 - Z4
        !Vector normal al plano
        vnx = v1y*v2z - v1z*v2y
        vny = v1z*v2x - v1x*v2z
        vnz = v1x*v2y - v1y*v2x
        !Se normaliza
        vn = dsqrt(vnx**2+vny**2+vnz**2)
        vnx=vnx/vn
        vny=vny/vn
        vnz=vnz/vn

        !Se calcula y normaliza 4-1 (en v1)
        v1x = X1 - X4
        v1y = Y1 - Y4
        v1z = Z1 - Z4
        v1 = dsqrt(v1x**2+v1y**2+v1z**2)
        v1x=v1x/v1
        v1y=v1y/v1
        v1z=v1z/v1

        !Se calcula el producto escalar de vn y v1
        im = v1x*vnx + v1y*vny + v1z*vnz
        im = acos(im)

        !El ángulo es el complementario del anterior
        im = PI/2.d0 - im

        return

    end function calc_improper
!
!======================================================================!
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
       end module genforcefield
!
!======================================================================!

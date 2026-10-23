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
                   refine_stereotypes,                                &
                   genffbonded,                                        &
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
               match = VF2_Equivalent(nat,nat,adj,adj,ideg,ideg,order,i,j)
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
! REFINE_STEREOTYPES - split topological classes by cis/trans environment
!
! Refine existing types using planar relations across rigid nonaromatic
! bonds. All substituents participate; graph refinement carries distinctions
! into their branches without measuring torsions about single bonds.
       subroutine refine_stereotypes(nat,adj,lrigid,laroma,coord,    &
                                     initial,eqv,changed,debug)
       implicit none
       integer,intent(in) :: nat,initial(nat)
       logical,intent(in) :: adj(nat,nat),lrigid(nat,nat)
       logical,intent(in) :: laroma(nat,nat),debug
       real(kind=8),intent(in) :: coord(3,nat) ! Angstrom
       logical,intent(inout) :: eqv(nat,nat)
       logical,intent(out) :: changed
       integer,allocatable :: relations(:,:),signature(:,:),tuples(:,:)
       integer :: degree(nat),colors(nat),nextcolors(nat),nrel(nat)
       integer :: i,j,a,d,k,q,nq,start,n,ntypes,oldtypes,iteration
       integer :: temporary(4),side
       logical :: blocked(nat),valid
       real(kind=8) :: axis(3),u(3),v(3),pu(3),pv(3),length,cosphi
! New cis/trans window in degrees; compatible with the existing OOP
! planarity criterion. The roundoff guard below is dimensionless.
       real(kind=8),parameter :: planar_window_deg = 25.d0
       real(kind=8),parameter :: pi = 4*atan(1.d0)
       real(kind=8),parameter :: roundoff_guard = 64*epsilon(1.d0)
       real(kind=8),parameter :: planar_cos =                     &
                                  cos(planar_window_deg*pi/180.d0)

       changed = .false.
       if (nat < 4) return
       degree = count(adj,dim=1)
! Each degree-2/3 center has at most two substituents per eligible bond:
! four cross-bond pairs, each stored in both directions (eight per edge).
       allocate(relations(5,4*count(adj)))
       nq = 0
       blocked = .false.
       do i = 1,nat-1
         do j = i+1,nat
           if (.not. adj(i,j) .or. .not. lrigid(i,j)) cycle
           if (laroma(i,j)) cycle
           if (degree(i)<2 .or. degree(i)>3) cycle
           if (degree(j)<2 .or. degree(j)>3) cycle
           start = nq
           axis = coord(:,j)-coord(:,i)
           length = sqrt(sum(axis**2))
           valid = length > tiny(1.d0)
           if (valid) axis = axis/length
           do a = 1,nat
             if (.not. adj(a,i) .or. a==j) cycle
             do d = 1,nat
               if (.not. adj(d,j) .or. d==i) cycle
               if (a==d .or. .not. valid) then
                 valid = .false.
                 cycle
               end if
               u = coord(:,a)-coord(:,i)
               v = coord(:,d)-coord(:,j)
               pu = u-dot_product(u,axis)*axis
               pv = v-dot_product(v,axis)*axis
               if (sum(pu**2)<=roundoff_guard*sum(u**2) .or.       &
                   sum(pv**2)<=roundoff_guard*sum(v**2) .or.       &
                   sum(u**2)<=tiny(1.d0) .or.                     &
                   sum(v**2)<=tiny(1.d0)) then
                 valid = .false.
                 cycle
               end if
               cosphi = dot_product(pu,pv)/                      &
                                  sqrt(sum(pu**2)*sum(pv**2))
               if (abs(cosphi)<planar_cos) then
                 valid = .false.
                 cycle
               end if
               side = 1
               if (cosphi<0.d0) side = 2
               relations(:,nq+1) = [a,i,j,d,side]
               relations(:,nq+2) = [d,j,i,a,side]
               nq = nq+2
             end do
           end do
           if (.not. valid) then
! Do not split a class using partial/ambiguous geometric information.
             nq = start
             do a = 1,nat
               if ((adj(a,i).and.a/=j).or.(adj(a,j).and.a/=i))     &
                 blocked(initial(a)) = .true.
             end do
             write(*,'(A,I0,A,I0,A)') 'WARNING: stereo typing at bond ', &
               i,'-',j,': ambiguous geometry; affected types preserved.'
           end if
         end do
       end do
       if (nq==0) return
       nrel = 0
       do q = 1,nq
         a = relations(1,q)
         nrel(a) = nrel(a)+1
       end do
       allocate(tuples(4,maxval(nrel)))
       allocate(signature(1+nat+4*maxval(nrel),nat))
       colors = initial
       oldtypes = maxval(colors)
       do iteration = 1,nat
         signature = 0
         do a = 1,nat
           signature(1,a) = colors(a)
           if (blocked(initial(a))) cycle
! First split only by stereochemistry. Later rounds propagate the split
! through connectivity and reevaluate stereo signatures until stable.
           if (iteration>1) then
             do i = 1,nat
               if (adj(a,i)) signature(1+colors(i),a) =           &
                                          signature(1+colors(i),a)+1
             end do
           end if
           n = 0
           do q = 1,nq
             if (relations(1,q)/=a) cycle
             n = n+1
             tuples(:,n) = [colors(relations(2,q)),              &
               colors(relations(3,q)),colors(relations(4,q)),    &
               relations(5,q)]
           end do
! Sort tuples lexicographically, retaining multiplicity and avoiding an
! arbitrary reference substituent or atom-index-dependent orientation.
           do i = 2,n
             temporary = tuples(:,i)
             j = i-1
             do while (j>=1)
               do k = 1,4
                 if (temporary(k)/=tuples(k,j)) exit
               end do
               if (k>4) exit
               if (temporary(k)>tuples(k,j)) exit
               tuples(:,j+1) = tuples(:,j)
               j = j-1
             end do
             tuples(:,j+1) = temporary
           end do
           if (n>0) signature(nat+2:nat+1+4*n,a) =              &
                                             reshape(tuples(:,:n),[4*n])
         end do
         ntypes = 0
         do i = 1,nat
           do j = 1,i-1
             if (all(signature(:,i)==signature(:,j))) exit
           end do
           if (j<i) then
             nextcolors(i) = nextcolors(j)
           else
             ntypes = ntypes+1
             nextcolors(i) = ntypes
           end if
         end do
         if (ntypes==oldtypes) exit
         changed = .true.
         colors = nextcolors
         oldtypes = ntypes
       end do
       if (.not. changed) return
! eqvatms may represent a class by a sparse star. Rebuild its complete
! partition, otherwise removing the star center could lose valid matches.
       do i = 1,nat
         do j = 1,nat
           eqv(i,j) = colors(i)==colors(j) .and. i/=j
         end do
       end do
       if (debug) write(*,'(A,I0,A,I0)') 'Stereo atom types: ',     &
                                     maxval(initial),' -> ',oldtypes
       end subroutine refine_stereotypes

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
       integer,dimension(nat)                            ::  nmol     !
       integer,dimension(nat)                            ::  nagg     !
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
       subroutine genffbonded(nat,idat,nidat,coord,wiberg,adj,ideg,    &
                              lcycle,lrigid,laroma,latar,lch3,ich3,    &
                              znum,bonded,dihed,iroute,fring,fhetero,  &
                              fdouble,fnitro,fideal,debug)
!
       use datatypes, only: grobonded,                                 &
                            dihedrals
!
       implicit none
!
! Input/output variables
!
       type(grobonded),intent(inout)                   ::  bonded    !
       type(dihedrals),intent(inout)                   ::  dihed     !
       real(kind=8),dimension(3,nat),intent(in)        ::  coord     !  Atomic coordinates
       real(kind=8),dimension(nat,nat),intent(in)      ::  wiberg    !  Wiberg bond indices
       integer,dimension(:,:),allocatable,intent(out)  ::  ich3      !
       integer,dimension(nat),intent(in)               ::  idat      !
       logical,dimension(nat,nat),intent(in)           ::  adj       !  Boolean adjacency
       logical,dimension(nat,nat),intent(in)           ::  lcycle    !  Bonds belonging to rings
       logical,dimension(nat,nat),intent(in)           ::  lrigid    !
       logical,dimension(nat,nat),intent(in)           ::  laroma    !  Aromatic bonds
       logical,dimension(nat),intent(in)               ::  latar     !  Aromatic atoms
       logical,dimension(:),allocatable,intent(out)    ::  lch3      !
!
       integer,dimension(nat),intent(in)               ::  znum      !
       integer,dimension(nat),intent(in)               ::  ideg      !
       integer,intent(in)                              ::  nidat     !
       integer,intent(in)                              ::  nat       !  Number of atoms
!
       integer,intent(in)                              ::  iroute    !
       logical,intent(in)                              ::  fring     !
       logical,intent(in)                              ::  fhetero   !
       logical,intent(in)                              ::  fdouble   !
       logical,intent(in)                              ::  fnitro    !
       logical,intent(in)                              ::  fideal    !
       logical,intent(in)                              ::  debug     !  Debug mode
!
! Local variables
!
       logical,dimension(:,:),allocatable              ::  adjbond   !
       logical,dimension(:,:),allocatable              ::  adjang    !
       integer,dimension(:,:),allocatable              ::  edgeang   !
!
! Generating bonded terms
! -----------------------
!
      call genbonds(nat,coord,adj,bonded%nbond,bonded%ibond,           &
                    bonded%bond,bonded%kbond,bonded%fbond,             &
                    adjbond,debug)
!
      if ( nat .ge. 3 ) then
        call genangles(nat,coord,bonded%nbond,bonded%ibond,adjbond,    &
                       bonded%nang,bonded%iang,bonded%ang,bonded%kang, &
                       bonded%fang,adjang,edgeang,debug)
       end if
!
       if ( nat .ge. 4 ) then
         call gendihe(nat,adj,ideg,idat,coord,wiberg,znum,lcycle,      &
                      lrigid,laroma,latar,lch3,ich3,bonded%nbond,      &
                      bonded%ibond,bonded%nang,bonded%iang,adjang,     &
                      edgeang,dihed,iroute,fring,fhetero,fdouble,      &
                      fnitro,fideal,debug)
       end if
!
! Diatomics have bonds but never allocate the angle representation.
       if ( allocated(adjbond) ) deallocate(adjbond)
       if ( allocated(adjang) ) deallocate(adjang)
       if ( allocated(edgeang) ) deallocate(edgeang)
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
       logical,dimension(nat,nat),intent(in)              ::  adj      !  Boolean adjacency
       logical,dimension(:,:),allocatable,intent(out)     ::  adjbond  !  Boolean adjacency
       real(kind=8),dimension(3,nat),intent(in)           ::  coord    !  Atomic coordinates
       integer,intent(in)                                 ::  nat      !  Number of atoms
!
       logical,intent(in)                                 ::  debug    !  Debug mode
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
!
         ang(i) = calc_angle(coord(:,iang(1,i)),coord(:,iang(2,i)),    &
                                                     coord(:,iang(3,i)))
         ang(i) = ang(i)*180.0d0/pi
!
         if ( ang(i) .gt. 170 ) ang(i) = 180.0d0 ! TODO: add ideal value option
!
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
       subroutine gendihe(nat,adj,ideg,idat,coord,wiberg,znum,lcycle,  &
                          lrigid,laroma,latar,lch3,ich3,nbond,ibond,   &
                          nang,iang,adjang,edgeang,dihed,iroute,fring, &
                          fhetero,fdouble,fnitro,fideal,debug)
!
       use datatypes, only: dihedrals
       use graphtools
!
       implicit none
!
! Input/output variables
!
       type(dihedrals),intent(out)                     ::  dihed    !  Equilibrium dihedral angle
       real(kind=8),dimension(3,nat),intent(in)        ::  coord    !  Atomic coordinates
       real(kind=8),dimension(nat,nat),intent(in)      ::  wiberg   !  Wiberg bond indices
       logical,dimension(nat,nat),intent(in)           ::  adj      !
       logical,dimension(nat,nat),intent(in)           ::  lcycle   !  Bonds belonging to rings
       logical,dimension(nat,nat),intent(in)           ::  lrigid   !  Rigid bonds
       logical,dimension(nat,nat),intent(in)           ::  laroma   !  Aromatic bonds
       logical,dimension(nat),intent(in)               ::  latar    !  Aromatic atoms
       logical,dimension(:),allocatable,intent(out)    ::  lch3     !
       integer,dimension(:,:),allocatable,intent(out)  ::  ich3     !
       integer,dimension(nat),intent(in)               ::  idat     !
       integer,dimension(nat),intent(in)               ::  znum     !
       integer,dimension(nat),intent(in)               ::  ideg     !
       integer,intent(in)                              ::  nat      !  Number of atoms
!
       logical,dimension(nang,nang),intent(in)         ::  adjang   !  Boolean adjacency
       integer,dimension(2,nang),intent(in)            ::  edgeang  !
       integer,dimension(2,nbond),intent(in)           ::  ibond    !
       integer,dimension(3,nang),intent(in)            ::  iang     !
       integer,intent(in)                              ::  nang     !  Number of angles
       integer,intent(in)                              ::  nbond    !  Number of bonds
!
       integer,intent(in)                              ::  iroute   !
       logical,intent(in)                              ::  fring    !
       logical,intent(in)                              ::  fhetero  !
       logical,intent(in)                              ::  fdouble  !
       logical,intent(in)                              ::  fnitro !
       logical,intent(in)                              ::  fideal   !
       logical,intent(in)                              ::  debug    !  Debug mode
!
! Local variables
!
       integer,dimension(:,:),allocatable              ::  edges
       integer                                         ::  ndihe    !  Number of dihedrals
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
       allocate(lch3(ndihe),ich3(3,ndihe))
!
       allocate(edges(2,ndihe))
!
       call adj2edge(nang,adjang,ndihe,edges)
!
! Classifying adjacent angles as proper or improper dihedrals
!
       call setdihe(nat,adj,ideg,idat,coord,wiberg,nbond,ibond,nang,  &
                    edgeang,iang,ndihe,dihed,edges,lcycle,lrigid,      &
                    laroma,latar,lch3,ich3,znum,iroute,fring,fhetero,  &
                    fdouble,fnitro,fideal,debug)
!
       deallocate(edges)
!
       return
       end subroutine gendihe
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
       subroutine setdihe(nat,adj,ideg,idat,coord,wiberg,nbond,ibond, &
                          nang,edgeang,iang,ndihe,dihed,edgedihe,      &
                          lcycle,lrigid,laroma,latar,lch3,ich3,znum,   &
                          iroute,fring,fhetero,fdouble,fnitro,      &
                          fideal,debug)
!
       use datatypes,   only: dihedrals
       use genfftools
       use quadruplets, only: genquad
!
       implicit none
!
! Input/output variables
!
       type(dihedrals),intent(inout)              ::  dihed     !
       real(kind=8),dimension(3,nat),intent(in)   ::  coord     !
       real(kind=8),dimension(nat,nat),intent(in) ::  wiberg    !  Wiberg bond indices
       logical,dimension(nat,nat),intent(in)      ::  adj       !
       logical,dimension(nat,nat),intent(in)      ::  lcycle    !  Bonds belonging to rings
       logical,dimension(nat,nat),intent(in)      ::  lrigid    !  Bonds belonging to rings
       logical,dimension(nat,nat),intent(in)      ::  laroma    !  Aromatic bonds
       logical,dimension(nat),intent(in)          ::  latar     !  Aromatic atoms
       logical,dimension(ndihe),intent(out)       ::  lch3      !
       integer,dimension(3,ndihe),intent(out)     ::  ich3      !
       integer,dimension(nat),intent(in)          ::  idat      !
       integer,dimension(nat),intent(in)          ::  znum      !
       integer,dimension(nat),intent(in)          ::  ideg      !
       integer,intent(in)                         ::  nat       !
       logical,intent(in)                         ::  debug     !
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
       logical,intent(in)                         ::  fring     !
       logical,intent(in)                         ::  fhetero   !
       logical,intent(in)                         ::  fdouble   !
       logical,intent(in)                         ::  fnitro    !
       logical,intent(in)                         ::  fideal    !
!
! Local variables
!
       real(kind=8),dimension(ndihe)              ::  dimpro    !
       integer,dimension(4,ndihe)                 ::  iimpro    !
       integer,dimension(ndihe)                   ::  fimpro    !
       integer,dimension(ndihe)                   ::  primpro   !
       integer,dimension(ndihe)                   ::  radeq     ! Reference adequecy
       real(kind=8),dimension(ndihe)              ::  rwiberg   ! Reference Wiberg
       integer                                    ::  nimpro    !
       real(kind=8),dimension(ndihe)              ::  dinv      !
       integer,dimension(4,ndihe)                 ::  iinv      !
       integer,dimension(ndihe)                   ::  finv      !
       integer                                    ::  ninv      !
       logical,dimension(ndihe)                   ::  torsion   !
       logical,dimension(4)                       ::  lcheck    !
       logical                                    ::  flag      !
       logical                                    ::  laromatic !
       logical                                    ::  lcarbonyl !
       logical                                    ::  lconjhet  !
       logical                                    ::  lalkene   !
       logical                                    ::  lnitrogen !
       logical                                    ::  lfound    !
       real(kind=8)                               ::  daux      !
       integer,dimension(4)                       ::  vaux      !
       integer,dimension(2,4)                     ::  bonds     !
       integer,dimension(2)                       ::  rbond     !
       integer,dimension(3)                       ::  ang1      !
       integer,dimension(3)                       ::  ang2      !
       integer,dimension(4)                       ::  isigna    ! Target improper signature
       integer                                    ::  iprio     ! Improper priority
       integer                                    ::  iradeq    !
       integer                                    ::  ncandi    ! Number of improper candidates
       integer                                    ::  irepcandi ! Representative candidate
       integer                                    ::  iselcandi ! Selected candidate
       integer                                    ::  icenter   ! Center atom
       integer                                    ::  itmp      !
       integer                                    ::  i,j,k     !  Indexes
!
       real(kind=8),parameter                     ::  pi =  4*atan(1.0_8)
       real(kind=8),parameter                     ::  linear_angle_thr = 170.0d0
!
! Finding atoms belonging to each dihedral
! ----------------------------------------
!
! Start by assuming every candidate is a proper torsion.  Geminal-angle
! candidates below are reclassified as OOP impropers or inversions and marked
! as non-torsional.
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
! OOP impropers and inversion coordinates are generated only for three-fold
! centers.  Other geminal angle combinations are ignored at this stage.
!
           if ( ideg(vaux(1)) .ne. 3 ) cycle
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
! Canonical terminal ordering makes the selected representative
! independent of adjacent-angle enumeration order
!
! The first two positions of an improper or inversion are an oriented
! center-reference bond.  setdeps relies on this invariant to associate
! aromatic rigid terms with the corresponding improper.
!
           if ( (idat(vaux(3)).gt.idat(vaux(4))) .or.                  &
                ((idat(vaux(3)).eq.idat(vaux(4))).and.                 &
                 (vaux(3).gt.vaux(4))) ) then
             itmp = vaux(3)
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
! Classify the chemical reason for adding a star-like improper.  The class
! priority is later used to choose one representative when several equivalent
! OOP candidates are possible for the same atom type.
!
             lconjhet = ((znum(vaux(2)).eq.7)                          &
                         .or. (znum(vaux(2)).eq.8)                     &
                         .or. (znum(vaux(2)).eq.16))                   &
                        .and. (lrigid(vaux(1),vaux(2))                 &
                             .or. (wiberg(vaux(1),vaux(2)).ge. 1.20d0))
!
             lcarbonyl = (znum(vaux(1)).eq.6)                          &
                          .and. (.NOT.lcycle(vaux(1),vaux(2)))         &
                          .and. lconjhet
!
             lalkene = (znum(vaux(1)).eq.6)                            &
                        .and. (.NOT.lcycle(vaux(1),vaux(2)))           &
                        .and. lrigid(vaux(1),vaux(2))                  &
                        .and. (znum(vaux(2)).eq.6)
!
             laromatic = (znum(vaux(1)).eq.6)                          &
                          .and. (.NOT.lcycle(vaux(1),vaux(2)))         &
                          .and. (.NOT.lrigid(vaux(1),vaux(2)))         &
                          .and. (lcycle(vaux(1),vaux(3))               &
                                 .or. laroma(vaux(1),vaux(3)))         &
                          .and. (lcycle(vaux(1),vaux(4))               &
                                 .or. laroma(vaux(1),vaux(4)))
!
             lnitrogen = (znum(vaux(1)).eq.7)                          &
                          .and. (lrigid(vaux(1),vaux(2))               &
                                 .or. lrigid(vaux(1),vaux(3))          &
                                 .or. lrigid(vaux(1),vaux(4))          &
                                 .or. latar(vaux(2))                   &
                                 .or. latar(vaux(3))                   &
                                 .or. latar(vaux(4))                   &
                                 .or. laroma(vaux(1),vaux(2))          &
                                 .or. laroma(vaux(1),vaux(3))          &
                                 .or. laroma(vaux(1),vaux(4)))
!
             flag   = .FALSE.
             iprio  = 0
             iradeq = 1
!
! fhetero, fdouble, fring and fnitro are user-facing switches.  They gate the
! individual OOP families without changing the underlying candidate search.
!
             if ( fhetero .and. lcarbonyl ) then
               flag = .TRUE.
               iprio = 1
               iradeq = 0
             else if ( fdouble .and. lalkene ) then
               flag = .TRUE.
               iprio = 2
               iradeq = 0
             else if ( fring .and. laromatic ) then
               flag = .TRUE.
               iprio = 3
               iradeq = 0
             else if ( fnitro .and. lnitrogen ) then
               flag = .TRUE.
               iprio = 4
!
! Prefer the conjugated N neighbour as the oriented reference.
! A candidate with another reference is retained only as a fallback
!
               if ( lrigid(vaux(1),vaux(2)) .or. latar(vaux(2))       &
                    .or. laroma(vaux(1),vaux(2)) ) iradeq = 0
             end if
!
             if ( flag ) then
!
! Store the candidate as an OOP improper.  The final representative selection
! happens after the full candidate list is known.
!
               dihed%nimpro = dihed%nimpro + 1
!
               dihed%iimpro(:,dihed%nimpro) = vaux(:)
               dihed%fimpro(dihed%nimpro)   = 2
!
               if ( fideal ) then
! Idealized stiff/OOP values use 0 degrees for planar impropers.
!
                 dihed%dimpro(dihed%nimpro) = 0.0d0
               else
! The --noideal route preserves the actual equilibrium geometry.
!
                 dihed%dimpro(dihed%nimpro) = daux
               end if
!
               primpro(dihed%nimpro) = iprio
               radeq(dihed%nimpro)   = iradeq
               rwiberg(dihed%nimpro) = wiberg(vaux(1),vaux(2))
!
             end if
!
! If dihedral angle deviates from planarity classify it as inversion
!
           else if ( znum(vaux(1)) .eq. 7 ) then
!
! Non-planar three-fold nitrogen centers are treated as inversion coordinates
! rather than planar OOP impropers.
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
! Keep one improper per physical center, with a common oriented chemical
! signature for all centers sharing the same atom type.
! The representative is selected with the existing chemical/Wiberg
! comparator.
! Each equivalent center must then reproduce its idat(1:4) signature
!
       ncandi = dihed%nimpro
       nimpro = 0
!
! Find the best representative for each central atom type
!
       do i = 1, ncandi
         flag = .FALSE.
! Skip atom-type classes that were already processed.  The first candidate of
! each central atom type starts a new representative search.
!
         do j = 1, i-1
           if ( idat(dihed%iimpro(1,i)) .eq.                           &
                                          idat(dihed%iimpro(1,j)) ) then
             flag = .TRUE.
             exit
           end if
         end do
         if ( flag ) cycle
!
         irepcandi = i
         do j = i+1, ncandi
           if ( idat(dihed%iimpro(1,j)) .eq.                           &
                                          idat(dihed%iimpro(1,i)) ) then
! Find the chemically best signature for this central atom type.  The selected
! signature will be required for every equivalent physical center.
!
             call improcandi(idat,znum,dihed%iimpro(:,j),primpro(j),   &
                             radeq(j),rwiberg(j),                      &
                             dihed%iimpro(:,irepcandi),                &
                             primpro(irepcandi),radeq(irepcandi),      &
                             rwiberg(irepcandi),lfound)
             if ( lfound ) then
               irepcandi = j
             end if
           end if
         end do

         isigna(:) = idat(dihed%iimpro(:,irepcandi))
!
! Select the matching signature for every physical center of this type
!
         do j = i, ncandi
           if ( idat(dihed%iimpro(1,j)) .ne.                         &
                idat(dihed%iimpro(1,i)) ) cycle

           icenter = dihed%iimpro(1,j)
           flag = .FALSE.
! Avoid writing two impropers for the same physical center within this atom-type
! class.
!
           do k = i, j-1
             if ( (idat(dihed%iimpro(1,k)).eq.idat(dihed%iimpro(1,i))) &
                             .and. (dihed%iimpro(1,k).eq.icenter) ) then
               flag = .TRUE.
               exit
             end if
           end do
           if ( flag ) cycle

           iselcandi = 0
           do k = 1, ncandi
             if ( dihed%iimpro(1,k) .ne. icenter ) cycle
             if ( .not.all(idat(dihed%iimpro(:,k)).eq.isigna(:)) ) cycle

             if ( iselcandi .eq. 0 ) then
               iselcandi = k
             else
! If several candidates have the required signature on the same center, choose
! the best-oriented one using the same comparator.
!
               call improcandi(idat,znum,dihed%iimpro(:,k),            &
                               primpro(k),radeq(k),rwiberg(k),         &
                               dihed%iimpro(:,iselcandi),              &
                               primpro(iselcandi),radeq(iselcandi),    &
                               rwiberg(iselcandi),lfound)
               if ( lfound ) iselcandi = k
             end if
           end do

           if ( iselcandi .eq. 0 ) then
! This is a genuine consistency error: an equivalent center lacks the chemical
! signature chosen as representative for its atom-type class.
!
             write(*,'(A,I0,A,4(1X,I0))')                              &
                  'ERROR in setdihe: no OOP candidate for center ',    &
                  icenter,' matches chemical signature',isigna(:)
             stop 1
           end if

           nimpro           = nimpro + 1
           iimpro(:,nimpro) = dihed%iimpro(:,iselcandi)
           dimpro(nimpro)   = dihed%dimpro(iselcandi)
           fimpro(nimpro)   = dihed%fimpro(iselcandi)
           primpro(nimpro)  = primpro(iselcandi)
           radeq(nimpro)    = radeq(iselcandi)
           rwiberg(nimpro)  = rwiberg(iselcandi)
         end do
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
! Multiple inversion candidates can be generated for the same center.  Keep the
! best reference atom so Joyce receives a single inversion coordinate.
!
             call invcandi(idat,znum,wiberg,dihed%iinv(:,i),           &
                           iinv(:,j),lfound)
!
             if ( lfound ) then
               iinv(:,j) = dihed%iinv(:,i)
               dinv(j)   = dihed%dinv(i)
               finv(j)   = dihed%finv(i)
             end if
!
             exit
           end if
         end do
         if ( .not. flag ) then
!
           ninv = ninv + 1
!
           iinv(:,ninv) = dihed%iinv(:,i)
           dinv(ninv)   = dihed%dinv(i)
           finv(ninv)   = dihed%finv(i)
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
! Remaining candidates are proper dihedrals.  They are reconstructed from two
! adjacent angles by finding the common central bond and the two terminal atoms.
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
!
! The duplicated bond among the four angle edges is the central bond of the
! proper dihedral.
!
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
!
! The two non-central bonds provide the terminal atoms of the torsion.
!
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
! Sorting dihedral indexes
!
! Orient the central bond deterministically by atom type.  This keeps labels,
! equivalence detection and dependency matching reproducible.
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
! Checking if there are 3 atoms forming an angle of 180 deg
!
! Skip near-linear torsions because their dihedral angle is numerically ill
! defined and not suitable as an internal coordinate.
!
           daux = calc_angle(coord(:,vaux(1)),coord(:,vaux(2)),        &
                                                       coord(:,vaux(3)))
           daux = daux*180.0d0/pi
           if ( abs(daux) .ge. linear_angle_thr ) cycle
!
           daux = calc_angle(coord(:,vaux(2)),coord(:,vaux(3)),        &
                                                       coord(:,vaux(4)))
           daux = daux*180.0d0/pi
           if ( abs(daux) .ge. linear_angle_thr ) cycle
!
! Computing equilibrium value
!
           call Diedro(coord(:,vaux(1)),coord(:,vaux(2)),              &
                       coord(:,vaux(3)),coord(:,vaux(4)),daux)
           daux = daux*180.0d0/pi
!
! If dihedral angle is rigid classify it as improper
!
           if ( lrigid(vaux(2),vaux(3)) ) then  ! TODO: Removing 3-member cycle dihedrals
!
! Rigid central bonds are stiff dihedrals.  They are included in Joyce Step1
! rather than relaxed-scan torsional fitting.
!
             dihed%nrigid = dihed%nrigid + 1
!
             dihed%irigid(:,dihed%nrigid) = vaux(:)
             dihed%frigid(dihed%nrigid)   = 2
             if ( fideal ) then
! Idealize stiff torsions to the nearest planar/trans value.
!
               if ( abs(daux) .le. 90.0d0 ) then
                 dihed%drigid(dihed%nrigid) = 0.0d0
               else
                 dihed%drigid(dihed%nrigid) = 180.0d0
               end if
             else
               dihed%drigid(dihed%nrigid) = daux
             end if
!
! If dihedral angle belongs to a cycle classify it as inversion
!
           else if ( lcycle(vaux(2),vaux(3)) ) then
!
! Ring torsions are not scanned independently; they are treated as inversion-
! like stiff coordinates because the ring geometry constrains them strongly.
!
             dihed%ninv = dihed%ninv + 1
!
             dihed%iinv(:,dihed%ninv) = vaux(:)
             dihed%dinv(dihed%ninv)   = daux  !  TODO: how to handle equilibrium value
             dihed%finv(dihed%ninv)   = 3
!
           else
!
! Only acyclic non-rigid central bonds become flexible torsions for relaxed QM
! scans and Fourier fitting.
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
! genquad chooses which flexible torsions are independent scan coordinates and
! records dependencies for additional ring-exocyclic terms.
!
       call genquad(nat,coord,adj,ideg,lcycle,lrigid,znum,dihed,ndihe, &
                    debug)
!
! genmethlist marks methyl rotors before Fourier terms are allocated so they
! can use a compact n=3 functional form.
!
       call genmethlist(nat,adj,znum,dihed,dihed%ndihe,lch3,ich3,debug)
!
! Generating Fourier series
!
       do i = 1, dihed%nflexi
!
         if ( lch3(i) ) then
!
! Methyl rotors are dominated by threefold periodicity; avoid adding a full
! seven-term Fourier series for these simple rotations.
!
           dihed%flexi(i)%ntor = 1
           allocate(dihed%flexi(i)%tor(dihed%flexi(i)%ntor))
!
           dihed%flexi(i)%itor(:) = dihed%iflexi(:,i)
           dihed%fflexi(i)        = 1
!
           dihed%flexi(i)%tor(1)%vtor  = 0.0d0
           dihed%flexi(i)%tor(1)%phase = dihed%dflexi(i)
           dihed%flexi(i)%tor(1)%multi = 3
!
         else
!
! General flexible torsions start with a complete n=0..6 Fourier expansion.
! The fitted topology or user editing can later remove unnecessary terms.
!
           dihed%flexi(i)%ntor = 7 ! TODO: choose number of functions in Fourier expansion
           allocate(dihed%flexi(i)%tor(dihed%flexi(i)%ntor))
!
           dihed%flexi(i)%itor(:) = dihed%iflexi(:,i)
           dihed%fflexi(i)        = 1
!
           do j = 1, dihed%flexi(i)%ntor
             dihed%flexi(i)%tor(j)%vtor  = 0.0d0
             dihed%flexi(i)%tor(j)%phase = dihed%dflexi(i) ! TODO: shift phase for even multiplicities
             dihed%flexi(i)%tor(j)%multi = j-1
           end do
!
         end if
!
       end do
!
       return
       end subroutine setdihe
!
!======================================================================!
!
       subroutine improcandi(idat,znum,candidate,cclass,cadeq,cwiberg, &
                             selected,sclass,sadeq,swiberg,better)
!
       implicit none
!
! Input/output variables
!
       integer,dimension(:),intent(in)  ::  idat       !
       integer,dimension(:),intent(in)  ::  znum       !
       integer,dimension(4),intent(in)  ::  candidate  !
       integer,dimension(4),intent(in)  ::  selected   !
       integer,intent(in)               ::  cclass     ! Candidate class
       integer,intent(in)               ::  sclass     ! Selected class
       integer,intent(in)               ::  cadeq      ! Candidate adequacy
       integer,intent(in)               ::  sadeq      ! Selected adequacy
       real(kind=8),intent(in)          ::  cwiberg    ! Candidate Wiberg
       real(kind=8),intent(in)          ::  swiberg    ! Selected Wiberg
       logical,intent(out)              ::  better     !
!
       better = .FALSE.
!
! Candidate ordering is deliberately hierarchical:
!  1) chemical class
!  2) adequacy of the oriented reference atom
!  3) Wiberg bond strength
!  4) heavy atom over H
!  5) element priority
!  6) deterministic atom-type/index signature
!
       if ( cclass .ne. sclass ) then
         better = cclass .lt. sclass
       else if ( cadeq .ne. sadeq ) then
         better = cadeq .lt. sadeq
       else if ( abs(cwiberg-swiberg) .gt. 1.0d-12 ) then
         better = cwiberg .gt. swiberg
       else if ( (znum(candidate(2)).ne.1) .neqv.                      &
                                         (znum(selected(2)).ne.1) ) then
         better = znum(candidate(2)) .ne. 1
       else if ( refeleprio(znum(candidate(2))) .ne. &
                                    refeleprio(znum(selected(2))) ) then
         better = refeleprio(znum(candidate(2))) .lt. &
                                           refeleprio(znum(selected(2)))
       else
         call checksign(idat,candidate,selected,better)
       end if
!
       return
       end subroutine improcandi
!
!======================================================================!
!
       subroutine invcandi(idat,znum,wiberg,candidate,selected,better)
!
       implicit none
!
       integer,dimension(:),intent(in)       ::  idat,znum
       real(kind=8),dimension(:,:),intent(in) ::  wiberg
       integer,dimension(4),intent(in)        ::  candidate,selected
       logical,intent(out)                     ::  better
!
       better = .FALSE.
!
! For an amine inversion, H is a fallback reference: choose the
! strongest heavy-atom bond before considering any N--H candidate
!
       if ( (znum(candidate(2)).ne.1) .neqv.                          &
                 (znum(selected(2)).ne.1) ) then
         better = znum(candidate(2)) .ne. 1
       else if ( abs(wiberg(candidate(1),candidate(2)) -              &
                    wiberg(selected(1),selected(2))) .gt. 1.0d-12 ) then
         better = wiberg(candidate(1),candidate(2)) .gt.              &
                  wiberg(selected(1),selected(2))
       else if ( refeleprio(znum(candidate(2))) .ne. &
                                    refeleprio(znum(selected(2))) ) then
         better = refeleprio(znum(candidate(2))) .lt. &
                                           refeleprio(znum(selected(2)))
       else
         call checksign(idat,candidate,selected,better)
       end if
!
       return
       end subroutine invcandi
!
!======================================================================!
!
! REFELEPRIO - REFerence ELEment PRIOrity
!
! Lower values are preferred as oriented reference atoms when all more specific
! criteria are tied.  This keeps O/S/N/C references ordered reproducibly.
!
       integer function refeleprio(atomic_number)
!
       implicit none
!
       integer,intent(in) :: atomic_number
!
       select case (atomic_number)
       case (8)
         refeleprio = 1
       case (16)
         refeleprio = 2
       case (7)
         refeleprio = 3
       case (6)
         refeleprio = 4
       case default
         refeleprio = 5
       end select
!
       return
       end function refeleprio
!
!======================================================================!
!
       subroutine checksign(idat,candidate,selected,greater)
!
       implicit none
!
       integer,dimension(:),intent(in)  ::  idat
       integer,dimension(4),intent(in)  ::  candidate,selected
       logical,intent(out)              ::  greater
       integer                          ::  i
!
       greater = .FALSE.
!
! Tie-break by atom-type signature first and by raw atom index second.  This is
! not a chemical criterion; it only makes otherwise equivalent choices stable.
!
       do i = 1, 4
         if ( idat(candidate(i)) .ne. idat(selected(i)) ) then
           greater = idat(candidate(i)) .gt. idat(selected(i))
           return
         end if
       end do
       do i = 1, 4
         if ( candidate(i) .ne. selected(i) ) then
           greater = candidate(i) .gt. selected(i)
           return
         end if
       end do
!
       return
       end subroutine checksign
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
       end module genforcefield
!
!======================================================================!

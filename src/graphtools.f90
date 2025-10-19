!======================================================================!
!
       module graphtools
!
       implicit none
!
       contains
!
!======================================================================!
!
! BUILDADJ - BUILD ADJacency matrix
!
! This subroutine 
!
       subroutine buildadj(nat,coord,adj)
!
       implicit none
!
! Input/output variables
!
       real(kind=8),dimension(3,nat),intent(in)  ::  coord  !  Atomic coordinates
       logical,dimension(nat,nat),intent(out)    ::  adj    !  Adjacency matrix
       integer,intent(in)                        ::  nat    !  Number of nodes
!
! Local variables
!
       real(kind=8),dimension(3)                 ::  r      !  Relative position
       real(kind=8)                              ::  dist   !  Interatomic distance
       real(kind=8)                              ::  bond   !  Bond threshold
       integer                                   ::  i,j    !  Indexes
!
! Building adjacency matrix according to geometric criteria
!
       bond = 1.6d0
!
       adj(:,:) = .FALSE.
!
       do i = 1, nat-1
         do j = i+1, nat
           r    = coord(:,i)-coord(:,j)
           dist = sqrt(dot_product(r,r))
           if ( dist .le. bond ) then
             adj(i,j) = .TRUE.
             adj(j,i) = .TRUE.
           end if
         end do
       end do
!
       return
       end subroutine buildadj
!
!======================================================================!
!
! WIBERG2ADJ - WIBERG bond indexes TO ADJacency matrix
!
! This subroutine 
!
       subroutine wiberg2adj(nat,wiberg,adj)
!
       implicit none
!
! Input/output variables
!
       real(kind=8),dimension(nat,nat),intent(in)  ::  wiberg  !  Atomic coordinates
       logical,dimension(nat,nat),intent(out)      ::  adj     !  Adjacency matrix
       integer,intent(in)                          ::  nat     !  Number of nodes
!
! Local variables
!
       real(kind=8)                                ::  thr     !  Bond threshold
       integer                                     ::  i,j     !  Indexes
!
! Building adjacency matrix according to Wiberg bond indexes
!
       thr = 0.75d0
!
       adj(:,:) = .FALSE.
!
       do i = 1, nat-1
         do j = i+1, nat
           if ( wiberg(j,i) .ge. thr ) then
             adj(i,j) = .TRUE.
             adj(j,i) = .TRUE.
           end if
         end do
       end do
!
       return
       end subroutine wiberg2adj
!
!======================================================================!
!
! BONDS2ADJ - BOND termS TO ADJacency matrix
!
! This subroutine 
!
       subroutine bonds2adj(nbond,ibond,nat,adj)
!
       implicit none
!
! Input/output variables
!
       logical,dimension(nat,nat),intent(out)  ::  adj    !  Adjacency matrix
       integer,dimension(2,nbond),intent(in)   ::  ibond  !  Atom indexes
       integer,intent(in)                      ::  nbond  !  Number of edges
       integer,intent(in)                      ::  nat    !  Number of nodes
!
! Local variables
!
       integer                                     ::  i  !  Index
!
! Building adjacency matrix according to Wiberg bond indexes
!
       adj(:,:) = .FALSE.
!
       do i = 1, nbond
         adj(ibond(1,i),ibond(2,i)) = .TRUE.
         adj(ibond(2,i),ibond(1,i)) = .TRUE.
       end do
!
       return
       end subroutine bonds2adj
!
!======================================================================!
!
! COUNTEDGES - COUNT EDGES
!
! This subroutine 
!
       subroutine countedges(nat,adj,nbond)
!
       implicit none
!
! Input/output variables
!
       logical,dimension(nat,nat),intent(in)  ::  adj    !  Adjacency matrix
       integer,intent(in)                     ::  nat    !  Number of nodes
       integer,intent(out)                    ::  nbond  !  Number of edges
!
! Local variables
!
       integer                                ::  i,j   !  Indexes
!
! Counting number of edges
!
       nbond = 0
       do i = 1, nat-1
         do j = i+1, nat
           if ( adj(j,i) ) nbond = nbond + 1
         end do
       end do
!
       return
       end subroutine countedges
!
!======================================================================!
!
! ADJ2EDGE - transform from ADJacency TO EDGE
!
! This subroutine 
!
       subroutine adj2edge(nat,adj,nbond,edge)
!
       implicit none
!
! Input/output variables
!
       logical,dimension(nat,nat),intent(in)     ::  adj      !  Adjacency matrix
       integer,dimension(2,nbond),intent(out)    ::  edge     !  
       integer,intent(in)                        ::  nat      !  Number of nodes
       integer,intent(in)                        ::  nbond    !  Number of edges
!
! Local variables
!
       integer                                   ::  iibond   !  Index
       integer                                   ::  i,j      !  Indexes
!
! Obtaining edges information
!
       iibond = 0
       do i = 1, nat-1
         do j = i+1, nat
           if ( adj(j,i) ) then
!
             iibond = iibond + 1
! 
             edge(1,iibond) = i
             edge(2,iibond) = j
!
           end if
         end do
       end do
!
       return
       end subroutine adj2edge
!
!======================================================================!
!
! ADJ2INC - transform representation from ADJacency TO INCidence
!
! This subroutine 
!
       subroutine adj2inc(nat,adj,nbond,incbond)
!
       implicit none
!
! Input/output variables
!
       logical,dimension(nat,nat),intent(in)     ::  adj      !  Adjacency matrix
       logical,dimension(nbond,nat),intent(out)  ::  incbond  !  Incidence matrix
       integer,intent(in)                        ::  nat      !  Number of nodes
       integer,intent(in)                        ::  nbond    !  Number of edges
!
! Local variables
!
       integer                                   ::  iibond   !  Index
       integer                                   ::  i,j      !  Indexes
!
! Obtaining edges information
!
       incbond(:,:) = .FALSE.
!
       iibond = 0
       do i = 1, nat-1
         do j = i+1, nat
           if ( adj(j,i) ) then
!
             iibond = iibond + 1
!
             incbond(iibond,i) = .TRUE.
             incbond(iibond,j) = .TRUE.
!
           end if
         end do
       end do
!
       return
       end subroutine adj2inc
!
!======================================================================!
!
! INC2ADJ - transform representation from INCidence TO ADJacency
!
! This subroutine 
!
       subroutine inc2adj(nat,nbond,inc,edge,adj)
!
       implicit none
!
! Input/output variables
!
       logical,dimension(nbond,nbond),intent(out)  ::  adj    !  Adjacency matrix
       logical,dimension(nbond,nat),intent(in)     ::  inc    !  Incidence matrix
       integer,dimension(2,nbond),intent(in)       ::  edge   !  
       integer,intent(in)                          ::  nat    !  Number of nodes
       integer,intent(in)                          ::  nbond  !  Number of edges
!
! Local variables
!
       integer                                     ::  i,j,k  !  Indexes
!
! Building adjacency matrix in one representation from the incidence ma-
!  trix of the reference representation
!
       adj(:,:) = .FALSE.
!
       do i = 1, nbond-1
         do j = 1, 2
           do k = i+1, nbond
              if ( inc(k,edge(j,i)) ) then
                adj(i,k) = .TRUE.
                adj(k,i) = .TRUE.
              end if
            end do
         end do
       end do
!
       return
       end subroutine inc2adj
!
!======================================================================!
!
! STRNORMALLABELS - STRing NORMalize LABELS
!
! This subroutine 
!
       subroutine strnormlabels(n,lab,ilab,nlab)
!
       use lengths, only: lenlab
!
       implicit none
!
! Input/output variables
! 
       character(len=lenlab),dimension(n),intent(in)   ::  lab    !  
       integer,dimension(n),intent(out)                ::  ilab   !  
       integer,intent(out)                             ::  nlab   !  
       integer,intent(in)                              ::  n      !  
!
! Local variables
!
       logical                                         ::  match  !
       integer                                         ::  i,j    !  Indexes
!
! Relabeling nodes in the dictionary of values (1,#different labels)
!
       nlab = 0
       ilab(:) = 0
!
       do i = 1, n
         match = .FALSE.
         do j = 1, i-1
           if ( lab(i) .eq. lab(j) ) then
             match = .TRUE.
             ilab(i) = ilab(j)
             exit
           end if
         end do
         if ( .not. match ) then
           nlab    = nlab + 1
           ilab(i) = nlab
         end if
       end do
!
       return
       end subroutine strnormlabels
!
!======================================================================!
!
! INORMLABELS - Integer NORMalize LABELS
!
! This subroutine 
!
       subroutine inormlabels(nlab,inlab,outlab,mlab)
!
       implicit none
!
! Input/output variables
!
       integer,dimension(nlab),intent(in)   ::  inlab   !  
       integer,dimension(nlab),intent(out)  ::  outlab  !  
       integer,intent(in)                   ::  nlab    !  
       integer,intent(out)                  ::  mlab    !  
!
! Local variables
!
       logical                              ::  match   !
       integer                              ::  i,j     !  Indexes
!
! Relabeling nodes in the dictionary of values (1,#different labels)
!
       mlab      = 0
       outlab(:) = 0
!
       do i = 1, nlab
         match = .FALSE.
         do j = 1, i-1
           if ( inlab(i) .eq. inlab(j) ) then
             match = .TRUE.
             outlab(i) = outlab(j)
             exit
           end if
         end do
         if ( .not. match ) then
           mlab      = mlab + 1
           outlab(i) = mlab
         end if
       end do
!
       return
       end subroutine inormlabels
!
!======================================================================!
!
! DNORMLABELS - Double precision NORMalize LABELS
!
! This subroutine 
!
       subroutine dnormlabels(nlab,inlab,outlab,mlab)
!
       implicit none
!
! Input/output variables
!
       real(kind=8),dimension(nlab),intent(in)   ::  inlab   !  
       integer,dimension(nlab),intent(out)       ::  outlab  !  
       integer,intent(in)                        ::  nlab    !  
       integer,intent(out)                       ::  mlab    !  
!
! Local variables
!
       logical                                   ::  match   !
       integer                                   ::  i,j     !  Indexes
!
       real(kind=8),parameter                    ::  thr = 1.0E-6
!
! Relabeling nodes in the dictionary of values (1,#different labels)
!
       mlab      = 0
       outlab(:) = 0
!
       do i = 1, nlab
         match = .FALSE.
         do j = 1, i-1
           if ( (inlab(i)-inlab(j)) .lt. thr ) then
             match = .TRUE.
             outlab(i) = outlab(j)
             exit
           end if
         end do
         if ( .not. match ) then
           mlab      = mlab + 1
           outlab(i) = mlab
         end if
       end do
!
       return
       end subroutine dnormlabels
!
!======================================================================!
!
! FINDCOMPUNDIR - FIND COMPonents UNDIRected
!
! This subroutine finds the connected components in an undirected 
!  unweighted graph of NNODE vertices given as an adjacency matrix 
!  representaion ADJ(NODE,NODE) using Breadth First Search.  
!
       subroutine findcompundir(nnode,adj,imol,iagg,itag,maxagg,nmol,  &
                                nagg)
!
       implicit none
!
! Input/output variables
!
       logical,dimension(nnode,nnode),intent(in)  ::  adj     !  Adjacency matrix
       integer,dimension(nnode),intent(out)       ::  imol    !  Molecules identifier
       integer,dimension(nnode),intent(out)       ::  iagg    !  Aggregates identifier
       integer,dimension(nnode),intent(out)       ::  itag    !  Aggregates size
       integer,dimension(nnode),intent(out)       ::  nmol    !  Number of aggregates of each size
       integer,intent(out)                        ::  nagg    !  Number of aggregates
       integer,intent(in)                         ::  nnode   !  Number of molecules
       integer,intent(out)                        ::  maxagg  !  Maximum aggregate size
!
! Local variables
!
       logical,dimension(nnode)                   ::  notvis  !  Nodes visited
       integer,dimension(nnode)                   ::  queue   !  Queue of connected nodes
       integer                                    ::  inode   !  Node index
       integer                                    ::  jnode   !  Node index
       integer                                    ::  knode   !  Node index
       integer                                    ::  iqueue  !  Queue index
       integer                                    ::  nqueue  !  Number of queue elements
       integer                                    ::  nnmol   !  Number of aggregates index
       integer                                    ::  ntag    !  Size of the aggregate
       integer                                    ::  nntag   !  Size of the aggregate index
       integer                                    ::  i       !  Indexes
!
! Performing Breadth First Search over the target molecules
! ---------------------------------------------------------
!
! Marking all the vertices as not visited
       notvis(:) = .TRUE.
! Initializing the molecules information
       imol(:) = 0   
       nmol(:) = 0
       nnmol   = 1
! Initializing the aggregates information
       nagg    = 0
       iagg(:) = 0
! Initializing the size information
       ntag    = 0
       nntag   = 0
       itag(:) = 0
! 
       maxagg  = 1
!
! Outer loop over each node
!
       do inode = 1, nnode
         if ( notvis(inode) ) then
! Marking head node as visited
           notvis(inode) = .FALSE.
! Updating the system information
           nagg        = nagg + 1
           iagg(nnmol) = nagg
!
           imol(nnmol) = inode
           nnmol       = nnmol + 1
!
           ntag        = 1
! Initializing queue
           queue(:) = 0
! Adding current node to the queue
           queue(1) = inode
! Initializing the queue counter
           iqueue = 1
! Setting the next position in the queue
           nqueue = 2
!
! Inner loop over the queue elements
!
           do while ( iqueue .lt. nqueue )
! Saving actual element in the queue
             knode = queue(iqueue)
! Checking the connection between actual queue element and the rest of nodes
             do jnode = inode + 1, nnode
! Checking if node j is connected to node k and has not been already visited
               if ( notvis(jnode) .and. adj(jnode,knode) ) then
! Updating the system information
                 iagg(nnmol)   = nagg
!
                 imol(nnmol)   = jnode
                 nnmol         = nnmol + 1
!
                 ntag          = ntag + 1
! Marking the node connected to node k as visited
                 notvis(jnode) = .FALSE.
! Adding to the queue the node connected to node k
                 queue(nqueue) = jnode
! Updating next position in the queue
                 nqueue        = nqueue + 1
               end if
             end do
! Updating the queue counter
             iqueue = iqueue + 1
           end do
! Saving the size of the aggregate found
           do i = nntag+1, nntag+ntag
             itag(i) = ntag
           end do
           nntag = nntag + ntag
! Update the number of aggregates of each size
           if ( ntag .gt. maxagg ) maxagg = ntag
           nmol(ntag)   = nmol(ntag)   + 1
         end if
       end do
!
       return
       end subroutine findcompundir
!
!======================================================================!
!
! BLOCKDIAG - BLOCK DIAGonalization
!
! This subroutine block diagonalizes an input adjacency matrix 
!  ADJ(NNODE,NNODE) of an undirected unweighted graph of NNODE vertices.
! The subroutine FINDCOMPUNDIR is employed to find the connected compo-
!  nents of the graph using BFS and the QUICKSHORT algorithm is used to
!  find the basis of nodes that block diagonalizes the adjacency matrix
!  sorting the blocks by size and identifier, and then the molecules
!  forming the aggregate according to their canonical order.
!
       subroutine blockdiag(nnode,adj,mol,tag,agg,nsize,nagg,iagg, &
                            nmol,imol,magg)
!
       use sorting,    only:  ivvqsort,ivqsort,iqsort
!
       implicit none
!
! Input/output variables
!
       logical,dimension(nnode,nnode),intent(in)   ::  adj       !  Adjacency matrix
       integer,dimension(nnode),intent(out)        ::  mol       !  Molecules identifier
       integer,dimension(nnode),intent(out)        ::  tag       !  Aggregates identifier
       integer,dimension(nnode),intent(out)        ::  agg       !  Aggregates size
       integer,dimension(nnode),intent(out)        ::  nagg      !  Number of aggregates of each size
       integer,dimension(nnode),intent(out)        ::  iagg      !
       integer,dimension(nnode),intent(out)        ::  nmol      !
       integer,dimension(nnode),intent(out)        ::  imol      !
       integer,intent(in)                          ::  nnode     !  Number of molecules
       integer,intent(out)                         ::  nsize     !  Maximum aggregate size
       integer,intent(out)                         ::  magg      !  Number of aggregates
!
! Local variables
! 
       integer                                     ::  i,j,k     !  Indexes 
!
! Block-diagonalizing the adjacency matrix 
!     
       call findcompundir(nnode,adj,mol,tag,agg,nsize,nagg,magg)
!
! Sorting the blocks of the adjacency matrix according to their size,
!  identifier, and canonical order of the constituent molecules
!
! Setting up iagg, imol and nmol arrays
!
       iagg(1) = 0
       imol(1) = 0
       do i = 1, nsize-1
         iagg(i+1) = iagg(i) + nagg(i)
         nmol(i)   = i*nagg(i) 
         imol(i+1) = imol(i) + nmol(i)
       end do
       nmol(nsize) = nnode - imol(nsize)
! 
! Sorting molecules and aggregate identifiers based on the size of 
!  the aggregates
!
       call ivvqsort(nnode,agg,tag,mol,1,nnode)
!
! Sorting molecules based on their aggregate identifier
!
       do i = 2, nsize-1
         if ( nagg(i) .gt. 1 ) then
           call ivqsort(nnode,tag,mol,imol(i)+1,imol(i+1))
         end if
       end do
!
       if ( nagg(nsize) .gt. 1 )                                       &
                         call ivqsort(nnode,tag,mol,imol(nsize)+1,nnode)
!
! Sorting molecules based on their canonical order
!
       do i = 2, nsize
         k = imol(i)
         do j = 1, nagg(i)
           call iqsort(nnode,mol,k+1,k+i)
           k = k + i
         end do
       end do
!
       return
       end subroutine blockdiag
!
!======================================================================!
!
! CALCDEGUNDIR - CALCulate DEGrees UNDIRected
!
! This function calculates the degree of the NNODE vertices in an undi-
!  rected unweighted graph given as an adjacency matrix representaion 
!  ADJ(NODE,NODE) and stores them in the array DEGREE(NNODE)
!
       function calcdegundir(nnode,adj) result(degree)
!
       implicit none
!
! Input/output variables
!
       logical,dimension(nnode,nnode),intent(in)  ::  adj     !
       integer,intent(in)                         ::  nnode   !
       integer,dimension(nnode)                   ::  degree  !
!
! Local variables
!
       integer                                    ::  inode   !
       integer                                    ::  jnode   !
!
! Finding the degree of the vertices in a graph
! ---------------------------------------------
!
       do inode = 1, nnode
         degree(inode) = 0
! Traversing through row/column of each vertex 
         do jnode = 1, nnode
! If a path from this vertex to other exists then increment the degree
           if ( adj(jnode,inode) ) degree(inode) = degree(inode) + 1
         end do
       end do
!    
       return
       end function calcdegundir
!
!======================================================================!
!
! DIJUNWUNDIR - DIJkstra UNWeighted-UNDIRrected graph
!
! This subroutine computes the minimum distance matrix in an unweighted-
!   undirected graph making use of Dijkstra's algorithm
!
       subroutine dijkunwundir(n,adj,mindis)
!
       implicit none
!
! Input/output variables
!
       integer,intent(in)                  ::  n       !  Number of nodes
       integer,intent(out),dimension(n,n)  ::  mindis  !  Minimum distance matrix
       logical,intent(in),dimension(n,n)   ::  adj     !  Adjacency matrix
!
! Local variables
!
       logical,dimension(n)                ::  visit   !  Visited nodes
       integer                             ::  dist    !  Distance to previous node
       integer                             ::  prev    !  Previous node
       integer                             ::  i,j,k   !  Indexes
       integer,parameter                   ::  INF = 2147483647
!
! Initializing variablies
!
       do k = 1, n
         do i = k+1, n
           if ( adj(i,k) ) then
             mindis(i,k) = 1
             mindis(k,i) = 1
           else
             mindis(i,k) = INF
             mindis(k,i) = INF
           end if
         end do
       end do
!
       do k = 1, n
         mindis(k,k) = 0
       end do
!
! Running Dijkstra algorithm starting from every node
!
       do k = 1, n
!
         visit(:) = .FALSE.
         visit(k) = .TRUE.
!
         do i = 1, n
!
           prev = -1
           dist = INF
!
           do j = 1, n
             if ( (.NOT.visit(j)) .AND. (mindis(j,k).lt.dist) ) then
               dist = mindis(j,k)
               prev = j
             end if
           end do
!
           if ( prev .eq. -1 ) exit  ! Exit if no unvisited node is found
!
           visit(prev) = .TRUE.
!
           do j = 1, n
             if ( .NOT. visit(j) ) then
               if ( adj(j,prev) .AND. ((mindis(prev,k)+1).lt.mindis(j,k)) ) then
                 mindis(j,k) = mindis(prev,k) + 1
!~                  mindis(k,j) = mindis(j,k)
               end if
             end if
           end do
!
         end do
!
       end do
!
!~        do k = 1, n
!~          do i = 1, n
!~            if ( .NOT. adj(i,k) ) then
!~              mindis(k,i) = mindis(i,k)
!~            end if
!~          end do
!~        end do
!
       return
       end subroutine dijkunwundir
!
!======================================================================!
!
!~   SUBROUTINE findcycles(n,  graph,cyclenumber, cycleSize, cycles) ! FIXME: doesnt work
!~     INTEGER, INTENT(IN) :: n
!~     INTEGER, INTENT(OUT) :: cyclenumber
!~     logical, DIMENSION(N, N), INTENT(IN) :: graph
!~     INTEGER, DIMENSION(N, N), INTENT(OUT) :: cycles
!~     INTEGER, DIMENSION(N), INTENT(OUT) :: cycleSize
!~ !
!~   INTEGER, DIMENSION(n) :: color, par

!~   ! Inicializar variables
!~   cyclenumber = 0
!~   color = 0
!~   par = 0
!~   cycles = 0
!~   cycleSize = 0  

!~   ! Llamar a DFS para marcar los ciclos
!~   CALL dfs_cycle(n,1, 0, color, par, cyclenumber, graph, cycles, cycleSize)

!~   END SUBROUTINE findcycles
!~ !
!~ !======================================================================!
!~ !
!~   RECURSIVE SUBROUTINE dfs_cycle(n,u, p, color, par, cyclenumber, graph, cycles, cycleSize)
!~     INTEGER, INTENT(IN) :: u, p,n
!~     INTEGER, DIMENSION(N), INTENT(INOUT) :: color, par
!~     INTEGER, INTENT(INOUT) :: cyclenumber
!~     logical, DIMENSION(N, N), INTENT(IN) :: graph
!~     INTEGER, DIMENSION(N, N), INTENT(INOUT) :: cycles
!~     INTEGER, DIMENSION(N), INTENT(INOUT) :: cycleSize
!~     INTEGER :: v, cur, cycleIndex

!~     IF (color(u) == 2) RETURN

!~     IF (color(u) == 1) THEN
!~       cyclenumber = cyclenumber + 1
!~       cycleIndex = 0
!~       cur = p
!~       DO
!~         cycleIndex = cycleIndex + 1
!~         cycles(cyclenumber, cycleIndex) = cur
!~         IF (cur == u) EXIT
!~         cur = par(cur)
!~       END DO
!~       cycleSize(cyclenumber) = cycleIndex
!~       RETURN
!~     END IF

!~     par(u) = p
!~     color(u) = 1

!~     DO v = 1, N
!~       IF (graph(u, v) .AND. v /= par(u)) THEN
!~         CALL dfs_cycle(n, v, u, color, par, cyclenumber, graph, cycles, cycleSize)
!~       END IF
!~     END DO

!~     color(u) = 2
!~   END SUBROUTINE dfs_cycle
!
!======================================================================!
!
! FINDCYCLE - FIND minimum CYCLE basis in undirected-unweighted graph
!
! This subroutine computes the minimum cycle basis in an unweighted-
!   undirected graph 
!
! Adapted from networkx
!
!
       subroutine findcycle(ne,n,r,adj,mcycle,ncycle,cycles)
!
       implicit none
!
! Input/output variables
!
       logical,dimension(n,n),intent(in)   ::  adj     !  Adjacency matrix
       integer,dimension(r,n),intent(out)  ::  cycles  !  Nodes belonging to cycles
       integer,dimension(r),intent(out)    ::  ncycle  !  Number of nodes per cycle
       integer,intent(out)                 ::  mcycle  !  Number of cycles
       integer,intent(in)                  ::  ne      !  Number of edges
       integer,intent(in)                  ::  n       !  Number of nodes
       integer,intent(in)                  ::  r       !  Rank of the graph
!
! Local variables
!
       logical,dimension(n,n)              ::  mstadj
       logical,dimension(n,n)              ::  chordsadj
       logical,dimension(n,n)              ::  set_orthadj
       logical,dimension(n,n)              ::  baseadj
       logical,dimension(n,n)              ::  minadj
       integer,dimension(2,ne)             ::  iedge   !  Adjacency list
       integer,dimension(2,n-1)            ::  mstedge  
       integer,dimension(2,ne)             ::  set_orthedge
       integer,dimension(2,ne)             ::  minedge
       integer,dimension(n)                ::  minpath
       integer                             ::  nset_orth
       integer                             ::  nbase
       integer                             ::  minlength
       integer                             ::  icount
       integer                             ::  u,v
       integer                             ::  i,j,k
!
! Initializing variables
!
       k = 0
       do i = 1, n-1
         do j = i+1, n
           if ( adj(i,j) ) then
             k = k + 1
             iedge(1,k) = i
             iedge(2,k) = j
           end if
         end do
       end do
!~ write(*,*) 'CYCLE RANK',r
!~ write(*,*)
!
       cycles(:,:) = 0
       ncycle(:)   = 0
       mcycle      = 0
!
! Printing adjacency matrix
!
!~ write(*,*) 'Adjacency matrix'
!~ write(*,*) '----------------'
!~ do i = 1, n
!~   write(*,*) (adj(i,j),j=1,n)
!~ end do
!~ write(*,*)
!
       call mst_kruskal(ne,iedge,n,mstedge,mstadj)
!
!~ write(*,*) 'Minimum spanning tree edges:'
!~ do i = 1, n-1
!~   write(*,*) mstedge(1, i), '-', mstedge(2, i)
!~ end do
!~ write(*,*)
!
       chordsadj(:,:) = adj(:,:)
       do i = 1, n-1
         chordsadj(mstedge(1,i),mstedge(2,i)) = .FALSE.
         chordsadj(mstedge(2,i),mstedge(1,i)) = .FALSE.
       end do
!
!~ write(*,*) 'Chords matrix'
!~ write(*,*) '-------------'
!~ do i = 1, n
!~   write(*,*) (chordsadj(i,j),j=1,n)
!~ end do
!~ write(*,*)
!
       nset_orth = 0
       do i = 1, n-1
         do j = i+1, n
           if ( chordsadj(i,j) ) then
             nset_orth = nset_orth + 1
             set_orthedge(1,nset_orth) = i
             set_orthedge(2,nset_orth) = j
           end if
         end do
       end do
!
       set_orthadj(:,:) = chordsadj(:,:)
!
!~ write(*,*) 'nset_orth    ',nset_orth
!~ write(*,*) 'Initial set_orth edges'
!~ write(*,*) '----------------------'
!~ do i = 1, nset_orth
!~   write(*,*) set_orthedge(1, i), '-', set_orthedge(2, i)
!~ end do
!~ write(*,*)
!~ write(*,*) 'Initial set_orth matrix'
!~ write(*,*) '-----------------------'
!~ do i = 1, n
!~   write(*,*) (set_orthadj(i,j),j=1,n)
!~ end do
!~ write(*,*)
!
       baseadj(:,:) = .FALSE.
       nbase = 0
!
! Main loop
!
       do while ( nset_orth .gt. 0 )
!
!~ write(*,*) 'STARTING NEW ITERATION'
!~ write(*,*) '----------------------'
!
         u = set_orthedge(1,nset_orth)
         v = set_orthedge(2,nset_orth)
!
         nbase = nbase + 1
!
         baseadj(u,v) = .TRUE.
         baseadj(v,u) = .TRUE.
!
         set_orthadj(u,v) = .FALSE.
         set_orthadj(v,u) = .FALSE.
!
         nset_orth = nset_orth - 1
!
!~ write(*,*) 'Updated base matrix'
!~ write(*,*) '-------------------'
!~ do i = 1, n
!~   write(*,*) (baseadj(i,j),j=1,n)
!~ end do
!~ write(*,*)
!
         call min_cycle(ne,iedge,n,baseadj,minpath,minedge,minadj,     &
                        minlength) 
!
!~ write(*,*) 'cycle ',minpath
!~ write(*,*) 'length',minlength 
!
         mcycle = mcycle + 1
         cycles(mcycle,:) = minpath(:)
         ncycle(mcycle)   = minlength
!
!~ write(*,*)
!~ write(*,*) 'cycles size'
!~ write(*,*) '-----------'
!~ write(*,*) ncycle(:)
!~ write(*,*) 'Updated cycles matrix'
!~ write(*,*) '---------------------'
!~ do i = 1, n
!~   write(*,*) (cycles(i,j),j=1,r)
!~ end do
!~ write(*,*)
!
! Counting the number of commom edges between set_orth and minedge
!
         icount = 0
         do i = 1, minlength
           if ( set_orthadj(minedge(1,i),minedge(2,i)) ) then
             icount = icount + 1
           end if
         end do
!
! If the number of common edges is odd, update set_orth
!
        if ( mod(icount,2) .eq. 0 ) then
          do i = 1, n-1
            do j = i+1, n
              if ( baseadj(i,j) .AND. set_orthadj(i,j) ) then
! If the edge is in base and set_orth, remove it from set_orth
                set_orthadj(i,j) = .FALSE.
                set_orthadj(j,i) = .FALSE.
              else if ( baseadj(i,j) .AND. set_orthadj(i,j) ) then
! If the edge is in base but not in set_orth, add it to set_orth
                set_orthadj(i,j) = .TRUE.
                set_orthadj(j,i) = .TRUE.
              end if
! If the edge is in set_orth but not in base, keep it in set_orth
            end do
          end do
        end if
!
         nset_orth = 0
         do i = 1, n-1
           do j = i+1, n
             if ( set_orthadj(i,j) ) then
               nset_orth = nset_orth + 1
               set_orthedge(1,nset_orth) = i
               set_orthedge(2,nset_orth) = j
             end if
           end do
         end do
!
!write(*,*) 'nset_orth',nset_orth
!write(*,*) 'Updated set_orth edges'
!write(*,*) '----------------------'
!do i = 1, nset_orth
!  write(*,*) set_orthedge(1, i), '-', set_orthedge(2, i)
!end do
!write(*,*)
!write(*,*) 'Updated set_orth matrix'
!write(*,*) '-----------------------'
!do i = 1, n
!  write(*,*) (set_orthadj(i,j),j=1,n)
!end do
!write(*,*)
!
       end do
!
       return
       end subroutine findcycle
!
!======================================================================!
!
       subroutine min_cycle(ne,iedge,n,orthadj,minpath,minedge,        &
                            minadj,minlength)
!
! Input/output variables
!
       logical,dimension(n,n),intent(in)    ::  orthadj    !
       logical,dimension(n,n),intent(out)   ::  minadj     !
       integer,dimension(2,ne),intent(in)   ::  iedge      !
       integer,dimension(2,ne),intent(out)  ::  minedge    !
       integer,dimension(n),intent(out)     ::  minpath    !
       integer,intent(out)                  ::  minlength  !
       integer,intent(in)                   ::  ne         !
       integer,intent(in)                   ::  n          !
!
! Local variables
!
       logical,dimension(n*2,n*2)           ::  auxadj     !
       integer,dimension(n*2,n*2)           ::  mindis     !
       integer,dimension(n)                 ::  lift       !
       integer,dimension(n*2)               ::  path       !
       integer                              ::  length     ! 
       integer                              ::  istart     ! 
       integer                              ::  iend       !
       integer                              ::  u,v        !
       integer                              ::  i,j        !
!
! Initializing variables
! ----------------------
!
       auxadj(:,:) = .FALSE.
!
! Adding 2 copies of each edge in G to Gi
!
       do i = 1, ne
!   
         u = iedge(1,i)
         v = iedge(2,i)
!
! If edge is in orth then add cross edge
!
         if ( orthadj(u,v) ) then
           auxadj(u,n+v) = .TRUE.
           auxadj(n+v,u) = .TRUE.
!      
           auxadj(v,n+u) = .TRUE.
           auxadj(n+u,v) = .TRUE.
         else
!
! Otherwise add in-plane edge
!
           auxadj(u,v) = .TRUE.
           auxadj(v,u) = .TRUE.
!
           auxadj(n+u,n+v) = .TRUE.
           auxadj(n+v,n+u) = .TRUE.
         end if
!
       end do
!
!~ write(*,*) 'Auxiliary matrix'
!~ write(*,*) '----------------'
!~ do i = 1, n*2
!~   write(*,*) (auxadj(i,j),j=1,n*2)
!~ end do
!~ write(*,*)
!
! Computing minimum distance matrix of Gi
!
       call dijkunwundir(n*2,auxadj,mindis)
!
!~ write(*,*) 'Minimum distance matrix'
!~ write(*,*) '-----------------------'
!~ do i = 1, n*2
!~   write(*,*) (mindis(i,j),j=1,n*2)
!~ end do
!~ write(*,*)
!
! Finding shortest length in Gi between n and its copy
!
       lift(:) = 0
       do j = 1, n
         lift(j) = mindis(j,n+j)
       end do
!
!~ write(*,*) 'lift',lift(:)
!~ write(*,*)
!
! Compute the shortest path in Gi between n and its copy,
!  which translates to a cycle in G
!
       istart = minloc(lift,DIM=1)
       iend   = n + istart
!
!~ write(*,*) 'istart',istart
!~ write(*,*) 'iend  ',iend
!~ write(*,*)
!
       call dijkunwundirpath(istart,iend,n*2,auxadj,path,length) 
!
!~ write(*,*) 'tmppath  ',path
!~ write(*,*) 'tmplength',length
!~ write(*,*)
!
! Re-mapping nodes in Gi to those in G
!
       do j = 1, length
         if ( path(j) .gt. n ) path(j) = path(j) - n
       end do
!
!~ write(*,*) 'tmppath after re-mapping',path
!~ write(*,*) 
!
! Removing the edges that occur two times
!
       minpath(:) = 0
!    
       minedge(:,:) = 0
       minadj(:,:) = .FALSE.
!
       minlength = 0
       do j = 1, length-1
         if ( .NOT. minadj(path(j),path(j+1)) ) then
!
           minadj(path(j),path(j+1)) = .TRUE.
           minadj(path(j+1),path(j)) = .TRUE.
!
           minlength = minlength + 1
!
           minedge(1,minlength) = path(j)
           minedge(2,minlength) = path(j+1)
!
           minpath(minlength) = path(j)
!
         end if
       end do
!
!~ write(*,*) 'minlength',minlength
!~ write(*,*) 'minpath:',minpath(:)
!~ write(*,*) 'minedge:'
!~ do i = 1, minlength
!~   write(*,*) minedge(1, i), '-', minedge(2, i)
!~ end do
!~ write(*,*)
!~ write(*,*) 'Minimum path adjacency matrix'
!~ write(*,*) '-----------------------------'
!~ do i = 1, n
!~   write(*,*) (minadj(i,j),j=1,n)
!~ end do
!~ write(*,*)
!
       return
       end subroutine min_cycle
!
!======================================================================!        
!    
       subroutine dijkunwundirpath(istart,iend,n,adj,path,length)
!
       implicit none
!
! Input/output variables
!
       logical,intent(in),dimension(n,n)  ::  adj     !  Adjacency matrix
       integer,intent(out),dimension(n)   ::  path    !  Minimum distance path 
       integer,intent(out)                ::  length  !  Minimum distnace
       integer,intent(in)                 ::  n       !  Number of nodes
       integer,intent(in)                 ::  istart  ! 
       integer,intent(in)                 ::  iend    ! 
!
! Local variables
!
       logical,dimension(n)               ::  visit   !  Visited nodes
       integer,dimension(n)               ::  mindis  !  Minimum distances
       integer,dimension(n)               ::  iprev   !  Previous node
       integer                            ::  dist    !  Distance to previous node
       integer                            ::  prev    !  Actual previous node
       integer                            ::  i,j     !  Indexes
       integer,parameter                  ::  INF = 2147483647
!
! Initializing variablies
!
       do i = 1, n
         if ( adj(istart,i) ) then
           mindis(i) = 1
         else
           mindis(i) = INF
         end if
       end do
       mindis(istart) = 1
!
       visit(:)      = .FALSE.
       visit(istart) = .TRUE.
!
       iprev(:) = -1
!
! Running Dijkstra algorithm starting from every node
!
       do i = 1, n
!
         prev = -1
         dist = INF
!
         do j = 1, n
           if ( (.NOT.visit(j)) .AND. (mindis(j).lt.dist) ) then
             dist = mindis(j)
             prev = j
           end if
         end do
!
         if ( prev .eq. -1 ) exit  ! Exit if no unvisited node is found
!
         visit(prev) = .TRUE.
!
         do j = 1, n
           if ( .NOT. visit(j) ) then
             if ( adj(j,prev) .AND. ((mindis(prev)+1).lt.mindis(j)) ) then
               mindis(j) = mindis(prev) + 1
               iprev(j) = prev
             end if
           end if
         end do
!
       end do
!
       path(:)  = 0
       length   = 0
!
       i = iend
       do while ( i .ne. -1 )
         path(length+1) = i
         length = length + 1
         i = iprev(i)
       end do
!
       if ( path(length) .ne. istart ) then
         length = length + 1
         path(length) = istart
       end if
!
       return
       end subroutine dijkunwundirpath
!
!======================================================================!
!
       subroutine mst_kruskal(ne,iedge,n,mst,mstadj)
!
       implicit none
!
! Input/output variables
!
       integer,dimension(2,n-1),intent(out)  ::  mst     !
       integer,dimension(2,ne),intent(out)   ::  iedge   !
       logical,dimension(n,n),intent(out)    ::  mstadj  !
       integer,intent(in)                    ::  ne
       integer,intent(in)                    ::  n
!
! Local variables
!
       integer,dimension(n)                  ::  parent  !
       integer                               ::  u,v     !
       integer                               ::  uroot   !
       integer                               ::  vroot   !
       integer                               ::  i,j     !
! 
!
!
       mstadj(:,:) = .FALSE.
!
       do i = 1, n
         parent(i) = i
       end do 

       j = 0
       do i = 1, ne
!
         u = iedge(1,i)
         v = iedge(2,i)
!
         uroot = find(u,n,parent)
         vroot = find(v,n,parent)
!
         if ( uroot .ne. vroot ) then
           call union(uroot,vroot,n,parent)
           j = j + 1
           mst(1,j) = u
           mst(2,j) = v
           mstadj(u,v) = .TRUE.
           mstadj(v,u) = .TRUE.
           if ( j .eq. (n-1) ) return
         end if
!
       end do
!
       return
       end subroutine mst_kruskal
!
!======================================================================!
!
       recursive function find(u,n,parent) result(root)
!
       implicit none
!
! Input/output variables
!
       integer,dimension(n),intent(inout)  ::  parent  !
       integer,intent(in)                  ::  n       !
       integer,intent(in)                  ::  u       !
       integer                             ::  root    !
!
!
!
       if ( parent(u) .ne. u ) then
         parent(u) = find(parent(u),n,parent)
       end if
       root = parent(u)
!
       return
       end function
!
!======================================================================!
!
       subroutine union(u,v,n,parent)
!
       implicit none
!
! Input/output variables
!
       integer,dimension(n),intent(inout)  ::  parent  !
       integer,intent(in)                  ::  u,v     !
       integer,intent(in)                  ::  n       !
!
! Local variables
!
       integer                             ::  uroot
       integer                             ::  vroot
!
!
!
       uroot = find(u,n,parent)
       vroot = find(v,n,parent)
       parent(uroot) = vroot
!
       return
       end subroutine
!
!======================================================================!
!
       end module graphtools
!
!======================================================================!

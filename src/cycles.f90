!======================================================================!
!
       module cycles
!
       implicit none
!
       contains
!
!======================================================================!
!
  SUBROUTINE find_minimum_cycle_basis()
  IMPLICIT NONE
  integer,dimension(2,num_nodes-1)       ::  mstedge
  integer,dimension(2,num_edges)         ::  set_orthedge
  integer,dimension(2,num_edges)         ::  minedge
  integer,dimension(num_nodes)           ::  minpath
  logical,dimension(num_nodes,num_nodes) ::  mstadj
  logical,dimension(num_nodes,num_nodes) ::  chordsadj
  logical,dimension(num_nodes,num_nodes) ::  set_orthadj
  logical,dimension(num_nodes,num_nodes) ::  minadj
  logical,dimension(num_nodes,num_nodes) ::  baseadj
  
  integer,dimension(:,:),allocatable     ::  cycles
  integer,dimension(:),allocatable       ::  ncycle
  integer                                ::  rank
  integer                                ::  mcycle
  
  integer                                ::  nset_orth
  integer                                ::  nbase
  integer                                ::  minlength
  integer                                ::  icount
  integer                                ::  u,v
  integer                                ::  i,j,k
!
  rank = num_edges - num_nodes + 1
!
write(*,*) 'CYCLE RANK',rank
write(*,*)
!
  allocate(cycles(num_nodes,rank),ncycle(rank))
!
  cycles(:,:) = 0
  ncycle(:)   = 0
  mcycle      = 0
!
! Printing adjacency matrix
!
write(*,*) 'Adjacency matrix'
write(*,*) '----------------'
do i = 1, num_nodes
  write(*,*) (adj(i,j),j=1,num_nodes)
end do
write(*,*)
!
  call mst_kruskal(mstedge,mstadj)
!
! Imprimir el resultado del Árbol de Expansión Mínima
!
  write(*,*) 'El Árbol de Expansión Mínima incluye las siguientes aristas:'
  do i = 1, num_nodes-1
      write(*,*) mstedge(1, i), '-', mstedge(2, i)
  end do
  write(*,*)
!
  chordsadj(:,:) = adj(:,:)
  do i = 1, num_nodes-1
    chordsadj(mstedge(1,i),mstedge(2,i)) = .FALSE.
    chordsadj(mstedge(2,i),mstedge(1,i)) = .FALSE.
  end do
!
write(*,*) 'Chords matrix'
write(*,*) '-------------'
do i = 1, num_nodes
  write(*,*) (chordsadj(i,j),j=1,num_nodes)
end do
write(*,*)
!
  nset_orth = 0
  do i = 1, num_nodes-1
    do j = i+1, num_nodes
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
write(*,*) 'nset_orth    ',nset_orth
write(*,*) 'Initial set_orth edges'
write(*,*) '----------------------'
do i = 1, nset_orth
  write(*,*) set_orthedge(1, i), '-', set_orthedge(2, i)
end do
write(*,*)
write(*,*) 'Initial set_orth matrix'
write(*,*) '-----------------------'
do i = 1, num_nodes
  write(*,*) (set_orthadj(i,j),j=1,num_nodes)
end do
write(*,*)
!
  baseadj(:,:) = .FALSE.
  nbase = 0
!
! Main loop
!
  do while ( nset_orth .gt. 0 )
!
write(*,*) 'STARTING NEW ITERATION'
write(*,*) '----------------------'
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
write(*,*) 'Updated base matrix'
write(*,*) '-------------------'
do i = 1, num_nodes
  write(*,*) (baseadj(i,j),j=1,num_nodes)
end do
write(*,*)
!
    call min_cycle(baseadj,minpath,minedge,minadj,minlength) 
!
write(*,*) 'cycle ',minpath
write(*,*) 'length',minlength 
!
    mcycle = mcycle + 1
    cycles(:,mcycle) = minpath(:)
    ncycle(mcycle)   = minlength
!
write(*,*)
write(*,*) 'cycles size'
write(*,*) '-----------'
write(*,*) ncycle(:)
write(*,*) 'Updated cycles matrix'
write(*,*) '---------------------'
do i = 1, num_nodes
  write(*,*) (cycles(i,j),j=1,rank)
end do
write(*,*)
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
     do i = 1, num_nodes-1
       do j = i+1, num_nodes
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
    do i = 1, num_nodes-1
      do j = i+1, num_nodes
        if ( set_orthadj(i,j) ) then
          nset_orth = nset_orth + 1
          set_orthedge(1,nset_orth) = i
          set_orthedge(2,nset_orth) = j
        end if
      end do
    end do
!
write(*,*) 'nset_orth',nset_orth
write(*,*) 'Updated set_orth edges'
write(*,*) '----------------------'
do i = 1, nset_orth
  write(*,*) set_orthedge(1, i), '-', set_orthedge(2, i)
end do
write(*,*)
write(*,*) 'Updated set_orth matrix'
write(*,*) '-----------------------'
do i = 1, num_nodes
  write(*,*) (set_orthadj(i,j),j=1,num_nodes)
end do
write(*,*)

!
  end do
!
  END SUBROUTINE find_minimum_cycle_basis
!
!======================================================================!
!
  SUBROUTINE min_cycle(orthadj,minpath,minedge,minadj,minlength)

  logical,dimension(num_nodes,num_nodes),intent(in)     ::  orthadj
  integer,dimension(num_nodes),intent(out)              ::  minpath
  logical,dimension(num_nodes,num_nodes),intent(out)    ::  minadj
  integer,dimension(2,num_edges),intent(out)            ::  minedge
  integer,intent(out)                                   ::  minlength

  logical,dimension(num_nodes*2,num_nodes*2)            ::  auxadj
  integer,dimension(num_nodes)                          ::  lift
  integer,dimension(num_nodes*2)                        ::  path
  integer                                               ::  length
  integer                                               ::  istart
  integer                                               ::  iend
  integer                                               ::  n
  integer                                               ::  nedge
  integer                                               ::  u,v
  integer                                               ::  i,j
!
  auxadj(:,:) = .FALSE.
!
  n = num_nodes
!
! Adding 2 copies of each edge in G to Gi
!
  do i = 1, num_edges
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
write(*,*) 'Auxiliary matrix'
write(*,*) '----------------'
do i = 1, num_nodes*2
  write(*,*) (auxadj(i,j),j=1,num_nodes*2)
end do
write(*,*)
!
! Computing minimum distance matrix of Gi
!
  CALL dijkunwundir(num_nodes*2,auxadj,mindis)
!
write(*,*) 'Minimum distance matrix'
write(*,*) '-----------------------'
do i = 1, num_nodes*2
  write(*,*) (mindis(i,j),j=1,num_nodes*2)
end do
write(*,*)
!
! Finding shortest length in Gi between n and its copy
!
  lift(:) = 0
  do j = 1, num_nodes
    lift(j) = mindis(j,n+j)
  end do
!
write(*,*) 'lift',lift(:)
write(*,*)
!
! Compute the shortest path in Gi between n and its copy,
!  which translates to a cycle in G
!
  istart = minloc(lift,DIM=1)
  iend   = n + istart
!
write(*,*) 'istart',istart
write(*,*) 'iend  ',iend
write(*,*)
!
  call dijkunwundirpath(istart,iend,n*2,auxadj,path,length) 
!
write(*,*) 'tmppath  ',path
write(*,*) 'tmplength',length
write(*,*)
!
! Re-mapping nodes in Gi to those in G
!
  do j = 1, length
    if ( path(j) .gt. n ) path(j) = path(j) - n
  end do
write(*,*) 'tmppath after re-mapping',path
write(*,*) 
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
write(*,*) 'minlength',minlength
write(*,*) 'minpath:',minpath(:)
write(*,*) 'minedge:'
do i = 1, minlength
  write(*,*) minedge(1, i), '-', minedge(2, i)
end do
write(*,*)
write(*,*) 'Minimum path adjacency matrix'
write(*,*) '-----------------------------'
do i = 1, num_nodes
  write(*,*) (minadj(i,j),j=1,num_nodes)
end do
write(*,*)
!
  END SUBROUTINE min_cycle
!
!======================================================================!        
!    
       subroutine dijkunwundirpath(istart,iend,n,adj,path,length)
!
       implicit none
!
! Input/output variables
!
       logical,intent(in),dimension(n,n)   ::  adj     !  Adjacency matrix
       integer,intent(out),dimension(n)    ::  path    !  Minimum distance path 
       integer,intent(out)                 ::  length  !  Minimum distnace
       integer,intent(in)                  ::  n       !  Number of nodes
       integer,intent(in)                  ::  istart  ! 
       integer,intent(in)                  ::  iend    ! 
!
! Local variables
!
       logical,dimension(n)                ::  visit   !  Visited nodes
       integer,dimension(n)                ::  mindis  !  Minimum distances
       integer,dimension(n)                ::  iprev   !  Previous node
       integer                             ::  dist    !  Distance to previous node
       integer                             ::  prev    !  Actual previous node
       integer                             ::  iaux    !  
       integer                             ::  i,j,k   !  Indexes
       integer,parameter                   ::  INF = 2147483647
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
  SUBROUTINE mst_kruskal(mst,mstadj)
  IMPLICIT NONE
  integer,dimension(2,num_nodes-1),intent(out)        ::  mst     !
  logical,dimension(num_nodes,num_nodes),intent(out)  ::  mstadj     !

  integer,dimension(num_nodes)                ::  parent  !
  integer                                     ::  u,v     !
  integer                                     ::  uroot   !
  integer                                     ::  vroot   !
  integer                                     ::  i,j     !
! 
  mstadj(:,:) = .FALSE.
!
  do i = 1, num_nodes
    parent(i) = i
  end do 

  j = 0
  do i = 1, num_edges
!
    u = iedge(1,i)
    v = iedge(2,i)
!
    uroot = find(u,parent)
    vroot = find(v,parent)
!
    if ( uroot .ne. vroot ) then
      call union(uroot,vroot,parent)
      j = j + 1
      mst(1,j) = u
      mst(2,j) = v
      mstadj(u,v) = .TRUE.
      mstadj(v,u) = .TRUE.
      if ( j .eq. (num_nodes-1) ) return
    end if
!
  end do
!
  END SUBROUTINE mst_kruskal
!
!======================================================================!
!
  RECURSIVE FUNCTION find(u,parent) RESULT(root)
  IMPLICIT NONE
  integer,dimension(num_nodes),intent(inout)  ::  parent  !
  integer,intent(in)                          ::  u       !
  integer                                     ::  root    !
!
  if ( parent(u) .ne. u ) then
    parent(u) = find(parent(u),parent)
  end if
  root = parent(u)
!
  END FUNCTION
!
!======================================================================!
!
  SUBROUTINE union(u,v,parent)
  IMPLICIT NONE
  integer,dimension(num_nodes),intent(inout)  ::  parent  !
  integer,intent(in)                          ::  u,v     !

  integer                                     ::  uroot
  integer                                     ::  vroot
!
  uroot = find(u,parent)
  vroot = find(v,parent)
  parent(uroot) = vroot
!
  END SUBROUTINE
!
!======================================================================!
!
       end module cycles
!
!======================================================================!

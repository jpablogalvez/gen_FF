!======================================================================!
!
       module isomorphism
!
       implicit none
!
       contains
!
!======================================================================!
!
! Función para comprobar el isomorfismo entre dos grafos
! ------------------------------------------------------
!
!~        logical function VF2_Isomorphism(N1,N2,AdjMatrix1,AdjMatrix2,   &
!~                                         NodeDegrees1,NodeDegrees2)
!~ !
!~        implicit none
!~ !
!~        integer,intent(in) :: N1, N2
!~        integer,intent(in) :: NodeDegrees1(N1)
!~        integer,intent(in) ::  NodeDegrees2(N2)
!~        logical,intent(in) :: AdjMatrix1(N1,N1)
!~        logical,intent(in) :: AdjMatrix2(N2,N2)
!~ !
!~ ! Variables locales
!~ !
!~        integer :: t1(N1) ! Añadido para almacenar el estado de mapeo
!~        integer :: t2(N2) ! Añadido para almacenar el estado de mapeo
!~        logical :: In1(N1)
!~        logical :: In2(N2)
!~        logical :: InMatch(N1)
!~        logical :: SolutionFound
!~ !
!~ ! Inicializar variables
!~ !
!~        t1 = 0
!~        t2 = 0
!~        In1 = .FALSE.
!~        In2 = .FALSE.
!~        InMatch = .FALSE.
!~        SolutionFound = .FALSE.
!~ !
!~ ! Llamar a la función recursiva para buscar el isomorfismo
!~ !
!~        call RecursiveSearch(N1,N2,AdjMatrix1,AdjMatrix2,NodeDegrees1,  &
!~                             NodeDegrees2,In1,In2,InMatch,              &
!~                             SolutionFound,t1,t2)
!~ !
!~ ! Devolver el resultado del isomorfismo
!~ !
!~        VF2_Isomorphism = SolutionFound
!~ !
!~        return
!~        end function VF2_Isomorphism
!
!======================================================================!
!
! Función para comprobar el isomorfismo entre dos grafos
! ------------------------------------------------------
!
       logical function VF2_Equivalent(N1,N2,AdjMatrix1,AdjMatrix2,    &
                                       NodeDegrees1,NodeDegrees2,      &
                                       order,v,w)
!
       implicit none
!
       integer,intent(in) :: N1, N2
       integer,intent(in) :: v,w
       integer,intent(in) :: Order(N1)
       integer,intent(in) :: NodeDegrees1(N1)
       integer,intent(in) ::  NodeDegrees2(N2)
       logical,intent(in) :: AdjMatrix1(N1,N1)
       logical,intent(in) :: AdjMatrix2(N2,N2)
!
! Variables locales
!
       integer :: t1(N1) ! Añadido para almacenar el estado de mapeo
       integer :: t2(N2) ! Añadido para almacenar el estado de mapeo
       logical :: In1(N1)
       logical :: In2(N2)
       logical :: InMatch(N1)
       logical :: SolutionFound
!
! Inicializar variables
!
       t1 = 0
       t2 = 0
       In1 = .FALSE.
       In2 = .FALSE.
       InMatch = .FALSE.
       SolutionFound = .FALSE.
!
       InMatch(v) = .TRUE.
       In1(v) = .TRUE.
       In2(w) = .TRUE.
       t1(v) = w
       t2(w) = v
!
! Llamar a la función recursiva para buscar el isomorfismo
!
       call RecursiveSearch(N1,N2,AdjMatrix1,AdjMatrix2,NodeDegrees1,  &
                            NodeDegrees2,In1,In2,InMatch,              &
                            SolutionFound,t1,t2,Order)
!
! Devolver el resultado del isomorfismo
!
       VF2_Equivalent = SolutionFound
!~ !
!~        if ( SolutionFound ) then
!~          do i = 1, N1
!~ !
!~            eqv(t1(i),i) = .TRUE.
!~            eqv(i,t1(i)) = .TRUE.
!~ !
!~            if ( t1(i) .gt. i ) check(t1(i)) = .FALSE.
!~ !
!~          end do
!~        end if      
!
       return
       end function VF2_Equivalent
!
!======================================================================!
!
! Subrutina para la búsqueda recursiva de isomorfismos
! ---------------------------------------------------- 
! 
       recursive subroutine RecursiveSearch(N1,N2,AdjMatrix1,          &
                                            AdjMatrix2,NodeDegrees1,   &
                                            NodeDegrees2,In1,In2,      &
                                            InMatch,SolutionFound,     &
                                            t1,t2,Order)
!
       implicit none
!
       integer,intent(in) :: N1, N2
       integer,intent(in) :: Order(N1)
       integer,intent(in) :: NodeDegrees1(N1)
       integer,intent(in) :: NodeDegrees2(N2)
       logical,intent(in) :: AdjMatrix1(N1,N1)
       logical,intent(in) :: AdjMatrix2(N2,N2)
       logical,intent(inout) :: In1(N1)
       logical,intent(inout) :: In2(N2)
       logical,intent(inout) :: InMatch(N1)
       logical,intent(inout) :: SolutionFound
       integer,intent(inout) :: t1(N1) ! Añadido para almacenar el estado de mapeo
       integer,intent(inout) :: t2(N2) ! Añadido para almacenar el estado de mapeo
!
! Variables locales
!
!~        logical :: CheckRules
       integer :: v, w, i
!
! Si se ha encontrado una solución, devolver
!
       if ( SolutionFound ) return
!
! Si todos los nodos están en el mapeo, se ha encontrado una solución
!
       if ( all(InMatch) ) then
         SolutionFound = .TRUE.
         return
       end if
!
! Seleccionar el siguiente candidato
!
       v = 0
       do i = 1, N1
         if ( .NOT. InMatch(order(i)) ) then
           v = order(i)
           exit
         end if
       end do
!~ write(*,*) 'next node for the mapping',v
!
! Probar todos los nodos posibles en el grafo 2
!
       do w = 1, N2
         if ( .not.In2(w) ) then
!
! Comprobar las reglas de corte y factibilidad
!
           if ( CheckRules(v,w,N1,N2,AdjMatrix1,AdjMatrix2,            &
                           NodeDegrees1,NodeDegrees2,t1,t2) ) then
!
! Añadir los nodos al mapeo
!
             InMatch(v) = .TRUE.
             In1(v) = .TRUE.
             In2(w) = .TRUE.
             t1(v) = w
             t2(w) = v
!
! Llamar a la función recursiva para el siguiente nivel de la búsqueda
!
             call RecursiveSearch(N1,N2,AdjMatrix1,AdjMatrix2,         &
                                  NodeDegrees1,NodeDegrees2,In1,In2,   &
                                  InMatch,SolutionFound,t1,t2,Order)
!
! Si se ha encontrado una solución, devolver
!
             if ( SolutionFound ) return
!
! Eliminar los nodos del mapeo
!
             InMatch(v) = .FALSE.
             In1(v) = .FALSE.
             In2(w) = .FALSE.
             t1(v) = 0
             t2(w) = 0
!
           end if
         end if
       end do
!
       return
       end subroutine RecursiveSearch
!
!======================================================================!
!
! Función para verificar las reglas de corte y factibilidad
! ---------------------------------------------------------
!
       logical function CheckRules(v,w,N1,N2,AdjMatrix1,AdjMatrix2,    &
                                   NodeDegrees1,NodeDegrees2,t1,t2)
!
       implicit none
!
       integer, intent(in) :: v,w
       integer, intent(in) :: N1,N2
       logical, intent(in) :: AdjMatrix1(N1,N1)
       logical, intent(in) :: AdjMatrix2(N2,N2)
       integer, intent(in) :: t1(N1)
       integer, intent(in) :: t2(N2)
       integer,intent(in) :: NodeDegrees1(N1)
       integer,intent(in) :: NodeDegrees2(N2)

!
! Variables locales
!
       integer :: i
!~ !
!~ !
!~ !
       CheckRules = .FALSE.
!
! Regla 1: Los nodos v y w deben tener el mismo grado
!
       if ( .NOT. (NodeDegrees1(v) == NodeDegrees2(w)) ) return
!
! Regla 2: Para cada nodo adyacente a v en el grafo 1, debe existir un nodo adyacente a w en el grafo 2
!
       do i = 1, N1
         if ( AdjMatrix1(v,i) .and. t1(i) > 0 ) then
!~ !write(*,*) 'BACKTRACE'
!~ !write(*,*) t1(i),t2(t1(i))
           if ( .NOT. AdjMatrix2(w,t1(i)) ) return
         end if
       end do
!
! Regla 3: Para cada nodo adyacente a w en el grafo 2, debe existir un nodo adyacente a v en el grafo 1
!
       do i = 1, N2
         if ( AdjMatrix2(w,i) .and. t2(i) > 0 ) then
           if ( .NOT. AdjMatrix1(v,t2(i)) ) return
         end if
       end do
!
! Regla 4: Los nodos v y w deben tener el mismo número de nodos adyacentes no mapeados
!
       if ( .NOT. (count(AdjMatrix1(v,:) .and.                         &
                  (t1 == 0)) == count(AdjMatrix2(w,:) .and.            &
                                                    (t2 == 0))) ) return
!
! Regla 5: Para cada nodo no mapeado adyacente a v en el grafo 1, debe existir un nodo no mapeado adyacente a w en el grafo 2
!
       do i = 1, N1
         if (AdjMatrix1(v,i) .and. t1(i) == 0) then
           if (.not. any(AdjMatrix2(w,:) .and. (t2 == 0))) return
         end if
       end do
!
! Regla 6: Para cada nodo no mapeado adyacente a w en el grafo 2, debe existir un nodo no mapeado adyacente a v en el grafo 1
!
       do i = 1, N2
         if (AdjMatrix2(w,i) .and. t2(i) == 0) then
           if (.not. any(AdjMatrix1(v,:) .and. (t1 == 0))) return
         end if
       end do
!
! Devolver el resultado de las reglas de corte y factibilidad
!
       CheckRules = .TRUE.
!
       return
       end function CheckRules
!
!======================================================================!
!
! Función para comprobar el isomorfismo entre dos grafos
! ------------------------------------------------------
!
       logical function VF2_EquivalentEVC(N1,N2,AdjMatrix1,AdjMatrix2, &
                                          NodeDegrees1,NodeDegrees2,   &
                                          NodeEVC1,NodeEVC2,Order,v,w)
!
       implicit none
!
       integer,intent(in) :: N1, N2
       integer,intent(in) :: v,w
       integer,intent(in) :: Order(N1)
       integer,intent(in) :: NodeDegrees1(N1)
       integer,intent(in) ::  NodeDegrees2(N2)
       real(kind=8),intent(in) :: NodeEVC1(N1)
       real(kind=8),intent(in) :: NodeEVC2(N2)
       logical,intent(in) :: AdjMatrix1(N1,N1)
       logical,intent(in) :: AdjMatrix2(N2,N2)
!
! Variables locales
!
       integer :: t1(N1) ! Añadido para almacenar el estado de mapeo
       integer :: t2(N2) ! Añadido para almacenar el estado de mapeo
       logical :: In1(N1)
       logical :: In2(N2)
       logical :: InMatch(N1)
       logical :: SolutionFound
!
! Inicializar variables
!
       t1 = 0
       t2 = 0
       In1 = .FALSE.
       In2 = .FALSE.
       InMatch = .FALSE.
       SolutionFound = .FALSE.
!
       InMatch(v) = .TRUE.
       In1(v) = .TRUE.
       In2(w) = .TRUE.
       t1(v) = w
       t2(w) = v   
!
! Llamar a la función recursiva para buscar el isomorfismo
!
       call RecursiveSearchEVC(N1,N2,AdjMatrix1,AdjMatrix2,NodeDegrees1,  &
                               NodeDegrees2,In1,In2,InMatch,              &
                               SolutionFound,t1,t2,NodeEVC1,NodeEVC2,Order)
!
! Devolver el resultado del isomorfismo
!
       VF2_EquivalentEVC = SolutionFound
!
       return
       end function VF2_EquivalentEVC
!
!======================================================================!
!
! Subrutina para la búsqueda recursiva de isomorfismos
! ---------------------------------------------------- 
! 
       recursive subroutine RecursiveSearchEVC(N1,N2,AdjMatrix1,       &
                                               AdjMatrix2,NodeDegrees1,&
                                               NodeDegrees2,In1,In2,   &
                                               InMatch,SolutionFound,  &
                                               t1,t2,NodeEVC1,NodeEVC2,&
                                               Order)
!
       implicit none
!
       integer,intent(in) :: N1, N2
       real(kind=8),intent(in) :: NodeEVC1(N1)
       real(kind=8),intent(in) :: NodeEVC2(N2)
       integer,intent(in) :: Order(N1)
       integer,intent(in) :: NodeDegrees1(N1)
       integer,intent(in) :: NodeDegrees2(N2)
       logical,intent(in) :: AdjMatrix1(N1,N1)
       logical,intent(in) :: AdjMatrix2(N2,N2)
       logical,intent(inout) :: In1(N1)
       logical,intent(inout) :: In2(N2)
       logical,intent(inout) :: InMatch(N1)
       logical,intent(inout) :: SolutionFound
       integer,intent(inout) :: t1(N1) ! Añadido para almacenar el estado de mapeo
       integer,intent(inout) :: t2(N2) ! Añadido para almacenar el estado de mapeo
!
! Variables locales
!
!~        logical :: CheckRules
       integer :: v, w, i
!
! Si se ha encontrado una solución, devolver
!
       if ( SolutionFound ) return
!
! Si todos los nodos están en el mapeo, se ha encontrado una solución
!
       if ( all(InMatch) ) then
         SolutionFound = .TRUE.
         return
       end if
!
! Seleccionar el siguiente candidato
!
       v = 0
       do i = 1, N1
         if ( .NOT. InMatch(order(i)) ) then
           v = order(i)
           exit
         end if
       end do
!
! Probar todos los nodos posibles en el grafo 2
!
       do w = 1, N2
         if ( .not. In2(w) ) then
!
! Comprobar las reglas de corte y factibilidad
!
           if ( CheckRulesEVC(v,w,N1,N2,AdjMatrix1,AdjMatrix2,         &
                              NodeDegrees1,NodeDegrees2,t1,t2,         &
                                         NodeEVC1,NodeEVC2)       ) then
!
! Añadir los nodos al mapeo
!
             InMatch(v) = .TRUE.
             In1(v) = .TRUE.
             In2(w) = .TRUE.
             t1(v) = w
             t2(w) = v
!
! Llamar a la función recursiva para el siguiente nivel de la búsqueda
!
             call RecursiveSearchEVC(N1,N2,AdjMatrix1,AdjMatrix2,      &
                                     NodeDegrees1,NodeDegrees2,In1,In2,&
                                     InMatch,SolutionFound,t1,t2,      &
                                     NodeEVC1,NodeEVC2,Order)
!
! Si se ha encontrado una solución, devolver
!
             if ( SolutionFound ) return
!
! Eliminar los nodos del mapeo
!
             InMatch(v) = .FALSE.
             In1(v) = .FALSE.
             In2(w) = .FALSE.
             t1(v) = 0
             t2(w) = 0
!
           end if
         end if
       end do
!
       return
       end subroutine RecursiveSearchEVC
!
!======================================================================!
!
! Función para verificar las reglas de corte y factibilidad
! ---------------------------------------------------------
!
       logical function CheckRulesEVC(v,w,N1,N2,AdjMatrix1,AdjMatrix2, &
                                      NodeDegrees1,NodeDegrees2,t1,t2, &
                                      NodeEVC1,NodeEVC2)
!
       implicit none
!
       integer, intent(in) :: v,w
       integer, intent(in) :: N1,N2
       logical, intent(in) :: AdjMatrix1(N1,N1)
       logical, intent(in) :: AdjMatrix2(N2,N2)
       integer, intent(in) :: t1(N1)
       integer, intent(in) :: t2(N2)
       integer,intent(in) :: NodeDegrees1(N1)
       integer,intent(in) :: NodeDegrees2(N2)
       real(kind=8),intent(in) :: NodeEVC1(N1)
       real(kind=8),intent(in) :: NodeEVC2(N2)

!
! Variables locales
!
       integer :: i
!~ !
!~ !
!~ !
       CheckRulesEVC = .FALSE.
!
! Regla 1: Los nodos v y w deben tener el mismo grado
!
       if ( .NOT. (NodeDegrees1(v) == NodeDegrees2(w)) ) return
       if ( .NOT. abs(NodeEVC1(v)-NodeEVC2(w)) .lt. 1.0E-6 ) return
!
! Regla 2: Para cada nodo adyacente a v en el grafo 1, debe existir un nodo adyacente a w en el grafo 2
!
       do i = 1, N1
         if ( AdjMatrix1(v,i) .and. t1(i) > 0 ) then
!~ !write(*,*) 'BACKTRACE'
!~ !write(*,*) t1(i),t2(t1(i))
           if ( .NOT. AdjMatrix2(w,t1(i)) ) return
         end if
       end do
!
! Regla 3: Para cada nodo adyacente a w en el grafo 2, debe existir un nodo adyacente a v en el grafo 1
!
       do i = 1, N2
         if ( AdjMatrix2(w,i) .and. t2(i) > 0 ) then
           if ( .NOT. AdjMatrix1(v,t2(i)) ) return
         end if
       end do
!
! Regla 4: Los nodos v y w deben tener el mismo número de nodos adyacentes no mapeados
!
       if ( .NOT. (count(AdjMatrix1(v,:) .and.                         &
                  (t1 == 0)) == count(AdjMatrix2(w,:) .and.            &
                                                    (t2 == 0))) ) return
!
! Regla 5: Para cada nodo no mapeado adyacente a v en el grafo 1, debe existir un nodo no mapeado adyacente a w en el grafo 2
!
       do i = 1, N1
         if (AdjMatrix1(v,i) .and. t1(i) == 0) then
           if (.not. any(AdjMatrix2(w,:) .and. (t2 == 0))) return
         end if
       end do
!
! Regla 6: Para cada nodo no mapeado adyacente a w en el grafo 2, debe existir un nodo no mapeado adyacente a v en el grafo 1
!
       do i = 1, N2
         if (AdjMatrix2(w,i) .and. t2(i) == 0) then
           if (.not. any(AdjMatrix1(v,:) .and. (t1 == 0))) return
         end if
       end do
!
! Devolver el resultado de las reglas de corte y factibilidad
!
       CheckRulesEVC = .TRUE.
!
       return
       end function CheckRulesEVC
!
!======================================================================!
!
       subroutine DFS(source,nat,adj,num)
!
! Input/output variables
!
       logical,dimension(nat,nat),intent(in)  ::  adj       !  Adjacency matrix
       integer,dimension(nat),intent(out)     ::  num       !
       integer,intent(in)                     ::  source    !  Source vertex
       integer,intent(in)                     ::  nat       !  Number of nodes
!
! Local variables
!
       logical,dimension(nat)                 ::  visited   !
       integer                                ::  inum      !
!
! Initialize variables
!
       num(:) = 0
       inum   = 0
!        
! Perform depth-first search
!
       visited(:)      = .FALSE.
       visited(source) = .TRUE.
!
       inum = inum + 1
       num(inum) = source
!
       call DFSrecursion(source,nat,adj,visited,num,inum)
!
!~        write(*,*) 'Order of vertices visited:'
!~        do i = 1, nat
!~          write(*,*) 'Vertex ',i,' was visited at order ',num(i)
!~        end do
!~        write(*,*)
!
       return
       end subroutine DFS
!
!======================================================================!
!
       recursive subroutine DFSrecursion(source,nat,adj,visited,       &
                                         num,inum)
!
! Input/output variables
!
       logical,dimension(nat,nat),intent(in)     ::  adj       !  Adjacency matrix
       logical,dimension(nat),intent(inout)      ::  visited   !
       integer,dimension(nat),intent(inout)      ::  num       !
       integer,intent(inout)                     ::  inum      !
       integer,intent(in)                        ::  source    !  Number of nodes
       integer,intent(in)                        ::  nat       !  Number of nodes
!
! Local variables
!
       integer                                   ::  i         !
! 
! DFS
!
!~ write(*,*) 'source vertex:',source
!~        if ( .NOT. visited(source) ) then
!~          visited(source) = .TRUE.
!~          inum = inum + 1
!~          num(source) = inum
!~        end if
!
       do i = 1, nat
         if ( adj(source,i) .and. (.NOT.visited(i)) ) then
!
           inum = inum + 1
           num(inum) = i
!
           visited(i) = .TRUE.
           call DFSrecursion(i,nat,adj,visited,num,inum)
!
         end if
       end do
!
! Alternative
!
!procedure DFS(G,v):
!label v as discovered
!for all edges from v to w in G.adjacentEdges(v) do
!{
!       if (v < w) add edge(v,w) to output edges
!       if vertex w is not labeled as discovered then
!          recursively call DFS(G,w)
!}
!
!~        visited(source) = .TRUE.
!~ !
!~        do i = 1, nat
!~          if ( adj(source,i) ) then
!~ !
!~ !
!~            if ( .NOT. visited(i) ) then
!~              inum = inum + 1
!~              num(i) = inum
!~            end if
!~ !
!~            if ( num(source) .lt. num(i) ) then
!~ write(*,*) 'adding vertex',source,i,'-->',num(source),num(i)!,'-->',min(num(source),num(i)),max(num(source),num(i))
!~ !
!~              iedge = iedge + 1
!~ !
!~              efrom(iedge) = num(source)
!~              eto(iedge)   = num(i) 
!~            end if
!~ !
!~            if ( .NOT. visited(i) ) then
!~              call DFSrecursion(i,nat,adj,visited,evisited,num,inum,      &
!~                                efrom,eto,order,iedge)
!~            end if
!~ !
!~          end if
!~        end do
!
! First trial
!
!~ ! Check if the current vertex has been visited
!~        if (.not. visited(current)) then
!~ ! Mark the current vertex as visited
!~          visited(current) = .true.  
!~                inum = inum + 1
!~                num(i) = inum 
!~                visited(i) = .TRUE.                              
!~        end if
!~ ! Push the unvisited neighbors onto the stack
!~          do i = 1, nat
!           if ( adj(current,i) .and. (.not.visited(i)) ) then
!~            if ( adj(current,i) ) then
!~ write(*,*) 'adding vertex',current,i
!~ ! Assign the order number to the current vertex
!~              if ( .NOT. visited(i) ) then

!~                call DFSrecursion(i,nat,adj,visited,num,inum,efrom,eto,order)
!~              end if
!~            end if
!~          end do
!
       return
       end subroutine DFSrecursion
!
!======================================================================!
!
       end module isomorphism
!
!======================================================================!

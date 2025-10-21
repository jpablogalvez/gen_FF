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
                   findarcycles,                                       &
                   bonded2dihe,                                        &
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
       subroutine genrigidlist(nat,r,wiberg,coord,cycles,ncycle,       &
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
       integer,intent(in)                          ::  nat     !  Number of nodes
       integer,intent(in)                          ::  r       !  Graph rank
       integer,intent(in)                          ::  mcycle  !  Number of cycles
!
! Local variables
! 
       real(kind=8)                                ::  daux    !
       logical                                     ::  lcheck  !
       integer,dimension(4)                        ::  ivaux   !
       integer                                     ::  i,j,k   !
       integer                                     ::  ii,jj   !
!
       real(kind=8),parameter                      ::  pi =  4*atan(1.0_8) 

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
!  if Wiberg index is greater than 1.5 then it is a double bond
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
           else if ( wiberg(j,i) .ge. 1.5d0 ) then
             lrigid(i,j) = .TRUE.
             lrigid(j,i) = .TRUE.
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
       integer                                         ::  i,j      !  Index
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
       end module genfftools
!
!======================================================================!

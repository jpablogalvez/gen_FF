!======================================================================!
!
       program evaluation
!
       use timings
       use lengths,       only:  leninp,lentag,lenlab
       use units,         only:  uniinp,uniic,unideps,                 &
                                 unitop,uninb,unicorr
!
       use datatypes
!
       use printings
!
       use gromacs_files
       use g16_files
       use graphtools
       use genforcefield
!
       implicit none
!
! Input variables
!
       character(len=leninp)                           ::  inp      !
       character(len=leninp)                           ::  ref      !
       character(len=leninp)                           ::  intop    !
       character(len=leninp)                           ::  topout   !
       character(len=leninp)                           ::  qmout    !
       character(len=lenlab)                           ::  resname  !
       integer                                         ::  knei     !  Maximum distance from source node
       integer                                         ::  iroute   !  
       logical                                         ::  fsymm    !
       logical                                         ::  debug    !
!
       character(len=leninp)                           ::  geo      !
       character(len=leninp)                           ::  xyz      !
       character(len=leninp)                           ::  topnb    !
       character(len=leninp)                           ::  topcorr  !
       character(len=leninp)                           ::  bas      !
       character(len=leninp)                           ::  topbas   !
       character(len=lentag)                           ::  ext      !
       character(len=lentag)                           ::  topext   !
       character(len=lentag)                           ::  formt    !  QC  output format
       logical                                         ::  fqmout   !
!
       character(len=20)                               ::  meth     !
       character(len=20)                               ::  basis    !
       character(len=30)                               ::  disp     !
       integer                                         ::  chrg     !
       integer                                         ::  mult     !
!
       real(kind=8),dimension(:,:),allocatable         ::  coord    !  Atomic coordinates
       real(kind=8),dimension(:),allocatable           ::  mass     !  Atomic masses
       character(len=lenlab),dimension(:),allocatable  ::  lab      !  Atomic labels
       integer,dimension(:),allocatable                ::  znum     !  Atomic number
       integer                                         ::  nat      !  Number of atoms
!
! Graph properties
!       
       real(kind=8),dimension(:,:),allocatable         ::  wiberg   !  Wiberg bond index matrix
       integer,dimension(:,:),allocatable              ::  mindis   !  Minimum distance matrix
       integer,dimension(:,:),allocatable              ::  cycles   !  Cycles information
       integer,dimension(:),allocatable                ::  ncycle   !  Number of atoms in each cycle
       integer,dimension(:),allocatable                ::  ideg     !  
       integer                                         ::  mcycle   !  Number of cycles
       integer                                         ::  rank     !  Maximum number of cycles
       integer                                         ::  nedge    !  Number of edges (bonds)
       logical,dimension(:,:),allocatable              ::  adj      !  Boolean adjacency matrix
       logical,dimension(:,:),allocatable              ::  lcycle   !  Bonds belonging to rings
       logical,dimension(:,:),allocatable              ::  lrigid   !  Double/triple bonds       
! 
! Topology information
!
        type(grotop)                                   ::  reftop   !  New topology
        type(grotop)                                   ::  top      !  New topology
        type(dihedrals)                                ::  dihe     !  Dihedrals
! 
! Atom types information
!
       character(len=lenlab),dimension(:),allocatable  ::  newlab   !  New atom type name
       logical,dimension(:,:),allocatable              ::  eqv      !  Equivalent atoms relationships
       integer,dimension(:),allocatable                ::  ilab     !  Integers associated to initial labels
       integer,dimension(:),allocatable                ::  eqlab    !  Integers associated to equivalent nuclei
       integer,dimension(:),allocatable                ::  idat     !  Integers associated to refined labels
       integer,dimension(:),allocatable                ::  itype    !  Integer associated to atomtypes
       integer                                         ::  nlab     !  Number of initial labels
       integer                                         ::  nidat    !  Number of refined labels
       integer                                         ::  ntype    !  Number of atom types
!
! Auxiliary variables
!
       logical                                         ::  match    !
       integer                                         ::  lin      !  
       integer                                         ::  lfin     ! 
       integer                                         ::  io       !  Status
       integer                                         ::  i,j,k    !  Indexes
!
! Declaration of time control variables
!
       integer                                         ::  t1,t2    !  Wall times
!
! Printing header
! ---------------
!
       write(*,*)
       call print_start()
!
! Initializing line lengths
!
       lin  = 35
       lfin = 70
!
! Initializing timings
!
       call system_clock(count_max=count_max,count_rate=count_rate)
       call system_clock(t1)  
!
       tcpu  = 0.0d0
!
! Reading command line options
!
       call command_line(inp,ref,intop,topout,qmout,knei,iroute,       &
                         formt,meth,basis,disp,chrg,mult,resname,      &
                         fsymm,debug)
!
! Defaults
!
      top%def%nbfunc  = 1
      top%def%comrule = 2
      top%def%genpair = 'yes'
      top%def%fudgelj = 0.5d0
      top%def%fudgeqq = 0.5d0/0.6d0
!
      top%mol%resname = resname
      top%mol%nrexcl  = 3
!
! Printing summary of the input information
!
       call print_title(6,1,'Input information','-')
       write(*,*)
       call line_str(6,2,'General input file name',lin,':',            &
                     trim(inp),lfin)
       call line_int(6,2,'K-Nearest neighbor',lin,':','I3',knei,lfin)
       write(*,*)
!
       call FLUSH()
!
! Reading input file
! ------------------
!
! Extracting basename from input file name
!
       if ( SCAN(inp,'.',BACK=.TRUE.) .ne. 0 ) then
         bas = inp(:SCAN(inp,'.',BACK=.TRUE.)-1)
         ext = inp(SCAN(inp,'.',BACK=.TRUE.)+1:)
       else
         bas = inp
         ext = ''
       end if
!
! Cleaning basename
!
       if ( bas(len_trim(bas)-2:) .eq. '_of' ) then
         bas = bas(:len_trim(bas)-3)
       end if
!
       xyz = trim(bas)//'.xyz'
!
! Extracting topology basename from topology input file name  
!
       if ( len_trim(topout) .gt. 0 ) then
!
         if ( SCAN(topout,'.',BACK=.TRUE.) .ne. 0 ) then
           topbas = topout(:SCAN(topout,'.',BACK=.TRUE.)-1)
           topext = topout(SCAN(topout,'.',BACK=.TRUE.)+1:)
         else
           topbas = topout
           topext = '.top'
         end if
!
       else
!
         topout = trim(bas)//'_GMX.top'
         topbas = bas
         topext = '.top'
!
       end if
!
       topnb   = trim(topbas)//'_nb.top'
       topcorr = trim(topbas)//'_template.top'
! 
! Checking if QC output is present in current file
!
       if ( len_trim(qmout) .gt. 0 ) then
!
         INQUIRE(FILE=trim(qmout),EXIST=fqmout)
!
       else
!
         qmout = trim(bas)//'_of'
!
         INQUIRE(FILE=trim(qmout)//'.log',EXIST=fqmout)
!
         if ( .NOT. fqmout ) then
           INQUIRE(FILE=trim(qmout)//'.out',EXIST=fqmout)
           qmout = trim(qmout)//'.out'
         else
           qmout = trim(qmout)//'.log'
         end if
!
       end if
!
       if ( fqmout ) then
         write(*,*)
         write(*,'(2X,68("*"))')
         write(*,'(3X,A)') 'NOTE:  QC output found in working directory'
         write(*,*)
         write(*,'(3X,A)') 'Force field will be generated from QM information'
         write(*,'(3X,A)') 'QC output name : '//trim(qmout)
         write(*,'(2X,68("*"))')
         write(*,*)  
       end if
!
! Checking if a reference topology is provided
!
       if ( len_trim(intop) .gt. 0 ) then
!
         call read_top(reftop,itype,intop,uniinp)
!
         nat = reftop%nat       
!
       else if ( .NOT. fqmout ) then 
!
! If only coordinates file is provided then generate QC input
!
         write(*,'(2X,68("*"))')
         write(*,'(3X,A)') 'NOTE: XYZ Files only contain coordinat'//  &
                                                        'es information'
         write(*,*) 
         write(*,'(3X,A)') '      Generating QC input in format '//    &
                                                             trim(formt)
         write(*,'(2X,68("*"))')
         write(*,*) 
!
         open(unit=uniinp,file=trim(xyz),action='read',status='old',   &
              iostat=io)
         if ( io .ne. 0 ) call print_missinp(xyz)
!
         read(uniinp,*) nat
         read(uniinp,*)
!
         allocate(coord(3,nat),lab(nat))
!
! Reading geometry and initial labels
!
         do i = 1, nat
           read(uniinp,*) lab(i),coord(:,i)
         end do
!
         close(uniinp)
!
! Generating QC input
!
         call geninp(nat,lab,coord,bas,formt,meth,basis,disp,chrg,mult)
         GO TO 1000
!
       end if
!
! Extracting molecular information
! 
       if ( .NOT. fqmout ) then 
!
! If only initial topology is provided then extract FF information
!
         call read_bond(reftop%bonded,intop,uniinp)
         call read_angle(reftop%bonded,intop,uniinp)
         call read_dihe(reftop%bonded,intop,uniinp)
!
       else
!
! If QM output is provided then extract molecular information
!
!~          call chk_qmout(qmout,qmformt) ! TODO: guess qm output format
!
         call chk_log(qmout,nat)
!
         allocate(coord(3,nat),lab(nat),znum(nat),mass(nat))
!
         call read_log(qmout,nat,coord,lab,znum,mass) !  TODO: check inconsistency between topology and qm output
!
         top%nat = nat
!
! Allocating memory
!
         allocate(wiberg(nat,nat))
         allocate(lcycle(nat,nat),lrigid(nat,nat))
!
!  Reading Wiberg bond index matrix  ! TODO: only supports gaussian QC output
!
         call readwiberg(qmout,nat,wiberg)  
!
       end if
!
!  If reference geometry is provided extract atomic positions and labels
!
       geo = qmout
!
       if ( len_trim(ref) .ne. 0 ) then
         geo = ref
!
         open(unit=uniinp,file=trim(ref),action='read',                &
              status='old',iostat=io)
!
         if ( io .ne. 0 ) call print_missinp(ref)
!
         read(uniinp,*) nat  !  TODO: check inconsistency between topology, reference input and qm output
         read(uniinp,*)
!
! Reading geometry and initial labels
!
         do i = 1, nat  
           read(uniinp,*) lab(i),coord(:,i)  !  TODO: check inconsistency between reference input and qm output
         end do
!
         close(uniinp) 
       end if
!
! Allocating general variables
!
       allocate(adj(nat,nat),eqv(nat,nat),idat(nat))
       allocate(ilab(nat),newlab(nat),eqlab(nat))
       allocate(mindis(nat,nat),ideg(nat))
!
! Generation of the adjacency matrix
! ----------------------------------  !  TODO: check bonds are consistent with input topology
!
       if ( .NOT. fqmout ) then 
!
! If only initial topology is provided then generate topology from bonds
!
         call bonds2adj(reftop%bonded%nbond,reftop%bonded%ibond,nat,adj)
!
         allocate(lab(nat))
!
         lab(:) = reftop%atom%attype(:)
!
       else
!
! If QM output is provided then generate topology from Wiberg matrix
!
         call wiberg2adj(nat,wiberg,adj)
!
       end if
!
! Computing graph properties
! -------------------------- 
!
!   Computing degrees
!
       ideg(:) = calcdegundir(nat,adj)
!
! Computing edges and cycle rank
!
       nedge = 0
       do i = 1, nat
         nedge = nedge + ideg(i)
       end do
       nedge = nedge/2
!
       rank = nedge - nat + 1 
!
       allocate(cycles(rank,nat),ncycle(rank)) 
!
!  Cycles identification
!
       call findcycle(nedge,nat,rank,adj,mcycle,ncycle,cycles)
!
       if ( debug ) then
         write(*,*) 'Cycle information'
         write(*,*) '-----------------'
         write(*,*) ' Number of cycles',mcycle
         do i = 1, mcycle
           write(*,*) 'Cycle Number', i, ':',(cycles(i,j),j=1,ncycle(i))
         end do
         write(*,*)
       end if
!
!  Computing minimum distance matrix  
!
       call dijkunwundir(nat,adj,mindis)
!
       if ( debug ) then
         write(*,*) 'Minimum distance matrix'
         write(*,*) '-----------------------'
         do i = 1, nat
           write(*,'(20(1X,I10))') (mindis(i,j),j=1,nat)
         end do
         write(*,*)
       end if
!
       if ( maxval(mindis) .le. knei ) then
         write(*,'(2X,68("="))')
         write(*,'(3X,A)') 'NOTE:  Largest minimum distance shorte'//  &
                                               'r or equal to k-nearest'
         write(*,'(3X,A)') '        neighbour distance'
         write(*,*)
         write(*,'(3X,A)') 'Result will converge to the Chemically'//  &
                                           ' Equivalent Atoms algorithm'
         write(*,'(2X,68("="))')
         write(*,*)
       end if  
!
       if ( fqmout ) then
!
!  Generating cyclelist with bonds belonging to cycles
!
         call gencyclelist(nat,rank,adj,mcycle,ncycle,cycles,lcycle)
!
!  Generating rigidlist with rigid bonds
!
         call genrigidlist(nat,wiberg,lcycle,lrigid)   
!
       end if
!
! Setting initial partitions  ! TODO: incorporate isotopes information through mass
! --------------------------
!             
       if ( len_trim(intop) .eq. 0 )  then
         call strnormlabels(nat,lab,ilab,nlab)         
         allocate(itype(nat))
         itype(:) = ilab(:)
       else
         ilab(:) = itype(:)
         nlab    = reftop%attype%ntype
       end if
!
! Atom typing strategy
! --------------------
!
       if ( knei .le. 0 ) then
!
! Chemically equivalent atoms strategy
!
         call eqvatms(nat,adj,ilab,nlab,eqv)  
!
       else
!
! k-Nearest neighbors algorithm
!
         call nearnei(nat,adj,eqv,nlab,ilab,knei,mindis,               &
                      rank,mcycle,ncycle,cycles,debug)
!
       end if
!
! Generating atom types labels
!
       call gentypes(nat,eqv,lab,newlab,ilab,nlab,eqlab,idat,nidat,    &
                     debug)
!
       top%nat = nat
!
       allocate(top%attype%atname(nat),top%attype%bond(nat),           &
                top%attype%ptype(nat),top%attype%atnum(nat),           &
                top%attype%mass(nat),top%attype%charge(nat),           &
                top%attype%sig(nat),top%attype%eps(nat))
!
       allocate(top%atom%attype(nat),top%atom%residue(nat),            &
                top%atom%atom(nat),top%atom%mass(nat),                 &
                top%atom%charge(nat),top%atom%cgnr(nat),               &
                top%atom%atnr(nat),top%atom%resnr(nat),                &
                top%atom%itype(nat))
!
       if ( len_trim(intop) .gt. 0 ) then
         top%def = reftop%def
         top%mol = reftop%mol
       end if
!
! Generating new atomtypes section
!
       top%attype%ntype     = nidat
       top%attype%ptype(:)  = 'A'
       top%attype%mass(:)   = 0.0d0
       top%attype%charge(:) = 0.0d0
       top%attype%sig(:)    = 0.0d0
       top%attype%eps(:)    = 0.0d0
       top%attype%atnum(:)  = 0
!
       k = 0
       do i = 1, nat
         match = .FALSE.
         do j = 1, k
           if ( newlab(i) .eq. top%attype%atname(j) ) then
             match = .TRUE.
             top%atom%itype(i) = j
             exit
           end if
         end do
         if ( .not. match ) then
!
           k = k + 1
!
           top%atom%itype(i) = k
!
           top%attype%atname(k) = newlab(i)
           top%attype%bond(k)   = newlab(i) 
!
           if ( len_trim(intop) .gt. 0 ) then
             top%attype%ptype(k)  = reftop%attype%ptype(itype(i))
             top%attype%mass(k)   = reftop%attype%mass(itype(i))
             top%attype%charge(k) = reftop%attype%charge(itype(i))
             top%attype%sig(k)    = reftop%attype%sig(itype(i))
             top%attype%eps(k)    = reftop%attype%eps(itype(i))
             top%attype%atnum(k)  = reftop%attype%atnum(itype(i))
           end if
!
         end if
       end do
!
       itype(:) = top%atom%itype(:)
!
! Generating new atoms section
!
       top%atom%nat = nat
!
       if ( len_trim(intop) .gt. 0 ) then
         top%atom = reftop%atom
       else
!
         top%mol%resname = resname
!
         do i = 1, nat
           top%atom%residue(i) = top%mol%resname
           top%atom%atom(i)    = lab(i)
           top%atom%mass(i)    = 0.0d0
           top%atom%charge(i)  = 0.0d0
           top%atom%cgnr(i)    = i
           top%atom%atnr(i)    = i
           top%atom%resnr(i)   = 1
         end do
       end if
!
       top%atom%attype(:) = newlab(:)
!
! Force Field generation
! ----------------------
!
       if ( fqmout ) then
!
         call genffbonded(nat,coord,adj,ideg,lcycle,lrigid,            &
                          znum,top%bonded,dihe,iroute,debug)
!
       else
!
         allocate(top%bonded%bond(reftop%bonded%nbond),                &
                  top%bonded%kbond(reftop%bonded%nbond),               &
                  top%bonded%fbond(reftop%bonded%nbond),               &
                  top%bonded%ibond(2,reftop%bonded%nbond))     
!
         allocate(top%bonded%ang(reftop%bonded%nang),                  &
                  top%bonded%kang(reftop%bonded%nang),                 &
                  top%bonded%fang(reftop%bonded%nang),                 &
                  top%bonded%iang(3,reftop%bonded%nang))
!
         allocate(top%bonded%dihe(reftop%bonded%ndihe),                &
                  top%bonded%kdihe(reftop%bonded%ndihe),               &
                  top%bonded%fdihe(reftop%bonded%ndihe),               &
                  top%bonded%multi(reftop%bonded%ndihe),               &
                  top%bonded%idihe(4,reftop%bonded%ndihe))
!
         allocate(top%bonded%c0(reftop%bonded%ndihe),                  &
                  top%bonded%c1(reftop%bonded%ndihe),                  &
                  top%bonded%c2(reftop%bonded%ndihe),                  &
                  top%bonded%c3(reftop%bonded%ndihe),                  &
                  top%bonded%c4(reftop%bonded%ndihe),                  &
                  top%bonded%c5(reftop%bonded%ndihe))
!
         top%bonded = reftop%bonded 
!
         call bonded2dihe(reftop%bonded%ndihe,top%bonded,dihe,nat,adj)
!
       end if
!
! Symmetrizing force field terms
! ------------------------------
!
       call symffbonded(nat,nidat,idat,top%bonded,dihe,fsymm,debug)
!
! Printing force field
! --------------------
!
       call print_top(unitop,nat,itype,mindis,top,dihe,bas,geo,intop,topout,debug)
!
! Printing JOYCE input files
! --------------------------
!
!~        call print_deps
!~        call print_ic
!~        call print_step1
!
! Printing summary of the input information
!
       call print_title(6,1,'Output information','-')
       write(*,*)
!
       write(*,*) 'Atomtyping'
       write(*,*) '----------'
       do i = 1, nat
         write(*,'(2(1X,I3),2(1X,A5,1X,I3))') i,ilab(i),lab(i),        &
                                              eqlab(i),newlab(i),idat(i)
       end do
       write(*,*)
!
       if ( fqmout ) then
!
         write(*,*) 'Bond-stretching terms'
         write(*,*) '---------------------'
         do i = 1, top%bonded%nbond
           write(*,'(1X,I3,3(1X,A,1X,I3),1X,A,1X,F6.4)') i,'=',top%bonded%ibond(1,i), &
                                 '-',top%bonded%ibond(2,i)!,':',top%bonded%idbond(i),'=',top%bonded%bond(i)
         end do
         write(*,*)
!

         write(*,*) 'Angle-bending terms'
         write(*,*) '-------------------'
         do i = 1, top%bonded%nang
           write(*,'(1X,I3,4(1X,A,1X,I3),1X,A,1X,F9.4)') i,'=',top%bonded%iang(1,i), &
                     '-',top%bonded%iang(2,i),'-',top%bonded%iang(3,i)!,':',top%bonded%idang(i),'=',top%bonded%ang(i)
         end do
         write(*,*)
!
         if ( dihe%ndihe .gt. 0 ) then
           write(*,*) 'Dihedral terms'
           write(*,*) '--------------'
         end if
!
         if ( dihe%nimpro .gt. 0 ) then
           write(*,*) '; Impropers o.o.p'
           do i = 1, dihe%nimpro
             write(*,'(1X,I3,5(1X,A,1X,I3),1X,A,1X,F9.4)') i,'=',dihe%iimpro(1,i), &
     '-',dihe%iimpro(2,i),'-',dihe%iimpro(3,i),'-',dihe%iimpro(4,i)!,':',dihe%idimpro(i),'=',dihe%dimpro(i)
           end do
           write(*,*)
         end if
!
         if ( dihe%ninv .gt. 0 ) then
           write(*,*) '; Inversion dihedrals'
           do i = 1, dihe%ninv
             write(*,'(1X,I3,5(1X,A,1X,I3),1X,A,1X,F9.4)') i,'=',dihe%iinv(1,i), &
     '-',dihe%iinv(2,i),'-',dihe%iinv(3,i),'-',dihe%iinv(4,i)!,':',dihe%idinv(i),'=',dihe%dinv(i)
           end do
           write(*,*)
         end if
!
         if ( dihe%nrigid .gt. 0 ) then
           write(*,*) '; Impropers on double bonds/aromatic cycles'
           do i = 1, dihe%nrigid
             write(*,'(1X,I3,5(1X,A,1X,I3),1X,A,1X,F9.4)') i,'=',dihe%irigid(1,i), &
     '-',dihe%irigid(2,i),'-',dihe%irigid(3,i),'-',dihe%irigid(4,i)!,':',dihe%idrigid(i),'=',dihe%drigid(i)
           end do
           write(*,*)
         end if
!
         if ( dihe%nflexi .gt. 0 ) then
           write(*,*) '; Flexible'
           do i = 1, dihe%nflexi
             write(*,'(1X,I3,5(1X,A,1X,I3),1X,A,1X,F9.4)') i,'=',dihe%iflexi(1,i),  &
      '-',dihe%iflexi(2,i),'-',dihe%iflexi(3,i),'-',dihe%iflexi(4,i)!,':',dihe%idflexi(i),'=',dihe%dflexi(i)
           end do
           write(*,*)
         end if
!        
       end if
!
! Deallocating memory
! -------------------
!
! Deallocating topology information
!
!       if ( fqmout ) then
!         deallocate(top%attype%atname,top%attype%bond,                 &
!                    top%attype%ptype,top%attype%atnum,                 &
!                    top%attype%mass,top%attype%charge,                 &
!                    top%attype%mass,top%attype%charge,                 &
!                    top%attype%sigma,top%attype%eps)
!
!         deallocate(top%atom%attype,top%atom%residue,top%atom%atom,    &
!                    top%atom%mass,top%atom%charge,top%atom%cgnr,       &
!                    top%atom%atnr,top%atom%resnr)
!
!         deallocate(top%bonded%bond,top%bonded%kbond,                  &
!                    top%bonded%ibond,top%bonded%ang,                 &
!                    top%bonded%kang,top%bonded%iang,               &
!                    top%bonded%dihe,top%bonded%kdihe,                  &
!                    top%bonded%idihe,top%bonded%multi)
!       end if
!
! Deallocating reference topology information
!
       if ( len_trim(intop) .gt. 0 ) then
!
         deallocate(reftop%attype%atname,reftop%attype%bond,           &
                    reftop%attype%ptype,reftop%attype%atnum,           &
                    reftop%attype%mass,reftop%attype%charge,           &
                    reftop%attype%sig,reftop%attype%eps)
!
         deallocate(reftop%atom%attype,reftop%atom%residue,            &
                    reftop%atom%atom,reftop%atom%mass,                 &
                    reftop%atom%charge,reftop%atom%cgnr,               &
                    reftop%atom%atnr,reftop%atom%resnr)
!
         if ( .NOT. fqmout ) then
!
           deallocate(reftop%bonded%bond,reftop%bonded%kbond,          &
                      reftop%bonded%ibond,reftop%bonded%ang,           &
                      reftop%bonded%kang,reftop%bonded%iang,           &
                      reftop%bonded%dihe,reftop%bonded%kdihe,          &
                      reftop%bonded%idihe,reftop%bonded%multi)
!
           deallocate(reftop%bonded%c0,reftop%bonded%c1,               &
                      reftop%bonded%c2,reftop%bonded%c3,               &
                      reftop%bonded%c4,reftop%bonded%c5)
!
         end if
!
       end if
!
! Deallocating local variables
!
!       deallocate(itype)
! 
!~        if ( fqmout ) then
!~          deallocate(mass,znum)
!~          deallocate(wiberg,lcycle,lrigid)
!~        end if
!~ !
!~        deallocate(adj,mindis,ideg)
!~        deallocate(cycles,ncycle)
!~        deallocate(ilab,eqlab,idat,newlab)
!~ !
!~ 1000   deallocate(coord,lab)
1000   continue
!
! Printing timings
!
       call system_clock(t2)
!
       tcpu = dble(t2-t1)/dble(count_rate)
!
       write(*,'(1X,90("="))')
       write(*,*)
       write(*,'(1X,A,3(X,I2,X,A))') 'Total CPU time        ',         &
                                      int(tcpu/(60*60)),'hours',       &
                                      mod(int(tcpu/60),60),'minutes',  &
                                      mod(int(tcpu),60),'seconds'  
       write(*,*)
!
! Printing finishing date 
!    
       call print_end()
!
       end program evaluation
!
!======================================================================!
!
       subroutine command_line(inp,ref,top,topout,qmout,knei,iroute,   &
                               formt,meth,basis,disp,chrg,mult,        & 
                               resname,fsymm,debug)
!
       use lengths, only: leninp,lencmd,lenarg,lentag,lenlab
       use printings
       use utils
!
       implicit none
!
! Input/output variables
!
       character(len=leninp),intent(out)  ::  inp      !  Input file name
       character(len=leninp),intent(out)  ::  top      !  Topology file name
       character(len=leninp),intent(out)  ::  topout   !  Output topology file name
       character(len=leninp),intent(out)  ::  qmout    !  QC input file name
       character(len=leninp),intent(out)  ::  ref      !  Reference xyz file name
       character(len=lentag),intent(out)  ::  formt    !
       character(len=lenlab),intent(out)  ::  resname  !
       integer,intent(out)                ::  knei     !  
       integer,intent(out)                ::  iroute   !  
       logical,intent(out)                ::  fsymm    !  
!
       character(len=20),intent(out)      ::  meth     !
       character(len=20),intent(out)      ::  basis    !
       character(len=30),intent(out)      ::  disp     !
       integer,intent(out)                ::  chrg     !
       integer,intent(out)                ::  mult     !       
!
       logical,intent(out)                ::  debug    !  Debug mode
!
! Local variables
!
       character(len=lencmd)              ::  cmd     !  Command executed
       character(len=lenarg)              ::  code    !  Executable name
       character(len=lenarg)              ::  arg     !  Argument read
       character(len=lenarg)              ::  next    !  Next argument to be read
       integer                            ::  io      !  Status
       integer                            ::  i       !  Index
!
! Setting defaults
!
       inp    = 'mol.xyz'
       ref    = ''
       top    = ''
       topout = ''
       qmout  = ''
!
       resname = 'MOL'
!
       formt  = 'gaussian'
       knei   = -1
       iroute = 1
!
       meth  = 'PBEPBE'
       basis = '6-31+G(D)'
       disp  = 'GD3BJ'
       chrg  = 0
       mult  = 1
!
       fsymm = .TRUE.
       debug = .FALSE.
!
! Reading command line
!
       call get_command_argument(0,code)
       call get_command(cmd)
! Reading command line options
       if ( command_argument_count().eq.0) return
!
       i = 1
       do
         call get_command_argument(i,arg)
         if ( len_trim(arg) == 0 ) exit
         i = i+1
         select case ( arg )
           case ('-f','-file','--file')
             call get_command_argument(i,inp,status=io)
             call check_arg(inp,io,arg,cmd)
             i = i + 1
!
           case ('-c','-r','-ref','-coord','-xyz','--coord','--xyz',   &
                                   '--ref','--coordinates','--xyz-file')
             call get_command_argument(i,ref,status=io)
             call check_arg(ref,io,arg,cmd)
             i = i + 1
!
           case ('-p','-reftop','-topref','--reftop','--topref',       &
                 '--input-topology')
             call get_command_argument(i,top,status=io)
             call check_arg(top,io,arg,cmd)
             i = i + 1
!
           case ('-t','-top','--top','-topout','-outtop','--topout',   &
                 '--outtop','--topology','--output-topology')
             call get_command_argument(i,topout,status=io)
             call check_arg(topout,io,arg,cmd)
             i = i + 1
!
           case ('-qc','-qm','--qc','--qm','-qmout','--qmout')
             call get_command_argument(i,qmout,status=io)
             call check_arg(qmout,io,arg,cmd)
             i = i + 1
!
           case ('-name','-resname','--resname')
             call get_command_argument(i,resname,status=io)
             call check_arg(resname,io,arg,cmd)
             i = i + 1
!
           case ('-fmt','--format')
             call get_command_argument(i,formt,status=io)
             call check_arg(formt,io,arg,cmd)
!
             formt = lowercase(formt)
             select case (trim(formt))
               case('g16','gaussian')
                 formt='gaussian'
               case('orca')
                 STOP 'ORCA format not implemented yet!'
               case default
                 write(*,*)
                 write(*,'(2X,68("="))')
                 write(*,'(3X,A)')    'ERROR:  Invalid value intro'//  &
                                             'duced for --format option'
                 write(*,*)
                 write(*,'(3X,2(A))') 'Unrecognised value     : ',     &
                                                             trim(formt)
                 write(*,*)
                 write(*,'(3X,A)')    'Please, to know the possibl'//  &
                                                     'e options execute'
                 write(*,*)
                 write(*,'(4X,2(A))')  code(1:len_trim(code)), ' -h'
                 write(*,'(2X,68("="))')
                 write(*,*)
                 call print_end()
             end select
!
             i = i + 1
!
           case ('-knei','--knei')
             call get_command_argument(i,next,status=io)
             read(next,*) knei
             i = i + 1
!
           case ('-iroute','--iroute')
             call get_command_argument(i,next,status=io)
             read(next,*) iroute
             i = i + 1
!
           case ('-meth','--method')
             call get_command_argument(i,meth,status=io)
             call check_arg(meth,io,arg,cmd)
             i = i + 1
!
           case ('-bas','--basis')
             call get_command_argument(i,basis,status=io)
             call check_arg(basis,io,arg,cmd)
             i = i + 1
!
           case ('-disp','--dispersion')
             call get_command_argument(i,disp,status=io)
             call check_arg(disp,io,arg,cmd)
             i = i + 1
!
           case ('-chrg','--charge')
             call get_command_argument(i,next,status=io)
             read(next,*) chrg
             i = i + 1
!
           case ('-mult','-multi','--multi','--multiplicity')
             call get_command_argument(i,next,status=io)
             read(next,*) mult
             i = i + 1
!
           case ('-s','-sym','-symm','--symmetrize','--symm','--sym')
             fsymm = .TRUE.
!
           case ('-nos','-nosym','-nosymm','--nosymmetrize',           &
                                                   '--nosymm','--nosym')
             fsymm = .FALSE.
!
           case ('-v','--debug','--verbose')
             debug = .TRUE.
           case ('-h','-help','--help')
             call print_help()
             call print_end()
           case default
             write(*,*)
             write(*,'(2X,68("="))')
             write(*,'(3X,A)')    'ERROR:  Unknown statements from'//  &
                                  ' command line'
             write(*,*)
             write(*,'(4X,A)')     trim(cmd)
             write(*,*)
             write(*,'(3X,2(A))') 'Unrecognised command-line option'// &
                                  '  :  ', arg
             write(*,'(3X,A)')    'Please, to request help execute'
             write(*,*)
             write(*,'(4X,2(A))')  code(1:len_trim(code)), ' -h'
             write(*,'(2X,68("="))')
             write(*,*)
             call print_end()
         end select
       end do
!
       return
       end subroutine command_line
!
!======================================================================!
!
       subroutine print_help()  ! TODO: check names from command_line routine
!
       implicit none
!
       write(*,'(1X,A)') 'Command-line options'
       write(*,'(1X,20("-"))')
       write(*,*)
       write(*,'(2X,A)') '-h,--help               Print usage info'//  &
                                                       'rmtion and exit'
       write(*,*)
       write(*,'(2X,A)') '-f,--file               Input file name'
       write(*,'(2X,A)') '-c,--coordinates        Input coordinate'//  &
                                                                's name'
       write(*,'(2X,A)') '-p,--input-topology     Input Gromacs to'//  &
                                                                'pology'
       write(*,'(2X,A)') '-t,--output-topology    Output Gromacs t'//  & 
                                                               'opology'
       write(*,'(2X,A)') '-qc,--qmout             QC input file name'
       write(*,*) 
       write(*,'(2X,A)') '-knei,--knei            K-nearest neighbor'
       write(*,'(2X,A)') '-iroute','--iroute      Functional form of the FF'
       write(*,*)
       write(*,'(2X,A)') '-fmt,--format           Format of QC input'
       write(*,'(2X,A)') '                         ( gaussian | orca )'
       write(*,'(2X,A)') '-meth,--method          QM method'
       write(*,'(2X,A)') '-bas,--basis            Basis set'
       write(*,'(2X,A)') '-disp,--dispersion      Dispersion correction'
       write(*,'(2X,A)') '-chrg,--charge          Molecular charge'
       write(*,'(2X,A)') '-mult,--multiplicity    Spin multiplicity'
       write(*,*)
       write(*,'(2X,A)') '-v,--verbose            Debug mode'
       write(*,*)
!
       return
       end subroutine print_help
!
!======================================================================!
!
! GENINP - GENerate INPut
!
! This subroutine 
!
       subroutine geninp(nat,lab,coord,bas,formt,meth,basis,disp,      &
                         chrg,mult)
!
       use lengths,  only: leninp,lentag,lenlab
       use units,    only: uniout
!
       implicit none
!
! Input/output variables
!
       character(len=leninp),intent(in)          ::  bas     !
       character(len=lenlab),dimension(nat)      ::  lab     !  Atomic labels
       real(kind=8),dimension(3,nat),intent(in)  ::  coord   ! 
       integer,intent(in)                        ::  nat     !
!
       character(len=lentag),intent(in)          ::  formt   !
       character(len=20)                         ::  meth    !
       character(len=20)                         ::  basis   !
       character(len=30)                         ::  disp    !
       integer,intent(in)                        ::  chrg    !
       integer,intent(in)                        ::  mult    !
!
! Local variables
!
       integer                                   ::  i       !
!
! Generating QC input
! -------------------
!
       select case (trim(formt))
         case('orca')
!
           open(unit=uniout,file=trim(bas)//'_of.inp',action='write')
!

!
           close(20)
!
         case('gaussian')
!
           if ( (trim(disp).eq.'none') .or. (len_trim(disp).eq.0) ) then
             disp = ''
           else
             disp = 'empiricaldispersion='//trim(disp)
           end if 
!
           open(unit=uniout,file=trim(bas)//'_of.com',action='write')
!
           write(uniout,'(A)') '%nprocshared=8'
           write(uniout,'(A)') '%chk='//trim(bas)//'_of.chk'
           write(uniout,'(A)') '#p opt freq=intmodes pop=NBORead g'//  &
                                                                 'finput'
           write(uniout,'(A)') trim(meth)//'/'//trim(basis)//' '//     &
                                                              trim(disp)
           write(uniout,*)
           write(uniout,'(A)') trim(bas)//' optimization and frequ'//  &
                                                                'encies'
           write(uniout,*)
           write(uniout,'(I1,1X,I1)') chrg,mult
           do i = 1, nat
             write(uniout,'(A5,3(1X,F12.6))') lab(i),coord(:,i)
           end do
           write(uniout,*)
           write(uniout,'(A)') '$nbo BNDIDX $end'
           write(uniout,*)
!
           close(uniout)
!
       end select
!
       return
       end subroutine geninp
!
!======================================================================!
!
! GENNBINP - GENerate NonBonded INPut
!
! This subroutine 
!
       subroutine gennbinp(nat,lab,coord,bas,top,formt,chrg,mult)
!
       use lengths,  only: leninp,lentag,lenlab
       use units,    only: uniout
!
       implicit none
!
! Input/output variables
!
       character(len=leninp),intent(in)          ::  bas     !
       character(len=leninp),intent(in)          ::  top     !
       character(len=lenlab),dimension(nat)      ::  lab     !  Atomic labels
       real(kind=8),dimension(3,nat),intent(in)  ::  coord   ! 
       integer,intent(in)                        ::  nat     !
!
       character(len=lentag),intent(in)          ::  formt   !
       integer,intent(in)                        ::  chrg    !
       integer,intent(in)                        ::  mult    !
!
! Local variables
!
       integer                                   ::  i       !
!
! Generating nb MM input
! ----------------------
!
       select case (trim(formt))
         case('orca')
!
           open(unit=uniout,file=trim(bas)//'_nb.inp',action='write')
!

!
           close(uniout)
!
         case('gaussian')
!
           open(unit=uniout,file=trim(bas)//'_nb.com',action='write')
!
           write(uniout,'(A)') '%nprocshared=8'
           write(uniout,'(A)') '%chk='//trim(bas)//'_nb.chk'
           write(uniout,'(A)') '# freq external="gromacs_link.sh -'//  &
                                                    'p '//trim(top)//'"'
           write(uniout,*)
           write(uniout,'(A)') trim(bas)//' nb frequencies'
           write(uniout,*)
           write(uniout,'(I1,1X,I1)') chrg,mult
           do i = 1, nat
             write(uniout,'(A5,3(1X,F12.6))') lab(i),coord(:,i)
           end do
           write(uniout,*)
!
           close(uniout)
!
       end select
!
       return
       end subroutine gennbinp
!
!======================================================================!
!
       subroutine chk_qmout(inp,lab)
!
       use lengths,    only: leninp,lentag,lenline
       use units,      only: uniinp
!
       use printings,  only: print_end,print_missinp
!
       implicit none
!
! Input/output variables
!
       character(len=leninp),intent(in)   ::  inp   !  Input file name
       character(len=lentag),intent(out)  ::  lab   !
!
! Local variables
!
       character(len=lenline)            ::  line  !
       integer                           ::  posi  !
       integer                           ::  io    !
!
! Reading information from QM output file
!
       open(unit=uniinp,file=trim(inp),action='read',                  &
            status='old',iostat=io)
!
       if ( io .ne. 0 ) call print_missinp(inp)
!
! Checking if QM output comes from a Gaussian16 calculation
!
!   Find the first non-blank line
       do
         read(uniinp,'(A)',iostat=io) line
         if ( len_trim(line) == 0 ) then
           cycle
         else
           line = adjustl(line)
           if ( line(:5) .ne. 'nohup' ) exit
         end if
       end do
!
       posi = scan(line,',')
       line = line(:posi-1)
       line = adjustl(line)
!
       if ( trim(line) .eq. 'Entering Gaussian System' ) then
         lab = 'g16'  
         close(uniinp)  
         return    
       end if
!
! Checking if QM output comes from a ORCA calculation
!
       rewind (uniinp)
!   Find the first non-blank line
       do
         read(uniinp,'(A)',iostat=io) line
         if ( len_trim(line) == 0 ) then
           cycle
         else
           line = adjustl(line)
           if ( line(:5) .ne. 'nohup' ) exit
         end if
       end do
!
       read(uniinp,'(A)') line
       line = adjustl(line)
!
       if ( trim(line) .eq. '* O   R   C   A *' ) then
         lab = 'orca'
         close(uniinp)  
         return    
       end if
!
! Error check
!
       write(*,*)
       write(*,'(2X,68("="))')
       write(*,'(3X,A)') 'ERROR:  Format of QM output not supported'
       write(*,*)
       write(*,'(3X,A)') 'Input file '//trim(inp)//' comes from an'//  &
                                                  ' unknown QM software'
       write(*,'(2X,68("="))')
       write(*,*)
       call print_end()       
!
       return
       end subroutine chk_qmout
!
!======================================================================!

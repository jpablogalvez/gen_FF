!======================================================================!
!
       module datatypes
! 
       use lengths
!
       implicit none
!
       type groinp
         character(len=32)                          ::  fname    !  Input file
         character(len=72)                          ::  title    !  Title
         integer                                    ::  nat      !  Number of atoms
         real(kind=8)                               ::  totm     !  Total mass 
         integer,dimension(:),allocatable           ::  renum    !  Residue number
         character(len=5),dimension(:),allocatable  ::  rename   !  Residue name
         character(len=5),dimension(:),allocatable  ::  atname   !  Atom name
         integer,dimension(:),allocatable           ::  atnum    !  Atom number
         real(kind=8),dimension(:),allocatable      ::  mass     !  Atom mass
         real(kind=8),dimension(:,:),allocatable    ::  coord    !  Coordinates
         real(kind=8),dimension(:,:),allocatable    ::  vel      !  Velocities
         real(kind=8),dimension(3)                  ::  latvec   !  Box vectors
       end type groinp  
!
       type grodefaults
         character(len=3)                           ::  genpair  !  Pairs generation flag
         real(kind=8)                               ::  fudgelj  !  Lennard-Jones fudge
         real(kind=8)                               ::  fudgeqq  !  Electrostatic fudge       
         integer                                    ::  nbfunc   !  Non-bonded function type
         integer                                    ::  comrule  !  Combination rule
       end type grodefaults 
!
       type groattype
         character(len=5),dimension(:),allocatable  ::  atname   !  Atom type
         character(len=5),dimension(:),allocatable  ::  bond     !  Bonded type
         character(len=1),dimension(:),allocatable  ::  ptype    !  Particle type
         real(kind=8),dimension(:),allocatable      ::  mass     !  Atom mass
         real(kind=8),dimension(:),allocatable      ::  charge   !  Atom charge
         real(kind=8),dimension(:),allocatable      ::  sig      !  Sigma parameter
         real(kind=8),dimension(:),allocatable      ::  eps      !  Epsilon parameter
         integer,dimension(:),allocatable           ::  atnum    !  Atomic number
         integer                                    ::  ntype    ! 
       end type groattype 
!
       type groatoms
         character(len=5),dimension(:),allocatable  ::  attype   !  Atom type
         character(len=5),dimension(:),allocatable  ::  residue  !  Residue name
         character(len=5),dimension(:),allocatable  ::  atom     !  Atom name
         real(kind=8),dimension(:),allocatable      ::  charge   !  Atom charge
         real(kind=8),dimension(:),allocatable      ::  mass     !  Atom mass
         integer,dimension(:),allocatable           ::  atnr     !  Atom number
         integer,dimension(:),allocatable           ::  resnr    !  Residue number
         integer,dimension(:),allocatable           ::  cgnr     !  Charge group number
         integer,dimension(:),allocatable           ::  itype    !  Charge group number
         integer                                    ::  nat      !  
       end type groatoms
!
       type gromolecule
         character(len=lenarg)                      ::  sysname  !
         character(len=lenlab)                      ::  resname  !
         integer                                    ::  nrexcl   !  Exclusions of nb interactions
         integer                                    ::  nmol     ! 
       end type gromolecule
!
       type grobonded
         character(len=50),dimension(:),allocatable ::  labbond  !
         character(len=50),dimension(:),allocatable ::  labang   !
         real(kind=8),dimension(:),allocatable      ::  bond     !  Equilibrium bond distance
         real(kind=8),dimension(:),allocatable      ::  ang      !  Equilibrium angle  
         real(kind=8),dimension(:),allocatable      ::  dihe     !  Equilibrium dihedral angle
         real(kind=8),dimension(:),allocatable      ::  kbond    !  Stretchings force constant
         real(kind=8),dimension(:),allocatable      ::  kang     !  Bending force constant
         real(kind=8),dimension(:),allocatable      ::  kdihe    !  Torsional barrier height
         real(kind=8),dimension(:),allocatable      ::  c0       !  Torsional series
         real(kind=8),dimension(:),allocatable      ::  c1       !  Torsional series
         real(kind=8),dimension(:),allocatable      ::  c2       !  Torsional series
         real(kind=8),dimension(:),allocatable      ::  c3       !  Torsional series
         real(kind=8),dimension(:),allocatable      ::  c4       !  Torsional series
         real(kind=8),dimension(:),allocatable      ::  c5       !  Torsional series
         integer,dimension(:,:),allocatable         ::  ibond    !  Bond index
         integer,dimension(:,:),allocatable         ::  iang     !  Angle index
         integer,dimension(:,:),allocatable         ::  idihe    !  Dihedral index
         integer,dimension(:),allocatable           ::  fbond    !  Bond type
         integer,dimension(:),allocatable           ::  fang     !  Angle type
         integer,dimension(:),allocatable           ::  fdihe    !  Dihedral type
         integer,dimension(:),allocatable           ::  sbond    !  Bond type
         integer,dimension(:),allocatable           ::  sang     !  Angle type
         integer,dimension(:),allocatable           ::  sdihe    !  Dihedral type
         integer,dimension(:),allocatable           ::  idbond   !  
         integer,dimension(:),allocatable           ::  idang    !  
         integer,dimension(:),allocatable           ::  multi    !  Torsional potential multiplicity
         integer                                    ::  nbond    !  Bonds number
         integer                                    ::  nang     !  Angles number
         integer                                    ::  ndihe    !  Torsions number
         integer                                    ::  nidbond  !  
         integer                                    ::  nidang   !  
         integer                                    ::  mdihe    !  Multiplicity
       end type grobonded
!
       type torsion
         real(kind=8)                               ::  vtor     !  Torsional potential
         real(kind=8)                               ::  phase    !  Phase
         integer                                    ::  stor     !  
         integer                                    ::  multi    !  Multiplicity
       end type torsion
!
       type flexible
         type(torsion),dimension(:),allocatable     ::  tor      !
         integer,dimension(4)                       ::  itor     !  
         integer                                    ::  ntor     !
       end type flexible
!
       type dihedrals
         type(flexible),dimension(:),allocatable    ::  flexi    !
         character(len=50),dimension(:),allocatable ::  labflexi !
         character(len=50),dimension(:),allocatable ::  labrigid !
         character(len=50),dimension(:),allocatable ::  labimpro !
         character(len=50),dimension(:),allocatable ::  labinv   !
         real(kind=8),dimension(:),allocatable      ::  dquad    !
         real(kind=8),dimension(:),allocatable      ::  dflexi   !
         real(kind=8),dimension(:),allocatable      ::  drigid   !
         real(kind=8),dimension(:),allocatable      ::  dimpro   !
         real(kind=8),dimension(:),allocatable      ::  dinv     !
         real(kind=8),dimension(:),allocatable      ::  kflexi   !
         real(kind=8),dimension(:),allocatable      ::  krigid   !
         real(kind=8),dimension(:),allocatable      ::  kimpro   !
         real(kind=8),dimension(:),allocatable      ::  kinv     !
         integer,dimension(:,:),allocatable         ::  iquad    !
         integer,dimension(:,:),allocatable         ::  iflexi   !
         integer,dimension(:,:),allocatable         ::  irigid   !
         integer,dimension(:,:),allocatable         ::  iimpro   !
         integer,dimension(:,:),allocatable         ::  iinv     !
         integer,dimension(:),allocatable           ::  fquad    !
         integer,dimension(:),allocatable           ::  fflexi   !
         integer,dimension(:),allocatable           ::  frigid   !
         integer,dimension(:),allocatable           ::  fimpro   !
         integer,dimension(:),allocatable           ::  finv     !
         integer,dimension(:),allocatable           ::  squad    !
         integer,dimension(:),allocatable           ::  sflexi   !
         integer,dimension(:),allocatable           ::  srigid   !
         integer,dimension(:),allocatable           ::  simpro   !
         integer,dimension(:),allocatable           ::  sinv     !
         integer,dimension(:),allocatable           ::  idflexi  !
         integer,dimension(:),allocatable           ::  idrigid  !
         integer,dimension(:),allocatable           ::  idimpro  !
         integer,dimension(:),allocatable           ::  idinv    !
         integer,dimension(:),allocatable           ::  mapquad  !
         integer                                    ::  nquad    !
         integer                                    ::  ndihe    !
         integer                                    ::  nflexi   !
         integer                                    ::  nrigid   !
         integer                                    ::  nimpro   !
         integer                                    ::  ninv     !
         integer                                    ::  nidflexi !
         integer                                    ::  nidrigid !
         integer                                    ::  nidimpro !
         integer                                    ::  nidinv   !
       end type dihedrals
!
       type grotop
         type(grodefaults)                          ::  def      !  Default parameters
         type(groattype)                            ::  attype   !  Atom types information
         type(groatoms)                             ::  atom     !  Definition of the molecule
         type(gromolecule)                          ::  mol      !  Definition of the molecule
         type(grobonded)                            ::  bonded   !  Bonded interactions
         integer                                    ::  nat      !
         integer                                    ::  nstiff   !
         integer                                    ::  nsoft    !
       end type grotop 
!
       end module datatypes
!
!======================================================================!

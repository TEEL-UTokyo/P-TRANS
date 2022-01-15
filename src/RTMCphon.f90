!>
!! @file RTMCphon.f90
!! @brief Ray-tracing Monte-Carlo phonon simulator.
!! @copyright (c) 2015 Junichiro Shiomi All Rights Reserved.
!! @version 1.02
!! Redistribution or modification of the source or binary file
!! is not permitted without permission of the copyright holder.

!!$!>
!!$!! @brief random numbers for DEBUG
!!$module MyRandom
!!$  use omp_lib
!!$  implicit none
!!$
!!$  integer, parameter :: max_threads = 64
!!$  integer(8) :: x(0:max_threads)
!!$
!!$  integer(8), parameter :: A = 1664525_8
!!$  integer(8), parameter :: C = 1013904223_8
!!$  integer(8), parameter :: M = 2147483647_8
!!$
!!$contains
!!$  subroutine random_seed( size, put )
!!$    integer, intent(out), optional :: size
!!$    integer, intent(in), optional :: put(*)
!!$    integer :: ithread
!!$
!!$    if( present(size) ) size=1
!!$    if( present(put) ) then
!!$      ithread = omp_get_thread_num()
!!$
!!$      x(ithread) = put(1)
!!$    end if
!!$  end subroutine random_seed
!!$
!!$  subroutine random_number( r )
!!$    real(8), intent(out) :: r
!!$    integer :: ithread
!!$
!!$    ithread = omp_get_thread_num()
!!$    x(ithread) = mod( x(ithread)*A + C, M )
!!$    r = dble(x(ithread))/M
!!$  end subroutine random_number
!!$end module MyRandom

!>
!! @brief module for ray-tracing monte-carlo phonon simulation.
module RTMC
#ifdef MPI
  use mpi
#endif
  !!$  use MyRandom ! for DEBUG
  implicit none

  integer :: mpi_ierr, mpi_size=1, mpi_rank=0
  !!$  integer, parameter :: fd_trace = 0 ! not output trace
  integer, parameter :: fd_trace = 1001 ! output trace for debug
  integer :: particle_trace = 1000 ! particles for output trace
  integer :: interval_trace ! interval of output trace

  !< mathematical constants.
  real(8), parameter :: M_PI = 3.14159265358979323846d0 !< pi
  real(8), parameter :: M_TWOPI = M_PI*2.0d0 !< pi*2
  real(8), parameter :: M_PI_2 = M_PI*0.5d0 !< pi/2

  !< physics constants.
  real(8), parameter :: hbar   = 1.054571d-34 !< plank constant
  real(8), parameter :: kB     = 1.3806503d-23 !< boltzman constant
  real(8), parameter :: SMALL_VAL     = 1.0d-5 !< boltzman constant

  !< control constants for shape of surface.
  integer, parameter :: ID_SHAPE_POLYGON  = 1 !< polygon defined by vertexes
  integer, parameter :: ID_SHAPE_SPHERE   = 2 !< sphere defined by center and radius
  integer, parameter :: ID_SHAPE_CYLINDER = 3 !< cylinder along x-axis defined by center and radius

  !< control constants for transmissivity
  integer, parameter :: ID_TRANS_VALUE = 1 !< transmissivity given by value
  integer, parameter :: ID_TRANS_GAMMA = 2 !< transmissivity given by gamma
  integer, parameter :: ID_TRANS_TBC   = 3 !< transmissivity given by TBC

  !---- verbose for output
  !! verbose >=1, basic outputs like nPhonon, nAngle, nthread, etc
  !! verbose >=2, nGrain, nSurface for each grain, facet indexs of each grain, materials of the grain, etc.
  !! verbose >=3, exception check and print the warning, e.g.out-of-the box or not, inside the correct grain or not.
  !! verbose >=4, debug-mode. Tracking the trajectory of one random particle and visuallize it.  
  integer :: verbose  = 1   !< 0, 1, 2, 3 different level of output information, debug purpose
  integer :: fd_stdout = 6

  !< contral flag for virtual interface
  logical, parameter :: FLAG_ADD_VIRTUAL_INTERFACE = .false. !< use virtual interface or not

  !>
  !! @brief structure for each surface.
  type Surface_type
    logical :: incident      !< incident surface or not.
    logical :: internal      !< internal surface or not.
    integer :: grain_inside  !< index of the grain inside of this surface.
    integer :: grain_outside !< index of the grain outside of this surface.
    integer :: nvertex       !< total number of vertexes belonged to this surface.
    real(8), pointer :: vertex(:,:) => null() !< all vertexes belonged to this surface.
    real(8) :: normal(3)     !< normal vector on this surface.
    real(8) :: tangent(3,3)  !< tangent vectors on this surface.
    real(8) :: radius        !< radius of Sphere or Cylinder.
    real(8) :: center(3)     !< for Sphere and Cylinder
    integer :: kindShape     !< type of shape of this surface, Polygon, Sphere, or Cylinder.
    integer :: kindTrans     !< type of transmissivity, VALUE, Gamma, or TBC.
    integer :: imaterial1    !< material index of this surface
    integer :: imaterial2    !< material index of this surface
    real(8) :: specularity   !< specularity of this surface.
    real(8), pointer :: transmissivity(:,:) => null() !< transsmissivity of this surface for (nbranch,nkpoint).
  end type Surface_type

  !>
  !! @brief structure for each grain.
  type Grain_type
    logical :: incident      !< have incident surface or not.
    real(8) :: generatrix(3) !< position of generatrix point.
    integer :: nsurface      !< number of surfaces belonged to this grain.
    type(Surface_type), pointer :: vsurface(:) => null()
  end type Grain_type
  integer :: ngrain !< total number of grains.
  type(Grain_type), allocatable, target :: vgrain(:) !< array of all grains.
  integer :: igrain_box  !< index of grain for the simulation box
  integer :: igrain_fill !< index of grain for the filling

  !< control constants for mean free path sampling.
  integer, parameter :: ID_MFP_SAMPLE_MATTIESSEN = 1 !< use Mattiessen's rule for sampling
  integer, parameter :: ID_MFP_SAMPLE_MFP        = 2 !< mean free path sampling

  !< control constants for grain type
  integer, parameter :: ID_GRAIN_VORONOI = 1 !< simulation box is filled by voronoi grains
  integer, parameter :: ID_GRAIN_BRICKS  = 2 !< left of bricks is filled by a single filling grain

  !>
  !! @brief structure for input parameters.
  !! SHAO, define a subroutine 'loadParams' to load all params from mcp file
  type Parameter_type
    real(8) :: boxx, boxy, boxz !< size of the simulation box.
    real(8) :: cellspec(6) ! specularity for surfaces on cell boundaries
    real(8) :: specfluc ! fluctuation in specularity on surface.
    real(8) :: temperature !< temperature of the simulation box.
    integer :: nphonon !< total number of sampling phononsl.
    integer :: nangle !< total number of incident polar angle.
    real(8), allocatable :: vangle(:) !< incident polar angle for (nangle).
    real(8) :: scopey, scopez !< range of the scope on the incident boundary surface.
    integer :: kindMFP !< type of mean free path sampling.
    integer :: nMeanFreePath !< total number of mean free pathes.
    real(8), allocatable :: vMeanFreePath(:) !< mean free path for (nMeanFreePath).
    real(8) :: MeanFreePath_current
    integer :: nmaterial !< total number of materials.
    character(len=32), allocatable :: vmaterial(:) !< name of material for (nmaterial).
    integer, allocatable :: kindTrans(:,:)
    real(8), allocatable :: transmissivity(:,:)
    real(8), allocatable :: specularity(:,:)
    integer :: kindGrain !< type of all grains.
  end type Parameter_type
  type(Parameter_type) :: param !< parameters.

  !>
  !! @brief structure for data loaded from Anphon file.
  type Anphon_type
    character(len=32)  :: material !< name of this material.
    real(8) :: Rtemp   !< temperature.
    real(8) :: vol     !< volume.
    integer :: nkpoint !< total number of kpoints recorded in anphon file.
    integer :: nbranch !< total number of branches recorded in anphon file.
    real(8), allocatable :: omega(:,:)  !< omega for (nbranch,nkpoint).
    real(8), allocatable :: tau(:,:)    !< tau for (nbranch,nkpoint).
    real(8), allocatable :: vel(:,:)    !< vel for (nbranch,nkpoint).
    real(8), allocatable :: MFP(:,:)    !< MFP for (nbranch,nkpoint).
    real(8), allocatable :: weight(:,:) !< weight for (nbranch,nkpoint).
    real(8) :: maxOmega !< maximum of omega among all (nbranch,nkpoint).
  end type Anphon_type
  type(Anphon_type) :: anphon !< Anphon data for the first material.

  real(8), allocatable :: transmit(:,:,:) !< transmit for (nbranch,nkpoint,nangle)

contains

  !>
  !! @brief calculates outer product of 3D vectors.
  !! @param [in] va a 3D vector
  !! @param [in] vb a 3D vector
  !! @retval outer product
  function cross_product( va, vb ) result(vc)
    real(8), intent(in) :: va(3)
    real(8), intent(in) :: vb(3)
    real(8) :: vc(3)

    vc(1) = va(2)*vb(3) - va(3)*vb(2)
    vc(2) = va(3)*vb(1) - va(1)*vb(3)
    vc(3) = va(1)*vb(2) - va(2)*vb(1)
  end function cross_product

  logical function isEqual(val, ref)
      real(8), intent(in) :: val
      real(8), intent(in) :: ref
      if (abs(val -ref) < SMALL_VAL) then
          isEqual = .true.
      else 
          isEqual = .false.
      end if
  end function isEqual


  

  !>
  !! @brief read lines from a file until a key found
  !! @param unit unit number of the input file
  !! @param key  a key string to find
  logical function findKey( unit, key )
    integer, intent(in) :: unit
    character(len=*), intent(in) :: key

    integer :: status
    character(len=256) :: line, string

    rewind(unit)
    findKey = .false.
    do
      read(unit,'(a)',iostat=status) line
      if( status /= 0 ) return
      if( len_trim(line)==0 ) cycle
      read(line,*) string
      if( trim(string) == key ) then
        findKey = .true.
        return
      end if
    end do
  end function findKey

  !>
  !! @brief find index of material from the array
  function findMaterial( vmaterial, material ) result(imaterial)
    character(len=32), intent(in) :: vmaterial(:)
    character(len=32), intent(in) :: material
    integer :: imaterial

    do imaterial=1, size(vmaterial)
      if( param%vmaterial(imaterial) == material ) then
        return
      end if
    end do ! imaterial

    imaterial = 0
    return
  end function findMaterial

  !>
  !! @brief load parameters and geometry from MCP file
  subroutine loadMCP( GeoFile, status )
    character(len=*), intent(in) :: GeoFile
    integer, intent(out) :: status

    call loadParameterBox( GeoFile, status )
    if( status /= 0 ) return

    call loadParameterPhonon( GeoFile, status )
    if( status /= 0 ) return

    call loadParameterMaterial( GeoFile, status )
    if( status /= 0 ) return

    call loadParameterTransmission( GeoFile, status )
    if( status /= 0 ) return

    call loadGeometry( GeoFile, status )
    if( status /= 0 ) return

    if( mpi_rank == 0 ) then
      call printSimulationInfo
    end if
  end subroutine loadMCP


  !>
  !! @brief load parameters about simulation box from MCP file
  subroutine loadParameterBox( GeoFile, status )
    character(len=*), intent(in) :: GeoFile
    integer, intent(out) :: status

    integer, parameter :: fd = 11
    character(len=256) :: line
    integer :: NXcell, NYcell, NZcell

    open(fd,file=trim(GeoFile),action='read',iostat=status)
    if( status /= 0 ) then
      write(fd_stdout,*) "Error: input file not found: ", trim(GeoFile)
      return
    end if

    try_block: do ! try block, it is not a loop

      !---- NXcell, number of cells along x-axis
      if( .not. findKey( 11, 'NXcell' ) ) then
        write(fd_stdout,*) "Parameter Error : not found NXcell"
        status = 1
        exit try_block
      else
        read(fd,'(a)') line
        read(line,*,iostat=status) NXcell
        if( status /= 0 ) then
          write(fd_stdout,*) "Parameter Error : broken NXcell"
          exit try_block
        end if
      end if

      !---- NYcell, number of cells along y-axis
      if( .not. findKey( 11, 'NYcell' ) ) then
        write(fd_stdout,*) "Parameter Error : not found NYcell"
        status = 1
        exit try_block
      else
        read(fd,'(a)') line
        read(line,*,iostat=status) NYcell
        if( status /= 0 ) then
          write(fd_stdout,*) "Parameter Error : broken NYcell"
          exit try_block
        end if
      end if

      !---- NZcell, number of cells along z-axis
      if( .not. findKey( 11, 'NZcell' ) ) then
        write(fd_stdout,*) "Parameter Error : not found NZcell"
        status = 1
        exit try_block
      else
        read(fd,'(a)') line
        read(line,*,iostat=status) NZcell
        if( status /= 0 ) then
          write(fd_stdout,*) "Parameter Error : broken NZcell"
          exit try_block
        end if
      end if

      !---- LXcell, length of a cell along x-axis
      if( .not. findKey( 11, 'LXcell' ) ) then
        write(fd_stdout,*) "Parameter Error : not found LXcell"
        status = 1
        exit try_block
      else
        read(fd,'(a)') line
        read(line,*,iostat=status) param%boxx
        if( status /= 0 ) then
          write(fd_stdout,*) "Parameter Error : broken LXcell"
          exit try_block
        end if
        param%boxx = param%boxx * NXcell
      end if

      !---- LYcell, length of a cell along y-axis
      if( .not. findKey( 11, 'LYcell' ) ) then
        write(fd_stdout,*) "Parameter Error : not found LYcell"
        status = 1
        exit try_block
      else
        read(fd,'(a)') line
        read(line,*,iostat=status) param%boxy
        if( status /= 0 ) then
          write(fd_stdout,*) "Parameter Error : broken LYcell"
          exit try_block
        end if
        param%boxy = param%boxy * NYcell
      end if

      !---- LZcell, length of a cell along z-axis
      if( .not. findKey( 11, 'LZcell' ) ) then
        write(fd_stdout,*) "Parameter Error : not found LZcell"
        status = 1
        exit try_block
      else
        read(fd,'(a)') line
        read(line,*,iostat=status) param%boxz
        if( status /= 0 ) then
          write(fd_stdout,*) "Parameter Error : broken LZcell"
          exit try_block
        end if
        param%boxz = param%boxz * NZcell
      end if

      !---- CellSpecularity
      if( .not. findKey( 11, 'CellSpecularity' ) ) then
        param%cellspec(1:6) = 1.0d0 ! default
      else
        read(fd,'(a)') line
        read(line,*,iostat=status) param%cellspec(1:6) ! try read six values
        if( status /= 0 ) then
          read(line,*,iostat=status) param%cellspec(1) ! try read single value
          if( status /= 0 ) then
            write(fd_stdout,*) "Parameter Error : broken CellSpecularity"
            exit try_block
          end if
          param%cellspec(1:6) = param%cellspec(1)
        end if
      end if

      !---- FluctuationSpecularity
      if( .not. findKey( 11, 'FluctuationSpecularity' ) ) then
        param%specfluc = 0.0d0 ! default
      else
        read(fd,'(a)') line
        read(line,*,iostat=status) param%specfluc
        if( status /= 0 ) then
          write(fd_stdout,*) "Parameter Error : broken FluctuationSpecularity"
          exit try_block
        end if
      end if

      exit try_block
    end do try_block
    !---------------------------------------------------------------------------

    close(11)

    ! catch block
    if( status /= 0 ) then
      write(fd_stdout,*) "Error: loadParameterBox"
    end if
  end subroutine loadParameterBox

  !>
  !! @brief load parameters about phonon from MCP file
  subroutine loadParameterPhonon( GeoFile, status )
    character(len=*), intent(in) :: GeoFile
    integer, intent(out) :: status

    integer, parameter :: fd = 11
    character(len=256) :: line, word
    real(8) :: mfp_min, mfp_max ! for MeanFreePath
    integer :: s

    open(fd,file=trim(GeoFile),action='read',iostat=status)
    if( status /= 0 ) then
      write(fd_stdout,*) "Error: input file not found: ", trim(GeoFile)
      return
    end if

    try_block: do ! try block, it is not a loop

      !---- temp, temperature for phonons
      if( .not. findKey( 11, 'temp' ) ) then
        write(fd_stdout,*) "Parameter Error : not found temp"
        status = 1
        exit try_block
      else
        read(fd,'(a)') line
        read(line,*,iostat=status) param%temperature
        if( status /= 0 ) then
          write(fd_stdout,*) "Parameter Error : broken temp"
          exit try_block
        end if
      end if

      !---- Np, number of sampling particles
      if( .not. findKey( 11, 'Np' ) ) then
        write(fd_stdout,*) "Parameter Error : not found Np"
        status = 1
        exit try_block
      else
        read(fd,'(a)') line
        read(line,*,iostat=status) param%nphonon
        if( status /= 0 ) then
          write(fd_stdout,*) "Parameter Error : broken Np"
          exit try_block
        end if
      end if

      !---- Np, number of sampling angles
      if( .not. findKey( 11, 'Nangle' ) ) then
        param%nangle = 91 ! default
      else
        read(fd,'(a)') line
        read(line,*,iostat=status) param%nangle
        if( status /= 0 ) then
          write(fd_stdout,*) "Parameter Error : broken Nangle"
          exit try_block
        end if
      end if
      allocate( param%vangle(param%nangle) )

      !---- ScopeY, range of the scope on the incident boundary surface.
      if( .not. findKey( 11, 'ScopeY' ) ) then
        param%scopey = 1.0d0 ! default
      else
        read(fd,'(a)') line
        read(line,*,iostat=status) param%scopey
        if( status /= 0 ) then
          write(fd_stdout,*) "Parameter Error : broken ScopeY"
          exit try_block
        end if
      end if

      !---- ScopeZ, range of the scope on the incident boundary surface.
      if( .not. findKey( 11, 'ScopeZ' ) ) then
        param%scopez = 1.0d0 ! default
      else
        read(fd,'(a)') line
        read(line,*,iostat=status) param%scopez
        if( status /= 0 ) then
          write(fd_stdout,*) "Parameter Error : broken ScopeZ"
          exit try_block
        end if
      end if

      !---- MeanFreePath
      if( .not. findKey( 11, 'MeanFreePath' ) ) then
        param%kindMFP = ID_MFP_SAMPLE_MATTIESSEN

        ! use default value
        param%nMeanFreePath = 1
        allocate(param%vMeanFreePath(param%nMeanFreePath))
        param%vMeanFreePath(1) = 1.0d10
      else
        param%kindMFP = ID_MFP_SAMPLE_MFP

        read(fd,'(a)') line
        ! try reading min max steps and scale
        read(line,*,iostat=status) mfp_min, mfp_max, param%nMeanFreePath, word
        if( status == 0 ) then
          ! setup param%vMeanFreePath
          if( mfp_min<0.0d0 .or. mfp_min>mfp_max .or. param%nMeanFreePath<=0 ) then
            write(fd_stdout,*) "Parameter Error : broken MeanFreePath"
            status = 1
            exit try_block
          end if
          allocate(param%vMeanFreePath(param%nMeanFreePath))

          select case( word )
          case( "linear" )
            do s=1, param%nMeanFreePath
              param%vMeanFreePath(s) = mfp_min &
                + (mfp_max-mfp_min)/(param%nMeanFreePath-1)*(s-1)
            end do
          case( "log" )
            do s=1, param%nMeanFreePath
              param%vMeanFreePath(s) = exp( &
                log(mfp_min) &
                + (log(mfp_max)-log(mfp_min))/(param%nMeanFreePath-1)*(s-1) )
            end do
          case default
            write(fd_stdout,*) "Parameter Error : broken MeanFreePath"
            status = 1
            exit try_block
          end select

        else
          ! try reading single value
          read(line,*,iostat=status) mfp_min
          if( status /= 0 ) then
            write(fd_stdout,*) "Parameter Error : broken MeanFreePath"
            exit try_block
          end if

          ! setup param%vMeanFreePath
          if( mfp_min<0.0d0 ) then
            write(fd_stdout,*) "Parameter Error : broken MeanFreePath"
            status = 1
            exit try_block
          end if
          param%nMeanFreePath = 1
          allocate(param%vMeanFreePath(param%nMeanFreePath))
          param%vMeanFreePath(1) = mfp_min
        end if
      end if

      !---- particle_trace
      if( .not. findKey( 11, 'Trace' ) ) then
        particle_trace = 1000 ! default
      else
        read(fd,'(a)') line
        read(line,*,iostat=status) particle_trace
        if( status /= 0 ) then
          write(fd_stdout,*) "Parameter Error : Trace"
          exit try_block
        end if
      end if

      exit try_block
    end do try_block
    !---------------------------------------------------------------------------

    close(11)

    ! catch block
    if( status /= 0 ) then
      write(fd_stdout,*) "Error: loadParameterPhonon"
    end if
  end subroutine loadParameterPhonon

  !>
  !! @brief load parameters about material from MCP file
  subroutine loadParameterMaterial( GeoFile, status )
    character(len=*), intent(in) :: GeoFile
    integer, intent(out) :: status

    integer, parameter :: fd = 11
    character(len=256) :: line

    integer :: im

    open(fd,file=trim(GeoFile),action='read',iostat=status)
    if( status /= 0 ) then
      write(fd_stdout,*) "Error: input file not found: ", trim(GeoFile)
      return
    end if

    try_block: do ! try block, it is not a loop

      !---- Nm, number of materials
      if( .not. findKey( 11, 'Nm' ) ) then
        write(fd_stdout,*) "Parameter Error : not found Nm"
        status = 1
        exit try_block
      else
        read(fd,'(a)') line
        read(line,*,iostat=status) param%nmaterial
        if( status /= 0 ) then
          write(fd_stdout,*) "Parameter Error : broken Nm"
          exit try_block
        end if
        allocate( param%vmaterial(param%nmaterial) )
        read(fd,*) ! skip "material contents"
        do im=1,param%nmaterial
          read(fd,'(a)') line
          read(line,*,iostat=status) param%vmaterial(im)
          if( status /= 0 ) then
            write(fd_stdout,*) "Parameter Error : broken Nm"
            exit try_block
          end if
          read(fd,*) ! skip a value
        end do ! im
      end if

      ! setup of anphon data
      !!SHAO: nmaterials will always equal 1 in this code, remove the else option
      if( param%nmaterial==1 ) then
        call loadAnphon( anphon, param%vmaterial(1), status )
        if( status /= 0 ) then
          write(fd_stdout,*) "Error: failed to load anphon."
          exit try_block
        end if
      else
        anphon%nkpoint = 1
        anphon%nbranch = 1
      end if

      exit try_block
    end do try_block
    !---------------------------------------------------------------------------
    close(11)

    ! catch block
    if( status /= 0 ) then
      write(fd_stdout,*) "Error: loadParameterMaterial"
    end if
  end subroutine loadParameterMaterial

  !>
  !! @brief load anphon data of material
  subroutine loadAnphon( anphon, material, status )
    type(Anphon_type), intent(out) :: anphon
    character(len=32), intent(in)  :: material
    integer, intent(out) :: status

    character(len=256) :: filename, word
    integer, parameter :: fd = 100
    integer :: ik, ib

    anphon%material = material
    filename = 'Materials/'//trim(anphon%material)//'.anphon'

    open( fd, file=trim(filename), action='read', iostat=status )
    if( status /= 0 ) then
      write(fd_stdout,*) "Error: anphon file not found: ", trim(filename)
      return
    end if

    read(fd,*)
    read(fd,*)
    read(fd,*,iostat=status) word, word, word, word, word, anphon%RTemp
    if( status /= 0 ) then
      write(fd_stdout,*) "Anphon Error : broken RTemp"
      close(fd); return
    end if

    read(fd,*,iostat=status) word, word, word, anphon%vol
    if( status /= 0 ) then
      write(fd_stdout,*) "Anphon Error : broken vol"
      close(fd); return
    end if

    read(fd,*,iostat=status) word, word, word, word, anphon%nkpoint
    if( status /= 0 ) then
      write(fd_stdout,*) "Anphon Error : broken nkpoint"
      close(fd); return
    end if

    read(fd,*,iostat=status) word, word, word, word, anphon%nbranch
    if( status /= 0 ) then
      write(fd_stdout,*) "Anphon Error : broken nBranch"
      close(fd); return
    end if

    read(fd,*)

    allocate( anphon%omega(anphon%nbranch,anphon%nkpoint) )
    allocate( anphon%tau(anphon%nbranch,anphon%nkpoint) )
    allocate( anphon%vel(anphon%nbranch,anphon%nkpoint) )
    allocate( anphon%MFP(anphon%nbranch,anphon%nkpoint) )
    allocate( anphon%weight(anphon%nbranch,anphon%nkpoint) )

    do ik=1,anphon%nkpoint
      do ib=1,anphon%nbranch
        read(fd,*,iostat=status) word, word, &
          anphon%omega(ib,ik), &
          anphon%tau(ib,ik), &
          anphon%vel(ib,ik), &
          anphon%MFP(ib,ik), &
          anphon%weight(ib,ik)
        if( status /= 0 ) then
          write(fd_stdout,*) "Anphon Error : broken omega"
          close(fd); return
        end if
      end do ! ll
    end do ! kk
    anphon%omega(:,:) = anphon%omega(:,:)*0.03d0*M_TWOPI*1.0d12
    anphon%maxOmega = maxval(anphon%omega(:,:))

    close(fd)

    status = 0
  end subroutine loadAnphon

  !>
  !! @brief free memory for anphon data
  subroutine clearAnphon( anphon )
    type(Anphon_type), intent(out) :: anphon
    
    deallocate( anphon%omega )
    deallocate( anphon%tau )
    deallocate( anphon%vel )
    deallocate( anphon%MFP )
    deallocate( anphon%weight )
  end subroutine clearAnphon

  !>
  !! @brief load parameters about transmission from MCP file
  subroutine loadParameterTransmission( GeoFile, status )
    character(len=*), intent(in) :: GeoFile
    integer, intent(out) :: status

    integer, parameter :: fd = 11
    character(len=256) :: line, word
    character(len=32) :: material1, material2
    integer :: imaterial1, imaterial2
    real(8) :: gamma, TBC
    integer :: Nt1, Nt2, Nt, it
    character(len= 5) :: kindTrans

    open(fd,file=trim(GeoFile),action='read',iostat=status)
    if( status /= 0 ) then
      write(fd_stdout,*) "Error: input file not found: ", trim(GeoFile)
      return
    end if

    try_block: do ! try block, it is not a loop

      !---- Nt, number of truns functions
      if( .not. findKey( 11, 'Nt' ) ) then
        write(fd_stdout,*) "Parameter Error : not found Nt"
        status = 1
        exit try_block
      else
        read(fd,'(a)') line
        read(line,*,iostat=status) Nt1, Nt2
        if( status /= 0 ) then
          write(fd_stdout,*) "Parameter Error : broken Nt"
          exit try_block
        end if
      end if

      !! SHAO: the following part is to construct the transmission/specularity tables betweern material ij
      !! Better to wrapping it into a subroutine/code-block
      !! Nt: trans function table
      ! setup of param%transmissivity
      Nt = Nt1*Nt2
      read(fd,*) ! skip "materials"
      allocate( param%transmissivity(Nt2+1,Nt1) )
      allocate( param%specularity(Nt2+1,Nt1) )
      allocate( param%kindTrans(Nt2+1,Nt1) )

      ! Set Boundary conditions
      param%transmissivity(Nt2+1,:) = 0.0d0
      param%specularity(Nt2+1,:) = 0.0d0
      param%kindTrans(Nt2+1,:) = ID_TRANS_VALUE

      !! by SHAO, subroutine, `read_transmission`
      !! read the transmissivity and specularity matrix.
      !! kindTrans records how the transmissivity is given: "Value", "TBC", et, et, etc
      !! The global one is saved in specularity, Mode, and  Transmissivity(:,:,:)

      do it=1,Nt
        read(fd,'(a)') line
        read(line,*,iostat=status) material1, material2, word, kindTrans
        if( status /= 0 ) then
          write(fd_stdout,*) "Parameter Error : broken Specular"
          exit try_block
        end if

        imaterial1 = findMaterial( param%vmaterial(1:Nt1), material1 )
        if( imaterial1 == 0 ) then
          write(fd_stdout,*) "Parameter Error : unknown material in Specular:", trim(material1)
          status = 1
          exit try_block
        end if

        imaterial2 = findMaterial( param%vmaterial(1:Nt2), material2 )
        if( imaterial2 == 0 ) then
          write(fd_stdout,*) "Parameter Error : unknown material in Specular:", trim(material2)
          status = 1
          exit try_block
        end if

        select case( kindTrans )
        case("Value")
          param%kindTrans(imaterial2,imaterial1) = ID_TRANS_VALUE
        case("Gamma")
          param%kindTrans(imaterial2,imaterial1) = ID_TRANS_GAMMA
        case("TBC")
          param%kindTrans(imaterial2,imaterial1) = ID_TRANS_TBC
        case default
          write(fd_stdout,*) "Error in reading file : Transmission mode is Mismatched"
          status = 1
          exit try_block
        end select

        select case( param%kindTrans(imaterial2,imaterial1) )
        case(ID_TRANS_VALUE)
          read(fd,'(a)') line
          read(line,*,iostat=status) param%transmissivity(imaterial2,imaterial1), param%specularity(imaterial2,imaterial1)
          if( status /= 0 ) then
            write(fd_stdout,*) "Parameter Error : broken Specular Value"
            status = 1
            exit try_block
          end if

        case(ID_TRANS_GAMMA)
          read(fd,'(a)') line
          read(line,*,iostat=status) param%transmissivity(imaterial2,imaterial1), param%specularity(imaterial2,imaterial1)
          if( status /= 0 ) then
            write(fd_stdout,*) "Parameter Error : broken Specular Gamma"
            exit try_block
          end if

        case(ID_TRANS_TBC)
          read(fd,'(a)') line
          read(line,*,iostat=status) TBC, param%specularity(imaterial2,imaterial1)
          if( status /= 0 ) then
            write(fd_stdout,*) "Parameter Error : broken Specular TBC"
            exit try_block
          end if

          !SHAO: giving the TBC, and calculate a effective **gamma**
          !SHAO: the wavevector-dependent transmission is then calculated from gamma
          TBC = TBC * 1.0d6
          call GammaDeterm( material1, TBC, gamma, 1 )
          param%transmissivity(imaterial2,imaterial1) = gamma
          if( gamma<0 ) then
            write(fd_stdout,*) "Parameter Error : Gamma cannot determine"
            status = 1
            exit try_block
          end if
        end select
      end do ! it
      !! by SHAO, end-of-subroutine, `read_transmission`

      allocate( transmit(anphon%nbranch,anphon%nkpoint,param%nangle) )

      exit try_block
    end do try_block
    !---------------------------------------------------------------------------
    close(11)

    ! catch block
    if( status /= 0 ) then
      write(fd_stdout,*) "Error: loadParameterTransmission"
      return
    end if
  end subroutine loadParameterTransmission

  !>
  !! @brief load geometry from MCP file
  subroutine loadGeometry( GeoFile, status )
    character(len=*), intent(in) :: GeoFile
    integer, intent(out) :: status

    integer, parameter :: fd = 11
    character(len=256) :: line, word
    integer :: igrain, isurface
    type(Surface_type), pointer :: surface
    integer :: kindShape !< type of shape of each surface.
    real(8) :: diameter ! for Sphere and Cylinder
    character(len=32) :: material1, material2
    integer :: imaterial1
    integer :: kk
    type(Grain_type), pointer :: grain

    open(fd,file=trim(GeoFile),action='read',iostat=status)
    if( status /= 0 ) then
      write(fd_stdout,*) "Error: input file not found: ", trim(GeoFile)
      return
    end if

    try_block: do ! try block, it is not a loop

      !---- Shape, shape of grain
      if( .not. findKey( 11, 'Shape' ) ) then
        write(fd_stdout,*) "Parameter Error : not found Shape"
        status = 1
        exit try_block
      else
        read(fd,'(a)') line
        read(line,*,iostat=status) word
        if( status /= 0 ) then
          write(fd_stdout,*) "Parameter Error : broken Shape"
          exit try_block
        end if

        select case( word )
        case("voronoi")
          param%kindGrain = ID_GRAIN_VORONOI
        case("bricks")
          param%kindGrain = ID_GRAIN_BRICKS
        case default
          write(fd_stdout,*) "Parameter Error : unknown grain Shape: ", trim(word)
          status = 1
          exit try_block
        end select
      end if

      !---- Ng, number of grains
      if( .not. findKey( 11, 'Ng' ) ) then
        write(fd_stdout,*) "Parameter Error : not found Ng"
        status = 1
        exit try_block
      else
        read(fd,'(a)') line
        read(line,*,iostat=status) ngrain
        if( status /= 0 ) then
          write(fd_stdout,*) "Parameter Error : broken Ng"
          exit try_block
        end if
      end if

      ! Coordinate of grains
      select case( param%kindGrain )
      case(ID_GRAIN_VORONOI)
        allocate( vgrain(0:ngrain) )
        igrain_box  = 0
        igrain_fill = -1
      case(ID_GRAIN_BRICKS)
        allocate( vgrain(0:ngrain+1) )
        igrain_box  = 0
        igrain_fill = ngrain+1
      end select

      call setupBoxGrain( vgrain(igrain_box) )

      do igrain=1, ngrain
        grain => vgrain(igrain)
        grain%incident = .false.

        read(fd,*) ! skip "## grain No.*"
        read(fd,*) ! skip "generator"
        read(fd,'(a)') line
        read(line,*,iostat=status) grain%generatrix(:) ! xyz
        if( status /= 0 ) then
          write(fd_stdout,*) "Parameter Error : broken generatrix"
          exit try_block
        end if

        read(fd,*,iostat=status) word !   box, sphere, cylinder, or polyhedron
        if( status /= 0 ) then
          write(fd_stdout,*) "Parameter Error : broken grain shape"
          exit try_block
        else
          if( param%kindGrain == ID_GRAIN_VORONOI ) then
            read(fd,*) ! skip lengthes
            kindShape = ID_SHAPE_POLYGON
          else
            select case( word )
            case("box","length")
              read(fd,*) ! skip lengthes
              kindShape = ID_SHAPE_POLYGON
              diameter = 0.0d0

            case("sphere")
              read(fd,'(a)') line
              read(line,*,iostat=status) diameter
              if( status /= 0 ) then
                write(fd_stdout,*) "Parameter Error : broken diameter"
                exit try_block
              end if
              kindShape = ID_SHAPE_SPHERE

            case("cylinder")
              read(fd,'(a)') line
              read(line,*,iostat=status) diameter
              if( status /= 0 ) then
                write(fd_stdout,*) "Parameter Error : broken diameter"
                exit try_block
              end if
              kindShape = ID_SHAPE_CYLINDER

            case("polyhedron")
              read(fd,*) ! skip lengthes
              write(fd_stdout,*) "load grain from STL file"
              kindShape = ID_SHAPE_POLYGON

            case default
              write(fd_stdout,*) "Parameter Error : unknown grain shape: ", trim(word)
              status = 1
              exit try_block
            end select
          end if
        end if
        !---- grain material
        read(fd,*) ! skip "material"
        read(fd,'(a)') line
        read(line,*,iostat=status) material1
        if( status /= 0 ) then
          write(fd_stdout,*) "Parameter Error : broken material"
          exit try_block
        else
          imaterial1 = findMaterial( param%vmaterial, material1 )
          if( imaterial1 == 0 ) then
            write(fd_stdout,*) "Parameter Error : unknown material:", trim(material1)
            status = 1
            exit try_block
          end if
        end if
        !---- load nsurface of grain
        if( kindShape == ID_SHAPE_POLYGON ) then
          read(fd,*) ! skip "# facets of grain No.*"
          read(fd,*) ! skip "Nf # number of facets"
          read(fd,'(a)') line
          read(line,*,iostat=status) grain%nsurface
          if( status /= 0 ) then
            write(fd_stdout,*) "Parameter Error : broken Nf"
            exit try_block
          end if
        else
          grain%nsurface = 1 ! sphere or cylinder
        end if

        allocate( grain%vsurface(grain%nsurface) )
        !---- load surfaces of current grain
        do isurface=1,grain%nsurface
          surface => grain%vsurface(isurface)
          surface%incident = .false.
          surface%internal = .true.
          surface%grain_inside = igrain ! record index of grain of this surface
          surface%imaterial1 = imaterial1
          surface%kindShape = kindShape
          surface%center(:) = grain%generatrix(:)
          surface%radius = diameter*0.5d0

          !---- load vertexs of current surface
          if( surface%kindShape == ID_SHAPE_POLYGON ) then

            ! facet
            read(fd,*) ! skip "#  facet No.*"
            read(fd,*) ! skip "Nv # number of vertex"
            read(fd,'(a)') line
            read(line,*,iostat=status) surface%nvertex
            if( status /= 0 ) then
              write(fd_stdout,*) "Parameter Error : broken Nv"
              exit try_block
            else  !---- load surface vertex
              allocate( surface%vertex(3,surface%nvertex+1) )
              do kk=1,surface%nvertex
                read(fd,*) ! skip "vertex *"
                read(fd,'(a)') line
                read(line,*,iostat=status) surface%vertex(:,kk)
                if( status /= 0 ) then
                  write(fd_stdout,*) "Parameter Error : broken vertex"
                  exit try_block
                end if
              end do ! kk

              ! add the begin vertex as the end vertex
              surface%vertex(:,surface%nvertex+1) = surface%vertex(:,1)

              call setupSurfaceGeometry( surface )
            end if
            
            !! write(*,*) "surface%normal = ", surface%normal
            !! write(*,*) "surface%center = ", surface%center
            !! if( .not. isInsideSurface(surface, surface%center)) then
            !!     write(*,*) "surface center is insider the current surface"
            !! end if

            if( isIncidentSurface( surface ) ) then
              grain%incident   = .true.
              surface%incident = .true.
              surface%internal = .false.
            end if

            if( .not. isInternalSurface( surface ) ) then
              surface%internal = .false.
            end if
          else
            ! sphere, cylinder
            surface%nvertex = 0
          end if
          !---- "name" of material_to
          read(fd,'(a)') line
          read(line,*,iostat=status) word, word, material2 ! material_to
          if( status /= 0 ) then
            write(fd_stdout,*) "Parameter Error : broken matrial_to"
            exit try_block
          end if
          !---- grain_to, index
          ! 0 means the outermost box
          read(fd,'(a)') line
          read(line,*,iostat=status) word, word, surface%grain_outside ! grain_to
          if( status /= 0 ) then
            write(fd_stdout,*) "Parameter Error : broken grain_to"
            exit try_block
          end if

          if( .not. surface%internal ) cycle ! do not set the following properties for this surface.
          
          !---- for "brick" grain, all material_to == "Internal", means connected to the fill_grain
          if( trim(material2)=='Internal' ) then
            surface%imaterial2 = 1 !! Material of the fill_grain
            surface%grain_outside = igrain_fill  !! SHAO: Last grain, the simulation box
          else
            surface%imaterial2 = findMaterial( param%vmaterial, material2 )
            if( surface%imaterial2 == 0 ) then
              write(fd_stdout,*) "Parameter Error : unknown material:", trim(material2)
              status = 1
              exit try_block
            end if
          end if
          
          !! Mode: transmission type of each surface,  
          !! can be "Gamma", "TBC", or "Value"
          surface%specularity = param%specularity(surface%imaterial2,surface%imaterial1)
          if( param%specfluc>0.0d0 ) then
            surface%specularity = surface%specularity + normal_distribution(0.0d0,param%specfluc)
            if( surface%specularity < 0.0d0 ) surface%specularity = 0.0d0
            if( surface%specularity > 1.0d0 ) surface%specularity = 1.0d0
            ! SHAO, check surface specularity
            write(fd_stdout,*) "surface speculartiy", surface%specularity
          end if
          surface%kindTrans   = param%kindTrans(surface%imaterial2,surface%imaterial1)

          allocate( surface%transmissivity(anphon%nbranch,anphon%nkpoint) )
          call setupSurfaceTransmissivity( surface, surface%kindTrans, &
            param%transmissivity(surface%imaterial2,surface%imaterial1) )

        end do ! isurface
        !! if(isInsideGrain(grain, grain%generatrix)) then
        !!     write(*,*) "generatrix is inside the currecnt grain"
        !! end if

      end do ! igrain
      !! SHAO: end-of-subroutine `count_surface`

      !! SHAO: copying all the surfaces into the **last** GRAIN
      if( param%kindGrain == ID_GRAIN_BRICKS ) then
        call setupFillGrain( vgrain(igrain_fill) )
      end if ! kindGrain
      !! end-of-subroutine `construct_simulation_box`
      exit try_block
    end do try_block
    !---------------------------------------------------------------------------
    close(11)
    ! catch block
    if( status /= 0 ) then
      write(fd_stdout,*) "Error: loadGeometry"
      return
    end if

  contains
    !>
    !! @brief calculate a random number that obeys normal distribution
    real(8) function normal_distribution( mean, deviation )
      real(8), intent(in) :: mean, deviation
      real(8) :: sum, rand
      integer :: i

      sum = 0.0d0
      do i=1, 12
        call random_number(rand)
        sum = sum + rand
      end do
      normal_distribution = (sum - 6.0d0)*deviation + mean
    end function normal_distribution
    

    !>
    !! @brief check this surface locates outside of the z- boundary of simulation box
    logical function isIncidentSurface( surface ) 
      type(Surface_type), intent(in) :: surface

      isIncidentSurface = all(surface%vertex(1,:)<=0.0d0)  ! <== Fixed a bug here
    end function isIncidentSurface
    

    !>
    !! @brief check this surface locates inside of the simulation box
    logical function isInternalSurface( surface )
      type(Surface_type), intent(in) :: surface

      if( all(surface%vertex(1,:)<=0.0d0) ) then
        isInternalSurface = .false.
        return
      end if

      if( all(surface%vertex(1,:)>=param%boxx) ) then
        isInternalSurface = .false.
        return
      end if

      if( all(surface%vertex(2,:)<=0.0d0) ) then
        isInternalSurface = .false.
        return
      end if

      if( all(surface%vertex(2,:)>=param%boxy) ) then
        isInternalSurface = .false.
        return
      end if

      if( all(surface%vertex(3,:)<=0.0d0       ) ) then
        isInternalSurface = .false.
        return
      end if

      if( all(surface%vertex(3,:)>=param%boxz  ) ) then
        isInternalSurface = .false.
        return
      end if

      isInternalSurface = .true.
      return
    end function isInternalSurface
  end subroutine loadGeometry

  !>
  !! @brief setup surfaces of the simulation box grain
  subroutine setupBoxGrain( grain )
    integer :: isurface
    type(Grain_type), intent(inout) :: grain
    type(Surface_type), pointer :: surface

    grain%incident      = .false.
    grain%nsurface      = 6
    allocate( grain%vsurface(grain%nsurface) )
    grain%generatrix(:) = [ 0.0d0, 0.0d0, 0.0d0 ]

    ! setup surfaces that form simulation box.
    do isurface=1, 6
      surface => grain%vsurface(isurface)

      surface%incident = .false.
      surface%internal = .true.
      surface%grain_inside = 0 ! index of point of surface
      surface%kindShape = ID_SHAPE_POLYGON
      !surface%center(:) = [ 0.0d0, 0.0d0, 0.0d0 ]
      surface%center(:) = [ 0.5*param%boxx, 0.5*param%boxy, 0.5*param%boxz ]
      surface%radius    = 0.0d0
      surface%grain_outside = 0

      surface%nvertex = 4
      allocate( surface%vertex(3,surface%nvertex+1) )
      select case( isurface )
      case(1) ! x-
        surface%vertex(:,1) = [ 0.0d0,      0.0d0,      0.0d0      ]
        surface%vertex(:,2) = [ 0.0d0,      param%boxy, 0.0d0      ]
        surface%vertex(:,3) = [ 0.0d0,      param%boxy, param%boxz ]
        surface%vertex(:,4) = [ 0.0d0,      0.0d0,      param%boxz ]
        
        surface%vertex(:,1) = [ 0.0d0,      0.0d0,      0.0d0      ]
        surface%vertex(:,2) = [ 0.0d0,      0.0d0,      param%boxz ]
        surface%vertex(:,3) = [ 0.0d0,      param%boxy, param%boxz ]
        surface%vertex(:,4) = [ 0.0d0,      param%boxy, 0.0d0      ]
      case(2) ! x+
        surface%vertex(:,1) = [ param%boxx, 0.0d0,      0.0d0      ]
        surface%vertex(:,2) = [ param%boxx, 0.0d0,      param%boxz ]
        surface%vertex(:,3) = [ param%boxx, param%boxy, param%boxz ]
        surface%vertex(:,4) = [ param%boxx, param%boxy, 0.0d0      ]
        
        surface%vertex(:,1) = [ param%boxx, 0.0d0,      0.0d0      ]
        surface%vertex(:,4) = [ param%boxx, 0.0d0,      param%boxz ]
        surface%vertex(:,3) = [ param%boxx, param%boxy, param%boxz ]
        surface%vertex(:,2) = [ param%boxx, param%boxy, 0.0d0      ]
      case(3) ! y-
        surface%vertex(:,1) = [ 0.0d0,      0.0d0,      0.0d0      ]
        surface%vertex(:,4) = [ 0.0d0,      0.0d0,      param%boxz ]
        surface%vertex(:,3) = [ param%boxx, 0.0d0,      param%boxz ]
        surface%vertex(:,2) = [ param%boxx, 0.0d0,      0.0d0      ]
      case(4) ! y+
        surface%vertex(:,1) = [ param%boxx, param%boxy, 0.0d0      ]
        surface%vertex(:,4) = [ param%boxx, param%boxy, param%boxz ]
        surface%vertex(:,3) = [ 0.0d0,      param%boxy, param%boxz ]
        surface%vertex(:,2) = [ 0.0d0,      param%boxy, 0.0d0      ]
      case(5) ! z-
        surface%vertex(:,1) = [ 0.0d0,      0.0d0,      0.0d0      ]
        surface%vertex(:,4) = [ param%boxx, 0.0d0,      0.0d0      ]
        surface%vertex(:,3) = [ param%boxx, param%boxy, 0.0d0      ]
        surface%vertex(:,2) = [ 0.0d0,      param%boxy, 0.0d0      ]
      case(6) ! z+
        surface%vertex(:,1) = [ 0.0d0,      param%boxy, param%boxz ]
        surface%vertex(:,4) = [ param%boxx, param%boxy, param%boxz ]
        surface%vertex(:,3) = [ param%boxx, 0.0d0,      param%boxz ]
        surface%vertex(:,2) = [ 0.0d0,      0.0d0,      param%boxz ]
      end select

      surface%vertex(:,5) = surface%vertex(:,1)

      call setupSurfaceGeometry( surface )

      surface%specularity = param%cellspec(isurface)
      surface%kindTrans   = ID_TRANS_VALUE
      allocate( surface%transmissivity(anphon%nbranch,anphon%nkpoint) )
      call setupSurfaceTransmissivity( surface, surface%kindTrans, 0.0d0 )
    end do ! isurface
  end subroutine setupBoxGrain

  !>
  !! @brief setup surfaces of the filling grain
  subroutine setupFillGrain( grain_fill )
    type(Grain_type), intent(inout) :: grain_fill

    type(Grain_type), pointer :: grain
    type(Surface_type), pointer :: surfacej
    type(Surface_type), pointer :: surface
    integer :: igrain, jsurf, isurface, ivertex

    ! count up all surfaces
    grain_fill%nsurface = 0
    do igrain=1, ngrain
      grain => vgrain(igrain)
      do jsurf=1, grain%nsurface
        surfacej => grain%vsurface(jsurf)
        if( .not. surfacej%internal ) cycle

        grain_fill%nsurface = grain_fill%nsurface + 1
      end do
    end do

    allocate( grain_fill%vsurface(grain_fill%nsurface) )

    isurface = 1
    do igrain=1, ngrain
      grain => vgrain(igrain)
      do jsurf=1, grain%nsurface
        surfacej => grain%vsurface(jsurf)
        if( .not. surfacej%internal ) cycle  !! Skips surfaces overlap with the simulation box

        surface  => grain_fill%vsurface(isurface)
        surface%incident = .false.
        surface%internal = .true.
        surface%grain_inside = igrain_fill
        surface%imaterial1 = 1
        surface%imaterial2 = surfacej%imaterial1
        surface%kindShape = surfacej%kindShape
        surface%center(:) = surfacej%center(:)  ! TODO, why we need surface%center?
        surface%radius = surfacej%radius

        surface%nvertex = surfacej%nvertex
        if( surface%nvertex>0 ) then
          allocate( surface%vertex(3,surface%nvertex+1) )
          do ivertex=1,surface%nvertex+1  !<== each surface has nvertex+1 point
              surface%vertex(:,ivertex) = surfacej%vertex(:,surface%nvertex+2-ivertex)
          end do
          !isurface%vertex(:,:) = surfacej%vertex(:,:)
        end if
        write(*,*) "grain=", igrain, "surface =", jsurf, "nvertex =", surface%nvertex
        call setupSurfaceGeometry( surface )

        surface%grain_outside = surfacej%grain_inside

        surface%specularity = param%specularity(surface%imaterial2,surface%imaterial1)
        surface%kindTrans   = param%kindTrans(surface%imaterial2,surface%imaterial1)

        allocate( surface%transmissivity(anphon%nbranch,anphon%nkpoint) )
        call setupSurfaceTransmissivity( surface, surface%kindTrans, &
          param%transmissivity(surface%imaterial2,surface%imaterial1) )

        isurface = isurface + 1
      end do ! jsurf
    end do ! igrain
  end subroutine setupFillGrain

  !>
  !! @brief setup surface transmissivity
  subroutine setupSurfaceTransmissivity( surface, kindTrans, gamma )
    type(Surface_type), intent(inout) :: surface
    integer, intent(in) :: kindTrans
    real(8), intent(in) :: gamma

    integer :: ik, ib

    select case( kindTrans )
      !! SHAO: "gray" mode for phonon transmission, wavevector and frequency independent
    case( ID_TRANS_VALUE )
      surface%transmissivity(:,:) = gamma

      !! SHAO: frequency-dependent transmission mode
      !! SHAO: param%transmissivity eq. real transmission or the gamma value
      !! SHAO: for TBC or GAMMA, transmission = 1/(gamma*omega + 1)
    case( ID_TRANS_GAMMA, ID_TRANS_TBC )
      if( param%nmaterial==1 ) then
        do ik=1,anphon%nkpoint
          do ib=1,anphon%nbranch
            surface%transmissivity(ib,ik) &
              = 1.0d0/(gamma*anphon%omega(ib,ik)/anphon%maxOmega + 1.0d0)
          end do ! ib
        end do ! ik
      else
        surface%transmissivity(:,:) = 1.0d0/(gamma+1.0d0)
      end if
    end select
  end subroutine setupSurfaceTransmissivity

  !>
  !! @brief setup normal and tangent vectors of a surface
  subroutine setupSurfaceGeometry( surface )
    type(Surface_type), intent(inout) :: surface
    integer :: iv

    real(8) :: crossab(3), vectora(3), vectorb(3)

    select case( surface%kindShape )
    case( ID_SHAPE_POLYGON )
      ! facet
      surface%center(:) = 0.0d0
      do iv=1,surface%nvertex
        surface%center(:) = surface%center(:) + surface%vertex(:,iv) 
      end do
      surface%center(:) = 1.0d0/surface%nvertex*surface%center(:)

      surface%normal(:) = 0.0d0
      do iv=1,surface%nvertex-1 
        vectora = surface%vertex(:,iv+1) - surface%vertex(:,iv)
        vectorb = surface%vertex(:,iv+2) - surface%vertex(:,iv)
        crossab(:) = cross_product( vectora, vectorb )
        surface%normal(:) = surface%normal(:) + crossab(:)/sqrt(dot_product(crossab,crossab))
      end do
      vectora = surface%vertex(:,surface%nvertex+1) - surface%vertex(:,surface%nvertex)
      vectorb = surface%vertex(:,2) - surface%vertex(:,surface%nvertex)
      crossab(:) = cross_product( vectora, vectorb )
      surface%normal(:) = surface%normal(:) + crossab(:)/sqrt(dot_product(crossab,crossab))
      surface%normal(:) = 1.0d0/surface%nvertex*surface%normal(:)

      !vectora = surface%vertex(:,2) - surface%vertex(:,1)
      !vectorb = surface%vertex(:,surface%nvertex) - surface%vertex(:,1)
      !crossab(:) = cross_product( vectora, vectorb )
      !surface%normal(:) = crossab(:)/sqrt(dot_product(crossab,crossab))

      surface%tangent(:,3) = surface%normal(:)

      if( abs(surface%normal(3))==1.0d0 ) then
        surface%tangent(:,2) = [ 0.0d0, 1.0d0, 0.0d0 ] 
        surface%tangent(:,1) = [ 1.0d0, 0.0d0, 0.0d0 ]
      else if( surface%normal(1)**2==0.0d0 .and. &
        surface%normal(2)**2==0.0d0 ) then
        surface%tangent(:,2) = [ 0.0d0, 1.0d0, 0.0d0 ]
        surface%tangent(:,1) = [ 1.0d0, 0.0d0, 0.0d0 ]
      else
        surface%tangent(:,2) &
          = [ &
          -surface%normal(2) &
          / sqrt(surface%normal(1)**2 + surface%normal(2)**2), &
          surface%normal(1) &
          / sqrt(surface%normal(1)**2 + surface%normal(2)**2), &
          0.0d0 ]

        surface%tangent(:,1) &
          = cross_product( surface%tangent(:,2), surface%tangent(:,3) )
      end if

    case( ID_SHAPE_SPHERE )
    case( ID_SHAPE_CYLINDER )
    end select
  end subroutine setupSurfaceGeometry

  !>
  !! @brief print informations
  subroutine printSimulationInfo
    use omp_lib

    write(fd_stdout,*) ""
    write(fd_stdout,*) "*** Simulation infomation ***"

    select case(param%kindMFP)
    case(ID_MFP_SAMPLE_MATTIESSEN)
      write(fd_stdout,*) "- Enabled Mattiessen's rule"
      write(fd_stdout,*) "- MeanFreePath ", real(param%vMeanFreePath(:))
    case(ID_MFP_SAMPLE_MFP)
      write(fd_stdout,*) "- Enabled MFP sampling"
      write(fd_stdout,*) "- vMeanFreePath ", real(param%vMeanFreePath(:))
    end select

    !write(fd_stdout,'(a, i4)') " No. of  MFP to sample:       ", param%nMeanFreePath
    write(fd_stdout,*) "- Sampling particle: ", param%nphonon
    write(fd_stdout,*) "- Sampling angle: ", param%nangle
    write(fd_stdout,*) "- Max threads: ", omp_get_max_threads()
    write(fd_stdout,*) "- Incident particle:", param%nphonon

    write(fd_stdout,*) ""
    write(fd_stdout,*) "*** Geometry information ***"
    select case( param%kindGrain )
    case(ID_GRAIN_VORONOI)
      write(fd_stdout,*) "- Share of grains: voronoi"
    case(ID_GRAIN_BRICKS)
      write(fd_stdout,*) "- Share of grains: bricks"
    end select
    !!$      write(fd_stdout,*) "- Facet Num: ", nfacet
    !!$      write(fd_stdout,*) "- nsurface: ", nsurface
    !!$      write(fd_stdout,*) "- Valid surface inside domain: ", nvalidsurface
    !!$      write(fd_stdout,*) "- ngrain: ", ngrain

    if (verbose >= 1) then
      write(fd_stdout,*) "- Num of grain: ", ngrain
      write(fd_stdout,*) '- Indexs for grain surfaces:'
    !!$        write(fd_stdout,*) "- Surface: grain inside: ", vsurface(:)%grain_inside
    !!$        write(fd_stdout,*) "- Surface: grain outside: ", vsurface(:)%grain_outside
      write(fd_stdout,*) "- Facets in grain (mcp):", vgrain(0:)%nsurface
      write(*,*) vgrain(0)%vsurface(0)%center(:)
      write(*,*) vgrain(0)%vsurface(1)%center(:)
      write(*,*) vgrain(0)%vsurface(2)%center(:)
      write(*,*) vgrain(0)%vsurface(3)%center(:)
    end if ! verbose >= 1
  end subroutine printSimulationInfo

  !>
  !! @brief calculate all transmission functions
  subroutine TransmissionFunction
    use omp_lib

    integer :: ik, ib, ia
    real(8) :: trans_mean

#ifdef MPI
    real(8), allocatable :: transmit_mpi(:,:,:)
#endif

    do ia=1,param%nangle
      param%vangle(ia) = 0.0d0 + dble(ia-1)/dble(param%nangle-1)* M_PI_2
    end do ! angle

    transmit(:,:,:) = 0.0d0

    !!$    call omp_set_num_threads(1)
    !!$    call omp_set_num_threads( omp_get_num_procs() )
    if( mpi_rank==0 ) then
      if( omp_get_max_threads()>1 ) then
        write(fd_stdout,'(i3," OpenMP threads share loops." )') omp_get_max_threads()
      end if

      write(fd_stdout,*)
      write(fd_stdout,'(" Start calculate transmission...")')
      write(fd_stdout,*)
      call flush(6)
    end if

    !$OMP PARALLEL DO PRIVATE(ik,ib,ia) SHARED(transmit) COLLAPSE(3)
    do ik=1, anphon%nkpoint
      do ib=1, anphon%nbranch
        do ia=1, param%nangle-1
        !!$          call init_random_seed() ! for DEBUG
          call calcTransmit( transmit(ib,ik,ia), ib, ik, ia )
        end do ! ia
      end do ! ib
    end do ! ik
    !$OMP END PARALLEL DO

#ifdef MPI
    allocate( transmit_mpi(anphon%nbranch,anphon%nkpoint,param%nangle) )
    transmit_mpi = transmit
    call MPI_Allreduce( transmit_mpi, transmit, size(transmit), MPI_DOUBLE_PRECISION, &
         MPI_SUM, MPI_COMM_WORLD, mpi_ierr )
    deallocate( transmit_mpi )
#endif

    if ( mpi_rank == 0 .and. verbose >= 1) then
       trans_mean = 0.0
       do ia=1, param%nangle-1
         trans_mean = trans_mean + transmit(1,1,ia)
         if (transmit(1,1,ia) > 1) then 
           write(fd_stdout,*) "Error: trans is too large", transmit(1,1,ia)
         end if
       end do
       trans_mean = trans_mean / (param%nangle -1)
       write(fd_stdout,*) "- Trans mean: ", trans_mean
    end if 

  end subroutine TransmissionFunction

  !>
  !! @brief calculate a transmission function
  subroutine calcTransmit( trans, cbranch, ckpoint, angle )
    real(8), intent(out) :: trans
    integer, intent(in)  :: cbranch, ckpoint, angle

    integer :: particle
    integer :: collision, Jreflect
    integer :: StartSide, EndSide !1,2,or 0
    real(8) :: Start_Position(3), Start_Direction(3)
    real(8) :: Position(3), Direction(3), dir_old(3), tmpNormal(3), tmpDisp(3)
    real(8) :: Distance, FreePath
    integer :: igrain, isurface, tmpigrain, tmpisurface
    integer :: whichV, tmpwhichV
    integer, save :: id = 0

    type(Surface_type), pointer :: surface

    trans = 0.0d0
    ! [MPI] parallilized by sampling particles.
    ! [MPI] resultant value "trans" should be mpi_allreduced.
    do particle=1+mpi_rank, param%nphonon, mpi_size
      !$OMP MASTER
      id=id+1
      !$OMP END MASTER

      call createRandomPosition( Start_Position )
      call createRandomDirectionPhi( Start_Direction, angle )

      if( Start_Position(1)==0.0d0 ) then
        StartSide = 1 ! it means inject a phonon from X- outside surface
        isurface  = 1 ! the 1st surface is X- outside surface
      else
        StartSide = 2 ! it means inject a phonon from X+ outside surface
        isurface  = 2 ! the 2nd surface is X+ outside surface
      end if

      !---- determine end side, seems not necessary now
      if( Start_Direction(1)<=0.0d0 .and. StartSide==1 ) then 
        ! going out from X- outside surface
        EndSide = 1 ! X- outside surface
      else if( Start_Direction(1)>=0.0d0 .and. StartSide==2 ) then 
        ! going out from X+ outside surface
        EndSide = 2 ! X+ outside surface
      else
        EndSide = 0 ! inside simulation box
      end if

      Position  = Start_Position
      Direction = Start_Direction
      Distance  = 0.0d0
      igrain    = igrain_box
      whichV    = igrain_box
      collision = 1

      !$OMP MASTER
      if ( fd_trace>0 .and. mpi_rank == 0 .and. mod(id,interval_trace)==1 ) then
        write(fd_trace,*) real(Position) ! check trace
      end if
      !$OMP END MASTER

      if( ngrain>0 .and. param%kindGrain == ID_GRAIN_VORONOI ) then
        call findNearestIncidentGrain( whichV, Start_Position )
      else if( ngrain>0 .and. param%kindGrain == ID_GRAIN_BRICKS ) then
        call findIncidentCollisionGrain( whichV, Start_Position )
      end if ! kindGrain
      
      if(verbose >=1) then
        call findInsideWhichGrain(tmpwhichV, Start_Position)
        if(tmpwhichV /= whichV) then !make sure the initial grain is correct
            write(*,*) ""
            write(*,*) "Error: initial which grain is not consistent"
            write(*,*) "whichV = ", whichV, "tmpwhichV = ", tmpwhichV
            write(*,*) "Start_Pos = ", Start_Position
            exit
        end if
      end if

     
      !---- Initialize the FreePath
      if( param%MeanFreePath_current==1.0d10 ) then
        FreePath = 1.0d10
      else
        call createRandomFreePath( FreePath )
      end if

      ! loop until the phonon has gone outside the simulation box
      do while (EndSide==0)
        !! Exception check, make sure the particle is inside the simulation box 
        if (verbose >= 1) then
          if( .not. isInsideBox( Position ) ) then
            write(fd_stdout,*) " "
            write(fd_stdout,*) "Error, particle is out of box"
            write(fd_stdout,*) 'grain = ',igrain, "isurface = ", isurface
            write(fd_stdout,*) Position
            exit ! skip this phonon
          end if
        end if !! verbose

        collision = collision + 1

        call findCollisionSurface( Distance, igrain, &
            isurface, whichV, Position(:), Direction(:) )
            
        ! ---- collide with an internal surface , debug purpose
        if(verbose >= 3) then
            surface => vgrain(igrain)%vsurface(isurface)
            call calcSurfaceNormal( tmpNormal(:), &
                Position(:) + (Distance - SMALL_VAL)*Direction(:), surface )
            if( dot_product(Direction(:), tmpNormal(:)) <= -0 ) then
                write(*,*) ""
                write(*,*) "ray direction is not the same as the surface normal direction"
                write(*,*) "igrain = ", igrain, "isurface = ", isurface, "kind =", surface%kindShape
                write(*,*)  isInsideSurface( surface, Position(:) + Distance*Direction(:))
                write(*,*) "Position old", Position(:) 
                write(*,*) "Position new", Position(:) + Distance*Direction(:)
                write(*,*) "direction = ", Direction(:)
                write(*,*) "surface%normal = ", tmpNormal(:)
                write(*,*) "dot product val: ", dot_product(Direction(:), tmpNormal(:))
                write(*,*) "fligh distance =", Distance
                write(*,*) "collision = ", collision
                call findInsideWhichGrain(tmpwhichV, Position - Distance*Direction(:))
                write(*,*) "inside which grain",  tmpwhichV
                exit
            end if
        end if

        Jreflect = 1
        !---- update the position and direction 
        if( FreePath<Distance ) then    !! pp scattering
          Distance = FreePath
          Position(:) = Position(:) + Distance*Direction(:)
          call createRandomDirectionAllSphere( Direction )
          
          ! TODO: why for p-p scattering, set isurface = 0 and igrain = 0 ???
          isurface = 0 ! no surface collides
          igrain   = 0 ! no grain collides
        else
          ! position where phonon collides with the surface
          Position(:) = Position(:) + Distance*Direction(:)
          
          ! ---- Check which surface the phonon hits
          if( igrain==0 .and. isurface==1 ) then
            ! ---- hix X- surface of the box
            EndSide = 1
          else if( igrain==0 .and. isurface==2 ) then
            ! ---- hix X- surface of the box
            EndSide = 2
          else if( igrain==-1 .and. isurface ==-1) then   
            !! TODO, will remove later, scatter_flag == "virtual surface" 
            ! collide with -1 surface means phonon going out from X+ surface.
            !??? MIZUHO: isurface == -1, it is out of the range of array vsurface.
            !! Shao: 'surfaceNew == -1' is a temporary feature that I added to the code
            !! it is a "virtual" surface, and the norm direction is [100]
            surface => vgrain(igrain)%vsurface(isurface)
            call reflectVirtual( Jreflect, Direction, surface )
          else
            ! ---- collide with an internal surface 
            surface => vgrain(igrain)%vsurface(isurface)

            if( determinSpecular(surface) ) then
              call reflectSpecular( Jreflect, Direction, Position, surface )
            else
              call reflectDiffusive( Jreflect, Direction, Position, surface )
            end if
          end if
        end if ! FreePath
        

        surface => vgrain(igrain)%vsurface(isurface)
        ! move forward by a little
        tmpDisp = 1.0e-4*(Direction(:) + (surface%center(:) - Position(:))/&
            sqrt(dot_product(surface%center(:) - Position(:),surface%center(:) - Position(:))) ) 
        tmpDisp = 1.0e-7*Direction(:)
        Position(:) = Position(:) + tmpDisp

        ! ---- Updates the FreePath 
        if( param%MeanFreePath_current==1.0d10 ) then   !! Ballistic transport
          FreePath = 1.0d10
        else
          FreePath = FreePath - Distance
          if (FreePath < 1.0d-10) then
            call createRandomFreePath( FreePath )
          end if 
        end if
        ! ---- Update whichV
        !if( ngrain>0 ) then
          if( Jreflect==0 ) then  !--- transmitted, update grain index
            surface => vgrain(igrain)%vsurface(isurface)
            whichV = surface%grain_outside
          end if
        !end if ! ngrain
        
        !---- Verifying the grain ID
        if (verbose >=3) then
            call findInsideWhichGrain(tmpwhichV, Position)
            if(tmpwhichV /= whichV) then
              write(*,*) ""
              write(*,*) "whichV=", whichV
              write(*,*) "after  move whichV=", tmpwhichV
              call findInsideWhichGrain(tmpwhichV, Position-tmpDisp)
              write(*,*) "before move whichV=", tmpwhichV
              write(*,*)  "collision =", collision
              write(*,*) "position = ", Position
              do tmpigrain=1, ngrain 
                write(*,*) "is inside grain",tmpigrain, &
                    isInsideGrain(vgrain(tmpigrain), Position)
              end do
              exit
            end if
        end if
        
        !$OMP MASTER
        if ( fd_trace>0 .and. mpi_rank == 0 .and. mod(id,interval_trace)==1) then
          write(fd_trace,*) real(Position) ! check trace
        end if
        !$OMP END MASTER

      end do ! while(EndSide==0)

      !$OMP MASTER
      if ( fd_trace>0 .and. mpi_rank == 0 .and. mod(id,interval_trace)==1 ) then
        write(fd_trace,*) ! check trace
        write(fd_trace,*) ! check trace
      end if
      !$OMP END MASTER

      if( StartSide /= EndSide ) then
        trans = trans + 1.0d0
      end if
    end do ! particle

    trans = trans/dble(param%nphonon)

    return

  contains
    !>
    !! @brief calculate a position by random
    subroutine createRandomPosition( position )
      real(8), intent(out) :: position(3)

      real(8) :: y, z

      call random_number(y)
      call random_number(z)
      position(:) = [ &
        0.0d0, &
        (y)*param%boxy*param%scopey, &
        (z)*param%boxz*param%scopez ]
    end subroutine createRandomPosition

    !>
    !! @brief calculate a direction by random
    subroutine createRandomDirectionPhi( direction, itheta )
      real(8), intent(out) :: direction(3)
      integer, intent(in) :: itheta

      real(8) :: theta, phi

      call random_number(phi)
      theta = param%vangle(itheta)
      phi   = phi*M_TWOPI

      direction(:) = [ cos(theta), sin(theta)*cos(phi), sin(theta)*sin(phi) ]
    end subroutine createRandomDirectionPhi

    
    logical function isInsideSurface( surface, position)
       real(8), intent(in) :: position(3)
       type(Surface_type), intent(in) :: surface
       
       isInsideSurface = .false.
       if (surface%kindShape == ID_SHAPE_POLYGON) then
          isInsideSurface = &
          dot_product(position(:) - surface%center(:), surface%normal(:)) < 0
       else if (surface%kindShape == ID_SHAPE_SPHERE) then
          isInsideSurface = &
              dot_product(position - surface%center, position - surface%center) < surface%radius**2
       else if (surface%kindShape == ID_SHAPE_CYLINDER) then
          isInsideSurface = &
              (position(1)-surface%center(1))**2 + (position(2)-surface%center(2))**2 < surface%radius**2
      else
          write(*,*) "Error: unknow surface type, isInsideSurface"
      end if
    end function isInsideSurface

    logical function isInsideGrain(grain, position)
      type(Grain_type), intent(in) :: grain
      real(8), intent(in) :: position(3)
      integer :: isurface
      type(Surface_type), pointer :: surface
      
      do isurface=1,grain%nsurface
          surface => grain%vsurface(isurface)
          if( .not. surface%internal ) cycle ! do not set the following properties for this surface.
          if( surface%incident ) cycle
          if (.not. isInsideSurface(surface, position) ) then
            isInsideGrain = .false.
            return
          end if
      end do !isurface
      isInsideGrain  = .true.
      return
    end function isInsideGrain

   
    subroutine findInsideWhichGrain(whichV, Position)
      real(8), intent(in) :: Position(3)
      integer, intent(inout) :: whichV
      type(Grain_type), pointer :: grain
      type(Surface_type), pointer :: surface
      integer :: ig, is, igmin, ismin
      real(8) :: dist, dist_min
      integer :: count 
      
      count = 0
      
      do ig=1,ngrain  ! grain0 is the box
        grain => vgrain(ig)
        if(isInsideGrain(grain, Position)) then
            whichV = ig
            count = count + 1
        end if
      end do
      if (count == 1) then
          return 
      else if (count == 0 .and. param%kindGrain == ID_GRAIN_BRICKS ) then
        whichV = ngrain + 1  ! fill grain
        return
      else if (count == 0 .and. param%kindGrain == ID_GRAIN_VORONOI ) then
        write(*,*) "Error: cannot find the inside which grain for Voronoi structure"
        !---- calc minimum distance to surface
        dist_min = 1.0d10
        do ig=0,ngrain
          grain => vgrain(ig)
          do is=1, grain%nsurface
            surface => grain%vsurface(is)
            call calcDistToSurface(surface,Position, dist)
            if (abs(dist) < dist_min) then
                dist_min = abs(dist)
                igmin = ig
                ismin = is
            end if
          end do !nsurface
        end do !ngrain
        write(*,*) "dist_min to surface", dist_min
        write(*,*) "igmin = ", igmin, "ismin = ", ismin
        write(*,*) Position
        !---- calc minimum distance to surface
        whichV = -1
        return
      else if (count >= 2) then
        write(*,*) "=================================="
        write(*, *) "find inside more than grain, count=", count
        do ig=1,ngrain  ! grain0 is the box
          grain => vgrain(ig)
          write(*,*) "check inside which grain, igrain=", ig, isInsideGrain(grain, Position)
          if(isInsideGrain(grain, Position) ) then
            do is=1, grain%nsurface
              surface => grain%vsurface(is)
              call calcDistToSurface(surface,Position, dist)
              write(*,*) " "
              write(*,*) "dist to surface ", is, dist
              write(*,*) "surface normal ", surface%normal
              write(*,*) "vector:", Position(:) - surface%vertex(:,1)
              write(*,*) "ref pos",  surface%vertex(:,1)
              write(*,*) "Position",  Position(:)
              write(*,*) "inside surface: ", isInsideSurface(surface, Position)
              
            end do
          end if
        end do ! ngrain
        write(*,*) "whichV = ", whichV
        write(*,*) "position = ", position
        return
      else
        write(*,*) "unknown situation"
      end if
    end subroutine findInsideWhichGrain
    
    !>
    !! @brief find a surface which this phonon collides
    subroutine findCollisionSurface( distance, igrain, isurface, &
            whichV, Position, Direction )
      real(8), intent(out) :: distance
      integer, intent(inout) :: igrain
      integer, intent(inout) :: isurface
      integer, intent(in) :: whichV
      real(8), intent(in) :: Position(3)
      real(8), intent(in) :: Direction(3)

      integer :: ig, is
      type(Grain_type), pointer :: grain
      type(Surface_type), pointer :: surface
      integer :: isurface_min, igrain_min
      real(8) :: distance_min, distancetocollision
      logical :: hit, hit_flag
      logical :: isInPolygon
      
      distance_min = 1.0d10
      isurface_min = isurface
      igrain_min   = igrain
     
      hit_flag = .false.
      
      do ig=0,whichV,whichV
        grain => vgrain(ig)
        surfaceloop: do is=1, grain%nsurface
          surface => grain%vsurface(is)
          
          if( .not. surface%internal ) cycle
          
          call calcSurfaceDistance( hit, distancetocollision, &
            Position(:), Direction(:), surface )

          if( hit .and. distance_min > distancetocollision ) then
            hit_flag = .true.
            distance_min = distancetocollision
            isurface_min  = is
            igrain_min    = ig
          end if
        end do surfaceloop

        !  !! Check the virtual interface for the last grain
        !!! it is an temperary extention for one preject, I will check it
        if (FLAG_ADD_VIRTUAL_INTERFACE .eqv. .true.) then
          if ( param%kindGrain == ID_GRAIN_BRICKS .and. whichV == igrain_fill )then
            call calcSurfaceDistanceVirtualInterface( &
              hit, distancetocollision, Position(:), Direction(:))

            if( hit .and. distance_min > distancetocollision ) then
              distance_min = distancetocollision
              isurface_min  = -1
              igrain_min    = -1
            end if
          end if
        end if ! Virtual interface
      end do ! ig=0,whichV,whichV
      
      ! for debug purpose, can be deleted
      if (.false.) then
        if(surface%kindShape == ID_SHAPE_POLYGON) then
          if(  dot_product(Direction(:), &
               vgrain(igrain_min)%vsurface(isurface_min)%normal(:)) <= -0) then
               write(*,*) "hit a surface at wrong direction"
               write(*,*) "igrain =", igrain_min, "isurf=", isurface_min
               write(*,*) "surface norm",  vgrain(igrain_min)%vsurface(isurface_min)%normal(:)
               do is=1, vgrain(igrain_min)%nsurface
                 surface => vgrain(igrain_min)%vsurface(is)
                 call calcSurfaceDistance( hit, distancetocollision, &
                 Position(:), Direction(:), surface )
                 
                 write(*,*) " "
                 write(*,*) "hit = ", hit
                 write(*,*) 'is=', is, "dist = ", distancetocollision
                 write(*,*) "norm", surface%normal
                 write(*,*) "hit pos =", Position + distancetocollision*Direction
                 write(*,*) "ig=", igrain_min, "is =", is, "nvertex=", surface%nvertex
                 isInPolygon = isInsidePolygon(surface%vertex, surface%nvertex, &
                     Position + distancetocollision*Direction)
                 write(*,*) isInPolygon
               end do
          end if
        end if
      end if
    
      if (hit_flag .eqv. .false.) then
            write(*,*) "Error, cannot find collition surface"
            write(*,*) "distance_min", distance_min 
        do ig=0,whichV,whichV
          grain => vgrain(ig)
          do is=1, grain%nsurface
             surface => grain%vsurface(is)
             if( .not. surface%internal ) cycle
             call calcSurfaceDistance( hit, distancetocollision, &
               Position(:), Direction(:), surface )
             write(*,*) "ig =", ig, "is =", is
             write(*,*) "surface norm = ", surface%normal
             write(*,*) "ray vector" , Direction
             write(*,*) "dis to collision = ", distancetocollision
             write(*,*) "hitted position", Position(:)  + distancetocollision * Direction
             write(*,*) "dot product val= ", dot_product( &
                 surface%vertex(:,1)- Position(:), surface%normal(:) )
          end do !is
        end do ! whichV
      end if
      distance = distance_min
      igrain   = igrain_min
      isurface = isurface_min
    end subroutine findCollisionSurface

    !>
    !! @brief find the nearest grain which has incident surface from this phonon, for voronoi
    subroutine findNearestIncidentGrain( igrain_nearest, position )
      integer, intent(out) :: igrain_nearest
      real(8), intent(in) :: position(3)

      real(8) :: distance2, distance2_min
      integer :: igrain
      type(Grain_type), pointer :: grain

      distance2_min = 1.0d10
      do igrain=1, ngrain
        grain => vgrain(igrain)
        if( .not. grain%incident ) cycle
        distance2 = sum( (position(:) - grain%generatrix(:))**2 )
        if( distance2_min>distance2 ) then
          distance2_min = distance2
          igrain_nearest = igrain
        end if
      end do ! igrain
    end subroutine findNearestIncidentGrain


    subroutine calcDistToSurface(surface, position, dist)
      type(Surface_type), intent(in) :: surface
      real(8), intent(in) :: position(3)
      real(8), intent(out) :: dist

      dist = (dot_product(position(:) - surface%vertex(:,1), surface%normal(:)))
    end subroutine calcDistToSurface

    !>
    !! @brief find a grain which has incident surface and this phonon attaches, for Bricks
    subroutine findIncidentCollisionGrain( igrain_incident, position )
      integer, intent(out) :: igrain_incident
      real(8), intent(in) :: position(3)

      integer :: igrain, isurface, ivertex
      type(Grain_type), pointer :: grain
      type(Surface_type), pointer :: surface
      real(8) :: vv(3), vp(3), outer(3)
      logical :: found_one, found_min

      found_min = .false.
      ! find an incident surface that involves Start_Position
      do igrain=1, ngrain
        grain => vgrain(igrain)
        if( .not. grain%incident ) cycle
        do isurface=1, grain%nsurface
          surface => grain%vsurface(isurface)
          if( .not. surface%incident ) cycle
          
          found_one =  isInsidePolygon(surface%vertex, surface%nvertex, position)
          ! if found an incident surface whichV is set to the grain of the surface
          if( found_one ) then
            igrain_incident = igrain
            found_min = .true.
            exit
          end if
        end do ! isurface
      end do ! igrain

      if( .not. found_min ) then
        igrain_incident = igrain_fill
      end if
    end subroutine findIncidentCollisionGrain

    !>
    !! @brief calculate a free path by random
    subroutine createRandomFreePath( FreePath )
      real(8), intent(out) :: FreePath
      real(8) :: rand

      call random_number(rand)
      FreePath = -param%MeanFreePath_current*log(rand)
    end subroutine createRandomFreePath

    !>
    !! @brief check this position locates inside the simulation box
    logical function isInsideBox( position )
      real(8), intent(in) :: position(3)

      isInsideBox = &
        position(1)>=0.0d0 - SMALL_VAL .and. &
        position(1)<=param%boxx + SMALL_VAL .and.&
        position(2)>=0.0d0 - SMALL_VAL .and. &
        position(2)<=param%boxy + SMALL_VAL .and.&
        position(3)>=0.0d0 - SMALL_VAL .and. &
        position(3)<=param%boxz + SMALL_VAL
    end function isInsideBox


    !>
    !! @brief calculate a direction in all sphere by random
    subroutine createRandomDirectionAllSphere( direction )
      real(8), intent(out) :: direction(3)

      real(8) :: costheta, sintheta, phi

      call random_number(costheta)
      call random_number(phi)
      costheta = costheta*2.0d0 - 1.0d0
      phi = phi*M_TWOPI

      sintheta = sqrt(1.0d0-costheta**2)
      direction(:) = [ costheta, sintheta*cos(phi), sintheta*sin(phi) ]
    end subroutine createRandomDirectionAllSphere

    !>
    !! @brief calculate a direction in semi-sphere by random
    subroutine createRandomDirectionSemiSphere( direction )
      real(8), intent(out) :: direction(3)

      real(8) :: sintheta, costheta, phi

      call random_number(sintheta)

      sintheta = sqrt(sintheta)
      costheta = sqrt(1.0d0-sintheta**2)

      call random_number(phi)
      phi = phi*M_TWOPI
      direction(:) = [ sintheta*cos(phi), sintheta*sin(phi), costheta ]
    end subroutine createRandomDirectionSemiSphere

    !>
    !! @brief determin reflection type is specular or not by random
    logical function determinSpecular( surface )
      type(Surface_type), intent(in) :: surface

      real(8) :: value

      call random_number(value)
      determinSpecular = ( value<=surface%specularity )
      !write(fd_stdout,*) "in determiSpec surface speculartiy", surface%specularity
      !write(fd_stdout,*) "specularity = ", determinSpecular

    end function determinSpecular

    !>
    !! @brief reflect phonon by viritual surface
    subroutine reflectVirtual( reflect, direction, surface )
      integer, intent(out) :: reflect
      real(8), intent(inout) :: direction(3)
      type(Surface_type), intent(in) :: surface

      real(8) :: tmpRotation(3,3), tmpDirection(3)
      real(8) :: rand

      ! for specular scatter
      call calcSurfaceRotationVirtualInterface( tmpRotation(:,:) )

      !! SHAO: subroutine `set_direct_diffuse_surface`
      call createRandomDirectionSemiSphere( tmpDirection )
      !! tmpRotation(:,:) is used to obtain the norm vector of the surface

      ! hit from backside
      !! here dot_product(DirectionPre, [100]))
      if( dot_product(direction(:),surface%normal(:)) > 0.0d0 ) then
        ! Direction should be opposite side against normal
        direction(:) = - matmul(tmpRotation(:,:),tmpDirection(:))
      else
        direction(:) = + matmul(tmpRotation(:,:),tmpDirection(:))
      end if
      !! end of subroutine set_direct_diffuse_surface

      call random_number(rand)
      if( rand<=0.5d0) then !!Pass
        direction(:) = -direction(:)
        reflect = 1  ! It is a virtual surface, so Jreflect is always 1
      else
        reflect = 1
      end if
    end subroutine reflectVirtual

    !>
    !! @brief reflect phonon as specular scattering
    subroutine reflectSpecular( reflect, direction, position, surface )
      integer, intent(out) :: reflect
      real(8), intent(inout) :: direction(3)
      real(8), intent(in) :: position(3)
      type(Surface_type), intent(in) :: surface

      real(8) :: tmpNormal(3)

      if( determinTransmit(surface) ) then
        ! not change direction(:)
        reflect = 0
      else
        call calcSurfaceNormal( tmpNormal(:), position(:), surface )
        !!$        write(200,*) real(position), real(tmpNormal) ! MIZUHO DEBUG check ray and normal
      
        !---- SHAO, for sphere and cylinder,  make sure the norm direction is correct 
        if( dot_product(tmpNormal(:) , direction) < 0) then
            tmpNormal(:) = -1.0d0*tmpNormal(:)
        end if
        
        if(verbose >=3) then
            if(.not. isEqual(dot_product(tmpNormal, tmpNormal), 1.0d0)) then
                write(*,*) "inside reflec specular, temp norm is not unit "
            end if
            if(.not. isEqual(dot_product(direction, direction), 1.0d0)) then
                write(*,*) "inside reflec specular, direction  is not unit before "
                write(*,*) dot_product(direction, direction)
            end if
        end if
        direction(:) = direction(:) - 2.0d0*tmpNormal(:) &
          * dot_product(tmpNormal(:),direction(:))
        reflect = 1
        if(verbose >=3) then
            if(.not. isEqual(dot_product(direction, direction), 1.0d0)) then
                write(*,*) "inside reflec specular, direction  is not unit "
                write(*,*) dot_product(direction, direction)
            end if
        end if
      end if
      ! SHAO: important here
      ! there is a bug before that the norm of direction is not equal 1
      if(verbose >=3) then
        if ( abs(sqrt(dot_product(direction,direction)) -1.0d0) > 1.0d-8) then
            write(*,*) "direction is not very accurate at reflectDiffusive"
            write(*,*) dot_product(direction,direction)
            direction = 1.0d0/sqrt(dot_product(direction,direction))*direction
            if ( abs(sqrt(dot_product(direction,direction)) -1.0d0) > 1.0d-8) then
              write(*,*) "direction is not very accurate", "reflectSpecular"
              write(*,*) dot_product(direction,direction)
              write(*,*) abs(sqrt(dot_product(direction,direction)) -1.0d0)
            end if
        end if
      end if
    end subroutine reflectSpecular

    !>
    !! @brief reflect phonon as diffusive scattering
    subroutine reflectDiffusive( reflect, direction, position, surface )
      integer, intent(out) :: reflect
      real(8), intent(inout) :: direction(3)
      real(8), intent(in) :: position(3)
      type(Surface_type), intent(in) :: surface

      real(8) :: tmpNormal(3), tmpRotation(3,3), tmpDirection(3)
      real(8) :: rand

      !! SHAO: diffusive reflected, calc the new direction based on "tmpRotation"
      call calcSurfaceNormal( tmpNormal(:), position(:), surface ) ! MIZUHO fixed
      !---- SHAO, make sure the norm direction is correct 
      if( dot_product(tmpNormal(:) , direction) < 0.0d0) then
            tmpNormal(:) = -1.0d0*tmpNormal(:)
      end if
      call calcSurfaceRotation( tmpRotation(:,:), position(:), surface )

      call createRandomDirectionSemiSphere( tmpDirection )
      if( dot_product(direction(:),tmpNormal(:)) > 0.0d0 ) then
        direction(:) = - matmul(tmpRotation(:,:),tmpDirection(:))
      else
        direction(:) = + matmul(tmpRotation(:,:),tmpDirection(:))
      end if
      
      !---- SHAO: ensure the direction is correct 
      if( dot_product(direction(:),tmpNormal(:)) > 0.0d0 ) then
          direction(:) = -1.0d0*direction(:)
      end if

      call random_number(rand)
      if( determinTransmit(surface,rand) ) then
        direction(:) = - direction(:)
        reflect = 0
      else
        reflect = 1
      end if
      
      if(verbose >=3) then
        if ( .not. isEqual(sqrt(dot_product(direction,direction)), 1.0d0)) then
            write(*,*) "direction is not very accurate ", "reflectDiffusive"
            direction = 1.0d0/sqrt(dot_product(direction,direction))*direction
            if ( .not. isEqual(sqrt(dot_product(direction,direction)), 1.0d0)) then
              write(*,*) "direction is not very accurate", "reflectDiffusive"
              write(*,*) dot_product(direction,direction)
              write(*,*) abs(sqrt(dot_product(direction,direction)) -1.0d0)
            end if
        end if
      end if
    end subroutine reflectDiffusive

    !>
    !! @brief determin this phonon transit a surface or reflect on a surface by random
    logical function determinTransmit( surface, rand )
      type(Surface_type), intent(in) :: surface
      real(8), optional, intent(in) :: rand

      real(8) :: threshold, value
      real(8) :: omega

      select case( surface%kindTrans )
      case( ID_TRANS_VALUE )
        threshold = surface%transmissivity(1,1)

      case( ID_TRANS_TBC, ID_TRANS_GAMMA )
        if( param%nmaterial==1 ) then
          threshold = surface%transmissivity(cbranch,ckpoint)
        else
          call random_number(omega)
          threshold = surface%transmissivity(1,1) &
            + omega*(1-surface%transmissivity(1,1))
        end if
      case default
        threshold = 0.0d0
      end select

      if( present(rand) ) then
        value = rand
      else
        call random_number(value)
      end if
      determinTransmit = ( value<=threshold )
    end function determinTransmit
  end subroutine calcTransmit

  !>
  !! @brief calculate effective mean free pathes
  subroutine effectiveMFP( EMFP )
    real(8), intent(out) :: EMFP(anphon%nkpoint,anphon%nbranch)

    integer :: ia, ib, ik
    real(8) :: IntensityTransmission, dS
    real(8) :: WT(param%nangle)

    dS = (param%vangle(2) - param%vangle(1)) / 3.0d0 * 3.0d0*0.5d0*param%boxx

    !$OMP PARALLEL DO PRIVATE(ik,ib,ia,IntensityTransmission) COLLAPSE(2)
    do ik=1, anphon%nkpoint
      do ib=1,anphon%nbranch
        do ia=1,param%nangle
          WT(ia) = sin(param%vangle(ia))*cos(param%vangle(ia)) * transmit(ib,ik,ia)
        end do ! ia

        IntensityTransmission = 0.0d0
        do ia=2,param%nangle-1,2
          IntensityTransmission = IntensityTransmission &
            + (WT(ia-1) + 4.0d0*WT(ia) + WT(ia+1))
        end do ! ia
        if(ik == 1 .and. ib == 1) then  ! SHAO, check the transmission
            write(*,*) " averaged tranmission = ", IntensityTransmission * &
                (param%vangle(2) - param%vangle(1)) / 3.0d0
        end if
        EMFP(ik,ib) = IntensityTransmission * dS
      end do ! ib
    end do ! ik
    !$OMP END PARALLEL DO

    if( mpi_rank == 0 ) then
      write(fd_stdout,*) '     RTMC Simulation has completed.'
    endif
  end subroutine effectiveMFP

  subroutine clearGeometry
    integer :: igrain, isurface
    type(Grain_type), pointer :: grain
    type(Surface_type), pointer :: surface

    select case( param%kindGrain )
    case(ID_GRAIN_VORONOI)
      do igrain=0, ngrain
        grain => vgrain(igrain)
        do isurface=1, grain%nsurface
          surface => grain%vsurface(isurface)
          if( associated(surface%vertex) ) deallocate(surface%vertex)
          if( associated(surface%transmissivity) ) deallocate(surface%transmissivity)
        end do
        deallocate(grain%vsurface)
      end do
      deallocate(vgrain)

    case(ID_GRAIN_BRICKS)
      do igrain=0, ngrain+1
        grain => vgrain(igrain)
        do isurface=1, grain%nsurface
          surface => grain%vsurface(isurface)
          if( associated(surface%vertex) ) deallocate(surface%vertex)
          if( associated(surface%transmissivity) ) deallocate(surface%transmissivity)
        end do
        deallocate(grain%vsurface)
      end do
      deallocate(vgrain)

    end select
  end subroutine clearGeometry

  subroutine clearArray

    if( allocated(param%vangle) ) deallocate(param%vangle)
    if( allocated(param%vMeanFreePath) ) deallocate(param%vMeanFreePath)
    if( allocated(param%vmaterial) ) deallocate(param%vmaterial)
    if( allocated(param%kindTrans) ) deallocate(param%kindTrans)
    if( allocated(param%transmissivity) ) deallocate(param%transmissivity)
    if( allocated(param%specularity) ) deallocate(param%specularity)

    if( allocated(transmit) ) deallocate(transmit)

  end subroutine clearArray

  subroutine GammaDeterm( Material, TBC, gamma, ModeNum )
    character(len=32), intent(in) :: Material
    real(8), intent(in) :: TBC
    real(8), intent(out) :: gamma
    integer, intent(in) :: ModeNum

    integer :: ii, jj, kk, nite, rloopend, recalc, nqibz, nqfbz
    real(8) :: eps, differential
    real(8) :: tmpX, tmpGamma, tmpTBC, tmpTBCa, tmpdTBC, tmpF, tmpdF, tmpH
    real(8), allocatable :: cap(:,:)

    real(8) :: tbc1, tbc2, testTrans, tmpCond

    type(Anphon_type) :: anphon
    integer :: status

    recalc = 0 ! 0:Recalculation off , 1:Recalculation on
    nite = 10000 ! maximum of iteration
    eps = 1.0d-6

    call loadAnphon( anphon, Material, status )

    rloopend = anphon%nkpoint * anphon%nbranch
    allocate( cap(anphon%nbranch,anphon%nkpoint) )
    cap = 0.0d0

    nqibz = rloopend
    nqfbz = int(sum(anphon%weight(:,:))/anphon%nbranch)

    if( ModeNum==1 ) then
      do ii=1,anphon%nkpoint
        do jj=1,anphon%nbranch
          tmpX = hbar*anphon%omega(jj,ii)/(kB*param%temperature)
          if( tmpX<=1.0d-10 ) then
            cap(jj,ii) = kB
          else
            cap(jj,ii) = kB * ((tmpX/2)/sinh(tmpX/2))**2
          end if
        end do ! jj
      end do ! ii

      tmpGamma = 1.0d-2
      do ii=1,nite
        tmpTBC = 0.0d0
        tmpTBCa = 0.0d0
        tmpdTBC = 0.0d0
        do jj=1,anphon%nkpoint
          do kk=1,anphon%nbranch
            tmpH = 1/(tmpGamma*(anphon%omega(kk,jj)/anphon%maxOmega) + 1)
            tmpTBC = tmpTBC &
              + anphon%weight(kk,jj)*cap(kk,jj)*anphon%vel(kk,jj)*tmpH/2/anphon%vol/nqfbz
          end do ! kk
        end do ! jj

        do jj=1,anphon%nkpoint
          do kk=1,anphon%nbranch
            tmpH = 1
            tmpTBCa = tmpTBCa &
              + anphon%weight(kk,jj)*cap(kk,jj)*anphon%vel(kk,jj)*tmpH/2/anphon%vol/nqfbz
          end do ! kk
        end do ! jj

        do jj=1,anphon%nkpoint
          do kk=1,anphon%nbranch
            tmpH = -(anphon%omega(kk,jj)/anphon%maxOmega) &
              /(tmpGamma*(anphon%omega(kk,jj)/anphon%maxOmega) + 1)**2
            tmpdTBC = tmpdTBC &
              + anphon%weight(kk,jj)*cap(kk,jj)*anphon%vel(kk,jj)*tmpH/2/anphon%vol/nqfbz
          end do ! kk
        end do ! jj

        tmpF = tmpTBC/(1-tmpTBC/tmpTBCa)-TBC
        tmpdF = tmpdTBC/(1-tmpTBC/tmpTBCa)**2
        differential = tmpF/tmpdF
        if( abs(differential)<=eps ) then
          exit
        end if
        tmpGamma = tmpGamma - differential
      end do ! ii
      gamma = tmpGamma
    end if ! ModeNum

    if( recalc>=4 .and. ModeNum==1 ) then
      testTrans = 0
      tbc1 = 0
      tbc2 = 0
      do ii=1,anphon%nkpoint
        do jj=1,anphon%nbranch
          testTrans = 1.0d0/(tmpGamma*anphon%omega(jj,ii)/anphon%maxOmega + 1)
          tbc1 = tbc1 + anphon%weight(jj,ii)*cap(jj,ii) &
            *anphon%vel(jj,ii)*testTrans/2/anphon%vol/nqfbz
          tbc2 = tbc2 + anphon%weight(jj,ii)*cap(jj,ii) &
            *anphon%vel(jj,ii)/2/anphon%vol/nqfbz
        end do ! jj
      end do ! ii
      tmpCond = tbc1/(1 - tbc1/tbc2)
      write(fd_stdout,*) gamma, tmpCond
    end if ! recalc

    call clearAnphon( anphon )
  end subroutine  GammaDeterm

  subroutine fileoutput( SEMFP, GeoFile )
    real(8), intent(inout) :: SEMFP(anphon%nkpoint,anphon%nbranch,param%nMeanFreePath)
    character(len=*), intent(in) :: GeoFile

    integer :: ii, jj, kk, s
    integer :: rloopend, nqibz, nqfbz, cfreq
    real(8) :: aveEMFP, interpolate
    real(8) :: tmpOmega, tmpTau, tmpMFP, tmpWeight, tmpEMFP
    real(8) :: tmpX, tmpfreq, tmpVel
    real(8) :: Thermal_conductivity(anphon%nbranch,anphon%nkpoint)
    real(8), allocatable :: cap(:,:), spectrumKappa(:,:)
    real(8), allocatable :: spectrumOmega(:)
    real(8), allocatable :: EMFP(:,:)
    character(len=256) :: FileName
    integer :: cPath

    type data_type
      integer :: ibranch
      integer :: iKpoint
      real(8) :: bulkMFP
      real(8) :: effectiveMFP
      real(8) :: conductivity
      real(8) :: conductivity_accum
    end type data_type

    type(data_type), allocatable :: vdata(:)
    type(data_type) :: tmpData

    cPath = index(GeoFile, "/", back=.true.)
    if( cPath==0 ) then
      cPath = index(GeoFile, "\\", back=.true.)
    end if

    write(FileName,"(a,'thermalproperty.out')") GeoFile(:cPath)
    write(fd_stdout,*) "SEMFP=", SEMFP(1,1,:)

    open(41,file=trim(FileName))
    write(41,"('# Datatype = thermalproperty')")
    write(41,"('# Incident_Phonon                               =   ',i6)") param%nphonon
    if( param%nmaterial==1 ) then
      aveEMFP = sum(SEMFP(:,:,:))/(anphon%nbranch*anphon%nkpoint*param%nMeanFreePath)
      write(41,"('# Averaged Effective_MeanFreePath [nm]          = ',f20.10)") aveEMFP
      write(fd_stdout,"('# Averaged Effective_MeanFreePath [nm]          = ',f20.10)") aveEMFP

      do ii=1,anphon%nkpoint
        do jj=1,anphon%nbranch
          anphon%tau(jj,ii) = anphon%tau(jj,ii)*1.0d-12*(anphon%RTemp/param%temperature)
          anphon%MFP(jj,ii) = anphon%tau(jj,ii)*anphon%vel(jj,ii)
        end do ! jj
      end do ! ii

      allocate(  cap(anphon%nbranch,anphon%nkpoint) )
      allocate( EMFP(anphon%nbranch,anphon%nkpoint) )
      cap = 0.0d0
      rloopend = anphon%nkpoint * anphon%nbranch
      nqibz = rloopend
      nqfbz = int(sum(anphon%weight(:,:))/anphon%nbranch)
      do ii=1,anphon%nbranch
        do jj=1,anphon%nkpoint-1
          do kk=jj+1,anphon%nkpoint
            if( anphon%omega(ii,jj)>anphon%omega(ii,kk) ) then

              tmpOmega = anphon%omega(ii,jj)
              anphon%omega(ii,jj) = anphon%omega(ii,kk)
              anphon%omega(ii,kk) = tmpOmega

              tmpWeight = anphon%weight(ii,jj)
              anphon%weight(ii,jj) = anphon%weight(ii,kk)
              anphon%weight(ii,kk) = tmpWeight

              tmpTau = anphon%tau(ii,jj)
              anphon%tau(ii,jj) = anphon%tau(ii,kk)
              anphon%tau(ii,kk) = tmpTau

              tmpVel = anphon%vel(ii,jj)
              anphon%vel(ii,jj) = anphon%vel(ii,kk)
              anphon%vel(ii,kk) = tmpVel

              tmpMFP = anphon%MFP(ii,jj)
              anphon%MFP(ii,jj) = anphon%MFP(ii,kk)
              anphon%MFP(ii,kk) = tmpMFP

              do s=1, param%nMeanFreePath
                tmpEMFP = SEMFP(jj,ii,s)
                SEMFP(jj,ii,s) = SEMFP(kk,ii,s)
                SEMFP(kk,ii,s) = tmpEMFP
              end do
            end if
          end do ! kk
        end do ! jj
      end do ! ii

      !$OMP PARALLEL DO PRIVATE(ii,jj) PRIVATE(tmpX) COLLAPSE(2)
      do ii=1,anphon%nbranch
        do jj=1,anphon%nkpoint
          tmpX = hbar*anphon%omega(ii,jj)/(kB*param%temperature)
          if( tmpX<=1.0d-10 ) then
            cap(ii,jj) = kB
          else
            cap(ii,jj) = kB * ((tmpX/2)/sinh(tmpX/2))**2
          end if
        end do ! jj
      end do ! ii
      !$OMP END PARALLEL DO

      !$OMP PARALLEL DO PRIVATE(ii,jj) COLLAPSE(2)
      do ii=1,anphon%nbranch
        do jj=1,anphon%nkpoint
          ! unit of MFP is [m], bulk MFP.
          ! unit of SEMFP is [nm], effective MFP sampled at some Mean Free Path.
          ! unit of vMeanFreePath is [nm], sample of Mean Free Path.
          ! unit of EMFP is [m], effective MFP interpolated by bulk MFP.
          if( param%kindMFP == ID_MFP_SAMPLE_MATTIESSEN ) then
            ! use Matthiessen's rule
            if( anphon%MFP(ii,jj)==0.0d0 ) then
              EMFP(ii,jj) = SEMFP(jj,ii,1)*1.0d-9
            else
              EMFP(ii,jj) = ((anphon%MFP(ii,jj))**(-1)+(SEMFP(jj,ii,1)*1.0d-9)**(-1))**(-1)
            end if
          else
            ! use MFP sampling method
            if( anphon%MFP(ii,jj)<param%vMeanFreePath(1)*1.0d-9 ) then
              ! under minimum of sampling vMeanFreePath 
              EMFP(ii,jj) = SEMFP(jj,ii,1)*1.0d-9 * anphon%MFP(ii,jj)/ (param%vMeanFreePath(1)*1.0d-9)
            else if( anphon%MFP(ii,jj)>param%vMeanFreePath(param%nMeanFreePath)*1.0d-9 ) then
              ! over maximum of sampling vMeanFreePath 
              EMFP(ii,jj) = SEMFP(jj,ii,param%nMeanFreePath)*1.0d-9 
            else
              ! find nearest sampling vMeanFreePath 
              do s=1, param%nMeanFreePath-1
                if( param%vMeanFreePath(s)*1.0d-9 <= anphon%MFP(ii,jj) .and. &
                  anphon%MFP(ii,jj) < param%vMeanFreePath(s+1)*1.0d-9 ) then

                  ! calculate ratio between nearest two samples
                  interpolate = &
                    (anphon%MFP(ii,jj) - param%vMeanFreePath(s)*1.0d-9) &
                    / ((param%vMeanFreePath(s+1) - param%vMeanFreePath(s)) * 1.0d-9)
                  ! interpolate effective MFP
                  EMFP(ii,jj) = SEMFP(jj,ii,s)*1.0d-9 * (1.0d0-interpolate) &
                    + SEMFP(jj,ii,s+1)*1.0d-9 * (interpolate)
                  exit
                end if
              end do
            end if
          end if

          Thermal_conductivity(ii,jj) &
            = cap(ii,jj)*anphon%weight(ii,jj) &
            *EMFP(ii,jj)*anphon%vel(ii,jj)/3/anphon%vol/nqfbz
        end do ! jj
      end do ! ii
      !$OMP END PARALLEL DO

      write(41,'("# Thermal_conductivity@",f0.6,"K [W m-1 K-1]  = ",f20.10)') &
        param%temperature ,sum(Thermal_conductivity(:,:))
      close(41)

      cfreq = anphon%nbranch*anphon%nkpoint
      allocate( vdata(cfreq) )

      ! copy values into a one-dimensional array.
      kk=0
      do ii=1,anphon%nbranch
        do jj=1,anphon%nkpoint
          kk=kk+1
          vdata(kk)%ibranch = ii
          vdata(kk)%iKpoint = jj
          vdata(kk)%bulkMFP = anphon%MFP(ii,jj)
          vdata(kk)%effectiveMFP = EMFP(ii,jj)
          vdata(kk)%conductivity = Thermal_conductivity(ii,jj)
          vdata(kk)%conductivity_accum = 0.0d0
        end do
      end do

      ! sort data by effectiveMFP in ascending order.
      do ii=1, cfreq
        do jj=ii+1, cfreq
          if( vdata(ii)%effectiveMFP > vdata(jj)%effectiveMFP ) then
            tmpData   = vdata(ii)
            vdata(ii) = vdata(jj)
            vdata(jj) = tmpData
          end if
        end do
      end do

      ! accumulate conductivity.
      vdata(1)%conductivity_accum = vdata(1)%conductivity
      do ii=2, cfreq
        vdata(ii)%conductivity_accum = &
          vdata(ii)%conductivity + vdata(ii-1)%conductivity_accum 
      end do

      ! output Effective MFP as a function of bulk MFP.
      write(FileName,"(a,'effectiveMFP.out')") GeoFile(:cPath)
      write(fd_stdout,"('output sorted effective MFP in ',a)") trim(FileName)
      open(102,file=trim(FileName))
      write(102,"('# Datatype = Effective Mean Free Path of phonon.')")
      write(102,"('# data are sorted by effective MFP in ascending order.')")
      write(102,"('#   total branches: ',i3)") anphon%nbranch
      write(102,"('#   total Kpoints:  ',i3)") anphon%nkpoint
      write(102,"('#')")
      write(102,"('# branch, kpoint, bulk MFP[nm], effective MFP[nm], ThermalConductivity, and its accumulated[W m-1 K-1].')")
      do ii=1, cfreq
        write(102,"(2i5,2e15.6,2e15.6)") &
          vdata(ii)%ibranch, vdata(ii)%iKpoint, &
          vdata(ii)%bulkMFP*1.0d9, vdata(ii)%effectiveMFP*1.0d9, &
          vdata(ii)%conductivity, &
          vdata(ii)%conductivity_accum
      end do
      close(102)

      deallocate( vdata )

      ! allocate arrays with the maximum size
      ! to prevent error occurs in the succesive code.
      allocate( spectrumOmega(cfreq),spectrumKappa(cfreq,anphon%nbranch) )

      cfreq = 0
      do ii=1,anphon%nbranch
        do jj=1,anphon%nkpoint
          cfreq = cfreq + 1
          spectrumOmega(cfreq) = anphon%omega(ii,jj) ! error occured
          do kk=1,cfreq-1
            if( spectrumOmega(cfreq)==spectrumOmega(kk) ) then
              cfreq = cfreq - 1
              exit
            end if
          end do ! kk
        end do ! jj
      end do ! ii

      do ii=1,cfreq-1
        do jj=ii+1,cfreq
          if( spectrumOmega(ii)>spectrumOmega(jj) ) then
            tmpfreq = spectrumOmega(ii)
            spectrumOmega(ii) = spectrumOmega(jj)
            spectrumOmega(jj) = tmpfreq
          end if
        end do ! jj
      end do ! ii

      spectrumKappa = 0.0d0

      !$OMP PARALLEL DO PRIVATE(ii,jj,kk) COLLAPSE(3) REDUCTION(+:spectrumKappa)
      do ii=1,cfreq-1
        do jj=1,anphon%nbranch
          do kk=1,anphon%nkpoint
            if( anphon%omega(jj,kk) == spectrumOmega(ii) ) then
              spectrumKappa(ii,jj) = spectrumKappa(ii,jj) + Thermal_conductivity(jj,kk)
            end if
          end do ! kk
        end do ! jj
      end do ! ii
      !$OMP END PARALLEL DO

      !! SHAO: phononconductivity.out -> spectralTC.out
      write(FileName,"(a,'phononconductivity.out')") GeoFile(:cPath)

      open(43,file=trim(FileName))
      write(43,"('# Datatype = phononconductivity')")
      write(43,"('# Nb = ',i6)") cfreq-1
      write(43,"('# Nbranch = ',i4)") anphon%nbranch
      write(43,"('# omega, fluxphononout(1:Nbranch)')")
      do ii=1,cfreq-1
        write(43,"(e20.10$)") spectrumOmega(ii)
        do jj=1,anphon%nbranch
          write(43,"(e20.10$)") spectrumKappa(ii,jj)
        end do ! jj
        write(43,*)
      end do ! ii
      close(43)

      deallocate( spectrumOmega, spectrumKappa )
      deallocate( EMFP )
      deallocate( cap )

    else if( param%nmaterial>1 ) then
      write(41,"('# Effective_MeanFreePath [nm]       = ',e20.10)") SEMFP(1,1,1)
    end if
    close(41)
  end subroutine fileoutput

  subroutine calcSurfaceDistance( hit, Lhit, Rstart, Ndir, surface )
    logical, intent(out) :: hit         ! hit or not
    real(8), intent(out) :: Lhit        ! length to the hitting position
    real(8), intent(in)  :: Rstart(3)   ! starting position of the ray
    real(8), intent(in)  :: Ndir(3)     ! normalized direction of the ray
    type(Surface_type), intent(in) :: surface

    select case( surface%kindShape )
    case(ID_SHAPE_POLYGON)
      call calcSurfaceDistancePolygon( &
        hit, Lhit, Rstart(:), Ndir(:), &
        surface%vertex(:,:), surface%normal(:), surface%nvertex )
    case(ID_SHAPE_SPHERE)
      call calcSurfaceDistanceSphere( &
        hit, Lhit, Rstart(:), Ndir(:), &
        surface%center(:), surface%radius )
    case(ID_SHAPE_CYLINDER)
      call calcSurfaceDistanceCylinder( &
        hit, Lhit, Rstart(:), Ndir(:), &
        surface%center(:), surface%radius )
    end select
  end subroutine calcSurfaceDistance

  ! for Sphere and Cylinder
  subroutine calcSurfaceDistancePolygon( hit, Lhit, Rstart, Ndir, Rsurf, Nsurf, npoly )
    logical, intent(out) :: hit         ! hit or not
    real(8), intent(out) :: Lhit        ! length to the hitting position
    real(8), intent(in)  :: Rstart(3)   ! starting position of the ray
    real(8), intent(in)  :: Ndir(3)     ! normalized direction of the ray
    real(8), intent(in)  :: Rsurf(3,1:npoly+1)  ! positions of polygon vertexes
    real(8), intent(in)  :: Nsurf(3)    ! normal of the polygon
    integer, intent(in)  :: npoly       ! number of vertexes

    integer :: i
    real(8) :: cosangle, Lint
    real(8) :: cross1(3), crossi(3), Ra(3), Rb(3)
    real(8) :: Rhit(3)
    
    Lhit = 0.0d0
    Rhit(:) = 0.0d0

    cosangle = dot_product( Ndir(:), Nsurf(:) )
    if( cosangle <= 0.0d0 ) then  ! wrong direction
      hit = .false.
      ! write(*,*) "parallel to surface"
      ! write(*,*) "cosangle=" , cosangle
      ! write(*,*) "Ndir = ", Ndir
      ! write(*,*) "Nsurf = ", Nsurf
      return
    end if

    Lhit = dot_product( Rsurf(:,1) - Rstart(:), Nsurf(:) )/cosangle

    if( Lhit < 0.0d0 ) then
      hit = .false.
      return
    end if

    Rhit(:) = Rstart(:) + Ndir(:)*Lhit

    do i=1, npoly
      Ra(:) = Rsurf(:,i+1) - Rsurf(:,i)
      Rb(:) = Rhit(:) - Rsurf(:,i)
      crossi(:) = cross_product( Ra(:), Rb(:) )
      if( i == 1 ) then
        cross1(:) = crossi(:)
      else
        if( dot_product(cross1(:),crossi(:))<0.0d0 ) then
          hit = .false.
          return
        end if
      end if
    end do ! i

    !Lhit = Lint
    hit = .true.
  end subroutine calcSurfaceDistancePolygon
  
  logical function isInsidePolygon(Rsurf, npoly, Rhit)
    real(8), intent(in)  :: Rsurf(3,1:npoly+1)  ! positions of polygon vertexes
    integer, intent(in)  :: npoly       ! number of vertexes
    real(8), intent(in)  :: Rhit(3)
    
    real(8) :: cross1(3), crossi(3), Ra(3), Rb(3)
    integer   :: i
    
    isInsidePolygon = .true. 
    do i=1, npoly
      Ra(:) = Rsurf(:,i+1) - Rsurf(:,i)
      Rb(:) = Rhit(:) - Rsurf(:,i)
      crossi(:) = cross_product( Ra(:), Rb(:) )
      if( i == 1 ) then
        cross1(:) = crossi(:)
      else
        if( dot_product(cross1(:),crossi(:))<0.0d0 ) then
          isInsidePolygon = .false.
          return
        end if
      end if
    end do ! i
  end function isInsidePolygon

  subroutine calcSurfaceDistanceSphere( hit, Lhit, Rstart, Ndir, Rcenter, radius )
    logical, intent(out) :: hit          ! hit or not
    real(8), intent(out) :: Lhit         ! length to the hitting position
    real(8), intent(in)  :: Rstart(3)    ! starting position of the ray
    real(8), intent(in)  :: Ndir(3)      ! normalized direction of the ray
    real(8), intent(in)  :: Rcenter(3)   ! center of the sphere
    real(8), intent(in)  :: radius       ! radius of the sphere

    real(8) :: Rsc(3), RD, RR, D, Lhit1, Lhit2
    real(8) :: Rhit(3)


    Rsc(:) = Rcenter(:) - Rstart(:)
    RD = dot_product(Rsc(:), NDir(:))

    RR = dot_product(Rsc(:),Rsc(:))
    if (RR > radius**2 - 1e-10) then
        if (RD <= 0.0d0) then
          hit = .false.
          return
        end if
    end if

    D = RD**2 + radius**2 - RR
    if( D<0.0d0 ) then
      hit = .false.
      return
    end if

    Lhit1 = RD - sqrt(D)
    Lhit2 = RD + sqrt(D)
    if( Lhit1 > 0.0d0 ) then
      Lhit = Lhit1
    else if( Lhit2 > 0.0d0 ) then ! Lhit1<=0
      Lhit = Lhit2
    else
      hit = .false.
      return
    end if

    Rhit(:) = Rstart(:) + Ndir(:)*Lhit
    if (abs(dot_product(Rhit(:) - Rcenter(:), Rhit(:) - Rcenter(:)) - radius**2) > 1e-8) then
        write(*,*) "hitted position is too far away from the sphere"
        write(*,*) "radius = ", radius
        write(*,*) "Rstart = ", Rstart(:) 
        write(*,*) "Rstart to center =", sqrt(dot_product(Rstart(:) - Rcenter(:), Rstart(:) - Rcenter(:)))
        write(*,*) "NDir(:)" , NDir(:)
        write(*,*) "NDir_norm", dot_product(NDir, NDir)
        write(*,*) "Rhit = ", Rhit
        write(*,*) "Rhit to center = ", sqrt(dot_product(Rhit(:) - Rcenter(:), Rhit(:) - Rcenter(:)) )
        write(*,*) "L_flight = ", Ndir(:)*Lhit
    end if 
    hit = .true.
  end subroutine calcSurfaceDistanceSphere

  subroutine calcSurfaceDistanceCylinder( hit, Lhit, Rstart, Ndir, Rcenter, radius )
    logical, intent(out) :: hit        ! hit or not
    real(8), intent(out) :: Lhit       ! length to the hitting position
    real(8), intent(in)  :: Rstart(3)  ! starting position of the ray
    real(8), intent(in)  :: Ndir(3)    ! normalized direction of the ray
    real(8), intent(in)  :: Rcenter(3) ! center of the cylinder
    real(8), intent(in)  :: radius     ! radius of the cylinder

    real(8) :: Rsc(3), RD, RR, D, Lhit1, Lhit2
    real(8) :: Rhit(3), temp_Dir(3), norm_val

    ! cylindrical axis is parallel to Z axis
    Rsc(3) = 0.0d0
    Rsc(1:2) = Rcenter(1:2) - Rstart(1:2)

    temp_Dir = Ndir
    temp_Dir(3) = 0.0d0
    norm_val = dot_product(temp_Dir(:), temp_Dir(:))
    temp_Dir(:) = temp_Dir(:)/sqrt(norm_val)

    RD = dot_product(Rsc(1:2),temp_Dir(1:2))
    RR = dot_product(Rsc(1:2),Rsc(1:2))

    D = RD**2 + radius**2 - RR
    if( D<0.0d0 ) then
      hit = .false.
      return
    end if

    Lhit1 = RD - sqrt(D)
    Lhit2 = RD + sqrt(D)
    if( Lhit1 > 0.0d0 ) then
      Lhit = Lhit1
    else if( Lhit2 > 0.0d0 ) then ! Lhit1<=0
      Lhit = Lhit2
    else
      hit = .false.
      return
    end if

    Lhit = Lhit / sqrt(Ndir(1) * Ndir(1) + Ndir(2) * Ndir(2))
    Rhit(:) = Rstart(:) + Ndir(:)*Lhit
    if( Rhit(3) < 0.0d0 .or. param%boxz < Rhit(3) ) then
      hit = .false.
      return
    end if

    hit = .true.
  end subroutine calcSurfaceDistanceCylinder

  subroutine calcSurfaceDistanceVirtualInterface(hit, Lhit, Rstart, Ndir)
    logical, intent(out) :: hit        ! hit or not
    real(8), intent(out) :: Lhit       ! length to the hitting position
    real(8), intent(in)  :: Rstart(3)  ! starting position of the ray
    real(8), intent(in)  :: Ndir(3)    ! normalized direction of the ray
    real(8) :: Rhit(3)

    if ((Rstart(1) - param%boxx ) * Ndir(1)  >= 0.0d0) then
      hit = .false.
      return
    end if

    Lhit =  abs(  (Rstart(1) - param%boxx) / Ndir(1) )
    Rhit(:) = Rstart(:) + Ndir(:)*Lhit

    if( Rhit(2) < 0.0d0 .or. param%boxy < Rhit(2) ) then
      hit = .false.
      return
    end if

    if( Rhit(3) < 0.0d0 .or. param%boxz < Rhit(3) ) then
      hit = .false.
      return
    end if

    hit = .true.
  end subroutine calcSurfaceDistanceVirtualInterface


  subroutine calcSurfaceNormal( Nsurf, Rsurf, surface )
    real(8), intent(out) :: Nsurf(3) ! normal at the hitting position
    real(8), intent(in) :: Rsurf(3) ! the hitting position
    type(Surface_type), intent(in) :: surface

    select case( surface%kindShape )
    case(ID_SHAPE_POLYGON)
      call calcSurfaceNormalPolygon( Nsurf(:), surface%normal(:) )
    case(ID_SHAPE_SPHERE)
      call calcSurfaceNormalSphere( Nsurf(:), Rsurf(:), surface%center(:), surface%radius )
    case(ID_SHAPE_CYLINDER)
      call calcSurfaceNormalCylinder( Nsurf(:), Rsurf(:), surface%center(:), surface%radius )
    end select
  end subroutine calcSurfaceNormal

  subroutine calcSurfaceNormalPolygon( Nsurf, Nsurf_polygon )
    real(8), intent(out) :: Nsurf(3) ! normal at the hitting position
    real(8), intent(in)  :: Nsurf_polygon(3) ! normal of the polygon

    Nsurf(:) = Nsurf_polygon(:)
  end subroutine calcSurfaceNormalPolygon

  subroutine calcSurfaceNormalVirtualInterface( Nsurf)
    real(8), intent(out) :: Nsurf(3) ! normal at the hitting position
    Nsurf(:) = [ 1.0d0, 0.0d0, 0.0d0 ]
  end subroutine calcSurfaceNormalVirtualInterface

  subroutine calcSurfaceNormalSphere( Nsurf, Rsurf, Rcenter, radius)
    real(8), intent(out) :: Nsurf(3) ! normal at the hitting poistion
    real(8), intent(in) :: Rsurf(3) ! the hitting position
    real(8), intent(in) :: Rcenter(3) ! center of the sphere
    real(8), intent(in) :: radius

    real(8) :: Rcs(3)

    Rcs(:) = Rsurf(:) - Rcenter(:)

    Nsurf(:) = Rcs(:)*(1.0d0/sqrt(dot_product(Rcs,Rcs)))
    !----  Outside the sphere
    if (dot_product(Rcs,Rcs) > radius**2) then 
        Nsurf = -1.0d0 *Nsurf
    end if 
  end subroutine calcSurfaceNormalSphere

  subroutine calcSurfaceNormalCylinder( Nsurf, Rsurf, Rcenter, radius )
    real(8), intent(out) :: Nsurf(3) ! normal at the hitting position
    real(8), intent(in) :: Rsurf(3) ! the hitting position
    real(8), intent(in) :: Rcenter(3) ! center of the cylinder
    real(8), intent(in) :: radius

    real(8) :: Rcs(3)

    ! cylindrical axis is parallel to Z axis
    Rcs(3) = 0.0d0
    Rcs(1:2) = Rsurf(1:2) - Rcenter(1:2)

    Nsurf(:) = Rcs(:)*(1.0d0/sqrt(dot_product(Rcs,Rcs)))
    if (dot_product(Rcs,Rcs) > radius**2) then 
        Nsurf = -1.0d0 *Nsurf
    end if
  end subroutine calcSurfaceNormalCylinder

  subroutine calcSurfaceRotation( Nrotation, Rsurf, surface )
    real(8), intent(out) :: Nrotation(3,3) ! rotation axes at the hitting position
    real(8), intent(in) :: Rsurf(3) ! the hitting position
    type(Surface_type), intent(in) :: surface

    select case( surface%kindShape )
    case(ID_SHAPE_POLYGON)
      call calcSurfaceRotationPolygon( Nrotation(:,:), &
        surface%tangent(:,:) )
    case(ID_SHAPE_SPHERE)
      call calcSurfaceRotationSphere( Nrotation(:,:), Rsurf(:), &
        surface%center(:) )
    case(ID_SHAPE_CYLINDER)
      call calcSurfaceRotationCylinder( Nrotation(:,:), Rsurf(:), &
        surface%center(:) )
    end select
  end subroutine calcSurfaceRotation

  subroutine calcSurfaceRotationPolygon( Nrotation, Nrotation_polygon )
    real(8), intent(out) :: Nrotation(3,3) ! rotation axes at the hitting position
    real(8), intent(in)  :: Nrotation_polygon(3,3)

    Nrotation(:,:) = Nrotation_polygon(:,:)
  end subroutine calcSurfaceRotationPolygon

  subroutine calcSurfaceRotationVirtualInterface( Nrotation)
    real(8), intent(out) :: Nrotation(3,3) ! rotation axes at the hitting position

    Nrotation(:,1) = [ 0.0d0, 0.0d0, 1.0d0 ]
    Nrotation(:,2) = [ 0.0d0, 1.0d0, 0.0d0 ]
    Nrotation(:,3) = [ 1.0d0, 0.0d0, 0.0d0 ]
  end subroutine calcSurfaceRotationVirtualInterface

  subroutine calcSurfaceRotationSphere( Nrotation, Rsurf, Rcenter )
    real(8), intent(out) :: Nrotation(3,3) ! rotation axes at the hitting position
    real(8), intent(in) :: Rsurf(3) ! the hitting position
    real(8), intent(in) :: Rcenter(3) ! center of the sphere

    real(8) :: Rcs(3)

    Rcs(:) = Rsurf(:) - Rcenter(:)
    Rcs = Rcs*(1.0d0/sqrt(dot_product(Rcs,Rcs)))
    Nrotation(:,3) = Rcs 
    if( abs(Rcs(3))==1.0d0 ) then
      Nrotation(:,2) = [ 0.0d0, 1.0d0, 0.0d0 ] 
      Nrotation(:,1) = [ 1.0d0, 0.0d0, 0.0d0 ]
    else if( Rcs(1)**2==0.0d0 .and. &
      Rcs(2)**2==0.0d0 ) then
      Nrotation(:,2) = [ 0.0d0, 1.0d0, 0.0d0 ]
      Nrotation(:,1) = [ 1.0d0, 0.0d0, 0.0d0 ]
    else
      Nrotation(:,2) &
        = [ &
        -Rcs(2) &
        / sqrt(Rcs(1)**2 + Rcs(2)**2), &
        Rcs(1) &
        / sqrt(Rcs(1)**2 + Rcs(2)**2), &
        0.0d0 ]

      Nrotation(:,1) &
        = cross_product( Nrotation(:,2), Nrotation(:,3) )
    end if

    ! Nrotation(:,3) = Rcs(:)*(1.0d0/sqrt(dot_product(Rcs,Rcs)))


    ! ! suppose Nsurf is not parallel to X axis
    ! Nrotation(:,1) = cross_product( Nrotation(:,3), [1.0d0, 0.0d0, 0.0d0] )
    ! if( dot_product(Nrotation(:,1),Nrotation(:,1)) > 0.0d0 ) then ! OK
    !   Nrotation(:,2) = cross_product( Nrotation(:,3), Nrotation(:,1) )
    !   return
    ! end if

    ! ! suppose Nsurf is not parallel to Y axis
    ! Nrotation(:,1) = cross_product( Nrotation(:,3), [0.0d0, 1.0d0, 0.0d0] )
    ! if( dot_product(Nrotation(:,1),Nrotation(:,1)) > 0.0d0 ) then ! OK
    !   Nrotation(:,2) = cross_product( Nrotation(:,3), Nrotation(:,1) )
    !   return
    ! end if

    ! ! suppose Nsurf is not parallel to Z axis
    ! Nrotation(:,1) = cross_product( Nrotation(:,3), [0.0d0, 0.0d0, 1.0d0] )
    ! if( dot_product(Nrotation(:,1),Nrotation(:,1)) > 0.0d0 ) then ! OK
    !   Nrotation(:,2) = cross_product( Nrotation(:,3), Nrotation(:,1) )
    !   return
    ! end if
  end subroutine calcSurfaceRotationSphere

  subroutine calcSurfaceRotationCylinder( Nrotation, Rsurf, Rcenter )
    real(8), intent(out) :: Nrotation(3,3) ! rotation axes at the hitting position
    real(8), intent(in) :: Rsurf(3) ! the hitting position
    real(8), intent(in) :: Rcenter(3) ! center of the cylinder

    real(8) :: Rcs(3), Nnormal(3)

    Rcs(:) = Rsurf(:) - Rcenter(:)

    Nnormal(:) = Rcs(:)
    Nnormal(3) =  0.0d0
    Nnormal(:) = Nnormal(:) * (1.0d0/sqrt(dot_product(Nnormal,Nnormal)))

    Nrotation(:,3) = Nnormal(:)
    if( abs(Nnormal(1))==1.0d0 ) then
      Nrotation(:,2) = [ 0.0d0, 1.0d0, 0.0d0 ]
      Nrotation(:,1) = [ 1.0d0, 0.0d0, 0.0d0 ]
    else if( Nnormal(1)**2==0.0d0 .and. &
      Nnormal(2)**2==0.0d0 ) then
      Nrotation(:,2) = [ 0.0d0, 1.0d0, 0.0d0 ]
      Nrotation(:,1) = [ 1.0d0, 0.0d0, 0.0d0 ]
    else
      Nrotation(:,2) &
        = [ &
        -Nnormal(2) &
        / sqrt(Nnormal(1)**2+Nnormal(2)**2), &
        Nnormal(1) &
        / sqrt(Nnormal(1)**2+Nnormal(2)**2), &
        0.0d0 ]

      Nrotation(:,1) = cross_product( Nrotation(:,2), Nrotation(:,3) )
    end if
  end subroutine calcSurfaceRotationCylinder

  subroutine print_log()
    write(fd_stdout, *)"  _____     _______ _____            _   _  _____    "
    write(fd_stdout, *)"  |  __ \   |__   __|  __ \     /\   | \ | |/ ____|  "
    write(fd_stdout, *)"  | |__) |_____| |  | |__) |   /  \  |  \| | (___    "
    write(fd_stdout, *)"  |  ___/______| |  |  _  /   / /\ \ | . ` |\___ \   "
    write(fd_stdout, *)"  | |          | |  | | \ \  / ____ \| |\  |____) |  "
    write(fd_stdout, *)"  |_|          |_|  |_|  \_\/_/    \_\_| \_|_____/   "
    write(fd_stdout,'("              C.Shao @ Shiomi-lab, UTokyo, ver. 1.0.1")')
    write(fd_stdout,*) "                                                           "
  end subroutine print_log

  subroutine init_random_seed()
    ! random_seed based on system time
    integer :: i, n, clock
    integer, dimension(:), allocatable :: seed

    call random_seed(size = n)
    allocate(seed(n))

    call system_clock(count=clock)

    !!$    seed = 0 ! for DEBUG
    seed = clock + 37 * (/ (i - 1, i = 1, n) /)
    write(fd_stdout,*) "seed =", seed
    call random_seed(put = seed)

    deallocate(seed)
  end subroutine init_random_seed

end module RTMC

subroutine RayTrace( GeoFile, status )
  use RTMC
  implicit none

  character(len=*), intent(in) :: GeoFile
  integer, intent(out) :: status

  integer :: s
  real(8), allocatable :: SEMFP(:,:,:)

  if (verbose >= 1)  then
    ! write to screen
    fd_stdout = 6
  else
    fd_stdout = 106
    open(unit=fd_stdout,file="standard.out",status="replace")
  end if

  if( mpi_rank == 0 ) then
    call print_log()
  end if

  call init_random_seed()
  !call random_seed()
  call flush(fd_stdout)
  status = 0
  try_block: do ! try block, it is not a loop
    call loadMCP( GeoFile, status )
    if( status /= 0 ) exit try_block

    if( mpi_rank == 0 ) then
      write(fd_stdout,*)
      write(fd_stdout,*) 'MeanFreePathes: ', param%nMeanFreePath
      write(fd_stdout,*) 'KPoints:        ', anphon%nkpoint
      write(fd_stdout,*) 'Branches:       ', anphon%nbranch
      write(fd_stdout,*) 'Angles:         ', param%nangle-1
      write(fd_stdout,*) 'Particles:      ', param%nphonon
    end if

    interval_trace = param%nMeanFreePath &
      *anphon%nkpoint*anphon%nbranch*(param%nangle-1) &
      *param%nphonon / (particle_trace*mpi_size)

    if( fd_trace>0 .and. mpi_rank == 0 ) then
      open(fd_trace,file="trace.out")
      write(fd_trace,"('# Datatype = trace')")
      write(fd_trace,"('# particle_trace = ',i8)") particle_trace
      write(fd_trace,"('# interval_trace = ',i8)") interval_trace
      write(fd_trace,*)
    end if

    allocate( SEMFP(anphon%nkpoint,anphon%nbranch,param%nMeanFreePath) )

    do s=1, param%nMeanFreePath
      param%MeanFreePath_current = param%vMeanFreePath(s)

      if( mpi_rank == 0 ) then
        write(fd_stdout,'(" ")')
        write(fd_stdout,'(" Run:", i4, ",  Total: ", i4)') s, param%nMeanFreePath
        write(fd_stdout,"(a,es10.3,a)") " Current MFP: ", param%MeanFreePath_current, " nm"
        call flush(fd_stdout)
      end if

      call TransmissionFunction
      call effectiveMFP( SEMFP(:,:,s) )
    end do

    if( mpi_rank == 0 ) then
      call fileoutput( SEMFP, GeoFile )
    end if
#ifdef MPI
    call MPI_Barrier( MPI_COMM_WORLD, mpi_ierr );
#endif

    deallocate( SEMFP )

    exit try_block
  end do try_block

  if( fd_trace>0 .and. mpi_rank == 0 ) then
    close(fd_trace)
  end if

  ! catch block
  if( status /= 0 ) then
    write(fd_stdout,*) "Error: RayTrace"
  end if

  call clearGeometry
  call clearArray
  call flush(fd_stdout)
  if(fd_stdout /= 6 ) close(fd_stdout)

  ! invert status flag
  if( status == 0 ) then
    status = 1 ! True means success
  else
    status = 0 ! False means error
  end if
end subroutine RayTrace

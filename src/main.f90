program main
  use RTMC
  implicit none
  character(len=256) :: GeoFile
  logical :: status

#ifdef MPI
  call MPI_Init(mpi_ierr)
  call MPI_Comm_rank( MPI_COMM_WORLD, mpi_rank, mpi_ierr )
  call MPI_Comm_size( MPI_COMM_WORLD, mpi_size, mpi_ierr )
#endif
  call getarg(1,GeoFile)

  verbose = 1 ! overwrite default value 0

  call initRand
  call RayTrace(GeoFile,status)

#ifdef MPI
  call MPI_Finalize(mpi_ierr)
#endif
endprogram main

subroutine initRand
  integer :: sizeseed
  integer, allocatable :: seed(:)

  call random_seed(size=sizeseed)
  allocate( seed(sizeseed) )
  seed(:) = 0
  call random_seed(put=seed)
  deallocate( seed )
end subroutine initRand

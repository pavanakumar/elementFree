program readGrid
  use connectivity
  use boundary
  implicit none
  character (len=*), parameter :: gridFile = "blcone.dat"
  character (len=*), parameter :: normalFile = "blcone_normal.dat"
  integer :: num_bcs, nngbr, bc_flag, iorder
  integer :: node_id, i, bc_node_id
  real (kind=8) :: tgt(3), norm(3)
  integer :: wall_count = 0, ext_count = 0, out_count = 0, sym_count = 0
  integer :: in_count = 0

  open( unit=111, file=gridFile, status='old' )
  open( unit=222, file=normalFile, status='old' )

  read(111,*) my_size, num_bcs
  allocate( xyz(3,size), xadj(size + 1) )

  !!! Read nodes from blcone and check their bc_flag 
  !!! and increment counter to get my_sizes
  xadj(1) = 1
  adj_size = 0
  do i = 1, my_size
    read(111,*) node_id, xyz(1:3,i), nngbr, bc_flag, iorder
    read(111,*) ! ignore neighbours
    xadj(i+1) = xadj(i) + nngbr
    adj_size = adj_size + nngbr
    if( bc_flag .ne. 0 ) then
      read(222,*) bc_node_id! ignore tangents and normals
      if( bc_node_id .ne. node_id ) then
        write(*,*) 'Make sure the blcone_normal.dat is sorted'
        write(*,*) 'mv blcone_normal.dat blcone_normal.old'
        write(*,*) 'sort -n blcone_normal.old > blcone_normal.dat'
        stop
      end if
      !!! Count wall bc points
      if( bc_flag.eq.1 ) wall_sz = wall_sz + 1
    end if
  end do

  ! re-read file
  allocate( adjncy(adj_size) )
  call AllocateBoundary

  rewind(111)
  rewind(222)

  read(111,*) ! ignore my_size header
  do i = 1, my_size
    read(111,*) node_id, xyz(1:3,i), nngbr, bc_flag, iorder
    read(111,*) adjncy( xadj(i) : xadj(i+1) - 1 )
     if( bc_flag .ne. 0 ) then
      read(222,*) bc_node_id, tgt(1:3), norm(1:3)
      !!! Wall boundary
      if( bc_flag.eq.1 ) then
        wall_count = wall_count + 1
        wall_pts(wall_count) = node_id + 1
        wall_nxyz(:,wall_count) = norm
      end if
    end if
  end do
  connectivityAllocated = .true.

  do i = 1, adj_size
    adjncy(i) = adjncy(i) + 1
  end do

  close(111)
  close(222)
end program readGrid


module Task
  use :: mpi
  implicit none
  contains
 
  subroutine GetMaxCoordinates(A, x1, y1, x2, y2)
    implicit none
    real(8), intent(in), dimension(:,:) :: A
    real(8), allocatable :: current_column(:)
    real(8) :: current_sum, max_sum, send_max_sum
    integer(4), dimension(4) :: send_max_coords
    integer(4), intent(out) :: x1, y1, x2, y2
    integer(4) :: i, n, L, R, Up, Down, m, tmp
    integer(4) :: mpiErr, mpiSize, mpiRank
    integer(4), dimension(MPI_STATUS_SIZE) :: status

    call mpi_comm_size(MPI_COMM_WORLD, mpiSize, mpiErr)
    call mpi_comm_rank(MPI_COMM_WORLD, mpiRank, mpiErr)
    
    m = size(A, dim=1) 
    n = size(A, dim=2) 
 
    allocate(current_column(m))

    x1=1
    y1=1
    x2=1
    y2=1
    max_sum = A(1,1)

    do L = mpiRank+1, n, mpiSize
      current_column = A(:, L)

      do R = L, n
        if (R > L) then
          current_column = current_column + A(:, R)
        endif 
         
        call FindMaxInArray(current_column, current_sum, Up, Down)

        if (current_sum > max_sum) then
          max_sum = current_sum
          x1 = Up
          x2 = Down
          y1 = L
          y2 = R
        endif
      end do
    end do

    if (mpiRank /= 0) then
      send_max_sum = max_sum
      send_max_coords(1) = x1
      send_max_coords(2) = x2
      send_max_coords(3) = y1
      send_max_coords(4) = y2
      call mpi_send(   send_max_sum, 1, MPI_REAL8,    0, 777+mpiRank, MPI_COMM_WORLD, mpiErr)
      call mpi_send(send_max_coords, 4, MPI_INTEGER4, 0, 444+mpiRank, MPI_COMM_WORLD, mpiErr)
    endif

    if (mpiRank == 0) then
      do i=1, mpiSize-1
        call mpi_recv(   send_max_sum, 1, MPI_REAL8,    i, 777+i, MPI_COMM_WORLD, status, mpiErr)
        call mpi_recv(send_max_coords, 4, MPI_INTEGER4, i, 444+i, MPI_COMM_WORLD, status, mpiErr)

        if (send_max_sum > max_sum) then
          max_sum = send_max_sum
          x1 = send_max_coords(1)
          x2 = send_max_coords(2)
          y1 = send_max_coords(3)
          y2 = send_max_coords(4)
        endif
      enddo

      send_max_sum = max_sum
      send_max_coords(1) = x1
      send_max_coords(2) = x2
      send_max_coords(3) = y1
      send_max_coords(4) = y2
    endif

    call mpi_bcast(send_max_sum,    1, MPI_REAL8,    0, MPI_COMM_WORLD, mpiErr)
    call mpi_bcast(send_max_coords, 4, MPI_INTEGER4, 0, MPI_COMM_WORLD, mpiErr)

    max_sum = send_max_sum
    x1 = send_max_coords(1)
    x2 = send_max_coords(2)
    y1 = send_max_coords(3)
    y2 = send_max_coords(4)

   deallocate(current_column)
  end subroutine

  subroutine FindMaxInArray(A, Summ, Up, Down)
    implicit none
    real(8), intent(in), dimension(:) :: A
    integer(4), intent(out) :: Up, Down
    real(8), intent(out) :: Summ
    real(8) :: cur_sum
    integer(4) :: minus_pos, i

    Summ = A(1)
    Up = 1
    Down = 1
    cur_sum = 0
    minus_pos = 0

    do i=1, size(A)
      cur_sum = cur_sum + A(i)
      if (cur_sum > Summ) then
        Summ = cur_sum
        Up = minus_pos + 1
        Down = i
      endif
     
      if (cur_sum < 0) then
        cur_sum = 0
        minus_pos = i
      endif
    enddo

  end subroutine
end module

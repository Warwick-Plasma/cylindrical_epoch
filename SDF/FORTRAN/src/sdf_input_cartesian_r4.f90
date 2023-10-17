!
! SDF (Self-Describing Format) Fortran Library
! Copyright (c) 2011-2016, SDF Development Team
!
! Distributed under the terms of the BSD 3-clause License.
! See the LICENSE file for details.
!

MODULE sdf_input_cartesian_r4

  USE sdf_input_cartesian_ru

  IMPLICIT NONE

CONTAINS

  ! Mesh loading functions

  SUBROUTINE read_plain_mesh_info_r4(h, geometry, dims, extents, dim_labels, &
      dim_units, dim_mults)

    TYPE(sdf_file_handle) :: h
    INTEGER, INTENT(OUT) :: geometry
    INTEGER, DIMENSION(:), INTENT(OUT) :: dims
    REAL(r4), DIMENSION(:), INTENT(OUT) :: extents
    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL :: dim_labels(:), dim_units(:)
    REAL(r4), DIMENSION(:), INTENT(OUT), OPTIONAL :: dim_mults
    INTEGER :: i, clen
    TYPE(sdf_block_type), POINTER :: b

    CALL read_plain_mesh_info_ru(h, geometry, dims)
    b => h%current_block

    extents(1:2*b%ndims) = REAL(b%extents(1:2*b%ndims),r4)
    IF (PRESENT(dim_mults)) dim_mults = REAL(b%dim_mults,r4)
    IF (PRESENT(dim_labels)) THEN
      DO i = 1,b%ndims
        clen = MIN(LEN(dim_labels(i)),INT(c_id_length))
        dim_labels(i)(1:clen) = b%dim_labels(i)(1:clen)
      END DO
    END IF
    IF (PRESENT(dim_units)) THEN
      DO i = 1,b%ndims
        clen = MIN(LEN(dim_units(i)),INT(c_id_length))
        dim_units(i)(1:clen) = b%dim_units(i)(1:clen)
      END DO
    END IF

  END SUBROUTINE read_plain_mesh_info_r4



  SUBROUTINE read_srl_1d_mesh_r4(h, x)

    TYPE(sdf_file_handle) :: h
    REAL(r4), DIMENSION(:), INTENT(OUT) :: x
    INTEGER :: errcode, intn
    TYPE(sdf_block_type), POINTER :: b

    IF (sdf_check_block_header(h)) RETURN

    b => h%current_block
    IF (.NOT. b%done_info) CALL read_plain_mesh_info_ru(h)

    h%current_location = b%data_location

    intn = b%dims(1)
    IF (h%rank == h%rank_master) THEN
      CALL MPI_FILE_SEEK(h%filehandle, h%current_location, MPI_SEEK_SET, &
          errcode)
      CALL MPI_FILE_READ(h%filehandle, x, intn, b%mpitype, &
          MPI_STATUS_IGNORE, errcode)
    END IF

    CALL MPI_BCAST(x, intn, b%mpitype, h%rank_master, h%comm, errcode)

    ! That should be it, so now skip to end of block
    h%current_location = b%next_block_location
    b%done_data = .TRUE.

  END SUBROUTINE read_srl_1d_mesh_r4



  SUBROUTINE read_srl_2d_mesh_r4(h, x, y)

    TYPE(sdf_file_handle) :: h
    REAL(r4), DIMENSION(:), INTENT(OUT) :: x, y
    INTEGER :: errcode, intn
    TYPE(sdf_block_type), POINTER :: b

    IF (sdf_check_block_header(h)) RETURN

    b => h%current_block
    IF (.NOT. b%done_info) CALL read_plain_mesh_info_ru(h)

    h%current_location = b%data_location

    intn = b%dims(1)
    IF (h%rank == h%rank_master) THEN
      CALL MPI_FILE_SEEK(h%filehandle, h%current_location, MPI_SEEK_SET, &
          errcode)
      CALL MPI_FILE_READ(h%filehandle, x, intn, b%mpitype, &
          MPI_STATUS_IGNORE, errcode)
    END IF

    CALL MPI_BCAST(x, intn, b%mpitype, h%rank_master, h%comm, errcode)

    intn = b%dims(2)
    IF (h%rank == h%rank_master) THEN
      CALL MPI_FILE_READ(h%filehandle, y, intn, b%mpitype, &
          MPI_STATUS_IGNORE, errcode)
    END IF

    CALL MPI_BCAST(y, intn, b%mpitype, h%rank_master, h%comm, errcode)

    ! That should be it, so now skip to end of block
    h%current_location = b%next_block_location
    b%done_data = .TRUE.

  END SUBROUTINE read_srl_2d_mesh_r4



  SUBROUTINE read_srl_3d_mesh_r4(h, x, y, z)

    TYPE(sdf_file_handle) :: h
    REAL(r4), DIMENSION(:), INTENT(OUT) :: x, y, z
    INTEGER :: errcode, intn
    TYPE(sdf_block_type), POINTER :: b

    IF (sdf_check_block_header(h)) RETURN

    b => h%current_block
    IF (.NOT. b%done_info) CALL read_plain_mesh_info_ru(h)

    h%current_location = b%data_location

    intn = b%dims(1)
    IF (h%rank == h%rank_master) THEN
      CALL MPI_FILE_SEEK(h%filehandle, h%current_location, MPI_SEEK_SET, &
          errcode)
      CALL MPI_FILE_READ(h%filehandle, x, intn, b%mpitype, &
          MPI_STATUS_IGNORE, errcode)
    END IF

    CALL MPI_BCAST(x, intn, b%mpitype, h%rank_master, h%comm, errcode)

    intn = b%dims(2)
    IF (h%rank == h%rank_master) THEN
      CALL MPI_FILE_READ(h%filehandle, y, intn, b%mpitype, &
          MPI_STATUS_IGNORE, errcode)
    END IF

    CALL MPI_BCAST(y, intn, b%mpitype, h%rank_master, h%comm, errcode)

    intn = b%dims(3)
    IF (h%rank == h%rank_master) THEN
      CALL MPI_FILE_READ(h%filehandle, z, intn, b%mpitype, &
          MPI_STATUS_IGNORE, errcode)
    END IF

    CALL MPI_BCAST(z, intn, b%mpitype, h%rank_master, h%comm, errcode)

    ! That should be it, so now skip to end of block
    h%current_location = b%next_block_location
    b%done_data = .TRUE.

  END SUBROUTINE read_srl_3d_mesh_r4



  SUBROUTINE read_1d_mesh_r4(h, x, distribution, subarray)

    TYPE(sdf_file_handle) :: h
    REAL(r4), DIMENSION(:), INTENT(OUT) :: x
    INTEGER, INTENT(IN) :: distribution(:)
    INTEGER, INTENT(IN), OPTIONAL :: subarray(:)
    INTEGER :: errcode
    TYPE(sdf_block_type), POINTER :: b

    IF (sdf_check_block_header(h)) RETURN

    b => h%current_block
    IF (.NOT. b%done_info) CALL read_plain_mesh_info_ru(h)

    ! Read the actual data

    h%current_location = b%data_location

    CALL MPI_FILE_SET_VIEW(h%filehandle, h%current_location, MPI_BYTE, &
        distribution(1), 'native', MPI_INFO_NULL, errcode)

    IF (PRESENT(subarray)) THEN
      CALL MPI_FILE_READ_ALL(h%filehandle, x, 1, subarray(1), &
          MPI_STATUS_IGNORE, errcode)
    ELSE
      CALL MPI_FILE_READ_ALL(h%filehandle, x, SIZE(x), b%mpitype, &
          MPI_STATUS_IGNORE, errcode)
    END IF

    CALL MPI_FILE_SET_VIEW(h%filehandle, c_off0, MPI_BYTE, MPI_BYTE, 'native', &
        MPI_INFO_NULL, errcode)

    ! That should be it, so now skip to end of block
    h%current_location = b%next_block_location
    b%done_data = .TRUE.

  END SUBROUTINE read_1d_mesh_r4



  SUBROUTINE read_2d_mesh_r4(h, x, y, distribution, subarray)

    TYPE(sdf_file_handle) :: h
    REAL(r4), DIMENSION(:), INTENT(OUT) :: x, y
    INTEGER, INTENT(IN) :: distribution(:)
    INTEGER, INTENT(IN), OPTIONAL :: subarray(:)
    INTEGER :: errcode
    TYPE(sdf_block_type), POINTER :: b

    IF (sdf_check_block_header(h)) RETURN

    b => h%current_block
    IF (.NOT. b%done_info) CALL read_plain_mesh_info_ru(h)

    ! Read the actual data

    h%current_location = b%data_location

    CALL MPI_FILE_SET_VIEW(h%filehandle, h%current_location, MPI_BYTE, &
        distribution(1), 'native', MPI_INFO_NULL, errcode)

    IF (PRESENT(subarray)) THEN
      CALL MPI_FILE_READ_ALL(h%filehandle, x, 1, subarray(1), &
          MPI_STATUS_IGNORE, errcode)
    ELSE
      CALL MPI_FILE_READ_ALL(h%filehandle, x, SIZE(x), b%mpitype, &
          MPI_STATUS_IGNORE, errcode)
    END IF

    h%current_location = h%current_location + b%dims(1) * b%type_size

    CALL MPI_FILE_SET_VIEW(h%filehandle, h%current_location, MPI_BYTE, &
        distribution(2), 'native', MPI_INFO_NULL, errcode)

    IF (PRESENT(subarray)) THEN
      CALL MPI_FILE_READ_ALL(h%filehandle, y, 1, subarray(2), &
          MPI_STATUS_IGNORE, errcode)
    ELSE
      CALL MPI_FILE_READ_ALL(h%filehandle, y, SIZE(y), b%mpitype, &
          MPI_STATUS_IGNORE, errcode)
    END IF

    CALL MPI_FILE_SET_VIEW(h%filehandle, c_off0, MPI_BYTE, MPI_BYTE, 'native', &
        MPI_INFO_NULL, errcode)

    ! That should be it, so now skip to end of block
    h%current_location = b%next_block_location
    b%done_data = .TRUE.

  END SUBROUTINE read_2d_mesh_r4



  SUBROUTINE read_3d_mesh_r4(h, x, y, z, distribution, subarray)

    TYPE(sdf_file_handle) :: h
    REAL(r4), DIMENSION(:), INTENT(OUT) :: x, y, z
    INTEGER, INTENT(IN) :: distribution(:)
    INTEGER, INTENT(IN), OPTIONAL :: subarray(:)
    INTEGER :: errcode
    TYPE(sdf_block_type), POINTER :: b

    IF (sdf_check_block_header(h)) RETURN

    b => h%current_block
    IF (.NOT. b%done_info) CALL read_plain_mesh_info_ru(h)

    ! Read the actual data

    h%current_location = b%data_location

    CALL MPI_FILE_SET_VIEW(h%filehandle, h%current_location, MPI_BYTE, &
        distribution(1), 'native', MPI_INFO_NULL, errcode)

    IF (PRESENT(subarray)) THEN
      CALL MPI_FILE_READ_ALL(h%filehandle, x, 1, subarray(1), &
          MPI_STATUS_IGNORE, errcode)
    ELSE
      CALL MPI_FILE_READ_ALL(h%filehandle, x, SIZE(x), b%mpitype, &
          MPI_STATUS_IGNORE, errcode)
    END IF

    h%current_location = h%current_location + b%dims(1) * b%type_size

    CALL MPI_FILE_SET_VIEW(h%filehandle, h%current_location, MPI_BYTE, &
        distribution(2), 'native', MPI_INFO_NULL, errcode)

    IF (PRESENT(subarray)) THEN
      CALL MPI_FILE_READ_ALL(h%filehandle, y, 1, subarray(2), &
          MPI_STATUS_IGNORE, errcode)
    ELSE
      CALL MPI_FILE_READ_ALL(h%filehandle, y, SIZE(y), b%mpitype, &
          MPI_STATUS_IGNORE, errcode)
    END IF

    h%current_location = h%current_location + b%dims(2) * b%type_size

    CALL MPI_FILE_SET_VIEW(h%filehandle, h%current_location, MPI_BYTE, &
        distribution(3), 'native', MPI_INFO_NULL, errcode)

    IF (PRESENT(subarray)) THEN
      CALL MPI_FILE_READ_ALL(h%filehandle, z, 1, subarray(3), &
          MPI_STATUS_IGNORE, errcode)
    ELSE
      CALL MPI_FILE_READ_ALL(h%filehandle, z, SIZE(z), b%mpitype, &
          MPI_STATUS_IGNORE, errcode)
    END IF

    CALL MPI_FILE_SET_VIEW(h%filehandle, c_off0, MPI_BYTE, MPI_BYTE, 'native', &
        MPI_INFO_NULL, errcode)

    ! That should be it, so now skip to end of block
    h%current_location = b%next_block_location
    b%done_data = .TRUE.

  END SUBROUTINE read_3d_mesh_r4



  SUBROUTINE read_1d_lag_mesh_r4(h, x, distribution, subarray)

    TYPE(sdf_file_handle) :: h
    REAL(r4), DIMENSION(:), INTENT(OUT) :: x
    INTEGER, INTENT(IN) :: distribution
    INTEGER, INTENT(IN), OPTIONAL :: subarray
    INTEGER :: errcode
    TYPE(sdf_block_type), POINTER :: b

    IF (sdf_check_block_header(h)) RETURN

    b => h%current_block
    IF (.NOT. b%done_info) CALL read_plain_mesh_info_ru(h)

    ! Read the actual data

    h%current_location = b%data_location

    CALL MPI_FILE_SET_VIEW(h%filehandle, h%current_location, MPI_BYTE, &
        distribution, 'native', MPI_INFO_NULL, errcode)

    IF (PRESENT(subarray)) THEN
      CALL MPI_FILE_READ_ALL(h%filehandle, x, 1, subarray, &
          MPI_STATUS_IGNORE, errcode)
    ELSE
      CALL MPI_FILE_READ_ALL(h%filehandle, x, SIZE(x), b%mpitype, &
          MPI_STATUS_IGNORE, errcode)
    END IF

    CALL MPI_FILE_SET_VIEW(h%filehandle, c_off0, MPI_BYTE, MPI_BYTE, 'native', &
        MPI_INFO_NULL, errcode)

    ! That should be it, so now skip to end of block
    h%current_location = b%next_block_location
    b%done_data = .TRUE.

  END SUBROUTINE read_1d_lag_mesh_r4



  SUBROUTINE read_2d_lag_mesh_r4(h, x, y, distribution, subarray)

    TYPE(sdf_file_handle) :: h
    REAL(r4), DIMENSION(:,:), INTENT(OUT) :: x, y
    INTEGER, INTENT(IN) :: distribution
    INTEGER, INTENT(IN), OPTIONAL :: subarray
    INTEGER :: errcode
    TYPE(sdf_block_type), POINTER :: b

    IF (sdf_check_block_header(h)) RETURN

    b => h%current_block
    IF (.NOT. b%done_info) CALL read_plain_mesh_info_ru(h)

    ! Read the actual data

    h%current_location = b%data_location

    CALL MPI_FILE_SET_VIEW(h%filehandle, h%current_location, MPI_BYTE, &
        distribution, 'native', MPI_INFO_NULL, errcode)

    IF (PRESENT(subarray)) THEN
      CALL MPI_FILE_READ_ALL(h%filehandle, x, 1, subarray, &
          MPI_STATUS_IGNORE, errcode)

      CALL MPI_FILE_READ_ALL(h%filehandle, y, 1, subarray, &
          MPI_STATUS_IGNORE, errcode)
    ELSE
      CALL MPI_FILE_READ_ALL(h%filehandle, x, SIZE(x), b%mpitype, &
          MPI_STATUS_IGNORE, errcode)

      CALL MPI_FILE_READ_ALL(h%filehandle, y, SIZE(y), b%mpitype, &
          MPI_STATUS_IGNORE, errcode)
    END IF

    CALL MPI_FILE_SET_VIEW(h%filehandle, c_off0, MPI_BYTE, MPI_BYTE, 'native', &
        MPI_INFO_NULL, errcode)

    ! That should be it, so now skip to end of block
    h%current_location = b%next_block_location
    b%done_data = .TRUE.

  END SUBROUTINE read_2d_lag_mesh_r4



  SUBROUTINE read_3d_lag_mesh_r4(h, x, y, z, distribution, subarray)

    TYPE(sdf_file_handle) :: h
    REAL(r4), DIMENSION(:,:,:), INTENT(OUT) :: x, y, z
    INTEGER, INTENT(IN) :: distribution
    INTEGER, INTENT(IN), OPTIONAL :: subarray
    INTEGER :: errcode
    TYPE(sdf_block_type), POINTER :: b

    IF (sdf_check_block_header(h)) RETURN

    b => h%current_block
    IF (.NOT. b%done_info) CALL read_plain_mesh_info_ru(h)

    ! Read the actual data

    h%current_location = b%data_location

    CALL MPI_FILE_SET_VIEW(h%filehandle, h%current_location, MPI_BYTE, &
        distribution, 'native', MPI_INFO_NULL, errcode)

    IF (PRESENT(subarray)) THEN
      CALL MPI_FILE_READ_ALL(h%filehandle, x, 1, subarray, &
          MPI_STATUS_IGNORE, errcode)

      CALL MPI_FILE_READ_ALL(h%filehandle, y, 1, subarray, &
          MPI_STATUS_IGNORE, errcode)

      CALL MPI_FILE_READ_ALL(h%filehandle, z, 1, subarray, &
          MPI_STATUS_IGNORE, errcode)
    ELSE
      CALL MPI_FILE_READ_ALL(h%filehandle, x, SIZE(x), b%mpitype, &
          MPI_STATUS_IGNORE, errcode)

      CALL MPI_FILE_READ_ALL(h%filehandle, y, SIZE(y), b%mpitype, &
          MPI_STATUS_IGNORE, errcode)

      CALL MPI_FILE_READ_ALL(h%filehandle, z, SIZE(z), b%mpitype, &
          MPI_STATUS_IGNORE, errcode)
    END IF

    CALL MPI_FILE_SET_VIEW(h%filehandle, c_off0, MPI_BYTE, MPI_BYTE, 'native', &
        MPI_INFO_NULL, errcode)

    ! That should be it, so now skip to end of block
    h%current_location = b%next_block_location
    b%done_data = .TRUE.

  END SUBROUTINE read_3d_lag_mesh_r4



  ! Variable loading functions

  SUBROUTINE read_plain_variable_info_r4(h, dims, units, mesh_id, stagger, mult)

    TYPE(sdf_file_handle) :: h
    INTEGER, DIMENSION(:), INTENT(OUT) :: dims
    CHARACTER(LEN=*), INTENT(OUT) :: units, mesh_id
    INTEGER, INTENT(OUT) :: stagger
    REAL(r4), INTENT(OUT) :: mult
    TYPE(sdf_block_type), POINTER :: b

    CALL read_plain_variable_info_ru(h, dims, units, mesh_id, stagger)
    b => h%current_block
    mult = REAL(b%mult,r4)

  END SUBROUTINE read_plain_variable_info_r4



  !----------------------------------------------------------------------------
  ! Code to read a 1D cartesian variable in parallel
  ! using the mpitype {distribution} for distribution of data
  !----------------------------------------------------------------------------

  SUBROUTINE read_1d_float_r4(h, variable, distribution, subarray)

    TYPE(sdf_file_handle) :: h
    REAL(r4), DIMENSION(:), INTENT(OUT) :: variable
    INTEGER, INTENT(IN) :: distribution, subarray

    CALL read_nd_float_r4(h, variable(1), distribution, subarray)

  END SUBROUTINE read_1d_float_r4



  !----------------------------------------------------------------------------
  ! Code to read a nD cartesian variable in parallel
  ! using the mpitype {distribution} for distribution of data
  !----------------------------------------------------------------------------

  SUBROUTINE read_nd_float_r4(h, variable, distribution, subarray)

    TYPE(sdf_file_handle) :: h
    REAL(r4), INTENT(OUT) :: variable
    INTEGER, INTENT(IN) :: distribution, subarray
    INTEGER :: errcode
    TYPE(sdf_block_type), POINTER :: b

    IF (sdf_check_block_header(h)) RETURN

    b => h%current_block
    IF (.NOT. b%done_info) CALL read_plain_variable_info_ru(h)

    ! Read the actual data

    h%current_location = b%data_location

    CALL MPI_FILE_SET_VIEW(h%filehandle, h%current_location, MPI_BYTE, &
        distribution, 'native', MPI_INFO_NULL, errcode)

    CALL MPI_FILE_READ_ALL(h%filehandle, variable, 1, subarray, &
        MPI_STATUS_IGNORE, errcode)

    CALL MPI_FILE_SET_VIEW(h%filehandle, c_off0, MPI_BYTE, MPI_BYTE, 'native', &
        MPI_INFO_NULL, errcode)

    h%current_location = b%data_location + b%data_length
    b%done_data = .TRUE.

  END SUBROUTINE read_nd_float_r4



  !----------------------------------------------------------------------------
  ! Code to read a 1D cartesian variable component (with material index first)
  ! in parallel using the mpitype {distribution} for distribution of data
  !----------------------------------------------------------------------------

  SUBROUTINE read_1d_var_first_r4(h, nm, dims, variable, idx, distribution, &
      subarray)

    TYPE(sdf_file_handle) :: h
    INTEGER, INTENT(IN) :: nm
    INTEGER(i4), INTENT(IN) :: dims(:)
    REAL(r4), INTENT(OUT) :: variable(:,:)
    INTEGER, INTENT(IN) :: idx, distribution, subarray

    CALL read_nd_float_r4(h, variable(idx,1), distribution, subarray)

  END SUBROUTINE read_1d_var_first_r4



  !----------------------------------------------------------------------------
  ! Code to read a 2D cartesian variable component (with material index first)
  ! in parallel using the mpitype {distribution} for distribution of data
  !----------------------------------------------------------------------------

  SUBROUTINE read_2d_var_first_r4(h, nm, dims, variable, idx, distribution, &
      subarray)

    TYPE(sdf_file_handle) :: h
    INTEGER, INTENT(IN) :: nm
    INTEGER(i4), INTENT(IN) :: dims(:)
    REAL(r4), INTENT(OUT) :: variable(:,:,:)
    INTEGER, INTENT(IN) :: idx, distribution, subarray

    CALL read_nd_float_r4(h, variable(idx,1,1), distribution, subarray)

  END SUBROUTINE read_2d_var_first_r4



  !----------------------------------------------------------------------------
  ! Code to read a 3D cartesian variable component (with material index first)
  ! in parallel using the mpitype {distribution} for distribution of data
  !----------------------------------------------------------------------------

  SUBROUTINE read_3d_var_first_r4(h, nm, dims, variable, idx, distribution, &
      subarray)

    TYPE(sdf_file_handle) :: h
    INTEGER, INTENT(IN) :: nm
    INTEGER(i4), INTENT(IN) :: dims(:)
    REAL(r4), INTENT(OUT) :: variable(:,:,:,:)
    INTEGER, INTENT(IN) :: idx, distribution, subarray

    CALL read_nd_float_r4(h, variable(idx,1,1,1), distribution, subarray)

  END SUBROUTINE read_3d_var_first_r4



  !----------------------------------------------------------------------------
  ! Code to read a 3D cartesian variable component (with material index last)
  ! in parallel using the mpitype {distribution} for distribution of data
  !----------------------------------------------------------------------------

  SUBROUTINE read_1d_var_last_r4(h, nm, dims, variable, idx, distribution, &
      subarray)

    TYPE(sdf_file_handle) :: h
    INTEGER, INTENT(IN) :: nm
    INTEGER(i4), INTENT(IN) :: dims(:)
    REAL(r4), INTENT(OUT) :: variable(:,:)
    INTEGER, INTENT(IN) :: idx, distribution, subarray

    CALL read_nd_float_r4(h, variable(1,idx), distribution, subarray)

  END SUBROUTINE read_1d_var_last_r4



  !----------------------------------------------------------------------------
  ! Code to read a 2D cartesian variable component (with material index last)
  ! in parallel using the mpitype {distribution} for distribution of data
  !----------------------------------------------------------------------------

  SUBROUTINE read_2d_var_last_r4(h, nm, dims, variable, idx, distribution, &
      subarray)

    TYPE(sdf_file_handle) :: h
    INTEGER, INTENT(IN) :: nm
    INTEGER(i4), INTENT(IN) :: dims(:)
    REAL(r4), INTENT(OUT) :: variable(:,:,:)
    INTEGER, INTENT(IN) :: idx, distribution, subarray

    CALL read_nd_float_r4(h, variable(1,1,idx), distribution, subarray)

  END SUBROUTINE read_2d_var_last_r4



  !----------------------------------------------------------------------------
  ! Code to read a 3D cartesian variable component (with material index last)
  ! in parallel using the mpitype {distribution} for distribution of data
  !----------------------------------------------------------------------------

  SUBROUTINE read_3d_var_last_r4(h, nm, dims, variable, idx, distribution, &
      subarray)

    TYPE(sdf_file_handle) :: h
    INTEGER, INTENT(IN) :: nm
    INTEGER(i4), INTENT(IN) :: dims(:)
    REAL(r4), INTENT(OUT) :: variable(:,:,:,:)
    INTEGER, INTENT(IN) :: idx, distribution, subarray

    CALL read_nd_float_r4(h, variable(1,1,1,idx), distribution, subarray)

  END SUBROUTINE read_3d_var_last_r4



  !----------------------------------------------------------------------------
  ! Code to read a 1D cartesian material mesh in parallel using the
  ! mpitype {distribution} for distribution of data
  !----------------------------------------------------------------------------

  SUBROUTINE read_1d_material_r4(h, variable, distribution, subarray, last_in)

    TYPE(sdf_file_handle) :: h
    REAL(r4), DIMENSION(:,:), INTENT(OUT) :: variable
    INTEGER, INTENT(IN) :: distribution, subarray
    LOGICAL, INTENT(IN), OPTIONAL :: last_in
    INTEGER :: i, j, nm, nsub
    LOGICAL :: last, found
    TYPE(sdf_block_type), POINTER :: cur, b

    IF (PRESENT(last_in)) THEN
      last = last_in
    ELSE
      last = .FALSE.
    END IF

    cur => h%current_block
    nm = cur%ndims
    nsub = 0
    DO i = 1,nm
      IF (cur%variable_ids(i) /= '') nsub = nsub + 1
    END DO

    j = 0
    IF (last) THEN
      DO i = 1,nm
        IF (cur%variable_ids(i) == '') CYCLE
        found = sdf_find_block(h, b, cur%variable_ids(i))
        IF (.NOT. found) RETURN
        h%current_block => b
        j = j + 1
        CALL read_1d_var_last_r4(h, nsub, b%dims, variable, j, distribution, &
            subarray)
      END DO
    ELSE
      DO i = 1,nm
        IF (cur%variable_ids(i) == '') CYCLE
        found = sdf_find_block(h, b, cur%variable_ids(i))
        IF (.NOT. found) RETURN
        h%current_block => b
        j = j + 1
        CALL read_1d_var_first_r4(h, nsub, b%dims, variable, j, distribution, &
            subarray)
      END DO
    END IF

    h%current_block => cur

  END SUBROUTINE read_1d_material_r4



  !----------------------------------------------------------------------------
  ! Code to read a 2D cartesian material mesh in parallel using the
  ! mpitype {distribution} for distribution of data
  !----------------------------------------------------------------------------

  SUBROUTINE read_2d_material_r4(h, variable, distribution, subarray, last_in)

    TYPE(sdf_file_handle) :: h
    REAL(r4), DIMENSION(:,:,:), INTENT(OUT) :: variable
    INTEGER, INTENT(IN) :: distribution, subarray
    LOGICAL, INTENT(IN), OPTIONAL :: last_in
    INTEGER :: i, j, nm, nsub
    LOGICAL :: last, found
    TYPE(sdf_block_type), POINTER :: cur, b

    IF (PRESENT(last_in)) THEN
      last = last_in
    ELSE
      last = .FALSE.
    END IF

    cur => h%current_block
    nm = cur%ndims
    nsub = 0
    DO i = 1,nm
      IF (cur%variable_ids(i) /= '') nsub = nsub + 1
    END DO

    j = 0
    IF (last) THEN
      DO i = 1,nm
        IF (cur%variable_ids(i) == '') CYCLE
        found = sdf_find_block(h, b, cur%variable_ids(i))
        IF (.NOT. found) RETURN
        h%current_block => b
        j = j + 1
        CALL read_2d_var_last_r4(h, nsub, b%dims, variable, j, distribution, &
            subarray)
      END DO
    ELSE
      DO i = 1,nm
        IF (cur%variable_ids(i) == '') CYCLE
        found = sdf_find_block(h, b, cur%variable_ids(i))
        IF (.NOT. found) RETURN
        h%current_block => b
        j = j + 1
        CALL read_2d_var_first_r4(h, nsub, b%dims, variable, j, distribution, &
            subarray)
      END DO
    END IF

    h%current_block => cur

  END SUBROUTINE read_2d_material_r4



  !----------------------------------------------------------------------------
  ! Code to read a 3D cartesian material mesh in parallel using the
  ! mpitype {distribution} for distribution of data
  !----------------------------------------------------------------------------

  SUBROUTINE read_3d_material_r4(h, variable, distribution, subarray, last_in)

    TYPE(sdf_file_handle) :: h
    REAL(r4), DIMENSION(:,:,:,:), INTENT(OUT) :: variable
    INTEGER, INTENT(IN) :: distribution, subarray
    LOGICAL, INTENT(IN), OPTIONAL :: last_in
    INTEGER :: i, j, nm, nsub
    LOGICAL :: last, found
    TYPE(sdf_block_type), POINTER :: cur, b

    IF (PRESENT(last_in)) THEN
      last = last_in
    ELSE
      last = .FALSE.
    END IF

    cur => h%current_block
    nm = cur%ndims
    nsub = 0
    DO i = 1,nm
      IF (cur%variable_ids(i) /= '') nsub = nsub + 1
    END DO

    j = 0
    IF (last) THEN
      DO i = 1,nm
        IF (cur%variable_ids(i) == '') CYCLE
        found = sdf_find_block(h, b, cur%variable_ids(i))
        IF (.NOT. found) RETURN
        h%current_block => b
        j = j + 1
        CALL read_3d_var_last_r4(h, nsub, b%dims, variable, j, distribution, &
            subarray)
      END DO
    ELSE
      DO i = 1,nm
        IF (cur%variable_ids(i) == '') CYCLE
        found = sdf_find_block(h, b, cur%variable_ids(i))
        IF (.NOT. found) RETURN
        h%current_block => b
        j = j + 1
        CALL read_3d_var_first_r4(h, nsub, b%dims, variable, j, distribution, &
            subarray)
      END DO
    END IF

    h%current_block => cur

  END SUBROUTINE read_3d_material_r4



  !----------------------------------------------------------------------------
  ! Code to read either 2D cartesian variable or 1D cartesian multi-material
  ! in parallel
  !----------------------------------------------------------------------------

  SUBROUTINE read_2d_variable_r4(h, variable, distribution, subarray, last_in)

    TYPE(sdf_file_handle) :: h
    REAL(r4), DIMENSION(:,:), INTENT(OUT) :: variable
    INTEGER, INTENT(IN) :: distribution, subarray
    LOGICAL, INTENT(IN), OPTIONAL :: last_in
    TYPE(sdf_block_type), POINTER :: b

    b => h%current_block
    IF (b%blocktype == c_blocktype_plain_variable) THEN
      CALL read_nd_float_r4(h, variable(1,1), distribution, subarray)
    ELSE
      CALL read_1d_material_r4(h, variable, distribution, subarray, last_in)
    END IF

  END SUBROUTINE read_2d_variable_r4



  !----------------------------------------------------------------------------
  ! Code to read either 3D cartesian variable or 2D cartesian multi-material
  ! in parallel
  !----------------------------------------------------------------------------

  SUBROUTINE read_3d_variable_r4(h, variable, distribution, subarray, last_in)

    TYPE(sdf_file_handle) :: h
    REAL(r4), DIMENSION(:,:,:), INTENT(OUT) :: variable
    INTEGER, INTENT(IN) :: distribution, subarray
    LOGICAL, INTENT(IN), OPTIONAL :: last_in
    TYPE(sdf_block_type), POINTER :: b

    b => h%current_block
    IF (b%blocktype == c_blocktype_plain_variable) THEN
      CALL read_nd_float_r4(h, variable(1,1,1), distribution, subarray)
    ELSE
      CALL read_2d_material_r4(h, variable, distribution, subarray, last_in)
    END IF

  END SUBROUTINE read_3d_variable_r4

END MODULE sdf_input_cartesian_r4

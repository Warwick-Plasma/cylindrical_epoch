!
! SDF (Self-Describing Format) Fortran Library
! Copyright (c) 2011-2016, SDF Development Team
!
! Distributed under the terms of the BSD 3-clause License.
! See the LICENSE file for details.
!

MODULE sdf_output_cartesian_r4

  USE sdf_output_cartesian_ru

  IMPLICIT NONE

  INTEGER, PARAMETER, PRIVATE :: sof = 4
  INTEGER, PARAMETER, PRIVATE :: datatype_real = c_datatype_real4
  INTEGER, PARAMETER, PRIVATE :: mpitype_real = MPI_REAL4

CONTAINS

  !----------------------------------------------------------------------------
  ! Code to write a 1D cartesian mesh in serial from the node with
  ! rank {rank_write}
  ! Serial operation, so no need to specify dims
  !----------------------------------------------------------------------------

  SUBROUTINE write_srl_1d_mesh_r4(h, id, name, x, convert_in, &
      dim_labels, dim_units, dim_mults, geometry, rank_write)

    INTEGER, PARAMETER :: ndims = 1
    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name
    REAL(r4), DIMENSION(:), INTENT(IN) :: x
    LOGICAL, INTENT(IN), OPTIONAL :: convert_in
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: dim_labels(:), dim_units(:)
    REAL(r4), DIMENSION(:), INTENT(IN), OPTIONAL :: dim_mults
    INTEGER, INTENT(IN), OPTIONAL :: geometry, rank_write
    REAL(r4), DIMENSION(:), ALLOCATABLE :: r4array
    INTEGER :: errcode, intn
    TYPE(sdf_block_type), POINTER :: b
    LOGICAL :: convert

    CALL sdf_get_next_block(h)
    b => h%current_block

    IF (PRESENT(convert_in)) THEN
      convert = convert_in
    ELSE
      convert = .FALSE.
    END IF

    IF (convert) THEN
      b%type_size = 4
      b%datatype = c_datatype_real4
      b%mpitype = MPI_REAL4
    ELSE
      b%type_size = sof
      b%datatype = datatype_real
      b%mpitype = mpitype_real
    END IF

    IF (PRESENT(geometry)) THEN
      b%geometry = geometry
    ELSE
      b%geometry = c_geometry_cartesian
    END IF
    b%ndims = ndims

    IF (PRESENT(rank_write)) h%rank_master = rank_write

    b%dims(1) = INT(SIZE(x),i4)

    IF (h%rank == h%rank_master) THEN
      b%extents(1) = REAL(MINVAL(x(1:b%dims(1))),r8)
      b%extents(ndims+1) = REAL(MAXVAL(x(1:b%dims(1))),r8)
    END IF

    ! Write header

    CALL write_mesh_meta_r4(h, id, name, dim_labels, dim_units, dim_mults)

    ! Write the actual data

    IF (h%rank == h%rank_master) THEN
      h%current_location = b%data_location

      CALL MPI_FILE_SEEK(h%filehandle, h%current_location, MPI_SEEK_SET, &
          errcode)

      IF (convert) THEN
        intn = b%dims(1)
        ALLOCATE(r4array(intn))

        intn = b%dims(1)
        r4array(1:intn) = REAL(x(1:intn),r4)
        CALL MPI_FILE_WRITE(h%filehandle, r4array, intn, b%mpitype, &
            MPI_STATUS_IGNORE, errcode)

        DEALLOCATE(r4array)
      ELSE
        intn = b%dims(1)
        CALL MPI_FILE_WRITE(h%filehandle, x, intn, b%mpitype, &
            MPI_STATUS_IGNORE, errcode)
      END IF
    END IF

    h%rank_master = h%default_rank
    h%current_location = b%data_location + b%data_length
    b%done_data = .TRUE.

  END SUBROUTINE write_srl_1d_mesh_r4



  !----------------------------------------------------------------------------
  ! Code to write a 2D cartesian mesh in serial from the node with
  ! rank {rank_write}
  ! Serial operation, so no need to specify dims
  !----------------------------------------------------------------------------

  SUBROUTINE write_srl_2d_mesh_r4(h, id, name, x, y, convert_in, &
      dim_labels, dim_units, dim_mults, geometry, rank_write)

    INTEGER, PARAMETER :: ndims = 2
    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name
    REAL(r4), DIMENSION(:), INTENT(IN) :: x, y
    LOGICAL, INTENT(IN), OPTIONAL :: convert_in
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: dim_labels(:), dim_units(:)
    REAL(r4), DIMENSION(:), INTENT(IN), OPTIONAL :: dim_mults
    INTEGER, INTENT(IN), OPTIONAL :: geometry, rank_write
    REAL(r4), DIMENSION(:), ALLOCATABLE :: r4array
    INTEGER :: errcode, intn
    TYPE(sdf_block_type), POINTER :: b
    LOGICAL :: convert

    CALL sdf_get_next_block(h)
    b => h%current_block

    IF (PRESENT(convert_in)) THEN
      convert = convert_in
    ELSE
      convert = .FALSE.
    END IF

    IF (convert) THEN
      b%type_size = 4
      b%datatype = c_datatype_real4
      b%mpitype = MPI_REAL4
    ELSE
      b%type_size = sof
      b%datatype = datatype_real
      b%mpitype = mpitype_real
    END IF

    IF (PRESENT(geometry)) THEN
      b%geometry = geometry
    ELSE
      b%geometry = c_geometry_cartesian
    END IF
    b%ndims = ndims

    IF (PRESENT(rank_write)) h%rank_master = rank_write

    b%dims(1) = INT(SIZE(x),i4)
    b%dims(2) = INT(SIZE(y),i4)

    IF (h%rank == h%rank_master) THEN
      b%extents(1) = REAL(MINVAL(x(1:b%dims(1))),r8)
      b%extents(2) = REAL(MINVAL(y(1:b%dims(2))),r8)
      b%extents(ndims+1) = REAL(MAXVAL(x(1:b%dims(1))),r8)
      b%extents(ndims+2) = REAL(MAXVAL(y(1:b%dims(2))),r8)
    END IF

    ! Write header

    CALL write_mesh_meta_r4(h, id, name, dim_labels, dim_units, dim_mults)

    ! Write the actual data

    IF (h%rank == h%rank_master) THEN
      h%current_location = b%data_location

      CALL MPI_FILE_SEEK(h%filehandle, h%current_location, MPI_SEEK_SET, &
          errcode)

      IF (convert) THEN
        intn = MAX(b%dims(1),b%dims(2))
        ALLOCATE(r4array(intn))

        intn = b%dims(1)
        r4array(1:intn) = REAL(x(1:intn),r4)
        CALL MPI_FILE_WRITE(h%filehandle, r4array, intn, b%mpitype, &
            MPI_STATUS_IGNORE, errcode)

        intn = b%dims(2)
        r4array(1:intn) = REAL(y(1:intn),r4)
        CALL MPI_FILE_WRITE(h%filehandle, r4array, intn, b%mpitype, &
            MPI_STATUS_IGNORE, errcode)

        DEALLOCATE(r4array)
      ELSE
        intn = b%dims(1)
        CALL MPI_FILE_WRITE(h%filehandle, x, intn, b%mpitype, &
            MPI_STATUS_IGNORE, errcode)

        intn = b%dims(2)
        CALL MPI_FILE_WRITE(h%filehandle, y, intn, b%mpitype, &
            MPI_STATUS_IGNORE, errcode)
      END IF
    END IF

    h%rank_master = h%default_rank
    h%current_location = b%data_location + b%data_length
    b%done_data = .TRUE.

  END SUBROUTINE write_srl_2d_mesh_r4



  !----------------------------------------------------------------------------
  ! Code to write a 3D cartesian mesh in serial from the node with
  ! rank {rank_write}
  ! Serial operation, so no need to specify dims
  !----------------------------------------------------------------------------

  SUBROUTINE write_srl_3d_mesh_r4(h, id, name, x, y, z, convert_in, &
      dim_labels, dim_units, dim_mults, geometry, rank_write)

    INTEGER, PARAMETER :: ndims = 3
    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name
    REAL(r4), DIMENSION(:), INTENT(IN) :: x, y, z
    LOGICAL, INTENT(IN), OPTIONAL :: convert_in
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: dim_labels(:), dim_units(:)
    REAL(r4), DIMENSION(:), INTENT(IN), OPTIONAL :: dim_mults
    INTEGER, INTENT(IN), OPTIONAL :: geometry, rank_write
    REAL(r4), DIMENSION(:), ALLOCATABLE :: r4array
    INTEGER :: errcode, intn
    TYPE(sdf_block_type), POINTER :: b
    LOGICAL :: convert

    CALL sdf_get_next_block(h)
    b => h%current_block

    IF (PRESENT(convert_in)) THEN
      convert = convert_in
    ELSE
      convert = .FALSE.
    END IF

    IF (convert) THEN
      b%type_size = 4
      b%datatype = c_datatype_real4
      b%mpitype = MPI_REAL4
    ELSE
      b%type_size = sof
      b%datatype = datatype_real
      b%mpitype = mpitype_real
    END IF

    IF (PRESENT(geometry)) THEN
      b%geometry = geometry
    ELSE
      b%geometry = c_geometry_cartesian
    END IF
    b%ndims = ndims

    IF (PRESENT(rank_write)) h%rank_master = rank_write

    b%dims(1) = INT(SIZE(x),i4)
    b%dims(2) = INT(SIZE(y),i4)
    b%dims(3) = INT(SIZE(z),i4)

    IF (h%rank == h%rank_master) THEN
      b%extents(1) = REAL(MINVAL(x(1:b%dims(1))),r8)
      b%extents(2) = REAL(MINVAL(y(1:b%dims(2))),r8)
      b%extents(3) = REAL(MINVAL(z(1:b%dims(3))),r8)
      b%extents(ndims+1) = REAL(MAXVAL(x(1:b%dims(1))),r8)
      b%extents(ndims+2) = REAL(MAXVAL(y(1:b%dims(2))),r8)
      b%extents(ndims+3) = REAL(MAXVAL(z(1:b%dims(3))),r8)
    END IF

    ! Write header

    CALL write_mesh_meta_r4(h, id, name, dim_labels, dim_units, dim_mults)

    ! Write the actual data

    IF (h%rank == h%rank_master) THEN
      h%current_location = b%data_location

      CALL MPI_FILE_SEEK(h%filehandle, h%current_location, MPI_SEEK_SET, &
          errcode)

      IF (convert) THEN
        intn = MAX(MAX(b%dims(1),b%dims(2)),b%dims(3))
        ALLOCATE(r4array(intn))

        intn = b%dims(1)
        r4array(1:intn) = REAL(x(1:intn),r4)
        CALL MPI_FILE_WRITE(h%filehandle, r4array, intn, b%mpitype, &
            MPI_STATUS_IGNORE, errcode)

        intn = b%dims(2)
        r4array(1:intn) = REAL(y(1:intn),r4)
        CALL MPI_FILE_WRITE(h%filehandle, r4array, intn, b%mpitype, &
            MPI_STATUS_IGNORE, errcode)

        intn = b%dims(3)
        r4array(1:intn) = REAL(z(1:intn),r4)
        CALL MPI_FILE_WRITE(h%filehandle, r4array, intn, b%mpitype, &
            MPI_STATUS_IGNORE, errcode)

        DEALLOCATE(r4array)
      ELSE
        intn = b%dims(1)
        CALL MPI_FILE_WRITE(h%filehandle, x, intn, b%mpitype, &
            MPI_STATUS_IGNORE, errcode)

        intn = b%dims(2)
        CALL MPI_FILE_WRITE(h%filehandle, y, intn, b%mpitype, &
            MPI_STATUS_IGNORE, errcode)

        intn = b%dims(3)
        CALL MPI_FILE_WRITE(h%filehandle, z, intn, b%mpitype, &
            MPI_STATUS_IGNORE, errcode)
      END IF
    END IF

    h%rank_master = h%default_rank
    h%current_location = b%data_location + b%data_length
    b%done_data = .TRUE.

  END SUBROUTINE write_srl_3d_mesh_r4



  !----------------------------------------------------------------------------
  ! Code to write a 1D cartesian mesh in parallel using the
  ! mpitype {distribution} for distribution of data
  ! It's up to the coder to design the distribution parallel operation, so
  ! need global dims
  !----------------------------------------------------------------------------

  SUBROUTINE write_1d_mesh_r4(h, id, name, x, dims, xmin, xmax, &
      distribution, subarray, convert_in, dim_labels, &
      dim_units, dim_mults, geometry)

    INTEGER, PARAMETER :: ndims = 1
    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name
    REAL(r4), DIMENSION(:), INTENT(IN) :: x
    INTEGER, DIMENSION(:), INTENT(IN) :: dims
    REAL(r4), INTENT(IN) :: xmin, xmax
    INTEGER, DIMENSION(:), INTENT(IN) :: distribution, subarray
    LOGICAL, INTENT(IN), OPTIONAL :: convert_in
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: dim_labels(:), dim_units(:)
    REAL(r4), DIMENSION(:), INTENT(IN), OPTIONAL :: dim_mults
    INTEGER, INTENT(IN), OPTIONAL :: geometry
    REAL(r8), DIMENSION(ndims) :: gmn, gmx
    REAL(r4), DIMENSION(:), ALLOCATABLE :: r4array
    INTEGER :: i, errcode, intn = 0
    TYPE(sdf_block_type), POINTER :: b
    LOGICAL :: convert

    CALL sdf_get_next_block(h)
    b => h%current_block

    IF (PRESENT(convert_in)) THEN
      convert = convert_in
    ELSE
      convert = .FALSE.
    END IF

    IF (convert) THEN
      b%type_size = 4
      b%datatype = c_datatype_real4
      b%mpitype = MPI_REAL4

      intn = INT(SIZE(x),i4)
      ALLOCATE(r4array(intn))
    ELSE
      b%type_size = sof
      b%datatype = datatype_real
      b%mpitype = mpitype_real
    END IF

    IF (PRESENT(geometry)) THEN
      b%geometry = geometry
    ELSE
      b%geometry = c_geometry_cartesian
    END IF
    b%ndims = ndims

    DO i = 1,ndims
      b%dims(i) = INT(dims(i),i4)
    END DO

    gmn(1) = xmin
    gmx(1) = xmax

    CALL MPI_ALLREDUCE(gmn, b%extents, ndims, MPI_REAL8, MPI_MIN, &
        h%comm, errcode)
    CALL MPI_ALLREDUCE(gmx, b%extents(ndims+1), ndims, MPI_REAL8, MPI_MAX, &
        h%comm, errcode)

    ! Write header

    CALL write_mesh_meta_r4(h, id, name, dim_labels, dim_units, dim_mults)

    ! Write the actual data

    CALL MPI_FILE_SET_VIEW(h%filehandle, h%current_location, MPI_BYTE, &
        distribution(1), 'native', MPI_INFO_NULL, errcode)
    IF (convert) THEN
      r4array(1:intn) = REAL(x(1:intn),r4)
      CALL MPI_FILE_WRITE_ALL(h%filehandle, r4array, 1, subarray(1), &
          MPI_STATUS_IGNORE, errcode)
      DEALLOCATE(r4array)
    ELSE
      CALL MPI_FILE_WRITE_ALL(h%filehandle, x, 1, subarray(1), &
          MPI_STATUS_IGNORE, errcode)
    END IF

    CALL MPI_FILE_SET_VIEW(h%filehandle, c_off0, MPI_BYTE, MPI_BYTE, 'native', &
        MPI_INFO_NULL, errcode)

    h%current_location = b%data_location + b%data_length
    b%done_data = .TRUE.

  END SUBROUTINE write_1d_mesh_r4



  !----------------------------------------------------------------------------
  ! Code to write a 2D cartesian mesh in parallel using the
  ! mpitype {distribution} for distribution of data
  ! It's up to the coder to design the distribution parallel operation, so
  ! need global dims
  !----------------------------------------------------------------------------

  SUBROUTINE write_2d_mesh_r4(h, id, name, x, y, dims, xmin, xmax, &
      ymin, ymax, distribution, subarray, convert_in, dim_labels, &
      dim_units, dim_mults, geometry)

    INTEGER, PARAMETER :: ndims = 2
    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name
    REAL(r4), DIMENSION(:), INTENT(IN) :: x, y
    INTEGER, DIMENSION(:), INTENT(IN) :: dims
    REAL(r4), INTENT(IN) :: xmin, xmax, ymin, ymax
    INTEGER, DIMENSION(:), INTENT(IN) :: distribution, subarray
    LOGICAL, INTENT(IN), OPTIONAL :: convert_in
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: dim_labels(:), dim_units(:)
    REAL(r4), DIMENSION(:), INTENT(IN), OPTIONAL :: dim_mults
    INTEGER, INTENT(IN), OPTIONAL :: geometry
    REAL(r8), DIMENSION(ndims) :: gmn, gmx
    REAL(r4), DIMENSION(:), ALLOCATABLE :: r4array
    INTEGER :: i, errcode, sz(ndims), intn
    TYPE(sdf_block_type), POINTER :: b
    LOGICAL :: convert

    CALL sdf_get_next_block(h)
    b => h%current_block

    IF (PRESENT(convert_in)) THEN
      convert = convert_in
    ELSE
      convert = .FALSE.
    END IF

    IF (convert) THEN
      b%type_size = 4
      b%datatype = c_datatype_real4
      b%mpitype = MPI_REAL4

      sz(1) = INT(SIZE(x),i4)
      sz(2) = INT(SIZE(y),i4)

      intn = MAX(sz(1),sz(2))
      ALLOCATE(r4array(intn))
    ELSE
      b%type_size = sof
      b%datatype = datatype_real
      b%mpitype = mpitype_real
    END IF

    IF (PRESENT(geometry)) THEN
      b%geometry = geometry
    ELSE
      b%geometry = c_geometry_cartesian
    END IF
    b%ndims = ndims

    DO i = 1,ndims
      b%dims(i) = INT(dims(i),i4)
    END DO

    gmn(1) = xmin
    gmx(1) = xmax
    gmn(2) = ymin
    gmx(2) = ymax

    CALL MPI_ALLREDUCE(gmn, b%extents, ndims, MPI_REAL8, MPI_MIN, &
        h%comm, errcode)
    CALL MPI_ALLREDUCE(gmx, b%extents(ndims+1), ndims, MPI_REAL8, MPI_MAX, &
        h%comm, errcode)

    ! Write header

    CALL write_mesh_meta_r4(h, id, name, dim_labels, dim_units, dim_mults)

    ! Write the actual data

    CALL MPI_FILE_SET_VIEW(h%filehandle, h%current_location, MPI_BYTE, &
        distribution(1), 'native', MPI_INFO_NULL, errcode)
    IF (convert) THEN
      intn = sz(1)
      r4array(1:intn) = REAL(x(1:intn),r4)
      CALL MPI_FILE_WRITE_ALL(h%filehandle, r4array, 1, subarray(1), &
          MPI_STATUS_IGNORE, errcode)
    ELSE
      CALL MPI_FILE_WRITE_ALL(h%filehandle, x, 1, subarray(1), &
          MPI_STATUS_IGNORE, errcode)
    END IF

    h%current_location = h%current_location + b%dims(1) * b%type_size

    CALL MPI_FILE_SET_VIEW(h%filehandle, h%current_location, MPI_BYTE, &
        distribution(2), 'native', MPI_INFO_NULL, errcode)
    IF (convert) THEN
      intn = sz(2)
      r4array(1:intn) = REAL(y(1:intn),r4)
      CALL MPI_FILE_WRITE_ALL(h%filehandle, r4array, 1, subarray(2), &
          MPI_STATUS_IGNORE, errcode)
      DEALLOCATE(r4array)
    ELSE
      CALL MPI_FILE_WRITE_ALL(h%filehandle, y, 1, subarray(2), &
          MPI_STATUS_IGNORE, errcode)
    END IF

    CALL MPI_FILE_SET_VIEW(h%filehandle, c_off0, MPI_BYTE, MPI_BYTE, 'native', &
        MPI_INFO_NULL, errcode)

    h%current_location = b%data_location + b%data_length
    b%done_data = .TRUE.

  END SUBROUTINE write_2d_mesh_r4



  !----------------------------------------------------------------------------
  ! Code to write a 3D cartesian mesh in parallel using the
  ! mpitype {distribution} for distribution of data
  ! It's up to the coder to design the distribution parallel operation, so
  ! need global dims
  !----------------------------------------------------------------------------

  SUBROUTINE write_3d_mesh_r4(h, id, name, x, y, z, dims, xmin, xmax, &
      ymin, ymax, zmin, zmax, distribution, subarray, convert_in, dim_labels, &
      dim_units, dim_mults, geometry)

    INTEGER, PARAMETER :: ndims = 3
    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name
    REAL(r4), DIMENSION(:), INTENT(IN) :: x, y, z
    INTEGER, DIMENSION(:), INTENT(IN) :: dims
    REAL(r4), INTENT(IN) :: xmin, xmax, ymin, ymax, zmin, zmax
    INTEGER, DIMENSION(:), INTENT(IN) :: distribution, subarray
    LOGICAL, INTENT(IN), OPTIONAL :: convert_in
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: dim_labels(:), dim_units(:)
    REAL(r4), DIMENSION(:), INTENT(IN), OPTIONAL :: dim_mults
    INTEGER, INTENT(IN), OPTIONAL :: geometry
    REAL(r8), DIMENSION(ndims) :: gmn, gmx
    REAL(r4), DIMENSION(:), ALLOCATABLE :: r4array
    INTEGER :: i, errcode, sz(ndims), intn
    TYPE(sdf_block_type), POINTER :: b
    LOGICAL :: convert

    CALL sdf_get_next_block(h)
    b => h%current_block

    IF (PRESENT(convert_in)) THEN
      convert = convert_in
    ELSE
      convert = .FALSE.
    END IF

    IF (convert) THEN
      b%type_size = 4
      b%datatype = c_datatype_real4
      b%mpitype = MPI_REAL4

      sz(1) = INT(SIZE(x),i4)
      sz(2) = INT(SIZE(y),i4)
      sz(3) = INT(SIZE(z),i4)

      intn = MAX(MAX(sz(1),sz(2)),sz(3))
      ALLOCATE(r4array(intn))
    ELSE
      b%type_size = sof
      b%datatype = datatype_real
      b%mpitype = mpitype_real
    END IF

    IF (PRESENT(geometry)) THEN
      b%geometry = geometry
    ELSE
      b%geometry = c_geometry_cartesian
    END IF
    b%ndims = ndims

    DO i = 1,ndims
      b%dims(i) = INT(dims(i),i4)
    END DO

    gmn(1) = xmin
    gmx(1) = xmax
    gmn(2) = ymin
    gmx(2) = ymax
    gmn(3) = zmin
    gmx(3) = zmax

    CALL MPI_ALLREDUCE(gmn, b%extents, ndims, MPI_REAL8, MPI_MIN, &
        h%comm, errcode)
    CALL MPI_ALLREDUCE(gmx, b%extents(ndims+1), ndims, MPI_REAL8, MPI_MAX, &
        h%comm, errcode)

    ! Write header

    CALL write_mesh_meta_r4(h, id, name, dim_labels, dim_units, dim_mults)

    ! Write the actual data

    CALL MPI_FILE_SET_VIEW(h%filehandle, h%current_location, MPI_BYTE, &
        distribution(1), 'native', MPI_INFO_NULL, errcode)
    IF (convert) THEN
      intn = sz(1)
      r4array(1:intn) = REAL(x(1:intn),r4)
      CALL MPI_FILE_WRITE_ALL(h%filehandle, r4array, 1, subarray(1), &
          MPI_STATUS_IGNORE, errcode)
    ELSE
      CALL MPI_FILE_WRITE_ALL(h%filehandle, x, 1, subarray(1), &
          MPI_STATUS_IGNORE, errcode)
    END IF

    h%current_location = h%current_location + b%dims(1) * b%type_size

    CALL MPI_FILE_SET_VIEW(h%filehandle, h%current_location, MPI_BYTE, &
        distribution(2), 'native', MPI_INFO_NULL, errcode)
    IF (convert) THEN
      intn = sz(2)
      r4array(1:intn) = REAL(y(1:intn),r4)
      CALL MPI_FILE_WRITE_ALL(h%filehandle, r4array, 1, subarray(2), &
          MPI_STATUS_IGNORE, errcode)
    ELSE
      CALL MPI_FILE_WRITE_ALL(h%filehandle, y, 1, subarray(2), &
          MPI_STATUS_IGNORE, errcode)
    END IF

    h%current_location = h%current_location + b%dims(2) * b%type_size

    CALL MPI_FILE_SET_VIEW(h%filehandle, h%current_location, MPI_BYTE, &
        distribution(3), 'native', MPI_INFO_NULL, errcode)
    IF (convert) THEN
      intn = sz(3)
      r4array(1:intn) = REAL(z(1:intn),r4)
      CALL MPI_FILE_WRITE_ALL(h%filehandle, r4array, 1, subarray(3), &
          MPI_STATUS_IGNORE, errcode)
      DEALLOCATE(r4array)
    ELSE
      CALL MPI_FILE_WRITE_ALL(h%filehandle, z, 1, subarray(3), &
          MPI_STATUS_IGNORE, errcode)
    END IF

    CALL MPI_FILE_SET_VIEW(h%filehandle, c_off0, MPI_BYTE, MPI_BYTE, 'native', &
        MPI_INFO_NULL, errcode)

    h%current_location = b%data_location + b%data_length
    b%done_data = .TRUE.

  END SUBROUTINE write_3d_mesh_r4



  !----------------------------------------------------------------------------
  ! Code to write a 2D lagrangian mesh in serial from the node with
  ! rank {rank_write}
  ! Serial operation, so no need to specify dims
  !----------------------------------------------------------------------------

  SUBROUTINE write_srl_2d_lag_mesh_r4(h, id, name, x, y, convert_in, &
      dim_labels, dim_units, dim_mults, geometry, rank_write)

    INTEGER, PARAMETER :: ndims = 2
    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name
    REAL(r4), DIMENSION(:,:), INTENT(IN) :: x, y
    LOGICAL, INTENT(IN), OPTIONAL :: convert_in
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: dim_labels(:), dim_units(:)
    REAL(r4), DIMENSION(:), INTENT(IN), OPTIONAL :: dim_mults
    INTEGER, INTENT(IN), OPTIONAL :: geometry, rank_write
    REAL(r4), DIMENSION(:,:), ALLOCATABLE :: r4array
    INTEGER :: i, errcode, intn
    TYPE(sdf_block_type), POINTER :: b
    LOGICAL :: convert

    CALL sdf_get_next_block(h)
    b => h%current_block

    intn = 1
    DO i = 1,ndims
      b%dims(i) = INT(SIZE(x,i),i4)
      intn = intn * b%dims(i)
    END DO

    IF (PRESENT(convert_in)) THEN
      convert = convert_in
    ELSE
      convert = .FALSE.
    END IF

    IF (convert) THEN
      b%type_size = 4
      b%datatype = c_datatype_real4
      b%mpitype = MPI_REAL4

      ALLOCATE(r4array(b%dims(1),b%dims(2)))
    ELSE
      b%type_size = sof
      b%datatype = datatype_real
      b%mpitype = mpitype_real
    END IF

    IF (PRESENT(geometry)) THEN
      b%geometry = geometry
    ELSE
      b%geometry = c_geometry_cartesian
    END IF
    b%ndims = ndims

    IF (h%rank == h%rank_master) THEN
      b%extents(1) = REAL(MINVAL(x),r8)
      b%extents(2) = REAL(MINVAL(y),r8)
      b%extents(ndims+1) = REAL(MAXVAL(x),r8)
      b%extents(ndims+2) = REAL(MAXVAL(y),r8)
    END IF

    ! Write header

    b%blocktype = c_blocktype_lagrangian_mesh
    CALL write_mesh_meta_r4(h, id, name, dim_labels, dim_units, dim_mults)

    ! Write the actual data

    IF (h%rank == h%rank_master) THEN
      h%current_location = b%data_location

      CALL MPI_FILE_SEEK(h%filehandle, h%current_location, MPI_SEEK_SET, &
          errcode)

      IF (convert) THEN
        r4array = REAL(x,r4)
        CALL MPI_FILE_WRITE(h%filehandle, r4array, intn, b%mpitype, &
            MPI_STATUS_IGNORE, errcode)

        r4array = REAL(y,r4)
        CALL MPI_FILE_WRITE(h%filehandle, r4array, intn, b%mpitype, &
            MPI_STATUS_IGNORE, errcode)

        DEALLOCATE(r4array)
      ELSE
        CALL MPI_FILE_WRITE(h%filehandle, x, intn, b%mpitype, &
            MPI_STATUS_IGNORE, errcode)

        CALL MPI_FILE_WRITE(h%filehandle, y, intn, b%mpitype, &
            MPI_STATUS_IGNORE, errcode)
      END IF
    END IF

    h%rank_master = h%default_rank
    h%current_location = b%data_location + b%data_length
    b%done_data = .TRUE.

  END SUBROUTINE write_srl_2d_lag_mesh_r4



  !----------------------------------------------------------------------------
  ! Code to write a 3D lagrangian mesh in serial from the node with
  ! rank {rank_write}
  ! Serial operation, so no need to specify dims
  !----------------------------------------------------------------------------

  SUBROUTINE write_srl_3d_lag_mesh_r4(h, id, name, x, y, z, convert_in, &
      dim_labels, dim_units, dim_mults, geometry, rank_write)

    INTEGER, PARAMETER :: ndims = 3
    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name
    REAL(r4), DIMENSION(:,:,:), INTENT(IN) :: x, y, z
    LOGICAL, INTENT(IN), OPTIONAL :: convert_in
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: dim_labels(:), dim_units(:)
    REAL(r4), DIMENSION(:), INTENT(IN), OPTIONAL :: dim_mults
    INTEGER, INTENT(IN), OPTIONAL :: geometry, rank_write
    REAL(r4), DIMENSION(:,:,:), ALLOCATABLE :: r4array
    INTEGER :: i, errcode, intn
    TYPE(sdf_block_type), POINTER :: b
    LOGICAL :: convert

    CALL sdf_get_next_block(h)
    b => h%current_block

    intn = 1
    DO i = 1,ndims
      b%dims(i) = INT(SIZE(x,i),i4)
      intn = intn * b%dims(i)
    END DO

    IF (PRESENT(convert_in)) THEN
      convert = convert_in
    ELSE
      convert = .FALSE.
    END IF

    IF (convert) THEN
      b%type_size = 4
      b%datatype = c_datatype_real4
      b%mpitype = MPI_REAL4

      ALLOCATE(r4array(b%dims(1),b%dims(2),b%dims(3)))
    ELSE
      b%type_size = sof
      b%datatype = datatype_real
      b%mpitype = mpitype_real
    END IF

    IF (PRESENT(geometry)) THEN
      b%geometry = geometry
    ELSE
      b%geometry = c_geometry_cartesian
    END IF
    b%ndims = ndims

    IF (h%rank == h%rank_master) THEN
      b%extents(1) = REAL(MINVAL(x),r8)
      b%extents(2) = REAL(MINVAL(y),r8)
      b%extents(3) = REAL(MINVAL(z),r8)
      b%extents(ndims+1) = REAL(MAXVAL(x),r8)
      b%extents(ndims+2) = REAL(MAXVAL(y),r8)
      b%extents(ndims+3) = REAL(MAXVAL(z),r8)
    END IF

    ! Write header

    b%blocktype = c_blocktype_lagrangian_mesh
    CALL write_mesh_meta_r4(h, id, name, dim_labels, dim_units, dim_mults)

    ! Write the actual data

    IF (h%rank == h%rank_master) THEN
      h%current_location = b%data_location

      CALL MPI_FILE_SEEK(h%filehandle, h%current_location, MPI_SEEK_SET, &
          errcode)

      IF (convert) THEN
        r4array = REAL(x,r4)
        CALL MPI_FILE_WRITE(h%filehandle, r4array, intn, b%mpitype, &
            MPI_STATUS_IGNORE, errcode)

        r4array = REAL(y,r4)
        CALL MPI_FILE_WRITE(h%filehandle, r4array, intn, b%mpitype, &
            MPI_STATUS_IGNORE, errcode)

        r4array = REAL(z,r4)
        CALL MPI_FILE_WRITE(h%filehandle, r4array, intn, b%mpitype, &
            MPI_STATUS_IGNORE, errcode)

        DEALLOCATE(r4array)
      ELSE
        CALL MPI_FILE_WRITE(h%filehandle, x, intn, b%mpitype, &
            MPI_STATUS_IGNORE, errcode)

        CALL MPI_FILE_WRITE(h%filehandle, y, intn, b%mpitype, &
            MPI_STATUS_IGNORE, errcode)

        CALL MPI_FILE_WRITE(h%filehandle, z, intn, b%mpitype, &
            MPI_STATUS_IGNORE, errcode)
      END IF
    END IF

    h%rank_master = h%default_rank
    h%current_location = b%data_location + b%data_length
    b%done_data = .TRUE.

  END SUBROUTINE write_srl_3d_lag_mesh_r4



  !----------------------------------------------------------------------------
  ! Code to write a 1D path mesh for a 2D domain in serial from the node with
  ! rank {rank_write}
  ! Serial operation, so no need to specify dims
  !----------------------------------------------------------------------------

  SUBROUTINE write_srl_2d_path_mesh_r4(h, id, name, x, y, convert_in, &
      dim_labels, dim_units, dim_mults, geometry, rank_write)

    INTEGER, PARAMETER :: ndims = 2
    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name
    REAL(r4), DIMENSION(:), INTENT(IN) :: x, y
    LOGICAL, INTENT(IN), OPTIONAL :: convert_in
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: dim_labels(:), dim_units(:)
    REAL(r4), DIMENSION(:), INTENT(IN), OPTIONAL :: dim_mults
    INTEGER, INTENT(IN), OPTIONAL :: geometry, rank_write
    REAL(r4), DIMENSION(:), ALLOCATABLE :: r4array
    INTEGER :: errcode, intn
    TYPE(sdf_block_type), POINTER :: b
    LOGICAL :: convert

    CALL sdf_get_next_block(h)
    b => h%current_block

    intn = INT(SIZE(x),i4)
    b%dims(:) = 1
    b%dims(1) = intn

    IF (PRESENT(convert_in)) THEN
      convert = convert_in
    ELSE
      convert = .FALSE.
    ENDIF

    IF (convert) THEN
      b%type_size = 4
      b%datatype = c_datatype_real4
      b%mpitype = MPI_REAL4

      ALLOCATE(r4array(b%dims(1)))
    ELSE
      b%type_size = sof
      b%datatype = datatype_real
      b%mpitype = mpitype_real
    ENDIF

    IF (PRESENT(geometry)) THEN
      b%geometry = geometry
    ELSE
      b%geometry = c_geometry_cartesian
    ENDIF
    b%ndims = ndims

    IF (h%rank == h%rank_master) THEN
      b%extents(1) = REAL(MINVAL(x(1:intn)),r8)
      b%extents(2) = REAL(MINVAL(y(1:intn)),r8)
      b%extents(ndims+1) = REAL(MAXVAL(x(1:intn)),r8)
      b%extents(ndims+2) = REAL(MAXVAL(y(1:intn)),r8)
    ENDIF

    ! Write header

    b%blocktype = c_blocktype_lagrangian_mesh
    CALL write_mesh_meta_r4(h, id, name, dim_labels, dim_units, dim_mults)

    ! Write the actual data

    IF (h%rank == h%rank_master) THEN
      h%current_location = b%data_location

      CALL MPI_FILE_SEEK(h%filehandle, h%current_location, MPI_SEEK_SET, &
          errcode)

      IF (convert) THEN
        r4array = REAL(x,r4)
        CALL MPI_FILE_WRITE(h%filehandle, r4array, intn, b%mpitype, &
            MPI_STATUS_IGNORE, errcode)

        r4array = REAL(y,r4)
        CALL MPI_FILE_WRITE(h%filehandle, r4array, intn, b%mpitype, &
            MPI_STATUS_IGNORE, errcode)

        DEALLOCATE(r4array)
      ELSE
        CALL MPI_FILE_WRITE(h%filehandle, x, intn, b%mpitype, &
            MPI_STATUS_IGNORE, errcode)

        CALL MPI_FILE_WRITE(h%filehandle, y, intn, b%mpitype, &
            MPI_STATUS_IGNORE, errcode)
      ENDIF
    ENDIF

    h%rank_master = h%default_rank
    h%current_location = b%data_location + b%data_length
    b%done_data = .TRUE.

  END SUBROUTINE write_srl_2d_path_mesh_r4



  !----------------------------------------------------------------------------
  ! Code to write a 1D path mesh for a 3D domain in serial from the node with
  ! rank {rank_write}
  ! Serial operation, so no need to specify dims
  !----------------------------------------------------------------------------

  SUBROUTINE write_srl_3d_path_mesh_r4(h, id, name, x, y, z, convert_in, &
      dim_labels, dim_units, dim_mults, geometry, rank_write)

    INTEGER, PARAMETER :: ndims = 3
    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name
    REAL(r4), DIMENSION(:), INTENT(IN) :: x, y, z
    LOGICAL, INTENT(IN), OPTIONAL :: convert_in
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: dim_labels(:), dim_units(:)
    REAL(r4), DIMENSION(:), INTENT(IN), OPTIONAL :: dim_mults
    INTEGER, INTENT(IN), OPTIONAL :: geometry, rank_write
    REAL(r4), DIMENSION(:), ALLOCATABLE :: r4array
    INTEGER :: errcode, intn
    TYPE(sdf_block_type), POINTER :: b
    LOGICAL :: convert

    CALL sdf_get_next_block(h)
    b => h%current_block

    intn = INT(SIZE(x),i4)
    b%dims(:) = 1
    b%dims(1) = intn

    IF (PRESENT(convert_in)) THEN
      convert = convert_in
    ELSE
      convert = .FALSE.
    ENDIF

    IF (convert) THEN
      b%type_size = 4
      b%datatype = c_datatype_real4
      b%mpitype = MPI_REAL4

      ALLOCATE(r4array(b%dims(1)))
    ELSE
      b%type_size = sof
      b%datatype = datatype_real
      b%mpitype = mpitype_real
    ENDIF

    IF (PRESENT(geometry)) THEN
      b%geometry = geometry
    ELSE
      b%geometry = c_geometry_cartesian
    ENDIF
    b%ndims = ndims

    IF (h%rank == h%rank_master) THEN
      b%extents(1) = REAL(MINVAL(x(1:intn)),r8)
      b%extents(2) = REAL(MINVAL(y(1:intn)),r8)
      b%extents(3) = REAL(MINVAL(z(1:intn)),r8)
      b%extents(ndims+1) = REAL(MAXVAL(x(1:intn)),r8)
      b%extents(ndims+2) = REAL(MAXVAL(y(1:intn)),r8)
      b%extents(ndims+3) = REAL(MAXVAL(z(1:intn)),r8)
    ENDIF

    ! Write header

    b%blocktype = c_blocktype_lagrangian_mesh
    CALL write_mesh_meta_r4(h, id, name, dim_labels, dim_units, dim_mults)

    ! Write the actual data

    IF (h%rank == h%rank_master) THEN
      h%current_location = b%data_location

      CALL MPI_FILE_SEEK(h%filehandle, h%current_location, MPI_SEEK_SET, &
          errcode)

      IF (convert) THEN
        r4array = REAL(x,r4)
        CALL MPI_FILE_WRITE(h%filehandle, r4array, intn, b%mpitype, &
            MPI_STATUS_IGNORE, errcode)

        r4array = REAL(y,r4)
        CALL MPI_FILE_WRITE(h%filehandle, r4array, intn, b%mpitype, &
            MPI_STATUS_IGNORE, errcode)

        r4array = REAL(z,r4)
        CALL MPI_FILE_WRITE(h%filehandle, r4array, intn, b%mpitype, &
            MPI_STATUS_IGNORE, errcode)

        DEALLOCATE(r4array)
      ELSE
        CALL MPI_FILE_WRITE(h%filehandle, x, intn, b%mpitype, &
            MPI_STATUS_IGNORE, errcode)

        CALL MPI_FILE_WRITE(h%filehandle, y, intn, b%mpitype, &
            MPI_STATUS_IGNORE, errcode)

        CALL MPI_FILE_WRITE(h%filehandle, z, intn, b%mpitype, &
            MPI_STATUS_IGNORE, errcode)
      ENDIF
    ENDIF

    h%rank_master = h%default_rank
    h%current_location = b%data_location + b%data_length
    b%done_data = .TRUE.

  END SUBROUTINE write_srl_3d_path_mesh_r4



  !----------------------------------------------------------------------------
  ! Code to write a 1D lagrangian mesh in parallel using the
  ! mpitype {distribution} for distribution of data
  ! It's up to the coder to design the distribution parallel operation, so
  ! need global dims
  !----------------------------------------------------------------------------

  SUBROUTINE write_1d_lag_mesh_r4(h, id, name, x, dims, xmin, xmax, &
      distribution, subarray, convert_in, dim_labels, &
      dim_units, dim_mults, geometry)

    INTEGER, PARAMETER :: ndims = 1
    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name
    REAL(r4), DIMENSION(:), INTENT(IN) :: x
    INTEGER, DIMENSION(:), INTENT(IN) :: dims
    REAL(r4), INTENT(IN) :: xmin, xmax
    INTEGER, INTENT(IN) :: distribution, subarray
    LOGICAL, INTENT(IN), OPTIONAL :: convert_in
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: dim_labels(:), dim_units(:)
    REAL(r4), DIMENSION(:), INTENT(IN), OPTIONAL :: dim_mults
    INTEGER, INTENT(IN), OPTIONAL :: geometry
    REAL(r8), DIMENSION(ndims) :: gmn, gmx
    REAL(r4), DIMENSION(:), ALLOCATABLE :: r4array
    INTEGER :: i, errcode, sz(ndims)
    TYPE(sdf_block_type), POINTER :: b
    LOGICAL :: convert

    CALL sdf_get_next_block(h)
    b => h%current_block

    IF (PRESENT(convert_in)) THEN
      convert = convert_in
    ELSE
      convert = .FALSE.
    END IF

    IF (convert) THEN
      b%type_size = 4
      b%datatype = c_datatype_real4
      b%mpitype = MPI_REAL4

      sz(1) = INT(SIZE(x,1),i4)

      ALLOCATE(r4array(sz(1)))
    ELSE
      b%type_size = sof
      b%datatype = datatype_real
      b%mpitype = mpitype_real
    END IF

    IF (PRESENT(geometry)) THEN
      b%geometry = geometry
    ELSE
      b%geometry = c_geometry_cartesian
    END IF
    b%ndims = ndims

    DO i = 1,ndims
      b%dims(i) = INT(dims(i),i4)
    END DO

    gmn(1) = xmin
    gmx(1) = xmax

    CALL MPI_ALLREDUCE(gmn, b%extents, ndims, MPI_REAL8, MPI_MIN, &
        h%comm, errcode)
    CALL MPI_ALLREDUCE(gmx, b%extents(ndims+1), ndims, MPI_REAL8, MPI_MAX, &
        h%comm, errcode)

    ! Write header

    b%blocktype = c_blocktype_lagrangian_mesh
    CALL write_mesh_meta_r4(h, id, name, dim_labels, dim_units, dim_mults)

    ! Write the actual data

    CALL MPI_FILE_SET_VIEW(h%filehandle, h%current_location, MPI_BYTE, &
        distribution, 'native', MPI_INFO_NULL, errcode)
    IF (convert) THEN
      r4array = REAL(x,r4)
      CALL MPI_FILE_WRITE_ALL(h%filehandle, r4array, 1, subarray, &
          MPI_STATUS_IGNORE, errcode)

      DEALLOCATE(r4array)
    ELSE
      CALL MPI_FILE_WRITE_ALL(h%filehandle, x, 1, subarray, &
          MPI_STATUS_IGNORE, errcode)
    END IF

    CALL MPI_FILE_SET_VIEW(h%filehandle, c_off0, MPI_BYTE, MPI_BYTE, 'native', &
        MPI_INFO_NULL, errcode)

    h%current_location = b%data_location + b%data_length
    b%done_data = .TRUE.

  END SUBROUTINE write_1d_lag_mesh_r4



  !----------------------------------------------------------------------------
  ! Code to write a 2D lagrangian mesh in parallel using the
  ! mpitype {distribution} for distribution of data
  ! It's up to the coder to design the distribution parallel operation, so
  ! need global dims
  !----------------------------------------------------------------------------

  SUBROUTINE write_2d_lag_mesh_r4(h, id, name, x, y, dims, xmin, xmax, &
      ymin, ymax, distribution, subarray, convert_in, dim_labels, &
      dim_units, dim_mults, geometry)

    INTEGER, PARAMETER :: ndims = 2
    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name
    REAL(r4), DIMENSION(:,:), INTENT(IN) :: x, y
    INTEGER, DIMENSION(:), INTENT(IN) :: dims
    REAL(r4), INTENT(IN) :: xmin, xmax, ymin, ymax
    INTEGER, INTENT(IN) :: distribution, subarray
    LOGICAL, INTENT(IN), OPTIONAL :: convert_in
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: dim_labels(:), dim_units(:)
    REAL(r4), DIMENSION(:), INTENT(IN), OPTIONAL :: dim_mults
    INTEGER, INTENT(IN), OPTIONAL :: geometry
    REAL(r8), DIMENSION(ndims) :: gmn, gmx
    REAL(r4), DIMENSION(:,:), ALLOCATABLE :: r4array
    INTEGER :: i, errcode, sz(ndims)
    TYPE(sdf_block_type), POINTER :: b
    LOGICAL :: convert

    CALL sdf_get_next_block(h)
    b => h%current_block

    IF (PRESENT(convert_in)) THEN
      convert = convert_in
    ELSE
      convert = .FALSE.
    END IF

    IF (convert) THEN
      b%type_size = 4
      b%datatype = c_datatype_real4
      b%mpitype = MPI_REAL4

      sz(1) = INT(SIZE(x,1),i4)
      sz(2) = INT(SIZE(x,2),i4)

      ALLOCATE(r4array(sz(1),sz(2)))
    ELSE
      b%type_size = sof
      b%datatype = datatype_real
      b%mpitype = mpitype_real
    END IF

    IF (PRESENT(geometry)) THEN
      b%geometry = geometry
    ELSE
      b%geometry = c_geometry_cartesian
    END IF
    b%ndims = ndims

    DO i = 1,ndims
      b%dims(i) = INT(dims(i),i4)
    END DO

    gmn(1) = xmin
    gmx(1) = xmax
    gmn(2) = ymin
    gmx(2) = ymax

    CALL MPI_ALLREDUCE(gmn, b%extents, ndims, MPI_REAL8, MPI_MIN, &
        h%comm, errcode)
    CALL MPI_ALLREDUCE(gmx, b%extents(ndims+1), ndims, MPI_REAL8, MPI_MAX, &
        h%comm, errcode)

    ! Write header

    b%blocktype = c_blocktype_lagrangian_mesh
    CALL write_mesh_meta_r4(h, id, name, dim_labels, dim_units, dim_mults)

    ! Write the actual data

    CALL MPI_FILE_SET_VIEW(h%filehandle, h%current_location, MPI_BYTE, &
        distribution, 'native', MPI_INFO_NULL, errcode)
    IF (convert) THEN
      r4array = REAL(x,r4)
      CALL MPI_FILE_WRITE_ALL(h%filehandle, r4array, 1, subarray, &
          MPI_STATUS_IGNORE, errcode)

      r4array = REAL(y,r4)
      CALL MPI_FILE_WRITE_ALL(h%filehandle, r4array, 1, subarray, &
          MPI_STATUS_IGNORE, errcode)

      DEALLOCATE(r4array)
    ELSE
      CALL MPI_FILE_WRITE_ALL(h%filehandle, x, 1, subarray, &
          MPI_STATUS_IGNORE, errcode)

      CALL MPI_FILE_WRITE_ALL(h%filehandle, y, 1, subarray, &
          MPI_STATUS_IGNORE, errcode)
    END IF

    CALL MPI_FILE_SET_VIEW(h%filehandle, c_off0, MPI_BYTE, MPI_BYTE, 'native', &
        MPI_INFO_NULL, errcode)

    h%current_location = b%data_location + b%data_length
    b%done_data = .TRUE.

  END SUBROUTINE write_2d_lag_mesh_r4



  !----------------------------------------------------------------------------
  ! Code to write a 3D lagrangian mesh in parallel using the
  ! mpitype {distribution} for distribution of data
  ! It's up to the coder to design the distribution parallel operation, so
  ! need global dims
  !----------------------------------------------------------------------------

  SUBROUTINE write_3d_lag_mesh_r4(h, id, name, x, y, z, dims, xmin, xmax, &
      ymin, ymax, zmin, zmax, distribution, subarray, convert_in, dim_labels, &
      dim_units, dim_mults, geometry)

    INTEGER, PARAMETER :: ndims = 3
    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name
    REAL(r4), DIMENSION(:,:,:), INTENT(IN) :: x, y, z
    INTEGER, DIMENSION(:), INTENT(IN) :: dims
    REAL(r4), INTENT(IN) :: xmin, xmax, ymin, ymax, zmin, zmax
    INTEGER, INTENT(IN) :: distribution, subarray
    LOGICAL, INTENT(IN), OPTIONAL :: convert_in
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: dim_labels(:), dim_units(:)
    REAL(r4), DIMENSION(:), INTENT(IN), OPTIONAL :: dim_mults
    INTEGER, INTENT(IN), OPTIONAL :: geometry
    REAL(r8), DIMENSION(ndims) :: gmn, gmx
    REAL(r4), DIMENSION(:,:,:), ALLOCATABLE :: r4array
    INTEGER :: i, errcode, sz(ndims)
    TYPE(sdf_block_type), POINTER :: b
    LOGICAL :: convert

    CALL sdf_get_next_block(h)
    b => h%current_block

    IF (PRESENT(convert_in)) THEN
      convert = convert_in
    ELSE
      convert = .FALSE.
    END IF

    IF (convert) THEN
      b%type_size = 4
      b%datatype = c_datatype_real4
      b%mpitype = MPI_REAL4

      sz(1) = INT(SIZE(x,1),i4)
      sz(2) = INT(SIZE(x,2),i4)
      sz(3) = INT(SIZE(x,3),i4)

      ALLOCATE(r4array(sz(1),sz(2),sz(3)))
    ELSE
      b%type_size = sof
      b%datatype = datatype_real
      b%mpitype = mpitype_real
    END IF

    IF (PRESENT(geometry)) THEN
      b%geometry = geometry
    ELSE
      b%geometry = c_geometry_cartesian
    END IF
    b%ndims = ndims

    DO i = 1,ndims
      b%dims(i) = INT(dims(i),i4)
    END DO

    gmn(1) = xmin
    gmx(1) = xmax
    gmn(2) = ymin
    gmx(2) = ymax
    gmn(3) = zmin
    gmx(3) = zmax

    CALL MPI_ALLREDUCE(gmn, b%extents, ndims, MPI_REAL8, MPI_MIN, &
        h%comm, errcode)
    CALL MPI_ALLREDUCE(gmx, b%extents(ndims+1), ndims, MPI_REAL8, MPI_MAX, &
        h%comm, errcode)

    ! Write header

    b%blocktype = c_blocktype_lagrangian_mesh
    CALL write_mesh_meta_r4(h, id, name, dim_labels, dim_units, dim_mults)

    ! Write the actual data

    CALL MPI_FILE_SET_VIEW(h%filehandle, h%current_location, MPI_BYTE, &
        distribution, 'native', MPI_INFO_NULL, errcode)
    IF (convert) THEN
      r4array = REAL(x,r4)
      CALL MPI_FILE_WRITE_ALL(h%filehandle, r4array, 1, subarray, &
          MPI_STATUS_IGNORE, errcode)

      r4array = REAL(y,r4)
      CALL MPI_FILE_WRITE_ALL(h%filehandle, r4array, 1, subarray, &
          MPI_STATUS_IGNORE, errcode)

      r4array = REAL(z,r4)
      CALL MPI_FILE_WRITE_ALL(h%filehandle, r4array, 1, subarray, &
          MPI_STATUS_IGNORE, errcode)

      DEALLOCATE(r4array)
    ELSE
      CALL MPI_FILE_WRITE_ALL(h%filehandle, x, 1, subarray, &
          MPI_STATUS_IGNORE, errcode)

      CALL MPI_FILE_WRITE_ALL(h%filehandle, y, 1, subarray, &
          MPI_STATUS_IGNORE, errcode)

      CALL MPI_FILE_WRITE_ALL(h%filehandle, z, 1, subarray, &
          MPI_STATUS_IGNORE, errcode)
    END IF

    CALL MPI_FILE_SET_VIEW(h%filehandle, c_off0, MPI_BYTE, MPI_BYTE, 'native', &
        MPI_INFO_NULL, errcode)

    h%current_location = b%data_location + b%data_length
    b%done_data = .TRUE.

  END SUBROUTINE write_3d_lag_mesh_r4



  !----------------------------------------------------------------------------
  ! Code to write a 1D cartesian variable in serial from the node with
  ! rank {rank_write}
  !----------------------------------------------------------------------------

  SUBROUTINE write_srl_1d_float_gen_r4(h, id, name, units, ndims, nm, dims, &
      sz, stagger, mesh_id, variable, subarray, convert_in, mult, rank_write)

    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name, units
    INTEGER, INTENT(IN) :: ndims, nm
    INTEGER, DIMENSION(:), INTENT(IN) :: dims, sz
    INTEGER(i4), INTENT(IN) :: stagger
    CHARACTER(LEN=*), INTENT(IN) :: mesh_id
    REAL(r4), INTENT(IN) :: variable(sz(1))
    INTEGER, OPTIONAL, INTENT(IN) :: subarray
    LOGICAL, OPTIONAL, INTENT(IN) :: convert_in
    REAL(r4), OPTIONAL, INTENT(IN) :: mult
    INTEGER, OPTIONAL, INTENT(IN) :: rank_write
    REAL(r4), DIMENSION(:), ALLOCATABLE :: r4array
    INTEGER :: i, errcode, nd, mpitype, nitems, nelements
    TYPE(sdf_block_type), POINTER :: b
    LOGICAL :: convert

    CALL sdf_get_next_block(h)
    b => h%current_block

    IF (PRESENT(convert_in)) THEN
      convert = convert_in
    ELSE
      convert = .FALSE.
    END IF

    IF (convert) THEN
      b%type_size = 4
      b%datatype = c_datatype_real4
      b%mpitype = MPI_REAL4
    ELSE
      b%type_size = sof
      b%datatype = datatype_real
      b%mpitype = mpitype_real
    END IF
    b%ndims = ndims
    b%stagger = stagger

    nd = ndims
    IF (nm > 1) THEN
      b%dims(nd) = nm
      nd = nd - 1
    END IF

    nelements = 1
    DO i = 1,nd
      b%dims(i) = INT(dims(i),i4)
      nelements = nelements * b%dims(i)
    END DO

    IF (PRESENT(rank_write)) h%rank_master = rank_write

    ! Write header

    CALL write_mesh_variable_meta_r4(h, id, name, units, mesh_id, mult)

    ! Write the actual data

    IF (h%rank == h%rank_master) THEN
      h%current_location = b%data_location

      CALL MPI_FILE_SEEK(h%filehandle, h%current_location, MPI_SEEK_SET, &
          errcode)

      IF (PRESENT(subarray)) THEN
        mpitype = subarray
        nitems = nm
      ELSE
        mpitype = b%mpitype
        nitems = nelements
      END IF

      IF (convert) THEN
        ALLOCATE(r4array(sz(1)))
        r4array = REAL(variable,r4)
        CALL MPI_FILE_WRITE(h%filehandle, r4array, nitems, mpitype, &
            MPI_STATUS_IGNORE, errcode)
        DEALLOCATE(r4array)
      ELSE
        CALL MPI_FILE_WRITE(h%filehandle, variable, nitems, mpitype, &
            MPI_STATUS_IGNORE, errcode)
      END IF
    END IF

    h%current_location = b%data_location + b%data_length
    b%done_data = .TRUE.

  END SUBROUTINE write_srl_1d_float_gen_r4



  !----------------------------------------------------------------------------
  ! Code to write a 2D cartesian variable in serial from the node with
  ! rank {rank_write}
  !----------------------------------------------------------------------------

  SUBROUTINE write_srl_2d_float_gen_r4(h, id, name, units, ndims, nm, dims, &
      sz, stagger, mesh_id, variable, subarray, convert_in, mult, rank_write)

    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name, units
    INTEGER, INTENT(IN) :: ndims, nm
    INTEGER, DIMENSION(:), INTENT(IN) :: dims, sz
    INTEGER(i4), INTENT(IN) :: stagger
    CHARACTER(LEN=*), INTENT(IN) :: mesh_id
    REAL(r4), INTENT(IN) :: variable(sz(1),sz(2))
    INTEGER, OPTIONAL, INTENT(IN) :: subarray
    LOGICAL, OPTIONAL, INTENT(IN) :: convert_in
    REAL(r4), OPTIONAL, INTENT(IN) :: mult
    INTEGER, OPTIONAL, INTENT(IN) :: rank_write
    REAL(r4), DIMENSION(:,:), ALLOCATABLE :: r4array
    INTEGER :: i, errcode, nd, mpitype, nitems, nelements
    TYPE(sdf_block_type), POINTER :: b
    LOGICAL :: convert

    CALL sdf_get_next_block(h)
    b => h%current_block

    IF (PRESENT(convert_in)) THEN
      convert = convert_in
    ELSE
      convert = .FALSE.
    END IF

    IF (convert) THEN
      b%type_size = 4
      b%datatype = c_datatype_real4
      b%mpitype = MPI_REAL4
    ELSE
      b%type_size = sof
      b%datatype = datatype_real
      b%mpitype = mpitype_real
    END IF
    b%ndims = ndims
    b%stagger = stagger

    nd = ndims
    IF (nm > 1) THEN
      b%dims(nd) = nm
      nd = nd - 1
    END IF

    nelements = 1
    DO i = 1,nd
      b%dims(i) = INT(dims(i),i4)
      nelements = nelements * b%dims(i)
    END DO

    IF (PRESENT(rank_write)) h%rank_master = rank_write

    ! Write header

    CALL write_mesh_variable_meta_r4(h, id, name, units, mesh_id, mult)

    ! Write the actual data

    IF (h%rank == h%rank_master) THEN
      h%current_location = b%data_location

      CALL MPI_FILE_SEEK(h%filehandle, h%current_location, MPI_SEEK_SET, &
          errcode)

      IF (PRESENT(subarray)) THEN
        mpitype = subarray
        nitems = nm
      ELSE
        mpitype = b%mpitype
        nitems = nelements
      END IF

      IF (convert) THEN
        ALLOCATE(r4array(sz(1),sz(2)))
        r4array = REAL(variable,r4)
        CALL MPI_FILE_WRITE(h%filehandle, r4array, nitems, mpitype, &
            MPI_STATUS_IGNORE, errcode)
        DEALLOCATE(r4array)
      ELSE
        CALL MPI_FILE_WRITE(h%filehandle, variable, nitems, mpitype, &
            MPI_STATUS_IGNORE, errcode)
      END IF
    END IF

    h%current_location = b%data_location + b%data_length
    b%done_data = .TRUE.

  END SUBROUTINE write_srl_2d_float_gen_r4



  !----------------------------------------------------------------------------
  ! Code to write a 3D cartesian variable in serial from the node with
  ! rank {rank_write}
  !----------------------------------------------------------------------------

  SUBROUTINE write_srl_3d_float_gen_r4(h, id, name, units, ndims, nm, dims, &
      sz, stagger, mesh_id, variable, subarray, convert_in, mult, rank_write)

    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name, units
    INTEGER, INTENT(IN) :: ndims, nm
    INTEGER, DIMENSION(:), INTENT(IN) :: dims, sz
    INTEGER(i4), INTENT(IN) :: stagger
    CHARACTER(LEN=*), INTENT(IN) :: mesh_id
    REAL(r4), INTENT(IN) :: variable(sz(1),sz(2),sz(3))
    INTEGER, OPTIONAL, INTENT(IN) :: subarray
    LOGICAL, OPTIONAL, INTENT(IN) :: convert_in
    REAL(r4), OPTIONAL, INTENT(IN) :: mult
    INTEGER, OPTIONAL, INTENT(IN) :: rank_write
    REAL(r4), DIMENSION(:,:,:), ALLOCATABLE :: r4array
    INTEGER :: i, errcode, nd, mpitype, nitems, nelements
    TYPE(sdf_block_type), POINTER :: b
    LOGICAL :: convert

    CALL sdf_get_next_block(h)
    b => h%current_block

    IF (PRESENT(convert_in)) THEN
      convert = convert_in
    ELSE
      convert = .FALSE.
    END IF

    IF (convert) THEN
      b%type_size = 4
      b%datatype = c_datatype_real4
      b%mpitype = MPI_REAL4
    ELSE
      b%type_size = sof
      b%datatype = datatype_real
      b%mpitype = mpitype_real
    END IF
    b%ndims = ndims
    b%stagger = stagger

    nd = ndims
    IF (nm > 1) THEN
      b%dims(nd) = nm
      nd = nd - 1
    END IF

    nelements = 1
    DO i = 1,nd
      b%dims(i) = INT(dims(i),i4)
      nelements = nelements * b%dims(i)
    END DO

    IF (PRESENT(rank_write)) h%rank_master = rank_write

    ! Write header

    CALL write_mesh_variable_meta_r4(h, id, name, units, mesh_id, mult)

    ! Write the actual data

    IF (h%rank == h%rank_master) THEN
      h%current_location = b%data_location

      CALL MPI_FILE_SEEK(h%filehandle, h%current_location, MPI_SEEK_SET, &
          errcode)

      IF (PRESENT(subarray)) THEN
        mpitype = subarray
        nitems = nm
      ELSE
        mpitype = b%mpitype
        nitems = nelements
      END IF

      IF (convert) THEN
        ALLOCATE(r4array(sz(1),sz(2),sz(3)))
        r4array = REAL(variable,r4)
        CALL MPI_FILE_WRITE(h%filehandle, r4array, nitems, mpitype, &
            MPI_STATUS_IGNORE, errcode)
        DEALLOCATE(r4array)
      ELSE
        CALL MPI_FILE_WRITE(h%filehandle, variable, nitems, mpitype, &
            MPI_STATUS_IGNORE, errcode)
      END IF
    END IF

    h%current_location = b%data_location + b%data_length
    b%done_data = .TRUE.

  END SUBROUTINE write_srl_3d_float_gen_r4



  !----------------------------------------------------------------------------
  ! Code to write a 4D cartesian variable in serial from the node with
  ! rank {rank_write}
  !----------------------------------------------------------------------------

  SUBROUTINE write_srl_4d_float_gen_r4(h, id, name, units, ndims, nm, dims, &
      sz, stagger, mesh_id, variable, subarray, convert_in, mult, rank_write)

    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name, units
    INTEGER, INTENT(IN) :: ndims, nm
    INTEGER, DIMENSION(:), INTENT(IN) :: dims, sz
    INTEGER(i4), INTENT(IN) :: stagger
    CHARACTER(LEN=*), INTENT(IN) :: mesh_id
    REAL(r4), INTENT(IN) :: variable(sz(1),sz(2),sz(3),sz(4))
    INTEGER, OPTIONAL, INTENT(IN) :: subarray
    LOGICAL, OPTIONAL, INTENT(IN) :: convert_in
    REAL(r4), OPTIONAL, INTENT(IN) :: mult
    INTEGER, OPTIONAL, INTENT(IN) :: rank_write
    REAL(r4), DIMENSION(:,:,:,:), ALLOCATABLE :: r4array
    INTEGER :: i, errcode, nd, mpitype, nitems, nelements
    TYPE(sdf_block_type), POINTER :: b
    LOGICAL :: convert

    CALL sdf_get_next_block(h)
    b => h%current_block

    IF (PRESENT(convert_in)) THEN
      convert = convert_in
    ELSE
      convert = .FALSE.
    END IF

    IF (convert) THEN
      b%type_size = 4
      b%datatype = c_datatype_real4
      b%mpitype = MPI_REAL4
    ELSE
      b%type_size = sof
      b%datatype = datatype_real
      b%mpitype = mpitype_real
    END IF
    b%ndims = ndims
    b%stagger = stagger

    nd = ndims
    IF (nm > 1) THEN
      b%dims(nd) = nm
      nd = nd - 1
    END IF

    nelements = 1
    DO i = 1,nd
      b%dims(i) = INT(dims(i),i4)
      nelements = nelements * b%dims(i)
    END DO

    IF (PRESENT(rank_write)) h%rank_master = rank_write

    ! Write header

    CALL write_mesh_variable_meta_r4(h, id, name, units, mesh_id, mult)

    ! Write the actual data

    IF (h%rank == h%rank_master) THEN
      h%current_location = b%data_location

      CALL MPI_FILE_SEEK(h%filehandle, h%current_location, MPI_SEEK_SET, &
          errcode)

      IF (PRESENT(subarray)) THEN
        mpitype = subarray
        nitems = nm
      ELSE
        mpitype = b%mpitype
        nitems = nelements
      END IF

      IF (convert) THEN
        ALLOCATE(r4array(sz(1),sz(2),sz(3),sz(4)))
        r4array = REAL(variable,r4)
        CALL MPI_FILE_WRITE(h%filehandle, r4array, nitems, mpitype, &
            MPI_STATUS_IGNORE, errcode)
        DEALLOCATE(r4array)
      ELSE
        CALL MPI_FILE_WRITE(h%filehandle, variable, nitems, mpitype, &
            MPI_STATUS_IGNORE, errcode)
      END IF
    END IF

    h%current_location = b%data_location + b%data_length
    b%done_data = .TRUE.

  END SUBROUTINE write_srl_4d_float_gen_r4



  !----------------------------------------------------------------------------
  ! Code to write a 1D cartesian variable in serial from the node with
  ! rank {rank_write}
  !----------------------------------------------------------------------------

  SUBROUTINE write_srl_1d_float_r4(h, id, name, units, dims, stagger, mesh_id, &
      variable, subarray, convert, mult, rank_write)

    INTEGER, PARAMETER :: ndims = 1
    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name, units
    INTEGER, DIMENSION(:), INTENT(IN) :: dims
    INTEGER(i4), INTENT(IN) :: stagger
    CHARACTER(LEN=*), INTENT(IN) :: mesh_id
    REAL(r4), DIMENSION(:), INTENT(IN) :: variable
    INTEGER, OPTIONAL, INTENT(IN) :: subarray
    LOGICAL, OPTIONAL, INTENT(IN) :: convert
    REAL(r4), OPTIONAL, INTENT(IN) :: mult
    INTEGER, OPTIONAL, INTENT(IN) :: rank_write
    INTEGER :: i, sz(ndims)

    DO i = 1,ndims
      sz(i) = SIZE(variable,i)
    END DO

    CALL write_srl_1d_float_gen_r4(h, id, name, units, ndims, 1, dims, sz, &
        stagger, mesh_id, variable, subarray, convert, mult, rank_write)

  END SUBROUTINE write_srl_1d_float_r4



  !----------------------------------------------------------------------------
  ! Code to write a 2D cartesian variable in serial from the node with
  ! rank {rank_write}
  !----------------------------------------------------------------------------

  SUBROUTINE write_srl_2d_float_r4(h, id, name, units, dims, stagger, mesh_id, &
      variable, subarray, convert, mult, rank_write)

    INTEGER, PARAMETER :: ndims = 2
    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name, units
    INTEGER, DIMENSION(:), INTENT(IN) :: dims
    INTEGER(i4), INTENT(IN) :: stagger
    CHARACTER(LEN=*), INTENT(IN) :: mesh_id
    REAL(r4), DIMENSION(:,:), INTENT(IN) :: variable
    INTEGER, OPTIONAL, INTENT(IN) :: subarray
    LOGICAL, OPTIONAL, INTENT(IN) :: convert
    REAL(r4), OPTIONAL, INTENT(IN) :: mult
    INTEGER, OPTIONAL, INTENT(IN) :: rank_write
    INTEGER :: i, sz(ndims)

    DO i = 1,ndims
      sz(i) = SIZE(variable,i)
    END DO

    CALL write_srl_2d_float_gen_r4(h, id, name, units, ndims, 1, dims, sz, &
        stagger, mesh_id, variable, subarray, convert, mult, rank_write)

  END SUBROUTINE write_srl_2d_float_r4



  !----------------------------------------------------------------------------
  ! Code to write a 3D cartesian variable in serial from the node with
  ! rank {rank_write}
  !----------------------------------------------------------------------------

  SUBROUTINE write_srl_3d_float_r4(h, id, name, units, dims, stagger, mesh_id, &
      variable, subarray, convert, mult, rank_write)

    INTEGER, PARAMETER :: ndims = 3
    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name, units
    INTEGER, DIMENSION(:), INTENT(IN) :: dims
    INTEGER(i4), INTENT(IN) :: stagger
    CHARACTER(LEN=*), INTENT(IN) :: mesh_id
    REAL(r4), DIMENSION(:,:,:), INTENT(IN) :: variable
    INTEGER, OPTIONAL, INTENT(IN) :: subarray
    LOGICAL, OPTIONAL, INTENT(IN) :: convert
    REAL(r4), OPTIONAL, INTENT(IN) :: mult
    INTEGER, OPTIONAL, INTENT(IN) :: rank_write
    INTEGER :: i, sz(ndims)

    DO i = 1,ndims
      sz(i) = SIZE(variable,i)
    END DO

    CALL write_srl_3d_float_gen_r4(h, id, name, units, ndims, 1, dims, sz, &
        stagger, mesh_id, variable, subarray, convert, mult, rank_write)

  END SUBROUTINE write_srl_3d_float_r4



  !----------------------------------------------------------------------------
  ! Code to write a 1D cartesian variable in serial from the node with
  ! rank {rank_write}
  !----------------------------------------------------------------------------

  SUBROUTINE write_srl_1d_float_num_r4(h, id, name, units, dims, stagger, &
      nm, mesh_id, variable, subarray, convert, mult, rank_write)

    INTEGER, PARAMETER :: ndims = 1
    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name, units
    INTEGER, DIMENSION(:), INTENT(IN) :: dims
    INTEGER(i4), INTENT(IN) :: stagger
    INTEGER, INTENT(IN) :: nm
    CHARACTER(LEN=*), INTENT(IN) :: mesh_id
    REAL(r4), DIMENSION(:), INTENT(IN) :: variable
    INTEGER, OPTIONAL, INTENT(IN) :: subarray
    LOGICAL, OPTIONAL, INTENT(IN) :: convert
    REAL(r4), OPTIONAL, INTENT(IN) :: mult
    INTEGER, OPTIONAL, INTENT(IN) :: rank_write
    INTEGER :: i, sz(ndims)

    DO i = 1,ndims
      sz(i) = SIZE(variable,i)
    END DO

    CALL write_srl_1d_float_gen_r4(h, id, name, units, ndims, nm, dims, sz, &
        stagger, mesh_id, variable, subarray, convert, mult, rank_write)

  END SUBROUTINE write_srl_1d_float_num_r4



  !----------------------------------------------------------------------------
  ! Code to write a 2D cartesian variable in serial from the node with
  ! rank {rank_write}
  !----------------------------------------------------------------------------

  SUBROUTINE write_srl_2d_float_num_r4(h, id, name, units, dims, stagger, &
      nm, mesh_id, variable, subarray, convert, mult, rank_write)

    INTEGER, PARAMETER :: ndims = 2
    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name, units
    INTEGER, DIMENSION(:), INTENT(IN) :: dims
    INTEGER(i4), INTENT(IN) :: stagger
    INTEGER, INTENT(IN) :: nm
    CHARACTER(LEN=*), INTENT(IN) :: mesh_id
    REAL(r4), DIMENSION(:,:), INTENT(IN) :: variable
    INTEGER, OPTIONAL, INTENT(IN) :: subarray
    LOGICAL, OPTIONAL, INTENT(IN) :: convert
    REAL(r4), OPTIONAL, INTENT(IN) :: mult
    INTEGER, OPTIONAL, INTENT(IN) :: rank_write
    INTEGER :: i, sz(ndims)

    DO i = 1,ndims
      sz(i) = SIZE(variable,i)
    END DO

    CALL write_srl_2d_float_gen_r4(h, id, name, units, ndims, nm, dims, sz, &
        stagger, mesh_id, variable, subarray, convert, mult, rank_write)

  END SUBROUTINE write_srl_2d_float_num_r4



  !----------------------------------------------------------------------------
  ! Code to write a 3D cartesian variable in serial from the node with
  ! rank {rank_write}
  !----------------------------------------------------------------------------

  SUBROUTINE write_srl_3d_float_num_r4(h, id, name, units, dims, stagger, &
      nm, mesh_id, variable, subarray, convert, mult, rank_write)

    INTEGER, PARAMETER :: ndims = 3
    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name, units
    INTEGER, DIMENSION(:), INTENT(IN) :: dims
    INTEGER(i4), INTENT(IN) :: stagger
    INTEGER, INTENT(IN) :: nm
    CHARACTER(LEN=*), INTENT(IN) :: mesh_id
    REAL(r4), DIMENSION(:,:,:), INTENT(IN) :: variable
    INTEGER, OPTIONAL, INTENT(IN) :: subarray
    LOGICAL, OPTIONAL, INTENT(IN) :: convert
    REAL(r4), OPTIONAL, INTENT(IN) :: mult
    INTEGER, OPTIONAL, INTENT(IN) :: rank_write
    INTEGER :: i, sz(ndims)

    DO i = 1,ndims
      sz(i) = SIZE(variable,i)
    END DO

    CALL write_srl_3d_float_gen_r4(h, id, name, units, ndims, nm, dims, sz, &
        stagger, mesh_id, variable, subarray, convert, mult, rank_write)

  END SUBROUTINE write_srl_3d_float_num_r4



  !----------------------------------------------------------------------------
  ! Code to write a 4D cartesian variable in serial from the node with
  ! rank {rank_write}
  !----------------------------------------------------------------------------

  SUBROUTINE write_srl_4d_float_num_r4(h, id, name, units, dims, stagger, &
      nm, mesh_id, variable, subarray, convert, mult, rank_write)

    INTEGER, PARAMETER :: ndims = 4
    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name, units
    INTEGER, DIMENSION(:), INTENT(IN) :: dims
    INTEGER(i4), INTENT(IN) :: stagger
    INTEGER, INTENT(IN) :: nm
    CHARACTER(LEN=*), INTENT(IN) :: mesh_id
    REAL(r4), DIMENSION(:,:,:,:), INTENT(IN) :: variable
    INTEGER, OPTIONAL, INTENT(IN) :: subarray
    LOGICAL, OPTIONAL, INTENT(IN) :: convert
    REAL(r4), OPTIONAL, INTENT(IN) :: mult
    INTEGER, OPTIONAL, INTENT(IN) :: rank_write
    INTEGER :: i, sz(ndims)

    DO i = 1,ndims
      sz(i) = SIZE(variable,i)
    END DO

    CALL write_srl_4d_float_gen_r4(h, id, name, units, ndims, nm, dims, sz, &
        stagger, mesh_id, variable, subarray, convert, mult, rank_write)

  END SUBROUTINE write_srl_4d_float_num_r4



  !----------------------------------------------------------------------------
  ! Code to write a 1D cartesian variable in parallel
  ! using the mpitype {distribution} for distribution of data
  ! It's up to the coder to design the distribution parallel operation, so
  ! need global dims
  !----------------------------------------------------------------------------

  SUBROUTINE write_1d_float_gen_r4(h, id, name, units, ndims, nm, dims, sz, &
      stagger, mesh_id, variable, distribution, subarray, convert_in, mult)

    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name, units
    INTEGER, INTENT(IN) :: ndims, nm
    INTEGER, DIMENSION(:), INTENT(IN) :: dims, sz
    INTEGER(i4), INTENT(IN) :: stagger
    CHARACTER(LEN=*), INTENT(IN) :: mesh_id
    REAL(r4), INTENT(IN) :: variable(sz(1))
    INTEGER, INTENT(IN) :: distribution, subarray
    LOGICAL, OPTIONAL, INTENT(IN) :: convert_in
    REAL(r4), OPTIONAL, INTENT(IN) :: mult
    REAL(r4), DIMENSION(:), ALLOCATABLE :: r4array
    INTEGER :: i, errcode, nd
    TYPE(sdf_block_type), POINTER :: b
    LOGICAL :: convert

    CALL sdf_get_next_block(h)
    b => h%current_block

    IF (PRESENT(convert_in)) THEN
      convert = convert_in
    ELSE
      convert = .FALSE.
    END IF

    IF (convert) THEN
      b%type_size = 4
      b%datatype = c_datatype_real4
      b%mpitype = MPI_REAL4
    ELSE
      b%type_size = sof
      b%datatype = datatype_real
      b%mpitype = mpitype_real
    END IF
    b%ndims = ndims
    b%stagger = stagger

    nd = ndims
    IF (nm > 1) THEN
      b%dims(nd) = nm
      nd = nd - 1
    END IF

    DO i = 1,nd
      b%dims(i) = INT(dims(i),i4)
    END DO

    ! Write header

    CALL write_mesh_variable_meta_r4(h, id, name, units, mesh_id, mult)

    ! Write the actual data

    h%current_location = b%data_location
    CALL MPI_FILE_SET_VIEW(h%filehandle, h%current_location, MPI_BYTE, &
        distribution, 'native', MPI_INFO_NULL, errcode)
    IF (convert) THEN
      ALLOCATE(r4array(sz(1)))
      r4array = REAL(variable,r4)
      CALL MPI_FILE_WRITE_ALL(h%filehandle, r4array, nm, subarray, &
          MPI_STATUS_IGNORE, errcode)
      DEALLOCATE(r4array)
    ELSE
      CALL MPI_FILE_WRITE_ALL(h%filehandle, variable, nm, subarray, &
          MPI_STATUS_IGNORE, errcode)
    END IF

    CALL MPI_FILE_SET_VIEW(h%filehandle, c_off0, MPI_BYTE, MPI_BYTE, 'native', &
        MPI_INFO_NULL, errcode)

    h%current_location = b%data_location + b%data_length
    b%done_data = .TRUE.

  END SUBROUTINE write_1d_float_gen_r4



  !----------------------------------------------------------------------------
  ! Code to write a 2D cartesian variable in parallel
  ! using the mpitype {distribution} for distribution of data
  ! It's up to the coder to design the distribution parallel operation, so
  ! need global dims
  !----------------------------------------------------------------------------

  SUBROUTINE write_2d_float_gen_r4(h, id, name, units, ndims, nm, dims, sz, &
      stagger, mesh_id, variable, distribution, subarray, convert_in, mult)

    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name, units
    INTEGER, INTENT(IN) :: ndims, nm
    INTEGER, DIMENSION(:), INTENT(IN) :: dims, sz
    INTEGER(i4), INTENT(IN) :: stagger
    CHARACTER(LEN=*), INTENT(IN) :: mesh_id
    REAL(r4), INTENT(IN) :: variable(sz(1),sz(2))
    INTEGER, INTENT(IN) :: distribution, subarray
    LOGICAL, OPTIONAL, INTENT(IN) :: convert_in
    REAL(r4), OPTIONAL, INTENT(IN) :: mult
    REAL(r4), DIMENSION(:,:), ALLOCATABLE :: r4array
    INTEGER :: i, errcode, nd
    TYPE(sdf_block_type), POINTER :: b
    LOGICAL :: convert

    CALL sdf_get_next_block(h)
    b => h%current_block

    IF (PRESENT(convert_in)) THEN
      convert = convert_in
    ELSE
      convert = .FALSE.
    END IF

    IF (convert) THEN
      b%type_size = 4
      b%datatype = c_datatype_real4
      b%mpitype = MPI_REAL4
    ELSE
      b%type_size = sof
      b%datatype = datatype_real
      b%mpitype = mpitype_real
    END IF
    b%ndims = ndims
    b%stagger = stagger

    nd = ndims
    IF (nm > 1) THEN
      b%dims(nd) = nm
      nd = nd - 1
    END IF

    DO i = 1,nd
      b%dims(i) = INT(dims(i),i4)
    END DO

    ! Write header

    CALL write_mesh_variable_meta_r4(h, id, name, units, mesh_id, mult)

    ! Write the actual data

    h%current_location = b%data_location
    CALL MPI_FILE_SET_VIEW(h%filehandle, h%current_location, MPI_BYTE, &
        distribution, 'native', MPI_INFO_NULL, errcode)
    IF (convert) THEN
      ALLOCATE(r4array(sz(1),sz(2)))
      r4array = REAL(variable,r4)
      CALL MPI_FILE_WRITE_ALL(h%filehandle, r4array, nm, subarray, &
          MPI_STATUS_IGNORE, errcode)
      DEALLOCATE(r4array)
    ELSE
      CALL MPI_FILE_WRITE_ALL(h%filehandle, variable, nm, subarray, &
          MPI_STATUS_IGNORE, errcode)
    END IF

    CALL MPI_FILE_SET_VIEW(h%filehandle, c_off0, MPI_BYTE, MPI_BYTE, 'native', &
        MPI_INFO_NULL, errcode)

    h%current_location = b%data_location + b%data_length
    b%done_data = .TRUE.

  END SUBROUTINE write_2d_float_gen_r4



  !----------------------------------------------------------------------------
  ! Code to write a 3D cartesian variable in parallel
  ! using the mpitype {distribution} for distribution of data
  ! It's up to the coder to design the distribution parallel operation, so
  ! need global dims
  !----------------------------------------------------------------------------

  SUBROUTINE write_3d_float_gen_r4(h, id, name, units, ndims, nm, dims, sz, &
      stagger, mesh_id, variable, distribution, subarray, convert_in, mult)

    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name, units
    INTEGER, INTENT(IN) :: ndims, nm
    INTEGER, DIMENSION(:), INTENT(IN) :: dims, sz
    INTEGER(i4), INTENT(IN) :: stagger
    CHARACTER(LEN=*), INTENT(IN) :: mesh_id
    REAL(r4), INTENT(IN) :: variable(sz(1),sz(2),sz(3))
    INTEGER, INTENT(IN) :: distribution, subarray
    LOGICAL, OPTIONAL, INTENT(IN) :: convert_in
    REAL(r4), OPTIONAL, INTENT(IN) :: mult
    REAL(r4), DIMENSION(:,:,:), ALLOCATABLE :: r4array
    INTEGER :: i, errcode, nd
    TYPE(sdf_block_type), POINTER :: b
    LOGICAL :: convert

    CALL sdf_get_next_block(h)
    b => h%current_block

    IF (PRESENT(convert_in)) THEN
      convert = convert_in
    ELSE
      convert = .FALSE.
    END IF

    IF (convert) THEN
      b%type_size = 4
      b%datatype = c_datatype_real4
      b%mpitype = MPI_REAL4
    ELSE
      b%type_size = sof
      b%datatype = datatype_real
      b%mpitype = mpitype_real
    END IF
    b%ndims = ndims
    b%stagger = stagger

    nd = ndims
    IF (nm > 1) THEN
      b%dims(nd) = nm
      nd = nd - 1
    END IF

    DO i = 1,nd
      b%dims(i) = INT(dims(i),i4)
    END DO

    ! Write header

    CALL write_mesh_variable_meta_r4(h, id, name, units, mesh_id, mult)

    ! Write the actual data

    h%current_location = b%data_location
    CALL MPI_FILE_SET_VIEW(h%filehandle, h%current_location, MPI_BYTE, &
        distribution, 'native', MPI_INFO_NULL, errcode)
    IF (convert) THEN
      ALLOCATE(r4array(sz(1),sz(2),sz(3)))
      r4array = REAL(variable,r4)
      CALL MPI_FILE_WRITE_ALL(h%filehandle, r4array, nm, subarray, &
          MPI_STATUS_IGNORE, errcode)
      DEALLOCATE(r4array)
    ELSE
      CALL MPI_FILE_WRITE_ALL(h%filehandle, variable, nm, subarray, &
          MPI_STATUS_IGNORE, errcode)
    END IF

    CALL MPI_FILE_SET_VIEW(h%filehandle, c_off0, MPI_BYTE, MPI_BYTE, 'native', &
        MPI_INFO_NULL, errcode)

    h%current_location = b%data_location + b%data_length
    b%done_data = .TRUE.

  END SUBROUTINE write_3d_float_gen_r4



  !----------------------------------------------------------------------------
  ! Code to write a 4D cartesian variable in parallel
  ! using the mpitype {distribution} for distribution of data
  ! It's up to the coder to design the distribution parallel operation, so
  ! need global dims
  !----------------------------------------------------------------------------

  SUBROUTINE write_4d_float_gen_r4(h, id, name, units, ndims, nm, dims, sz, &
      stagger, mesh_id, variable, distribution, subarray, convert_in, mult)

    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name, units
    INTEGER, INTENT(IN) :: ndims, nm
    INTEGER, DIMENSION(:), INTENT(IN) :: dims, sz
    INTEGER(i4), INTENT(IN) :: stagger
    CHARACTER(LEN=*), INTENT(IN) :: mesh_id
    REAL(r4), INTENT(IN) :: variable(sz(1),sz(2),sz(3),sz(4))
    INTEGER, INTENT(IN) :: distribution, subarray
    LOGICAL, OPTIONAL, INTENT(IN) :: convert_in
    REAL(r4), OPTIONAL, INTENT(IN) :: mult
    REAL(r4), DIMENSION(:,:,:,:), ALLOCATABLE :: r4array
    INTEGER :: i, errcode, nd
    TYPE(sdf_block_type), POINTER :: b
    LOGICAL :: convert

    CALL sdf_get_next_block(h)
    b => h%current_block

    IF (PRESENT(convert_in)) THEN
      convert = convert_in
    ELSE
      convert = .FALSE.
    END IF

    IF (convert) THEN
      b%type_size = 4
      b%datatype = c_datatype_real4
      b%mpitype = MPI_REAL4
    ELSE
      b%type_size = sof
      b%datatype = datatype_real
      b%mpitype = mpitype_real
    END IF
    b%ndims = ndims
    b%stagger = stagger

    nd = ndims
    IF (nm > 1) THEN
      b%dims(nd) = nm
      nd = nd - 1
    END IF

    DO i = 1,nd
      b%dims(i) = INT(dims(i),i4)
    END DO

    ! Write header

    CALL write_mesh_variable_meta_r4(h, id, name, units, mesh_id, mult)

    ! Write the actual data

    h%current_location = b%data_location
    CALL MPI_FILE_SET_VIEW(h%filehandle, h%current_location, MPI_BYTE, &
        distribution, 'native', MPI_INFO_NULL, errcode)
    IF (convert) THEN
      ALLOCATE(r4array(sz(1),sz(2),sz(3),sz(4)))
      r4array = REAL(variable,r4)
      CALL MPI_FILE_WRITE_ALL(h%filehandle, r4array, nm, subarray, &
          MPI_STATUS_IGNORE, errcode)
      DEALLOCATE(r4array)
    ELSE
      CALL MPI_FILE_WRITE_ALL(h%filehandle, variable, nm, subarray, &
          MPI_STATUS_IGNORE, errcode)
    END IF

    CALL MPI_FILE_SET_VIEW(h%filehandle, c_off0, MPI_BYTE, MPI_BYTE, 'native', &
        MPI_INFO_NULL, errcode)

    h%current_location = b%data_location + b%data_length
    b%done_data = .TRUE.

  END SUBROUTINE write_4d_float_gen_r4



  !----------------------------------------------------------------------------
  ! Code to write a 1D cartesian variable in parallel
  ! using the mpitype {distribution} for distribution of data
  ! It's up to the coder to design the distribution parallel operation, so
  ! need global dims
  !----------------------------------------------------------------------------

  SUBROUTINE write_1d_float_r4(h, id, name, units, dims, stagger, mesh_id, &
      variable, distribution, subarray, convert, mult)

    INTEGER, PARAMETER :: ndims = 1
    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name, units
    INTEGER, DIMENSION(:), INTENT(IN) :: dims
    INTEGER(i4), INTENT(IN) :: stagger
    CHARACTER(LEN=*), INTENT(IN) :: mesh_id
    REAL(r4), DIMENSION(:), INTENT(IN) :: variable
    INTEGER, INTENT(IN) :: distribution, subarray
    LOGICAL, OPTIONAL, INTENT(IN) :: convert
    REAL(r4), OPTIONAL, INTENT(IN) :: mult
    INTEGER :: i, sz(ndims)

    DO i = 1,ndims
      sz(i) = SIZE(variable,i)
    END DO

    CALL write_1d_float_gen_r4(h, id, name, units, ndims, 1, dims, sz, &
        stagger, mesh_id, variable, distribution, subarray, convert, mult)

  END SUBROUTINE write_1d_float_r4



  !----------------------------------------------------------------------------
  ! Code to write a 2D cartesian variable in parallel
  ! using the mpitype {distribution} for distribution of data
  ! It's up to the coder to design the distribution parallel operation, so
  ! need global dims
  !----------------------------------------------------------------------------

  SUBROUTINE write_2d_float_r4(h, id, name, units, dims, stagger, mesh_id, &
      variable, distribution, subarray, convert, mult)

    INTEGER, PARAMETER :: ndims = 2
    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name, units
    INTEGER, DIMENSION(:), INTENT(IN) :: dims
    INTEGER(i4), INTENT(IN) :: stagger
    CHARACTER(LEN=*), INTENT(IN) :: mesh_id
    REAL(r4), DIMENSION(:,:), INTENT(IN) :: variable
    INTEGER, INTENT(IN) :: distribution, subarray
    LOGICAL, OPTIONAL, INTENT(IN) :: convert
    REAL(r4), OPTIONAL, INTENT(IN) :: mult
    INTEGER :: i, sz(ndims)

    DO i = 1,ndims
      sz(i) = SIZE(variable,i)
    END DO

    CALL write_2d_float_gen_r4(h, id, name, units, ndims, 1, dims, sz, &
        stagger, mesh_id, variable, distribution, subarray, convert, mult)

  END SUBROUTINE write_2d_float_r4



  !----------------------------------------------------------------------------
  ! Code to write a 3D cartesian variable in parallel
  ! using the mpitype {distribution} for distribution of data
  ! It's up to the coder to design the distribution parallel operation, so
  ! need global dims
  !----------------------------------------------------------------------------

  SUBROUTINE write_3d_float_r4(h, id, name, units, dims, stagger, mesh_id, &
      variable, distribution, subarray, convert, mult)

    INTEGER, PARAMETER :: ndims = 3
    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name, units
    INTEGER, DIMENSION(:), INTENT(IN) :: dims
    INTEGER(i4), INTENT(IN) :: stagger
    CHARACTER(LEN=*), INTENT(IN) :: mesh_id
    REAL(r4), DIMENSION(:,:,:), INTENT(IN) :: variable
    INTEGER, INTENT(IN) :: distribution, subarray
    LOGICAL, OPTIONAL, INTENT(IN) :: convert
    REAL(r4), OPTIONAL, INTENT(IN) :: mult
    INTEGER :: i, sz(ndims)

    DO i = 1,ndims
      sz(i) = SIZE(variable,i)
    END DO

    CALL write_3d_float_gen_r4(h, id, name, units, ndims, 1, dims, sz, &
        stagger, mesh_id, variable, distribution, subarray, convert, mult)

  END SUBROUTINE write_3d_float_r4



  !----------------------------------------------------------------------------
  ! Code to write a 1D cartesian variable in parallel
  ! using the mpitype {distribution} for distribution of data
  ! It's up to the coder to design the distribution parallel operation, so
  ! need global dims
  !----------------------------------------------------------------------------

  SUBROUTINE write_1d_float_num_r4(h, id, name, units, dims, nm, stagger, &
      mesh_id, variable, distribution, subarray, convert, mult)

    INTEGER, PARAMETER :: ndims = 1
    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name, units
    INTEGER, DIMENSION(:), INTENT(IN) :: dims
    INTEGER, INTENT(IN) :: nm
    INTEGER(i4), INTENT(IN) :: stagger
    CHARACTER(LEN=*), INTENT(IN) :: mesh_id
    REAL(r4), DIMENSION(:), INTENT(IN) :: variable
    INTEGER, INTENT(IN) :: distribution, subarray
    LOGICAL, OPTIONAL, INTENT(IN) :: convert
    REAL(r4), OPTIONAL, INTENT(IN) :: mult
    INTEGER :: i, sz(ndims)

    DO i = 1,ndims
      sz(i) = SIZE(variable,i)
    END DO

    CALL write_1d_float_gen_r4(h, id, name, units, ndims, nm, dims, sz, &
        stagger, mesh_id, variable, distribution, subarray, convert, mult)

  END SUBROUTINE write_1d_float_num_r4



  !----------------------------------------------------------------------------
  ! Code to write a 2D cartesian variable in parallel
  ! using the mpitype {distribution} for distribution of data
  ! It's up to the coder to design the distribution parallel operation, so
  ! need global dims
  !----------------------------------------------------------------------------

  SUBROUTINE write_2d_float_num_r4(h, id, name, units, dims, nm, stagger, &
      mesh_id, variable, distribution, subarray, convert, mult)

    INTEGER, PARAMETER :: ndims = 2
    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name, units
    INTEGER, DIMENSION(:), INTENT(IN) :: dims
    INTEGER, INTENT(IN) :: nm
    INTEGER(i4), INTENT(IN) :: stagger
    CHARACTER(LEN=*), INTENT(IN) :: mesh_id
    REAL(r4), DIMENSION(:,:), INTENT(IN) :: variable
    INTEGER, INTENT(IN) :: distribution, subarray
    LOGICAL, OPTIONAL, INTENT(IN) :: convert
    REAL(r4), OPTIONAL, INTENT(IN) :: mult
    INTEGER :: i, sz(ndims)

    DO i = 1,ndims
      sz(i) = SIZE(variable,i)
    END DO

    CALL write_2d_float_gen_r4(h, id, name, units, ndims, nm, dims, sz, &
        stagger, mesh_id, variable, distribution, subarray, convert, mult)

  END SUBROUTINE write_2d_float_num_r4



  !----------------------------------------------------------------------------
  ! Code to write a 3D cartesian variable in parallel
  ! using the mpitype {distribution} for distribution of data
  ! It's up to the coder to design the distribution parallel operation, so
  ! need global dims
  !----------------------------------------------------------------------------

  SUBROUTINE write_3d_float_num_r4(h, id, name, units, dims, nm, stagger, &
      mesh_id, variable, distribution, subarray, convert, mult)

    INTEGER, PARAMETER :: ndims = 3
    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name, units
    INTEGER, DIMENSION(:), INTENT(IN) :: dims
    INTEGER, INTENT(IN) :: nm
    INTEGER(i4), INTENT(IN) :: stagger
    CHARACTER(LEN=*), INTENT(IN) :: mesh_id
    REAL(r4), DIMENSION(:,:,:), INTENT(IN) :: variable
    INTEGER, INTENT(IN) :: distribution, subarray
    LOGICAL, OPTIONAL, INTENT(IN) :: convert
    REAL(r4), OPTIONAL, INTENT(IN) :: mult
    INTEGER :: i, sz(ndims)

    DO i = 1,ndims
      sz(i) = SIZE(variable,i)
    END DO

    CALL write_3d_float_gen_r4(h, id, name, units, ndims, nm, dims, sz, &
        stagger, mesh_id, variable, distribution, subarray, convert, mult)

  END SUBROUTINE write_3d_float_num_r4



  !----------------------------------------------------------------------------
  ! Code to write a 4D cartesian variable in parallel
  ! using the mpitype {distribution} for distribution of data
  ! It's up to the coder to design the distribution parallel operation, so
  ! need global dims
  !----------------------------------------------------------------------------

  SUBROUTINE write_4d_float_num_r4(h, id, name, units, dims, nm, stagger, &
      mesh_id, variable, distribution, subarray, convert, mult)

    INTEGER, PARAMETER :: ndims = 4
    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name, units
    INTEGER, DIMENSION(:), INTENT(IN) :: dims
    INTEGER, INTENT(IN) :: nm
    INTEGER(i4), INTENT(IN) :: stagger
    CHARACTER(LEN=*), INTENT(IN) :: mesh_id
    REAL(r4), DIMENSION(:,:,:,:), INTENT(IN) :: variable
    INTEGER, INTENT(IN) :: distribution, subarray
    LOGICAL, OPTIONAL, INTENT(IN) :: convert
    REAL(r4), OPTIONAL, INTENT(IN) :: mult
    INTEGER :: i, sz(ndims)

    DO i = 1,ndims
      sz(i) = SIZE(variable,i)
    END DO

    CALL write_4d_float_gen_r4(h, id, name, units, ndims, nm, dims, sz, &
        stagger, mesh_id, variable, distribution, subarray, convert, mult)

  END SUBROUTINE write_4d_float_num_r4



  !----------------------------------------------------------------------------
  ! Code to write a 1D cartesian variable component (with material index first)
  ! in parallel using the
  ! mpitype {distribution} for distribution of data
  ! It's up to the coder to design the distribution parallel operation, so
  ! need global dims
  !----------------------------------------------------------------------------

  SUBROUTINE write_1d_var_first_r4(h, id, name, units, dims, sz, stagger, &
      mesh_id, variable, idx, distribution, subarray, convert, mult)

    INTEGER, PARAMETER :: ndims = 1
    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name, units
    INTEGER, DIMENSION(:), INTENT(IN) :: dims, sz
    INTEGER(i4), INTENT(IN) :: stagger
    CHARACTER(LEN=*), INTENT(IN) :: mesh_id
    REAL(r4), INTENT(IN) :: variable(sz(1),sz(2))
    INTEGER, INTENT(IN) :: idx, distribution, subarray
    LOGICAL, OPTIONAL, INTENT(IN) :: convert
    REAL(r4), OPTIONAL, INTENT(IN) :: mult

    CALL write_2d_float_gen_r4(h, id, name, units, ndims, 1, dims, sz, &
        stagger, mesh_id, variable(idx,1), distribution, subarray, &
        convert, mult)

  END SUBROUTINE write_1d_var_first_r4



  !----------------------------------------------------------------------------
  ! Code to write a 2D cartesian variable component (with material index first)
  ! in parallel using the
  ! mpitype {distribution} for distribution of data
  ! It's up to the coder to design the distribution parallel operation, so
  ! need global dims
  !----------------------------------------------------------------------------

  SUBROUTINE write_2d_var_first_r4(h, id, name, units, dims, sz, stagger, &
      mesh_id, variable, idx, distribution, subarray, convert, mult)

    INTEGER, PARAMETER :: ndims = 2
    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name, units
    INTEGER, DIMENSION(:), INTENT(IN) :: dims, sz
    INTEGER(i4), INTENT(IN) :: stagger
    CHARACTER(LEN=*), INTENT(IN) :: mesh_id
    REAL(r4), INTENT(IN) :: variable(sz(1),sz(2),sz(3))
    INTEGER, INTENT(IN) :: idx, distribution, subarray
    LOGICAL, OPTIONAL, INTENT(IN) :: convert
    REAL(r4), OPTIONAL, INTENT(IN) :: mult

    CALL write_3d_float_gen_r4(h, id, name, units, ndims, 1, dims, sz, &
        stagger, mesh_id, variable(idx,1,1), distribution, subarray, &
        convert, mult)

  END SUBROUTINE write_2d_var_first_r4



  !----------------------------------------------------------------------------
  ! Code to write a 3D cartesian variable component (with material index first)
  ! in parallel using the
  ! mpitype {distribution} for distribution of data
  ! It's up to the coder to design the distribution parallel operation, so
  ! need global dims
  !----------------------------------------------------------------------------

  SUBROUTINE write_3d_var_first_r4(h, id, name, units, dims, sz, stagger, &
      mesh_id, variable, idx, distribution, subarray, convert, mult)

    INTEGER, PARAMETER :: ndims = 3
    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name, units
    INTEGER, DIMENSION(:), INTENT(IN) :: dims, sz
    INTEGER(i4), INTENT(IN) :: stagger
    CHARACTER(LEN=*), INTENT(IN) :: mesh_id
    REAL(r4), INTENT(IN) :: variable(sz(1),sz(2),sz(3),sz(4))
    INTEGER, INTENT(IN) :: idx, distribution, subarray
    LOGICAL, OPTIONAL, INTENT(IN) :: convert
    REAL(r4), OPTIONAL, INTENT(IN) :: mult

    CALL write_4d_float_gen_r4(h, id, name, units, ndims, 1, dims, sz, &
        stagger, mesh_id, variable(idx,1,1,1), distribution, subarray, &
        convert, mult)

  END SUBROUTINE write_3d_var_first_r4



  !----------------------------------------------------------------------------
  ! Code to write a 1D cartesian variable component (with material index last)
  ! in parallel using the
  ! mpitype {distribution} for distribution of data
  ! It's up to the coder to design the distribution parallel operation, so
  ! need global dims
  !----------------------------------------------------------------------------

  SUBROUTINE write_1d_var_last_r4(h, id, name, units, dims, sz, stagger, &
      mesh_id, variable, idx, distribution, subarray, convert, mult)

    INTEGER, PARAMETER :: ndims = 1
    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name, units
    INTEGER, DIMENSION(:), INTENT(IN) :: dims, sz
    INTEGER(i4), INTENT(IN) :: stagger
    CHARACTER(LEN=*), INTENT(IN) :: mesh_id
    REAL(r4), INTENT(IN) :: variable(sz(1),sz(2))
    INTEGER, INTENT(IN) :: idx, distribution, subarray
    LOGICAL, OPTIONAL, INTENT(IN) :: convert
    REAL(r4), OPTIONAL, INTENT(IN) :: mult

    CALL write_2d_float_gen_r4(h, id, name, units, ndims, 1, dims, sz, &
        stagger, mesh_id, variable(1,idx), distribution, subarray, &
        convert, mult)

  END SUBROUTINE write_1d_var_last_r4



  !----------------------------------------------------------------------------
  ! Code to write a 2D cartesian variable component (with material index last)
  ! in parallel using the
  ! mpitype {distribution} for distribution of data
  ! It's up to the coder to design the distribution parallel operation, so
  ! need global dims
  !----------------------------------------------------------------------------

  SUBROUTINE write_2d_var_last_r4(h, id, name, units, dims, sz, stagger, &
      mesh_id, variable, idx, distribution, subarray, convert, mult)

    INTEGER, PARAMETER :: ndims = 2
    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name, units
    INTEGER, DIMENSION(:), INTENT(IN) :: dims, sz
    INTEGER(i4), INTENT(IN) :: stagger
    CHARACTER(LEN=*), INTENT(IN) :: mesh_id
    REAL(r4), INTENT(IN) :: variable(sz(1),sz(2),sz(3))
    INTEGER, INTENT(IN) :: idx, distribution, subarray
    LOGICAL, OPTIONAL, INTENT(IN) :: convert
    REAL(r4), OPTIONAL, INTENT(IN) :: mult

    CALL write_3d_float_gen_r4(h, id, name, units, ndims, 1, dims, sz, &
        stagger, mesh_id, variable(1,1,idx), distribution, subarray, &
        convert, mult)

  END SUBROUTINE write_2d_var_last_r4



  !----------------------------------------------------------------------------
  ! Code to write a 3D cartesian variable component (with material index last)
  ! in parallel using the
  ! mpitype {distribution} for distribution of data
  ! It's up to the coder to design the distribution parallel operation, so
  ! need global dims
  !----------------------------------------------------------------------------

  SUBROUTINE write_3d_var_last_r4(h, id, name, units, dims, sz, stagger, &
      mesh_id, variable, idx, distribution, subarray, convert, mult)

    INTEGER, PARAMETER :: ndims = 3
    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name, units
    INTEGER, DIMENSION(:), INTENT(IN) :: dims, sz
    INTEGER(i4), INTENT(IN) :: stagger
    CHARACTER(LEN=*), INTENT(IN) :: mesh_id
    REAL(r4), INTENT(IN) :: variable(sz(1),sz(2),sz(3),sz(4))
    INTEGER, INTENT(IN) :: idx, distribution, subarray
    LOGICAL, OPTIONAL, INTENT(IN) :: convert
    REAL(r4), OPTIONAL, INTENT(IN) :: mult

    CALL write_4d_float_gen_r4(h, id, name, units, ndims, 1, dims, sz, &
        stagger, mesh_id, variable(1,1,1,idx), distribution, subarray, &
        convert, mult)

  END SUBROUTINE write_3d_var_last_r4



  !----------------------------------------------------------------------------
  ! Code to write a 1D cartesian material mesh in parallel using the
  ! mpitype {distribution} for distribution of data
  ! It's up to the coder to design the distribution parallel operation, so
  ! need global dims
  !----------------------------------------------------------------------------

  SUBROUTINE write_1d_material_r4(h, id, name, units, dims, nmat, stagger, &
      mesh_id, material_names, variable, distribution, subarray, convert, &
      mult, last_in)

    INTEGER, PARAMETER :: ndims = 1
    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name, units
    INTEGER, DIMENSION(:), INTENT(IN) :: dims
    INTEGER, INTENT(IN) :: nmat
    INTEGER(i4), INTENT(IN) :: stagger
    CHARACTER(LEN=*), INTENT(IN) :: mesh_id
    CHARACTER(LEN=*), INTENT(IN) :: material_names(:)
    REAL(r4), DIMENSION(:,:), INTENT(IN) :: variable
    INTEGER, INTENT(IN) :: distribution, subarray
    LOGICAL, OPTIONAL, INTENT(IN) :: convert
    REAL(r4), OPTIONAL, INTENT(IN) :: mult
    LOGICAL, OPTIONAL, INTENT(IN) :: last_in
    INTEGER :: i, idx
    INTEGER(i8) :: nsize, data_length
    INTEGER :: sz(ndims+1)
    LOGICAL :: last
    INTEGER, PARAMETER :: maxstring = 512
    CHARACTER(LEN=maxstring) :: temp_name
    CHARACTER(LEN=c_id_length), DIMENSION(:), ALLOCATABLE :: variable_ids

    IF (PRESENT(last_in)) THEN
      last = last_in
    ELSE
      last = .FALSE.
    END IF

    DO i = 1,ndims+1
      sz(i) = INT(SIZE(variable,i),i4)
    END DO

    ALLOCATE(variable_ids(nmat))

    nsize = 1
    DO i = 1,ndims
      nsize = nsize * dims(i)
    END DO
    nsize = sof * nsize

    data_length = 0
    DO i = 1,nmat
      IF (LEN_TRIM(material_names(i)) == 0) THEN
        variable_ids(i) = ''
      ELSE
        CALL sdf_safe_string_composite(h, id, &
            sdf_string_lowercase(material_names(i)), variable_ids(i))
        data_length = data_length + nsize
      END IF
    END DO

    h%datatype = datatype_real

    CALL sdf_write_stitched_material(h, id, name, mesh_id, stagger, &
        material_names, variable_ids, nmat, data_length)

    h%data_location = h%current_block%data_location
    idx = 0
    DO i = 1,nmat
      IF (LEN_TRIM(material_names(i)) == 0) CYCLE
      idx = idx + 1
      CALL sdf_safe_string_composite(h, name, material_names(i), temp_name)
      IF (last) THEN
        CALL write_1d_var_last_r4(h, variable_ids(i), temp_name, units, dims, &
            sz, stagger, mesh_id, variable, idx, distribution, subarray, &
            convert, mult)
      ELSE
        CALL write_1d_var_first_r4(h, variable_ids(i), temp_name, units, dims, &
            sz, stagger, mesh_id, variable, idx, distribution, subarray, &
            convert, mult)
      END IF
      h%data_location = h%data_location + nsize
    END DO

    h%data_location = 0

    DEALLOCATE(variable_ids)

  END SUBROUTINE write_1d_material_r4



  !----------------------------------------------------------------------------
  ! Code to write a 2D cartesian material mesh in parallel using the
  ! mpitype {distribution} for distribution of data
  ! It's up to the coder to design the distribution parallel operation, so
  ! need global dims
  !----------------------------------------------------------------------------

  SUBROUTINE write_2d_material_r4(h, id, name, units, dims, nmat, stagger, &
      mesh_id, material_names, variable, distribution, subarray, convert, &
      mult, last_in)

    INTEGER, PARAMETER :: ndims = 2
    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name, units
    INTEGER, DIMENSION(:), INTENT(IN) :: dims
    INTEGER, INTENT(IN) :: nmat
    INTEGER(i4), INTENT(IN) :: stagger
    CHARACTER(LEN=*), INTENT(IN) :: mesh_id
    CHARACTER(LEN=*), INTENT(IN) :: material_names(:)
    REAL(r4), DIMENSION(:,:,:), INTENT(IN) :: variable
    INTEGER, INTENT(IN) :: distribution, subarray
    LOGICAL, OPTIONAL, INTENT(IN) :: convert
    REAL(r4), OPTIONAL, INTENT(IN) :: mult
    LOGICAL, OPTIONAL, INTENT(IN) :: last_in
    INTEGER :: i, idx
    INTEGER(i8) :: nsize, data_length
    INTEGER :: sz(ndims+1)
    LOGICAL :: last
    INTEGER, PARAMETER :: maxstring = 512
    CHARACTER(LEN=maxstring) :: temp_name
    CHARACTER(LEN=c_id_length), DIMENSION(:), ALLOCATABLE :: variable_ids

    IF (PRESENT(last_in)) THEN
      last = last_in
    ELSE
      last = .FALSE.
    END IF

    DO i = 1,ndims+1
      sz(i) = INT(SIZE(variable,i),i4)
    END DO

    ALLOCATE(variable_ids(nmat))

    nsize = 1
    DO i = 1,ndims
      nsize = nsize * dims(i)
    END DO
    nsize = sof * nsize

    data_length = 0
    DO i = 1,nmat
      IF (LEN_TRIM(material_names(i)) == 0) THEN
        variable_ids(i) = ''
      ELSE
        CALL sdf_safe_string_composite(h, id, &
            sdf_string_lowercase(material_names(i)), variable_ids(i))
        data_length = data_length + nsize
      END IF
    END DO

    h%datatype = datatype_real

    CALL sdf_write_stitched_material(h, id, name, mesh_id, stagger, &
        material_names, variable_ids, nmat, data_length)

    h%data_location = h%current_block%data_location
    idx = 0
    DO i = 1,nmat
      IF (LEN_TRIM(material_names(i)) == 0) CYCLE
      idx = idx + 1
      CALL sdf_safe_string_composite(h, name, material_names(i), temp_name)
      IF (last) THEN
        CALL write_2d_var_last_r4(h, variable_ids(i), temp_name, units, dims, &
            sz, stagger, mesh_id, variable, idx, distribution, subarray, &
            convert, mult)
      ELSE
        CALL write_2d_var_first_r4(h, variable_ids(i), temp_name, units, dims, &
            sz, stagger, mesh_id, variable, idx, distribution, subarray, &
            convert, mult)
      END IF
      h%data_location = h%data_location + nsize
    END DO

    h%data_location = 0

    DEALLOCATE(variable_ids)

  END SUBROUTINE write_2d_material_r4



  !----------------------------------------------------------------------------
  ! Code to write a 3D cartesian material mesh in parallel using the
  ! mpitype {distribution} for distribution of data
  ! It's up to the coder to design the distribution parallel operation, so
  ! need global dims
  !----------------------------------------------------------------------------

  SUBROUTINE write_3d_material_r4(h, id, name, units, dims, nmat, stagger, &
      mesh_id, material_names, variable, distribution, subarray, convert, &
      mult, last_in)

    INTEGER, PARAMETER :: ndims = 3
    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name, units
    INTEGER, DIMENSION(:), INTENT(IN) :: dims
    INTEGER, INTENT(IN) :: nmat
    INTEGER(i4), INTENT(IN) :: stagger
    CHARACTER(LEN=*), INTENT(IN) :: mesh_id
    CHARACTER(LEN=*), INTENT(IN) :: material_names(:)
    REAL(r4), DIMENSION(:,:,:,:), INTENT(IN) :: variable
    INTEGER, INTENT(IN) :: distribution, subarray
    LOGICAL, OPTIONAL, INTENT(IN) :: convert
    REAL(r4), OPTIONAL, INTENT(IN) :: mult
    LOGICAL, OPTIONAL, INTENT(IN) :: last_in
    INTEGER :: i, idx
    INTEGER(i8) :: nsize, data_length
    INTEGER :: sz(ndims+1)
    LOGICAL :: last
    INTEGER, PARAMETER :: maxstring = 512
    CHARACTER(LEN=maxstring) :: temp_name
    CHARACTER(LEN=c_id_length), DIMENSION(:), ALLOCATABLE :: variable_ids

    IF (PRESENT(last_in)) THEN
      last = last_in
    ELSE
      last = .FALSE.
    END IF

    DO i = 1,ndims+1
      sz(i) = INT(SIZE(variable,i),i4)
    END DO

    ALLOCATE(variable_ids(nmat))

    nsize = 1
    DO i = 1,ndims
      nsize = nsize * dims(i)
    END DO
    nsize = sof * nsize

    data_length = 0
    DO i = 1,nmat
      IF (LEN_TRIM(material_names(i)) == 0) THEN
        variable_ids(i) = ''
      ELSE
        CALL sdf_safe_string_composite(h, id, &
            sdf_string_lowercase(material_names(i)), variable_ids(i))
        data_length = data_length + nsize
      END IF
    END DO

    h%datatype = datatype_real

    CALL sdf_write_stitched_material(h, id, name, mesh_id, stagger, &
        material_names, variable_ids, nmat, data_length)

    h%data_location = h%current_block%data_location
    idx = 0
    DO i = 1,nmat
      IF (LEN_TRIM(material_names(i)) == 0) CYCLE
      idx = idx + 1
      CALL sdf_safe_string_composite(h, name, material_names(i), temp_name)
      IF (last) THEN
        CALL write_3d_var_last_r4(h, variable_ids(i), temp_name, units, dims, &
            sz, stagger, mesh_id, variable, idx, distribution, subarray, &
            convert, mult)
      ELSE
        CALL write_3d_var_first_r4(h, variable_ids(i), temp_name, units, dims, &
            sz, stagger, mesh_id, variable, idx, distribution, subarray, &
            convert, mult)
      END IF
      h%data_location = h%data_location + nsize
    END DO

    h%data_location = 0

    DEALLOCATE(variable_ids)

  END SUBROUTINE write_3d_material_r4



  !----------------------------------------------------------------------------
  ! Code to write a 1D cartesian material variable in parallel using the
  ! mpitype {distribution} for distribution of data
  ! It's up to the coder to design the distribution parallel operation, so
  ! need global dims
  !----------------------------------------------------------------------------

  SUBROUTINE write_1d_matvar_r4(h, id, name, units, dims, nmat, stagger, &
      mesh_id, material_id, material_names, variable, distribution, &
      subarray, convert, mult, last_in)

    INTEGER, PARAMETER :: ndims = 1
    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name, units
    INTEGER, DIMENSION(:), INTENT(IN) :: dims
    INTEGER, INTENT(IN) :: nmat
    INTEGER(i4), INTENT(IN) :: stagger
    CHARACTER(LEN=*), INTENT(IN) :: mesh_id, material_id
    CHARACTER(LEN=*), INTENT(IN) :: material_names(:)
    REAL(r4), DIMENSION(:,:), INTENT(IN) :: variable
    INTEGER, INTENT(IN) :: distribution, subarray
    LOGICAL, OPTIONAL, INTENT(IN) :: convert
    REAL(r4), OPTIONAL, INTENT(IN) :: mult
    LOGICAL, OPTIONAL, INTENT(IN) :: last_in
    INTEGER :: i, idx
    INTEGER(i8) :: nsize, data_length
    INTEGER :: sz(ndims+1)
    LOGICAL :: last
    INTEGER, PARAMETER :: maxstring = 512
    CHARACTER(LEN=maxstring) :: temp_name
    CHARACTER(LEN=c_id_length), DIMENSION(:), ALLOCATABLE :: variable_ids

    IF (PRESENT(last_in)) THEN
      last = last_in
    ELSE
      last = .FALSE.
    END IF

    DO i = 1,ndims+1
      sz(i) = INT(SIZE(variable,i),i4)
    END DO

    ALLOCATE(variable_ids(nmat))

    nsize = 1
    DO i = 1,ndims
      nsize = nsize * dims(i)
    END DO
    nsize = sof * nsize

    data_length = 0
    DO i = 1,nmat
      IF (LEN_TRIM(material_names(i)) == 0) THEN
        variable_ids(i) = ''
      ELSE
        CALL sdf_safe_string_composite(h, id, &
            sdf_string_lowercase(material_names(i)), variable_ids(i))
        data_length = data_length + nsize
      END IF
    END DO

    h%datatype = datatype_real

    CALL sdf_write_stitched_matvar(h, id, name, mesh_id, stagger, &
        material_id, variable_ids, nmat, data_length)

    h%data_location = h%current_block%data_location
    idx = 0
    DO i = 1,nmat
      IF (LEN_TRIM(material_names(i)) == 0) CYCLE
      idx = idx + 1
      CALL sdf_safe_string_composite(h, name, material_names(i), temp_name)
      IF (last) THEN
        CALL write_1d_var_last_r4(h, variable_ids(i), temp_name, units, dims, &
            sz, stagger, mesh_id, variable, idx, distribution, subarray, &
            convert, mult)
      ELSE
        CALL write_1d_var_first_r4(h, variable_ids(i), temp_name, units, dims, &
            sz, stagger, mesh_id, variable, idx, distribution, subarray, &
            convert, mult)
      END IF
      h%data_location = h%data_location + nsize
    END DO

    h%data_location = 0

    DEALLOCATE(variable_ids)

  END SUBROUTINE write_1d_matvar_r4



  !----------------------------------------------------------------------------
  ! Code to write a 2D cartesian material variable in parallel using the
  ! mpitype {distribution} for distribution of data
  ! It's up to the coder to design the distribution parallel operation, so
  ! need global dims
  !----------------------------------------------------------------------------

  SUBROUTINE write_2d_matvar_r4(h, id, name, units, dims, nmat, stagger, &
      mesh_id, material_id, material_names, variable, distribution, &
      subarray, convert, mult, last_in)

    INTEGER, PARAMETER :: ndims = 2
    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name, units
    INTEGER, DIMENSION(:), INTENT(IN) :: dims
    INTEGER, INTENT(IN) :: nmat
    INTEGER(i4), INTENT(IN) :: stagger
    CHARACTER(LEN=*), INTENT(IN) :: mesh_id, material_id
    CHARACTER(LEN=*), INTENT(IN) :: material_names(:)
    REAL(r4), DIMENSION(:,:,:), INTENT(IN) :: variable
    INTEGER, INTENT(IN) :: distribution, subarray
    LOGICAL, OPTIONAL, INTENT(IN) :: convert
    REAL(r4), OPTIONAL, INTENT(IN) :: mult
    LOGICAL, OPTIONAL, INTENT(IN) :: last_in
    INTEGER :: i, idx
    INTEGER(i8) :: nsize, data_length
    INTEGER :: sz(ndims+1)
    LOGICAL :: last
    INTEGER, PARAMETER :: maxstring = 512
    CHARACTER(LEN=maxstring) :: temp_name
    CHARACTER(LEN=c_id_length), DIMENSION(:), ALLOCATABLE :: variable_ids

    IF (PRESENT(last_in)) THEN
      last = last_in
    ELSE
      last = .FALSE.
    END IF

    DO i = 1,ndims+1
      sz(i) = INT(SIZE(variable,i),i4)
    END DO

    ALLOCATE(variable_ids(nmat))

    nsize = 1
    DO i = 1,ndims
      nsize = nsize * dims(i)
    END DO
    nsize = sof * nsize

    data_length = 0
    DO i = 1,nmat
      IF (LEN_TRIM(material_names(i)) == 0) THEN
        variable_ids(i) = ''
      ELSE
        CALL sdf_safe_string_composite(h, id, &
            sdf_string_lowercase(material_names(i)), variable_ids(i))
        data_length = data_length + nsize
      END IF
    END DO

    h%datatype = datatype_real

    CALL sdf_write_stitched_matvar(h, id, name, mesh_id, stagger, &
        material_id, variable_ids, nmat, data_length)

    h%data_location = h%current_block%data_location
    idx = 0
    DO i = 1,nmat
      IF (LEN_TRIM(material_names(i)) == 0) CYCLE
      idx = idx + 1
      CALL sdf_safe_string_composite(h, name, material_names(i), temp_name)
      IF (last) THEN
        CALL write_2d_var_last_r4(h, variable_ids(i), temp_name, units, dims, &
            sz, stagger, mesh_id, variable, idx, distribution, subarray, &
            convert, mult)
      ELSE
        CALL write_2d_var_first_r4(h, variable_ids(i), temp_name, units, dims, &
            sz, stagger, mesh_id, variable, idx, distribution, subarray, &
            convert, mult)
      END IF
      h%data_location = h%data_location + nsize
    END DO

    h%data_location = 0

    DEALLOCATE(variable_ids)

  END SUBROUTINE write_2d_matvar_r4



  !----------------------------------------------------------------------------
  ! Code to write a 3D cartesian material variable in parallel using the
  ! mpitype {distribution} for distribution of data
  ! It's up to the coder to design the distribution parallel operation, so
  ! need global dims
  !----------------------------------------------------------------------------

  SUBROUTINE write_3d_matvar_r4(h, id, name, units, dims, nmat, stagger, &
      mesh_id, material_id, material_names, variable, distribution, &
      subarray, convert, mult, last_in)

    INTEGER, PARAMETER :: ndims = 3
    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name, units
    INTEGER, DIMENSION(:), INTENT(IN) :: dims
    INTEGER, INTENT(IN) :: nmat
    INTEGER(i4), INTENT(IN) :: stagger
    CHARACTER(LEN=*), INTENT(IN) :: mesh_id, material_id
    CHARACTER(LEN=*), INTENT(IN) :: material_names(:)
    REAL(r4), DIMENSION(:,:,:,:), INTENT(IN) :: variable
    INTEGER, INTENT(IN) :: distribution, subarray
    LOGICAL, OPTIONAL, INTENT(IN) :: convert
    REAL(r4), OPTIONAL, INTENT(IN) :: mult
    LOGICAL, OPTIONAL, INTENT(IN) :: last_in
    INTEGER :: i, idx
    INTEGER(i8) :: nsize, data_length
    INTEGER :: sz(ndims+1)
    LOGICAL :: last
    INTEGER, PARAMETER :: maxstring = 512
    CHARACTER(LEN=maxstring) :: temp_name
    CHARACTER(LEN=c_id_length), DIMENSION(:), ALLOCATABLE :: variable_ids

    IF (PRESENT(last_in)) THEN
      last = last_in
    ELSE
      last = .FALSE.
    END IF

    DO i = 1,ndims+1
      sz(i) = INT(SIZE(variable,i),i4)
    END DO

    ALLOCATE(variable_ids(nmat))

    nsize = 1
    DO i = 1,ndims
      nsize = nsize * dims(i)
    END DO
    nsize = sof * nsize

    data_length = 0
    DO i = 1,nmat
      IF (LEN_TRIM(material_names(i)) == 0) THEN
        variable_ids(i) = ''
      ELSE
        CALL sdf_safe_string_composite(h, id, &
            sdf_string_lowercase(material_names(i)), variable_ids(i))
        data_length = data_length + nsize
      END IF
    END DO

    h%datatype = datatype_real

    CALL sdf_write_stitched_matvar(h, id, name, mesh_id, stagger, &
        material_id, variable_ids, nmat, data_length)

    h%data_location = h%current_block%data_location
    idx = 0
    DO i = 1,nmat
      IF (LEN_TRIM(material_names(i)) == 0) CYCLE
      idx = idx + 1
      CALL sdf_safe_string_composite(h, name, material_names(i), temp_name)
      IF (last) THEN
        CALL write_3d_var_last_r4(h, variable_ids(i), temp_name, units, dims, &
            sz, stagger, mesh_id, variable, idx, distribution, subarray, &
            convert, mult)
      ELSE
        CALL write_3d_var_first_r4(h, variable_ids(i), temp_name, units, dims, &
            sz, stagger, mesh_id, variable, idx, distribution, subarray, &
            convert, mult)
      END IF
      h%data_location = h%data_location + nsize
    END DO

    h%data_location = 0

    DEALLOCATE(variable_ids)

  END SUBROUTINE write_3d_matvar_r4



  !----------------------------------------------------------------------------
  ! Code to write a 1D cartesian material species in parallel using the
  ! mpitype {distribution} for distribution of data
  ! It's up to the coder to design the distribution parallel operation, so
  ! need global dims
  !----------------------------------------------------------------------------

  SUBROUTINE write_1d_species_r4(h, id, name, units, dims, nmat, stagger, &
      mesh_id, material_id, material_name, specnames, variable, distribution, &
      subarray, convert, mult, last_in)

    INTEGER, PARAMETER :: ndims = 1
    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name, units
    INTEGER, DIMENSION(:), INTENT(IN) :: dims
    INTEGER, INTENT(IN) :: nmat
    INTEGER(i4), INTENT(IN) :: stagger
    CHARACTER(LEN=*), INTENT(IN) :: mesh_id, material_id
    CHARACTER(LEN=*), INTENT(IN) :: material_name
    CHARACTER(LEN=*), INTENT(IN) :: specnames(:)
    REAL(r4), DIMENSION(:,:), INTENT(IN) :: variable
    INTEGER, INTENT(IN) :: distribution, subarray
    LOGICAL, OPTIONAL, INTENT(IN) :: convert
    REAL(r4), OPTIONAL, INTENT(IN) :: mult
    LOGICAL, OPTIONAL, INTENT(IN) :: last_in
    INTEGER :: i, idx
    INTEGER(i8) :: nsize, data_length
    INTEGER :: sz(ndims+1)
    LOGICAL :: last
    INTEGER, PARAMETER :: maxstring = 512
    CHARACTER(LEN=maxstring) :: temp_name
    CHARACTER(LEN=c_id_length), DIMENSION(:), ALLOCATABLE :: variable_ids

    IF (PRESENT(last_in)) THEN
      last = last_in
    ELSE
      last = .FALSE.
    END IF

    DO i = 1,ndims+1
      sz(i) = INT(SIZE(variable,i),i4)
    END DO

    ALLOCATE(variable_ids(nmat))

    nsize = 1
    DO i = 1,ndims
      nsize = nsize * dims(i)
    END DO
    nsize = sof * nsize

    data_length = 0
    DO i = 1,nmat
      IF (LEN_TRIM(specnames(i)) == 0) THEN
        variable_ids(i) = ''
      ELSE
        CALL sdf_safe_string_composite(h, id, &
            sdf_string_lowercase(specnames(i)), variable_ids(i))
        data_length = data_length + nsize
      END IF
    END DO

    h%datatype = datatype_real

    CALL sdf_write_stitched_species(h, id, name, mesh_id, stagger, &
        material_id, material_name, specnames, variable_ids, nmat, data_length)

    h%data_location = h%current_block%data_location
    idx = 0
    DO i = 1,nmat
      IF (LEN_TRIM(specnames(i)) == 0) CYCLE
      idx = idx + 1
      CALL sdf_safe_string_composite(h, name, specnames(i), temp_name)
      IF (last) THEN
        CALL write_1d_var_last_r4(h, variable_ids(i), temp_name, units, dims, &
            sz, stagger, mesh_id, variable, idx, distribution, subarray, &
            convert, mult)
      ELSE
        CALL write_1d_var_first_r4(h, variable_ids(i), temp_name, units, dims, &
            sz, stagger, mesh_id, variable, idx, distribution, subarray, &
            convert, mult)
      END IF
      h%data_location = h%data_location + nsize
    END DO

    h%data_location = 0

    DEALLOCATE(variable_ids)

  END SUBROUTINE write_1d_species_r4



  !----------------------------------------------------------------------------
  ! Code to write a 2D cartesian material species in parallel using the
  ! mpitype {distribution} for distribution of data
  ! It's up to the coder to design the distribution parallel operation, so
  ! need global dims
  !----------------------------------------------------------------------------

  SUBROUTINE write_2d_species_r4(h, id, name, units, dims, nmat, stagger, &
      mesh_id, material_id, material_name, specnames, variable, distribution, &
      subarray, convert, mult, last_in)

    INTEGER, PARAMETER :: ndims = 2
    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name, units
    INTEGER, DIMENSION(:), INTENT(IN) :: dims
    INTEGER, INTENT(IN) :: nmat
    INTEGER(i4), INTENT(IN) :: stagger
    CHARACTER(LEN=*), INTENT(IN) :: mesh_id, material_id
    CHARACTER(LEN=*), INTENT(IN) :: material_name
    CHARACTER(LEN=*), INTENT(IN) :: specnames(:)
    REAL(r4), DIMENSION(:,:,:), INTENT(IN) :: variable
    INTEGER, INTENT(IN) :: distribution, subarray
    LOGICAL, OPTIONAL, INTENT(IN) :: convert
    REAL(r4), OPTIONAL, INTENT(IN) :: mult
    LOGICAL, OPTIONAL, INTENT(IN) :: last_in
    INTEGER :: i, idx
    INTEGER(i8) :: nsize, data_length
    INTEGER :: sz(ndims+1)
    LOGICAL :: last
    INTEGER, PARAMETER :: maxstring = 512
    CHARACTER(LEN=maxstring) :: temp_name
    CHARACTER(LEN=c_id_length), DIMENSION(:), ALLOCATABLE :: variable_ids

    IF (PRESENT(last_in)) THEN
      last = last_in
    ELSE
      last = .FALSE.
    END IF

    DO i = 1,ndims+1
      sz(i) = INT(SIZE(variable,i),i4)
    END DO

    ALLOCATE(variable_ids(nmat))

    nsize = 1
    DO i = 1,ndims
      nsize = nsize * dims(i)
    END DO
    nsize = sof * nsize

    data_length = 0
    DO i = 1,nmat
      IF (LEN_TRIM(specnames(i)) == 0) THEN
        variable_ids(i) = ''
      ELSE
        CALL sdf_safe_string_composite(h, id, &
            sdf_string_lowercase(specnames(i)), variable_ids(i))
        data_length = data_length + nsize
      END IF
    END DO

    h%datatype = datatype_real

    CALL sdf_write_stitched_species(h, id, name, mesh_id, stagger, &
        material_id, material_name, specnames, variable_ids, nmat, data_length)

    h%data_location = h%current_block%data_location
    idx = 0
    DO i = 1,nmat
      IF (LEN_TRIM(specnames(i)) == 0) CYCLE
      idx = idx + 1
      CALL sdf_safe_string_composite(h, name, specnames(i), temp_name)
      IF (last) THEN
        CALL write_2d_var_last_r4(h, variable_ids(i), temp_name, units, dims, &
            sz, stagger, mesh_id, variable, idx, distribution, subarray, &
            convert, mult)
      ELSE
        CALL write_2d_var_first_r4(h, variable_ids(i), temp_name, units, dims, &
            sz, stagger, mesh_id, variable, idx, distribution, subarray, &
            convert, mult)
      END IF
      h%data_location = h%data_location + nsize
    END DO

    h%data_location = 0

    DEALLOCATE(variable_ids)

  END SUBROUTINE write_2d_species_r4



  !----------------------------------------------------------------------------
  ! Code to write a 3D cartesian material species in parallel using the
  ! mpitype {distribution} for distribution of data
  ! It's up to the coder to design the distribution parallel operation, so
  ! need global dims
  !----------------------------------------------------------------------------

  SUBROUTINE write_3d_species_r4(h, id, name, units, dims, nmat, stagger, &
      mesh_id, material_id, material_name, specnames, variable, distribution, &
      subarray, convert, mult, last_in)

    INTEGER, PARAMETER :: ndims = 3
    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name, units
    INTEGER, DIMENSION(:), INTENT(IN) :: dims
    INTEGER, INTENT(IN) :: nmat
    INTEGER(i4), INTENT(IN) :: stagger
    CHARACTER(LEN=*), INTENT(IN) :: mesh_id, material_id
    CHARACTER(LEN=*), INTENT(IN) :: material_name
    CHARACTER(LEN=*), INTENT(IN) :: specnames(:)
    REAL(r4), DIMENSION(:,:,:,:), INTENT(IN) :: variable
    INTEGER, INTENT(IN) :: distribution, subarray
    LOGICAL, OPTIONAL, INTENT(IN) :: convert
    REAL(r4), OPTIONAL, INTENT(IN) :: mult
    LOGICAL, OPTIONAL, INTENT(IN) :: last_in
    INTEGER :: i, idx
    INTEGER(i8) :: nsize, data_length
    INTEGER :: sz(ndims+1)
    LOGICAL :: last
    INTEGER, PARAMETER :: maxstring = 512
    CHARACTER(LEN=maxstring) :: temp_name
    CHARACTER(LEN=c_id_length), DIMENSION(:), ALLOCATABLE :: variable_ids

    IF (PRESENT(last_in)) THEN
      last = last_in
    ELSE
      last = .FALSE.
    END IF

    DO i = 1,ndims+1
      sz(i) = INT(SIZE(variable,i),i4)
    END DO

    ALLOCATE(variable_ids(nmat))

    nsize = 1
    DO i = 1,ndims
      nsize = nsize * dims(i)
    END DO
    nsize = sof * nsize

    data_length = 0
    DO i = 1,nmat
      IF (LEN_TRIM(specnames(i)) == 0) THEN
        variable_ids(i) = ''
      ELSE
        CALL sdf_safe_string_composite(h, id, &
            sdf_string_lowercase(specnames(i)), variable_ids(i))
        data_length = data_length + nsize
      END IF
    END DO

    h%datatype = datatype_real

    CALL sdf_write_stitched_species(h, id, name, mesh_id, stagger, &
        material_id, material_name, specnames, variable_ids, nmat, data_length)

    h%data_location = h%current_block%data_location
    idx = 0
    DO i = 1,nmat
      IF (LEN_TRIM(specnames(i)) == 0) CYCLE
      idx = idx + 1
      CALL sdf_safe_string_composite(h, name, specnames(i), temp_name)
      IF (last) THEN
        CALL write_3d_var_last_r4(h, variable_ids(i), temp_name, units, dims, &
            sz, stagger, mesh_id, variable, idx, distribution, subarray, &
            convert, mult)
      ELSE
        CALL write_3d_var_first_r4(h, variable_ids(i), temp_name, units, dims, &
            sz, stagger, mesh_id, variable, idx, distribution, subarray, &
            convert, mult)
      END IF
      h%data_location = h%data_location + nsize
    END DO

    h%data_location = 0

    DEALLOCATE(variable_ids)

  END SUBROUTINE write_3d_species_r4



  !----------------------------------------------------------------------------
  ! Code to write a 1D cartesian stitched variable in parallel using the
  ! mpitype {distribution} for distribution of data
  ! It's up to the coder to design the distribution parallel operation, so
  ! need global dims
  !----------------------------------------------------------------------------

  SUBROUTINE write_1d_stitched_var_r4(h, id, name, units, dims, nmat, stagger, &
      mesh_id, material_names, variable, distribution, subarray, convert, &
      mult, last_in)

    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name, units
    INTEGER, DIMENSION(:), INTENT(IN) :: dims
    INTEGER, INTENT(IN) :: nmat
    INTEGER(i4), INTENT(IN) :: stagger
    CHARACTER(LEN=*), INTENT(IN) :: mesh_id
    CHARACTER(LEN=*), INTENT(IN) :: material_names(:)
    REAL(r4), DIMENSION(:,:), INTENT(IN) :: variable
    INTEGER, INTENT(IN) :: distribution, subarray
    LOGICAL, OPTIONAL, INTENT(IN) :: convert
    REAL(r4), OPTIONAL, INTENT(IN) :: mult
    LOGICAL, OPTIONAL, INTENT(IN) :: last_in

    h%blocktype = c_blocktype_contiguous
    CALL write_1d_material_r4(h, id, name, units, dims, nmat, stagger, &
        mesh_id, material_names, variable, distribution, subarray, convert, &
        mult, last_in)
    h%blocktype = 0

  END SUBROUTINE write_1d_stitched_var_r4



  !----------------------------------------------------------------------------
  ! Code to write a 2D cartesian stitched variable in parallel using the
  ! mpitype {distribution} for distribution of data
  ! It's up to the coder to design the distribution parallel operation, so
  ! need global dims
  !----------------------------------------------------------------------------

  SUBROUTINE write_2d_stitched_var_r4(h, id, name, units, dims, nmat, stagger, &
      mesh_id, material_names, variable, distribution, subarray, convert, &
      mult, last_in)

    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name, units
    INTEGER, DIMENSION(:), INTENT(IN) :: dims
    INTEGER, INTENT(IN) :: nmat
    INTEGER(i4), INTENT(IN) :: stagger
    CHARACTER(LEN=*), INTENT(IN) :: mesh_id
    CHARACTER(LEN=*), INTENT(IN) :: material_names(:)
    REAL(r4), DIMENSION(:,:,:), INTENT(IN) :: variable
    INTEGER, INTENT(IN) :: distribution, subarray
    LOGICAL, OPTIONAL, INTENT(IN) :: convert
    REAL(r4), OPTIONAL, INTENT(IN) :: mult
    LOGICAL, OPTIONAL, INTENT(IN) :: last_in

    h%blocktype = c_blocktype_contiguous
    CALL write_2d_material_r4(h, id, name, units, dims, nmat, stagger, &
        mesh_id, material_names, variable, distribution, subarray, convert, &
        mult, last_in)
    h%blocktype = 0

  END SUBROUTINE write_2d_stitched_var_r4



  !----------------------------------------------------------------------------
  ! Code to write a 3D cartesian stitched variable in parallel using the
  ! mpitype {distribution} for distribution of data
  ! It's up to the coder to design the distribution parallel operation, so
  ! need global dims
  !----------------------------------------------------------------------------

  SUBROUTINE write_3d_stitched_var_r4(h, id, name, units, dims, nmat, stagger, &
      mesh_id, material_names, variable, distribution, subarray, convert, &
      mult, last_in)

    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name, units
    INTEGER, DIMENSION(:), INTENT(IN) :: dims
    INTEGER, INTENT(IN) :: nmat
    INTEGER(i4), INTENT(IN) :: stagger
    CHARACTER(LEN=*), INTENT(IN) :: mesh_id
    CHARACTER(LEN=*), INTENT(IN) :: material_names(:)
    REAL(r4), DIMENSION(:,:,:,:), INTENT(IN) :: variable
    INTEGER, INTENT(IN) :: distribution, subarray
    LOGICAL, OPTIONAL, INTENT(IN) :: convert
    REAL(r4), OPTIONAL, INTENT(IN) :: mult
    LOGICAL, OPTIONAL, INTENT(IN) :: last_in

    h%blocktype = c_blocktype_contiguous
    CALL write_3d_material_r4(h, id, name, units, dims, nmat, stagger, &
        mesh_id, material_names, variable, distribution, subarray, convert, &
        mult, last_in)
    h%blocktype = 0

  END SUBROUTINE write_3d_stitched_var_r4

END MODULE sdf_output_cartesian_r4

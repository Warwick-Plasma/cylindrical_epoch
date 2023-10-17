!
! SDF (Self-Describing Format) Fortran Library
! Copyright (c) 2011-2016, SDF Development Team
!
! Distributed under the terms of the BSD 3-clause License.
! See the LICENSE file for details.
!

MODULE sdf_output_cartesian_ru

  USE sdf_output_ru

  IMPLICIT NONE

CONTAINS

  SUBROUTINE write_mesh_meta_r8(h, id, name, dim_labels, dim_units, dim_mults)

    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: id, name
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: dim_labels(:), dim_units(:)
    REAL(r8), DIMENSION(:), INTENT(IN), OPTIONAL :: dim_mults
    INTEGER :: ndims
    TYPE(sdf_block_type), POINTER :: b
    INTEGER :: i, errcode

    b => h%current_block

    ndims = b%ndims
    IF (b%blocktype == c_blocktype_lagrangian_mesh) THEN
      b%nelements = ndims
      DO i = 1,ndims
        b%nelements = b%nelements * b%dims(i)
      END DO
    ELSE
      b%blocktype = c_blocktype_plain_mesh
      b%nelements = 0
      DO i = 1,ndims
        b%nelements = b%nelements + b%dims(i)
      END DO
    END IF

    ! Metadata is
    ! - mults     REAL(r8), DIMENSION(ndims)
    ! - labels    CHARACTER(id_length), DIMENSION(ndims)
    ! - units     CHARACTER(id_length), DIMENSION(ndims)
    ! - geometry  INTEGER(i4)
    ! - minval    REAL(r8), DIMENSION(ndims)
    ! - maxval    REAL(r8), DIMENSION(ndims)
    ! - dims      INTEGER(i4), DIMENSION(ndims)

    b%info_length = h%block_header_length + (ndims + 1) * soi4 &
        + (3 * ndims) * sof8 + 2 * ndims * c_id_length
    b%data_length = b%nelements * b%type_size

    ! Write header

    IF (PRESENT(id)) THEN
      ALLOCATE(b%dim_labels(ndims), b%dim_units(ndims), b%dim_mults(ndims))

      IF (PRESENT(dim_labels)) THEN
        DO i = 1,ndims
          CALL sdf_safe_copy_id(h, dim_labels(i), b%dim_labels(i))
        END DO
      ELSE
        IF (ndims >= 1) CALL sdf_safe_copy_id(h, 'X', b%dim_labels(1))
        IF (ndims >= 2) CALL sdf_safe_copy_id(h, 'Y', b%dim_labels(2))
        IF (ndims >= 3) CALL sdf_safe_copy_id(h, 'Z', b%dim_labels(3))
      END IF

      IF (PRESENT(dim_units)) THEN
        DO i = 1,ndims
          CALL sdf_safe_copy_id(h, dim_units(i), b%dim_units(i))
        END DO
      ELSE
        DO i = 1,ndims
          CALL sdf_safe_copy_id(h, 'm', b%dim_units(i))
        END DO
      END IF

      IF (PRESENT(dim_mults)) THEN
        DO i = 1,ndims
          b%dim_mults(i) = REAL(dim_mults(i),r8)
        END DO
      ELSE
        DO i = 1,ndims
          b%dim_mults(i) = 1.d0
        END DO
      END IF

      CALL sdf_write_block_header(h, id, name)
    ELSE
      CALL write_block_header(h)
    END IF

    IF (h%rank == h%rank_master) THEN
      CALL MPI_FILE_WRITE(h%filehandle, b%dim_mults, ndims, MPI_REAL8, &
          MPI_STATUS_IGNORE, errcode)

      DO i = 1,ndims
        CALL sdf_safe_write_id(h, b%dim_labels(i))
      END DO

      DO i = 1,ndims
        CALL sdf_safe_write_id(h, b%dim_units(i))
      END DO

      CALL MPI_FILE_WRITE(h%filehandle, b%geometry, 1, MPI_INTEGER4, &
          MPI_STATUS_IGNORE, errcode)

      CALL MPI_FILE_WRITE(h%filehandle, b%extents, 2 * ndims, MPI_REAL8, &
          MPI_STATUS_IGNORE, errcode)

      CALL MPI_FILE_WRITE(h%filehandle, b%dims, ndims, MPI_INTEGER4, &
          MPI_STATUS_IGNORE, errcode)
    END IF

    h%current_location = b%block_start + b%info_length
    b%done_info = .TRUE.

  END SUBROUTINE write_mesh_meta_r8



  SUBROUTINE write_mesh_meta_r4(h, id, name, dim_labels, dim_units, dim_mults)

    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: id, name
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: dim_labels(:), dim_units(:)
    REAL(r4), DIMENSION(:), INTENT(IN), OPTIONAL :: dim_mults
    REAL(r8), DIMENSION(c_maxdims) :: dim_mults8
    TYPE(sdf_block_type), POINTER :: b
    INTEGER :: i

    IF (PRESENT(dim_mults)) THEN
      b => h%current_block
      DO i = 1,b%ndims
        dim_mults8(i) = REAL(dim_mults(i),r8)
      END DO

      CALL write_mesh_meta_r8(h, id, name, dim_labels, dim_units, dim_mults8)
    ELSE
      CALL write_mesh_meta_r8(h, id, name, dim_labels, dim_units)
    END IF

  END SUBROUTINE write_mesh_meta_r4



  SUBROUTINE write_mesh_variable_meta_r8(h, id, name, units, mesh_id, mult)

    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: id, name, units, mesh_id
    REAL(r8), INTENT(IN), OPTIONAL :: mult
    INTEGER :: ndims
    TYPE(sdf_block_type), POINTER :: b
    INTEGER :: i, errcode

    b => h%current_block

    b%blocktype = c_blocktype_plain_variable
    ndims = b%ndims

    b%nelements = 1
    DO i = 1,ndims
      b%nelements = b%nelements * b%dims(i)
    END DO

    ! Metadata is
    ! - mult      REAL(r8)
    ! - units     CHARACTER(id_length)
    ! - meshid    CHARACTER(id_length)
    ! - dims      INTEGER(i4), DIMENSION(ndims)
    ! - stagger   INTEGER(i4)

    b%info_length = h%block_header_length + (ndims + 1) * soi4 + sof8 &
        + 2 * c_id_length
    b%data_length = b%nelements * b%type_size

    ! Write header

    IF (PRESENT(id)) THEN
      CALL sdf_safe_copy_id(h, units, b%units)
      CALL sdf_safe_copy_id(h, mesh_id, b%mesh_id)

      IF (PRESENT(mult)) THEN
        b%mult = REAL(mult,r8)
      ELSE
        b%mult = 1.d0
      END IF

      CALL sdf_write_block_header(h, id, name)
    ELSE
      CALL write_block_header(h)
    END IF

    IF (h%rank == h%rank_master) THEN
      CALL MPI_FILE_WRITE(h%filehandle, b%mult, 1, MPI_REAL8, &
          MPI_STATUS_IGNORE, errcode)

      CALL sdf_safe_write_id(h, b%units)

      CALL sdf_safe_write_id(h, b%mesh_id)

      CALL MPI_FILE_WRITE(h%filehandle, b%dims, ndims, MPI_INTEGER4, &
          MPI_STATUS_IGNORE, errcode)

      CALL MPI_FILE_WRITE(h%filehandle, b%stagger, 1, MPI_INTEGER4, &
          MPI_STATUS_IGNORE, errcode)
    END IF

    h%current_location = b%block_start + b%info_length
    b%done_info = .TRUE.

  END SUBROUTINE write_mesh_variable_meta_r8



  SUBROUTINE write_mesh_variable_meta_r4(h, id, name, units, mesh_id, mult)

    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: id, name, units, mesh_id
    REAL(r4), INTENT(IN), OPTIONAL :: mult

    IF (PRESENT(mult)) THEN
      CALL write_mesh_variable_meta_r8(h, id, name, units, mesh_id, &
          REAL(mult,r8))
    ELSE
      CALL write_mesh_variable_meta_r8(h, id, name, units, mesh_id)
    END IF

  END SUBROUTINE write_mesh_variable_meta_r4



  !----------------------------------------------------------------------------
  ! Code to write a 1D cartesian integer variable in serial from the node
  ! with rank {rank_write}
  !----------------------------------------------------------------------------

  SUBROUTINE write_srl_1d_integer_i4_r8(h, id, name, units, dims, stagger, &
      mesh_id, variable, subarray, mult, rank_write)

    INTEGER, PARAMETER :: ndims = 1
    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name, units
    INTEGER, DIMENSION(:), INTENT(IN) :: dims
    INTEGER(i4), INTENT(IN) :: stagger
    CHARACTER(LEN=*), INTENT(IN) :: mesh_id
    INTEGER(i4), DIMENSION(:), INTENT(IN) :: variable
    INTEGER, OPTIONAL, INTENT(IN) :: subarray
    REAL(r8), OPTIONAL, INTENT(IN) :: mult
    INTEGER, OPTIONAL, INTENT(IN) :: rank_write
    INTEGER :: i, errcode, mpitype, nitems, nelements
    TYPE(sdf_block_type), POINTER :: b

    CALL sdf_get_next_block(h)
    b => h%current_block

    b%type_size = 4
    b%datatype = c_datatype_integer4
    b%mpitype = MPI_INTEGER4
    b%ndims = ndims
    b%stagger = stagger

    IF (PRESENT(mult)) THEN
      b%mult = REAL(mult,r8)
    ELSE
      b%mult = 1.d0
    END IF

    nelements = 1
    DO i = 1,ndims
      b%dims(i) = INT(dims(i),i4)
      nelements = nelements * b%dims(i)
    END DO

    IF (PRESENT(rank_write)) h%rank_master = rank_write

    ! Write header

    CALL write_mesh_variable_meta_r8(h, id, name, units, mesh_id, mult)

    ! Write the actual data

    IF (h%rank == h%rank_master) THEN
      h%current_location = b%data_location

      CALL MPI_FILE_SEEK(h%filehandle, h%current_location, MPI_SEEK_SET, &
          errcode)

      IF (PRESENT(subarray)) THEN
        mpitype = subarray
        nitems = 1
      ELSE
        mpitype = b%mpitype
        nitems = nelements
      END IF

      CALL MPI_FILE_WRITE(h%filehandle, variable, nitems, mpitype, &
          MPI_STATUS_IGNORE, errcode)
    END IF

    h%current_location = b%data_location + b%data_length
    b%done_data = .TRUE.

  END SUBROUTINE write_srl_1d_integer_i4_r8



  !----------------------------------------------------------------------------
  ! Code to write a 2D cartesian integer variable in serial from the node
  ! with rank {rank_write}
  !----------------------------------------------------------------------------

  SUBROUTINE write_srl_2d_integer_i4_r8(h, id, name, units, dims, stagger, &
      mesh_id, variable, subarray, mult, rank_write)

    INTEGER, PARAMETER :: ndims = 2
    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name, units
    INTEGER, DIMENSION(:), INTENT(IN) :: dims
    INTEGER(i4), INTENT(IN) :: stagger
    CHARACTER(LEN=*), INTENT(IN) :: mesh_id
    INTEGER(i4), DIMENSION(:,:), INTENT(IN) :: variable
    INTEGER, OPTIONAL, INTENT(IN) :: subarray
    REAL(r8), OPTIONAL, INTENT(IN) :: mult
    INTEGER, OPTIONAL, INTENT(IN) :: rank_write
    INTEGER :: i, errcode, mpitype, nitems, nelements
    TYPE(sdf_block_type), POINTER :: b

    CALL sdf_get_next_block(h)
    b => h%current_block

    b%type_size = 4
    b%datatype = c_datatype_integer4
    b%mpitype = MPI_INTEGER4
    b%ndims = ndims
    b%stagger = stagger

    IF (PRESENT(mult)) THEN
      b%mult = REAL(mult,r8)
    ELSE
      b%mult = 1.d0
    END IF

    nelements = 1
    DO i = 1,ndims
      b%dims(i) = INT(dims(i),i4)
      nelements = nelements * b%dims(i)
    END DO

    IF (PRESENT(rank_write)) h%rank_master = rank_write

    ! Write header

    CALL write_mesh_variable_meta_r8(h, id, name, units, mesh_id, mult)

    ! Write the actual data

    IF (h%rank == h%rank_master) THEN
      h%current_location = b%data_location

      CALL MPI_FILE_SEEK(h%filehandle, h%current_location, MPI_SEEK_SET, &
          errcode)

      IF (PRESENT(subarray)) THEN
        mpitype = subarray
        nitems = 1
      ELSE
        mpitype = b%mpitype
        nitems = nelements
      END IF

      CALL MPI_FILE_WRITE(h%filehandle, variable, nitems, mpitype, &
          MPI_STATUS_IGNORE, errcode)
    END IF

    h%current_location = b%data_location + b%data_length
    b%done_data = .TRUE.

  END SUBROUTINE write_srl_2d_integer_i4_r8



  !----------------------------------------------------------------------------
  ! Code to write a 3D cartesian integer variable in serial from the node
  ! with rank {rank_write}
  !----------------------------------------------------------------------------

  SUBROUTINE write_srl_3d_integer_i4_r8(h, id, name, units, dims, stagger, &
      mesh_id, variable, subarray, mult, rank_write)

    INTEGER, PARAMETER :: ndims = 3
    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name, units
    INTEGER, DIMENSION(:), INTENT(IN) :: dims
    INTEGER(i4), INTENT(IN) :: stagger
    CHARACTER(LEN=*), INTENT(IN) :: mesh_id
    INTEGER(i4), DIMENSION(:,:,:), INTENT(IN) :: variable
    INTEGER, OPTIONAL, INTENT(IN) :: subarray
    REAL(r8), OPTIONAL, INTENT(IN) :: mult
    INTEGER, OPTIONAL, INTENT(IN) :: rank_write
    INTEGER :: i, errcode, mpitype, nitems, nelements
    TYPE(sdf_block_type), POINTER :: b

    CALL sdf_get_next_block(h)
    b => h%current_block

    b%type_size = 4
    b%datatype = c_datatype_integer4
    b%mpitype = MPI_INTEGER4
    b%ndims = ndims
    b%stagger = stagger

    IF (PRESENT(mult)) THEN
      b%mult = REAL(mult,r8)
    ELSE
      b%mult = 1.d0
    END IF

    nelements = 1
    DO i = 1,ndims
      b%dims(i) = INT(dims(i),i4)
      nelements = nelements * b%dims(i)
    END DO

    IF (PRESENT(rank_write)) h%rank_master = rank_write

    ! Write header

    CALL write_mesh_variable_meta_r8(h, id, name, units, mesh_id, mult)

    ! Write the actual data

    IF (h%rank == h%rank_master) THEN
      h%current_location = b%data_location

      CALL MPI_FILE_SEEK(h%filehandle, h%current_location, MPI_SEEK_SET, &
          errcode)

      IF (PRESENT(subarray)) THEN
        mpitype = subarray
        nitems = 1
      ELSE
        mpitype = b%mpitype
        nitems = nelements
      END IF

      CALL MPI_FILE_WRITE(h%filehandle, variable, nitems, mpitype, &
          MPI_STATUS_IGNORE, errcode)
    END IF

    h%current_location = b%data_location + b%data_length
    b%done_data = .TRUE.

  END SUBROUTINE write_srl_3d_integer_i4_r8



  SUBROUTINE write_srl_1d_integer_i4_r4(h, id, name, units, dims, stagger, &
      mesh_id, variable, subarray, mult, rank_write)

    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name, units
    INTEGER, DIMENSION(:), INTENT(IN) :: dims
    INTEGER(i4), INTENT(IN) :: stagger
    CHARACTER(LEN=*), INTENT(IN) :: mesh_id
    INTEGER(i4), DIMENSION(:), INTENT(IN) :: variable
    INTEGER, INTENT(IN) :: subarray
    REAL(r4), INTENT(IN) :: mult
    INTEGER, OPTIONAL, INTENT(IN) :: rank_write

    CALL write_srl_1d_integer_i4_r8(h, id, name, units, dims, stagger, &
        mesh_id, variable, subarray, REAL(mult,r8), rank_write)

  END SUBROUTINE write_srl_1d_integer_i4_r4



  SUBROUTINE write_srl_2d_integer_i4_r4(h, id, name, units, dims, stagger, &
      mesh_id, variable, subarray, mult, rank_write)

    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name, units
    INTEGER, DIMENSION(:), INTENT(IN) :: dims
    INTEGER(i4), INTENT(IN) :: stagger
    CHARACTER(LEN=*), INTENT(IN) :: mesh_id
    INTEGER(i4), DIMENSION(:,:), INTENT(IN) :: variable
    INTEGER, INTENT(IN) :: subarray
    REAL(r4), INTENT(IN) :: mult
    INTEGER, OPTIONAL, INTENT(IN) :: rank_write

    CALL write_srl_2d_integer_i4_r8(h, id, name, units, dims, stagger, &
        mesh_id, variable, subarray, REAL(mult,r8), rank_write)

  END SUBROUTINE write_srl_2d_integer_i4_r4



  SUBROUTINE write_srl_3d_integer_i4_r4(h, id, name, units, dims, stagger, &
      mesh_id, variable, subarray, mult, rank_write)

    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name, units
    INTEGER, DIMENSION(:), INTENT(IN) :: dims
    INTEGER(i4), INTENT(IN) :: stagger
    CHARACTER(LEN=*), INTENT(IN) :: mesh_id
    INTEGER(i4), DIMENSION(:,:,:), INTENT(IN) :: variable
    INTEGER, INTENT(IN) :: subarray
    REAL(r4), INTENT(IN) :: mult
    INTEGER, OPTIONAL, INTENT(IN) :: rank_write

    CALL write_srl_3d_integer_i4_r8(h, id, name, units, dims, stagger, &
        mesh_id, variable, subarray, REAL(mult,r8), rank_write)

  END SUBROUTINE write_srl_3d_integer_i4_r4



  !----------------------------------------------------------------------------
  ! Code to write a 1D cartesian integer variable in serial from the node
  ! with rank {rank_write}
  !----------------------------------------------------------------------------

  SUBROUTINE write_srl_1d_integer_i8_r8(h, id, name, units, dims, stagger, &
      mesh_id, variable, subarray, mult, rank_write)

    INTEGER, PARAMETER :: ndims = 1
    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name, units
    INTEGER, DIMENSION(:), INTENT(IN) :: dims
    INTEGER(i4), INTENT(IN) :: stagger
    CHARACTER(LEN=*), INTENT(IN) :: mesh_id
    INTEGER(i8), DIMENSION(:), INTENT(IN) :: variable
    INTEGER, OPTIONAL, INTENT(IN) :: subarray
    REAL(r8), OPTIONAL, INTENT(IN) :: mult
    INTEGER, OPTIONAL, INTENT(IN) :: rank_write
    INTEGER :: i, errcode, mpitype, nitems, nelements
    TYPE(sdf_block_type), POINTER :: b

    CALL sdf_get_next_block(h)
    b => h%current_block

    b%type_size = 8
    b%datatype = c_datatype_integer8
    b%mpitype = MPI_INTEGER8
    b%ndims = ndims
    b%stagger = stagger

    IF (PRESENT(mult)) THEN
      b%mult = REAL(mult,r8)
    ELSE
      b%mult = 1.d0
    END IF

    nelements = 1
    DO i = 1,ndims
      b%dims(i) = INT(dims(i),i4)
      nelements = nelements * b%dims(i)
    END DO

    IF (PRESENT(rank_write)) h%rank_master = rank_write

    ! Write header

    CALL write_mesh_variable_meta_r8(h, id, name, units, mesh_id, mult)

    ! Write the actual data

    IF (h%rank == h%rank_master) THEN
      h%current_location = b%data_location

      CALL MPI_FILE_SEEK(h%filehandle, h%current_location, MPI_SEEK_SET, &
          errcode)

      IF (PRESENT(subarray)) THEN
        mpitype = subarray
        nitems = 1
      ELSE
        mpitype = b%mpitype
        nitems = nelements
      END IF

      CALL MPI_FILE_WRITE(h%filehandle, variable, nitems, mpitype, &
          MPI_STATUS_IGNORE, errcode)
    END IF

    h%current_location = b%data_location + b%data_length
    b%done_data = .TRUE.

  END SUBROUTINE write_srl_1d_integer_i8_r8



  !----------------------------------------------------------------------------
  ! Code to write a 2D cartesian integer variable in serial from the node
  ! with rank {rank_write}
  !----------------------------------------------------------------------------

  SUBROUTINE write_srl_2d_integer_i8_r8(h, id, name, units, dims, stagger, &
      mesh_id, variable, subarray, mult, rank_write)

    INTEGER, PARAMETER :: ndims = 2
    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name, units
    INTEGER, DIMENSION(:), INTENT(IN) :: dims
    INTEGER(i4), INTENT(IN) :: stagger
    CHARACTER(LEN=*), INTENT(IN) :: mesh_id
    INTEGER(i8), DIMENSION(:,:), INTENT(IN) :: variable
    INTEGER, OPTIONAL, INTENT(IN) :: subarray
    REAL(r8), OPTIONAL, INTENT(IN) :: mult
    INTEGER, OPTIONAL, INTENT(IN) :: rank_write
    INTEGER :: i, errcode, mpitype, nitems, nelements
    TYPE(sdf_block_type), POINTER :: b

    CALL sdf_get_next_block(h)
    b => h%current_block

    b%type_size = 8
    b%datatype = c_datatype_integer8
    b%mpitype = MPI_INTEGER8
    b%ndims = ndims
    b%stagger = stagger

    IF (PRESENT(mult)) THEN
      b%mult = REAL(mult,r8)
    ELSE
      b%mult = 1.d0
    END IF

    nelements = 1
    DO i = 1,ndims
      b%dims(i) = INT(dims(i),i4)
      nelements = nelements * b%dims(i)
    END DO

    IF (PRESENT(rank_write)) h%rank_master = rank_write

    ! Write header

    CALL write_mesh_variable_meta_r8(h, id, name, units, mesh_id, mult)

    ! Write the actual data

    IF (h%rank == h%rank_master) THEN
      h%current_location = b%data_location

      CALL MPI_FILE_SEEK(h%filehandle, h%current_location, MPI_SEEK_SET, &
          errcode)

      IF (PRESENT(subarray)) THEN
        mpitype = subarray
        nitems = 1
      ELSE
        mpitype = b%mpitype
        nitems = nelements
      END IF

      CALL MPI_FILE_WRITE(h%filehandle, variable, nitems, mpitype, &
          MPI_STATUS_IGNORE, errcode)
    END IF

    h%current_location = b%data_location + b%data_length
    b%done_data = .TRUE.

  END SUBROUTINE write_srl_2d_integer_i8_r8



  !----------------------------------------------------------------------------
  ! Code to write a 3D cartesian integer variable in serial from the node
  ! with rank {rank_write}
  !----------------------------------------------------------------------------

  SUBROUTINE write_srl_3d_integer_i8_r8(h, id, name, units, dims, stagger, &
      mesh_id, variable, subarray, mult, rank_write)

    INTEGER, PARAMETER :: ndims = 3
    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name, units
    INTEGER, DIMENSION(:), INTENT(IN) :: dims
    INTEGER(i4), INTENT(IN) :: stagger
    CHARACTER(LEN=*), INTENT(IN) :: mesh_id
    INTEGER(i8), DIMENSION(:,:,:), INTENT(IN) :: variable
    INTEGER, OPTIONAL, INTENT(IN) :: subarray
    REAL(r8), OPTIONAL, INTENT(IN) :: mult
    INTEGER, OPTIONAL, INTENT(IN) :: rank_write
    INTEGER :: i, errcode, mpitype, nitems, nelements
    TYPE(sdf_block_type), POINTER :: b

    CALL sdf_get_next_block(h)
    b => h%current_block

    b%type_size = 8
    b%datatype = c_datatype_integer8
    b%mpitype = MPI_INTEGER8
    b%ndims = ndims
    b%stagger = stagger

    IF (PRESENT(mult)) THEN
      b%mult = REAL(mult,r8)
    ELSE
      b%mult = 1.d0
    END IF

    nelements = 1
    DO i = 1,ndims
      b%dims(i) = INT(dims(i),i4)
      nelements = nelements * b%dims(i)
    END DO

    IF (PRESENT(rank_write)) h%rank_master = rank_write

    ! Write header

    CALL write_mesh_variable_meta_r8(h, id, name, units, mesh_id, mult)

    ! Write the actual data

    IF (h%rank == h%rank_master) THEN
      h%current_location = b%data_location

      CALL MPI_FILE_SEEK(h%filehandle, h%current_location, MPI_SEEK_SET, &
          errcode)

      IF (PRESENT(subarray)) THEN
        mpitype = subarray
        nitems = 1
      ELSE
        mpitype = b%mpitype
        nitems = nelements
      END IF

      CALL MPI_FILE_WRITE(h%filehandle, variable, nitems, mpitype, &
          MPI_STATUS_IGNORE, errcode)
    END IF

    h%current_location = b%data_location + b%data_length
    b%done_data = .TRUE.

  END SUBROUTINE write_srl_3d_integer_i8_r8



  SUBROUTINE write_srl_1d_integer_i8_r4(h, id, name, units, dims, stagger, &
      mesh_id, variable, subarray, mult, rank_write)

    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name, units
    INTEGER, DIMENSION(:), INTENT(IN) :: dims
    INTEGER(i4), INTENT(IN) :: stagger
    CHARACTER(LEN=*), INTENT(IN) :: mesh_id
    INTEGER(i8), DIMENSION(:), INTENT(IN) :: variable
    INTEGER, INTENT(IN) :: subarray
    REAL(r4), INTENT(IN) :: mult
    INTEGER, OPTIONAL, INTENT(IN) :: rank_write

    CALL write_srl_1d_integer_i8_r8(h, id, name, units, dims, stagger, &
        mesh_id, variable, subarray, REAL(mult,r8), rank_write)

  END SUBROUTINE write_srl_1d_integer_i8_r4



  SUBROUTINE write_srl_2d_integer_i8_r4(h, id, name, units, dims, stagger, &
      mesh_id, variable, subarray, mult, rank_write)

    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name, units
    INTEGER, DIMENSION(:), INTENT(IN) :: dims
    INTEGER(i4), INTENT(IN) :: stagger
    CHARACTER(LEN=*), INTENT(IN) :: mesh_id
    INTEGER(i8), DIMENSION(:,:), INTENT(IN) :: variable
    INTEGER, INTENT(IN) :: subarray
    REAL(r4), INTENT(IN) :: mult
    INTEGER, OPTIONAL, INTENT(IN) :: rank_write

    CALL write_srl_2d_integer_i8_r8(h, id, name, units, dims, stagger, &
        mesh_id, variable, subarray, REAL(mult,r8), rank_write)

  END SUBROUTINE write_srl_2d_integer_i8_r4



  SUBROUTINE write_srl_3d_integer_i8_r4(h, id, name, units, dims, stagger, &
      mesh_id, variable, subarray, mult, rank_write)

    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name, units
    INTEGER, DIMENSION(:), INTENT(IN) :: dims
    INTEGER(i4), INTENT(IN) :: stagger
    CHARACTER(LEN=*), INTENT(IN) :: mesh_id
    INTEGER(i8), DIMENSION(:,:,:), INTENT(IN) :: variable
    INTEGER, INTENT(IN) :: subarray
    REAL(r4), INTENT(IN) :: mult
    INTEGER, OPTIONAL, INTENT(IN) :: rank_write

    CALL write_srl_3d_integer_i8_r8(h, id, name, units, dims, stagger, &
        mesh_id, variable, subarray, REAL(mult,r8), rank_write)

  END SUBROUTINE write_srl_3d_integer_i8_r4



  !----------------------------------------------------------------------------
  ! Code to write a 1D cartesian integer variable in parallel
  ! using the mpitype {distribution} for distribution of data
  ! It's up to the coder to design the distribution parallel operation, so
  ! need global dims
  !----------------------------------------------------------------------------

  SUBROUTINE write_1d_integer_i4_r8(h, id, name, units, dims, stagger, &
      mesh_id, variable, distribution, subarray, convert_in, mult)

    INTEGER, PARAMETER :: ndims = 1
    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name, units
    INTEGER, DIMENSION(:), INTENT(IN) :: dims
    INTEGER(i4), INTENT(IN) :: stagger
    CHARACTER(LEN=*), INTENT(IN) :: mesh_id
    INTEGER(i4), DIMENSION(:), INTENT(IN) :: variable
    INTEGER, INTENT(IN) :: distribution, subarray
    LOGICAL, OPTIONAL, INTENT(IN) :: convert_in
    REAL(r8), OPTIONAL, INTENT(IN) :: mult
    INTEGER :: i, errcode
    TYPE(sdf_block_type), POINTER :: b

    CALL sdf_get_next_block(h)
    b => h%current_block

    b%type_size = 4
    b%datatype = c_datatype_integer4
    b%mpitype = MPI_INTEGER4
    b%ndims = ndims
    b%stagger = stagger

    IF (PRESENT(mult)) THEN
      b%mult = REAL(mult,r8)
    ELSE
      b%mult = 1.d0
    END IF

    DO i = 1,ndims
      b%dims(i) = INT(dims(i),i4)
    END DO

    ! Write header

    CALL write_mesh_variable_meta_r8(h, id, name, units, mesh_id, mult)

    ! Write the actual data

    CALL MPI_FILE_SET_VIEW(h%filehandle, h%current_location, MPI_BYTE, &
        distribution, 'native', MPI_INFO_NULL, errcode)
    CALL MPI_FILE_WRITE_ALL(h%filehandle, variable, 1, subarray, &
        MPI_STATUS_IGNORE, errcode)

    CALL MPI_FILE_SET_VIEW(h%filehandle, c_off0, MPI_BYTE, MPI_BYTE, 'native', &
        MPI_INFO_NULL, errcode)

    h%current_location = b%data_location + b%data_length
    b%done_data = .TRUE.

  END SUBROUTINE write_1d_integer_i4_r8



  !----------------------------------------------------------------------------
  ! Code to write a 2D cartesian integer variable in parallel
  ! using the mpitype {distribution} for distribution of data
  ! It's up to the coder to design the distribution parallel operation, so
  ! need global dims
  !----------------------------------------------------------------------------

  SUBROUTINE write_2d_integer_i4_r8(h, id, name, units, dims, stagger, &
      mesh_id, variable, distribution, subarray, convert_in, mult)

    INTEGER, PARAMETER :: ndims = 2
    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name, units
    INTEGER, DIMENSION(:), INTENT(IN) :: dims
    INTEGER(i4), INTENT(IN) :: stagger
    CHARACTER(LEN=*), INTENT(IN) :: mesh_id
    INTEGER(i4), DIMENSION(:,:), INTENT(IN) :: variable
    INTEGER, INTENT(IN) :: distribution, subarray
    LOGICAL, OPTIONAL, INTENT(IN) :: convert_in
    REAL(r8), OPTIONAL, INTENT(IN) :: mult
    INTEGER :: i, errcode
    TYPE(sdf_block_type), POINTER :: b

    CALL sdf_get_next_block(h)
    b => h%current_block

    b%type_size = 4
    b%datatype = c_datatype_integer4
    b%mpitype = MPI_INTEGER4
    b%ndims = ndims
    b%stagger = stagger

    IF (PRESENT(mult)) THEN
      b%mult = REAL(mult,r8)
    ELSE
      b%mult = 1.d0
    END IF

    DO i = 1,ndims
      b%dims(i) = INT(dims(i),i4)
    END DO

    ! Write header

    CALL write_mesh_variable_meta_r8(h, id, name, units, mesh_id, mult)

    ! Write the actual data

    CALL MPI_FILE_SET_VIEW(h%filehandle, h%current_location, MPI_BYTE, &
        distribution, 'native', MPI_INFO_NULL, errcode)
    CALL MPI_FILE_WRITE_ALL(h%filehandle, variable, 1, subarray, &
        MPI_STATUS_IGNORE, errcode)

    CALL MPI_FILE_SET_VIEW(h%filehandle, c_off0, MPI_BYTE, MPI_BYTE, 'native', &
        MPI_INFO_NULL, errcode)

    h%current_location = b%data_location + b%data_length
    b%done_data = .TRUE.

  END SUBROUTINE write_2d_integer_i4_r8



  !----------------------------------------------------------------------------
  ! Code to write a 3D cartesian integer variable in parallel
  ! using the mpitype {distribution} for distribution of data
  ! It's up to the coder to design the distribution parallel operation, so
  ! need global dims
  !----------------------------------------------------------------------------

  SUBROUTINE write_3d_integer_i4_r8(h, id, name, units, dims, stagger, &
      mesh_id, variable, distribution, subarray, convert_in, mult)

    INTEGER, PARAMETER :: ndims = 3
    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name, units
    INTEGER, DIMENSION(:), INTENT(IN) :: dims
    INTEGER(i4), INTENT(IN) :: stagger
    CHARACTER(LEN=*), INTENT(IN) :: mesh_id
    INTEGER(i4), DIMENSION(:,:,:), INTENT(IN) :: variable
    INTEGER, INTENT(IN) :: distribution, subarray
    LOGICAL, OPTIONAL, INTENT(IN) :: convert_in
    REAL(r8), OPTIONAL, INTENT(IN) :: mult
    INTEGER :: i, errcode
    TYPE(sdf_block_type), POINTER :: b

    CALL sdf_get_next_block(h)
    b => h%current_block

    b%type_size = 4
    b%datatype = c_datatype_integer4
    b%mpitype = MPI_INTEGER4
    b%ndims = ndims
    b%stagger = stagger

    IF (PRESENT(mult)) THEN
      b%mult = REAL(mult,r8)
    ELSE
      b%mult = 1.d0
    END IF

    DO i = 1,ndims
      b%dims(i) = INT(dims(i),i4)
    END DO

    ! Write header

    CALL write_mesh_variable_meta_r8(h, id, name, units, mesh_id, mult)

    ! Write the actual data

    CALL MPI_FILE_SET_VIEW(h%filehandle, h%current_location, MPI_BYTE, &
        distribution, 'native', MPI_INFO_NULL, errcode)
    CALL MPI_FILE_WRITE_ALL(h%filehandle, variable, 1, subarray, &
        MPI_STATUS_IGNORE, errcode)

    CALL MPI_FILE_SET_VIEW(h%filehandle, c_off0, MPI_BYTE, MPI_BYTE, 'native', &
        MPI_INFO_NULL, errcode)

    h%current_location = b%data_location + b%data_length
    b%done_data = .TRUE.

  END SUBROUTINE write_3d_integer_i4_r8



  SUBROUTINE write_1d_integer_i4_r4(h, id, name, units, dims, stagger, &
      mesh_id, variable, distribution, subarray, convert_in, mult)

    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name, units
    INTEGER, DIMENSION(:), INTENT(IN) :: dims
    INTEGER(i4), INTENT(IN) :: stagger
    CHARACTER(LEN=*), INTENT(IN) :: mesh_id
    INTEGER(i4), DIMENSION(:), INTENT(IN) :: variable
    INTEGER, INTENT(IN) :: distribution, subarray
    LOGICAL, INTENT(IN) :: convert_in
    REAL(r4), INTENT(IN) :: mult

    CALL write_1d_integer_i4_r8(h, id, name, units, dims, stagger, &
        mesh_id, variable, distribution, subarray, convert_in, REAL(mult,r8))

  END SUBROUTINE write_1d_integer_i4_r4



  SUBROUTINE write_2d_integer_i4_r4(h, id, name, units, dims, stagger, &
      mesh_id, variable, distribution, subarray, convert_in, mult)

    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name, units
    INTEGER, DIMENSION(:), INTENT(IN) :: dims
    INTEGER(i4), INTENT(IN) :: stagger
    CHARACTER(LEN=*), INTENT(IN) :: mesh_id
    INTEGER(i4), DIMENSION(:,:), INTENT(IN) :: variable
    INTEGER, INTENT(IN) :: distribution, subarray
    LOGICAL, INTENT(IN) :: convert_in
    REAL(r4), INTENT(IN) :: mult

    CALL write_2d_integer_i4_r8(h, id, name, units, dims, stagger, &
        mesh_id, variable, distribution, subarray, convert_in, REAL(mult,r8))

  END SUBROUTINE write_2d_integer_i4_r4



  SUBROUTINE write_3d_integer_i4_r4(h, id, name, units, dims, stagger, &
      mesh_id, variable, distribution, subarray, convert_in, mult)

    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name, units
    INTEGER, DIMENSION(:), INTENT(IN) :: dims
    INTEGER(i4), INTENT(IN) :: stagger
    CHARACTER(LEN=*), INTENT(IN) :: mesh_id
    INTEGER(i4), DIMENSION(:,:,:), INTENT(IN) :: variable
    INTEGER, INTENT(IN) :: distribution, subarray
    LOGICAL, INTENT(IN) :: convert_in
    REAL(r4), INTENT(IN) :: mult

    CALL write_3d_integer_i4_r8(h, id, name, units, dims, stagger, &
        mesh_id, variable, distribution, subarray, convert_in, REAL(mult,r8))

  END SUBROUTINE write_3d_integer_i4_r4



  !----------------------------------------------------------------------------
  ! Code to write a 1D cartesian integer variable in parallel
  ! using the mpitype {distribution} for distribution of data
  ! It's up to the coder to design the distribution parallel operation, so
  ! need global dims
  !----------------------------------------------------------------------------

  SUBROUTINE write_1d_integer_i8_r8(h, id, name, units, dims, stagger, &
      mesh_id, variable, distribution, subarray, convert_in, mult)

    INTEGER, PARAMETER :: ndims = 1
    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name, units
    INTEGER, DIMENSION(:), INTENT(IN) :: dims
    INTEGER(i4), INTENT(IN) :: stagger
    CHARACTER(LEN=*), INTENT(IN) :: mesh_id
    INTEGER(i8), DIMENSION(:), INTENT(IN) :: variable
    INTEGER, INTENT(IN) :: distribution, subarray
    LOGICAL, OPTIONAL, INTENT(IN) :: convert_in
    REAL(r8), OPTIONAL, INTENT(IN) :: mult
    INTEGER(i4), DIMENSION(:), ALLOCATABLE :: i4array
    INTEGER :: i, errcode
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
      b%datatype = c_datatype_integer4
      b%mpitype = MPI_INTEGER4
    ELSE
      b%type_size = 8
      b%datatype = c_datatype_integer8
      b%mpitype = MPI_INTEGER8
    END IF
    b%ndims = ndims
    b%stagger = stagger

    IF (PRESENT(mult)) THEN
      b%mult = REAL(mult,r8)
    ELSE
      b%mult = 1.d0
    END IF

    DO i = 1,ndims
      b%dims(i) = INT(dims(i),i4)
    END DO

    ! Write header

    CALL write_mesh_variable_meta_r8(h, id, name, units, mesh_id, mult)

    ! Write the actual data

    CALL MPI_FILE_SET_VIEW(h%filehandle, h%current_location, MPI_BYTE, &
        distribution, 'native', MPI_INFO_NULL, errcode)
    IF (convert) THEN
      ALLOCATE(i4array(b%dims(1)))
      i4array = INT(variable,i4)
      CALL MPI_FILE_WRITE_ALL(h%filehandle, i4array, 1, subarray, &
          MPI_STATUS_IGNORE, errcode)
      DEALLOCATE(i4array)
    ELSE
      CALL MPI_FILE_WRITE_ALL(h%filehandle, variable, 1, subarray, &
          MPI_STATUS_IGNORE, errcode)
    END IF

    CALL MPI_FILE_SET_VIEW(h%filehandle, c_off0, MPI_BYTE, MPI_BYTE, 'native', &
        MPI_INFO_NULL, errcode)

    h%current_location = b%data_location + b%data_length
    b%done_data = .TRUE.

  END SUBROUTINE write_1d_integer_i8_r8



  !----------------------------------------------------------------------------
  ! Code to write a 2D cartesian integer variable in parallel
  ! using the mpitype {distribution} for distribution of data
  ! It's up to the coder to design the distribution parallel operation, so
  ! need global dims
  !----------------------------------------------------------------------------

  SUBROUTINE write_2d_integer_i8_r8(h, id, name, units, dims, stagger, &
      mesh_id, variable, distribution, subarray, convert_in, mult)

    INTEGER, PARAMETER :: ndims = 2
    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name, units
    INTEGER, DIMENSION(:), INTENT(IN) :: dims
    INTEGER(i4), INTENT(IN) :: stagger
    CHARACTER(LEN=*), INTENT(IN) :: mesh_id
    INTEGER(i8), DIMENSION(:,:), INTENT(IN) :: variable
    INTEGER, INTENT(IN) :: distribution, subarray
    LOGICAL, OPTIONAL, INTENT(IN) :: convert_in
    REAL(r8), OPTIONAL, INTENT(IN) :: mult
    INTEGER(i4), DIMENSION(:,:), ALLOCATABLE :: i4array
    INTEGER :: i, errcode
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
      b%datatype = c_datatype_integer4
      b%mpitype = MPI_INTEGER4
    ELSE
      b%type_size = 8
      b%datatype = c_datatype_integer8
      b%mpitype = MPI_INTEGER8
    END IF
    b%ndims = ndims
    b%stagger = stagger

    IF (PRESENT(mult)) THEN
      b%mult = REAL(mult,r8)
    ELSE
      b%mult = 1.d0
    END IF

    DO i = 1,ndims
      b%dims(i) = INT(dims(i),i4)
    END DO

    ! Write header

    CALL write_mesh_variable_meta_r8(h, id, name, units, mesh_id, mult)

    ! Write the actual data

    CALL MPI_FILE_SET_VIEW(h%filehandle, h%current_location, MPI_BYTE, &
        distribution, 'native', MPI_INFO_NULL, errcode)
    IF (convert) THEN
      ALLOCATE(i4array(b%dims(1),b%dims(2)))
      i4array = INT(variable,i4)
      CALL MPI_FILE_WRITE_ALL(h%filehandle, i4array, 1, subarray, &
          MPI_STATUS_IGNORE, errcode)
      DEALLOCATE(i4array)
    ELSE
      CALL MPI_FILE_WRITE_ALL(h%filehandle, variable, 1, subarray, &
          MPI_STATUS_IGNORE, errcode)
    END IF

    CALL MPI_FILE_SET_VIEW(h%filehandle, c_off0, MPI_BYTE, MPI_BYTE, 'native', &
        MPI_INFO_NULL, errcode)

    h%current_location = b%data_location + b%data_length
    b%done_data = .TRUE.

  END SUBROUTINE write_2d_integer_i8_r8



  !----------------------------------------------------------------------------
  ! Code to write a 3D cartesian integer variable in parallel
  ! using the mpitype {distribution} for distribution of data
  ! It's up to the coder to design the distribution parallel operation, so
  ! need global dims
  !----------------------------------------------------------------------------

  SUBROUTINE write_3d_integer_i8_r8(h, id, name, units, dims, stagger, &
      mesh_id, variable, distribution, subarray, convert_in, mult)

    INTEGER, PARAMETER :: ndims = 3
    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name, units
    INTEGER, DIMENSION(:), INTENT(IN) :: dims
    INTEGER(i4), INTENT(IN) :: stagger
    CHARACTER(LEN=*), INTENT(IN) :: mesh_id
    INTEGER(i8), DIMENSION(:,:,:), INTENT(IN) :: variable
    INTEGER, INTENT(IN) :: distribution, subarray
    LOGICAL, OPTIONAL, INTENT(IN) :: convert_in
    REAL(r8), OPTIONAL, INTENT(IN) :: mult
    INTEGER(i4), DIMENSION(:,:,:), ALLOCATABLE :: i4array
    INTEGER :: i, errcode
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
      b%datatype = c_datatype_integer4
      b%mpitype = MPI_INTEGER4
    ELSE
      b%type_size = 8
      b%datatype = c_datatype_integer8
      b%mpitype = MPI_INTEGER8
    END IF
    b%ndims = ndims
    b%stagger = stagger

    IF (PRESENT(mult)) THEN
      b%mult = REAL(mult,r8)
    ELSE
      b%mult = 1.d0
    END IF

    DO i = 1,ndims
      b%dims(i) = INT(dims(i),i4)
    END DO

    ! Write header

    CALL write_mesh_variable_meta_r8(h, id, name, units, mesh_id, mult)

    ! Write the actual data

    CALL MPI_FILE_SET_VIEW(h%filehandle, h%current_location, MPI_BYTE, &
        distribution, 'native', MPI_INFO_NULL, errcode)
    IF (convert) THEN
      ALLOCATE(i4array(b%dims(1),b%dims(2),b%dims(3)))
      i4array = INT(variable,i4)
      CALL MPI_FILE_WRITE_ALL(h%filehandle, i4array, 1, subarray, &
          MPI_STATUS_IGNORE, errcode)
      DEALLOCATE(i4array)
    ELSE
      CALL MPI_FILE_WRITE_ALL(h%filehandle, variable, 1, subarray, &
          MPI_STATUS_IGNORE, errcode)
    END IF

    CALL MPI_FILE_SET_VIEW(h%filehandle, c_off0, MPI_BYTE, MPI_BYTE, 'native', &
        MPI_INFO_NULL, errcode)

    h%current_location = b%data_location + b%data_length
    b%done_data = .TRUE.

  END SUBROUTINE write_3d_integer_i8_r8



  SUBROUTINE write_1d_integer_i8_r4(h, id, name, units, dims, stagger, &
      mesh_id, variable, distribution, subarray, convert_in, mult)

    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name, units
    INTEGER, DIMENSION(:), INTENT(IN) :: dims
    INTEGER(i4), INTENT(IN) :: stagger
    CHARACTER(LEN=*), INTENT(IN) :: mesh_id
    INTEGER(i8), DIMENSION(:), INTENT(IN) :: variable
    INTEGER, INTENT(IN) :: distribution, subarray
    LOGICAL, INTENT(IN) :: convert_in
    REAL(r4), INTENT(IN) :: mult

    CALL write_1d_integer_i8_r8(h, id, name, units, dims, stagger, &
        mesh_id, variable, distribution, subarray, convert_in, REAL(mult,r8))

  END SUBROUTINE write_1d_integer_i8_r4



  SUBROUTINE write_2d_integer_i8_r4(h, id, name, units, dims, stagger, &
      mesh_id, variable, distribution, subarray, convert_in, mult)

    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name, units
    INTEGER, DIMENSION(:), INTENT(IN) :: dims
    INTEGER(i4), INTENT(IN) :: stagger
    CHARACTER(LEN=*), INTENT(IN) :: mesh_id
    INTEGER(i8), DIMENSION(:,:), INTENT(IN) :: variable
    INTEGER, INTENT(IN) :: distribution, subarray
    LOGICAL, INTENT(IN) :: convert_in
    REAL(r4), INTENT(IN) :: mult

    CALL write_2d_integer_i8_r8(h, id, name, units, dims, stagger, &
        mesh_id, variable, distribution, subarray, convert_in, REAL(mult,r8))

  END SUBROUTINE write_2d_integer_i8_r4



  SUBROUTINE write_3d_integer_i8_r4(h, id, name, units, dims, stagger, &
      mesh_id, variable, distribution, subarray, convert_in, mult)

    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name, units
    INTEGER, DIMENSION(:), INTENT(IN) :: dims
    INTEGER(i4), INTENT(IN) :: stagger
    CHARACTER(LEN=*), INTENT(IN) :: mesh_id
    INTEGER(i8), DIMENSION(:,:,:), INTENT(IN) :: variable
    INTEGER, INTENT(IN) :: distribution, subarray
    LOGICAL, INTENT(IN) :: convert_in
    REAL(r4), INTENT(IN) :: mult

    CALL write_3d_integer_i8_r8(h, id, name, units, dims, stagger, &
        mesh_id, variable, distribution, subarray, convert_in, REAL(mult,r8))

  END SUBROUTINE write_3d_integer_i8_r4



  !----------------------------------------------------------------------------
  ! Code to write a 1D cartesian character variable in parallel
  ! using the mpitype {distribution} for distribution of data
  ! It's up to the coder to design the distribution parallel operation, so
  ! need global dims
  !----------------------------------------------------------------------------

  SUBROUTINE write_1d_character_r8(h, id, name, units, dims, stagger, mesh_id, &
      variable, distribution, subarray, mult)

    INTEGER, PARAMETER :: ndims = 1
    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name, units
    INTEGER, DIMENSION(:), INTENT(IN) :: dims
    INTEGER(i4), INTENT(IN) :: stagger
    CHARACTER(LEN=*), INTENT(IN) :: mesh_id
    CHARACTER(LEN=1), DIMENSION(:), INTENT(IN) :: variable
    INTEGER, INTENT(IN) :: distribution, subarray
    REAL(r8), OPTIONAL, INTENT(IN) :: mult
    INTEGER :: i, errcode
    TYPE(sdf_block_type), POINTER :: b

    CALL sdf_get_next_block(h)
    b => h%current_block

    b%type_size = 1
    b%datatype = c_datatype_character
    b%mpitype = MPI_CHARACTER
    b%ndims = ndims
    b%stagger = stagger

    IF (PRESENT(mult)) THEN
      b%mult = REAL(mult,r8)
    ELSE
      b%mult = 1.d0
    END IF

    DO i = 1,ndims
      b%dims(i) = INT(dims(i),i4)
    END DO

    ! Write header

    CALL write_mesh_variable_meta_r8(h, id, name, units, mesh_id, mult)

    ! Write the actual data

    CALL MPI_FILE_SET_VIEW(h%filehandle, h%current_location, MPI_BYTE, &
        distribution, 'native', MPI_INFO_NULL, errcode)
    CALL MPI_FILE_WRITE_ALL(h%filehandle, variable, 1, subarray, &
        MPI_STATUS_IGNORE, errcode)

    CALL MPI_FILE_SET_VIEW(h%filehandle, c_off0, MPI_BYTE, MPI_BYTE, 'native', &
        MPI_INFO_NULL, errcode)

    h%current_location = b%data_location + b%data_length
    b%done_data = .TRUE.

  END SUBROUTINE write_1d_character_r8



  !----------------------------------------------------------------------------
  ! Code to write a 2D cartesian character variable in parallel
  ! using the mpitype {distribution} for distribution of data
  ! It's up to the coder to design the distribution parallel operation, so
  ! need global dims
  !----------------------------------------------------------------------------

  SUBROUTINE write_2d_character_r8(h, id, name, units, dims, stagger, mesh_id, &
      variable, distribution, subarray, mult)

    INTEGER, PARAMETER :: ndims = 2
    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name, units
    INTEGER, DIMENSION(:), INTENT(IN) :: dims
    INTEGER(i4), INTENT(IN) :: stagger
    CHARACTER(LEN=*), INTENT(IN) :: mesh_id
    CHARACTER(LEN=1), DIMENSION(:,:), INTENT(IN) :: variable
    INTEGER, INTENT(IN) :: distribution, subarray
    REAL(r8), OPTIONAL, INTENT(IN) :: mult
    INTEGER :: i, errcode
    TYPE(sdf_block_type), POINTER :: b

    CALL sdf_get_next_block(h)
    b => h%current_block

    b%type_size = 1
    b%datatype = c_datatype_character
    b%mpitype = MPI_CHARACTER
    b%ndims = ndims
    b%stagger = stagger

    IF (PRESENT(mult)) THEN
      b%mult = REAL(mult,r8)
    ELSE
      b%mult = 1.d0
    END IF

    DO i = 1,ndims
      b%dims(i) = INT(dims(i),i4)
    END DO

    ! Write header

    CALL write_mesh_variable_meta_r8(h, id, name, units, mesh_id, mult)

    ! Write the actual data

    CALL MPI_FILE_SET_VIEW(h%filehandle, h%current_location, MPI_BYTE, &
        distribution, 'native', MPI_INFO_NULL, errcode)
    CALL MPI_FILE_WRITE_ALL(h%filehandle, variable, 1, subarray, &
        MPI_STATUS_IGNORE, errcode)

    CALL MPI_FILE_SET_VIEW(h%filehandle, c_off0, MPI_BYTE, MPI_BYTE, 'native', &
        MPI_INFO_NULL, errcode)

    h%current_location = b%data_location + b%data_length
    b%done_data = .TRUE.

  END SUBROUTINE write_2d_character_r8



  !----------------------------------------------------------------------------
  ! Code to write a 3D cartesian character variable in parallel
  ! using the mpitype {distribution} for distribution of data
  ! It's up to the coder to design the distribution parallel operation, so
  ! need global dims
  !----------------------------------------------------------------------------

  SUBROUTINE write_3d_character_r8(h, id, name, units, dims, stagger, mesh_id, &
      variable, distribution, subarray, mult)

    INTEGER, PARAMETER :: ndims = 3
    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name, units
    INTEGER, DIMENSION(:), INTENT(IN) :: dims
    INTEGER(i4), INTENT(IN) :: stagger
    CHARACTER(LEN=*), INTENT(IN) :: mesh_id
    CHARACTER(LEN=1), DIMENSION(:,:,:), INTENT(IN) :: variable
    INTEGER, INTENT(IN) :: distribution, subarray
    REAL(r8), OPTIONAL, INTENT(IN) :: mult
    INTEGER :: i, errcode
    TYPE(sdf_block_type), POINTER :: b

    CALL sdf_get_next_block(h)
    b => h%current_block

    b%type_size = 1
    b%datatype = c_datatype_character
    b%mpitype = MPI_CHARACTER
    b%ndims = ndims
    b%stagger = stagger

    IF (PRESENT(mult)) THEN
      b%mult = REAL(mult,r8)
    ELSE
      b%mult = 1.d0
    END IF

    DO i = 1,ndims
      b%dims(i) = INT(dims(i),i4)
    END DO

    ! Write header

    CALL write_mesh_variable_meta_r8(h, id, name, units, mesh_id, mult)

    ! Write the actual data

    CALL MPI_FILE_SET_VIEW(h%filehandle, h%current_location, MPI_BYTE, &
        distribution, 'native', MPI_INFO_NULL, errcode)
    CALL MPI_FILE_WRITE_ALL(h%filehandle, variable, 1, subarray, &
        MPI_STATUS_IGNORE, errcode)

    CALL MPI_FILE_SET_VIEW(h%filehandle, c_off0, MPI_BYTE, MPI_BYTE, 'native', &
        MPI_INFO_NULL, errcode)

    h%current_location = b%data_location + b%data_length
    b%done_data = .TRUE.

  END SUBROUTINE write_3d_character_r8



  SUBROUTINE write_1d_character_r4(h, id, name, units, dims, stagger, mesh_id, &
      variable, distribution, subarray, mult)

    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name, units
    INTEGER, DIMENSION(:), INTENT(IN) :: dims
    INTEGER(i4), INTENT(IN) :: stagger
    CHARACTER(LEN=*), INTENT(IN) :: mesh_id
    CHARACTER(LEN=1), DIMENSION(:), INTENT(IN) :: variable
    INTEGER, INTENT(IN) :: distribution, subarray
    REAL(r4), INTENT(IN) :: mult

    CALL write_1d_character_r8(h, id, name, units, dims, stagger, mesh_id, &
        variable, distribution, subarray, REAL(mult,r8))

  END SUBROUTINE write_1d_character_r4



  SUBROUTINE write_2d_character_r4(h, id, name, units, dims, stagger, mesh_id, &
      variable, distribution, subarray, mult)

    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name, units
    INTEGER, DIMENSION(:), INTENT(IN) :: dims
    INTEGER(i4), INTENT(IN) :: stagger
    CHARACTER(LEN=*), INTENT(IN) :: mesh_id
    CHARACTER(LEN=1), DIMENSION(:,:), INTENT(IN) :: variable
    INTEGER, INTENT(IN) :: distribution, subarray
    REAL(r4), INTENT(IN) :: mult

    CALL write_2d_character_r8(h, id, name, units, dims, stagger, mesh_id, &
        variable, distribution, subarray, REAL(mult,r8))

  END SUBROUTINE write_2d_character_r4



  SUBROUTINE write_3d_character_r4(h, id, name, units, dims, stagger, mesh_id, &
      variable, distribution, subarray, mult)

    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name, units
    INTEGER, DIMENSION(:), INTENT(IN) :: dims
    INTEGER(i4), INTENT(IN) :: stagger
    CHARACTER(LEN=*), INTENT(IN) :: mesh_id
    CHARACTER(LEN=1), DIMENSION(:,:,:), INTENT(IN) :: variable
    INTEGER, INTENT(IN) :: distribution, subarray
    REAL(r4), INTENT(IN) :: mult

    CALL write_3d_character_r8(h, id, name, units, dims, stagger, mesh_id, &
        variable, distribution, subarray, REAL(mult,r8))

  END SUBROUTINE write_3d_character_r4



  !----------------------------------------------------------------------------
  ! Code to write a 1D cartesian character variable in serial from the node
  ! with rank {rank_write}
  !----------------------------------------------------------------------------

  SUBROUTINE write_srl_1d_character_r8(h, id, name, units, dims, stagger, &
      mesh_id, variable, subarray, mult, rank_write)

    INTEGER, PARAMETER :: ndims = 1
    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name, units
    INTEGER, DIMENSION(:), INTENT(IN) :: dims
    INTEGER(i4), INTENT(IN) :: stagger
    CHARACTER(LEN=*), INTENT(IN) :: mesh_id
    CHARACTER(LEN=1), DIMENSION(:), INTENT(IN) :: variable
    INTEGER, OPTIONAL, INTENT(IN) :: subarray
    REAL(r8), OPTIONAL, INTENT(IN) :: mult
    INTEGER, OPTIONAL, INTENT(IN) :: rank_write
    INTEGER :: i, errcode, mpitype, nitems, nelements
    TYPE(sdf_block_type), POINTER :: b

    CALL sdf_get_next_block(h)
    b => h%current_block

    b%type_size = 1
    b%datatype = c_datatype_character
    b%mpitype = MPI_CHARACTER
    b%ndims = ndims
    b%stagger = stagger

    IF (PRESENT(mult)) THEN
      b%mult = REAL(mult,r8)
    ELSE
      b%mult = 1.d0
    END IF

    nelements = 1
    DO i = 1,ndims
      b%dims(i) = INT(dims(i),i4)
      nelements = nelements * b%dims(i)
    END DO

    IF (PRESENT(rank_write)) h%rank_master = rank_write

    ! Write header

    CALL write_mesh_variable_meta_r8(h, id, name, units, mesh_id, mult)

    ! Write the actual data

    IF (h%rank == h%rank_master) THEN
      h%current_location = b%data_location

      CALL MPI_FILE_SEEK(h%filehandle, h%current_location, MPI_SEEK_SET, &
          errcode)

      IF (PRESENT(subarray)) THEN
        mpitype = subarray
        nitems = 1
      ELSE
        mpitype = b%mpitype
        nitems = nelements
      END IF

      CALL MPI_FILE_WRITE(h%filehandle, variable, nitems, mpitype, &
          MPI_STATUS_IGNORE, errcode)
    END IF

    h%current_location = b%data_location + b%data_length
    b%done_data = .TRUE.

  END SUBROUTINE write_srl_1d_character_r8



  !----------------------------------------------------------------------------
  ! Code to write a 2D cartesian character variable in serial from the node
  ! with rank {rank_write}
  !----------------------------------------------------------------------------

  SUBROUTINE write_srl_2d_character_r8(h, id, name, units, dims, stagger, &
      mesh_id, variable, subarray, mult, rank_write)

    INTEGER, PARAMETER :: ndims = 2
    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name, units
    INTEGER, DIMENSION(:), INTENT(IN) :: dims
    INTEGER(i4), INTENT(IN) :: stagger
    CHARACTER(LEN=*), INTENT(IN) :: mesh_id
    CHARACTER(LEN=1), DIMENSION(:,:), INTENT(IN) :: variable
    INTEGER, OPTIONAL, INTENT(IN) :: subarray
    REAL(r8), OPTIONAL, INTENT(IN) :: mult
    INTEGER, OPTIONAL, INTENT(IN) :: rank_write
    INTEGER :: i, errcode, mpitype, nitems, nelements
    TYPE(sdf_block_type), POINTER :: b

    CALL sdf_get_next_block(h)
    b => h%current_block

    b%type_size = 1
    b%datatype = c_datatype_character
    b%mpitype = MPI_CHARACTER
    b%ndims = ndims
    b%stagger = stagger

    IF (PRESENT(mult)) THEN
      b%mult = REAL(mult,r8)
    ELSE
      b%mult = 1.d0
    END IF

    nelements = 1
    DO i = 1,ndims
      b%dims(i) = INT(dims(i),i4)
      nelements = nelements * b%dims(i)
    END DO

    IF (PRESENT(rank_write)) h%rank_master = rank_write

    ! Write header

    CALL write_mesh_variable_meta_r8(h, id, name, units, mesh_id, mult)

    ! Write the actual data

    IF (h%rank == h%rank_master) THEN
      h%current_location = b%data_location

      CALL MPI_FILE_SEEK(h%filehandle, h%current_location, MPI_SEEK_SET, &
          errcode)

      IF (PRESENT(subarray)) THEN
        mpitype = subarray
        nitems = 1
      ELSE
        mpitype = b%mpitype
        nitems = nelements
      END IF

      CALL MPI_FILE_WRITE(h%filehandle, variable, nitems, mpitype, &
          MPI_STATUS_IGNORE, errcode)
    END IF

    h%current_location = b%data_location + b%data_length
    b%done_data = .TRUE.

  END SUBROUTINE write_srl_2d_character_r8



  !----------------------------------------------------------------------------
  ! Code to write a 3D cartesian character variable in serial from the node
  ! with rank {rank_write}
  !----------------------------------------------------------------------------

  SUBROUTINE write_srl_3d_character_r8(h, id, name, units, dims, stagger, &
      mesh_id, variable, subarray, mult, rank_write)

    INTEGER, PARAMETER :: ndims = 3
    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name, units
    INTEGER, DIMENSION(:), INTENT(IN) :: dims
    INTEGER(i4), INTENT(IN) :: stagger
    CHARACTER(LEN=*), INTENT(IN) :: mesh_id
    CHARACTER(LEN=1), DIMENSION(:,:,:), INTENT(IN) :: variable
    INTEGER, OPTIONAL, INTENT(IN) :: subarray
    REAL(r8), OPTIONAL, INTENT(IN) :: mult
    INTEGER, OPTIONAL, INTENT(IN) :: rank_write
    INTEGER :: i, errcode, mpitype, nitems, nelements
    TYPE(sdf_block_type), POINTER :: b

    CALL sdf_get_next_block(h)
    b => h%current_block

    b%type_size = 1
    b%datatype = c_datatype_character
    b%mpitype = MPI_CHARACTER
    b%ndims = ndims
    b%stagger = stagger

    IF (PRESENT(mult)) THEN
      b%mult = REAL(mult,r8)
    ELSE
      b%mult = 1.d0
    END IF

    nelements = 1
    DO i = 1,ndims
      b%dims(i) = INT(dims(i),i4)
      nelements = nelements * b%dims(i)
    END DO

    IF (PRESENT(rank_write)) h%rank_master = rank_write

    ! Write header

    CALL write_mesh_variable_meta_r8(h, id, name, units, mesh_id, mult)

    ! Write the actual data

    IF (h%rank == h%rank_master) THEN
      h%current_location = b%data_location

      CALL MPI_FILE_SEEK(h%filehandle, h%current_location, MPI_SEEK_SET, &
          errcode)

      IF (PRESENT(subarray)) THEN
        mpitype = subarray
        nitems = 1
      ELSE
        mpitype = b%mpitype
        nitems = nelements
      END IF

      CALL MPI_FILE_WRITE(h%filehandle, variable, nitems, mpitype, &
          MPI_STATUS_IGNORE, errcode)
    END IF

    h%current_location = b%data_location + b%data_length
    b%done_data = .TRUE.

  END SUBROUTINE write_srl_3d_character_r8



  SUBROUTINE write_srl_1d_character_r4(h, id, name, units, dims, stagger, &
      mesh_id, variable, subarray, mult, rank_write)

    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name, units
    INTEGER, DIMENSION(:), INTENT(IN) :: dims
    INTEGER(i4), INTENT(IN) :: stagger
    CHARACTER(LEN=*), INTENT(IN) :: mesh_id
    CHARACTER(LEN=1), DIMENSION(:), INTENT(IN) :: variable
    INTEGER, INTENT(IN) :: subarray
    REAL(r4), INTENT(IN) :: mult
    INTEGER, OPTIONAL, INTENT(IN) :: rank_write

    CALL write_srl_1d_character_r8(h, id, name, units, dims, stagger, &
        mesh_id, variable, subarray, REAL(mult,r8), rank_write)

  END SUBROUTINE write_srl_1d_character_r4



  SUBROUTINE write_srl_2d_character_r4(h, id, name, units, dims, stagger, &
      mesh_id, variable, subarray, mult, rank_write)

    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name, units
    INTEGER, DIMENSION(:), INTENT(IN) :: dims
    INTEGER(i4), INTENT(IN) :: stagger
    CHARACTER(LEN=*), INTENT(IN) :: mesh_id
    CHARACTER(LEN=1), DIMENSION(:,:), INTENT(IN) :: variable
    INTEGER, INTENT(IN) :: subarray
    REAL(r4), INTENT(IN) :: mult
    INTEGER, OPTIONAL, INTENT(IN) :: rank_write

    CALL write_srl_2d_character_r8(h, id, name, units, dims, stagger, &
        mesh_id, variable, subarray, REAL(mult,r8), rank_write)

  END SUBROUTINE write_srl_2d_character_r4



  SUBROUTINE write_srl_3d_character_r4(h, id, name, units, dims, stagger, &
      mesh_id, variable, subarray, mult, rank_write)

    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id, name, units
    INTEGER, DIMENSION(:), INTENT(IN) :: dims
    INTEGER(i4), INTENT(IN) :: stagger
    CHARACTER(LEN=*), INTENT(IN) :: mesh_id
    CHARACTER(LEN=1), DIMENSION(:,:,:), INTENT(IN) :: variable
    INTEGER, INTENT(IN) :: subarray
    REAL(r4), INTENT(IN) :: mult
    INTEGER, OPTIONAL, INTENT(IN) :: rank_write

    CALL write_srl_3d_character_r8(h, id, name, units, dims, stagger, &
        mesh_id, variable, subarray, REAL(mult,r8), rank_write)

  END SUBROUTINE write_srl_3d_character_r4

END MODULE sdf_output_cartesian_ru

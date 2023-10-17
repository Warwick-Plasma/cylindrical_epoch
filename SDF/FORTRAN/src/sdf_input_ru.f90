!
! SDF (Self-Describing Format) Fortran Library
! Copyright (c) 2011-2016, SDF Development Team
!
! Distributed under the terms of the BSD 3-clause License.
! See the LICENSE file for details.
!

MODULE sdf_input_ru

  USE sdf_common

  IMPLICIT NONE

CONTAINS

  SUBROUTINE read_header(h)

    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=4) :: sdf
    INTEGER :: errcode

    h%current_location = 0
    h%start_location = 0

    ! Read the header
    CALL read_entry_stringlen(h, sdf, 4)

    ! If this isn't SDF_MAGIC then this isn't an SDF file
    IF (sdf /= c_sdf_magic) THEN
      CALL MPI_FILE_CLOSE(h%filehandle, errcode)
      h%error_code = c_err_unsupported_file + 64 * h%nblocks
      h%handled_error = .TRUE.
      RETURN
    END IF

    CALL read_entry_int4(h, h%endianness)

    CALL read_entry_int4(h, h%file_version)

    IF (h%file_version > sdf_version) THEN
      CALL MPI_FILE_CLOSE(h%filehandle, errcode)
      h%error_code = c_err_unsupported_file + 64 * h%nblocks
      h%handled_error = .TRUE.
      RETURN
    END IF

    CALL read_entry_int4(h, h%file_revision)

    CALL read_entry_id(h, h%code_name)

    CALL read_entry_int8(h, h%first_block_location)

    CALL read_entry_int8(h, h%summary_location)

    CALL read_entry_int4(h, h%summary_size)

    CALL read_entry_int4(h, h%nblocks)

    CALL read_entry_int4(h, h%block_header_length)

    CALL read_entry_int4(h, h%step)

    CALL read_entry_real8(h, h%time)

    CALL read_entry_int4(h, h%jobid%start_seconds)

    CALL read_entry_int4(h, h%jobid%start_milliseconds)

    CALL read_entry_int4(h, h%string_length)

    CALL read_entry_int4(h, h%code_io_version)

    CALL read_entry_logical(h, h%restart_flag)

    CALL read_entry_logical(h, h%other_domains)

    CALL read_entry_logical(h, h%station_file)

  END SUBROUTINE read_header



  SUBROUTINE read_header_ru(h, step)

    TYPE(sdf_file_handle) :: h
    INTEGER, INTENT(OUT), OPTIONAL :: step
    INTEGER :: errcode

    IF (h%done_header) THEN
      IF (h%print_warnings .AND. h%rank == h%rank_master) THEN
        PRINT*,'*** WARNING ***'
        PRINT*,'SDF header already read. Ignoring extra call.'
      END IF
      RETURN
    END IF

    ALLOCATE(h%buffer(c_header_length))

    h%current_location = 0
    CALL MPI_FILE_SEEK(h%filehandle, h%current_location, MPI_SEEK_SET, errcode)
    IF (h%rank == h%rank_master) THEN
      CALL MPI_FILE_READ(h%filehandle, h%buffer, c_header_length, &
          MPI_CHARACTER, MPI_STATUS_IGNORE, errcode)
    END IF

    CALL MPI_BCAST(h%buffer, c_header_length, MPI_CHARACTER, h%rank_master, &
        h%comm, errcode)

    ! Read the header

    CALL read_header(h)

    DEALLOCATE(h%buffer)
    NULLIFY(h%buffer)

    IF (h%file_revision > sdf_revision) THEN
      IF (h%rank == h%rank_master) &
          PRINT *, 'Revision number of file is too high. Writing disabled'
      h%writing = .FALSE.
    END IF

    h%current_location = h%first_block_location
    h%done_header = .TRUE.

    IF (PRESENT(step)) step = h%step

  END SUBROUTINE read_header_ru



  SUBROUTINE read_block_header(h)

    TYPE(sdf_file_handle) :: h
    TYPE(sdf_block_type), POINTER :: b
    INTEGER(i4) :: info_length
    INTEGER :: errcode

    IF (.NOT. ASSOCIATED(h%current_block)) THEN
      IF (h%print_warnings .AND. h%rank == h%rank_master) THEN
        PRINT*,'*** WARNING ***'
        PRINT*,'SDF block not initialised. Ignoring call.'
      END IF
      RETURN
    END IF

    b => h%current_block

    IF (b%done_header) THEN
      h%current_location = b%block_start + h%block_header_length
      RETURN
    END IF

    h%current_location = b%block_start
    IF (.NOT. ASSOCIATED(h%buffer)) THEN
      CALL MPI_FILE_SEEK(h%filehandle, h%current_location, MPI_SEEK_SET, &
          errcode)
    END IF

    CALL read_entry_int8(h, b%next_block_location)

    CALL read_entry_int8(h, b%data_location)

    CALL read_entry_id(h, b%id)
    CALL add_to_hash_table(h, b)

    CALL read_entry_int8(h, b%data_length)

    CALL read_entry_int4(h, b%blocktype)

    CALL read_entry_int4(h, b%datatype)

    CALL read_entry_int4(h, b%ndims)

    CALL read_entry_string(h, b%name)

    CALL read_entry_int4(h, info_length)
    b%info_length = INT(info_length,i8)

    b%done_header = .TRUE.

    !print*,'block header: b%id: ', TRIM(b%id)
    !print*,'   b%name: ', TRIM(b%name)
    !print*,'   b%blocktype: ', b%blocktype
    !print*,'   b%next_block_location: ', b%next_block_location
    !print*,'   b%data_location: ', b%data_location
    !print*,'   b%datatype: ', b%datatype
    !print*,'   b%ndims: ', b%ndims
    !print*,'   b%data_length: ', b%data_length

  END SUBROUTINE read_block_header



  SUBROUTINE sdf_read_block_header(h, id, name, blocktype, ndims, datatype)

    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL :: id, name
    INTEGER, INTENT(OUT), OPTIONAL :: blocktype, ndims, datatype
    INTEGER :: errcode
    TYPE(sdf_block_type), POINTER :: b

    IF (sdf_check_header(h)) RETURN

    b => h%current_block

    IF (.NOT. b%done_header) THEN
      h%current_location = b%block_start

      IF (PRESENT(id)) THEN
        CALL MPI_FILE_SEEK(h%filehandle, h%current_location, MPI_SEEK_SET, &
            errcode)
      END IF

      CALL read_block_header(h)

      IF (b%datatype == c_datatype_real4) THEN
        b%type_size = 4
        b%mpitype = MPI_REAL4
      ELSE IF (b%datatype == c_datatype_real8) THEN
        b%type_size = 8
        b%mpitype = MPI_REAL8
      ELSE IF (b%datatype == c_datatype_integer4) THEN
        b%type_size = 4
        b%mpitype = MPI_INTEGER4
      ELSE IF (b%datatype == c_datatype_integer8) THEN
        b%type_size = 8
        b%mpitype = MPI_INTEGER8
      ELSE IF (b%datatype == c_datatype_character) THEN
        b%type_size = 1
        b%mpitype = MPI_CHARACTER
      ELSE IF (b%datatype == c_datatype_logical) THEN
        b%type_size = 1
        b%mpitype = MPI_CHARACTER
      END IF

      b%done_header = .TRUE.
    END IF

    IF (PRESENT(id)) CALL sdf_safe_copy_string(b%id, id)
    IF (PRESENT(name)) CALL sdf_safe_copy_string(b%name, name)
    IF (PRESENT(blocktype)) blocktype = b%blocktype
    IF (PRESENT(ndims)) ndims = b%ndims
    IF (PRESENT(datatype)) datatype = b%datatype

    h%current_location = b%block_start + h%block_header_length

  END SUBROUTINE sdf_read_block_header



  SUBROUTINE sdf_read_next_block_header(h, id, name, blocktype, ndims, datatype)

    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL :: id, name
    INTEGER, INTENT(OUT), OPTIONAL :: blocktype, ndims, datatype

    IF (sdf_check_header(h)) RETURN

    CALL sdf_get_next_block(h)

    CALL sdf_read_block_header(h, id, name, blocktype, ndims, datatype)

  END SUBROUTINE sdf_read_next_block_header



  SUBROUTINE read_run_info_minor(h, version, revision, minor_rev, commit_id, &
      sha1sum, compile_machine, compile_flags, defines, compile_date, &
      run_date, io_date)

    TYPE(sdf_file_handle) :: h
    INTEGER(i4), INTENT(OUT) :: version, revision, minor_rev
    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL :: commit_id, sha1sum
    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL :: compile_machine, compile_flags
    INTEGER(i8), INTENT(OUT), OPTIONAL :: defines
    INTEGER(i4), INTENT(OUT), OPTIONAL :: compile_date, run_date, io_date
    TYPE(sdf_block_type), POINTER :: b

    IF (sdf_info_init(h)) RETURN

    b => h%current_block

    IF (.NOT. b%done_info) THEN
      ! Metadata is
      ! - version   INTEGER(i4)
      ! - revision  INTEGER(i4)
      ! - commit_id CHARACTER(string_length)
      ! - sha1sum   CHARACTER(string_length)
      ! - compmac   CHARACTER(string_length)
      ! - compflag  CHARACTER(string_length)
      ! - defines   INTEGER(i8)
      ! - compdate  INTEGER(i4)
      ! - rundate   INTEGER(i4)
      ! - iodate    INTEGER(i4)
      ! - minor_rev INTEGER(i4)

      IF (.NOT. ASSOCIATED(b%run)) ALLOCATE(b%run)

      CALL read_entry_int4(h, b%run%version)

      CALL read_entry_int4(h, b%run%revision)

      CALL read_entry_string(h, b%run%commit_id)

      CALL read_entry_string(h, b%run%sha1sum)

      CALL read_entry_string(h, b%run%compile_machine)

      CALL read_entry_string(h, b%run%compile_flags)

      CALL read_entry_int8(h, b%run%defines)

      CALL read_entry_int4(h, b%run%compile_date)

      CALL read_entry_int4(h, b%run%run_date)

      CALL read_entry_int4(h, b%run%io_date)

      IF (h%file_version == 1 .AND. h%file_revision < 2) THEN
        b%run%minor_rev = 0
      ELSE
        CALL read_entry_int4(h, b%run%minor_rev)
      END IF
    END IF

    version = b%run%version
    revision = b%run%revision
    minor_rev = b%run%minor_rev
    IF (PRESENT(commit_id)) &
        CALL sdf_safe_copy_string(b%run%commit_id, commit_id)
    IF (PRESENT(sha1sum)) CALL sdf_safe_copy_string(b%run%sha1sum, sha1sum)
    IF (PRESENT(compile_machine)) &
        CALL sdf_safe_copy_string(b%run%compile_machine, compile_machine)
    IF (PRESENT(compile_flags)) &
        CALL sdf_safe_copy_string(b%run%compile_flags, compile_flags)
    IF (PRESENT(defines)) defines = b%run%defines
    IF (PRESENT(compile_date)) compile_date = b%run%compile_date
    IF (PRESENT(run_date)) run_date = b%run%run_date
    IF (PRESENT(io_date)) io_date = b%run%io_date

    b%done_info = .TRUE.
    b%done_data = .TRUE.

  END SUBROUTINE read_run_info_minor



  SUBROUTINE read_run_info_old(h, version, revision, commit_id, &
      sha1sum, compile_machine, compile_flags, defines, compile_date, &
      run_date, io_date)

    TYPE(sdf_file_handle) :: h
    INTEGER(i4), INTENT(OUT) :: version, revision
    CHARACTER(LEN=*), INTENT(OUT) :: commit_id
    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL :: sha1sum
    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL :: compile_machine, compile_flags
    INTEGER(i8), INTENT(OUT), OPTIONAL :: defines
    INTEGER(i4), INTENT(OUT), OPTIONAL :: compile_date, run_date, io_date
    INTEGER(i4) :: minor_rev

    CALL read_run_info_minor(h, version, revision, minor_rev, commit_id, &
        sha1sum, compile_machine, compile_flags, defines, compile_date, &
        run_date, io_date)

  END SUBROUTINE read_run_info_old



  SUBROUTINE read_run_info(h)

    TYPE(sdf_file_handle) :: h
    INTEGER(i4) :: version, revision, minor_rev

    CALL read_run_info_minor(h, version, revision, minor_rev)

  END SUBROUTINE read_run_info



  SUBROUTINE sdf_read_constant(h)

    TYPE(sdf_file_handle) :: h
    INTEGER :: errcode
    INTEGER(i4) :: int4
    INTEGER(i8) :: int8
    REAL(r4) :: real4
    REAL(r8) :: real8
    LOGICAL :: logic
    TYPE(sdf_block_type), POINTER :: b

    IF (sdf_check_block_header(h)) RETURN

    b => h%current_block
    IF (b%done_data) RETURN

    CALL read_block_header(h)

    h%current_location = b%block_start + h%block_header_length

    IF (.NOT. ASSOCIATED(h%buffer)) THEN
      CALL MPI_FILE_SEEK(h%filehandle, h%current_location, MPI_SEEK_SET, &
          errcode)
    END IF

    IF (b%datatype == c_datatype_integer4) THEN
      CALL read_entry_int4(h, int4)
      b%const_value(1:soi4) = TRANSFER(int4, b%const_value(1:soi4))
    ELSE IF (b%datatype == c_datatype_integer8) THEN
      CALL read_entry_int8(h, int8)
      b%const_value(1:soi8) = TRANSFER(int8, b%const_value(1:soi8))
    ELSE IF (b%datatype == c_datatype_real4) THEN
      CALL read_entry_real4(h, real4)
      b%const_value(1:sof4) = TRANSFER(real4, b%const_value(1:sof4))
    ELSE IF (b%datatype == c_datatype_real8) THEN
      CALL read_entry_real8(h, real8)
      b%const_value(1:sof8) = TRANSFER(real8, b%const_value(1:sof8))
    ELSE IF (b%datatype == c_datatype_logical) THEN
      CALL read_entry_logical(h, logic)
      b%const_value(1:sol) = TRANSFER(logic, b%const_value(1:sol))
    END IF

    b%done_info = .TRUE.
    b%done_data = .TRUE.

  END SUBROUTINE sdf_read_constant



  SUBROUTINE read_constant_integer(h, value)

    TYPE(sdf_file_handle) :: h
    INTEGER, INTENT(OUT) :: value
    INTEGER(i4) :: integer4
    INTEGER(i8) :: integer8
    TYPE(sdf_block_type), POINTER :: b

    CALL sdf_read_constant(h)

    b => h%current_block

    IF (b%datatype == c_datatype_integer4) THEN
      integer4 = TRANSFER(b%const_value, integer4)
      value = INT(integer4)
    ELSE IF (b%datatype == c_datatype_integer8) THEN
      integer8 = TRANSFER(b%const_value, integer8)
      value = INT(integer8)
    END IF

  END SUBROUTINE read_constant_integer



  SUBROUTINE read_constant_logical(h, value)

    TYPE(sdf_file_handle) :: h
    LOGICAL, INTENT(OUT) :: value
    TYPE(sdf_block_type), POINTER :: b

    CALL sdf_read_constant(h)

    b => h%current_block

    value = TRANSFER(b%const_value, value)

  END SUBROUTINE read_constant_logical



  SUBROUTINE sdf_read_array_info(h, dims)

    TYPE(sdf_file_handle) :: h
    INTEGER, DIMENSION(:), INTENT(OUT), OPTIONAL :: dims
    TYPE(sdf_block_type), POINTER :: b

    IF (sdf_info_init(h)) RETURN

    ! Metadata is
    ! - dims      ndims*INTEGER(i4)

    b => h%current_block
    IF (b%done_info) THEN
      IF (PRESENT(dims)) dims(1:b%ndims) = b%dims(1:b%ndims)
      RETURN
    END IF

    ! Size of array
    CALL read_entry_array_int4(h, b%dims, INT(b%ndims))

    IF (PRESENT(dims)) dims(1:b%ndims) = b%dims(1:b%ndims)

    b%done_info = .TRUE.

  END SUBROUTINE sdf_read_array_info



  SUBROUTINE read_1d_array_integer(h, values)

    TYPE(sdf_file_handle) :: h
    INTEGER, DIMENSION(:), INTENT(OUT) :: values
    INTEGER, DIMENSION(c_maxdims) :: dims
    INTEGER :: errcode, n1
    TYPE(sdf_block_type), POINTER :: b

    IF (sdf_check_block_header(h)) RETURN

    b => h%current_block
    IF (.NOT. b%done_info) CALL sdf_read_array_info(h, dims)

    n1 = b%dims(1)

    h%current_location = b%data_location

    IF (h%rank == h%rank_master) THEN
      CALL MPI_FILE_READ_AT(h%filehandle, h%current_location, values, n1, &
          b%mpitype, MPI_STATUS_IGNORE, errcode)
    END IF

    CALL MPI_BCAST(values, n1, b%mpitype, h%rank_master, h%comm, errcode)

    h%current_location = b%next_block_location
    b%done_data = .TRUE.

  END SUBROUTINE read_1d_array_integer



  SUBROUTINE read_2d_array_integer(h, values)

    TYPE(sdf_file_handle) :: h
    INTEGER, DIMENSION(:,:), INTENT(OUT) :: values
    INTEGER, DIMENSION(c_maxdims) :: dims
    INTEGER :: errcode, n1
    TYPE(sdf_block_type), POINTER :: b

    IF (sdf_check_block_header(h)) RETURN

    b => h%current_block
    IF (.NOT. b%done_info) CALL sdf_read_array_info(h, dims)

    n1 = b%dims(1) * b%dims(2)

    h%current_location = b%data_location

    IF (h%rank == h%rank_master) THEN
      CALL MPI_FILE_READ_AT(h%filehandle, h%current_location, values, n1, &
          b%mpitype, MPI_STATUS_IGNORE, errcode)
    END IF

    CALL MPI_BCAST(values, n1, b%mpitype, h%rank_master, h%comm, errcode)

    h%current_location = b%next_block_location
    b%done_data = .TRUE.

  END SUBROUTINE read_2d_array_integer



  SUBROUTINE read_1d_array_integer8(h, values)

    TYPE(sdf_file_handle) :: h
    INTEGER(i8), DIMENSION(:), INTENT(OUT) :: values
    INTEGER, DIMENSION(c_maxdims) :: dims
    INTEGER :: errcode, n1
    TYPE(sdf_block_type), POINTER :: b

    IF (sdf_check_block_header(h)) RETURN

    b => h%current_block
    IF (.NOT. b%done_info) CALL sdf_read_array_info(h, dims)

    n1 = b%dims(1)

    h%current_location = b%data_location

    IF (h%rank == h%rank_master) THEN
      CALL MPI_FILE_READ_AT(h%filehandle, h%current_location, values, n1, &
          b%mpitype, MPI_STATUS_IGNORE, errcode)
    END IF

    CALL MPI_BCAST(values, n1, b%mpitype, h%rank_master, h%comm, errcode)

    h%current_location = b%next_block_location
    b%done_data = .TRUE.

  END SUBROUTINE read_1d_array_integer8



  SUBROUTINE read_1d_array_logical(h, values)

    TYPE(sdf_file_handle) :: h
    LOGICAL, DIMENSION(:), INTENT(OUT) :: values
    CHARACTER(LEN=1), DIMENSION(:), ALLOCATABLE :: cvalues
    INTEGER, DIMENSION(c_maxdims) :: dims
    INTEGER :: errcode, i, n1
    TYPE(sdf_block_type), POINTER :: b

    IF (sdf_check_block_header(h)) RETURN

    b => h%current_block
    IF (.NOT. b%done_info) CALL sdf_read_array_info(h, dims)

    n1 = b%dims(1)

    ALLOCATE(cvalues(n1))

    h%current_location = b%data_location

    IF (h%rank == h%rank_master) THEN
      CALL MPI_FILE_READ_AT(h%filehandle, h%current_location, cvalues, n1, &
          b%mpitype, MPI_STATUS_IGNORE, errcode)
    END IF

    CALL MPI_BCAST(cvalues, n1, b%mpitype, h%rank_master, h%comm, errcode)

    DO i = 1,n1
      IF (cvalues(i) == ACHAR(0)) THEN
        values(i) = .FALSE.
      ELSE
        values(i) = .TRUE.
      END IF
    END DO

    DEALLOCATE(cvalues)

    h%current_location = b%next_block_location
    b%done_data = .TRUE.

  END SUBROUTINE read_1d_array_logical



  SUBROUTINE read_2d_array_character(h, values)

    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), DIMENSION(:), INTENT(OUT) :: values
    INTEGER, DIMENSION(c_maxdims) :: dims
    INTEGER :: errcode, i, n1, n2
    TYPE(sdf_block_type), POINTER :: b

    IF (sdf_check_block_header(h)) RETURN

    b => h%current_block
    IF (.NOT. b%done_info) CALL sdf_read_array_info(h, dims)

    n1 = b%dims(1)
    n2 = b%dims(2)

    h%current_location = b%data_location

    IF (h%rank == h%rank_master) THEN
      CALL MPI_FILE_SEEK(h%filehandle, h%current_location, MPI_SEEK_SET, &
          errcode)

      DO i = 1,n2
        CALL MPI_FILE_READ(h%filehandle, values(i), n1, b%mpitype, &
            MPI_STATUS_IGNORE, errcode)
      END DO
    END IF

    DO i = 1,n2
      CALL MPI_BCAST(values(i), n1, b%mpitype, h%rank_master, h%comm, errcode)
    END DO

    h%current_location = b%next_block_location
    b%done_data = .TRUE.

  END SUBROUTINE read_2d_array_character



  SUBROUTINE sdf_read_cpu_split_info(h, dims, geometry)

    TYPE(sdf_file_handle) :: h
    INTEGER, INTENT(OUT), OPTIONAL :: dims(:), geometry
    TYPE(sdf_block_type), POINTER :: b

    IF (sdf_info_init(h)) RETURN

    b => h%current_block
    IF (.NOT.b%done_info) THEN
      CALL read_entry_int4(h, b%geometry)
      CALL read_entry_array_int4(h, b%dims, INT(b%ndims))
    END IF

    IF (PRESENT(dims)) dims(1:b%ndims) = b%dims(1:b%ndims)
    IF (PRESENT(geometry)) geometry = b%geometry

    b%done_info = .TRUE.

  END SUBROUTINE sdf_read_cpu_split_info



  SUBROUTINE read_srl_cpu_split_part(h, part)

    TYPE(sdf_file_handle) :: h
    INTEGER(i8), INTENT(OUT) :: part(:)
    INTEGER :: errcode, n1
    TYPE(sdf_block_type), POINTER :: b

    IF (sdf_check_block_header(h)) RETURN

    b => h%current_block
    IF (.NOT. b%done_info) CALL sdf_read_cpu_split_info(h)

    n1 = b%dims(1)

    h%current_location = b%data_location

    IF (h%rank == h%rank_master) THEN
      CALL MPI_FILE_READ_AT(h%filehandle, h%current_location, part, n1, &
          b%mpitype, MPI_STATUS_IGNORE, errcode)
    END IF

    CALL MPI_BCAST(part, n1, b%mpitype, h%rank_master, h%comm, errcode)

    h%current_location = b%next_block_location
    b%done_data = .TRUE.

  END SUBROUTINE read_srl_cpu_split_part



  SUBROUTINE read_srl_cpu_split(h, x, y, z)

    TYPE(sdf_file_handle) :: h
    INTEGER, INTENT(OUT) :: x(:)
    INTEGER, INTENT(OUT), OPTIONAL :: y(:), z(:)
    INTEGER :: errcode, n1
    TYPE(sdf_block_type), POINTER :: b

    IF (sdf_check_block_header(h)) RETURN

    b => h%current_block
    IF (.NOT. b%done_info) CALL sdf_read_cpu_split_info(h)

    n1 = b%dims(1)

    h%current_location = b%data_location

    IF (h%rank == h%rank_master) THEN
      CALL MPI_FILE_SEEK(h%filehandle, h%current_location, MPI_SEEK_SET, &
          errcode)
      CALL MPI_FILE_READ(h%filehandle, x, n1, b%mpitype, &
          MPI_STATUS_IGNORE, errcode)
    END IF

    CALL MPI_BCAST(x, n1, b%mpitype, h%rank_master, h%comm, errcode)

    IF (PRESENT(y)) THEN
      n1 = b%dims(2)

      IF (h%rank == h%rank_master) THEN
        CALL MPI_FILE_READ(h%filehandle, y, n1, b%mpitype, &
            MPI_STATUS_IGNORE, errcode)
      END IF

      CALL MPI_BCAST(y, n1, b%mpitype, h%rank_master, h%comm, errcode)
    END IF

    IF (PRESENT(z)) THEN
      n1 = b%dims(3)

      IF (h%rank == h%rank_master) THEN
        CALL MPI_FILE_READ(h%filehandle, z, n1, b%mpitype, &
            MPI_STATUS_IGNORE, errcode)
      END IF

      CALL MPI_BCAST(z, n1, b%mpitype, h%rank_master, h%comm, errcode)
    END IF

    h%current_location = b%next_block_location
    b%done_data = .TRUE.

  END SUBROUTINE read_srl_cpu_split



  SUBROUTINE sdf_read_stitched(h)

    TYPE(sdf_file_handle) :: h
    INTEGER :: iloop, errcode
    TYPE(sdf_block_type), POINTER :: b

    IF (sdf_check_block_header(h)) RETURN

    b => h%current_block
    IF (b%done_data) RETURN

    CALL read_block_header(h)

    IF (.NOT. ASSOCIATED(h%buffer)) THEN
      CALL MPI_FILE_SEEK(h%filehandle, h%current_location, MPI_SEEK_SET, &
          errcode)
    END IF

    ! Metadata is
    ! - stagger   INTEGER(i4)
    ! - meshid    CHARACTER(id_length)
    ! - varids    ndims*CHARACTER(id_length)

    CALL read_entry_int4(h, b%stagger)

    CALL read_entry_id(h, b%mesh_id)

    ALLOCATE(b%variable_ids(b%ndims))
    DO iloop = 1, b%ndims
      CALL read_entry_id(h, b%variable_ids(iloop))
    END DO

    b%done_info = .TRUE.
    b%done_data = .TRUE.

  END SUBROUTINE sdf_read_stitched



  SUBROUTINE sdf_read_stitched_material(h)

    TYPE(sdf_file_handle) :: h
    INTEGER :: iloop, errcode
    TYPE(sdf_block_type), POINTER :: b

    IF (sdf_check_block_header(h)) RETURN

    b => h%current_block
    IF (b%done_data) RETURN

    CALL read_block_header(h)

    IF (.NOT. ASSOCIATED(h%buffer)) THEN
      CALL MPI_FILE_SEEK(h%filehandle, h%current_location, MPI_SEEK_SET, &
          errcode)
    END IF

    ! Metadata is
    ! - stagger   INTEGER(i4)
    ! - meshid    CHARACTER(id_length)
    ! - material_names ndims*CHARACTER(string_length)
    ! - varids    ndims*CHARACTER(id_length)

    CALL read_entry_int4(h, b%stagger)

    CALL read_entry_id(h, b%mesh_id)

    ALLOCATE(b%material_names(b%ndims))
    DO iloop = 1, b%ndims
      CALL read_entry_string(h, b%material_names(iloop))
    END DO

    ALLOCATE(b%variable_ids(b%ndims))
    DO iloop = 1, b%ndims
      CALL read_entry_id(h, b%variable_ids(iloop))
    END DO

    b%done_info = .TRUE.
    b%done_data = .TRUE.

  END SUBROUTINE sdf_read_stitched_material



  SUBROUTINE sdf_read_stitched_matvar(h)

    TYPE(sdf_file_handle) :: h
    INTEGER :: iloop, errcode
    TYPE(sdf_block_type), POINTER :: b

    IF (sdf_check_block_header(h)) RETURN

    b => h%current_block
    IF (b%done_data) RETURN

    CALL read_block_header(h)

    IF (.NOT. ASSOCIATED(h%buffer)) THEN
      CALL MPI_FILE_SEEK(h%filehandle, h%current_location, MPI_SEEK_SET, &
          errcode)
    END IF

    ! Metadata is
    ! - stagger   INTEGER(i4)
    ! - meshid    CHARACTER(id_length)
    ! - matid     CHARACTER(id_length)
    ! - varids    ndims*CHARACTER(id_length)

    CALL read_entry_int4(h, b%stagger)

    CALL read_entry_id(h, b%mesh_id)

    CALL read_entry_id(h, b%material_id)

    ALLOCATE(b%variable_ids(b%ndims))
    DO iloop = 1, b%ndims
      CALL read_entry_id(h, b%variable_ids(iloop))
    END DO

    b%done_info = .TRUE.
    b%done_data = .TRUE.

  END SUBROUTINE sdf_read_stitched_matvar



  SUBROUTINE sdf_read_stitched_species(h)

    TYPE(sdf_file_handle) :: h
    INTEGER :: iloop, errcode
    TYPE(sdf_block_type), POINTER :: b

    IF (sdf_check_block_header(h)) RETURN

    b => h%current_block
    IF (b%done_data) RETURN

    CALL read_block_header(h)

    IF (.NOT. ASSOCIATED(h%buffer)) THEN
      CALL MPI_FILE_SEEK(h%filehandle, h%current_location, MPI_SEEK_SET, &
          errcode)
    END IF

    ! Metadata is
    ! - stagger   INTEGER(i4)
    ! - meshid    CHARACTER(id_length)
    ! - matid     CHARACTER(id_length)
    ! - matname   CHARACTER(string_length)
    ! - specnames ndims*CHARACTER(string_length)
    ! - varids    ndims*CHARACTER(id_length)

    CALL read_entry_int4(h, b%stagger)

    CALL read_entry_id(h, b%mesh_id)

    CALL read_entry_id(h, b%material_id)

    CALL read_entry_string(h, b%material_name)

    ALLOCATE(b%material_names(b%ndims))
    DO iloop = 1, b%ndims
      CALL read_entry_string(h, b%material_names(iloop))
    END DO

    ALLOCATE(b%variable_ids(b%ndims))
    DO iloop = 1, b%ndims
      CALL read_entry_id(h, b%variable_ids(iloop))
    END DO

    b%done_info = .TRUE.
    b%done_data = .TRUE.

  END SUBROUTINE sdf_read_stitched_species



  SUBROUTINE sdf_read_stitched_obstacle_group(h)

    TYPE(sdf_file_handle) :: h
    INTEGER :: iloop
    TYPE(sdf_block_type), POINTER :: b

    IF (sdf_info_init(h)) RETURN

    b => h%current_block
    IF (b%done_data) RETURN

    ! Metadata is
    ! - stagger   INTEGER(i4)
    ! - obstacle_id    CHARACTER(id_length)
    ! - vfm_id         CHARACTER(id_length)
    ! - obstacle_names ndims*CHARACTER(string_length)

    CALL read_entry_int4(h, b%stagger)

    CALL read_entry_id(h, b%obstacle_id)

    CALL read_entry_id(h, b%vfm_id)

    ALLOCATE(b%material_names(b%ndims))
    DO iloop = 1, b%ndims
      CALL read_entry_string(h, b%material_names(iloop))
    END DO

    b%done_info = .TRUE.
    b%done_data = .TRUE.

  END SUBROUTINE sdf_read_stitched_obstacle_group



  SUBROUTINE sdf_read_obstacle_group_info(h, material_names)

    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), DIMENSION(:), INTENT(OUT) :: material_names
    INTEGER :: iloop
    TYPE(sdf_block_type), POINTER :: b

    IF (sdf_check_block_header(h)) RETURN

    CALL sdf_read_stitched_obstacle_group(h)
    b => h%current_block

    DO iloop = 1, b%ndims
      CALL sdf_safe_copy_string(b%material_names(iloop), material_names(iloop))
    END DO

  END SUBROUTINE sdf_read_obstacle_group_info



  SUBROUTINE read_entry_int4(h, value)

    INTEGER, PARAMETER :: n = soi4
    TYPE(sdf_file_handle) :: h
    INTEGER(i4), INTENT(OUT) :: value
    INTEGER(i8) :: i
    INTEGER :: j, errcode, mpitype = MPI_INTEGER4
    CHARACTER(LEN=n) :: buf

    IF (ASSOCIATED(h%buffer)) THEN
      i = h%current_location - h%start_location + 1
      DO j = 1,n
        buf(j:j) = h%buffer(i)
        i = i + 1
      END DO
      value = TRANSFER(buf, value)
    ELSE
      IF (h%rank == h%rank_master) THEN
        CALL MPI_FILE_READ(h%filehandle, value, 1, mpitype, &
            MPI_STATUS_IGNORE, errcode)
      END IF
      CALL MPI_BCAST(value, 1, mpitype, h%rank_master, h%comm, errcode)
    END IF

    h%current_location = h%current_location + n

  END SUBROUTINE read_entry_int4



  SUBROUTINE read_entry_int8(h, value)

    INTEGER, PARAMETER :: n = soi8
    TYPE(sdf_file_handle) :: h
    INTEGER(i8), INTENT(OUT) :: value
    INTEGER(i8) :: i
    INTEGER :: j, errcode, mpitype = MPI_INTEGER8
    CHARACTER(LEN=n) :: buf

    IF (ASSOCIATED(h%buffer)) THEN
      i = h%current_location - h%start_location + 1
      DO j = 1,n
        buf(j:j) = h%buffer(i)
        i = i + 1
      END DO
      value = TRANSFER(buf, value)
    ELSE
      IF (h%rank == h%rank_master) THEN
        CALL MPI_FILE_READ(h%filehandle, value, 1, mpitype, &
            MPI_STATUS_IGNORE, errcode)
      END IF
      CALL MPI_BCAST(value, 1, mpitype, h%rank_master, h%comm, errcode)
    END IF

    h%current_location = h%current_location + n

  END SUBROUTINE read_entry_int8



  SUBROUTINE read_entry_real4(h, value)

    INTEGER, PARAMETER :: n = sof4
    TYPE(sdf_file_handle) :: h
    REAL(r4), INTENT(OUT) :: value
    INTEGER(i8) :: i
    INTEGER :: j, errcode, mpitype = MPI_REAL4
    CHARACTER(LEN=n) :: buf

    IF (ASSOCIATED(h%buffer)) THEN
      i = h%current_location - h%start_location + 1
      DO j = 1,n
        buf(j:j) = h%buffer(i)
        i = i + 1
      END DO
      value = TRANSFER(buf, value)
    ELSE
      IF (h%rank == h%rank_master) THEN
        CALL MPI_FILE_READ(h%filehandle, value, 1, mpitype, &
            MPI_STATUS_IGNORE, errcode)
      END IF
      CALL MPI_BCAST(value, 1, mpitype, h%rank_master, h%comm, errcode)
    END IF

    h%current_location = h%current_location + n

  END SUBROUTINE read_entry_real4



  SUBROUTINE read_entry_real8(h, value)

    INTEGER, PARAMETER :: n = sof8
    TYPE(sdf_file_handle) :: h
    REAL(r8), INTENT(OUT) :: value
    INTEGER(i8) :: i
    INTEGER :: j, errcode, mpitype = MPI_REAL8
    CHARACTER(LEN=n) :: buf

    IF (ASSOCIATED(h%buffer)) THEN
      i = h%current_location - h%start_location + 1
      DO j = 1,n
        buf(j:j) = h%buffer(i)
        i = i + 1
      END DO
      value = TRANSFER(buf, value)
    ELSE
      IF (h%rank == h%rank_master) THEN
        CALL MPI_FILE_READ(h%filehandle, value, 1, mpitype, &
            MPI_STATUS_IGNORE, errcode)
      END IF
      CALL MPI_BCAST(value, 1, mpitype, h%rank_master, h%comm, errcode)
    END IF

    h%current_location = h%current_location + n

  END SUBROUTINE read_entry_real8



  SUBROUTINE read_entry_logical(h, ovalue)

    INTEGER, PARAMETER :: n = 1
    TYPE(sdf_file_handle) :: h
    LOGICAL, INTENT(OUT) :: ovalue
    INTEGER(i8) :: i
    INTEGER :: errcode, mpitype = MPI_CHARACTER
    CHARACTER(LEN=n) :: value

    IF (ASSOCIATED(h%buffer)) THEN
      i = h%current_location - h%start_location + 1
      value(1:1) = h%buffer(i)
    ELSE
      IF (h%rank == h%rank_master) THEN
        CALL MPI_FILE_READ(h%filehandle, value, 1, mpitype, &
            MPI_STATUS_IGNORE, errcode)
      END IF
      CALL MPI_BCAST(value, 1, mpitype, h%rank_master, h%comm, errcode)
    END IF

    IF (value(1:1) == ACHAR(1)) THEN
      ovalue = .TRUE.
    ELSE
      ovalue = .FALSE.
    END IF

    h%current_location = h%current_location + n

  END SUBROUTINE read_entry_logical



  SUBROUTINE read_entry_stringlen(h, value, n)

    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(OUT) :: value
    INTEGER, INTENT(IN) :: n
    INTEGER(i8) :: i
    INTEGER :: j, idx, errcode, mpitype = MPI_CHARACTER

    idx = 1

    IF (ASSOCIATED(h%buffer)) THEN
      i = h%current_location - h%start_location + 1
      DO j = 1,n
        value(j:j) = h%buffer(i+j-1)
        IF (value(j:j) == ACHAR(0)) EXIT
        idx = idx + 1
      END DO
    ELSE
      IF (h%rank == h%rank_master) THEN
        CALL MPI_FILE_READ(h%filehandle, value, n, mpitype, &
            MPI_STATUS_IGNORE, errcode)
      END IF
      CALL MPI_BCAST(value, n, mpitype, h%rank_master, h%comm, errcode)
      DO j = 1,n
        IF (value(j:j) == ACHAR(0)) EXIT
        idx = idx + 1
      END DO
    END IF

    DO j = idx,n
      value(j:j) = ' '
    END DO

    h%current_location = h%current_location + n

  END SUBROUTINE read_entry_stringlen



  SUBROUTINE read_entry_string(h, value)

    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(OUT) :: value

    CALL read_entry_stringlen(h, value, INT(h%string_length))

  END SUBROUTINE read_entry_string



  SUBROUTINE read_entry_id(h, value)

    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(OUT) :: value

    CALL read_entry_stringlen(h, value, INT(c_id_length))

  END SUBROUTINE read_entry_id



  SUBROUTINE read_entry_array_int4(h, value, nentries)

    INTEGER, PARAMETER :: n = soi4
    TYPE(sdf_file_handle) :: h
    INTEGER(i4), INTENT(OUT) :: value(:)
    INTEGER, INTENT(IN) :: nentries
    INTEGER(i8) :: i
    INTEGER :: j, k, errcode, mpitype = MPI_INTEGER4
    CHARACTER(LEN=n) :: buf

    IF (ASSOCIATED(h%buffer)) THEN
      i = h%current_location - h%start_location + 1
      DO j = 1, nentries
        DO k = 1,n
          buf(k:k) = h%buffer(i)
          i = i + 1
        END DO
        value(j) = TRANSFER(buf, value(1))
      END DO
    ELSE
      IF (h%rank == h%rank_master) THEN
        CALL MPI_FILE_READ(h%filehandle, value, nentries, mpitype, &
            MPI_STATUS_IGNORE, errcode)
      END IF
      CALL MPI_BCAST(value, nentries, mpitype, h%rank_master, h%comm, errcode)
    END IF

    h%current_location = h%current_location + n * nentries

  END SUBROUTINE read_entry_array_int4



  SUBROUTINE read_entry_array_int8(h, value, nentries)

    INTEGER, PARAMETER :: n = soi8
    TYPE(sdf_file_handle) :: h
    INTEGER(i8), INTENT(OUT) :: value(:)
    INTEGER, INTENT(IN) :: nentries
    INTEGER(i8) :: i
    INTEGER :: j, k, errcode, mpitype = MPI_INTEGER8
    CHARACTER(LEN=n) :: buf

    IF (ASSOCIATED(h%buffer)) THEN
      i = h%current_location - h%start_location + 1
      DO j = 1, nentries
        DO k = 1,n
          buf(k:k) = h%buffer(i)
          i = i + 1
        END DO
        value(j) = TRANSFER(buf, value(1))
      END DO
    ELSE
      IF (h%rank == h%rank_master) THEN
        CALL MPI_FILE_READ(h%filehandle, value, nentries, mpitype, &
            MPI_STATUS_IGNORE, errcode)
      END IF
      CALL MPI_BCAST(value, nentries, mpitype, h%rank_master, h%comm, errcode)
    END IF

    h%current_location = h%current_location + n * nentries

  END SUBROUTINE read_entry_array_int8



  SUBROUTINE read_entry_array_real4(h, value, nentries)

    INTEGER, PARAMETER :: n = sof4
    TYPE(sdf_file_handle) :: h
    REAL(r4), INTENT(OUT) :: value(:)
    INTEGER, INTENT(IN) :: nentries
    INTEGER(i8) :: i
    INTEGER :: j, k, errcode, mpitype = MPI_REAL4
    CHARACTER(LEN=n) :: buf

    IF (ASSOCIATED(h%buffer)) THEN
      i = h%current_location - h%start_location + 1
      DO j = 1, nentries
        DO k = 1,n
          buf(k:k) = h%buffer(i)
          i = i + 1
        END DO
        value(j) = TRANSFER(buf, value(1))
      END DO
    ELSE
      IF (h%rank == h%rank_master) THEN
        CALL MPI_FILE_READ(h%filehandle, value, nentries, mpitype, &
            MPI_STATUS_IGNORE, errcode)
      END IF
      CALL MPI_BCAST(value, nentries, mpitype, h%rank_master, h%comm, errcode)
    END IF

    h%current_location = h%current_location + n * nentries

  END SUBROUTINE read_entry_array_real4



  SUBROUTINE read_entry_array_real8(h, value, nentries)

    INTEGER, PARAMETER :: n = sof8
    TYPE(sdf_file_handle) :: h
    REAL(r8), INTENT(OUT) :: value(:)
    INTEGER, INTENT(IN) :: nentries
    INTEGER(i8) :: i
    INTEGER :: j, k, errcode, mpitype = MPI_REAL8
    CHARACTER(LEN=n) :: buf

    IF (ASSOCIATED(h%buffer)) THEN
      i = h%current_location - h%start_location + 1
      DO j = 1, nentries
        DO k = 1,n
          buf(k:k) = h%buffer(i)
          i = i + 1
        END DO
        value(j) = TRANSFER(buf, value(1))
      END DO
    ELSE
      IF (h%rank == h%rank_master) THEN
        CALL MPI_FILE_READ(h%filehandle, value, nentries, mpitype, &
            MPI_STATUS_IGNORE, errcode)
      END IF
      CALL MPI_BCAST(value, nentries, mpitype, h%rank_master, h%comm, errcode)
    END IF

    h%current_location = h%current_location + n * nentries

  END SUBROUTINE read_entry_array_real8



  SUBROUTINE read_entry_array_logical_c(h, value, nentries)

    INTEGER, PARAMETER :: n = 1
    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=n), INTENT(OUT) :: value(:)
    INTEGER, INTENT(IN) :: nentries
    INTEGER(i8) :: i
    INTEGER :: j, errcode, mpitype = MPI_CHARACTER

    IF (ASSOCIATED(h%buffer)) THEN
      i = h%current_location - h%start_location + 1
      DO j = 1, nentries
        value(j) = h%buffer(i)
        i = i + 1
      END DO
    ELSE
      IF (h%rank == h%rank_master) THEN
        CALL MPI_FILE_READ(h%filehandle, value, nentries, mpitype, &
            MPI_STATUS_IGNORE, errcode)
      END IF
      CALL MPI_BCAST(value, nentries, mpitype, h%rank_master, h%comm, errcode)
    END IF

    h%current_location = h%current_location + n * nentries

  END SUBROUTINE read_entry_array_logical_c



  SUBROUTINE read_entry_array_string(h, value, nentries)

    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(OUT) :: value(:)
    INTEGER, INTENT(IN) :: nentries
    INTEGER :: i

    DO i = 1, nentries
      CALL read_entry_string(h, value(i))
    END DO

  END SUBROUTINE read_entry_array_string



  SUBROUTINE sdf_safe_read_string(h, string)

    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(OUT) :: string

    CALL sdf_safe_read_string_len(h, string, INT(h%string_length))

  END SUBROUTINE sdf_safe_read_string



  SUBROUTINE sdf_safe_skip_string(h)

    TYPE(sdf_file_handle) :: h
    INTEGER(KIND=MPI_OFFSET_KIND) :: offset
    INTEGER :: errcode

    offset = h%string_length
    CALL MPI_FILE_SEEK_SHARED(h%filehandle, offset, MPI_SEEK_CUR, errcode)

  END SUBROUTINE sdf_safe_skip_string



  SUBROUTINE sdf_safe_read_id(h, string)

    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(OUT) :: string

    CALL sdf_safe_read_string_len(h, string, INT(c_id_length))

  END SUBROUTINE sdf_safe_read_id



  SUBROUTINE sdf_safe_read_string_len(h, string, length)

    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(OUT) :: string
    INTEGER, INTENT(IN) :: length
    CHARACTER(LEN=length) :: string_l
    INTEGER :: string_len, errcode, mpitype = MPI_CHARACTER

    string_len = LEN(string)

    IF (h%rank == h%rank_master) THEN
      CALL MPI_FILE_READ(h%filehandle, string_l, length, mpitype, &
          MPI_STATUS_IGNORE, errcode)
    END IF
    CALL MPI_BCAST(string_l, length, mpitype, h%rank_master, h%comm, errcode)

    string = ' '
    string = string_l(1:MIN(string_len, length))

  END SUBROUTINE sdf_safe_read_string_len



  FUNCTION sdf_check_header(h) RESULT(error)

    LOGICAL :: error
    TYPE(sdf_file_handle) :: h

    IF (.NOT. h%done_header) THEN
      IF (h%print_errors .AND. h%rank == h%rank_master) THEN
        PRINT*,'*** ERROR ***'
        PRINT*,'SDF header has not been read. Unable to read block.'
      END IF
      h%error_code = c_err_sdf
      error = .TRUE.
      RETURN
    END IF

    error = .FALSE.

  END FUNCTION sdf_check_header



  FUNCTION sdf_check_block_header(h) RESULT(error)

    LOGICAL :: error
    TYPE(sdf_file_handle) :: h

    IF (.NOT. ASSOCIATED(h%current_block)) THEN
      IF (h%print_errors .AND. h%rank == h%rank_master) THEN
        PRINT*,'*** ERROR ***'
        PRINT*,'SDF block header has not been read. Ignoring call.'
      END IF
      h%error_code = c_err_sdf
      error = .TRUE.
      RETURN
    END IF

    error = .FALSE.

  END FUNCTION sdf_check_block_header



  FUNCTION sdf_info_init(h) RESULT(error)

    LOGICAL :: error
    TYPE(sdf_file_handle) :: h
    INTEGER :: errcode
    TYPE(sdf_block_type), POINTER :: b

    error = sdf_check_block_header(h)
    IF (error) RETURN

    b => h%current_block
    h%current_location = b%block_start + h%block_header_length

    IF (b%done_info) RETURN

    IF (.NOT. ASSOCIATED(h%buffer)) THEN
      CALL MPI_FILE_SEEK(h%filehandle, h%current_location, MPI_SEEK_SET, &
          errcode)
    END IF

  END FUNCTION sdf_info_init



  SUBROUTINE read_namevalue(h, names)

    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL :: names(:)
    INTEGER :: i
    TYPE(sdf_block_type), POINTER :: b

    IF (sdf_info_init(h)) RETURN

    b => h%current_block
    IF (.NOT. b%done_info) THEN
      ALLOCATE(b%material_names(b%ndims))

      DO i = 1,b%ndims
        CALL read_entry_string(h, b%material_names(i))
      END DO

      IF (b%datatype == c_datatype_integer4) THEN
        ALLOCATE(b%i4_array(b%ndims))
        CALL read_entry_array_int4(h, b%i4_array, INT(b%ndims))

      ELSE IF (b%datatype == c_datatype_integer8) THEN
        ALLOCATE(b%i8_array(b%ndims))
        CALL read_entry_array_int8(h, b%i8_array, INT(b%ndims))

      ELSE IF (b%datatype == c_datatype_real4) THEN
        ALLOCATE(b%r4_array(b%ndims))
        CALL read_entry_array_real4(h, b%r4_array, INT(b%ndims))

      ELSE IF (b%datatype == c_datatype_real8) THEN
        ALLOCATE(b%r8_array(b%ndims))
        CALL read_entry_array_real8(h, b%r8_array, INT(b%ndims))

      ELSE IF (b%datatype == c_datatype_logical) THEN
        ALLOCATE(b%logical_array(b%ndims))
        CALL read_entry_array_logical_c(h, b%logical_array, INT(b%ndims))

      ELSE IF (b%datatype == c_datatype_character) THEN
        ALLOCATE(b%string_array(b%ndims))
        CALL read_entry_array_string(h, b%string_array, INT(b%ndims))

      END IF
    END IF

    IF (PRESENT(names)) THEN
      DO i = 1,b%ndims
        CALL sdf_safe_copy_string(b%material_names(i), names(i))
      END DO
    END IF

    h%current_location = b%data_location + b%data_length
    b%done_info = .TRUE.
    b%done_data = .TRUE.

  END SUBROUTINE read_namevalue



  SUBROUTINE read_namevalue_i4(h, names, values)

    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(OUT) :: names(:)
    INTEGER(i4), INTENT(OUT) :: values(:)
    INTEGER :: i
    TYPE(sdf_block_type), POINTER :: b

    CALL read_namevalue(h, names)

    b => h%current_block

    DO i = 1,b%ndims
      values(i) = b%i4_array(i)
    END DO

  END SUBROUTINE read_namevalue_i4



  SUBROUTINE read_namevalue_i8(h, names, values)

    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(OUT) :: names(:)
    INTEGER(i8), INTENT(OUT) :: values(:)
    INTEGER :: i
    TYPE(sdf_block_type), POINTER :: b

    CALL read_namevalue(h, names)

    b => h%current_block

    DO i = 1,b%ndims
      values(i) = b%i8_array(i)
    END DO

  END SUBROUTINE read_namevalue_i8



  SUBROUTINE read_namevalue_r4(h, names, values)

    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(OUT) :: names(:)
    REAL(r4), INTENT(OUT) :: values(:)
    INTEGER :: i
    TYPE(sdf_block_type), POINTER :: b

    CALL read_namevalue(h, names)

    b => h%current_block

    DO i = 1,b%ndims
      values(i) = b%r4_array(i)
    END DO

  END SUBROUTINE read_namevalue_r4



  SUBROUTINE read_namevalue_r8(h, names, values)

    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(OUT) :: names(:)
    REAL(r8), INTENT(OUT) :: values(:)
    INTEGER :: i
    TYPE(sdf_block_type), POINTER :: b

    CALL read_namevalue(h, names)

    b => h%current_block

    DO i = 1,b%ndims
      values(i) = b%r8_array(i)
    END DO

  END SUBROUTINE read_namevalue_r8



  SUBROUTINE read_namevalue_logical(h, names, values)

    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(OUT) :: names(:)
    LOGICAL, INTENT(OUT) :: values(:)
    INTEGER :: i
    TYPE(sdf_block_type), POINTER :: b

    CALL read_namevalue(h, names)

    b => h%current_block

    DO i = 1,b%ndims
      IF (b%logical_array(i) == ACHAR(0)) THEN
        values(i) = .FALSE.
      ELSE
        values(i) = .TRUE.
      END IF
    END DO

  END SUBROUTINE read_namevalue_logical



  SUBROUTINE read_namevalue_string(h, names, values)

    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(OUT) :: names(:)
    CHARACTER(LEN=*), INTENT(OUT) :: values(:)
    INTEGER :: i
    TYPE(sdf_block_type), POINTER :: b

    CALL read_namevalue(h, names)

    b => h%current_block

    DO i = 1,b%ndims
      CALL sdf_safe_copy_string(b%string_array(i), values(i))
    END DO

  END SUBROUTINE read_namevalue_string

END MODULE sdf_input_ru

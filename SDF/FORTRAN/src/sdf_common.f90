!
! SDF (Self-Describing Format) Fortran Library
! Copyright (c) 2010-2016, SDF Development Team
!
! Distributed under the terms of the BSD 3-clause License.
! See the LICENSE file for details.
!

MODULE sdf_common

  USE mpi
  USE sdf_job_info

  IMPLICIT NONE

  INTEGER, PARAMETER :: i4  = SELECTED_INT_KIND(9)  ! 4-byte 2^31 ~ 10^9
  INTEGER, PARAMETER :: i8  = SELECTED_INT_KIND(18) ! 8-byte 2^63 ~ 10^18

  INTEGER, PARAMETER :: r4  = SELECTED_REAL_KIND(r=30)
  INTEGER, PARAMETER :: r8  = SELECTED_REAL_KIND(r=300)
  INTEGER, PARAMETER :: r16 = SELECTED_REAL_KIND(r=3000)

  INTEGER(i8), PARAMETER :: hash_size = 2039_i8
  INTEGER, PARAMETER :: c_maxdims = 4
  INTEGER(i4), PARAMETER :: c_id_length = 32
  INTEGER(i4), PARAMETER :: c_long_id_length = 256
  INTEGER(i4), PARAMETER :: c_max_string_length = 128
  INTEGER(i8) :: npoint_per_iteration = 10000
  CHARACTER(LEN=4), PARAMETER :: c_sdf_magic = 'SDF1'
  REAL(r8), PARAMETER :: c_tiny = TINY(1.0_r8)

  LOGICAL :: print_errors   = .TRUE.
  LOGICAL :: print_warnings = .TRUE.
  LOGICAL :: exit_on_error  = .TRUE.

  TYPE sdf_run_type
    INTEGER(i4) :: version, revision, minor_rev, compile_date, run_date, io_date
    INTEGER(i8) :: defines
    CHARACTER(LEN=c_max_string_length) :: commit_id, sha1sum
    CHARACTER(LEN=c_max_string_length) :: compile_machine, compile_flags
  END TYPE sdf_run_type

  TYPE sdf_hash_list
    TYPE(sdf_block_type), POINTER :: block
    TYPE(sdf_hash_list), POINTER :: next
  END TYPE sdf_hash_list

  TYPE sdf_block_type
    REAL(r4), POINTER :: r4_array(:)
    REAL(r8), DIMENSION(2*c_maxdims) :: extents
    REAL(r8) :: mult, time, time_increment
    REAL(r8), DIMENSION(:), POINTER :: dim_mults
    REAL(r8), DIMENSION(:,:), POINTER :: station_grid
    REAL(r8), POINTER :: r8_array(:)
    INTEGER(KIND=MPI_OFFSET_KIND) :: block_start
    INTEGER(i8) :: next_block_location, data_location
    INTEGER(i8) :: nelements, npoints, data_length, info_length
    INTEGER(i8) :: id_hash
    INTEGER(i8), POINTER :: i8_array(:)
    INTEGER(i4) :: ndims, geometry, datatype, blocktype
    INTEGER(i4) :: mpitype, type_size, stagger
    INTEGER(i4) :: nstations, nvariables, step, step_increment
    INTEGER(i4) :: padding
    INTEGER(i4), DIMENSION(c_maxdims) :: dims
    INTEGER(i4), POINTER :: station_nvars(:), station_move(:), variable_types(:)
    INTEGER(i4), POINTER :: station_index(:)
    INTEGER(i4), POINTER :: i4_array(:)
    CHARACTER(LEN=1), POINTER :: logical_array(:)
    CHARACTER(LEN=8) :: const_value
    CHARACTER(LEN=c_id_length) :: id, units, mesh_id, material_id
    CHARACTER(LEN=c_id_length) :: vfm_id, obstacle_id, species_id
    CHARACTER(LEN=c_id_length) :: mimetype, checksum_type
    CHARACTER(LEN=c_long_id_length) :: long_id
    CHARACTER(LEN=c_id_length), POINTER :: station_ids(:), variable_ids(:)
    CHARACTER(LEN=c_id_length), POINTER :: dim_labels(:), dim_units(:)
    CHARACTER(LEN=c_max_string_length) :: name, material_name, checksum
    CHARACTER(LEN=c_max_string_length), POINTER :: station_names(:)
    CHARACTER(LEN=c_max_string_length), POINTER :: material_names(:)
    CHARACTER(LEN=c_max_string_length), POINTER :: string_array(:)
    LOGICAL :: done_header, done_info, done_data, truncated_id, use_mult
    TYPE(sdf_run_type), POINTER :: run
    TYPE(sdf_block_type), POINTER :: next_block
  END TYPE sdf_block_type

  TYPE sdf_file_handle
    INTEGER(KIND=MPI_OFFSET_KIND) :: current_location
    REAL(r8) :: time, time_wrote
    INTEGER(i8) :: first_block_location, summary_location, start_location
    INTEGER(i8) :: soi ! large integer to prevent overflow in calculations
    INTEGER(i8) :: data_location, summary_location_wrote
    INTEGER(i4) :: endianness, summary_size
    INTEGER(i4) :: block_header_length, string_length, nblocks, error_code
    INTEGER(i4) :: file_version, file_revision, code_io_version, step
    INTEGER(i4) :: datatype_integer, mpitype_integer
    INTEGER(i4) :: blocktype, summary_size_wrote, nblocks_wrote, step_wrote
    INTEGER(i4) :: datatype
    INTEGER :: filehandle, comm, rank, rank_master, default_rank, mode
    INTEGER :: errhandler, old_errhandler, nstations
    LOGICAL :: done_header, restart_flag, other_domains, writing, handled_error
    LOGICAL :: station_file, first, print_errors, print_warnings, exit_on_error
    LOGICAL :: station_file_wrote, writing_summary
    CHARACTER(LEN=1), POINTER :: buffer(:)
    CHARACTER(LEN=c_id_length) :: code_name
    CHARACTER(LEN=c_id_length), POINTER :: station_ids(:)
    CHARACTER(LEN=c_long_id_length) :: filename
    TYPE(jobid_type) :: jobid
    TYPE(sdf_block_type), POINTER :: blocklist, current_block
    TYPE(sdf_hash_list) :: hash_table(hash_size)
  END TYPE sdf_file_handle

  TYPE sdf_handle_type
    INTEGER :: filehandle
    TYPE(sdf_file_handle), POINTER :: handle
  END TYPE sdf_handle_type
  INTEGER, PARAMETER :: max_handles = 64
  TYPE(sdf_handle_type) :: sdf_handles(max_handles)

  INTEGER, SAVE :: errhandler_handle = MPI_ERRHANDLER_NULL
  INTEGER, SAVE :: open_handles = 0

  INTEGER, PARAMETER :: c_sdf_read = 0
  INTEGER, PARAMETER :: c_sdf_write = 1
  INTEGER, PARAMETER :: c_sdf_append = 3

  INTEGER(i4), PARAMETER :: c_blocktype_scrubbed = -1
  INTEGER(i4), PARAMETER :: c_blocktype_null = 0
  INTEGER(i4), PARAMETER :: c_blocktype_plain_mesh = 1
  INTEGER(i4), PARAMETER :: c_blocktype_point_mesh = 2
  INTEGER(i4), PARAMETER :: c_blocktype_plain_variable = 3
  INTEGER(i4), PARAMETER :: c_blocktype_point_variable = 4
  INTEGER(i4), PARAMETER :: c_blocktype_constant = 5
  INTEGER(i4), PARAMETER :: c_blocktype_array = 6
  INTEGER(i4), PARAMETER :: c_blocktype_run_info = 7
  INTEGER(i4), PARAMETER :: c_blocktype_source = 8
  INTEGER(i4), PARAMETER :: c_blocktype_stitched_tensor = 9
  INTEGER(i4), PARAMETER :: c_blocktype_stitched_material = 10
  INTEGER(i4), PARAMETER :: c_blocktype_stitched_matvar = 11
  INTEGER(i4), PARAMETER :: c_blocktype_stitched_species = 12
  INTEGER(i4), PARAMETER :: c_blocktype_species = 13
  INTEGER(i4), PARAMETER :: c_blocktype_plain_derived = 14
  INTEGER(i4), PARAMETER :: c_blocktype_point_derived = 15
  INTEGER(i4), PARAMETER :: c_blocktype_contiguous_tensor = 16
  INTEGER(i4), PARAMETER :: c_blocktype_contiguous_material = 17
  INTEGER(i4), PARAMETER :: c_blocktype_contiguous_matvar = 18
  INTEGER(i4), PARAMETER :: c_blocktype_contiguous_species = 19
  INTEGER(i4), PARAMETER :: c_blocktype_cpu_split = 20
  INTEGER(i4), PARAMETER :: c_blocktype_stitched_obstacle_group = 21
  INTEGER(i4), PARAMETER :: c_blocktype_unstructured_mesh = 22
  INTEGER(i4), PARAMETER :: c_blocktype_stitched = 23
  INTEGER(i4), PARAMETER :: c_blocktype_contiguous = 24
  INTEGER(i4), PARAMETER :: c_blocktype_lagrangian_mesh = 25
  INTEGER(i4), PARAMETER :: c_blocktype_station = 26
  INTEGER(i4), PARAMETER :: c_blocktype_station_derived = 27
  INTEGER(i4), PARAMETER :: c_blocktype_datablock = 28
  INTEGER(i4), PARAMETER :: c_blocktype_namevalue = 29
  INTEGER(i4), PARAMETER :: c_blocktype_max = 29

  INTEGER(i4), PARAMETER :: c_datatype_null = 0
  INTEGER(i4), PARAMETER :: c_datatype_integer4 = 1
  INTEGER(i4), PARAMETER :: c_datatype_integer8 = 2
  INTEGER(i4), PARAMETER :: c_datatype_real4 = 3
  INTEGER(i4), PARAMETER :: c_datatype_real8 = 4
  INTEGER(i4), PARAMETER :: c_datatype_real16 = 5
  INTEGER(i4), PARAMETER :: c_datatype_character = 6
  INTEGER(i4), PARAMETER :: c_datatype_logical = 7
  INTEGER(i4), PARAMETER :: c_datatype_other = 8
  INTEGER(i4), PARAMETER :: c_datatype_max = 8
  INTEGER, PARAMETER :: c_type_sizes(0:c_datatype_max) = &
      (/ 0, 4, 8, 4, 8, 16, 1, 1, 0 /)

  INTEGER(i4), PARAMETER :: c_geometry_null = 0
  INTEGER(i4), PARAMETER :: c_geometry_cartesian = 1
  INTEGER(i4), PARAMETER :: c_geometry_cylindrical = 2
  INTEGER(i4), PARAMETER :: c_geometry_spherical = 3

  ! c_dimension_irrelevant is used where the dimensionality isn't needed, as
  ! with point variables still keep dimensionality as a common quantity
  ! because other than this, they really are very alike
  INTEGER(i4), PARAMETER :: c_dimension_irrelevant = 0
  INTEGER(i4), PARAMETER :: c_dimension_1d = 1
  INTEGER(i4), PARAMETER :: c_dimension_2d = 2
  INTEGER(i4), PARAMETER :: c_dimension_3d = 3

  INTEGER(i4), PARAMETER :: c_stagger_cell_centre = 0
  INTEGER(i4), PARAMETER :: c_stagger_face_x = 1
  INTEGER(i4), PARAMETER :: c_stagger_face_y = 2
  INTEGER(i4), PARAMETER :: c_stagger_face_z = 4
  INTEGER(i4), PARAMETER :: c_stagger_edge_x = &
      c_stagger_face_y + c_stagger_face_z
  INTEGER(i4), PARAMETER :: c_stagger_edge_y = &
      c_stagger_face_x + c_stagger_face_z
  INTEGER(i4), PARAMETER :: c_stagger_edge_z = &
      c_stagger_face_x + c_stagger_face_y
  INTEGER(i4), PARAMETER :: c_stagger_vertex = &
      c_stagger_face_x + c_stagger_face_y + c_stagger_face_z

  CHARACTER(LEN=*), PARAMETER :: c_checksum_null = ''
  CHARACTER(LEN=*), PARAMETER :: c_checksum_md5 = 'md5'
  CHARACTER(LEN=*), PARAMETER :: c_checksum_sha1 = 'sha1'
  CHARACTER(LEN=*), PARAMETER :: c_checksum_sha256 = 'sha256'

  INTEGER(i4), PARAMETER :: sdf_version = 1, sdf_revision = 4
  INTEGER(i4), PARAMETER :: sdf_minor_rev = 0

  INTEGER(i4), PARAMETER :: soi4 = 4 ! Size of 4-byte integer
  INTEGER(i4), PARAMETER :: soi8 = 8 ! Size of 8-byte integer
  INTEGER(i4), PARAMETER :: sof4 = 4 ! Size of 4-byte real
  INTEGER(i4), PARAMETER :: sof8 = 8 ! Size of 8-byte real
  INTEGER(i4), PARAMETER :: sol  = 4 ! Size of logical

  ! header length (including padding) - must be updated if sdf_write_header
  ! changes
  INTEGER, PARAMETER :: c_header_length = 11 * soi4 + 2 * soi8 + sof8 + 12 &
      + c_id_length

  ! summary offset - must be updated if sdf_write_header changes
  INTEGER(i4), PARAMETER :: c_summary_offset = 4 + 3 * soi4 + c_id_length + soi8

  INTEGER(i4), PARAMETER :: c_endianness = 16911887

  INTEGER(KIND=MPI_OFFSET_KIND), PARAMETER :: c_off0 = 0
  INTEGER, PARAMETER :: max_mpi_error_codes = 21
  INTEGER, PARAMETER :: mpi_error_codes(max_mpi_error_codes) = (/ &
      MPI_ERR_ACCESS, MPI_ERR_AMODE, MPI_ERR_BAD_FILE, MPI_ERR_CONVERSION, &
      MPI_ERR_DUP_DATAREP, MPI_ERR_FILE, MPI_ERR_FILE_EXISTS, &
      MPI_ERR_FILE_IN_USE, MPI_ERR_INFO, MPI_ERR_INFO_KEY, MPI_ERR_INFO_NOKEY, &
      MPI_ERR_INFO_VALUE, MPI_ERR_IO, MPI_ERR_NOT_SAME, MPI_ERR_NO_SPACE, &
      MPI_ERR_NO_SUCH_FILE, MPI_ERR_QUOTA, MPI_ERR_READ_ONLY, &
      MPI_ERR_UNSUPPORTED_DATAREP, MPI_ERR_UNSUPPORTED_OPERATION, &
      MPI_ERR_TYPE /)

  ! SDF errors
  INTEGER, PARAMETER :: c_err_success = 0
  INTEGER, PARAMETER :: c_err_unknown = 1
  INTEGER, PARAMETER :: c_err_unsupported_file = 2
  INTEGER, PARAMETER :: c_err_sdf = 3

  ! MPI errors
  INTEGER, PARAMETER :: c_mpi_error_start = 31
  INTEGER, PARAMETER :: c_err_access = 31
  INTEGER, PARAMETER :: c_err_amode = 32
  INTEGER, PARAMETER :: c_err_bad_file = 33
  INTEGER, PARAMETER :: c_err_conversion = 34
  INTEGER, PARAMETER :: c_err_dup_datarep = 35
  INTEGER, PARAMETER :: c_err_file = 36
  INTEGER, PARAMETER :: c_err_file_exists = 37
  INTEGER, PARAMETER :: c_err_file_in_use = 38
  INTEGER, PARAMETER :: c_err_info = 39
  INTEGER, PARAMETER :: c_err_info_key = 40
  INTEGER, PARAMETER :: c_err_info_nokey = 41
  INTEGER, PARAMETER :: c_err_info_value = 42
  INTEGER, PARAMETER :: c_err_io = 43
  INTEGER, PARAMETER :: c_err_not_same = 44
  INTEGER, PARAMETER :: c_err_no_space = 45
  INTEGER, PARAMETER :: c_err_no_such_file = 46
  INTEGER, PARAMETER :: c_err_quota = 47
  INTEGER, PARAMETER :: c_err_read_only = 48
  INTEGER, PARAMETER :: c_err_unsupported_datarep = 49
  INTEGER, PARAMETER :: c_err_unsupported_operation = 50
  INTEGER, PARAMETER :: c_err_type = 51
  INTEGER, PARAMETER :: c_err_max = 51

  CHARACTER(LEN=*), PARAMETER :: c_blocktypes_char(-1:c_blocktype_max) = (/ &
      'SDF_BLOCKTYPE_SCRUBBED               ', &
      'SDF_BLOCKTYPE_NULL                   ', &
      'SDF_BLOCKTYPE_PLAIN_MESH             ', &
      'SDF_BLOCKTYPE_POINT_MESH             ', &
      'SDF_BLOCKTYPE_PLAIN_VARIABLE         ', &
      'SDF_BLOCKTYPE_POINT_VARIABLE         ', &
      'SDF_BLOCKTYPE_CONSTANT               ', &
      'SDF_BLOCKTYPE_ARRAY                  ', &
      'SDF_BLOCKTYPE_RUN_INFO               ', &
      'SDF_BLOCKTYPE_SOURCE                 ', &
      'SDF_BLOCKTYPE_STITCHED_TENSOR        ', &
      'SDF_BLOCKTYPE_STITCHED_MATERIAL      ', &
      'SDF_BLOCKTYPE_STITCHED_MATVAR        ', &
      'SDF_BLOCKTYPE_STITCHED_SPECIES       ', &
      'SDF_BLOCKTYPE_SPECIES                ', &
      'SDF_BLOCKTYPE_PLAIN_DERIVED          ', &
      'SDF_BLOCKTYPE_POINT_DERIVED          ', &
      'SDF_BLOCKTYPE_MULTI_TENSOR           ', &
      'SDF_BLOCKTYPE_MULTI_MATERIAL         ', &
      'SDF_BLOCKTYPE_MULTI_MATVAR           ', &
      'SDF_BLOCKTYPE_MULTI_SPECIES          ', &
      'SDF_BLOCKTYPE_CPU_SPLIT              ', &
      'SDF_BLOCKTYPE_STITCHED_OBSTACLE_GROUP', &
      'SDF_BLOCKTYPE_UNSTRUCTURED_MESH      ', &
      'SDF_BLOCKTYPE_STITCHED               ', &
      'SDF_BLOCKTYPE_CONTIGUOUS             ', &
      'SDF_BLOCKTYPE_LAGRANGIAN_MESH        ', &
      'SDF_BLOCKTYPE_STATION                ', &
      'SDF_BLOCKTYPE_STATION_DERIVED        ', &
      'SDF_BLOCKTYPE_DATABLOCK              ', &
      'SDF_BLOCKTYPE_NAMEVALUE              ' /)

  CHARACTER(LEN=*), PARAMETER :: c_datatypes_char(0:c_datatype_max) = (/ &
      'SDF_DATATYPE_NULL     ', &
      'SDF_DATATYPE_INTEGER4 ', &
      'SDF_DATATYPE_INTEGER8 ', &
      'SDF_DATATYPE_REAL4    ', &
      'SDF_DATATYPE_REAL8    ', &
      'SDF_DATATYPE_REAL16   ', &
      'SDF_DATATYPE_CHARACTER', &
      'SDF_DATATYPE_LOGICAL  ', &
      'SDF_DATATYPE_OTHER    ' /)

  CHARACTER(LEN=*), PARAMETER :: c_errcodes_char(0:c_err_max) = (/ &
      'SDF_ERR_SUCCESS              ', 'SDF_ERR_UNKNOWN              ', &
      'SDF_ERR_UNSUPPORTED_FILE     ', 'SDF_ERR_SDF                  ', &
      '                             ', '                             ', &
      '                             ', '                             ', &
      '                             ', '                             ', &
      '                             ', '                             ', &
      '                             ', '                             ', &
      '                             ', '                             ', &
      '                             ', '                             ', &
      '                             ', '                             ', &
      '                             ', '                             ', &
      '                             ', '                             ', &
      '                             ', '                             ', &
      '                             ', '                             ', &
      '                             ', '                             ', &
      '                             ', 'MPI_ERR_ACCESS               ', &
      'MPI_ERR_AMODE                ', 'MPI_ERR_BAD_FILE             ', &
      'MPI_ERR_CONVERSION           ', 'MPI_ERR_DUP_DATAREP          ', &
      'MPI_ERR_FILE                 ', 'MPI_ERR_FILE_EXISTS          ', &
      'MPI_ERR_FILE_IN_USE          ', 'MPI_ERR_INFO                 ', &
      'MPI_ERR_INFO_KEY             ', 'MPI_ERR_INFO_NOKEY           ', &
      'MPI_ERR_INFO_VALUE           ', 'MPI_ERR_IO                   ', &
      'MPI_ERR_NOT_SAME             ', 'MPI_ERR_NO_SPACE             ', &
      'MPI_ERR_NO_SUCH_FILE         ', 'MPI_ERR_QUOTA                ', &
      'MPI_ERR_READ_ONLY            ', 'MPI_ERR_UNSUPPORTED_DATAREP  ', &
      'MPI_ERR_UNSUPPORTED_OPERATION', 'MPI_ERR_TYPE                 ' /)

  INTERFACE sdf_set_point_array_size
    MODULE PROCEDURE &
        set_point_array_size_i8, &
        set_point_array_size_i4
  END INTERFACE sdf_set_point_array_size

CONTAINS

  SUBROUTINE sdf_get_next_block(h)

    TYPE(sdf_file_handle) :: h
    TYPE(sdf_block_type), POINTER :: next

    IF (ASSOCIATED(h%blocklist)) THEN
      IF (.NOT. ASSOCIATED(h%current_block)) THEN
        h%current_block => h%blocklist
        RETURN
      ELSE IF (h%first .AND. ASSOCIATED(h%current_block, h%blocklist)) THEN
        h%first = .FALSE.
        RETURN
      ELSE IF (ASSOCIATED(h%current_block%next_block)) THEN
        h%current_block => h%current_block%next_block
        RETURN
      ELSE
        ALLOCATE(h%current_block%next_block)
        CALL initialise_block_type(h%current_block%next_block)
        next => h%current_block%next_block
        next%block_start = h%current_block%next_block_location
      END IF
    ELSE
      ALLOCATE(h%blocklist)
      CALL initialise_block_type(h%blocklist)
      next => h%blocklist
      next%block_start = h%summary_location
    END IF

    h%first = .FALSE.
    next%done_header = .FALSE.
    next%done_info = .FALSE.
    next%done_data = .FALSE.
    NULLIFY(next%run)
    NULLIFY(next%next_block)
    h%current_block => next

  END SUBROUTINE sdf_get_next_block



  SUBROUTINE sdf_seek_start(h)

    TYPE(sdf_file_handle) :: h

    h%current_block => h%blocklist
    h%first = .TRUE.

  END SUBROUTINE sdf_seek_start



  FUNCTION sdf_find_block(h, b, block_id) RESULT(found)

    TYPE(sdf_file_handle) :: h
    TYPE(sdf_block_type), POINTER :: b
    CHARACTER(LEN=*), INTENT(IN) :: block_id
    INTEGER :: i
    INTEGER(i8) :: id_hash
    LOGICAL :: found, use_truncated
    TYPE(sdf_hash_list), POINTER :: hash_item

    id_hash = sdf_hash_function(TRIM(block_id))
    i = INT(MOD(ABS(id_hash), hash_size)) + 1

    b => h%hash_table(i)%block
    IF (.NOT.ASSOCIATED(b)) THEN
      found = .FALSE.
      NULLIFY(b)
      RETURN
    END IF

    use_truncated = (LEN_TRIM(block_id) > c_id_length)
    found = .TRUE.

    IF (b%id_hash == id_hash) THEN
      IF (use_truncated .AND. b%truncated_id) THEN
        IF (sdf_string_equal(block_id, b%long_id)) RETURN
      ELSE
        IF (sdf_string_equal(block_id, b%id)) RETURN
      END IF
    END IF

    hash_item => h%hash_table(i)%next
    DO WHILE (ASSOCIATED(hash_item))
      b => hash_item%block
      IF (b%id_hash == id_hash) THEN
        IF (use_truncated .AND. b%truncated_id) THEN
          IF (sdf_string_equal(block_id, b%long_id)) RETURN
        ELSE
          IF (sdf_string_equal(block_id, b%id)) RETURN
        END IF
      END IF
      hash_item => hash_item%next
    END DO

    found = .FALSE.
    NULLIFY(b)

  END FUNCTION sdf_find_block



  FUNCTION sdf_find_block_by_id(h, block_id) RESULT(found)

    TYPE(sdf_file_handle) :: h
    TYPE(sdf_block_type), POINTER :: b
    CHARACTER(LEN=*), INTENT(IN) :: block_id
    LOGICAL :: found

    found = sdf_find_block(h, b, block_id)
    IF (found) h%current_block => b

  END FUNCTION sdf_find_block_by_id



  FUNCTION sdf_get_block_id(h, long_id, block_id) RESULT(found)

    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: long_id
    CHARACTER(LEN=*), INTENT(OUT) :: block_id
    TYPE(sdf_block_type), POINTER :: b
    LOGICAL :: found

    found = sdf_find_block(h, b, long_id)
    IF (found) THEN
      CALL sdf_safe_copy_string(b%id, block_id)
    ELSE
      CALL sdf_safe_copy_string(long_id, block_id)
    END IF

  END FUNCTION sdf_get_block_id



  FUNCTION sdf_seek_block(h, block_id) RESULT(found)

    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: block_id
    TYPE(sdf_block_type), POINTER :: b
    LOGICAL :: found

    found = sdf_find_block(h, b, block_id)
    IF (found) h%current_block => b

  END FUNCTION sdf_seek_block



  FUNCTION sdf_get_data_location(h) RESULT(data_location)

    TYPE(sdf_file_handle) :: h
    INTEGER(i8) :: data_location

    data_location = h%current_block%data_location

  END FUNCTION sdf_get_data_location



  SUBROUTINE sdf_set_data_location(h, data_location)

    TYPE(sdf_file_handle) :: h
    INTEGER(i8), INTENT(IN) :: data_location

    h%data_location = data_location

  END SUBROUTINE sdf_set_data_location



  FUNCTION sdf_string_equal(str1, str2) RESULT(equal)

    CHARACTER(LEN=*), INTENT(IN) :: str1, str2
    INTEGER :: len1, len2
    LOGICAL :: equal

    len1 = LEN_TRIM(str1)
    len2 = LEN_TRIM(str2)

    IF (len1 > 0) THEN
      IF (IACHAR(str1(len1:len1)) == 0) len1 = len1 - 1
    END IF
    IF (len2 > 0) THEN
      IF (IACHAR(str2(len2:len2)) == 0) len2 = len2 - 1
    END IF

    IF (len1 /= len2) THEN
      equal = .FALSE.
      RETURN
    END IF

    equal = (str1(1:len1) == str2(1:len1))

  END FUNCTION sdf_string_equal



  SUBROUTINE sdf_safe_copy_string(s1, s2)

    CHARACTER(LEN=*), INTENT(IN) :: s1
    CHARACTER(LEN=*), INTENT(OUT) :: s2
    INTEGER :: len1, len2, olen, i

    len1 = LEN_TRIM(s1)
    len2 = LEN(s2)
    olen = MIN(len1,len2)
    IF (olen > 0) THEN
      s2(1:olen) = s1(1:olen)
      DO i = olen+1,len2
        s2(i:i) = ' '
      END DO
    ELSE
      DO i = 1,len2
        s2(i:i) = ' '
      END DO
    END IF

  END SUBROUTINE sdf_safe_copy_string



  SUBROUTINE sdf_safe_copy_id(h, id, new_id)

    TYPE(sdf_file_handle) :: h
    CHARACTER(LEN=*), INTENT(IN) :: id
    CHARACTER(LEN=c_id_length), INTENT(OUT) :: new_id

    IF (LEN_TRIM(id) > c_id_length) THEN
      IF (h%print_warnings .AND. h%rank == h%rank_master) THEN
        PRINT*, '*** WARNING ***'
        PRINT*, 'SDF ID string "' // TRIM(id) // '" was truncated.'
      END IF
    END IF

    CALL sdf_safe_copy_string(id, new_id)

  END SUBROUTINE sdf_safe_copy_id



  SUBROUTINE sdf_safe_copy_unique_id(h, b, id)

    TYPE(sdf_file_handle) :: h
    TYPE(sdf_block_type), POINTER :: b
    CHARACTER(LEN=*), INTENT(IN) :: id
    LOGICAL :: found
    CHARACTER(LEN=*), PARAMETER :: numbers = '0123456789'
    INTEGER :: i, n, num, pos, old_len
    TYPE(sdf_block_type), POINTER :: tmp

    IF (LEN_TRIM(id) > c_id_length) THEN
      b%truncated_id = .TRUE.
      CALL sdf_safe_copy_string(id, b%long_id)
      IF (LEN_TRIM(id) > c_long_id_length) THEN
        IF (h%print_warnings .AND. h%rank == h%rank_master) THEN
          PRINT*, '*** WARNING ***'
          PRINT*, 'SDF ID string "' // TRIM(id) // '" was truncated.'
        END IF
      END IF
    END IF

    CALL sdf_safe_copy_string(id, b%id)
    old_len = LEN_TRIM(b%id)
    found = sdf_find_block(h, tmp, b%id)
    i = 1
    DO WHILE(found)
      num = i

      ! Count digits
      pos = 1
      n = MOD(num,10)
      num = num / 10
      DO WHILE(num > 0)
        pos = pos + 1
        n = MOD(num,10)
        num = num / 10
      END DO

      num = i

      ! Generate a new ID by adding ASCII digits to the end.
      pos = MIN(c_id_length, old_len+pos)
      n = MOD(num,10)
      b%id(pos:pos) = numbers(n+1:n+1)
      num = num / 10
      DO WHILE(num > 0)
        pos = pos - 1
        n = MOD(num,10)
        b%id(pos:pos) = numbers(n+1:n+1)
        num = num / 10
      END DO

      found = sdf_find_block(h, tmp, b%id)
      i = i + 1
    END DO

  END SUBROUTINE sdf_safe_copy_unique_id



  SUBROUTINE initialise_block_type(var)

    TYPE(sdf_block_type) :: var

    NULLIFY(var%dim_mults)
    NULLIFY(var%variable_ids)
    NULLIFY(var%dim_labels)
    NULLIFY(var%dim_units)
    NULLIFY(var%material_names)
    NULLIFY(var%run)
    NULLIFY(var%next_block)
    NULLIFY(var%station_ids)
    NULLIFY(var%station_names)
    NULLIFY(var%station_nvars)
    NULLIFY(var%station_move)
    NULLIFY(var%station_grid)
    NULLIFY(var%station_index)
    NULLIFY(var%variable_types)
    NULLIFY(var%i4_array)
    NULLIFY(var%i8_array)
    NULLIFY(var%r4_array)
    NULLIFY(var%r8_array)
    NULLIFY(var%logical_array)
    NULLIFY(var%string_array)
    var%done_header = .FALSE.
    var%done_info = .FALSE.
    var%done_data = .FALSE.
    var%truncated_id = .FALSE.
    var%use_mult = .FALSE.
    var%data_location = 0
    var%blocktype = c_blocktype_null
    var%datatype = 0
    var%step = 0
    var%step_increment = 0
    var%time = 0
    var%time_increment = 0
    var%padding = 0
    var%id = ''
    var%name = ''

  END SUBROUTINE initialise_block_type



  SUBROUTINE deallocate_block_type(var)

    TYPE(sdf_block_type) :: var

    IF (ASSOCIATED(var%dim_mults))      DEALLOCATE(var%dim_mults)
    IF (ASSOCIATED(var%variable_ids))   DEALLOCATE(var%variable_ids)
    IF (ASSOCIATED(var%dim_labels))     DEALLOCATE(var%dim_labels)
    IF (ASSOCIATED(var%dim_units))      DEALLOCATE(var%dim_units)
    IF (ASSOCIATED(var%material_names)) DEALLOCATE(var%material_names)
    IF (ASSOCIATED(var%run))            DEALLOCATE(var%run)
    IF (ASSOCIATED(var%station_ids))    DEALLOCATE(var%station_ids)
    IF (ASSOCIATED(var%station_names))  DEALLOCATE(var%station_names)
    IF (ASSOCIATED(var%station_nvars))  DEALLOCATE(var%station_nvars)
    IF (ASSOCIATED(var%station_move))   DEALLOCATE(var%station_move)
    IF (ASSOCIATED(var%station_grid))   DEALLOCATE(var%station_grid)
    IF (ASSOCIATED(var%station_index))  DEALLOCATE(var%station_index)
    IF (ASSOCIATED(var%variable_types)) DEALLOCATE(var%variable_types)
    IF (ASSOCIATED(var%i4_array))       DEALLOCATE(var%i4_array)
    IF (ASSOCIATED(var%i8_array))       DEALLOCATE(var%i8_array)
    IF (ASSOCIATED(var%r4_array))       DEALLOCATE(var%r4_array)
    IF (ASSOCIATED(var%r8_array))       DEALLOCATE(var%r8_array)
    IF (ASSOCIATED(var%logical_array))  DEALLOCATE(var%logical_array)
    IF (ASSOCIATED(var%string_array))   DEALLOCATE(var%string_array)

    CALL initialise_block_type(var)

  END SUBROUTINE deallocate_block_type



  SUBROUTINE initialise_file_handle(var, set_handler)

    TYPE(sdf_file_handle) :: var
    LOGICAL, INTENT(IN), OPTIONAL :: set_handler
    LOGICAL :: set_err_handler
    INTEGER :: ierr, i

    NULLIFY(var%buffer)
    NULLIFY(var%blocklist)
    NULLIFY(var%current_block)
    NULLIFY(var%station_ids)
    ! Set filehandle to -1 to show that the file is closed
    var%filehandle = -1
    var%string_length = c_max_string_length
    var%done_header = .FALSE.
    var%restart_flag = .FALSE.
    var%other_domains = .FALSE.
    var%writing = .FALSE.
    var%handled_error = .FALSE.
    var%station_file = .FALSE.
    var%first = .TRUE.
    var%writing_summary = .FALSE.
    var%print_errors = print_errors
    var%print_warnings = print_warnings
    var%exit_on_error = exit_on_error
    var%nblocks = 0
    var%error_code = 0
    var%errhandler = MPI_ERRHANDLER_NULL
    var%old_errhandler = MPI_ERRHANDLER_NULL
    var%comm = 0

    var%summary_location = 0
    var%summary_size = 0
    var%step = 0
    var%time = 0
    DO i = 1, hash_size
      NULLIFY(var%hash_table(i)%block)
      NULLIFY(var%hash_table(i)%next)
    END DO

    var%summary_location_wrote = var%summary_location
    var%summary_size_wrote = var%summary_size
    var%nblocks_wrote = var%nblocks
    var%step_wrote = var%step
    var%time_wrote = var%time
    var%station_file_wrote = var%station_file

    IF (PRESENT(set_handler)) THEN
      set_err_handler = set_handler
    ELSE
      set_err_handler = .TRUE.
    END IF

    IF (set_err_handler) THEN
      CALL MPI_FILE_GET_ERRHANDLER(MPI_FILE_NULL, var%old_errhandler, ierr)
      IF (errhandler_handle == MPI_ERRHANDLER_NULL) THEN
        CALL MPI_FILE_CREATE_ERRHANDLER(error_handler, errhandler_handle, ierr)
      END IF
      var%errhandler = errhandler_handle
      CALL MPI_FILE_SET_ERRHANDLER(MPI_FILE_NULL, var%errhandler, ierr)
    END IF

  END SUBROUTINE initialise_file_handle



  SUBROUTINE deallocate_file_handle(var)

    TYPE(sdf_file_handle) :: var
    INTEGER :: errcode, i
    TYPE(sdf_hash_list), POINTER :: hash_item, hash_item_next

    IF (ASSOCIATED(var%buffer)) DEALLOCATE(var%buffer)
    IF (ASSOCIATED(var%station_ids)) DEALLOCATE(var%station_ids)

    DO i = 1, hash_size
      hash_item => var%hash_table(i)%next
      DO WHILE (ASSOCIATED(hash_item))
        hash_item_next => hash_item%next
        DEALLOCATE(hash_item)
        hash_item => hash_item_next
      END DO
      NULLIFY(var%hash_table(i)%block)
      NULLIFY(var%hash_table(i)%next)
    END DO

    var%errhandler = MPI_ERRHANDLER_NULL

    IF (var%old_errhandler /= MPI_ERRHANDLER_NULL) THEN
      CALL MPI_FILE_SET_ERRHANDLER(MPI_FILE_NULL, var%old_errhandler, errcode)
      var%old_errhandler = MPI_ERRHANDLER_NULL
    END IF

    IF (var%comm /= 0) CALL MPI_COMM_FREE(var%comm, errcode)

    DO i = 1, max_handles
      IF (sdf_handles(i)%filehandle == var%filehandle) THEN
        sdf_handles(i)%filehandle = 0
        open_handles = open_handles - 1
        EXIT
      END IF
    END DO

    CALL initialise_file_handle(var, set_handler=.FALSE.)

    IF (open_handles == 0 .AND. errhandler_handle /= MPI_ERRHANDLER_NULL) THEN
      CALL MPI_ERRHANDLER_FREE(errhandler_handle, errcode)
      errhandler_handle = MPI_ERRHANDLER_NULL
    END IF

  END SUBROUTINE deallocate_file_handle



  FUNCTION map_error_code(error_code) RESULT(errcode)

    INTEGER, INTENT(IN) :: error_code
    INTEGER :: errcode, i

    errcode = c_err_unknown
    DO i = 1, max_mpi_error_codes
      IF (error_code == mpi_error_codes(i)) THEN
        errcode = i + c_mpi_error_start - 1
        RETURN
      END IF
    END DO

  END FUNCTION map_error_code



  SUBROUTINE error_handler(filehandle, error_code)

    INTEGER :: filehandle, error_code
    TYPE(sdf_file_handle), POINTER :: h
    CHARACTER(LEN=MPI_MAX_ERROR_STRING) :: message
    CHARACTER(LEN=MPI_MAX_INFO_KEY) :: key
    CHARACTER(LEN=MPI_MAX_INFO_VAL) :: info_value
    INTEGER(KIND=MPI_OFFSET_KIND) :: filepos
    INTEGER :: sdf_error, message_len, info, nkeys, i, ierr
    LOGICAL :: found, print_error, do_abort
    REAL :: zz

    found = .FALSE.
    IF (filehandle > 0) THEN
      DO i = 1, max_handles
        IF (sdf_handles(i)%filehandle == filehandle) THEN
          h => sdf_handles(i)%handle
          found = .TRUE.
          EXIT
        END IF
      END DO
    END IF

    sdf_error = map_error_code(error_code)

    print_error = print_errors
    do_abort = exit_on_error

    IF (found) THEN
      print_error = .FALSE.
      IF (filehandle > 0) THEN
        do_abort = h%exit_on_error
        IF (.NOT.h%handled_error) THEN
          h%error_code = sdf_error + 64 * h%nblocks
          h%handled_error = .TRUE.
          print_error = h%print_errors
        END IF
      END IF
    END IF

    IF (print_error) THEN
      CALL MPI_ERROR_STRING(error_code, message, message_len, ierr)

      WRITE(0,*) 'An MPI-I/O error has occurred'
      IF (found) THEN
        WRITE(0,*) 'Process:     ', h%rank
        WRITE(0,*) 'Filename:    ' // TRIM(h%filename)
      END IF
      WRITE(0,*) 'File handle: ', filehandle
      WRITE(0,*) 'Error code:  ', error_code
      WRITE(0,*) 'SDF error:   ' // TRIM(c_errcodes_char(sdf_error))
      WRITE(0,*) 'Message:     ' // TRIM(message)
      IF (filehandle > 0) THEN
        CALL MPI_FILE_GET_POSITION(filehandle, filepos, ierr)
        WRITE(0,*) 'Position:    ', filepos
      END IF

      IF (filehandle > 0) THEN
        CALL MPI_FILE_GET_INFO(filehandle, info, ierr)
        IF (ierr == 0) THEN
          CALL MPI_INFO_GET_NKEYS(info, nkeys, ierr)
          IF (nkeys > 0) THEN
            WRITE(0,*) 'Info:'
            DO i = 0,nkeys-1
              CALL MPI_INFO_GET_NTHKEY(info, i, key, ierr)
              CALL MPI_INFO_GET(info, key, MPI_MAX_INFO_VAL, info_value, &
                                found, ierr)
              WRITE(0,'(10X,A,": ",A)') TRIM(key), TRIM(info_value)
            END DO
          END IF
          CALL MPI_INFO_FREE(info, ierr)
        END IF
      END IF
    END IF

    IF (do_abort) THEN
      ! First try to generate a floating-point error.
      ! This sometimes allows us to get a backtrace.
      zz = 0.0
      zz = 1.0 / zz
      CALL MPI_ABORT(MPI_COMM_WORLD, 10, ierr)
      STOP
    END IF

  END SUBROUTINE error_handler



  SUBROUTINE set_point_array_size_i8(value)

    INTEGER(i8) :: value

    npoint_per_iteration = value

  END SUBROUTINE set_point_array_size_i8



  SUBROUTINE set_point_array_size_i4(value)

    INTEGER(i4) :: value

    npoint_per_iteration = INT(value,i8)

  END SUBROUTINE set_point_array_size_i4



  FUNCTION sdf_get_point_array_size() RESULT(value)

    INTEGER(i8) :: value

    value = npoint_per_iteration

  END FUNCTION sdf_get_point_array_size



  SUBROUTINE sdf_set_exit_on_error(exit_on_error_val)

    LOGICAL, INTENT(IN) :: exit_on_error_val

    exit_on_error = exit_on_error_val

  END SUBROUTINE sdf_set_exit_on_error



  FUNCTION sdf_get_exit_on_error()

    LOGICAL :: sdf_get_exit_on_error

    sdf_get_exit_on_error = exit_on_error

  END FUNCTION sdf_get_exit_on_error



  FUNCTION sdf_hash_function(str) RESULT(hash)

    CHARACTER(LEN=*), INTENT(IN) :: str
    INTEGER(i8) :: hash
    INTEGER :: i

    hash = 5381

    DO i = 1, LEN(str)
      hash = (ISHFT(hash,5) + hash) + ICHAR(str(i:i))
    END DO

  END FUNCTION sdf_hash_function



  SUBROUTINE add_to_hash_table(h, b)

    TYPE(sdf_file_handle), INTENT(INOUT) :: h
    TYPE(sdf_block_type), POINTER , INTENT(INOUT):: b
    TYPE(sdf_hash_list), POINTER :: hash_item
    INTEGER :: i

    IF (h%writing_summary) RETURN

    b%id_hash = sdf_hash_function(TRIM(b%id))
    i = INT(MOD(ABS(b%id_hash), hash_size)) + 1

    IF (ASSOCIATED(h%hash_table(i)%block)) THEN
      IF (ASSOCIATED(h%hash_table(i)%next)) THEN
        hash_item => h%hash_table(i)%next
        DO WHILE (ASSOCIATED(hash_item%next))
          hash_item => hash_item%next
        END DO
        ALLOCATE(hash_item%next)
        hash_item => hash_item%next
      ELSE
        ALLOCATE(h%hash_table(i)%next)
        hash_item => h%hash_table(i)%next
      END IF
      hash_item%block => b
      NULLIFY(hash_item%next)
    ELSE
      h%hash_table(i)%block => b
      NULLIFY(h%hash_table(i)%next)
    END IF

  END SUBROUTINE add_to_hash_table

END MODULE sdf_common

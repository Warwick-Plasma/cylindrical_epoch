!
! SDF (Self-Describing Format) Fortran Library
! Copyright (c) 2010-2016, SDF Development Team
!
! Distributed under the terms of the BSD 3-clause License.
! See the LICENSE file for details.
!

MODULE sdf_output_util

  USE sdf_output_cartesian_ru
  USE sdf_output_point_ru

  IMPLICIT NONE

CONTAINS

  SUBROUTINE sdf_write_summary(h)

    TYPE(sdf_file_handle) :: h
    INTEGER :: i, errcode
    TYPE(sdf_block_type), POINTER :: b

    IF (h%rank /= h%rank_master) RETURN
    IF (.NOT.ASSOCIATED(h%current_block)) RETURN

    b => h%current_block

    h%writing_summary = .TRUE.
    h%summary_location = b%next_block_location
    h%current_location = h%summary_location
    CALL MPI_FILE_SEEK(h%filehandle, h%current_location, MPI_SEEK_SET, &
        errcode)

    b => h%blocklist
    b%block_start = h%current_location
    b%next_block_location = b%block_start + b%info_length
    b%done_header = .FALSE.
    h%current_block => b

    CALL sdf_write_block_info(h)

    DO i = 2,h%nblocks
      h%current_location = b%next_block_location
      b => b%next_block
      b%block_start = h%current_location
      b%next_block_location = b%block_start + b%info_length
      b%done_header = .FALSE.
      h%current_block => b

      CALL sdf_write_block_info(h)
    END DO

    h%summary_size = INT(h%current_location - h%summary_location,i4)
    h%writing_summary = .FALSE.

  END SUBROUTINE sdf_write_summary



  SUBROUTINE sdf_write_block_info(h)

    TYPE(sdf_file_handle) :: h
    TYPE(sdf_block_type), POINTER :: b

    b => h%current_block

    IF (b%blocktype == c_blocktype_plain_mesh &
        .OR. b%blocktype == c_blocktype_lagrangian_mesh) THEN
      CALL write_mesh_meta_r8(h)
    ELSE IF (b%blocktype == c_blocktype_point_mesh) THEN
      CALL write_point_mesh_meta_r8(h)
    ELSE IF (b%blocktype == c_blocktype_plain_variable) THEN
      CALL write_mesh_variable_meta_r8(h)
    ELSE IF (b%blocktype == c_blocktype_point_variable) THEN
      CALL write_point_variable_meta_r8(h)
    ELSE IF (b%blocktype == c_blocktype_constant) THEN
      CALL write_constant_meta(h)
    ELSE IF (b%blocktype == c_blocktype_array) THEN
      CALL write_array_meta(h)
    ELSE IF (b%blocktype == c_blocktype_cpu_split) THEN
      CALL write_cpu_split_meta(h)
    ELSE IF (b%blocktype == c_blocktype_run_info) THEN
      CALL write_run_info_meta(h)
    ELSE IF (b%blocktype == c_blocktype_source) THEN
      CALL write_block_header(h)
    ELSE IF (b%blocktype == c_blocktype_datablock) THEN
      CALL write_datablock_meta(h)
    ELSE IF (b%blocktype == c_blocktype_stitched &
        .OR. b%blocktype == c_blocktype_contiguous &
        .OR. b%blocktype == c_blocktype_stitched_tensor &
        .OR. b%blocktype == c_blocktype_contiguous_tensor) THEN
      CALL sdf_write_stitched(h)
    ELSE IF (b%blocktype == c_blocktype_stitched_material &
        .OR. b%blocktype == c_blocktype_contiguous_material) THEN
      CALL sdf_write_stitched_material(h)
    ELSE IF (b%blocktype == c_blocktype_stitched_matvar &
        .OR. b%blocktype == c_blocktype_contiguous_matvar) THEN
      CALL sdf_write_stitched_matvar(h)
    ELSE IF (b%blocktype == c_blocktype_stitched_species &
        .OR. b%blocktype == c_blocktype_contiguous_species) THEN
      CALL sdf_write_stitched_species(h)
    ELSE IF (b%blocktype == c_blocktype_stitched_obstacle_group) THEN
      CALL sdf_write_stitched_obstacle_group(h)
    ELSE IF (b%blocktype == c_blocktype_namevalue) THEN
      CALL write_namevalue_meta(h)
    ELSE
      CALL write_block_header(h)
    END IF

  END SUBROUTINE sdf_write_block_info

END MODULE sdf_output_util

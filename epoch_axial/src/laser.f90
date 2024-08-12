! Copyright (C) 2009-2019 University of Warwick
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.

MODULE laser

  USE custom_laser
  USE evaluator

  IMPLICIT NONE

CONTAINS

  SUBROUTINE init_laser(boundary, laser)

    INTEGER, INTENT(IN) :: boundary
    TYPE(laser_block), INTENT(INOUT) :: laser

    laser%boundary = boundary
    laser%id = -1
    laser%use_time_function = .FALSE.
    laser%use_phase_function = .FALSE.
    laser%use_profile_function = .FALSE.
    laser%use_omega_function = .FALSE.
    laser%amp = -1.0_num
    laser%omega = -1.0_num
    laser%pol_angle = 0.0_num
    laser%t_start = 0.0_num
    laser%t_end = t_end
    laser%current_integral_phase = 0.0_num
    NULLIFY(laser%profile)
    NULLIFY(laser%phase)
    NULLIFY(laser%next)

    CALL allocate_with_boundary(laser%profile, boundary)
    CALL allocate_with_boundary(laser%phase, boundary)
    laser%profile = 1.0_num
    laser%phase = 0.0_num

  END SUBROUTINE init_laser



  SUBROUTINE setup_laser_phases(laser_init, phases)

    TYPE(laser_block), POINTER :: laser_init
    REAL(num), DIMENSION(:), INTENT(IN) :: phases
    TYPE(laser_block), POINTER :: laser
    INTEGER :: ilas

    ilas = 1
    laser => laser_init
    DO WHILE(ASSOCIATED(laser))
      laser%current_integral_phase = phases(ilas)
      ilas = ilas + 1
      laser => laser%next
    END DO

  END SUBROUTINE setup_laser_phases



  SUBROUTINE deallocate_laser(laser)

    TYPE(laser_block), POINTER :: laser

    IF (ASSOCIATED(laser%profile)) DEALLOCATE(laser%profile)
    IF (ASSOCIATED(laser%phase)) DEALLOCATE(laser%phase)
    IF (laser%use_profile_function) &
        CALL deallocate_stack(laser%profile_function)
    IF (laser%use_phase_function) &
        CALL deallocate_stack(laser%phase_function)
    IF (laser%use_time_function) &
        CALL deallocate_stack(laser%time_function)
    IF (laser%use_omega_function) &
        CALL deallocate_stack(laser%omega_function)
    DEALLOCATE(laser)

  END SUBROUTINE deallocate_laser



  SUBROUTINE deallocate_lasers

    TYPE(laser_block), POINTER :: current, next

    current => laser_x_min
    DO WHILE(ASSOCIATED(current))
      next => current%next
      CALL deallocate_laser(current)
      current => next
    END DO

    current => laser_x_max
    DO WHILE(ASSOCIATED(current))
      next => current%next
      CALL deallocate_laser(current)
      current => next
    END DO

    current => laser_y_min
    DO WHILE(ASSOCIATED(current))
      next => current%next
      CALL deallocate_laser(current)
      current => next
    END DO

    current => laser_y_max
    DO WHILE(ASSOCIATED(current))
      next => current%next
      CALL deallocate_laser(current)
      current => next
    END DO

  END SUBROUTINE deallocate_lasers



  ! Subroutine to attach a created laser object to the correct boundary
  SUBROUTINE attach_laser(laser)

    INTEGER :: boundary
    TYPE(laser_block), POINTER :: laser

    boundary = laser%boundary

    IF (boundary == c_bd_x_min) THEN
      n_laser_x_min = n_laser_x_min + 1
      CALL attach_laser_to_list(laser_x_min, laser)
    ELSE IF (boundary == c_bd_x_max) THEN
      n_laser_x_max = n_laser_x_max + 1
      CALL attach_laser_to_list(laser_x_max, laser)
    ELSE IF (boundary == c_bd_y_min) THEN
      n_laser_y_min = n_laser_y_min + 1
      CALL attach_laser_to_list(laser_y_min, laser)
    ELSE IF (boundary == c_bd_y_max) THEN
      n_laser_y_max = n_laser_y_max + 1
      CALL attach_laser_to_list(laser_y_max, laser)
    END IF

  END SUBROUTINE attach_laser



  ! This routine populates the constant elements of a parameter pack
  ! from a laser

  SUBROUTINE populate_pack_from_laser(laser, parameters)

    TYPE(laser_block), POINTER :: laser
    TYPE(parameter_pack), INTENT(INOUT) :: parameters

    parameters%pack_ix = 0
    parameters%pack_iy = 0

    SELECT CASE(laser%boundary)
      CASE(c_bd_x_min)
        parameters%pack_ix = 0
      CASE(c_bd_x_max)
        parameters%pack_ix = nx
      CASE(c_bd_y_min)
        parameters%pack_iy = 0
      CASE(c_bd_y_max)
        parameters%pack_iy = ny
    END SELECT

  END SUBROUTINE populate_pack_from_laser



  FUNCTION laser_time_profile(laser)

    TYPE(laser_block), POINTER :: laser
    REAL(num) :: laser_time_profile
    INTEGER :: err
    TYPE(parameter_pack) :: parameters

    err = 0
    CALL populate_pack_from_laser(laser, parameters)
    IF (laser%use_time_function) THEN
      laser_time_profile = evaluate_with_parameters(laser%time_function, &
          parameters, err)
      RETURN
    END IF

    laser_time_profile = custom_laser_time_profile(laser)

  END FUNCTION laser_time_profile



  SUBROUTINE laser_update_phase(laser)

    TYPE(laser_block), POINTER :: laser
    INTEGER :: i, err
    TYPE(parameter_pack) :: parameters

    err = 0
    CALL populate_pack_from_laser(laser, parameters)
    SELECT CASE(laser%boundary)
      CASE(c_bd_x_min, c_bd_x_max)
        DO i = 0,ny
          parameters%pack_iy = i
          laser%phase(i) = &
              evaluate_with_parameters(laser%phase_function, parameters, err)
        END DO
      CASE(c_bd_y_min, c_bd_y_max)
        DO i = 0,nx
          parameters%pack_ix = i
          laser%phase(i) = &
              evaluate_with_parameters(laser%phase_function, parameters, err)
        END DO
    END SELECT

  END SUBROUTINE laser_update_phase



  SUBROUTINE laser_update_profile(laser)

    TYPE(laser_block), POINTER :: laser
    INTEGER :: i, err
    TYPE(parameter_pack) :: parameters

    err = 0
    CALL populate_pack_from_laser(laser, parameters)
    SELECT CASE(laser%boundary)
      CASE(c_bd_x_min, c_bd_x_max)
        DO i = 0,ny
          parameters%pack_iy = i
          laser%profile(i) = &
              evaluate_with_parameters(laser%profile_function, parameters, err)
        END DO
      CASE(c_bd_y_min, c_bd_y_max)
        DO i = 0,nx
          parameters%pack_ix = i
          laser%profile(i) = &
              evaluate_with_parameters(laser%profile_function, parameters, err)
        END DO
    END SELECT

  END SUBROUTINE laser_update_profile



  SUBROUTINE laser_update_omega(laser)

    TYPE(laser_block), POINTER :: laser
    INTEGER :: err
    TYPE(parameter_pack) :: parameters

    err = 0
    CALL populate_pack_from_laser(laser, parameters)
    laser%omega = &
        evaluate_with_parameters(laser%omega_function, parameters, err)
    IF (laser%omega_func_type == c_of_freq) &
        laser%omega = 2.0_num * pi * laser%omega
    IF (laser%omega_func_type == c_of_lambda) &
        laser%omega = 2.0_num * pi * c / laser%omega

  END SUBROUTINE laser_update_omega



  SUBROUTINE update_laser_omegas

    TYPE(laser_block), POINTER :: current

    current => laser_x_min
    DO WHILE(ASSOCIATED(current))
      IF (current%use_omega_function) THEN
        CALL laser_update_omega(current)
        current%current_integral_phase = current%current_integral_phase &
            + current%omega * dt
      ELSE
        current%current_integral_phase = current%omega * time
      END IF
      current => current%next
    END DO

    current => laser_x_max
    DO WHILE(ASSOCIATED(current))
      IF (current%use_omega_function) THEN
        CALL laser_update_omega(current)
        current%current_integral_phase = current%current_integral_phase &
            + current%omega * dt
      ELSE
        current%current_integral_phase = current%omega * time
      END IF
      current => current%next
    END DO

    current => laser_y_min
    DO WHILE(ASSOCIATED(current))
      IF (current%use_omega_function) THEN
        CALL laser_update_omega(current)
        current%current_integral_phase = current%current_integral_phase &
            + current%omega * dt
      ELSE
        current%current_integral_phase = current%omega * time
      END IF
      current => current%next
    END DO

    current => laser_y_max
    DO WHILE(ASSOCIATED(current))
      IF (current%use_omega_function) THEN
        CALL laser_update_omega(current)
        current%current_integral_phase = current%current_integral_phase &
            + current%omega * dt
      ELSE
        current%current_integral_phase = current%omega * time
      END IF
      current => current%next
    END DO

  END SUBROUTINE update_laser_omegas



  ! Actually does the attaching of the laser to the correct list
  SUBROUTINE attach_laser_to_list(list, laser)

    TYPE(laser_block), POINTER :: list
    TYPE(laser_block), POINTER :: laser
    TYPE(laser_block), POINTER :: current

    IF (ASSOCIATED(list)) THEN
      current => list
      DO WHILE(ASSOCIATED(current%next))
        current => current%next
      END DO
      current%next => laser
    ELSE
      list => laser
    END IF

  END SUBROUTINE attach_laser_to_list



  SUBROUTINE allocate_with_boundary(array, boundary)

    REAL(num), DIMENSION(:), POINTER :: array
    INTEGER, INTENT(IN) :: boundary

    IF (boundary == c_bd_x_min .OR. boundary == c_bd_x_max) THEN
      ALLOCATE(array(1-ng:ny+ng))
    ELSE IF (boundary == c_bd_y_min .OR. boundary == c_bd_y_max) THEN
      ALLOCATE(array(1-ng:nx+ng))
    END IF

  END SUBROUTINE allocate_with_boundary



  SUBROUTINE set_laser_dt

    REAL(num) :: dt_local
    TYPE(laser_block), POINTER :: current

    dt_laser = HUGE(1.0_num)

    current => laser_x_min
    DO WHILE(ASSOCIATED(current))
      dt_local = 2.0_num * pi / current%omega
      dt_laser = MIN(dt_laser, dt_local)
      current => current%next
    END DO

    current => laser_x_max
    DO WHILE(ASSOCIATED(current))
      dt_local = 2.0_num * pi / current%omega
      dt_laser = MIN(dt_laser, dt_local)
      current => current%next
    END DO

    current => laser_y_min
    DO WHILE(ASSOCIATED(current))
      dt_local = 2.0_num * pi / current%omega
      dt_laser = MIN(dt_laser, dt_local)
      current => current%next
    END DO

    current => laser_y_max
    DO WHILE(ASSOCIATED(current))
      dt_local = 2.0_num * pi / current%omega
      dt_laser = MIN(dt_laser, dt_local)
      current => current%next
    END DO

    ! Need at least two iterations per laser period
    ! (Nyquist)
    dt_laser = dt_laser / 2.0_num

  END SUBROUTINE set_laser_dt



  SUBROUTINE outflow_bcs_x_min

    ! Calculates field corrections for outflow and laser boundaries on x_min

    REAL(num) :: t_env
    REAL(num) :: dtc2, lx, lr, sum, diff, dt_eps, base
    REAL(num), DIMENSION(:), ALLOCATABLE :: source1, source2
    COMPLEX(num), DIMENSION(:), ALLOCATABLE :: source_r, source_t
    REAL(num), DIMENSION(:), ALLOCATABLE :: r_d_vals
    INTEGER :: n, ir, im, ir_l, ir_h
    TYPE(laser_block), POINTER :: current

    n = c_bd_x_min

    dtc2 = dt * c**2
    lx = dtc2 / dx
    lr = dtc2 / dy
    sum = 1.0_num / (lx + c)
    diff = lx - c
    dt_eps = dt / epsilon0

    ! In cartesian co-ordinates, the x-propagating plane waves described by
    ! sources are:
    ! - source1: Ey, Bz plane wave
    ! - source2: Ez, By plane wave
    ALLOCATE(source1(0:ny))
    ALLOCATE(source2(0:ny))
    source1 = 0.0_num
    source2 = 0.0_num

    bxm(0, 0:ny, :) = bxm_x_min(0:ny, :)

    ! Evaulate Ey and Ez (source1, source2) on this laser boundary
    IF (add_laser(n)) THEN
      current => laser_x_min
      DO WHILE(ASSOCIATED(current))
        ! evaluate the temporal evolution of the laser
        IF (time >= current%t_start .AND. time <= current%t_end) THEN
          IF (current%use_phase_function) CALL laser_update_phase(current)
          IF (current%use_profile_function) CALL laser_update_profile(current)
          t_env = laser_time_profile(current) * current%amp
          DO ir = 0,ny
            base = t_env * current%profile(ir) &
              * SIN(current%current_integral_phase + current%phase(ir))
            source1(ir) = source1(ir) + base * COS(current%pol_angle)
            source2(ir) = source2(ir) + base * SIN(current%pol_angle)
          END DO
        END IF
        current => current%next
      END DO
    END IF

    ! Get cell-centred radial positions of each ir cell
    ALLOCATE(r_d_vals(0:ny))
    DO ir = 0,ny
      r_d_vals(ir) = ABS(REAL(ir-1, num) * dy + y_grid_min_local)
    END DO

    ! Allocate arrays for source terms - we convert Cartesian laser fields to
    ! cylindrical field modes
    ALLOCATE(source_t(0:ny))
    ALLOCATE(source_r(0:ny))
    source_t = 0.0_num
    source_r = 0.0_num
    ! Loop over all azimuthal modes
    DO im = 0, n_mode-1

      ! Currently, laser.f90 is only set to do axially symmetric plane waves.
      ! The source terms for these waves only exist for mode m=1
      ! source1 = 0.5 * (Ey + c*Bz) = Ey for a plane wave, source2 = Ez
      ! See eq. (1-3) in Lifschitz (2009) for Ey, Ez, Er1, Et1 conversion
      IF (im == 1) THEN
        source_t = source1 + imagi * source2
        source_r = -imagi * source1 + source2
      ELSE
        source_t = 0.0_num
        source_r = 0.0_num
      END IF

      btm(1,1:ny,im) = sum * ( 4.0_num * source_t(1:ny) &
          + 2.0_num * (erm_x_min(1:ny,im) + c * btm_x_min(1:ny,im)) &
          - 2.0_num * erm(1,1:ny,im) &
          + imagi * im * c**2 * dt * bxm(1,1:ny,im) / r_d_vals &
          + dt_eps * jrm(1,1:ny,im) &
          + diff * btm(2,1:ny,im))

      ! Set the low ir limit for the Brm field
      IF (y_min_boundary) THEN
        ! On these ranks, ir=0 is the r=0, where Br1 is always 0
        ir_l = 1
      ELSE
        ir_l = 0
      END IF
      IF (y_max_boundary) THEN 
        ir_h = ny-1
      ELSE
        ir_h = ny
      END IF
      brm(1,ir_l:ir_h,im) = sum * (-4.0_num * source_r(ir_l:ir_h) &
          - 2.0_num * (etm_x_min(ir_l:ir_h,im) + c * brm_x_min(ir_l:ir_h,im)) &
          + 2.0_num * etm(1,ir_l:ir_h,im) &
          - lr * (bxm(1,ir_l+1:ir_h+1,im) - bxm(1,ir_l:ir_h,im)) &
          - dt_eps * jtm(1,ir_l:ir_h,im) &
          + diff * brm(2,ir_l:ir_h,im))

    END DO

    DEALLOCATE(source1, source2, source_t, source_r, r_d_vals)

  END SUBROUTINE outflow_bcs_x_min



  SUBROUTINE outflow_bcs_x_max

    ! Calculates outflow and laser corrections to fields on the x_max boundary

    REAL(num) :: t_env
    REAL(num) :: dtc2, lx, lr, sum, diff, dt_eps, base
    REAL(num), DIMENSION(:), ALLOCATABLE :: source1, source2
    COMPLEX(num), DIMENSION(:), ALLOCATABLE :: source_t, source_r
    REAL(num), DIMENSION(:), ALLOCATABLE :: r_d_vals
    INTEGER :: laserpos, n, ir, im, ir_l, ir_h
    TYPE(laser_block), POINTER :: current

    n = c_bd_x_max

    dtc2 = dt * c**2
    lx = dtc2 / dx
    lr = dtc2 / dy
    sum = 1.0_num / (lx + c)
    diff = lx - c
    dt_eps = dt / epsilon0

    ! In cartesian co-ordinates, the x-propagating plane waves described by
    ! sources are:
    ! - source1: Ey, Bz plane wave
    ! - source2: Ez, By plane wave
    ALLOCATE(source1(0:ny))
    ALLOCATE(source2(0:ny))
    source1 = 0.0_num
    source2 = 0.0_num

    bxm(nx, 0:ny, :) = bxm_x_max(0:ny, :)

    ! Evaulate Ey and Ez (source1, source2) on this laser boundary
    IF (add_laser(n)) THEN
      current => laser_x_max
      DO WHILE(ASSOCIATED(current))
        ! evaluate the temporal evolution of the laser
        IF (time >= current%t_start .AND. time <= current%t_end) THEN
          IF (current%use_phase_function) CALL laser_update_phase(current)
          IF (current%use_profile_function) CALL laser_update_profile(current)
          t_env = laser_time_profile(current) * current%amp
          DO ir = 0,ny
            base = t_env * current%profile(ir) &
              * SIN(current%current_integral_phase + current%phase(ir))
            source1(ir) = source1(ir) + base * COS(current%pol_angle)
            source2(ir) = source2(ir) + base * SIN(current%pol_angle)
          END DO
        END IF
        current => current%next
      END DO
    END IF

    ! Get cell-centred radial positions of each ir cell
    ALLOCATE(r_d_vals(0:ny))
    DO ir = 0,ny
      r_d_vals(ir) = ABS(REAL(ir-1, num) * dy + y_grid_min_local)
    END DO

    ! Allocate arrays for source terms - we convert Cartesian laser fields to
    ! cylindrical field modes
    ALLOCATE(source_t(0:ny))
    ALLOCATE(source_r(0:ny))
    source_t = 0.0_num
    source_r = 0.0_num

    ! Loop over all azimuthal modes
    DO im = 0, n_mode-1

      ! Currently, laser.f90 is only set to do axially symmetric plane waves.
      ! The source terms for these waves only exist for mode m=1
      ! source1 = 0.5 * (Ey + c*Bz) = Ey for a plane wave, source2 = Ez
      ! See eq. (1-3) in Lifschitz (2009) for Ey, Ez, Er1, Et1 conversion
      IF (im == 1) THEN
        source_t = source1 + imagi * source2
        source_r = -imagi * source1 + source2
      ELSE
        source_t = 0.0_num
        source_r = 0.0_num
      END IF

      btm(nx,1:ny,im) = sum * ( -4.0_num * source_t &
          - 2.0_num * (erm_x_max(1:ny,im) + c * btm_x_max(1:ny,im)) &
          + 2.0_num * erm(nx-1,1:ny,im) &
          - imagi * im * c**2 * dt * bxm(nx-1,1:ny,im) / r_d_vals &
          - dt_eps * jrm(nx-1,1:ny,im) &
          + diff * btm(nx-1,1:ny,im))

      ! Set the low ir limit for the Brm field
      IF (y_min_boundary) THEN
        ! On these ranks, ir=1 is the x-axis, where Br1 is always 0
        ir_l = 1
      ELSE
        ir_l = 0
      END IF
      IF (y_max_boundary) THEN 
        ir_h = ny-1
      ELSE
        ir_h = ny
      END IF
      brm(nx,ir_l:ir_h,im) = sum * ( 4.0_num * source_r(ir_l:ir_h) &
          + 2.0_num * (etm_x_max(ir_l:ir_h,im) + c * brm_x_max(ir_l:ir_h,im)) &
          - 2.0_num * etm(nx-1,ir_l:ir_h,im) &
          + lr * (bxm(nx-1,ir_l+1:ir_h+1,im) - bxm(nx-1,ir_l:ir_h,im)) &
          + dt_eps * jtm(nx-1,ir_l:ir_h,im) &
          + diff * brm(nx-1,ir_l:ir_h,im))
    END DO

    DEALLOCATE(source1, source2, source_r, source_t, r_d_vals)

  END SUBROUTINE outflow_bcs_x_max



  SUBROUTINE outflow_bcs_r_max

    REAL(num) :: t_env
    REAL(num) :: dtc2, inv_r, dtc2_4r, icdt_2r, lx, ly, sum_x, sum_t, dt_2eps
    REAL(num), DIMENSION(:), ALLOCATABLE :: source1, source2
    INTEGER :: i, im, ix_l, ix_h
    TYPE(laser_block), POINTER :: current

    dtc2 = dt * c**2
    inv_r = 1.0_num / (REAL(ny-1.5, num) * dy + y_grid_min_local)
    dtc2_4r = 0.25_num * dtc2 * inv_r
    icdt_2r = 0.5_num * imagi * c * dt * inv_r
    lx = dtc2 / dx
    ly = dtc2 / dy
    sum_x = 1.0_num / (ly + c)
    sum_t = 1.0_num / (ly + c + dtc2_4r)
    dt_2eps = 0.5_num * dt / epsilon0

    ! Ignore cells in the high-r simulation corners. Leave these for x-bcs
    IF (x_min_boundary) THEN
      ix_l = 1
    ELSE
      ix_l = 0
    END IF

    IF (x_max_boundary) THEN
      ix_h = nx-1
    ELSE
      ix_h = nx
    END IF

    DO im = 0, n_mode-1
      bxm(ix_l:ix_h,ny,im) = sum_x * (-bxm(ix_l:ix_h,ny-1,im) * (c - ly) &
          -bxm_old(ix_l:ix_h,ny,im) * (-c + ly) &
          -bxm_old(ix_l:ix_h,ny-1,im) * (-c - ly) &
          -c * dt * inv_r * etm(ix_l:ix_h,ny-1,im) &
          + 0.5_num * lx * (brm(ix_l+1:ix_h+1,ny-1,im) &
          - brm(ix_l:ix_h,ny-1,im) &
          + brm_old(ix_l+1:ix_h+1,ny-1,im) - brm_old(ix_l:ix_h,ny-1,im)) &
          - icdt_2r * REAL(im, num) * (erm(ix_l:ix_h,ny,im) &
          + erm(ix_l:ix_h,ny-1,im)) &
          - dt_2eps * (jtm(ix_l:ix_h,ny-1,im) + jtm_old(ix_l:ix_h,ny-1,im)))

      btm(1:nx,ny,im) = sum_t * (-btm(1:nx,ny-1,im) * (c - ly + dtc2_4r) &
          -btm_old(1:nx,ny,im) * (-c + ly + dtc2_4r) &
          -btm_old(1:nx,ny-1,im) * (-c - ly + dtc2_4r) &
          - 0.5_num * lx / c * (erm(1:nx,ny,im) + erm(1:nx,ny-1,im) &
          - erm(0:nx-1,ny,im) - erm(0:nx-1,ny-1,im)) &
          - icdt_2r * REAL(im, num) * c * (brm(1:nx,ny-1,im) &
          + brm_old(1:nx,ny-1,im)) &
          + dt_2eps * (jxm(1:nx,ny-1,im) + jxm_old(1:nx,ny-1,im)))
    END DO

  END SUBROUTINE outflow_bcs_r_max



  SUBROUTINE calc_absorption(bd, lasers)

    TYPE(laser_block), POINTER, OPTIONAL :: lasers
    INTEGER, INTENT(IN) :: bd
    TYPE(laser_block), POINTER :: current
    REAL(num) :: t_env, dir, dd, factor, lfactor, laser_inject_sum
    REAL(num), DIMENSION(:), ALLOCATABLE :: e1, e2, b1, b2
    INTEGER :: mm, ibc, icell

    ! Note: ideally e1, e2, b1, b2 should be face-centred. However, this is not
    ! possible with 'open' boundaries since E-fields are not defined in the
    ! ghost cell, so we use the cell-centred quantities in the first cell.

    dir = 1.0_num
    mm = 1

    SELECT CASE(bd)
      CASE(c_bd_x_min, c_bd_x_max)
        dd = dy
        mm = ny
        ALLOCATE(e1(mm), e2(mm), b1(mm), b2(mm))

        ibc = 1
        IF (bd == c_bd_x_max) THEN
          dir = -1.0_num
          ibc = nx
        END IF

        e1 = 0.5_num  * (ey(ibc  , 0:ny-1) + ey(ibc, 1:ny  ))
        e2 = ez(ibc, 1:ny)
        b1 = 0.25_num * (bz(ibc-1, 0:ny-1) + bz(ibc, 0:ny-1) &
                       + bz(ibc-1, 1:ny  ) + bz(ibc, 1:ny  ))
        b2 = 0.5_num  * (by(ibc-1, 1:ny  ) + by(ibc, 1:ny  ))

      CASE(c_bd_y_min, c_bd_y_max)
        dd = dx
        mm = nx
        ALLOCATE(e1(mm), e2(mm), b1(mm), b2(mm))

        ibc = 1
        IF (bd == c_bd_y_max) THEN
          dir = -1.0_num
          ibc = ny
        END IF

        e1 = ez(1:nx, ibc)
        e2 = 0.5_num  * (ex(0:nx-1, ibc  ) + ex(1:nx  , ibc))
        b1 = 0.5_num  * (bx(1:nx  , ibc-1) + bx(1:nx  , ibc))
        b2 = 0.25_num * (bz(0:nx-1, ibc-1) + bz(0:nx-1, ibc) &
                       + bz(1:nx  , ibc-1) + bz(1:nx  , ibc))

      CASE DEFAULT
        dd = 0.0_num
        ALLOCATE(e1(mm), e2(mm), b1(mm), b2(mm))

        e1 = 0.0_num
        e2 = 0.0_num
        b1 = 0.0_num
        b2 = 0.0_num
    END SELECT

    factor = dt * dd * dir
    laser_absorb_local = laser_absorb_local &
        + (factor / mu0) * SUM(e1 * b1 - e2 * b2)

    IF (PRESENT(lasers)) THEN
      current => lasers
      DO WHILE(ASSOCIATED(current))
        laser_inject_sum = 0.0_num
        DO icell = 1, mm
          laser_inject_sum = laser_inject_sum + current%profile(icell)**2
        END DO
        t_env = laser_time_profile(current)
        lfactor = 0.5_num * epsilon0 * c * factor * (t_env * current%amp)**2
        laser_inject_local = laser_inject_local + lfactor * laser_inject_sum
        current => current%next
      END DO
    END IF

    DEALLOCATE(e1, e2, b1, b2)

  END SUBROUTINE calc_absorption

END MODULE laser

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
!
! Field stagger in the bottom-left cell in the simulation grid
!
!   ^ Increasing r               |
!   |                            |
!   |----------------------------|----
!   |                            |
!   |                            |
!  Bxm(0,1,:)                    |
!  Jrm(0,1,:)   Btm(1,1,:)       |
!  Erm(0,1,:)                    |
!   |                            |
!   |           Brm(1,0,:)       |
!  Jtm(0,0,:)---Jxm(1,0,:)-------|---- > Increasing x
!  Etm(0,0,:)   Exm(1,0,:)
!
! Field stagger in the top-right cell in the simulation grid
!
!              Brm(nx,ny,:)              
!   |--------- Jxm(nx,ny,:)---Jtm(nx,ny,:)
!   |          Exm(nx,ny,:)   Etm(nx,ny,:)
!   |                            | 
!   |                            |
!   |                         Bxm(nx,ny,:)                    
!   |          Btm(nx,ny,:)   Jrm(nx,ny,:)
!   |                         Erm(nx,ny,:)                    
!   |                            |
!   |                            |
!   |----------------------------|

MODULE fields

  USE boundary

  IMPLICIT NONE

CONTAINS

  SUBROUTINE update_e_field

    ! Updates the electric field modes of the Fourier decomposition by half a
    ! time-step, as discussed in Sections 2 and 3 of Lifschitz (2009). These
    ! equations have been converted to SI units here
    
    INTEGER :: ix, ir, im, ir_min
    REAL(num) :: r_p, r_d, fac_x, mode_sign
    COMPLEX(num) :: im_fac_x, im_fac_r

    ir_min = 0
    IF (y_min_boundary) ir_min = 1

    ! Loop over field modes, m
    DO im = 0, n_mode-1

      ! Loop over radial cells
      DO ir = ir_min, ny

        ! 'Dual' r refers to the radial distance from axis to cell centre
        ! 'Primal' r refers to the r distance from axis to high-r cell-edge
        ! r_grid_min_local is the radial distance from axis to cell centre
        ! in the cell with the local index ir=1
        r_d = ABS(REAL(ir - 1, num) * dy + y_grid_min_local)
        r_p = r_d + 0.5_num * dy
        fac_x = c**2 / r_p
        im_fac_x = imagi * im * fac_x
        im_fac_r = imagi * im * c**2 / r_d

        ! Loop over x cells, along cylindrical axis
        DO ix = 0, nx

          ! exm(x,r)m
          ! (11) - Lifschitz (2009)
          exm(ix, ir, im) = exm(ix, ir, im) &
              + (fac_x * 0.5_num * (btm(ix, ir+1, im) + btm(ix, ir, im)) &
              + im_fac_x * brm(ix, ir, im) &
              + c**2 * (btm(ix, ir+1, im) - btm(ix, ir, im)) / dy &
              - jxm(ix, ir, im) / epsilon0) * 0.5_num * dt

          ! erm(x,r)m
          ! (9) - Lifschitz (2009)
          erm(ix, ir, im) = erm(ix, ir, im) &
              + (-im_fac_r * bxm(ix, ir, im) &
              - c**2 * (btm(ix+1, ir, im) - btm(ix, ir, im)) / dx &
              - jrm(ix, ir, im) / epsilon0) * 0.5_num * dt

          ! etm(x,r)m
          ! (10) - Lifschitz (2009)
          etm(ix, ir, im) = etm(ix, ir, im) &
              + (c**2 * (brm(ix+1, ir, im) - brm(ix, ir, im)) / dx &
              - c**2 * (bxm(ix, ir+1, im) - bxm(ix, ir, im)) / dy &
              - jtm(ix, ir, im) / epsilon0) * 0.5_num * dt
        END DO
      END DO
    END DO

    ! r_min boundary conditions and ghost cells for ranks with r=0
    ! Below r=0, Exm switches sign for odd m, no switch for even
    !            Erm, Etm switch sign for even m, no switch for odd
    ! Multi-value: Ex has no theta variation for m=0, E_perp has no theta 
    !              variation for m=1. All other modes vary with theta, and
    !              imply multiple-values at r=0, which is forbidden. Set to 0.
    IF (y_min_boundary) THEN

      ! m = 0 conditions:

      ! (20) - Lifschitz (2009)
      exm(1-ng:nx+ng, 0, 0) = exm(1-ng:nx+ng, 0, 0) &
          + (4.0_num * c**2 / dy * btm(1-ng:nx+ng, 1, 0) &
          - jxm(1-ng:nx+ng, 0, 0) / epsilon0) * 0.5_num * dt

      ! Multi-value forbidden
      etm(1-ng:nx+ng,0,0) = 0.0_num
      erm(1-ng:nx+ng,0,0) = -erm(1-ng:nx+ng,1,0)

      ! Physical reflection below r=0
      DO ir = 1-ng,-1
        etm(1-ng:nx+ng,ir,0) = - etm(1-ng:nx+ng,-ir,0)
        erm(1-ng:nx+ng,ir,0) = - erm(1-ng:nx+ng,-ir+1,0)
        exm(1-ng:nx+ng,ir,0) = + exm(1-ng:nx+ng,-ir,0)
      END DO
      
      ! m = 1 conditions:
      IF (n_mode > 1) THEN 
        
        ! Multi-value forbidden
        exm(1-ng:nx+ng,0,1) = 0.0_num

        ! dEy/dtheta = 0 on r=0 implies Erm = i*Etm/m
        ! Also use div(E) = rho/epsilon0, where Ex=rho=0 for m=1, r=0
        erm(1-ng:nx+ng,0,1) = 2.0_num * imagi * etm(1-ng:nx+ng,0,1) &
            - erm(1-ng:nx+ng,1,1)
        etm(1-ng:nx+ng,0,1) = -imagi / 8.0_num * &
            (9.0_num*erm(1-ng:nx+ng,1,1) - erm(1-ng:nx+ng,2,1))

        ! Physical reflection below r=0
        DO ir = 1-ng,-1
          etm(1-ng:nx+ng,ir,1) = + etm(1-ng:nx+ng,-ir,1)
          erm(1-ng:nx+ng,ir,1) = + erm(1-ng:nx+ng,-ir+1,1)
          exm(1-ng:nx+ng,ir,1) = - exm(1-ng:nx+ng,-ir,1)
        END DO
      END IF

      ! m > 1 conditions:
      IF (n_mode > 2) THEN
        mode_sign = 1.0_num
        DO im = 2, n_mode-1

          ! Multi-value forbidden
          exm(1-ng:nx+ng,0,im) = 0.0_num 
          etm(1-ng:nx+ng,0,im) = 0.0_num
          erm(1-ng:nx+ng,0,im) = -erm(1-ng:nx+ng,1,im)

          ! dEy/dtheta = 0, div(E) = rho/epsilon0, etm=0 on r=0 for m>1
          erm(1-ng:nx+ng,1,im) = erm(1-ng:nx+ng,2,im) / 9.0_num

          ! Physical reflection below r=0
          DO ir = 1-ng,-1
            etm(1-ng:nx+ng,ir,im) = - mode_sign * etm(1-ng:nx+ng,-ir,im)
            erm(1-ng:nx+ng,ir,im) = - mode_sign * erm(1-ng:nx+ng,-ir+1,im)
            exm(1-ng:nx+ng,ir,im) = + mode_sign * exm(1-ng:nx+ng,-ir,im)
          END DO

          mode_sign = -mode_sign
        END DO
      END IF
    END IF

  END SUBROUTINE update_e_field



  SUBROUTINE update_b_field

    ! Updates the magnetic field modes of the Fourier decomposition by half a
    ! time-step, as discussed in Sections 2 and 3 of Lifschitz (2009). These
    ! equations have been converted to SI units here

    INTEGER :: ix, ir, im, ir_min, ir_max
    REAL(num) :: r_p, r_d, mode_sign
    COMPLEX(num) :: im_fac_x, im_fac_r
    COMPLEX(num) :: dbxm_dt_i(3), dbxm_dt

    ir_min = 0
    IF (y_min_boundary) ir_min = 1
    ir_max = ny 
    IF (y_max_boundary) ir_max = ny-1

    ! Loop over field modes, m
    DO im = 0, n_mode-1

      ! Loop over radial positions
      DO ir = ir_min, ir_max

        ! 'Dual' r refers to the radial distance from axis to cell centre
        ! 'Primal' r refers to the radial distance from axis to high-r cell-edge
        ! r_grid_min_local is the radial distance from axis to cell centre in
        ! the cell with the local index ir=1
        r_d = ABS(REAL(ir - 1, num) * dy + y_grid_min_local)
        r_p = r_d + 0.5_num * dy

        im_fac_x = imagi * im / r_d
        im_fac_r = imagi * im / r_p

        ! Loop over x, along cylindrical axis
        DO ix = 0, nx

          ! bxm(x,r)m
          ! (8) - Lifschitz (2009)
          bxm(ix, ir, im) = bxm(ix, ir, im) &
              - (im_fac_x * erm(ix, ir, im) &
              + 0.5_num * (etm(ix, ir, im) + etm(ix, ir-1, im)) / r_d &
              + (etm(ix, ir, im) - etm(ix, ir-1, im)) / dy) * 0.5_num * dt

          ! brm(x,r)m
          ! (6) - Lifschitz (2009)
          brm(ix, ir, im) = brm(ix, ir, im) &
              + (im_fac_r * exm(ix, ir, im) &
              + (etm(ix, ir, im) - etm(ix-1, ir, im)) / dx) * 0.5_num * dt

          ! btm(x,r)m
          ! (7) - Lifschitz (2009)
          btm(ix, ir, im) = btm(ix, ir, im) &
              + (-(erm(ix, ir, im) - erm(ix-1, ir, im)) / dx &
              + (exm(ix, ir, im) - exm(ix, ir-1, im) ) / dy) * 0.5_num * dt
        END DO
      END DO
    END DO

    ! r_min boundary conditions and ghost cells for ranks with r=0
    ! Below r=0, Bxm switches sign for odd m, no switch for even
    !            Brm, Btm switch sign for even m, no switch for odd
    ! Multi-value: Bx has no theta variation for m=0, B_perp has no theta 
    !              variation for m=1. All other modes vary with theta, and
    !              imply multiple-values at r=0, which is forbidden. Set to 0.
    IF (y_min_boundary) THEN

      ! m = 0 conditions:

      ! Multi-value forbidden
      brm(1-ng:nx+ng,0,0) = 0.0_num

      ! Physical reflection below r=0
      bxm(1-ng:nx+ng,0,0) = + bxm(1-ng:nx+ng,1,0)
      btm(1-ng:nx+ng,0,0) = - btm(1-ng:nx+ng,1,0)
      DO ir = 1-ng,-1
        btm(1-ng:nx+ng,ir,0) = - btm(1-ng:nx+ng,-ir+1,0)
        brm(1-ng:nx+ng,ir,0) = - brm(1-ng:nx+ng,-ir,0)
        bxm(1-ng:nx+ng,ir,0) = + bxm(1-ng:nx+ng,-ir+1,0)
      END DO
      
      ! m = 1 conditions:
      IF (n_mode > 1) THEN 
        
        ! Multi-value forbidden
        bxm(1-ng:nx+ng,0,1) = -bxm(1-ng:nx+ng,1,1)

        ! (22) - Lifschitz (2009)
        brm(2-ng:nx+ng,0,1) = brm(2-ng:nx+ng,0,1) &
            + (imagi/dy * exm(2-ng:nx+ng,0,1) &
            + (etm(2-ng:nx+ng,0,1) - etm(1-ng:nx+ng-1,0,1))/dx) * 0.5_num * dt

        ! dBy/dtheta = 0 on r=0, implies Btm = -i*m*Brm
        btm(1-ng:nx+ng,0,1) = -2 * imagi * brm(1-ng:nx+ng,0,1) &
            - btm(1-ng:nx+ng,1,1)

        ! Physical reflection below r=0
        DO ir = 1-ng,-1
          btm(1-ng:nx+ng,ir,1) = + btm(1-ng:nx+ng,-ir+1,1)
          brm(1-ng:nx+ng,ir,1) = + brm(1-ng:nx+ng,-ir,1)
          bxm(1-ng:nx+ng,ir,1) = - bxm(1-ng:nx+ng,-ir+1,1)
        END DO
      END IF

      ! m > 1 conditions:
      IF (n_mode > 2) THEN
        mode_sign = 1.0_num
        DO im = 2, n_mode-1

          ! Multi-value forbidden
          bxm(1-ng:nx+ng,0,im) = -bxm(1-ng:nx+ng,1,im)
          brm(1-ng:nx+ng,0,im) = 0.0_num

          ! dBy/dtheta = 0, implies Btm=-i*m*Brm on r=0, but btm=0 here for m>1
          btm(1-ng:nx+ng,0,im) = -btm(1-ng:nx+ng,1,im)

          ! Physical reflection below r=0
          DO ir = 1-ng,-1
            btm(1-ng:nx+ng,ir,im) = - mode_sign * btm(1-ng:nx+ng,-ir+1,im)
            brm(1-ng:nx+ng,ir,im) = - mode_sign * brm(1-ng:nx+ng,-ir,im)
            bxm(1-ng:nx+ng,ir,im) = + mode_sign * bxm(1-ng:nx+ng,-ir+1,im)
          END DO

          mode_sign = -mode_sign
        END DO
      END IF
    END IF

  END SUBROUTINE update_b_field



  SUBROUTINE update_eb_fields_half

    ! Update E field to t+dt/2
    CALL update_e_field

    ! Now have E(t+dt/2), do boundary conditions on E
    CALL efield_bcs

    ! Save B field at t
    ! Update B field to t+dt/2 using E(t+dt/2)
    bxm_old = bxm 
    brm_old = brm 
    btm_old = btm
    CALL update_b_field

    ! Now have B field at t+dt/2. Do boundary conditions on B
    CALL bfield_bcs(.TRUE.)

    ! Now have E&B fields at t = t+dt/2
    ! Move to particle pusher

  END SUBROUTINE update_eb_fields_half



  SUBROUTINE update_eb_fields_final

    CALL update_b_field

    ! Do boundary conditions considering lasers and outflow boundaries
    ! We need B at t and t+dt for the r_max boundary
    CALL bfield_final_bcs

    CALL update_e_field

    CALL efield_bcs

  END SUBROUTINE update_eb_fields_final

END MODULE fields

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

MODULE particles

  USE boundary
  USE partlist
#ifdef PREFETCH
  USE prefetch
#endif

  IMPLICIT NONE

CONTAINS

  SUBROUTINE push_particles

    ! 2nd order accurate particle pusher using parabolic weighting
    ! on and off the grid. The calculation of J looks rather odd
    ! Since it works by solving d(rho)/dt = div(J) and doing a 1st order
    ! Estimate of rho(t+1.5*dt) rather than calculating J directly
    ! This gives exact charge conservation on the grid

    ! Contains the integer cell position of the particle in x, y, z
    INTEGER :: cell_x1, cell_x2, cell_x3
    INTEGER :: cell_y1, cell_y2, cell_y3

    ! Xi (space factor see page 38 in manual)
    ! The code now uses gx and hx instead of xi0 and xi1

    ! J from a given particle, can be spread over up to 3 cells in
    ! Each direction due to parabolic weighting. We allocate 4 or 5
    ! Cells because the position of the particle at t = t+1.5dt is not
    ! known until later. This part of the algorithm could probably be
    ! Improved, but at the moment, this is just a straight copy of
    ! The core of the PSC algorithm
    INTEGER, PARAMETER :: sf0 = sf_min, sf1 = sf_max
  !  REAL(num) :: jxh, jzh
    COMPLEX(num), DIMENSION(sf0-1:sf1+1) :: jyh
    COMPLEX(num) :: jxh, jzh

    ! Properties of the current particle. Copy out of particle arrays for speed
    REAL(num) :: part_x, part_y, part_z
    REAL(num) :: part_ux, part_uy, part_uz
    REAL(num) :: part_q, part_mc, ipart_mc, part_weight, part_m
#ifdef HC_PUSH
    REAL(num) :: beta_x, beta_y, beta_z, beta2, beta_dot_u, alpha, sigma
#endif 

    ! Additional variables for cylindrical co-ordinates
    REAL(num) :: part_r
    COMPLEX(num) :: exp_min_itheta, exp_min_imtheta, exp_min_itheta_05
    COMPLEX(num) :: exp_itheta_05, exp_itheta_10, exp_itheta_15
    REAL(num) :: part_x_local, part_r_local
    REAL(num) :: part_vy, part_vz, part_vt

    ! Used for particle probes (to see of probe conditions are satisfied)
#ifndef NO_PARTICLE_PROBES
    REAL(num) :: init_part_x, final_part_x
    REAL(num) :: init_part_y, final_part_y
    REAL(num) :: init_part_z, final_part_z
    TYPE(particle_probe), POINTER :: current_probe
    TYPE(particle), POINTER :: particle_copy
    REAL(num) :: d_init, d_final
    REAL(num) :: probe_energy, part_mc2
#endif

    ! Contains the floating point version of the cell number (never actually
    ! used)
    REAL(num) :: cell_x_r, cell_y_r

    ! The fraction of a cell between the particle position and the cell boundary
    REAL(num) :: cell_frac_x, cell_frac_y

    ! Weighting factors as Eqn 4.77 page 25 of manual
    ! Eqn 4.77 would be written as
    ! F(j-1) * gmx + F(j) * g0x + F(j+1) * gpx
    ! Defined at the particle position
    REAL(num), DIMENSION(sf_min-1:sf_max+1) :: gx, gy

    ! Defined at the particle position - 0.5 grid cell in each direction
    ! This is to deal with the grid stagger
    REAL(num), DIMENSION(sf_min-1:sf_max+1) :: hx, hy

    ! Fields at particle location
    INTEGER :: im
    REAL(num) :: ex_part, ey_part, ez_part, er_part, et_part 
    REAL(num) :: bx_part, by_part, bz_part, br_part, bt_part

    ! P+, P- and Tau variables from Boris1970, page27 of manual
    REAL(num) :: uxp, uxm, uyp, uym, uzp, uzm
    REAL(num) :: tau, taux, tauy, tauz, taux2, tauy2, tauz2

    ! charge to mass ratio modified by normalisation
    REAL(num) :: cmratio, ccmratio

    ! Used by J update
    INTEGER :: xmin, xmax, ymin, ymax
    REAL(num), ALLOCATABLE :: inv_area_rt(:), inv_area_xt(:), inv_volume(:)
    REAL(num), ALLOCATABLE :: ratio_area_xt(:)
    REAL(num) :: r_low 
    REAL(num) :: theta_05, theta_15, dtheta, mdth, m2dth2, inv_mdth, inv_m2dth2
    LOGICAL :: small_angle_change
    COMPLEX(num) :: w_xr, w_xt, w_rt
    COMPLEX(num) :: exp_idtheta, exp_imdtheta, exp_imtheta0
    COMPLEX(num) :: m_fac_1, m_fac_2, m_fac_3, m_fac_4
    COMPLEX(num) :: ym_fac_1

    ! Temporary variables
    REAL(num) :: idx, idy
    REAL(num) :: idty, idtx, idxy
    REAL(num) :: idt, dto2, dtco2
    REAL(num) :: q_weight_fac
    REAL(num) :: fcx, fcy, fcz, fjx, fjy, fjz
    REAL(num) :: root, dtfac, gamma_rel, part_u2, third, sixth, igamma
    REAL(num) :: delta_x, delta_y, delta_z
    REAL(num) :: hy_iy, xfac1, yfac1, yfac2
    INTEGER :: ispecies, ix, iy, dcellx, dcelly, cx, cy
    INTEGER(i8) :: ipart
#ifdef WORK_DONE_INTEGRATED
    REAL(num) :: tmp_x, tmp_y, tmp_z
    REAL(num) :: work_x, work_y, work_z
#endif
#ifndef NO_PARTICLE_PROBES
    LOGICAL :: probes_for_species
    REAL(num) :: gamma_rel_m1
    REAL(num) :: path_frac
#endif
#ifndef NO_TRACER_PARTICLES
    LOGICAL :: not_zero_current_species
#endif
    ! Particle weighting multiplication factor
#ifdef PARTICLE_SHAPE_BSPLINE3
    REAL(num) :: cf2
    REAL(num), PARAMETER :: fac = (1.0_num / 24.0_num)**c_ndims
#elif  PARTICLE_SHAPE_TOPHAT
    REAL(num), PARAMETER :: fac = 1.0_num
#else
    REAL(num) :: cf2
    REAL(num), PARAMETER :: fac = (0.5_num)**c_ndims
#endif
#ifdef DELTAF_METHOD
    REAL(num) :: weight_back
#endif

    TYPE(particle), POINTER :: current, next
#ifdef PREFETCH
    CALL prefetch_particle(species_list(1)%attached_list%head)
#endif

    jxm_old = jxm 
    jrm_old = jrm 
    jtm_old = jtm

    jxm = 0.0_num
    jrm = 0.0_num
    jtm = 0.0_num

    gx = 0.0_num
    gy = 0.0_num

    ! Unvarying multiplication factors

    idx = 1.0_num / dx
    idy = 1.0_num / dy
    idt = 1.0_num / dt
    dto2 = dt / 2.0_num
    dtco2 = c * dto2
    dtfac = 0.5_num * dt * fac
    third = 1.0_num / 3.0_num
    sixth = 0.5_num * third

    idty = idt * idy * fac
    idtx = idt * idx * fac
    idxy = idx * idy * fac

    ! Arrays containing 1.0/[areas of staggered cells] and 1.0/[volume]
    ALLOCATE(inv_area_rt(1-jng:ny+jng), inv_area_xt(0-jng:ny+jng))
    ALLOCATE(inv_volume(1-jng:ny+jng))
    
    ! Area of low-r xt face for the first cell
    r_low = y_grid_min_local - jng * dy
    inv_area_xt(0-jng) = 1.0_num / (2.0_num * pi * ABS(r_low) * dx)
    
    DO iy = 1-jng,ny+jng        
      ! Area of rt face of cell
      IF (NINT(2.0_num * r_low / dy) == -1) THEN 
        ! If r_low = -0.5*dy, r_high = +0.5*dy - circle of radius 0.5*dy
        inv_area_rt(iy) = 1.0_num / (pi * (0.5_num * dy)**2)
      ELSE
        ! Otherwise, area is ring between two radii  
        inv_area_rt(iy) = 1.0_num / (pi * ABS((r_low+dy)**2 - r_low**2))
      END IF

      ! Area of high-r xt face of the cell
      inv_area_xt(iy) = 1.0_num / (2.0_num * pi * ABS(r_low + dy) * dx)

      r_low = r_low + dy       
    END DO

    inv_volume = inv_area_rt / dx

    ! Ratio of low-r xt face area to high-r xt face area
    ALLOCATE(ratio_area_xt(1-jng:ny+jng))
    ratio_area_xt = inv_area_xt(1-jng:ny+jng) / inv_area_xt(0-jng:ny+jng-1)

    DO ispecies = 1, n_species
      current => species_list(ispecies)%attached_list%head
      IF (species_list(ispecies)%immobile) CYCLE
      IF (species_list(ispecies)%ignore_fields) THEN 
        CALL no_field_push(ispecies)
        CYCLE 
      END IF
      IF (species_list(ispecies)%species_type == c_species_id_photon) THEN
#ifdef BREMSSTRAHLUNG
        IF (ispecies == bremsstrahlung_photon_species) THEN
          IF (bremsstrahlung_photon_dynamics) THEN
            CALL push_photons(ispecies)
          ELSE
            CYCLE
          END IF
        END IF
#endif
#ifdef PHOTONS
        IF (photon_dynamics) CALL push_photons(ispecies)
#endif
        CYCLE
      END IF
#ifndef NO_PARTICLE_PROBES
      current_probe => species_list(ispecies)%attached_probes
      probes_for_species = ASSOCIATED(current_probe)
#if defined(PARTICLE_ID) || defined(PARTICLE_ID4) 
      IF (probes_for_species) THEN 
        CALL generate_particle_ids(species_list(ispecies)%attached_list)
      END IF
#endif
#endif
#ifndef NO_TRACER_PARTICLES
      not_zero_current_species = .NOT. species_list(ispecies)%zero_current
#endif

#ifdef PER_SPECIES_WEIGHT
      part_weight = species_list(ispecies)%weight
      fcx = idty * part_weight
      fcy = idtx * part_weight
      fcz = idxy * part_weight
#endif
#ifndef PER_PARTICLE_CHARGE_MASS
      part_q   = species_list(ispecies)%charge
      part_m   = species_list(ispecies)%mass
      part_mc  = c * species_list(ispecies)%mass
      ipart_mc = 1.0_num / part_mc
      cmratio  = part_q * dtfac * ipart_mc
      ccmratio = c * cmratio
#ifndef NO_PARTICLE_PROBES
      part_mc2 = c * part_mc
#endif
#endif
      !DEC$ VECTOR ALWAYS
      DO ipart = 1, species_list(ispecies)%attached_list%count

        next => current%next
#ifdef PREFETCH
        CALL prefetch_particle(next)
#endif
#ifndef PER_SPECIES_WEIGHT
        part_weight = current%weight
        fcx = idty * part_weight
        fcy = idtx * part_weight
        fcz = idxy * part_weight
#endif
#ifndef NO_PARTICLE_PROBES
        init_part_x = current%part_pos(1)
        init_part_y = current%part_pos(2)
        init_part_z = current%part_pos(3)
#endif
#ifdef PER_PARTICLE_CHARGE_MASS
        part_q   = current%charge
        part_m   = current%mass
        part_mc  = c * current%mass
        ipart_mc = 1.0_num / part_mc
        cmratio  = part_q * dtfac * ipart_mc
        ccmratio = c * cmratio
#ifndef NO_PARTICLE_PROBES
        part_mc2 = c * part_mc
#endif
#endif
        ! Copy the particle properties out for speed
        part_x  = current%part_pos(1)
        part_y  = current%part_pos(2)
        part_z  = current%part_pos(3)
        part_ux = current%part_p(1) * ipart_mc
        part_uy = current%part_p(2) * ipart_mc
        part_uz = current%part_p(3) * ipart_mc 

        ! Calculate v(t) from p(t)
        ! See PSC manual page (25-27)
        gamma_rel = SQRT(part_ux**2 + part_uy**2 + part_uz**2 + 1.0_num)
        root = dtco2 / gamma_rel

        ! Move particles to half timestep position to first order
        part_x = part_x + part_ux * root
        part_y = part_y + part_uy * root
        part_z = part_z + part_uz * root

        ! Shift the origin of part_x and part_r for interpolation
        part_x_local = part_x - x_grid_min_local 
        part_r = SQRT(part_y**2 + part_z**2)
        part_r_local = part_r - y_grid_min_local
        
        ! exp(-i*theta) for field interpolation
        exp_min_itheta = (part_y - imagi * part_z) / part_r
        theta_05 = ATAN2(part_z, part_y)
        exp_min_itheta_05 = exp_min_itheta
        exp_itheta_05 = 1.0_num / exp_min_itheta_05

#ifdef WORK_DONE_INTEGRATED
        ! This is the actual total work done by the fields: Results correspond
        ! with the electron's gamma factor
        root = cmratio / gamma_rel

        tmp_x = part_ux * root
        tmp_y = part_uy * root
        tmp_z = part_uz * root
#endif

        ! Grid cell position as a fraction.
#ifdef PARTICLE_SHAPE_TOPHAT
        cell_x_r = part_x_local * idx - 0.5_num
        cell_y_r = part_r_local * idy - 0.5_num
#else
        cell_x_r = part_x_local * idx
        cell_y_r = part_r_local * idy
#endif
        ! Round cell position to nearest cell
        cell_x1 = FLOOR(cell_x_r + 0.5_num)
        ! Calculate fraction of cell between nearest cell boundary and particle
        cell_frac_x = REAL(cell_x1, num) - cell_x_r
        cell_x1 = cell_x1 + 1

        cell_y1 = FLOOR(cell_y_r + 0.5_num)
        cell_frac_y = REAL(cell_y1, num) - cell_y_r
        cell_y1 = cell_y1 + 1

        ! Particle weight factors as described in the manual, page25
        ! These weight grid properties onto particles
        ! Also used to weight particle properties onto grid, used later
        ! to calculate J
        ! NOTE: These weights require an additional multiplication factor!
#ifdef PARTICLE_SHAPE_BSPLINE3
#include "bspline3/gx.inc"
#elif  PARTICLE_SHAPE_TOPHAT
#include "tophat/gx.inc"
#else
#include "triangle/gx.inc"
#endif

        ! Now redo shifted by half a cell due to grid stagger.
        ! Use shifted version for ex in X, ey in Y, ez in Z
        ! And in Y&Z for bx, X&Z for by, X&Y for bz
        cell_x2 = FLOOR(cell_x_r)
        cell_frac_x = REAL(cell_x2, num) - cell_x_r + 0.5_num
        cell_x2 = cell_x2 + 1

        cell_y2 = FLOOR(cell_y_r)
        cell_frac_y = REAL(cell_y2, num) - cell_y_r + 0.5_num
        cell_y2 = cell_y2 + 1

        dcellx = 0
        dcelly = 0
        hx = 0.0_num 
        hy = 0.0_num
        ! NOTE: These weights require an additional multiplication factor!
#ifdef PARTICLE_SHAPE_BSPLINE3
#include "bspline3/hx_dcell.inc"
#elif  PARTICLE_SHAPE_TOPHAT
#include "tophat/hx_dcell.inc"
#else
#include "triangle/hx_dcell.inc"
#endif

        ! These are the electric and magnetic fields interpolated to the
        ! particle position. They have been checked and are correct.
        ! Actually checking this is messy.
#ifdef PARTICLE_SHAPE_BSPLINE3
#include "bspline3/e_part.inc"
#include "bspline3/b_part.inc"
#elif  PARTICLE_SHAPE_TOPHAT
#include "tophat/e_part.inc"
#include "tophat/b_part.inc"
#else
#include "triangle/e_part.inc"
#include "triangle/b_part.inc"
#endif

        ! update particle momenta using weighted fields
        uxm = part_ux + cmratio * ex_part
        uym = part_uy + cmratio * ey_part
        uzm = part_uz + cmratio * ez_part

#ifdef HC_PUSH
        ! Half timestep, then use Higuera-Cary push see
        ! https://aip.scitation.org/doi/10.1063/1.4979989
        gamma_rel = uxm**2 + uym**2 + uzm**2 + 1.0_num
        alpha = 0.5_num * part_q * dt / part_m
        beta_x = alpha * bx_part
        beta_y = alpha * by_part
        beta_z = alpha * bz_part
        beta2 = beta_x**2 + beta_y**2 + beta_z**2
        sigma = gamma_rel - beta2
        beta_dot_u = beta_x * uxm + beta_y * uym + beta_z * uzm
        gamma_rel = sigma + SQRT(sigma**2 + 4.0_num * (beta2 + beta_dot_u**2))
        gamma_rel = SQRT(0.5_num * gamma_rel)
#else
        ! Half timestep, then use Boris1970 rotation, see Birdsall and Langdon
        gamma_rel = SQRT(uxm**2 + uym**2 + uzm**2 + 1.0_num)
#endif
        root = ccmratio / gamma_rel

        taux = bx_part * root
        tauy = by_part * root
        tauz = bz_part * root

        taux2 = taux**2
        tauy2 = tauy**2
        tauz2 = tauz**2

        tau = 1.0_num / (1.0_num + taux2 + tauy2 + tauz2)

        uxp = ((1.0_num + taux2 - tauy2 - tauz2) * uxm &
            + 2.0_num * ((taux * tauy + tauz) * uym &
            + (taux * tauz - tauy) * uzm)) * tau
        uyp = ((1.0_num - taux2 + tauy2 - tauz2) * uym &
            + 2.0_num * ((tauy * tauz + taux) * uzm &
            + (tauy * taux - tauz) * uxm)) * tau
        uzp = ((1.0_num - taux2 - tauy2 + tauz2) * uzm &
            + 2.0_num * ((tauz * taux + tauy) * uxm &
            + (tauz * tauy - taux) * uym)) * tau

        ! Rotation over, go to full timestep
        part_ux = uxp + cmratio * ex_part
        part_uy = uyp + cmratio * ey_part
        part_uz = uzp + cmratio * ez_part

        ! Calculate particle velocity from particle momentum
        part_u2 = part_ux**2 + part_uy**2 + part_uz**2
        gamma_rel = SQRT(part_u2 + 1.0_num)
        igamma = 1.0_num / gamma_rel
        root = dtco2 * igamma

        delta_x = part_ux * root
        delta_y = part_uy * root
        delta_z = part_uz * root

        ! Move particles to end of time step at 2nd order accuracy
        part_x = part_x + delta_x
        part_y = part_y + delta_y
        part_z = part_z + delta_z

        ! particle has now finished move to end of timestep, so copy back
        ! into particle array
        current%part_pos = (/ part_x, part_y, part_z /)
        current%part_p   = part_mc * (/ part_ux, part_uy, part_uz /)

        ! Calculate v_theta with the position at time t+dt
        part_vy = part_uy * c * igamma
        part_vz = part_uz * c * igamma
        part_r = SQRT(part_y**2 + part_z**2)
        exp_itheta_10 = (part_y + imagi * part_z) / part_r
        part_vt = -part_vy * AIMAG(exp_itheta_10) &
            + part_vz * REAL(exp_itheta_10,num)

#ifdef WORK_DONE_INTEGRATED
        ! This is the actual total work done by the fields: Results correspond
        ! with the electron's gamma factor
        root = cmratio / gamma_rel

        work_x = ex_part * (tmp_x + part_ux * root)
        work_y = ey_part * (tmp_y + part_uy * root)
        work_z = ez_part * (tmp_z + part_uz * root)

        current%work_x = work_x
        current%work_y = work_y
        current%work_z = work_z
        current%work_x_total = current%work_x_total + work_x
        current%work_y_total = current%work_y_total + work_y
        current%work_z_total = current%work_z_total + work_z
#endif

#ifndef NO_PARTICLE_PROBES
        final_part_x = current%part_pos(1)
        final_part_y = current%part_pos(2)
        final_part_z = current%part_pos(3)
#endif
        ! Original code calculates densities of electrons, ions and neutrals
        ! here. This has been removed to reduce memory footprint

        ! If the code is compiled with zero-current particle support then put in
        ! an IF statement so that the current is not calculated for this species
#ifndef NO_TRACER_PARTICLES
        IF (not_zero_current_species) THEN
#endif
          ! Now advance to t+1.5dt to calculate current. This is detailed in
          ! the manual between pages 37 and 41. The version coded up looks
          ! completely different to that in the manual, but is equivalent.
          ! Use t+1.5 dt so that can update J to t+dt at 2nd order
          part_x_local = part_x + delta_x - x_grid_min_local
          part_y = part_y + delta_y
          part_z = part_z + delta_z
          part_r = SQRT(part_y**2 + part_z**2) 
          part_r_local = part_r - y_grid_min_local

          exp_itheta_15 = (part_y + imagi * part_z) / part_r 
          theta_15 = ATAN2(part_z, part_y)
          exp_idtheta = exp_itheta_15 * exp_min_itheta_05
          dtheta = theta_15 - theta_05

          ! Delta-f calcuation: subtract background from
          ! calculated current.
#ifdef DELTAF_METHOD
          weight_back = current%pvol * f0(ispecies, part_mc / c, current%part_p)
          fcx = idty * (part_weight - weight_back)
          fcy = idtx * (part_weight - weight_back)
          fcz = idxy * (part_weight - weight_back)
#endif

          ! Save old weightings on staggered grid to the g variables
          ! g now represents staggered grid weights at t+0.5dt
          gx = hx
          gy = hy

#ifdef PARTICLE_SHAPE_TOPHAT
          cell_x_r = part_x_local * idx - 0.5_num
          cell_y_r = part_r_local * idy - 0.5_num
#else
          cell_x_r = part_x_local * idx
          cell_y_r = part_r_local * idy
#endif

          cell_x3 = FLOOR(cell_x_r)
          cell_frac_x = REAL(cell_x3, num) - cell_x_r + 0.5_num
          cell_x3 = cell_x3 + 1

          cell_y3 = FLOOR(cell_y_r)
          cell_frac_y = REAL(cell_y3, num) - cell_y_r + 0.5_num
          cell_y3 = cell_y3 + 1

          hx = 0.0_num
          hy = 0.0_num
          dcellx = cell_x3 - cell_x2
          dcelly = cell_y3 - cell_y2
          ! NOTE: These weights require an additional multiplication factor!
#ifdef PARTICLE_SHAPE_BSPLINE3
#include "bspline3/hx_dcell.inc"
#elif  PARTICLE_SHAPE_TOPHAT
#include "tophat/hx_dcell.inc"
#else
#include "triangle/hx_dcell.inc"
#endif

          ! hx, hy currently represent cell weights on the staggered grid at 
          ! t+1.5dt - change these to represent the cell-weight-change between 
          ! t+0.5dt and t+1.5dt
          hx = hx - gx
          hy = hy - gy

          ! Remember that due to CFL condition particle can never cross more
          ! than one gridcell in one timestep

          xmin = sf_min + (dcellx - 1) / 2
          xmax = sf_max + (dcellx + 1) / 2

          ymin = sf_min + (dcelly - 1) / 2
          ymax = sf_max + (dcelly + 1) / 2

          q_weight_fac = part_q * part_weight * fac
          fcx = q_weight_fac * idt
          fcy = q_weight_fac * idt
          fcz = q_weight_fac * part_vt

          DO im = 0, n_mode-1
            
            mdth = REAL(im, num) * dtheta
            m2dth2 = mdth**2
            small_angle_change = ABS(mdth) < 1.0e-4_num
            IF (.NOT. small_angle_change .AND. im > 0) THEN 
              inv_mdth = 1.0_num / mdth 
              inv_m2dth2 = inv_mdth**2
            END IF

            ! Pre-calculate terms which depend on:
            ! ! m * delta_theta
            ! ! exp(i*m*theta)
            ! ! exp(i*m*delta_theta)
            IF (im == 0) THEN
              exp_imtheta0 = 1.0_num
              exp_imdtheta = 1.0_num
            ELSE 
              exp_imtheta0 = exp_imtheta0 * exp_itheta_05
              exp_imdtheta = exp_imdtheta * exp_idtheta

              IF (small_angle_change) THEN 
                m_fac_1 = 2.0_num * exp_imtheta0
                m_fac_2 = m_fac_1 * &
                    (1.0_num + 0.5_num * imagi * mdth - sixth * m2dth2)
                m_fac_3 = m_fac_1 * &
                    (0.5_num + third * imagi * mdth - 0.125_num * m2dth2)  
                m_fac_4 = m_fac_1 * &
                    (third + 0.25_num * imagi * mdth - 0.1_num * m2dth2) 
              ELSE
                m_fac_1 = 2.0_num * inv_mdth * exp_imtheta0 
                m_fac_2 = m_fac_1 * (-imagi * (exp_imdtheta - 1.0_num))
                m_fac_3 = m_fac_1 * (inv_mdth &
                    * (exp_imdtheta * (1.0_num - imagi * mdth) - 1.0_num))
                m_fac_4 = m_fac_1 * (imagi * inv_m2dth2 * (exp_imdtheta &
                    * (-m2dth2 - 2.0_num * imagi * mdth + 2.0_num) - 2.0_num))
              END IF
            END IF  
             
            jyh = 0.0_num
            DO iy = ymin, ymax
              cy = cell_y2 + iy

              ! Pre-calculate terms which vary with grid y
              IF (im == 0) THEN 
                w_rt = gy(iy) + 0.5_num * hy(iy)
                ym_fac_1 = (0.5_num * gy(iy) + third * hy(iy)) 
              ELSE 
                w_rt = m_fac_2 * gy(iy) + m_fac_3 * hy(iy)
                ym_fac_1 = m_fac_3 * gy(iy) + m_fac_4 * hy(iy)
              END IF
                
              fjx = fcx * inv_area_rt(cy)
              fjy = fcy * hy(iy) * inv_area_xt(cy)
              fjz = fcz * inv_volume(cy)

              jxh = 0.0_num
              DO ix = xmin, xmax
                cx = cell_x2 + ix

                IF (im == 0) THEN 
                  w_xt = gx(ix) + 0.5_num * hx(ix) 
                ELSE 
                  w_xt = m_fac_2 * gx(ix) + m_fac_3 * hx(ix)
                END IF
                w_xr = gx(ix) * w_rt + hx(ix) * ym_fac_1

                ! This is the bit that actually solves d(rho)/dt = -div(J)
                jxh = jxh - fjx * hx(ix) * w_rt
                jyh(ix) = jyh(ix) * ratio_area_xt(cy) - fjy * w_xt
                jzh = fjz * w_xr
                jxm(cx+1, cy, im) = jxm(cx+1, cy, im) + jxh
                jrm(cx, cy+1, im) = jrm(cx, cy+1, im) + jyh(ix)
                jtm(cx, cy, im) = jtm(cx, cy, im) + jzh
              END DO
            END DO
          END DO
#ifndef NO_TRACER_PARTICLES
        END IF
#endif
#if !defined(NO_PARTICLE_PROBES) && !defined(NO_IO)
        IF (probes_for_species) THEN
          ! Compare the current particle with the parameters of any probes in
          ! the system. These particles are copied into a separate part of the
          ! output file.

          gamma_rel_m1 = part_u2 / (gamma_rel + 1.0_num)

          current_probe => species_list(ispecies)%attached_probes

          ! Cycle through probes
          DO WHILE(ASSOCIATED(current_probe))
            ! Note that this is the energy of a single REAL particle in the
            ! pseudoparticle, NOT the energy of the pseudoparticle
            probe_energy = gamma_rel_m1 * part_mc2

            ! Unidirectional probe
            IF (probe_energy > current_probe%ek_min) THEN
              IF (probe_energy < current_probe%ek_max) THEN

                d_init  = SUM(current_probe%normal * (current_probe%point &
                    - (/init_part_x, init_part_y, init_part_z/)))
                d_final = SUM(current_probe%normal * (current_probe%point &
                    - (/final_part_x, final_part_y, final_part_z/)))
                IF (d_final < 0.0_num .AND. d_init >= 0.0_num) THEN
                  ! this particle is wanted so copy it to the list associated
                  ! with this probe
                  ALLOCATE(particle_copy)
                  particle_copy = current
                  ! Fraction of step until particle hits probe
                  path_frac = SUM(current_probe%normal*(current_probe%point - &
                      (/init_part_x,init_part_y,init_part_z/))) &
                      /SUM(current_probe%normal * (/final_part_x-init_part_x, &
                      final_part_y-init_part_y, final_part_z-init_part_z/))
                  ! Position of particle on probe
                  particle_copy%part_pos = path_frac * &
                      (/final_part_x-init_part_x, &
                      final_part_y-init_part_y, &
                      final_part_z-init_part_z/) &
                      + (/init_part_x, init_part_y, init_part_z/)
#ifdef PROBE_TIME
                  ! Note: time variable corresponds to (time at x_init)+0.5*dt
                  particle_copy%probe_time = time + dt * (path_frac - 0.5_num)
#endif
                  CALL add_particle_to_partlist(&
                      current_probe%sampled_particles, particle_copy)
                  NULLIFY(particle_copy)
                END IF

              END IF
            END IF
            current_probe => current_probe%next
          END DO
        END IF
#endif
        current => next
      END DO
      CALL current_bcs(species=ispecies)
    END DO

    CALL current_bcs_r_min_final
    CALL particle_bcs

    DEALLOCATE(inv_area_rt, inv_area_xt, inv_volume, ratio_area_xt)

  END SUBROUTINE push_particles



  ! Background distribution function used for delta-f calculations.
  ! Specialise to a drifting (tri)-Maxwellian to simplify and ensure
  ! zero density/current divergence.
  ! Can effectively switch off deltaf method by setting zero background density.

  FUNCTION f0(ispecies, mass, p)

    INTEGER, INTENT(IN) :: ispecies
    REAL(num), INTENT(IN) :: mass
    REAL(num), DIMENSION(:), INTENT(IN) :: p
    REAL(num) :: f0
    REAL(num) :: Tx, Ty, Tz, driftx, drifty, driftz, density
    REAL(num) :: f0_exponent, norm, two_kb_mass, two_pi_kb_mass3
    TYPE(particle_species), POINTER :: species

    species => species_list(ispecies)

    IF (ABS(species%initial_conditions%density_back) > c_tiny) THEN
       two_kb_mass = 2.0_num * kb * mass
       two_pi_kb_mass3 = (pi * two_kb_mass)**3

       Tx = species%initial_conditions%temp_back(1)
       Ty = species%initial_conditions%temp_back(2)
       Tz = species%initial_conditions%temp_back(3)
       driftx  = species%initial_conditions%drift_back(1)
       drifty  = species%initial_conditions%drift_back(2)
       driftz  = species%initial_conditions%drift_back(3)
       density = species%initial_conditions%density_back
       f0_exponent = ((p(1) - driftx)**2 / Tx &
                    + (p(2) - drifty)**2 / Ty &
                    + (p(3) - driftz)**2 / Tz) / two_kb_mass
       norm = density / SQRT(two_pi_kb_mass3 * Tx * Ty * Tz)
       f0 = norm * EXP(-f0_exponent)
    ELSE
       f0 = 0.0_num
    END IF

  END FUNCTION f0



#if defined(PHOTONS) || defined(BREMSSTRAHLUNG)
  SUBROUTINE push_photons(ispecies)

    ! Very simple photon pusher
    ! Properties of the current particle. Copy out of particle arrays for speed
    REAL(num) :: delta_x, delta_y, delta_z
    INTEGER,INTENT(IN) :: ispecies
    TYPE(particle), POINTER :: current
    REAL(num) :: current_energy, dtfac, fac

    ! Used for particle probes (to see of probe conditions are satisfied)
#ifndef NO_PARTICLE_PROBES
    REAL(num) :: init_part_x, final_part_x
    REAL(num) :: init_part_y, final_part_y
    REAL(num) :: init_part_z, final_part_z
    TYPE(particle_probe), POINTER :: current_probe
    TYPE(particle), POINTER :: particle_copy
    REAL(num) :: d_init, d_final
    LOGICAL :: probes_for_species
    REAL(num) :: path_frac
#endif

#ifndef NO_PARTICLE_PROBES
    current_probe => species_list(ispecies)%attached_probes
    probes_for_species = ASSOCIATED(current_probe)
#if defined(PARTICLE_ID) || defined(PARTICLE_ID4) 
    IF (probes_for_species) THEN 
      CALL generate_particle_ids(species_list(ispecies)%attached_list)
    END IF
#endif
#endif
    dtfac = dt * c**2

    ! set current to point to head of list
    current => species_list(ispecies)%attached_list%head
    ! loop over photons
    DO WHILE(ASSOCIATED(current))
      ! Note that this is the energy of a single REAL particle in the
      ! pseudoparticle, NOT the energy of the pseudoparticle
      current_energy = current%particle_energy

      fac = dtfac / current_energy
      delta_x = current%part_p(1) * fac
      delta_y = current%part_p(2) * fac
      delta_z = current%part_p(3) * fac
#ifndef NO_PARTICLE_PROBES
      init_part_x = current%part_pos(1)
      init_part_y = current%part_pos(2)
      init_part_z = current%part_pos(3)
#endif
      current%part_pos = current%part_pos + (/delta_x, delta_y, delta_z/)
#ifndef NO_PARTICLE_PROBES
      final_part_x = current%part_pos(1)
      final_part_y = current%part_pos(2)
      final_part_z = current%part_pos(3)
#endif

#ifndef NO_PARTICLE_PROBES
      IF (probes_for_species) THEN
        ! Compare the current particle with the parameters of any probes in
        ! the system. These particles are copied into a separate part of the
        ! output file.

        current_probe => species_list(ispecies)%attached_probes

        ! Cycle through probes
        DO WHILE(ASSOCIATED(current_probe))
          ! Unidirectional probe
          IF (current_energy > current_probe%ek_min) THEN
            IF (current_energy < current_probe%ek_max) THEN

              d_init  = SUM(current_probe%normal * (current_probe%point &
                  - (/init_part_x, init_part_y, init_part_z/)))
              d_final = SUM(current_probe%normal * (current_probe%point &
                  - (/final_part_x, final_part_y, final_part_z/)))
              IF (d_final < 0.0_num .AND. d_init >= 0.0_num) THEN
                ! this particle is wanted so copy it to the list associated
                ! with this probe
                ALLOCATE(particle_copy)
                particle_copy = current
                ! Fraction of step until particle hits probe
                path_frac = SUM(current_probe%normal*(current_probe%point - &
                    (/init_part_x,init_part_y,init_part_z/))) &
                    /SUM(current_probe%normal * (/final_part_x-init_part_x, &
                    final_part_y-init_part_y, &
                    final_part_z-init_part_z/))
                ! Position of particle on probe
                particle_copy%part_pos = path_frac * &
                    (/final_part_x-init_part_x, &
                    final_part_y-init_part_y, &
                    final_part_z-init_part_z/) &
                    + (/init_part_x,init_part_y,init_part_z/)
#ifdef PROBE_TIME
                ! Note: time variable corresponds to (time at x_init)+0.5*dt
                particle_copy%probe_time = time + dt * (path_frac - 0.5_num)
#endif
                CALL add_particle_to_partlist(&
                    current_probe%sampled_particles, particle_copy)
                NULLIFY(particle_copy)
              END IF

            END IF
          END IF
          current_probe => current_probe%next
        END DO
      END IF
#endif

      current => current%next
    END DO

  END SUBROUTINE push_photons
#endif



SUBROUTINE no_field_push(ispecies)

  ! Particle pusher using momentum only.
  ! Does not contribute to fields or current density
  REAL(num) :: delta_x, delta_y, delta_z
  REAL(num) :: part_ux, part_uy, part_uz
  REAL(num) :: part_x, part_y, part_z
  INTEGER,INTENT(IN) :: ispecies
  TYPE(particle), POINTER :: current
  REAL(num) :: cdt, ipart_mc, gamma_rel, root

  ! Used for particle probes (to see of probe conditions are satisfied)
#ifndef NO_PARTICLE_PROBES
  REAL(num) :: init_part_x, final_part_x
  REAL(num) :: init_part_y, final_part_y
  REAL(num) :: init_part_z, final_part_z
  TYPE(particle_probe), POINTER :: current_probe
  TYPE(particle), POINTER :: particle_copy
  REAL(num) :: d_init, d_final, current_energy
  LOGICAL :: probes_for_species
  REAL(num) :: path_frac
#endif

  cdt = c * dt

#ifndef NO_PARTICLE_PROBES
  current_probe => species_list(ispecies)%attached_probes
  probes_for_species = ASSOCIATED(current_probe)
#if defined(PARTICLE_ID) || defined(PARTICLE_ID4) 
  IF (probes_for_species) THEN 
    CALL generate_particle_ids(species_list(ispecies)%attached_list)
  END IF
#endif
#endif

#ifndef PER_PARTICLE_CHARGE_MASS
  ipart_mc = 1.0_num / (species_list(ispecies)%mass * c)
#endif

  ! set current to point to head of list
  current => species_list(ispecies)%attached_list%head
  ! loop over photons
  DO WHILE(ASSOCIATED(current))
#ifdef PER_PARTICLE_CHARGE_MASS
    ipart_mc = 1.0_num / (current%mass * c)
#endif
#ifndef NO_PARTICLE_PROBES
    init_part_x = current%part_pos(1)
    init_part_y = current%part_pos(2)
    init_part_z = current%part_pos(3)
#endif

    part_ux = current%part_p(1) * ipart_mc
    part_uy = current%part_p(2) * ipart_mc
    part_uz = current%part_p(3) * ipart_mc 

    ! Calculate v(t) from p(t)
    ! See PSC manual page (25-27)
    gamma_rel = SQRT(part_ux**2 + part_uy**2 + part_uz**2 + 1.0_num)
    root = cdt / gamma_rel

    ! Move particles by a full time-step
    part_x = part_x + part_ux * root
    part_y = part_y + part_uy * root
    part_z = part_z + part_uz * root

    current%part_pos = current%part_pos + (/delta_x, delta_y, delta_z/)
#ifndef NO_PARTICLE_PROBES
    final_part_x = current%part_pos(1)
    final_part_y = current%part_pos(2)
    final_part_z = current%part_pos(3)
#endif

#ifndef NO_PARTICLE_PROBES
    IF (probes_for_species) THEN
      ! Compare the current particle with the parameters of any probes in
      ! the system. These particles are copied into a separate part of the
      ! output file.

      current_probe => species_list(ispecies)%attached_probes
      current_energy = (gamma_rel - 1.0_num) * c / ipart_mc

      ! Cycle through probes
      DO WHILE(ASSOCIATED(current_probe))
        ! Unidirectional probe
        IF (current_energy > current_probe%ek_min) THEN
          IF (current_energy < current_probe%ek_max) THEN

            d_init  = SUM(current_probe%normal * (current_probe%point &
                - (/init_part_x, init_part_y, init_part_z/)))
            d_final = SUM(current_probe%normal * (current_probe%point &
                - (/final_part_x, final_part_y, final_part_z/)))
            IF (d_final < 0.0_num .AND. d_init >= 0.0_num) THEN
              ! this particle is wanted so copy it to the list associated
              ! with this probe
              ALLOCATE(particle_copy)
              particle_copy = current
              ! Fraction of step until particle hits probe
              path_frac = SUM(current_probe%normal*(current_probe%point - &
                  (/init_part_x,init_part_y,init_part_z/))) &
                  /SUM(current_probe%normal * (/final_part_x-init_part_x, &
                  final_part_y-init_part_y, &
                  final_part_z-init_part_z/))
              ! Position of particle on probe
              particle_copy%part_pos = path_frac * &
                  (/final_part_x-init_part_x, &
                  final_part_y-init_part_y, &
                  final_part_z-init_part_z/) &
                  + (/init_part_x,init_part_y,init_part_z/)
#ifdef PROBE_TIME
              ! Note: time variable corresponds to (time at x_init)+0.5*dt
              particle_copy%probe_time = time + dt * (path_frac - 0.5_num)
#endif
              CALL add_particle_to_partlist(&
                  current_probe%sampled_particles, particle_copy)
              NULLIFY(particle_copy)
            END IF

          END IF
        END IF
        current_probe => current_probe%next
      END DO
    END IF
#endif

    current => current%next
  END DO

END SUBROUTINE no_field_push

END MODULE particles

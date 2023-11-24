exp_min_imtheta = 1.0_num
part_var = 0.0_num

DO im = 0, n_mode-1
  part_var = part_var + REAL( exp_min_imtheta * (&
        gy(-1) * (gx(-1) * grid_var(cell_x1-1,cell_y1-1,im) &
      +           gx( 0) * grid_var(cell_x1  ,cell_y1-1,im) &
      +           gx( 1) * grid_var(cell_x1+1,cell_y1-1,im)) &
      + gy( 0) * (gx(-1) * grid_var(cell_x1-1,cell_y1  ,im) &
      +           gx( 0) * grid_var(cell_x1  ,cell_y1  ,im) &
      +           gx( 1) * grid_var(cell_x1+1,cell_y1  ,im)) &
      + gy( 1) * (gx(-1) * grid_var(cell_x1-1,cell_y1+1,im) &
      +           gx( 0) * grid_var(cell_x1  ,cell_y1+1,im) &
      +           gx( 1) * grid_var(cell_x1+1,cell_y1+1,im))), num)

  exp_min_imtheta = exp_min_imtheta * exp_min_itheta
END DO

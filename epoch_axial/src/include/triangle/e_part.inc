exp_min_imtheta = 1.0_num
ex_part = 0.0_num 
er_part = 0.0_num 
et_part = 0.0_num

DO im = 0, n_mode-1
  ex_part = ex_part + REAL( exp_min_imtheta * (&
        hy(-1) * (gx(-1) * exm(cell_x1-1,cell_y2-1,im) &
      +           gx( 0) * exm(cell_x1  ,cell_y2-1,im) &
      +           gx( 1) * exm(cell_x1+1,cell_y2-1,im)) &
      + hy( 0) * (gx(-1) * exm(cell_x1-1,cell_y2  ,im) &
      +           gx( 0) * exm(cell_x1  ,cell_y2  ,im) &
      +           gx( 1) * exm(cell_x1+1,cell_y2  ,im)) &
      + hy( 1) * (gx(-1) * exm(cell_x1-1,cell_y2+1,im) &
      +           gx( 0) * exm(cell_x1  ,cell_y2+1,im) &
      +           gx( 1) * exm(cell_x1+1,cell_y2+1,im))), num)

  er_part = er_part + REAL( exp_min_imtheta * (&
        gy(-1) * (hx(-1) * erm(cell_x2-1,cell_y1-1,im) &
      +           hx( 0) * erm(cell_x2  ,cell_y1-1,im) &
      +           hx( 1) * erm(cell_x2+1,cell_y1-1,im)) &
      + gy( 0) * (hx(-1) * erm(cell_x2-1,cell_y1  ,im) &
      +           hx( 0) * erm(cell_x2  ,cell_y1  ,im) &
      +           hx( 1) * erm(cell_x2+1,cell_y1  ,im)) &
      + gy( 1) * (hx(-1) * erm(cell_x2-1,cell_y1+1,im) &
      +           hx( 0) * erm(cell_x2  ,cell_y1+1,im) &
      +           hx( 1) * erm(cell_x2+1,cell_y1+1,im))), num)

  et_part = et_part + REAL( exp_min_imtheta * (&
        hy(-1) * (hx(-1) * etm(cell_x2-1,cell_y2-1,im) &
      +           hx( 0) * etm(cell_x2  ,cell_y2-1,im) &
      +           hx( 1) * etm(cell_x2+1,cell_y2-1,im)) &
      + hy( 0) * (hx(-1) * etm(cell_x2-1,cell_y2  ,im) &
      +           hx( 0) * etm(cell_x2  ,cell_y2  ,im) &
      +           hx( 1) * etm(cell_x2+1,cell_y2  ,im)) &
      + hy( 1) * (hx(-1) * etm(cell_x2-1,cell_y2+1,im) &
      +           hx( 0) * etm(cell_x2  ,cell_y2+1,im) &
      +           hx( 1) * etm(cell_x2+1,cell_y2+1,im))), num)

  exp_min_imtheta = exp_min_imtheta * exp_min_itheta
END DO

ey_part =  er_part * REAL(exp_min_itheta) + et_part * AIMAG(exp_min_itheta)
ez_part = -er_part * AIMAG(exp_min_itheta) + et_part * REAL(exp_min_itheta)
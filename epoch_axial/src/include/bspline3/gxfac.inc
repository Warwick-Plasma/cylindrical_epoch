cx2 = cell_frac_x**2
gx(-2) = fac1 * (0.5_num + cell_frac_x)**4
gx(-1) = fac2 * (1.1875_num + 2.75_num * cell_frac_x &
    + cx2 * (1.5_num - cell_frac_x - cx2))
gx( 0) = 0.25_num * (fac3 + cx2 * (cx2 - 2.5_num))
gx( 1) = fac2 * (1.1875_num - 2.75_num * cell_frac_x &
    + cx2 * (1.5_num + cell_frac_x - cx2))
gx( 2) = fac1 * (0.5_num - cell_frac_x)**4

cy2 = cell_frac_y**2
gy(-2) = fac1 * (0.5_num + cell_frac_y)**4
gy(-1) = fac2 * (1.1875_num + 2.75_num * cell_frac_y &
    + cy2 * (1.5_num - cell_frac_y - cy2))
gy( 0) = 0.25_num * (fac3 + cy2 * (cy2 - 2.5_num))
gy( 1) = fac2 * (1.1875_num - 2.75_num * cell_frac_y &
    + cy2 * (1.5_num + cell_frac_y - cy2))
gy( 2) = fac1 * (0.5_num - cell_frac_y)**4

IF (part_r < 2.0_num * dy) THEN
  IF (part_r < dy) THEN 
    gy(0) = gy(0) + gy(-1) 
    gy(1) = gy(1) + gy(-2)
    gy(-1) = 0.0_num 
    gy(-2) = 0.0_num
  ELSE 
    gy(-1) = gy(-1) + gy(-2) 
    gy(-2) = 0.0_num
  END IF
END IF
   
 
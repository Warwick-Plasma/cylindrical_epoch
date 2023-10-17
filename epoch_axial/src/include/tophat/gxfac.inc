gx( 0) = 0.5_num + cell_frac_x
gx( 1) = 0.5_num - cell_frac_x

gy( 0) = 0.5_num + cell_frac_y
gy( 1) = 0.5_num - cell_frac_y

IF (part_r < 0.5_num * dy) THEN
  gy(1) = 1.0_num
  gy(0) = 0.0_num
END IF
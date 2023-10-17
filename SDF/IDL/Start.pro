;
; SDF (Self-Describing Format) IDL reader
; Copyright (c) 2008-2016, SDF Development Team
;
; Distributed under the terms of the BSD 3-clause License.
; See the LICENSE file for details.
;


.r IsoPlot.pro
.r LoadSDF.pro
.r SDFHelp.pro
.r StartPIC.pro
.r widget.pro
.r ReadNameVal.pro

init_StartPIC

q0 = 1.602176565d-19 ; elementary charge [C]
m0 = 9.10938291d-31  ; electron mass [kg]
v0 = 2.99792458d8    ; speed of light [m/s]
kb = 1.3806488d-23   ; Boltzmann's constant [J/K]
mu0 = 4.0d-7 * !dpi  ; Magnetic constant [N/A^2]
epsilon0 = 8.8541878176203899d-12 ; Electric constant [F/m]
h_planck = 6.62606957d-34 ; Planck constant [J s]

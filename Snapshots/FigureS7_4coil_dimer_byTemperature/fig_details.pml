bg_color white

color black, all
show cell

alter resn ILE, ss='H'
alter resn ALA, ss='H'
alter resn LYS, ss='H'
alter resn GLU, ss='H'

set ray_trace_mode=1

set cartoon_cylindrical_helices, 1
as cartoon

cartoon tube, resname GLY


set_color olive, [0.809, 0.719, 0.484]
color olive, resn ILE
color olive, resn ALA
color olive, resn LYS
color olive, resn GLU

color black, resn GLY

set antialias,2


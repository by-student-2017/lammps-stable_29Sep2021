# small W unit cell for use with VASP

variable        mode index file
 
if "${mode} == file" then &
  "message client md file tmp.couple" &
elif "${mode} == zmq" &
  "message client md zmq localhost:5555" &

variable	x index 1
variable	y index 1
variable	z index 1

units		metal
atom_style	atomic
atom_modify     sort 0 0.0 map yes

read_data       data.W
mass		1 183.85

replicate       $x $y $z

velocity	all create 300.0 87287 loop geom

#neighbor        0.3 bin # 0.3 is lj unit. 2.0 [Angstrom] is metal or real units
#neigh_modify    delay 0 every 10 check no
comm_modify cutoff 2.0 # = comm_modify mode single group all cutoff 2.0 vel yes

fix		1 all nve
fix             2 all client/md
fix_modify      2 energy yes

thermo          1
run		3

message         quit

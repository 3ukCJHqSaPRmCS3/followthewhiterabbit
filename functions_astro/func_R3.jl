###################################################################
# function func_eci2ecef(ECI, PropagationTimeSec, t_start, dt_sec)
###################################################################
#
# Purpose:
# Creates a rotation matrix about the 3rd-axis (or the Z-axis)
#
# INPUT:
# x              = Rotation angle      [1x1] (rad)
#
# OUTPUT:
# A              =  rotation matrix about the 3rd-axis (or the Z-axis)      [3x3] (-)
#
###################################################################
# Written by:
# Carlos Deccia
# July 12th, 2022
###################################################################

function func_R3(x)

cx = cos(x)
sx = sin(x)
A = [cx sx 0; -sx cx 0; 0 0 1]
return A

 end

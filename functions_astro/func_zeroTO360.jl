###################################################################
# function func_zeroTO360(x)
###################################################################
#
# Purpose: This function reduces an angle to the range of 0 - 360 degrees or 0 - 2*pi radians.
#
# INPUT:
# z              = Angle to be reduced, may be an array of angles
#
# OUTPUT:
# y              = Reduced angle
#
###################################################################
# Written by:
# Carlos Deccia
# July 12th, 2022
###################################################################

function func_zeroTO360(x)

degg = 2*pi

if (x >= 2*pi)
     x = x - trunc(x/degg)*degg
elseif (x < 0)
     x = x - (trunc(x/degg) - 1)*degg
end
y = x

return y

end

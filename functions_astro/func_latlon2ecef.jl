###################################################################
# function func_latlon2ecef(Lat_rad, Lon_rad, R, Alt)
###################################################################
#
# Purpose:
# Creates a rotation matrix about the 3rd-axis (or the Z-axis)
#
# INPUT:
# Lat_deg           = Latitude   [1x1] (deg) ############################################################# NOTE GEODETIC or GEOCENTRI ??? I think geocentric. BUT CHECK
# Lon_deg           = Longitude  [1x1] (deg) ############################################################# NOTE GEODETIC or GEOCENTRI ??? MAYBE VALLADO has code to check ?
# Rearth_km         = Radius of Earth                 [1x1] (km)
# Alt_km            = Altitude                        [1x1] (km)
#
# OUTPUT:
# ECEF              = ECEF position                   [3x1] (km)
#
###################################################################
# Written by:
# Carlos Deccia
# July 13th, 2022
###################################################################
function func_latlon2ecef(Lat_deg,Lon_deg,R,Alt)

R_ECEF_tmp = zeros(3)

R_ECEF_tmp[1] = (R+Alt)*cosd(Lat_deg)*cosd(Lon_deg)
R_ECEF_tmp[2] = (R+Alt)*cosd(Lat_deg)*sind(Lon_deg)
R_ECEF_tmp[3] = (R+Alt)*sind(Lat_deg)

return R_ECEF_tmp

end

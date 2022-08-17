###################################################################
# function func_eci2ecef(ECI, PropagationTimeSec, t_start, dt_sec)
###################################################################
#
# INPUT:
# RV_ECI              = ECI [position velocity] vector      [6x1] (km;km/s)
# PropagationTimeSec  = Time of propagation      [1x1] (s)
# dt_sec              = Time step of propagation [1x1] (s)
#
# OUTPUT:
# RV_ECI              = ECEF [position velocity] vector      [6x1] (km;km/s)
#
###################################################################
# Written by:
# Carlos Deccia
# July 12th, 2022
###################################################################

function func_eci2ecef(ECI, PropagationTimeSec, t_start, dt_sec, w_earth_rad_s)

# GMSTo = 1.7494710664781898   # computed with http://neoprogrammics.com/sidereal_time_calculator/index.php using 2003-01-01-00:00 as input, but it shouldn't really matter because the computation should apply for any dates
GMSTo    = 0.0
time_vec = collect(t_start:dt_sec:PropagationTimeSec)
col_len  = length(time_vec)
GMST_rad = zeros(col_len)

ECEF     = [zeros(6) for _ in 1:length(ECI)]

for jj = 1:1:col_len
    GMST_rad[jj]  = func_zeroTO360(GMSTo + w_earth_rad_s*time_vec[jj]);  # GMST in radians
    ECEF[jj][1:3] = func_R3(GMST_rad[jj])*ECI[jj][1:3]
    ECEF[jj][4:6] = func_R3(GMST_rad[jj])*ECI[jj][4:6]
end

return ECEF

end

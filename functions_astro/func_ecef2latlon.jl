###################################################################
# function func_ecef2latlon(ECEF, Rearth_km)
###################################################################
#
# Purpose:
# Creates a rotation matrix about the 3rd-axis (or the Z-axis)
#
# INPUT:
# ECEF              = ECEF position and velocity      [6x1] (km km/s)
# Rearth_km         = Radius of Earth                 [1x1] (km)
#
# OUTPUT:
# lat_geod_deg      =
# lon_geod_deg      =
# lat_geoc_deg      =
# long_geoc_deg     =
#
###################################################################
# Written by:
# Carlos Deccia
# July 12th, 2022
###################################################################

function func_ecef2latlon(ECEF, Rearth_km)

lat_geod_deg = zeros(length(ECEF))
lon_geod_deg = zeros(length(ECEF))
lat_geoc_deg = zeros(length(ECEF))
lon_geoc_deg = zeros(length(ECEF))


for i=1:length(ECEF)

    # Rearth_km = 6378.137 #km
    eEarth = 0.081819221456 ######################### NOTE Earth eccentricity ??? probably the (f1-f2)/a equation {provide SOURCE !!!!}

    r_dsat = sqrt(ECEF[i][1]^2+ECEF[i][2]^2)

    sinA = ECEF[i][2]/r_dsat
    cosA = ECEF[i][1]/r_dsat
    lambda = atand(sinA,cosA) # LONG

    delta = atand(ECEF[i][3]/r_dsat)

    global psi_old = delta
    tol = 0.0001

    CEarth = Rearth_km/sqrt(1-eEarth^2*sind(psi_old)^2)
    global psi = atand( (ECEF[i][3] + CEarth*eEarth^2*sind(psi_old) )/r_dsat )

    while abs(psi - psi_old)>tol
        global psi_old = psi
        CEarth = Rearth_km / sqrt( 1-eEarth^2*sind(psi)^2 )
        global psi = atand( (ECEF[i][3] + CEarth*eEarth^2*sind(psi) )/r_dsat )
    end

    lat_geod_deg[i] = psi               # geoDetic LAT [deg]
    lon_geod_deg[i] = lambda            # geoDetic LON [deg]

    lat_geoc_deg[i] = atand((1-eEarth^2)*tand(psi))              # geoCentric LAT [deg]
    lon_geoc_deg[i] = lambda                                     # geoCentric LON [deg]

end

return lat_geod_deg, lon_geod_deg, lat_geoc_deg, lon_geoc_deg

end

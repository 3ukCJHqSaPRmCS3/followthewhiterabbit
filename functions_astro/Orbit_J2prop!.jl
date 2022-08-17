###################################################################
# function Orbit_J2prop!(dstate, state, params, t)
###################################################################
#
# INPUT:
# dstate = [rdot vdot] [3X1 3X1]
# state = [r v]        [3X1 3X1]
# params = [mu_km3s2 J2 Rearth_km] [1X1 1X1 1X1]
#
# OUTPUT:
# none
#
###################################################################
# Written by:
# Carlos Deccia
# July 12th, 2022
###################################################################

function Orbit_J2prop!(dstate, state, param, t)
                    # (du,      u,     p,     t)
    #For use with "DifferentialEquations.jl"

    mu     = param[1]     # 398600.4415 #[km/s^3]
    J2     = param[2]     # 1.08264E-3 [-]
    Rearth = param[3]     # 6378.1363 #[km]

    x = state[1]
    y = state[2]
    z = state[3]

    rmag = sqrt(x^2 + y^2 + z^2)

    dstate[1] = state[4]
    dstate[2] = state[5]
    dstate[3] = state[6]
    dstate[4] = -mu*x/rmag^3 -J2*mu*Rearth^2*( (-15*z^2*x)/(2*rmag^7) + (3*x)/(2*rmag^5) )
    dstate[5] = -mu*y/rmag^3 -J2*mu*Rearth^2*( (-15*z^2*y)/(2*rmag^7) + (3*y)/(2*rmag^5) )
    dstate[6] = -mu*z/rmag^3 -J2*mu*Rearth^2*( (-15*z^3)/(2*rmag^7)   + (9*z)/(2*rmag^5) )

end

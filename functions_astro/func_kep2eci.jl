###################################################################
# function func_eci2ecef(ECI, PropagationTimeSec, t_start, dt_sec)
###################################################################
#
# INPUT:
# mu              = G*M                     [1x1] (km^3/s^2)
# a               = Semi-major axis         [1x1] (km)
# e               = Eccentricity            [1x1] (-)
# i               = Inclination             [1x1] (deg)
# W               = RAAN                    [1x1] (deg)
# w               = Argument of pericenter  [1x1] (deg)
# M               = True anomaly            [1x1] (deg)
#
# OUTPUT:
# RV_ECI              = ECEF [position velocity] vector      [6x1] (km;km/s)
#
###################################################################
# Written by:
# Carlos Deccia
# July 12th, 2022
###################################################################

function func_kep2eci(mu,a,e,i,W,w,M)
    # mu = 398600.4418;
    maxiter=100
    count=0
    tol=1E-9

    #Convert angles into radians
    i=i*pi/180
    W=W*pi/180
    w=w*pi/180
    M=M*pi/180

    #Calculate n
    # a=(mu/n^2)^(1/3);    # T=(2*pi)/n
    n = sqrt(mu/a^3)

    #Calculate Initial Position
        To=M/n
        E=M
        delta=(M-E+e*sin(E))/(1-e*cos(E))

        while abs(delta) > tol && count < maxiter
            count = count + 1
            delta=(M-E+e*sin(E))/(1-e*cos(E))
            E=E+delta
        end

         v=2*atan((sqrt((1+e)/(1-e))*tan(E/2)));

             # R,V = Kep2ECI(a,e,i,W,w,v);
             # mu = 398600.4418;
             p=a*(1-e^2)

             # Position Vector
             rx=p*cos(v)/(1+e*cos(v))
             ry=p*sin(v)/(1+e*cos(v))
             rz=0
             Rp=[rx;ry;rz]

             # Velocity Vector
             vx = -sqrt(mu/p)*sin(v)
             vy = sqrt(mu/p)*(e+cos(v))
             vz=0
             Vp=[vx;vy;vz]

             R=[cos(W) -sin(W) 0;sin(W) cos(W) 0;0 0 1]*[1 0 0;0 cos(i) -sin(i); 0 sin(i) cos(i)]*[cos(w) -sin(w) 0; sin(w) cos(w) 0; 0 0 1]*Rp
             V=[cos(W) -sin(W) 0;sin(W) cos(W) 0;0 0 1]*[1 0 0;0 cos(i) -sin(i); 0 sin(i) cos(i)]*[cos(w) -sin(w) 0; sin(w) cos(w) 0; 0 0 1]*Vp
             # R,V = Kep2ECI(a,e,i,W,w,v)

    return R, V
end

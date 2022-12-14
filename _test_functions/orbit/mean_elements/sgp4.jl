# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==============================================================================
#
#   Functions to convert osculating elements to mean elements using SGP4.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
# ==============================================================================
#
#   [1] Vallado, D. A; Crawford, P (2008). SGP4 orbit determination. AIAA/AAS
#       Astrodynamics Specialist Conference, Honoulu, HI.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

export rv_to_mean_elements_sgp4, rv_to_tle

"""
    rv_to_mean_elements_sgp4(vJD::AbstractVector{T}, vr::AbstractVector{Tv}, vv::AbstractVector{Tv}, W = I; estimate_bstar::Bool = true, mean_elements_epoch::Symbol = :end, max_it::Int = 50, sgp4_gc = sgp4_gc_wgs84, atol::Number = 2e-4, rtol::Number = 2e-4) where {T,Tv<:AbstractVector}

Compute the mean elements for SGP4 based on the position `vr` and velocity
vectors `vr` represented in TEME reference frame. The epoch of those
measurements [Julian Day] must be in `vJD`.

The matrix `W` defined the weights for the least-square algorithm.

# Keywords

* `estimate_bstar`: If `true`, then the BSTAR parameters of the TLE will be
                    estimated.
* `mean_elements_epoch`: If it is  `:end`, the epoch of the mean elements will
                         be equal to the last value in `vJD`. Otherwise, if it
                         is `:begin`, the epoch will be selected as the first
                         value in `vJD`.
* `max_it`: The maximum allowed number of iterations.
* `sgp4_gc`: SPG4 constants (see `SGP4_GravCte`).
* `atol`: The tolerance for the absolute value of the residue. If, at any
          iteration, the residue is lower than `atol`, then the iterations stop.
* `rtol`: The tolerance for the relative difference between the residues. If, at
          any iteration, the relative difference between the residues in two
          consecutive iterations is lower than `rtol`, then the iterations stop.

# Returns

* The epoch of the elements [Julian Day].
* The mean elements for SGP4 algorithm:
    * Semi-major axis [m];
    * Eccentricity [ ];
    * Inclination [rad];
    * Right ascension of the ascending node [rad];
    * Argument of perigee [rad];
    * True anomaly [rad];
    * BSTAR (0 if `estimate_bstar` is `false`).
* The covariance matrix of the mean elements estimation.

"""
function rv_to_mean_elements_sgp4(vJD::AbstractVector{T},
                                  vr::AbstractVector{Tv},
                                  vv::AbstractVector{Tv},
                                  W = I;
                                  estimate_bstar::Bool = true,
                                  mean_elements_epoch::Symbol = :end,
                                  max_it::Int = 50,
                                  sgp4_gc = sgp4_gc_wgs84,
                                  atol::Number = 2e-4,
                                  rtol::Number = 2e-4) where
    {T,Tv<:AbstractVector}

    # Number of measurements.
    num_meas = length(vr)

    # Check if the orbit epoch must be the first or the last element.
    if mean_elements_epoch == :end
        r??? = last(vr)
        v??? = last(vv)
        epoch = last(vJD)
    else
        r??? = first(vr)
        v??? = first(vv)
        epoch = first(vJD)
    end

    # Initial guess of the mean elements.
    #
    # NOTE: x??? is the previous estimate and x??? is the current estimate.
    x??? = estimate_bstar ?
        SVector{7,T}(r???[1], r???[2], r???[3], v???[1], v???[2], v???[3], T(0.00001)) :
        SVector{7,T}(r???[1], r???[2], r???[3], v???[1], v???[2], v???[3], T(0))
    x??? = x???

    # Number of states in the input vector.
    num_states = length(x???)

    # Covariance matrix.
    P = SMatrix{num_states, num_states, T}(I)

    # Variable to store the last residue.
    ??_i?????? = T(NaN)

    # Variable to store how many iterations the residue increased. This is used
    # to account for divergence.
    ??d = 0

    # Header
    @printf("          %10s %20s %20s\n", "Iter.", "Residue", "Res. variation")

    # Loop until the maximum allowed iteration.
    @inbounds for it = 1:max_it
        x??? = x???

        # Variables to store the summations to compute the least square fitting
        # algorithm.
        ??A???WA = @SMatrix zeros(num_states, num_states)
        ??A???Wb = @SVector zeros(num_states)

        # Variable to store the residue in this iteration.
        ??_i = T(0)

        @views for k = 1:num_meas
            # Obtain the measured ephemerides.
            y = vcat(vr[k], vv[k])

            # Obtain the computed ephemerides considering the current estimate
            # of the mean elements.
            ??t = (vJD[k] - epoch)*86400

            r??, v??, ~ = estimate_bstar ?
                _sgp4_sv(??t, sgp4_gc, epoch, x???[1], x???[2], x???[3], x???[4], x???[5], x???[6], x???[7]) :
                _sgp4_sv(??t, sgp4_gc, epoch, x???[1], x???[2], x???[3], x???[4], x???[5], x???[6])

            y?? = vcat(r??, v??)

            # Compute the error.
            b = y - y??

            # Compute the Jacobian.
            A = _sgp4_jacobian(??t, epoch, x???, y??; estimate_bstar = estimate_bstar)

            # Accumulation.
            ??A???WA += A'*W*A
            ??A???Wb += A'*W*b
            ??_i   += sum(W*(b.^2))
        end

        # Normalize the residue.
        ??_i /= num_meas

        # Update the estimate.
        P  = SMatrix{num_states, num_states, T}(pinv(??A???WA))
        ??x = P*??A???Wb

        # Limit the correction to avoid divergence.
        for i = 1:num_states
            abs(??x[i] / x???[i]) > 0.01 &&
                (??x = setindex(??x, 0.01 * abs(x???[i]) * sign(??x[i]), i))
        end

        x??? = x??? + ??x

        # Compute the residue variation.
        ??_p = (??_i - ??_i??????)/??_i??????

        if it != 1
            @printf("PROGRESS: %10d %20g %20g %%\n", it, ??_i, 100??_p)
        else
            @printf("PROGRESS: %10d %20g %20s\n", it, ??_i, "---")
        end

        # Check if the residue is increasing.
        if ??_i < ??_i??????
            ??d = 0
        else
            ??d += 1
        end

        # If the residue increased by three iterations and the residue is higher
        # than 5e8, then we abort because the iterations are diverging.
        ( (??d ??? 3) && (??_i > 5e11) ) && error("The iterations diverged!")

        # Check if the condition to stop has been reached.
        ((abs(??_p) < rtol) || (??_i < atol) || (it ??? max_it)) && break

        ??_i?????? = ??_i
    end

    # Convert the state vector to mean elements.
    orb = rv_to_kepler(x???[1], x???[2], x???[3], x???[4], x???[5], x???[6])
    bstar = x???[7]

    # Assemble the output vector.
    xo = @SVector [orb.a, orb.e, orb.i, orb.??, orb.??, orb.f, bstar]

    # Return the mean elements for SGP4 and the covariance matrix.
    return epoch, xo, P
end

"""
    rv_to_tle(args...; name::String = "UNDEFINED", sat_num::Int = 9999, classification::Char = 'U', int_designator = "999999", elem_set_number::Int = 0, rev_num, kwargs...)

Convert a set of position and velocity vectors represented in TEME reference
frame to a TLE. The arguments `args` and keywords `kwargs` are the same as those
described in the function `rv_to_mean_elements_sgp4`.

Additionally, the user can specify some parameters of the generated TLE.

This function returns the TLE and the covariance of the estimated elements
(state vector).

"""
function rv_to_tle(args...;
                   name::String = "UNDEFINED",
                   sat_num::Int = 9999,
                   classification::Char = 'U',
                   int_designator = "999999",
                   elem_set_number::Int = 0,
                   rev_num::Int = 0,
                   sgp4_gc::SGP4_GravCte = sgp4_gc_wgs84,
                   kwargs...)

    # Convert the position and velocity vectors to mean elements.
    JD, x, P = rv_to_mean_elements_sgp4(args...; kwargs...)

    # Compute the data as required by the TLE format.
    dt  = jd_to_date(DateTime, JD)
    dt??? = jd_to_date(DateTime, date_to_jd(year(dt), 1, 1, 0, 0, 0))

    dt_year    = year(dt)
    epoch_year = dt_year < 1980 ? dt_year - 1900 : dt_year - 2000
    epoch_day  = (dt - dt???).value/1000/86400 + 1

    # Obtain the Keplerian elements with the right units.
    a??? = x[1]/(1000sgp4_gc.R0)
    e??? = x[2]
    i??? = rad2deg(x[3])
    ????? = rad2deg(x[4])
    ????? = rad2deg(x[5])
    M??? = rad2deg(f_to_M(e???, x[6]))

    # Obtain the mean motion [rad/min].
    n??? = sgp4_gc.XKE / sqrt(a??? * a??? * a???)

    # Construct the TLE.
    tle = TLE(name,
              sat_num,
              classification,
              int_designator,
              epoch_year,
              epoch_day,
              JD,
              0.0,
              0.0,
              x[7],
              elem_set_number,
              0,
              i???,
              ?????,
              e???,
              ?????,
              M???,
              720n???/??,
              rev_num,
              0)

    # Return the TLE string.
    return tle, P
end


################################################################################
#                              Private functions
################################################################################

# Compute the SGP4 algorithm considering all variables in a state vector.
function _sgp4_sv(??t::Number,
                  sgp4_gc::SGP4_GravCte,
                  epoch::Number,
                  rx_TEME::Number,
                  ry_TEME::Number,
                  rz_TEME::Number,
                  vx_TEME::Number,
                  vy_TEME::Number,
                  vz_TEME::Number,
                  bstar::Number = 0)

    r_TEME = @SVector [rx_TEME, ry_TEME, rz_TEME]
    v_TEME = @SVector [vx_TEME, vy_TEME, vz_TEME]

    orb_TEME = rv_to_kepler(r_TEME, v_TEME, epoch)

    # Obtain the required mean elements to initialize the SGP4.
    a??? = orb_TEME.a/(1000sgp4_gc.R0) # .................... Semi-major axis [ER]
    e??? = orb_TEME.e                  # ........................ Eccentricity [ ]
    i??? = orb_TEME.i                  # ....................... Inclination [rad]
    ????? = orb_TEME.??                  # .............................. RAAN [rad]
    ????? = orb_TEME.??                  # ................... Arg. of perigee [rad]
    M??? = f_to_M(e???, orb_TEME.f)      # ...................... Mean anomaly [rad]

    # Obtain the mean motion [rad/min].
    n??? = sgp4_gc.XKE / sqrt(a??? * a??? * a???)

    r, v, sgp4d = sgp4(??t/60, sgp4_gc, epoch, n???, e???, i???, ?????, ?????, M???, bstar)

    # Return the elements using SI units.
    return 1000r, 1000v, sgp4d
end

function _sgp4_jacobian(??t::T,
                        epoch::T,
                        x???::SVector{NS, T},
                        y???::SVector{NO, T};
                        estimate_bstar::Bool = true,
                        pert::T = 1e-3,
                        perttol::T = 1e-5,
                        sgp4_gc::SGP4_GravCte = sgp4_gc_wgs84) where {NS, NO, T}

    num_states = NS
    dim_obs    = NO
    M          = zeros(T, dim_obs, num_states)

    # Auxiliary variables.
    x??? = copy(x???)

    # If B* does not need to be estimated, then does not compute the last
    # column.
    J = num_states - !estimate_bstar

    @inbounds for j = 1:J
        # State that will be perturbed.
        ?? = x???[j]

        # Obtain the perturbation, taking care to avoid small values.
        ?? = T(0)
        pert_i = pert

        for _ = 1:5
            ?? = ?? * pert_i
            abs(??) > perttol && break
            pert_i *= 1.4
        end

        ?? += ??

        # Avoid division by zero in cases that ?? is very small. In this
        # situation, we force |??| = perttol.
        abs(??) < perttol && (?? = sign(??) * perttol)

        # Modify the perturbed state.
        x??? = setindex(x???, ??, j)

        # Obtain the Jacobian by finite differentiation.
        r, v, ~ = _sgp4_sv(??t, sgp4_gc, epoch, x???...)
        y???      = vcat(r,v)

        M[:,j] .= (y??? .- y???)./??

        # Restore the value of the perturbed state for the next iteration.
        x??? = setindex(x???, x???[j], j)
    end

    return M
end

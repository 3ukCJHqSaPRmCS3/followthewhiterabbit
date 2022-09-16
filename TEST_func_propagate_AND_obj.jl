# Test:
#     Given Sat population [a,e,i,RAAN,w,M][6xN] @ t=t0
#     Propagator for t0->tn
#         output: ECEF, lat, lon
#
#     obj_fun1: Earth global coverage
#     obj_fun2: Earth local coverage = global*weight vector
#     obj_fun3: Quality of coverage
#                 a) Unobstructed time over point(s) of interest on Earth
#                 b)
#
#     Create test show cases
#         a.) city M, P, RM, NY, LA
#             import / create dictionary of city names with lat lon description
#         b.) Pick ship routes
#             Look for ship routes information [lat lon] e.g. what is the busiest global ship route? OR Shanghai <-> LA
##########################################################################################################################
dir_main = @__DIR__
dir_functions = "\\functions_astro"
push!(LOAD_PATH, string(dir_main, dir_functions))

using BenchmarkTools
using OrdinaryDiffEq    # used in: Propagation module
using Astrofunctions
############################################# NOTE TO MOVE TO main function #############################################################################
using NearestNeighbors
using NetCDF

Rearth_km 	 = 6378.137;    # [km]
Alt_kd_km 	 = 0;
foldername 	 = string(dir_main,"\\_spatial_grids\\");
filename_JPL = "JPL_MSCNv01_PLACEMENT.nc"; # source: https://podaac.jpl.nasa.gov/dataset/TELLUS_MASCON_GRID_POINT_PLACEMENT_V1
# ncinfo(string(foldername, filename_JPL))
JPL_mascon_radius_km      =    ncread(string(foldername, filename_JPL), "mascon_rad");          # Float64: radius  1xN [km]
JPL_mascon_lat_center_deg =    ncread(string(foldername, filename_JPL), "mascon_lat");          # Float64: center of latitude 1xN [deg]
JPL_mascon_lon_center_deg =    ncread(string(foldername, filename_JPL), "mascon_lon");          # Float64: center of Longitude 1xN [deg]
JPL_mascon_id             =    ncread(string(foldername, filename_JPL), "mascon_id");           # integer: id of spatial cell 1xN [-]
Nr_mascons                =    trunc(Int, JPL_mascon_id[end]);                                  # integer: number of spatial cell 1x1 [-]
############################ KD Tree ##############################
JPL_mascon_lon_center_deg_180 = JPL_mascon_lon_center_deg .- 180;
JPL_mascon_ECEF =[NaN*ones(3) for _ in 1:length(JPL_mascon_lat_center_deg)];

for i=1:length(JPL_mascon_lat_center_deg)
	JPL_mascon_ECEF[i] = func_latlon2ecef(JPL_mascon_lat_center_deg[i],JPL_mascon_lon_center_deg_180[i],Rearth_km,Alt_kd_km)
	# JPL_mascon_ECEF[i] = func_latlon2ecef(JPL_mascon_lat_center_deg[i],JPL_mascon_lon_center_deg[i],Rearth_km,Alt_kd_km)
end
kdtree = KDTree(hcat(JPL_mascon_ECEF...)); # build a KD tree


##########################################################################################################################
mu_km3s2  = 398600.4418;   # [km/s^3] mu of Earth
J2        = 0.0010826269;  # value of J2
Rearth_km = 6378.137;      # [km]
param 	  = [mu_km3s2, J2, Rearth_km];
w_earth_rad_s  = 7.2921158553E-5  # rad/s  Rotation rate of Earth
days2sec = 86400 #[s]

##########################################################################################################################


a_km = Rearth_km+500; ecc = 0; i_deg = 90; RAAN_deg = 45; w_deg = 5; M_deg=5; # [km] [deg]
const_TEST  = [a_km ecc i_deg    RAAN_deg w_deg M_deg]
const_TEST2 = [a_km ecc i_deg-30 RAAN_deg w_deg M_deg]
const_TEST3 = [a_km ecc i_deg-60 RAAN_deg w_deg M_deg]

npop = 5; nvar = 6;

population_kep = [repeat(const_TEST,outer=npop) repeat(const_TEST2,outer=npop) repeat(const_TEST3,outer=npop)]

nsat = Int(size(population_kep)[2]/6)

R_tot = [NaN*ones(nsat,3) for _ in 1:npop]
V_tot = [NaN*ones(nsat,3) for _ in 1:npop]
for ipop=1:npop
for isat = 1:nsat
	ind_i = nvar*(isat+1)-nvar-nvar+1
	a = population_kep[ipop,ind_i];      e = population_kep[ipop,ind_i+1];
	i = population_kep[ipop,ind_i+2];    W = population_kep[ipop,ind_i+3];
	w = population_kep[ipop,ind_i+4];    M = population_kep[ipop,ind_i+5];
	R, V = func_kep2eci(mu_km3s2,a,e,i,W,w,M)
	R_tot[ipop][isat,1:3] = R
	V_tot[ipop][isat,1:3] = V
end
end
# 0.000170 seconds (2.05 k allocations: 108.156 KiB)

t_start = 5; dt_sec = 5; tfinal_day = 1;
PropagationTimeSec = tfinal_day*days2sec
time_vec = t_start:dt_sec:PropagationTimeSec
time_steps = length(time_vec)


ECEF_tot = [NaN*ones(nsat,time_steps,6) for _ in 1:npop]
latlon_tot  = [NaN*ones(nsat,time_steps,2) for _ in 1:npop]

for ipop = 1:npop
for isat = 1:nsat
    Rtmp = R_tot[ipop][isat,:] #R[isat]
    Vtmp = V_tot[ipop][isat,:] #V[isat]
    u0 = [Rtmp; Vtmp]
    t = (Float64(t_start), PropagationTimeSec)
    prob = ODEProblem(Orbit_J2prop!, u0, t, param) # Setup Propagator
    sol = solve(prob,Vern8(),saveat = dt_sec) #Propagate
    ECI = sol.u
    ECEF = func_eci2ecef(ECI, PropagationTimeSec, t_start, dt_sec, w_earth_rad_s)
	lat_geod_deg, lon_geod_deg, lat_geoc_deg, lon_geoc_deg = func_ecef2latlon(ECEF, Rearth_km)
	ECEF_tot[ipop][isat,:,:] = hcat(ECEF...)'
	latlon_tot[ipop][isat,:,:] = hcat([lat_geod_deg, lon_geod_deg]...)
end
end


flag_plot_globe = false
flag_plot_Robinson = false
	if flag_plot_globe
		using GMT
		GMT.coast(region=[0 360 -90 90], proj=(name=:laea, center=(300,30)), frame=:g, res=:crude, land=:navy, figsize=6)
		GMT.plot!(lon_geod_deg, lat_geod_deg, lw=0.1, lc=:red, fmt=:png, marker=:circle, markeredgecolor=0, size=0.0005, markerfacecolor=:red, show=true)
	end
	if flag_plot_Robinson
		using GMT
		GMT.coast(region=:d, proj=:Robinson, frame=:g, res=:crude, area=10000, land="#CEB9A0", water=:lightblue, figsize=6)
		GMT.plot!(lon_geod_deg, lat_geod_deg, lw=0.1, lc=:red, fmt=:png, marker=:circle, markeredgecolor=0, size=0.0005, markerfacecolor=:red, show=true)
	end

# Get a list of cities with [Latitude, Longitude] information && get a ship route (easy [PortA<->B], medium [portA<->B<->C], hard [A..G])
# convert cities[Latitude,Longitude] -> ECEF
# convert ship routes[Lat,Lon] -> ECEF
cities_TEST_latlon_deg = [55.7558 37.6173;
						39.9042 116.4074;
						40.3399 127.5101;
						35.7219 51.3347] 
						#### Moskow,Beijing,Pyongyang,Teheran ################################################################################# NOTE to be changed #########################################################################
ntarget = size(cities_TEST_latlon_deg)[1]

flag_plot_Robinson_p_cities = false
if flag_plot_Robinson_p_cities # plot cities on MAP
	using GMT
	GMT.coast(region=:d, proj=:Robinson, frame=:g, res=:crude, area=10000, land="#CEB9A0", water=:lightblue, figsize=6)
	GMT.plot!(lon_geod_deg, lat_geod_deg, lw=0.1, lc=:red, fmt=:png, marker=:circle, markeredgecolor=0, size=0.0005, markerfacecolor=:red, show=true)
	GMT.plot!(cities_TEST_lonlat_deg[:,1], cities_TEST_lonlat_deg[:,2], lw=0.1, lc=:red, fmt=:png, marker=:circle, markeredgecolor=0, size=2, markerfacecolor=:red, show=true)
end

# end ##collect the ECEF + lat lon of all satellites per consetellation
##collect the ECEF of each sat of the constellation


## which location on Earth is the satellite(t=ti) the closest? && which ones are closest to it?
grid_idx_visited 	= zeros(npop,nsat,time_steps); grid_dis_visited = zeros(npop,nsat,time_steps);
grid_idx_visited10 	= zeros(npop,nsat,time_steps,10); grid_dis_visited10 = zeros(npop,nsat,time_steps,10);
grid_idx_target  	= zeros(ntarget);
grid_dis_target  	= zeros(ntarget);

# @sync @distributed for i=1:time_steps_Nsat_tot

co_obs = 1 # closest_outputs @ observation
for iobs=1:time_steps
for isat=1:nsat
for ipop=1:npop
	grid_idx_visited[ipop,isat,iobs] = NearestNeighbors.knn(kdtree, ECEF_tot[ipop][isat,iobs,1:3],co_obs)[1][1]
	grid_dis_visited[ipop,isat,iobs] = NearestNeighbors.knn(kdtree, ECEF_tot[ipop][isat,iobs,1:3],co_obs)[2][1]
end
end
end

co_obs = 1 # closest_outputs @ observation
for iobs=1:time_steps
for isat=1:nsat
for ipop=1:npop
	grid_idx_visited10[ipop,isat,iobs,:] = NearestNeighbors.knn(kdtree, ECEF_tot[ipop][isat,iobs,1:3],co_obs)[1]
	grid_dis_visited10[ipop,isat,iobs,:] = NearestNeighbors.knn(kdtree, ECEF_tot[ipop][isat,iobs,1:3],co_obs)[2]
end
end
end

co_tar = 1 # closest_outputs @ target
for itarget=1:ntarget
	lat, lon = cities_TEST_latlon_deg[itarget,:]
	ECEF_target = func_latlon2ecef(lat,lon,Rearth_km,Alt_kd_km)
	grid_idx_target[itarget] = NearestNeighbors.knn(kdtree, ECEF_target,co_tar)[1][1]
	grid_dis_target[itarget] = NearestNeighbors.knn(kdtree, ECEF_target,co_tar)[2][1]
end

grid_idx_visited_target = NaN*zeros(npop,nsat,time_steps)

for iobs=1:time_steps
for isat=1:nsat
for ipop=1:npop
	if sum(grid_idx_visited[ipop,isat,iobs] .== grid_idx_target)>=1 #observation has YES been performed at target
		grid_idx_visited_target[ipop,isat,iobs] = grid_idx_visited[ipop,isat,iobs]
	else #observation has NOT been performed at target -> do NOT apply changes
	end
end
end
end


target_cell_lat_lon = zeros(ntarget,2)
for itarget=1:ntarget
	target_cell_lat_lon[itarget,:] = [JPL_mascon_lat_center_deg[Int(grid_idx_target[itarget])] , JPL_mascon_lon_center_deg_180[Int(grid_idx_target[itarget])]]
end

grid_idx_visited_lat_lon = NaN*zeros(npop,nsat,time_steps,2)
grid_idx_visited_target_lat_lon = NaN*zeros(npop,nsat,time_steps,2)
for iobs=1:time_steps
for isat=1:nsat
for ipop=1:npop
	grid_idx_visited_lat_lon[ipop,isat,iobs,:] = [JPL_mascon_lat_center_deg[Int(grid_idx_visited[ipop,isat,iobs])] , JPL_mascon_lon_center_deg_180[Int(grid_idx_visited[ipop,isat,iobs])]]
	if isnan(grid_idx_visited_target[ipop,isat,iobs]) == false
		grid_idx_visited_target_lat_lon[ipop,isat,iobs,:] = [JPL_mascon_lat_center_deg[Int(grid_idx_visited_target[ipop,isat,iobs])] , JPL_mascon_lon_center_deg_180[Int(grid_idx_visited_target[ipop,isat,iobs])]]
	end
end
end
end


# plot the observed targets
using GMT
GMT.coast(region=:d, proj=:Robinson, frame=:g, res=:crude, area=10000, land="#CEB9A0", water=:lightblue, figsize=6)
GMT.plot!(target_cell_lat_lon[:,2], target_cell_lat_lon[:,1], fmt=:png, marker=:diamod, markeredgecolor=0, size=0.05, markerfacecolor=:red, show=true)
# plot all ground track observations : blue
GMT.coast(region=:d, proj=:Robinson, frame=:g, res=:crude, area=10000, land="#CEB9A0", water=:lightblue, figsize=6)

ipop=1
lat_tmp = reshape(grid_idx_visited_lat_lon[ipop,:,:,1],(nsat*time_steps,1))
lon_tmp = reshape(grid_idx_visited_lat_lon[ipop,:,:,2],(nsat*time_steps,1))
GMT.coast(region=:d, proj=:Robinson, frame=:g, res=:crude, area=10000, land="#CEB9A0", water=:lightblue, figsize=6)
GMT.plot!(lon_tmp, lat_tmp, fmt=:png, marker=:circle, markeredgecolor=0, size=0.05, markerfacecolor=:red, show=true)


ipop=1
lat_tmp = reshape(grid_idx_visited_target_lat_lon[ipop,:,:,1],(nsat*time_steps,1))
lon_tmp = reshape(grid_idx_visited_target_lat_lon[ipop,:,:,2],(nsat*time_steps,1))
GMT.coast(region=:d, proj=:Robinson, frame=:g, res=:crude, area=10000, land="#CEB9A0", water=:lightblue, figsize=6)
GMT.plot!(lon_tmp, lat_tmp, lw=0.1, lc=:red, fmt=:png, marker=:circle, markeredgecolor=0, size=0.05, markerfacecolor=:red, show=true)


# GMT.plot!(cities_TEST_lonlat_deg[:,1], cities_TEST_lonlat_deg[:,2], lw=0.1, lc=:red, fmt=:png, marker=:circle, markeredgecolor=0, size=2, markerfacecolor=:red, show=true)

# MAP (lat lon):

# plot the ground track (immediately before over + after the target) : red


ipop=1
lat_tmp = reshape(grid_idx_visited_target_lat_lon[ipop,:,:,1],(nsat*time_steps,1))
lon_tmp = reshape(grid_idx_visited_target_lat_lon[ipop,:,:,2],(nsat*time_steps,1))
GMT.coast(region=:d, proj=:Robinson, frame=:g, res=:crude, area=10000, land="#CEB9A0", water=:lightblue, figsize=6)
GMT.plot!(lon_tmp, lat_tmp, lw=0.1, lc=:red, fmt=:png, marker=:circle, markeredgecolor=0, size=0.05, markerfacecolor=:red, show=true)


# ONLY save index interested instead of creating a new variable!
# create animation
# ground track trail with dissapearing line
# highlight observations over target

#
# Given target location on Earth & observations
# determine quantity and "quality" of observation








# findall(isequal(maximum(grid_dis_visited)),grid_dis_visited)
# a = findall(x -> x >=600 , grid_dis_visited)
# grid_dis_visited[a]
# grid_idx_visited[1,2,12]

# test = reshape(grid_dis_visited,size(grid_dis_visited)[1]*size(grid_dis_visited)[2]*size(grid_dis_visited)[3],1)
# 17280*nsat*npop
# unique(sort(test,dims=1))

# using Plots
# plot(1:length(unique(sort(test,dims=1))),unique(sort(test,dims=1)),label=false)

## check for observation values for that location + surrounding locations
# ipop=1
# isat=2
# # lon_geod_deg = latlon_tot[ipop][isat,:,2]
# # lat_geod_deg = latlon_tot[ipop][isat,:,1]
# lon_geod_deg_tmp = zeros(length(a))
# lat_geod_deg_tmp = zeros(length(a))
# for i=1:length(a)
# 	lon_geod_deg_tmp[i] = latlon_tot[a[i][1]][a[i][2],a[i][3],2]
# 	lat_geod_deg_tmp[i] = latlon_tot[a[i][1]][a[i][2],a[i][3],1]
# end
#
# using GMT
# GMT.coast(region=:d, proj=:Robinson, frame=:g, res=:crude, area=10000, land="#CEB9A0", water=:lightblue, figsize=6)
# GMT.plot!(lon_geod_deg_tmp, lat_geod_deg_tmp, lw=0.1, lc=:red, fmt=:png, marker=:circle, markeredgecolor=0, size=0.0005, markerfacecolor=:red, show=true)



## OBJ function
#Spatial quality of observation
# LOS_km = sqrt( (ECEF[1]-POI[1])^2 + (ECEF[2]-POI[2])^2 + (ECEF[3]-POI[3])^2 )
# beta_deg =

# define range(LOS or max distance , angle from the horizon due to obstructions)
	# in Line of sight (LOS) + 5 degrees above horizon to account for obstructions
	# what is the observational range of current on-board instrumentation? Search for sources. Otherwise assume LOS (aka as long as it is above the horizon + 5deg)
	# to compute the local angle transform ECEF to local position (see HMW 1-2 6080 StatOD)
# if satellite is in range
	# calculate the beta angle between horizon and satellite
	# if above horizon
		# observation is accepted


#Temporal quality of observation
	# Duration: time length of observation -> length of observation vector
	# Uniformity: of seperate observation instances -> gini
	#


# create animation about how the ground track changes and adapts to the cities/ship routes chosen
# gen 1..5 etc

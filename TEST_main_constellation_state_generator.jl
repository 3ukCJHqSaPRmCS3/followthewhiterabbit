# load JPL mascon files
# load NGSA-II
# initialize constellations

# dir_save_data = "_GA_data\\"
dir_main = @__DIR__
dir_functions = "\\functions_astro"
push!(LOAD_PATH, string(dir_main, dir_functions))

using NetCDF            # Needed for JPL mascons
using NearestNeighbors
using Astrofunctions



############################ CONSTANTS ########################################
mu_km3s2  = 398600.4418; # [km/s^3] mu of the Earth
we = 7.2921158553E-5;       # Rotation rate of the Earth
mu = 398600.4418;           # mu of the Earth
J2 = 0.0010826269;          # value of J2
Rearth_km = 6378.137;    # [km]
daysinsec = 86400;       # [s]
params 	  = [mu_km3s2, J2, Rearth_km];
Alt_kd_km = 0;
###############################################################################

############################ Inputs ########################################
a_0_km = 500 #km
# a_vec_km = [ 400 300 ]
# a_elevator_km = a_vec_km[end]:25:250
Targets_lat_lon = [55 37] # 55°45′21″N 37°37′2″E Moskow
weights = [100]
###############################################################################

######################### READ IN MASCONS of choice ###########################
# Call JPL Mascons
foldername = string(dir_main,"\\_spatial_grids\\");
filename_JPL = "JPL_MSCNv01_PLACEMENT.nc"; # source: https://podaac.jpl.nasa.gov/dataset/TELLUS_MASCON_GRID_POINT_PLACEMENT_V1
# ncinfo(string(foldername, filename_JPL))
JPL_mascon_radius_km      =    ncread(string(foldername, filename_JPL), "mascon_rad");          # Float64: radius  1xN [km]
JPL_mascon_lat_center_deg =    ncread(string(foldername, filename_JPL), "mascon_lat");          # Float64: center of latitude 1xN [deg]
JPL_mascon_lon_center_deg =    ncread(string(foldername, filename_JPL), "mascon_lon");          # Float64: center of Longitude 1xN [deg]
JPL_mascon_id             =    ncread(string(foldername, filename_JPL), "mascon_id");           # integer: id of spatial cell 1xN [-]
Nr_mascons                =    trunc(Int, JPL_mascon_id[end]);                                  # integer: number of spatial cell 1x1 [-]
############################ KD Tree ##############################
# JPL_mascon_ECEF = [NaN * ones(3) for _ = 1:length(JPL_mascon_lat_center_deg)];
JPL_mascon_lon_center_deg_180 = JPL_mascon_lon_center_deg .- 180;



######################### READ IN MASCONS of choice ###########################
###############################################################################
##################################NOTE TODO ###################################
###############################################################################
######################### READ IN MASCONS of choice ###########################

# create KD tree
JPL_mascon_ECEF =[NaN*ones(3) for _ in 1:length(JPL_mascon_lat_center_deg)];
JPL_mascon_lon_center_deg_180 = JPL_mascon_lon_center_deg.-180;

for i=1:length(JPL_mascon_lat_center_deg)
	JPL_mascon_ECEF[i] = func_latlon2ecef(JPL_mascon_lat_center_deg[i],JPL_mascon_lon_center_deg_180[i],Rearth_km,Alt_kd_km)
end
kdtree = KDTree(hcat(JPL_mascon_ECEF...)); # build a KD tree

# data = [JPL_mascon_lat_center_deg, JPL_mascon_lon_center_deg.-180];







##################### USER INPUTS #####################
################# IDENTIFY constellation constants #################
PropagationTime_days = 1
halt_km = 500
a = Rearth_km + halt_km    # semimajor axis (km)
e = 0           # eccentricity
w = 90          # argument of perigee (deg)
nsat = 3
TC = 16
RP = 29 # NOTE or free?###########################################################################################################################################################################################################################
t_start = 5
dt_sec  = 5                                                     # Sampling time [s]

####################################################################
# explanation: [inc, RAAN, M, RD]
nbits_vec    = [5, 5, 5, 0]

hi_inc   = 179; lo_inc   = 1;
hi_RAAN  = 360; lo_RAAN  = trunc(Int,hi_RAAN/2^nbits_vec[2])
hi_M     = 360; lo_M  	 = trunc(Int,hi_M/2^nbits_vec[3])

lim_vector  = [lo_inc, hi_inc, lo_RAAN, hi_RAAN, lo_M, hi_M];

PropagationTimeSec 	   = round(PropagationTime_days * daysinsec);
Propagation_time_vec_sec   = collect(t_start:dt_sec:PropagationTimeSec);
Propagation_steps	   = length(Propagation_time_vec_sec);
Propagation_steps_vec_Nsat = collect(0:Propagation_steps:nsat*Propagation_steps);


############################ GEN alg SETUP ############################

popsize 			= 100 											  # Size of the population per generation NOTE needs to be an EVEN number
run_max				= 1												  # Number of iterations
gen_max				= 40 											  # Number of generations per iteration
mutation_rate 			= 0.5
tol_SC  			= 1E-3
tol_TC 				= 1E-3
crossover_method 		= "n_point_uniform_parents"
npar    = 4 #size(nbits_vec)[1]
obj_var = 2
ind_tmp01 = (nsat*npar);                # V
ind_tmp02 = (nsat*npar)+obj_var;        # V+M
#########################################################################

# Staffetta_constants = [we, mu_km3s2, J2, Rearth_km, disp_info_RVgenerator];
# Staffetta_constants = [we, mu_km3s2, J2, Rearth_km];
# Staffetta_sizes     = [Nsat, npar, obj_var, popsize];
# Steffetta_sat_const = [a, e, w, sep_req, Nsat_chain_tot];
# Staffetta_grid      = [kdtree, data, Nr_mascons];
# Staffetta_time      = [t_start, dt_sec, PropagationTimeSec, Propagation_steps, Propagation_steps_vec_Nsat, N_time_mascon, TM_N];
# Staffetta_ga	    = [crossover_method, mutation_rate];
# Staffetta_lim	    = [lim_vector_i; lim_vector_R; lim_vector_M; lim_vector_RD; nbits_vec];
# Staffetta_setup	    = [PropagationTimeDays,popsize,run_max,gen_max,Nsat,Staffetta_lim,Staffetta_ga];
# Staffetta_setup_ALL = [Staffetta_constants, Staffetta_sizes, Steffetta_sat_const, Staffetta_grid, Staffetta_time, Staffetta_ga];


# setup initial values and data_save matrix
# population_SAVE = NaN*ones(popsize,ind_tmp02+3,gen_max+1, run_max);
# values_init 	= NaN*ones(popsize,nsat*npar);


###########################################################################################################################################################################################################################

# create random population
# for isat = 1:Nsat
# 	if isat == 1
# 		pol_inc_val, pol_RAAN_val, pol_M_val, pol_RD_val = func_value_generation_v3(isat, Staffetta_lim, popsize)
# 		values_init[1:popsize,1:npar] = [pol_inc_val pol_RAAN_val pol_M_val pol_RD_val]
#
# 	else
# 		lim1_tmp = npar*(isat+1)-npar-npar+1
# 		lim2_tmp = npar*(isat+1)-npar
#
# 		inc_val, RAAN_val, M_val, RD_val = func_value_generation_v3(isat, Staffetta_lim, popsize)
# 		values_init[1:popsize,lim1_tmp:lim2_tmp] = [ inc_val RAAN_val M_val RD_val ]
# 	end
# end

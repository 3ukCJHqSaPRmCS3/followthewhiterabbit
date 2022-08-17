# develop equations for spherical harmonics time-varying gravity coefficients
using SatelliteToolbox

# change this to create coefs(t)
# coefs,t_vec = load_gravity_model(ESA(t_start,t_end))
# coefs = CS[degree x ord x time]
# t_vec = t_start:dt:t_end

#################################### START: INPUT #####################################
lat_deg = -22 # [deg]
lon_deg = -45 # [deg]
h_km    = 700 # [km]
deg_max = 20 # [-]
ord_max = deg_max # usually = deg_max
time_vector = YYYYMMDDHH
#################################### END: INPUT #####################################
coefs = load_gravity_model(EGM96())
g1 = compute_g(coefs, geodetic_to_ecef(lat_deg*pi/180, lon_deg*pi/180, h_km), deg_max, ord_max)  # [m/s²]

# # LOAD the saved data
using JLD
dir_main = @__DIR__
dir_folder = "ESM_models"
dir_models = string(dir_main,"\\",dir_folder)
setup_name = "mtm3h_2006_AOHIS"
d = load(string(dir_models,"\\",setup_name,".jld"))
# C_mat = d["C_mat"];
# S_mat = d["S_mat"];
# time_vec = d["time_vec"];
# time_step = 1
# itime     = time_vec[time_step]

#######################


fileName = "ESM_models\\mtm3h_2006_A\\01\\test.360"
####################### read in time vector from file name #######################

new_model_tmp = create_gravity_model_coefs(parse_icgem(string(dir_main,"\\",fileName)))
g2 = compute_g(new_model_tmp, geodetic_to_ecef(lat_deg*pi/180, lon_deg*pi/180, h_km), deg_max, ord_max)  # [m/s²]


################################ TODO ################################
## change header from all AOHIS files
# # from this format
# product type            anomalous gravity potential
# modelname               Improved Mass Transport Model
# model content           AOHIS
# version                 1.0
# earth_gravity_constant  0.39860050000000D+15
# radius                  6378137.0000
# max_degree              360
# error                   no
# norm                    fully_normalized
# tide_system             does-not-apply
# end_of_head
# #to this format
# product_type                gravity_field
# modelname                   EGM96
# earth_gravity_constant      0.3986004415E+15
# radius                      0.6378136300E+07
# max_degree                     360
# errors                      no
# tide_system                 tide_free
# end_of_head

# Use sed D/E bash script for quick changes
# tide_system             does-not-apply
# tide_system                 tide_free

# earth_gravity_constant  0.39860050000000D+15
# earth_gravity_constant      0.3986004415E+15

# modelname               Improved Mass Transport Model
# TO modelname            ESM_AOHIS
# modelname                   EGM96

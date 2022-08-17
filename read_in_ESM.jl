# Read in ESM CS coefficients.
# Name: module_codeDeposit_tools
# Author: CMA Deccia
# Creation date: June 15 2022
# Modifications: ...

dir_main = @__DIR__
dir_function = joinpath(@__DIR__,"function_general\\")
push!(LOAD_PATH, dir_main)
push!(LOAD_PATH, dir_function)

using module_codeDeposit_tools       # LOAD general functions
using JLD
using ProgressMeter

flag_load_models = false
if flag_load_models
dir_folder = "ESM_models"
for imodel = 1:6
    if imodel==1
        dir_model = "mtm3h_2006_A"
    elseif imodel==2
        dir_model = "mtm3h_2006_O"
    elseif imodel==3
        dir_model = "mtm3h_2006_H"
    elseif imodel==4
        dir_model = "mtm3h_2006_I"
    elseif imodel==5
        dir_model = "mtm3h_2006_S"
    elseif imodel==6
        dir_model = "mtm3h_2006_AOHIS"
    end
    dir_name_month = readdir(string(dir_main,"\\",dir_folder,"\\",dir_model))
    # initialize N steps
    dt = 3
    files_p_day = Int(24/dt)
    tot_time_steps = 365*files_p_day-1 # sum([248 224 248 240 248 240 248 248 240 248 240 247]) # 2919 time steps
    itime = 1
    # setup matrix size
    nmax = 360
    C_mat = zeros(nmax+1,nmax+1,tot_time_steps);
    S_mat = zeros(nmax+1,nmax+1,tot_time_steps);
    time_vec = zeros(tot_time_steps);
    @showprogress 1 "Computing..." for ifolder in 1:length(dir_name_month)
        file_name_vector = readdir(string(dir_main,"\\",dir_folder,"\\",dir_model,"\\",dir_name_month[ifolder]))
        for ifile = 1:length(file_name_vector)
            data_tmp = readlines(string(string(dir_main,"\\",dir_folder,"\\",dir_model,"\\",dir_name_month[ifolder]),"\\",file_name_vector[ifile]))
            indx = findall( x -> occursin("gfc", x), data_tmp)
            x_tmp = permutedims(hcat(split.(data_tmp[indx])...))
            if length(dir_model) == 12
                time_vec[itime] = parse(Int64, join([file_name_vector[ifile][10:10+7], file_name_vector[ifile][10+9:10+10]]) ) # date format 2006010100 => 2006 01 01 00 YYYY MM DD HH NOTE: HH = [03,06,09,12,15,18,21]
            elseif length(dir_model) == 16
                time_vec[itime] = parse(Int64, join([file_name_vector[ifile][14:14+7], file_name_vector[ifile][14+9:14+10]]) ) # date format 2006010100 => 2006 01 01 00 YYYY MM DD HH NOTE: HH = [03,06,09,12,15,18,21]
            end
            for iline = 1:size(x_tmp)[1]
                deg_tmp = parse(Int,x_tmp[iline,2])  #deg
                ord_tmp = parse(Int,x_tmp[iline,3]) #ord
                C_tmp   = parse(Float64,func_convert_D2E(x_tmp[iline,4])) #C_val
                S_tmp   = parse(Float64,func_convert_D2E(x_tmp[iline,5])) #S_val
                C_mat[deg_tmp+1,ord_tmp+1,itime] = C_tmp
                S_mat[deg_tmp+1,ord_tmp+1,itime] = S_tmp
            end
            itime+=1
        end
    end
    # save variables in jld file
    dir_models = string(dir_main,"\\",dir_folder)
    setup_name = dir_model
    save(string(dir_models,"\\",setup_name,".jld"),"C_mat",C_mat,"S_mat",S_mat,"time_vec",Int.(time_vec),"model_name",dir_model)
end

# output:
# Int.(time_vec), C_mat, S_mat, nmax, model_name = "mtm3h_2006_AOHIS"

end

##
# # LOAD the saved data
using JLD
dir_main = @__DIR__
dir_folder = "ESM_models"
dir_models = string(dir_main,"\\",dir_folder)

setup_name = "mtm3h_2006_A"
d = load(string(dir_models,"\\",setup_name,".jld"))
C_A_mat = d["C_mat"];
S_A_mat = d["S_mat"];
time_A_vec = d["time_vec"];

setup_name = "mtm3h_2006_O"
d = load(string(dir_models,"\\",setup_name,".jld"))
C_O_mat = d["C_mat"];
S_O_mat = d["S_mat"];
time_O_vec = d["time_vec"];

setup_name = "mtm3h_2006_H"
d = load(string(dir_models,"\\",setup_name,".jld"))
C_H_mat = d["C_mat"];
S_H_mat = d["S_mat"];
time_H_vec = d["time_vec"];

setup_name = "mtm3h_2006_I"
d = load(string(dir_models,"\\",setup_name,".jld"))
C_I_mat = d["C_mat"];
S_I_mat = d["S_mat"];
time_I_vec = d["time_vec"];

setup_name = "mtm3h_2006_S"
d = load(string(dir_models,"\\",setup_name,".jld"))
C_S_mat = d["C_mat"];
S_S_mat = d["S_mat"];
time_S_vec = d["time_vec"];

setup_name = "mtm3h_2006_AOHIS"
d = load(string(dir_models,"\\",setup_name,".jld"))
C_AOHIS_mat = d["C_mat"];
S_AOHIS_mat = d["S_mat"];
time_AOHIS_vec = d["time_vec"];


C_AOHIS2_mat = C_A_mat + C_O_mat + C_H_mat + C_I_mat + C_S_mat
C_AOHI_mat = C_A_mat + C_O_mat + C_H_mat + C_I_mat

test = C_AOHIS2_mat - C_AOHIS_mat
sum(C_AOHI_mat[:,:,1:240],dims=3)

# create a way for SatelliteToolbox to read in AOHIS files instead of the given ones
# can you modify the file that reads in the text files??? if you take it out of the module?
# call the jl files directly OR the text files
# call it outside the toolbox
# if yes
# copy all jl files required: outside the toolbox
# create new module to load all files
# with the modified jl readin file + text files for the NEW CS ESM models

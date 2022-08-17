###################################################################
# Written by:
# Carlos Deccia
# July 12th, 2022
###################################################################

##Declare module name:
module Astrofunctions

	## Modules used by func_XXX
	# e.g. using ProgressMeter

	## List of functions exported by this module:
	export Orbit_J2prop!
	export func_kep2eci
	export func_eci2ecef
	export func_zeroTO360
	export func_R3
	export func_ecef2latlon
	export func_latlon2ecef

	## List of Julia scripts included in this module:
	include("Orbit_J2prop!.jl")
	include("func_kep2eci.jl")
	include("func_eci2ecef.jl")
	include("func_zeroTO360.jl")
	include("func_R3.jl")
	include("func_ecef2latlon.jl")
	include("func_latlon2ecef.jl")



end

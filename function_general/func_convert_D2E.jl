# Name: module_codeDeposit_tools
# Written by:
# Carlos Deccia
# Morpheus Space
# Last Modified: June 15 2022

"""
    func_convert_D2E(str)

# input: string
#
# function: Converts "D" to "E" in a given string 'str'.
#
# output: string

# Example:
```julia-repl
julia> func_convert_D2E("-0.87763439035909D-10")
-0.87763439035909E-10
```
"""


function func_convert_D2E(str)
    subs = Dict("D" => "E")
    subbed = replace(str, r"D" => s -> subs[s])
end

using KitAMR
include("airfoil/airfoil_udf.jl")
include("convergence_ps/convergence_ps_udf.jl")
include("cylinder/cylinder_udf.jl")
include("sphere/sphere_udf.jl")

KitAMR.boundary_result2csv("RESULT_FOLDER","NAME")   
## Run all the greenfields for fusion and no_fusion

# Fusion cases

for i in 2001:2021
    
    include("fusion\\$i\\Run.jl") 

       
    include("no_fusion\\$i\\Run.jl") 


end

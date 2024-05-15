using LinearAlgebra

G_values = [-5.92, -5.942, -5.97, -5.948, -5.937, -5.94] #Amorphous unsure
deltaG = [G_values[j] - G_values[i] for j in 1:length(G_values), i in 1:length(G_values)]
Ea_constants = [0.00 1.00 0.01 0.01 1.00 1.00;  #set the activation energy from amorphous to kappa and beta to be small
                1.00 0.00 1.00 1.00 1.00 1.00;
                0.01 1.00 0.00 0.01 1.00 1.00;
                0.01 1.00 0.01 0.00 1.00 1.00;
                1.00 1.00 1.00 1.00 0.00 1.00;
                1.00 1.00 1.00 1.00 1.00 0.00]
Ea_constants = Symmetric(Ea_constants)
Ea_values = [
    if i == j # Diagonal elements are 0
        0
    else
        if deltaG[i, j] > 0 # deltaG>0 means Ea = deltaG + Ea
            deltaG[i, j] + Ea_constants[i, j]
        else # deltaG<0 means Ea = Ea
            Ea_constants[i, j]
        end
    end
    for j in 1:length(G_values), i in 1:length(G_values)
]
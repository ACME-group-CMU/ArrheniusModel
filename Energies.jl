G_values = [-5.92, -5.942, -5.97, -5.948, -5.937, -5.94] #Amorphous unsure
deltaG = [G_values[j] - G_values[i] for j in 1:length(G_values), i in 1:length(G_values)]
Ea_constants = [0.00 1.00 0.01 0.01 1.00 1.00;  #set the activation energy from amorphous to kappa and beta to be small
                1.00 0.00 1.00 1.00 1.00 1.00;
                0.01 1.00 0.00 0.01 1.00 1.00;
                0.01 1.00 0.01 0.00 1.00 1.00;
                1.00 1.00 1.00 1.00 0.00 1.00;
                1.00 1.00 1.00 1.00 1.00 0.00]
Ea_values = [i == j ? 0 : (deltaG[i, j] > 0 ? deltaG[i, j] + Ea_constants[i, j] : Ea_constants[i, j])
                for j in 1:length(G_values), i in 1:length(G_values)]
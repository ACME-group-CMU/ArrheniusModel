G_values = [-5.93, -5.942, -5.97-0.3, -5.948-0.2, -5.937, -5.94] #Amorphous unsure
deltaG = [G_values[j] - G_values[i] for j in 1:length(G_values), i in 1:length(G_values)]
Ea_constants = [0.000 0.001 0.001 0.001 0.001 0.001;  #set the activation energy from amorphous to small
                0.001 0.000 0.002 0.004 0.005 0.006;
                0.001 0.002 0.000 0.002 0.002 0.002;
                0.001 0.004 0.002 0.000 0.004 0.004;
                0.001 0.005 0.002 0.004 0.000 0.006;
                0.001 0.006 0.002 0.004 0.006 0.000]

Ea_values = [i == j ? 0 : (deltaG[i, j] > 0 ? deltaG[i, j] + Ea_constants[i, j] : Ea_constants[i, j])
                for j in 1:length(G_values), i in 1:length(G_values)]
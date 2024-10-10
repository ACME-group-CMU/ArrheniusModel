using Enzyme

function sens_arrhenius_rate(b, T)
    T_ = [T]
    dT = [0.0]
    db = zero(b)
    K = zero(b)
    dK = zero(b)
    sens_b = [zero(b) for _ in 1:size(dK, 1), _ in 1:size(dK, 2)]
    sens_T = [zero(dT) for _ in 1:size(dK, 1), _ in 1:size(dK, 2)]
    for i in 1:size(dK, 1)
        for j in 1:size(dK, 2)
            dK .= 0.0  # Reset dK to zero
            dK[i, j] = 1.0
            db = zero(b)  # Reinitialize db
            dT = [0.0]    # Reinitialize dT
            Enzyme.autodiff(ReverseWithPrimal, ar_matrixT!, Duplicated(b, db), Duplicated(K, dK), Duplicated(T_, dT))
            sens_b[i, j] = copy(db)  # Store the entire db matrix
            sens_T[i, j] = copy(dT)  # Store the entire dT vector
        end
    end
    return sens_b, sens_T
end

function sens_simulation(flow_rate, T, pe, para_sim)
    num_steps, num_layers, dt = para_sim
    f_ = [flow_rate]
    T_ = [T]
    df = [0.0]
    dT = [0.0]
    b = copy(pe.barriers)
    db = zero(b)
    s = zeros(num_layers, n_phases(pe))
    ds = zero(s)
    sens_f = [zero(f_) for _ in 1:size(ds, 1), _ in 1:size(ds, 2)]
    sens_T = [zero(T_) for _ in 1:size(ds, 1), _ in 1:size(ds, 2)]
    sens_b = [zero(b) for _ in 1:size(ds, 1), _ in 1:size(ds, 2)]
    for i in 1:size(ds, 1)
        for j in 1:size(ds, 2)
            df = [0.0]  # Reinitialize df
            dT = [0.0]  # Reinitialize dT
            db = zero(b)  # Reinitialize db
            ds[i,j] = 1.0
            Enzyme.autodiff(ReverseWithPrimal, simulate_deposition!, Duplicated(s, ds), Duplicated(f_, df), Duplicated(T_, dT), Duplicated(b, db), Const(para_sim))
            sens_f[i, j] = copy(df)  # Store the entire df vector
            sens_T[i, j] = copy(dT)  # Store the entire dT vector
            sens_b[i, j] = copy(db)  # Store the entire db matrix
        end
    end
    return sens_f, sens_T, sens_b
end
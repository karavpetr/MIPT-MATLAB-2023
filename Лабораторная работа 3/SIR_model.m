function dXdt = SIR_model(t, X, params)
    alpha = params(1);
    beta_1 = params(2);
    beta_2 = params(3);
    k = params(4);
    gamma = params(5);
    delta = params(6);
    omega = params(7);
    c = params(8);
    
    S_1 = X(1);
    S_2 = X(2);
    I = X(3);
    R = X(4);
    u = X(5);
    
    V_1 = c / (delta + beta_1 * I);
    V_2 = c / (delta + beta_2 * I);
    
    dS_1dt = - beta_1 * I * S_1 + k * S_2;
    dS_2dt = - beta_2 * I * S_2 + alpha * R - k * S_2;
    dIdt = beta_1 * I * S_1 + beta_2 * I * S_2 - gamma * I;
    dRdt = gamma * I - delta * R - alpha * R;
    dudt = omega * (V_2 - V_1) * u;
    
    dXdt = [dS_1dt; dS_2dt; dIdt; dRdt; dudt]; 
end

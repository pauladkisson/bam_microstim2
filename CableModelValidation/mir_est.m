function Vmir = mir_est(z, I_el, thresh_cor, depol)
    rattay_z_constants(z) %overwrite z value
    load("rattay_constants.mat")
    r = sqrt(x.^2 + z^2); %distance from electrode to node (cm)
    V_e = rho_e*I_el ./ (4*pi*r);
    Vmir = ((mean(V_e) - V_e(n==0))*1e-3 * thresh_cor) + depol*1e-3;
end
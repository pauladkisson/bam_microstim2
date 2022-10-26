%Paul Adkisson
%4/12/21
%Purpose: Simulates Rattay Cable Equation
%   Inputs:
%       Electrode Current I_el,
%       tspan,
%       initial conditions y0(k) = [V0(k), n0(k), m0(k), h0(k)]
%   Outputs:
%       time t, 3D matrix of state variables y = [V(k), ]
function [t, y, V_e] = rattayrun(y0, I_el, I_int, tstart_int, tend_int, ...
                           tstart_ext, tflip_ext, tend_ext, tmax, tstep, show_progress)
        load("rattay_constants.mat")
        ramp_dur = 10*1e-3; %10us ramp
        t = (0:tstep:tmax)';
        V_e = rho_e*I_el ./ (4*pi*r);
        y = zeros(size(t, 1), N, size(y0, 2));
        y(1, :, :) = y0;
        
        if show_progress
            wait = waitbar(1 / length(t), compose("rattayrun progress: %f %%", 1 / length(t)));
        end
        for i = 2:length(t)
            x = reshape(y(i-1, :, :), N, 4);
            if t(i) >= tstart_int && t(i) <= tend_int
                I_int_ = I_int;
            else
                I_int_ = 0;
            end
            if t(i) >= tstart_ext && t(i)<tflip_ext && t(i) <= tend_ext
                %ramp = min((t(i) - tstart_ext) / ramp_dur, 1); %ramps 0 --> 1
                ramp = 1; %no ramp
                V_e_ = V_e*ramp;
            elseif t(i) >= tstart_ext && t(i)>=tflip_ext && t(i) <= tend_ext
                %ramp = min((t(i) - (tflip_ext+ramp_dur)) / ramp_dur, 1); %ramps -1 --> 1
                ramp = 1;
                V_e_ = -V_e * ramp;
            elseif t(i) > tend_ext
                %ramp = 1 - min((t(i) - tend_ext) / ramp_dur, 1); %ramps 1 --> 0
                ramp = 0;
                V_e_ = -V_e * ramp;
            else
                V_e_ = zeros(size(V_e));
            end
            %RK4 Method
            k1 = rattay_odes(x, V_e_, I_int_);
            k2 = rattay_odes(x+tstep*k1/2, V_e_, I_int_);
            k3 = rattay_odes(x+tstep*k2/2, V_e_, I_int_);
            k4 = rattay_odes(x+tstep*k3, V_e_, I_int_);
            dydt = 1/6*(k1 + 2*k2 + 2*k3 + k4);
            dydt = reshape(dydt, size(y(i, :, :)));
            y(i, :, :) = y(i-1, :, :) + tstep*dydt;
            
            if show_progress
                waitbar(i / length(t), wait, compose("rattayrun progress: %f %%", i / length(t)));
            end
        end
        if show_progress
            close(wait)
        end
end
function [Pin,FQ] = total_power_in(R_star,T_star,r_orb,R_grain,fnum)
%%
% total_power_in
%
% Calculate the total power absorbed by the debris ring.
%
% INPUT
% R_star - [double] body radius of star [R_sun]
% T_star - [double] temperature of star [K]
% r_orb  - [double] orbital radius of secondary object [AU]
% R_grain - [double] average radius of grain [micron]
%
% OUTPUT
% Pin - [double] power absorbed
%
% AUTHOR
% Chun-Yi Wu

global Rdebs lambdaQs Qabs vQs

%% Calculation
% get flux density
Fvs = flux_density_mks(R_star,T_star,r_orb,lambdaQs,0);

% determine with set of Qabs values to use
lambdaQi = 22;
for ( i = 1 : 21 )
    if ( R_grain*1e-6 <= Rdebs(i) )
        lambdaQi = i;
        break;
    end
end

if ( lambdaQi < 22 )
    Qvs = Qabs(:,lambdaQi);
else
    Qvs = ones(size(vQs));
end

FQ = Fvs .* Qvs;

Pin = pi * (R_grain*1e-6)^2 * trapz(vQs,FQ);

%% Plotting (if requested)
if ( fnum > 0 )
    figure(fnum); clf(fnum);
    plot(lambdaQs*1e6,Fvs*1e26,'-');
    grid on; hold on;
    plot(lambdaQs*1e6,FQ*1e26,'-');
    xlabel('wavelength [\mum]');
    ylabel('flux density [Jy]');
    title('dust ring input');
    legend('flux from star','energy absorbed');
    axis([0,3.5,0,max(Fvs*1e26)*1.1]);
end
end
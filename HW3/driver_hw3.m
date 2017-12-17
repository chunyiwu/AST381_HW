clear; clc; close all;

global Rdebs lambdaQs Qabs vQs
readQtable()


R_star = 1.842;
T_star = 8590;
lambdas = linspace(0.1,3.5,1001) * 1e-6;
rs_orb = [10,130];
Rs_grain = [1e-1,1,10,1000];

R_sun = 6.957e8;   % solar radius [m]
M_sun = 1.98855e30; % solar mass [kg];
AU = 1.496e+11; % astronomical unit [m]
k = 1.38064852e-23; % Boltzmann constant [J/K]
h = 6.62607004e-34; % Planck constant [J-s]
c = 2.99792458e+08; % speed of light in vacuum [m/s]
G = 6.67384e-11; % universal gravitational constant

M = 1.92 * M_sun;

%% part 1
fprintf('\n\n\nPart 1\n');
figure(1); clf(1); grid on;
for ( ir = 1 : 2 )
    [Fvs] = flux_density_mks(R_star,T_star,rs_orb(ir),lambdas,0);
    semilogy(lambdas*1e6,Fvs*1e26,'-'); hold on;
end
xlabel('wavelength [\mum]');
ylabel('flux density [Jy]');
title('Stellar output at distance');
legend('10 AU','130 AU');

%% part 2
fprintf('\n\n\nPart 2\n');
for ( ir = 1 : 2 )
    fprintf('orbital radius: %g AU\n',rs_orb(ir));
    for ( iR = 1 : 4 )
        Pin = total_power_in(R_star,T_star,rs_orb(ir),Rs_grain(iR),0);
        fprintf('\tGrain size: %5.0e micron, P = %e W\n', Rs_grain(iR), Pin);
    end
end

%% part 3
fprintf('\n\n\nPart 3\n');
for ( ir = 1 : 2 )
    figure(300+ir); clf(300+ir); grid on;
    leg = cell(1,4);
    fprintf('orbital radius: %g AU\n',rs_orb(ir));
    for ( iR = 1 : 4 )
        Pin = total_power_in(R_star,T_star,rs_orb(ir),Rs_grain(iR),0);
        [Teq,BQ] = dust_spectrum(Pin,R_star,rs_orb(ir),Rs_grain(iR),0);
        fprintf('\tGrain size: %5.0e micron, Teq = %f K\n', Rs_grain(iR), Teq);
        
        loglog(lambdaQs*1e6,BQ*1e26,'-'); hold on;
        leg{iR} = sprintf('%5.0e um - %.1f K',Rs_grain(iR),Teq);
    end
    xlabel('wavelength [\mum]');
    ylabel('flux density [Jy]');
    title(sprintf('dust ring output - %g AU',rs_orb(ir)));
    legend(leg,'Location','Best');
    axis([1,1e3,1e-5,1e15]);
end

%% part 4
fprintf('\n\n\nPart 4\n');
for ( ir = 1 : 2 )
    fprintf('orbital radius: %g AU\n',rs_orb(ir));
    for ( iR = 1 : 4 )
        Pin = total_power_in(R_star,T_star,rs_orb(ir),Rs_grain(iR),0);
        fprintf('\tGrain size: %5.0e micron\n', Rs_grain(iR));
        
        F_RP = Pin/c;
        F_PR = Pin * sqrt(G*M/(rs_orb(ir)*AU)) / c^2;
        
        fprintf('\t\tRad press: %.3e N\n',F_RP);
        fprintf('\t\tPoynt-Rob: %.3e N\n',F_PR);
        
        Fnet = F_RP - F_PR;
    end
end

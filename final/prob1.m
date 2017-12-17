clear; clc; close all;

as = logspace(8,12,1001);
Ms_rv = 2e22 * as;
Ms_tr = ones(size(as)) * 3e23;
Ms_am = 3e44 ./ as;

figure(1); clf(1); grid on;
loglog(as,Ms_rv); hold on;
loglog(as,Ms_tr);
loglog(as,Ms_am);
legend('RV','transit','astrometry','Location','Best');
xlabel('a [m]');
ylabel('M [kg]');
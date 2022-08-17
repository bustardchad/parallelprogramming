%CS759, Chad Bustard, Final Project
%Make plot of times for sequential and parallel cases

%sequential:
N = [10 20 50 100 200 500];
numIters = [18  41  123 207 276 328];
time = [0.409   1.22    22.7    159.4  843 6307];

loglog(N,time,'bo','LineWidth',3);
ylabel('Time (ms)','FontSize',24);
xlabel('Nx, Ny','FontSize',24);
legend('Sequential');
set(gca,'FontSize',20);
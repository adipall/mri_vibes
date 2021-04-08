density = [1400 3100 5150 7200 9200 11000];
disp = [2.1425 2.1153 2.065 2.017 1.9723 1.934];
disp = disp*1e-4;
figure(); clf;
plot(density,disp,'ko-')
xlabel('density kg/m^3')
ylabel('displacement (m)')
%%
beam_rad=[5 10 10 5 2 5];
beam_rad=beam_rad*1e-3;
ball_rad=[6.5 8 6.5 8 6.5 4];
ball_rad=ball_rad*1e-3;
disp = [2.170 2.125 2.14 2.146 2.175 2.187];
disp = disp*1e-4;
figure(); clf;
plot3(beam_rad,ball_rad,disp,'r*')
hold on
plot3([10e-3 2e-3],[4e-3 8e-3],[min(disp) min(disp)],'k--')
plot3([2e-3 10e-3],[8e-3 4e-3],[max(disp) max(disp)],'k--')
% zline(max(disp))
xlabel('beam radius (m)')
ylabel('mass radius (m)')
zlabel('displacement (m)')

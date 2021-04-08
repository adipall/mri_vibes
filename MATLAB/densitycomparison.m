density = [1400 3100 5150 7200 9200 11000];
disp = [2.1425 2.1153 2.065 2.017 1.9723 1.934];
disp = disp*1e-4;
figure(); clf;
plot(density,disp,'ko-')
xlabel('density kg/m^3')
ylabel('displacement (m)')
% Simulation Documentation Script for Demo Day

close all
model_NR = [20 30 40 50 60;
            -0.0394297 -0.0396946 -0.0397873 -0.0398303 -0.0398536];
    
model_OG = [20 30 40 50 60;
            -0.0376659 -0.0389107 -0.0393464 -0.039548 -0.0396576];
    
model_OP = [20 30 40 50 60;
            -0.0385851 -0.0389889 -0.0393904 -0.0395762 -0.0396771];
        
Freq_20 = [ 1 2 3;
    -0.0394297  -0.0376659 -0.0385851];

Freq_40 = [ 1 2 3;
           -0.0397873 -0.0393464 -0.0393904];

        
%% Plot

figure; clf
plot(model_NR(1,:), model_NR(2,:), 'k*-')
hold on
plot(model_OG(1,:), model_OG(2,:), 'r*-')
plot(model_OP(1,:), model_OP(2,:), 'b*-')
grid on
legend('Model Without Resonator', 'Iteration-1 Resonator', 'Iteration-2 Resonator', 'interpreter', 'latex')
xlabel('Applied Frequency (Hz)', 'interpreter', 'latex')
ylabel({'Negative Z-Displacement at', 'Top Surface of Model (mm)'}, 'interpreter', 'latex')
title('Top Surface Displacement at Varying Applied Frequencies: Three Models', 'interpreter', 'latex')
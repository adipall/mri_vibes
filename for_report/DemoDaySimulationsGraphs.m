% Simulation Documentation Script for Demo Day

close all



model_NR = [20 30 40 50 60;
            1-0.0394297 1-0.0396946 1-0.0397873 1-0.0398303 1-0.0398536];
        
model_NR = [model_NR(1,:);
            model_NR(2,:)/24.3];
        
model_OG = [20 30 40 50 60;
            1-0.0376659 1-0.0389107 1-0.0393464 1-0.039548 1-0.0396576];

model_OG = [model_OG(1,:);
            model_OG(2,:)/24.3];        
    
model_OP = [20 30 40 50 60;
            1-0.0385851 1-0.0389889 1-0.0393904 1-0.0395762 1-0.0396771];
        
model_OP = [model_OP(1,:);
            model_OP(2,:)/24.3];
        

        
%% Plot

figure; clf
plot(model_NR(1,:), model_NR(2,:), 'k*-')
hold on
plot(model_OG(1,:), model_OG(2,:), 'r*-')
plot(model_OP(1,:), model_OP(2,:), 'b*-')
grid on
xlim([0 60])
legend('Model Without Resonator', 'Iteration-1 Resonator', 'Iteration-2 Resonator', 'interpreter', 'latex')
xlabel('Applied Frequency (Hz)', 'interpreter', 'latex')
ylabel({'Negative Z-Displacement at', 'Top Surface of Model (mm)'}, 'interpreter', 'latex')
title({'Top Surface Displacement (mm) at Varying', 'Applied Frequencies: Three Models'}, 'interpreter', 'latex')
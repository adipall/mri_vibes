%% physical parameters (Cello A-string)
tension = 133.447; % N (30 lbf)
len = 0.69; % m
freq = 220; % Hz
density = tension/(2*len*freq)^2; % kg / m
c = sqrt(tension/density); % transverse wave speed

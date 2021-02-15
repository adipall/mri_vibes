% psd_driver
fx = 'simple-frfx.frq';
fy = 'simple-frfx.frq'; % = fx
fz = 'simple-frfx.frq'; % = fy
Sxx = [10 .001;1e4 .001];
Syy = [10 0 ;1e4 0];
Szz = [10 0;1e4 0];
block_id  = 201;
relative_disp_psd(fx,fy,fz,Sxx,Syy,Szz,block_id);
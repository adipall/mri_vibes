% function [xrms,yrms,zrms,xrotrms,yrotrms,zrotrms] = relative_disp_psd(frfx,frfy,frfz,Sxx,Syy,Szz,block_id)
function [xrms,yrms,zrms,xrotrms,yrotrms,zrotrms] = relative_disp_psd(frfx,frfy,frfz,Sxx,Syy,Szz,block_id)
    %RELATIVE_DISP_PSD Summary of this function goes here
    %   Detailed explanation goes here

    fx = exo_rd(frfx);
    fy = exo_rd(frfy);
    fz = exo_rd(frfz);

    freq = fx.Time;  % "time" = frequency,  size( freq ) = number time steps = 400

    [node_index_list] = group_j2g_nodes(fx,block_id);

    [numpairs,~] = size(node_index_list);  % There are 4 pairs.  What is a pair?
    % The conmass block, 4th block, id 500,  has 4 bars
    [frfx_reldispx,frfx_reldispy,frfx_reldispz,frfx_reldispxrot,frfx_reldispyrot,frfx_reldispzrot] = get_left_frf_minus_right_frf(fx,node_index_list);
    [frfy_reldispx,frfy_reldispy,frfy_reldispz,frfy_reldispxrot,frfy_reldispyrot,frfy_reldispzrot] = get_left_frf_minus_right_frf(fy,node_index_list);
    [frfz_reldispx,frfz_reldispy,frfz_reldispz,frfz_reldispxrot,frfz_reldispyrot,frfz_reldispzrot] = get_left_frf_minus_right_frf(fz,node_index_list);
    Sxxi = log_interp(Sxx(:,1),Sxx(:,2),freq);
    Syyi = log_interp(Syy(:,1),Syy(:,2),freq);
    Szzi = log_interp(Szz(:,1),Szz(:,2),freq);
    Smat = zeros(6,6,length(freq)); %3 translational and 3 rotational input DOFs at seismic mass
    % Is Smat diagonal?
    Smat(1,1,:) = Sxxi; % uncorrelated inputs for the moment...could be generalized
    Smat(2,2,:) = Syyi;
    Smat(3,3,:) = Szzi; % doesn't consider rotational inputs at the moment
    for j = 1:numpairs
        for i=1:length(freq);
            % This is a 6 x 6 matrix
            Htf = [frfx_reldispx(j,i),    frfy_reldispx(j,i),    frfz_reldispx(j,i), zeros(1,3);
                   frfx_reldispy(j,i),    frfy_reldispy(j,i),    frfz_reldispy(j,i), zeros(1,3);
                   frfx_reldispz(j,i),    frfy_reldispz(j,i),    frfz_reldispz(j,i), zeros(1,3);
                   frfx_reldispxrot(j,i), frfy_reldispxrot(j,i), frfz_reldispxrot(j,i), zeros(1,3);
                   frfx_reldispyrot(j,i), frfy_reldispyrot(j,i), frfz_reldispyrot(j,i), zeros(1,3);
                   frfx_reldispzrot(j,i), frfy_reldispzrot(j,i), frfz_reldispzrot(j,i), zeros(1,3)];
            G(:,:,i) = Htf*Smat(:,:,i)*Htf';
        end
        xrms(j) = calc_grms(freq,real(G(1,1,:)));
        yrms(j) = calc_grms(freq,real(G(2,2,:)));
        zrms(j) = calc_grms(freq,real(G(3,3,:)));
        xrotrms(j) = calc_grms(freq,real(G(4,4,:)));
        yrotrms(j) = calc_grms(freq,real(G(5,5,:)));
        zrotrms(j) = calc_grms(freq,real(G(6,6,:)));
    end
    disp('**** Relative displacement at joint2g elements ****');
     for j=1:numpairs

        fprintf('#%i: Coord: (%d, %d, %d)\n',j,fx.Nodes.Coordinates(node_index_list(j,1),1),fx.Nodes.Coordinates(node_index_list(j,1),2),fx.Nodes.Coordinates(node_index_list(j,1),3));
        fprintf('\t X RMS: %5.3e\n', xrms(j));
        fprintf('\t Y RMS: %5.3e\n', yrms(j));
        fprintf('\t Z RMS: %5.3e\n', zrms(j));
%         fprintf('\t X ROT RMS: %5.3e\n', xrotrms(j));
%         fprintf('\t Y ROT RMS: %5.3e\n', yrotrms(j));
%         fprintf('\t Z ROT RMS: %5.3e\n', zrotrms(j));
        fprintf('\n\n');
    end

end

function [node_index_list] = group_j2g_nodes(f,block_id)
    idx = f.Blocks.id2idx(block_id);
    node_index_list = f.Blocks(idx).Connectivity;
  
end

function [frf_reldispx,frf_reldispy,frf_reldispz,frf_reldisprotx,frf_reldisproty,frf_reldisprotz] = get_left_frf_minus_right_frf(f,node_index_list)
    [nr,~] = size(node_index_list);
    left = 1;
    right = 2;
    for i=1:nr
        [frf_dispx1,frf_dispy1,frf_dispz1,frf_disprotx1,frf_disproty1,frf_disprotz1] = extract_either_left_or_right_frf(f,node_index_list(i,left));
        [frf_dispx2,frf_dispy2,frf_dispz2,frf_disprotx2,frf_disproty2,frf_disprotz2] = extract_either_left_or_right_frf(f,node_index_list(i,right));
            frf_reldispx(i,:) = frf_dispx1 -frf_dispx2;     % hold to your hats
            frf_reldispy(i,:) = frf_dispy1 -frf_dispy2;
            frf_reldispz(i,:) = frf_dispz1 -frf_dispz2;
            frf_reldisprotx(i,:) = frf_disprotx1 -frf_disprotx2;
            frf_reldisproty(i,:) = frf_disproty1 -frf_disproty2;
            frf_reldisprotz(i,:) = frf_disprotz1 -frf_disprotz2;
    end
end

function fieldIndex = isNodalVar(  mesh,  aName )

  numberFields = size( mesh.NodalVars, 2);
  for i = 1:numberFields,
      
      
  end
end


function [frf_dispx,frf_dispy,frf_dispz,frf_disprotx,frf_disproty,frf_disprotz] = extract_either_left_or_right_frf(f,node_index)

    left = 1;
    right = 2;
    
    

    vid1 = f.NodalVars.isNodalVar('DispX');
    vid2 = f.NodalVars.isNodalVar('ImagDispX');
    frf_dispx = f.NodalVars(vid1).Data(node_index,:) + 1j.*f.NodalVars(vid2).Data(node_index,:);


    vid1 = f.NodalVars.isNodalVar('DispY');
    vid2 = f.NodalVars.isNodalVar('ImagDispY');
    frf_dispy = f.NodalVars(vid1).Data(node_index,:) + 1j.*f.NodalVars(vid2).Data(node_index,:);


    vid1 = f.NodalVars.isNodalVar('DispZ');
    vid2 = f.NodalVars.isNodalVar('ImagDispZ');
    frf_dispz = f.NodalVars(vid1).Data(node_index,:) + 1j.*f.NodalVars(vid2).Data(node_index,:);

    
    vid1 = f.NodalVars.isNodalVar('RotX');
    vid2 = f.NodalVars.isNodalVar('ImagRotX');
    frf_disprotx = f.NodalVars(vid1).Data(node_index,:) + 1j.*f.NodalVars(vid2).Data(node_index,:);


    vid1 = f.NodalVars.isNodalVar('RotY');
    vid2 = f.NodalVars.isNodalVar('ImagRotY');
    frf_disproty = f.NodalVars(vid1).Data(node_index,:) + 1j.*f.NodalVars(vid2).Data(node_index,:);


    vid1 = f.NodalVars.isNodalVar('RotZ');
    vid2 = f.NodalVars.isNodalVar('ImagRotZ');
    frf_disprotz = f.NodalVars(vid1).Data(node_index,:) + 1j.*f.NodalVars(vid2).Data(node_index,:);
    
end
function rmsval = calc_grms(freq,psd)
    numFrequencies = length(freq);
    num_locations = size(psd,1);
    for j=1:num_locations    
        val = 0;
        for i=1:numFrequencies-1
            psdl = psd(j,i);
            psdh = psd(j,i+1);
            fl = freq(i);
            fh = freq(i+1);
            numoct = log10(fh/fl)/log10(2);
            db = 10*log10(psdh/psdl);
            mo = 10*log10(2);
            m = db/numoct;
            if (m == -mo)
                temp = psdl * fl *log(fh/fl);
            else
                temp = mo*(psdh/(mo+m))*(    fh-fl*(fl/fh)^(m/mo)    );
            end
            val = val + temp;
        end
        rmsval(j) = sqrt(val);
    end
end
function val_new = log_interp(f,val,fnew)  % f = [initial,final],  val=[v_initial,v_final]
    val_log = 10*log10(val);   %    y = 10 log(10 x),   x = 10^(y/10)
    val_log_new = interp1(f,val_log,fnew);  % linear interpolation
    val_new = 10.^(val_log_new/10);
end

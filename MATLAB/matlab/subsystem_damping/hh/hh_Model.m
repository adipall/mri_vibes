function [K,M,V,freq,vi_resp,vi_inp,dof_inp,dof_resp]=hh_Model()
    %if exist('OnehexModel.mat','file')
     %   load OnehexModel
    %else
        ak=Kssr; %Kaa;
        am=Mssr; %Maa;
        M=am+am'-diag(diag(am));
        K=ak+ak'-diag(diag(ak));
        fexo=exo_get('hh-out.exo');
%         V=zeros(8*size(fexo.Nodes.Coordinates,1),length(fexo.Time));
%         V(1:3:end,:)=fexo.NodalVars(1).Data;
%         V(2:3:end,:)=fexo.NodalVars(2).Data;
%         V(3:3:end,:)=fexo.NodalVars(3).Data;
        %V(:,1) = hh_disp_a1;
        freq=fexo.Time;
        for i = 1:length(freq)
            eval(sprintf('V(:,%d)=hh_disp_a%d',i,i));
        end
        %freq(1:6)=0; %%clean up the rigid modes 
        % I=dsearchn(fexo.Nodes.Coordinates,[-5 -.05 -.5]);
        % II=find(abs(fexo.Nodes.Coordinates(:,1)+5)<1e-10);
        % vi_resp=I;
        % vi_inp=II
        %% Figure out DOF used in model after boundary conditions.
        fetimap_a = FetiMap_a;
        nodes_v = 1:length(fexo.Nodes.Coordinates);
        count_bd = 0;
        for i = 1:length(fexo.Nodes.Coordinates)
            if (fexo.Nodes.Coordinates(i,1) == 0.0)
                count_bd = count_bd + 1;
                bound_nd(count_bd) = i;
            end
        end
        nodes_v(bound_nd) = [];
        count_vi =0;
        for i = 1:length(fexo.Nodes.Coordinates)
            if (fexo.Nodes.Coordinates(i,1) == 2.0)
                count_vi = count_vi + 1;
                node_inp(count_vi) = i; % node for input
            end
        end
        % input DOF
        for i = 1:length(node_inp)
            vi_inp(i)=find(nodes_v==node_inp(i));
            dof_inp(i) = fetimap_a(node_inp(i)*8-6); % -6 puts the load in the y-direction.
        end
        %vi_inp = 1;
       node_resp = 63; % node to look at for comparison
       vi_resp = find(nodes_v==node_resp);
       dof_resp = fetimap_a(node_resp*8-6); % -6 looks at the dof in the y-direction
       
          
        % save OnehexModel K M V freq vi_resp vi_inp
    %end
    
    
end
    
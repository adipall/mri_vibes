function [K,M,V,freq,vi_resp,vi_inp]=Beam_Model()
    if exist('BeamModel.mat','file')
        load BeamModel
    else
        ak=Kaa;
        am=Maa;
        M=am+am'-diag(diag(am));
        K=ak+ak'-diag(diag(ak));
        fexo=exo_rd('beam-out.g');
        V=zeros(3*size(fexo.Nodes.Coordinates,1),length(fexo.Time));
        V(1:3:end,:)=fexo.NodalVars(1).Data;
        V(2:3:end,:)=fexo.NodalVars(2).Data;
        V(3:3:end,:)=fexo.NodalVars(3).Data;
        freq=fexo.Time;
        freq(1:6)=0; %%clean up the rigid modes 
        I=dsearchn(fexo.Nodes.Coordinates,[-5 -.05 -.5]);
        II=find(abs(fexo.Nodes.Coordinates(:,1)+5)<1e-10);
        vi_resp=I;
        vi_inp=II;
        save BeamModel K M V freq vi_resp vi_inp
    end
end
    
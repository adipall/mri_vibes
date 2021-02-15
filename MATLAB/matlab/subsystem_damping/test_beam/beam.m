%function [t,x,v,a,iter,force]=beam()
    %
    %  [t,x,v,a,iter,params]=beam()
    %
    dbg=0;  % 1 for just making matrices
    %%%
    % Read the model
    %%%
    [K,M,V,freq,vi_resp,vi_inp]=Beam_Model;
    dir=2;  %% drive in the y direction
    
    dr = '/home/mross/projects/salinas.d/joints_dissipation.d/subsystem_damping_matlab';
    path(path,dr);
    
    %%%%%
    % Set up the Joint Like Distributed Damping
    %%%%%
    if 0, %1=viscous damper, 0=Iwan,
        %% viscous damper
        disp('Just Viscous Damper')
        for i=1:length(freq),
            params(i).type='lindamp';
            params(i).data=[2*freq(i)*2*pi*0.02 0];
        end
    else
        disp('Iwan')
        %% Iwan damper
        Kt=0.1*freq.*freq*4*pi*pi;   %  10% of modal stiffness
        Fs=0.1*ones(length(freq),1);  % about 10 x max seen in viscous run
        beta=10*ones(length(freq),1);  % just a nice number
        chi=2.6*ones(length(freq),1)-3;  % 2.6 from theory for material damping
        %%% derived quantities
        phi_max=(Fs./Kt).*(1+beta)./(beta+(chi+1)./(chi+2));
        R=Fs.*(chi+1)./(phi_max.^(chi+2))./(beta+(chi+1)./(chi+2));
        S=(Fs./phi_max).*beta./(beta+(chi+1)./(chi+2));
        
        %% "Viscous" for the rigid modes
        for i=1:6,
            params(i).type='lindamp';
            params(i).data=[2*freq(i)*2*pi*0.02 0];
        end
        
        %% Iwan for the flexible modes
        for i=7:length(freq),
            params(i).type='iwan';
            params(i).data=[chi(i) phi_max(i) R(i) S(i)];
        end
    end
    n=size(M,1);
    %Cm=sparse(n,n,.01);
    %Cm=M*V*diag(2*0.02*freq.*2*pi)*V'*M;
    
    %%%
    % Set up some time integration parameters
    %%%
    nsteps=2000;
    dt=1e-4;
    num_iter=10;
    tolerance=1e-6;
    g=100;
    
    %%%
    % Allocate the displacement vectors
    %%%
    x=zeros(length(vi_resp),nsteps);
    iter=zeros(1,nsteps);
    v=zeros(length(vi_resp),nsteps);
    a=zeros(length(vi_resp),nsteps);
    
    
    % initial conditions are zero
    x0=zeros(n,1);x(:,1)=x0(3*(vi_resp-1)+dir);
    v0=zeros(n,1);v(:,1)=v0(3*(vi_resp-1)+dir);
    
    %%
    % Set up the impulse function
    %%
    func=[0 0;3e-4 1;6e-3 0;1e6 0];
    t=0:dt:dt*(nsteps-1);
    f=g*interp1q(func(:,1),func(:,2),t')';
    
    F=spalloc(n,length(t),length(vi_inp)*length(f));
    F(3*(vi_inp-1)+dir,:)=repmat(f,length(vi_inp),1); % z direction input
    FintM=zeros(length(freq),length(t));
    if dbg==1,
        return;
    end
    alpha=[0 .25 0 .5]; %% this is an undamped integrator.  Replace .9 below with alpha
    
    %%%%
    % instatiate the integrator
    %%%%%
    Tobj=ImplicitIntegrator(x0,v0,dt,.9);
    %%%%%
    % instatiate the distributed damping
    %%%%%
    JLDDobj=DistributedDamping(params,V,freq,M*V);
    
    %%%
    % Provide entertainment for the analyst while the simulation is running
    %%%
    i10=1;
    h=waitbar(0,sprintf('Performing Time Integration... Step %d/%d',i10,nsteps),'Name','Time Integration Waitbar');
    
    %% loop over timesteps, i=1 is the initial condition
    for i=2:nsteps,
        %%% Estimate the displacements for the timestep
        Tobj=Tobj.estimateDispl;
        %% iterate until it converges
        for j=1:num_iter,
            %% get a current estimate of the force
            [JLDDobj,force]=JLDDobj.updateForce(Tobj.X,Tobj.Vel);
            %% integrate the displacements
            [Tobj]=Tobj.GamIntegrator(M,sparse(0),K,full(F(:,i)),force+K*Tobj.X);
            %% check for convergence
            if abs(Tobj.Tol)<tolerance,
                break;
            end
        end
        %% save the number of iterations
        iter(i)=j;
        %% if too many iterations were used, stop.
        if j==num_iter,
            disp(['Num_iter exceeded ',num2str(i)]);
            break;
        end
        %% update states in integrator and distributed damping
        Tobj=Tobj.updateStates;
        JLDDobj=JLDDobj.updateStates;
        FintM(:,i)=V'*M*force;
        %% save the measurement location
        x(:,i)=Tobj.X(3*(vi_resp-1)+dir);
        v(:,i)=Tobj.Vel(3*(vi_resp-1)+dir);
        a(:,i)=Tobj.Accel(3*(vi_resp-1)+dir);
        %% Keep the analyst entertained
        waitbar(i/nsteps,h,sprintf('Performing Time Integration... Step %d/%d',i10,nsteps))
        if ~mod(i,100)
            i10=i;
        end
    end
    close(h)
    %% head home.
    
%end

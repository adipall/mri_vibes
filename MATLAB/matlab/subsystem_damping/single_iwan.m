function [t,x,v,a,iter,FintM,xm,params]=single_iwan()
    %
    %  [t,x,v,a,iter,params]=beam()
    %
    xm=[];
    dbg=0;  % 1 for just making matrices
    %%%
    % Read the model
    %%%
    M=100;
    
    %%%%%
    % Set up the Joint Like Distributed Damping
    %%%%%
    if 0,
        %% viscous damper
        params.type='lindamp';
        params.data=[2*freq*2*pi*0.02 0];
    end
    if 1,
        params.type='iwan';
        %params.data=[chi(i) phi_max(i) R(i) S(i)];
        %params.data=[-.6 7.1e-5 4.9e7 8.8e6];
        params.data=[-0.802 7.6e-8 1.01e7 9.88e5];
    end
    
    %%%
    % Set up some time integration parameters
    %%%
    nsteps=10000;
    dt=1e-5;
    num_iter=2000;
    tolerance=1e-6;
    
    %%%
    % Allocate the displacement vectors
    %%%
    x=zeros(1,nsteps);
    iter=zeros(1,nsteps);
    v=zeros(1,nsteps);
    a=zeros(1,nsteps);
    FintM=zeros(1,nsteps);
    
    % initial conditions are zero
    x0=zeros(1,1);x(:,1)=x0;
    v0=zeros(1,1);v(:,1)=v0;
    
    %%
    % Set up the impulse function
    %%
    g=.1;
    func=[0 0;3e-4 1;6e-3 0;1e6 0];
    t=0:dt:dt*(nsteps-1);
    F=g*interp1q(func(:,1),func(:,2),t')';
    %F=sin(2*pi*200*t);
    if dbg==1,
        return;
    end
    alpha=[0 .25 0 .5]; %% this is an undamped integrator.  Replace .9 below with alpha
   
    %%%%
    % instatiate the integrator
    %%%%%
    Tobj=ImplicitIntegrator(x0,v0,dt,alpha);
    %%%%%
    % instatiate the Iwan element
    %%%%%
    Iwobj=Iwan(params.data);
    Iwobj.K0
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
            [Iwobj,force]=Iwobj.updateForce(Tobj.X,Tobj.Vel);
            %% integrate the displacements
            [Tobj]=Tobj.GamIntegrator(M,sparse(0),sparse(0),full(F(:,i)),force);
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
        Iwobj=Iwobj.updateStates;
        %% save the measurement location
        x(:,i)=Tobj.X;
        v(:,i)=Tobj.Vel;
        a(:,i)=Tobj.Accel;
        FintM(:,i)=force;
        %% Keep the analyst entertained
        waitbar(i/nsteps,h,sprintf('Performing Time Integration... Step %d/%d',i10,nsteps))
        if ~mod(i,100)
            i10=i;
        end
    end
    close(h)
    %% head home.
    datestr(now)
end

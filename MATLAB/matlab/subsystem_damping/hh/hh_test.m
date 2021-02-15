%function [t,x,v,a,iter,force]=hh_test()
    %
    %  [t,x,v,a,iter,params]=beam()
    %
    clear all
    close all
    
    dr = '/home/mross/projects/salinas.d/joints_dissipation.d/subsystem_damping_matlab';
    path(path,dr);
    dr = '/home/mross/matlab.d/toolbox.d/exodus';
    path(path,dr);
    
    %% ---- Parameteres ----------------------------------------------------
    dir=2;  %% drive in the y direction
    
    
    %%%
    % Set up some time integration parameters
    %%%
    nsteps=2000;
    dt=5e-6;
    num_iter=20;
    tolerance=1e-6;
    g=1e3; %  scale factor for the force
    
    %% --------------------------------------------------------------------
    
    dbg=0;  % 1 for just making matrices
    %%%
    % Read the model
    %%%
    [K,M,V,freq,vi_resp,vi_inp,dof_inp,dof_resp]=hh_Model;
    
    %%%%%
    % Set up the Joint Like Distributed Damping
    %%%%%
    if 0, %1=viscous damper, 0=Iwan,
        %% viscous damper
        disp('Just Viscous Damper')
        for i=1:length(freq),
            params(i).type='lindamp';
            %params(i).data=[2*freq(i)*2*pi*0.02 0];
            params(i).data=[10000 0];
        end
    else
        disp('Iwan')
        %% Iwan damper
        %Kt=0.1*freq.*freq*4*pi*pi;   %  10% of modal stiffness
        Kt(1:10,1) = 1e7; %5e6; %[1e4;1e4];
        Fs=1e2*ones(length(freq),1); %0.1e-5*ones(length(freq),1);  % about 10 x max seen in viscous run
        beta=0.5*ones(length(freq),1); %10*ones(length(freq),1);  % just a nice number
        chi=2.6*ones(length(freq),1)-3;  % 2.6 from theory for material damping
        %%% derived quantities
        phi_max=(Fs./Kt).*(1+beta)./(beta+(chi+1)./(chi+2));
        R=Fs.*(chi+1)./(phi_max.^(chi+2))./(beta+(chi+1)./(chi+2));
        S=(Fs./phi_max).*beta./(beta+(chi+1)./(chi+2));
        
%         Kt=2.97e6*ones(length(freq),1);   %  10% of modal stiffness
%         Fs=0.1e-5*ones(length(freq),1);  % about 10 x max seen in viscous run
%         beta=0.498*ones(length(freq),1);  % just a nice number
%         chi=-0.802*ones(length(freq),1);  % 2.6 from theory for material damping
%         %%% derived quantities
%         phi_max=(7.6e-8)*ones(length(freq),1);
%         R=1.01e7*ones(length(freq),1);
%         S=9.88e5*ones(length(freq),1);
        
        %% Set mode type
        
        % "Viscous" for the rigid modes
%         for i=1:6,
%             params(i).type='lindamp';
%             params(i).data=[2*freq(i)*2*pi*0.02 0];
%         end
        
%         %% Iwan for the flexible modes
%         for i=1:length(freq), %7:length(freq),
%             params(i).type='iwan';
%             params(i).data=[chi(i) phi_max(i) R(i) S(i)];
%         end
    
        % Only set Mode 1 has iwan parameter
        params(1).type = 'iwan';
        params(1).data=[chi(1) phi_max(1) R(1) S(1)];

%         for i = 2:length(freq)
%             params(i).type='lindamp';
%             params(i).data=[0 0];
%         end 
        
    end
    n=size(M,1);
    %Cm=sparse(n,n,.01);
    %Cm=M*V*diag(2*0.02*freq.*2*pi)*V'*M;
      
    %%%
    % Allocate the displacement vectors
    %%%
%     x=zeros(length(vi_resp),nsteps);
%     iter=zeros(1,nsteps);
%     v=zeros(length(vi_resp),nsteps);
%     a=zeros(length(vi_resp),nsteps);
    x=zeros(length(dof_resp),nsteps);
    iter=zeros(1,nsteps);
    v=zeros(length(dof_resp),nsteps);
    a=zeros(length(dof_resp),nsteps);
%     x2 = zeros(3*length(vi_resp),nsteps);
%     v2 = zeros(3*length(vi_resp),nsteps);
    
    
    % initial conditions are zero
    x0=zeros(n,1);
    %x(:,1)=x0(3*(vi_resp-1)+dir);
    x(:,1)=x0(dof_resp);
    v0=zeros(n,1);
    %v(:,1)=v0(3*(vi_resp-1)+dir);
    v(:,1)=v0(dof_resp);
    force = zeros(size(K,1),1);
    
    %%
    % Set up the impulse function
    %%
    func=[0 0;1e-4 1.0;2e-4 0;1.0 0];
    t=0:dt:dt*(nsteps-1);
    f=g*interp1q(func(:,1),func(:,2),t')';
    
    F=spalloc(n,length(t),length(vi_inp)*length(f));
    %F(3*(vi_inp-1)+dir,:)=repmat(f,length(vi_inp),1); % direction input
    F(dof_inp,:)=repmat(f,length(dof_inp),1); % direction input
    FintM=zeros(length(freq),length(t));
    if dbg==1,
        return;
    end
    alpha=[0 .25 0 .5]; %% this is an undamped integrator.  Replace .9 below with alpha
    %alpha = 0.0; % This matches the Salinas default.
    
    %%%%
    % instatiate the integrator
    %%%%%
    Tobj=ImplicitIntegrator(x0,v0,dt,alpha);
    %%%%%
    % instatiate the distributed damping
    %%%%%
    JLDDobj=DistributedDamping(params,V,freq,M*V);
    
    %%% UNIT TEST PORTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     x2 = dlmread('x2_new.txt');
%     v = zeros(12,1);
%     [JLDDobj,force]=JLDDobj.updateForce(x2,v);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
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
            Force_store(:,i) = force; % + K*Tobj.X;
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
        FintM(:,i)=V'*force;
        % FintM(:,i)=modalforce;
        xm(:,i)=V'*M*Tobj.X;
        % xm(:,i)=xq;
        % FintM(:,i)=V'*force;
        % xm(:,i)=V'*Tobj.X;
        %% save the measurement location
%         x(:,i)=Tobj.X(3*(vi_resp-1)+dir);
%         v(:,i)=Tobj.Vel(3*(vi_resp-1)+dir);
%         a(:,i)=Tobj.Accel(3*(vi_resp-1)+dir);
        x(:,i)=Tobj.X(dof_resp);
        v(:,i)=Tobj.Vel(dof_resp);
        a(:,i)=Tobj.Accel(dof_resp);
%         x2(:,i) = Tobj.X;
%         v2(:,i) = Tobj.Vel;
        %% Keep the analyst entertained
        waitbar(i/nsteps,h,sprintf('Performing Time Integration... Step %d/%d',i10,nsteps))
        if ~mod(i,100)
            i10=i;
        end
    end
    close(h)
    %% head home.
    
    figure(1)
    plot(t,x);
    
    figure(2)
    plot(x,Force_store(dof_resp,:));
    
    figure(3)
    plot(t,v)
    
    %% Compare with Salinas Results
    load hh-two
    figure(4)
    plot(t,x,time,nvar02(63,:),'LineWidth',2)
    legend('Simmermacher','Salinas')
    title('Linear Damper','Fontsize',20)
    xlabel('Time (sec)','Fontsize',18)
    ylabel('Displacement','Fontsize',18)
    set(gca,'Fontsize',16)
    set(gcf,'position',[50 60 800 600])
    
%     diffx = x-nvar02(63,:);
%     figure(5)
%     plot(t,diffx,'LineWidth',2)
%     title('Difference Linear Damper','Fontsize',20)
%     xlabel('Time (sec)','Fontsize',18)
%     ylabel('Difference Salinas-Simmermacher','Fontsize',18)
%     set(gca,'Fontsize',16)
%     set(gcf,'position',[50 60 800 600])  

    figure(6)
    plot(xm(1,:),FintM(1,:))
    
%end
%%

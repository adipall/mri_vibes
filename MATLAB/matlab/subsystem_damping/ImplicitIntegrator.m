classdef ImplicitIntegrator
    %
    % This is the integrator that salinas uses
    %
    properties (Access=public)
        X     % Current Displacement
        Xhat   % estimated displacement used for nonlinear iterations
        Vel     % Current Velocity
        Accel     % Current Acceleration
        Fint      % Current internal force
        Fext   % Current external force
        Dt    % current timestep
        Tol   % convergence tolerance
    end
    properties (Access=private)
        XPrevious      % previous displacement
        VPrevious      % previous velocity
        APrevious      % previous acceleration
        FintPrevious   % previous internal force
        FextPrevious   % previous external force
        %%
        % constants affiliated with the integrator
        %%
        Am     = 0.
        Bn     = 0.25
        Af     = 0.
        Gamman = 0.5
    end
    
    methods
        function obj=ImplicitIntegrator(x0,v0,dt,rho)
            % obj.ImplicitIntegrator(x0,v0,dt,rho)
            %  or
            % obj.ImplicitIntegrator(x0,v0,dt,alpha)
            % where
            %
            % alpha=[am bn af gamman]
            %
            % rho=> constant used to calculate the alphas (same as salinas)
            %
            if nargin>0
                if length(rho)==1,
                    alpha=zeros(4,1);
                    alpha(1)=rho/(1+rho);
                    alpha(2)=(2*rho-1)/(1+rho);
                    alpha(3)=(1-alpha(2)+alpha(1))*(1-alpha(2)+alpha(1))/4;
                    alpha(4)=0.5-alpha(2)+alpha(1);
                else
                    alpha=rho;
                end
                %%
                obj.XPrevious=x0;
                obj.VPrevious=v0;
                obj.APrevious=zeros(size(v0));
                obj.FintPrevious=zeros(size(v0));
                obj.FextPrevious=zeros(size(v0));
                obj.Dt=dt;
                obj.Am=alpha(1);
                obj.Bn=alpha(2);
                obj.Af=alpha(3);
                obj.Gamman=alpha(4);
            end
        end
        function [obj]=GamIntegrator(obj,M,C,K,Fext,Fint)
            %
            obj.Fext=Fext;
            obj.Fint=Fint;
            alpha=[obj.Am obj.Bn obj.Af obj.Gamman]; % define a vector with the constants for speed
            dt=obj.Dt;
            %
            
            A=M*((1-alpha(1))/alpha(2)/dt/dt) + C*((1-alpha(3))*alpha(4)/alpha(2)/dt) ...
                + K*((1-alpha(3)));
            
            Fext_af=(1-alpha(3))*Fext + alpha(3)*obj.FextPrevious;
            %Fint=Fintn+K*(dhat-xn);
            %%%
            obj.Accel=(obj.X-obj.XPrevious-obj.VPrevious*dt)/(alpha(2)*dt*dt) - ...
                (1-2*alpha(2))*obj.APrevious/(2*alpha(2));
            %
            obj.Vel=obj.VPrevious+dt*((1-alpha(4))*obj.APrevious+alpha(4)*obj.Accel);
            %%%
            res=Fext_af-(1-alpha(3))*Fint-alpha(3)*obj.FintPrevious ...
                -C*((1-alpha(3))*obj.Vel+alpha(3)*obj.VPrevious) ...
                -M*((1-alpha(1))*obj.Accel+alpha(1)*obj.APrevious);
            
            if issparse(A),
                obj.Tol=A\res;
            else
                obj.Tol=sparse(A)\res;
            end
            
            
            obj.X=obj.X+obj.Tol;
            obj.Accel=(obj.X-obj.XPrevious-obj.VPrevious*dt)/(alpha(2)*dt*dt) - ...
                ((1-2*alpha(2))/2/alpha(2))*obj.APrevious;
            obj.Vel=obj.VPrevious+dt*((1-alpha(4))*obj.APrevious+alpha(4)*obj.Accel);
        end
        function obj=estimateDispl(obj)
            alpha=[obj.Am obj.Bn obj.Af obj.Gamman]; % define a vector with the constants for speed
            dt=obj.Dt;
            obj.X=obj.XPrevious+dt*obj.VPrevious+dt*dt/2*(1-2*alpha(2))*obj.APrevious;
            %
            obj.Accel=(obj.X-obj.XPrevious-obj.VPrevious*dt)/alpha(2)/dt/dt - ...
                (1-2*alpha(2))*obj.APrevious/2/alpha(2);
            %
            obj.Vel=obj.VPrevious+dt*((1-alpha(4))*obj.APrevious+alpha(4)*obj.Accel);
        end
        function obj=updateStates(obj)
            obj.XPrevious = obj.X;     
            obj.VPrevious = obj.Vel;     
            obj.APrevious = obj.Accel;    
            obj.FintPrevious = obj.Fint; 
            obj.FextPrevious =  obj.Fext;
        end
    end
end
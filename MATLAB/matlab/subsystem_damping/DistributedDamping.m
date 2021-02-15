classdef DistributedDamping
    properties
        Modeshapes=0    % mass normalized
        InvModeshapes
        Frequencies=0   % in Hz
        ModalDampingObj     % vector of nonlinear element objects
        NModes      % number of modes/damping elements
        NDOF        % number of dofs in the subsystem (size(Modeshapes,1))
        Force
    end
    
    methods
        function obj=DistributedDamping(params,V,freq,Vinv)
            %
            %   obj=DistributedDamping(params,V,freq,Vinv)
            %
            %
            % params(i).type ='iwan' or 'lindamp';
            % params(i).data=[chi,phi_max,R,S] for iwan
            %            =[zeta Kt] for lindamp
            %   where i is the ith mode.  so length(params)==size(V,2)
            %
            
            %
            obj.Modeshapes=V
            if nargin<4,
                obj.InvModeshapes=pinv(V)';
            else
                obj.InvModeshapes=Vinv; % Vinv = M*V;
            end
            obj.Frequencies=freq;
            if length(params)~=size(V,2),
                error('number of nonlinear elements and number of modes are not equal');
            end
            obj.NModes=size(V,2);
            obj.NDOF=size(V,1);
            for i=1:length(params),
                switch lower(params(i).type)
                    case 'iwan'
                        nlobj=Iwan(params(i).data);
                    case 'lindamp'
                        nlobj=LinDamp(params(i).data);
                    otherwise
                        error('Unknown nonlinear element')
                end
                obj.ModalDampingObj{i}=nlobj;
            end
        end
        function [obj,force]=updateForce(obj,x,v)
            modalforce=zeros(obj.NModes,1);
            % Transform to modal coordinates
            xq=obj.InvModeshapes'*x;
            vq=obj.InvModeshapes'*v;
            % update the modal forces and subtract off no load stiffness
            % effects
            for i=1:obj.NModes,
                % Stupid overloader screws me up everytime. This update
                % force is with LinDamp.
                [obj.ModalDampingObj{i},modalforce(i)]=obj.ModalDampingObj{i}.updateForce(xq(i),vq(i));
                modalforce(i)=modalforce(i)-obj.ModalDampingObj{i}.K0*xq(i);
            end
            % transform back to physical coordinates
            force=obj.InvModeshapes*modalforce;
            obj.Force=force;
        end
        function [obj]=updateStates(obj)
            for i=1:obj.NModes,
                obj.ModalDampingObj{i}=obj.ModalDampingObj{i}.updateStates;
            end
        end
    end
end
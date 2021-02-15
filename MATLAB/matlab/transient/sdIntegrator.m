%SDINTEGRATOR sd time integration algorithm
classdef sdIntegrator
   properties
      alpha=0;       % mass proportional damping
      beta=0;        % stiffness proportional damping
      rho=1;
      alphaF=0;
      alphaM=0;
      newmarkBeta=0;
      newmarkGamma=0;
      rampRate=0;         % a load
      impulseDuration=0;  % another ...
      impulseMagnitude=0; % ... load
   end

   methods
      function obj = sdIntegrator(mDamp,kDamp,rho)
         if nargin == 3
            if isnumeric(mDamp) && isnumeric(kDamp) && isnumeric(rho)
               obj.alpha = mDamp;
               obj.beta  = kDamp;
               obj.rho   = rho;
               obj.alphaF=rho/(rho+1);
               obj.alphaM=(2*rho-1)/(rho+1);
               x=1-obj.alphaM+obj.alphaF;
               obj.newmarkBeta=0.25*(1-obj.alphaM+obj.alphaF)*x;
               obj.newmarkGamma=0.5-obj.alphaM+obj.alphaF;
            else
               error('Values must be numeric')
            end
         else
            fprintf('ctor %d args\n',nargin);
         end
      end
      function rhs = getRhs(obj,K,M,deltaH,accel,veloc,disp)
          xx = (2*obj.newmarkBeta);
          fCoeffOfAccel=(1-obj.alphaF)*(deltaH*(1-obj.newmarkGamma)+...
              -obj.newmarkGamma*deltaH*(1-2*obj.newmarkBeta)/xx);
          fCoeffOfVel=1-(1-obj.alphaF)*obj.newmarkGamma/obj.newmarkBeta;
          bh = (obj.newmarkBeta*deltaH);
          fCoeffOfD=(1-obj.alphaF)*(-obj.newmarkGamma/bh);
          fcoeff=[fCoeffOfAccel;fCoeffOfVel;fCoeffOfD];
          vector=[accel,veloc,disp]*fcoeff;
          C = obj.alpha*M + obj.beta*K;
          proportionalDampingTerms=C*vector;
          dVec = disp*obj.alphaF;
          Kterms = K*dVec;
          coeffOfAccel=obj.alphaM-(1-obj.alphaM)*(1-2*obj.newmarkBeta)/xx;
          coeffOfVel=-(1-obj.alphaM)/(obj.newmarkBeta*deltaH);
          deltaSq=deltaH*deltaH;
          coeffOfD=-(1-obj.alphaM)/(obj.newmarkBeta*deltaSq);
          coeff = [coeffOfAccel;coeffOfVel;coeffOfD];
          vec= [accel,veloc,disp]*coeff;
          Mterms = M*vec;
          rhs = -Kterms -Mterms -proportionalDampingTerms;
      end % getRhs function
      function dPlus = linearSystem(obj,K,M,deltaH,rhs)
          deltaSq=deltaH*deltaH;
          ms1 = (1-obj.alphaM)/(obj.newmarkBeta*deltaSq);
          bh = (obj.newmarkBeta*deltaH);
          mScale = ms1 + obj.alpha*(1-obj.alphaF)*obj.newmarkGamma/bh;
          K_coeff = 1-obj.alphaF+obj.beta*(1-obj.alphaF)*obj.newmarkGamma/bh;
          C_coeff =(1-obj.alphaF)*obj.newmarkGamma/(obj.newmarkBeta*deltaH);
          nscale = K_coeff/C_coeff;
          C = 0; % What is C?
          conv = C + nscale*K;
          conv=conv*C_coeff;
          sigma=-1/mScale;
          matrix = M - sigma * conv;
          sol=rhs/matrix;
          dPlus = sol/mScale;
      end % linearSystem
      function [D_n,Vel_n,Accel_n]=updateState(obj,deltaH,accel,veloc,disp,dPlus)
        coeffOfAccelForVel=deltaH*((1-obj.newmarkGamma)-obj.newmarkGamma*(1-2*obj.newmarkBeta)/(2*obj.newmarkBeta));
        coeffOfVelForVel=1-obj.newmarkGamma/(obj.newmarkBeta);
        coeffOfDispForVel=obj.newmarkGamma/(obj.newmarkBeta*deltaH);
        coefForVel=[coeffOfAccelForVel;coeffOfVelForVel;coeffOfDispForVel];
        deltaSq = deltaH * deltaH;
        coeffOfDispForAccel=1/(obj.newmarkBeta*deltaSq);
        coeffOfVelForAccel=-1/(obj.newmarkBeta*deltaH);
        coeffOfAccelForAccel=-(1-2*obj.newmarkBeta)/(2*obj.newmarkBeta);
        coefForAccel=[coeffOfAccelForAccel;coeffOfVelForAccel;coeffOfDispForAccel];
        Accel_n = [accel,veloc,dPlus-disp]*coefForAccel;
        Velplus_temp=[accel,veloc,dPlus-disp]*coefForVel;
        Vel_n = Velplus_temp;
        D_n = dPlus;
      end
      function history = integrator(obj,K,M,deltaH,levels,uo,vo)
        disp = uo; veloc=vo; accel =0;
        history(1,:) = [disp,veloc,accel];
        for level = 2:levels,
            rhs = obj.getRhs(K,M,deltaH,accel,veloc,disp);
            rhs=rhs+obj.force(K,level,deltaH);
            dPlus=obj.linearSystem(K,M,deltaH,rhs);
            [Dn,Vely,Acceln]=obj.updateState(deltaH,accel,veloc,disp,dPlus);
            disp = Dn;  veloc=Vely; accel = Acceln;
            history(level,:) =[disp,veloc,accel];
        end
      end % integrator
      function load=force(obj,K,level,deltaH)
        time=level*deltaH;
        load=time*K*obj.rampRate;
        if time<obj.impulseDuration,
            load=load+obj.impulseMagnitude;
        end
      end % integrator
      function disp(td)
         fprintf(1,'a: %g\nb %g\nr %1.5g\naF %g\naM %g\nB %g\nG %g\n',...
            td.alpha,td.beta,td.rho,td.alphaF,td.alphaM,td.newmarkBeta,...
            td.newmarkGamma);
      end % disp
   end % methods
end % class

% I. time integration algorithms
%   sdIntegrator               - sd time integration algorithm
%   generalizedAlpha           - function [d,v,a]=generalizedAlpha
%                                     (n,tf,m,c,k,rampRate,rho,uo,vo)
%   generalizedAlphaParameters - function [alphaF,alphaM,beta,gamma]=
%                                    generalizedAlphaParameters(rho)
% II. tests, examples, scripts to make figures, helper functions
%   testSDSetMatrix            - regression test SD time integration
%   testSDSecondStep           - regression test integrator level 2
%   testSDInitialStep          - regression test integrator level 1
%   testConstantForce          - test solutionConstantForce
%   testPiecewiseConstF        - remove
%   testImpulse                - test solutionImpulse
%   testRampedForce            - test solutionRampedForce
%   testIvp                    - test solutionIvp
%   testAlgoChung              - Chung&Hulbert's Generalized Alpha
%   testTrapezoidal            - test generalizedA(rho=1)
%   testAlphaParameters        - test generalizedAlphaParameters
%   sdTestSuite                - test initial value problems, f=0
%   oneSdTest                  - function [history,u,v,a,uh,vh,ah,time]
%                                    oneSdTest(id)
%   getSuiteParameters         - function [uo,vo,rho,mDamp,kDamp,name,number]
%                                    = getSuiteParameters(id)
%   exampleRamped              - example zero initial conditions
%                                    and ramped force
%   exampleChung               - example Generalized Alpha method
%   figGood                    - trapezoid works sometimes
%   figAo.m                    - initial acceleration
%   figLinAccel.m              - accerations are lousy if rho < 1
%   figRamp.m                  - time dependent forces are lousy
%   getRandomModel             - function [k,c,m] = getRandomModel()
%   getKandM                   - function [k,m] = getKandM(),defualt
%   getRhoH                    - function [rho,h] = getRhoH(),defualt
%   getProportionalDamping     - function [alpha,beta]=
%                                    getProportaionalDamping()
% III. sdof const coef, homogeneous initial conditions
%   solutionConstantForce      - function [u,v] =
%                                    solutionConstantForce(m,c,k,uf,n,tf)
%   solutionPiecewiseConstF    - remove
%   solutionImpulse            - function u = solutionImpulse
%                                    (m,c,k,uforce,T,num,tfinal)
%   solutionRampedForce        - function [u,v,a]=
%                                    solutionRampedForce(m,c,k,v,n,ti,tf)
%   solutionIvp                - function [u,v,a] =
%                                    solutionIvp(m,c,k,uo,vo,n,tf)
%   rampedForce                - function f=
%                                    rampedForce(k,v,tf) f=v*k*time

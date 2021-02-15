% MATHIEU_FUNCTIONS_TOOLBOX
%
% Files
%   Cpm                   - function correlationFactors = Cpm(category,mv,mvv,nmax)
%   dHpm1                 - function y = dHpm1(category,u,q,mv,nmax)
%   dHpm2                 - function y = dHpm2(KF,u,q,mv,nmax)
%   dJpm                  - function y = dJpm(KF,u,q,mv,nmax)
%   dSpm                  - functiony = dSpm(KF,v,mv,nmax)
%   dYpm                  - function y = dYpm(KF,u,q,mv,nmax)
%   eig_Spm               - function [va,mv,vt]=eig_Spm(category,ellipticalParameter)
%   example1_Jpm          - example evaluate radial Mathieu function Jpm at nmax and orders
%   example1_Spm          - example angular Mathieu function at given: q > 0, nmax, and order n
%   example2_Jpm          - example radial Mathieu function Jpm at n, nmax, elliptical parameter q
%   example2_Spm          - example angular Mathieu function Spm at order n, nmax, q
%   extract_one_column    - function vec=extract_one_column(KF,t,mvncoeffs)
%   extract_one_value     - function value=extract_one_value(KF,t,vncoeffs)
%   gpm                   - function joiningFactors = gpm(KF,parameter,mv,nmax)
%   Hpm1                  - function values = Hpm1(category,u,q,mv,nmax)
%   Hpm2                  - function y = Hpm2(KF,u,q,mv,nmax)
%   Jpm                   - function values = Jpm(category,u,q,mv,nmax)
%   Npm                   - function normalizationFactor= Npm(category,coefficients)
%   Spm                   - function y = Spm(KF,v,mv,nmax)
%   Ypm                   - function y = Ypm(KF,u,q,mv,nmax)

%   getMatrix             - function matrix=getMatrix(category,parameter,numberCoefficients)
%   getVt                 - function coefficientIndices = getVt(catetory, numberCoefficients)
%   getNumberCoefficients - function ncoeffs = getNumberCoefficients()
%   getIdList              - function ik = getIdList(category,ncoeffs)
%   getTimeIndex           - function index =getTimeIndex(KF,time)
%   safeNumberCoefficients - function boolean = safeNumberCoefficients( n, nmax, ncoeffs )
%   symmetricEig           - function [V,E] = symmetricEig(J,KF)
%   b2nTest                 - Gutierrez-Vega, Table 4.1. 
%   bEvenModes              - q=25, b_{2n+2}(q),  n = 1:8 from Gutierrez-Vega and Leeb (1977).
%   ce2nTest                - test Fourier coefficients of angular Mathieu functions
%   coefficientCe10         - q = .1,  ce_10(eta,q) = sum A_2n cos(2n eta)


%   exampleMathieu          - elliptic cylinder functions
%   first3Modes             - q=1, k=1,2,3,4, modes [a0,a2,a4,a1,a3,a5,b2,b4,b6,b1,b3,b5]'
%   first4Coefficients      - q=10, k=1,2,3,4, coefficients ce0(0:2:6),ce1(1:2:7),se2(2:2:8),se1(1:2:7)
%   firstCoefficientsTest   - test Fourier coefficients of angular Mathieu functions
%   firstModeTest           - 


%   getCartesian            - function [x,y] = getCartesian(xi,eta,focalLength) Elliptic -> Cartesian
%   getElliptic             - function [xi,eta]=getElliptic(x,y,focalLength) Cartesian>Elliptic



%   parameterizedOde        - 
%   solveIvp2               - 
%   evaluateMathieu        - function [axial,radial] = evaluateMathieu(math,angles,radii)
%   exSurf                 - elliptic cylinder function
%   radialFunctionData     - function [u,q,value] = radialFunctionData(category)
%   testSeries             - test series
%   testRadial             - test radial Mathieu functions
%   testCoordinates        - test change of coordinates
%   series                 - function value=series(angles,coefficient,indices,trigonometry)

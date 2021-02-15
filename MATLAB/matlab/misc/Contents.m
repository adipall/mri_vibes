% MISC
%
% Files
%   createFunctionInput       - utility for SierraSD text input files
%   dofmap                    - function map=dofmap(basename,proc)
%   parallelmap               - function [localMap,globalRow]=parallelmap(name,numProc)
%   isVolumetic               - function check = isVolumetic(fmap)
%   exParallelMap             - example 
%   expandRmodel              - function [dispgr,nodes]=expandRmodel( cbmap, OTM, OutMap, vr )
%   getdispg                  - function disp_gset=getdispg(FetiMap_a,dispAset)
%   getdispr                  - function dispr=getdispr(cbmap,disp_gset)
%   getdispr_g                - function disp_gset=getdispr_g(cbmap,OTM,OutMap,fetimapdim,vr)
%   mk_netcdf_se              - generate a netcdf formatted text file of a superelement.
%   mkphi6                    - function phi=mkphi6(nsteps,nvar01,nvar02,nvar03,nvar04,nvar05,nvar06)
%   mkphi3                    - assemble the phi matrix from nvar01 through nvar03
%   getDofPerNode             - function eight = getDofPerNode(); DOF_PER_NODE
%   triplet                   - function map=triplet(gid,fmap) [equation nodeIndex local_dof]
%   exTriplet                 - example triplet
%   loadGids                  - function C = loadGids(root,p)
%   overlaps                  - example loadGids
%   loadMaps                  - function aCell = loadMaps(p)
%   loadDisp                  - function aCell = loadDisp(p)
%   loadConstraintMatrix      - function aCell = loadConstraintMatrix(p)
%   loadMass                  - function aCell = loadMass(p)
%   loadMatrix                - function aCell = loadMatrix(p)
%   getSharedNodes            - function [domains,numberSharedNodes,sharedGids,currentLids,adjacentLids]= getSharedNodes(gids)
%   getLastGlobalRow          - function last = getLastGlobalRow(globalRow)
%   casaMap                   - 1. Generate analysis and structural set maps
%   casaDisp                  - 2. convert mode shapes from structural set to analysis set
%   casaConsistencyTest       - 3. check consistency of matrices and vectors
%   casaModel                 - 4. Compute force = K*rbm, which should vanish
%   casaWorstConstrainedNodes - 5. Find the worst constraints and corresponding node numbers
%   casaConstraint            - 5.A. load constraint matrices
%   casaMassToo               - optional, e.g. check R'MR=eye
%   casaStudy                 - 6.






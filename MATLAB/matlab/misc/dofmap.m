%DOFMAP function map=dofmap(basename,proc)
% For parallel use parallelmap.
function map=dofmap(basename,varargin)
% function map=dofmap(basename,proc)
% returns a map of the dof versus GID and CID
% first column is DOF(1:n), next is GID(1:N) and last is CID(1:8)
%
% requires that 'basename_gid.m' and FetiMap_a.m exist. Both are
% written by Salinas.
procstr='';
nargin
if nargin==2
   np=varargin{1};
   procstr=['_' int2str(np)];
   disp(['the process string is: ' procstr]);
end

% get gids
fgidname=[basename '_gid' procstr '.m'];
if exist(fgidname,'file')~=2
   error(['file "' fgidname '" must exist']);
end
fgid=[basename '_gid' procstr ';'];
gid = [];
eval(fgid);
if  exist('gid','var')~=1
   error(['unable to execute "' fgid '"']);
end

% get feti map
mapname=['FetiMap_a' procstr '.m'];
if exist(mapname,'file')~=2
    error(['file "' mapname '" must exist.']);
end
fmap = [];
mapcmd=['fmap=FetiMap_a' procstr ';'];
eval(mapcmd);
if exist('fmap','var')~=1
    error(['Unable get feti map using: "' mapcmd '"']);
end


% check dimensions
nnodes=max(size(gid));
if ( nnodes*8 ~= size(fmap,1) )
    error(['mismatched map dimensions: ' int2str(nnodes) '*8 should equal ' int2str(size(fmap,1))]);
end

ndofs=max(fmap);
map=zeros(ndofs,3);
cdof=1;
for node=1:nnodes;
   for d=1:8
       k=(node-1)*8+d;
       if fmap(k)>0
           map(cdof,:)=[fmap(k) gid(node) d];
           cdof = cdof+1;
        end
    end
end


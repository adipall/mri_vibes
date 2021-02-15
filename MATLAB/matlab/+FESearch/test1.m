nodexyz=[0 0 1 % block 1:    nodexyz[1:27],   3x3x3 mesh on [0,1]^3
   .5 0 1
   1 0 1
   0 0 .5
   .5 0 .5
   1 0 .5
   0 0 0
   .5 0 0
   1 0 0
   0 .5 1
   .5 .5 1
   1 .5 1
   0 .5 .5
   .5 .5 .5    % This, the only interior node in block 1 will be skipped.
   1 .5 .5
   0 .5 0
   .5 .5 0
   1 .5 0
   0 1 1
   .5 1 1
   1 1 1
   0 1 .5
   .5 1 .5
   1 1 .5
   0 1 0
   .5 1 0
   1 1 0
   0 1 1   %%% block 2,  nodexyz(28:54,:)   3x3x3 mesh on [0,1] x [1,2] x [0,1]
   .5 1 1
   1 1 1
   0 1 .5
   .5 1 .5
   1 1 .5
   0 1 0
   .5 1 0
   1 1 0
   0 1.5 1
   .5 1.5 1
   1 1.5 1
   0 1.5 .5
   .5 1.5 .5  % This the only interior node in block 2 will also be dropped
   1 1.5 .5
   0 1.5 0
   .5 1.5 0
   1 1.5 0
   0 2 1
   .5 2 1
   1 2 1
   0 2 .5
   .5 2 .5
   1 2 .5
   0 2 0
   .5 2 0
   1 2 0];

idx=[1:13 15:40 42:54]; % omit interior nodes
disp=zeros(size(nodexyz));
disp(1:27,2)=5.5e-4*ones(27,1);

%% the nodes of the masters surface are based upon the idx node numbers
%  not the original 54 nodes
ms=[1 2 11 10    % ms(1:12,:) is a surface mesh of block 1
   2 3 12 11
   10 11 19 18
   11 12 20 19
   18 19 22 21
   19 20 23 22
   22 23 26 25
   21 22 25 24
   3 6 14 12
   6 9 17 14
   14 17 26 23
   12 14 23 20
   27 28 37 36  % ms(13:24,:) is a surface mesh of block 2
   28 29 38 37
   36 37 45 44
   37 38 46 45
   29 32 40 38
   32 35 43 40
   38 40 49 46
   40 43 52 49
   44 45 48 47
   45 46 49 48
   48 49 52 51
   47 48 51 50];

obj=FESearch.Search;
inp=obj.ParticleSearch(nodexyz(idx,:),ms,nodexyz(idx,:),max(abs(disp(idx,:)))); % disp could be dt*v
% For each surfaceElement, inp{surfaceElement,:} is the nodes closest to ms(surfaceElement,:)

numElement = size(ms,1);
for element=1:numElement,
    neighbor = setdiff(  inp{element,:},  ms(element,:));
    assert( size( inp{element,:},1 ) - size( neighbor,1 ) == 4 );
end

%CASASTUDY 6.
nodes = [40656 40657 52257;...
         41330 41487 52417;...
         52665 52666 52450];
load('casamesh.mat');% blk**, x0,y0,z0
% blocks 11,12,16,17
% [a,e] = find(blk11==nodes(1,1)); % (2,515)
% [a,e] = find(blk11==nodes(1,2)); % (6,515)
% [a,e] = find(blk16==nodes(1,3)),
% [a,e] = find(blk12==nodes(2,1)); % (5,393)
% [a,e] = find(blk12==nodes(2,2)); % (1,393)
% [a,e] = find(blk16==nodes(2,3)),
% [a,e] = find(blk17==nodes(3,1)); % (1,294)
% [a,e] = find(blk17==nodes(3,2)); % (2,294)
% [a,e] = find(blk16==nodes(3,3)),
% x0(nodes) = 17.7
for i = 1:3,
    liney = y0(nodes(1:2,i));
    linez = z0(nodes(1:2,i));
    pointy = y0(nodes(3,i));
    pointz = z0(nodes(3,i));
    subplot(1,3,i);
    plot(liney,linez,'k',pointy,pointz,'+');
end





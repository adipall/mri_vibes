%function plotPoly( element, connectivity, coords, color )
function plotPoly( element, connectivity, coords, color )
nGon = size( connectivity, 2); 
for pt = 1:nGon,
    i = connectivity(element,pt);
    if pt < nGon,
       next = pt + 1;
    else 
       next = 1;
    end
    f = connectivity(element,next);
    plot(  [coords(i,1), coords(f,1)], [coords(i,2), coords(f,2)], color  );  
    hold on;
end

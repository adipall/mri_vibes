%function plotQuad( ids, coords, flag)
function plotQuad( ids, coords, flag )
ids(5) = ids(1);
x = coords(ids,1);
y = coords(ids,2);
z = coords(ids,3);
if flag == 0,
    plot3(x,y,z,'b');
end
if flag == 1,
    plot3(x,y,z,'k');
end
if flag == 2,
    plot3(x,y,z,'--b');
end
if flag == 3,
    plot3(x,y,z,'--k');
end
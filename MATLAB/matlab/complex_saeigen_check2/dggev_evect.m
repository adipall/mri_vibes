function u=dggev_evect(ai,v)
% function u=dggev_evect(ai,v)
% convert vector v in dggev format to standard form

n=max(size(ai));

i=sqrt(-1);
j=1;
while ( j<=n )
    if ( ai(j) == 0 )
        u(:,j)=v(:,j);
        j=j+1;
    else
        u(:,j)=v(:,j)+i*v(:,j+1);
        u(:,j+1)=v(:,j)-i*v(:,j+1);
        j=j+2;
    end
end


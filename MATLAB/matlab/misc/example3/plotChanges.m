figure(1);
subdAd = getSubdAd();
for subdomain = 1:4,   % subdomain = rank
    [forceVector, nodalSum, displacementVector] = getSubdomainFeVectors(subdomain);

    % nodal sum step.   gather to owner
    delta =  nodalSum - forceVector;
    [matchListSize,matchList,owner]= getSubdomainAdjacency(subdomain);
    numberNeighbors = size(matchListSize,1) - 1;
    filteredDelta = delta;
    for i = 1:numberNeighbors,
        neighbor = subdAd{subdomain}(i);  % neighbor sends to subdomain
        if( neighbor < subdomain ),
            start = matchListSize(i) + 1;
            stop = matchListSize(i+1);
            numberShared = stop - start + 1;
            sharedIds = matchList( start:stop ) + int32( ones(numberShared,1) );
            filteredDelta( sharedIds ) =  zeros(numberShared,1);
        end
    end
    mustBeZero(4) = norm( filteredDelta );
    % scatter from owner to shared dof
    change = displacementVector - nodalSum;
    filteredChange = change;  
    for i = 1:numberNeighbors,
        neighbor = subdAd{subdomain}(i);
        if( neighbor > subdomain ),  % neighbor > receives from subdomain
            start = matchListSize(i) + 1;
            stop = matchListSize(i+1);
            numberShared = stop - start + 1;
            sharedIds = matchList( start: stop ) + int32( ones(numberShared,1) );
            filteredChange( sharedIds ) = zeros(numberShared,1);
        end
    end
    mustBeZero(5) = norm( filteredChange );
    subplot(4,1,subdomain);
    plot(delta,'gs'); hold on;
    plot(change,'k+'); hold off;
    assert( sum( mustBeZero ) == 0 );
end
disp('Passed');  % all assertions hold




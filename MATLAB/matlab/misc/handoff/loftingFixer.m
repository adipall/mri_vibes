% hypothesis: tripetTest passes
clear
load('simple_loft.mat');
filename='simple_loft.txt';
threshold = 1.e-12;
[mpcVector,mpcOrder]=readMpc(filename);
numMpc = size(mpcOrder,1)/3;
newMpcOrder = zeros(numMpc,1);
residual = zeros(numMpc,1);
resid2 = zeros(numMpc,1);
updateSize = zeros(numMpc,1);
start=0;
newstart=0;
coords = [x0,y0,z0] + [nvar01(:,nsteps),nvar02(:,nsteps),nvar03(:,nsteps)];
for k=1:numMpc,
    n = mpcOrder(3*k-2);
    newMpcOrder(k)=n;
    serialNode = zeros(n,1);
    for j = 1:n,
        serialNode(j) = find( node_num_map == mpcVector(start+j,1));
    end
    weight = mpcVector(start+1:start+n,3);
    start = start + 3*n;
    e = ones(n,1);
    R = [e, coords(serialNode,:)];
    residual(k)= norm(R'*weight,1);%what the code does but a bit
    centroid = e'*R(:,2:4)/n;% of an unfair comparison
    R(:,2:4)=R(:,2:4) - e*centroid; % I really care abut the svd.
    n2 = n;
    if residual(k) > threshold,
        [active,~]=find(abs(weight)>threshold);
        numActive = size(active,1);
        [U,S,V] = svd(R); % todo use QR instead of SVD
        singularValues = diag(S);
        if singularValues(4) < .1,
            maxWeight = norm(weight,inf);
            slaveId = find(abs(weight)==maxWeight);
            if size(slaveId,1) == 2,
                masterIds = slaveId(2);
                slaveId = slaveId(1);
            else
                e = ones(mpcOrder(3*k-2),1);
                masterIds=find((abs(weight)>e*threshold) & (abs(weight)<e));
            end
            assert( size(slaveId,1) == 1 );
            slaveNode = serialNode( slaveId );
            masterNodes = serialNode(masterIds);
            matchingBlock = 0;
            act = 'elementblock=';
            root = 'blk0';
            assert( 0<nblks && nblks < 100 );
            block = 1;
            match = false;
            while ((~match) && (block<=nblks)),
                if  9 < block
                    root = 'blk';
                end
               eval( [act,root,int2str(block),';']);
               if isempty(find(elementblock == masterNodes(1),1)),
                   block=block+1;
               else
                   match=true;
               end
            end
            assert( block<=nblks);
                            % fixme
                % assignment of one element to a constraint
                % an element contains a node.  add to the weight of the
                % element the absolut value of the coefficient.
                % Choose one of the heaviest elements.  
            [row1,col1] = find( elementblock == masterNodes(1) );
            if numActive > 2,
                [row2,col2] = find( elementblock == masterNodes(2) );
                element = intersect(col1,col2);%slow
            else
                element = col1;
            end
            if size(element,1) > 1 && numActive <= 3,
               element = element(1);  
            end
            if numActive > 4,
                firstpair = element;
                [row3,col3] = find( elementblock == masterNodes(3) );
                [row4,col4] = find( elementblock == masterNodes(3) );
                secondpair = intersect( col3, col4 );%slow
                element = intersect( firstpair, secondpair );%slow
            end
            assert( size(element,1) == 1 );% fixme numActive > 5
            if numActive > 4,
                clear firstpair row3 col3 row4 col4 seconpair
            end
            n2=9;
            newMpcOrder(k)=n2;
            newWeight = zeros(n2,1);
            newWeight(1) = weight(slaveId);
            slaveId = 1;
            serialNode = [slaveNode;elementblock(:,element)];
            retainedMasterNodes = intersect(serialNode,masterNodes);%slow
            carryOver = size(retainedMasterNodes,1);
            newId = zeros(carryOver,1);
            for i=1:carryOver,
                newId(i) = find(serialNode == masterNodes(i));
            end
            for i=1:carryOver,
               newWeight(newId(i)) = weight(masterIds(i));
            end
            if carryOver < numActive-1,
                disp('Warning: dropping some active nodes from constraint');
                gap = norm(weight,1)-norm( newWeight,1 );
                disp(gap);
            end
            weight = newWeight;
            e = ones(n2,1);
            R = [e, coords(serialNode,:)];
            centroid = e'*R(:,2:4)/n;
            R(:,2:4)=R(:,2:4) - e*centroid;
            [U,S,V] = svd(R);
            singularValues = diag(S);  
           
        
             
        end
        assert( singularValues(4)>= .1);% fixme: replace w/?
        U = U(:,1:4);
        w = weight - U*(U'*weight);
        
        updateSize(k) = norm(U'*weight);
        resid2(k) = norm(R'*w,1);
        newMpcVector(newstart+1:newstart+n2,2)=w;
    else
        newMpcVector(newstart+1:newstart+n2,2)=weight;
    end % residual is large
    newMpcVector(newstart+1:newstart+n2,1)=node_num_map(serialNode);
    newstart = newstart + n2;
end % k = 1: nummpc
d = linspace(1,numMpc,numMpc)';
figure(1);
subplot(2,1,1);
[~,tightOrder]=tightMpc(mpcVector,mpcOrder);
handle=plot(d,mpcOrder(1:3:end),'ko',d,tightOrder,'r+',d,newMpcOrder,'bs');
legend(handle,'order','filtered','new','location','Best');
subplot(2,1,2);
handle = semilogy(d,residual,'ko',d,resid2,'r+',d,updateSize,'bx');
legend(handle,'resid','resid2','update','location','Best');
assert( norm( resid2,inf) < threshold ),

if strcmp(filename,'simple_loft.txt'),
    simpleloftFigure(x0,y0,z0,mpcOrder,mpcVector,newMpcOrder,newMpcVector)
end

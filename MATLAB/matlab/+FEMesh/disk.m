classdef disk
    properties (SetAccess = private)
    end
    methods
        function obj = disk()
        end % ctor
    end  % methods
    methods (Static = true )
        function [connect,coords] = constructor(object,n)
            [connect,coords] = object.setQuarterMesh(object,n);
            [connect2,coords2,top2] = object.semicircleMesh(object,n,connect,coords);
            connect = connect2;
            coords = coords2;
            top = top2;
            clear connect2 coords2 top2;
            [connect2,coords2] = object.diskMesh(coords, connect,top, n);
            connect = connect2;
            coords = coords2;
            clear connect2 coords2 top;

        end
        % function  = plotQuadMesh(conn,coords,fig)
        function plotQuadMesh(conn,coords,fig)
            assert( size(conn,2) == 4 );
            figure(fig);
            m = size(conn,1);
            for i=1:m,
               x = coords(conn(i,:),1);
               y = coords(conn(i,:),2);
               x(5) = x(1);
               y(5) = y(1);
               plot(x,y); hold on;
            end
            hold off;
        end
        %function [connect,z4] = setQuarterMesh(n)
        function [connect,z4] = setQuarterMesh(object,n)
            w = object.setSquareDomain(n);
            z = exp(-.5*w.*w);
            rmin = z(n*n);
            r = abs(z);
            dir = z./r;
            z2 = (z - dir*rmin)/(1-rmin);
            z3 = zeros(n*n,1);
            top = (1:n)*n;
            angles = imag(log(z2));
            radius = abs(z2);
            for i = 0:n-1,
                for j = 0:n-1,
                   k = (1+j+i*n);
                   if angles(k) >= 0,
                       base=angles(top(i+1));
                   else
                       base=-angles(n*(n-1)+j+1);
                   end
                   if base > 0
                       newArg = pi*angles(k)/(4*base);
                   else
                       newArg = angles(k);
                   end
                   r = radius(k);
                   z3(k) = r*exp(1i*newArg);
                end
            end
            [z4,connect] = object.finishQuarterMesh(object,n,z3);
        end
        %function [connect2,coords2,top2] = semicircleMesh(n,connect,coords)
        function [connect2,coords2,top2] = semicircleMesh(object,n,connect,coords)
            % The numbering of the flat edges is from top to bottom,
            numNode = size(coords,1);% top(1:n),  bottom(1:n),   top(n) = bottom(1).
            [top,bottom] = object.getFlatEdges(object,n,numNode);
            x = real(coords);% Make the new nodes. (x,y) -> (y,x)
            y = imag(coords);
            newOrTop = complex( y, x );% includes the top
            [newNodeId,selectNew] = setdiff( 1:numNode, top);
            newCoords = newOrTop(newNodeId);
            assert( size(top,2) == n);
            numNewNode = numNode - size(top,2);% number the new nodes
            oldToNew = zeros(1,numNode);
            oldToNew( selectNew ) = (1:numNewNode) + ones(1,numNewNode)*numNode;
            % The next stop is to reflect across the new top nodes,
            % the union of the bottom nodes and the new bottom nodes.
            justBottom = bottom(2:n);
            newBottom = oldToNew( justBottom );
            top2 = [newBottom,bottom];
            % Make new elements
            tol = 1.e-9;
            numElement = size( connect,1);
            nskip = 0;
            adj = zeros(numElement,4);
            for element = 1:numElement,
                ids = connect(element,:);
                z = coords(ids);
                x = real(z);
                y = imag(z);
                distanceToMirror = min(abs(y-x));
                reverse = ids(4:-1:1); % reflections reverse orientation
                if distanceToMirror > tol,
                    adj(element,:) = oldToNew( reverse );
                else
                    g = abs(y-x);
                    g = g(4:-1:1);
                    mixed = reverse;
                    for i=1:4,
                        if g(i) > tol,
                           mixed(i) = oldToNew( mixed(i) );
                        end
                    end
                    adj(element,:) = mixed;
                    nskip = nskip + 1;
                end
            end
            assert( nskip == n-1 );
            coords2 = [coords;newCoords];
            connect2 = [connect;adj];
        end
        %function [connect2,coords2] = diskMesh(coords, connect,top, n)
        function [connect2,coords2] = diskMesh(coords, connect,top, n)
            x = real(coords);% (x,y) -> (-y,-x)
            y = imag(coords);
            newOrTop = complex( -y, -x );% includes the top
            numNode = size(coords,1);
            [newNodeId,selectNew] = setdiff( 1:numNode, top);
            newCoords = newOrTop(newNodeId);
            assert( size(top,2) == 2*n-1);
            numNewNode = numNode - size(top,2);  % number the new nodes
            oldToNew = zeros(1,numNode);
            oldToNew( selectNew ) = (1:numNewNode)+ones(1,numNewNode)*numNode;
            tol = 1.e-9;% Make new elements
            numElement = size( connect,1);
            nskip = 0;
            adj = zeros(numElement,4);
            for element = 1:numElement,
                ids = connect(element,:);
                z = coords(ids);
                x = real(z);
                y = imag(z);
                  g = abs(y+x);
                distanceToMirror = min(g);
                reverse = ids(4:-1:1); % reflections reverse orientation
                if distanceToMirror > tol,
                    adj(element,:) = oldToNew( reverse );
                else
                    g = g(4:-1:1);
                    mixed = reverse;
                    for i=1:4,
                        if g(i) > tol,
                           mixed(i) = oldToNew( mixed(i) );
                        end
                    end
                    adj(element,:) = mixed;
                    nskip = nskip + 1;
                end
            end
            assert( nskip == 2*n-2 );
            coords2 = [coords;newCoords];
            connect2 = [connect;adj];
        end
        %
        %
        function w = setSquareDomain(n)
            parameter = .5;
            K = ellipke(parameter);
            h = 1/(n-1);
            wo = (1+1i)*K*.5;
            w = zeros(n*n,1);
            for i = 0:n-1,
                for j = 0:n-1,
                   w(1+j+i*n)=wo*i*h + wo'*j*h;
                end
            end
        end
        % Step 5.1: elide nodes (0,1) and (1,0),
        % Step 5.2: generate element to node adjacency
        % with respect to the new node numbering.
        function [coords,connect] = finishQuarterMesh(object,n,w)
            k = 0;
            numNode = n*n -2;
            coords = zeros(numNode,1);
            for i = 0:n-1,
                for j = 0:n-1,
                   if i+j ~= 1,
                       k = k + 1;
                       coords(k) = w(1+j+i*n);
                       assert(k == object.newNodeIdAfterStep5(i,j,n));
                   end
                end
            end
            numElement = (n-1)*(n-1) -1;
            connect= zeros( numElement ,4);
            element = 0;
            for i = 0:n-2,
                for j = 0:n-2,
                    if i + j > 0, % skip (0,0)
                        element = element + 1;
                        if i+j > 1,
                            a = object.newNodeIdAfterStep5(i,j,n);
                        else % (1,0),(0,1)
                            a = object.newNodeIdAfterStep5(0,0,n);
                        end
                        b = object.newNodeIdAfterStep5(i,j+1,n);
                        c = object.newNodeIdAfterStep5(i+1,j+1,n);
                        d = object.newNodeIdAfterStep5(i+1,j,n);
                        connect(element,:) = [a,b,c,d];
                    end
                end
            end
        end
        %function k = newNodeIdAfterStep5(i,j,n)
        function k = newNodeIdAfterStep5(i,j,n)
            if i == 0,
               if j == 1,
                   k = 0; %error
               else
                   if j == 0,
                       k = 1;
                   else
                       k = j;
                   end
               end
            else
               if i==1 && j == 0
                   k=0;  %error
               else
                   k= j + i*n-1; % 1+j+i*n-2
               end
            end
        end
        function [top,bottom] = getFlatEdges(object,n,numNode)
            % The numbering of the flat edges is from top to bottom,
            % top(1:n),  bottom(1:n),   top(n) = bottom(1).
            bottom = numNode-n+1 : numNode;
            bottom = bottom(n:-1:1);
            top = zeros(1,n);
            k = 0;
            for i=0:n-1,
                j = n-1;
                index = object.newNodeIdAfterStep5(i,j,n);
                if index > 0
                    k = k + 1;
                    top(k) = index;
                end
            end
        end
    end % methods
end % class

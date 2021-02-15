classdef icosahedron
    properties (SetAccess = private)
        coordinates = [];
        connectivity = [];
    end
    methods
        function obj = icosahedron()
            obj.coordinates = obj.getCoordinates();
            n = size(obj.coordinates,1);
            edges = obj.getEdges(obj.coordinates); % beam mesh
            obj.connectivity = obj.getFaces(edges);
            twenty = size(obj.connectivity,1);
            for face = 1:twenty,
                i = obj.connectivity(face,1);
                j = obj.connectivity(face,2);
                e1 = obj.coordinates(j,:) - obj.coordinates(i,:);
                k = obj.connectivity(face,3);
                e2 = obj.coordinates(k,:) - obj.coordinates(j,:);
                normal = cross(e1,e2)';
                orientation = obj.coordinates(i,:)*normal;
                if orientation < 0
                    obj.connectivity(face,2) = k;
                    obj.connectivity(face,3) = j;
                end
            end %face
        end % ctor
    end  % methods 
    methods (Static = true )
        function coords = getCoordinates()
            p = .5*(1+sqrt(5));
            sigh = 1;
            sp = 1;
            coords = zeros(12,3);
            k = 0;
            for i = 1:2,
                for j = 1:2
                    point = [0,sigh,sp*p];
                    coords(k+1,:) = point;
                    forward = [3 1 2];
                    coords(k+2,:) = point(forward);
                    ff = [2 3 1];
                    coords(k+3,:) = point(ff);
                    sigh = -sigh;
                    k=k+3;
                end
                sp = -sp;
            end
        end
        function edges = getEdges(coordinates)
            edges = zeros(12,5);
            for i = 1:12, 
                  d = coordinates - ones(12,1)* coordinates(i,:);
                  r2 = sum( d.*d , 2);
                  neighbors = find(r2 == 4);
                  assert( size(neighbors,1) == 5 );
                  edges(i,:) = neighbors';
            end
        end
        function faces = getFaces(edges)
            ten = size(edges,1) - 2;
            d = size(edges,2);
            n = 0;
            faces = zeros(20,3);
            for vertex = 1:ten,
                for i = 1:d-1,
                    e = edges(vertex,i);
                    if e > vertex,
                        adj = edges(e,:);
                        for j = i+1:d,
                            neighbor = edges(vertex,j);
                            location = find( neighbor == adj );
                            assert( vertex < neighbor ) % skip
                            if( size(location,2) == 1 )
                                % adj(location) == edges(1,j)
                                triplet = [vertex,e,neighbor];
                                n = n + 1;
                                faces(n,:) = triplet;
                            end
                        end
                    end
                end
            end % for vertex
        end  %  getFaces(edges)
    end % methods (Status=ture)
end % icosahedron class

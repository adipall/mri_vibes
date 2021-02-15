classdef cubedSphere
    properties (SetAccess = private)
    end
    methods
        function obj = cubedSphere()
        end % ctor
    end  % methods
    methods (Static = true )
        function [X,Y,Z] = cartesianCoordinates(sphericalCoordinates)
            latitude  = sphericalCoordinates(:,1);
            longitude = sphericalCoordinates(:,2);
            X = cos(latitude).*cos(longitude);
            Y = sin(latitude).*cos(longitude);
            Z = sin(longitude);
        end
%function connectivity = getConnectivity(n)
function connectivity = getConnectivity(n)
    n = round(n);
    if n < 2,
      n = 2;
    end
    if n == 2,
      connectivity = [1 2 3 4; 2 6 7 3; 6 5 8 7 ;5 1 4 8; 4 3 7 8 ; 5 6 2 1];
      return;
    end
    m = n-2;
    edgeOffset = 6*m^2; % number of internal nodes
    cornerOffset = edgeOffset + 12*m;
    numberFaces = 6*(n-1)^2;
    connectivity = zeros(numberFaces,4);
    cornerIds = [1 2 3 4; 2 6 7 3; 6 5 8 7 ;5 1 4 8; 4 3 7 8 ; 5 6 2 1];
    edgeIds = [1 2 3 4; 9 6 10 2; 5 8 7 6;12 4 11 8;3 10 7 11;5 9 1 12];
    % [adjacentElement,adjacentLocalEdge,orientation]=getAdjacency();
    k = 0;
    for cubeFace = 1:6,
        localEdge = 1;
        for locationInEdge=1:m+1,
            if locationInEdge==1,
                lowLeft = getCornerId(n,cubeFace,localEdge);
                lowRight = getEdgeId(n,cubeFace,localEdge,locationInEdge);
                highRight = m*m*(cubeFace-1) + m*(localEdge-1)+locationInEdge;
                previous = minusBase1(localEdge,4);
                highLeft = getEdgeId(n,cubeFace,previous,locationInEdge);
            else
                if locationInEdge <= m,
                    lowLeft   = connectivity(k,2);
                    highLeft  = connectivity(k,3);
                    lowRight = getEdgeId(n,cubeFace,localEdge,locationInEdge);
                    highRight = highLeft + 1;
                else
                    % locationInEdge=n-1;
                    previous = locationInEdge-1;%m
                    lowLeft = getEdgeId(n,cubeFace,localEdge,previous);
                    next = plusBase1(localEdge,4);
                    lowRight = getCornerId(n,cubeFace,next);
                    ei1 = edgeIds(cubeFace,2);
                    if ei1 ~= 9,
                        e1 = edgeOffset + m*(ei1-1)+1;
                    else
                        e1 = edgeOffset + m*ei1;
                    end
                    highRight = e1;
                    highLeft = m*m*(cubeFace-1) + m*(localEdge-1) + m;
                end
            end
            k = k + 1;
            connectivity(k,:) = [lowLeft,lowRight,highRight,highLeft];
        end
        for localEdge=2:m, %m=2,  m=3
            for i=1:m+1,   % 5:8, 9:12
                if i == 1,
                    locationInEdge=localEdge-1;
                    previous = 4;
                    lowLeft = getEdgeId(n,cubeFace,previous,locationInEdge);
                    lowRight = m*m*(cubeFace-1)+m*(localEdge-2)+1;
                    locationInEdge=localEdge;
                    highLeft = getEdgeId(n,cubeFace,previous,locationInEdge);
                    highRight = lowRight + m;
                else
                    lowLeft = lowRight;
                    highLeft = highRight;
                    if i <= m,
                        lowRight = lowLeft + 1;
                        highRight  = highLeft  + 1;
                    else
                        locationInEdge=localEdge-1;
                        next = 2;
                        lowRight = getEdgeId(n,cubeFace,next,locationInEdge);
                        locationInEdge=localEdge;
                        highRight = getEdgeId(n,cubeFace,next,locationInEdge);
                    end
                end
                k = k + 1;
                connectivity(k,:) = [lowLeft,lowRight,highRight,highLeft];
            end
        end
        localEdgeofFace = n-1;
            for i=1:n-2,
                if i == 1,
                    elementEdge = 4;
                    locationInEdge=localEdgeofFace-1;
                    lowLeft = getEdgeId(n,cubeFace,elementEdge,locationInEdge);
                    lowRight = m*m*(cubeFace-1) + m*(m-1)+ 1;
                    elementEdge = 3;
                    highRight = getEdgeId(n,cubeFace,elementEdge,i);
                    highLeft = getCornerId(n,cubeFace,previous);
                 else
                    lowLeft = lowRight;
                    highLeft = highRight;
                    lowRight = lowLeft + 1;
                    elementEdge = 3;
                    highRight = getEdgeId(n,cubeFace,elementEdge,i);
                 end
                k = k + 1;
                connectivity(k,:) = [lowLeft,lowRight,highRight,highLeft];
            end
            %i=n-1; % face, edge 2, corner 3, edge 3
            lowLeft = lowRight;
            highLeft = highRight;
            ei1 = edgeIds(cubeFace,2);
            if ei1 ~= 9,
                lowRight = edgeOffset + m*ei1;
            else
                 lowRight = edgeOffset + m*(ei1-1)+1;
            end
            highRight = cornerOffset + cornerIds(cubeFace,3);
            k = k + 1;
            connectivity(k,:) = [lowLeft,lowRight,highRight,highLeft];
    end
end
%function coordinates = getCoordinates(n)
function coordinates = getCoordinates(n)
%                _________ % North Pole - 5
%                |       |
%                |   5   |
%        ________|_______|_______________
%        |       |       |       |       |
%        |   4   |   1   |   2   |   3   |
%        |_______|_______|_______|_______|
%                |       |
%                |   6   |
%                |_______| %  South Pole - 6
    n = round(n);
    if n < 2,
      n = 2;
    end
    if n == 2,
      coordinates = getCornerAtlas();
      return;
    end
%4 3:local corner order
%1 2
%                8_______7:global corner order
%                |       |
%                |   5   |
%        8_______4_______3_______7_______8
%        |       |       |       |       |
%        |   4   |   1   |   2   |   3   |
%        5_______1_______2_______6_______5
%                |       |
%                |   6   |
%                5_______6   compare cornerIds
%
%
%  3:local edge order
%4    2
%  1
%                ____7____global edge order
% edge 11 flips  |       |
% like edge 9    11      10
%        ___11___|___3___|__10_______7___
%        |       |       |       |       |
%        8       4       2       6       8
%        |__12___|___1___|___9___|___5___|
%                |       |
%               12       9
% cf edgeIds     |___5___|   edge 9 orientation changes between elements
% re-write to use spherical coorinates.
    m = n-2;
    alpha = linspace(-pi/4, pi/4, n);
    t = tan(alpha)';  % -1 = t(1)<t(n) = 1
    clear alpha;
    gnomonicCoordinates = zeros(m*m,2);
    eine = ones(m,1);
    for j =2:n-1,
        ki = 1 + (j-2)*m;
        kf = (j-1)*m;
        gnomonicCoordinates(ki:kf,:) = [t(2:n-1),eine*t(j)];
    end
    clear t j
    coordinates = zeros(6*m*m,2);
    for cubeFace = 1:6,
        ki = 1 + (cubeFace-1)*m*m;
        kf = cubeFace*m*m;
        coordinates(ki:kf,:) =  localFrame(gnomonicCoordinates,cubeFace);
    end
    clear ki kf cubeFace gnomonicCoordinates
    numberCorners = 8;
    cornerAtlas = getCornerAtlas();
    numberEdges = 12;
    edgeAtlas =  getEdgeAtlas(n);
    coordinates = [coordinates;edgeAtlas;cornerAtlas];
    numberVertices = numberCorners + numberEdges*m + 6*m*m;
    assert( size(coordinates,1) == numberVertices );
end
% function [cornerAtlas,copyCornerAtlas]=getCornerAtlas()
function [cornerAtlas,copyCornerAtlas]=getCornerAtlas()
%4 3:local corner order
%1 2
%                8_______7:global corner order
%                |       |
%                |   5   |
%        8_______4_______3_______7_______8
%        |       |       |       |       |
%        |   4   |   1   |   2   |   3   |
%        5_______1_______2_______6_______5
%                |       |
%                |   6   |
%                5_______6   compare cornerIds
    localCorners = [-1,-1;1,-1;1,1;-1,1];
    cornerIds = [1 2 3 4; 2 6 7 3; 6 5 8 7 ;5 1 4 8; 4 3 7 8 ; 5 6 2 1];
    cornerAtlas = zeros(8,2);
    pair = [1,3];
    for jj = 1:2,
       cubeFace = pair(jj);
       cornerCoordinates = localFrame( localCorners, cubeFace );
       cc = cornerIds(cubeFace,:);
       for ii = 1:4,
           gg = cc(ii);
           cornerAtlas(gg,:) = cornerCoordinates(ii,:);
       end
    end
    copyCornerAtlas = zeros(8,2);
    % pair = [2,4];
    pair =[5,6];
    for jj = 1:2,
       cubeFace = pair(jj);
       cornerCoordinates = localFrame( localCorners, cubeFace );
       cc = cornerIds(cubeFace,:);
       for ii = 1:4,
           gg = cc(ii);
           copyCornerAtlas(gg,:) = cornerCoordinates(ii,:);
       end
    end
end
% function id = getCornerId(n,cubeFace,localEdge)
function id = getCornerId(n,cubeFace,localEdge)
    cornerIds = [1 2 3 4; 2 6 7 3; 6 5 8 7 ;5 1 4 8; 4 3 7 8 ; 5 6 2 1];
    m = n-2;
    edgeOffset = 6*m^2; % number of internal nodes
    cornerOffset = edgeOffset + 12*m; % internal + edges
    id = cornerOffset + cornerIds(cubeFace,localEdge);
end
%function edgeAtlas = getEdgeAtlas(n)
function edgeAtlas =  getEdgeAtlas(n)
    n = round(n);
    if n < 2,
      n = 2;
    end
    if n == 2,
       edgeAtlas= [];
       return;
    end
%  3:local  order
%4    2
%  1
%                ____7____global  order
% edge 11 flips  |       |
% like edge 9   11      10
%        ___11___|___3___|__10_______7___
%        |       |       |       |       |
%        8       4       2       6       8
%        |__12___|___1___|___9___|___5___|
%                |       |
%               12       9
% cf edgeIds     |___5___|    9 orientation changes between elements
    t = tan( linspace(-pi/4, pi/4, n))'; % -1 = t(1)<t(n) = 1
    m = n-2;
    eine = ones(m,1);
    t = t(2:n-1);
    coordinates = [t,-eine;eine,t;t,eine;-eine,t];
    maskedEdgeIds = [1 2 3 4; 9 6 10 0; 5 8 7 0;12 0 11 0;0  0 0  0;0 0 0 0];
    %edgeIds = [1 2 3 4; 9 6 10 2; 5 8 7 6;12 4 11 8;3 10 7 11;5 9 1 12];
    edgeAtlas = zeros(2*m,2);
    for cubeFace = 1:6,
       edgeCoordinates = localFrame( coordinates, cubeFace );
       ee = maskedEdgeIds(cubeFace,:);
       for localId = 1:4,
           i = 1 + (localId-1)*m;
           f = localId*m;
           gg = ee(localId);
           if gg > 0,
               first = 1 + (gg-1)*m;
               last = gg*m;
               edgeAtlas(first:last,:) = edgeCoordinates(i:f,:);
           end
       end
    end
end
% function id = getEdgeId(n,cubeFace,localEdge,locationInEdge)
function id = getEdgeId(n,cubeFace,localEdge,locationInEdge)
    m = n-2;
    edgeOffset = 6*m^2; % number of internal nodes
    %
    edgeIds = [1 2 3 4; 9 6 10 2; 5 8 7 6;12 4 11 8;3 10 7 11;5 9 1 12];
    ei0 = edgeIds(cubeFace,localEdge);
    %
    orientation = ones(6,4);
    orientation(5,3:4) = [-1,-1];
    orientation(6,1:2) = [-1,-1]; % of my edge relative to the reference edge
    edgeOrientation = orientation(cubeFace,localEdge);
    if edgeOrientation == 1,
        id = edgeOffset + m*(ei0-1)+locationInEdge;
    else
        id = edgeOffset + m*ei0-locationInEdge+1;
    end
end
% function data = elementDiagnostic(x,type)
function data = elementDiagnostic(x,type)
    assert( size(x,1) == 4);
    assert( size(x,2) == 3);
    if type == 1, %
        forward = [2 3 4 1];
        y = x(forward,:);
        d = y-x;
        edgeLengths = sqrt( sum( d.*d,2));
        data = edgeLengths;
    else
        if type == 2,
            forward = [2 3 4 1];
            y = x(forward,:);
            d = y-x;
            n = zeros(2,3);
            n(1,:) = cross(d(1,:),d(2,:));
            n(2,:) = cross(d(3,:),d(4,:));
            norms = sqrt( sum( n.*n,2));
            n = diag(norms)\n;
    %         normalized =  sum( n.*n,2);
            cosine = sum( n(1,:) .* n(2,:) );
            angle = acos( cosine );
            data = angle;
        else
           distance = zeros(4,4);
           for i = 1:4,
             for j = 1:4,
                distance(i,j) = norm( x(i,:) - x(j,:));
             end
           end
           data = distance;
        end
    end
end
%function sphericalCoordinates =  localFrame(localCoordinates,cubeFace)
function sphericalCoordinates =  localFrame(localCoordinates,cubeFace)
    n = size(localCoordinates,1);
    if n < 1,
        sphericalCoordinates = [];
        return
    end
    assert( size(localCoordinates,2)== 2);
    R = sqrt(ones(n,1)+localCoordinates(:,1).^2 + localCoordinates(:,2).^2 );
    % locals: (:,1),  (n,:), (:,n),  (1,:)
      switch cubeFace,
         case 1,
           Latitude = atan2( localCoordinates(:,1) ,1.0);
           Longitude = asin( localCoordinates(:,2)./R );
         case 2
           Latitude = atan2(1.0,-localCoordinates(:,1));
           Longitude = asin( localCoordinates(:,2)./R );
         case 3
           %Latitude = atan2(-localCoordinates(:,1),-1.0);
           Latitude = atan2( localCoordinates(:,1) ,1.0) + ones(n,1)*pi;
           Longitude = asin( localCoordinates(:,2) ./R );
         case 4
           Latitude = atan2(-1.0, localCoordinates(:,1));
           Longitude = asin( localCoordinates(:,2) ./R );
         case 5
            Latitude = atan2( localCoordinates(:,1), -localCoordinates(:,2)) ;
            Longitude = asin( 1./R );
         case 6
            Latitude = atan2( localCoordinates(:,1) , localCoordinates(:,2));
            Longitude = asin( -1./R );
      end
    sphericalCoordinates = [ Latitude, Longitude ];
end
%function previous = minusBase1(value, modulus )
function previous = minusBase1(value, modulus )
    if value == 1,
       previous = modulus;
    else
      previous = value - 1;
    end
end
%function next = plusBase1(value, modulus )
function next = plusBase1(value, modulus )
    if value == modulus,
       next = 1;
    else
      next = value + 1;
    end
end
    end % methods (Status=ture)
end % cubedSphere class

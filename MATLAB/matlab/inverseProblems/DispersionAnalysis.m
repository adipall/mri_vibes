%
% Kevin Manktelow
% 10/25/2016
% 
% Dispersion analysis calculation using Salinas as a driver.  A mesh/unit
% cell is created in Cubit.  I use Salinas to write out the mass and
% stiffness matrices (obviously, this only works in Serial).  With this 
% information, I calculate the dispersion relation and plot it.  
%
% This script was written quickly in a short night without any checking, so
% still to-do is some simple verification.  But it seems to be doing
% something that looks correct.
%
function DispersionAnalysis()
    % Set some input parameters
    dir = './mfile_cube';
    mesh = 'cube.g';
    nFreq = 20;

    % Read in Mass & Stiffness matrix output from Salinas
    [K, M] = ImportSalinasMatrices(dir);
    
    % Read in exodus file
    f = exo_rd(fullfile(dir,mesh));
    C = f.Nodes.Coordinates;
    
    % By my definition, matching nodeset pairs are:
    %  Primary  Secondary
    %  -------------------------------
    %     1       2    (neg. X, pos. X)
    %     3       4    (neg. Y, pos. Y)
    %     5       6    (neg. Z, pos. Z)
    %  -------------------------------
    nsIdx1 = f.Nodesets(f.Nodesets.id2idx(1)).Nodes;
    nsIdx2 = f.Nodesets(f.Nodesets.id2idx(3)).Nodes;
    nsIdx3 = f.Nodesets(f.Nodesets.id2idx(5)).Nodes;
    edgeNodesXY = intersect(nsIdx1,nsIdx2);
    edgeNodesYZ = intersect(nsIdx2,nsIdx3);
    edgeNodesXZ = intersect(nsIdx1,nsIdx3);
    cornerXYZ = intersect(edgeNodesXZ,edgeNodesXY);
    
    % The last few arguments here are defining an "offset" unit cell length
    % that is used to match nodes via coordinates.  Probably a better way
    % to do this.
    %
    matchIdx1 = findMatchingNodes(nsIdx1, C, 1, 0, 0);
    matchIdx2 = findMatchingNodes(nsIdx2, C, 0, 1, 0);
    matchIdx3 = findMatchingNodes(nsIdx3, C, 0, 0, 1);
    matchIdxXY = findMatchingNodes(edgeNodesXY, C, 1, 1, 0);
    matchIdxYZ = findMatchingNodes(edgeNodesYZ, C, 0, 1, 1);
    matchIdxXZ = findMatchingNodes(edgeNodesXZ, C, 1, 0, 1);
    matchIdxXYZ = findMatchingNodes(cornerXYZ, C, 1, 1, 1);
    
    % Find matching nodes for Bloch theorem
    numNodes = length(f.Nodes.NodeNumMap);
    iKeep = setdiff(1:numNodes, [matchIdx1; matchIdx2; matchIdx3]);
    
    iX = setdiff(matchIdx1, [matchIdx2; matchIdx3]);
    iY = setdiff(matchIdx2, [matchIdx1; matchIdx3]);
    iZ = setdiff(matchIdx3, [matchIdx1; matchIdx2]);
    iXY = setdiff(matchIdxXY, matchIdx3);
    iYZ = setdiff(intersect(matchIdx2, matchIdx3), matchIdx1);
    iXZ = setdiff(intersect(matchIdx1, matchIdx3), matchIdx2);
    iXYZ = intersect(intersect(matchIdx1, matchIdx2), matchIdx3);
    
    % Sanity check to make sure each node is indexed only once.
    numReduce = length(iX)+length(iY)+length(iZ) ...
           + length(iXY)+length(iYZ)+length(iXZ)+length(iXYZ);
    assert( numReduce == length(unique([matchIdx1;matchIdx2;matchIdx3])));
           
    % Transfer matrices for nodes on:
    %   the X, Y, or Z plane (TX, TY, TZ)
    %   edges shared by planes (TXY, TYZ, TXZ)
    %   edges shared by all planes (TXYZ)
    T = zeros(numNodes, numNodes - numReduce);
    TX = T; TY = T; TZ = T;
    TXY = T; TYZ = T; TXZ = T;
    TXYZ = T;
    
    % Retained indices
    reducedIdx = 1:length(iKeep);
    Tidx = sub2ind(size(T), iKeep, reducedIdx);
    T(Tidx) = 1;
    
    % Matching X
    [~,~,idx] = intersect(iX, matchIdx1);
    [~,~,idx] = intersect(nsIdx1(idx), iKeep);
    Tidx = sub2ind(size(T), iX', reducedIdx(idx));
    TX(Tidx) = 1;

    % Matching Y
    [~,~,idx] = intersect(iY, matchIdx2);
    [~,~,idx] = intersect(nsIdx2(idx), iKeep);
    Tidx = sub2ind(size(T), iY', reducedIdx(idx));
    TY(Tidx) = 1;

    % Matching Z
    [~,~,idx] = intersect(iZ, matchIdx3);
    [~,~,idx] = intersect(nsIdx3(idx), iKeep);
    Tidx = sub2ind(size(T), iZ', reducedIdx(idx));
    TZ(Tidx) = 1;
    
    % Matching XY
    [~,~,idx] = intersect(iXY, matchIdxXY);
    [~,~,idx] = intersect(edgeNodesXY(idx), iKeep);
    Tidx = sub2ind(size(T), iXY', reducedIdx(idx));
    TXY(Tidx) = 1;
    
    % Matching YZ
    [~,~,idx] = intersect(iYZ, matchIdxYZ);
    [~,~,idx] = intersect(edgeNodesYZ(idx), iKeep);
    Tidx = sub2ind(size(T), iYZ', reducedIdx(idx));
    TYZ(Tidx) = 1;
    
    % Matching XZ
    [~,~,idx] = intersect(iXZ, matchIdxXZ);
    [~,~,idx] = intersect(edgeNodesXZ(idx), iKeep);
    Tidx = sub2ind(size(T), iXZ', reducedIdx(idx));
    TXZ(Tidx) = 1;

    % Matching XYZ
    [~,~,idx] = intersect(iXYZ, matchIdxXYZ);
    [~,~,idx] = intersect(cornerXYZ(idx), iKeep);
    Tidx = sub2ind(size(T), iXYZ', reducedIdx(idx));
    TXYZ(Tidx) = 1;
    
    T_all = zeros(size(T,1), size(T,2), 8);
    T_all(:,:,1) = T;
    T_all(:,:,2) = TX;
    T_all(:,:,3) = TY;
    T_all(:,:,4) = TZ;
    T_all(:,:,5) = TXY;
    T_all(:,:,6) = TYZ;
    T_all(:,:,7) = TXZ;
    T_all(:,:,8) = TXYZ;
    
    % Reduce the matrices by applying the Floquet-Bloch theorem
    Kr = ReduceMatrix( K, T_all );
    Mr = ReduceMatrix( M, T_all );
    
    % Set a number of wave-vector points to evaluate.  Here I'm only
    % varying muX for expediency of writing this script, but the general
    % approach is to vary the wave vector over the first brillouin zone.
    %
    muX = linspace(0.01, pi, 25);
    omega = zeros(length(muX), nFreq);
    
    %
    % Call MATLAB's Hermitian sparse eigenvalue solver once for each
    % wave-vector.  The result should be a real-valued eigenvalue, but may
    % have some small numerical imaginary part from the calculation.
    %
    disp('Solving eigenvalue problems ...')
    opts.tol = 1.0e-3;
    for i=1:length(muX)
        muY = 0;
        muZ = 0;
        
        [~,D] = eigs(Kr(muX(i), muY, muZ), Mr(muX(i), muY, muZ), nFreq, 'sm', opts);
        omega(i, :) = sqrt(real(diag(D)));
        omega(i, :) = sort(omega(i,:), 'ascend');
        
        disp(['  ', num2str(i), ' of ', num2str(length(muX)), ' completed ...']);
    end
    
    %
    % Plot the results
    %
    figure(1);
    clf;
    plot(muX, omega);
    xlabel('Wave vector, muX (rad)');
    ylabel('Frequency, \omega (rad/s)');
    xlim([0 pi]);
    
    % Estimated wave speed in LWL
    disp(['Estimated LWL wave speed: ', num2str(omega(1,1)/muX(1))]);
end

%
% Find nodes on a 'matching' plane, edge, or vertex.  I.e., locate the
% nodes that can be related to another node inside the unit cell via the
% Bloch theorem.
%
% I've chosen to define the lengths Lx, Ly, and Lz as the side-lengths of a
% rectangular prism unit cell.  Other unit cell shapes can be considered
% but the code is not generalized for that case.  
%
function MatchingNodes = findMatchingNodes( PrimaryNodes, Coords, Lx, Ly, Lz )
    N = length(PrimaryNodes);
    MatchingNodes = zeros(N,1);
    MatchCoords = Coords;
    
    MatchCoords(:,1) = Coords(:,1) - Lx;
    MatchCoords(:,2) = Coords(:,2) - Ly;
    MatchCoords(:,3) = Coords(:,3) - Lz;
        
    tol = max([Lx,Ly,Lz])*1e-6;
    for i=1:N
        Cprimary = Coords(PrimaryNodes(i),:);
        Cmatch = MatchCoords - repmat(Cprimary, size(Coords,1), 1);        
        idx = find( abs(Cmatch(:,1)) < tol & ...
                    abs(Cmatch(:,2)) < tol & ...
                    abs(Cmatch(:,3)) < tol );

        % there should be only one matching node
        assert(length(idx) == 1);
        
        MatchingNodes(i) = idx;
    end
end

%
% Applying the Floquet-Bloch theorem for a 3D array of unit cells.  Would
% need to be slightly modified for something that is only a 2D tesselation.
%
function Ar = ReduceMatrix( A, T_all )
    % Expand T_all
    d1 = size(T_all,1);   % total number of unit cell nodes
    d2 = size(T_all,2);   % number of unique unit cell nodes

    num_dof = 3;

    % Expand 'node-matching matrix' to DOF transfer matrix
    % Note: this assumes DOF ordering as
    %  {U} = [x1, y1, z1, x2, y2, z2, ..., xN, yN, zN] 
    % for node 1, 2, ..., N
    %
    T_all_tmp = zeros(d1*3, d2*3, 8);
    T_all_tmp(1:num_dof:(3*d1), 1:num_dof:(3*d2), :) = T_all;
    T_all_tmp(2:num_dof:(3*d1), 2:num_dof:(3*d2), :) = T_all;
    T_all_tmp(3:num_dof:(3*d1), 3:num_dof:(3*d2), :) = T_all;
    T_all = T_all_tmp;

    T =@(muX, muY, muZ) T_all(:,:,1) + ...
        T_all(:,:,2)*exp(1i*muX) + ...
        T_all(:,:,3)*exp(1i*muY) + ...
        T_all(:,:,4)*exp(1i*muZ) + ...
        T_all(:,:,5)*exp(1i*muX)*exp(1i*muY) + ...
        T_all(:,:,6)*exp(1i*muY)*exp(1i*muZ) + ...
        T_all(:,:,7)*exp(1i*muX)*exp(1i*muZ) + ...
        T_all(:,:,8)*exp(1i*muX)*exp(1i*muY)*exp(1.i*muZ);
    
    Ar =@(muX,muY,muZ) T(muX,muY,muZ)'*A*T(muX,muY,muZ);
end

%
% Read in Salinas' Mass() and Stiff() files.  This is an attempt to avoid
% m-file caching that is suspect for strange results that are not
% consistent with what the Mass.m and Stiff.m files should be outputting.
%
function [K, M] = ImportSalinasMatrices(dir)
    stiff_file = fullfile(dir, 'Stiff.m');
    mass_file = fullfile(dir, 'Mass.m');
    
    %======================================================
    % Read MASS matrix file
    %======================================================
    
    % Define sparse matrix data size
    num_data_rows = CountSparseMatrixRows(mass_file);
    array_size = [3, num_data_rows];
    data_fmt = ' %f %f %f';
    
    fid = fopen(mass_file, 'r');
    fgetl(fid);         % throw away 1st header line
    line = fgetl(fid);  % get number of entries from 2nd line
    fgetl(fid);         % throw away 3rd header line
    matrix_dim = str2double(strtok(line, 'NumRows = '));
    
    % Read data in and populate matrix
    sparse_dat = fscanf(fid, data_fmt, array_size)';
    
    % Close file
    fclose(fid);
    
    % Populate M and make symmetric
    M = sparse(sparse_dat(:,1), sparse_dat(:,2), sparse_dat(:,3), matrix_dim, matrix_dim);
    M = M + M' - diag(diag(M));
    
    %======================================================
    % Read STIFFNESS matrix file
    %======================================================
    % Define sparse matrix data size
    num_data_rows = CountSparseMatrixRows(stiff_file);
    array_size = [3, num_data_rows];
    data_fmt = ' %f %f %f';
    
    fid = fopen(stiff_file, 'r');
    fgetl(fid);         % throw away 1st header line
    line = fgetl(fid);  % get number of entries from 2nd line
    fgetl(fid);         % throw away 3rd header line
    matrix_dim = str2double(strtok(line, 'NumRows = '));
    
    % Read data in and populate matrix
    sparse_dat = fscanf(fid, data_fmt, array_size)';
    
    % Populate K and make symmetric
    K = sparse(sparse_dat(:,1), sparse_dat(:,2), sparse_dat(:,3), matrix_dim, matrix_dim);
    K = K + K' - diag(diag(K));
end

%
%   Open a file and count number of lines, then subtract # of header and
%   footer lines to calculate number of actual data lines for sparse
%   matrices.
%   
function num_data_rows = CountSparseMatrixRows(fname)
    fid = fopen(fname, 'r');    
    
    % Count number of data rows
    num_header = 3;
    num_footer = 3;
    num_rows = 0;
    
    while ~feof(fid)
        fgetl(fid);
        num_rows = num_rows + 1;
    end
    num_data_rows = num_rows - num_header - num_footer;
    
    fclose(fid);
end

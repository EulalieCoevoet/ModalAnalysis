function [V,F,U,q] = computeModalBasis(filename, varargin)
    % Input:
    %   filename  holding volume or surface mesh (in m)
    %             (.stl, .obj, .msh)
    %             if a volume mesh is provided, the corresponding surface mesh with the same
    %             name should be located next to it
    % Output: 
    %   .msh      hodling the generated volume mesh if a surface mesh is provided (generated by tetgen)
    %   .stl      hodling the computed surface mesh if a volume mesh is provided 
    %   .javabin  holding the modal basis (k modes) 
    %
    % Options:
    %   'Young'   followed by Young's modulus (Pa) of object (default 1e4)
    %   'Density' followed by density (kg/m3) of object (default 5000)
    %   'NbModes' followed by number of modes to compute (default 10)
    %   'Indices' followed by list of indices of fixed nodes (default none)
    %   'Boxes'   followed by list of axis aligned boxes ROI for fixed nodes (default none)
    %   'Draw'    followed by true or false to visualize the modes
    %
    % Requirement: 
    % 1- Install tetgen: http://wias-berlin.de/software/tetgen/ 
    
    if nargin<1
        error('Filename required. Should be a surface mesh ".stl" or ".obj".');
    end
    
    %%%% OPTIONIAL PARAMS %%%%
    %% Default values
    k = 10;
    indices = [];
    boxes = [];
    density = 5000;
    young = 1e4;
    draw = false;
    
    %% Map of parameter names to variable names
    params = containers.Map( ...
    {'NbModes','Boxes','Indices','Density','Young','Draw'},{'k','boxes','indices','density','young','draw'});
    v = 1;
    while v <= numel(varargin)
        param_name = varargin{v};
        if isKey(params, param_name)
              assert(v+1<=numel(varargin));
              v = v+1;
              % Trick: use feval on anonymous function to use assignin to this workspace 
              feval(@()assignin('caller',params(param_name),varargin{v}));
        else
             error('Unsupported parameter: %s',varargin{v});
        end
        v=v+1;
    end
    
    %% LOAD or, GENENERATE and EXPORT VOLUME MESH 
    [filepath,name,ext] = fileparts(filename);
    if ext == ".msh"
        [V,T,F] = readMSH(filename);
        if isfile(filepath+'/'+name+'.obj')
            [SV,SF] = loadMesh(filepath+'/'+name+'.obj');
        elseif isfile(filepath+'/'+name+'.stl')
            [SV,SF] = readSTL(filepath+'/'+name+'.stl','JoinCorners',1);
        else
            error('Surface mesh with name '+name+' not found in '+filepath);
            return;
        end        
    else
        if ext == ".stl"
            [SV,SF] = readSTL(filename,'JoinCorners',1); % Option is to remove duplicate vertices
        else
            [SV,SF] = loadMesh(filename);
        end
        filename = filepath+"/"+name+".msh";
        [V,T,F] = tetgen(SV, SF);
        writeMSH(filename,V,T,F);
    end

    %% SOLVE EIGEN VALUE PROBLEM 
    [K,M] = computeFEMLEMatrices( V, T, young, 0.48, density );
    K = 0.5*(K'+K);
    M = diag(M);
    MM = zeros(size(K,1),1);
    MM(1:3:end) = M;
    MM(2:3:end) = M;
    MM(3:3:end) = M;
    M = eye(size(K,1)).*MM;
    
    if ~isempty(boxes) % ROI for fixed points
        indices=boxROI(V,F,boxes); % Note this will find duplicate indices
    end
    
    if ~isempty(indices) % Fixed points
        ix = setdiff( 1:size(V,1)*3, [ indices*3-2, indices*3-1, indices*3 ] ); % Make sure to remove duplicated indices
        KF = K(ix,ix); % Remove fixed points from K and M
        MF = M(ix,ix);
        
        [UF, q] = eigs(KF, MF, k, 'smallestabs'); % Gets the k smallest eigen values and modes

        U = zeros(size(V,1)*3,k); % Set zeros on fixed points location
        U(ix,:) = UF;
    else % No boundary conditions 
        [U, q] = eigs(K, M, 6+k, 'smallestabs'); % Gets the k+6 smallest eigen values and modes
        % The first six eigenvalues are zero and the modes correspond to rigid translations and infinitesimal rotations
        % We discard these modes
        q = q(7:6+k,7:6+k);
        U = U(:,7:6+k);  
    end
    
    %% EXPORT 
    for i = 1:size(U,2)
        U(:,i) = U(:,i) ./ norm(U(:,i),2); % Normalize U
    end
    
    % Compute Center of Mass
    com = zeros(1, 3);
    totalMass = 0.0;
    for i = 1:size(V, 1)
        mi = M(3*i, 3*i);
        totalMass = totalMass + mi;
        com = com + mi * V(i, :);
    end
    com = com / totalMass;
    
    % Compute Rotation Inertia
    J = zeros(3, 3);
    for i = 1:size(V, 1)
        mi = M(3*i, 3*i);
        r = V(i, :) - com;
        J = J + mi * (r*r'*eye(3) - r'*r);
    end
    
    g = zeros(size(M, 1), 1);
    g(2:3:end) = -1;
    
    Mi = diag(U'*M*U);
    Ki = diag(U'*K*U); 
    U = U.*Ki'./Mi'/2.; % Scale columns U with w2/2
    UTMg = U'*M*g;
    Mi = diag(U'*M*U);
    Ki = diag(U'*K*U);
    SU = U(1:size(SV,1)*3,:);
    fos = java.io.FileOutputStream(strcat(filepath,'/',name,'.javabin'));
    oos = java.io.ObjectOutputStream(fos);
    
    oos.writeObject(SV);
    oos.writeObject(cast(SF,'int32'));
    oos.writeObject(Mi);
    oos.writeObject(totalMass);
    oos.writeObject([J(1,:) J(2,:) J(3,:)]); % Matrix3d in Java has a constructor which reads a vector[9]
    oos.writeObject(com);
    oos.writeObject(UTMg);
    oos.writeObject(Ki);
    oos.writeObject(SU);
    oos.close;
    fos.close;
    
    %% VISUALIZE MODES
    if draw
        figure;
        for i=1:k
            visualizeMode(V,F,U,q,i);
        end
    end
    
end

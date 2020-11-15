function [V,F,N] = readSTL(filename,varargin)
  % READSTL read a triangle mesh from an .stl file.
  %
  % [V,F] = readSTL(filename)
  % [V,F,N] = readSTL(filename,'ParameterName',ParameterValue, ...)
  %
  % Inputs:
  %   filename  path to .stl file
  %   Optional:
  %     'JoinCorners' followed by whether to join perfectly match corners to
  %     form a combinatorially connected mesh {false}.
  % Outputs:
  %   V  #V by 3 list of vertex positions
  %   F  #F by 3 list of faces
  %   N  #F by 3 list of face normals
  % 
  % Copyright https://github.com/alecjacobson/gptoolbox/blob/master/mesh/readSTL.m
  
  % default values
  join_corners = false;
  variable_name2 = [1,2,3];
  % Map of parameter names to variable names
  params_to_variables = containers.Map( {'JoinCorners'}, {'join_corners'});
  v = 1;
  while v <= numel(varargin)
    param_name = varargin{v};
    if isKey(params_to_variables,param_name)
      assert(v+1<=numel(varargin));
      v = v+1;
      % Trick: use feval on anonymous function to use assignin to this workspace
      feval(@()assignin('caller',params_to_variables(param_name),varargin{v}));
    else
      error('Unsupported parameter: %s',varargin{v});
    end
    v=v+1;
  end 


  is_ascii = false;
  fid = fopen(filename, 'r');
  header =fread(fid,80,'uchar=>schar'); % Read file title
  header = char(header(:)');
  is_ascii = startsWith(lower(header),'solid');
  fclose(fid);

  if is_ascii
    fid = fopen(filename, 'r');
    % discard header line
    fgets(fid);
    % The prefixing space is important here.
    D = fscanf(fid,' facet normal %f %f %f outer loop vertex %f %f %f vertex %f %f %f vertex %f %f %f endloop endfacet ');
    D = reshape(D,12,[])';
    N = D(:,1:3);
    V = reshape(D(:,4:12)',3,[])';
    F = reshape(1:size(V,1),3,size(V,1)/3)';
    fclose(fid);
  else
    [FVX,FVY,FVZ] = localstlread(filename);
    V = [FVX(:) FVY(:) FVZ(:)];
    F = reshape(1:size(V,1),3,size(V,1)/3)';
    N = [];
  end

  if join_corners
    [V,~,J] = remove_duplicate_vertices(V,0);
    % remap faces
    F = J(F);
  end

end

function [x, y, z, varargout] = localstlread(filename)
    % This function reads an STL file in binary format into matrixes X, Y and
    % Z, and C.  C is optional and contains color rgb data in 5 bits.  
    %
    % USAGE: [x, y, z, c] = stlread(filename);
    %
    % To plot use patch(x,y,z,c), or patch(x,y,z)
    %
    % Written by Doron Harlev
    % code copied from : https://github.com/alecjacobson/gptoolbox/

    if nargout>4
        error('Too many output arguments')
    end
    use_color=(nargout==4);

    fid=fopen(filename, 'r'); %Open the file, assumes STL Binary format.
    if fid == -1 
        error('File could not be opened, check name or path.')
    end

    ftitle=fread(fid,80,'uchar=>schar'); % Read file title
    num_facet=fread(fid,1,'int32'); % Read number of Facets

    %fprintf('\nTitle: %s\n', char(ftitle'));
    %fprintf('Num Facets: %d\n', num_facet);

    % Preallocate memory to save running time
    x=zeros(3,num_facet); y=zeros(3,num_facet); z=zeros(3,num_facet);
    if use_color
        c=uint8(zeros(3,num_facet));
    end

    %h = waitbar(0,'Please wait...');
    for i=1:num_facet,
        norm=fread(fid,3,'float32'); % normal coordinates, ignored for now
        ver1=fread(fid,3,'float32'); % vertex 1
        ver2=fread(fid,3,'float32'); % vertex 2
        ver3=fread(fid,3,'float32'); % vertex 3
        col=fread(fid,1,'uint16'); % color bytes
        if (bitget(col,16)==1 & use_color)
            r=bitshift(bitand(2^16-1, col),-10);
            g=bitshift(bitand(2^11-1, col),-5);
            b=bitand(2^6-1, col);
            c(:,i)=[r; g; b];
        end
        x(:,i)=[ver1(1); ver2(1); ver3(1)]; % convert to matlab "patch" compatible format
        y(:,i)=[ver1(2); ver2(2); ver3(2)];
        z(:,i)=[ver1(3); ver2(3); ver3(3)];
        if mod(i,floor(num_facet/10))==0
            %waitbar(i/num_facet,h);
        end
    end
    if use_color
        varargout(1)={c};
    end
    fclose(fid);
    %close(h);

    % For more information http://rpdrc.ic.polyu.edu.hk/old_files/stl_binary_format.htm

end

function [SV,SVI,SVJ] = remove_duplicate_vertices(V,epsilon,varargin)
  % REMOVE_DUPLICATE_VERTICES Remove duplicate vertices upto a uniqueness
  % tolerance (epsilon)
  %
  % SV = remove_duplicate_vertices(V,epsilon)
  % [SV,SVI,SVJ] = ...
  %   remove_duplicate_vertices(V,epsilon,'ParameterName',ParameterValue,...)
  %
  % Inputs:
  %   V  #V by dim list of vertex positions
  %   epsilon  uniqueness tolerance (significant digit)
  %   Optional:
  %     'WhiteList' Only merge vertices from the following selection (not
  %     working correctly, yet)
  % Outputs:
  %   SV  #SV by dim new list of vertex positions
  %   SVI #SV by 1 list of indices so SV = V(SVI,:) 
  %   SVJ #V by 1 list of indices so V = SV(SVJ,:)
  %
  % Example:
  %   % Mesh in (V,F)
  %   [SV,SVI,SVJ] = remove_duplicate_vertices(V,1e-7);
  %   % remap faces
  %   SF = SVJ(F);
  %
  % Copyright https://github.com/alecjacobson/gptoolbox/

  % default values
  whitelist = [];
  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'WhiteList'}, ...
    {'whitelist'});
  v = 1;
  while v <= numel(varargin)
    param_name = varargin{v};
    if isKey(params_to_variables,param_name)
      assert(v+1<=numel(varargin));
      v = v+1;
      % Trick: use feval on anonymous function to use assignin to this workspace 
      feval(@()assignin('caller',params_to_variables(param_name),varargin{v}));
    else
      error('Unsupported parameter: %s',varargin{v});
    end
    v=v+1;
  end

  if isempty(whitelist)
    assert(nargin==1 || epsilon >= 0);
    if nargin==1 || epsilon == 0
      [SV,SVI,SVJ] = unique(V,'rows','stable');
    else
      [~,SVI,SVJ] = unique(round(V/(10*epsilon)),'rows','stable');
      SV = V(SVI,:);
    end
  else
    error('not implemented correctly')
    VW = V(whitelist,:);
    J = 1:size(V,1);
    JW = J(whitelist);
    [SVW,SVIW,SVJW] = remove_duplicate_vertices(VW,epsilon);
    SJ = 1:size(V,1);
    SJ(whitelist) = JW(SVJ);
  end
end

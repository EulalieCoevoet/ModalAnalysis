function visualizeMode(V,F,U,q,index)
    % Visualizes a mode U(index) 
    % Input:
    %   V  list of vertices
    %   F  list of tetrahedra indices
    %   U  modal basis
    %   q  eigen values
    %   index  index of the mode to vizualize
    
    t=0;
    L=zeros(size(U,2),1);
    N=size(V,1);
    while t<1
        cla;
        L(index)=sin(sqrt(q(index,index))*t);
        UX = U*L;
        UX = [UX(1:3:N*3),UX(2:3:N*3),UX(3:3:N*3)] ;
        X=V+UX;
        plot_mesh(X,F);
        pause(0.1);
        t = t + 0.1;
    end
end
% Subroutine to generate meshes apriori and store them. This way redundant
% computations can be eliminated.

% all meshes are on the element (-1,1).
meshes = {'mesh1.mat','mesh2.mat','mesh3.mat','mesh4.mat',...
           'mesh5.mat','mesh6.mat'};
stimas = {'stima1.mat','stima2.mat','stima3.mat',...
             'stima4.mat','stima5.mat','stima6.mat'};

ref_fin = 6;
for ref = 1:ref_fin
    dx = 2^(-ref);
    vertex = -1:dx:1;
    center = 0.5*(vertex(1:end-1) + vertex(2:end))';
    N = length(vertex);
    cell_v = [[1:N-1]' [2:N]'];
    nodes_b = [1;N];
    area = [vertex(2:N) - vertex(1:N-1)]';
    ncell = N-1;
    nvert = N;
    
    filename =  meshes{ref};
    save(filename,'area','vertex','cell_v','nodes_b','center',...
                  'ncell','nvert');
end

for ref = 1:ref_fin
    
    loadmesh = ['load ' meshes{ref}];
	disp(loadmesh);
	eval(loadmesh);
	disp('mesh loaded');

    stima_grad = zeros(4*ncell,1);
    stima_mass = zeros(4*ncell,1);
    
     i_rows = zeros(4*ncell,1);
     i_cols = zeros(4*ncell,1);
    

   for i = 1:ncell
       stima_grad(4*i-3:4*i,1) = (area(i)^-1)*[1;-1;-1;1];
       stima_mass(4*i-3:4*i,1) = area(i)*[1/3;1/6;1/6;1/3];
       cur_nodes = cell_v(i,:);
       i_rows(4*i-3:4*i,1) = [cur_nodes(1) cur_nodes(1)...
                              cur_nodes(2) cur_nodes(2)]';
       i_cols(4*i-3:4*i,1) = [cur_nodes(1) cur_nodes(2)...
                              cur_nodes(1) cur_nodes(2)]';
   end
   
    filename =  stimas{ref};
    save(filename,'stima_grad','stima_mass','i_rows','i_cols');
end
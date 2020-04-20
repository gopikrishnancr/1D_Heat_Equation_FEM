% Time-Implicit FEM code for one dimensional heat equation
%
%   d_t u - div(a(x) grad(u)) = f(t,x),  
% with  Dirichlet boundary conditions
%               u(t,a) = ul(t),  u(t,b) = ur(t)
% and initial condition 
%               u(t
%------------------------------------------------------------------------
%   FEM discretisation:
%
%  (U_{n+1},v) - (U_{n},v) + dt*(grad(U_{n+1}),grad(v)) = (f(t_{n+1},.),v)
%

clearvars;
close all;

cfl = 0.5;  ref_fin = 6;  


% Meshes and stima matrices: all predefiend...assembly is using optimised
% 'sparse commands'.
meshes = {'mesh1.mat','mesh2.mat','mesh3.mat','mesh4.mat',...
           'mesh5.mat','mesh6.mat'};
stimas = {'stima1.mat','stima2.mat','stima3.mat',...
             'stima4.mat','stima5.mat','stima6.mat'};

% data functions
fex = @(t,x) exp(-t^2/2)*cos(x) - t*exp(-t^2/2)*cos(x);    % source
uex = @(t,x) exp((-t^2)/2).*cos(x);         % exact solution
a = @(x) 1;                                 % coefficcient function
duex = @(t,x) -exp(-t^2/2).*sin(x);         % spatial derivative of  exact sol.

for ref = 1:ref_fin        % main loop running over refinements
    
fprintf('Refinement %d started ....\n',ref);

    loadmesh = ['load ' meshes{ref}];
	disp(loadmesh);
	eval(loadmesh);
    
    loadmesh = ['load ' stimas{ref}];
	disp(loadmesh);
	eval(loadmesh);
  
    
    % CFL condition:   max(abs(a(x))*dt/dx <= cfl. Since we are using an
    % implicit scheme, CFL condition is not necessary. However, keep it so
    % just as an additional stability imposition.
    dx = max(area);        % spatial disc. factor
    dt = cfl*dx^2;         % dt based of CFL condition.
    tdis = 0:dt:1;  TN = length(tdis);    % temporal discretisation
    % storage matrix for the approximate solution. 'nvert' is the number of
    % vertices and 'TN' is the number of time nodes. 
    u_hd = zeros(nvert,TN);  
    u_hd(:,1) = uex(tdis(1),vertex)';     % initial condition
    
    %% Assebly of the coefficient matrix;
      C1 = zeros(4*ncell,1);
    for i = 1:ncell
        aval = a(center(i));
        C1(4*i - 3:4*i,1) = aval;
    end
    % A_grad  - \int a(x) (grad(u),grad(v))
    % A_mass -  \int (u,v)
    A_grad = sparse(i_rows,i_cols,C1.*stima_grad,nvert,nvert);
    A_mass = sparse(i_rows,i_cols,stima_mass,nvert,nvert);
  
    
    snodes = setdiff(1:nvert,nodes_b); 

for n = 1:TN-1     % main time loop
    
 % assembly of load function. We use a Boole's quadrature to compute the
 % integrals.
fload = zeros(nvert,1);      
    for i = 1:ncell
        curnodes = cell_v(i,:);
% Booles integration:
 f1 = int_bool(fex,vertex(curnodes(1)),vertex(curnodes(2)),1,tdis(n+1));
 f2 = int_bool(fex,vertex(curnodes(1)),vertex(curnodes(2)),2,tdis(n+1));
 fload(curnodes,:) = fload(curnodes,:) + dt*[f1;f2];
    end

 % imposition of the boundary data.
    for i = 1:length(nodes_b)
        u_hd(nodes_b(i),n+1) = uex(tdis(n+1),vertex(nodes_b(i)));
    end
 
 % implicit matrix operator
CM = dt*A_grad + A_mass;
%-------------------------------------
% load matrix adjusted for explicit terms and boundary values.
fload = fload + A_mass*u_hd(:,n);
fload = fload - CM*u_hd(:,n+1);
u_hd(snodes,n+1) = CM(snodes,snodes)\fload(snodes);
    
end

figure(ref);    
uval = uex(tdis(end),vertex);

fprintf('plotting ...\n');
subplot(1,2,1)
plot(vertex,uval,'b-');
title('Exact solution');
subplot(1,2,2)
plot(vertex,u_hd(:,end),'r-');
title('Approximate solution');

fprintf('Refinement %d completed ....\n\n',ref);
[L2(ref),H1(ref)] = erros(u_hd(:,end),uex,duex,...
                    vertex,area,cell_v,ncell,center,tdis(end));
DX(ref) = max(abs(area));
end


if ref_fin > 1
    for ref = 2:ref_fin
        ocl2(ref-1) = log(L2(ref)/L2(ref-1))/log(DX(ref)/DX(ref-1));
        och1(ref-1) = log(H1(ref)/H1(ref-1))/log(DX(ref)/DX(ref-1));
    end
end


L2
H1
DX
ocl2
och1








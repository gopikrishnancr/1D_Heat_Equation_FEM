function [L2,H1] = erros(u_app,u,du,vertex,area,cell_v,ncell,center,t)

L2e = 0; H1e = 0;
for i = 1:ncell
    curnodes = cell_v(i,:);
    vert = vertex(curnodes);
    
ua = linspace(u_app(curnodes(1)),u_app(curnodes(2)),5); 
ue = u(t,linspace(vert(1),vert(2),5));
ud = dot([7,32,12,32,7],(ua - ue).^2);

   L2e = L2e + ud*(2*area(i)/45);
%    L2e = L2e + area(i)*((u_app(curnodes(1)) - u(t,vert(1)))^2 +...
%         (u_app(curnodes(2)) - u(t,vert(2)))^2 +...
%            4*(0.5*(u_app(curnodes(2)) + u_app(curnodes(2))) - u(t,0.5*(vert(1) + vert(2))))^2)/3;
    
%     L2e = L2e + ((u_app(curnodes(1)) - u(t,vert(1)))^2 +...
%         (u_app(curnodes(2)) - u(t,vert(2)))^2)*area(i)/2;
%     

    duh = (u_app(curnodes(2)) - u_app(curnodes(1)))/area(i);
    du_diff = (du(t,linspace(vert(1),vert(2),5)) - duh).^2;
    udd = dot([7,32,12,32,7],du_diff);
    H1e = H1e + udd*(2*area(i)/45);
  %  H1e = H1e + area(i)*(du(t,center(i)) - (u_app(curnodes(2)) - u_app(curnodes(1)))/area(i))^2;
end
L2 = sqrt(L2e);
H1 = sqrt(H1e);

function [grad] =laplacian(nx,ny,dx,dy)

format long;

nxny=nx*ny;

r=zeros(1,nx);
r(1:2)=[2,-1];
T=toeplitz(r);

E=speye(nx);

grad=-(kron(T,E)+kron(E,T));

%-- for periodic boundaries

for i=1:nx
ii=(i-1)*ny+1; %left boundary
jj=ii+ny-1; %right boundary
%grad(ii,jj)=1.0;
%grad(jj,ii)=1.0;
%grad(jj,jj-1) = 1.0;Why have I written this? - It already was 1 - would not have made a difference
%grad(ii,ii+1) = 1.0; Why have I written this??
kk=nxny-nx+i;%bottom boundary
grad(i,kk)= 1.0;
grad(kk,i)= 1.0;
grad(jj,jj) = -3.0; %Condition you get after puttting Neumann boundary condition for phi on right boundary. You do not get do for the left boundary as I do not update phi at the left boundary.
grad(ii,ii) = -3.0;%neumann on left
end

grad = grad /(dx*dy);

end %endfunction


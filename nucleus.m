function [phi,tempr] = nucleus(Nx,Ny,inittemp)

format long;

for i=1:Nx
for j=1:Ny

tempr(i,j) = inittemp;

end
end
phi = zeros(Nx,Ny);
Nx2 = Nx/2;
phi((Nx2-2):(Nx2+2),(Nx2-2):(Nx2+2)) = 1;
phi(1:Nx,1) = 1;
phi(1:Nx,Ny) = 1;
phi(1,1:Ny) = 1;
phi(Nx,1:Ny) = 1;
end %endfunction

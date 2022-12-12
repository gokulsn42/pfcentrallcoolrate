%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    PHASE-FIELD FINITE-DIFFRENCE CODE FOR      %
%             DENDRITIC SOLIDIFICATION          %
%           (OPTIMIZED FOR MATLAB/OCTAVE)       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pkg load io;
%== get intial wall time:
time0=clock();
format long;
noisemag = 0.1;
%-- Simulation cell parameters:
%coolrate = 1000000; %unit is K/s : 1e5 is fair acc to literature as well
coolgrad = 100000; %unit is K/m
coolgradright = 100000; %unit is K/m
coolgradtop = 100000;
coolgradbottom = 100000;
Nx = 100;
Ny = 100;
NxNy= Nx*Ny;
A = sparse(NxNy,NxNy);
Id = speye(NxNy,NxNy);
B = sparse(NxNy,1);
dx = 1.5e-8;
dy = 1.5e-8;
inthickness = 1.5e-8
x = zeros(NxNy,1);
%--- Time integration parameters:
nstep  =  54000;
nprint = 1800;
dtime  = 1e-11;
%--- Material specific parameters:
W = 1.52e8;
tau = 1/42.948;
epsilonb = 0.0001311;
mu = 1.0;
kappa = 439;
delta = 0.03;
aniso = 4.0;
alpha = 72.13;
gamma = 146.5;
teq = 933.15;
theta0 = 0.2;
seed = 10.0;
td = 9.6e-5 %Thermal diffusivity
%
pix=4.0*atan(1.0);
init = 1;
if(init==1)
%--- Initiliaze and introduce initial nuclei:
inittemp = 850;
[phi,tempr] = nucleus(Nx,Ny,inittemp);
endif
if(init != 1)
  phiname = ["time_" num2str(init) ".xlsx"];
  phi = xlsread(phiname);
  tempname = ["temperature_" num2str(init) ".xlsx"];
  tempr = xlsread(tempname);
endif
%--- Laplacian templet:
[laplacian] =laplacian(Nx,Ny,dx,dy);
%--- Evolution
[A,lftbdry,rightbdry,topbdry,bottombdry] = thermalinit(A,Nx,Ny,dx,dy,dtime,coolgrad,coolgradright,coolgradtop,coolgradbottom,td);
fprintf('Hello World');
for istep = init:nstep
phiold = phi;
%---
% calculate the laplacians and epsilon:
%---

phi2 = reshape(phi',NxNy,1);

lap_phi2 = laplacian*phi2;

[lap_phi]=vec2matx(lap_phi2,Nx);

%--gradients of phi:

[phidy,phidx]=gradient_mat_phi(phi,Nx,Ny,dx,dy);
%The gradient_mat function is only used to calculate the gradient of phi.
%-- calculate angle:

theta =atan2(phidy,phidx);

%--- epsilon and its derivative:

epsilon = epsilonb*(1.0+delta*cos(aniso*(theta-theta0)));

epsilon_deriv = -epsilonb*aniso*delta*sin(aniso.*(theta-theta0));

%--- first term:

dummyx =epsilon.*epsilon_deriv.*phidx;

[term1,dummy] =gradient_mat(dummyx,Nx,Ny,dx,dy);

%--- second term:

dummyy =-epsilon.*epsilon_deriv.*phidy;

[dummy,term2] = gradient_mat(dummyy,Nx,Ny,dx,dy);

%--- factor m:

m = (0.9/3.1415)*atan(gamma*(teq-tempr)/teq);%This is slight cheating....
if(max(abs(m))>0.5)
    disp('istep');
    disp(istep);
    return;
endif
psi = -0.5 + rand(Nx,Ny);
%-- Time integration:

phi = phi +(dtime/tau) *(term1 +term2 + epsilon.^2 .* lap_phi + ...
			 W*phiold.*(1.0-phiold).*(phiold - 0.5 + m) +...
			W*noisemag*phiold.*(1.0-phiold).*psi);

%-- evolve temperature:
% DO A MATRIX UNROLLING HERE
phioldvec = reshape(phiold',[NxNy,1]);
phivec = reshape(phi',[NxNy,1]);
temprvec = reshape(tempr',[NxNy,1]);
fprintf('hi');
B = (kappa/dtime)*(phivec-phioldvec) + temprvec/dtime + lftbdry + rightbdry + topbdry + bottombdry;
x = A\B;
tempr = reshape(x,[Ny,Nx]);
tempr = tempr';
fprintf('done step: %5d\n',istep);
%---- print results

if(mod(istep,nprint) == 0  || istep == 20 || istep == 1)
  fname = ["time_" num2str(istep) ".jpg"];
  f1name = ["temperature_" num2str(istep) ".jpg"];
  f = figure('visible','off');
  imagesc(phi);
  colorbar;
  saveas(f,fname);
  f1 = figure('visible','off');
  imagesc(tempr);
  colorbar;
  saveas(f1,f1name);
  f4name = ["time_" num2str(istep) ".xlsx"];
  f5name = ["temperature_" num2str(istep) ".xlsx"];
  f6name = ["laplacian_phi_" num2str(istep) ".xlsx"];
  f7name = ["phidx_" num2str(istep) ".xlsx"];
  f8name = ["phidy_" num2str(istep) ".xlsx"];
  xlswrite(f4name,phi);
  xlswrite(f5name,tempr);
  xlswrite(f6name,lap_phi);
  xlswrite(f7name,phidx);
  xlswrite(f8name,phidy);
end %if
end %istep
%--- calculate compute time:
compute_time = etime(clock(),time0);
fprintf('Compute Time: %10d\n',compute_time);

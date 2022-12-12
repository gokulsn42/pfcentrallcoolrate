function [A,lftbdry,rightbdry,topbdry,bottombdry] = thermalinit(A,Nx,Ny,dx,dy,dtime,coolgrad,coolgradright,coolgradtop,coolgradbottom,td)
  format long;
  A = sparse(NxNy+2,NxNy+2);
  for i= 0:Nx-1
  for j= 1:Ny
    ii = Ny*i+j;
      im = Ny*(i-1)+j;
      ip = Ny*(i+1)+j;
      jm = Ny*i+j-1;
      jp = Ny*i+j+1;
    A(ii,ii) = (1/dtime) + (4*td/(dx*dx));
    A(ii,im) = -(td/(dx*dx));
    A(ii,ip) = -(td/(dx*dx));
    A(ii,jp) = -(td/(dx*dx));
    A(ii,jm) = -(td/(dx*dx));
    if(i==0)
      A(ii,ii) = (1/dtime) + (3*td/(dx*dx));
      if(j==1) %top left corner
        A(ii,ii) = (1/dtime) + (2*td/(dx*dx));
        lftbdry(ii) = -coolgrad*td/dx;
        A(ii,ip) = -(td/(dx*dx));
        A(ii,jp) = -(td/(dx*dx));
        A(ii,jm) = 0;
        A(ii,im) = 0;
      endif
      if(j==Ny) %top right corner
        A(ii,ii) = (1/dtime) + (2*td/(dx*dx));
        rightbdry(ii) = -coolgradright*td/dx;
        A(ii,ip) = -(td/(dx*dx));
        A(ii,jm) = -(td/(dx*dx));
        A(ii,jp) = 0;
        A(ii,im) = 0;
      endif
      topbdry(ii) = -coolgradtop*td/dx;
    endif
    if(i==(Nx-1))
    bottombdry(ii) = -coolgradbottom*td/dx;
    A(ii,ii) = (1/dtime) + (3*td/(dx*dx));
    A(ii,jp) = -(td/(dx*dx));
    A(ii,jm) = -(td/(dx*dx));
    A(ii,im) = -(td/(dx*dx));
      if(j==Ny) %bottom right corner
        A(ii,ii) = (1/dtime) + (2*td/(dx*dx));
        rightbdry(ii) = -coolgradright*td/dx;
        A(ii,im) = -(td/(dx*dx));
        A(ii,jm) = -(td/(dx*dx));
        A(ii,jp) = 0;
        A(ii,ip) = 0;
      endif
      if(j==1) %bottom left corner
        A(ii,ii) = (1/dtime) + (2*td/(dx*dx));
        A(ii,jp) = -(td/(dx*dx));
        A(ii,im) = -(td/(dx*dx));
        A(ii,ip) = 0;
        A(ii,jm) = 0;
        lftbdry(ii) = -coolgrad*td/dx;
      endif
    endif
    if(j==Ny && i!=0 && i!=(Nx-1)) %Right boundary. Corners have been handled above and excluded here
    A(ii,ii) = (1/dtime) + (3*td/(dx*dx));
    A(ii,im) = -(td/(dx*dx));
    A(ii,ip) = -(td/(dx*dx));
    A(ii,jm) = -(td/(dx*dx));
    rightbdry(ii) = -coolgradright*td/dx;%be mindful of the sign here!!!
    endif
    if(j==1 && i!=0 && i!=(Nx-1)) %Left boundary. Corners have been handled above and excluded here
      A(ii,ii) = (1/dtime) + (3*td/(dx*dx));
      A(ii,im) = -(td/(dx*dx));
      A(ii,ip) = -(td/(dx*dx));
      A(ii,jp) = -(td/(dx*dx));
      lftbdry(ii) = -coolgrad*td/dx;
    endif
  end
  end
end %function

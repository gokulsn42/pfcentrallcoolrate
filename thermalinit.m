function [A,lftbdry,rightbdry,topbdry,bottombdry] = thermalinit(A,Nx,Ny,dx,dy,dtime,coolgrad,coolgradright,coolgradtop,coolgradbottom,td)
  format long;
  NxNy = Nx*Ny;
  lftbdry = sparse(NxNy,1);
  rightbdry = sparse(NxNy,1);
  topbdry = sparse(NxNy,1);
  bottombdry = sparse(NxNy,1);
  for i= 0:Nx-1
  for j= 1:Ny
    ii = Ny*i+j;
    if(i!=0)
      im = Ny*(i-1)+j;
    endif
    if(i!=(Nx-1))
      ip = Ny*(i+1)+j;
    endif
    if(j!=1)
      jm = Ny*i+j-1;
    endif
    if(j!=Ny)
      jp = Ny*i+j+1;
    endif
    if(j!=Ny && j!=1 && i!=0 && i!=(Nx-1))
    A(ii,ii) = (1/dtime) + (4*td/(dx*dx));
    A(ii,im) = -(td/(dx*dx));
    A(ii,ip) = -(td/(dx*dx));
    A(ii,jp) = -(td/(dx*dx));
    A(ii,jm) = -(td/(dx*dx));
    endif
    if(i==0)
      if(j==1) %top left corner
        A(ii,ii) = (1/dtime) + (2*td/(dx*dx));
        lftbdry(ii) = -coolgrad*td/dx;
        A(ii,ip) = -(td/(dx*dx));
        A(ii,jp) = -(td/(dx*dx));
      endif
      if(j==Ny) %top right corner
        A(ii,ii) = (1/dtime) + (2*td/(dx*dx));
        rightbdry(ii) = -coolgradright*td/dx;
        A(ii,ip) = -(td/(dx*dx));
        A(ii,jm) = -(td/(dx*dx));
      endif
      if(j!=1 && j!=Ny)
      A(ii,ii) = (1/dtime) + (3*td/(dx*dx));
      A(ii,ip) = -(td/(dx*dx));
      A(ii,jp) = -(td/(dx*dx));
      A(ii,jm) = -(td/(dx*dx));
      endif
      topbdry(ii) = -coolgradtop*td/dx;
    endif
    if(i==(Nx-1))
    bottombdry(ii) = -coolgradbottom*td/dx;
      if(j==Ny) %bottom right corner
        A(ii,ii) = (1/dtime) + (2*td/(dx*dx));
        rightbdry(ii) = -coolgradright*td/dx;
        A(ii,im) = -(td/(dx*dx));
        A(ii,jm) = -(td/(dx*dx));
      endif
      if(j==1) %bottom left corner
        A(ii,ii) = (1/dtime) + (2*td/(dx*dx));
        A(ii,jp) = -(td/(dx*dx));
        A(ii,im) = -(td/(dx*dx));
        lftbdry(ii) = -coolgrad*td/dx;
      endif
      if(j!=1 && j!=Ny)
        A(ii,ii) = (1/dtime) + (3*td/(dx*dx));
        A(ii,jp) = -(td/(dx*dx));
        A(ii,jm) = -(td/(dx*dx));
        A(ii,im) = -(td/(dx*dx));
      endif
    endif
    if(j==Ny && i!=0 && i!=(Nx-1))
    A(ii,ii) = (1/dtime) + (3*td/(dx*dx));
    A(ii,im) = -(td/(dx*dx));
    A(ii,ip) = -(td/(dx*dx));
    A(ii,jm) = -(td/(dx*dx));
    rightbdry(ii) = -coolgradright*td/dx;%be mindful of the sign here!!!
    endif
    if(j==1 && i!=0 && i!=(Nx-1))
      A(ii,ii) = (1/dtime) + (3*td/(dx*dx));
      A(ii,im) = -(td/(dx*dx));
      A(ii,ip) = -(td/(dx*dx));
      A(ii,jp) = -(td/(dx*dx));
      lftbdry(ii) = -coolgrad*td/dx;
    endif
  end
  end
end %function

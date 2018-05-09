%Final Project Helmhotz 2-D:
%d^2u/dx^2+d^2u/dy^2+Lambda*u=F(x,y)
%Over domain of ax<x<bx and ay<y<by
%Using 2 methods of Gauss elimination, Gauss-Seidel, 
%or Successive Over Relaxation

clc
clear;

%Given variables
ucoef1=1/2;             %Given
ax=0;                   %Given
ay=0;                   %Given
Pi=4*atan(1);           %Given
bx=2*Pi;                %Given
by=2*Pi;                %Given
v=0;                    %Given (du/dy @y=by = 0)

%Initial Problem Setup
iterate=20000;
Nx=12;                  %Initial number of points in x
Ny=Nx;                  %Initial number of points in y
Lx=bx;                  %Length of x
Ly=by;                  %Length of y
deltax=Lx/(Nx+1);           %Step size in x
deltay=Ly/(Ny+1);           %Step size in y


%Constants
fbay = (by-ay)*(by-ay)*cos(Pi*ay/by);
gbay = ay*(by-ay)*(by-ay);
cons1 = bx-ax;
Cons2=deltay*deltay;
Cons3=deltax*deltax;
Cons4=deltax*deltax*deltay*deltay;
Cons5=-2*deltay*deltay-2*deltax*deltax;
Cons5a=Cons5+Cons4*ucoef1;

%Solution Matrix
U1 = zeros(Nx+2,Ny+2);	% Create solution matix with initial guess of ZERO

%Boundary Conditions
c=0;                    %counter
for j = 1:Ny+2                      %u(x=ax,y)=fb(y) boundary condition
    U1(1,j)=(by-deltay*c)*(by-deltay*c)*cos(Pi*deltay*c/by);
    c=c+1;
end
c=0;

for j=1:Ny+2                        %u(x=bx,y)=gb(y) boundary condition
    U1(Ny+2,j)=deltay*c*(by-deltay*c)*(by-deltay*c);
    c=c+1;
end
c=0;

for i=1:Nx+2                        %u(x,y=ay)  boundary conditoin
    U1(i,1)=fbay+(deltax*c-ax)/cons1*(gbay-fbay);
    c=c+1;
end
c=0;

for i=1:Nx+2                        %du/dy @y=by = 0 boundary condition
    U1(i,Nx+2)=2*deltay*c*v;
    c=c+1;
end

%Part 1 Using Gauss Siedel
F1=zeros(Nx+2,Ny+2);

for i=2:Nx+2
    for j=2:Ny+2
        F1(i,j)=cos(Pi/2*(2*((i-1)-ax)/(bx-ax)+1))*sin(Pi*((j-1)-ay)/(by-ay));
    end
end

for z=1:50000
for i=2:Nx+1
    for j=2:Ny+1
        %Checking U w/Poisson
        U1(i,j)=(1/4)*(U1(i-1,j)+U1(i+1,j)+U1(i,j-1)+U1(i,j+1))-((deltax*deltax/4)*F1(i,j)-ucoef1*deltax*deltax);
        %U1(i,j)=(Cons2.*U1(i+1,j)+Cons2.*U1(i-1,j)+Cons3.*U1(i,j+1)+Cons3.*U1(i,j-1)-Cons4.*F1(i,j))/(2*deltay*deltay+2*deltax*deltax-deltax*deltax*deltay*deltay*ucoef1);
    end
end
end

%Part 2
%Boundary Conditions
ucoef2=0;
Cons5b=Cons5+Cons4*ucoef2;
U2 = zeros(Nx+2,Ny+2);


c=0;                    %counter
for j = 1:Ny+2                      %u(x=ax,y)=fb(y) boundary condition
    U2(1,j)=(by-deltay*c)*(by-deltay*c)*cos(Pi*deltay*c/by);
    c=c+1;
end
c=0;

for j=1:Ny+2                        %u(x=bx,y)=gb(y) boundary condition
    U2(Ny+2,j)=deltay*c*(by-deltay*c)*(by-deltay*c);
    c=c+1;
end
c=0;

for i=1:Nx+2                        %u(x,y=ay)  boundary conditoin
    U2(i,1)=fbay+(deltax*c-ax)/cons1*(gbay-fbay);
    c=c+1;
end
c=0;

for i=1:Nx+2                        %du/dy @y=by = 0 boundary condition
    U2(i,Nx+2)=2*deltay*c*v;
    c=c+1;
end

F2=zeros(Nx+2,Ny+2);

for i=2:Nx+1
    for j=2:Ny+1
       U2(i,j)=(F2(i,j)-(deltay*deltay).*U2(i-1,j)-(deltay*deltay).*U2(i+1,j)-(deltax*deltax).*U2(i,j-1)-(deltax*deltax).*U2(i,j+1))/(-2*deltay*deltay-2*deltax*deltax+deltax*deltax*deltay*deltay*ucoef2);
    end
end
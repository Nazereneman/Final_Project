%Final Project Poisson 2-D:
%d^2u/dx^2+d^2u/dy^2=F(x,y)
%%%Old for Helmholtz%%%
%d^2u/dx^2+d^2u/dy^2+Lambda*u=F(x,y)
%%%%%%%%%%%%%%%%%%%%%%%
%Over domain of ax<x<bx and ay<y<by
%Using Gauss-Seidel

clc
clear;
%close all;

if exist('checkpt_GS.mat','file')     %If a checkpoint file exists, open it
    load('checkpt_GS.mat')
end

%Given variables
%lamda1=0.5;            %Given when using Helhmohtz, not needed for Poisson
ax=0;                   %Given
ay=0;                   %Given
Pi=4*atan(1);           %Given
bx=2*Pi;                %Given
by=2*Pi;                %Given
v=0;                    %Given (du/dy @y=by = 0)

%Initial Problem Setup 
Lx=bx-ax;                  %Length of x
Ly=by-ay;                  %Length of y
Nx=60;                  %Initial number of points in x
Ny=Nx;                %Initial number of points in y
deltax=Lx/(Nx+1);           %Step size in x
deltay=Ly/(Ny+1);           %Step size in y
    

%Constants 
fbay = (by-ay)*(by-ay)*cos(Pi*ay/by);
gbay = ay*(by-ay)*(by-ay);
cons1 = bx-ax;
dy2=deltay*deltay;
dx2=deltax*deltax;
dx2dy2=deltax*deltax*deltay*deltay;
denomin=-2*deltay*deltay-2*deltax*deltax;
Pidy=Pi*deltay;
gfbay=(gbay-fbay)/cons1;
halfPi=Pi/2;

%Initialize x and y to graph
x=0:deltax:2*pi;
y=0:deltay:2*pi;

%Error info 
tol=10^-8;                                      %Pick tolerance required


%Solution Matrix
U1 = zeros(Nx+2,Ny+2);	% Preallocate solution matix with initial guess of ZERO

%Boundary Conditions
for j = 1:Ny+2                      
    U1(1,j)=(by-deltay*(j-1))*(by-deltay*(j-1))*cos(Pidy*(j-1)/by);     %u(x=ax,y)=fb(y) boundary condition
    U1(Nx+2,j)=(deltay*(j-1))*(by-deltay*(j-1))*(by-deltay*(j-1));      %u(x=bx,y)=gb(y) boundary condition
    U1(j,1)=fbay+(deltax*(j-1)-ax)*gfbay;                               %u(x,y=ay)  boundary condition
end

%Error
err1=1000;                          %Initialize error as 1000
iter1=0;                            %Set initial iteration to 0
Uprev1=U1;                          %Set Uprev = U, for iteration error
    
%Part 1 Using Gauss Siedel
F1=zeros(Nx+2,Ny+2);                %Preallocate F1

for i=1:Nx+2
    for j=1:Ny+2
        F1(i,j)=cos(halfPi*(2*(deltax*(i-1)-ax)/(Lx)+1))*sin(Pi*(deltay*(j-1)-ay)/(Ly));
    end
end


while err1>tol
    for i=2:Nx+1
        for j=2:Ny+1
            %Solving for U
            U1(i,j)=(dx2dy2*F1(i,j)-dy2*U1(i-1,j)-dy2*U1(i+1,j)-dx2*U1(i,j-1)-dx2*U1(i,j+1))/(denomin);
        end
        %Neumann Condition
        U1(i,Ny+2)=(dx2dy2*F1(i,Ny+2)-dy2*U1(i-1,Ny+2)-dy2*U1(i+1,Ny+2)-dx2*U1(i,Ny+1)-dx2*U1(i,Ny+1))/(denomin);
    end
    err1=max(max(abs(Uprev1-U1)./abs(Uprev1)));
    Uprev1=U1;
    iter1=iter1+1;
    if mod(iter1,1500)==0                               %Save checkpoint file every 1500 iterations
        save('checkpt_GS.mat');                   %Saving the file
    end                                                 %Ending if loop
end

Uaverage=sum(sum(U1))/((Nx+2)*(Ny+2));
U1T=transpose(U1);      %Transpose the matrix so x and y axes are correct

figure()
h1=surf(x,y,U1T);
xlabel('x')
ylabel('y')
set(h1,'linestyle','none')
zlabel('U values');
title('Surface Plot for GS (Nx=Ny=60) w/1st Order Neumann');
colorbar;

figure()
contour(x,y,U1T);
xlabel('x')
ylabel('y')
set(h1,'linestyle','none')
zlabel('U values');
title('Contour Plot for GS (Nx=Ny=60) w/1st Order Neumann');
colorbar;

%%
delete('checkpt_GS.mat');             %Delete checkpoint file once complete
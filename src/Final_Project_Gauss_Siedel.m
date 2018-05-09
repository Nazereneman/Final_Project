%Final Project Poisson 2-D:
%d^2u/dx^2+d^2u/dy^2=F(x,y)
%%%Old for Helmholtz%%%
%d^2u/dx^2+d^2u/dy^2+Lambda*u=F(x,y)
%%%%%%%%%%%%%%%%%%%%%%%
%Over domain of ax<x<bx and ay<y<by
%Using Gauss-Seidel

clc
clear;

%Given variables
%lamda1=0.5;            %Given when using Helhmohtz, not needed for Poisson
ax=0;                   %Given
ay=0;                   %Given
Pi=4*atan(1);           %Given
bx=2*Pi;                %Given
by=2*Pi;                %Given
v=0;                    %Given (du/dy @y=by = 0)

%Initial Problem Setup that doesn't depend on # of points (ie size of Nx,Ny)
Lx=bx;                  %Length of x
Ly=by;                  %Length of y
Nx1=6;                  %Initial number of points in x
Ny1=Nx1;                %Initial number of points in y

%Constants that do not depend on the # of points (ie size of Nx,Ny)
fbay = (by-ay)*(by-ay)*cos(Pi*ay/by);
gbay = ay*(by-ay)*(by-ay);
cons1 = bx-ax;

%Error info that does not depend on the #of points (ie size of Nx,Ny)
k=0;                                            %Multiple of 2 being multiplied by Nx and Ny
grid_error=100;                                 %Set initial grid error greater than tol
tol=10^-8;                                      %Pick tolerance required
Uold1=zeros(Nx1+2,Ny1+2);                          %Allocate Uold

while grid_error > tol
    %Variables that change based on the number of points (ie size of Nx,Ny)
    if k==0
    Nx=2^(k)*Nx1;              %Double # of points in x for each iteration
    Ny=2^(k)*Ny1;              %Double # of points in y for each iteration
    else
    Nx=Nx*2+2
    Ny=Ny*2+2
    end
    OldNx=.5*Nx;
    OldNy=.5*Ny;
    deltax=Lx/(Nx+1);           %Step size in x
    deltay=Ly/(Ny+1);           %Step size in y
    
    %Constants that depend on the number of points (ie size of Nx,Ny)
    dy2=deltay*deltay;
    dx2=deltax*deltax;
    dx2dy2=deltax*deltax*deltay*deltay;
    denomin=-2*deltay*deltay-2*deltax*deltax;
    %Initialize x and y to graph
    x=0:deltax:2*pi;
    y=0:deltay:2*pi;

    %Solution Matrix
    U1 = zeros(Nx+2,Ny+2);	% Preallocate solution matix with initial guess of ZERO

    %Boundary Conditions
    for j = 1:Ny+2                      %u(x=ax,y)=fb(y) boundary condition
        U1(1,j)=(by-deltay*(j-1))*(by-deltay*(j-1))*cos(Pi*deltay*(j-1)/by);
    end

    for j=1:Ny+2                        %u(x=bx,y)=gb(y) boundary condition
        U1(Nx+2,j)=(deltay*(j-1))*(by-deltay*(j-1))*(by-deltay*(j-1));
    end

    for i=1:Nx+2                        %u(x,y=ay)  boundary conditoin
        U1(i,1)=fbay+(deltax*(i-1)-ax)/cons1*(gbay-fbay);
    end

    %Error
    err1=1000;                          %Initialize error as 1000
    iter1=0;                            %Set initial iteration to 0
    Uprev1=U1;                          %Set Uprev = U, for iteration error
    
    %Part 1 Using Gauss Siedel
    F1=zeros(Nx+2,Ny+2);                %Preallocate F1

    for i=1:Nx+2
        for j=1:Ny+2
            F1(i,j)=cos(Pi/2*(2*(deltax*(i-1)-ax)/(bx-ax)+1))*sin(Pi*(deltay*(j-1)-ay)/(by-ay));
        end
    end


    while err1>tol
        for i=2:Nx+1
            for j=2:Ny+1
                %Solving for U
                U1(i,j)=(dx2dy2*F1(i,j)-dy2*U1(i-1,j)-dy2*U1(i+1,j)-dx2*U1(i,j-1)-dy2*U1(i,j+1))/(denomin);
            end
            %Neumann Condition
            U1(i,Ny+2)=(dx2dy2*F1(i,Ny+2)-dy2*U1(i-1,Ny+2)-dy2*U1(i+1,Ny+1)-dx2*U1(i,Ny+1)-dy2*U1(i,Ny+1))/(denomin);
        end
        err1=max(max(((abs(U1)-abs(Uprev1))./abs(U1))*100));
        Uprev1=U1;
        iter1=iter1+1;
    end
    
    griderror = zeros(OldNx+2,OldNy+2);                 %Preallocate Grid Error for each node in previous U
    
        for i=2:OldNx+1
            for j=2:OldNy+1
                griderror(i,j)=abs(Uold1(i,j)-U1(2*i,2*j))/abs(U1(2*i,2*j));
            end
        end
        g1=max(griderror);
        grid_error = max(max(griderror));
        
        Uold1=U1;
        k=k+1;
end
    figure()
    surf(x,y,U1);
%%
% % %Part 2
% %Boundary Conditions
% ucoef2=0.0;
% U2 = zeros(Nx+2,Ny+2);
% err2=1000;
% k2=0;
% 
% U2(1,:)=U1(1,:);
% U2(Nx+2,:)=U1(Nx+2,:);
% U2(:,1)=U1(:,1);
% 
% F2=zeros(Nx+2,Ny+2);
% Uprev2 = U2;
% while err2>tol
% %for z=1:40000
%     for i=2:Ny+1
%         for j=2:Nx+1
%             U2(i,j)=(dy2*U2(i-1,j)+dy2*U2(i+1,j)+dx2*U2(i,j-1)+dx2*U2(i,j+1)-dx2dy2*F2(i,j))/(denomb);
%         end
%         U2(i,Ny+2)=(dy2*U2(i-1,Ny+2)+dy2*U2(i+1,Ny+2)+dx2*U2(i,Ny+1)+dx2*U2(i,Ny+1)-dx2dy2*F2(i,Ny+2))/(denomb);
%     end
%     err2=max(max(((abs(U2)-abs(Uprev2))./abs(U2))*100));
%     Uprev2=U2;
%     k2=k2+1;
% end
% 
% figure()
% surf(x,y,U2);

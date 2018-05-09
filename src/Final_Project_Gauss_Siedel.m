%Final Project Helmhotz 2-D:
%d^2u/dx^2+d^2u/dy^2+Lambda*u=F(x,y)
%Over domain of ax<x<bx and ay<y<by
%Using Gauss-Seidel

clc
clear;

%Given variables
lamda1=0.5;             %Given
ax=0;                   %Given
ay=0;                   %Given
Pi=4*atan(1);           %Given
bx=2*Pi;                %Given
by=2*Pi;                %Given
v=0;                    %Given (du/dy @y=by = 0)

%Initial Problem Setup
Nx=20;                  %Initial number of points in x
Ny=Nx;                  %Initial number of points in y
Lx=bx;                  %Length of x
Ly=by;                  %Length of y
deltax=Lx/(Nx+1);           %Step size in x
deltay=Ly/(Ny+1);           %Step size in y


%Constants
fbay = (by-ay)*(by-ay)*cos(Pi*ay/by);
gbay = ay*(by-ay)*(by-ay);
cons1 = bx-ax;
dy2=deltay*deltay;
dx2=deltax*deltax;
dx2dy2=deltax*deltax*deltay*deltay;
denominator=2*deltay*deltay+2*deltax*deltax;
denoma=denominator-dx2dy2*lamda1;
x=0:deltax:2*pi;
y=0:deltay:2*pi;

%Solution Matrix
U1 = zeros(Nx+2,Ny+2);	% Create solution matix with initial guess of ZERO

%Boundary Conditions
for j = 1:Ny+2                      %u(x=ax,y)=fb(y) boundary condition
    U1(1,j)=(by-deltay*(j-1))*(by-deltay*(j-1))*cos(Pi*deltay*(j-1)/by);
end

for j=1:Ny+2                        %u(x=bx,y)=gb(y) boundary condition
    U1(Ny+2,j)=(deltay*(j-1))*(by-deltay*(j-1))*(by-deltay*(j-1));
end

for i=1:Nx+2                        %u(x,y=ay)  boundary conditoin
    U1(i,1)=fbay+(deltax*(i-1)-ax)/cons1*(gbay-fbay);
end

% % %Wrong
% for i=1:Nx+2                        %du/dy @y=by = 0 boundary condition
%    U2(i,Ny+2)=2*deltay*c*v;
%    c=c+1;
% end

%Part 1 Using Gauss Siedel
F1=zeros(Nx+2,Ny+2);

for i=1:Nx+2
    for j=1:Ny+2
        F1(i,j)=cos(Pi/2*(2*(deltax*(i-1)-ax)/(bx-ax)+1))*sin(Pi*(deltay*(j-1)-ay)/(by-ay));
    end
end

for z=1:100
    for i=2:Nx+1
        for j=2:Ny+1
            %%%%% Different Equations for U if wanted %%%%%%%%%%%%%%%%%%%%
            U1(i,j)=(deltax*deltax*F1(i,j)-(U1(i-1,j)+U1(i+1,j)+U1(i,j-1)+U1(i,j+1)))/(-4+deltax*deltax*lamda1);
            %U1a(i,j)=(Cons4*F1(i,j)-dy2*U1(i+1,j)-dy2*U1(i-1,j)-dx2*U1(i,j+1)-dx2*U1(i,j-1))/(Cons5a);
            %U1(i,j)=((U1(i-1,j)+U1(i+1,j)+U1(i,j-1)+U1(i,j+1))-deltax*deltax*F1(i,j))/(4-deltax*deltax*lamda1);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %U1(i,j)=(dy2*U1(i-1,j)+dy2*U1(i+1,j)+dx2*U1(i,j-1)+dx2*U1(i,j+1))/(denoma)-(dx2dy2*F1(i,j))/(denoma);
        end
        %%%%%%%%%%%% Different Neumann Conditions if wanted %%%%%%%%%%%%%%
        U1(i,Ny+2)=(deltax*deltax*F1(i,Ny+2)-(U1(i-1,Ny+2)+U1(i+1,j)+U1(i,Ny+1)+U1(i,Ny+1)))/(-4+deltax*deltax*lamda1);
        %U1a(i,Ny+2)=(Cons4*F1(i,Ny+2)-dy2*U1(i+1,Ny+2)-dy2*U1(i-1,Ny+2)-dx2*U1(i,Ny+1)-dx2*U1(i,Ny+1))/(Cons5a);
        %U1(i,Ny+2)=((U1(i-1,Ny+2)+U1(i+1,Ny+2)+U1(i,Ny+1)+U1(i,Ny+1))-deltax*deltax*F1(i,Ny+2))/(4-deltax*deltax*lamda1);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %U1(i,Ny+2)=((dy2)*U1(i-1,Ny+2)+(dy2)*U1(i+1,Ny+2)+(dx2)*U1(i,Ny+1)+(dx2)*U1(i,Ny+1))/(denoma)-(dx2dy2*F1(i,Ny+2))/(denoma);
    end
    
end

figure()
surf(x,y,U1);

% %%
% % %Part 2
% %Boundary Conditions
% ucoef2=0.0;
% denomb=denominator-dx2dy2*ucoef2;
% U2 = zeros(Nx+2,Ny+2);
% 
% U2(1,:)=U1(1,:);
% U2(Nx+2,:)=U1(Nx+2,:);
% U2(:,1)=U1(:,1);
% 
% F2=zeros(Nx+2,Ny+2);
% 
% for z=1:1000
%     for i=2:Ny+1
%         for j=2:Nx+1
%             U2(i,j)=(dy2*U2(i-1,j)+dy2*U2(i+1,j)+dx2*U2(i,j-1)+dx2*U2(i,j+1)-dx2dy2*F2(i,j))/(denomb);
%         end
%         U2(i,Ny+2)=(dy2*U2(i-1,Ny+2)+dy2*U2(i+1,Ny+2)+dx2*U2(i,Ny+1)+dx2*U2(i,Ny+1)-dx2dy2*F2(i,Ny+2))/(denomb);
%     end
% end
% 
% figure()
% surf(x,y,U2);

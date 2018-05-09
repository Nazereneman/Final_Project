clc
clear;

if exist('checkpt_SOR_Manuf.mat','file')     %If a checkpoint file exists, open it
    load('checkpt_SOR_Manuf.mat')
end

%Given variables
%lamda1=0.5;            %Given for Helmholtz, not needed for Poisson
ax=-pi;                   %Given lower x bound
ay=-pi;                   %Given lower x bound
Pi=4*atan(1);           %Known value
bx=pi;                %Given upper x bound
by=pi;                %Given upper y bound
v=0;                    %Given (du/dy @y=by = 0)

%% Initial Problem Setup Nodes and Constants
%Nodes
Nx=320;                         %Initial number of points along x-axis
Ny=Nx;                          %Initial number of points along y-axis
Lx=bx-ax;                       %Length of x-axis
Ly=by-ay;                       %Length of y-axis
deltax=Lx/(Nx+1);               %Step size in x
deltay=Ly/(Ny+1);               %Step size in y
k=2;                            %Chosen constant value for Manufactured Solution
h=2;                            %Chosen constant value for Manufactured Solution
w=1.9;    %Relaxation constant
w_1=(1-w);

%Constants that will be used inside the loop
%Put here to make code run faster
fbay = (by-ay)*(by-ay)*cos(Pi*ay/by);
gbay = ay*(by-ay)*(by-ay);
cons1 = bx-ax;
dy2=deltay*deltay;
dx2=deltax*deltax;
dx2dy2=deltax*deltax*deltay*deltay;
denomin=-2*deltay*deltay-2*deltax*deltax;
k2=k*k;
h2=h*h;
kh2=-k2-h2;
k_ax=k*ax;
k_bx=k*bx;
h_ay=h*ay;
h_by=h*by;
k_dx=k*deltax;
h_dy=h*deltay;

%% Solution Matrix
U = zeros(Nx+2,Ny+2);       %Preallocate solution matix with initial guess of ZERO
Uprev = zeros(Nx+2,Ny+2);   %Preallocate dummy matrix

%Boundary Conditions
for j = 1:Ny+2                      
    U(1,j)=cos(k_ax)*cos(h_dy*(j-1));       %u(x=ax,y)=fb(y) boundary condition
    U(Nx+2,j)=cos(k_bx)*cos(h_dy*(j-1));    %u(x=bx,y)=gb(y) boundary condition
    U(j,1)=cos(k_dx*(j-1))*cos(h_ay);       %u(x,y=ay)  boundary condition
    U(j,Ny+2)=cos(k_dx*(j-1))*cos(h_by);    %u(x,y=by)  boundary condition
end

%Part 1 Using Gauss Siedel
F1=zeros(Nx+2,Ny+2);    %Preallocate F matrix outside since it doesn't change based on number of iterations

%Create F matrix
for i=1:Nx+2
    for j=1:Ny+2
        F1(i,j)=kh2*cos(k_dx*(i-1))*cos(h_dy*(j-1));
    end
end
%% Gauss-Seidel Loop

%Error Initialization
err1=1;                          %Initialize error as 1000 (# greater than tol)
iter1=0;                            %Set initial iteration to 0
tol=10^-8;                          %Define tolerance

while err1>tol
%for z=1:10000
    for j=2:Ny+1                    %loop for all interior x nodes
        for i=2:Nx+1                %loop for all interior y nodes
            %Solving for U
            U(i,j)=(w_1)*U(i,j)+w*((dx2dy2*F1(i,j)-dy2*U(i-1,j)-dy2*U(i+1,j)-dx2*U(i,j-1)-dx2*U(i,j+1))/(denomin));
        end
    end
    %er=abs(U1-Uexact)./abs(Uexact);        %Useful if trying to see where error is
    %err=max(er);
    err1=max(max(abs(Uprev-U)./abs(Uprev)));            %Find overall maximum error
    iter1=iter1+1;                                      %Increase the iteration
    Uprev=U;                                            %Seeing when to stop iterations
    if mod(iter1,100)==0                               %Save checkpoint file every 1500 iterations
        save('checkpt_SOR_Manuf.mat');                   %Saving the file
    end                                                 %Ending if loop
end                                                     %Ending while loop
%% Plotting

x=-pi:deltax:pi;    %Discretize the x axis
y=-pi:deltay:pi;    %Discretize the y axis
Ut=transpose(U);    %Transpose U so that x and y axes are on correct sides

%Exact U
Uexact=zeros(Nx+2,Ny+2);            %Preallocate Uexact
for i=1:Nx+2                        %Set Uexact
    for j=1:Ny+2
        Uexact(i,j)=cos(k*deltax*(i-1))*cos(h*deltay*(j-1));
    end
end

%Set Uref (if wanted for relative error)
Uref=1/(Nx*Ny)*sum(sum(abs(Uexact)));

%Overall error matrix
Error=zeros(Nx+2,Ny+2);
for j=1:Ny+2
    for i =1:Nx+2
        Error(i,j)=abs(U(i,j)-Uexact(i,j));
    end
end

Totalerr=sum(sum(Error));
L1err=Totalerr/(Nx*Ny);                     %L1error
L1errrel=L1err/Uref;                        %Relative L1 error
L2err=sqrt(sum(sum((Error).^2))/(Nx*Ny));   %L2error
L2errrel=L2err/Uref;                        %Relative L2 error
logL2=log10(L2err);
logdeltax=log10(deltax);

figure()
surf(x,y,U);

figure()
surf(x,y,Uexact);

figure()
plot(logL2,logdeltax);

figure()
h=surf(x,y,Error.^2);                       %Create surface plot
ylabel('y');                                %Label x-axis
xlabel('x');                                %Label y-axis
set(h,'linestyle','none');                  %Removing gridlines
%%
delete('checkpt_SOR_Manuf.mat');             %Delete checkpoint file once complete
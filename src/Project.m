%Housekeeping
clc 
clear 
tic
% Domain Initialization
% Domain: -pi<X<pi   -pi<y<pi
Ax = -pi;
Ay = -pi;
Bx = pi;
By = pi;

% Starting by Discretizing FEB 5th PAGE 2 BOT
% L: number of points, h: interval, N: number of points
Lx = Bx-Ax
Ly = By-Ay
Nx = 30;
Ny = 30;
Bt=10;
hx = Lx/(Nx-1)
hy = Ly/(Ny-1)



%Discretly 

x = Ax:hx:Bx
y = Ay:hy:By



%Fuck Yeah Boundary Conditions
% LB: Left Boundary
GLB = ((Bx-Ax)^2)*cos(((pi*Ax)/Bx));
FLB = Ax*((Bx-Ax)^2);
ULB = GLB + ((y-Ay)./(By-Ay)).*(FLB-GLB);


% TB: Top Boundary
FTB = x.*((Bx-x).^2);
UTB = FTB;

% BB: Bottom Boundary
GBB = ((Bx-x).^2).*cos(((pi.*x)./Bx));
UBB = GBB;

ULBT = ULB'

%Let make a big ass matrix
 U = [UTB ;ULBT(2:Ny-1), zeros(Ny-2,Nx-1); UBB];
 
 
 % Assuming 
 D = 1;
 
 % Von Neumann Stability Method
%  ht = (hx^2)/(D*4)
ht = ((hx^2)*(hy^2))/(2*D*((hx^2)+(hy^2)))
 Nt=round(Bt/ht)-1;
 
 
% It Begins Explicit 

v=0

UnE = zeros(Ny,Nx,Nt);
UnE(:,:,1) = U;

for k = 0:ht:Bt
   v=v+1 
for i = 2:Ny-1
    for j = 2:Nx
 
        if j == Nx
            UnE(i,j,v+1) = (UnE(i+1,j,v)-2*UnE(i,j,v)+UnE(i-1,j,v))*((ht*D)/(hx^2)) + (-2*UnE(i,j,v)+2*UnE(i,j-1,v))*((ht*D)/(hy^2))+ UnE(i,j,v);   
     
        else
        UnE(i,j,v+1) = (UnE(i+1,j,v)-2*UnE(i,j,v)+UnE(i-1,j,v))*((ht*D)/(hx^2)) + (UnE(i,j+1,v)-2*UnE(i,j,v)+UnE(i,j-1,v))*((ht*D)/(hy^2))+ UnE(i,j,v);

        end
    end
end
%Adding BC
UnE(1,1:Nx,v+1) = UTB;
UnE(Ny,1:Nx,v+1) = UBB;
UnE(:,1,v+1) = ULB;

        h= surf(x,y,UnE(:,:,v+1));
        set(h,'edgecolor','none')
            drawnow;
            refreshdata(h)
end

 UnE;
toc


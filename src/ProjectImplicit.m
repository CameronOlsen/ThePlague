%Housekeeping
clc 
clear 
close
tic
% Domain Initialization
% Domain: -pi<X<pi   -pi<y<pi
Ax = -pi;
Ay = -pi;
Bx = pi;
By = pi;

% Starting by Discretizing FEB 5th PAGE 2 BOT
% L: number of points, h: interval, N: number of points
Lx = Bx-Ax;
Ly = By-Ay;
Nx = 7;
Ny = 7;
hx = Lx/(Nx-1);
hy = Ly/(Ny-1);
ht = 1;


%Discretly 
x = Ax:hx:Bx;
y = Ay:hy:By;




%Fuck Yeah Boundary Conditions
% LB: Left Boundary
GLB = ((Bx-Ax)^2)*cos(((pi*Ax)/Bx));
FLB = Ax*((Bx-Ax)^2);
ULB = GLB + ((y-Ay)./(By-Ay)).*(FLB-GLB);

% RB
URB = 0;

% TB: Top Boundary
FTB = x.*((Bx-x).^2);
UTB = FTB;

% BB: Bottom Boundary
GBB = ((Bx-x).^2).*cos(((pi.*x)./Bx));
UBB = GBB;
ULBT = ULB';


Nt=20;

%Let make a big ass matrix
 U = [UTB ;ULBT(2:Ny-1), zeros(Ny-2,Nx-1); UBB];
 
 
% Implicit Method

Nxy= Nx*Ny;
ht = .01;
D = .1;


Lambx = (D*ht)/(hx^2);
Lamby = (D*ht)/(hy^2);


%Matrix of Coefficeints 
A = zeros(Nxy,Nxy);
for i = 2:Nxy
    for j = 2:Nxy-1
       A(i,i) = -Lamby;
       A(i,i+(Ny-1)) = -Lambx;
       A(i,i+Ny) = 2*Lambx+2*Lamby+1;
       A(i,i+(Ny+1)) = -Lambx;
       A(i,i+2*Ny) = -Lamby ;
    end
end

A=A(1:Nxy,1:Nxy);




ULBs = ULB(2:Nx-1);
UTBs = UTB(2:Nx-1);
UBBs = UBB(2:Nx-1);





%Initial U
UNEW = U;




for k = 1:Nt       
   p = reshape(U,[Nxy,1]);
   q = A\p;
   UNEW = reshape(q,[Ny,Nx]);
   %BC
   UNEW(1,:) = UTB;
   UNEW(Ny,:) = UBB;
   UNEW(:,1) = ULB;
   UNEW(Nx,:) = UNEW(Nx-1,:);


    
    
            h= surf(x,y,UNEW);
        set(h,'edgecolor','none');
            drawnow;
            refreshdata(h);
            
            
             colormap gray
                 title('Implicit Scheme','fontsize',40);
                 xlabel('X','fontsize',50) ;
                 ylabel('Y','fontsize',50) ;
                 zlabel('U','fontsize',50);              
                 colorbar
            
            U = UNEW;
            
            
            
    end
toc


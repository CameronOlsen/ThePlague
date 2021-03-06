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
Nx = 10;
Ny = 10;
Bt=20;
hx = Lx/(Nx-1);
hy = Ly/(Ny-1);


%Checkpointing

check='check.mat';


%Discretly 

x = Ax:hx:Bx;
y = Ay:hy:By;



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

ULBT = ULB';

%Let make a big ass matrix
 U = [UTB ;ULBT(2:Ny-1), zeros(Ny-2,Nx-1); UBB];
 
 
 % Assuming 
 D = 1;
 
 % Von Neumann Stability Method
%  ht = (hx^2)/(D*4)
ht = ((hx^2)*(hy^2))/(2*D*((hx^2)+(hy^2)));
 Nt=round(Bt/ht)-1;
 

 
 
 
 
% It Begins Explicit 
v=0;
c=1;
UnE = zeros(Ny,Nx,Nt);
UnE(:,:,1) = U;
Desired_Error = 10^-70;
Error=1;
    %Initializing Time Step
    for k = 0:ht:Bt
       v=v+1;
       %Space Step X
        for j = 2:Nx
        %Space Step Y
            for i = 2:Ny-1

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
                
                colormap gray
                 title('Explicit Scheme','fontsize',40)
                 xlabel('X','fontsize',50) 
                 ylabel('Y','fontsize',50) 
                 zlabel('U','fontsize',50)              
                 colorbar
                 
                if max(max(abs(UnE(:,:,v+1)-UnE(:,:,v))))< Desired_Error
                    Error = max(max(abs(UnE(:,:,v+1)-UnE(:,:,v))));
                    break
                    
                end
                if mod(v,100) == 0
                   save(check);
                end

  UnE(round(Ny/2),round(Nx/2),v+1);
    end

 UnE;
toc


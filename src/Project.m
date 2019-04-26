%Housekeeping
clc 
clear 

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
Nx = 50;
Ny = 50;
hx = Lx/(Nx-1)
hy = Ly/(Ny-1)
ht = .001

%Hello Ghosts Nodes
% Startx = Ax-hx
% Starty = Ay-hy
% Endx = Bx+hx;
% Endy = By+hy;

Startx = Ax;
Starty = Ay;
Endx = Bx;
Endy = By;



%Discretly 
m=1
for n = Startx:hx:Endx
    
    x(m) =n;
    m = m + 1;

end
m=1
for n = Starty:hy:Endy
    
    y(m) = n;
    m = m + 1;

end



%Fuck Yeah Boundary Conditions
% LB: Left Boundary
GLB = ((Bx-Ax)^2)*cos(((pi*Ax)/Bx))
FLB = Ax*((Bx-Ax)^2)
ULB = GLB + ((y-Ay)./(By-Ay)).*(FLB-GLB)

% Right Boundary Condition U = Constant
% Using the Upper or Lower Boundary Conidtions I can solve for the right side
% RB: Right Boundary
GRB = ((Bx-Bx)^2)*cos(((pi*Bx)/Bx))
FRB = Bx*((Bx-Bx)^2)

% Using Bottom BC U(x,ay) = GB(X) on bottom right corner but the whole right
% % side is constant so the corner is the same as the rest
% URB = GRB
%                UNSURE
% % Verification Using Top BC U(x,by) = FB(X)
% URB_ver = FRB

% Verified U(bx,y) = 0

% TB: Top Boundary
FTB = x.*((Bx-x).^2)
UTB = FTB

% BB: Bottom Boundary
GBB = ((Bx-x).^2).*cos(((pi.*x)./Bx));
UBB = GBB;

ULBT = ULB'

%Let make a big ass matrix
 U = [UTB ;ULBT(2:Ny-1), zeros(Ny-2,Nx-1); UBB]
 
 
% It Begins 

syms D
Un = U
for i = 1:Nx
    for j = Ny
 Un = (U(i+1,j)-2*U(i,j)+U(i-1,j))*((ht*D)/(hx^2)) + (U(i,j+1)-2*U(i,j)+U(i,j-1))*((ht*D)/(hy^2))+ U(i,j)

    end
end
% lamb = 

% for n = Ax:hx:Lx
%    
%     U(n) = ((n-h)-2n+(n+h))/(h^2) 
%     
%     
%     
%     
% end


%The Goal
% U = a + int(f(T)dt 0->t APRIL 17 PAGE 1 TOP
% For this problem a =  
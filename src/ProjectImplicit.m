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
Nx = 10;
Ny = 10;
hx = Lx/(Nx-1)
hy = Ly/(Ny-1)
ht = 1


%Discretly 
m=1
for n = Ax:hx:Bx
    
    x(m) =n;
    m = m + 1;

end
m=1
for n = Ay:hy:By
    
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
 
 
 % Setup only temporary
 D = .1
 Nt=4
 
 
% It Begins Implicit 

v=0

UnE = zeros(Ny,Nx,Nt)
UnE(:,:,1) = U


v=1
Nt=12
UnI = zeros(Ny,Nx,Nt)
UnI(:,:,1) = U
UnI(:,:,2) = U

for k = 0:ht:10
   v=v+1 
for i = 2:Ny-1
    for j = 2:Nx
        if j == Nx
         UnI(i,j,v) = (UnI(i+1,j,v)-2*UnI(i,j,v)+UnI(i-1,j,v))*((ht*D)/(hx^2)) + (-2*UnI(i,j,v)+2*UnI(i,j-1,v))*((ht*D)/(hy^2))+ UnI(i,j,v-1);    
        else
        UnI(i,j,v) = (UnI(i+1,j,v)-2*UnI(i,j,v)+UnI(i-1,j,v))*((ht*D)/(hx^2)) + (UnI(i,j+1,v)-2*UnI(i,j,v)+UnI(i,j-1,v))*((ht*D)/(hy^2))+ UnI(i,j,v-1);
        end
    end
end
%Adding BC
% Un(:,:,v+1) = [Un, zeros(Ny-1,1);zeros(1,Nx)]
UnE(1,1:Nx,v+1) = UTB;
UnE(Ny,1:Nx,v+1) = UBB;

UnE(:,1,v+1) = ULB;
% UnE(:,Nx+1,v+1) = UnE(:,Nx-1,v+1)


end
 
UnE
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
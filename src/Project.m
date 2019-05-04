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
Nx = 20;
Ny = 20;
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
GLB = ((Bx-Ax)^2)*cos(((pi*Ax)/Bx));
FLB = Ax*((Bx-Ax)^2);
ULB = GLB + ((y-Ay)./(By-Ay)).*(FLB-GLB);

% Right Boundary Condition U = Constant
% Using the Upper or Lower Boundary Conidtions I can solve for the right side
% RB: Right Boundary
GRB = ((Bx-Bx)^2)*cos(((pi*Bx)/Bx));
FRB = Bx*((Bx-Bx)^2);

% Using Bottom BC U(x,ay) = GB(X) on bottom right corner but the whole right
% % side is constant so the corner is the same as the rest
% URB = GRB
%                UNSURE
% % Verification Using Top BC U(x,by) = FB(X)
% URB_ver = FRB

% Verified U(bx,y) = 0

% TB: Top Boundary
FTB = x.*((Bx-x).^2);
UTB = FTB;

% BB: Bottom Boundary
GBB = ((Bx-x).^2).*cos(((pi.*x)./Bx));
UBB = GBB;

ULBT = ULB'

%Let make a big ass matrix
 U = [UTB ;ULBT(2:Ny-1), zeros(Ny-2,Nx-1); UBB];
 
 
 % Setup only temporary
 D = .01;
 Nt=300;

 
 
% It Begins Explicit 

v=0

UnE = zeros(Ny,Nx,Nt);
UnE(:,:,1) = U;
 % Lambda 1 and 2 must be less than 1/2 for this to be stable
 if D*ht/(hx*hx)<=.5
     if D*ht/(hx*hx)<=.5
for k = 0:ht:300
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

        h= surf(x,y,UnE(:,:,k+1));
            drawnow;
            refreshdata(h)
end

 UnE
     else
         fprintf("Will not Converge")
     end
 else
      fprintf("Will not Converge")
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
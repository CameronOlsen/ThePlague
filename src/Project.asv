clc 
clear
% Domain Initialization
% Domain: -pi<X<pi   -pi<y<pi
Sx = -pi;
Sy = -pi;
Ex = pi;
Ey = pi;
% Starting by Discretizing FEB 5th PAGE 2 BOT
% L: number of points, h: interval, N: number of points
Lx = Ex-Sx
Ly = Ey-Sy
Nx = 50;
Ny = 50;
hx = Lx/(Nx-1)
hy = Ly/(Ny-1)


%Hello Ghosts Nodes
Startx = Sx-hx
Starty = Sy-hy
Endx = Ex+hx;
Endy = Ey+hy;

%Discretly 
m=1
for n = Startx:hx:Endx
    
    Ux(m) =n;
    m = m + 1;

end
m=1
for n = Starty:hy:Endy
    
    Uy(m) = n;
    m = m + 1;

end


% It Begins 

% lamb = 

% for n = Sx:hx:Lx
%    
%     U(n) = ((n-h)-2n+(n+h))/(h^2) 
%     
%     
%     
%     
% end
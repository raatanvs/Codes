% 2D Steady State Heat Conduction

clear
clc
tic

nodesCoordinates = ([[0,0];[0,1];[0,2];[1,0];[1,1/2];[1,3/2];[1,2];[2,0];[2,1];[2,2];[3,0];[3,1/2];[3,3/2];[3,2]]);


elementNodes = ([[1,2,5];[2,3,6];[6,5,2];[4,1,5];[6,3,7];[10,6,7];[6,10,9];[6,9,5];[5,9,8];[8,4,5];[12,11,8];[8,9,12];[12,9,13];[9,10,13];[10,14,13]]);

%% Looping through to find out the local and global stiffness matrices
for ii=1:length(elementNodes)  
    nodeList = elementNodes(ii,:); 
    % Gets the local stiffness matrix in 3 x 3 form
    localK(:,:,ii)  = localStiffness(nodesCoordinates(nodeList(1),:),nodesCoordinates(nodeList(2),:),nodesCoordinates(nodeList(3),:));
    
    % Gets the stiffness matrix in 14 x 14 form for superposition
    globalKall(:,:,ii) = globalStiffness(localK(:,:,ii),nodeList);
      
end
% Getting the global stiffness matrices by superposition of nodes
% across elements
globalK = sum(globalKall,3);   %adding across third dimension

%applying boundary conditions to find the load vectors at corner nodes
kUse=globalK(4:10,4:10);
k1Use = globalK(4:10,1:3);     % Matrix Formed by T1,T2,T3
k2Use = globalK(4:10,11:14);   % Matrix Formed by T11 - T14

% Conversion of Temperature from Celcius to Kelvin
tempUse = convtemp([200,200,200,25,25,25,25],'C','K')';

Mat1 = zeros(7,1)-k1Use*tempUse(1:3)-k2Use*tempUse(4:end);

TempRest = linsolve(kUse,Mat1);
toc

temp = [tempUse(1:3);TempRest;tempUse(4:end)]



%% Formulating local stiffness matrix given coordinates of nodes

function[K] =  localStiffness(co1,co2,co3)
xa=co1(1);
ya=co1(2);

xb=co2(1);
yb=co2(2);

xc=co3(1);
yc=co3(2);

Ae_matrix = [1,xa,ya;1,xb,yb;1,xc,yc];

Ae = abs(0.5*det(Ae_matrix));

a1 = xb*yc - xc*yb;
b1= yb-yc;
c1 = xc-xb;

a2 =  xc*ya-xa*yc;
b2=yc-ya;
c2 = xa-xc;

a3 = xa*yb -xb*ya;
b3 = ya-yb;
c3 = xb-xa;

kxx=45;
kyy=45;
K = 1/(4*Ae)* [kxx*b1^2+kyy*c1^2  ,kxx*b1*b2+kyy*c1*c2,kxx*b1*b3+kyy*c1*c3;...
               kxx*b1*b2+kyy*c1*c2,kxx*b2^2+kyy*c2^2  ,kxx*b2*b3+kyy*c2*c3;...
               kxx*b1*b3+kyy*c1*c3,kxx*b2*b3+kyy*c2*c3,kxx*b3^2+kyy*c3^2];      

end
%% Assembly of Global Stiffness matrix for each element
function [globalK] = globalStiffness(localStiffness,nodesArray)
   
globalK = zeros(14,14);

[x,y]=meshgrid(nodesArray);
X=x(:);
Y=y(:);
localK = localStiffness(:);

for ii = 1:length(localK)
globalK(X(ii),Y(ii)) = localK(ii);
end
  
end




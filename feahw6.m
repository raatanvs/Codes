%% 2 Dimensional Finite Element Structural Solver for the custom mesh specified
clear all
clc

% Indicating the number of nodes
Nodes = 1:10;   

% Indicating the coordinates of nodes
Coordinates =[0,0;
    3.33,0;
    6.66,0;
    10,0;
    0,2.5;
    3.33,2.5;
    6.66,2.5;
    0,5;
    3.33,5;
    10,5];

% Element Mapping
Elements = [6 5 1 2;
    7 6 2 3;
    10 7 3 4;
    9 8 5 6;
    10 9 6 7];



%% Constructing Shape Functions
syms zeta;
syms neta;

shapeFunction(1) = (1+zeta)*(1+neta)/4;
shapeFunction(2) = (1-zeta)*(1+neta)/4;
shapeFunction(3) = (1-zeta)*(1-neta)/4;
shapeFunction(4) = (1+zeta)*(1-neta)/4;

for jj=1:5

elemCoordinates = Coordinates(Elements(jj,:)',:);

% Performing Geometric Transformation
syms x;
syms y;

x= shapeFunction*elemCoordinates(:,1);
y= shapeFunction*elemCoordinates(:,2);

J(:,:,jj) = [diff(x,zeta),diff(y,zeta);
    diff(x,neta),diff(y,neta)];

B1 = inv(J(:,:,jj)) * [diff(shapeFunction,zeta);diff(shapeFunction,neta)]


trial = [B1(1,:);zeros(1,numel(shapeFunction))];
B(1,:,jj) = reshape(trial,[1,8]);
clear trial
trial = [zeros(1,numel(shapeFunction));B1(2,:)];
B(2,:,jj) = reshape(trial,[1,8]);
clear trial
trial = [B1(2,:);B1(1,:)];
B(3,:,jj) = reshape(trial,[1,8]);



end
shapeDiff=B;




%% Applying Plane Strain Condition
e = 1e3;    % Young's Modulus
g = 0.3   % Poisson ratio of the material


E  = e/((1+g)*(1-2*g))*[1-g g 0;g 1-g 0;0 0 0.5*(1-2*g)];


%% Gauss Quadrature Integration, since weights are 1 for 2x2 rule, they are added straight away 
k=zeros(8,8,5);
for kk=1:5
I(:,:,kk) = B(:,:,kk)'*E*B(:,:,kk)*det(J(:,:,kk));
a = 1/sqrt(3);
k(:,:,kk)=double(k(:,:,kk)+subs(I(:,:,kk),[zeta,neta],[a,a])+subs(I(:,:,kk),[zeta,neta],[a,-a])+subs(I(:,:,kk),[zeta,neta],[-a,-a])+subs(I(:,:,kk),[zeta,neta],[-a,a]));
sanityCheck(kk) = issymmetric(k(:,:,kk));
end
clear B;
%% Expanding local to Global Size

K_Global = zeros(20,20);

for ii = 1:5
A = [Elements(ii,:)];
I_exact = k(:,:,ii);

    for m = 1:8
        B = [2*A(1,1)-1 2*A(1,1) 2*A(1,2)-1 2*A(1,2) 2*A(1,3)-1 2*A(1,3) 2*A(1,4)-1 2*A(1,4)];
        if rem(m,2) == 1
            pos = A(1,(m-1)/2 + 1);
            for l = 1: 8
                K_Global(B(l),2*pos-1) = K_Global(B(l),2*pos-1) + I_exact(l,m);
            end
        elseif rem(m,2) == 0
            pos = A(1,(m)/2);
            for h = 1:8
                K_Global(B(h),2*pos) = K_Global(B(h),2*pos) + I_exact(h,m);
            end
        end
    end
end
globalStiffness = K_Global;
R = zeros(10,2);   % 10 rows each nodes and 2 columns for dofs

R(8,2) = 33.33/2;
R(9,2) = (33.33+66.66)/2;
R(10,2)= 66.66/2;
R(10,1)= 4.2;
R(4,1) = 1.8;

R = R.*1e3;
R = reshape(R',[20,1]);
u=ones(1,10);
v=ones(1,10);
v(1:4)=0;
u(1)=0;
u(7)=0;
u(8)=0;

d = [u;v];

D = reshape(d,[20,1]);

for kk=1:numel(D)
if (D(kk)==0)
    disp(['Yes Boss',num2str(kk)])
globalStiffness(kk,:) = 0;
globalStiffness(kk,kk)=1;
end

end

nodalDisplacements = reshape(linsolve(globalStiffness,R),[2,10])

% stress Calculation - using plane strain condition

%for element A - Element 3 

nodesA = Elements(3,:);
loadVecA = reshape(nodalDisplacements(:,nodesA),[8,1]);

strain = shapeDiff(:,:,3)*loadVecA;

stress = E*strain;

% Finding the Stress Values at the gauss points

disp(['Gauss Points: (',num2str(a),num2str(-a),')'])
S_GP(:,1) = subs(stress,[zeta,neta],[a,a]);
S_GP(:,2) = subs(stress,[zeta,neta],[-a,a]);
S_GP(:,3) = subs(stress,[zeta,neta],[-a,-a]);
S_GP(:,4) = subs(stress,[zeta,neta],[a,-a]);

SGP = double(S_GP)









clear; %clears workspace data and variables
clc; %clears output in command
close all; %close all figures and windows related to program

%Translate 3D trilateration to 2D trilateration using known z-coordinate
%reading. Z reading is from barometer, with expected error of 0.10m/10cm.

%values are in metres.

% generate the position and distances of nodes
% rows correspond to anchors, columns coorespond to dimension/axes x,y,z
A1 = [0; 0; 0]; %For simplicity, start as close to origin as possible
A2 = [10; 0; 0]; %extreme end for positive x and y axis, use for z-axis check
A3 = [0; 10; 0]; % positive x, positive y
A4 = [10; 10; 0]; %positive x, negative y
A5 = [-5; 5; 0]; % negative x, positive y
A6 = [-10; -10; 0]; %extreme end for negative x and y axis, use for z-axis check
% Anchor = [A1 A2 A3 A4]; %Unique 2D solution
Anchor = [A1 A2 A3 A4 A5 A6]; %Unique 3D solution 

%n--number of nodes
n = length(Anchor);
% n = 4;

% Area_cap = 10; %Set maximum coordinate aka size of floor
% %Generate 2D random coordinates, then add in a fixed z-coordinate. Assume
% %all anchors are on the same level.
% Anchor = Area_cap*randn(2,n); 
% %Impose ground level assumption for all anchors first. Height variation of
% %anchors can come later
% for i = 1:n
%     Anchor(3,n) = 0;
% end
% %randn(n,3) generates a nx3 matrix of normally distributed random values
% %from 0 to 1

%Define degree of random measurement error from ranging sensor here
random_error = 0.2;
baro_error = 0.1; %barometer has error up to 10 cm i.e. 0.1 m

%define height reading from barometer, h
%typical storey height is 3m.
h = 18; %- baro_error + 2*baro_error*rand;

% target position, let this be known as tag.
Tag = [3;4;h];

di = zeros(1,n); %d(i) is true distance between anchor and tag; i = 1 to n
Radius = zeros(1,n); %Radius is matrix of measured distance estimate between anchor and target coordinate
for i=1:n
     di(i) = sqrt(((Tag(1) - Anchor(1,i))^2) + ...
                  ((Tag(2) - Anchor(2,i))^2) + ...
                  ((Tag(3) - Anchor(3,i))^2));
     %If Radius is ideally same as measured distance di, then we can add
     %random error to di to create a realistic Radius.
%      Radius(i) = di(i) - random_error + 2*random_error*rand;  %with random error
     Radius(i) = di(i); %no error, base algorithm check 
end

%d(i)(r) is true distance between anchor and reference anchor; i = 2 to n, r is reference
d_ir = zeros(1,n-1); 
for i = 1:n-1
    d_ir(i) = sqrt(((Anchor(1,1) - Anchor(1,i+1))^2) + ...
                   ((Anchor(2,1) - Anchor(2,i+1))^2) + ...
                   ((Anchor(3,1) - Anchor(3,i+1))^2));
end

%Set up Linear Least Squares model
% Ax = b

%In traditional LSQ evaluation (A' is tranpose in MATLAB terms)
% If A'A is non-singular and well-conditioned, A'Ax = A'b
% then x = inv(A'*A)*A'*b;
% If poorly conditioned, perform Orthogonal-triangular decomposition qr on
% A, A = QR, Rx = Q'b, then x = inv(R)*Q'*b
%For streamlining of code, native MATLAB function lsqr is used here.

%A is state matrix of distances between anchor and reference anchor
%A = [transpose(A2-A1); ...; transpose(An-A1)];
A = zeros(n-1,3);
for i = 2:n %first element corresponds to 2nd anchor A2, not reference anchor A1
    for k = 1:3 
    %k is number of axes/dimensions
    %i is coordinate in axis corresponding to the anchor number
    %A in this case is transpose of (Ai - A1)
    A(i-1,k) = Anchor(k,i) - Anchor(k,1);
    end
end

b = zeros(n-1,1); %b is matrix of measured distance estimate from target coordinate
for i = 2:n %first element corresponds to 2nd anchor A2, not reference anchor A1
    %obtain approximation by cosine rule
    %bi1 = 0.5*((r1)^2 - (ri)^2 + (di1)^2)
    b(i-1,1) = 0.5*(Radius(1)^2 - Radius(i)^2 + d_ir(i-1)^2);
end

%Implement MATLAB "lqsr" function to find LSQ coordinates (without
%reference tag compensation)
%Sum of square error for LSQ has to be found separately because of this
%compensation
maxit = 100; %maximum iterations
tol = 10^-4; %tolerance, setting accuracy up to mm.
x_raw = lsqr(A,b,tol,maxit);

%Actual coordinates of tag by LSQ method
%Because x in this case is really [x-x(A1); y-y(A1); z-z(A1)], need to add
%back A1's coordinates to find actual tag coordinates by LSQ method
x_lsq = x_raw + Anchor(:,1);
%LSQ end

%Let there be horizontal projection from ranging called p
p = zeros(1,n);
%Let angle of range vector/measurement from horizontal porjection be theta
theta = zeros(1,n); %Actual angle from anchor to tag
meas_theta = zeros(1,n); %measured angle from anchor to tag, different from actual theta due to error
%Let there be 2D trilateration measurement r comprising of x and y
%coordinates only.
true_r = zeros(1,n);
for i=1:n
    %As range measurement d is known, height h is known, p can be found by
    %Pythogoras theorem. Assume height is zeroed to ground level or
    %reference anchor.
    p(i) = sqrt(Radius(i)^2 - h^2);
    theta(i) = asind(h / di(i));
    meas_theta(i) = asind(h / Radius(i));
    true_r(i) = sqrt((Tag(1) - Anchor(1,i))^2 + (Tag(2) - Anchor(2,i))^2);
end

%Objective function is error between 2D trilateration distance and
%horizontal projection
% r = @(x)[sqrt((x(1) - Anchor(1,1))^2 + (x(2) - Anchor(2,1))^2) - p(1);
%          sqrt((x(1) - Anchor(1,2))^2 + (x(2) - Anchor(2,2))^2) - p(2); 
%          sqrt((x(1) - Anchor(1,3))^2 + (x(2) - Anchor(2,3))^2) - p(3); 
%          sqrt((x(1) - Anchor(1,4))^2 + (x(2) - Anchor(2,4))^2) - p(4)];

r = @(x)[sqrt((x(1) - Anchor(1,1))^2 + (x(2) - Anchor(2,1))^2) - p(1);
         sqrt((x(1) - Anchor(1,2))^2 + (x(2) - Anchor(2,2))^2) - p(2); 
         sqrt((x(1) - Anchor(1,3))^2 + (x(2) - Anchor(2,3))^2) - p(3); 
         sqrt((x(1) - Anchor(1,4))^2 + (x(2) - Anchor(2,4))^2) - p(4);
         sqrt((x(1) - Anchor(1,5))^2 + (x(2) - Anchor(2,5))^2) - p(5);
         sqrt((x(1) - Anchor(1,6))^2 + (x(2) - Anchor(2,6))^2) - p(6)];     
     
     
%Objective function is r-p = 0 i.e. the 2D range measurement is equal to
%the "scalar horizontal projection" p. Solve for system of non-linear
%equations.
%If not, use minimum sum of sqaure error; objective function is min (r-p)^2
options = optimoptions('lsqnonlin','Display','iter');
[Target,sqerror_sum] = lsqnonlin(r,x_lsq,[],[],options);
%Substitute the barometer reading into the calculated 2D coordinate
Target(3) = h;

%Find the sum of square error from this method in 3D.
final_3Derror = 0;
threed_e = zeros(1,n); %create matrix of error of each anchor
for i = 1:n
    threed_e(1,i) = sqrt((Target(1) - Anchor(1,i))^2 + ...
                         (Target(2) - Anchor(2,i))^2 + ...
                         (Target(3) - Anchor(3,i))^2) - Radius(i);
end

%Compile sum of square error for NLSQ method
final_3Dsqerror = 0;
for i = 1:n
   final_3Dsqerror = final_3Dsqerror + threed_e(i)^2;
end

%%%%%%%%%%%%%%%%%
%%%  Display  %%%
%%%%%%%%%%%%%%%%%
fprintf('LSQ estimate is (%g, %g, %g) \n',x_lsq(1),x_lsq(2),x_lsq(3));
fprintf('Calaulcated Target is (%g, %g, %g) \n',Target(1),Target(2),Target(3));
fprintf('Actual Target is (%g, %g, %g) \n',Tag(1),Tag(2),Tag(3));
fprintf('Sum of square error in 2D is %g \n', sqerror_sum);
fprintf('Sum of square error in 3D is %g \n', final_3Dsqerror);
%Horizontal distnace comparison
fprintf('Horizontal projection is'); disp(p);
fprintf('True 2D distance is'); disp(true_r);
%Angle comparison
fprintf('Angle of anchor to Target is '); disp(theta);
fprintf('Measured Angle of anchor to Target is '); disp(meas_theta);
   
%%%%%%%%%%%%%%%%%%
%%%%   plot   %%%%
%%%%%%%%%%%%%%%%%%

% %%% figure of found 2D 'target' point  %%%
% figure('Name','Trend of coordinates through iterations of NLSQ');
% hold on;
% plot(Plot_point(1,:),Plot_point(2,:),'b*');
% grid on;
% 
% %%% Target from LSQ in 2D %%%
% plot(x_lsq(1),x_lsq(2),'r.','MarkerSize',20);
% 
% %%% Final Target from NLSQ  %%%
% plot(Target(1),Target(2),'kx','MarkerSize',20,'LineWidth',2);
% title('Trend of 2D Sum of square error through iterations of NLSQ');
% xlabel('x' );
% ylabel('y');
% hold off;
% 
% %%%  plot the 2D error trend  %%%
% figure('Name','Trend of Sum of square error through iterations of NLSQ');
% plot(NLSQ_sqe_trend);
% title('Trend of Sum of square error through iterations of NLSQ');
% xlabel('iteration(n)' );
% ylabel('NLSQ(\epsilon)','Fontname', 'Times New Roman');
% 
% %%%  plot the anchor point and the target point  %%%
figure('Name','The position of set Anchor point and target point');
hold on;
% plot(Plot_point(1,:),Plot_point(2,:),'b*');
grid on;

%%% Target from LSQ  %%%
plot3(x_lsq(1),x_lsq(2),x_lsq(3),'r.','MarkerSize',20);

%%% Final Target from NLSQ  %%%
plot3(Target(1),Target(2),Target(3),'kx','MarkerSize',20,'LineWidth',2);
plot3(Anchor(1,:),Anchor(2,:),Anchor(3,:),'mo');
title('The position of set Anchor point and target point');
xlabel('x');
ylabel('y');
zlabel('z');
hold off;
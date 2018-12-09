clear; %clears workspace data and variables
clc; %clears output in command
close all; %close all figures and windows related to program

%Translate 3D trilateration to 2D trilateration using known z-coordinate
%reading. Z reading is from barometer, with expected error of 0.10m/10cm.

%values are in metres.

%Let reference anchor be A1

A1 = [0; 0; 0]; %For simplicity, start as close to origin as possible
A2 = [10; 0; 0]; %extreme end for positive x and y axis, use for z-axis check
A3 = [0; 10; 10]; % positive x, positive y
A4 = [10; 10; 0]; %positive x, negative y
A5 = [-5; 5; 0]; % negative x, positive y
A6 = [-10; -10; 0]; %extreme end for negative x and y axis, use for z-axis check
% A7 = [-3; -7; 0]; %negative x, negative y, not to extreme end
% A8 = [6; 8; 3]; %wild card, as long as not far away from extreme end

%Create matrix system of anchor coordinates.
%As a result, each row correspond to a dimension/axis i.e. x,y,z
%Each column corresponds to the specific anchor
Anchor = [A1 A2 A3 A4]; %Unique 2D solution
% Anchor = [A1 A2 A3 A4 A5 A6]; %Unique 3D solution 

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
h = 15; %- baro_error + 2*baro_error*rand;

% target position, let this be known as tag.
Tag = [8; 3; h];

%tolerance for NLSQ method
tol = 10^-12;

%%%%%%%%%%%%%%%%%%
%Range measurement
%%%%%%%%%%%%%%%%%%
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Horizontal Projection and other statistics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    p(i) = sqrt(Radius(i)^2 - (h- Anchor(3,i))^2);  %% FIXED THE BUG
    theta(i) = asind(h / di(i));
    meas_theta(i) = asind(h / Radius(i));
    true_r(i) = sqrt((Tag(1) - Anchor(1,i))^2 + (Tag(2) - Anchor(2,i))^2);
end


disp("p")
disp(p)
disp("true_r")
disp(true_r)
disp("Tag")
disp(Tag)
p = true_r
%%%%
%LSQ
%%%%

% Ax = b
% In traditional LSQ evaluation (A' is tranpose in MATLAB terms)
% If A'A is non-singular and well-conditioned, A'Ax = A'b
% then x = inv(A'*A)*A'*b;
% If poorly conditioned, perform Orthogonal-triangular decomposition qr on
% A, A = QR, Rx = Q'b, then x = inv(R)*Q'*b

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

% %Implement MATLAB "lqsr" function to find LSQ coordinates (without
% %reference tag compensation)
% %Sum of square error for LSQ has to be found separately because of this
% %compensation
% maxit = 100; %maximum iterations
% tol = 10^-4; %tolerance, setting accuracy up to mm.
% x_raw = lsqr(A,b,tol,maxit);

%In traditional LSQ evaluation (A' is transpose in MATLAB terms)
% If A'A is non-singular and well-conditioned,
% then x = inv(A' * A) * A' * b;
% If poorly conditioned, perform Orthogonal-triangular decomposition qr on
% A, A = QR, Rx = Q'b, then x = inv(R)*Q'*b

sing_cond = cond(A'*A,inf);
if sing_cond < 10 %find condition number and determine if matrix is well conditioned
    x_raw = inv(A'*A)*A'*b; %assume A'A is non-singular and well-conditioned
else
    [Q,R] = qr(A);
    %Because A is not necessarily sqaure, use Moore-Penrose inverse.
    x_raw = pinv(R)*Q'*b;
end

%Actual coordinates of tag by LSQ method
%Because x in this case is really [x-x(A1); y-y(A1); z-z(A1)], need to add
%back A1's coordinates to find actual tag coordinates by LSQ method
x_lsq = x_raw + Anchor(:,1);
%Anchor(:,1) refers to the 1st column of Anchor i.e. Anchor 1 coordinates

%Let there be error variable e(i), where e(i) is difference between 
%actual 2D distance from anchor to tag and horizontal projection.
%ei = sqrt((x-xi)^2+(y-yi)^2) - p(i), i = 1 to n
lsq_e = zeros(1,n); %create matrix of error of each anchor, 1xn matrix
for i = 1:n
    lsq_e(1,i) = sqrt((x_lsq(1) - Anchor(1,i))^2 + ...
                      (x_lsq(2) - Anchor(2,i))^2 + ...
                      (x_lsq(3) - Anchor(3,i))^2) - Radius(i);
end

%Compile sum of square error for LSQ method
LSQ_square_error = 0;
for i = 1:n
    LSQ_square_error = LSQ_square_error + lsq_e(i)^2;
end

%%%%%%%%%%%%%
%NLSQ Model
%%%%%%%%%%%%%



%Aim of NLSQ is to minimise Error function F = sum(ei^2)

%Objective function is r-p = 0 i.e. the 2D range measurement is equal to
%the "scalar horizontal projection" p. Solve for system of non-linear
%equations.
%If not, use minimum sum of sqaure error; objective function is min (r-p)^2

%As we are constraining a 3D to a 2D problem, the coordinate going through
%NLSQ is in 2 dimensions instead of 3.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Newton Ralphton loop for NLSQ
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Use Newton Ralphton method to reach minimum error function

%Get corrected coordinates changed by Newton's method. Let this be C.
T = [x_lsq(1); x_lsq(2)]; %Set LSQ result as initial guess

%This will be the final sum of square error after Newton's method, 
%initial value is the error from LSQ
NEsq = LSQ_square_error; 
%This will be final coordinates after Newton's method,
%initial value is C, the variable to be corrected iteratively
%Initialise variable to store C at iteration k for iteration k+1
previous_T = zeros(2,1);
Target = T; 

%Impose cap on iterations run.
iter = 0;
iter_cap = 10000;

%initialise variable for plotting points over NLSQ method
Plot_point = [x_lsq(1); x_lsq(2)]; %initialise with initial LSQ coordinate

disp("Anchor")
disp(Anchor)

while iter < iter_cap
    %Let there be Jacobian of tag coordinates J(q).
    %Let there be variable JqE aka transpose of J(q)*e
    xe = 0; ye = 0;
    %Let there be square matrix g = J(q)'*J(q)
    Jx = 0; Jy = 0;
    Jxy = 0; 
    
    %Let there be distance between anchor and tag d(i), i = 0 to n
    d_i = zeros(1,n);
    for i = 1:n
    %get actual 2D distance from one anchor to the target
        d_i(1,i) = sqrt((T(1) - Anchor(1,i))^2 +...
                        (T(2) - Anchor(2,i))^2);
    end
    
    for i = 1:n
        xe = xe + (T(1) - Anchor(1,i)) * (d_i(i) - p(i)) / d_i(i);
        ye = ye + (T(2) - Anchor(2,i)) * (d_i(i) - p(i)) / d_i(i);
    end

    JqE = [xe; ye];
    
    for i =1:n
        Jx = Jx + (T(1) - Anchor(1,i)) ^ 2 / (d_i(i) ^ 2);
        Jy = Jy + (T(2) - Anchor(2,i)) ^ 2 / (d_i(i) ^ 2);
        
        Jxy = Jxy + (T(1) - Anchor(1,i)) * (T(2) - Anchor(2,i)) / (d_i(i) ^ 2);
    end
    
    g = [Jx Jxy; Jxy Jy];
    
    %Get corrected coordinates
    T = T - (inv(g) * JqE); %corrected coordinates here
        
    %Update calculated 2D points in the list
    Plot_point = [Plot_point T];
    
    %Find new error
    nlsq_e = zeros(1,n); %create matrix of error of each anchor
    for i = 1:n
        nlsq_e(1,i) = sqrt((T(1) - Anchor(1,i))^2 + ...
                           (T(2) - Anchor(2,i))^2) - p(i);
    end
    
    %Compile sum of square error for NLSQ method
    NLSQ_square_error = 0;
    for i = 1:n
        NLSQ_square_error = NLSQ_square_error + nlsq_e(i)^2;
    end
    
    %NLSQ_sqe_trend is collection of NLSQ Sum of square error at each iteration over
    %time, p is data point number
    %At every iteration, paste in the NLSQ Sum of square error at iteration
    %p
    NLSQ_sqe_trend(iter+1) = NLSQ_square_error;
    
    %Update iteration
    Target = T;
    NEsq = NLSQ_square_error;
    
    %If criterion reached before iteration cap, terminate Newton method
    %Criterion: Change in coordinates between current and previous
    %iteration is sufficiently small (ideally 0)
    if abs(Target(1) - previous_T(1)) < tol && ...
       abs(Target(2) - previous_T(2)) < tol
        break
    end
    
    %update iterations done so far, store current iteration k for next
    %cycle k+1
    iter = iter + 1;
    previous_T = T; %store previous iteration of T here
end

%Substitute back in known z-coordinate into calculations
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

%%%%%%%%%%%%%%%%%%
%%% True Error %%% 
%%%$$%%%%%%%%%%%%%
%find true error of lsq, nlsq, and lsqnonlin coordinates
te_lsq = zeros(1,3);
for i=1:3
    te_lsq(i) = Tag(i) - x_lsq(i);
end

te_nlsq = zeros(1,3);
for i=1:3
    te_nlsq(i) = Tag(i) - Target(i);
end

%Matrix of actual distance between anchor and tag is found earlier as di.

%find true error from tag to each anchor (true distance vs measured distance)
tag_e = zeros(1,n); %create matrix of error of each anchor
for i = 1:n
    tag_e(1,i) = di(i) - Radius(i);
end
%Compile sum of square error for true coordinates
True_square_error = 0;
for i = 1:n
    True_square_error = True_square_error + tag_e(i)^2;
end

%%%%%%%%%%%%%%%%%
%%%  Display  %%%
%%%%%%%%%%%%%%%%%

%Display true error between tag coordinates and calculated coordinates 
% R Result
disp('Set Radius and True target point: ');
fprintf('Measured Distance between Tag and anchor is:\n'); disp(Radius);
fprintf('Actual Distance between Tag and anchor is:\n'); disp(di);
fprintf('True measurement error of each anchor is:\n'); disp(tag_e); 
fprintf('Sum of true square error is %g \n\n', True_square_error);
fprintf('True coordinates is (%g, %g, %g) \n\n\n', ...
        Tag(1), Tag(2), Tag(3));

%Horizontal distnace comparison
fprintf('Horizontal projection is'); disp(p);
fprintf('True 2D distance is'); disp(true_r);
%Angle comparison
fprintf('Angle of anchor to Target is '); disp(theta);
fprintf('Measured Angle of anchor to Target is '); disp(meas_theta);

%LSQ Result
disp('LSQ Result:');
fprintf('LSQ coordinates were (%g, %g, %g) \n\n', x_lsq(1), x_lsq(2), x_lsq(3));
disp('Error of each anchor from LSQ is:'); disp(lsq_e);
fprintf('Sum of square error from LSQ is: %g \n \n', LSQ_square_error);
fprintf('True error of LSQ coordinates is: (%g, %g, %g) \n\n', ...
        te_lsq(1), te_lsq(2), te_lsq(3));
fprintf('Sum of square true error from LSQ is: %.5e \n\n\n', sum(te_lsq.^2)); 
    
%NLSQ results
disp('NLSQ Result:');
fprintf('NLSQ method converged after iteration %d \n\n', iter);
fprintf('NLSQ coordinates were: (%g, %g, %g) \n\n', Target(1), Target(2), Target(3));
disp('Error of each anchor from NLSQ is:'); disp(nlsq_e);
fprintf('Sum of 2D square error is %g \n', NEsq);
fprintf('Sum of 3D square error is %g \n\n', final_3Dsqerror);
fprintf('True error of NLSQ coordinates is: (%g, %g, %g) \n\n', ...
        te_nlsq(1), te_nlsq(2), te_nlsq(3));
fprintf('Sum of square true error from NLSQ is: %.5e \n\n\n', sum(te_nlsq.^2)); 
    
%%%%%%%%%%%%%%%%%%
%%%%   plot   %%%%
%%%%%%%%%%%%%%%%%%

%%% figure of found 2D 'target' point  %%%
figure('Name','Trend of coordinates through iterations of NLSQ');
hold on;
plot(Plot_point(1,:),Plot_point(2,:),'b*');
grid on;

%%% Target from LSQ in 2D %%%
plot(x_lsq(1),x_lsq(2),'r.','MarkerSize',20);

%%% Final Target from NLSQ  %%%
plot(Target(1),Target(2),'kx','MarkerSize',20,'LineWidth',2);
title('Trend of 2D Sum of square error through iterations of NLSQ');
xlabel('x' );
ylabel('y');
hold off;

%%%  plot the 2D error trend  %%%
figure('Name','Trend of Sum of square error through iterations of NLSQ');
plot(NLSQ_sqe_trend);
title('Trend of Sum of square error through iterations of NLSQ');
xlabel('iteration(n)' );
ylabel('NLSQ(\epsilon)','Fontname', 'Times New Roman');

%%%  plot the anchor point and the target point  %%%
figure('Name','The position of set Anchor point and target point');
hold on;
plot(Plot_point(1,:),Plot_point(2,:),'b*');
grid on;

%%% Target from LSQ  %%%
plot3(x_lsq(1),x_lsq(2),x_lsq(3),'r.','MarkerSize',20);

%%% Final Target from NLSQ  %%%
plot3(Target(1),Target(2),Target(3),'kx','MarkerSize',20,'LineWidth',2);
plot3(Anchor(1,:),Anchor(2,:),Anchor(3,:),'mo');
title('The position of set Anchor point and target point');
xlabel('x' );
ylabel('y');
zlabel('z');
hold off;
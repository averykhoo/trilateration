%30.505 Computational Science and Engineering
%Project: Trilateration with numerical methods

%Trilateration using Non Linear Least Squares Model
%Extension from Linear Least Squares code attempted

%Code written by: Alvin Goh
%Acknowledgements: Avery Khoo for code checking and advice

%LSQ is to get a linear model relation between measured distance and actual
%distance between anchor and tag using one of the anchors as reference
%point, and then finding the coordinate based on this relationship.
%For NLSQ, use LSQ as initial guess.
%Aim of NLSQ is to reduce the error between actual distance and measured 
%distance (from anchor to tag) to the minimum. 

%Source: "Trilateration: The Mathematics behind a Local Positioning System"
%by Willy Hereman, Turgut Ozal University, Department of Computer
%Engineering, 21st June 2011

%Difference from previous versions
%1. Changed calculation of distances Ax=b to be more general. Still
%requires manual input at the front, but is able to handle otherwise.
%2. Utilising MATLAB's native functions such as lsqr, fmincon, lsqnonlin
%instead of hardcoding the matrix calculations.

%All values are in metres.
%Accruacy is up to mm i.e. 3 decimal places if expressed in m.

clc;
clear;

%%%%%%%%%%%%%%%%%%%
%Initial variables
%%%%%%%%%%%%%%%%%%%

%Coordinates of anchors are [x; y; z]
%Let reference anchor be A1
A1 = [0; 0; 0]; %For simplicity, start as close to origin as possible
A2 = [9; 9; 5]; %extreme end for positive x and y axis, use for z-axis check
A3 = [3; 5; 0]; % positive x, positive y y
A4 = [5; -3; 0]; %positive x, negative 
A5 = [-5; 5; 0]; % negative x, positive y
A6 = [-9; -9; 5]; %extreme end for negative x and y axis, use for z-axis check
% A7 = [-3; -7; 0]; %negative x, negative y, not to extreme end
% A8 = [6; 8; 3]; %wild card, as long as not far away from extreme end

%Create matrix system of anchor coordinates.
%As a result, each row correspond to a dimension/axis i.e. x,y,z
%Each column corresponds to the specific anchor
%Anchor = [A1 A2 A3 A4]; %Unique 2D solution
Anchor = [A1 A2 A3 A4 A5 A6]; %Unique 3D solution
% Anchor = [A1 A2 A3 A4 A5 A6 A7 A8]; %Examine accuracy effects from
%overdetermined system

n = 6; %sorry, manual input. n is total number of anchors.

%Put in intended tag coordinates here, this is not meant for calculations
%but for comparison
Tag = [5; 5; 0];

%Set radius measurement

%If you want to check base accuracy, distances should ideally be close to
%the actual distances between anchor and tag. Subsequently, put in
%reasonable figures to see the limits of the algorithm

%Input distance from anchor to tag here
%Reminder: Accruacy is up to mm i.e. 3 decimal places if expressed in m.
r1 = 7; %r1 is measured distance between reference anchor and tag
r2 = 7;
r3 = 2;
r4 = 8;
r5 = 10;
r6 = 20;
% r7 = 17;
% r8 = 4;

%Create matrix system of anchor measurements.
%R = [r1 r2 r3 r4]; %Underdetermined for 3D, unique solution for 2D
R = [r1 r2 r3 r4 r5 r6]; %Unique 3D solution
%R = [r1 r2 r3 r4 r5 r6 r7 r8]; %Examine accuracy effects from
%overdetermined system

%Add measurement noise
%Decawave TREK1000 has accruacy <0.2m
acc_upper = 0.2; %Accuracy upper bound
acc_lower = -0.2; %Accuracy lower bound
%add noise to measurements using random number generator
% (acc_lower + (acc_upper - acc_lower).*rand)

%%%%%%%%%%%%
%LSQ portion
%%%%%%%%%%%%
%Use LSQ to get a initial guess

%Create matrix of actual distance between anchor and reference anchor
%d(i)(r) is distance between anchor and reference anchor; i = 1 to n, r is reference 
d_ir = zeros(1,n-1);
for i = 2:n
    %get actual distance from one anchor to the reference anchor
    d_ir(1,i-1) = sqrt((Anchor(1,i)-Anchor(1,1))^2 +...
                    (Anchor(2,i)-Anchor(2,1))^2 +...
                    (Anchor(3,i)-Anchor(3,1))^2);
end

%Set up Linear Least Squares model
%Ax = b, assume A'A is non-singular and well-conditioned; A' is tranpose in MATLAB terms
%A is state matrix of distances between anchor and reference anchor
%A = [transpose(A2-A1); transpose(A3-A1); transpose(A4-A1)];
A = zeros(n-1,3);
for i = 2:n
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
    b(i-1,1) = 0.5*(R(1)^2 - R(i)^2 + d_ir(i-1)^2);
end

%Implement MATLAB "lqsr" function to find LSQ coordinates (without
%reference tag compensation)
%Sum of square error for LSQ has to be found separately because of this
%compensation
maxit = 100; %maximum iterations
tol = 10^-4; %tolerance, setting accuracy up to cm.
x = lsqr(A,b,tol,maxit);

%Actual coordinates of tag by LSQ method
%Because x in this case is really [x-x(A1); y-y(A1); z-z(A1)], need to add
%back A1's coordinates to find actual tag coordinates by LSQ method
x_lsq = x + A1;

%Let there be error variable e(i), where e(i) is difference between actual 
%distance and measured distance between tag and anchor.
%ei = sqrt((x-xi)^2+(y-yi)^2+(z-zi)^2) - r(i), i = 1 to n
lsq_e = zeros(1,n); %create matrix of error of each anchor, 1xn matrix
for i = 1:n
    lsq_e(1,i) = sqrt((x_lsq(1) - Anchor(1,i))^2 + ...
                      (x_lsq(2) - Anchor(2,i))^2 + ...
                      (x_lsq(3) - Anchor(3,i))^2) - R(i);
end
%Compile sum of square error for LSQ method
LSQ_square_error = 0;
for i = 1:n
    LSQ_square_error = LSQ_square_error + lsq_e(i)^2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Non-linear squares portion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Use lsqnonlin function from MATLAB
%function is error_sum, which is the error found between measured distance
%and actual distance between each anchor and the tag
%initial guess x0 is obtained from the LSQ function of MATLAB
%[x_nlsq,resnorm] = lsqnonlin(@error_sum,x_lsq);

%fmincon function start%
%f here is for implementation with fmincon function
%Manual input of error function; to find better way to implement for
%lsqnonlin function
f = @(x) (sqrt((x(1) - Anchor(1,1))^2 + (x(2) - Anchor(2,1))^2 + (x(3) - Anchor(3,1))^2) - R(1))^2 + ...
         (sqrt((x(1) - Anchor(1,2))^2 + (x(2) - Anchor(2,2))^2 + (x(3) - Anchor(3,2))^2) - R(2))^2 + ...
         (sqrt((x(1) - Anchor(1,3))^2 + (x(2) - Anchor(2,3))^2 + (x(3) - Anchor(3,3))^2) - R(3))^2 + ...
         (sqrt((x(1) - Anchor(1,4))^2 + (x(2) - Anchor(2,4))^2 + (x(3) - Anchor(3,4))^2) - R(4))^2 + ...
         (sqrt((x(1) - Anchor(1,5))^2 + (x(2) - Anchor(2,5))^2 + (x(3) - Anchor(3,5))^2) - R(5))^2 + ...
         (sqrt((x(1) - Anchor(1,6))^2 + (x(2) - Anchor(2,6))^2 + (x(3) - Anchor(3,6))^2) - R(6))^2;
% f = @(x) (sqrt((x(1) - Anchor(1,i))^2 + (x(2) - Anchor(2,i))^2 + (x(3) - Anchor(3,i))^2) - R(i))^2;
% F = symsum(f,k,1,n);

%Display options
options = optimoptions(@fmincon,...
     'Display','iter','Algorithm','interior-point');
%Use fmincon to find the NLSQ coordinates and the minimum sum of square
%errors
[x_nlsq,resnorm] = fmincon(f,x_lsq,[],[],[],[],[],[],[],options);
% [x_nlsq,resnorm] = fmincon(F,x_lsq,[],[],[],[],[],[],[],options);

%fmincon function end%

%Let there be error variable e(i), where e(i) is difference between actual 
%distance and measured distance between tag and anchor.
%ei = sqrt((x-xi)^2+(y-yi)^2+(z-zi)^2) - r(i), i = 1 to n
nlsq_e = zeros(1,n); %create matrix of error of each anchor
for i = 1:n
    nlsq_e(1,i) = sqrt((x_nlsq(1) - Anchor(1,i))^2 + ...
                       (x_nlsq(2) - Anchor(2,i))^2 + ...
                       (x_nlsq(3) - Anchor(3,i))^2) - R(i);
end
%Compile sum of square error for NLSQ method
NLSQ_square_error = 0;
for i = 1:n
    NLSQ_square_error = NLSQ_square_error + nlsq_e(i)^2;
end

%%%%%%%%%%%
%True Error
%%%%%%%%%%%
%find true error of lsq and nlsq coordinates
te_lsq = zeros(1,3);
for i=1:3
    te_lsq(i) = Tag(i) - x_lsq(i);
end

te_nlsq = zeros(1,3);
for i=1:3
    te_nlsq(i) = Tag(i) - x_nlsq(i);
end

%Find actual distance between anchor and tag
D = zeros(1,n);
for i = 1:n
    %get actual distance from one anchor to the reference anchor
    D(1,i) = sqrt((Tag(1) - Anchor(1,i))^2 +...
                  (Tag(2) - Anchor(2,i))^2 +...
                  (Tag(3) - Anchor(3,i))^2);
end

%find true error from tag to each anchor
tag_e = zeros(1,n); %create matrix of error of each anchor
for i = 1:n
    tag_e(1,i) = sqrt((Tag(1) - Anchor(1,i))^2 + ...
                      (Tag(2) - Anchor(2,i))^2 + ...
                      (Tag(3) - Anchor(3,i))^2) - R(i);
end
%Compile sum of square error for true coordinates
True_square_error = 0;
for i = 1:n
    True_square_error = True_square_error + tag_e(i)^2;
end

%%%%%%%%%%%%%%%%%%
%Display Results
%%%%%%%%%%%%%%%%%%
%LSQ results
fprintf('Raw LSQ coordinates were (%.4f, %.4f, %.4f) \n', ...
        x(1), x(2), x(3));
fprintf('Actual LSQ coordinates were (%.4f, %.4f, %.4f) \n', ...
        x_lsq(1), x_lsq(2), x_lsq(3));
disp('LSQ measurement error of each anchor is'); disp(lsq_e);
fprintf('Sum of square error from LSQ is %g \n', LSQ_square_error);

%NLSQ results
fprintf('NLSQ coordinates were (%.4f, %.4f, %.4f) \n', ...
        x_nlsq(1), x_nlsq(2), x_nlsq(3));
disp('NLSQ measurement error of each anchor is'); disp(nlsq_e);
fprintf('Sum of square error from NLSQ is %g \n', NLSQ_square_error);

%Comparing MATLAB residual squre error from MATLAB Method
fprintf('MATLAB Sum of square error from NLSQ is %g \n', resnorm);

%Display true error between tag coordinates and calculated coordinates
fprintf('True error of LSQ coordinates is (%g, %g, %g) \n', ...
        te_lsq(1), te_lsq(2), te_lsq(3));
fprintf('True error of NLSQ coordinates is (%g, %g, %g) \n', ...
        te_nlsq(1), te_nlsq(2), te_nlsq(3));
disp('Measured Distance between Tag and anchor is'); disp(R);
disp('Actual Distance between Tag and anchor is'); disp(D);
disp('True measurement error of each anchor is'); disp(tag_e);
fprintf('Sum of true square error is %g \n', True_square_error);
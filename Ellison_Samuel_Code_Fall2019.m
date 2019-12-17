%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% HORIZONTAL PLATFORM MECHANISM %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% SAMUEL ELLISON %%%%%%%%%%%%%
%%%%%%%%%%%%%% 204977052 %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
clear all;
close all;

%% PART 1: SOLVING FOR CENTER OF MASS AND ANGULAR DISPLACEMENTS OF TABLE

% Chebyshev dimensions
r1 = 810;
r2 = 1.25*r1;
r3 = 0.5*r1;
r4 = r2;
% Ideal dimensions
r5 = 0.5*r2;
r6 = r3;
r26 = 0.5*r3;
r15 = 0.5*r2;
table_height = 175+30;
% Link Properties (keep mm)
rho = 2698.9/(1000^3);
r = [r2 r3 r4 r5 r6];
V = r*900;
m = rho*V;
g = 9807;
I = zeros(1,5);
for i = 1:1:5
    I(i) = (m(i)/12)*(30^2+r(i)^2);
end

% Theta2 range
theta2_min = 36.87*pi/180;
theta2_max = 96*pi/180;

% Initialize angle arrays, a complete cycle for theta 2
theta21 = linspace(theta2_min, theta2_max, 1000);
theta22 = linspace(theta2_max, theta2_min, 1000);
theta2 = [theta21 theta22];
theta2(1000) = [];
theta3 = zeros(1,numel(theta2));
theta4 = zeros(1,numel(theta2));
theta5 = zeros(1,numel(theta2));
theta6 = zeros(1,numel(theta2));

% Initial guesses
theta3(1) = pi/1.5;
theta4(1) = pi/2;
theta5(1) = 220*pi/180;
theta6(1) = 0;

% Newton-Raphson tolerance
eps = 10^-7;

% NR Vector Loop #1, Position
for i = 1:1:numel(theta2)
    
    t3 = theta3(i);
    t4 = theta4(i);
    
    J1 = zeros(2,2);
    B1 = zeros(2,1);
    C1 = zeros(2,1);
    
    dt3 = Inf;
    dt4 = Inf;
    
    while abs(dt3) > eps || abs(dt4) > eps
        
        % evaluate functions
        f1 = r2*cos(theta2(i))+r3*cos(t3)-r4*cos(t4)-r1;
        f2 = r2*sin(theta2(i))+r3*sin(t3)-r4*sin(t4);
        
        % evaluate Jacobian
        J1(1,1) = -r3*sin(t3);
        J1(1,2) = r4*sin(t4);
        J1(2,1) = r3*cos(t3);
        J1(2,2) = -r4*cos(t4);
        
        % Fill in matrices/vectors
        A1 = -inv(J1);
        B1(1,1) = f1;
        B1(2,1) = f2;
        C1(1,1) = t3;
        C1(2,1) = t4;
        
        % Solve for new solutions
        X1 = A1*B1 + C1;
        
        % update dt
        dt3 = X1(1) - t3;
        dt4 = X1(2) - t4;
        
        % update current solutions and index
        t3 = X1(1);
        t4 = X1(2);
        
    end
    
    % Fill current index with solution
    theta3(i) = t3;
    theta4(i) = t4;
    
    % Update guess
    if i < numel(theta2)
        theta3(i+1) = theta3(i);
        theta4(i+1) = theta4(i);
    end
    
end

% NR Vector Loop #2, Position
for i = 1:1:numel(theta2)
    
    t5 = theta5(i);
    t6 = theta6(i);
    
    J2 = zeros(2,2);
    B2 = zeros(2,1);
    C2 = zeros(2,1);
    
    dt5 = Inf;
    dt6 = Inf;
    
    while abs(dt5) > eps || abs(dt6) > eps
        
        % evaluate functions
        f3 = r2*cos(theta2(i))+r26*cos(theta3(i))+r6*cos(t6)+r5*cos(t5)-r15*cos(theta4(i))-r1;
        f4 = r2*sin(theta2(i))+r26*sin(theta3(i))+r6*sin(t6)+r5*sin(t5)-r15*sin(theta4(i));
        
        % evaluate Jacobian
        J2(1,1) = -r5*sin(t5);
        J2(1,2) = -r6*sin(t6);
        J2(2,1) = r5*cos(t5);
        J2(2,2) = r6*cos(t6);
        
        % Fill in matrices/vectors
        A2 = -inv(J2);
        B2(1,1) = f3;
        B2(2,1) = f4;
        C2(1,1) = t5;
        C2(2,1) = t6;
        
        % Solve for new solutions
        X2 = A2*B2 + C2;
        
        % update dt
        dt5 = X2(1) - t5;
        dt6 = X2(2) - t6;
        
        % update current solutions and index
        t5 = X2(1);
        t6 = X2(2);
        
    end
    
    % Fill current index with solution
    theta5(i) = t5;
    theta6(i) = t6;
    
    % Update guess
    if i < numel(theta2)
        theta5(i+1) = theta5(i);
        theta6(i+1) = theta6(i);
    end
    
end

% Center of Mass vectors wrt center of coupler
xcm = zeros(1,numel(theta2));
ycm = zeros(1,numel(theta2));
% Center of Mass positioning wrt center of coupler
xcm0 = 202.5;
ycm0 = 129.2994;
D = sqrt(xcm0^2+ycm0^2); % diagonal
tcm0 = atan(ycm0/xcm0); % angle

% Position of COM over entire theta2 interval
for i = 1:1:numel(theta2)
    xcm(i) = r2*cos(theta2(i))+r26*cos(theta3(i))+D*cos(tcm0+theta6(i));
    ycm(i) = r2*sin(theta2(i))+r26*sin(theta3(i))+D*sin(tcm0+theta6(i));
    % Displacement from original position
    xcmd(i) = xcm(i) - xcm(1);
    ycmd(i) = ycm(i) - ycm(1);
end

% Maximum Displacements of horizontal/vertical/angular COM
xcm_disp = max(xcm)-min(xcm);
ycm_disp = max(ycm)-min(ycm);
tdisp = (max(theta6)-min(theta6))*180/pi;

% Print COM data
fprintf('Horizontal Movement: %5.3f\n',xcm_disp);
fprintf('Vertical Movement: %5.3f\n',ycm_disp);
fprintf('Angular Movement: %5.3e\n',tdisp);
if xcm_disp > 850 && ycm_disp < 4 && tdisp < 1
    fprintf('\nAll conditions are met!\n')
else
    fprintf('\nAll conditions NOT met.\n')
end

% Convert to degrees for plotting purposes
theta2 = theta2*180/pi;
theta6 = theta6*180/pi;

% Plots
% Horizontal
figure(1)
hold on
grid on
plot(theta2,xcmd,'r-','LineWidth',2)
xlabel('Theta 2 [degrees]')
ylabel('Horizontal Displacement [mm]')
title('Displacement of Table Center of Mass, Horizontal')
% Vertical
figure(2)
hold on
grid on
plot(theta2,ycmd,'b-','LineWidth',2)
xlabel('Theta 2 [degrees]')
ylabel('Vertical Displacement [mm]')
title('Displacement of Table Center of Mass, Vertical')
% Angular
figure(3)
hold on
grid on
plot(theta2,theta6,'k-','LineWidth',0.5)
xlabel('Theta 2 [degrees]')
ylabel('Angular Position of Table [degrees]')
title('Angular Rotation of Table')
% Position
figure(4)
hold on
grid on
axis([100 1100 0 1000])
plot(xcm,ycm,'g-','LineWidth',2)
xlabel('Horizontal Position of Table COM [mm]')
ylabel('Vertical Position of Table COM [mm]')
title('Position of Table COM in Cartesian Space')

%% PART 2: SOLVING FOR THE INPUT TORQUE AS A FUNCTION OF TIME

% Constants in theta2 equation
w = pi/180;
A = ((theta2_max-theta2_min)/2);
B = (theta2_max+theta2_min)/2;

% Create a new theta2 array based on crank (3 cycles = 6pi/w)
t = linspace(0,6*pi/w,1000);
theta2 = A*sin(w*t)+B;
omega2 = w*A*cos(w*t);
alpha2 = -w^2*A*sin(w*t);

% New angle arrays
theta3 = zeros(1,numel(theta2));
theta4 = zeros(1,numel(theta2));
theta5 = zeros(1,numel(theta2));
theta6 = zeros(1,numel(theta2));

% Initial guesses
theta3(1) = pi/1.3;
theta4(1) = pi/2;
theta5(1) = 225*pi/180;
theta6(1) = 0;

% Recalculate angles over new range of theta2
% NR Vector Loop #1
for i = 1:1:numel(theta2)
    
    t3 = theta3(i);
    t4 = theta4(i);
    
    J = zeros(2,2);
    B = zeros(2,1);
    C = zeros(2,1);
    
    dt3 = Inf;
    dt4 = Inf;
    
    while abs(dt3) > eps || abs(dt4) > eps
        
        % evaluate functions
        f1 = r2*cos(theta2(i))+r3*cos(t3)-r4*cos(t4)-r1;
        f2 = r2*sin(theta2(i))+r3*sin(t3)-r4*sin(t4);
        
        % evaluate Jacobian
        J(1,1) = -r3*sin(t3);
        J(1,2) = r4*sin(t4);
        J(2,1) = r3*cos(t3);
        J(2,2) = -r4*cos(t4);
        
        % Fill in matrices/vectors
        A = -inv(J);
        B(1,1) = f1;
        B(2,1) = f2;
        C(1,1) = t3;
        C(2,1) = t4;
        
        % Solve for new solutions
        X = A*B + C;
        
        % update dt
        dt3 = X(1) - t3;
        dt4 = X(2) - t4;
        
        % update current solutions and index
        t3 = X(1);
        t4 = X(2);
        
    end
    
    theta3(i) = t3;
    theta4(i) = t4;
    
    if i < numel(theta2)
        theta3(i+1) = theta3(i);
        theta4(i+1) = theta4(i);
    end
    
end

% NR Vector Loop #2
for i = 1:1:numel(theta2)
    
    t5 = theta5(i);
    t6 = theta6(i);
    
    J = zeros(2,2);
    B = zeros(2,1);
    C = zeros(2,1);
    
    dt5 = Inf;
    dt6 = Inf;
    
    while abs(dt5) > eps || abs(dt6) > eps
        
        % evaluate functions
        f3 = r2*cos(theta2(i))+r26*cos(theta3(i))+r6*cos(t6)+r5*cos(t5)-r15*cos(theta4(i))-r1;
        f4 = r2*sin(theta2(i))+r26*sin(theta3(i))+r6*sin(t6)+r5*sin(t5)-r15*sin(theta4(i));
        
        % evaluate Jacobian
        J(1,1) = -r5*sin(t5);
        J(1,2) = -r6*sin(t6);
        J(2,1) = r5*cos(t5);
        J(2,2) = r6*cos(t6);
        
        % Fill in matrices/vectors
        A = -inv(J);
        B(1,1) = f3;
        B(2,1) = f4;
        C(1,1) = t5;
        C(2,1) = t6;
        
        % Solve for new solutions
        X = A*B + C;
        
        % update dt
        dt5 = X(1) - t5;
        dt6 = X(2) - t6;
        
        % update current solutions and index
        t5 = X(1);
        t6 = X(2);
        
    end
    
    theta5(i) = t5;
    theta6(i) = t6;
    
    if i < numel(theta2)
        theta5(i+1) = theta5(i);
        theta6(i+1) = theta6(i);
    end
    
end

% First & Second Order Kinematic Coefficients, Vector Loop 1
t3p = zeros(1,numel(theta2));
t4p = zeros(1,numel(theta2));
t3pp = zeros(1,numel(theta2));
t4pp = zeros(1,numel(theta2));

for i = 1:1:numel(theta2)
    J1 = [-r3*sin(theta3(i)) r4*sin(theta4(i)); r3*cos(theta3(i)) -r4*cos(theta4(i))];
    b1 = [r2*sin(theta2(i)); -r2*cos(theta2(i))];
    Xp = J1\b1;
    t3p(i) = Xp(1);
    t4p(i) = Xp(2);
    b2 = [r2*cos(theta2(i))+r3*t3p(i)^2*cos(theta3(i))-r4*t4p(i)^2*cos(theta4(i)); r2*sin(theta2(i))+r3*t3p(i)^2*sin(theta3(i))-r4*t4p(i)^2*sin(theta4(i))];
    Xpp = J1\b2;
    t3pp(i) = Xpp(1);
    t4pp(i) = Xpp(2);
end

% First & Second Order Kinematic Coefficients, Vector Loop 2
t5p = zeros(1,numel(theta2));
t6p = zeros(1,numel(theta2));
t5pp = zeros(1,numel(theta2));
t6pp = zeros(1,numel(theta2));

for i = 1:1:numel(theta2)
    J2 = [-r5*sin(theta5(i)) -r6*sin(theta6(i)); r5*cos(theta5(i)) r6*cos(theta6(i))];
    b3 = [r2*sin(theta2(i))+r26*t3p(i)*sin(theta3(i))-r15*t4p(i)*sin(theta4(i)); -r2*cos(theta2(i))-r26*t3p(i)*cos(theta3(i))+r15*t4p(i)*cos(theta4(i))];
    Yp = J2\b3;
    t5p(i) = Yp(1);
    t6p(i) = Yp(2);
    b4 = [r2*cos(theta2(i))+r26*t3pp(i)*sin(theta3(i))+r26*t3p(i)^2*cos(theta3(i))+r6*t6p(i)^2*cos(theta6(i))+r5*t5p(i)^2*cos(theta5(i))-r15*t4pp(i)*sin(theta4(i))-r15*t4p(i)^2*cos(theta4(i)); r2*sin(theta2(i))-r26*t3pp(i)*cos(theta3(i))+r26*t3p(i)^2*sin(theta3(i))+r6*t6p(i)^2*sin(theta6(i))+r5*t5p(i)^2*sin(theta5(i))+r15*t4pp(i)*cos(theta4(i))-r15*t4p(i)^2*sin(theta4(i))];
    Ypp = J2\b4;
    t5pp(i) = Ypp(1);
    t6pp(i) = Ypp(2);
end

% Link COM positions, velocities, and accelerations wrt theta 2
% Link 2
xG2 = 0.5*r2*cos(theta2);
yG2 = 0.5*r2*sin(theta2);
xG2p = -0.5*r2*sin(theta2);
yG2p = 0.5*r2*cos(theta2);
xG2pp = -0.5*r2*cos(theta2);
yG2pp = -0.5*r2*sin(theta2);
% Link 3
xG3 = r2*cos(theta2)+0.5*r3*cos(theta3);
yG3 = r2*sin(theta2)+0.5*r3*sin(theta3);
xG3p = zeros(1,numel(theta2));
yG3p = zeros(1,numel(theta2));
xG3pp = zeros(1,numel(theta2));
yG3pp = zeros(1,numel(theta2));
for i = 1:1:numel(theta2)
    xG3p(i) = -r2*sin(theta2(i))-0.5*r3*t3p(i)*sin(theta3(i));
    yG3p(i) = r2*cos(theta2(i))+0.5*r3*t3p(i)*cos(theta3(i));
    xG3pp(i) = -r2*cos(theta2(i))-0.5*r3*t3pp(i)*sin(theta3(i))-0.5*r3*t3p(i)^2*cos(theta3(i));
    yG3pp(i) = -r2*sin(theta2(i))+0.5*r3*t3pp(i)*cos(theta3(i))-0.5*r3*t3p(i)^2*sin(theta3(i));    
end
% Link 4
xG4 = r2*cos(theta2)+r3*cos(theta3)-0.5*r4*cos(theta4);
yG4 = r2*sin(theta2)+r3*sin(theta3)-0.5*r4*sin(theta4);
xG4p = zeros(1,numel(theta2));
yG4p = zeros(1,numel(theta2));
xG4pp = zeros(1,numel(theta2));
yG4pp = zeros(1,numel(theta2));
for i = 1:1:numel(theta2)
    xG4p(i) = -r2*sin(theta2(i))-r3*t3p(i)*sin(theta3(i))+0.5*r4*t4p(i)*sin(theta4(i));
    yG4p(i) = r2*cos(theta2(i))+r3*t3p(i)*cos(theta3(i))-0.5*r4*t4p(i)*cos(theta4(i));
    xG4pp(i) = -r2*cos(theta2(i))-r3*t3pp(i)*sin(theta3(i))-r3*t3p(i)^2*cos(theta3(i))+0.5*r4*t4pp(i)*sin(theta4(i))+0.5*r4*t4p(i)^2*cos(theta4(i));
    yG4pp(i) = -r2*sin(theta2(i))+r3*t3pp(i)*cos(theta3(i))-r3*t3p(i)^2*sin(theta3(i))-0.5*r4*t4pp(i)*cos(theta4(i))+0.5*r4*t4p(i)^2*sin(theta4(i));
end
% Link 5
xG5 = xG3+r6*cos(theta6)+0.5*r5*cos(theta5);
yG5 = yG3+r6*sin(theta6)+0.5*r5*sin(theta5);
xG5p = zeros(1,numel(theta2));
yG5p = zeros(1,numel(theta2));
xG5pp = zeros(1,numel(theta2));
yG5pp = zeros(1,numel(theta2));
for i = 1:1:numel(theta2)
    xG5p(i) = xG3p(i)-r6*t6p(i)*sin(theta6(i))-0.5*r5*t5p(i)*sin(theta5(i));
    yG5p(i) = yG3p(i)+r6*t6p(i)*cos(theta6(i))+0.5*r5*t5p(i)*cos(theta5(i));
    xG5pp(i) = xG3pp(i)-r6*t6pp(i)*sin(theta6(i))-r6*t6p(i)^2*cos(theta6(i))-0.5*r5*t5pp(i)*sin(theta5(i))-0.5*r5*t5p(i)^2*cos(theta5(i));
    yG5pp(i) = yG3pp(i)+r6*t6pp(i)*cos(theta6(i))-r6*t6p(i)^2*sin(theta6(i))+0.5*r5*t5pp(i)*cos(theta5(i))-0.5*r5*t5p(i)^2*sin(theta5(i));
end
% Link 6 (Table)
xG6 = xG3+D*cos(theta6+tcm0);
yG6 = yG3+D*sin(theta6+tcm0);
xG6p = zeros(1,numel(theta2));
yG6p = zeros(1,numel(theta2));
xG6pp = zeros(1,numel(theta2));
yG6pp = zeros(1,numel(theta2));
for i = 1:1:numel(theta2)
    xG6p(i) = xG3p(i)-D*t6p(i)*sin(theta6(i)+tcm0);
    yG6p(i) = yG3p(i)+D*t6p(i)*cos(theta6(i)+tcm0);
    xG6pp(i) = xG3pp(i)-D*t6pp(i)*sin(theta6(i)+tcm0)-D*t6p(i)^2*cos(theta6(i)+tcm0);
    yG6pp(i) = yG3pp(i)+D*t6pp(i)*cos(theta6(i)+tcm0)-D*t6p(i)^2*sin(theta6(i)+tcm0);
end
% Put Link information and kinematic coefficients in matrix form
xG = [xG2; xG3; xG4; xG5; xG6];
xGp = [xG2p; xG3p; xG4p; xG5p; xG6p];
xGpp = [xG2pp; xG3pp; xG4pp; xG5pp; xG6pp];
yG = [yG2; yG3; yG4; yG5; yG6];
yGp = [yG2p; yG3p; yG4p; yG5p; yG6p];
yGpp = [yG2pp; yG3pp; yG4pp; yG5pp; yG6pp];
tp = [ones(1,numel(theta2)); t3p; t4p; t5p; t6p];
tpp = [zeros(1,numel(theta2)); t3pp; t4pp; t5pp; t6pp];

% Torque Calculations
Torque = zeros(1,numel(theta2));
for i = 1:1:numel(theta2)
    A = zeros(1,5);
    B = zeros(1,5);
    W = zeros(1,5);
    for j = 1:1:5
        A(j) = m(j)*(xGp(j,i)^2+yGp(j,i)^2)+I(j)*tp(j,i)^2;
        B(j) = m(j)*(xGp(j,i)*xGpp(j,i)+yGp(j,i)*yGpp(j,i))+I(j)*tp(j,i)*tpp(j,i);
        W(j) = m(j)*g*yGp(j,i);
    end
    Torque(i) = sum(A)*alpha2(i)+sum(B)*omega2(i)^2+sum(W);
end

% Convert to N.m.
Torque = Torque/(1000^2);

% Plot Torque as a function of time over 3 cycles
figure(5)
hold on
grid on
plot(t,Torque,'m-','LineWidth',2)
xlabel('Time [s]')
ylabel('Torque [Nm]')
title('Input Torque of the Motor Over 3 Cycles')

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% Made by: Matas Sabaliauskas     %
% Student ID: 20199952471         %
% Korea University, Robotics 2020 %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %


syms m a g theta(t) r omega_0 x theta_0 theta_t0

eqn = m*a == -m*g*sin(theta);

eqn = subs(eqn,a,r*diff(theta,2));

eqn = isolate(eqn,diff(theta,2));

eqn = subs(eqn,g/r,omega_0^2);

approx = taylor(sin(x),x,'Order',2);
approx = subs(approx,x,theta(t));

eqnLinear = subs(eqn,sin(theta(t)),approx);

theta_t = diff(theta);
cond = [theta(0) == theta_0, theta_t(0) == theta_t0];
assume(omega_0,'real')
thetaSol(t) = dsolve(eqnLinear,cond);

gValue = 9.81;
rValue = 0.5;     %Since pendulum length is 50 cm = 0.50 m
omega_0Value = sqrt(gValue/rValue);
T = 2*pi/omega_0Value;

%Note: Solution only valid for small angles.
theta_0Value  = pi/18;  % theta_0 = 10 (deg) = pi/18 (rad)
theta_t0Value = 0;      % Initially at rest.

vars   = [omega_0      theta_0      theta_t0];
values = [omega_0Value theta_0Value theta_t0Value];
thetaSolPlot = subs(thetaSol,vars,values);



%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Video #1: theta(t) vs. time
figure(1);
t = 0:0.1:5;
y = thetaSolPlot(t*T)/pi;

myVideo = VideoWriter('Angular Displacement Video 1'); %open video file
myVideo.FrameRate = 10;
open(myVideo)

for i=1:1:length(t)
    plot(t(1:i), y(1:i), 'LineWidth', 2);
    grid on;
    title('Angular Motion');
    xlabel('Time t (s)');
    ylabel('Angular Displacement \theta (rad)');
    ylim([-1, 1])
    xlim([0, 5])
    pause(0.01) %Pause and grab frame
    frame = getframe(gcf); %get frame
    writeVideo(myVideo, frame);
end
close(myVideo)

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% This Second part solve differential equations%
% To create a horizontal and vertical displacement graphs %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %

syms theta(t) theta_t(t) omega_0
eqs = [diff(theta)   == theta_t; diff(theta_t) == -omega_0^2*sin(theta)];
   
   
eqs  = subs(eqs,omega_0,omega_0Value);
vars = [theta, theta_t];

[M,F] = massMatrixForm(eqs,vars);

f = M\F;

f = odeFunction(f, vars);

x0 = [0; 1.99*omega_0Value];

tInit  = 0;
tFinal = 10;

sols = ode45(f,[tInit tFinal],x0);

%x_pos and y_pos find solutions to the differential equation
x_pos = @(t) sin(deval(sols,t,1));
y_pos = @(t) -cos(deval(sols,t,1));



%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Video #2: X(t) and Y(t) vs. time
figure(2);
t = 0:0.01:5;
%t = 0:0.05:5;
%t = 0:0.1:5;

y = x_pos(t);
z = y_pos(t);

myVideo = VideoWriter('X(t) and Y(t) Video 2'); %open video file
myVideo.FrameRate = 10;
open(myVideo)

for i=1:1:length(t)
    yyaxis left; % X Displacement
    plot(t(1:i), y(1:i), 'LineWidth', 2)
    ylabel('Horizontal Displacement X(t)');

    yyaxis right; % Y Displacement
    plot(t(1:i), z(1:i), 'LineWidth', 2)
    ylabel('Vertical Displacement Y(t)');
    
    grid on;
    title('Cartpole Motion');
    xlabel('Time t (s)');
    
    ylim([-1, 1])
    xlim([0, 5])
    pause(0.01) %Pause and grab frame
    frame = getframe(gcf); %get frame
    writeVideo(myVideo, frame);
end
close(myVideo)


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Video #3: X(t) and Y(t)
figure(3);
t = 0:0.01:5;
y = x_pos(t);
z = y_pos(t);

myVideo = VideoWriter('X(t) and Y(t) Video 3'); %open video file
myVideo.FrameRate = 10;
open(myVideo)

for i=1:1:(length(t)/10)
    plot(x_pos(t), y_pos(t), 'LineWidth', 2);
    grid on;
    title('Y(t) vs. X(t) Plot');
    xlabel('X(t)');
    ylabel('Y(t)');
    ylim([-5, 5])
    xlim([-5, 5])
    
    pause(0.01) %Pause and grab frame
    frame = getframe(gcf); %get frame
    writeVideo(myVideo, frame);
end
close(myVideo)


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Video #4: Cartpole Motion Simulation video 4

figure(4);
fanimator(@(t) plot(x_pos(t),y_pos(t),'ko','MarkerFaceColor','k'));
hold on;
fanimator(@(t) plot([0 x_pos(t)],[0 y_pos(t)],'k-'));
fanimator(@(t) text(-0.3,1.5,"Timer: "+num2str(t,2)+" s"));
playAnimation

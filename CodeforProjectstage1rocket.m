function CodeforProjectstage1rocket
clc
clear
close all
% givens from problem statement and first term calculations for basic
% variables in use later in code
g0 = 9.807; 
Re = 6378e3; 
ho = 7500; 
rho0 = 1.225; 
CD = 0.5; 
Ttokg = 1450*907;
mPL = 5000;
m0 = Ttokg + mPL;
TtoN = 9806.65*1746;
% T = 17122.41; %kn
T = TtoN/1000; %kn
tb0 = 240;
Isp = 310;
d = 4.6; 
A = pi/4*(d)^2;
mo = 13631243.5 + mPL;
go = 9.807/1000;
deltat = tb0-0;
mf = (-T*deltat)/(Isp*go)+ mo;
n = mo/mf;
mfinal = m0/n;
% Thrust = 17122410; %n
Thrust = TtoN; %n
m_dot = (mo-mf)/tb0;
mp = m0 - mfinal;
tburn = 240; 
hturn = 130;
deg = pi/180; 
% Initial conditions for function
t0 = 0;
tf = tburn; 
tspan = [t0,tf]; 
% initial values for integration
v0 = 0; 
gamma0 = deg2rad(89.85); 
x0 = 0; 
h0 = 0;
vD0 = 0; 
vG0 = 0; 
%Initial conditions vector:
Initialv = [v0; gamma0; x0; h0; vD0; vG0];
% ode45 function solves the differential equations of motion using initial
% conditions and integrates over time until burnout
[t1,z] = ode45(@Rocketeq2flight, tspan, Initialv);
%  Rocketeq2flight is an embedded function using the differential eqs with 
% z giving the final values based on the time interval of integration
v     =  z(:,1)/1000;
gamma =  z(:,2)/deg; 
x     =  z(:,3)/1000; 
h     =  z(:,4)/1000; 
vD    = -z(:,5)/1000;
vG    = -z(:,6)/1000; 
% for loop that solves for the dynamic pressure over time and as a function
% of height as the launch vehicle continues its ascent. Using atmosisa
% function that returns atmospheric values based on given altitudes to
% solve for mach number at each time increment
for i = 1:length(t1)
   Rho = rho0 * exp(-h(i)*1000/ho); 
   q(i) = 1/2*Rho*(v(i)*1000)^2; 
   [f r(i) f f] = atmosisa(h(i)*1000); 
   M(i) = 1000*v(i)/r(i); 
end
% Max dynamic pressure and time of occurrence solved for using values from for loop
[maxQ,imax] = max(q); 
tQ = t1(imax); 
vQ = v(imax); 
hQ = h(imax); 
[f rQ f f] = atmosisa(h(imax)*1000); 
MQ = 1000*vQ/rQ;
Answers
return
function dzdt = Rocketeq2flight(t,z)
% Function calculates rate of change in each variable for the given
% equations of motion for gravity turn path
% initial  vector is defined
dzdt = zeros(6,1);
v     = z(1); 
gamma = z(2);
x     = z(3);
h     = z(4);
vD    = z(5); 
vG    = z(6);
% for loop used to calculate mass over burntime and the resulting thrust at
% each condition of launch vehicle mass
if t < tburn
   m = m0 - m_dot*t; 
   T = Thrust;
else
   m = m0 - m_dot*tburn;
   T = 0; 
end
% Density vs altitude, gravitational change with altitude and drag at given
% velocity calculated for later steps
g = g0/(1 + h/Re)^2; 
rho = rho0*exp(-h/ho); 
D = 0.5*rho*v^2*A*CD; 
% First if statement for vertical flight using differential equations of motion up
% until altitude is above 130m. Then gravity assist begins in else
% statement
if h <= hturn
   gamma_dot = 0;
   v_dot  = T/m - D/m - g;
   x_dot  = 0;
   h_dot  = v;
   vG_dot = -g;
else
   % when flight reaches above 130m, gravity turn begins using differential equations of
   % motion to calculate values for the vehicle as it continues it ascent
   v_dot = T/m - D/m - g*sin(gamma);
   gamma_dot = -1/v*(g - v^2/(Re + h))*cos(gamma);
   x_dot     = Re/(Re + h)*v*cos(gamma); 
   h_dot     = v*sin(gamma); 
   vG_dot    = -g*sin(gamma); 
end                            
vD_dot = -D/m;                 
% Values for the equations of motion calculated stored into matrix
dzdt(1) = v_dot;
dzdt(2) = gamma_dot;
dzdt(3) = x_dot;
dzdt(4) = h_dot;
dzdt(5) = vD_dot;
dzdt(6) = vG_dot;
end 
function Answers
    % printing of important values of solution
fprintf('\n Speed at Burnout        = %f km/s',v(end))
fprintf('\n Altitude at Burnout       = %f km ',h(end))
fprintf('\n Downrange distance        = %f km ',x(end))
fprintf('\n Mass of Propellent Expended = %f kg',mp)
figure (1)
plot(x,h)
title('Altitude vs Downrange Distance (Cape Canaveral)')
xlabel('Downrange Distance (km)')
ylabel('Altitude (km)')
grid
end 
end

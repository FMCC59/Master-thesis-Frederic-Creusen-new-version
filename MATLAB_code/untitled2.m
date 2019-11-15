Fibre_pos


xp = 1
yp = 1
 t = linspace(0,2*pi);              %linearly spaced vector (2pi)
 X = a*cos(t);                      %Definition of X axis of circumpherence
 Y = b*sin(t);                      %Definition of Y axis of circumpherence
 x = xp + X; %Definition of X axis of circumpherence
 y = yp + Y; %Definition of Y axis of circumpherence
 plot(x,y,'r-')
 axis equal
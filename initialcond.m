vw=10;  %Average Wind speed Magnitude
th=pi/3;  %Average angle of wind vector (on the xy plane)
vwx0=vw*cos(th);  %x-component of wind speed
vwy0=vw*sin(th);  %y-component of wind speed
v1=3;   %Lower Interval of wind speed magnitude
v2=25;   %Higher Interval of wind speed magnitude
th1=pi/18;  %Lower interval of angle of wind speed vector
th2=pi/2.5;  %Higher interval of angle of wind speed vector

vsx0=2;  %Initial velocity of sensor nodes
vsy0=4;
vsz0=10;

m=0.05; %Mass of sensor nodes
kw=0.4;  
kar=0.015;

ppx=20; %Initial position of nodes when released
ppy=20;
ppz=1000;

mex=0;
mey=0;

sigmax=2.58; 
sigmay=2.58;

g=3.69; %Mars gravitational acceleration

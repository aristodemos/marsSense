clear all
clc
%Calculation of Probability Density Function for Sensor Node Placement
initialcond   %Call values for Initial Conditions

%Calculate Wind Speed Probability Density Functions (for x and y). 
%'IF' conditions apply in regards to the values for the resultants of
%initial wind speeds (vwx and vwy)   
%x-direction
cvw=1/((th2-th1)*(v2-v1));  %Calculate constants for the Wind Speed(x-direction) PDF 
cvwx1=cvw*log((1-sqrt(1-(cos(th1))^2))/cos(th1));
cvwx2=cvw*log((1-sqrt(1-(cos(th2))^2))/cos(th2));

if (vwx0>(v1*cos(th2))&&vwx0<(v1*cos(th1)));
    vwx=v1*cos(th2):0.001:v2*cos(th2);
    vwxpdf=cvw*log((v1-sqrt((v1^2)-(vwx.^2)))./vwx)-cvwx2;
   
elseif (vwx0>(v1*cos(th1))&&vwx0<(v2*cos(th2)));
    vwx=v1*cos(th1):0.001:v2*cos(th2);
    vwxpdf=cvwx1-cvwx2;
   
elseif (vwx0>(v2*cos(th2))&&vwx0<(v2*cos(th1)));
    vwx=v2*cos(th2):0.001:v2*cos(th1);
    vwxpdf=cvwx1-cvw*log((v2-sqrt((v2^2)-(vwx.^2)))./vwx);
else vwxpdf=0;
end

%y-direction
cvwy1=cvw*log((1-sqrt(1-(sin(th1))^2))/sin(th1));  %Calculate constants for the Wind Speed(y-direction) PDF
cvwy2=cvw*log((1-sqrt(1-(sin(th2))^2))/sin(th2));
if (vwy0>(v1*sin(th1))&&vwy0<(v1*sin(th2)));
    vwy=v1*sin(th1):0.001:v1*sin(th2);
    vwypdf=cvw*log((v1-sqrt((v1^2)-(vwy.^2)))./vwy)-cvwy1;
elseif (vwy0>(v1*sin(th2))&&vwy0<(v2*sin(th1)));
    vwy=v1*sin(th2):0.001:v2*sin(th1);
    vwypdf=cvwy2-cvwy1;
elseif (vwy0>(v2*sin(th1))&&vwy0<(v2*sin(th2)));
    vwy=v2*sin(th1):0.001:v2*sin(th2);
    vwypdf=cvwy2-cvw*log((v2-sqrt((v2^2)-(vwy.^2)))./vwy);
else vwypdf=0;
end

%Calculate velocity of sensors (x, y, z)
cvs=kw/kar; %Constants
cvsx=cvs*vwx0+vsx0;
cvsy=cvs*vwy0+vsy0;
cvsz=-((m*g*cvs)/kw)+vsz0;
T=20;
t=0:0.001:T;
for i=1:length(t)
    vsx(i)=cvsx*exp(-(kar.*t(i))/m)-cvs*vwx0;
    vsy(i)=cvsy*exp(-(kar.*t(i))/m)-cvs*vwy0;
    vsz(i)=cvsz*exp(-(kar.*t(i))/m)+m*g*cvs/kw;
end

%Calculate positions of sensors on ground
cp=kw*T/kar+m*(kw/(kar^2))*(exp(-kar*T/m)-1); %Constants
cpx=ppx-m*(vsx0/kar)*(exp(-kar*T/m)-1);
cpy=ppy-m*(vsy0/kar)*(exp(-kar*T/m)-1);
cpz=ppz-m*(vsz0/kar)*(exp(-kar*T/m)-1);

psxfin=cp.*vwx+cpx;
psyfin=cp.*vwy+cpy;


%Calculation of the Sensor Ground Position PDF 
%in regards to the Wind Speed PDF
%introduce psxfin1 instead of psxfin, to account for the mathematical model 
if (vwx0>(v1*cos(th2))&&vwx0<(v1*cos(th1)));
    psxfin1=(cpx+cp*v1*cos(th2)):0.01:(cpx+cp*v2*cos(th1));
    psxfinpdf=(1/cp)*cvw*log((v1-sqrt((v1^2)-(((psxfin1-cpx)/cp).^2)))./((psxfin1-cpx)/cp))-cvwx2;
elseif (vwx0>(v2*cos(th2))&&vwx0<(v2*cos(th1)));
    psxfin1=(cpx+cp*v1*cos(th2)):0.01:(cpx+cp*v2*cos(th1));
    psxfinpdf=(1/cp)*(cvwx1-cvw*log((v2-sqrt((v2^2)-(((psxfin1-cpx)/cp).^2)))./((psxfin1-cpx)/cp)));
elseif (vwx0>(v1*cos(th1))&&vwx0<(v2*cos(th2)));
    psxfinpdf=(1/cp)*(cvwx1-cvwx2);
else psxfinpdf=0;   
end

if (vwy0>=(v1*sin(th1))&&vwy0<=(v1*sin(th2)));
    psyfin1=(cpy+cp*v1*sin(th1)):0.01:(cpy+cp*v2*sin(th2));
    psyfinpdf=(1/cp)*cvw*log((v1-sqrt((v1^2)-(((psyfin1-cpy)/cp).^2)))./((psyfin1-cpy)/cp))-cvwy1;
elseif (vwy0>=(v2*sin(th1))&&vwy0<=(v2*sin(th2)));
    psyfin1=(cpy+cp*v1*sin(th1)):0.01:(cpy+cp*v2*sin(th2));
    psyfinpdf=(1/cp)*(cvwy2-cvw*log((v2-sqrt((v2^2)-(((psyfin1-cpy)/cp).^2)))./((psyfin1-cpy)/cp)));
elseif (vwy0>(v1*sin(th2))&&vwy0<(v2*sin(th1)));
    psyfinpdf=(1/cp)*(cvwy2-cvwy1);
else psyfinpdf=0;
end


%Calculate PDF for actual ground position (x and y).
%Multiply the two PDFs above and integrate
%x-direction
psxactual=-400:1:400;
px=-400:1:400;
for i=1:length(psxactual)
    for j=1:length(px)

if (vwx0>(v1*cos(th2))&&vwx0<(v1*cos(th1)));
    psxactualdist(j,i)=(((1/cp)*cvw*log((v1-sqrt((v1^2)-(((px(j)-cpx)/cp).^2)))./((px(j)-cpx)/cp))-cvwx2)*((1/(sqrt(2*pi)*sigmax)).*exp(-((psxactual(i)-px(j)-mex).^2)/(2*sigmax^2))));
elseif (vwx0>(v2*cos(th2))&&vwx0<(v2*cos(th1)));
    psxactualdist(j,i)=(((1/cp)*(cvwx1-cvw*log((v2-sqrt((v2^2)-(((px(j)-cpx)/cp).^2)))./((px(j)-cpx)/cp))))*((1/(sqrt(2*pi)*sigmax)).*exp(-((psxactual(i)-px(j)-mex).^2)/(2*sigmax^2))));
elseif (vwx0>(v1*cos(th1))&&vwx0<(v2*cos(th2)));
        psxactualdist(j,i)=(((1/cp)*(cvwx1-cvwx2))*((1/(sqrt(2*pi)*sigmax)).*exp(-((psxactual(i)-px(j)-mex).^2)/(2*sigmax^2))));
else    psxactualdist(j,i)=(0*((1/(sqrt(2*pi)*sigmax)).*exp(-((psxactual(i)-px(j)-mex).^2)/(2*sigmax^2))));   
end
psxactualpdf=trapz(psxactualdist);
    end
end

%y-direction
psyactual=-400:1:400;
py=-400:1:400;

for i=1:length(psyactual)
    for j=1:length(py)
if (vwy0>=(v1*sin(th1))&&vwy0<=(v1*sin(th2)));
    psyactualdist(j,i)=(((1/cp)*cvw*log((v1-sqrt((v1^2)-(((py(j)-cpy)/cp).^2)))./((py(j)-cpy)/cp))-cvwy1)*((1/(sqrt(2*pi)*sigmay)).*exp(-((psyactual(i)-py(j)-mey).^2)/(2*sigmay^2))));
elseif (vwy0>=(v2*sin(th1))&&vwy0<=(v2*sin(th2)));
    psyactualdist(j,i)=((1/cp)*(cvwy2-cvw*log((v2-sqrt((v2^2)-(((py(j)-cpy)/cp).^2)))./((py(j)-cpy)/cp))))*((1/(sqrt(2*pi)*sigmay)).*exp(-((psyactual(i)-py(j)-mey).^2)/(2*sigmay^2)));
elseif (vwy0>(v1*sin(th2))&&vwy0<(v2*sin(th1)));
        psyactualdist(j,i)=(((1/cp)*(cvwy2-cvwy1))*((1/(sqrt(2*pi)*sigmay)).*exp(-((psyactual(i)-py(j)-mey).^2)/(2*sigmay^2))));
else    psyactualdist(j,i)=(0*((1/(sqrt(2*pi)*sigmay)).*exp(-((psyactual(i)-py(j)-mey).^2)/(2*sigmay^2))));   
end
psyactualpdf=trapz(psyactualdist);
    end
end


    


































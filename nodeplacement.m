clear all
clc
%Calculation of Probability Density Function for Sensor Node Placement
initialcond
cvw=1/((th2-th1)*(v2-v1));
cvwx1=cvw*log((1-sqrt(1-(cos(th1))^2))/cos(th1));
cvwx2=cvw*log((1-sqrt(1-(cos(th2))^2))/cos(th2));
if (vwx0>(v1*cos(th2))&&vwx0<(v1*cos(th1)));
    vwx=v1*cos(th2):0.001:v2*cos(th2);
    vwxpdf=cvw*log((v1-sqrt((v1^2)-(vwx.^2)))./vwx)-cvwx2;
   
elseif (vwx0>(v1*cos(th1))&&vwx0<(v2*cos(th2)));
    vwx=v2*cos(th2):0.001:v1*cos(th1);
    vwxpdf=cvwx1-cvwx2;
   
elseif (vwx0>(v2*cos(th2))&&vwx0<(v2*cos(th1)));
    vwx=v1*cos(th1):0.001:v2*cos(th1);
    vwxpdf=cvwx1-cvw*log((v2-sqrt((v2^2)-(vwx.^2)))./vwx);
else vwxpdf=0;
end


cvwy1=cvw*log((1-sqrt(1-(sin(th1))^2))/sin(th1));
cvwy2=cvw*log((1-sqrt(1-(sin(th2))^2))/sin(th2));
if (vwy0>(v1*sin(th1))&&vwy0<(v1*sin(th2)));
    vwy=v1*sin(th1):0.001:v1*sin(th2)
    vwypdf=cvw*log((v1-sqrt((v1^2)-(vwy.^2)))./vwy)-cvwy1;
elseif (vwy0>(v1*sin(th2))&&vwy0<(v2*sin(th1)));
    vwy=v1*sin(th2):0.001:v2*sin(th1);
    vwypdf=cvwy2-cvwy1;
elseif (vwy0>(v2*sin(th1))&&vwy0<(v2*sin(th2)));
    vwy=v2*sin(th1):0.001:v2*sin(th2);
    vwypdf=cvwy2-cvw*log((v2-sqrt((v2^2)-(vwy.^2)))./vwy);
else vwypdf=0;
end


cvs=kw/kar;
cvsx=cvs*vwx0+vsx0;
cvsy=cvs*vwy0+vsy0;
cvsz=-((m*g*cvs)/kw)+vsz0;
T=-log((cvs*vwx0)/cvsx)*m/kar;
for t=0:T
    vsx=cvsx*exp(-(kar*t)/m)-cvs*vwx0;
    vsy=cvsy*exp(-(kar*t)/m)-cvs*vwy0;
    vsz=cvsz*exp(-(kar*t)/m)+m*g*cvs/kw;
end
cp=kw*T/kar+m*(kw/(kar^2))*(exp(-kar*T/m)-1);
cpx=ppx-m*(vsx0/kar)*(exp(-kar*T/m)-1);
cpy=ppy-m*(vsy0/kar)*(exp(-kar*T/m)-1);
cpz=ppz-m*(vsz0/kar)*(exp(-kar*T/m)-1);
psxfin=cp.*vwx+cpx;
psyfin=cp.*vwy+cpy;

%introduce psxfin1 instead of psxfin, to account for the mathematical model 
if (vwx0>(v1*cos(th2))&&vwx0<(v1*cos(th1)));
    psxfin1=(cpx+cp*v1*cos(th2)):0.0001:(cpx+cp*v2*cos(th1));
    psxfinpdf=(1/cp)*cvw*log((v1-sqrt((v1^2)-(((psxfin1-cpx)/cp).^2)))./((psxfin1-cpx)/cp))-cvwx2;
elseif (vwx0>(v2*cos(th2))&&vwx0<(v2*cos(th1)));
    psxfin1=(cpx+cp*v1*cos(th2)):0.0001:(cpx+cp*v2*cos(th1));
    psxfinpdf=(1/cp)*(cvwx1-cvw*log((v2-sqrt((v2^2)-(((psxfin1-cpx)/cp).^2)))./((psxfin1-cpx)/cp)));
end

if (vwy0>=(v1*sin(th1))&&vwy0<=(v1*sin(th2)));
    psyfin1=(cpy+cp*v1*sin(th1)):0.0001:(cpy+cp*v2*sin(th2));
    psyfinpdf=(1/cp)*cvw*log((v1-sqrt((v1^2)-(((psyfin1-cpy)/cp).^2)))./((psyfin1-cpy)/cp))-cvwy1;
elseif (vwy0>=(v2*sin(th1))&&vwy0<=(v2*sin(th2)));
    psyfin1=(cpy+cp*v1*sin(th1)):0.0001:(cpy+cp*v2*sin(th2));
    psyfinpdf=(1/cp)*(cvwy2-cvw*log((v2-sqrt((v2^2)-(((psyfin-cpy)/cp).^2)))./((psyfin-cpy)/cp)));
end

psxactual=3
syms px 

    pnrpxpdf=(1/(sqrt(2*pi)*sigmax)).*exp(-(((psxactual-px)-mex).^2)/(2*sigmax^2));
    psxfinpdf=((1/cp).*(cvwx1-cvw.*log((v2-sqrt((v2^2)-(((px-cpx)/cp).^2)))./((px-cpx)/cp))));
multipl=-(5571089318940243*((23158249616816542938140033119749*log(-((64 - (562949953421312/5460031289076775*px - 750852422226929408/5460031289076775)^2)^(1/2) - 8)/((562949953421312*px)/5460031289076775 - 750852422226929408/5460031289076775)))/1298074214633706907132624082305024 + 30983901282646560938906088050025/5192296858534827628530496329220096))/(36028797018963968*exp((1250*(psxactual - px)^2)/16641))


   sol=int(multipl,px,1,10)


    


































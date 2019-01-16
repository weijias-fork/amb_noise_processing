function [ierr,trig,ftrig]=trigger(grvel,om,nf,tresh)
% test dispersion curve for jumps

ierr    = 0;
ftrig   = zeros(nf,1);
trig    = zeros(nf,1);


for i =1:nf-2
    trig(i+1) = 0;
    hh1 = om(i+1)-om(i);
    hh2 = om(i+2)-om(i+1);
    hh3 = hh1+hh2;
    r = (grvel(i,1)/hh1-(1/hh1+1/hh2)*grvel(i+1,1)+ ...
        grvel(i+2,1)/hh2)*hh3/4*100;
    ftrig(i+1) = r;
    if abs(r)>tresh
        trig(i+1) = sign(r);
        ierr = 1;
    end
end
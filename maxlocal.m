function [maxloc,indmax]=maxlocal(y)

ind=1;
maxloc(1)=0;
indmax(1)=length(y);
for i=2:length(y)-1
    if y(i)>y(i-1) && y(i)>=y(i+1)
      maxloc(ind)=y(i);
      indmax(ind)=i;
      ind=ind+1;
    end
end
function ChrIF=DistConv(ChrIF,s,alpha,dmean)
%Distance conversion and update

for i=1:s
ChrIF(i,4)=1/(ChrIF(i,3)^alpha);
end

beta=dmean/mean(ChrIF(:,4));

ChrIF(:,4)=ChrIF(:,4)*beta;
 

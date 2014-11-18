function lnL_a=GS_alpha(ChrIF,s,weight,Ensemble,SegLen,omega,Dim,dmean)
%Grid search for alpha

%The range of parameter and step length should be chosen appropriately
alpha=[0.1:0.05:1.2];

loop=size(alpha);loop=loop(2);

lnL_a=zeros(loop,1);

for j=1:loop

ChrIF=DistConv(ChrIF,s,alpha(j),dmean);

lnL_a(j)=lnL_overall(weight,Ensemble,ChrIF,SegLen,omega,Dim);

end

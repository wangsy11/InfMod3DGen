function lnL_o=GS_omega(ChrIF,weight,Ensemble,SegLen,Dim)
%Grid search for omega

omega=[0.5:0.25:6];

loop=size(omega);loop=loop(2);

lnL_o=zeros(loop,1);

for j=1:loop

lnL_o(j)=lnL_overall(weight,Ensemble,ChrIF,SegLen,omega(j),Dim);

end

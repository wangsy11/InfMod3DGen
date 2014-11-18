function lnL_all=lnL_overall(weight,Ensemble,ChrIF,SegLen,omega,Dim)
%calculate overall likelihood (logarithm) for whole ensemble

vol=size(weight);vol=vol(1);

Ptemp=zeros(Dim,3);

lnL_all=0;

for i=1:vol
    if weight(i)~=0
        
        Ptemp(:,:)=Ensemble(i,:,:);
        lnL_all=lnL_all+weight(i)*likelihood(Ptemp,ChrIF,SegLen,omega);
        
    end
  
end

function D=RevDist(Ensemble,n,ChrIF)
%Computed value of distances between gene loci pairs given structure ensemble 

s=size(ChrIF);s=s(1);

D=zeros(n,s);

for i=1:n
for j=1:s
  
  D(i,j)=sqrt(sum((Ensemble(i,ChrIF(j,1),:) - Ensemble(i,ChrIF(j,2),:)) .^ 2)); 

end
end


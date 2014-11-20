%---------------------------------------------------------------------%
%Inferential modeling of 3D chromatin structure. Author: Siyu Wang, De-
%partment of Automation, Tsinghua University. Date:18-11-2014
%---------------------------------------------------------------------%

function ChrMod_main(ChrNo,n)

%Process the original data to a simple format of interaction frequency information that can be further processed
ChrIF=PreProcess(ChrNo);

%Parameters setting for starting points of the modeling
alpha=0.5;   %initial alpha
omega=4;     %initial value for omega

%Get the dimension of chromosome (number of model segments) and size of datasets 
Dim=max(ChrIF(:,2));s=size(ChrIF);s=s(1);

%Estimate the size of chromosome as the starting point for model iteration
dmean=350; 
ChrDat=DataInput(ChrNo);
ChrSize=(max(ChrDat(:,4))/12100000* 0.5236)^(1/3)*2000;

%Estimate the physical length of each segments
EstLen=SegLen(ChrNo);

%Number for total iteration steps, large enough to convergent
loop=25;

%Initial for record some values during modeling procedure
Correlation=zeros(loop,1);        %correlation between reuslts and input data for each iteration step
OverallLikelihood=zeros(loop,1);  %overall likelihood
AlphaRec=zeros(loop,1);           %value of parameter alpha
OmegaRec=zeros(loop,1);           %omega

%Random initialize a ensemble of structure as start of modeling 
Ensemble=zeros(n,Dim,3);

for i=1:n
 Ensemble(i,:,:)=StrucInitial(Dim,ChrSize);
end
%------------------------------------------------------------------------------------------------------------------

%Iteration for EM-like algorithm
for k=1:loop

%Distance conversion from interaction frequency
%updated before each iteration as value of parameter changed
ChrIF=DistConv(ChrIF,s,alpha,dmean);

%E step 
%Optimization for structures in ensemble
for i=1:n
  t_struc=zeros(Dim,3);
  t_struc(:,:)=Ensemble(i,:,:);
  Ensemble(i,:,:)=StrucOpt(ChrIF,t_struc,Dim,EstLen,omega);
end

%Get the corresponding weight for each structure in ensemble
weight=zeros(n,1);  
  for i=1:n
    t_struc=zeros(Dim,3);
    t_struc(:,:)=Ensemble(i,:,:);
    weight(i)=exp(likelihood(t_struc,ChrIF,EstLen,omega)/10000);%amplify so that won't be zero
  end    

%Nomarlization of weights
wsum=sum(weight);
for i=1:n
   weight(i)=weight(i)/wsum;
end

%Computed value of distances between gene loci pairs from structure ensemble 
D=RevDist(Ensemble,n,ChrIF);
D=D';
D_excp=D*weight;   %Averaged with weights
%Computed Overall likelihood and correaltion between original data and
%results
Correlation(k)=corr(D_excp,ChrIF(:,4));
OverallLikelihood(k)=lnL_overall(weight,Ensemble,ChrIF,EstLen,omega,Dim);

%--------------------------------------------
if k==loop  %End of modeling, save result
    name=ChrNo*10000+dmean;
    save(num2str(name));
    %vairable 'Ensemble' is the 3D structure we want
%--------------------------------------------

else    
    %M step
    %Grid search for alpha 
    lnL_al=GS_alpha(ChrIF,s,weight,Ensemble,EstLen,omega,Dim,dmean);
    a=find(lnL_al==max(lnL_al));
    alpha=0.05*a+0.05;%value of alpha updated
    AlphaRec(k)=alpha;
    
    %Expected 3D distance update
    ChrIF=DistConv(ChrIF,s,alpha,dmean);
    
    %Grid search for omega 
    lnL_om=GS_omega(ChrIF,weight,Ensemble,EstLen,Dim);
    omega=find(lnL_om==max(lnL_om));
    omega=0.25+0.25*omega;
    OmegaRec(k)=omega;
    
end

end







function ChrData=DataInput(j)
%This function is used for read input data file and process it for further
%usage.
%We use yeast genome Hi-C data from Duan et al. (Nature. 2010 May 20;465(7296):363-7. doi: 10.1038/nature08973)
Data=textread('interactions_HindIII_fdr0.0001_intra.txt');

%Get the length for data of each chromosome.
seg=Data(:,2);

totalsize=size(seg,1);

length=zeros(17,1);

length(1)=0;

f=2;

for i=2:totalsize-1
 if seg(i)<seg(i-1)
    length(f)=i-1;
    f=f+1;
 end
end

length(f)=totalsize;

%Pick up the data of the chromosome that we will model 
ChrData=Data(length(j)+1:length(j+1),:);

end

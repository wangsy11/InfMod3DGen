function IF=PreProcess(ChrNo)
%This function is to process the data into a simpler format for our covience for further usage.

ChrData=DataInput(ChrNo);

seg=[ChrData(:,2);ChrData(:,4)];

dim=size(ChrData);dim=dim(1);

seg=sort(seg);

t=0;f=0;

%Correspond position of enzyme cutting sites to number
for i=1:2*dim
 if seg(i)~=t
 f=f+1;
 rep(f,1)=f;
 rep(f,2)=seg(i);
 t=seg(i);
 end
end

%We use sequencial number to take the plcae of position of enzyme cutting
%sites in original dataset.
for i=1:dim
 for j=1:f
     if ChrData(i,2)==rep(j,2)
         ChrData(i,2)=rep(j,1);
     end
     if ChrData(i,4)==rep(j,2)
         ChrData(i,4)=rep(j,1);
     end  
 end
end

IF=[ChrData(:,2) ChrData(:,4) ChrData(:,5)];

end
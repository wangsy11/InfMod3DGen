function EstLen=SegLen(num)
%This function is for estimation of real length for each segments of
%chromosome

chro=DataInput(num);

seg=[chro(:,2);chro(:,4)];

dim=size(chro);dim=dim(1);

seg=sort(seg);

t=0;f=0;

for i=1:2*dim
 if seg(i)~=t
 f=f+1;
 SegSeqLen(f,1)=f;
 SegSeqLen(f,2)=seg(i)-t;
 t=seg(i);
 end
end

SegSeqLen=SegSeqLen(:,2);

%130bp/nm
EstLen=SegSeqLen/130;
function t_struc=StrucOpt(ChrIF,t_struc,Dim,EstLen,omega)
%This is for optimization for each single structure in ensemble

%Number for iteration steps,should be large enough
loop=Dim*5;

%Record the likelihood change of each step, later used for stop judgement 
q=zeros(loop+1,1);%
f=0;

%Starting likelihood value
q(1)=likelihood(t_struc,ChrIF,EstLen,omega);

%-------------------------------------------------------------------------
for i=1:loop
    
j=mod(i,Dim)+1;
t_s=t_struc;
lnlt=-1e20;xp=0;yp=0;zp=0;

%Search for finding directions that likelihood rise fastest 
 for x=-15:3:15
    t_s(j,1)=t_struc(j,1)+x;
    lnltp=l_search(t_s,ChrIF,j,EstLen,omega);
    if lnltp>lnlt
        lnlt=lnltp;
        xp=x;
    end
 end
 t_s(j,1)=t_struc(j,1)+xp;
 
  for y=-15:3:15
    t_s(j,2)=t_struc(j,2)+y;
    lnltp=l_search(t_s,ChrIF,j,EstLen,omega);
    if lnltp>lnlt
        lnlt=lnltp;
        yp=y;
    end
 end
 t_s(j,2)=t_struc(j,2)+yp;
 
  for z=-15:3:15
    t_s(j,3)=t_struc(j,3)+z;
    lnltp=l_search(t_s,ChrIF,j,EstLen,omega);
    if lnltp>lnlt
        lnlt=lnltp;
        zp=z;
    end
 end
 t_s(j,3)=t_struc(j,3)+zp;
%End of optimizing search


 t_struc=t_s;

 q(i+1)=likelihood(t_struc,ChrIF,EstLen,omega);
 
%Judge if already convergent 
if i>Dim+2
  if q(i+1)-q(i-Dim)<abs(0.0001*q(1))
     f=f+1;
  end
end
 
 if f==1
     break
 end
 
end






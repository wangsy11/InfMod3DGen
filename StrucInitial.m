function P=StrucInitial(Dim,chrosize)
%Initialization for a single structure in 3D space

P1=zeros(Dim,3);
step=chrosize/Dim*30;

%Randomly generated as a Brownian motion
 for i=1:Dim-1
 P1(i+1,:)=P1(i,:)+step*(rand(1,3)-[0.5,0.5,0.5]);   
 end
 
P=P1;
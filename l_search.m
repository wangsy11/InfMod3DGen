function lnLp=l_search(t_s,ChrIF,j,EstLen,omega)
%This function is similar to function 'likelihood' but we clear away some
%unnecessary steps for optimization search so that save time for
%computation

w=1*10^-(omega);


k1=2.3*10^-21;
k2=1.7*10^-20;
k3=4.1*10^-21;
sigma=30;

s=size(ChrIF);s=s(1);Dim=max(ChrIF(:,2));

sqsum=0;

for i=1:s
    if j==ChrIF(i,1)||j==ChrIF(i,2)
    dist=sqrt(sum((t_s(ChrIF(i,1),:) - t_s(ChrIF(i,2),:)) .^ 2));          
    sqsum=(dist-ChrIF(i,4))^2+sqsum;   
    end   
end

%Energy change, branch structure used here in order to deal with occasion at the end chromosome 
E1=0;E2=0;E3=0;
for i=j-1:j
    if i~=Dim&&i~=0        
      E1=0.5*k1*(sqrt(sum((t_s(i+1,:) - t_s(i,:)) .^ 2))-EstLen(i)).^ 2+E1;
    end
end

if j==1
  theta1=acos(dot((t_s(j+1,:)-t_s(j,:)),(t_s(j+2,:)-t_s(j+1,:)))/(norm(t_s(j+1,:)-t_s(j,:))*norm(t_s(j+1,:)-t_s(j+2,:))));
  theta2=0;theta3=0;
   elseif j==2
    theta1=acos(dot((t_s(j+1,:)-t_s(j,:)),(t_s(j+2,:)-t_s(j+1,:)))/(norm(t_s(j+1,:)-t_s(j,:))*norm(t_s(j+1,:)-t_s(j+2,:))));
    theta2=acos(dot((t_s(j,:)-t_s(j-1,:)),(t_s(j+1,:)-t_s(j,:)))/(norm(t_s(j,:)-t_s(j-1,:))*norm(t_s(j+1,:)-t_s(j,:))));
    theta3=0;
     elseif j==Dim-1
      theta2=acos(dot((t_s(j,:)-t_s(j-1,:)),(t_s(j+1,:)-t_s(j,:)))/(norm(t_s(j,:)-t_s(j-1,:))*norm(t_s(j+1,:)-t_s(j,:))));
      theta3=acos(dot((t_s(j-1,:)-t_s(j-2,:)),(t_s(j,:)-t_s(j-1,:)))/(norm(t_s(j-1,:)-t_s(j-2,:))*norm(t_s(j,:)-t_s(j-1,:))));   
      theta1=0;
       elseif j==Dim
        theta3=acos(dot((t_s(j-1,:)-t_s(j-2,:)),(t_s(j,:)-t_s(j-1,:)))/(norm(t_s(j-1,:)-t_s(j-2,:))*norm(t_s(j,:)-t_s(j-1,:))));       
        theta1=0;theta2=0;
         else
           theta1=acos(dot((t_s(j+1,:)-t_s(j,:)),(t_s(j+2,:)-t_s(j+1,:)))/(norm(t_s(j+1,:)-t_s(j,:))*norm(t_s(j+1,:)-t_s(j+2,:))));
           theta2=acos(dot((t_s(j,:)-t_s(j-1,:)),(t_s(j+1,:)-t_s(j,:)))/(norm(t_s(j,:)-t_s(j-1,:))*norm(t_s(j+1,:)-t_s(j,:))));
           theta3=acos(dot((t_s(j-1,:)-t_s(j-2,:)),(t_s(j,:)-t_s(j-1,:)))/(norm(t_s(j-1,:)-t_s(j-2,:))*norm(t_s(j,:)-t_s(j-1,:))));   
end
E2=E2+0.5*k2*(theta1*theta1+theta2*theta2+theta3*theta3);

for j=1:Dim
   if j~=j-1&&j~=j&&j~=j+1
     if sqrt(sum((t_s(j,:) - t_s(j,:)) .^ 2))<1.1225*sigma%
         d=sqrt(sum((t_s(j,:) - t_s(j,:)) .^ 2));
        E3=E3+4*k3*((sigma/d).^12-(sigma/d).^6+0.25); 
     end
   end
 end

E=E1+E2+E3;

i1=E/k3;

i2=w*sqsum;

lnLp=-i1-i2;








function lnL=likelihood(t_struc,ChrIF,EstLen,omega)
%Calculate the likelihood for a given structure

w=1*10^-(omega);

%parameters for energy determination
k1=2.3*10^-21;
k2=1.7*10^-20;
k3=4.1*10^-21;

%parameter for LJ formula
sigma=30;

s=size(ChrIF);s=s(1);Dim=max(ChrIF(:,2));

DataViolence=0;

i=1;

%Summation of data violence
for j=1:(Dim-1)
for k=j+1:Dim
    if j==ChrIF(i,1)&&k==ChrIF(i,2)
    dist=sqrt(sum((t_struc(ChrIF(i,1),:) - t_struc(ChrIF(i,2),:)) .^ 2));          
    DataViolence=(dist-ChrIF(i,4))^2+DataViolence;
    i=i+1;
     if i==s
         i=s-1;
     end   
    end   
end
end




%energy
E1=0;E2=0;E3=0;

%Energy of stretching
for i=1:(Dim-1)
E1=0.5*k1*(sqrt(sum((t_struc(i+1,:) - t_struc(i,:)) .^ 2))-EstLen(i)).^ 2+E1;
end

%Energy of bending
for j=1:(Dim-2)
theta=acos(dot((t_struc(j+1,:)-t_struc(j,:)),(t_struc(j+2,:)-t_struc(j+1,:)))/(norm(t_struc(j+1,:)-t_struc(j,:))*norm(t_struc(j+1,:)-t_struc(j+2,:))));
E2=E2+0.5*k2*theta*theta;
end

%Energy of excluding
for j=1:Dim
 for k=j+2:Dim; 
     if sqrt(sum((t_struc(j,:) - t_struc(k,:)) .^ 2))<1.1225*sigma
         d=sqrt(sum((t_struc(j,:) - t_struc(k,:)) .^ 2));
        E3=E3+4*k3*((sigma/d).^12-(sigma/d).^6+0.25); 
     end
 end
end
E=E1+E2+E3;


L1=E/k3;

L2=w*DataViolence;

L3=s*sqrt(1/(2*w));

%logarithm of total likelihood
lnL=-L1-L2-L3;







clear
clc
N=4;
g=1;
P=zeros(N,N);
s=ones(1,N)*100/N
SOld=zeros(N,N);
SNew=SOld;
T=24;
Iterations=[1:T];
MaxDif=zeros(1,T);
Nwght=0.2

for i=1:N
    for j=i+1:i+g
        if j<=N
%            P(i,j)=i/(j+i);
            P(i,j)=(1+rand())/2;
            P(j,i)=1-P(i,j);
        else
%            P(i,j-N)=i/(i+j-N);
            P(i,j-N)=rand()/2;
            P(j-N,i)=1-P(i,j-N);
        end
    end
end
P
sNew=ones(1,N);
for t=1:T
    for i=1:N
        for j=1:N
            if P(i,j)>0
                SNew(i,j)=(P(i,j)/P(j,i))*s(j);
                sNew(i)=sNew(i)*(P(i,j)/P(j,i))*s(j);
            end
            if j==N
                    sNew(i)=sNew(i)^(1/(2*g));
            end
        end
    end
    sNew=sNew*100/(sNew*ones(N,1));
    MaxDif(t)=max(abs(sNew-s));
    s=(1-Nwght)*s+Nwght*sNew
    sNew=ones(1,N);
end
MaxDif
plot(Iterations,MaxDif)
PP=P;
for i=1:N
    for j=1:N
        if P(i,j)>0
            PP(i,j)=s(i)/(s(i)+s(j));
        end
    end
end
P
PP
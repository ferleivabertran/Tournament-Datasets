Abar=32;
Alpha=zeros(Abar,1);
Alpha(1,1)=0.0001;
Alpha(2,1)=0.0002;
Alpha(3,1)=0.0005;
Alpha(4,1)=0.001;
Alpha(5,1)=0.002;
Alpha(6,1)=0.005;
Alpha(7,1)=0.01;
Alpha(8,1)=0.02;
Alpha(9,1)=0.05;
Alpha(10,1)=0.10;
Alpha(11,1)=0.15;
Alpha(12,1)=0.20;
Alpha(13,1)=0.25;
Alpha(14,1)=0.30;
Alpha(15,1)=0.32;
Alpha(16,1)=0.34;
Alpha(17,1)=0.36;
Alpha(18,1)=0.37;
Alpha(19,1)=0.38;
Alpha(20,1)=0.39;
Alpha(21,1)=0.40;
Alpha(22,1)=0.41;
Alpha(23,1)=0.42;
Alpha(24,1)=0.43;
Alpha(25,1)=0.44;
Alpha(26,1)=0.45;
Alpha(27,1)=0.46;
Alpha(28,1)=0.47;
Alpha(29,1)=0.48;
Alpha(30,1)=0.49;
Alpha(31,1)=0.5;
Alpha(32,1)=0.6;
Alphas=Alpha';
Types=8;
TypesMin=4;
TypesMax=10;
BigType=20;
GoBig=0;
TypesTot=TypesMax-TypesMin+1;
AlphaTest=zeros(TypesTot+1,Abar);
AlphaTest(1,:)=Alphas;
g=11;
wT=(g-1)/g;
% wT=10/11;


for Types=TypesMin:TypesMax
     if Types==TypesMax
         GoBig=1;
     end
    if GoBig==1
        Typesold=Types;
        Types=BigType;
    end
    T=floor((Types+1)/2);
    odd=mod(Types,2);
    VT=ones(T,Abar)/2;
    VTnew=ones(T,Abar)/2;
    iter=1;
    Maxiter=100;
    while iter<Maxiter
        for k=1:Abar
            for t=1:T
                if t==1
                    VTnew(t,k)=Alpha(k)+(1-Alpha(k))*VT(t+1,k);
                elseif t==T
                    if odd==1
                        VTnew(t,k)=0.5;
                    end
                    if odd==0
                        VTnew(t,k)=Alpha(k)*wT+(1-Alpha(k))*(VT(t-1,k)*(1-wT)+(1-VT(t,k))*wT);
                    end
                else
                    VTnew(t,k)=Alpha(k)*(g-1)/g+(1-Alpha(k))*(VT(t-1,k)/g+VT(t+1,k)*(g-1)/g);
                end
            end
        end
        iterold=iter;
        iter=iter+1;
        MaxDiff=max(max(abs(VT-VTnew)));
        if MaxDiff<0.0000001
            iter=Maxiter;
        end
        VTold=VT;
        VT=VTnew;
    end
    if Types>TypesMax
        AlphaTest(Typesold-TypesMin+2,:)=2*(VT(1,:).*(ones(1,Abar)-Alphas));
    else
        AlphaTest(Types-TypesMin+2,:)=2*(VT(1,:).*(ones(1,Abar)-Alphas));
    end
end
AlphaTest;
figure
plot(AlphaTest(1,:),AlphaTest(2,:),AlphaTest(1,:),AlphaTest(3,:),AlphaTest(1,:),AlphaTest(4,:),AlphaTest(1,:),AlphaTest(5,:),AlphaTest(1,:),AlphaTest(6,:),AlphaTest(1,:),AlphaTest(7,:),AlphaTest(1,:),AlphaTest(8,:),AlphaTest(1,:),ones(1,Abar),ones(1,Abar)/2,[0:Abar-1]/(Abar-1))
title('Fairness criterion 4: Lower bound for minimum alpha')
xlabel('alpha')
ylabel('2 (1-alpha) v1')
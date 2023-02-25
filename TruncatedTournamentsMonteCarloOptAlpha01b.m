% SIMULATION PROGRAM for Scoring System for Incomplete Tournaments
% Difference with non-b version is that here we fix Str for all rep
% instead of changing every time.

clear
clc

%PARAMETERS AND OTHER FIXED COMPONENTS

%Conf=2;      % Total conferences - Select only even numbers
%Teams=4;    % Total teams in conference - Select only even numbers
%N=Conf*Teams; 
N=130; %Number of teams
repmin = -2; %If set to 1, no optimization of Alphas occurs
OptWght = 0.63; %Weight that new min/max threshold Alphas have in next iteration
reptot=20; %Number of times this program runs (diff Str and W each time)
correlatedsims=1; %If =0, then uniformly random simulations
                 %If =1, then correlated simulations
strordered=1;     %If =0 then strength is randomly distributed
                 %If =1 then teams are ordered by strength
%Thus, for both =1, correlated simultations are going to match strong teams                 
%with other strong teams and weak teams with other weak teams.

%Create Team Strengths using Pareto Distribution with curvature 'beta'
%Then subtract (1-gamma) so that the min is not 1 but gamma instead
beta=1.8; % Benchmark = 1.8
gamma=0.01;

% Create vector of std devs of the win percentage in the complete tournament

sigmaW = zeros(1,reptot);
% First, create the Alphas to consider
Abar=11; %Number of Alphas to consider. Each one for a diff scoring method
Alow=4;                    %Number of decimals for lowest alpha
%Note: This creates alphas that increase tenfold each time until we 
%reach 'Alow/10'. Then it increases in increments of '1/10' until 
% we reach 1, which represents win percentage.
OverrideAlpha=1;
%Override the Alphas to fit what I want;
AlphaAux=zeros(1,Abar);
    AlphaAux(1,1)= 0.0001  ;
    AlphaAux(2,1)= 0.001  ;
    AlphaAux(3,1)= 0.01  ;
    AlphaAux(4,1)= 0.10  ;
    AlphaAux(5,1)= 0.20  ;
    AlphaAux(6,1)= 0.30  ;
    AlphaAux(7,1)= 0.50  ;
    AlphaAux(8,1)= 0.70  ;
    AlphaAux(9,1)= 0.90  ;
    AlphaAux(10,1)=0.95  ;
    AlphaAux(11,1)=0.9999  ;
    
 %Comment out after use   
%     AlphaAux(1,1)= 0.03  ;
%     AlphaAux(2,1)= 0.034  ;
%     AlphaAux(3,1)= 0.038 ;
%     AlphaAux(4,1)= 0.042  ;
%     AlphaAux(5,1)= 0.046  ;
%     AlphaAux(6,1)= 0.05  ;
%     AlphaAux(7,1)= 0.054  ;
%     AlphaAux(8,1)= 0.058 ;
%     AlphaAux(9,1)= 0.062  ;
%     AlphaAux(10,1)=0.066 ;
%     AlphaAux(11,1)=0.07  ;

DiffThreshold=0.000001; %Convergence threshold for scoring method

onegamemax=1; %If =1, then 2 teams can play each other only once at most
sbar=100; %Number of Monte Carlo simulations (fixed Str, W; diff WT)
%Generate Incomplete Tournament with average 'gbar' games per team:
gbar=11;      %Average number of games played per team
g=ones(N,1);  %Number of Games each team plays
g=gbar*g;     %Assumes all teams play the same number of games
gmax=max(g);  %Equals gbar if all play the same number of games

rho = 0;
cwght = 4 * rho ;  %Correlation Weight: 0=Uniform;1=twice the neighbor;2=4times 
          %Note: At cwght=4 you are almost surely playing the closest neighbors

AggAccuracy=zeros(reptot+1,Abar); %These two are the final reports on the 
AggSDev=zeros(reptot+1,Abar); %Accuracy and SDev of each scoring method
VTNDifSDevCSAvRepTot = ones(reptot+1,Abar); %Output for Optimal Alpha

GTdistMeanMC = zeros(reptot,sbar); %Average distance between teams that play each other
GTdistMean = zeros(1,reptot); %Average distance between teams that play each other
GTdistSDev = zeros(1,reptot); %SDev of distance between teams that play each other

%Create Complete Tournament Fixture (where a row is a round of the tournament)

Local=zeros(N-1,N/2); 
Visita = zeros(N-1,N/2);
k=1;

for i=1:N-1
    Ind=(-1)^i;
    for j=1:N/2
        if j<N/2
            Local(i,j)=k;
        elseif Ind==1
            Visita(i,j)=k;
            Local(i,j)=N;
        else
            Local(i,j)=k;
            Visita(i,j)=N;
        end
        if k==N-1
            k=1;
        else
            k=k+1;
        end
    end
end

k=N-1;
for i=1:N-1
    for j=1:N/2-1
      Visita(i,j)=k;
      if k==1
          k=N-1;
      else
          k=k-1;
      end
    end
end

tic

%Create Team Strengths using Pareto Distribution with curvature 'beta'
%Then subtract (1-gamma) so that the min is not 1 but gamma instead
Str=zeros(N,1);
rdmv=rand(N,1);
AvgTopTeamStr=(1/(1-(2*N-1)/(2*N)))^(1/beta);
AvgNextTeamStr=(1/(1-(2*N-3)/(2*N)))^(1/beta);
for i=1:N
    Str(i)=(1/(1-rdmv(i)))^(1/beta)-(1-gamma);    
end
if strordered==1           %Teams numbers are assigned by strength
    Str=sortrows(Str,-1);  %from highest to lowest
end

% Placing the rep-loop here will guarantee that the Str vector does not change
for rep = repmin:reptot

%Create Win matrix of complete tournament (one round-robin)

W=zeros(N,N);            %Win Matrix
for i=1:N-1
    for j=i+1:N
        rdmb=rand();
        if rdmb<Str(i)/(Str(i)+Str(j))
            W(i,j)=1;
        else
            W(j,i)=1;
        end
    end
end

Wins=W*ones(N,1);               %Vector of total wins by each team
GamesPlayed=W.'*ones(N,1)+Wins; %Vector of total Games played by each team
            
%Order Teams by wins

[WWWins,index]=sortrows([W,Wins,Str,GamesPlayed],[-(N+1) -(N+2)]);%Order...
WWins=WWWins; %... all rows by wins (and then by strength)
for j=1:N                            %Orders first N columns by wins
    WWins(:,j)=WWWins(:,index(j));   %so that wins matrix is ordered
end                                  %by wins on both rows and columns
W=WWins(:,1:N);               
Wins=WWins(:,N+1);            %Now W, Wins, Losses, 
Losses=W'*ones(N,1);          %StrAux and GamesPlayed 
StrAux=WWins(:,N+2);          %are all ordered
GamesPlayed=WWins(:,N+3);     %by number of wins

%Obtain true scores: total wins and win-percentage of complete tournament

Winperc=zeros(N,1);
for i=1:N
    Winperc(i)=Wins(i)/GamesPlayed(i);
end

if rep>0
    sigmaW(rep) = std(Winperc);
end
%OUTPUT: True scores vs. Strength scores Comparison (Not necessary)

StrNorm=N*StrAux/(2*(StrAux.'*ones(N,1))); %normalize strengths so that avg str = 1/2
%StrNorm.'*ones(N,1)/N;

%OUTPUT: Obtain total points and scores under scoring system for multiple
%parameter values (alpha) of the scoring system.

% First, create the Alphas to consider

Alpha=zeros(Abar,1);       %Alphas to be created
%Note: This creates alphas that increase tenfold each time until we 
%reach 'Alow/10'. Then it increases in increments of '1/10' until 
% we reach 1, which represents win percentage.
for i=1:Alow
    Alpha(i)=Alow*10^(-Alow+(i-2));
end

for i=Alow+1:Abar
    Alpha(i)=(i-1)/(Abar-1);
end

%Create a better set of alphas to search for the optimal one
if rep<2
   if rep>repmin
      if ii-2>0
         AlphaMin = OptWght*Alphas(ii-2) + (1-OptWght)*AlphaMinOld
      elseif ii-1>0
         AlphaMin = OptWght*Alphas(ii-1) + (1-OptWght)*AlphaMinOld
      else
         AlphaMin = Alphas(ii)/2
      end
      if ii+2<=Abar
         AlphaMax = OptWght*Alphas(ii+2) + (1-OptWght)*AlphaMaxOld
      else
         AlphaMax = (Alphas(Abar)+Alphas(ii))/2
      end
      for i=1:Abar
         AlphaAux(i) = AlphaMin + (AlphaMax-AlphaMin)*(i-1)/(Abar-1); 
      end
   end
end

if OverrideAlpha==1
   'Override the Alphas to fit what I want';
   i=1;
   while i<=Abar 
     Alpha(i,1)= AlphaAux(i,1);
     i=i+1;
   end
end

if rep<2
    AlphaMinOld = Alpha(1);
    AlphaMaxOld = Alpha(Abar);
end
% Then create the Points and Scores matrices
P=ones(N,Abar)/2;  %Points Matrix 
Pold=P;
V=P;               %Scores Matrix (initially 1/2 for all i,alpha)
Vold=V;
iter=1;                 %Counts iterations required for convergence
done=0;                 %Turns into 1 when convergence is met
while done==0
  for k=1:Abar      %For all alphas
    for i=1:N       %For all teams
        sum=0;      
        for j=1:N   %For all opponents
            sum=sum+(W(i,j)+W(j,i))*Vold(j,k);
        end
        P(i,k)=Alpha(k)*Wins(i)+(1-Alpha(k))*sum; %SCORING METHOD
        V(i,k)=P(i,k)/GamesPlayed(i);             %AT WORK
    end
  end
  DiffV=abs(V-Vold);
  MaxDiff=max(max(DiffV));
  iter=iter+1;
  if MaxDiff<DiffThreshold
      done=1;
  end
  Pold=P;
  Vold=V;
  vbar=ones(1,N)*Vold/N;
end
iter;    %Number of iterations until convergence

%Control that normalization is working
 VN=V;
 PN=P;
 for j=1:Abar
   VN(:,j)=(V(:,j)-(1/2)*((1-Alpha(j))*N/(N-Alpha(j)))*ones(N,1))*((N-Alpha(j))/((N-1)*Alpha(j)));
   PN(:,j)=(P(:,j)-(1/2)*((1-Alpha(j))*N*(N-1)/(N-Alpha(j)))*ones(N,1))*((N-Alpha(j))/((N-1)*Alpha(j)));
 end

Vvector=V(:);  %Stacks all columns of V into one vector 
Pvector=P(:);  %Stacks all columns of P into one vector

%BEGIN ITERATION FOR MONTE CARLO SIMULATION
 
mciter=0;    %NOTE: Two teams may play twice but it would be very rare
VTmatrix=zeros(length(Vvector),sbar); %Matrix of all simulation scores
PTmatrix=zeros(length(Pvector),sbar); %Matrix of all simulation points
VTNmatrix=zeros(length(Vvector),sbar); %Matrix of all normalized simulation scores
PTNmatrix=zeros(length(Pvector),sbar); %Matrix of all normalized simulation point

Alphas=Alpha.';
SDevCS=zeros(sbar+1,Abar+1);
SDevCSN=SDevCS;
SDevCS(1,2:Abar+1)=Alphas;
SDevCSN(1,2:Abar+1)=Alphas;


%Generate Incomplete Tournament with average 'gbar' games per team:
if onegamemax==1
    if gmax>N-1
        'Too many games'
    end
    doubles=zeros(1,sbar);
end

%OUTPUT TO WRITE BEFORE MONTE CARLO ON FIRST REP (rep=1)
if rep==1
 N;
 %Local
 %Visita
 AvgTopTeamStr;
 AvgNextTeamStr;
 minN=min(N,10);
 StrWinLoss=[StrAux(1:minN,:),Wins(1:minN,:),Losses(1:minN,:),StrAux(N+1-minN:N,:),Wins(N+1-minN:N,:),Losses(N+1-minN:N,:)];
 %Wins
 %GamesPlayed
 %Winperc
 %Abar
 %P
 %V
 sbar;
 gbar;
end



%Generate Incomplete Tournament with correlated games matrix:
if correlatedsims==1
  FixedCorrelMatrix=zeros(N,N);  % The weight on the correlation matrix is given
  for i=1:N                 % by the distance between the teams (=abs(j-i))
    for j=1:N             % to a fixed power (=cwght). Thus, if cwght = 0
        if i==j           % then we have a uniform distribution and if 
        else              % cwght->inf we have a 50% for each neighbor.
            FixedCorrelMatrix(i,j)=abs(j-i)^(-cwght);
        end
    end
  end
  SumCorrel=FixedCorrelMatrix*ones(N,1);
  for i=1:N
    FixedCorrelMatrix(i,:)=FixedCorrelMatrix(i,:)/SumCorrel(i);
  end 
end

%START MONTE CARLO ITERATION
for s=1:sbar
    
% 1.UNIFORMLY RANDOM GAME ASSIGNMENT
% Randomly assign teams to fixture table teams
if correlatedsims==0
  TL=(1:N).';
  Teamlist=TL(randperm(length(TL)));

% Randomly select "g" dates of the fixture table
  GL=(1:N-1).';
  Gamelistfull=GL(randperm(length(GL)));
  Gamelist=Gamelistfull(1:gmax);

% Use team assignments and fixture dates to create games played matrix 
% G, where g_ij = g_ji = 1 if a game is played and zero otherwise.
  GT=zeros(N,N);
  for i=1:gmax             %Run through all fixture rounds
    for j=1:N/2          %Run through all games within fixture round 
        loc=Teamlist(Local(i,j));
        vis=Teamlist(Visita(i,j));
        GT(loc,vis)=1;
        GT(vis,loc)=1;
    end
  end
  
% 2. CORRELATED GAME ASSIGNMENT
elseif correlatedsims==1
CorrelMatrix=FixedCorrelMatrix;
GT=zeros(N,N);
Nord=(1:N).';
NorderAux=Nord(randperm(length(Nord))); %Order of teams to randomize from
Norder=[NorderAux,Nord]; %Vector of team numbers (Col1) and their index (Col2)
Pending=g; %Pending uses the indexes of teams, not the team numbers
gtot=ones(1,N)*g/2;
  for k=1:gtot
    [maxp, index] = max(Pending);
    i=Norder(index,1); %This is the first team to be assigned a game
    r=rand();
    sum=0;
    j=0;     %This will be the opponent team to be assigned a game
    while r>sum
        sum=sum+CorrelMatrix(i,j+1);
        j=j+1;
        if j==N        %No oppontent has been found, so this will generate
            sum=1.1;   %the need to assign the team an opponent it has 
            for h=1:N  %already played against. This should not occur often
                if Pending(h)>0
                    if h==index
                    else
                        jindex=h;
                    end
                end
            end
            if onegamemax==1
                doubles(s)=doubles(s)+1;
            end
            j=Norder(jindex,1);
        end
    end      %At this point opponent j has been found
    if rep>0
        GTdistMeanMC(rep,s) = GTdistMeanMC(rep,s) + abs(j-i);
    end
    NorderAux=sortrows(Norder);
    jindex=NorderAux(j,2); %The index for j (in 'Pending') has been found
    if onegamemax==1
        CorrelMatrix(i,j)=0;
        CorrelMatrix(j,i)=0;
        SumCorrel=CorrelMatrix*ones(N,1);
        for h=1:N
            if SumCorrel(h)>0
                CorrelMatrix(h,:)=CorrelMatrix(h,:)/SumCorrel(h);
            end
        end
    end
    Pending(index)=Pending(index)-1; %Adjusts pending games (by index of i)
    Pending(jindex)=Pending(jindex)-1;%Adjusts pending games (by index of j)
    if Pending(index)==0
        CorrelMatrix(:,i)=zeros(N,1);
        SumCorrel=CorrelMatrix*ones(N,1);
        for h=1:N
            if SumCorrel(h)>0
                CorrelMatrix(h,:)=CorrelMatrix(h,:)/SumCorrel(h);
            end
        end
    end
    if Pending(jindex)==0
        CorrelMatrix(:,j)=zeros(N,1);
        SumCorrel=CorrelMatrix*ones(N,1);
        for h=1:N
            if SumCorrel(h)>0
                CorrelMatrix(h,:)=CorrelMatrix(h,:)/SumCorrel(h);
            end
        end
    end
    GT(i,j)=GT(i,j)+1;
    GT(j,i)=GT(j,i)+1;
  end
  if rep>0
      GTdistMeanMC(rep,s) = GTdistMeanMC(rep,s)/gtot;
  end
 end
% Generate truncated win/loss matrix that includes only games played
% according to games played matrix (by multplying g_ij * w_ij for all ij)

  WT=zeros(N,N);     %Matrix of wins in truncated tournament
  for i=1:N
    for j=1:N
        WT(i,j)=GT(i,j)*W(i,j);
    end
  end

% Create Win Total and Win Percentage vector

  WinsT=WT*ones(N,1); %Vector of total wins by each team
  GamesPlayedT=WT.'*ones(N,1)+WinsT; %Vector of total Games played by each team
  WinpercT=zeros(N,1);
  for i=1:N
    WinpercT(i)=WinsT(i)/GamesPlayedT(i);
  end



%OUTPUT: Obtain total points and scores under scoring system for different
%parameter values (alpha).

 PT=ones(N,Abar)/2;  %Points Matrix
 PoldT=PT;
 VT=PT;              %Scores Matrix
 VoldT=VT;
 iter=1;
 done=0;
 while done==0
  for k=1:Abar
    for i=1:N
        sum=0;
        for j=1:N
            sum=sum+(WT(i,j)+WT(j,i))*VoldT(j,k);
        end
        PT(i,k)=Alpha(k)*WinsT(i)+(1-Alpha(k))*sum;
        VT(i,k)=PT(i,k)/GamesPlayedT(i);
    end
  end
  DiffV=abs(VT-VoldT);
  MaxDiff=max(max(DiffV));
  iter=iter+1;
  if MaxDiff<DiffThreshold
      done=1;
  end
  PoldT=PT;
  VoldT=VT;
  vbarT=ones(1,N)*VoldT/N;
 end
 mciter=(s-1)*mciter/s+iter/s;

 VTN=VT;
 PTN=PT;
 for j=1:Abar
   VTN(:,j)=(VT(:,j)-(1/2)*((1-Alpha(j))*N/(N-Alpha(j)))*ones(N,1))*((N-Alpha(j))/((N-1)*Alpha(j)));
   PTN(:,j)=(PT(:,j)-(1/2)*((1-Alpha(j))*N/(N-Alpha(j)))*GamesPlayedT)*((N-Alpha(j))/((N-1)*Alpha(j)));
 end
 
 SDevCS(s+1,:)=std([WinpercT,VT]);
 SDevCSN(s+1,:)=std([WinpercT,VTN]);
 VTvector=VT(:);
 PTvector=PT(:);
 VTNvector=VTN(:);
 PTNvector=PTN(:);
 VTmatrix(:,s)=VTvector;
 PTmatrix(:,s)=PTvector;
 
 VTNmatrix(:,s)=VTNvector;
 PTNmatrix(:,s)=PTNvector;
 
 % OUTPUT TO WRITE DURING MONTE CARLO
 %Teamlist
 %Gamelist
 %WT
 %WinsT
 %GamesPlayedT
 %WinpercT
 %PT
 %VT
 ss=floor(s/10);
 ssbar=floor(sbar/10);
 for h=1:ssbar
     if s==10*h
         s;
     end
 end
end            % End of Iteration over all 's'

if rep>0
    GTdistMean(1,rep) = mean(GTdistMeanMC(rep,:));
    GTdistSDev(1,rep) = std(GTdistMeanMC(rep,:));
end
  
mciter;
if onegamemax==1
    AvgDoubles=doubles*ones(sbar,1)/sbar;
end

%OUTPUT: Calculate Mean and Standard Deviation of recorded points and
%scores 

VTmatrixmean=VTmatrix*ones(sbar,1)/sbar;
VTNmatrixmean=VTNmatrix*ones(sbar,1)/sbar;
VTmean=reshape(VTmatrixmean,[N,Abar]);
VTNmean=reshape(VTNmatrixmean,[N,Abar]);
VTmatrixstd=std(VTmatrix,0,2);
VTNmatrixstd=std(VTNmatrix,0,2);
VTstd=reshape(VTmatrixstd,[N,Abar]);
VTNstd=reshape(VTNmatrixstd,[N,Abar]);

%PLOT: Q1- Comparison of Mean of Recorded Score vs True Score for all teams

MeanValuesbyAlphaVsTrueWinperc=[9.9999,Alpha.';Winperc,VTNmean];
WinpercAux=repmat(Winperc,1,Abar);
MeanValuesbyAlphaMinusTrueWinperc=[9.9999,Alpha.';Winperc,VTNmean-WinpercAux];

%      Q2- Standard Deviation across teams (ordered by Score) 

SDevValuesbyAlphaVsWinperc=[9.9999,Alpha.';Winperc,VTstd];
SDevNValuesbyAlphaVsWinperc=[9.9999,Alpha.';Winperc,VTNstd];

%      Q3- Recorded Standard Deviation as a function of parameter alpha

%Also, if we iterate over 'g' we can assess the evolution of the Standard
%Deviation as a function of 'g',
%for fixed values of parameter alpha.

% How close is each Alpha-method to the true values?
WinPercMatrix=zeros(N,Abar);
for k=1:Abar
    WinPercMatrix(:,k)=Winperc;
end


AggAccuracy(1,:)=Alpha.';
if rep>0
  AggAccuracy(rep+1,:)=ones(1,N)*(VTNmean-WinPercMatrix)/N;
end

AggSDevAux=zeros(N,Abar);
for i=1:N
    for k=1:Abar
        AggSDevAux(i,k)=(VTNmean(i,k)-Winperc(i))^2;
    end
end

AggSDev(1,:)=Alpha.';
if rep>0
AggSDev(rep+1,:)=sqrt(ones(1,N)*AggSDevAux/N);
end
if rep==1
    toc
    if correlatedsims==1
     cwght;
    end
end

% Create Output for Optimal Alpha
WinPercRepMatrix = repmat(Winperc,Abar,sbar);
VTNmatrixDifSq = [VTNmatrix - WinPercRepMatrix].^2;
VTNmatrixDifSDev = (VTNmatrixDifSq * ones(sbar,1) / (max(1,sbar-1))).^(1/2);
VTNDifSDev = reshape(VTNmatrixDifSDev,[N,Abar]);
VTNDifSDevCSAv = ones(1,N) * VTNDifSDev / N;
VTNDifSDevCSAvRepTot(1,:) = Alphas;
if rep>0
  VTNDifSDevCSAvRepTot(rep+1,:) = VTNDifSDevCSAv;
end
if rep<1
[MinSDev,ii] = min(VTNDifSDevCSAv);
ii
OptimalAlphaGuess = Alphas(ii)
end

end    % End of Iteration over all 'rep'
toc

WinpercSDev = std(Winperc);
SDevCS;
SDevCSN;
AggAccuracy;        
AggSDev;
VTNDifSDevCSAvRepTot;
figure
plot(Alphas,VTNDifSDevCSAvRepTot(2:reptot+1,:))
txt = ['N = ' num2str(N) ',  g = ' num2str(mean(g)) ',  rho = ' num2str(rho) ',  \sigma_{W} = ' num2str(mean(sigmaW))];
xloc = Alphas(1);
yloc = max(max(VTNDifSDevCSAvRepTot(2:reptot+1,:)));
text(xloc,yloc,txt)

[MinSDevs,I] = min(VTNDifSDevCSAvRepTot(2:reptot+1,:).');
OptimalAlphas = zeros(1,reptot);
for repp = 1:reptot
    OptimalAlphas(1,repp) = Alphas(I(repp));
end

GTdistMean
GTdistSDev
[mean(OptimalAlphas),median(OptimalAlphas),std(OptimalAlphas)]


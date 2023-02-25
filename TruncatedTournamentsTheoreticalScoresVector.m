% Generate truncated win/loss matrix that includes only games played
% according to games played matrix 

WTexp=zeros(N,N);            %Expected Win Matrix
for i=1:N
    for j=1:N
        WTexp(i,j)=Str(i)*GT(i,j)/(Str(i)+Str(j));
    end
end


%   WT=zeros(N,N);     %Matrix of wins in truncated tournament
%   for i=1:N
%     for j=1:N
%         WT(i,j)=GT(i,j)*W(i,j);
%     end
%   end

% Create Win Total and Win Percentage vector

WinsTexp=WTexp*ones(N,1); %Vector of total wins by each team
GamesPlayedTexp=WTexp.'*ones(N,1)+WinsTexp; %Vector of total Games played by each team
WinpercTexp=zeros(N,1);
 for i=1:N
   WinpercTexp(i)=WinsTexp(i)/GamesPlayedTexp(i);
 end


%OUTPUT: Obtain total points and scores under scoring system for different
%parameter values (alpha).
 tic
 PTexp=ones(N,Abar)/2;  %Points Matrix
 PoldTexp=PTexp;
 VTexp=PTexp;              %Scores Matrix
 VoldTexp=VTexp;
 iter=1;
 done=0;
 while done==0
  for k=1:Abar
    for i=1:N
        sum=0;
        for j=1:N
            sum=sum+(WTexp(i,j)+WTexp(j,i))*VoldTexp(j,k);
        end
        PTexp(i,k)=Alpha(k)*WinsTexp(i)+(1-Alpha(k))*sum;
        VTexp(i,k)=PTexp(i,k)/GamesPlayedTexp(i);
    end
  end
  DiffV=abs(VTexp-VoldTexp);
  MaxDiff=max(max(DiffV));
  iter=iter+1;
  if MaxDiff<DiffThreshold
      done=1;
  end
  PoldTexp=PTexp;
  VoldTexp=VTexp;
  vbarTexp=ones(1,N)*VoldTexp/N;
 end
  toc
 % mciter=(s-1)*mciter/s+iter/s;
 % Generate Output:
 VToldexp=VTexp;
 PToldexp=PTexp;
 VTNoldexp=VTexp;
 PTNoldexp=PTexp;
 for j=1:Abar
   VTNoldexp(:,j)=(VTexp(:,j)-(1/2)*((1-Alpha(j))*N/(N-Alpha(j)))*ones(N,1))*((N-Alpha(j))/((N-1)*Alpha(j)));
   PTNoldexp(:,j)=(PTexp(:,j)-(1/2)*((1-Alpha(j))*N/(N-Alpha(j)))*GamesPlayedTexp)*((N-Alpha(j))/((N-1)*Alpha(j)));
 end
 
% SDevCS(s+1,:)=std([WinpercT,VT]);
% SDevCSN(s+1,:)=std([WinpercT,VTN]);
 VToldvectorexp=VToldexp(:);
 PToldvectorexp=PToldexp(:);
 VTNoldvectorexp=VTNoldexp(:);
 PTNoldvectorexp=PTNoldexp(:);
% VToldmatrix(:,rep)=VTvector;
% PToldmatrix(:,rep)=PTvector;
% VTNoldmatrix(:,rep)=VTNoldvector;
% PTNoldmatrix(:,rep)=PTNvector;

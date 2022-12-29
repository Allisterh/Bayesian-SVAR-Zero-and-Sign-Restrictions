%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
path(path,'U:\LouisRonchi\EmpMacroPaper')
path(path,'U:\Luca\matprog\')
path(path,'U:\Luca\matprog\stats')
%options = optimset;
%options = optimset(options,'Display','off');
%warning('off','all')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
HorizonForSignRestrictions=2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Intercept='Y';
LagOrder=4;
NumberOfDraws=40000000;
Horizon=5*12;
Stationary='N';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                       Loading time series:
X=xlsread('U:\LouisRonchi\EmpMacroPaper\SwissData.xlsx','Data','A36:H243');
%X=xlsread('/Users/louis/Desktop/ALL/programming/matlab/EmpMacroPaper/EmpMacroPaper/SwissData.xlsx','Data','A36:H243');
[T,N]=size(X);
Time = X(:,1);
LogGDPPerCapita=log(X(:,2)./ X(:,6));
LogConsumptionPerCapita=log(X(:,5)./X(:,6));
LogGDPDeflator=log(X(:,3));
ShortRate=X(:,4)/100;
UnemploymentRate=X(:,7)/100;
Y=[LogConsumptionPerCapita LogGDPPerCapita LogGDPDeflator UnemploymentRate ShortRate];
[T,N]=size(Y);
%
PositionConsumption=1;
PositionGDP=2;
PositionGDPDeflator=3;
PositionUnemployment=4;
PositionShortRate=5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NumberOfCointegrationVectors=1;
BurnInNN=5000; % Number of burn-in iterations
SamplingInterval=25;
ErgodicNN=10000*SamplingInterval; % Number of iterations of the ergodic distribution
NumberOfPermanentGDPShocks=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DFILE=['U:\LouisRonchi\EmpMacroPaper\PosteriorDistributionReducedFormVECM'];
varname(1,:)=['SIGMA'];
varname(2,:)=['ALPHA'];
varname(3,:)=['BETAA'];
varname(4,:)=['BBBBB'];
% Y=Y*100;
SavedY=Y;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Choosing the lag order:
X=diff(Y);
MaxLagOrder=8;
% Estimating the VAR:
LagOrderSIC=varplagselect(Y,'Y',MaxLagOrder,'SIC',1);
LagOrderHannanQuinn=varplagselect(Y,'Y',MaxLagOrder,'HQQ',1);
LagOrder=max([LagOrderSIC LagOrderHannanQuinn 2]')-1; % The lag order has to be at least 2 in levels
LagOrderSaved=LagOrder;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                   The estimated VECM:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %  The stationary VAR we will bootstrap in order to get the distribution of the trace statistic under the null of no cointegration:
% The cointegration vector for Y=[LogRealGDP LogRealConsumption];
CV=GetStockWatson1993DOLS(Y(:,PositionGDP),Y(:,PositionConsumption),0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %  The stationary VAR we will bootstrap in order to get the distribution of the trace statistic under the null of no cointegration:
% The matrix of the cointegration vectors:
CointegrationVectors=zeros(N,NumberOfCointegrationVectors);
CointegrationVectors(PositionGDP,1)=1;
CointegrationVectors(PositionConsumption,1)=CV(2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %  The stationary VAR we will bootstrap in order to get the distribution of the trace statistic under the null of no cointegration:
% We impose the cointegration vectors and we estimate the other parameters of the VECM via OLS
%
% This is the structure we need to estimate (from the code BootstrapCointegratedVARP.m):
% xx(:,tt)=Alpha+CSI*vec(fliplr(xx(:,tt-LagOrder:tt-1)))+CSI0*yy(:,tt-1)+BootstrappedE(tt,:)';
% that is: we regress X=diff(Y) on an intercept, lags of itself, and the
% level of Y lagged one period.
[Alpha,A,B,U,Omega,CSI,CSI0,~,Y,Z,LevelsLagged1Period]=EstimateVECMConditionalOnCointegrationVector(Y,CointegrationVectors,LagOrder);
%
MU=Alpha;
BB=CSI;
GG=CSI0;
UU=U';
VAR=Omega;
%
CointegrationVectors=A;
Loadings=B;
%
SavedMU=MU;
SavedBB=BB;
SavedGG=GG;
SavedVAR=VAR;
SavedCointegrationVectors=CointegrationVectors;
SavedLoadings=Loadings;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
OrthCompLoadings=orthcomp(Loadings);
OrthCompCointegrationVectors=orthcomp(CointegrationVectors);
C=OrthCompCointegrationVectors*MyInverse(OrthCompLoadings'*(eye(N)-sum(reshape(CSI,N,N,LagOrder),3))*OrthCompCointegrationVectors)*OrthCompLoadings'; %see lecture 5 slide 29
[VV,DD]=eig(VAR);
A0Start=VV*DD.^0.5;
LRIStart=C*A0Start;
[~,P]=Household(LRIStart(1,:)','U');
LRI=LRIStart*P';
A0=A0Start*P';
IRFs=GetIRFsCointegratedVECM(A0,BB,GG,LagOrder,N,Horizon,LRI,1);
[~,FEVs]=GetFractionsOfVarianceCointegratedVAR(A0,BB,GG,LagOrder,N,Horizon,1);
% subplot(1,2,1)
% plot(IRFs(:,[PositionConsumption PositionGDP ]))
% subplot(1,2,2)
% plot(FEVs(:,[PositionConsumption PositionGDP ]))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%[mumu,bbbb]=GetVARInLevelsFromVECM(MU,reshape(BB,N,N,LagOrder),GG);
%AbsRootsVARInLevels=sort(abs(varroots(LagOrder+1,N,[mumu reshape(bbbb,N,N*(LagOrder+1))])))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Here we turn everything into the notation of Geweke (1996):
Y_Geweke=Y';
Z_Geweke=Z(1:LagOrder*N+1,:)';
% A=[MU BB]';
% Theta=GG';
X_Geweke=LevelsLagged1Period';
%
X1_Geweke=X_Geweke(:,1:NumberOfCointegrationVectors);
X2_Geweke=X_Geweke(:,NumberOfCointegrationVectors+1:size(X_Geweke,2));
Theta21Hat=(X2_Geweke'*X2_Geweke)\(X2_Geweke'*X1_Geweke);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                       These are the priors:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Tau_Geweke=sqrt(10);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DiffY=diff(SavedY);
LaggedY=SavedY(1:T-1,:);
%
Y_KoopEtAl2010=DiffY;
X_KoopEtAl2010=LaggedY;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Hg=zeros(N,NumberOfCointegrationVectors);
Hg(PositionGDP,1)=1;
Hg(PositionConsumption,1)=-1;
%
H=Hg*sqrtm(MyInverse(Hg'*Hg));
OrthcompH=orthcomp(H);
%
TauBeta=1;  % Scalar between 0 and 1: 0 = informative, 1 = uninformative
PofTau=GetPofTau(H,OrthcompH,TauBeta);
Pof1OverTau=GetPofTau(H,OrthcompH,1/TauBeta);
InvPofTau=MyInverse(PofTau);
%
Alpha=Loadings;
Beta=CointegrationVectors;
%
% The parameters of the transformed specification:
Kappa=sqrtm(Alpha'*Alpha);
A=Alpha*MyInverse(Kappa);
B=Beta*Alpha'*A;
%
Nu=10; % Informative
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                       Here we do Gibbs-sampling
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alpha=Alpha;
beta=Beta;
sigma=VAR;
%
gg=alpha*beta';
theta=gg';
%
aa=[MU BB]';
%
% Draw=1;
% while Draw<=BurnInNN
%     if (Draw/1000)==fix(Draw/1000)
%         BurnInNN-Draw
%     end
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     U=((Y_Geweke'-(X_Geweke*theta+Z_Geweke*aa)'))'; % OK
%     sigma=iwishrnd(U'*U+Nu*alpha*MyInverse(beta'*Pof1OverTau*beta)*alpha',T+NumberOfCointegrationVectors); % OK
%     InvSigma=MyInverse(sigma); % OK
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     OmegaAlphaBar=MyInverse(kron(beta'*X_KoopEtAl2010'*X_KoopEtAl2010*beta,InvSigma)+(1/Nu)*kron(beta'*InvPofTau*beta,InvSigma));
%     AlphaBar=OmegaAlphaBar*vec(InvSigma*Y_KoopEtAl2010'*X_KoopEtAl2010*beta);
%     VecAlphaStar=MultivariateNormalRandomDraw(AlphaBar,OmegaAlphaBar,1)'; % OK
%     AlphaStar=reshape(VecAlphaStar,N,NumberOfCointegrationVectors); % OK
%     Astar=AlphaStar*MyInverse(sqrtm(AlphaStar'*AlphaStar));
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     AstarPrimeInvSigmaAstar=Astar'*InvSigma*Astar;
%     OmegaBBar=MyInverse(kron(AstarPrimeInvSigmaAstar,X_KoopEtAl2010'*X_KoopEtAl2010)+kron(AstarPrimeInvSigmaAstar,(1/Nu)*InvPofTau));
%     BBar=OmegaBBar*vec(X_KoopEtAl2010'*Y_KoopEtAl2010*InvSigma*Astar);
%     VecB=MultivariateNormalRandomDraw(BBar,OmegaBBar,1)'; % OK
%     B=reshape(VecB,N,NumberOfCointegrationVectors); % OK
%     beta=B*MyInverse(sqrtm(B'*B));
%     alpha=Astar*sqrtm(B'*B);
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     gg=alpha*beta';
%     theta=gg';
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     AHat=(Z_Geweke'*Z_Geweke)\(Z_Geweke'*(Y_Geweke-X_Geweke*theta)); % OK
%     InvSigmaKronZprimeZ=kron(InvSigma,Z_Geweke'*Z_Geweke); % OK
%     CovarianceMatrix=MyInverse(InvSigmaKronZprimeZ+eye((LagOrder*N+1)*N)*Tau_Geweke^2); % OK
%     Mean=(CovarianceMatrix*InvSigmaKronZprimeZ)*AHat(:); % OK
%     VecA=MultivariateNormalRandomDraw(Mean,CovarianceMatrix,1)';
%     aa=reshape(VecA,1+N*LagOrder,N); % OK
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     aaaa=aa';
%     mu=aaaa(:,1);
%     bb=aaaa(:,2:size(aaaa,2));
%     [mumu,bbbb]=GetVARInLevelsFromVECM(mu,reshape(bb,N,N,LagOrder),gg);
%     %
%     AbsRootsVARInLevels=sort(abs(varroots(LagOrder+1,N,[mumu reshape(bbbb,N,N*(LagOrder+1))])));
%     if sum(AbsRootsVARInLevels>(1+1.0e-10))==0 & sum(abs(AbsRootsVARInLevels-1)<1.0e-10)==(N-NumberOfCointegrationVectors)
%         Draw=Draw+1;
%     end
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NN=ErgodicNN/SamplingInterval;
% SIGMA=zeros(N,N,NN); % OK
% ALPHA=zeros(size(alpha,1),size(alpha,2),NN);
% BETAA=zeros(size(beta,1),size(beta,2),NN);
% BBBBB=zeros(N,1+N*LagOrder,NN);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% xx=0;
% Draw=1;
% while Draw<=ErgodicNN
%     if (Draw/1000)==fix(Draw/1000)
%         ErgodicNN-Draw
%     end
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     U=((Y_Geweke'-(X_Geweke*theta+Z_Geweke*aa)'))'; % OK
%     sigma=iwishrnd(U'*U+Nu*alpha*MyInverse(beta'*Pof1OverTau*beta)*alpha',T+NumberOfCointegrationVectors); % OK
%     InvSigma=MyInverse(sigma); % OK
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     OmegaAlphaBar=MyInverse(kron(beta'*X_KoopEtAl2010'*X_KoopEtAl2010*beta,InvSigma)+(1/Nu)*kron(beta'*InvPofTau*beta,InvSigma));
%     AlphaBar=OmegaAlphaBar*vec(InvSigma*Y_KoopEtAl2010'*X_KoopEtAl2010*beta);
%     VecAlphaStar=MultivariateNormalRandomDraw(AlphaBar,OmegaAlphaBar,1)'; % OK
%     AlphaStar=reshape(VecAlphaStar,N,NumberOfCointegrationVectors); % OK
%     Astar=AlphaStar*MyInverse(sqrtm(AlphaStar'*AlphaStar));
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     AstarPrimeInvSigmaAstar=Astar'*InvSigma*Astar;
%     OmegaBBar=MyInverse(kron(AstarPrimeInvSigmaAstar,X_KoopEtAl2010'*X_KoopEtAl2010)+kron(AstarPrimeInvSigmaAstar,(1/Nu)*InvPofTau));
%     BBar=OmegaBBar*vec(X_KoopEtAl2010'*Y_KoopEtAl2010*InvSigma*Astar);
%     VecB=MultivariateNormalRandomDraw(BBar,OmegaBBar,1)'; % OK
%     B=reshape(VecB,N,NumberOfCointegrationVectors); % OK
%     beta=B*MyInverse(sqrtm(B'*B));
%     alpha=Astar*sqrtm(B'*B);
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     gg=alpha*beta';
%     theta=gg';
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     AHat=(Z_Geweke'*Z_Geweke)\(Z_Geweke'*(Y_Geweke-X_Geweke*theta)); % OK
%     InvSigmaKronZprimeZ=kron(InvSigma,Z_Geweke'*Z_Geweke); % OK
%     CovarianceMatrix=MyInverse(InvSigmaKronZprimeZ+eye((LagOrder*N+1)*N)*Tau_Geweke^2); % OK
%     Mean=(CovarianceMatrix*InvSigmaKronZprimeZ)*AHat(:); % OK
%     VecA=MultivariateNormalRandomDraw(Mean,CovarianceMatrix,1)';
%     aa=reshape(VecA,1+N*LagOrder,N); % OK
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     aaaa=aa';
%     mu=aaaa(:,1);
%     bb=aaaa(:,2:size(aaaa,2));
%     [mumu,bbbb]=GetVARInLevelsFromVECM(mu,reshape(bb,N,N,LagOrder),gg);
%     AbsRootsVARInLevels=sort(abs(varroots(LagOrder+1,N,[mumu reshape(bbbb,N,N*(LagOrder+1))])));
%     if sum(AbsRootsVARInLevels>(1+1.0e-10))==0 & sum(abs(AbsRootsVARInLevels-1)<1.0e-10)==(N-NumberOfCointegrationVectors)
%         if (Draw/SamplingInterval)==fix(Draw/SamplingInterval)
%             xx=xx+1;
%             SIGMA(:,:,xx)=sigma; % OK
%             ALPHA(:,:,xx)=alpha; % OK
%             BETAA(:,:,xx)=beta; % OK
%             BBBBB(:,:,xx)=[mu bb]; % OK
%         end
%         Draw=Draw+1;
%     end
% end
% save(DFILE,varname(1,:),varname(2,:),varname(3,:),varname(4,:))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('PosteriorDistributionReducedFormVECM')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Here we get the reduced-form parameters of the VECM:
NN=size(BBBBB,3);
%
MUMU=zeros(N,1,NN);
BBBB=zeros(N,N*LagOrder,NN);
GGGG=zeros(N,N,NN);
%
VVVV=SIGMA;
Alpha=ALPHA;
Beta=BETAA;
%
Draw=1;
while Draw<=NN
    bb=BBBBB(:,:,Draw);
    mu=bb(:,1);
    bb=bb(:,2:size(bb,2));
    MUMU(:,:,Draw)=mu;
    BBBB(:,:,Draw)=bb;
    GGGG(:,:,Draw)=Alpha(:,:,Draw)*Beta(:,:,Draw)';
    Draw=Draw+1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                            Here we impose the zero and sign restrictions
clear A0 IN MU BB GG UU CC
%
NumberOfRandomRotationMatrices=5000;
Horizon=15*4;
NumberOfGibbsSamplingIterations=5;
SignRestrictionsHorizon=0;
%
A0=zeros(N,N,NumberOfRandomRotationMatrices,NN);
IN=zeros(NN,1);
MU=zeros(N,1,NN);
BB=zeros(N,N*LagOrder,NN);
GG=zeros(N,N,NN);
CC=zeros(N,N,NN);
%
zz=0;
Draw=1;
%while Draw<=NN
%    mu=MUMU(:,:,Draw);
%    bb=BBBB(:,:,Draw);
%    gg=GGGG(:,:,Draw);
%    var=VVVV(:,:,Draw);
%    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    [vv,dd]=eig(var);
%    A0Start=vv*dd.^0.5;
%    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    AlphaOrthComp=orthcomp(Alpha(:,:,Draw));
%    BetaOrthComp=orthcomp(Beta(:,:,Draw));
%    C=BetaOrthComp*MyInverse(AlphaOrthComp'*(eye(N)-sum(reshape(bb,N,N,LagOrder),3))*BetaOrthComp)*AlphaOrthComp';
%    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    [a0,IND]=getA0_paper(bb,gg,var,C,N,LagOrder,NumberOfRandomRotationMatrices,SignRestrictionsHorizon);
%    if sum(IND)>0>0 & DetectInf(a0)==0
%        zz=zz+1;
%        Index=sum(IND);
%        A0(:,:,1:Index,zz)=a0;
%        IN(zz,1)=Index;
%        MU(:,:,zz)=mu;
%        BB(:,:,zz)=bb;
%        GG(:,:,zz)=gg;
%        CC(:,:,zz)=C;
%        sum(IN)
%        disp('Expected number of successful draws:')
%        (sum(IN)/Draw)*NN
%    end
%    Draw=Draw+1;
%end
%
%save('VariablesAfterZerosAndSignRestrictions.mat','A0','IN', 'MU', 'BB', 'GG','CC')
%save('A0Variable.mat', 'A0', '-v7.3') % need to save it like that because the size of the variable is too large
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('VariablesAfterZerosAndSignRestrictions.mat')
load('A0Variable.mat')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                           Here we re-organize the successful draws
clear AA00 MUMU BBBB GGGG CCCC UUUU
%
NNN=sum(IN);
NN=length(IN);
%
AA00=zeros(N,N,NNN);
MUMU=zeros(N,1,NNN);
BBBB=zeros(N,N*LagOrder,NNN);
GGGG=zeros(N,N,NNN);
CCCC=zeros(N,N,NNN);
%
xx=0;
Draw=1;
while Draw<=NN
    Index=IN(Draw);
    zz=1;
    while zz<=Index
        xx=xx+1;
        MUMU(:,:,xx)=MU(:,:,Draw);
        BBBB(:,:,xx)=BB(:,:,Draw);
        GGGG(:,:,xx)=GG(:,:,Draw);
        CCCC(:,:,xx)=CC(:,:,Draw);
        AA00(:,:,xx)=A0(:,:,zz,Draw);
        zz=zz+1;
    end
    Draw=Draw+1;
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %                                                           The fractions of FEV:

NumberOfPermanentGDPShocks=1;
NNN=size(AA00,3);
% %
FEVs=zeros(Horizon+1,N,NumberOfPermanentGDPShocks,NNN);
% %
Draw=1;
while Draw<=NNN
    mu=MUMU(:,:,Draw);
    bb=BBBB(:,:,Draw);
    gg=GGGG(:,:,Draw);
    cc=CCCC(:,:,Draw);
    a0=AA00(:,:,Draw);
    [~,FEVs(:,:,:,Draw)]=GetFractionsOfVarianceCointegratedVAR(a0,bb,gg,LagOrder,N,Horizon,NumberOfPermanentGDPShocks);
    Draw=Draw+1;
end
%
Hor=(0:1:Horizon)';
%

Shock=1;
while Shock<=NumberOfPermanentGDPShocks
    Variable=1;
    while Variable<=N
        Fevs=ExtractPercentiles(sort(squeeze(FEVs(:,Variable,Shock,:)),2)',[0.5 0.16 0.84 0.05 0.95]')';
        figure(3)
        subplot(1,N,(Shock-1)*(N)+Variable)
        plot(Hor,Fevs(:,1),'k',Hor,Fevs(:,2:3),'r:',Hor,Fevs(:,4:5),'r','LineWidth',2)
        axis([0 15*4 0 1])
        xlabel('Horizon (quarters ahead)')
        if Variable==1
            if Shock==1
                ylabel('GDP shock')
            end
        end
        if Shock==1
            if Variable==PositionConsumption
                title('Consumption')
            elseif Variable==PositionGDP
                title('GDP')
            elseif Variable==PositionGDPDeflator
                title('GDP Deflator')
            elseif Variable==PositionUnemployment
                title('Unemployment Rate')
            elseif Variable==PositionShortRate
                title('Short rate')
            else
            end
        end
        Variable=Variable+1;
    end
    Shock=Shock+1;
end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Fevs=ExtractPercentiles(sort(squeeze(sum(FEVs(:,PositionGDP,1:2,:),3))'),[0.5 0.16 0.84 0.05 0.95]')';
%figure(3)
%subplot(1,N+1,N+1)
%plot(Hor,Fevs(:,1),'k',Hor,Fevs(:,2:3),'r:',Hor,Fevs(:,4:5),'r','LineWidth',2)
%axis([0 15*4 0 1])
%title('Fraction of the FEV of GDP jointly explained by BQ and H shock')
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %                                                             The IRFs:
NNN=size(AA00,3);
% %
Horizon=15*4;
%Horizon=0;
IRFs=zeros(Horizon+1,N,NumberOfPermanentGDPShocks,NNN);
Draw=1;
while Draw<=NNN
    mu=MUMU(:,:,Draw);
    bb=BBBB(:,:,Draw);
    gg=GGGG(:,:,Draw);
    cc=CCCC(:,:,Draw);
    a0=AA00(:,:,Draw);
    lri=cc*a0;
    irfs=GetIRFsCointegratedVECM(a0,bb,gg,LagOrder,N,Horizon,eye(N),1);
    IRFs(:,:,:,Draw)=irfs*sign(lri(1,1)); %normalization
    LRI(:,:,Draw)=lri*sign(lri(1,1)); %normalization
    Draw=Draw+1;
end
Index=find(squeeze(LRI(PositionGDP,1,:))>0);
IRFs=IRFs(:,:,:,Index);
IRFs(:,:,1,:)=IRFs(:,:,1,:)/median(squeeze(IRFs(Horizon+1,PositionGDP,1,:)));
%IRFs(:,:,2,:)=IRFs(:,:,2,:)/median(squeeze(IRFs(Horizon+1,PositionGDP,2,:)));
%
Hor=(0:1:Horizon)';
%
Shock=1;
while Shock<=NumberOfPermanentGDPShocks
    Variable=1;
    while Variable<=N
        Irfs=ExtractPercentiles(sort(squeeze(IRFs(:,Variable,Shock,:)),2)',[0.5 0.16 0.84 0.05 0.95]')';
        figure(2)
        subplot(1,N,(Shock-1)*N+Variable)
        plot(Hor,Irfs(:,1),'k',Hor,Irfs(:,2:3),'r:',Hor,Irfs(:,4:5),'r',Hor,zeros(size(Hor)),'b:','LineWidth',2)
        xlim([0 Horizon])
        xlabel('Horizon (quarters after shock)')
        if Variable==1
            if Shock==1
                ylabel('GDP shock')
            end
        end
        if Shock==1
            if Variable==PositionConsumption
                title('Consumption')
            elseif Variable==PositionGDP
                title('GDP')
            elseif Variable==PositionGDPDeflator
                title('GDP Deflator')
            elseif Variable==PositionUnemployment
                title('Unemployment Rate')
            elseif Variable==PositionShortRate
                title('Short rate')
            else
            end
        end
        Variable=Variable+1;
    end
    Shock=Shock+1;
end





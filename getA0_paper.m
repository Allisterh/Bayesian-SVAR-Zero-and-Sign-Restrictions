function [AA00,A0IN]=getA0_paper(bb,gg,var,cc,N,LagOrder,MaximumNumberOfRandomRotationMatrices,SignRestrictionsHorizon);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PositionConsumption=1;
PositionGDP=2;
PositionGDPDeflator=3;
PositionUnemployment=4;
PositionShortRate=5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
IndicesRestrictions=[PositionGDP PositionGDPDeflator PositionUnemployment PositionShortRate]';
SignsRestrictions=[1 -1 1 -1; 1 1 -1 1; -1 1 1 1; -1 -1 1 1]';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N=size(var,1);
[vv,dd]=eig(var);
A0Start=vv*dd.^0.5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
LRIStart=cc*A0Start;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% I: The number of zero restrictions:
NumberOfZeroImpacts=1;
Z=zeros(NumberOfZeroImpacts,N);
Z(1,1)=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
AA00=zeros(N,N,MaximumNumberOfRandomRotationMatrices);
A0IN=zeros(MaximumNumberOfRandomRotationMatrices,1);
%
zz=1;
%one permanent shock one GDP and sign restrictions on the others 

while zz<=MaximumNumberOfRandomRotationMatrices
    % (1) Here we randomly draw, one column a time, the random rotation matrix
    %     that is going to jointly impose the zero and sign restrictions
    %     
    Q=[];
    jj=1;
    while jj<=N
        if jj<=N-1
            RjofA0ApPlus=[Z*LRIStart; Q'];
        else
            RjofA0ApPlus=[Q'];
        end
        NN=null(RjofA0ApPlus);
        xj=randn(N,1);
        qj=NN*((NN'*xj)/norm(NN'*xj,2));
        Q=[Q qj];
        jj=jj+1;
    end
    % Check that Q is indeed a rotation matrix: Q*Q'
    % 
    A0=fliplr(A0Start*Q);
    LRI=cc*A0;
    A0=A0*sign(LRI(1,1));
    LRI=cc*A0;
    % Here we check that the sign restricitons are in fact satisfied:
    SignImpactAt0=sign(A0(IndicesRestrictions,2:N));
    if any(SignImpactAt0(:)-SignsRestrictions(:))==0
        IND=1;
    elseif any(-SignImpactAt0(:)-SignsRestrictions(:))==0
        A0=-A0;
        IND=1;
    else
        IND=0;
    end
    if IND==1
        AA00(:,:,zz)=A0;
        A0IN(zz)=1;
    end
    %
    zz=zz+1;
end
%
Index=find(A0IN);
%
AA00=AA00(:,:,Index);
A0IN=A0IN(Index);
% 



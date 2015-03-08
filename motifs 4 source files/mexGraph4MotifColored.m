% function Out=mexGraph4MotifColored(G,NodeColors2D)
% Input:
% G is a Graph object
% NodeColors2D is a Nx2 matrix of [NodeID ColorID]
% Output:
function Out=mexGraph4MotifColored(G,NodeColors2D)
m=mex4ColorMotifs(G,NodeColors2D);
n=size(m,1);
K4=zeros(n,4);
nColors=max(NodeColors2D(:,2));
for i=1:n
    K4(i,1:4)=Num2Bits(m(i,2),4,nColors);
end
%% Shrinking the result by mulitplicity of motives
ClrPower=(nColors.^[0:3])';
FullM=zeros(n,3);
load('SelfPerms4.mat');
for i=1:n
    MotifIdx=m(i,1);
    CurClrs=K4(i,:);
    CurPerms=SelfPerms{MotifIdx};
    nPerms=numel(CurPerms);
    CurMinPermScr=zeros(1,nPerms);
    for j=1:nPerms
        PermedClrs=CurClrs(CurPerms{j});
        CurMinPermScr(j)=PermedClrs*ClrPower;
    end
    MinClr=min(CurMinPermScr);
    FullM(i,1:3)=[MotifIdx MinClr m(i,3)];
end
[U , ~, IB]=unique(FullM(:,1:2),'rows');
nU=size(U,1);
SMat=zeros(nU,3);
for i=1:size(U,1)
    SMat(i,1:3)=[U(i,1) U(i,2) sum(FullM(IB==i,3))];
end
%% Output
SK4=zeros(size(SMat,1),4);
for i=1:size(SMat,1)
    SK4(i,1:4)=Num2Bits(SMat(i,2),4,nColors);
end
Out=[SMat(:,[1 3]) SK4];
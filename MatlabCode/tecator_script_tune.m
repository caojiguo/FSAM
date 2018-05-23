% This script analyze the tecator data. 


tecator = dlmread('U:\fda\fplr\FPLR_Paper2016Dec07\r code\tecator.txt','', 150, 0);
tec=mat2vec(tecator');

all_samples=NaN(240,125);
for i=1:240
    all_samples(i,:)=tec((i-1)*125+1:i*125);    
end


allPCs=all_samples(:,101:122);
Allcontent=all_samples(:,123:end);

train=all_samples(1:129,:);
monit=all_samples(130:172,:);
test=all_samples(173:215,:);
E1=all_samples(216:223,:);
E2=all_samples(224:240,:);

s = RandStream('mcg16807','Seed',18265);
RandStream.setDefaultStream(s);

A=randsample(size(all_samples,1),185);
list=zeros(size(all_samples,1),1);
list(A)=1;

train_spec=all_samples(logical(list),1:100);
test_spec=all_samples(~logical(list),1:100);
train_response=all_samples(logical(list),123:end);
test_response=all_samples(~logical(list),123:end);

addpath('U:\fda\fplr\FPLR_Paper2016Dec07\r code\matlab code\PACE_release2.11\release2.11\PACE');

trainY=train_response(:,3); % protein
testY=test_response(:,3);
nm=(851:2:1050)';


[ntrain,numgrid]=size(train_spec);
ntest=size(test_spec,1);
train_cell=mat2cell(train_spec,ones(1,ntrain),numgrid)';
t_cell=mat2cell(repmat(nm',ntrain,1),ones(1,ntrain),numgrid)';


param_X = setOptions('regular',2,'selection_k',20,'corrPlot',0,'rho',-1,'ngrid',55); 
trainRes= FPCA(train_cell, t_cell, param_X);   %perform PCA on x


%save('FPCres_tune.mat','trainY','testY','train_spec','test_spec','trainRes');


numBasis= getVal(trainRes,'no_opt');
trainPCscore = getVal(trainRes,'xi_est'); 
Phihat=getVal(trainRes,'phi');
hatsig2_x=getVal(trainRes,'sigma');
lamhat=getVal(trainRes,'lambda');
Corx_hat=getVal(trainRes,'xcorr');
Mu_x=getVal(trainRes,'mu');

nm=(851:2:1050)';
[ntest,numgrid]=size(test_spec);
ntrain=size(train_spec,1);
train_cell=mat2cell(train_spec,ones(1,ntrain),numgrid)';
t_cell=mat2cell(repmat(nm',ntrain,1),ones(1,ntrain),numgrid)';
test_cell=mat2cell(test_spec,ones(1,ntest),numgrid)';
[testhat, testPCscore, testPC_var]=FPCApred(trainRes,test_cell,t_cell, 2);






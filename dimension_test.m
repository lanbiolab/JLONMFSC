warning off
clear
clc
close all
clear memory;
addpath('fea')

%dataset = {'BreastTumor','Forebrain','GMvsHEK','GMvsHL'};
%dataset = {'Darmanis','Mouse2','Haber','Chu'};
%dataset = {'Pollen','Mouse2','Chu','Haber'};
dataset = {'Goolam'};
collect_result_ACC=[];
collect_result_NMI=[];
collect_result_AMI=[];
collect_result_ARI=[];
collect_time=[];
r0=8e-3;
r1=3e-2;  
r2=10;
r3=5;
r=10;
root=100;
for ii = 1:length(dataset)
    disp(dataset{ii});  
    for ss = 1:r
        %d=root*ss;
        d=100;
        path_data = ['E:/single_data/',dataset{ii},'_ready.txt'];
        path_label = ['E:/single_data/',dataset{ii},'_ready_label.txt'];
        %% hyper paramters setting 
        fea=load(path_data);  
        gnd=load(path_label);
        selected_class =length(unique(gnd));
        fea = double(fea);  %% 
        nnClass = length(unique(gnd));     % The number of classes
        select_sample = [];
        select_gnd    = [];
        for i = 1:selected_class  %% 
            idx = find(gnd == i);
            idx_sample    = fea(idx,:);
            select_sample = [select_sample;idx_sample];
            select_gnd    = [select_gnd;gnd(idx)];
        end
        fea = select_sample';
        fea = abs(fea);
        fea = fea./repmat(sqrt(sum(fea.^2)),[size(fea,1) 1]);%归一化
        gnd = select_gnd;
        c   = selected_class;
        X = fea;  %% genenum*cellnum
        clear fea select_gnd select_sample idx
        %
        tic()
        % ---------- initilization Z -------- %
        options = [];
        options.NeighborMode = 'KNN';
        options.k = 10;
        % options.WeightMode = 'Binary';      % Binary  HeatKernel
        options.WeightMode = 'Cosine';
        Z = constructW(X',options);
        Z_ini = Z;
   
        max_iter=250;
        Ctg = inv((X')*X+2*eye(size(X,2)));
        [ng, nc]=size(X);%% ng--number of genes; nc--number of cells
    
        XX=load(path_data);
        XX = abs(XX);
        XX = XX';%%%%%%%Rows are genes, columns are cell sample
        XX(all(XX == 0,2),:)=[];
        %==============Constructing a weight matrix==============
        %Preset value before constructing weight matrix
        options1 = [];
        options1.Metric = 'Cosine';
        options1.NeighborMode = 'KNN';%KNN
        options1.k =10;%5 nearest neighbors
        options1.WeightMode = 'Cosine';%Weights are 0 or 1, it can eplace with 'HeatKernel', 'Euclidean' 
        W = constructW(XX',options1);
        clear options1;
        options1 = [];
        %[m,n]=size(XX);
        n=nnClass;
        k1=n+d; 
        k2=n;
        [Z_final,A_final,B_final, H_ini] = DR_source(XX, k1, k2, W, options,r2,r3);
        %[W, H_ini] = nmf(XX, k2, 250);
        %disp(H_ini);
        % ---------- obtain similarity matrix ------- %
        [Z,S,U,E,obj_error] = LRR(X,Z_ini,H_ini,r0,r1,max_iter,Ctg);
        clear W_ini H_ini Ctg
        similarity=(abs(Z)+abs(Z'))/2;
        [result_label, kerNS]= SpectralClustering(similarity, nnClass);  %
        time_nmflrr=toc();
        % ---------- evaluation ------- %
        NMI=Cal_NMI(gnd, result_label);
        AMI=Cal_AMI(gnd, result_label);
        ARI=Cal_ARI(gnd,max(gnd),result_label,max(result_label));
        ACC=ACC_ClusteringMeasure(gnd, result_label);
    
        %fprintf(['NMI_for_ ' dataset{ii} ' is %f\n'],NMI)
        %fprintf(['AMI_for_ ' dataset{ii} ' is %f\n'],AMI)
        %fprintf(['ARI_for_ ' dataset{ii} ' is %f\n'],ARI)
        %fprintf(['ACC_for_ ' dataset{ii} ' is %f\n'],ACC)
        %fprintf(['Time_for_ ' dataset{ii} ' is %f\n'],time_nmflrr)
    
        fprintf('%f\t',NMI);
        fprintf('%f\t',AMI);
        fprintf('%f\t',ARI);
        fprintf('%f\t',ACC);
        fprintf('%f\t',time_nmflrr);
        fprintf('\n'); 
    
        collect_result_NMI(ii)=NMI;
        collect_result_AMI(ii)=AMI;
        collect_result_ARI(ii)=ARI;
        collect_result_ACC(ii)=ACC;
        collect_time(ii)=time_nmflrr;
        clear LZ DZ Z fea
    end
    
end


% save collect_result_ACC.mat collect_result_ACC
% save collect_result_NMI.mat collect_result_NMI

% xlswrite(['NMI_collect_result_9data5' '.xlsx'],collect_result_NMI)
% xlswrite(['AMI_collect_result_9data5' '.xlsx'],collect_result_AMI)
% xlswrite(['ARI_collect_result_9data5' '.xlsx'],collect_result_ARI)
% xlswrite(['ACC_collect_result_9data5' '.xlsx'],collect_result_ACC)
% xlswrite(['time_collect_result_9data5' '.xlsx'],collect_time)



function [localX,coverage] = localize( C )
%C is the coefficient matrix
%[tmp,ind] = sort(C,1,'descend');
[m,n]=size(C);
localX=C;
coverage=zeros(1,n);
for i=1:n
    thr=C(i,i)/2; %%为何取1.5?
    localX(localX(:,i)<thr,i)=0;  %% 如果C(i,j)小于对角线C（i，i）/1.5的值，那么设置为0
    coverage(1,i)=mean(C(i,i)./localX(localX(:,i)>thr,i));
end
end



function [groups, kerNS] = SpectralClustering(CKSym,n)
%--------------------------------------------------------------------------
% This function takes an adjacency matrix of a graph and computes the
% clustering of the nodes using the spectral clustering algorithm of
% Ng, Jordan and Weiss.
% CMat: NxN adjacency matrix
% n: number of groups for clustering
% groups: N-dimensional vector containing the memberships of the N points
% to the n groups obtained by spectral clustering
%--------------------------------------------------------------------------
% Copyright @ Ehsan Elhamifar, 2012
% Modified @ Chong You, 2015
%--------------------------------------------------------------------------
warning off;
N = size(CKSym,1);
% Normalized spectral clustering according to Ng & Jordan & Weiss
% using Normalized Symmetric Laplacian L = I - D^{-1/2} W D^{-1/2}
DN = diag(1./sqrt(sum(CKSym)+eps) );   %eps=2.2204e-16
LapN = speye(N) - DN * CKSym * DN;  % speye(N)生成N*N对角元素为1的矩阵  构建拉普拉斯矩阵
[~,~,vN] = svd(LapN); %% 奇异值分解
kerN = vN(:,N-n+1:N);
%kerN = vN(:,N-12:N);
normN = sum(kerN .^2, 2) .^.5;%% normalize the matrix U by L2-Norm
kerNS = bsxfun(@rdivide, kerN, normN + eps);  %% 矩阵kerN 每一行除以normN+eps
%-------------
%Y = pdist(kerNS);
%Z = linkage(Y);
%groups = cluster(Z,'maxclust',n);
MAXiter = 1000; % Maximum number of iterations for KMeans
REPlic = 100; % Number of replications for KMeans
groups = kmeans(kerNS,n,'maxiter',MAXiter,'replicates',REPlic,'EmptyAction','singleton');  %% 调用kmeans聚类

end



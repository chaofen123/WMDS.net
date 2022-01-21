% % % % % % % At the whole level
% % % % % % %  BLCA BRCA COAD HNSC KICH KIRC KIRP LUAD LUSC LIHC PRAD STAD THCA UCEC % % % % % % %     
clear
load('Gene_Net.mat');
node=unique(edge);
tumor=importdata('TCGA_BLCAtumor.txt');
express_gene=tumor.textdata(2:end,1);
tumor_expression=tumor.data;
In=ismember(express_gene,node);
tumor_expression(In==0,:)=[];
normal=importdata('TCGA_BLCAnormal.txt');
normal_expression=normal.data;
normal_expression(In==0,:)=[];
N=40
Normal_num=min(N,length(normal_expression(1,:)))   
Tumor_num=N
rand('seed',7);
t=randperm(length(tumor_expression(1,:)));
n=randperm(length(normal_expression(1,:)));
Tumor_data=tumor_expression(:,t(1:Tumor_num));
Normal_data=normal_expression(:,n(1:Normal_num)); 
[~,location1]=ismember(edge(:,1),node);
[~,location2]=ismember(edge(:,2),node);
l=location1.*location2;
Interact=[location1 location2];
Interact(l==0,:)=[];
N1=length(node);
N2=size(Interact,1);
Net=zeros(N1,N1);
for i=1:N2  
   Net(Interact(i,2),Interact(i,1))=1;
   Net(Interact(i,1),Interact(i,2))=1; 
   
end
Net=Net-diag(diag(Net));
p=Fisher(Tumor_data,Normal_data);
p(p==0)=eps;
P7=1./p;
C7=P7.*Net; 
C7=log(C7);
C7(isinf(C7))=0;
sum_p=sum(C7,2);
weight=1./sum_p;
weight(isinf(weight))=eps;
p(p>=0.05)=0;
p(p~=0)=1; 
C=p.*Net;   
C=C-diag(diag(C));
C=C+eye(N1);
[e,~]=find(C);
table=tabulate(e);
[~,degree_rank]=sort(table(:,2));
degree_top=degree_rank(length(degree_rank)-20:length(degree_rank));
hub=degree_top;
hub1=hub;
 for j=1:12
    hub=hub1;
    for i=1:length(hub)
        [h1,h2]=find(C(hub(i),:)==1);
        hub1=union(h2,hub1);
    end
    if length(hub1)==length(hub)
        break
    end
 end
 f = weight;
 intcon = 1:N1;
 A = -C;
 b = -ones(N1,1);
 lb = zeros(N1,1);
 ub = ones(N1,1);
 options = optimoptions('intlinprog','Display','off');
 x = intlinprog(f,intcon,A,b,[],[],lb,ub,options);
[x1,~]=find(x);       
x2=intersect(x1,hub);
gene=node(x2);




%%
% % % % % % % At the Personal level
% % % % % % %  BLCA BRCA COAD HNSC KICH KIRC KIRP LUAD LUSC LIHC PRAD STAD THCA UCEC % % % % % % %     
clear
load('Gene_Net.mat');
node=unique(edge);
tumor=importdata('TCGA_BLCAtumor.txt');
express_gene=tumor.textdata(2:end,1);
tumor_expression=tumor.data;
In=ismember(express_gene,node);
tumor_expression(In==0,:)=[];
normal=importdata('TCGA_BLCAnormal.txt');
normal_expression=normal.data;
normal_expression(In==0,:)=[];
rand('seed',7);
Normal_num=min(30,length(normal_expression(1,:)));   
n=randperm(length(normal_expression(1,:)));
Tumor_data=tumor_expression;
Normal_data=normal_expression(:,n(1:Normal_num)); 
[~,location1]=ismember(edge(:,1),node);
[~,location2]=ismember(edge(:,2),node);
l=location1.*location2;
Interact=[location1 location2];
Interact(l==0,:)=[];
N1=length(node);
N2=size(Interact,1);
Net=zeros(N1,N1);
for i=1:N2  
   Net(Interact(i,2),Interact(i,1))=1;
   Net(Interact(i,1),Interact(i,2))=1; 
   end
Net=Net-diag(diag(Net));
for m=1:30
    m
 Tumor_sample_data=Tumor_data(:,m);
p=SSN(Tumor_sample_data,Normal_data);
p(p==0)=eps;
P7=1./p;
C7=P7.*Net; 
C7=log(C7);
C7(isinf(C7))=0;
sum_p=sum(C7,2);
weight=1./sum_p;
weight(isinf(weight))=eps;
p(p>=0.05)=0;
p(p~=0)=1; 
C=p.*Net;   
C=C-diag(diag(C));
C=C+eye(N1);
[e,~]=find(C);
table=tabulate(e);
[~,degree_rank]=sort(table(:,2));
degree_top=degree_rank(length(degree_rank)-20:length(degree_rank));
hub=degree_top;
hub1=hub;
 for j=1:12
    hub=hub1;
    for i=1:length(hub)
        [h1,h2]=find(C(hub(i),:)==1);
        hub1=union(h2,hub1);
    end
    if length(hub1)==length(hub)
        break
    end
 end
 f = weight;
 intcon = 1:N1;
 A = -C;
 b = -ones(N1,1);
 lb = zeros(N1,1);
 ub = ones(N1,1);
 options = optimoptions('intlinprog','Display','off');
 x = intlinprog(f,intcon,A,b,[],[],lb,ub,options);
[x1,~]=find(x);       
x2=intersect(x1,hub);
gene{m}=node(x2);
 clear A C C7 p P7 R0
end
classdef Bombing < handle

% 算法具体步骤
% 1，	计算全局高斯核密度（混合k邻近图） k1=10，k2=2
% 2，	    选取Gaus密度最高的点，进行可达区域搜索计算，将所有可达性较高的点纳入该点所属的簇，
%       对于遍历到的可达性较低的点随着搜索层数增加设置二级边界点。
% 	        同时将遍历到其他簇的边界点数量分别记录下来，选择接壤数量最多的簇作为目标簇。
%       根据相对隶属度（簇的数目除以接壤数）判断，如果相对隶属度较低则将该簇作为目标簇的边界点。此时簇的数目不会增加。
% 3，	选取其余点中高斯密度最高且未被可达性搜索访问到的点，重复2的计算，直到寻找到k个簇。
%       若是无参k的算法，会进行到所有数据点均被访问
% 4，	对于所有离群点或未被可达性搜索访问到的点，寻找比自己Gaus密度大且距离最近的点（若k近邻图中不连通，则取欧氏距离）

    %
    % Properties:
    % ----------
    %   data            -   数据
    %
    %   nPoints         -   数据点的个数
    %
    %   nDims           -   数据点的维度
    %   
    %   K               -   聚类的数目，如果K大于0，说明需要分为K个簇，K<=0则是无参算法
    %                         
    %   G               -   混合k近邻图的邻接矩阵，G(i,j)表示数据点i到j的欧氏距离，若为inf则说明不存在边                     
    %
    %   knn_neigh       -   每个数据点的k近邻关系，分别从近到远存储数据点的编号
    %
    %   knn_dist        -   每个数据点的k近邻关系，分别从近到远存储数据点的欧氏距离
    %                       
    %   dist            -   dist(i,j)存储数据点i到j的在混合k近邻图的测地距离
    %                       
    %   visited         -   广度优先遍历中被访问的数据点visited(i)=1
    %
    %   must_notreach   -   广度优先遍历中边界数据点must_notreach(i)>=1
    %
    %   reach_index     -   广度优先遍历中，访问到的数据点的混合密度如果达到了reach_index，则作为可达区域
    %
    %   tocluster_index -   相对隶属度，若大于该值，说明可以独立成簇（是论文中T2的倒数）
    %
    %   dc              -   截断距离，仅仅影响Gaus密度，由于数据进行了归一化，dc=0.1(默认值)
    %
    %   in_cluster;     -   分类的标签
    %
    %
    %
    %
    % Methods:
    % -------
    %   Bombing(data)           -传入data数据，初始化算法
    %   
    %   set_dc_ri_ti(dc,ri,ti)  -设置截断距离，可达阈值，相对隶属度阈值,不设置则采用默认值
    %   
    %   find_cluster(K)         -调用算法,传入期望聚类数K，若K<=0,则是无参算法

    % Written by feizexuan
    % 15/8/2023
    
    properties
        data
        K=0
        nPoints
        nDims
        %与图有关的变量
        G
        knn_neigh
        knn_dist
        dist
        visited
        must_notreach;%must_notreach=1说明是边界点
        outlier_incluster;%outlier_incluster=1表示边界点属于某个簇
        reach_index=0.6
        tocluster_index=5
        outnode_level=3;
        is_tempcluster=0;
        %与密度有关的变量
        dc=0.1
        density_Gaus
        density_Gaus_line
        density_sorted
        W=0.5%两次修正Gaus密度与density_neigh加权
        %与分类有关
        in_cluster;
        noise_label;
        noise_rate=0;

        clusterid=0;%当前分配簇的编号
        clusterid0=0;
        clusterid_outlier={};%簇号对于它的边界点
        temp_clusterid_outlier={}
        temp_clusterid=0
        cluster_record={}
        startnode=0;
        %与去噪有关
        is_findnoise=0
        must_notreach_notnoise={}
        noisenum
        %控制相关
        is_show_watermelon=0;%展示算法每次迭代结果
    end
    
    methods
        
        function self = Bombing( data )
            %数据归一化
            y=size(data,2);
            min_dim=min(data);
            max_dim=max(data);
            data(:,1:y)=data(:,1:y)-min_dim(1:y);
            data=data./(max_dim-min_dim);

            self.data = data;
            self.nPoints = size( data,1 );
            self.nDims = size( data,2 );
            [self.G,self.knn_neigh,self.knn_dist]=self.makeknnG(data,10,3);%测试无参
            self.dist = self.G_dist(data,self.G);
            self.visited=zeros(self.nPoints,1);
            
            self.must_notreach=zeros(self.nPoints,1);
            self.outlier_incluster=zeros(self.nPoints,1);
            self.in_cluster=zeros(self.nPoints,1);
        end
        function self = set_dc_ri_ti(self, dc ,ri , ti)
            self.dc=dc;%默认0.1
            self.reach_index=ri;%默认0.6
            self.tocluster_index=ti;%默认5.5
        end

        function noise = find_noise(self,noise_rate)
            %去噪声所需要的参数
            self.noise_rate=noise_rate;
            self.noisenum=self.nPoints*self.noise_rate;
            self.is_findnoise=1;
            self.outnode_level=1;
            self.tocluster_index=50;
            %算法执行
            self.run_alg();
            noise=self.noise_label;
        end
        function idx = find_cluster(self,K)%K大于0则有类别，K<=0说明不设置K
            self.K=K;
            self.run_alg();
            idx=self.in_cluster;
        end
    end
    % public methods
    
    
    %% private methods
    methods(Access=private)
        %% 获取k近邻图G，与k近邻关系knn_neigh,knn_dist
        function [G,knn_neigh,knn_dist] = makeknnG(self,data,k1,k2)
            %makeknnG(data,k1,k2),data是数据，k1是互为k近邻图，k2是全连接k近邻图
            %G为生成的k近邻图的邻接矩阵
            %knn_neigh（n*k）存储n个数据点的k个近邻下标
            %knn_dist（n*k）存储n个数据点分别到k个近邻的欧氏距离

            KDT=KDTreeSearcher(data);%使用k近邻树计算k近邻
            n=self.nPoints;%n为数据点的个数
            G=zeros(n,n);

            %k近邻全连接
            %只要是任何一方的k近邻都会存在边
            % knn_neigh=zeros(n,k);
            % knn_dist=zeros(n,k);
            % for point_i=1:n
            %     [neigh,dist]=knnsearch(KDT,data(point_i,:),'k',k);%neigh存储k个近邻的下标，dist存储到这些近邻的距离
            %     knn_neigh(point_i,:)=neigh;
            %     knn_dist(point_i,:)=dist;
            %     G(point_i,neigh)=dist;
            %     G(neigh,point_i)=dist;
            % end

            %互为k近邻连接
            %G用于构成k近邻图，G_count计算每条边被记录几次，如果小于2次说明不是互为k近邻
            %knn_neigh存储i号节点的k个近邻下标，knn_dist记录对应距离,并且已经按照距离排好序
            G_count=G;
            knn_neigh=zeros(n,k1);
            knn_dist=zeros(n,k1);
            for point_i=1:n
                [neigh,dist]=knnsearch(KDT,data(point_i,:),'k',k1);%neigh存储k个近邻的下标，dist存储到这些近邻的距离
                knn_neigh(point_i,:)=neigh;
                knn_dist(point_i,:)=dist;

                G(point_i,neigh)=dist;
                G(neigh,point_i)=dist;
                G_count(point_i,neigh)=G_count(point_i,neigh)+1;
                G_count(neigh,point_i)=G_count(neigh,point_i)+1;
            end
            for x=1:n
                for y=1:n
                    if G_count(x,y)==1
                        a=find(knn_neigh(x,:)==y);
                        b=find(knn_neigh(y,:)==x);
                        if size(a,2)==0 && b<=k2
                            G(x,y)=knn_dist(y,b);
                        elseif size(b,2)==0 && a<=k2
                            G(x,y)=knn_dist(x,a);
                        else
                            G(x,y)=0;
                        end
                    elseif G_count(x,y)==0
                        G(x,y)=0;
                    end
                end
            end
            % if size(data,2)==2%如果是二维数据集则展示k邻近图
            %     figure('Name','knnDPC');
            %     for i=1:n
            %         for j=1:n
            %             if G(i,j)>0 
            %                 line([data(i,1),data(j,1)],[data(i,2),data(j,2)],'color','k');
            %                 %line([起点横坐标，终点横坐标],[起点纵坐标，终点纵坐标])
            %                 hold on;
            %             end
            %         end
            %     end
            %     plot(data(:,1),data(:,2),'linestyle',"none",'marker','x','color','k');
            %     hold on;
            % end

        end
        %% dist(i,j)存储的是i到j单源最短路径的长度
        function [dist] = G_dist(self,data,G)
            G = sparse(G) ;
            x=size(data,1);
            A = zeros(x,x);
            for i = 1 : x
                [dd , ~] = dijkstra_sp(G,i);%计算出真实图的单元最短路径
                A(i,:) = dd;
            end
            dist = A ;%d(i,j)存储的是i到j的单元最短路径的长度
        end
        %% 得到密度
        function [density_Gaus,density_Gaus_line,density_sorted,dc_percent] = get_density_Gaus(self,data,dc,dist)
            %计算k近邻图上的Gaus密度，Gaus密度归一化，以及该密度的排序
            ND=self.nPoints;
            dc_percent=dc/(max(data(:,1))-min(data(:,1)));
            % 密度计算
            for i=1:ND
                Density(i)=0.;
            end
            % k近邻图上的高斯核函数
            for i=1:ND-1
                for j=i+1:ND
                    %i结点依次加上其余所有点即j
                    %其余所有点j也均加上i点
                    Density(i)=Density(i)+exp(-(dist(i,j)/dc)*(dist(i,j)/dc));
                    Density(j)=Density(j)+exp(-(dist(i,j)/dc)*(dist(i,j)/dc));
                end
            end
            [~,density_sorted]=sort(Density,'descend');%density_sorted记录密度排序后的结点编号
            density_Gaus=Density;
            %将密度进行归一化
            Density=Density';
            y=size(Density,2);
            min_dim=min(Density);
            max_dim=max(Density);
            Density(:,1:y)=Density(:,1:y)-min_dim(1:y);
            density_Gaus_line=Density./(max_dim-min_dim);
        end
        %% 算法执行
        function self = run_alg( self )
            [self.density_Gaus,self.density_Gaus_line,self.density_sorted,dc_percent]=self.get_density_Gaus(self.data,self.dc,self.dist);
            while true
                if size(self.startnode,2)==0
                    break;
                end
                self.clusterid=self.clusterid+1;
                if self.clusterid0+1==self.clusterid
                    self.clusterid_outlier=[self.clusterid_outlier;{self.clusterid,[]}];
                    self.clusterid0=self.clusterid0+1;
                else
                    self.clusterid_outlier{self.clusterid,2}=[];
                end

                %   聚类第一阶段
                %   1，选取startnode核心点
                %   2，计算局部近邻密度
                %   3，计算第一次修正Gaus密度
                ND=self.nPoints;
                must_reach=zeros(self.nPoints,1);
                self.startnode=[];
                %选出Gaus密度最高的且没归入任何簇的点，并设置局部近邻半径
                for i=1:ND
                    if self.in_cluster(self.density_sorted(i))==0 && self.must_notreach(self.density_sorted(i))==0
                        range=1.1*self.knn_dist(self.density_sorted(i),size(self.knn_neigh,2));
                        self.startnode=self.density_sorted(i);
                        break
                    end
                end
                if size(self.startnode,1)==0%如果startnode未选中，也就是所有点都选完了
                    density_neigh=[];
                    Density_sigmod=[];
                    Density_final=[];
                else
                    %计算局部近邻密度
                    density_neigh=zeros(ND,1);
                    for i=1:ND
                        a=find(self.knn_dist(i,:)>=range);
                        if size(a,2)~=0
                            density_neigh(i)=a(1);
                        else
                            %如果找不到比range大的值则取10
                            density_neigh(i)=10;
                        end
                    end

                    %求高斯核的第一次修正的结果
                    max_Gaus=self.density_Gaus(self.startnode);
                    %     min_Gaus=min(density_Gaus(:));
                    min_Gaus=0;
                    slope=10;
                    movedist=(max_Gaus-min_Gaus)/2+min_Gaus;
                    Density_sigmod= 1 ./ (1 + exp(-(slope/max_Gaus)*(self.density_Gaus-movedist)));
                    Density_sigmod=Density_sigmod';


                    %根据Gaus-sigmod函数高于0.9设为must-reach，核心区域密度必可达
                    for i=1:size(self.in_cluster,1)
                        if Density_sigmod(i)>=0.9 && self.in_cluster(i)==0
                            must_reach(i)=1;
                        end
                    end
                    self.in_cluster(self.startnode)=self.clusterid;%密度最高点设置为簇编号

%                     %测试只有一次sigmod函数与局部近邻加权的结果
%                     max_dbscan=max(density_neigh(self.in_cluster==0));
%                     Density_neigh_linear= density_neigh / max_dbscan;
%                     Density_final=self.W1*Density_sigmod+(1-self.W1)*Density_neigh_linear;
                end

                %   聚类第二阶段
                %   1，计算二次修正Gaus密度
                %   2，计算混合密度
                %   3，广度优先遍历形成簇
                %   4，求未归入簇也没有must-reach的其余点的最大Gaus密度
                max_Gaus=[];
                min_Gaus=[];
                flag=1;
                for i=1:size(self.density_Gaus,2)
                    if self.in_cluster(i)>0 || must_reach(i)==1 || self.must_notreach(i)==1
                        continue;
                    end
                    if flag
                        max_Gaus=self.density_Gaus(i);
                        flag=0;
                    end
                    max_Gaus=max(max_Gaus,self.density_Gaus(i));
                end

                if size(max_Gaus,2)~=0%若得不到maxGaus说明剩余点都被标记访问了
                    slope=10;
                    %slope是斜率，density_Gaus是所有数据点高速密度值，max_Gaus是最大Gaus密度的值
                    %     min_Gaus=min(density_Gaus(:));
                    min_Gaus=0;
                    movedist=(max_Gaus-min_Gaus)/2+min_Gaus;
                    density_sigmod= 1 ./ (1 + exp(-(slope/max_Gaus)*(self.density_Gaus-movedist)));
                    density_sigmod=density_sigmod';

                    %寻找除去归入簇和must-reach点其余点的max-dbscan
                    flag=1;
                    max_dbscan=[];
                    for i=1:size(density_neigh,1)
                        if (self.in_cluster(i)>0) || must_reach(i)==1 || self.must_notreach(i)==1
                            continue;
                        end
                        if flag
                            max_dbscan=density_neigh(i);
                            flag=0;
                        end
                        max_dbscan=max(max_dbscan,density_neigh(i));
                    end

                    if size(max_dbscan,1)~=0
                        %只有当局部近邻密度和二次修正Gaus密度都存在的时候才计算混合密度
                        %否则混合密度为空
                        Density_neigh_linear= density_neigh / max_dbscan;
                        %只观察两次sigmod函数的结果
                        density_mix=self.W*density_sigmod+(1-self.W)*Density_neigh_linear;
                    else
                        density_sigmod=[];
                        density_mix=[];

                    end
                else
                    density_sigmod=[];
                    density_mix=[];
                end
                [self.visited,self.in_cluster,self.must_notreach,self.clusterid_outlier,self.outlier_incluster]=self.BFS_density_reach(self.G,density_mix,self.startnode,self.visited,self.in_cluster,must_reach,self.must_notreach,self.clusterid_outlier,self.clusterid,self.outlier_incluster);


                %     展示每次迭代完毕的效果图
                if self.is_show_watermelon
                    figure();
                    cmap = hsv(self.K+1); %// define your colormap here
                    reach_outlier=zeros(self.nPoints,1);
                    reach_outlier(self.in_cluster>0)=1;
                    for i=1:self.nPoints
                        if reach_outlier(i)==0 && self.must_notreach(i)~=0
                            reach_outlier(i)=self.must_notreach(i)+1;
                        end
                    end
                    figure;
                    cmap(1,:)=[0 0 0];%设置颜色
                    cmap(2,:)=[1 0 0];
                    cmap(3,:)=[0 0.69 0.31];
                    cmap(4,:)=[0.6 1 0.6];
                    cmap(5,:)=[0.5 0.5 0.5];
                    gscatter(self.data(:,1), self.data(:,2), reach_outlier, cmap);
                    hold on;
                    plot(self.data(self.startnode,1),self.data(self.startnode,2),'linestyle',"none",'marker','o','MarkerEdgeColor', 'k','MarkerFaceColor', 'y','Markersize',5);
                end
                %     展示Gaus_line,Gaus_sigmoid,Gaus_sigmoid2
                %     figure();
                %     subplot(2,3,1);
                %     h=makehotmap(data,density_Gaus);
                %     title('densityGaus');
                %     subplot(2,3,4);
                %     Linechart(density_Gaus_line,0.01,'density_Gaus_line');
                %
                %     subplot(2,3,2);
                %     h=makehotmap(data,density_sigmod);
                %     title('densitysigmod');
                %     subplot(2,3,5);
                %     Linechart(density_sigmod,0.01,'density_sigmod');
                %
                %     subplot(2,3,3);
                %     h=makehotmap(data,density_sigmod2);
                %     title('densitysigmod2')
                %     subplot(2,3,6);
                %     Linechart(density_sigmod2,0.01,'density_sigmod2');

                if self.clusterid==self.K || isempty(find(self.visited(:)==0, 1))
                    break
                end
            end
            if self.clusterid==1
                %如果只发现了一个簇，那么用K=2进行簇的归入
                self.K=2;
            end
            %若找到的簇的数量小于聚类数目K，则需要根据相对隶属度来选择最有可能独立的簇
            if self.clusterid<self.K 
                if size(self.cluster_record,2)~=0
                    self.cluster_record=sortrows(self.cluster_record,1,"descend");
                    j=1;
                    while(j<=size(self.cluster_record,1)) && self.clusterid<self.K
                        self.clusterid=self.clusterid+1;
                        self.in_cluster(self.cluster_record{j,2})=self.clusterid;
                        j=j+1;
                    end
                end
            end
            %噪声修改
            if self.is_findnoise==1
                if size(self.cluster_record,2)~=0
                    self.cluster_record=sortrows(self.cluster_record,1,"descend");
                    j=1;
                    visited_outlier_incluster=zeros(size(self.clusterid_outlier,1),1);
                    while(j<=size(self.cluster_record,1)) && (length(find(self.must_notreach(:)==0))/self.nPoints)<=0.8
                        self.must_notreach(self.cluster_record{j,2})=0;%簇中的数据标记
                        self.must_notreach(self.cluster_record{j,3})=0;%和其他簇接壤的数据纪录
%                         to_clusterid=self.cluster_record{j,4};
%                         outliernum=size(self.clusterid_outlier{to_clusterid,2},1);%接壤的簇的边界点总数
%                         adjnum=size(self.cluster_record{j,3},2);%本次的接壤数
%                         visited_outlier_incluster(to_clusterid)=visited_outlier_incluster(to_clusterid)+adjnum/outliernum;
%                         if visited_outlier_incluster(to_clusterid)>0.5;
%                             %                             self.must_notreach_notnoise=[self.must_notreach_notnoise;{[]}];
%                             self.must_notreach(self.clusterid_outlier{to_clusterid,2})=0;
%                         end
                        j=j+1;
                    end
                    %此时选出噪声点的数量应该是多于应该选出的
                    self.noise_label=zeros(self.nPoints,1);
                    self.noise_label(self.must_notreach>0)=1;
                    self.noise_label(self.visited==0)=1;
                    for i=1:self.nPoints
                        a=length(find(self.noise_label(:)==1));
                        if self.noise_label(self.density_sorted(i))==1 && length(find(self.noise_label(:)==1))>=self.noisenum
                            self.noise_label(self.density_sorted(i))=0
                        end
                    end
                else
                    self.noise_label=zeros(self.nPoints,1);
                end
            end

            %算法迭代结束，将边界点和未访问点归入密度较高最近的簇中
            if self.is_findnoise==0
                self.in_cluster=self.Outlier2cluster(self.data,self.dist,self.density_Gaus,self.density_sorted,self.knn_neigh,self.knn_dist,self.in_cluster);
            end
        end
        %% 广度优先遍历算法

        function [visited,in_cluster,must_notreach,clusterid_outlier,outlier_incluster] = BFS_density_reach(self,G,density_mix,startnode,visited,in_cluster,must_reach,must_notreach,clusterid_outlier,clusterid,outlier_incluster)
            if size(startnode,1)~=0%如果选取了核心点
                must_notreach0=must_notreach;
                visited0=visited;
                count_visited=1;%计算访问到的新点的数量,也就是该簇的总数目包括边界点，初始为1因为核心点算一个
                Maxn=size(G,1);
                queue=[];
                visited(startnode)=1;
                queue=[queue,startnode];
                visit_outlier_incluster=zeros(clusterid,1);
                visit_outlier_incluster_temp=zeros(size(self.cluster_record,1),1);
                visit_outlier_incluster_id=cell(clusterid,1);
                visit_outlier_incluster_id_temp=cell(size(self.cluster_record,1),1);
                while(size(queue,2)~=0)
                    node=queue(1);
                    queue(1)=[];
                    node_neigh=find(G(node,:)~=0);
                    for i=1:size(node_neigh,2)
                        if(visited(node_neigh(i))==0 && size(density_mix,1)~=0 && must_notreach(node)<self.outnode_level)%如果没被访问过并且密度达标

                            if (density_mix(node_neigh(i))>=self.reach_index || must_reach(node_neigh(i))==1) && must_notreach(node)==0
                                in_cluster(node_neigh(i))=in_cluster(startnode);
                                visited(node_neigh(i))=1;
                                count_visited=count_visited+1;
                                queue=[queue,node_neigh(i)];
                            else
                                %开始设置一级边界点
                                must_notreach(node_neigh(i))=max(must_notreach(node_neigh(i)),must_notreach(node)+1);
                                visited(node_neigh(i))=1;
                                count_visited=count_visited+1;
                                queue=[queue,node_neigh(i)];
                            end
                        elseif visited(node_neigh(i))==0 && must_notreach(node)<self.outnode_level
                            %开始设置多级边界点
                            if must_reach(node_neigh(i))==1
                                in_cluster(node_neigh(i))=in_cluster(startnode);
                                visited(node_neigh(i))=1;
                                count_visited=count_visited+1;
                                queue=[queue,node_neigh(i)];
                            end
                        else
                            %接壤情况计算,需要知道must_notreach属于哪个簇
                            if outlier_incluster(node_neigh(i))==1%如果访问到了有主人的边界点
                                %查找是哪个主人
                                for j=1:clusterid-1
                                    if(ismember(node_neigh(i),clusterid_outlier{j,2}))
                                        visit_outlier_incluster(j)=visit_outlier_incluster(j)+1;
                                        %new
                                        visit_outlier_incluster_id{j}=[visit_outlier_incluster_id{j},node_neigh(i)];
                                    end
                                end
                                for j=1:size(self.cluster_record,1)
                                    if(ismember(node_neigh(i),self.cluster_record{j,6}))
                                        %访问到广义簇
                                        visit_outlier_incluster_temp(j)=visit_outlier_incluster_temp(j)+1;
                                        %访问到非广义簇
                                        visit_outlier_incluster_id_temp{j}=[visit_outlier_incluster_id_temp{j},node_neigh(i)];
                                    end
                                end
                                if max(visit_outlier_incluster)<max(visit_outlier_incluster_temp)
                                    self.is_tempcluster=1;
                                else
                                    self.is_tempcluster=0;
                                end
                            end

                        end
                    end
                end
                %广度优先遍历结束，开始判断当前簇是否需要归入其他簇
                %1，筛选接壤数量最多的簇，作为可能归入的对象
                %2，如果不需要归入，则将边界点作为该簇的边界点，否则都归入另一个簇作为其边界点
                index_to_cluster=inf;
                if self.is_tempcluster==0
                    [adj_num,to_cluster_id]=max(visit_outlier_incluster);
                else
                    [adj_num,to_cluster_id]=max(visit_outlier_incluster_temp);
                end

                %稍作修改
                if adj_num==0
                    adj_num=1;
                end
                index_to_cluster=count_visited/adj_num;


                if self.is_findnoise==1
                    index_to_cluster=count_visited^1.5/adj_num;
                elseif self.K<=0 && self.tocluster_index==5
                    self.tocluster_index=4.5;
                end

                %step2
                %获取当前簇的边界点
                outlier_new=[];
                for i=1:size(must_notreach,1)
                    if must_notreach0(i)==0 && must_notreach(i)>0
                        outlier_new=[outlier_new;i];
                        outlier_incluster(i)=1;
                    end
                end
                if index_to_cluster>self.tocluster_index%如果大于相对隶属度，则独立成簇
                    %作为正式簇
                    clusterid_outlier{clusterid,2}=outlier_new;
                else%若小于则作为其他簇的边界点,并且所有数据点被记录
                    %self.cluster_record用于记录当前簇的各种信息
                    %1，相对隶属度；2，与目标簇接壤的边界点；3，目标簇的簇号；4，目标簇是否是临时簇，5，当前簇的边界点
                    self.cluster_record=[self.cluster_record;{index_to_cluster,[],[],[],[],[]}];
                    reach_temp=[];
%                     must_notreach_temp=[];
                    for i=1:size(visited,1)
                        if visited0(i)==0 && visited(i)~=0%如果是本次访问的数据点，即当前簇
                            if(in_cluster(i)>0)
                                reach_temp=[reach_temp,i];
                            end
%                             if(must_notreach(i)>0)
%                                 must_notreach_temp=[must_notreach_temp,i];
%                             end
                            must_notreach(i)=self.outnode_level;%属于边界点
                            outlier_incluster(i)=1;%属于某簇的边界点
                            clusterid_outlier{to_cluster_id,2}=[clusterid_outlier{to_cluster_id,2};i];%记录属于哪个簇
                            in_cluster(i)=0;%不属于任何簇
                        end
                    end
                    %第一列记录相对隶属度，临时簇的簇编号隐含在行号
                    self.cluster_record{size(self.cluster_record,1),2}=reach_temp;%临时簇的可达区域
                    self.cluster_record{size(self.cluster_record,1),3}=[visit_outlier_incluster_id{to_cluster_id}];%临时簇与目标簇接壤的边界点
                    self.cluster_record{size(self.cluster_record,1),4}=to_cluster_id;%临时簇需要归入的簇
                    self.cluster_record{size(self.cluster_record,1),5}=self.is_tempcluster;%它的目标簇是不是临时簇
                    self.cluster_record{size(self.cluster_record,1),6}=outlier_new;%临时簇的边界点
                    self.clusterid=self.clusterid-1;
                end

            end
        end
        %% 边界点归入距离最近且Gaus密度比自身高的点中
        function [in_cluster] = Outlier2cluster(self,data,dist,density_Gaus,density_sorted,knn_neigh,knn_dist,in_cluster)
            outlier=find(in_cluster==0);
            outlier_belong=[outlier,zeros(size(outlier,1),1)];
            for i=1:size(outlier,1)
                point_i=outlier(i);

                flag=0;
                for j=1:size(knn_neigh,2)
                    if density_Gaus(knn_neigh(point_i,j))>density_Gaus(outlier(i))
                        flag=1;
                        outlier_belong(i,2)=knn_neigh(point_i,j);
                        break;
                    end
                end
                if flag==0%说明周围没有密度比自身高的结点
                    if dist(point_i,density_sorted(1))~=inf
                        min_dist=dist(point_i,density_sorted(1))
                    else
                        min_dist=pdist2(data(point_i,:),data(density_sorted(1),:))
                    end
                    outlier_belong(i,2)=density_sorted(1);
                    loc=find(density_sorted==point_i);
                    for j=2:loc-1
                        if dist(point_i,density_sorted(j))~=inf
                            dist_ij=dist(point_i,density_sorted(j));
                        else
                            dist_ij=pdist2(data(point_i,:),data(density_sorted(j),:));
                        end
                        if min_dist>dist_ij
                            min_dist=dist_ij;
                            outlier_belong(i,2)=density_sorted(j);
                        end
                    end
                end
            end
            for i=1:size(outlier_belong,1)
                point=outlier_belong(i,1);
                near_point=point;
                if in_cluster(point)==0
                    while 1
                        near_point=outlier_belong(find(outlier_belong(:,1)==near_point),2);
                        if in_cluster(near_point)~=0
                            in_cluster(point)=in_cluster(near_point);
                            break;
                        end
                    end
                end
            end
        end
    end
end
    

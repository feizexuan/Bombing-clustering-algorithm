classdef Bombing < handle

% 算法具体步骤
% 1.1	Calculate the k nearest neighboring points and Construct the joint k-nearest neighbor graph
%           we use function makeknnG to return an adjacency matrix of jNNG
%
% 1.2   Calculate the Gaussian kernel density of each point
%           we use function get_density_Gaus to get Density_Gaus
% 
% 2.1   Determine the core point as startnode which Denisty_Gaus is the max
% of all the point remained
%             
% 2.2   Calculate the modified Gaussian kernel density and the local neighbor density
%       and add up both to obtain mixture density
%
% 2.3   utilize breadth-first search to determine the reachable region,boundary region and number of adjacent points
%           we use function BFS_density_reach to get these areas
%
%       repeat 2.1~2.3 until getting K clusters or every point is visited
% 
% Finally:
% Assign the data points within the boundary region and unvisited data points to the clusters of the nearest data point with larger density value

    % Main Properties:
    % ----------
    %   data            -   Data
    %
    %   nPoints         -   Number of data points
    %
    %   nDims           -   Dimensions of the data point
    %   
    %   K               -   If K is greater than 0, it means that it needs to be divided into K clusters. If K<=0, it is a parameterless algorithm
    %                         
    %   G               -   Adjacency matrix of joint K-nearest neighbor graphs                   
    %
    %   knn_neigh       -   The K-nearest neighbor relationship of each data point stores the number of the data point from near to far
    %
    %   knn_dist        -   The K-nearest neighbor relationship of each data point stores the Euclidean distance of the data point from near to far
    %                       
    %   dist            -   dist(i,j) stores the geodesic distance of data points i to j in a mixed joint k-nearest neighbor map
    %                       
    %   visited         -   =1 means data points visited in BFS
    %
    %   must_notreach   -   >0 means the first or second level boundary point
    %
    %   reach_index     -   In BFS，mixture density >= reach_index means reachable region
    %
    %   tocluster_index -   Threshold of DA(degree of adjacent).The default value is 5 which is the reciprocal of T2(0.2) in the paper
    %
    %   dc              -   distance-cut value.dc=0.1(default)
    %
    %   in_cluster;     -   label
    %
    %
    %
    %
    % Main Methods:
    % -------
    %   Bombing(data)           -initialize the algorithm
    %   
    %   find_cluster(K)         -Call the algorithm, pass in the expected number of clusters K, if K<=0, then it is a parameterless algorithm

    % Written by feizexuan
    % 15/8/2023
    
    properties
        data
        K=0
        nPoints
        nDims
        %about graph
        G
        knn_neigh
        knn_dist
        dist
        visited
        must_notreach;%must_notreach=1说明是边界点
        outlier_incluster;%=1 means the boundary point belong to a cluster 
        reach_index=0.6
        tocluster_index=5
        outnode_level=3;
        is_tempcluster=0;
        %与密度有关的变量
        dc=0.1
        density_Gaus
        density_Gaus_line
        density_sorted
        W=0.5
        %与分类有关
        in_cluster;

        clusterid=0;%cluster id
        clusterid0=0;
        clusterid_outlier={};
        temp_clusterid_outlier={}
        temp_clusterid=0
        cluster_record={}
        startnode=0;
        %Presentation process
        is_show_watermelon=0;%=1 to show the results of each iteration of the algorithm
        time_BFS=0
        time_init_Gaus=0
        time_neigh=0
    end
    
    methods
        
        function self = Bombing( data )
            %数据归一化
            data=(data-min(data))./(max(data)-min(data));

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
            self.dc=dc;%0.1 as default
            self.reach_index=ri;%0.6 as default
            self.tocluster_index=ti;%5 as default
        end
        function idx = find_cluster(self,K)
            self.K=K;
            self.run_alg();
            idx=self.in_cluster;
        end
    end
    % public methods
    
    
    %% private methods
    methods(Access=private)
        %% Get the adjacency matrix of jNNG
        function [G,knn_neigh,knn_dist] = makeknnG(self,data,k1,k2)
            n=self.nPoints;
            G=zeros(n,n);
            %G is the adjacency matrix of jNNG
            [knn_neigh,knn_dist] = knnsearch(data,data,'k',k1);
            for point_i=1:n
                neigh=knn_neigh(point_i,:);
                dist=knn_dist(point_i,:);
                index=[];
                for j=k2+1:size(neigh,2)
                    neigh2=knn_neigh(neigh(j),:);
                    if isempty(find(neigh2==point_i,1))
                        index=[index,j];
                    end
                end
                dist(index)=0;
                G(point_i,neigh)=dist;
                G(neigh,point_i)=dist;
            end
              %show graph if it is two-dim dataset
%             if size(data,2)==2
%                 figure('Name','knnDPC');
%                 for i=1:n
%                     for j=1:n
%                         if G(i,j)>0
%                             line([data(i,1),data(j,1)],[data(i,2),data(j,2)],'color','k');
%                             %line([起点横坐标，终点横坐标],[起点纵坐标，终点纵坐标])
%                             hold on;
%                         end
%                     end
%                 end
%                 plot(data(:,1),data(:,2),'linestyle',"none",'marker','x','color','k');
%                 hold on;
%             end

        end
        %% dist(i,j) is the shortest path between Vi and Vj in the graph
        function [dist] = G_dist(self,data,G)
            G = sparse(G) ;
            [dist] = graphallshortestpaths(G);
        end
        %% get density
        function [density_Gaus,density_Gaus_line,density_sorted,dc_percent] = get_density_Gaus(self,data,dc,dist)
            ND=self.nPoints;
            dc_percent=dc/(max(data(:,1))-min(data(:,1)));
            t1=clock;
%           Calculate the Gaussian kernel density of each point
            dist2=exp(-(dist./dc).^2);
            Density=sum(dist2',1)-1;
            [~,density_sorted]=sort(Density,'descend');%density_sorted记录密度排序后的结点编号
            density_Gaus=Density;
            t2=clock;
            self.time_init_Gaus=self.time_init_Gaus+etime(t2,t1);
            density_Gaus_line=[];
        end
        %% Select core point and Modify density and BFS
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

                %   First stage
                %   1，Select startnode(The maximum Gaussian density of the remaining points)
                %   2，Calculate the modified Gaussian kernel density and the local neighbor density
                %   3，Calculate the first modified Gaus density
                ND=self.nPoints;
                must_reach=zeros(self.nPoints,1);
                self.startnode=[];
                %Select startnode and get local neighbor radius
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
                    
                    %Calculate the local neighbor density
                    t1=clock;
                    density_neigh=sum(self.knn_dist<=range,2);
                    density_neigh(density_neigh<10)=density_neigh(density_neigh<10)+1;
                    t2=clock;
                    self.time_neigh=self.time_neigh+etime(t2,t1);

                    %First modification to density
                    max_Gaus=self.density_Gaus(self.startnode);
                    %     min_Gaus=min(density_Gaus(:));
                    min_Gaus=0;
                    slope=10;
                    movedist=(max_Gaus-min_Gaus)/2+min_Gaus;
                    Density_sigmod= 1 ./ (1 + exp(-(slope/max_Gaus)*(self.density_Gaus-movedist)));
                    Density_sigmod=Density_sigmod';
                    must_reach=(Density_sigmod>=0.9)&(self.in_cluster==0);
                    self.in_cluster(self.startnode)=self.clusterid;%密度最高点设置为簇编号

                end

                %   Second stage
                %   1，Calculate the secondary modified Gaus density
                %   2，Calculate the mixture density
                %   3，BFS to get generalized cluster
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

                if size(max_Gaus,2)~=0
%                   Second modification to density
                    slope=10;
                    %     min_Gaus=min(density_Gaus(:));
                    min_Gaus=0;
                    movedist=(max_Gaus-min_Gaus)/2+min_Gaus;
                    density_sigmod= 1 ./ (1 + exp(-(slope/max_Gaus)*(self.density_Gaus-movedist)));
                    density_sigmod=density_sigmod';

                    %refresh max of  local neighbor radius
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
                        Density_neigh_linear= density_neigh / max_dbscan;
                        %Mixture Density
                        density_mix=self.W*density_sigmod+(1-self.W)*Density_neigh_linear;
                    else
                        density_sigmod=[];
                        density_mix=[];

                    end
                else
                    density_sigmod=[];
                    density_mix=[];
                end
                t1=clock;
                %BFS
                [self.visited,self.in_cluster,self.must_notreach,self.clusterid_outlier,self.outlier_incluster]=self.BFS_density_reach(self.G,density_mix,self.startnode,self.visited,self.in_cluster,must_reach,self.must_notreach,self.clusterid_outlier,self.clusterid,self.outlier_incluster);
                t2=clock;
                self.time_BFS=self.time_BFS+etime(t2,t1);
                %set is_show_watermelon=1 to show the effect of each algorithm iteration
                if self.is_show_watermelon
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
%               Stop when getting K clusters or each point is visited
                if self.clusterid==self.K || isempty(find(self.visited(:)==0, 1))
                    break
                end
            end
            if self.clusterid==1
                self.K=2;
            end
            %if num of cluster is smaller than K ,we choose some undetermined clusters
            %that are more likely to be real clusters according to the number of adjacent points
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

            %Assign the data points remained to the points with larger density value and shortest distance
            self.in_cluster=self.Outlier2cluster(self.data,self.dist,self.density_Gaus,self.density_sorted,self.knn_neigh,self.knn_dist,self.in_cluster);

        end
        %% BFS
        function [visited,in_cluster,must_notreach,clusterid_outlier,outlier_incluster] = BFS_density_reach(self,G,density_mix,startnode,visited,in_cluster,must_reach,must_notreach,clusterid_outlier,clusterid,outlier_incluster)
            if size(startnode,1)~=0
                must_notreach0=must_notreach;
                visited0=visited;
                count_visited=1;
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
                            %boundary point
                            if must_reach(node_neigh(i))==1
                                in_cluster(node_neigh(i))=in_cluster(startnode);
                                visited(node_neigh(i))=1;
                                count_visited=count_visited+1;
                                queue=[queue,node_neigh(i)];
                            end
                        else
                            if outlier_incluster(node_neigh(i))==1%If a boundary point with an owner is accessed
                                %find whose boundary point
                                for j=1:clusterid-1
                                    if(ismember(node_neigh(i),clusterid_outlier{j,2}))
                                        visit_outlier_incluster(j)=visit_outlier_incluster(j)+1;
                                        visit_outlier_incluster_id{j}=[visit_outlier_incluster_id{j},node_neigh(i)];
                                    end
                                end
                                for j=1:size(self.cluster_record,1)
                                    if(ismember(node_neigh(i),self.cluster_record{j,6}))
                                        visit_outlier_incluster_temp(j)=visit_outlier_incluster_temp(j)+1;
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
                %the end of BFS,then begin to determine whether the current cluster needs to be subdivided into other clusters
                index_to_cluster=inf;
                if self.is_tempcluster==0
                    [adj_num,to_cluster_id]=max(visit_outlier_incluster);
                else
                    [adj_num,to_cluster_id]=max(visit_outlier_incluster_temp);
                end
                if adj_num==0
                    adj_num=1;
                end
                index_to_cluster=count_visited/adj_num;
                if self.K<=0 && self.tocluster_index==5
                    self.tocluster_index=4.5;
                end
                %get boundary point
                outlier_new=[];
                for i=1:size(must_notreach,1)
                    if must_notreach0(i)==0 && must_notreach(i)>0
                        outlier_new=[outlier_new;i];
                        outlier_incluster(i)=1;
                    end
                end
                if index_to_cluster>self.tocluster_index%If it is greater than tocluster_index, it is clustered independently
                    clusterid_outlier{clusterid,2}=outlier_new;
                else
                    %self.cluster_record is used to record various information about the current cluster
                   
                    self.cluster_record=[self.cluster_record;{index_to_cluster,[],[],[],[],[]}];
                    reach_temp=[];
%                     must_notreach_temp=[];
                    for i=1:size(visited,1)
                        if visited0(i)==0 && visited(i)~=0
                            if(in_cluster(i)>0)
                                reach_temp=[reach_temp,i];
                            end
%                             if(must_notreach(i)>0)
%                                 must_notreach_temp=[must_notreach_temp,i];
%                             end
                            must_notreach(i)=self.outnode_level;
                            outlier_incluster(i)=1;
                            clusterid_outlier{to_cluster_id,2}=[clusterid_outlier{to_cluster_id,2};i];%记录属于哪个簇
                            in_cluster(i)=0;
                        end
                    end
                    self.cluster_record{size(self.cluster_record,1),2}=reach_temp;
                    self.cluster_record{size(self.cluster_record,1),3}=[visit_outlier_incluster_id{to_cluster_id}];
                    self.cluster_record{size(self.cluster_record,1),4}=to_cluster_id;
                    self.cluster_record{size(self.cluster_record,1),5}=self.is_tempcluster;
                    self.cluster_record{size(self.cluster_record,1),6}=outlier_new;
                    self.clusterid=self.clusterid-1;
                end

            end
        end
        %% The boundary points are grouped into the nearest points with a higher Gaus density than themselves
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
                if flag==0
                    if dist(point_i,density_sorted(1))~=inf
                        min_dist=dist(point_i,density_sorted(1));
                    else
                        min_dist=pdist2(data(point_i,:),data(density_sorted(1),:));
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
    

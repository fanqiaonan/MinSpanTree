%-------------------------用LDA进行数据分类完整程序------------------------%
clear
clc
P = 15; %被试个数
Fs = 200; %采样频率
channelnum = 62; %通道个数

tic
for k = 1:P
%     %------------------数据预处理：滤波并去公共平均------------------%
%     str_k = num2str(k);
%     filename = strcat('E:\Class\ADdata\signal',str_k,'.mat');
%     load (filename); %导入原始数据signal
%     f_signal = eegfilt(signal',256,0,47);
%     f_signal = eegfilt(f_signal,256,2,0); %滤波到2~47Hz
%     signalF = f_signal';
%     avg = mean(signalF(:,1:15),2); %15通道列平均
%     [a b] = size(signalF);
%     signalF(:,1:15) = signalF(:,1:15) - repmat(avg,1,b-1); %去公共平均
    
    str_k = num2str(k);
    filename = strcat('F:\脑纹识别\身份识别后续研究\脑网络\subject',str_k,'.mat');
    load (filename); %导入预处理后的数据signalF
    
    [B,A] = butter(2,[2/(Fs/2) 47/(Fs/2)]);
    signalF = filter(B,A,motionSignal(:,1:channelnum)); %滤波到2~47Hz
    avg = mean(signalF(:,1:channelnum),2); %22通道列平均
    [a b] = size(signalF);
    signalF(:,1:channelnum) = signalF(:,1:channelnum) - repmat(avg,1,b); %去公共平均
    
    %-------------------------滤波到所需频段------------------------%
%     [B,A] = butter(2,[4/(Fs/2) 8/(Fs/2)]); %滤波到theta频段
%     signal_filter = filter(B,A,signalF(:,1:channelnum));
%     [B,A] = butter(2,[8/(Fs/2) 13/(Fs/2)]); %滤波到alpha频段
%     signal_filter = filter(B,A,signalF(:,1:channelnum));
    [B,A] = butter(2,[13/(Fs/2) 30/(Fs/2)]); %滤波到beta频段
    signal_filter = filter(B,A,signalF(:,1:channelnum));
%     [B,A] = butter(2,[30/(Fs/2) 40/(Fs/2)]); %滤波到gamma频段
%     signal_filter = filter(B,A,signalF(:,1:channelnum));
    
    %-----------------------计算两两通道的PLV值----------------------%
    S = 16; %样本个数
    len = 30; %30s平均
    PLVM = zeros(channelnum,channelnum,S);
    PLVM_all = zeros(channelnum,channelnum,S);
    
    for j = 1:S
        for i = len*(j-1)+1:len*j
            allCH = signal_filter((i-1)*Fs+1:i*Fs,1:channelnum); %一秒钟数据计算一个PLV值
            [samples channels] = size(allCH);
            for ch1 = 1:channels
                for ch2 = ch1:channels
                    tmpCh1 = allCH(:,ch1);
                    tmpCh2 = allCH(:,ch2);
                    PLVch(ch1,ch2) = plv_hilbert(tmpCh1,tmpCh2); %plv_hilbert函数计算两通道间的相同步指数
                end
            end
            PLVch_all= PLVch+ PLVch';
            diacIdx = 1:(channelnum+1):channelnum*channelnum;
            PLVch_all(diacIdx) = 1;
            PLVM_all(:,:,j) = PLVM_all(:,:,j) + PLVch_all(:,:); %30s的PLV值叠加
            PLVch=[];PLVch_all=[];
        end
    end
    
    for j = 1:S
        PLVMavg_all(:,:,j) = PLVM_all(:,:,j)/len; %30s求平均
%          %画plv图
%           figure
%           subplot;
%           imagesc(PLVMavg_all);
%           caxis([0,1]);
%           colorbar;
%           set(gca, 'XTick', [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15]); %X坐标轴刻度数据点位置
%           set(gca,'XTickLabel',{'FPz','AF3','AF4','F3','Fz','F4','T7','C3','Cz','C4','T8','P3','Pz','P4','Oz'}); %X坐标轴刻度处显示的字符
%           set(gca, 'YTick', [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15]); %X坐标轴刻度数据点位置
%           set(gca,'YTickLabel',{'FPz','AF3','AF4','F3','Fz','F4','T7','C3','Cz','C4','T8','P3','Pz','P4','Oz'});
    end
    
    eval(['plv_',num2str(k),'=PLVMavg_all;']);
    eval(['plv=plv_',num2str(k),';']);
    trialNum=size(plv,3); %样本个数
    
    netPro.LF = []; % RnetPro.LF2 = []; 
    netPro.BC3 = [];
    netPro.TreeHierarchy = []; % netPro.TreeHierarchy2 = [];
    netPro.Ecc = [];
    netPro.Diam = []; 
    netPro.D = [];
    netPro.distance=[];
    
    for curTrial = 1:trialNum
        tmpAdj1 = plv(:,:,curTrial);%一个样本的PLV矩阵
        tmpAdj = Primf(tmpAdj1);%由prim算法得到最小生成树
        
        % 计算最小生成树属性-全局指标(Global measures)
        % 幸存率（Survival rate）
        
        % 分歧率（Divergence rate）
       
        % 连接熵（连接熵（Connection Entropy）
        
        % 叶分数（Leaf fraction）
        Leaf = []; Edge = [];
        netPro.Leafnum = []; netPro.Edgesnum = [];
        netPro.Leafnum(curTrial,:) = leaf_nodes(tmpAdj);
        netPro.Edgesnum(curTrial,:) = numedges(tmpAdj);
        Leaf = netPro.Leafnum(curTrial,:);
        [m,Leaf_num] = size(Leaf);
        Edge = netPro.Edgesnum(curTrial,:);
        netPro.LF(curTrial,:) = Leaf_num/Edge;
        % 偏心性（Eccentricity） 
        netPro.Ecc(curTrial,:) = vertex_eccentricity(tmpAdj);
        % 树的直径（Diameter）
        netPro.Diam(curTrial,:) = diameter(tmpAdj); 
        % 度（Degree）
        netPro.D(curTrial,:) = degrees(tmpAdj);
        % 平均路径长度
        netPro.distance(curTrial,:) = ave_path_length(tmpAdj);
        
        plv_net(curTrial,:)=[netPro.LF(curTrial,:),...
        netPro.Ecc(curTrial,:),netPro.Diam(curTrial,:),...
        netPro.D(curTrial,:),netPro.distance(curTrial,:)];%放入 plv_net，得到16*15，样本个数*样本特征
    end    
    eval(['netPro_',num2str(k),'=netPro; clear netPro;']);
    eval(['plv_net_',num2str(k),'=plv_net; clear plv_net;']);
end
toc

%-------------------------用LDA进行数据分类------------------------%
Tr = S/2; %训练样本
Te = S - Tr; %测试样本
threshold = 0.95;

%真实标记
label_real = [];
for i = 1:P
    label_real = [label_real i*ones(1,Te)];
end

for t = 1:10  %共进行10次分类
    nk = 0;
    nd = 0;
    for i = 1:P
        eval(['All=plv_net_',num2str(i),'; ']);   
        All=All'; % 15*16，16个样本，每个样本有15个特征
        temp = randperm(S);
        %获取训练数据
        for j = 1:Tr
            b = All(:,temp(j));%获取PLVall中的第i个被试的第temp（j）列
            nk = nk + 1;
            trainData(:,nk) = b;
        end
        %获取测试数据
        for j = Tr+1:S
            b = All(:,temp(j));
            nd = nd + 1;
            testData(:,nd) = b;
        end
    end
    
    [Wopt,BASE,PCAcorr4TR] = LDA_TRAIN1(trainData',threshold,P); %分类训练
    [LDAcorr4TR,LDAcorr4TE,LABELS_TEST] = LDA_TEST1(testData',Wopt,BASE,PCAcorr4TR,P); %分类测试
    
    compare = (label_real == LABELS_TEST);
    LDA_accuracy(t) = sum(compare)/length(compare); %分类正确率
end
%-------------------------��LDA�������ݷ�����������------------------------%
clear
clc
P = 15; %���Ը���
Fs = 200; %����Ƶ��
channelnum = 62; %ͨ������

tic
for k = 1:P
%     %------------------����Ԥ�����˲���ȥ����ƽ��------------------%
%     str_k = num2str(k);
%     filename = strcat('E:\Class\ADdata\signal',str_k,'.mat');
%     load (filename); %����ԭʼ����signal
%     f_signal = eegfilt(signal',256,0,47);
%     f_signal = eegfilt(f_signal,256,2,0); %�˲���2~47Hz
%     signalF = f_signal';
%     avg = mean(signalF(:,1:15),2); %15ͨ����ƽ��
%     [a b] = size(signalF);
%     signalF(:,1:15) = signalF(:,1:15) - repmat(avg,1,b-1); %ȥ����ƽ��
    
    str_k = num2str(k);
    filename = strcat('F:\����ʶ��\���ʶ������о�\������\subject',str_k,'.mat');
    load (filename); %����Ԥ����������signalF
    
    [B,A] = butter(2,[2/(Fs/2) 47/(Fs/2)]);
    signalF = filter(B,A,motionSignal(:,1:channelnum)); %�˲���2~47Hz
    avg = mean(signalF(:,1:channelnum),2); %22ͨ����ƽ��
    [a b] = size(signalF);
    signalF(:,1:channelnum) = signalF(:,1:channelnum) - repmat(avg,1,b); %ȥ����ƽ��
    
    %-------------------------�˲�������Ƶ��------------------------%
%     [B,A] = butter(2,[4/(Fs/2) 8/(Fs/2)]); %�˲���thetaƵ��
%     signal_filter = filter(B,A,signalF(:,1:channelnum));
%     [B,A] = butter(2,[8/(Fs/2) 13/(Fs/2)]); %�˲���alphaƵ��
%     signal_filter = filter(B,A,signalF(:,1:channelnum));
    [B,A] = butter(2,[13/(Fs/2) 30/(Fs/2)]); %�˲���betaƵ��
    signal_filter = filter(B,A,signalF(:,1:channelnum));
%     [B,A] = butter(2,[30/(Fs/2) 40/(Fs/2)]); %�˲���gammaƵ��
%     signal_filter = filter(B,A,signalF(:,1:channelnum));
    
    %-----------------------��������ͨ����PLVֵ----------------------%
    S = 16; %��������
    len = 30; %30sƽ��
    PLVM = zeros(channelnum,channelnum,S);
    PLVM_all = zeros(channelnum,channelnum,S);
    
    for j = 1:S
        for i = len*(j-1)+1:len*j
            allCH = signal_filter((i-1)*Fs+1:i*Fs,1:channelnum); %һ�������ݼ���һ��PLVֵ
            [samples channels] = size(allCH);
            for ch1 = 1:channels
                for ch2 = ch1:channels
                    tmpCh1 = allCH(:,ch1);
                    tmpCh2 = allCH(:,ch2);
                    PLVch(ch1,ch2) = plv_hilbert(tmpCh1,tmpCh2); %plv_hilbert����������ͨ�������ͬ��ָ��
                end
            end
            PLVch_all= PLVch+ PLVch';
            diacIdx = 1:(channelnum+1):channelnum*channelnum;
            PLVch_all(diacIdx) = 1;
            PLVM_all(:,:,j) = PLVM_all(:,:,j) + PLVch_all(:,:); %30s��PLVֵ����
            PLVch=[];PLVch_all=[];
        end
    end
    
    for j = 1:S
        PLVMavg_all(:,:,j) = PLVM_all(:,:,j)/len; %30s��ƽ��
%          %��plvͼ
%           figure
%           subplot;
%           imagesc(PLVMavg_all);
%           caxis([0,1]);
%           colorbar;
%           set(gca, 'XTick', [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15]); %X������̶����ݵ�λ��
%           set(gca,'XTickLabel',{'FPz','AF3','AF4','F3','Fz','F4','T7','C3','Cz','C4','T8','P3','Pz','P4','Oz'}); %X������̶ȴ���ʾ���ַ�
%           set(gca, 'YTick', [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15]); %X������̶����ݵ�λ��
%           set(gca,'YTickLabel',{'FPz','AF3','AF4','F3','Fz','F4','T7','C3','Cz','C4','T8','P3','Pz','P4','Oz'});
    end
    
    eval(['plv_',num2str(k),'=PLVMavg_all;']);
    eval(['plv=plv_',num2str(k),';']);
    trialNum=size(plv,3); %��������
    
    netPro.LF = []; % RnetPro.LF2 = []; 
    netPro.BC3 = [];
    netPro.TreeHierarchy = []; % netPro.TreeHierarchy2 = [];
    netPro.Ecc = [];
    netPro.Diam = []; 
    netPro.D = [];
    netPro.distance=[];
    
    for curTrial = 1:trialNum
        tmpAdj1 = plv(:,:,curTrial);%һ��������PLV����
        tmpAdj = Primf(tmpAdj1);%��prim�㷨�õ���С������
        
        % ������С����������-ȫ��ָ��(Global measures)
        % �Ҵ��ʣ�Survival rate��
        
        % �����ʣ�Divergence rate��
       
        % �����أ������أ�Connection Entropy��
        
        % Ҷ������Leaf fraction��
        Leaf = []; Edge = [];
        netPro.Leafnum = []; netPro.Edgesnum = [];
        netPro.Leafnum(curTrial,:) = leaf_nodes(tmpAdj);
        netPro.Edgesnum(curTrial,:) = numedges(tmpAdj);
        Leaf = netPro.Leafnum(curTrial,:);
        [m,Leaf_num] = size(Leaf);
        Edge = netPro.Edgesnum(curTrial,:);
        netPro.LF(curTrial,:) = Leaf_num/Edge;
        % ƫ���ԣ�Eccentricity�� 
        netPro.Ecc(curTrial,:) = vertex_eccentricity(tmpAdj);
        % ����ֱ����Diameter��
        netPro.Diam(curTrial,:) = diameter(tmpAdj); 
        % �ȣ�Degree��
        netPro.D(curTrial,:) = degrees(tmpAdj);
        % ƽ��·������
        netPro.distance(curTrial,:) = ave_path_length(tmpAdj);
        
        plv_net(curTrial,:)=[netPro.LF(curTrial,:),...
        netPro.Ecc(curTrial,:),netPro.Diam(curTrial,:),...
        netPro.D(curTrial,:),netPro.distance(curTrial,:)];%���� plv_net���õ�16*15����������*��������
    end    
    eval(['netPro_',num2str(k),'=netPro; clear netPro;']);
    eval(['plv_net_',num2str(k),'=plv_net; clear plv_net;']);
end
toc

%-------------------------��LDA�������ݷ���------------------------%
Tr = S/2; %ѵ������
Te = S - Tr; %��������
threshold = 0.95;

%��ʵ���
label_real = [];
for i = 1:P
    label_real = [label_real i*ones(1,Te)];
end

for t = 1:10  %������10�η���
    nk = 0;
    nd = 0;
    for i = 1:P
        eval(['All=plv_net_',num2str(i),'; ']);   
        All=All'; % 15*16��16��������ÿ��������15������
        temp = randperm(S);
        %��ȡѵ������
        for j = 1:Tr
            b = All(:,temp(j));%��ȡPLVall�еĵ�i�����Եĵ�temp��j����
            nk = nk + 1;
            trainData(:,nk) = b;
        end
        %��ȡ��������
        for j = Tr+1:S
            b = All(:,temp(j));
            nd = nd + 1;
            testData(:,nd) = b;
        end
    end
    
    [Wopt,BASE,PCAcorr4TR] = LDA_TRAIN1(trainData',threshold,P); %����ѵ��
    [LDAcorr4TR,LDAcorr4TE,LABELS_TEST] = LDA_TEST1(testData',Wopt,BASE,PCAcorr4TR,P); %�������
    
    compare = (label_real == LABELS_TEST);
    LDA_accuracy(t) = sum(compare)/length(compare); %������ȷ��
end
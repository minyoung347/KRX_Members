
%% 
% mean trading
% sort domestic & foreign list by corr with ifd, year
% for each decile
% output: [mean_vol_member, ca, cb, cc, domestic_sort, foreign_sort]
% 'ca' is the correlation domestic member(41) with ifd(3)
% 'cb' is the correlation foreign member(21) with ifd(3)
% 'cc' is the correlation institutions(8) with ifd(3)
% 'domestic_sort' is the index of domestic members according to ifd
% 'foreign_sort' is the index of foreign members according to ifd

% 기존에 상위 121개 기업으로 설정
% 지금은 코스피 또는 코스닥
ind_target = ind_kosdaq;
num_firm = length(ind_target);


price_ori = zeros(2722,num_firm);
for i = 1:2722
    firm_i = 1;
    for ii = 1:num_firm
        price_ori(i,firm_i) = price{ind_target(ii),1}(i,1);
        firm_i = firm_i + 1;
    end
end

% mean volume of member firm
% this matrix will be used in the inventory correlation with ifd
for i = 1:size(firm_net,1)
    sum_member = sum(firm_net{i,2}(1:2722,ind_target).*price_ori,2);
    ind = find(sum_member~=0); % 0인것 제외하고 평균
    mean_vol_member(i,1) = mean(sum_member(ind));
end
clear ind sum_member price_ori

figure
set(gcf,'color','w')
bar(sqrt(mean_vol_member([domestic,foreign])))
set(gca,'xtick',1:62,'xticklabel',inv_mem_name([domestic,foreign]),'fontsize',9)
xlabel('member')
ylabel('mean volume')
xtickangle(90)


% 'ca', 'cb', 'cc' will be used in the inventory correlation with ifd

% 'ca' is the correlation domestic member(41) with ifd(3)
a1 = zeros(num_firm,11);
a2 = zeros(num_firm,11);
a3 = zeros(num_firm,11);
for i = 1:41
    firm_i = 1;
    for f = 1:num_firm
        for y = 1:11
            if length(find(inv_mem_net{ind_target(f),2}(year_day{y,1},domestic(i))>0)) > length(year_day{y,1})/2
                a1(firm_i,y) = corr(inv_mem_net{ind_target(f),1}(year_day{y,1},1),inv_mem_net{ind_target(f),1}(year_day{y,1},domestic(i)));
                a2(firm_i,y) = corr(inv_mem_net{ind_target(f),1}(year_day{y,1},2),inv_mem_net{ind_target(f),1}(year_day{y,1},domestic(i)));
                a3(firm_i,y) = corr(inv_mem_net{ind_target(f),1}(year_day{y,1},3),inv_mem_net{ind_target(f),1}(year_day{y,1},domestic(i)));
            end
        end
        firm_i = firm_i + 1;
    end
    ca(i,1) = mean(delniz(a1(:),[1],0));
    ca(i,2) = mean(delniz(a2(:),[1],0));
    ca(i,3) = mean(delniz(a3(:),[1],0));
end
clear a1 a2 a3


% 'cb' is the correlation foreign member(21) with ifd(3)
a1 = zeros(num_firm,11);
a2 = zeros(num_firm,11);
a3 = zeros(num_firm,11);
for i = 1:21
    firm_i = 1;
    for f = 1:num_firm
        for y = 1:11
            if length(find(inv_mem_net{ind_target(f),2}(year_day{y,1},foreign(i))>0)) > length(year_day{y,1})/2
                a1(firm_i,y) = corr(inv_mem_net{ind_target(f),1}(year_day{y,1},1),inv_mem_net{ind_target(f),1}(year_day{y,1},foreign(i)));
                a2(firm_i,y) = corr(inv_mem_net{ind_target(f),1}(year_day{y,1},2),inv_mem_net{ind_target(f),1}(year_day{y,1},foreign(i)));
                a3(firm_i,y) = corr(inv_mem_net{ind_target(f),1}(year_day{y,1},3),inv_mem_net{ind_target(f),1}(year_day{y,1},foreign(i)));
            end
        end
        firm_i = firm_i + 1;
    end
    cb(i,1) = mean(delniz(a1(:),[1],0));
    cb(i,2) = mean(delniz(a2(:),[1],0));
    cb(i,3) = mean(delniz(a3(:),[1],0));
end
clear a1 a2 a3


% 'cc' is the correlation institutions(8) with ifd(3)
a1 = zeros(num_firm,11);
a2 = zeros(num_firm,11);
a3 = zeros(num_firm,11);
ins = 4:11;
for i = 1:8
    firm_i = 1;
    for f = 1:num_firm
        for y = 1:11
            if length(find(inv_mem_net{ind_target(f),2}(year_day{y,1},domestic(i))>0)) > length(year_day{y,1})/2
                a1(firm_i,y) = corr(inv_mem_net{ind_target(f),1}(year_day{y,1},1),inv_mem_net{ind_target(f),1}(year_day{y,1},ins(i)));
                a2(firm_i,y) = corr(inv_mem_net{ind_target(f),1}(year_day{y,1},2),inv_mem_net{ind_target(f),1}(year_day{y,1},ins(i)));
                a3(firm_i,y) = corr(inv_mem_net{ind_target(f),1}(year_day{y,1},3),inv_mem_net{ind_target(f),1}(year_day{y,1},ins(i)));
            end
        end
        firm_i = firm_i + 1;
    end
    cc(i,1) = mean(delniz(a1(:),[1],0));
    cc(i,2) = mean(delniz(a2(:),[1],0));
    cc(i,3) = mean(delniz(a3(:),[1],0));
end
clear a1 a2 a3


% sort by correlation with inv_mem_net of fore
[~,f1] = sort(cb(:,1),'descend');
[~,d1] = sort(ca(:,1),'descend');
[~,f2] = sort(cb(:,2),'descend');
[~,d2] = sort(ca(:,2),'descend');
[~,f3] = sort(cb(:,3),'descend');
[~,d3] = sort(ca(:,3),'descend');
domestic_sort(:,1) = domestic(d1);
domestic_sort(:,2) = domestic(d2);
domestic_sort(:,3) = domestic(d3);
foreign_sort(:,1) = foreign(f1);
foreign_sort(:,2) = foreign(f2);
foreign_sort(:,3) = foreign(f3);
clear f1 d1 f2 d2 f3 d3

%% PLOT Inventory correlation with ifd using 'ca', 'cb', 'cc'
% 'ca' is the correlation domestic member(41) with ifd(3)
% 'cb' is the correlation foreign member(21) with ifd(3)
% 'cc' is the correlation institutions(8) with ifd(3)

addpath('/home/minyoung/data_minyoung/Research/Matlab_function/linspecer')
color_3 = linspecer(3);

figure1 = figure;
set(gcf,'color','w')
subplot(1,2,1)
scatter(ca(:,1),ca(:,3),sqrt(mean_vol_member(domestic)/10^9*5),'filled','markerfacecolor',color_3(2,:),'MarkerEdgeColor',[0 0 0])
hold on
plot(ca(:,1),ca(:,3),'.','color',[0 0 0])
hold on
plot([-.15 .6],[0 0],'k--')
hold on
plot([0 0],[-.25 .2],'k--')
for i = 1:41
    if ismember(domestic(i),domestic([11,15,34,32,6,13,30,33,29,8,18,9,36,35,1,7,27,28,20,24,17,23,14,41,17,39,40]))
        
    elseif domestic(i) == domestic(31) % KTB
        text(ca(i,1)-0.015,ca(i,3)+0.015,inv_mem_name{domestic(i)},'fontsize', 12)
    elseif domestic(i) == domestic(19) % Samsung
        text(ca(i,1)-0.09,ca(i,3)+0.002,inv_mem_name{domestic(i)},'fontsize', 12)
    elseif domestic(i) == domestic(2) % Shinhan invest
        text(ca(i,1)-0.13,ca(i,3),inv_mem_name{domestic(i)},'fontsize', 12)
    elseif domestic(i) == domestic(38) % Hanwha invest
        text(ca(i,1)-0.05,ca(i,3)+0.01,inv_mem_name{domestic(i)},'fontsize', 12)
    elseif domestic(i) == domestic(5) % Mirae Daewoo
        text(ca(i,1)-0.135,ca(i,3),inv_mem_name{domestic(i)},'fontsize', 12)    
    elseif domestic(i) == domestic(10) % NH invest
        text(ca(i,1)-0.04,ca(i,3)+0.006,inv_mem_name{domestic(i)},'fontsize', 12)    
    elseif domestic(i) == domestic(12) % KB
        text(ca(i,1)-0.01,ca(i,3)-0.006,inv_mem_name{domestic(i)},'fontsize', 12)   
    elseif domestic(i) == domestic(3) % Korea invest
        text(ca(i,1)-0.125,ca(i,3),inv_mem_name{domestic(i)},'fontsize', 12)    
    elseif domestic(i) == domestic(21) % HI Invest
        text(ca(i,1)+.022,ca(i,3),inv_mem_name{domestic(i)},'fontsize', 12)
    elseif domestic(i) == domestic(26) % Korea Asset
        text(ca(i,1)+.015,ca(i,3),inv_mem_name{domestic(i)},'fontsize', 12)
    elseif domestic(i) == domestic(22) % Kiwoom
        text(ca(i,1)-0.03,ca(i,3)+0.025,inv_mem_name{domestic(i)},'fontsize', 12)   
    else
        text(ca(i,1)+.005,ca(i,3),inv_mem_name{domestic(i)},'fontsize', 12)
    end
end
text(-0.13,1.08,'(a)','Units', 'Normalized', 'VerticalAlignment', 'Top','fontsize',14)

xlabel('\rho_{ID,IV}','fontsize',16)
ylabel('\rho_{IS,IV}','fontsize',16)
%title('Domestic members','fontsize',16)
xlim([-.15 .6])
ylim([-.32 .21])

% % right upper
% a_ru = [11,15,34,32,6,13,30,33,29,8,18,9,36];
% axes('Position',[.65*.5 .65 .3*.5 .25])
% scatter(ca(a_ru,1),ca(a_ru,3),sqrt(mean_vol_member(domestic(a_ru))/10^6*5),'filled','markerfacecolor',color_3(2,:),'MarkerEdgeColor',[0 0 0])
% hold on
% plot(ca(a_ru,1),ca(a_ru,3),'k.')
% for i = 1:41
%     if sum(a_ru == i) > 0
%         if domestic(i) == domestic(9) % Meritz
%             text(ca(i,1)+.005,ca(i,3)-0.001,inv_mem_name{domestic(i)},'fontsize', 12)
%         elseif domestic(i) == domestic(18) % Golden Bridge
%             text(ca(i,1)+.005,ca(i,3)+0.001,inv_mem_name{domestic(i)},'fontsize', 12)
%         elseif domestic(i) == domestic(15) % Yuhwa
%             text(ca(i,1)-.022,ca(i,3)+0.001,inv_mem_name{domestic(i)},'fontsize', 12)
%         elseif domestic(i) == domestic(34) % IM
%             text(ca(i,1)+.005,ca(i,3)-0.001,inv_mem_name{domestic(i)},'fontsize', 12)
%         elseif domestic(i) == domestic(33) % BNK
%             text(ca(i,1)-.017,ca(i,3),inv_mem_name{domestic(i)},'fontsize', 12)
%         elseif domestic(i) == domestic(13) % Hanwha
%             text(ca(i,1)-.013,ca(i,3)+0.007,inv_mem_name{domestic(i)},'fontsize', 12)
%         elseif domestic(i) == domestic(36) % KB Invest
%             text(ca(i,1)+.008,ca(i,3)+0.001,inv_mem_name{domestic(i)},'fontsize', 12)
%         elseif domestic(i) == domestic(8) % Hanyang
%             text(ca(i,1),ca(i,3)-0.005,inv_mem_name{domestic(i)},'fontsize', 12)
%         else
%             text(ca(i,1)+.005,ca(i,3),inv_mem_name{domestic(i)},'fontsize', 12)
%         end
%     end
% end
% box on
% ylim([.105, .17])
% xlim([-.09, .024])
% 
% % left lower
% a_ll = [35,1,7,27,28,20,24,17,23,14,41,17,39,40];
% axes('Position',[.33*.5 .155 .35*.5 .27])
% scatter(ca(a_ll,1),ca(a_ll,3),sqrt(mean_vol_member(domestic(a_ll))/10^6*5),'filled','markerfacecolor',color_3(2,:),'MarkerEdgeColor',[0 0 0])
% hold on
% plot(ca(a_ll,1),ca(a_ll,3),'k.')
% for i = 1:41
%     if sum(a_ll == i) > 0
%         if domestic(i) == domestic(20) % DB Financial
%             text(ca(i,1)-.007,ca(i,3)+0.004,inv_mem_name{domestic(i)},'fontsize', 12)
%         elseif domestic(i) == domestic(20) % DB Financial
%             text(ca(i,1)-.007,ca(i,3)+0.004,inv_mem_name{domestic(i)},'fontsize', 12)
%         elseif domestic(i) == domestic(24) % Hana Financial
%             text(ca(i,1)-.022,ca(i,3)+0.007,inv_mem_name{domestic(i)},'fontsize', 12)
%         elseif domestic(i) == domestic(41) % Hanmang
%             text(ca(i,1)-0.008,ca(i,3)-0.003,inv_mem_name{domestic(i)},'fontsize', 12)
%         elseif domestic(i) == domestic(17) % SK
%             text(ca(i,1)+.008,ca(i,3),inv_mem_name{domestic(i)},'fontsize', 12)
%         elseif domestic(i) == domestic(1) % Kobyo
%             text(ca(i,1)-0.008,ca(i,3)-0.005,inv_mem_name{domestic(i)},'fontsize', 12)
%         elseif domestic(i) == domestic(35) % NH (Nonghyup)
%             text(ca(i,1)+.005,ca(i,3)+0.002,inv_mem_name{domestic(i)},'fontsize', 12)
%         else
%             text(ca(i,1)+.005,ca(i,3),inv_mem_name{domestic(i)},'fontsize', 12)
%         end
%     end
% end
% box on
% xlim([.018, .141])

subplot(1,2,2)
scatter(cb(:,2),cb(:,3),sqrt(mean_vol_member(foreign)/10^9*5),'filled','markerfacecolor',color_3(3,:),'MarkerEdgeColor',[0 0 0])
hold on
plot(cb(:,2),cb(:,3),'.k')
hold on
plot([.12, .295], [0, 0],'k--')
xlim([.12, .295])
ylim([-.08, .01])
for i = 1:21
    if foreign(i) == foreign(20)
        text(cb(i,2)-.01,cb(i,3)-0.003,inv_mem_name{foreign(i)},'fontsize', 12)
    elseif foreign(i) == foreign(21)
        text(cb(i,2)-.01,cb(i,3)+0.003,inv_mem_name{foreign(i)},'fontsize', 12)
    elseif foreign(i) == foreign(16) % Standard Chartered
        text(cb(i,2)+0.002,cb(i,3),inv_mem_name{foreign(i)},'fontsize', 12)
    elseif foreign(i) == foreign(17) % CIMB
        text(cb(i,2)+.002,cb(i,3)+0.0019,inv_mem_name{foreign(i)},'fontsize', 12)
    elseif foreign(i) == foreign(18) % RBS
        text(cb(i,2)-.01,cb(i,3)+0.001,inv_mem_name{foreign(i)},'fontsize', 12)
    elseif foreign(i) == foreign(19) % Newedge
        text(cb(i,2)-0.01,cb(i,3)-0.002,inv_mem_name{foreign(i)},'fontsize', 12)
    elseif foreign(i) == foreign(2) % Macquarie
        text(cb(i,2)-.024,cb(i,3),inv_mem_name{foreign(i)},'fontsize', 12)
    elseif foreign(i) == foreign(1) % JP Morgan
        text(cb(i,2)-0.01,cb(i,3)+.003,inv_mem_name{foreign(i)},'fontsize', 12)
    elseif foreign(i) == foreign(14) % Daiwa
        text(cb(i,2),cb(i,3)-0.0017,inv_mem_name{foreign(i)},'fontsize', 12)
    else
        text(cb(i,2)+.005,cb(i,3),inv_mem_name{foreign(i)},'fontsize', 12)
    end
end
text(-0.13,1.08,'(b)','Units', 'Normalized', 'VerticalAlignment', 'Top','fontsize',14)
xlabel('\rho_{FR,IV}','fontsize',16)
ylabel('\rho_{IS,IV}','fontsize',16)
%title('Foreign members','fontsize',16)


% axes('Position',[.83 .4 .15/1.3 .3/1.3])
% scatter(ca(:,2),ca(:,3),mean_vol_member(domestic)/10^9*5/5,'filled','markerfacecolor',color_3(2,:),'MarkerEdgeColor',[0 0 0])
% hold on
% scatter(cb(:,2),cb(:,3),mean_vol_member(foreign)/10^9*5/5,'filled','markerfacecolor',color_3(3,:),'MarkerEdgeColor',[0 0 0])
% hold on
% plot([0, 0],[-.4, .23],'k--')
% hold on
% plot([-.45, .4],[0, 0],'k--')
% xlim([-.4, .38])
% ylim([-.35, .23])
% box on


%% Correlation distribution
% domestic-domestic, foreign-foreign, domestic-foreign
% output: [corr_matrix]

corr_matrix = cell(10,1);
for f = 1:1
    corr_matrix{f,1} = cell(75,75);
    for i = 1+(f-1)*121:f*121
        sprintf('%d, %d',f, i)
        for y = 1:11
            for ii = 1:75
                for iii = 1:75
                    if sum(sign(inv_mem_net{i,2}(year_day{y,1},ii))) > length(year_day{y,1})/1.5 ...
                            && sum(sign(inv_mem_net{i,2}(year_day{y,1},iii))) > length(year_day{y,1})/1.5
                        corr_matrix{f,1}{ii,iii} = [corr_matrix{f,1}{ii,iii}; corr(inv_mem_net{i,1}(year_day{y,1},ii),inv_mem_net{i,1}(year_day{y,1},iii))];
                    end
                end
            end
        end
    end
end

d_d = [];
for i = 1:41
    for ii = 1:41
        if i ~= ii
            d_d = [d_d; corr_matrix{1,1}{domestic(i),domestic(ii)}];
        end
    end
end

f_f = [];
for i = 1:21
    for ii = 1:21
        if i ~= ii
            f_f = [f_f; corr_matrix{1,1}{foreign(i),foreign(ii)}];
        end
    end
end

d_f = [];
for i = 1:41
    for ii = 1:21
        d_f = [d_f; corr_matrix{1,1}{domestic(i),foreign(ii)}];
    end
end

[d_d_1, d_d_2, d_d_3, d_d_4] = pcdf(d_d,10000);
[f_f_1, f_f_2, f_f_3, f_f_4] = pcdf(f_f,10000);
[d_f_1, d_f_2, d_f_3, d_f_4] = pcdf(d_f,10000);
[d_d_1_m, d_d_2_m, d_d_3_m, d_d_4_m] = pcdf_minus(d_d,10000);
[f_f_1_m, f_f_2_m, f_f_3_m, f_f_4_m] = pcdf_minus(f_f,10000);
[d_f_1_m, d_f_2_m, d_f_3_m, d_f_4_m] = pcdf_minus(d_f,10000);

figure
set(gcf,'color','w')
plot(d_d_1,d_d_4,'g')
hold on
plot(f_f_1,f_f_4,'b')
hold on
plot(d_f_1,d_f_4,'k')
hold on
plot(d_d_1_m,d_d_4_m,'g--')
hold on
plot(f_f_1_m,f_f_4_m,'b--')
hold on
plot(d_f_1_m,d_f_4_m,'k--')
hold on
plot([0 0],[10^(-5) 1],'k:')
grid on
legend('domestic-domestic','foreign-foreign','domestic-foreign')
set(gca,'yscale','log')


%% Granger Causality test between the inventory variation (+return)
% output: [stat_test]

% causality - yearly
stat_test = cell(10,2);
for f = 1:1
    stat_test{f,1} = zeros(75,75);
    stat_test{f,2} = zeros(75,75);
    for i = 1+(f-1)*121:f*121
        sprintf('%d, %d',f, i)
        for y = 1:11
            for ii = 1:75
                for iii = 1:75
                    if sum(sign(inv_mem_net{i,2}(year_day{y,1},ii))) > length(year_day{y,1})/1.5 ...
                            && sum(sign(inv_mem_net{i,2}(year_day{y,1},iii))) > length(year_day{y,1})/1.5
                        [F,c_v] = granger_cause(inv_mem_net{i,1}(year_day{y,1},ii),inv_mem_net{i,1}(year_day{y,1},iii),.05,1);
                        if F > c_v
                            stat_test{f,1}(ii,iii) = stat_test{f,1}(ii,iii) + 1;
                        end
                        stat_test{f,2}(ii,iii) = stat_test{f,2}(ii,iii) + 1;
                        clear F c_v
                    end
                end
            end
        end
    end
end

% causality - whole period
stat_test = cell(10,2);
for f = 1:1
    stat_test{f,1} = zeros(75,75);
    stat_test{f,2} = zeros(75,75);
    for i = 1+(f-1)*121:f*121
        sprintf('%d, %d',f, i)
        for ii = 1:75
            for iii = 1:75
                if sum(sign(inv_mem_net{i,2}(:,ii))) > 2722/1.5 ...
                        && sum(sign(inv_mem_net{i,2}(:,iii))) > 2722/1.5
                    [F,c_v] = granger_cause(inv_mem_net{i,1}(:,ii),inv_mem_net{i,1}(:,iii),.05,5);
                    if F > c_v
                        stat_test{f,1}(ii,iii) = stat_test{f,1}(ii,iii) + 1;
                    end
                    stat_test{f,2}(ii,iii) = stat_test{f,2}(ii,iii) + 1;
                    clear F c_v
                end
            end
        end
    end
end

% price - inventory variation causality
stat_test = cell(10,2);
for f = 1:1
    stat_test{f,1} = zeros(75,1);
    for i = 1+(f-1)*121:f*121
        sprintf('%d, %d',f, i)
        for iii = 1:75
            if sum(sign(inv_mem_net{i,2}(:,iii))) > 2722/1.5
                [F,c_v] = granger_cause(price_ori(:,i),inv_mem_net{i,1}(1:2722,iii),.05,1);
                if F > c_v
                    stat_test{f,1}(iii,1) = stat_test{f,1}(iii,1) + 1;
                end
                clear F c_v
            end
        end
    end
end

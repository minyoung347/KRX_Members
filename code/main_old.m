%% from dynamics_of_member.m
% forked 2020.01.05 17:00

% import functions
addpath(genpath('/home/minyoung/data_minyoung/Research/Matlab_function/panel/paneldata'))
addpath('/home/minyoung/data_minyoung/Research/Matlab_function/')
addpath('/home/minyoung/data_minyoung/Research/Matlab_function/linspecer')

% load korean market data
load('/home/minyoung/data_minyoung/Research/korean_market/data_200106.mat');

% % add Nation data to domestic institutions
% for firm = 1:1216
%     inv_mem_net{firm,1}(1:1378,3) = inv_mem_net{firm,1}(1:1378,3) + inv_mem_net{firm,1}(1:1378,11);
%     inv_mem_net{firm,2}(1:1378,3) = inv_mem_net{firm,2}(1:1378,3) + inv_mem_net{firm,2}(1:1378,11);
% end

%% KOSPI/KOSDAQ
% output: [name_kospi(kosdaq), code_kospi(kosdaq), ind_kospi(kosdaq)]

% KOSPI/KOSDAQ code
code_temp = table2cell(readtable('stock_code.txt','ReadVariableNames',false));
code_kospi_temp = code_temp(1:1448);
code_kosdaq_temp = code_temp(1449:end);
clear code_temp

% label stocks to KOSPI/KOSDAQ by temp
temp = zeros(1210,1);
for i = 1:1210
    for ii = 1:length(code_kospi_temp)
        if code_firm{i,1} == code_kospi_temp{ii,1}
            temp(i,1) = 1;
        end
    end
    for ii = 1:length(code_kosdaq_temp)
        if code_firm{i,1} == code_kosdaq_temp{ii,1}
            temp(i,1) = 2;
        end
    end
end

% name, code, index of KOSPI/KOSDAQ
name_kospi = name_firm(temp == 1);
name_kosdaq = name_firm(temp == 2);
code_kospi = code_firm(temp == 1);
code_kosdaq = code_firm(temp == 2);
ind_kospi = find(temp == 1);
ind_kosdaq = find(temp == 2);

clear temp code_temp code_kospi_temp code_kosdaq_temp


%% Load Korea 10 years Bond
% output: [koreabond]

koreabond_raw = csvread('korea_10year_bond_yield_processed.csv');
koreabond_raw = flipud(koreabond_raw);

koreabond = zeros(2801,2);

for i = 1:2801
    for ii = 1:length(koreabond_raw)
        if day{1,1}(i,1) == koreabond_raw(ii,1)
            koreabond(i,1) = koreabond_raw(ii,2);
            koreabond(i,2) = koreabond_raw(ii,6);
        end
    end
end

empty_ind = find(koreabond(:,1) == 0);
double_empty_ind = empty_ind(find(empty_ind(2:end)-empty_ind(1:end-1)==1));
koreabond(double_empty_ind,1) = (koreabond(double_empty_ind-1,1) + koreabond(double_empty_ind+2,1))/2;
koreabond(double_empty_ind+1,1) = (koreabond(double_empty_ind-1,1) + koreabond(double_empty_ind+2,1))/2;
empty_ind(empty_ind==double_empty_ind) = [];
empty_ind(empty_ind==double_empty_ind+1) = [];
koreabond(empty_ind,1) = (koreabond(empty_ind-1) + koreabond(empty_ind+1))/2;
koreabond(empty_ind,2) = log(koreabond(empty_ind,1)) - log(koreabond(empty_ind-1,1));

koreabond(:,2) = koreabond(:,2)/100;

clear koreabond_raw empty_ind double_empty_ind


%% Load KOSDAQ Index
% output: [koreabond]

kosdaq_raw = csvread('kosdaq_linux.csv');
kosdaq_raw = flipud(kosdaq_raw);

kosdaq = zeros(2801,1);
kosdaq_return = zeros(2801,1);

for i = 1:2801
    for ii = 1:length(kosdaq_raw)
        if day{1,1}(i,1) == kosdaq_raw(ii,1)
            kosdaq(i,1) = kosdaq_raw(ii,2);
        end
    end
end

kosdaq_return(2:end,1) = log(kosdaq(2:end,1)) - log(kosdaq(1:end-1,1));
kosdaq_return(1,1) = kosdaq_return(2,1);

clear kosdaq_raw

%% 모든 주식종목을 합해서 회원사 레벨
% output: [firm_net_member, firm_net_member_mat]
% firm_net_member is the daily trading amount of each member
% [firm_net_member_mat] is the matrix form of [firm_net_member]

firm_net_member = cell(75,1);
firm_net_member_mat = zeros(2801,1);
for i = 1:75
    firm_net_member{i,1} = zeros(2801,1);
    for ii = 1:121
        firm_net_member{i,1} = firm_net_member{i,1} + firm_net{i,1}(:,ii) .* price{ii,2};
    end
    firm_net_member{i,1} = firm_net_member{i,1} / 121;
    firm_net_member_mat(:,i) = firm_net_member{i,1};
end

figure
set(gcf,'color','w')
for i = 1:21
    subplot(3,7,i)
    plot(cumsum(firm_net_member{foreign(i),1}))
    title(inv_mem_name{foreign(i),1},'fontsize',7)
end

figure
set(gcf,'color','w')
for i = 1:41
    subplot(6,7,i)
    plot(cumsum(firm_net_member{domestic(i),1}))
    title(inv_mem_name{domestic(i),1},'fontsize',7)
end


%% 외국인 비중 데이터 생성 (계산 불필요. 데이터 파일에 반영 됨)
% output: [fr, fr_ratio, fr_high, fr_low, fr_high_lmc]
% fr_ratio: 외국인 평균 비중
% fr_high: 외국인 평균 비중이 높은 주식 1 decile
% fr_low: 외국인 평균 비중이 낮은 주식 10 decile
% fr_high_lmc: 시가총액이 낮은(lmc) 기업 중 외국인 평균 비중이 높은(high) 주식

fr = cell(1216,1);
for i = 1:1216
    fr{i,1} = zeros(2801,2);
    if isfile(sprintf('./foreign/foreign_%s.txt',code_firm{i,1}))
        fid = fopen(sprintf('./foreign/foreign_%s.txt',code_firm{i,1}));
        temp = textscan(fid, '%s %f %f %f','Delimiter',',');
        fclose(fid);
        for ii = 1:2801
            ind = find(temp{1,2} == day{1,1}(ii,1));
            if length(ind) > 0
                fr{i,1}(ii,1) = temp{1,3}(ind);
                fr{i,1}(ii,2) = temp{1,4}(ind);
            end
        end
    end
end
clear fid temp

% Averaged foreign ratio
for i = 1:1216
    fr_ratio(i,1) = mean(delniz(fr{i,1}(:,1)./fr{i,1}(:,2),[1,2],0));
end


fr_high = find(fr_ratio>.23);
fr_low = find(fr_ratio<.0047);

fr_temp= find(fr_ratio>.14855);
fr_high_lmc = fr_temp(find(fr_temp>121));



%% Causality test between the inventory variation
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

%% member trading on/off
% output: [start_end, start_end_ind]
% start_end: start_day, end_day, # of days
% start_end_ind: index of start_end

total_volume = cell(1210,1);
agg_total_volume = zeros(2722,75);
for f = 1:1210
    f
    total_volume{f,1} = zeros(2722,75);
    for i = 1:2722
        for ii = 1:75
            total_volume{f,1}(i,ii) = sign(firm_net{ii,2}(i,f));
        end
    end
    
    agg_total_volume = sign(agg_total_volume + total_volume{f,1});
end


f = 1;
figure
set(gcf,'color','w')
subplot(1,2,1)
imagesc(agg_total_volume(:,[domestic,foreign]))
set(gca,'xtick',1:62,'xticklabel',inv_mem_name([domestic,foreign]),'fontsize',10)
xtickangle(90)
title('All firm')
subplot(1,2,2)
imagesc(total_volume{f,1}(:,[domestic,foreign]))
set(gca,'xtick',1:62,'xticklabel',inv_mem_name([domestic,foreign]),'fontsize',10)
xtickangle(90)
title(sprintf('%s',name_firm{f,1}))

figure
set(gcf,'color','w')
for i = 1:40
    subplot(4,10,i)
    imagesc(total_volume{i+0,1}(:,[domestic,foreign]))
    title(sprintf('%s',name_firm{i+0,1}))
end

start_end = zeros(75,3);
start_end_ind = zeros(75,2);
for i = 1:75
    start_end(i,1) = day{1,1}(min(find(agg_total_volume(:,i) == 1)));
    start_end(i,2) = day{1,1}(max(find(agg_total_volume(:,i) == 1)));
    start_end(i,3) = length(find(agg_total_volume(:,i) == 1));
    
    start_end_ind(i,1) = min(find(agg_total_volume(:,i) == 1));
    start_end_ind(i,2) = max(find(agg_total_volume(:,i) == 1));
end

clear total_volume agg_total_volume

%% statistics for member firms
% Market capitalization 1 decile stocks에 대해서, member trading 통계
% output: [total_pvol, active_day, active_stock]
% total_pvol: 전체기간에 대해서 회원사가 1 decile stocks을 거래한 pvolume
% active_day: 회원사가 활동한 기간
% active_stock: 회원사가 거래한 주식 수

total_pvol = zeros(1,75);
active_day = zeros(1,75);
active_stock = zeros(1,75);

for ii = 1:75
    active_temp = zeros(2722,121);
    for i = 1 : 121
        pvol = inv_mem_net{i,2}(1:2722,ii) .* price{i,2}(1:2722);
        total_pvol(1,ii) = total_pvol(1,ii) + sum(pvol);
        
        active_temp(:,i) = inv_mem_net{i,2}(1:2722,ii)>0;
    end
    active_stock(1,ii) = length(find(sum(active_temp,1) > 0));
    active_day(1,ii) = length(find(sum(active_temp,2) > 0));
end

figure
set(gcf,'color','w')
subplot(1,2,1)
bar(total_pvol(domestic))
ylim([0 8*10^15])
set(gca,'xtick',1:41,'xticklabel',inv_mem_name(domestic),'fontsize',9)
xtickangle(90)
subplot(1,2,2)
bar(total_pvol(foreign))
ylim([0 8*10^15])
set(gca,'xtick',1:21,'xticklabel',inv_mem_name(foreign),'fontsize',9)
xtickangle(90)

figure
set(gcf,'color','w')
subplot(1,2,1)
bar(total_pvol(domestic)./active_day(domestic))
ylim([0 3*10^12])
set(gca,'xtick',1:41,'xticklabel',inv_mem_name(domestic),'fontsize',9)
xtickangle(90)
subplot(1,2,2)
bar(total_pvol(foreign)./active_day(foreign))
ylim([0 3*10^12])
set(gca,'xtick',1:21,'xticklabel',inv_mem_name(foreign),'fontsize',9)
xtickangle(90)

figure
set(gcf,'color','w')
subplot(1,2,1)
bar(total_pvol(domestic)./active_day(domestic)./active_stock(domestic))
ylim([0 2.5*10^9])
set(gca,'xtick',1:41,'xticklabel',inv_mem_name(domestic),'fontsize',9)
xtickangle(90)
subplot(1,2,2)
bar(total_pvol(foreign)./active_day(foreign)./active_stock(foreign))
ylim([0 2.5*10^9])
set(gca,'xtick',1:21,'xticklabel',inv_mem_name(foreign),'fontsize',9)
xtickangle(90)

a2 = inv_mem_name([domestic,foreign]);
a3 = start_end([domestic,foreign],1:2);
a3p = active_day([domestic,foreign]);
a4 = active_stock([domestic,foreign]);
a5 = total_pvol([domestic,foreign]);
a6 = total_pvol([domestic,foreign])./active_day([domestic,foreign])./active_stock([domestic,foreign]);

fid = fopen('stat_member.csv','w');
for i = 1:62
    start_ym = num2str(a3(i,1));
    end_ym = num2str(a3(i,2));
    fprintf(fid,'%.0f,%s,%s-%s (%.0f),%.0f,%.2f,%.0f\n',i,a2{i,1},start_ym(3:6),end_ym(3:6),a3p(i),a4(i),a5(i)/10^12,a6(i)/10^6);
end
fclose(fid);


%% correlation distribution
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


%% correlation between inventory variation & return
% delay에 따라서 의미 달라짐 (delay = -1, 0, 1)

% cir: Correlation b/w Iventory & Return
cir_1 = cell(10,75);
cir_2 = cell(10,75);
cir_3 = cell(10,75);

for ii = 1:75
    sprintf('%d',ii)
    for f = 1:10
        temp_1 = [];
        temp_2 = [];
        temp_3 = [];
        for i = 1+121*(f-1):121*f
            for y = 1:11
                if length(inv_mem_net{i,1}(year_day{y,1},ii)) > length(year_day{y,1})/1.5
                    temp_1 = [temp_1; corr(inv_mem_net{i,1}(year_day{y,1},ii), price_return{i,1}(year_day{y,1},1))];
                    temp_2 = [temp_2; corr(inv_mem_net{i,1}(year_day{y,1}(1:end-1),ii), price_return{i,1}(year_day{y,1}(2:end),1))];
                    temp_3 = [temp_3; corr(inv_mem_net{i,1}(year_day{y,1}(2:end),ii), price_return{i,1}(year_day{y,1}(1:end-1),1))];
                end
            end
        end
        cir_1{f,ii} = temp_1;
        cir_2{f,ii} = temp_2;
        cir_3{f,ii} = temp_3;
    end
end

% reversed: cir_1{f,ii}, cc_1(ii,f)
for ii = 1:75
    for f = 1:10
        cc_1(ii,f) = mean(delniz(cir_1{f,ii},[1,2],0));
        cc_2(ii,f) = mean(delniz(cir_2{f,ii},[1,2],0));
        cc_3(ii,f) = mean(delniz(cir_3{f,ii},[1,2],0));
    end
end


figure
set(gcf,'color','w')
for i = 1:41
    subplot(6,7,i)
%     subplot(1,3,i)
    histfit(cir_1{1,domestic_sort(i)}) % histogram of cir for 1 decile and domestic members
    hold on
    plot([0 0],[0 100],'k','linewidth',3)
    xlim([-.5 .5])
    title(sprintf('%s, %.2d',inv_mem_name{domestic_sort(i,1),1},mean(delniz(cir_1{1,domestic_sort(i,1)},[1,2],0))))
end


figure
set(gcf,'color','w')
for i = 1:21
    subplot(3,7,i)
    histfit(cir_1{1,foreign_sort(i,2)}) % histogram of cir for 1 decile and foreign members
    hold on
    plot([0 0],[0 100],'k','linewidth',3)
    xlim([-.5 .5])
    title(sprintf('%s, %.2d',inv_mem_name{foreign_sort(i,2),1},mean(delniz(cir_1{1,foreign_sort(i,2)},[1,2],0))))
end


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


ind_target = ind_kospi;
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


%% mean correlation matrix for network construction
% 논문에서 사용한 네트워크 구성 방법
% 1 decile 주식들에 대해서 회원사들 사이의 correlation 합
% ct: 거래가 있는 investor(75)
% cr: correlation 합
% ct_m: ct에서 회원사만 (+정렬)
% cr_m: cr에서 회원사만 (+정렬)
% cr_m에서 거래가 적은 회원사 제거
% 네트워크 node: mean_vol_member
% 네트워크 edge: cr_m

ct = zeros(75,75);
cr = zeros(75,75);

f = 1;

for i = 1+(f-1)*121:121*f
    for y = 1:11
        temp = corr(inv_mem_net{i,1}(year_day{y,1},:));
        temp(isnan(temp)) = 0;
        temp(isinf(temp)) = 0;
        
        no_member_temp = [];
        for ii = 1:75
            if length(find(inv_mem_net{i,2}(year_day{y,1},ii)>0)) < length(year_day{y,1})/3
                no_member_temp = [no_member_temp; ii];
            end
        end
        
        temp(no_member_temp,:) = 0;
        temp(:,no_member_temp) = 0;
        
        ct = ct + abs(sign(temp));
        cr = cr + temp;
    end
end

for ii = 1:75
    cr(ii,ii) = 0;
end

% 75 investor 중에서 회원사 추출 및 정렬
df = [domestic_sort(:,1)', foreign];
ct_m = ct(df,df);
cr_m = cr(df,df) ./ ct(df,df);

% 거래가 적은 회원사 제거
cr_m(find(mean(ct_m,2)<360),:) = 0;
cr_m(:,find(mean(ct_m,2)<360)) = 0;
cr_m(isnan(cr_m)) = 0;
cr_m_valid = find(mean(cr_m)~=0);



figure
set(gcf,'color','w')
imagesc(cr_m(cr_m_valid,cr_m_valid),[-.2,.2])
set(gca,'xtick',1:length(cr_m_valid),'xticklabel',inv_mem_name(df(cr_m_valid)),'fontsize',9)
set(gca,'ytick',1:length(cr_m_valid),'yticklabel',inv_mem_name(df(cr_m_valid)),'fontsize',9)
xtickangle(90)
colorbar
colormap(jet)


sqrt(mean_vol_member([domestic,foreign]))
% save node
fid = fopen('./save_network/inv_node_mean.csv','w');
fprintf(fid,'id,label,color,size\n');
for i = 1:size(domestic,2)
    fprintf(fid,'%d,%s,%d,%d\n',domestic(i),inv_mem_name{domestic(i),1},1,2*sqrt(mean_vol_member(domestic(i)))/10^5);
end
for i = 1:size(foreign,2)
    fprintf(fid,'%d,%s,%d,%d\n',foreign(i),inv_mem_name{foreign(i),1},2,2*sqrt(mean_vol_member(foreign(i)))/10^5);
end
fclose(fid);


% save edge
fid = fopen('./save_network/inv_edge_mean.csv','w');
fprintf(fid,'source,target,weight\n');
% member = [domestic,foreign];
% corr_mean = cr([domestic,foreign],[domestic,foreign])./ct([domestic,foreign],[domestic,foreign]);
% corr_mean(isnan(corr_mean)) = 0;
for i = 1:62-1
    for ii = i+1:62
        fprintf(fid,'%d,%d,%d\n',df(i),df(ii),cr_m(i,ii));
    end
end
fclose(fid);




%% Herding of members


num_firm = 1210;

% corr_fr = zeros(num_firm,length(year_day));
% leval = zeros(num_firm,length(year_day));
% lambda_max_save = zeros(num_firm,length(year_day)-1);

herd_adf = cell(num_firm,3); % for herding distribution, herding all, domestic, foreign
% 'pnh' is the existence of herding
% 'pnhd' is the direction of herding
pnh = cell(num_firm,length(year_day)-1);
pnhd = cell(num_firm,length(year_day)-1);
pnh_d = cell(num_firm,length(year_day)-1);
pnhd_d = cell(num_firm,length(year_day)-1);
pnh_f = cell(num_firm,length(year_day)-1);
pnhd_f = cell(num_firm,length(year_day)-1);

pnh_d_id = cell(num_firm,length(year_day)-1);
pnh_d_is = cell(num_firm,length(year_day)-1);

pnh_raw = cell(num_firm,length(year_day)-1);
pnh_raw_d = cell(num_firm,length(year_day)-1);
pnh_raw_f = cell(num_firm,length(year_day)-1);
pnh_raw_d_id = cell(num_firm,length(year_day)-1);
pnh_raw_d_is = cell(num_firm,length(year_day)-1);

color_3 = linspecer(3);


for f = 1:num_firm
    f
    for y = 1:length(year_day)-1
        onoff_trading = find(vol{f,1}(year_day{y,1}) > 0);
        
        active_day = sum(sign(inv_mem_net{f,2}(year_day{y,1},:)),1);
        active_inv = find(active_day > length(onoff_trading)/1.5);
        clear active day onoff_trading
        
        domestic_act = intersect(domestic, active_inv);
        foreign_act = intersect(foreign, active_inv);
        member_act = [domestic_act, foreign_act];
        
        domestic_act_is = intersect(domestic(ca(:,1)<0), active_inv);
        domestic_act_id = intersect(domestic(ca(:,1)>0), active_inv);
        clear active_inv
        
        if length(domestic_act) + length(foreign_act) > 0
%             corr_df = corr(inv_mem_net{f,1}(year_day{y,1},member_act));
% 
%             [eig_vec, eig_val] = eig(corr_df);
%             [eig_val_sort, ind_val_sort] = sort(on_diagonal(eig_val),'descend');
%             largest_eig_val = eig_val_sort(1);
%             largest_eig_vec = eig_vec(:,ind_val_sort(1));
%             
%             Q = size(year_day{y,1},1) / size(member_act,2);
%             lambda_min = 1 + 1/Q - 2*sqrt(1/Q);
%             lambda_max = 1 + 1/Q + 2*sqrt(1/Q);
%             lambda_max_save(f,y) = lambda_max;
            
%             figure1 = figure;
%             set(gcf,'color','w')
%             subplot(1,2,1)
%             x=linspace(lambda_min,lambda_max,10000); 
%             histo = histogram(eig_val_sort,30);
%             histo.Visible = 'off';
%             bar((histo.BinEdges(1:end-1) + histo.BinEdges(2:end))/2,histo.BinCounts/8)
%             hold on
%             plot(x,Q/(2*pi)*sqrt((lambda_max-x).*(x-lambda_min))./x,'k--','Linewidth',1.5)
%             xlabel('\lambda','fontsize',16)
%             ylabel('\rho(\lambda)','fontsize',16)
%             legend('dataset','random matrix')
%             legend boxoff
%             text(2.8, 0.3,'largest eigenvalue(\lambda_1)')
%             ylim([0 .92])
%             text(-0.14,1.08,'(a)','Units', 'Normalized', 'VerticalAlignment', 'Top')
            
            
%             if sum(member_act==46) > 0
            
%             end
                 
            % 'x1' is the sum of direction in daily
            x1 = sum(sign(inv_mem_net{f,1}(year_day{y,1},member_act)),2);
            % 'x1_1' is the sum of positive direction in daily
            x1_1 = sum(sign(inv_mem_net{f,1}(year_day{y,1},member_act))>0,2);
            % 'x2' is the number of member_act in daily
            x2 = sum(sign(inv_mem_net{f,2}(year_day{y,1},member_act)),2);
            % 'x3' is the price return in daily
            x3 = price_return{f,1}(year_day{y,1});

            % 'x1_d' is the sum of direction of [domestic member] in daily
            x1_d = sum(sign(inv_mem_net{f,1}(year_day{y,1},domestic_act)),2);
            % 'x1_d_1' is the sum of positive direction of [domestic member] in daily
            x1_d_1 = sum(sign(inv_mem_net{f,1}(year_day{y,1},domestic_act))>0,2);
            % 'x2_d' is the number of member_act of [domestic member] in daily
            x2_d = sum(sign(inv_mem_net{f,2}(year_day{y,1},domestic_act)),2);

            % 'x1_f' is the sum of direction of [foreign member] in daily
            x1_f = sum(sign(inv_mem_net{f,1}(year_day{y,1},foreign_act)),2);
            % 'x1_f_1' is the sum of positive direction of [foreign member] in daily
            x1_f_1 = sum(sign(inv_mem_net{f,1}(year_day{y,1},foreign_act))>0,2);
            % 'x2_f' is the number of member_act of [foreign member] in daily
            x2_f = sum(sign(inv_mem_net{f,2}(year_day{y,1},foreign_act)),2);
            
            % 'x1_d' is the sum of direction of [domestic member similar to institution] in daily
            x1_d_is = sum(sign(inv_mem_net{f,1}(year_day{y,1},domestic_act_is)),2);
            % 'x1_d_1' is the sum of positive direction of [domestic member similar to institution] in daily
            x1_d_1_is = sum(sign(inv_mem_net{f,1}(year_day{y,1},domestic_act_is))>0,2);
            % 'x2_d' is the number of member_act of [domestic member similar to institution] in daily
            x2_d_is = sum(sign(inv_mem_net{f,2}(year_day{y,1},domestic_act_is)),2);
            
            % 'x1_d' is the sum of direction of [domestic member similar to individual] in daily
            x1_d_id = sum(sign(inv_mem_net{f,1}(year_day{y,1},domestic_act_id)),2);
            % 'x1_d_1' is the sum of positive direction of [domestic member similar to individual] in daily
            x1_d_1_id = sum(sign(inv_mem_net{f,1}(year_day{y,1},domestic_act_id))>0,2);
            % 'x2_d' is the number of member_act of [domestic member similar to individual] in daily
            x2_d_id = sum(sign(inv_mem_net{f,2}(year_day{y,1},domestic_act_id)),2);

            pnh_raw{f,y} = (x1_1 + (x1_1 - x2)) ./ x2;
            pnh_raw_d{f,y} = (x1_d_1 + (x1_d_1 - x2_d)) ./ x2_d;
            pnh_raw_f{f,y} = (x1_f_1 + (x1_f_1 - x2_f)) ./ x2_f;
            
            pnh_raw_d_is{f,y} = (x1_d_1_is + (x1_d_1_is - x2_d_is)) ./ x2_d_is;
            pnh_raw_d_id{f,y} = (x1_d_1_id + (x1_d_1_id - x2_d_id)) ./ x2_d_id;

            herd_adf{f,1} = [herd_adf{f,1}; x1./x2];
            herd_adf{f,2} = [herd_adf{f,2}; x1_d./x2_d];
            herd_adf{f,3} = [herd_adf{f,3}; x1_f./x2_f];

            % herding using binomial null hypothesis
            % 'pnh' is the existence of herding
            % 'pnhd' is the direction of herding
            for bn = 1:size(year_day{y,1},1)
                if nchoosek(x2(bn,1),x1_1(bn,1))*(1/2)^(x2(bn,1)) <= 0.05
                    pnh{f,y}(bn,1) = 1;
                    if x1_1(bn,1) < x2(bn,1)/2
                        pnhd{f,y}(bn,1) = -1;
                    else
                        pnhd{f,y}(bn,1) = 1;
                    end
                else
                    pnh{f,y}(bn,1) = 0;
                    pnhd{f,y}(bn,1) = 0;
                end

                if nchoosek(x2_f(bn,1),x1_f_1(bn,1))*(1/2)^(x2_f(bn,1)) <= 0.05
                    pnh_f{f,y}(bn,1) = 1;
                    if x1_f_1(bn,1) < x2_f(bn,1)/2
                        pnhd_f{f,y}(bn,1) = -1;
                    else
                        pnhd_f{f,y}(bn,1) = 1;
                    end
                else
                    pnh_f{f,y}(bn,1) = 0;
                    pnhd_f{f,y}(bn,1) = 0;
                end

                if nchoosek(x2_d(bn,1),x1_d_1(bn,1))*(1/2)^(x2_d(bn,1)) <= 0.05
                    pnh_d{f,y}(bn,1) = 1;
                    if x1_d_1(bn,1) < x2_d(bn,1)/2
                        pnhd_d{f,y}(bn,1) = -1;
                    else
                        pnhd_d{f,y}(bn,1) = 1;
                    end
                else
                    pnh_d{f,y}(bn,1) = 0;
                    pnhd_d{f,y}(bn,1) = 0;
                end
                
                if nchoosek(x2_d_is(bn,1),x1_d_1_is(bn,1))*(1/2)^(x2_d_is(bn,1)) <= 0.05
                    pnh_d_is{f,y}(bn,1) = 1;
                    if x1_d_1_is(bn,1) < x2_d_is(bn,1)/2
                        pnhd_d_is{f,y}(bn,1) = -1;
                    else
                        pnhd_d_is{f,y}(bn,1) = 1;
                    end
                else
                    pnh_d_is{f,y}(bn,1) = 0;
                    pnhd_d_is{f,y}(bn,1) = 0;
                end
                
                if nchoosek(x2_d_id(bn,1),x1_d_1_id(bn,1))*(1/2)^(x2_d_id(bn,1)) <= 0.05
                    pnh_d_id{f,y}(bn,1) = 1;
                    if x1_d_1_id(bn,1) < x2_d_id(bn,1)/2
                        pnhd_d_id{f,y}(bn,1) = -1;
                    else
                        pnhd_d_id{f,y}(bn,1) = 1;
                    end
                else
                    pnh_d_id{f,y}(bn,1) = 0;
                    pnhd_d_id{f,y}(bn,1) = 0;
                end
            end
        end
    end
end


% 'pnh' is the existence of herding
% 'pnhd' is the direction of herding

% 'pnhh' is number of day of herding
% 'pnhhp' is number of day of buy herding
% 'pnhhn' is number of day of sell herding
% 'pnhh_d' is number of day of domestic member herding
% 'pnhh_f' is number of day of foreign member herding
for f = 1:num_firm
    for y = 1:size(year_day,1)-1
        pnhh(f,y) = sum(pnh{f,y}) ./ size(year_day{y,1},1);
        pnhdp(f,y) = sum(pnhd{f,y}>0) ./  size(year_day{y,1},1);
        pnhdn(f,y) = sum(pnhd{f,y}<0) ./  size(year_day{y,1},1);
        
        pnhh_d(f,y) = sum(pnh_d{f,y}) ./ size(year_day{y,1},1);
        pnhh_f(f,y) = sum(pnh_f{f,y}) ./ size(year_day{y,1},1);
        
        pnhh_d_id(f,y) = sum(pnh_d_id{f,y}) ./ size(year_day{y,1},1);
        pnhh_d_is(f,y) = sum(pnh_d_is{f,y}) ./ size(year_day{y,1},1);
    end
end


% herding, price
for f = 1:num_firm
    for y = 1:size(year_day,1)-1
        if ~isempty(pnhd{f,y})
            cdr(f,y) = corr(pnh_raw{f,y},price_return{f,1}(year_day{y,1}));
        end
        if ~isempty(pnhd_d{f,y})
            cdr_d(f,y) = corr(pnh_raw_d{f,y},price_return{f,1}(year_day{y,1}));
        end
        if ~isempty(pnhd_f{f,y})
            cdr_f(f,y) = corr(pnh_raw_f{f,y},price_return{f,1}(year_day{y,1}));
        end
    end
end


for f = 1:10
    pnhh_mean(f,1) = mean(mean(pnhh(1+(f-1)*120:f*120,:)));
    pnhh_d_mean(f,1) = mean(mean(pnhh_d(1+(f-1)*120:f*120,:)));
    pnhh_f_mean(f,1) = mean(mean(pnhh_f(1+(f-1)*120:f*120,:)));
    
    pnhh_d_is_mean(f,1) = mean(mean(pnhh_d_is(1+(f-1)*120:f*120,:)));
    pnhh_d_id_mean(f,1) = mean(mean(pnhh_d_id(1+(f-1)*120:f*120,:)));
    
    temp = [];
    temp_d = [];
    temp_f = [];
    for i = 1+(f-1)*120:f*120
        temp = [temp; delniz(cdr(i,:),[1,2],0)];
        temp_d = [temp_d; delniz(cdr_d(i,:),[1,2],0)];
        temp_f = [temp_f; delniz(cdr_f(i,:),[1,2],0)];
    end
    cdr_mean(f,1) = mean(temp);
    cdr_d_mean(f,1) = mean(temp_d);
    cdr_f_mean(f,1) = mean(temp_f);
end



% first figure in the paper
figure1 = figure;
set(gcf,'color','w')
% subplot(1,2,1)
plot(pnhh_d_is_mean,'color','k','linewidth',2)
hold on
plot(pnhh_d_id_mean,'color','b','linewidth',2)
hold on
plot(pnhh_mean,'color',color_3(1,:),'linewidth',2)
hold on
plot(pnhh_d_mean,'color',color_3(2,:),'linewidth',2)
hold on
plot(pnhh_f_mean,'color',color_3(3,:),'linewidth',2)
xlim([1 10])
xlabel('Market capitalization decile')
ylabel('<H>')
legend('Domestic (Institutions)','Domestic (Individuals)','All members','Domestic members','Foreign members')
legend boxoff
text(-0.14,1.1,'(a)','Units', 'Normalized', 'VerticalAlignment', 'Top')

subplot(1,2,2)
plot(cdr_mean,'color',color_3(1,:),'linewidth',2)
hold on
plot(cdr_d_mean,'color',color_3(2,:),'linewidth',2)
hold on
plot(cdr_f_mean,'color',color_3(3,:),'linewidth',2)
hold on
plot([1 10],[0 0],'--k')
xlim([1 10])
xlabel('Market capitalization decile')
ylabel('<DH>')
legend('All members','Domestic members','Foreign members')
legend boxoff
text(-0.14,1.1,'(b)','Units', 'Normalized', 'VerticalAlignment', 'Top')

dx=0.06;
dy=0.06;
x=(1-4*dx)/1.9;
y=(1-4*dy)/1.05;
dxx=0.01;
dyy=0.01;
AxesHandle=findobj(figure1,'Type','axes');
set(AxesHandle(2),'Position',[dx+0.01,3*dy,x,y]);
set(AxesHandle(1),'Position',[x+2.6*dx,3*dy,x,y]);
set(gcf,'Position',[100, 100, 800, 400])


%% panel regression all (herding)

y = [];
id = [];
day_tot = [];
day_each = (1:2722)';
year_reg = zeros(2722,1);
ind_start = 1;
ind_end = length(year_day{1,1});
year_temp = 2007;
for i = 1:11
    year_reg(ind_start:ind_end) = year_temp;
    ind_start = ind_start + length(year_day{i,1});
    ind_end = ind_end + length(year_day{i+1,1});
    year_temp = year_temp + 1;
end

x0 = [];
x1 = [];
x2 = [];
x3 = [];
x_year = [];
const = [];
for i = 1:121
    ind = 0;
    
    for year_ind = 1:11
        if length(pnh_raw_d{i,year_ind}) + length(pnh_raw_f{i,year_ind}) > 0
            ind = ind + 1;
        end
    end
    
    if ind == 11
        for year_ind = 1:11
            temp_d_is = pnh_raw_d_is{i,year_ind};
            temp_d_is(isnan(temp_d_is)) = 0;
            temp_d_is(isinf(temp_d_is)) = 0;
            x1 = [x1; temp_d_is];
            
            temp_d_id = pnh_raw_d_id{i,year_ind};
            temp_d_id(isnan(temp_d_id)) = 0;
            temp_d_id(isinf(temp_d_id)) = 0;
            x2 = [x2; temp_d_id];

            temp_f = pnh_raw_f{i,year_ind};
            temp_f(isnan(temp_f)) = 0;
            temp_f(isinf(temp_f)) = 0;
            x3 = [x3; temp_f];
        end
    
        y = [y; price_return{i,1}(1:2722)-koreabond(1:2722,2)];
        day_tot = [day_tot; day_each];
        id = [id; i*ones(2722,1)];
        x_year = [x_year; year_reg];
        const = [const; ones(2722,1)];

        x0 = [x0; kospi_return(1:2722)-koreabond(1:2722,2)];
    end
end

const = randn(length(x0),1);
X = [x_year, x0, x1, x2, x3, const];


csvwrite('panel_data_herding_decile_1.csv',[id,day_tot,X,y])



%% Calculating Trending 
% Corr of ratio & net_buy with price_return

member_corr_net = cell(10,75);
member_corr_ratio = cell(10,75);

for f = 1:10
    sprintf('%d',f)
    for i = (f-1)*121+1 : f*121
        for y = 1:11
            for ii = 1:75
                if sum(sign(inv_mem_net{i,2}(year_day{y,1},ii))) > size(year_day{y,1},1)/1.5
                    % for net buy
                    member_corr_net{f,ii} = [member_corr_net{f,ii}; corr(price_return{i,1}(year_day{y,1}), inv_mem_net{i,1}(year_day{y,1},ii))];

                    % for ratio
                    temp = inv_mem_net{i,1}(year_day{y,1},ii)./inv_mem_net{i,2}(year_day{y,1},ii);
%                     [res, ind] = delniz(temp,[1,2],1);
%                     if length(res) > 3
%                         member_corr_ratio{f,ii} = [member_corr_ratio{f,ii}; corr(price_return{i,1}(year_day{y,1}(ind)), res)];
%                     end
                    temp(isnan(temp)) = 0;
                    temp(isinf(temp)) = 0;
                    member_corr_ratio{f,ii} = [member_corr_ratio{f,ii}; corr(price_return{i,1}(year_day{y,1}), temp)];
                end
            end
        end
    end
end

%% Correlation between inventory according to trending

% mean & length over 121 firms and 11 years
for i = 1:41
    domestic_corr_mean(i,1) = mean(member_corr_ratio{1,domestic(i)});
    domestic_corr_length(i,1) = length(member_corr_ratio{1,domestic(i)});
end
for i = 1:21
    foreign_corr_mean(i,1) = mean(member_corr_ratio{1,foreign(i)});
    foreign_corr_length(i,1) = length(member_corr_ratio{1,foreign(i)});
end

domestic_remain = domestic(domestic_corr_length >= 100);
foreign_remain = foreign(foreign_corr_length >= 100);


[df_sort, df_sort_i] = sort([domestic_corr_mean;foreign_corr_mean]);
df = [domestic,foreign];
df_sort_ind = df(df_sort_i);

df_sort_remain = [];
for i = 1:size(df_sort_ind,2)
    if sum(domestic_remain == df_sort_ind(i)) || sum(foreign_remain == df_sort_ind(i))
        df_sort_remain = [df_sort_remain; df_sort_ind(i)];
    end
end

% 'a': the averaged correlation over 1 decile and 11 years
a = zeros(75,75);
for i = 1:121
    for y = 1:11
        temp = corr(inv_mem_net{i,1}(year_day{y,1},:));
        temp(isnan(temp)) = 0;
        temp(isinf(temp)) = 0;
        a = a + temp;
    end
end

for i = 1:75
    a(i,i) = 0;
end

figure
set(gcf,'color','w')
imagesc(a(df_sort_remain,df_sort_remain)/(121*11))
colormap(jet)
set(gca,'xtick',1:length(df_sort_remain),'xticklabel',inv_mem_name(df_sort_remain),'fontsize',9)
xtickangle(90)
set(gca,'ytick',1:length(df_sort_remain),'yticklabel',inv_mem_name(df_sort_remain),'fontsize',9)
colorbar
title('Domestic + Foreign')





[d_sort, d_sort_i] = sort(domestic_corr_mean);
d_sort_ind = domestic(d_sort_i);

d_sort_remain = [];
for i = 1:size(d_sort_ind,2)
    if sum(domestic_remain == d_sort_ind(i))
        d_sort_remain = [d_sort_remain; d_sort_ind(i)];
    end
end

[f_sort, f_sort_i] = sort(foreign_corr_mean);
f_sort_ind = foreign(f_sort_i);

f_sort_remain = [];
for i = 1:size(f_sort_ind,2)
    if sum(foreign_remain == f_sort_ind(i))
        f_sort_remain = [f_sort_remain; f_sort_ind(i)];
    end
end



figure
set(gcf,'color','w')

subplot(1,2,1)
imagesc(a(d_sort_remain,d_sort_remain)/(121*11))
colormap(jet)
set(gca,'xtick',1:length(d_sort_remain),'xticklabel',inv_mem_name(d_sort_remain),'fontsize',9)
xtickangle(90)
set(gca,'ytick',1:length(d_sort_remain),'yticklabel',inv_mem_name(d_sort_remain),'fontsize',9)
caxis([-.14 .2])
colorbar
title('Domestic')

subplot(1,2,2)
imagesc(a(f_sort_remain,f_sort_remain)/(121*11))
colormap(jet)
set(gca,'xtick',1:length(f_sort_remain),'xticklabel',inv_mem_name(f_sort_remain),'fontsize',9)
xtickangle(90)
set(gca,'ytick',1:length(f_sort_remain),'yticklabel',inv_mem_name(f_sort_remain),'fontsize',9)
caxis([-.14 .2])
colorbar
title('Foreign')

%% calculating directionality

ratio_dist = cell(10,75);
for f = 1:10
    for i = (f-1)*121+1 : f*121
        for ii = 1:75
            temp = firm_net{ii,1}(:,i)./firm_net{ii,2}(:,i);
            temp = delniz(temp,[1,2],0);
            ratio_dist{f,ii} = [ratio_dist{f,ii}; temp];
        end
    end
end



threshold = .1:.1:.9;
% lr_indv = zeros(10,length(threshold));
% lr_fore = zeros(10,length(threshold));
% lr_inst = zeros(10,length(threshold));
directionality_member = cell(75,1);
for i = 1:75
    directionality_member{i,1} = zeros(10,length(threshold));
end
for f = 1:10
    for th = 1:length(threshold)
%         for ii = 1:3
%             if ii == 1
%                 lr_indv(f,th) = length(find(abs(ratio_dist{f,ii})>=threshold(th))) / length(ratio_dist{f,ii});
%             elseif ii == 2
%                 lr_fore(f,th) = length(find(abs(ratio_dist{f,ii})>=threshold(th))) / length(ratio_dist{f,ii});
%             elseif ii == 3
%                 lr_inst(f,th) = length(find(abs(ratio_dist{f,ii})>=threshold(th))) / length(ratio_dist{f,ii});
%             end
%         end
        for ii = 1:75
            directionality_member{ii,1}(f,th) = length(find(abs(ratio_dist{f,ii})>=threshold(th))) / length(ratio_dist{f,ii});
        end
    end
end
clear ratio_dist


figure
set(gcf,'color','w')
subplot(1,3,1)
plot(directionality_member{1,1})
ylim([0 1])
xlabel('percentile')
ylabel('portion of ratio > th')
title('indv')
subplot(1,3,2)
plot(directionality_member{2,1})
ylim([0 1])
xlabel('percentile')
ylabel('portion of ratio > th')
title('fore')
subplot(1,3,3)
plot(directionality_member{3,1})
ylim([0 1])
xlabel('percentile')
ylabel('portion of ratio > th')
title('inst')


%% directionality & trending for members

for i = 1:10
    for ii = 1:75
        trending(i,ii) = mean(delniz(member_corr_ratio{i,ii},[1,2],0));
%         trending(i,ii) = mean(delniz(member_corr_net{i,ii},[1,2],0));
    end
end

k = 1;

addpath('/home/minyoung/data_minyoung/Research/Matlab_function/linspecer')
color_3 = linspecer(3);



% directionality & trending for 41 domestic and 21 foreign
figure2 = figure;
set(gcf,'color','w')
subplot(1,2,1)
for ii = 1:41
    plot(trending(:,domestic(ii)),directionality_member{domestic(ii),1}(:,k),'--k','linewidth',1.5)
    hold on
end
plot([0 0],[0 1],'k--')
hold on
plot([-1 1],[.5 .5],'k--')
hold on
% xlim([-.5 .5])
% ylim([.4 1])
for ii = 1:41
    for i = 1:10
        scatter(trending(i,domestic(ii)),directionality_member{domestic(ii),1}(i,k),70,'filled','MarkerFaceColor',[i/12 i/12 i/12])
        hold on
    end
end

subplot(1,2,2)
for ii = 1:21
    plot(trending(:,foreign(ii)),directionality_member{foreign(ii),1}(:,k),'b','linewidth',1.5)
    hold on
end
plot([0 0],[0 1],'k--')
hold on
plot([-1 1],[.5 .5],'k--')
hold on
for ii = 1:21
    for i = 1:10
        scatter(trending(i,foreign(ii)),directionality_member{foreign(ii),1}(i,k),70,'filled','MarkerFaceColor',[i/12 i/12 i/12])
%         scatter(trending(i,foreign(ii)),directionality_member{foreign(ii),1}(i,k),70,'filled','MarkerFaceColor',[0 0 1])
        hold on
    end
end
% text(-0.13,1.08,'(a)','Units', 'Normalized', 'VerticalAlignment', 'Top')
xlim([-.6 .35])
xlabel('{\it <T>}','fontsize',12)
ylabel('{\it <D>}','fontsize',12)
lgd1 = legend('Individual','Foreigner','Institution');
lgd1.FontSize = 12;
legend boxoff



%% figure in the paper

figure1 = figure;
set(gcf,'color','w')

subplot(1,3,1)
plot(trending(:,1),directionality_member{1,1}(:,k),'-','color',color_3(:,3),'linewidth',1.5)
hold on
plot(trending(:,2),directionality_member{2,1}(:,k),'--','color',color_3(:,2),'linewidth',1.5)
hold on
plot(trending(:,3),directionality_member{3,1}(:,k),':','color',color_3(:,1),'linewidth',1.5)
hold on
plot([0 0],[0 1],'k--')
hold on
plot([-.6 .35],[.5 .5],'k--')
for i = 1:10
    hold on
    scatter(trending(i,1),directionality_member{1,1}(i,k),70,'filled','MarkerFaceColor',[i/12 i/12 i/12])
    hold on
    scatter(trending(i,2),directionality_member{2,1}(i,k),70,'filled','MarkerFaceColor',[i/12 i/12 i/12])
    hold on
    scatter(trending(i,3),directionality_member{3,1}(i,k),70,'filled','MarkerFaceColor',[i/12 i/12 i/12])
end 
text(-0.13,1.08,'(a)','Units', 'Normalized', 'VerticalAlignment', 'Top')
xlim([-.6 .35])
xlabel('{\it <T>}','fontsize',12)
ylabel('{\it <D>}','fontsize',12)
lgd1 = legend('Individual','Foreigner','Institution');
lgd1.FontSize = 12;
legend boxoff

subplot(1,3,2)
for ii = 1:41
    if ca(ii,3) < -.00 && mean_vol_member(domestic(ii)) > mean(mean_vol_member(domestic))/10
        hold on
        for i = 1:1
            scatter(trending(i,domestic(ii)),directionality_member{domestic(ii),1}(i,k),70,'filled','MarkerFaceColor',color_3(3,:))
            hold on
        end
        text(trending(1,domestic(ii))+.01,directionality_member{domestic(ii),1}(1,k),inv_mem_name(domestic(ii)),'fontsize',7)
    elseif ca(ii,3) > .00 && mean_vol_member(domestic(ii)) > mean(mean_vol_member(domestic))/10
        hold on
        for i = 1:1
            scatter(trending(i,domestic(ii)),directionality_member{domestic(ii),1}(i,k),70,'filled','MarkerFaceColor',color_3(1,:))
            hold on
        end
        text(trending(1,domestic(ii))+.01,directionality_member{domestic(ii),1}(1,k),inv_mem_name(domestic(ii)),'fontsize',7)
    end
end

for ii = 1:21
    if cb(ii,2) > .1 && mean_vol_member(foreign(ii)) > mean(mean_vol_member(foreign))/10
%         plot(trending(:,foreign(ii)),directionality_member{foreign(ii),1}(:,k),'--k','linewidth',1.5)
        hold on
        for i = 1:1
            scatter(trending(i,foreign(ii)),directionality_member{foreign(ii),1}(i,k),70,'filled','MarkerFaceColor',color_3(2,:))
            hold on
        end
        text(trending(1,foreign(ii))+.01,directionality_member{foreign(ii),1}(1,k),inv_mem_name(foreign(ii)),'fontsize',7)
    end
end
plot([0 0],[0 1],'k--')
hold on
plot([-1 1],[.5 .5],'k--')
hold on
text(-0.13,1.08,'(b)','Units', 'Normalized', 'VerticalAlignment', 'Top')
xlim([-.45 .28])
ylim([.43 1])
xlabel('{\it <T>}','fontsize',12)
ylabel('{\it <D>}','fontsize',12)
box on



subplot(1,3,3)
for ii = 1:41
    if ca(ii,3) < -.02 && mean_vol_member(domestic(ii)) > mean(mean_vol_member(domestic))/10
        hold on
        for i = 10:10
            scatter(trending(i,domestic(ii)),directionality_member{domestic(ii),1}(i,k),70,'filled','MarkerFaceColor',color_3(3,:))
            hold on
        end
        text(trending(10,domestic(ii))+.01,directionality_member{domestic(ii),1}(10,k),inv_mem_name(domestic(ii)),'fontsize',7)
    elseif ca(ii,3) > .02 && mean_vol_member(domestic(ii)) > mean(mean_vol_member(domestic))/10
        hold on
        for i = 10:10
            scatter(trending(i,domestic(ii)),directionality_member{domestic(ii),1}(i,k),70,'filled','MarkerFaceColor',color_3(1,:))
            hold on
        end
        text(trending(10,domestic(ii))+.01,directionality_member{domestic(ii),1}(10,k),inv_mem_name(domestic(ii)),'fontsize',7)
    end
end

for ii = 1:21
    if cb(ii,2) > .1 && mean_vol_member(foreign(ii)) > mean(mean_vol_member(foreign))/10
%         plot(trending(:,foreign(ii)),directionality_member{foreign(ii),1}(:,k),'--k','linewidth',1.5)
        hold on
        for i = 10:10
            scatter(trending(i,foreign(ii)),directionality_member{foreign(ii),1}(i,k),70,'filled','MarkerFaceColor',color_3(2,:))
            hold on
        end
        text(trending(10,foreign(ii))+.01,directionality_member{foreign(ii),1}(10,k),inv_mem_name(foreign(ii)),'fontsize',7)
    end
end
plot([0 0],[0 1],'k--')
hold on
plot([-1 1],[.5 .5],'k--')
hold on
text(-0.13,1.08,'(c)','Units', 'Normalized', 'VerticalAlignment', 'Top')
xlim([-.45 .28])
ylim([.43 1])
xlabel('{\it <T>}','fontsize',12)
ylabel('{\it <D>}','fontsize',12)
box on



dx=0.06;
dy=0.06;
x=(1-4*dx)/2.9;
y=(1-4*dy)/1;
dxx=0.01;
dyy=0.01;
AxesHandle=findobj(figure1,'Type','axes');
set(AxesHandle(3),'Position',[2*dx-dx,3*dy,x,y]);
set(AxesHandle(2),'Position',[x+2.6*dx-dx/2,3*dy,x,y]);
set(AxesHandle(1),'Position',[2*x+3.5*dx-dx/3,3*dy,x,y]);
set(gcf,'Position',[100, 100, 1800, 600])



%% Feedback & Price impact

a = zeros(50,75);
percent1 = 0;

for i = 1:121 %1210-121+1:1210
    
    temp1 = find(price_return{i,1}(1:2722) > log((100+percent1)/100));

    if ~isempty(temp1)
        for ii = 1:50
            delay = ii-25;
            temp_d = temp1(find(temp1+delay>0 & temp1+delay<=2722));
            if ~isempty(temp_d)
                a(ii,:) = a(ii,:) + mean(inv_mem_net{i,1}(temp_d+delay,:).*price{i,2}(temp_d+delay,:));
%                 a(ii,:) = a(ii,:) + mean(inv_mem_net{i,1}(temp_d+delay,:));
            end
        end
    end
end


% plot for individuals, foreigners, institutions
figure
set(gcf,'color','w')
for i = 1:3
    subplot(1,3,i)
    plot(-24:1:25,a(:,i),'k')
    hold on
    scatter(0,a(25,i),'k','filled')
    hold on
    scatter(-1,a(24,i),'k')
    hold on
    scatter(-2,a(23,i),'k')
    hold on
    scatter(1,a(26,i),'k')
    hold on
    scatter(2,a(27,i),'k')
    title(inv_mem_name{i,1},'fontsize',18)
    xlim([-25 25])
    xlabel('Time lag (day)','fontsize',16)
    ylabel('Net price (Won)','fontsize',16)
end


% plot mean & std
figure0 = figure;
set(gcf,'color','w')
subplot(1,3,1)
ind = [];
for i = 1:41
    if ca(i,1) < -.0
        ind = [ind; i];
    end
end
plot(-24:1:25,mean(a(:,domestic(ind)),2),'k','linewidth',1)
hold on
scatter(-1:1,mean(a(24:26,domestic(ind)),2),'k')
hold on
ciplot(mean(a(:,domestic(ind)),2)-std(a(:,domestic(ind))')',mean(a(:,domestic(ind)),2)+std(a(:,domestic(ind))')',-24:1:25,[color_3(1,1) color_3(2,1) color_3(3,1)])
alpha(.5)
hold on
text(-0.13,1.08,'(a)','Units', 'Normalized', 'VerticalAlignment', 'Top')
xlabel('lag (day)','fontsize',14)
ylabel('E(B-S|r>0)','fontsize',14)
xlim([-25 25])
legend('Domestic members (institution dominent)')
legend boxoff


subplot(1,3,2)
ind = [];
for i = 1:41
    if ca(i,1) > .05
        ind = [ind; i];
    end
end
plot(-24:1:25,mean(a(:,domestic(ind)),2),'k')
hold on
scatter(-1:1,mean(a(24:26,domestic(ind)),2),'k')
hold on
ciplot(mean(a(:,domestic(ind)),2)-std(a(:,domestic(ind))')',mean(a(:,domestic(ind)),2)+std(a(:,domestic(ind))')',-24:1:25,[color_3(1,3) color_3(2,3) color_3(3,3)])
alpha(.5)
hold on
text(-0.13,1.08,'(b)','Units', 'Normalized', 'VerticalAlignment', 'Top')
xlabel('lag (day)','fontsize',14)
ylabel('E(B-S|r>0)','fontsize',14)
xlim([-25 25])
legend('Domestic members (individual dominent)')
legend boxoff

subplot(1,3,3)
plot(-24:1:25,mean(a(:,foreign),2),'k')
hold on
scatter(-1:1,mean(a(24:26,foreign),2),'k')
hold on
ciplot(mean(a(:,foreign),2)-std(a(:,foreign)')',mean(a(:,foreign),2)+std(a(:,foreign)')',-24:1:25,[color_3(1,2) color_3(2,2) color_3(3,2)])
alpha(.5)
hold on
text(-0.13,1.08,'(c)','Units', 'Normalized', 'VerticalAlignment', 'Top')
xlabel('lag (day)','fontsize',14)
ylabel('E(B-S|r>0)','fontsize',14)
xlim([-25 25])
legend('Foreign members')
legend boxoff

dx=0.06;
dy=0.06;
x=(1-4*dx)/2.9;
y=(1-4*dy)/1;
dxx=0.01;
dyy=0.01;
AxesHandle=findobj(figure0,'Type','axes');
set(AxesHandle(3),'Position',[2*dx-dx,3*dy,x,y]);
set(AxesHandle(2),'Position',[x+2.6*dx-dx/2,3*dy,x,y]);
set(AxesHandle(1),'Position',[2*x+3.5*dx-dx/3,3*dy,x,y]);
set(gcf,'Position',[100, 100, 1200, 420])



% plot the mean of all members
figure
set(gcf,'color','w')
subplot(1,3,1)
for i = 1:41
    if ca(i,3) < -.05 && mean_vol_member(domestic(i)) > mean(mean_vol_member(domestic))/5
        plot(-24:1:25,a(:,domestic(i)),'k')
        hold on
    end
end
subplot(1,3,2)
for i = 1:41
    if ca(i,3) > .1 && mean_vol_member(domestic(i)) > mean(mean_vol_member(domestic))/5
        plot(-24:1:25,a(:,domestic(i)),'k')
        hold on
    end
end
subplot(1,3,3)
for i = 1:21
    plot(-24:1:25,a(:,foreign(i)),'k')
    hold on
end


%% Feedback & Price impact
% simultaneously 1 dcile & 10 decile

a = zeros(50,75);
percent1 = 0;

for i = 1:121
    
    temp1 = find(price_return{i,1}(1:2722) > log((100+percent1)/100));

    if ~isempty(temp1)
        for ii = 1:50
            delay = ii-25;
            temp_d = temp1(find(temp1+delay>0 & temp1+delay<=2722));
            if ~isempty(temp_d)
                a(ii,:) = a(ii,:) + mean(inv_mem_net{i,1}(temp_d+delay,:).*price{i,2}(temp_d+delay,:));
%                 a(ii,:) = a(ii,:) + mean(inv_mem_net{i,1}(temp_d+delay,:));
            end
        end
    end
end





figure0 = figure;
set(gcf,'color','w')


subplot(2,3,1)
ind = [];
for i = 1:41
    if ca(i,1) < -.0
        ind = [ind; i];
    end
end
plot(-24:1:25,mean(a(:,domestic(ind)),2),'k','linewidth',1)
hold on
scatter(-1:1,mean(a(24:26,domestic(ind)),2),'k')
hold on
ciplot(mean(a(:,domestic(ind)),2)-std(a(:,domestic(ind))')',mean(a(:,domestic(ind)),2)+std(a(:,domestic(ind))')',-24:1:25,[color_3(1,1) color_3(2,1) color_3(3,1)])
alpha(.5)
hold on
text(-0.13,1.08,'(a)','Units', 'Normalized', 'VerticalAlignment', 'Top')
xlabel('lag (day)','fontsize',14)
ylabel('E(B-S|r>0)','fontsize',14)
xlim([-25 25])
legend('Domestic members (institution dominent)')
legend boxoff


subplot(2,3,2)
ind = [];
for i = 1:41
    if ca(i,1) > .05
        ind = [ind; i];
    end
end
plot(-24:1:25,mean(a(:,domestic(ind)),2),'k')
hold on
scatter(-1:1,mean(a(24:26,domestic(ind)),2),'k')
hold on
ciplot(mean(a(:,domestic(ind)),2)-std(a(:,domestic(ind))')',mean(a(:,domestic(ind)),2)+std(a(:,domestic(ind))')',-24:1:25,[color_3(1,3) color_3(2,3) color_3(3,3)])
alpha(.5)
hold on
text(-0.13,1.08,'(b)','Units', 'Normalized', 'VerticalAlignment', 'Top')
xlabel('lag (day)','fontsize',14)
ylabel('E(B-S|r>0)','fontsize',14)
xlim([-25 25])
legend('Domestic members (individual dominent)')
legend boxoff

subplot(2,3,3)
plot(-24:1:25,mean(a(:,foreign),2),'k')
hold on
scatter(-1:1,mean(a(24:26,foreign),2),'k')
hold on
ciplot(mean(a(:,foreign),2)-std(a(:,foreign)')',mean(a(:,foreign),2)+std(a(:,foreign)')',-24:1:25,[color_3(1,2) color_3(2,2) color_3(3,2)])
alpha(.5)
hold on
text(-0.13,1.08,'(c)','Units', 'Normalized', 'VerticalAlignment', 'Top')
xlabel('lag (day)','fontsize',14)
ylabel('E(B-S|r>0)','fontsize',14)
xlim([-25 25])
legend('Foreign members')
legend boxoff



a = zeros(50,75);
percent1 = 0;

for i = 1210-121+1:1210
    
    temp1 = find(price_return{i,1}(1:2722) > log((100+percent1)/100));

    if ~isempty(temp1)
        for ii = 1:50
            delay = ii-25;
            temp_d = temp1(find(temp1+delay>0 & temp1+delay<=2722));
            if ~isempty(temp_d)
                a(ii,:) = a(ii,:) + mean(inv_mem_net{i,1}(temp_d+delay,:).*price{i,2}(temp_d+delay,:));
%                 a(ii,:) = a(ii,:) + mean(inv_mem_net{i,1}(temp_d+delay,:));
            end
        end
    end
end


subplot(2,3,4)
ind = [];
for i = 1:41
    if ca(i,1) < -.0
        ind = [ind; i];
    end
end
plot(-24:1:25,mean(a(:,domestic(ind)),2),'k','linewidth',1)
hold on
scatter(-1:1,mean(a(24:26,domestic(ind)),2),'k')
hold on
ciplot(mean(a(:,domestic(ind)),2)-std(a(:,domestic(ind))')',mean(a(:,domestic(ind)),2)+std(a(:,domestic(ind))')',-24:1:25,[color_3(1,1) color_3(2,1) color_3(3,1)])
alpha(.5)
hold on
text(-0.13,1.08,'(d)','Units', 'Normalized', 'VerticalAlignment', 'Top')
xlabel('lag (day)','fontsize',14)
ylabel('E(B-S|r>0)','fontsize',14)
xlim([-25 25])
legend('Domestic members (institution dominent)')
legend boxoff


subplot(2,3,5)
ind = [];
for i = 1:41
    if ca(i,1) > .05
        ind = [ind; i];
    end
end
plot(-24:1:25,mean(a(:,domestic(ind)),2),'k')
hold on
scatter(-1:1,mean(a(24:26,domestic(ind)),2),'k')
hold on
ciplot(mean(a(:,domestic(ind)),2)-std(a(:,domestic(ind))')',mean(a(:,domestic(ind)),2)+std(a(:,domestic(ind))')',-24:1:25,[color_3(1,3) color_3(2,3) color_3(3,3)])
alpha(.5)
hold on
text(-0.13,1.08,'(e)','Units', 'Normalized', 'VerticalAlignment', 'Top')
xlabel('lag (day)','fontsize',14)
ylabel('E(B-S|r>0)','fontsize',14)
xlim([-25 25])
legend('Domestic members (individual dominent)')
legend boxoff

subplot(2,3,6)
plot(-24:1:25,mean(a(:,foreign),2),'k')
hold on
scatter(-1:1,mean(a(24:26,foreign),2),'k')
hold on
ciplot(mean(a(:,foreign),2)-std(a(:,foreign)')',mean(a(:,foreign),2)+std(a(:,foreign)')',-24:1:25,[color_3(1,2) color_3(2,2) color_3(3,2)])
alpha(.5)
hold on
text(-0.13,1.08,'(f)','Units', 'Normalized', 'VerticalAlignment', 'Top')
xlabel('lag (day)','fontsize',14)
ylabel('E(B-S|r>0)','fontsize',14)
xlim([-25 25])
legend('Foreign members')
legend boxoff


dx=0.05;
dy=0.05;
x=(1-5*dx)/2.9;
y=(1-5*dy)/1.9;
dxx=0.01;
dyy=0.01;
AxesHandle=findobj(figure0,'Type','axes');
set(AxesHandle(6),'Position',[dx,y+3.45*dy,x,y]);
set(AxesHandle(5),'Position',[x+2.5*dx,y+3.45*dy,x,y]);
set(AxesHandle(4),'Position',[2*x+4*dx,y+3.45*dy,x,y]);
set(AxesHandle(3),'Position',[dx,dy+dy/2,x,y]);
set(AxesHandle(2),'Position',[x+2.5*dx,dy+dy/2,x,y]);
set(AxesHandle(1),'Position',[2*x+4*dx,dy+dy/2,x,y]);
set(gcf,'Position',[100, 100, 1300, 750])



%% RMT for correlation matrix

num_firm = 1210;

corr_fr = zeros(num_firm,length(year_day));
corr_fr_1 = zeros(num_firm,length(year_day));
corr_fr_2 = zeros(num_firm,length(year_day));
corr_fr_3 = zeros(num_firm,length(year_day));
herd_adf = cell(num_firm,3);
leval = zeros(num_firm,length(year_day));
vote = zeros(75,1);
vote_count = zeros(75,1);
RMT_valid = zeros(num_firm,length(year_day));
lambda_max_save = zeros(num_firm,length(year_day)-1);

ki_not = 0; % Kiwoom is not present
for f = 1:num_firm
    f
    for y = 1:length(year_day)-1
        onoff_trading = find(vol{f,1}(year_day{y,1}) > 0);
        
        active_day = sum(sign(inv_mem_net{f,2}(year_day{y,1},:)),1);
        active_inv = find(active_day > length(onoff_trading)/1.5);
        clear active day onoff_trading
        
        domestic_act = intersect(domestic, active_inv);
        foreign_act = intersect(foreign, active_inv);
        member_act = [domestic_act, foreign_act];
        clear active_inv
        
        if length(domestic_act) + length(foreign_act) > 0
            corr_df = corr(inv_mem_net{f,1}(year_day{y,1},member_act));

            [eig_vec, eig_val] = eig(corr_df);
            [eig_val_sort, ind_val_sort] = sort(on_diagonal(eig_val),'descend');
            largest_eig_val = eig_val_sort(1);
            largest_eig_vec = eig_vec(:,ind_val_sort(1));
            
            Q = size(year_day{y,1},1) / size(member_act,2);
            lambda_min = 1 + 1/Q - 2*sqrt(1/Q);
            lambda_max = 1 + 1/Q + 2*sqrt(1/Q);
            lambda_max_save(f,y) = lambda_max;
            
%             figure1 = figure;
%             set(gcf,'color','w')
%             subplot(1,2,1)
%             x=linspace(lambda_min,lambda_max,10000); 
%             histo = histogram(eig_val_sort,30);
%             histo.Visible = 'off';
%             bar((histo.BinEdges(1:end-1) + histo.BinEdges(2:end))/2,histo.BinCounts/8)
%             hold on
%             plot(x,Q/(2*pi)*sqrt((lambda_max-x).*(x-lambda_min))./x,'k--','Linewidth',1.5)
%             xlabel('\lambda','fontsize',16)
%             ylabel('\rho(\lambda)','fontsize',16)
%             legend('dataset','random matrix')
%             legend boxoff
%             text(2.8, 0.3,'largest eigenvalue(\lambda_1)')
%             ylim([0 .92])
%             text(-0.14,1.08,'(a)','Units', 'Normalized', 'VerticalAlignment', 'Top')
            
            
            if sum(member_act==46) > 0
                if largest_eig_vec(find(member_act==46)) <= 0
                    factor = inv_mem_net{f,1}(year_day{y,1},member_act) * largest_eig_vec;
                    vote(member_act) = vote(member_act) + largest_eig_vec;
                elseif largest_eig_vec(find(member_act==46)) > 0
                    factor = inv_mem_net{f,1}(year_day{y,1},member_act) * (-largest_eig_vec);
                    vote(member_act) = vote(member_act) + (-largest_eig_vec);
                end

                vote_count(member_act) = vote_count(member_act) + 1;

                q = length(year_day{y,1}) / length(member_act);
                if largest_eig_val > 1+ 1/q + 2*sqrt(1/q)
                    RMT_valid(f,y) = 1;
                end

                leval(f,y) = largest_eig_val;

                % 'corr_fr' is the correlation coefficient between factor and price_return
                corr_fr(f,y) = corr(factor,price_return{f,1}(year_day{y,1}));

                % 'corr_fr_1' is the correlation coefficient between individual and price_return
                % 'corr_fr_2' is the correlation coefficient between foreigner and price_return
                % 'corr_fr_3' is the correlation coefficient between institution and price_return
                corr_fr_1(f,y) = corr(inv_mem_net{f,1}(year_day{y,1},1),price_return{f,1}(year_day{y,1}));
                corr_fr_2(f,y) = corr(inv_mem_net{f,1}(year_day{y,1},2),price_return{f,1}(year_day{y,1}));
                corr_fr_3(f,y) = corr(inv_mem_net{f,1}(year_day{y,1},3),price_return{f,1}(year_day{y,1}));
            else
                
                ki_not = ki_not + 1;
                
                ki_not_list(ki_not,1) = f;
                ki_not_list(ki_not,2) = y;
            end
            
        end
    end
end



% check which direction for kiwoom to move compared with price return
clear aa
for f = 1:121
    for y = 1:11
        aa(f,y) = corr(inv_mem_net{f,1}(year_day{y,1},46),price_return{f,1}(year_day{y,1}));
    end
end
sprintf('direction of Kiwoom: negative %.2f %%',length(find(delniz(aa(:),[1,2,3],0)<0))/length(aa(:))*100)


% timeseries (yearly): factor(individual, foreigner, institution) vs. return
fs = 1;
fe = 121;
clear a b
for i = 1:11
    a(i,1) = mean(delniz(corr_fr(fs:fe,i),[1,2,3],0));
    a(i,2) = mean(delniz(corr_fr_1(fs:fe,i),[1,2,3],0));
    a(i,3) = mean(delniz(corr_fr_2(fs:fe,i),[1,2,3],0));
    a(i,4) = mean(delniz(corr_fr_3(fs:fe,i),[1,2,3],0));
end

figure
plot(a(:,1),'k--')
hold on
plot(a(:,2),'r')
hold on
plot(a(:,3),'g')
hold on
plot(a(:,4),'b')
legend('factor','individual','foreigner','institution')


figure
set(gcf,'color','w')
for f = 1:10
    l1(f,1) = mean(mean(lambda_max_save(1+(f-1)*121:f*121,:)));
    l2(f,1) = mean(mean(leval(1+(f-1)*121:f*121,:)));
end
plot(l1)
hold on
plot(l2,'--')

%% correlation vs. largest eigenvalue

figure
set(gcf,'color','w')

f = 10;
leval_temp = leval((f-1)*121+1 : f*121,:);
corr_fr_temp = corr_fr((f-1)*121+1 : f*121,:);

leval_121 = leval_temp(:);
corr_fr_121 = corr_fr_temp(:);

hold on
scatter(corr_fr_121,leval_121,'b.')
xlabel('corr. factor & return')
ylabel('largest eigen value')

f = 1;
leval_temp = leval((f-1)*121+1 : f*121,:);
corr_fr_temp = corr_fr((f-1)*121+1 : f*121,:);

leval_121 = leval_temp(:);
corr_fr_121 = corr_fr_temp(:);

hold on
scatter(corr_fr_121,leval_121,'r')





figure1 = figure;
set(gcf,'color','w')
subplot(1,2,1)
for ff = 1:10
%     f = 11 - ff;
	f = ff;
    for i = (f-1)*121+1 : f*121
        if f == 1 || f == 10
            scatter(mean(corr_fr(i,1:11)),mean(leval(i,1:11)),'filled','MarkerFaceColor',[f/14 f/14 f/14],'MarkerEdgeColor',[f/14 f/14 f/14])
            hold on
        end
    end
end
text(-0.14,1.08,'(a)','Units', 'Normalized', 'VerticalAlignment', 'Top')
xlabel('Corr. factor and return','Fontsize',12)
ylabel('Largest eigenvalue','Fontsize',12)

subplot(1,2,2)
for f = 1:10
    scatter(mean(mean(corr_fr((f-1)*121+1 : f*121,1:11))),mean(mean(leval((f-1)*121+1 : f*121,1:11))),'MarkerFaceColor',[f/14 f/14 f/14],'MarkerEdgeColor',[f/14 f/14 f/14])
    text(mean(mean(corr_fr((f-1)*121+1 : f*121,1:11)))+.01,mean(mean(leval((f-1)*121+1 : f*121,1:11))),sprintf('%.0f',f))
    hold on
    
    cc(f,1) = mean(mean(corr_fr((f-1)*121+1 : f*121,1:11)));
    cc(f,2) = mean(mean(leval((f-1)*121+1 : f*121,1:11)));
end

hold on
[b, ~, ~, ~,stats] = regress(cc(:,2),[ones(size(cc(:,1))), cc(:,1)]);
pt = plot([-.05 .4],[b(1)+b(2)*(-.05) b(1)+b(2)*(.4)],'k--');
legend(pt(1),sprintf('y=%.2f+%.2fx, R^2=%.2f',b(2),b(1),stats(1)))
legend boxoff
text(-0.14,1.08,'(b)','Units', 'Normalized', 'VerticalAlignment', 'Top')
xlabel('Corr. factor and return','Fontsize',12)
ylabel('Largest eigenvalue(\lambda_1)','Fontsize',12)
xlim([-.04 .4])
ylim([2.4, 3.7])

dx=0.08;
dy=0.06;
x=(1-4*dx)/1.9;
y=(1-4*dy)/1.05;
dxx=0.01;
dyy=0.01;
AxesHandle=findobj(figure1,'Type','axes');
set(AxesHandle(2),'Position',[dx+0.01,3*dy,x,y]);
set(AxesHandle(1),'Position',[x+2.6*dx,3*dy,x,y]);
set(gcf,'Position',[100, 100, 1000, 450])



%% Foreign herding vs. kospi

win_size = 20;

foreign_herd = zeros(2801-win_size+1,1);

for i = 1:2801-win_size+1
    
    for ii = 1:75
        active_day(ii,1) = length(find(abs(firm_net_member{ii,1}(i:i+win_size-1))>0));
    end
    active_firm = find(active_day > length(win_size)/1.5);
    clear active day onoff_trading

    domestic_act = intersect(domestic, active_firm);
    foreign_act = intersect(foreign, active_firm);
    member_act = [domestic_act; foreign_act];
    clear active_firm
    
    temp = zeros(win_size,length(domestic_act));
    for ii = 1:length(domestic_act)
        temp(:,ii) = firm_net_member{domestic_act(ii,1),1}(i:i+win_size-1);
    end
    
    foreign_herd(i,1) = mean(off_diagonal(corr(temp)));
    
end

kospi_win = kospi(win_size:end);
kospi_return_win = kospi_return(win_size:end);
    
fitlm(foreign_herd,kospi_return_win)
    
    
%% stock network

threshold = .05; % for ratio

% figure1 = figure;
% set(gcf,'color','w')

% figure2 = figure;
% set(gcf,'color','w')

for year_ind = 1:11
    start_day = year_day{year_ind,1}(1);
    end_day = year_day{year_ind,1}(end);
    len_day = end_day - start_day + 1;
    start_firm = 1;
    end_firm = 121;%1216
    firm_anl = end_firm - start_firm + 1; % firms analyzed

    price_ori = zeros(len_day,num_firm);
    day_i = 1;
    for i = start_day:end_day
        firm_i = 1;
        for ii = start_firm:end_firm
            price_ori(day_i,firm_i) = price{firm_i,1}(day_i,1);
            firm_i = firm_i + 1;
        end
        day_i = day_i + 1;
    end

    % mean volume of member firm
    % this matrix will be used in the inventory correlation with ifd
    member_cons = [];
    for i = 1:size(firm_net,1)
        sum_member = sum(firm_net{i,2}(start_day:end_day,start_firm:end_firm).*price_ori,2);
        if length(find(sum_member ~= 0)) >= len_day/1.5
            member_cons = [member_cons; i];
        end
        ind = find(sum_member~=0); % 0인것 제외하고 평균
        mean_vol_member(i,1) = mean(sum_member(ind));
    end

    domestic_cons = intersect(member_cons,domestic)';
    foreign_cons = intersect(member_cons,foreign)';

    clear ind sum_member price_ori member_cons

    i_name = inv_mem_name([domestic_cons,foreign_cons]);


    stock_corr = zeros(firm_anl,firm_anl);
    for i = 1:size(price_return(start_firm:end_firm,1),1)
        for ii = 1:size(price_return(start_firm:end_firm,1),1)
            stock_corr(i,ii) = corr(price_return{i,1}(start_day:end_day),price_return{ii,1}(start_day:end_day));
        end
    end

    ratio = zeros(size(inv_mem_net(start_firm:end_firm,:),1),size(inv_mem_name,1));

    for i = 1:size(inv_mem_net(start_firm:end_firm,:),1)
        for ii = 1:size(inv_mem_name,1)
            ratio(i,ii) = sum(inv_mem_net{i,1}(start_day:end_day,ii).*price{i,2}(start_day:end_day)) / sum(inv_mem_net{i,2}(start_day:end_day,ii).*price{i,2}(start_day:end_day));
            if isnan(ratio(i,ii))
                ratio(i,ii) = 0;
            end

            if ratio(i,ii) > threshold
                state(i,ii) = 1;
            elseif ratio(i,ii) < -threshold
                state(i,ii) = -1;
            else
                state(i,ii) = 0;
            end
        end
    end

    %     inv_mem_corr = similarity(ratio(:,[4,5,6,7,8,9,10,11,12,foreign]),'cosine');
    inv_mem_corr = similarity(ratio(:,[domestic_cons,foreign_cons]),'cosine');
    inv_ind = find(~isnan(inv_mem_corr(:,1)));
    inv_mem_corr = inv_mem_corr(~isnan(inv_mem_corr(:,1)),~isnan(inv_mem_corr(:,1)));

    asc(year_ind,1) = mean(stock_corr(:));
    asc(year_ind,2) = std(stock_corr(:));
    aic(year_ind,1) = mean(inv_mem_corr(:));
    aic(year_ind,2) = std(inv_mem_corr(:));

    
    figure(1)
    subplot(2,6,year_ind)
    hold on
    imagesc(stock_corr);
    colorbar
    lim = caxis;
    caxis([-.5 1])
    set(gca,'ytick',start_firm:end_firm,'yticklabel',code_firm(start_firm:end_firm),'fontsize',9)
    set(gca,'xtick',start_firm:end_firm,'xticklabel',code_firm(start_firm:end_firm),'fontsize',9)
    xlabel('larger --- Stock (sorted by market capital) --- small','FontSize',30)
    ylabel('larger --- Stock (sorted by market capital) --- small','FontSize',30)
    xtickangle(90)
    title('Stock correlation','FontSize',35)

    figure(2)
    hold on
    subplot(2,6,year_ind)
    hold on
    imagesc(inv_mem_corr);
    hold on
    plot([length(domestic_cons)+.5 length(domestic_cons)+.5],[.5 length(domestic_cons)+.5+size(foreign_cons,2)],'k')
    hold on
    plot([.5 length(domestic_cons)+.5+size(foreign_cons,2)],[length(domestic_cons)+.5 length(domestic_cons)+.5],'k')
    set(gca,'ytick',1:length(i_name),'yticklabel',i_name(inv_ind),'fontsize',14)
    set(gca,'xtick',1:length(i_name),'xticklabel',i_name(inv_ind),'fontsize',14)
    xtickangle(90)
    colorbar
    lim = caxis;
    caxis([-1 1])
    title('Similarity of investor','FontSize',35)

    % save node
    fid = fopen(sprintf('./save_network/inv_node_%d.csv',year_ind),'w');
    fprintf(fid,'id,label,color,size\n');
    for i = 1:size(domestic_cons,2)
        fprintf(fid,'%d,%s,%d,%d\n',domestic_cons(i),inv_mem_name{domestic_cons(i),1},1,2*sqrt(mean_vol_member(domestic_cons(i)))/10^5);
    end
    for i = 1:size(foreign_cons,2)
        fprintf(fid,'%d,%s,%d,%d\n',foreign_cons(i),inv_mem_name{foreign_cons(i),1},2,2*sqrt(mean_vol_member(foreign_cons(i)))/10^5);
    end
    fclose(fid);
    
    % save edge
    fid = fopen(sprintf('./save_network/inv_edge_%d.csv',year_ind),'w');
    fprintf(fid,'Source,Target,Weight\n');
    for i=2:length(domestic_cons)+length(foreign_cons)
        for ii=1:i-1
            if ~isnan(inv_mem_corr(i,ii))
                if i <= length(domestic_cons)
                    i_ind = domestic_cons(i);
                else
                    i_ind = foreign_cons(i-length(domestic_cons));
                end
                if ii <= length(domestic_cons)
                    ii_ind = domestic_cons(ii);
                else
                    ii_ind = foreign_cons(ii-length(domestic_cons));
                end
                fprintf(fid,'%d,%d,%d\n',i_ind,ii_ind,inv_mem_corr(i,ii));
            end
        end
    end
    fclose(fid);

end


%% Portfolio Asymmetry

ind_target = ind_kospi;
 
% domestic_ind = find(ca(:,3)<0 & ca(:,1)>0);
% domestic_ins = find(ca(:,3)>=0 & ca(:,1)<0);
domestic_ind = find(ca(:,3)<0);
domestic_ins = find(ca(:,3)>=0);

asym_domestic = zeros(2801,41);
asym_foreign = zeros(2801,21);
asym_domestic_ind = zeros(2801,length(domestic_ind));
asym_domestic_ins = zeros(2801,length(domestic_ins));

for i = 1:41
    for ii = 1:2801
        asym_domestic(ii,i) = length(find(sign(firm_net{domestic(i),1}(ii,ind_target))>0)) ./ length(find(sign(firm_net{domestic(i),1}(ii,ind_target))~=0));
    end
end

for i = 1:21
    for ii = 1:2801
        asym_foreign(ii,i) = length(find(sign(firm_net{foreign(i),1}(ii,ind_target))>0)) ./ length(find(sign(firm_net{foreign(i),1}(ii,ind_target))~=0));
    end
end

for i = 1:length(domestic_ind)
    for ii = 1:2801
        asym_domestic_ind(ii,i) = length(find(sign(firm_net{domestic(domestic_ind(i)),1}(ii,ind_target))>0)) ./ length(find(sign(firm_net{domestic(domestic_ind(i)),1}(ii,ind_target))~=0));
    end
end

for i = 1:length(domestic_ins)
    for ii = 1:2801
        asym_domestic_ins(ii,i) = length(find(sign(firm_net{domestic(domestic_ins(i)),1}(ii,ind_target))>0)) ./ length(find(sign(firm_net{domestic(domestic_ins(i)),1}(ii,ind_target))~=0));
    end
end
        


%%

th=0.00;

for i = 1:50

    asym_d1(i,1) = length(find(delniz(asym_domestic_ind(:),[1,2],0) < 0.5 -th | delniz(asym_domestic_ind(:),[1,2],0) > 0.5 +th))/length(delniz(asym_domestic_ind(:),[1,2],0));
    asym_d2(i,1) = length(find(delniz(asym_domestic_ins(:),[1,2],0) < 0.5 -th | delniz(asym_domestic_ins(:),[1,2],0) > 0.5 +th))/length(delniz(asym_domestic_ins(:),[1,2],0));
    asym_f(i,1) = length(find(delniz(asym_foreign(:),[1,2],0) < 0.5 -th | delniz(asym_foreign(:),[1,2],0) > 0.5 +th))/length(delniz(asym_foreign(:),[1,2],0));
    asym_d(i,1) = length(find(delniz(asym_domestic(:),[1,2],0) < 0.5 -th | delniz(asym_domestic(:),[1,2],0) > 0.5 +th))/length(delniz(asym_domestic(:),[1,2],0));

    th = th + 0.01;

end

figure
set(gcf,'color','w')
plot(asym_f,'k')
hold on
plot(asym_d2)
hold on
plot(asym_d1)
hold on
plot(asym_d)
% set(gca,'yscale','log')
legend('f','ins','ind','d')
xticks([10 20 30 40 50])
xticklabels({'0.1','0.2','0.3','0.4','0.5'})
legend('Foreign','Domestic(Institutions)','Domestic(Individuals)','Domestic')
xlabel('threshold')
ylabel('asymmetry')

%%

th = 0.00;
cr = zeros(50,3);

for t = 1:50

    for i = 1:2801
        d(i,1) = length(find(delniz(asym_domestic(i,:),[1,2],0) < 0.5 -th | delniz(asym_domestic(i,:),[1,2],0) > 0.5 +th)) / length(delniz(asym_domestic(i,:),[1,2],0));
        di(i,1) = length(find(delniz(asym_domestic_ind(i,:),[1,2],0) < 0.5 -th | delniz(asym_domestic_ind(i,:),[1,2],0) > 0.5 +th)) / length(delniz(asym_domestic_ind(i,:),[1,2],0));
        ds(i,1) = length(find(delniz(asym_domestic_ins(i,:),[1,2],0) < 0.5 -th | delniz(asym_domestic_ins(i,:),[1,2],0) > 0.5 +th)) / length(delniz(asym_domestic_ins(i,:),[1,2],0));
        fo(i,1) = length(find(delniz(asym_foreign(i,:),[1,2],0) < 0.5 -th | delniz(asym_foreign(i,:),[1,2],0) > 0.5 +th)) / length(delniz(asym_foreign(i,:),[1,2],0));
    end

    cr(t,1) = corr(di,abs(kospi_return));
    cr(t,2) = corr(ds,abs(kospi_return));
    cr(t,3) = corr(d,abs(kospi_return));
    cr(t,4) = corr(fo,abs(kospi_return));

    th = th + 0.01;
    
    aaa = fitlm([di ds fo],abs(kospi_return));
    bbb(t,1) = aaa.Rsquared.Adjusted;
    
end

% figure
% set(gcf,'color','w')
% plot(cr,'linewidth',1.5)
% xticks([10 20 30 40 50])
% xticklabels({'0.1','0.2','0.3','0.4','0.5'})
% legend('Domestic(Individuals)','Domestic(Institutions)','Domestic','Foreign')
% xlabel('threshold')
% ylabel('Correlation with KOSPI return')

[ta tb] = max(bbb);
th = (tb-1)/100;
for i = 1:2801
    d(i,1) = length(find(delniz(asym_domestic(i,:),[1,2],0) < 0.5 -th | delniz(asym_domestic(i,:),[1,2],0) > 0.5 +th)) / length(delniz(asym_domestic(i,:),[1,2],0));
    di(i,1) = length(find(delniz(asym_domestic_ind(i,:),[1,2],0) < 0.5 -th | delniz(asym_domestic_ind(i,:),[1,2],0) > 0.5 +th)) / length(delniz(asym_domestic_ind(i,:),[1,2],0));
    ds(i,1) = length(find(delniz(asym_domestic_ins(i,:),[1,2],0) < 0.5 -th | delniz(asym_domestic_ins(i,:),[1,2],0) > 0.5 +th)) / length(delniz(asym_domestic_ins(i,:),[1,2],0));
    fo(i,1) = length(find(delniz(asym_foreign(i,:),[1,2],0) < 0.5 -th | delniz(asym_foreign(i,:),[1,2],0) > 0.5 +th)) / length(delniz(asym_foreign(i,:),[1,2],0));
end

fitlm([di ds fo],abs(kospi_return))

%% Portfolio Asymmetry


asym_domestic = cell(10,1);
asym_foreign = cell(10,1);
asym_domestic_ind = cell(10,1);
asym_domestic_ins = cell(10,1);

for f = 1:10
    firm_start = (f-1)*121+1;
    firm_end = f*121;
    
    domestic_ind = find(ca(:,3)<0 & ca(:,1)>0);
    domestic_ins = find(ca(:,3)>=0 & ca(:,1)<0);
    
    asym_domestic{f,1} = zeros(2801,41);
    asym_foreign{f,1} = zeros(2801,21);
    asym_domestic_ind{f,1} = zeros(2801,length(domestic_ind));
    asym_domestic_ins{f,1} = zeros(2801,length(domestic_ins));

    for i = 1:41
        for ii = 1:2801
            asym_domestic{f,1}(ii,i) = length(find(sign(firm_net{domestic(i),1}(ii,firm_start:firm_end))>0)) ./ length(find(sign(firm_net{domestic(i),1}(ii,firm_start:firm_end))~=0));
        end
    end

    for i = 1:21
        for ii = 1:2801
            asym_foreign{f,1}(ii,i) = length(find(sign(firm_net{foreign(i),1}(ii,firm_start:firm_end))>0)) ./ length(find(sign(firm_net{foreign(i),1}(ii,firm_start:firm_end))~=0));
        end
    end

    for i = 1:length(domestic_ind)
        for ii = 1:2801
            asym_domestic_ind{f,1}(ii,i) = length(find(sign(firm_net{domestic(domestic_ind(i)),1}(ii,firm_start:firm_end))>0)) ./ length(find(sign(firm_net{domestic(domestic_ind(i)),1}(ii,firm_start:firm_end))~=0));
        end
    end

    for i = 1:length(domestic_ins)
        for ii = 1:2801
            asym_domestic_ins{f,1}(ii,i) = length(find(sign(firm_net{domestic(domestic_ins(i)),1}(ii,firm_start:firm_end))>0)) ./ length(find(sign(firm_net{domestic(domestic_ins(i)),1}(ii,firm_start:firm_end))~=0));
        end
    end
        
end

%%

figure
set(gcf,'color','w')

for f = 1:10

    th=0.00;

    for i = 1:50

        asym_d1(i,1) = length(find(delniz(asym_domestic_ind{f,1}(:),[1,2],0) < 0.5 -th | delniz(asym_domestic_ind{f,1}(:),[1,2],0) > 0.5 +th))/length(delniz(asym_domestic_ind{f,1}(:),[1,2],0));
        asym_d2(i,1) = length(find(delniz(asym_domestic_ins{f,1}(:),[1,2],0) < 0.5 -th | delniz(asym_domestic_ins{f,1}(:),[1,2],0) > 0.5 +th))/length(delniz(asym_domestic_ins{f,1}(:),[1,2],0));
        asym_f(i,1) = length(find(delniz(asym_foreign{f,1}(:),[1,2],0) < 0.5 -th | delniz(asym_foreign{f,1}(:),[1,2],0) > 0.5 +th))/length(delniz(asym_foreign{f,1}(:),[1,2],0));

        th = th + 0.01;

    end

    subplot(2,5,f)
    plot(asym_f)
    hold on
    plot(asym_d2)
    hold on
    plot(asym_d1)
    set(gca,'yscale','log')
    legend('f','ins','ind')
    
end


%%

th = 0.10;
f = 1;

for i = 1:2801
    di(i,1) = length(find(delniz(asym_domestic_ind{f,1}(i,:),[1,2],0) < 0.5 -th | delniz(asym_domestic_ind{f,1}(i,:),[1,2],0) > 0.5 +th)) / length(delniz(asym_domestic_ind{f,1}(i,:),[1,2],0));
    ds(i,1) = length(find(delniz(asym_domestic_ins{f,1}(i,:),[1,2],0) < 0.5 -th | delniz(asym_domestic_ins{f,1}(i,:),[1,2],0) > 0.5 +th)) / length(delniz(asym_domestic_ins{f,1}(i,:),[1,2],0));
    fo(i,1) = length(find(delniz(asym_foreign{f,1}(i,:),[1,2],0) < 0.5 -th | delniz(asym_foreign{f,1}(i,:),[1,2],0) > 0.5 +th)) / length(delniz(asym_foreign{f,1}(i,:),[1,2],0));
end


figure
set(gcf,'color','w')
subplot(1,3,1)
scatter(di,abs(kospi_return),'.')
title(sprintf('%.3f',corr(di,abs(kospi_return))))
set(gca,'yscale','log')
subplot(1,3,2)
scatter(ds,abs(kospi_return),'.')
title(sprintf('%.3f',corr(ds,abs(kospi_return))))
set(gca,'yscale','log')
subplot(1,3,3)
scatter(fo,abs(kospi_return),'.')
title(sprintf('%.3f',corr(fo,abs(kospi_return))))
set(gca,'yscale','log')


figure
set(gcf,'color','w')
scatter(di,abs(kospi_return),'.r')
hold on
scatter(ds,abs(kospi_return),'.g')
hold on
scatter(fo,abs(kospi_return),'.k')
set(gca,'yscale','log')








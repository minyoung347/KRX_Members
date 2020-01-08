
%% Add Nation data to domestic institutions
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


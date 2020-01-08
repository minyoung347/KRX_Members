
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



    
%% Stock Network (Portfolio)

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







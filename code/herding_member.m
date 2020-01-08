
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

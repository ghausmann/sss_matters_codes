function moments = quarterly_moments_5p(myt_sim_exp)

treplications=size(myt_sim_exp,3);
tn_quarters=floor(size(myt_sim_exp,2)/3);

Y_quarterly_sum_aggregation=NaN(treplications,tn_quarters);
C_quarterly_sum_aggregation=NaN(treplications,tn_quarters);
I_quarterly_sum_aggregation=NaN(treplications,tn_quarters);

for ii=1:tn_quarters
    %FGRU aggregation by sum over months
    Y_quarterly_sum_aggregation(:,ii)=squeeze(sum(myt_sim_exp(5,(ii-1)*3+1:ii*3,:)));
    %Y_level_quarterly_sum_aggregation(:,ii)=squeeze(sum(exp(simulated_series(strmatch('Y',M_.endo_names,'exact'),(ii-1)*3+1:ii*3,:))));
    C_quarterly_sum_aggregation(:,ii)=squeeze(sum(myt_sim_exp(1,(ii-1)*3+1:ii*3,:)));
    I_quarterly_sum_aggregation(:,ii)=squeeze(sum(myt_sim_exp(7,(ii-1)*3+1:ii*3,:)));
    %         NX_quarterly_sum_aggregation(:,ii)=Y_quarterly_sum_aggregation(:,ii)-C_quarterly_sum_aggregation(:,ii) - I_quarterly_sum_aggregation(:,ii);
    %         NX_Y_quarterly_sum_aggregation(:,ii)=NX_quarterly_sum_aggregation(:,ii)./Y_quarterly_sum_aggregation(:,ii); %take value third month for quarter as it was defined backward-looking
end

Y_quarterly_sum_aggregation_log = log(Y_quarterly_sum_aggregation);
C_quarterly_sum_aggregation_log = log(C_quarterly_sum_aggregation);
I_quarterly_sum_aggregation_log = log(I_quarterly_sum_aggregation);
%NX_quarterly_sum_aggregation_log = log(NX_quarterly_sum_aggregation);


%keyboard;

% compute cyclical components
%FGRU Aggregation
[~, Y_cyc_FGRU]=sample_hp_filter(Y_quarterly_sum_aggregation_log',1600);
[~, C_cyc_FGRU]=sample_hp_filter(C_quarterly_sum_aggregation_log',1600);
[~, I_cyc_FGRU]=sample_hp_filter(I_quarterly_sum_aggregation_log',1600);

% [Y_cyc_FGRU]=(Y_quarterly_sum_aggregation_log');
% [C_cyc_FGRU]=(C_quarterly_sum_aggregation_log');
% [I_cyc_FGRU]=(I_quarterly_sum_aggregation_log');



moments(1,1)=(std(Y_cyc_FGRU));
moments(2,1)=(std(C_cyc_FGRU))/moments(1,1);
moments(3,1)=(std(I_cyc_FGRU))/moments(1,1);
      

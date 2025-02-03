function [PCS,PGS,ANS,cpu_time,ci,ci_k]=DKT_1_main(n0,nb,k,alpha,vari,delta,trial,delta_choose,eta_all,good_num,delta_choose_good)
clear samplenumber

time_1 = tic;

BestSystem = zeros(trial,1);
SampleSize = zeros(trial,1);

MU = [];
VAR = [];

%---mean&variance

% [delta(10) ,0.5delta(good_num) &0(good_num)] ,-delta -2delta...平均

best_num=10;
good_layer_num=2;

layer=10-(good_layer_num+1); %10-1(best)-2(good)
remaining=k-10-good_layer_num*good_num;
num_per_layer=floor(remaining/layer);

MU = -layer*delta*ones(1,k);
MU(1:best_num)=delta;
for numb=1:good_layer_num
    MU( (best_num+1)+good_num*(numb-1) : (best_num+1)+good_num*numb )= (1/good_layer_num)*delta*(good_layer_num-numb);
end

for numb=1:layer
    MU( (best_num+good_layer_num*good_num+1)+num_per_layer*(numb-1) : (best_num+good_layer_num*good_num)+num_per_layer*numb )= delta*(-numb);
end


VAR = vari*ones(1,k);

%---------
Best = find(MU == max(MU));
GOOD = find(MU >= max(MU) - delta_choose_good);

FileName =  ['DKT_1(0,',num2str(vari^(1/2)),')  , k = ',num2str(k),' , trial = ',num2str(trial),'.mat'];


%-----循環測試------
for t = 1:trial
    disp(t)
    clear obser_sum 
    samplenumber = n0*k;%樣本數累積
    eta_num = k;
    
    k_system = k;%目前剩餘系統數
    SystemLeft = 1:k;

    %---初始抽樣
    obser_sum = sum(mvnrnd(MU,diag(VAR),n0),1);
    %---------
    while k_system > 1
        %---計算 N 這邊指的是論文中的S
        obser_bar = (1/k_system) * sum(obser_sum) ;
        obser_sum_bar = obser_sum - obser_bar ;
        obser2 = obser_sum_bar.^2;
        N = (1/vari)*sum(obser2,2);
        %---------

        %---刪除條件
        eta = eta_all(eta_num,k);
        delta_I = (delta_choose^2);
        

        condition = ((vari^(1/2))*eta/sqrt(delta_I))^2;
        %---

        %若滿足條件式 刪除系統 回到56重新計算 若不滿足則進入81 else再抽取一個樣本
        if N >= condition
            [a,idx]=min(obser_sum);

            SystemLeft(idx)=0;
            temp = find(SystemLeft~=0);
            SystemLeft = SystemLeft(temp);
            obser_sum(:,idx)=[];

            k_system = length(SystemLeft);
            eta_num = length(SystemLeft);
            
        else
            %---後續抽樣
            obser_sum = obser_sum+mvnrnd(MU(SystemLeft),diag(VAR(SystemLeft)),1);
            samplenumber = samplenumber + nb * length(SystemLeft);
            %---

        end        
    end
    
    BestSystem(t)=SystemLeft;
    SampleSize(t)=samplenumber;

end
PCS = mean(ismember(BestSystem, Best));
PGS = mean(ismember(BestSystem, GOOD));
ANS = mean(SampleSize)/k;

[ci]= Confidence_interval_95percent(SampleSize);
sample_temp =SampleSize/k;
[ci_k]= Confidence_interval_95percent(sample_temp);

cpu_time = toc(time_1);
save(FileName)
end

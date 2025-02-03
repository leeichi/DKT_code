function [PCS,PGS,ANS,cpu_time,ci,ci_k]=DKT_2_main(n0,nb,k,vari,delta,trial,delta_choose,eta_all,B_z,good_num,delta_choose_good)
clear samplenumber

time_1 = tic;


BestSystem = zeros(trial,1);
SampleSize = zeros(trial,1);

MU = [];
VAR = [];

%---mean&variance
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



VAR = vari*((1+3*((k-(1:k))./(k-1))).^2); %decrease


%---------


Best = find(MU == max(MU));
GOOD = find(MU >= max(MU) - delta_choose_good);
FileName = ['DKT_2(0,', num2str(vari^(1/2)), ')_k_', num2str(k), '_trial_', num2str(trial), '_delta_', num2str(delta), '.mat'];
%-----循環測試------
for t = 1:trial
    disp(t)
    samplenumber = n0*k;%樣本數累積
    eta_num = k;
    
    k_system = k;%目前剩餘系統數
    SystemLeft = 1:k;
    n_i = n0*ones(1,k);

    %---初始抽樣
    obser = mvnrnd(MU,diag(VAR),n0);
    obser_squ_sum = sum(obser.^2,1);%更新變異數用
    obser_sum = sum(obser,1)./n_i;
    obser_var = (obser_squ_sum./n_i) - (obser_sum.^2);
    
    %---------
    while k_system > 1
        %---計算 N 這邊指的是論文中的S
        lambda_2 = ( 1/sum(n_i) )*sum(obser_var);

        obser_bar = (1/k_system) * sum(obser_sum) ;
        obser_sum_bar = obser_sum - obser_bar ;
        obser2 = obser_sum_bar.^2;
        N = (1/lambda_2)*sum(obser2,2);
        %---------

        %---刪除條件
        eta = eta_all(eta_num,k);
        delta_I = (delta_choose^2) ;
        
        condition = ((lambda_2^(1/2))*eta/sqrt(delta_I))^2;
        %---

        %若滿足條件式 刪除系統 若不滿足則else再抽取一個樣本
        if N >= condition
            [a,idx]=min(obser_sum);

            SystemLeft(idx)=0;
            temp = find(SystemLeft~=0);
            SystemLeft = SystemLeft(temp);

            obser_sum(:,idx)=[];
            obser_var(:,idx)=[];
            obser_squ_sum(:,idx)=[];
            
            n_i(:,idx)=[];
            

            k_system = length(SystemLeft);
            eta_num = length(SystemLeft);
            
        else
            %---後續抽樣
            [a,idx]=min(n_i./obser_var);
            temp_cal = (n_i(:,idx)+B_z)/obser_var(:,idx);
            
            for i = 1:k_system
                Lambda_i = ceil((obser_var(i)*temp_cal)-1e-10); %
                add_num = Lambda_i - n_i(i);
               
                if add_num > 0

                   number_of_system = SystemLeft(i);
                   temp_ober = normrnd(MU(number_of_system),(VAR(number_of_system))^(1/2),[1,add_num]);%normrnd(mu,sigma,數量）mvnrnd(mu,covariance,數量）
                   temp_sum = sum(temp_ober);
                   temp_squ_sum = sum(temp_ober.^2);
                   
                   obser_squ_sum(i) = obser_squ_sum(i) + temp_squ_sum;
                   obser_sum(i) = ((obser_sum(i)*n_i(i)) + temp_sum)/ (n_i(i)+add_num);
                   obser_var(i) = obser_squ_sum(i)/(n_i(i)+add_num) - obser_sum(i)^2;%
                   n_i(i)= Lambda_i;
                   samplenumber = samplenumber +add_num;
                end
            end
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

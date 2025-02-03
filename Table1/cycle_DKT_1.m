clear

trail = 1000;%重複測試次數

k = [5,10,20,100,500,1000];%系統數量

alpha = 0.1;%pcs
delta = 0.5;% mean用
delta_choose = delta;%比較用

eta=open_eta(alpha);
%

n0 = 1;%初始抽樣
nb = 1;%每次抽樣


var = 100;
result = [];
count = 1;

for f = 1:length(delta_choose)
    for z = 1:length(k)
        for i = 1:length(n0)
            for j = 1:length(nb)
                [PCS,ANS,cpu_time,ci,ci_k]=DKT_1_main(n0(i),nb(j),k(z),var,delta,trail,delta_choose(f),eta);
                result(count,1) = k(z);
                result(count,2) = n0(i);
                result(count,3) = nb(j);
                result(count,4) = PCS;
                result(count,5) = ANS;
                result(count,6) = cpu_time;
                result(count,7) = ci;
                result(count,8) = ci_k;
                result(count,9) = delta_choose(f);
                count = count + 1;
            end
        end
    end
end
clear

trail = 1000;%重複測試次數
k = [5,10,20,100,500,1000];%系統數量

alpha =0.1;%pcs
delta = 0.5;% mean用
delta_choose = delta;


eta=open_eta(alpha);
%

n0 = 30;%初始抽樣
nb = 1; %沒用到這個參數
B_z=1;

var = 25;
result = [];
count = 1;

for z = 1:length(k)
    for i = 1:length(n0)

        [PCS,ANS,cpu_time,ci,ci_k]=DKT_2_main(n0(i),k(z),var,delta,trail,delta_choose,eta,B_z);
        result(count,1) = k(z);
        result(count,2) = n0(i);
        result(count,3) = PCS;
        result(count,4) = ANS;
        result(count,5) = cpu_time;
        result(count,6) = ci;
        result(count,7) = ci_k;

        count = count + 1;

    end
end

function[ci] = Confidence_interval_95percent(SampleSize)
n = length(SampleSize);
t = tinv(0.975,n-1);%雙尾檢定
sample_std = std(SampleSize);

ci = t * (sample_std / sqrt(n));
end

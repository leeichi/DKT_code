function [list_eta]=open_eta(alpha)
% alpha=0.1;
one_minus_alpha = 1-alpha;


str1="eta_";
str2=".mat";
use_alpha=string(one_minus_alpha);

read_name = strcat(str1,use_alpha,str2);
load(read_name)
% eta=list_eta;
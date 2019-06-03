%This code was used as part of my master's dissertation to model the value
%of an interest rate swap over time. This valuation was then used to test
%techniques of dealing with credit-risk (unfortunately I cannot share the code 
%for this until my dissertation has been made available in the university's
%library)

clear;
clc;

n = 20;

k_1 = 0.1;
k_2 = 0.25;
sigma_1 = 0.025;
sigma_2 = 0.02;
ro = -0.6;

k = 0.1;
theta = 0.09;
sigma = 0.04;
r0 = 0.07;

T = 5;
delta = 1/24;
t = 0:delta:T;

A = (1-exp(-k*t))/k;
D = (theta-sigma^2/(2*k^2)).*(A-t)-sigma^2*A.^2/(4*k);
initial_bond_prices = exp(-A*r0+D);

initial_spot_rates = -log(initial_bond_prices)./t;

var_phi = initial_spot_rates + sigma_1^2/(2*k_1^2)*(1-exp(-k_1*t)).^2+...
  +sigma_2^2/(2*k_2^2)*(1-exp(-k_2*t)).^2+ro*sigma_1*sigma_2/(k_1*k_2)*...
  (1-exp(-k_1*t)).*(1-exp(-k_2*t));
var_phi(1) = r0;

A1 = @(t1,t2) 1/k_1*(1-exp(-k_1*(t2-t1)));
A2 = @(t1,t2) 1/k_2*(1-exp(-k_2*(t2-t1)));

var_y = @(t1,t2) sigma_1^2/(k_1^2)*((t2-t1)-A1(t1,t2)-k_1/2*A1(t1,t2).^2)+...
  sigma_2^2/(k_2^2)*((t2-t1)-A2(t1,t2)-k_2/2*A2(t1,t2).^2)+...
  2*ro*sigma_1*sigma_2/(k_1*k_2)*((t2-t1)-A1(t1,t2)-A2(t1,t2)-...
  (exp(-(k_1+k_2)*(t2-t1))-1)/(k_1+k_2));
 
mean_y = @(t1,t2,x1,x2) (1-exp(-k_1*(t2-t1)))/k_1*x1+...
  (1-exp(-k_2*(t2-t1)))/k_2*x2;
  
mean_x1 = @(t1,t2,x1) exp(-k_1*(t2-t1))*x1;
mean_x2 = @(t1,t2,x2) exp(-k_2*(t2-t1))*x2;
 
var_x1 = @(t1,t2) sigma_1^2/(2*k_1)*(1-exp(-2*k_1*(t2-t1)));
var_x2 = @(t1,t2) sigma_2^2/(2*k_2)*(1-exp(-2*k_2*(t2-t1)));

cov_x1_y = @(t1,t2) ro*sigma_1*sigma_2/(k_1+k_2)*(A2(t1,t2)+...
  1/k_2*(exp(-(k_1+k_2)*(t2-t1))-exp(-k_1*(t2-t1))))+...
  sigma_1^2/(2*k_1)*(A1(t1,t2)+1/k_1*(exp(-2*k_1*(t2-t1))-exp(-k_1*(t2-t1))));
  
cov_x2_y = @(t1,t2)  ro*sigma_1*sigma_2/(k_1+k_2)*(A1(t1,t2)+...
  1/k_1*(exp(-(k_1+k_2)*(t2-t1))-exp(-k_2*(t2-t1))))+...
  sigma_2^2/(2*k_2)*(A2(t1,t2)+1/k_2*(exp(-2*k_2*(t2-t1))-exp(-k_2*(t2-t1))));

cov_x1_x2 = @(t1,t2) ro*sigma_1*sigma_2/(k_1+k_2)*(1-exp(-(k_1+k_2)*(t2-t1)));

sigma_generate = @(t1,t2) [1 cov_x1_x2(t1,t2)./sqrt(var_x1(t1,t2)*var_x2(t1,t2)) ...
  cov_x1_y(t1,t2)./sqrt(var_x1(t1,t2)*var_y(t1,t2)); cov_x1_x2(t1,t2)./sqrt(var_x1(t1,t2)*var_x2(t1,t2)) ...
  1 cov_x2_y(t1,t2)./sqrt(var_x2(t1,t2)*var_y(t1,t2)); cov_x1_y(t1,t2)./sqrt(var_x1(t1,t2)*var_y(t1,t2)) ...
  cov_x2_y(t1,t2)./sqrt(var_x2(t1,t2)*var_y(t1,t2)) 1];

nominal = 100;
p_freq = 1/4;

payment_times = p_freq:p_freq:T;
reset_times = 0:p_freq:(T-p_freq);
var_phi_int = cumtrapz(t,var_phi);
var_phi_int = var_phi_int(1+p_freq/delta:p_freq/delta:length(t));
bond_prices = exp(-var_phi_int-mean_y(0,payment_times,0,0)+1/2*var_y(0,payment_times));
forward_bond_prices = bond_prices(2:end)./bond_prices(1:end-1);
forward_bond_prices = [bond_prices(1),forward_bond_prices];
forward_rates = (1./forward_bond_prices-1)/(p_freq);

floating_payments_value = sum(nominal*forward_rates*p_freq.*bond_prices);
fixed_rate = floating_payments_value/(nominal*p_freq*sum(bond_prices));
fixed_pay = forward_rates(1);
  
sim_x1 = zeros(1,length(t));
sim_x2 = zeros(1,length(t));
sim_y = zeros(1,length(t));

irs_value = zeros(1,length(t));
irs_value(1) = nominal*fixed_rate*p_freq*sum(bond_prices)-nominal*p_freq*sum(forward_rates.*bond_prices);

for i = 2:length(t)-1
  inter_t2 = t(i);
  inter_t1 = t(i-1);
  Sigma = sigma_generate(inter_t1,inter_t2);
  L = chol(Sigma,"lower");
  Z = randn(3,1);
  corr_Z = L*Z;
  
  sim_x1(i) = mean_x1(inter_t1,inter_t2,sim_x1(i-1))+sqrt(var_x1(inter_t1,inter_t2))*corr_Z(1,1);
  sim_x2(i) = mean_x2(inter_t1,inter_t2,sim_x2(i-1))+sqrt(var_x2(inter_t1,inter_t2))*corr_Z(2,1);
  sim_y(i) = mean_y(inter_t1,inter_t2,sim_x1(i-1),sim_x2(i-1))+...
    sqrt(var_y(inter_t1,inter_t2))*corr_Z(3,1);
 
   
  next_time = (inter_t2 <= payment_times);
  next_time = find(next_time,1);
  next_time = payment_times(next_time);
  var_phi_int = cumtrapz(t(i:end),var_phi(i:end));
  var_phi_int = var_phi_int(find(t(i:end)==next_time):p_freq/delta:length(t(i:end)));
  test_index = floor((i-2)/(p_freq/delta))+1;
  bond_prices = exp(-var_phi_int-mean_y(inter_t2,payment_times(test_index:end),sim_x1(i),sim_x2(i))+...
      1/2*var_y(inter_t2,payment_times(test_index:end)));
   
 is_reset = (inter_t2 == reset_times);
 
 forward_bond_prices = bond_prices(2:end)./bond_prices(1:end-1);
 forward_rates = (1./forward_bond_prices-1)/(p_freq);
 floating_payments_value = sum(nominal*forward_rates*p_freq.*bond_prices(2:end));

 if (max(is_reset == 1))
    fixed_pay = (1./bond_prices(2)-1)/(p_freq);
    fixed_payments_value = nominal*p_freq*sum(bond_prices(2:end))*fixed_rate;
 else
    floating_payments_value = floating_payments_value + bond_prices(1)*...
        nominal*p_freq*fixed_pay;
    fixed_payments_value = nominal*p_freq*sum(bond_prices)*fixed_rate;
 end
 
 irs_value(i) = floating_payments_value-fixed_payments_value;
 end  

payoff = sort(V_t*notional,2);
PFE = payoff(:,(1-alpha)*n);



  
  
 
 
 
 
 
 
 
 


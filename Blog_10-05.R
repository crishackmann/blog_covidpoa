# Código-fonte do Blog dia 05/10
# Estimativas de óbitos/COVID-19 atualizadas até 05/10/20
# https://www.ufrgs.br/covidpoa/?p=751

#### CSV PMPA UTI
# https://infografico-covid.procempa.com.br

require(deSolve)
require(ggplot2)
require(reshape2)
library(tidyverse)
library(lubridate)
require(grid)
options(stringsAsFactors = FALSE)
###########
arq_UTI_pmpa = "1_UTI_obt.csv"
arq_Hosp_pmpa = "2_Hosp_obt.csv"
arq_Death_pmpa ="3_Death_obt.csv"
arq_Death_new_pmpa ="6_Death_New_obt.csv"

caminho = "Dados_PMPA/"

UTI_pmpa = read.csv(file = paste(caminho, arq_UTI_pmpa, sep = ""))
hosp_pmpa = read.csv(file = paste(caminho, arq_Hosp_pmpa, sep = ""))
death_pmpa = read.csv(file = paste(caminho, arq_Death_pmpa, sep = ""), encoding = 'latin1')
colnames(death_pmpa)[2] = "Obitos.acumulado"
death_new_pmpa = read.csv(file = paste(caminho, arq_Death_new_pmpa, sep = ""))

index_iqual_days_pmpa = which(UTI_pmpa[,1] %in% hosp_pmpa[,1])

number_reported_full = UTI_pmpa[index_iqual_days_pmpa,2] + hosp_pmpa[,2]
number_reported = number_reported_full[-c(1:24)]
time = c(1:(length(number_reported)))
reported_data = data_frame(time, number_reported)

# Death
death_pmpa_iqual_days_model = death_pmpa[-c(1:11),]
death_pmpa_iqual_days_model$category = mdy(death_pmpa_iqual_days_model$category)

# Death new
death_new_pmpa_iqual_days_model = death_new_pmpa[-c(1:11),]
death_new_pmpa_iqual_days_model$category = mdy(death_new_pmpa_iqual_days_model$category)

#### Dates
date_ini = dmy("25/Abr/2020") 
num_days = 178 
series_days = date_ini + days(0:num_days)
times  = seq(0, num_days, by = 1)

num_days_hospi = length(number_reported)
num_days_death = length(death_pmpa_iqual_days_model$category)
date_ini_hospi = dmy("25/Abr/2020") 
###
# INPUT
sigma  = 1/5.1

gammah = 1/4.86 
gammad = 0 
gammar = 1/5.6

etad   = 1/11.98 
etar   = 1/20.22 

Theta  = 0.007 
Lambda = 0.195 

N = 1014009
I = N * 0.0025
E = 1.5 * I
H = 49
D = 12
R = 0
S = N - I - H - E - D - R

initial_state_values = c(S = S, E = E, I = I, H = H, D = D, R = R, incid_acum_I = I, exit_I = 0)

times = seq(from = 0, to = num_days, by = 1)

# SEIDR 
seidr_model = function(times, state, parameters) {  
  
  with(as.list(c(state, parameters)), {
    
    incid_I = (sigma * E)
    exit_I  = (Theta * gammah * I) + 
              ((1 - Theta) * (1 - Lambda) * gammar * I) + 
              ((1 - Theta) * Lambda *  gammad * I)
    
    dS = - (beta * I) * (S / N)
    
    dE =   (beta * I) * (S / N) - (sigma * E)
    
    dI =   incid_I - exit_I
    
    dH =   (Theta * gammah * I) - (Lambda * etad * H) - 
      ((1 - Lambda) * etar * H)
    
    dD =   ((1 - Theta) * Lambda *  gammad * I) + 
      (Lambda * etad * H) 
    
    dR =   ((1 - Theta) * (1 - Lambda) * gammar * I) + 
      ((1 - Lambda) * etar * H)
    
    return(list(c(dS, dE, dI, dH, dD, dR, incid_I, exit_I))) 
  })
}

# Distance Function
loglikelihood_fun = function(parameters, dat) { 
  
  beta = parameters[1]
  
  output = as.data.frame(ode(y = initial_state_values, 
                             times = times, 
                             func = seidr_model,
                             parms = c(beta = beta)))  
  
  # Calculate log-likelihood 
  LL = -sum(dpois(x = dat$number_reported, lambda = output$D[output$time %in% dat$time], log = TRUE))
  
  return(LL) 
}

# Optimization
param_otimim = optim(par = c(0.4),
                     fn = loglikelihood_fun,
                     dat = reported_data, 
                     control = list(fnscale = 1),
                     method = 'Brent',
                     lower = 0,
                     upper = 2) 

beta_estimated = round(param_otimim$par[1], 4)
######################
parameters = c(beta   = beta_estimated, 
               sigma  = sigma,
               gammah = gammah, 
               gammad = gammad, 
               gammar = gammar,
               etad   = etad,
               etar   = etar, 
               Theta  = Theta,
               Lambda = Lambda)

output = as.data.frame(ode( y     = initial_state_values, 
                            times = times, 
                            func  = seidr_model,
                            parms = parameters))

# New Column incidence D
output$incid_D = c(0,diff(output$D, lag = 1))

### Dates
# Change time index
output$time = series_days
time = date_ini_hospi + days(1:num_days_hospi - 1)
reported_data = data.frame(time, number_reported)

last_day_hosp = tail(time, n = 1)
d = day(last_day_hosp)
m = months(last_day_hosp, abbreviate = T)
y = year(last_day_hosp)
day_month_year_last_day_hosp = paste(d,"/",m,"/",y)

##### Resíduos Óbitos
out_res_death = output[output$time %in% death_pmpa_iqual_days_model$category, output_column(c('time','D'))]
res_death = death_pmpa_iqual_days_model$Obitos.acumulado - out_res_death$D
res_quad_death = res_death^2
sum_res_death = sum(res_quad_death)/num_days_death
rmse_death = round(sqrt(sum_res_death),2)

# Plot the model fit D New
base_D_new = ggplot() +
  geom_bar(stat = "identity", data = output, aes(x = time, y = incid_D), colour = '#119299', size = 2) +
  geom_bar(stat = "identity", data = death_new_pmpa_iqual_days_model, aes(x = category, y = Obitos.por.dia, colour = "Número observado de Novos Casos/COVID-19 por Dia"), size = 2) +
  xlab("Data") +
  ylab("Número previsto de novos óbitos/COVID-19 por dia") +                                 
  labs(title = paste("Modelo Calibrado até",day_month_year_last_day_hosp,"."), colour = "") +
  theme(axis.title.x = element_text(size = 16), axis.text.x = element_text(angle = 30, vjust = 0.5, size = 13)) +
  theme(axis.title.y = element_text(size = 16), axis.text.y = element_text(angle = 0, vjust = 0.5, size = 13)) +
  theme(legend.position = "bottom") +
  theme(text = element_text(size = 18))+
  coord_cartesian(ylim = c(0, 20))

 base_D_new + scale_x_date(date_breaks = '2 week', date_labels = "%d %b")
# End D new

# Plot the model fit Death
base_death = ggplot() +
  geom_line(data = output, aes(x = time, y = D), colour = '#119299', size = 2) +
  geom_point(data = death_pmpa_iqual_days_model, aes(x = category, y = Obitos.acumulado, colour = "Número observado de óbitos acumulados/COVID-19"), size = 2) +
  xlab("Data") +
  ylab("Número previsto de óbitos acumulados/COVID-19") +                                 
  labs(title = paste("Modelo Calibrado até",day_month_year_last_day_hosp,"."), colour = "") +
  theme(axis.title.x = element_text(size = 16), axis.text.x = element_text(angle = 30, vjust = 0.5, size = 13)) +
  theme(axis.title.y = element_text(size = 16), axis.text.y = element_text(angle = 0, vjust = 0.5, size = 13)) +
  theme(legend.position = "bottom") +
  theme(text = element_text(size = 18))

base_death + scale_x_date(date_breaks = '2 week', date_labels = "%d %b")

# End Deaths

# R0 (next-generation matrix) - Ini

F1 = quote(beta * S * I / N)
F2 = 0
F3 = 0

###################################################
Vm1 = quote(sigma * E)
Vm2 = quote(Theta * gammah * I + (1 - Theta) * 
              (1 - Lambda) * gammar * I + (1 - Theta) * 
              Lambda * gammad * I)
Vm3 = quote(Lambda * etad * H + (1 - Lambda) * etar * H)

###################################################
Vp1 = 0
Vp2 = quote(sigma * E)
Vp3 = quote(Theta * gammah * I)

###################################################
V1 = substitute(a - b, list(a = Vm1, b = Vp1))
V2 = substitute(a - b, list(a = Vm2, b = Vp2))
V3 = substitute(a - b, list(a = Vm3, b = Vp3))

###################################################
f11 = D(F1, "E"); f12 = D(F1, "I"); f13 = D(F1, "H")
f21 = D(F2, "E"); f22 = D(F2, "I"); f23 = D(F2, "H")
f31 = D(F3, "E"); f32 = D(F3, "I"); f33 = D(F3, "H")

v11 = D(V1, "E"); v12 = D(V1, "I"); v13 = D(V1, "H")
v21 = D(V2, "E"); v22 = D(V2, "I"); v23 = D(V2, "H")
v31 = D(V3, "E"); v32 = D(V3, "I"); v33 = D(V3, "H")

###################################################

paras = list(N = N, S = N, E = 0, I = 0, 
             H = 0, D = 0, R = 0, 
             beta = beta_estimated,
             sigma = sigma, Theta = Theta, 
             Lambda = Lambda, gammah = gammah,
             gammad = gammad, gammar = gammar, 
             etad = etad, etar = etar)

f = with(paras, 
         matrix(c(eval(f11), eval(f12), eval(f13),
                  eval(f21), eval(f22), eval(f23),
                  eval(f31), eval(f32), eval(f33)),
                nrow = 3, byrow = T))

v = with(paras, 
         matrix(c(eval(v11), eval(v12), eval(v13),
                  eval(v21), eval(v22), eval(v23),
                  eval(v31), eval(v32), eval(v33)),
                nrow = 3, byrow = T))
###################################################
R0 = max(eigen(f %*% solve(v))$values)
# R0 (next-generation matrix) - End
###################################################
# R effective
Re = R0 * (output$S[output$time == last_day_hosp]/N)

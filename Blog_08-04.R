# Código-fonte do post do dia 04 de agosto de 2020
# Estimativas atualizadas até o dia 03/08/2020 e Ponto de Inflexão


#### CSV PMPA UTI
# https://infografico-covid.procempa.com.br/#section_1

require(deSolve)
require(ggplot2)
require(reshape2)
library(tidyverse)
library(lubridate)
require(grid)
options(stringsAsFactors = FALSE)
###########
arq_UTI_pmpa = "qqtc4yor.csv"
arq_Hosp_pmpa = "6o99iyrh.csv"

caminho = "Dados_PMPA/"

UTI_pmpa = read.csv(file = paste(caminho, arq_UTI_pmpa, sep = ""))
hosp_pmpa = read.csv(file = paste(caminho, arq_Hosp_pmpa, sep = ""))

index_iqual_days_pmpa = which(UTI_pmpa[,1] %in% hosp_pmpa[,1])

number_reported_full = UTI_pmpa[index_iqual_days_pmpa,2] + hosp_pmpa[,2]
number_reported = number_reported_full[-c(1:24)]
time = c(1: (length(number_reported)))
reported_data = data_frame(time, number_reported)

### Leitos
total_leitos_covid = 383
prop_UTI_H = 0.415
#### Dates
date_ini = dmy("25/Abr/2020") 
num_days = 300 # mais o primeiro dia
series_days = date_ini + days(0:num_days)
times  = seq(0, num_days, by = 1)

num_days_hospi = length(number_reported)
date_ini_hospi = dmy("25/Abr/2020") 
###
# INPUT
sigma  = 1/5.1
gammah = 1/4.86 
gammad = 0 
gammar = 1/5.6
etad   = 1/12.98 
etar   = 1/10.22  
Theta  = 0.01 
Lambda = 0.175

########
# Comentar/Descomentar para as estimativas do
# cenário 1 ou 2
########
# População adulta de Porto Alegre - Cenário 1
N = 1014009
# População de Porto Alegre - Cenário 2
#N = 1483771
  
I = N * 0.0013 
E = 1.5 * I 
H = 49
S = N - I - H - E
D = 12
R = 0
initial_state_values = c(S = S, E = E, I = I, H = H, D = D, R = R)

times = seq(from = 0, to = num_days, by = 1)

# SEIDR 
seidr_model = function(times, state, parameters) {  
  
  with(as.list(c(state, parameters)), {
    
    
    dS = - (beta * I) * (S / N)
    
    dE =   (beta * I) * (S / N) - (sigma * E)
    
    dI =   (sigma * E) - (Theta * gammah * I) - 
      ((1 - Theta) * (1 - Lambda) * gammar * I) - 
      ((1 - Theta) * Lambda *  gammad * I)
    
    dH =   (Theta * gammah * I) - (Lambda * etad * H) - 
      ((1 - Lambda) * etar * H)
    
    dD =   ((1 - Theta) * Lambda *  gammad * I) + 
      (Lambda * etad * H) 
    
    dR =   ((1 - Theta) * (1 - Lambda) * gammar * I) + 
      ((1 - Lambda) * etar * H)
    
    return(list(c(dS, dE, dI, dH, dD, dR))) 
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
  LL = -sum(dpois(x = dat$number_reported, lambda = output$H[output$time %in% dat$time], log = TRUE))
  
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

## Inflection Points
acceleration = c(0,0,0,diff(sign(diff(diff(output$H)))))
inflection = data.frame(series_days,acceleration)
inflection_days_index = which(inflection$acceleration != 0)

inflection_days = data.frame(output[inflection_days_index,])
keeps = c("time","H")
inflection_days = inflection_days[keeps]
inflection_days = inflection_days[1,]
###

##### UTI
UTI_df = data_frame(UTI_pmpa[-c(1:37),])   
date_ini_UTI_PMPA = dmy("25/Abr/2020")
num_days_UTI_pmpa = length(UTI_df$category)
UTI_days_PMPA = date_ini_UTI_PMPA + days(1:num_days_UTI_pmpa - 1)
UTI_df$category = UTI_days_PMPA

last_day_UTI_PMPA = tail(UTI_days_PMPA, n = 1)
d = day(last_day_UTI_PMPA)
m = months(last_day_UTI_PMPA, abbreviate = T)
y = year(last_day_UTI_PMPA)
day_month_year_last_day_UTI_PMPA = paste(d,"/",m,"/",y)

UTI_SEIHDR_prop_UTI_H = output$H * prop_UTI_H
UTI_SEIHDR = data_frame(series_days,UTI_SEIHDR_prop_UTI_H)

# Closest day = max UTI
max_UTI_index = which.max(UTI_SEIHDR$UTI_SEIHDR_prop_UTI_H)
closest_index_max_UTI = which.min(abs(UTI_SEIHDR$UTI_SEIHDR_prop_UTI_H[1:max_UTI_index] - total_leitos_covid))
closest_max_UTI = output[closest_index_max_UTI,]
closest_day_max_UTI = closest_max_UTI[1,1]

d = day(closest_day_max_UTI)
m = months(closest_day_max_UTI, abbreviate = T)
y = year(closest_day_max_UTI)

day_format_closest = paste(d,"/",m,"/",y)
day_format_closest = dmy(day_format_closest)

base_UTI = ggplot() +
  geom_line(data = UTI_SEIHDR, aes(x = series_days, y = UTI_SEIHDR_prop_UTI_H), colour = '#56B4E9', size = 2) +
  geom_point(data = UTI_df, aes(x = category, y = Leitos.UTI, colour = "Número observado de pacientes em UTIs/COVID-19"), size = 2) +
  xlab("Data") +
  ylab("Número previsto de UTIs/COVID-19") +                                 
  labs(title = paste("Modelo Calibrado até",day_month_year_last_day_hosp,".\n População Suscetível =",N), colour = "") +
  theme(axis.title.x = element_text(size = 16), axis.text.x = element_text(angle = 30, vjust = 0.5, size = 13)) +
  theme(axis.title.y = element_text(size = 16), axis.text.y = element_text(angle = 0, vjust = 0.5, size = 13)) +
  theme(legend.position = "bottom") +
  theme(text = element_text(size = 18))+
  geom_hline(yintercept = total_leitos_covid, colour = "orange", size = 1.5, linetype = "dashed")+
  geom_vline(xintercept = day_format_closest, colour = "green", size = 1.5, linetype = "dashed")+
  annotate("text", x = dmy("10/jun/2020"), y = total_leitos_covid + 15, label =  paste("Máximo de leitos UTI Porto Alegre/RS,",total_leitos_covid), colour = "orange")+
  geom_text()

base_UTI + scale_x_date(date_breaks = '2 week', date_labels = "%d %b")
# Plot the model fit
base_1 = ggplot() +
  geom_line(data = output, aes(x = time, y = H), colour = '#119299', size = 2) +
  geom_point(data = reported_data, aes(x = time, y = number_reported, colour = "Número observado de pacientes COVID-19 hospitalizados"), size = 2) +
  xlab("Data") +
  ylab("Número previsto de hospitalizações/COVID-19") +                                 
  labs(title = paste("Modelo Calibrado até",day_month_year_last_day_hosp,".\n População Suscetível =",N), colour = "") +
  theme(axis.title.x = element_text(size = 16), axis.text.x = element_text(angle = 30, vjust = 0.5, size = 13)) +
  theme(axis.title.y = element_text(size = 16), axis.text.y = element_text(angle = 0, vjust = 0.5, size = 13)) +
  theme(legend.position = "bottom") +
  theme(text = element_text(size = 18)) +
  geom_point(data = inflection_days, aes(x = time, y = H), size = 2) +
  annotate("text", x = inflection_days$time - days(30), y = inflection_days$H + 220, label = "Ponto de Inflexão" , colour = "black", size = 5)+
  geom_segment(aes(x = inflection_days$time - days(20), y = inflection_days$H + 200, xend = inflection_days$time - days(2), yend = inflection_days$H + 20),
               arrow = arrow(length = unit(0.5, "cm")))

base_1 + scale_x_date(date_breaks = '2 week', date_labels = "%d %b")
######################
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

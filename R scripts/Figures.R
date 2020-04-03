library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis)


# read data ##########################################################################
## data from cows
setwd("/Users/kristinaceres/Desktop/CT_HMM/data")
data = read.csv("hmm_ready_data.csv")

#model results
# 2 state model data 
seed = 1585510210

scale_out2 <- read.csv(paste("scale_out2",seed, ".csv", sep=""))
shape_out2 <- read.csv(paste("shape_out2",seed, ".csv", sep=""))
qmat_out2 <- read.csv(paste("qmat_out2",seed, ".csv", sep=""))
aic2 <- read.csv(paste("aic_2",seed, ".csv", sep=""))
ip2 <- read.csv(paste("ip_2",seed, ".csv", sep=""))
ip_out2 <- read.csv(paste("ip_out2",seed, ".csv", sep=""))
scale_in2 <- read.csv(paste("scale_2",seed, ".csv", sep=""))
shape_in2 <- read.csv(paste("shape_2",seed, ".csv", sep=""))
qmat_in2 <- read.csv(paste("qmat_2",seed, ".csv", sep=""))

# 3 state model data
scale_out3 <- read.csv(paste("scale_out3", seed, ".csv", sep=""))
shape_out3 <- read.csv(paste("shape_out3", seed, ".csv", sep=""))
qmat_out3 <- read.csv(paste("qmat_out3", seed, ".csv", sep=""))
aic3 <- read.csv(paste("aic_3", seed, ".csv", sep=""))
ip3 <- read.csv(paste("ip_3", seed, ".csv", sep=""))
ip_out3 <- read.csv(paste("ip_out3", seed, ".csv", sep=""))
scale_in3 <- read.csv(paste("scale_3", seed, ".csv", sep="")) 
shape_in3 <- read.csv(paste("shape_3", seed, ".csv", sep="")) 
qmat_in3 <- read.csv(paste("qmat_3", seed, ".csv", sep="")) 


# format data 
#seed 500 runs state 2 1583360998
get_result_df = function(df){
  df = as.matrix(df)
  df = t(df)
  colnames(df) = df[1,]
  df = df[-1, ] 
  df = as_tibble(df)
  return(df)
}

#transform estimated parameters
# 2 states 
scale_out2 = get_result_df(scale_out2) 
shape_out2 = get_result_df(shape_out2) 
qmat_out2 = get_result_df(qmat_out2)
aic2 = get_result_df(aic2)
ip_out2 = get_result_df(ip_out2)

outputs2 = bind_cols(ip_out2, shape_out2, scale_out2, qmat_out2, aic2)
colnames(outputs2) = c("ip_out0", "ip_out1", "shape_out0", "shape_out1",
                       "scale_out0", "scale_out1", "qout0", "qout1", "qout2", "qout3", "aic")

#transform input parameters
shape_in2 = get_result_df(shape_in2) 
scale_in2 = get_result_df(scale_in2)
ip2 = get_result_df(ip2)
qmat_in2 = get_result_df(qmat_in2)

inputs2 = bind_cols(ip2, shape_in2, scale_in2, qmat_in2)
colnames(inputs2) = c("ip0", "ip1", "shape0", "shape1", "scale0", "scale1", 
                      "qin0", "qin1", "qin2", "qin3")

#all parameters
allparam2 = bind_cols(inputs2, outputs2)

#################################################################################################
# 3 states
scale_out3 = get_result_df(scale_out3) 
shape_out3 = get_result_df(shape_out3) 
qmat_out3 = get_result_df(qmat_out3)
aic3 = get_result_df(aic3)
ip_out3 = get_result_df(ip_out3)

outputs3 = bind_cols(ip_out3, shape_out3, scale_out3, qmat_out3, aic3)
colnames(outputs3) = c("ip_out0", "ip_out1", "ip_out2", "shape_out0", "shape_out1", "shape_out2", "scale_out0", "scale_out1", "scale_out2", "qout0", "qout1", "qout2", "qout3", "qout4","qout5","qout6", "qout7", "qout8", "aic")

#transform input parameters
shape_in3 = get_result_df(shape_in3) 
scale_in3 = get_result_df(scale_in3)
ip3 = get_result_df(ip3)
qmat_in3 = get_result_df(qmat_in3)

inputs3 = bind_cols(ip3, shape_in3, scale_in3, qmat_in3)
colnames(inputs3) = c("ip0", "ip1", "ip2", "shape0", "shape1", "shape2", "scale0", "scale1", "scale2", 
                      "qin0", "qin1", "qin2", "qin3","qin4","qin5","qin6","qin7","qin8")

#all parameters
allparam3 = bind_cols(inputs3, outputs3)

###########################################################################################################
## data summary plots ######################################################################################
data$SampDate = as.POSIXct(strptime(data$SampDate, format="%Y-%m-%d"))
p1data = data %>% group_by(CombinedID) %>% mutate(firstdate = min(SampDate)) %>% 
  mutate(time2 =  difftime((SampDate), (firstdate), units = "weeks")) %>% mutate(SampTime = row_number())

idlist= c("A_1181", "C_106", "C_114", "A_1180", "A_1179")
subset = p1data %>% filter(CombinedID %in% idlist)
p1 = ggplot(data = subset) + 
  geom_col(aes(x=SampTime, y = cor_totCFU, group=CombinedID, fill=CombinedID), position=position_dodge2(width=.2, preserve = "single"))+
  scale_fill_viridis(discrete=T)+labs(x= "Sampling Interval", y = "log(CFU) MAP", fill="Cow ID")
setwd("../Figures")
ggsave("Fig1.tiff", p1, dpi=500)

# example distribution plots 
shape3<- c(19.130472109595726, 16.774795300202804, 198.4831018750443)
scale3<- c(0.12568632084200168, 0.2815762121436809, 0.04040960031804022)
d0_ex =rgamma(n = 10000, shape=shape3[1], scale=scale3[1])
d1_ex =rgamma(n = 10000, shape=shape3[2], scale=scale3[2])
d2_ex =rgamma(n = 10000, shape=shape3[3], scale=scale3[3])
gam_ex  = tibble(d0_ex, d1_ex, d2_ex)
gam_ex = gam_ex %>% gather(key = "distribution", value = "density")
p2 = ggplot()+
  geom_density(data = gam_ex, aes(density, fill = distribution), alpha = .5) +
  labs(x = "log(CFU) Value")+ scale_fill_viridis(discrete=T, labels=c("State 0", "State 1", "State 2"))+
  theme(legend.position="bottom")
p2
ggsave("Example_plot.png", p2, dpi = 500, width = 6.6, height = 4)


###############################################################################################
## qmat graphs
n=nrow(data)
outputs2 = outputs2 %>% filter(aic>0) %>% arrange(aic)
outputs3 = outputs3 %>% filter(aic > 0) %>% arrange(aic)
ip2 = c(outputs2[1,]$ip_out0, outputs2[1,]$ip_out1)
qmat2 = matrix(c(outputs2[1,]$qout0,outputs2[1,]$qout1,
                 outputs2[1,]$qout2,outputs2[1,]$qout3), nrow=2, byrow=T)


ip3 = c(outputs3[1,]$ip_out0, outputs3[1,]$ip_out1, outputs3[1,]$ip_out2)
qmat3 = matrix(c(outputs3[1,]$qout0,outputs3[1,]$qout1,outputs3[1,]$qout2,
                 outputs3[1,]$qout3,outputs3[1,]$qout4,outputs3[1,]$qout5,
                 outputs3[1,]$qout6,outputs3[1,]$qout7,outputs3[1,]$qout8), nrow=3, byrow=T)

n = 300
# make data frames 
q2 = tibble(x = seq(1,n, length=n), "0" = rep(0, n), "1" = rep(0,n), group = rep("2 State", n))
q3 = tibble(x = seq(1,n, length=n), "0" = rep(0, n), "1"= rep(0,n), "2" = rep(0,n), group = rep("3 State",n))

# get state probabilities for each position in the chain
for (i in 1:n){
  temp2 = ip2 %*% expm(qmat2 * q2$x[i])
  q2$`0`[i]=temp2[1]
  q2$`1`[i]=temp2[2]
  
  temp3 = ip3 %*% expm(qmat3 * q3$x[i])
  q3$`0`[i]=temp3[1]
  q3$`1`[i]=temp3[2]
  q3$`2`[i]=temp3[3]
  
}


q2 = gather(q2, -c(x, group), key= "State", value="Probability")
q3 = gather(q3, -c(x, group), key= "State", value="Probability")

q = bind_rows(q2, q3, q4, q5)
q = bind_rows(q2, q3)

p3 = ggplot(data = q)+geom_line(aes(x=x, y=Probability, color=State))+facet_wrap(~group, nrow=1)+
  labs(x = "t", y =expression(pi*e^(Qt)))+
  scale_color_manual(values=c("#874E9A","#50AAA3","#BAC438"))

ggsave("Fig3.tiff", p3, dpi=500, width= 8, height = 4)

## gamma distribution plots ################################################################################
n = nrow(data)
shape2 <- c(outputs2[1,]$shape_out0, outputs2[1,]$shape_out1)
scale2 <- c(outputs2[1,]$scale_out0, outputs2[1,]$scale_out1)
d02 =rgamma(n = round(q2$dist[q2$State == 0 & q2$x == 250]*n), shape=shape2[1], scale=scale2[1])
d12 =rgamma(n = round(q2$dist[q2$State == 1 & q2$x == 250]*n), shape=shape2[2], scale=scale2[2])

shape3 <- c(outputs3[1,]$shape_out0, outputs3[1,]$shape_out1, outputs3[1,]$shape_out2)
scale3 <- c(outputs3[1,]$scale_out0, outputs3[1,]$scale_out1, outputs3[1,]$scale_out2)
d03 =rgamma(n = round(q3$dist[q3$State == 0 & q3$x == 250] * n), shape=shape3[1], scale=scale3[1])
d13 =rgamma(n = round(q3$dist[q3$State == 1 & q3$x == 250] * n), shape=shape3[2], scale=scale3[2])
d23 =rgamma(n = round(q3$dist[q3$State == 2 & q3$x == 250] * n), shape=shape3[3], scale=scale3[3])

shape4 <- c(outputs4[1,]$shape_out0, outputs4[1,]$shape_out1, outputs4[1,]$shape_out2, outputs4[1,]$shape_out3)
scale4 <- c(outputs4[1,]$scale_out0, outputs4[1,]$scale_out1, outputs4[1,]$scale_out2, outputs4[1,]$scale_out3)
d04 =rgamma(n = round(q4$Probability[q4$State == 0 & q4$x == 250] * n), shape=shape4[1], scale=scale4[1])
d14 =rgamma(n = round(q4$Probability[q4$State == 1 & q4$x == 250] * n), shape=shape4[2], scale=scale4[2])
d24 = rgamma(n = round(q4$Probability[q4$State == 2 & q4$x == 250] * n), shape=shape4[3], scale=scale4[3])
d34 = rgamma(n = round(q4$Probability[q4$State == 3 & q4$x == 250] * n), shape=shape4[4], scale=scale4[4])

shape5 <- c(outputs5[1,]$shape_out0, outputs5[1,]$shape_out1, outputs5[1,]$shape_out2, outputs5[1,]$shape_out3, outputs5[1,]$shape_out4)
scale5 <- c(outputs5[1,]$scale_out0, outputs5[1,]$scale_out1, outputs5[1,]$scale_out2, outputs5[1,]$scale_out3, outputs5[1,]$scale_out4)
d05 =rgamma(n = round(q5$Probability[q5$State == 0 & q5$x == 250] * n), shape=shape5[1], scale=scale5[1])
d15 =rgamma(n = round(q5$Probability[q5$State == 1 & q5$x == 250] * n), shape=shape5[2], scale=scale5[2])
d25 = rgamma(n = round(q5$Probability[q5$State == 2 & q5$x == 250] * n), shape=shape5[3], scale=scale5[3])
d35 = rgamma(n = round(q5$Probability[q5$State == 3 & q5$x == 250] * n), shape=shape5[4], scale=scale5[4])
d45 = rgamma(n = round(q5$Probability[q5$State == 4 & q5$x == 250] * n), shape=shape5[5], scale=scale5[5])

gammas20 = tibble(dist = d02, State = rep("0", length(d02)), group = rep("2 State", length(d02)))
gammas21 = tibble(dist = d12, State = rep("1", length(d12)), group = rep("2 State", length(d12)))

gammas30 = tibble(dist = d03, State = rep("0", length(d03)), group = rep("3 State", length(d03)))
gammas31 = tibble(dist = d13, State = rep("1", length(d13)), group = rep("3 State", length(d13)))
gammas32 = tibble(dist = d23, State = rep("2", length(d23)), group = rep("3 State", length(d23)))

gammas40 = tibble(dist = d04, State = rep("0", length(d04)), group = rep("4 State", length(d04)))
gammas41 = tibble(dist = d14, State = rep("1", length(d14)), group = rep("4 State", length(d14)))
gammas42 = tibble(dist = d24, State = rep("2", length(d24)), group = rep("4 State", length(d24)))
gammas43 = tibble(dist = d34, State = rep("3", length(d34)), group = rep("4 State", length(d34)))

gammas50 = tibble(dist = d05, State = rep("0", length(d05)), group = rep("5 State", length(d05)))
gammas51 = tibble(dist = d15, State = rep("1", length(d15)), group = rep("5 State", length(d15)))
gammas52 = tibble(dist = d25, State = rep("2", length(d25)), group = rep("5 State", length(d25)))
gammas53 = tibble(dist = d35, State = rep("3", length(d35)), group = rep("5 State", length(d35)))
gammas54 = tibble(dist = d45, State = rep("4", length(d45)), group = rep("5 State", length(d45)))

gammas = bind_rows(gammas20, gammas21, gammas30, gammas31, gammas32, gammas40, gammas41, gammas42, 
                   gammas43, gammas50, gammas51, gammas52, gammas53, gammas54)

q$type = rep("stationary", nrow(q))
gammas$type = rep("emiss", nrow(gammas))

gammas = bind_rows(gammas20, gammas21, gammas30, gammas31, gammas32)
p4 = ggplot()+
    geom_histogram(data = data, aes(x=cor_totCFU), binwidth=.2, alpha = .7)+
    geom_histogram(data = gammas, aes(x=dist, fill = State), binwidth = .2, alpha = .7)+
    facet_wrap(~group, nrow=1)+
    scale_fill_manual(values=c("#874E9A","#50AAA3","#BAC438"))+
    labs(x= "log(CFU) MAP", y = "Number of observations")
p4

ggsave("Fig4.tiff", p4, dpi=500, width = 8, height= 4)

############################################################################################################
## posterior probabilities 

time_data = data %>% mutate(time = ifelse(is.na(time) == T, 0, time)) %>% 
  group_by(CombinedID) %>% mutate(times = cumsum(time)) %>% select(CombinedID, times) %>%
  rename(CowID = CombinedID) %>% arrange(CowID)

# 2 state 

max_post_plot = function(filename, time_data, nstates){
  posterior = read.csv(filename)
  colnames(posterior) = c("CowID", 1,2,3,4,5,6,7,8,9,10,11)
  posterior = gather(posterior, key="SampleTime", value = "prob", -c(CowID)) %>% filter(is.na(prob)==F)
  posterior = posterior %>% arrange(CowID) 
  p = bind_cols(posterior, time_data)
  p$nstates=as.factor(rep(nstates, length=nrow(p)))
  
  return(p)
}
setwd("../data/")

p0_2 = read_posterior("posterior0_df_2.csv", time_data, 2, 0)
p1_2 = read_posterior("posterior1_df_2.csv", time_data, 2, 1)
p2 = bind_rows(p0_2, p1_2)
p2 = p2  %>% group_by(CowID, state)  %>%
  mutate(Endstate = case_when(
    state==0 & last(prob) > 0.75 ~ "0",
    state== 1 & last(prob) > 0.75 ~ "1",
    state== 0 & last(prob) < 0.75 ~ "1",
    state== 1 & last(prob) < 0.75 ~ "0"
  ))

p0_3 = read_posterior("posterior0_df_3.csv", time_data, 3, 0)
p1_3 = read_posterior("posterior1_df_3.csv", time_data, 3, 1)
p2_3 = read_posterior("posterior2_df_3.csv", time_data, 3, 2)
p3 = bind_rows(p0_3, p1_3, p2_3)
p3 = p3 %>%  group_by(CowID, state) %>% mutate(max_last_prob = max(last(prob))) %>% 
    group_by(CowID) %>% mutate(max_end_prob = max(max_last_prob)) %>% 
    group_by(CowID, state)  %>%
      mutate(Endstate = case_when(
        state==0 & last(prob) == max_end_prob ~ "0",
        state==1 & last(prob) == max_end_prob ~ "1",
        state==2 & last(prob) == max_end_prob ~ "2")) %>% 
      group_by(CowID) %>% 
      mutate(Endstate = case_when(
        "0" %in% Endstate ~ "0",
        "1" %in% Endstate ~ "1",
        "2" %in% Endstate ~ "2"
      )) 

post= bind_rows(p2,p3)

post = post %>% group_by(CowID, nstates, times) %>% mutate(max_prob = ifelse(prob == max(prob), state, -1)) 
post = post %>% group_by(CowID, nstates, times) %>% 
  mutate(max_prob = ifelse(max_prob == max(max_prob), max_prob, max(max_prob)))
post = post %>% group_by(CowID, nstates)%>% mutate(max_prob = as.numeric(max_prob)+runif(1, -.05,.05)) %>%
  mutate(nstates2 = ifelse(nstates == "2" , "2 State", "3 State"))


p5 = ggplot(post)+geom_line(aes(x=times, y = max_prob, group=CowID, color=Endstate)) +
  facet_grid(vars(nstates2))+ scale_color_manual(values=c("#874E9A","#50AAA3","#BAC438"))+
  labs(color = "End State", x = "Time", y = "Maximum posterior probability")
p5

setwd("../Figures")
ggsave("Fig4.tiff", p5, dpi = 500, width = 6, height = 5)


## sojourn times
-1/qmat2[1,1]
-1/qmat2[2,2]

-1/qmat3[1,1]
-1/qmat3[2,2]





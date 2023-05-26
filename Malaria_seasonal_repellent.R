# This ODE model simulates the transmission of malaria between mosquito and human populations. 
# 
# Seasonality impacts 1) carrying capacity for mosquito population (K)
#                     2) incubation period of malaria (xi_m)
#                     3) death rate of mosquito (mu_m)
# Outputs include plots of seasonal dynamics of malaria cases and mosquito population

# by Hannah Meredith
# last updated October 4, 2019


# libraries

library("pracma")
library("deSolve")
library("ggplot2")
library("readxl")
library("grid")
library("gridExtra")
library("reshape2")
library("ggpubr")
library("DescTools")

### Baseline function: This function simulates what would happen if no controls (bednets, repellents, etc.) are present

baseline <- function(y0, P) {
  
  P = c(
    P,
    c_D = 0,  # no coverage in baseline
    hl_D = 1  # halflife should not have an effect here, but useful to have a place holder for troubleshooting
  )

  step = 1  # step size = 1 day
  t = seq(0, P["t_ss"] + P["treatment_pd"], by = step) #duration of simulation = t_ss + treatment period
  
  base <-
    rk(
      y = y0,
      times = t,
      func = biting,
      parms = P,
      method = "ode45",
      atol = 1e-10,
      rtol = 1e-10
    )

  total_infected_h <- base[ , "I_h"] + base[,"A_h"] # calculate total infected = infected + asymptomatic
  time <- base [, 1]
  
  baseline <- list(time, total_infected_h)
}

### Dosing function : This function introduces control methods at set dose efficies (dose_eff) and intervals (period_D)

dosing <- function(y0, P) {
  
  # initialize system. Let it run to steady state before introducing control methods
  t_steadyState <- P[["t_ss"]]         # length of time to achieve steady state
  treatment_pd <- P[["treatment_pd"]]  # treatment period (years)
  step <- 1
  times = seq(0, t_steadyState, by = step)
  out <-                                         # run ODE model for until steady state to let system reach equilibrium before introducing control methods
    rk(
      y = y0,  # y0
      times = times,
      func = biting,
      parms = P,
      method = "ode45",
      atol = 1e-10,
      rtol = 1e-10
    )
  
  # After system reaches steady state, introduce control methods
  dose_eff <- P[["dose_eff"]]    
  period_D <- P[["period_D"]]
  
  if (period_D <= step) {
    step <- period_D/2
  }
  
  # Define matrix for storage and new time windows. 
  overall.df <- as.data.frame(out)
  t_startD = tail(overall.df$time,1)
  t_endD = tail(overall.df$time,1) + period_D
  
  # Update initial conditions for next ODE solver. Add first round of repellent
  y0 <- tail(overall.df[,-1],1)
  y0$D <- dose_eff 
  y0 = unlist(y0)  # change y0 from data.frame to named num (necessary for desolve)
  
  out <-
    rk(
      y = y0,
      times = seq(t_startD, t_endD, by = step),
      func = biting,
      parms = P,
      method = "ode45",
      atol = 1e-10,
      rtol = 1e-10
    )
  
  # Update timers and dataframe
  overall.df <- rbind(overall.df, as.data.frame(out))
  t_startD = tail(overall.df$time, 1)
  t_endD = tail(overall.df$time, 1) + period_D
 
  # Update counters
  doseNum = 1
  toggle = 0
  
  # Run simulation for set treatment period (here, 1 year starting after steady state)
  while (tail(overall.df$time, 1) < t_steadyState + treatment_pd &
         toggle == 0) {
    
    # Ensures simulation cuts off at 1 year
    if (t_endD >= t_steadyState + treatment_pd) {
      t_endD = t_steadyState + treatment_pd
      toggle = 1
    } 
     if (abs(t_endD - t_startD) < 0.001){   # ends simulation if t_start and t_end for last round are too close
       dosing <- list(overall.df$time, overall.df$TI_h, overall.df$D)
       break
     }
    
      y0 <- tail(overall.df[,-1],1)
      y0$D <- dose_eff
      y0 = unlist(y0)         # change y0 from data.frame to named num (necessary for desolve)
      
      out <-
        rk(
          y = y0,
          times = seq(t_startD, t_endD, by = step),
          func = biting,
          parms = P,
          method = "ode45",
          atol = 1e-10,
          rtol = 1e-10
        )
      
      overall.df <- rbind(overall.df, as.data.frame(out))
      t_startD = tail(overall.df$time, 1)
      t_endD = tail(overall.df$time, 1) + period_D
      
      # Update counters
      doseNum = doseNum + 1
      
      # ggplot(subset(overall.df, time > t_steadyState), aes(x = time))+
      #   geom_line(aes(y = D)) +
      #   geom_line(aes(y = I_h))
      
  }
  
  overall.df$TI_h <- overall.df$I_h + overall.df$A_h
  dosing <- list(overall.df$time, overall.df$TI_h, overall.df$D)
}



# Biting function: ODE model for malaria transmission with control methods (LLINs and Systemic insecticides) -----

biting <- function(t, y, parms) {
  with(as.list(c(parms, y)), {
    
   if(D < 0){   # in case efficacy falls below 0
      y[["D"]] <- 0
      D <- 0
    }
    
    # EQs and ODEs
    V <- S_m + E_m + I_m             # total mosquitoes
    H <- S_h + E_h + I_h + R_h + A_h # total humans
    
    b_h <-  a * b * (1 - c_D * D)    #  proportion of bites resulting in human infections 
    b_m <-  a * c * (1 - c_D * D)     #  proportion of bites resulting in mosquito infections   
    
    K = K_0 * H * (1 + delta * cos(2 * pi * (t/365 + omega)))   # sinusoidal seasonal pattern, need to convert t is in years
    xi_m = xi_m_0 / ( 1 - delta * cos(2 * pi * (t/365 + omega)))
    mu_m = mu_m_0 / ( 1 + delta * cos(2 * pi * (t/365 + omega)))
    
    k_d <- log(2)/hl_D      # degradation rate of repellent efficacy [1/day]
    
    
    # humans
    dS_h <- r * (E_h + I_h + A_h) + q2 * R_h - b_h * m * S_h * I_m        # susceptible humans
    dE_h <- b_h * m * S_h * I_m  - (xi_h + r) * E_h                       # exposed humans
    dI_h <- xi_h * E_h - (q1 + r) * I_h                                   # infected humans
    dR_h <- q1 * (I_h + A_h) - (theta * b_h * m * I_m + q2) * R_h          # recovered/immune humans
    dA_h <- theta * b_h * m * R_h * I_m - (q1 + r) * A_h                  # asymptomatically infected humans
    
    # mosquitoes
    dO <- beta_m * V - d_o * O - mu_o * (1 + (O + L)/ K) * O                    # eggs and early larval instars
    dL <- d_o * O - d_l * L - mu_l * (1 + gam * (O + L)/ K) * L                 # late larval instars
    dP <- d_l * L - d_p * P - mu_p * P                                           # pupae
    dS_m <- 1/2 * d_p * P - b_m * (I_h + sigma * A_h) * S_m - mu_m * S_m         # susceptible mosquitoes
    dE_m <- b_m * (I_h + sigma * A_h) * S_m - (xi_m + mu_m) * E_m                # exposed mosquitoes  
    dI_m <- xi_m * E_m - mu_m * I_m                                              # infected mosquitoes
    
    # intervention
    dD <- -k_d * D
    
    # output
    res <- c(dS_h, dE_h, dI_h, dR_h, dA_h, dO, dL, dP, dS_m, dE_m, dI_m, dD)
    list(res)
  })
}


# define parameters and  ranges-----

# Parameters
 P <- c(
   beta_m = 21.19, # 1. number of eggs a mosquito lays per day [1/day]
   d_o = 0.15, # 2. rate at which early larval instars mature into late larval instars [1/day]
   d_l = 0.27, # 3. rate at which late larval instars mature into pupae [1/day]
   d_p = 1.56, # 4. rate at which pupae matures into mosquitoes [1/day]
   a = 0.2,    # 5. biting frequency [1/day]
   b = 0.5,    #6.  proportion of bites that produce infection in humans [unitless]
   c = 0.5,    # 7. proportion of bites that produce infection in mosquitoes [unitless]
   sigma = 0.25, # 8. adjustment factor for asymptomatic infection transmissibility to vector [unitless]
   xi_m_0 = 1/10, # 9. average rate of P. falciparum maturation in mosquito [1/day]
   xi_h = 1/21, # 10. rate of P. falciparum maturation in human [1/day]
   q1 = 1/200, # 11. rate of immunity acquisition [1/day]
   q2 = 1/1000, # 12. rate of immunity loss [1/day]
   theta = 0.5, # 13. level of reduced susceptibility to secondary infection [unitless]
   r = 0.01, # 14. rate of recovery [1/day]
   mu_o = 0.034, # 15. death rate of early larval instars [1/day]
   mu_l = 0.035, # 16. death rate of late larval instars [1/day]
   K_0 = 3,#100, # 17. average carrying capacity of environment assuming 50 mm water over the past 4.5 days and 150 mosquitoes captured
   gam = 13.25, # 18. adjustment factor to correct for different density dependence of late vs early instars
   mu_p = 0.25, # 19. death rate of pupae [1/day]
   mu_m_0 = 0.12, # 20. average death rate of mosquitoes [1/day]
   m = 10, # 21. ratio of mosquitoes to humans
   delta = 0.15, # 22. amplitude of seasonal forcing
   omega = 0, # 23. phase of seasonal forcing
   t_ss = 3*365,   # 24. time allowed to reach steady state [days]
   treatment_pd = 1 * 365   # 25. treatment period [days]
   )

# define initial conditions
y0 <- c(
  S_h = 0.96,
  E_h = 0.01,
  I_h = 0.01,
  R_h = 0.01,
  A_h = 0.01,
  O = 0.01,
  L = 0.01,
  P = 0.01,
  S_m = 0.01,
  E_m = 0.01,
  I_m = 0.01,
  D = 0
)

# define ranges for period, coverage, half life, etc
period_D_range <- c(1/6, 1/4, 1/3, 1/2, 1, 2, 3, 4, 5, 6, 7)      # frequency of repellent application [1/day] (every 4 hours, 6 hours...)
c_D_range <- seq(0.1, 1, .1)                                      # repellent coverage 
efficacy_range <- seq(0.25, 1, 0.25)                              # Instead of giving a specific concentration for a dose, I've coded things in terms of repellent efficacy (i.e. a dose could repel 25% of mosquitoes vs 100%)   *** UPDATE THIS WITH LAB DATA ****
halflife_D_range <- c( 6/24, 12/24, 1, 2)                         # half life of repellent **** UPDATE THIS WITH LAB DATA ****

# Simulate time courses, capture malaria prevalence at end of treatment window, and plot contour maps of
# dosing frequencies needed for each combination.
w <- size(c_D_range, 2)
x <- size(period_D_range, 2)
y <- size(efficacy_range, 2)
z <- size(halflife_D_range, 2)

someData <- rep(0, w*x*y*z)
# periodUse <- array(someData, c( x, y, z))
cases.avoided <- array(someData, c( w, x, y, z))

# calculate prevalence of malaria with repellent ----

# calculate number of infected cases without any control methods (baseline)
base <- baseline(y0, P[1:25])
time_base <- base[[1]]
infected_base <- base[[2]]
base.df <- cbind.data.frame(time_base, infected_base)
base.df$t.months <- (base.df$time - P["t_ss"])/30
base.df <- base.df[base.df$t.months >= 0 & base.df$t.months <= 12, ][-1]
tot_I_base <- AUC(base.df$t.months, base.df$infected_base)# add up proportion of infected people for time window (1 year)

# calculate number of infected cases for different repellent scenarios
# (changing the coverage, period between reappliance, efficacy, and half-life)
for (i in 1:size(c_D_range, 2)){
  for (j in 1:size(period_D_range, 2)){
    for (k in 1:size(efficacy_range, 2)){
      for (l in 1:size(halflife_D_range, 2)){
        # i = 1
        # j = 3
        # k = 1
        # l = 1
        
        print(c(i,j,k,l))  # indicator for which loop you're on
        
        P = c(
          P[1:25],
          c_D = c_D_range[i],
          period_D = period_D_range[j],
          dose_eff = efficacy_range[k], 
          hl_D = halflife_D_range[l]
        )
        
        # calculate number of infected cases with bednets and systemic insecticides 
        # compute the relative prevalence of malaria cases of avoided with addition of SI
        repellent <- dosing(y0, P)
        time_rep <- repellent[[1]]
        infected_rep <- repellent[[2]]
        drug <- repellent[[3]]
        
        repellent.df <- cbind.data.frame(time_rep, infected_rep)
        repellent.df$t.months <- (repellent.df$time_rep - P["t_ss"])/30
        repellent.df <- repellent.df[repellent.df$t.months >= 0 & repellent.df$t.months <= 12, ][-1]
        tot_I_repellent <- AUC(repellent.df$t.months, repellent.df$infected_rep)#sum(repellent.df$I_h)
        print(tot_I_repellent)
        
        
        names(repellent.df) <- c("I_h","time")
        names(base.df) <- c("I_h","time")
        newData <- melt(list(repellent.df = repellent.df, base.df = base.df),
                        id.vars = "time")
        
        ## plot to compare temporal dynamics of no treatment vs treatment   
        ## uncomment the plot below if you want to look at temporal dynamics within a certain loop. It will slow things down, so don't leave uncommented for running whole thing.
        
        # ggplot(newData, aes(time, value, color = L1, fill = L1))+
        #   geom_line(size = 2)+
        #   geom_area(aes(fill = L1, group = L1),
        #             alpha = 0.5, position = 'identity')+
        #   lims(x = c(0,12), y = c(0,1))+
        #   scale_color_manual(name = "Treatment", 
        #                      limits = c("base.df", "repellent.df"),
        #                      labels = c("None", "Repellent"),
        #                      values=c("black","red"))+
        #   scale_fill_manual(name = "Treatment", 
        #                      limits = c("base.df", "repellent.df"),
        #                      labels = c("None", "Repellent"),
        #                      values=c("black","red"))+
        #   theme(aspect.ratio = 1) +
        #   theme(
        #     panel.grid.major = element_blank(),
        #     panel.grid.minor = element_blank(),
        #     panel.background = element_blank(),
        #     axis.line = element_line(colour = "black")
        #   ) +
        #   xlab("Time (months)") +
        #   theme(
        #     axis.text.x = element_text(colour = 'black', size = 14),
        #     axis.title.x = element_text(
        #       size = 16,
        #       hjust = 0.5,
        #       vjust = 0.5
        #     )
        #   ) +
        #   ylab("Proportion of humans infected") +
        #   theme(
        #     axis.text.y = element_text(colour = 'black', size = 14),
        #     axis.title.y = element_text(
        #       size = 16,
        #       hjust = 0.5,
        #       vjust = 0.2
        #     )
        #   ) +
        #   coord_cartesian(xlim = c(0, 12), ylim = c(0, 1), expand=FALSE)
        
        cases.avoided[i,j,k,l] <- (tot_I_base - tot_I_repellent)/tot_I_base * 100
        
        print((c(c_D_range[i],
                 period_D_range[j],
                 efficacy_range[k], 
                 halflife_D_range[l], 
                 cases.avoided[i,j,k,l])))
        
      }
    }
  }
}

# turn array of cases avoided for different scenarios into long-form dataframe
cases.avoided <- provideDimnames(cases.avoided, sep= "_", base = list (as.character(c_D_range), as.character(period_D_range), as.character(efficacy_range), as.character(halflife_D_range)))
cases.avoided.long <- melt(cases.avoided, value.name = "cases avoided")
colnames(cases.avoided.long) <- c('coverage', 'period', 'efficacy', 'half.life', 'cases.avoided')
# define labels for facet_grid
half.life.labs <- c('Half-life = 6 h', '12 h', '1 d', '2 d')
names(half.life.labs) <- c(0.25, 0.5, 1, 2)

efficacy.labs <- c('Max efficacy = 25%', '50%', '75%', '100%')
names(efficacy.labs) <- c(0.25, 0.5, 0.75, 1)
# plot cases avoided as a function of period, coverage, half-life, and drug efficacy
ggplot(cases.avoided.long) + 
  aes(x = coverage, y = as.factor(period), z = cases.avoided, fill = cases.avoided) + 
  geom_tile() + 
  scale_fill_gradient(low = "blue", high = "red",
                      limits = c(0, 80)) + 
  theme_bw()+
  facet_grid(efficacy ~ half.life,
             labeller = labeller(half.life = half.life.labs,
                                 efficacy = efficacy.labs))+
  theme(aspect.ratio = 1,
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(colour = 'black', size = 10),
        axis.title.x = element_text(
          size = 12,
          hjust = 0.5,
          vjust = 0.5),
        axis.text.y = element_text(colour = 'black', size = 10),
        axis.title.y = element_text(
          size = 12,
          hjust = 0.5,
          vjust = 0.2)
  ) +
  xlab("Coverage (proportion)")+ 
  scale_x_continuous(breaks=c(0,0.5,1))+
  ylab("Time between applications")+
  scale_y_discrete(breaks = c(4/24, 6/24, 8/24, 12/24, 1, 2, 3, 4, 5, 6, 7),
                   labels = c('4 h', '6 h', '8 h', '12 h', '1 d', '2 d', '3 d', '4 d', '5 d', '6 d', '7 d'))+
  labs(fill = "Cases avoided (%)")

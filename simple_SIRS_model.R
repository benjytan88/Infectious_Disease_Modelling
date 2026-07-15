# Very simplified SIRS model using syphilis as an example

# -----------------------------------------------------------------------------
# IMPORTANT ASSUMPTIONS!!!
# (1) Single-stage disease without treatment
# (2) Closed population ignoring population dynamics
#     -- birth/death/new susceptibles are not considered i.e. fixed population size 
# -----------------------------------------------------------------------------

# Load libraries
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(deSolve))

# Set initial population state
init_state <- c(S = 9999,
                I = 1,
                R = 0)

# ODE parameters
# -- Use average data from Garnett et al. (1997) paper (10.1097/00007435-199704000-00002)
# -- Also referred Feldman & Mishra (2019) paper (10.1016/j.idm.2019.09.002)

# Beta = Transmission rate
# -- Transmission probability per partner is around 60%

# Gamma = Recovery rate
# -- Mean infectious period of syphillis is 154 days

# Omega = Reinfection rate among recovered
# -- Mean duration to relapse is 108
params <- c(beta = 0.627,
            gamma = 0.0065,  # 1/154
            omega = 0.0092)   # 1/108

# Time length (days)
time <- seq(from = 1,
            to = 365,
            by = 1)

# SIRS model function
SIRS <- function(time, state, parameters){
  with(as.list(c(state, parameters)),{
    N = S + I + R
    dS = (-beta * S * I / N) + (omega * R)
    dI = (beta * S * I / N) - (gamma * I)
    dR = (gamma * I) - (omega * R)
    return(list(c(dS, dI, dR)))
  })
}

# Solve the SIR model
output_SIR <- as.data.frame(ode(y = init_state, 
                                times = time,
                                func = SIRS,
                                parms = params))

# Convert to long format for plotting
output_SIR_long <- reshape2::melt(output_SIR, id = "time")
colnames(output_SIR_long) <- c("time", "compartment", "population")

# Plot populations over time
ggplot(data = output_SIR_long,          
       aes(x = time, y = population/10000, colour = compartment, group = compartment)) +  
  geom_line(size = 1.2) +
  labs(x = "Time (Days)", y = "Proportion of the population", title = "Simple SIRS Model") +
  scale_color_manual(values = c("steelblue", "firebrick", "darkgreen")) +
  theme_classic() +
  theme(legend.box.background = element_rect(fill = "transparent"),
        legend.background = element_rect(fill = "transparent"),
        axis.title = element_text(size = 14))


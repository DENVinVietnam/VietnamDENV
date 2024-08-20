##############################################################################################################
## This file contains:
## SEIRS model
## - explore seasonality
## - generate state names
## - define one-step transitions
## - extract coefficients from RO1 and apply to adjust for beta
## - SEIRS model function with initial conditions, parameters and simulation time points
##
## Visualization
## - Plot 1: Total population
## - Plot 2: Yearly Incidence
## - Plot 3: Number of people in specif compartment over the simulation period
## - Plot 4: Aggregated infection trends by serotype over years
## - Plot 5: Aggregated SEIR plot over years
## - Plot 6: Infection by Previous Infection History
##############################################################################################################

library(pacman)
p_load(dplyr, ggplot2, deSolve, stringr)

##################################
## Section 1: Model Development###
##################################


##############################################################################################################
## Explore and visualize the seasonality, rainy season occurs in the later half of the year
##############################################################################################################
t <- seq(1, 365, 1)
seas <- 1 + 0.2 * cos(2 * pi * (t / 365 + 0.3))
plot(t, seas, type = "l", xaxt = "n", xlab = "Month", ylab = "Seasonality Factor", 
     main = "Seasonality Over the Year")
months <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
axis(1, at = seq(1, 365, by = 365/12), labels = months)


##############################################################################################################
## Generate all the possible states for Susceptible, Exposed, Infected and Recovered
##############################################################################################################
PRNT50_states <- c("N", "L", "H") # antibody levels
# Each serotype can has "N", "L", "H" antibody levels
state_combinations <- expand.grid(sero1 = PRNT50_states, sero2 = PRNT50_states, 
                                  sero3 = PRNT50_states, sero4 = PRNT50_states)

## Susceptible states
## 81 susceptible states in total e.g S_NLLH, S_LNNL
## Note that susceptible states without "N" in the index do not exist, as people gain immunities against all
## serotypes, so he/she is no longer susceptible. For convenience, we include them but assign no transitions
## from them and no initial people in these states.
susceptible_states <- transform(state_combinations,
                 State_Name = apply(state_combinations, 1, function(x) paste0("S_", paste(x, collapse = ""))))

## Exposed States
# The numerical identifier matches the serotype exposed or infected
# It should be on the corresponding index position
# i.e exposed by serotype 1 -> "E_1NNN" "E_1LNN" "E_1HNN"... 1 can only be on the first position
# Initialize a vector to store all exposed states
exposed_states <- vector("list", length = 4)
# Loop through each serotype
for (serotype in 1:4) {
  # Generate all combinations of antibody levels for other serotypes
  other_states <- expand.grid(rep(list(PRNT50_states), 4))
  
  # Form the exposed state names, e.g E_1NLH, E_HLL4
  exposed_states[[serotype]] <- apply(other_states, 1, function(states) {
    # Only those individuals never infected by the serotype before can be infected (life-long immunity)
    if(states[serotype] == "N"){
      # Set the infected serotype with its numeric identifier, 
      # Note that the numeric identifier matches the index position
      states[serotype] <- as.character(serotype)  
      paste0("E_", paste(states, collapse = ""))
    }
  })
}

exposed_states <- unlist(exposed_states)

## Infectious States, same logic as exposed states
infectious_states <- vector("list", length = 4)
for (serotype in 1:4) {
  other_states <- expand.grid(rep(list(PRNT50_states), 4))

  infectious_states[[serotype]] <- apply(other_states, 1, function(states) {
    if(states[serotype] == "N"){
    states[serotype] <- as.character(serotype)  
    paste0("I_", paste(states, collapse = ""))
    }
  })
}

infectious_states <- unlist(infectious_states)

## Recovered states, same logic as susceptible states
## Note that R_NNNN state does not exist. For convenience, we include it but assign no initial people in it,
## and no transition to it.
recovered_states <- transform(state_combinations,
              State_Name = apply(state_combinations, 1, function(x) paste0("R_", paste(x, collapse = ""))))

# All the states
all_SEIR_states <- c(susceptible_states[, "State_Name"], exposed_states, 
                     infectious_states, recovered_states[, "State_Name"])


##############################################################################################################
## Function to generate all the one step movement from any states
## Input (character): state - the state to transit from
## Transition Rules:
## - S to E: Transition occurs only if the individual is exposed by the serotype which he/she has not been 
##           previously infected (without pre-existing antibody level)
## - E to I: Transition occurs if the individual become infectious with the exposed serotype
## - I to R: Transition occurs if the individual recovered with either a high or low antibody level against 
##           the infected serotype, and become immune to all the serotypes (TCI)
## - R to S: Transition occurs if the individual lose TCI and goes back to the corresponding susceptible states 
##           with the updated antibody level
## Output (character): next_states: all possible one step states to transit to
##############################################################################################################
generate_next_states_SEIR <- function(state) {
  parts <- unlist(strsplit(sub("^(.)_", "", state), ""))
  state_type <- substr(state, 1, 1)
  next_states <- c()
  
  if (state_type == "S") {
    # From Susceptible to Exposed transitions
    for (i in 1:length(parts)) {
      if (parts[i] == "N") {  # Transition only if no previous infection
        exposed_parts <- parts
        exposed_parts[i] <- as.character(i)  
        next_states <- c(next_states, paste0("E_", paste0(exposed_parts, collapse = "")))
      }
    }
  } else if (state_type == "E") {
    # From Exposed to Infected transitions
    for (i in 1:length(parts)) {
      if (parts[i] %in% c("1", "2", "3", "4")) {
        infected_parts <- parts
        infected_parts[i] <- as.character(i)  # Mark as infected by this serotype
        next_states <- c(next_states, paste0("I_", paste0(infected_parts, collapse = "")))
      }
    }
  } else if (state_type == "I") {
    # From Infected to Recovered transitions
    serotype_infected <- which(parts %in% c("1", "2", "3", "4"))
    for (serotype in serotype_infected) {
      high_parts <- parts
      low_parts <- parts
      high_parts[serotype] <- "H"
      low_parts[serotype] <- "L"
      next_states <- c(next_states, paste0("R_", paste0(high_parts, collapse = "")))
      next_states <- c(next_states, paste0("R_", paste0(low_parts, collapse = "")))
    }
  } else if (state_type == "R") {
    any_high <- any(parts == "H")
    any_none <- any(parts == "N")
    
    # All serotypes are infected, so no transition back to susceptible states
    if (!any_none) {
      next_states <- NULL
    } else{ # Not infected by all the serotypes
      next_states <- c(next_states, paste0("S_", paste0(parts, collapse = "")))
       }
  }
  return(unique(next_states))
}

## Checking
print(generate_next_states_SEIR("R_HNNH"))

# Initialize a list with all states having empty vectors as default for storing transitions
state_transitions <- setNames(lapply(seq_along(all_SEIR_states), 
                                     function(x) vector("list", 0)), all_SEIR_states)

# Apply the function to fill in the transition list for each state, 
# If no transitions to the state, remain the default empty vector
for (state in all_SEIR_states) {
  transitions <- generate_next_states_SEIR(state)
  if (length(transitions) > 0) {
    state_transitions[[state]] <- transitions
  }
}



##############################################################################################################
## Adjust beta base on the results from multinomial regression
##############################################################################################################
# Extract the coefficients from the model output in S2 Code
coef_matrix <- coef(multi_model1_low_high)

# Fill in the coefficients list, re-assign the values to D3 infection due to the small sample size
# Scale down the estimates for high PRNT DENV4 against D2 infection due to the large variability
coefficients <- list(
  D1 = c(Intercept = coef_matrix["D1", "(Intercept)"], 
         Past_PRNT50_D1high = coef_matrix["D1", "Past_PRNT50_D1high"], 
         Past_PRNT50_D1low = coef_matrix["D1", "Past_PRNT50_D1low"],
         Past_PRNT50_D2high = coef_matrix["D1", "Past_PRNT50_D2high"], 
         Past_PRNT50_D2low = coef_matrix["D1", "Past_PRNT50_D2low"],
         Past_PRNT50_D3high = coef_matrix["D1", "Past_PRNT50_D3high"], 
         Past_PRNT50_D3low = coef_matrix["D1", "Past_PRNT50_D3low"],
         Past_PRNT50_D4high = coef_matrix["D1", "Past_PRNT50_D4high"], 
         Past_PRNT50_D4low = coef_matrix["D1", "Past_PRNT50_D4low"]),
  
  D2 = c(Intercept = coef_matrix["D2", "(Intercept)"], 
         Past_PRNT50_D1high = coef_matrix["D2", "Past_PRNT50_D1high"], 
         Past_PRNT50_D1low = coef_matrix["D2", "Past_PRNT50_D1low"],
         Past_PRNT50_D2high = coef_matrix["D2", "Past_PRNT50_D2high"], 
         Past_PRNT50_D2low = coef_matrix["D2", "Past_PRNT50_D2low"],
         Past_PRNT50_D3high = coef_matrix["D2", "Past_PRNT50_D3high"], 
         Past_PRNT50_D3low = coef_matrix["D2", "Past_PRNT50_D3low"],
         Past_PRNT50_D4high = 5, 
         Past_PRNT50_D4low = coef_matrix["D2", "Past_PRNT50_D4low"]),
  
  D3 = c(Intercept = 0, 
         Past_PRNT50_D1high = 0, 
         Past_PRNT50_D1low = 0,
         Past_PRNT50_D2high = 0, 
         Past_PRNT50_D2low = 0,
         Past_PRNT50_D3high = 0, 
         Past_PRNT50_D3low = 0,
         Past_PRNT50_D4high = 0, 
         Past_PRNT50_D4low = 0),
  
  D4 = c(Intercept = coef_matrix["D4", "(Intercept)"], 
         Past_PRNT50_D1high = coef_matrix["D4", "Past_PRNT50_D1high"], 
         Past_PRNT50_D1low = coef_matrix["D4", "Past_PRNT50_D1low"],
         Past_PRNT50_D2high = coef_matrix["D4", "Past_PRNT50_D2high"], 
         Past_PRNT50_D2low = coef_matrix["D4", "Past_PRNT50_D2low"],
         Past_PRNT50_D3high = coef_matrix["D4", "Past_PRNT50_D3high"], 
         Past_PRNT50_D3low = coef_matrix["D4", "Past_PRNT50_D3low"],
         Past_PRNT50_D4high = coef_matrix["D4", "Past_PRNT50_D4high"], 
         Past_PRNT50_D4low = coef_matrix["D4", "Past_PRNT50_D4low"])
)


##############################################################################################################
## adjust_beta function: Adjusts the base transmission rate (beta) based on immune levels
## Inputs:
##  - state (character): The current state of the individual
##  - serotype (integer): The numeric identifier of the exposed serotype 
##  - base_betas (vector): The base transmission rates for each serotype
##  - coefficients (list): The coefficients for high and low antibody levels against each serotype
## Output: beta (numeric) - The adjusted transmission rate based on the antibody levels
##############################################################################################################
adjust_beta <- function(state, serotype, base_betas, coefficients) {
  antibody_levels <- strsplit(substr(state, 3, nchar(state)), "")[[1]]
  # Initialize beta with the base transmission rate for the given serotype
  beta <- base_betas[serotype]
  # Extract the coefficients for the given serotype
  coef <- coefficients[[paste0("D", serotype)]]
  
  # Adjust beta
  for (i in seq_along(antibody_levels)) {
    if (antibody_levels[i] == "H") {
      beta <- beta * exp(coef[paste0("Past_PRNT50_D", i, "high")])
    } else if (antibody_levels[i] == "L") {
      beta <- beta * exp(coef[paste0("Past_PRNT50_D", i, "low")])
    }
    # "N doesn't change the beta value
  }
  
  return(beta)
}

##############################################################################################################
## This is the function of the SEIRS compartment model with seasonality, TCI, birth and death rate
## Inputs:
##  - t (numeric): The current time point
##  - y (vector): The current state vector of the number of individuals in each compartment
##  - params (list): A list of parameters used in the differential equation as defined below
## Output:
##  - dy (vector): The rate of change for each compartment
##############################################################################################################
DENV_SEIRS_bd <- function(t, y, params) {
  # Initialize rate of change for each state with default value 0
  dy <- numeric(length(y))
  
  # Names of each state
  names <- names(y)
  names(dy) <- names
  # Total population, excluding CInc
  P <- sum(y[-379])
  
  # Birth rate is applied to the whole population and goes into the fully naive cohort
  ## mu*P (dS_NNNN)
  dy["S_NNNN"] <- dy["S_NNNN"] + params$mu * P
  
  # Iterate through each state to apply transition rules to find their possible next states 
  for (state in names(y)[-379]) {
    current_index <- which(names(y) == state)
    
    # Get transitions for the current state
    transitions <- state_transitions[[state]]
    
    # Update each possible transition from the current state
    if (!is.null(transitions)) {
      for (next_state in transitions) {
        next_index <- which(names(y) == next_state)
        
        # Determine the type of transition based on state names
        # Susceptible to Exposed transition
        if (grepl("^S_", state) && grepl("^E_", next_state)) {
          # Extracting exposed serotype numeric identifier
          serotype_exposed <- as.numeric(gsub("[^0-9]", "", sub("^[^_]*_", "", next_state)))
          
          # Calculate total infected individuals for this serotype since all the people 
          # infected by that serotype contributes to the infectiousness regardless of 
          # their pre-existing antibody levels
          serotype_infected_states <- grep(paste0("I_"), names(y), value = TRUE)
          # Find those states that are infected with the extracted exposed serotype
          serotype_infected_states <- serotype_infected_states[sapply(serotype_infected_states,
           function(x) any(unlist(strsplit(substr(x, 3, nchar(x)), "")) == as.character(serotype_exposed)))]
          total_infected <- sum(y[serotype_infected_states])
          
          # Apply adjust_beta function to get the dynamic beta value
          beta <- adjust_beta(state, serotype_exposed, params$betas, coefficients)
          
          ## Add seasonality
          seas <- 1 + params$amp * cos(2 * pi * (t / 365 - params$phi))
          lam <- seas * beta * total_infected / P
          
          # Calculate transitions from S to E for cumulative incidence
          S_to_E_transition_rate <- lam * y[current_index]
          dy["CInc"] <- dy["CInc"] + S_to_E_transition_rate
          
          ## -lam*S (dS)
          dy[current_index] <- dy[current_index] - lam * y[current_index] 
          ## lam*S (dE)
          dy[next_index] <- dy[next_index] + lam * y[current_index] 
          
          # Exposed to Infected transition
        } else if (grepl("^E_", state) && grepl("^I_", next_state)) {
          # Extract infected serotype number
          serotype_exposed <- as.numeric(gsub("[^0-9]", "", sub("^[^_]*_", "", next_state))) 
          
          ## -alpha*E (dE)
          dy[current_index] <- dy[current_index] - params$alpha * y[current_index]
          # alpha*E (dI)
          dy[next_index] <- dy[next_index] + params$alpha * y[current_index] 
          
          # Infected to Recovered transition
          # Assign a proportion of p to be recovered with high antibody level, (1-p) recovered with low level
        } else if (grepl("^I_", state) && grepl("^R_", next_state)) {
          # Extracting infected serotype number
          serotype_infected <- as.numeric(gsub("[^0-9]", "", sub("^[^_]*_", "", state))) 
          
          # "+2" means we are looking at the string starting from the 3rd character in the state name 
          # the first 2 characters are the state initial and "_"
          # Recovered with a high antibody level
          if (substr(next_state, serotype_infected + 2, serotype_infected + 2) == "H") {
            # -p*gamma*I (dI)
            dy[current_index] <- dy[current_index] - params$p * params$gamma * y[current_index] 
            # p*gamma*I (dR_H)
            dy[next_index] <- dy[next_index] + params$p * params$gamma * y[current_index] 
            
            # Recovered with a low antibody level
          } else if (substr(next_state, serotype_infected + 2, serotype_infected + 2) == "L") {
            # -(1-p)*gamma*I (dI)
            dy[current_index] <- dy[current_index] - (1 - params$p) * params$gamma * y[current_index] 
            # (1-p)*gamma*I (dR_L)
            dy[next_index] <- dy[next_index] + (1 - params$p) * params$gamma * y[current_index]
          }
          # Recovered to Susceptible transition
        } else if (grepl("^R_", state) && grepl("^S_", next_state)) {
          # -sigma*R (dR)
          dy[current_index] <- dy[current_index] - params$sigma * y[current_index] 
          # sigma*R (dS)
          dy[next_index] <- dy[next_index] + params$sigma * y[current_index]
        } 
      }
    }
    
    # Apply death rate to all compartments
    dy[current_index] <- dy[current_index] - params$mu * y[current_index]
  }
  
  return(list(dy))
}



##############################################################################################################
## Initial condition
##############################################################################################################
# Initialize all states with 0 people
SEIRS_initial_conditions <- rep(0, length(all_SEIR_states) + 1)
names(SEIRS_initial_conditions) <- c(all_SEIR_states, "CInc")

# Total population in Vietnam, 
# assuming half of the total population at the risk of dengue  (Jones et al., 2020)
total_population_target <- 98000000 * 0.5 

# Set 10 people in each recovery state, except for "R_NNNN"
recovery_states <- grep("^R_", all_SEIR_states, value = TRUE)
SEIRS_initial_conditions[recovery_states] <- 10000
# No transitions go into this state, not exist this state
SEIRS_initial_conditions["R_NNNN"] <- 0

# Prevalence of each serotype estimated from the data 
prevalence_D1 <- 0.0320
prevalence_D2 <- 0.0240
prevalence_D3 <- 0.00267
prevalence_D4 <- 0.0307

# The number of individual in E and I compartments depends on the prevalence of each serotype from the data
# 27 states for each infectious or exposed serotype
for (state in all_SEIR_states) {
  if (grepl("^E_", state) || grepl("^I_", state)) {
      if (grepl("1", state)) {
        SEIRS_initial_conditions[[state]] <- (prevalence_D1 * total_population_target)/27/2
      } else if (grepl("2", state)) {
        SEIRS_initial_conditions[[state]] <- (prevalence_D2 * total_population_target)/27/2
      } else if (grepl("3", state)){
        SEIRS_initial_conditions[[state]] <- (prevalence_D3 * total_population_target)/27/2
      } else if (grepl("4", state)){
        SEIRS_initial_conditions[[state]] <- (prevalence_D4 * total_population_target)/27/2
      }
  }
}


# Calculate the total number of people already assigned
assigned_population <- sum(SEIRS_initial_conditions)
# Calculate the remaining population for susceptible states
remaining_population <- total_population_target - assigned_population

# Identify susceptible states
susceptible_states <- grep("^S_", all_SEIR_states, value = TRUE)

# Allocate 80% of the remaining population to S_NNNN
SEIRS_initial_conditions["S_NNNN"] <- 0.8 * remaining_population
remaining_population <- 0.2 * remaining_population

# Distribute the remaining population among the other susceptible states based on the number of 'N's,
susceptible_states <- setdiff(susceptible_states, "S_NNNN")
num_Ns <- sapply(susceptible_states, function(state) sum(strsplit(state, "")[[1]] == "N"))
distribution_weights <- num_Ns / sum(num_Ns) 
# S states with measurable antibody level against all serotype do not exist, since those people
# are immune to all the serotypes due to the life long immunity assumption against homologous serotypes
# so 0 people assigned to these states (e.g S_HLLH)
SEIRS_initial_conditions[susceptible_states] <- distribution_weights * remaining_population

# Ensure integer values for the initial conditions
SEIRS_initial_conditions <- round(SEIRS_initial_conditions)
SEIRS_initial_conditions
##############################################################################################################
## Simulation Period
##############################################################################################################
# Time points
SEIRS_bd_times <- seq(0, 365*50, by = 1) 

##############################################################################################################
## Model Parameters
##############################################################################################################
SEIRS_params_bd2 <- list(
  mu = (1/(70*52*7)),
  #betas = c(0.2, 0.2, 0.2, 0.2), # homogeneous base transmission rate
  betas = c(0.35, 0.22, 0.1, 0.2), # heterogeneous base transmission rate
  p = 1/10,           # proportion of people gain high antibody level after recovery from I
  alpha = (1/7),      # incubation period, rate of onset of infectiousness 
  gamma = (1/10),     # recovery rate
  sigma = (1/150),    # rate of loss of immunity, R to S
  amp = 0.2,          # amplitude of the seasonality (magnitude of variation in the transmission rate)
  phi = -0.3          # phase shift of seasonality (peak time)
)

# Run the model
data1_40_160_heter<- ode(y = SEIRS_initial_conditions, times = SEIRS_bd_times, 
                        func = DENV_SEIRS_bd, parms = SEIRS_params_bd2)

###################################
## Section 2: Visualization########
###################################
# Convert the output to a data frame, change the output name accordingly 
SEIRS_bd_results_df <- as.data.frame(data1_40_160_heter)

##############################################################################################################
## Plot 1: Total population
##############################################################################################################
# Remove the last column (CInc) when calculate the population
plot(round(rowSums(SEIRS_bd_results_df[,-c(1,380)])), xlab = "Time(years)", ylab = "Total population")


##############################################################################################################
## Plot 2: Yearly Incidence
##############################################################################################################
# Calculate the incidence from cummulative incidence
data_inc <- as_tibble(SEIRS_bd_results_df) %>%
  mutate(Inc = c(rep(0,364), diff(CInc, lag = 364))) %>%
  pivot_longer(names_to = "variable",cols = !1)

# Select only the incidence at the end of the year
data_CInc_yearly <- data_inc %>%
  filter(variable %in% c("Inc")) %>%
  mutate(year = time %/% 364) %>%
  filter (time %% 364 == 0 & time != 0) %>%
  group_by(year) %>%
  summarize(sum_value = sum(value))

# Exclude the incidence in the first year
plotdata_CInc_yearly <- data_CInc_yearly %>%
  filter(!(year %in% c(1)))


ggplot() +
  geom_line(data = plotdata_CInc_yearly, aes(x = year, y = sum_value), color = "blue", size = 1) +
  labs(title = "Yearly Dengue Incidence and Actual Dengue Cases",
       x = "Year",
       y = "Number of Cases",
       color = "Legend") +
  #ylim(c(0, 1*10^7)) +
  scale_color_manual(values = c("Modelled Incidence" = "blue", "Actual Cases" = "red")) +
  theme(legend.position = "none")


##############################################################################################################
## Plot 3: Visualize the trend of specific compartments, choose the interested compartments to geenrate plot
##############################################################################################################
# Specify the compartments of interest
ggplot(SEIRS_bd_results_df, aes(x = time)) +
  geom_line(aes(y = S_NNNN, colour = "S_NNNN")) +
  geom_line(aes(y = I_1LNN, colour = "I_1LNN")) +
  geom_line(aes(y = I_N2NN, colour = "I_N2NN")) +
  geom_line(aes(y = I_NN3N, colour = "I_LN3N")) +
  geom_line(aes(y = I_NNN4, colour = "I_NNN4")) +
  geom_line(aes(y = R_HHNN, colour = "R_HHNN")) +
  labs(title = "Dengue Infection Dynamics", y = "Number of Individuals", x = "Time (days)") +
  scale_color_manual(values = c("S_NNNN" = "blue", "I_1LNN" = "red", "I_N2NN" = "purple",
                                "I_LN3N" = "cornflowerblue", "I_NNN4" = "black", "R_HHNN" = "green")) 

##############################################################################################################
## Plot 4: Aggregate Total Infections for each serotype
##############################################################################################################
# Extract all the states and sum the total infections for each serotype at each time point
total_infected_S1 <- rowSums(SEIRS_bd_results_df[, grep("^I_1", names(SEIRS_bd_results_df))])
total_infected_S2 <- rowSums(SEIRS_bd_results_df[, grep("^I_.2", names(SEIRS_bd_results_df))])
total_infected_S3 <- rowSums(SEIRS_bd_results_df[, grep("^I_..3", names(SEIRS_bd_results_df))])
total_infected_S4 <- rowSums(SEIRS_bd_results_df[, grep("^I_...4", names(SEIRS_bd_results_df))])

# Add these infections to the dataframe for plotting
SEIRS_bd_results_df$total_infected_S1 <- total_infected_S1
SEIRS_bd_results_df$total_infected_S2 <- total_infected_S2
SEIRS_bd_results_df$total_infected_S3 <- total_infected_S3
SEIRS_bd_results_df$total_infected_S4 <- total_infected_S4

# Plot the aggregated total infections by serotype
ggplot(SEIRS_bd_results_df, aes(x = time)) +
  geom_line(aes(y = total_infected_S1, colour = "Serotype 1")) +
  geom_line(aes(y = total_infected_S2, colour = "Serotype 2")) +
  geom_line(aes(y = total_infected_S3, colour = "Serotype 3")) +
  geom_line(aes(y = total_infected_S4, colour = "Serotype 4")) +
  labs(title = "Total Infected Individuals by Serotype", 
       y = "Number of Infected Individuals", 
       x = "Time (years)") +
  ylim(c(0, 5*10^5)) +
  theme_minimal() +
  scale_color_manual(values = c("Serotype 1" = "red", "Serotype 2" = "blue", 
                                "Serotype 3" = "green", "Serotype 4" = "purple"))+
  scale_x_continuous(breaks = seq(0, 3650*5, by = 365*5), labels = paste("Year", seq(0,50,5)))

##############################################################################################################
# Plot 5: Aggregated SEIR plot
##############################################################################################################
# Calculate the total number of susceptible, exposed, infected, and recovered individuals
total_susceptible <- rowSums(SEIRS_bd_results_df[, grep("^S_", names(SEIRS_bd_results_df))])
total_exposed <- rowSums(SEIRS_bd_results_df[, grep("^E_", names(SEIRS_bd_results_df))])
total_infected <- rowSums(SEIRS_bd_results_df[, grep("^I_", names(SEIRS_bd_results_df))])
total_recovered <- rowSums(SEIRS_bd_results_df[, grep("^R_", names(SEIRS_bd_results_df))])

# Create a data frame for plotting
plot_data_SEIR_bd <- data.frame(
  time = SEIRS_bd_results_df$time,
  Susceptible = total_susceptible,
  Exposed = total_exposed,
  Infected = total_infected,
  Recovered = total_recovered
)

ggplot(plot_data_SEIR_bd, aes(x = time)) +
  geom_line(aes(y = Susceptible, colour = "Susceptible")) +
  geom_line(aes(y = Exposed, colour = "Exposed")) +
  geom_line(aes(y = Infected, colour = "Infected")) +
  geom_line(aes(y = Recovered, colour = "Recovered")) +
  labs(title = "Aggregate SEIR Compartment Dynamics", 
       y = "Number of Individuals", x = "Time (years)") +
  scale_color_manual(values = c("Susceptible" = "blue", "Exposed" = "orange", 
                                "Infected" = "red", "Recovered" = "green"))+
  scale_x_continuous(breaks = seq(0, 3650*5, by = 365*5), labels = paste("Year", seq(0,50,5))) +
  theme_minimal()

##############################################################################################################
# Plot 6: Number of people in infectious state by their previous infected number
##############################################################################################################
# Count the number of "N" states to represent the number of uninfected serotypes 
count_inf_sero <- function(state) {
  return(nchar(gsub("[^N]", "", state)))
}

# Create a data frame to store the previous infected number at each time point
infected_number_data <- data.frame(time = SEIRS_bd_results_df$time)

# Calculate the sum of people in states 'i' for each time point
# i is the number of serotypes that the individuals is susceptible to, (3-i) is the infected serotypes
for (i in 0:3) {
  infected_number_data[[paste0("Infected_", (3-i), "_times")]] <- 
    rowSums(SEIRS_bd_results_df[infectious_states][sapply(infectious_states, count_inf_sero) == i])
}

# Plotting
ggplot(infected_number_data, aes(x = time)) +
  geom_line(aes(y = Infected_0_times, color = "Infected 0 times")) +
  geom_line(aes(y = Infected_1_times, color = "Infected 1 time")) +
  geom_line(aes(y = Infected_2_times, color = "Infected 2 times")) +
  geom_line(aes(y = Infected_3_times, color = "Infected 3 times")) +
  labs(title = "Number of People in Infectious States by Previous Infected Number",
       x = "Year",
       y = "Number of People in Infectious states",
       color = "Previous Infected Number") +
  ylim(c(0, 2*10^5)) +
  scale_x_continuous(breaks = seq(0, 3650*5, by = 365*5), labels = paste("Year", seq(0,50,5))) +
  theme_minimal()


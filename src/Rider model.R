# ################################################################# #
#### LOAD LIBRARY AND DEFINE CORE SETTINGS                       ####
# ################################################################# #

# Clear memory to avoid conflicts with previous variables
timestamp()
rm(list = ls())

# Load the Apollo library for discrete choice modeling
library(apollo)

# Initialize Apollo settings
apollo_initialise()

# Set core model controls
apollo_control = list(
  modelName       = "Rider_mix_hsk",  # Model name (used for saving outputs)
  modelDescr      = "Mixed model (TT, money, info & error component) in preference space", 
  indivID         = "ID",  # Column identifying individuals in panel data
  nCores          = 11,  # Number of CPU cores to use for estimation
  outputDirectory = "out_mix_hsk"  # Directory where results will be saved
)

# ################################################################# #
#### LOAD DATA AND APPLY ANY TRANSFORMATIONS                     ####
# ################################################################# #

# Load dataset from CSV file
database = read.csv("F_use_070224.csv", header = TRUE)

# ################################################################# #
#### DEFINE MODEL PARAMETERS                                     ####
# ################################################################# #

# Define model parameters (estimated during model fitting)
apollo_beta = c(
  asc_sco               = 0,  # Alternative-specific constant for 'scoober'
  asc_alt               = 0,  # Alternative-specific constant for 'alt'
  b_ageL34              = 0,  # Coefficient for age group (1-3)
  b_nonEUcitizen        = 0,  # Coefficient for non-EU citizenship
  b_timepres_yes        = 0,  # Coefficient for time presence (binary)
  b_expL1               = 0,  # Coefficient for an experience of less than one year
  mu_log_b_tt           = -3, # Mean of log-normal distribution for travel time
  sigma_log_b_tt        = 0,  # Standard deviation for travel time
  mu_log_b_info         = -3, # Mean of log-normal distribution for information
  sigma_log_b_info      = 0,  # Standard deviation for information
  mu_log_b_money        = -3, # Mean of log-normal distribution for money
  sigma_log_b_money     = 0,  # Standard deviation for money
  sigma_hsk             = 0   # Standard deviation for unobserved heterogeneity
)

# Specify parameters that will be held fixed during estimation
apollo_fixed = c("asc_alt")  # 'asc_alt' remains unchanged

# ################################################################# #
#### SPECIFY RANDOM COMPONENTS & DRAW SETTINGS                  ####
# ################################################################# #

# Define Monte Carlo draws for capturing random taste heterogeneity
apollo_draws = list(
  interDrawsType = "halton",  # Use Halton sequences for inter-individual draws
  interNDraws    = 80000,     # Number of draws per individual
  interUnifDraws = c(),       # No uniform draws
  interNormDraws = c("draws_tt", "draws_info", "draws_money", "draws_hsk"),
  intraDrawsType = "halton",  # Halton draws for intra-individual variations
  intraNDraws    = 0,         # No intra-individual draws
  intraUnifDraws = c(),
  intraNormDraws = c()
)

# Function defining how random coefficients are generated
apollo_randCoeff = function(apollo_beta, apollo_inputs){
  randcoeff = list()
  randcoeff[["b_tt"]]    = -exp(mu_log_b_tt + sigma_log_b_tt * draws_tt) # Travel time
  randcoeff[["b_info"]]  = exp(mu_log_b_info + sigma_log_b_info * draws_info) # Information
  randcoeff[["b_money"]] = exp(mu_log_b_money + sigma_log_b_money * draws_money) # Money
  randcoeff[["hsk_sco"]] = sigma_hsk * draws_hsk  # Heterogeneity component
  return(randcoeff)
}

# Validate all model inputs
apollo_inputs = apollo_validateInputs()

# ################################################################# #
#### DEFINE MODEL AND LIKELIHOOD FUNCTION                        ####
# ################################################################# #

# Function to compute choice probabilities
apollo_probabilities = function(apollo_beta, apollo_inputs, functionality = "estimate"){
  apollo_attach(apollo_beta, apollo_inputs)
  on.exit(apollo_detach(apollo_beta, apollo_inputs))
  
  # Initialize probabilities list
  P = list()
  
  # Define utility functions for each alternative
  V = list()
  V[["scoober"]] = asc_sco + b_tt * time_sco + b_info * (saf_info / 10) + 
    b_nonEUcitizen * (citizen == 3) + 
    b_ageL34 * (age %in% c(1,2,3)) + 
    b_timepres_yes * (time_pres == 2) + 
    b_expL1 * (job_exp == 1) + 
    b_money * (mon_sco * 10) + hsk_sco
  
  V[["alt"]] = asc_alt + b_tt * time_alt
  
  # Define settings for the multinomial logit (MNL) model
  mnl_settings = list(
    alternatives  = c(scoober = 1, alt = 2),
    avail         = list(scoober = 1, alt = 1),
    choiceVar     = choice,
    utilities     = V
  )
  
  # Compute MNL probabilities
  P[["model"]] = apollo_mnl(mnl_settings, functionality)
  
  # Aggregate probabilities across panel observations
  P = apollo_panelProd(P, apollo_inputs, functionality)
  P = apollo_avgInterDraws(P, apollo_inputs, functionality)
  P = apollo_prepareProb(P, apollo_inputs, functionality)
  return(P)
}

# ################################################################# #
#### MODEL ESTIMATION AND OUTPUT                                 ####
# ################################################################# #

# Load previous estimation results to use as starting values
prevmodel = apollo_loadModel("F_com_mix_hsk1_final")

# Estimate model parameters using maximum likelihood
model = apollo_estimate(prevmodel$estimate, apollo_fixed, apollo_probabilities, apollo_inputs)

# Display and save model output
apollo_modelOutput(model, list(printPVal = 2))
apollo_saveOutput(model, list(printPVal = 2))

# ################################################################# #
##### POST-PROCESSING AND INTERPRETATION                        ####
# ################################################################# #

# Compute conditional and unconditional parameter distributions
unconditionals = apollo_unconditionals(model, apollo_probabilities, apollo_inputs)
conditionals   = apollo_conditionals(model, apollo_probabilities, apollo_inputs)

# Compute and summarize key measures
VRR = (unconditionals[["b_info"]] / unconditionals[["b_tt"]])
WTA = (unconditionals[["b_tt"]] / unconditionals[["b_money"]])

cat("Mean VRR:", mean(VRR), " SD VRR:", sd(VRR), "\n")
cat("Mean WTA:", mean(WTA), " SD WTA:", sd(WTA), "\n")

# Print summaries of conditional distributions
summary(conditionals[["b_info"]] / conditionals[["b_tt"]])
summary(conditionals[["b_tt"]] / conditionals[["b_money"]])

# Close all open sinks (output redirection)
apollo_sink()

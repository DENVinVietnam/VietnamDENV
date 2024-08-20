##############################################################################################################
## This file write out a csv file using the following process to generate a fake data
## Note that although the process of generating fake data follows the general rules in the original data, 
## the whole dynamic and underlying distributions may not be fully captured, 
## so the results are different from those presented in the dissertation.
##############################################################################################################

# Load necessary libraries
library(pacman)
p_load(dplyr, lubridate, stringr)
# lubridate: for dates formatting

##############################################################################################################
## Data Generating
##############################################################################################################
# For reproducibility
set.seed(123) 

# Number of records
n <- 4856

# Generate MAC-ELISA and GAC-ELISA results
MAC_ELISA <- sample(c("Positive", "Negative", "Equivocal"), n, replace = TRUE, prob = c(0.18,0.8,0.02))
GAC_ELISA <- sample(c("Positive", "Negative", "Equivocal"), n, replace = TRUE, prob = c(0.6,0.3,0.01))

# Generate PRNT samples column based on MAC-ELISA and GAC-ELISA results
prnt_samples <- !(MAC_ELISA == "Negative" & GAC_ELISA == "Negative")

# Define PRNT50 values and their probabilities
prnt50_values <- c(20, 40, 80, 160, 320, 640, 1280, 2560, 5120, 10240, 20480, 40960, 81920, 163840, NA)
prnt50_probs <- c(0.01, 0.01, 0.03, 0.04, 0.03, 0.01, 0.01, 0.005, 
                    0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.70)

# Define probabilities for D1 and D2 to have higher PRNT50 values (around 70% NA)
prnt50_probs_high <- c(0.01, 0.01, 0.05, 0.10, 0.05, 0.03, 0.03, 
                         0.03, 0.03, 0.05, 0.05, 0.03, 0.03, 0.03, 0.70)

# Generate PRNT50 columns
prnt50_D1 <- ifelse(prnt_samples, sample(prnt50_values, n, replace = TRUE, prob = prnt50_probs_high), NA)
prnt50_D2 <- ifelse(prnt_samples, sample(prnt50_values, n, replace = TRUE, prob = prnt50_probs_high), NA)
prnt50_D3 <- ifelse(prnt_samples, sample(prnt50_values, n, replace = TRUE, prob = prnt50_probs), NA)
prnt50_D4 <- ifelse(prnt_samples, sample(prnt50_values, n, replace = TRUE, prob = prnt50_probs), NA)

# Generate PRNT50 serotype column
# This column contains all the serotypes with a detectable antibody level for each individual, separated by ","
prnt50 <- sapply(1:n, function(i) {
  if (prnt_samples[i]) {
    serotypes <- c(ifelse(!is.na(prnt50_D1[i]), "D1", NA),
                   ifelse(!is.na(prnt50_D2[i]), "D2", NA),
                   ifelse(!is.na(prnt50_D3[i]), "D3", NA),
                   ifelse(!is.na(prnt50_D4[i]), "D4", NA))
    serotypes <- na.omit(serotypes)
    if (length(serotypes) > 0) {
      paste0(serotypes, collapse = ",")
    } else {
      "Neg"
    }
  } else {
    "Neg"
  }
})

# Generate Serotype column, which stores the PCR results
pcr_values <- c("D1", "D2", "D3", "D4", "Neg")
pcr_prob <- c(0.12, 0.18, 0.01, 0.09, 0.60)
pcr <- sample(pcr_values, n, replace = TRUE, prob = pcr_prob)

# Generate onset dates with the all the possible formats occurred in the original data
onset_date_raw <- sample(seq.Date(as.Date("2020-09-01"), as.Date("2022-07-31"), by = "day"), 
                         n, replace = TRUE)
onset_date <- sapply(onset_date_raw, function(x) {
  formats <- c("%Y-%m-%d", "%d/%m/%Y", "%d %b %y")
  format_choice <- sample(formats, 1)
  format(x, format = format_choice)
})

# Generate examination dates with all the possible formats occurred in the original data
exam_date_raw <- onset_date_raw + sample(0:10, n, replace = TRUE)
exam_date <- sapply(exam_date_raw, function(x) {
  formats <- c("%Y-%m-%d", "%d/%m/%Y", "%d %b %y")
  format_choice <- sample(formats, 1)
  format(x, format = format_choice)
})

# Ensure the majority have a difference within 5 days between the onset and examination dates
exam_date_raw <- as.Date(ifelse(runif(n) < 0.8, 
                                onset_date_raw + sample(0:5, n, replace = TRUE), exam_date_raw))
exam_date <- sapply(exam_date_raw, function(x) {
  formats <- c("%Y-%m-%d", "%d/%m/%Y", "%d %b %y")
  format_choice <- sample(formats, 1)
  format(x, format = format_choice)
})

# Create the ID-NIHE column, the unique index for each individual
index <- 1:n
ID_NIHE <- paste0("VT", substr(as.character(exam_date_raw), 3, 4), ".", sprintf("%04d", index))

# Range of birth year in the data
birthyear <- 1924:2022

# Higher probabilities for years after 2012 (age under 10), moderate for 1993-2012 (age under 30), 
# and lower for 1924-1992 (older than 30)
prob_dist <- c(
  rep(0.01, length(1924:1952)), # Rare probability for 1924-1952
  rep(0.04, length(1953:1992)),  # Lower probability for 1953-1992
  rep(0.20, length(1993:2002)),  # Middle probability for 1993-2002
  rep(0.45, length(2003:2012)),   # Higher probability for 2003-2012
  rep(0.30, length(2013:2022))   # Higher probability for 2013-2022
)

# Generate birth year
birth_year <- sample(birthyear, n, replace = TRUE, prob = prob_dist)

## Generate gender
gender <- sample(c("Male", "Female"), n, replace = TRUE)

# Create the data frame, follow the name in the orginal data
DENV_raw <- data.frame(
  `ID-NIHE` = ID_NIHE,
  `MAC-ELISA- Final result` = MAC_ELISA,
  `GAC-ELISA-Final result` = GAC_ELISA,
  Serotype = pcr,
  `PRNT samples` = prnt_samples,
  PRNT50 = prnt50,
  `PRNT50 DENV1` = prnt50_D1,
  `PRNT50 DENV2` = prnt50_D2,
  `PRNT50 DENV3` = prnt50_D3,
  `PRNT50 DENV4` = prnt50_D4,
  `Date of birth ` = birth_year,
  Gender = gender,
  `Date of onset` = onset_date,
  `Examination date` = exam_date,
  check.names = FALSE
)

# Extract the index from ID-NIHE
ID_index <- as.numeric(substring(DENV_raw$`ID-NIHE`, 
                                 nchar(DENV_raw$`ID-NIHE`) - 3, nchar(DENV_raw$`ID-NIHE`)))

# Order by the index
DENV_raw <- DENV_raw[order(ID_index), ]

##############################################################################################################
## Initial data cleaning by extracting the onset year and age for analysis
## Write out as a xlsx file
##############################################################################################################
# Extract onset years and months
DENV_raw <- DENV_raw %>%
  mutate(
    # Using 'Date of onset' and handle various formats, if it's NA then use `Examination date`
    OnsetDate = if_else(
      is.na(`Date of onset`) | `Date of onset` == "",
      case_when(
        str_detect(`Examination date`, "-") ~
          as.Date(`Examination date`, format = "%Y-%m-%d"), # e.g 19-10-2020
        str_detect(`Examination date`, "/") ~
          as.Date(`Examination date`, format = "%d/%m/%Y"), # e.g 19/10/2020
        str_detect(`Examination date`, "\\d{1,2} [A-Za-z]{3} \\d{2}") ~
          as.Date(`Examination date`, format = "%d %b %y") #e.g 19 Oct 20
      ),
      case_when(
        str_detect(`Date of onset`, "-") ~ as.Date(`Date of onset`, format = "%Y-%m-%d"),
        str_detect(`Date of onset`, "/") ~ as.Date(`Date of onset`, format = "%d/%m/%Y"),
        str_detect(`Date of onset`, "\\d{1,2} [A-Za-z]{3} \\d{2}") ~
                    as.Date(`Date of onset`, format = "%d %b %y")
      )
    ),
    # Extract Year and Month from the parsed date
    OnsetYear = as.integer(format(OnsetDate, "%Y")),
    OnsetMonth = format(OnsetDate, "%m")
  ) %>%
  # Use the year from Examination date if not between 2020 and 2022
  mutate(
    OnsetYear = if_else(
      OnsetYear >= 2020 & OnsetYear <= 2022,
      OnsetYear,
      as.integer(format(case_when(
        str_detect(`Examination date`, "-") ~ as.Date(`Examination date`, format = "%Y-%m-%d"),
        str_detect(`Examination date`, "/") ~ as.Date(`Examination date`, format = "%d/%m/%Y"),
        str_detect(`Examination date`, "\\d{1,2} [A-Za-z]{3} \\d{2}") ~
                                       as.Date(`Examination date`, format = "%d %b %y")
      ), "%Y"))
    )
  )

# Extract birth year from "Date of birth"
DENV_raw <- DENV_raw %>%
  mutate(
    # Parse Date of birth in various formats
    BirthDate = case_when(
      str_detect(`Date of birth `, "/") ~ as.Date(`Date of birth `, format = "%d/%m/%Y"),
      str_detect(`Date of birth `, "-") ~ as.Date(`Date of birth `, format = "%Y-%m-%d"),
      # Add a month and day if the entry is year-only to unify the foramt
      str_detect(`Date of birth `, "^\\d{4}$") ~
                    as.Date(paste(`Date of birth `, "01", "01", sep = "-"), format = "%Y-%m-%d"),
      TRUE ~ as.Date(NA)
    ),
    # Extract Year 
    BirthYear = if_else(
      str_detect(`Date of birth `, "^\\d{4}$"), # Directly use the year if only years 
      as.integer(`Date of birth `), # Convert to integer to match the type
      as.integer(format(BirthDate, "%Y")) # If it is date, extract the year from the birth date
    )
  )
#
# # Calculate age of the individuals
DENV_raw$Age_clean <- ifelse(DENV_raw$OnsetYear - as.numeric(DENV_raw$BirthYear) < 0, NA,
                                    DENV_raw$OnsetYear - as.numeric(DENV_raw$BirthYear))

# Write out as in csv file
write.csv(DENV_raw, "S1 Data. Vietnam Dengue Fake Data.csv", row.names = FALSE)
##############################################################################################################


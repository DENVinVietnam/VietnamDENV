##############################################################################################################
## This file contains:
## Section 1: Basic Data Preparation and Data Cleaning
## Section 2: Update PRNT50 columns based on the thresholds
## Section 3: Differentiate current and pat infections
## Section 4: Data Preparation for multinomial regression
##
## Analysis 1: Multinomial regression using PCR+ cases v.s PCR- & IgM
## Analysis 2: Multinomial regression using data excluding S9
## Analysis 3: Multinomial regression using all the data (s1-s9)
## Analysis 4: Use all the data from S1-S8 with 20% individuals from S10 using EM algorithm
##
## Section 5: Data Visualization
## - Plot 1a (bar): visualizes age distribution by serotype
## - Plot 1b (bar): visualizes gender distribution by serotype
## - Plot 2 (line + bar): visualizes yearly dengue infection counts and percentages by serotype
## - Plot 3 (violin): visualizes the distribution of PRNT50 values by serotype
## - Plot 4 (forest): visualizes the odds ratios with confidence intervals
## - Plot 5 (bar): Distribution between examination date and onset date
## - Plot 6 (line): Cumulative distribution plots for PRNT50 values 
##############################################################################################################

library(pacman)
p_load(dplyr, lubridate, stringr, tidyr, ggplot2, nnet, RColorBrewer, matrixStats)

# Read enciphered data (only available on Windows, for processing the original data), no need for the fake one
#library("excel.link") 

###################################################
## Section 1: Data Preparation and Data Cleaning###
###################################################

## Initial Data Cleaning######################################################################################
## This process is already done in S1 Code for the fake data
## Initial data cleaning by extracting the onset year and age for analysis
## Write out as a csv file
##############################################################################################################
# ### Import and read enciphered data
# DENV_raw <- xl.read.file("PW_VungTau_Dengue-Result-Summary_240621updated.xlsx", 
#             xl.sheet = 'DS', password = "20240621", write.res.password="pass")
# 
# # Rename duplicated column names
# names(DENV_raw)[4] <- "Index.Value.1"
# names(DENV_raw)[7] <- "Index.Value.2"
# 
# # Extract onset years and months
# DENV_raw <- DENV_raw %>%
#   mutate(
#     # Using 'Date of onset', if it's NA then use `Examination date`
#     OnsetDate = if_else(
#       is.na(`Date of onset`) | `Date of onset` == "",
#       case_when(
#         str_detect(`Examination date`, "-") ~ 
#           as.Date(`Examination date`, format = "%Y-%m-%d"), # e.g 19-10-2020
#         str_detect(`Examination date`, "/") ~ 
#           as.Date(`Examination date`, format = "%d/%m/%Y"), # e.g 19/10/2020
#         str_detect(`Examination date`, "\\d{1,2} [A-Za-z]{3} \\d{2}") ~ 
#           as.Date(`Examination date`, format = "%d %b %y") #e.g 19 Oct 20
#       ),
#       case_when(
#         str_detect(`Date of onset`, "-") ~ as.Date(`Date of onset`, format = "%Y-%m-%d"),
#         str_detect(`Date of onset`, "/") ~ as.Date(`Date of onset`, format = "%d/%m/%Y"),
#         str_detect(`Date of onset`, "\\d{1,2} [A-Za-z]{3} \\d{2}") ~ 
#                     as.Date(`Date of onset`, format = "%d %b %y")
#       )
#     ),
#     # Extract Year and Month from the OnsetDate
#     OnsetYear = as.integer(format(OnsetDate, "%Y")),
#     OnsetMonth = format(OnsetDate, "%m")
#   ) %>%
#   # Use the year from Examination date if not between 2020 and 2022
#   mutate(
#     OnsetYear = if_else(
#       OnsetYear >= 2020 & OnsetYear <= 2022,
#       OnsetYear,
#       as.integer(format(case_when(
#         str_detect(`Examination date`, "-") ~ as.Date(`Examination date`, format = "%Y-%m-%d"),
#         str_detect(`Examination date`, "/") ~ as.Date(`Examination date`, format = "%d/%m/%Y"),
#         str_detect(`Examination date`, "\\d{1,2} [A-Za-z]{3} \\d{2}") ~ 
#                                        as.Date(`Examination date`, format = "%d %b %y")
#       ), "%Y"))
#     )
#   )
# 
# # Extract birth year from "Date of birth"
# DENV_raw <- DENV_raw %>%
#   mutate(
#     # Extract birth year from "Date of birth"
#     BirthDate = case_when(
#       str_detect(`Date of birth `, "/") ~ as.Date(`Date of birth `, format = "%d/%m/%Y"),
#       str_detect(`Date of birth `, "-") ~ as.Date(`Date of birth `, format = "%Y-%m-%d"),
#       # Add a month and day if the entry is year-only to unify the foramt
#       str_detect(`Date of birth `, "^\\d{4}$") ~ 
#                     as.Date(paste(`Date of birth `, "01", "01", sep = "-"), format = "%Y-%m-%d"), 
#       TRUE ~ as.Date(NA)
#     ),
#     # Extract Year from the birth date or directly from the year-only entry
#     BirthYear = if_else(
#       str_detect(`Date of birth `, "^\\d{4}$"), # Directly use the year if the entry is year-only
#       `Date of birth `,
#       format(BirthDate, "%Y") 
#     )
#   )
# 
# # Calculate age of the individuals
# DENV_raw$Age_clean <- ifelse(DENV_raw$OnsetYear - as.numeric(DENV_raw$BirthYear) < 0, NA, 
# DENV_raw$OnsetYear - as.numeric(DENV_raw$BirthYear))
# DENV_raw <- DENV_raw %>%
#   mutate(
#     Age = as.numeric(as.character(`Age `)), 
#     Age_clean = as.numeric(as.character(Age_clean)),  
#     Age_clean = coalesce(Age_clean, `Age `)  # fill NA values with the values in age column
#   )
# DENV_raw <- DENV_raw %>% select(-Age) 

# Write out in csv file
# write.csv(DENV_raw,"DENV version 5-age cleaned.csv")
##############################################################################################################


## Data Cleaning 2############################################################################################
## For the original data set
## Further data cleaning to the excel generated above for handling missing data, 
##############################################################################################################
# 72:4927 is the range of data where we have complete test results for analysis, 
# which is from Sep 2020 to Jul 2022
# DENV <- read.csv("DENV version 5-age cleaned.csv", check.names = FALSE)[c(72:4927),] # Original Data
# # Re-assign the row names to start from 1
# row.names(DENV) <- NULL
# 
# 
# ## Data cleaning
# DENV_cleaned <- DENV %>%
#  # Convert factors into character
#   mutate_if(is.factor, as.character) %>%
#   # Remove extra whitespace from character columns
#   mutate_if(is.character, function(x) str_squish(x)) %>% 
#   # PRNT50 column stores the serotype with measurable PRNT50 values, e.g D1,D3,D4
#   # For PRNT50 DENV1 values, we have 2 types of test results in the data set (D1-VN and D1),
#   # as the data side suggested, we choose to use D1 for a consistent standard
#   mutate(PRNT50 = if_else(
#     # For those only have D1-VN, treat them as Neg against all serotypes
#     PRNT50 %in% c("D1-VN","D1-VN,"), "Neg", 
#     # Otherwise, remove D1-VN in PRNT50 column
#     if_else(grepl("D1-VN", PRNT50), str_remove_all(PRNT50, "D1-VN,?|,?D1-VN"), PRNT50) 
#   )) %>% 
#   # filter out the unexpected Serotype rows
#   filter(Serotype != "" & !is.na(Serotype) & Serotype != "chưa làm"
#          & Serotype != "No sample" & Serotype != "Re-check" & Serotype != "ND"
#          & Serotype != "Out of sample" & Serotype != "Repeat")
# 
# ## Replace NA or empty values in PRNT50 by "empty" since NA is not comparable using "=="
# DENV_cleaned$PRNT50 <- ifelse(is.na(DENV_cleaned$PRNT50) | DENV_cleaned$PRNT50 == "", 
#                               "empty", DENV_cleaned$PRNT50)
# # Remove unexpected PRNT50 values: "Not enough sample"
# DENV_cleaned <- DENV_cleaned[DENV_cleaned$PRNT50 != "Not enough sample", ]
# 
# 
# DENV_cleaned <- DENV_cleaned %>%
#   # Convert `PRNT samples` to logical values 
#   # PRNT samples has values TRUE/FALSE/NA,indicating whether PRNT50 tests is performed or not,
#   # NA rows are associated with PRNT50 value equals "Not enough sample" which are already removed above
#   mutate(`PRNT samples` = as.logical(`PRNT samples`)) %>%
#   # Remove those PRNT50 samples is TRUE but no PRNT50 entries
#   filter(!(`PRNT samples` & (PRNT50 == "empty"))) %>% 
#   # Replace PRNT50 with "Neg" if both `MAC-ELISA- Final result` and `GAC-ELISA-Final result` are "Negative",
#   # since PRNT50 test is not performed when both IgG and IgM are negative and 
#   # we assume they negative against all serotypes in our study
#   mutate(PRNT50 = if_else(`MAC-ELISA- Final result` == "Negative"
#                           & `GAC-ELISA-Final result` == "Negative", "Neg", PRNT50)) %>%
#   # Filter out the data within the correct range of sample collecting years
#   filter(OnsetYear %in% c(2020, 2021, 2022))
# 
# 
# ## Unify the IgG and IgM formats
# DENV_cleaned <- DENV_cleaned %>%
#   mutate(
#     `MAC-ELISA- Final result` = case_when(
#       `MAC-ELISA- Final result` %in% c("NEG", "Negative") ~ "Negative",
#       `MAC-ELISA- Final result` %in% c("POS", "Positive") ~ "Positive",
#       `MAC-ELISA- Final result` %in% c("EQUIVOCAL", "equivocal") ~ "Equivocal",
#       TRUE ~ `MAC-ELISA- Final result` 
#     )
#   ) %>%
#   mutate(
#     `GAC-ELISA-Final result` = case_when(
#       `GAC-ELISA-Final result` %in% c("NEG", "Negative") ~ "Negative",
#       `GAC-ELISA-Final result` %in% c("POS", "Positive") ~ "Positive",
#       `GAC-ELISA-Final result` %in% c("EQUIVOCAL", "equivocal") ~ "Equivocal",
#       TRUE ~ `GAC-ELISA-Final result` 
#     )
#   )

##############################################################################################################
# Import the generated fake data set (keep the column name as in the original data)
## unifying column formats and updating incorrect PRNT50 values
DENV_cleaned <- read.csv("S1 Data. Vietnam Dengue Fake Data.csv", check.names = FALSE)

DENV_cleaned_comp <- DENV_cleaned %>%
  filter(!is.na(Age_clean) & Age_clean != "",   # remove rows where age is NA or empty
         !is.na(Gender) & Gender != "",  # remove rows where gender is NA or empty
         !is.na(OnsetYear) & OnsetYear != "") # remove rows where OnsetYear is NA or empty 

# Convert age into numeric values
DENV_cleaned_comp$`PRNT samples` <- as.logical(DENV_cleaned_comp$`PRNT samples`)
DENV_cleaned_comp$Age_clean <- as.numeric(DENV_cleaned_comp$Age_clean)
##############################################################################################################


## Function:clean_prnt50######################################################################################
## Function to unify the format of the PRNT50 column
## Input: prnt50 (string) - Characters in the PRNT50 column which were used to clean 
## Output (string): A cleaned PRNT50 characters separated by comma and no extra white space
##############################################################################################################
clean_prnt50 <- function(prnt50) {
  if (prnt50 == "Neg" || is.na(prnt50)) {
    return(prnt50)
  }
  
  # remove spaces at the begining or end
  prnt50 <- str_trim(prnt50)
  # remove commas that appear at the beginning or the end of a string.
  prnt50 <- str_remove_all(prnt50, "^,|,$")
  
  # Split the string in PRNT50 column, then connect by ","
  prnt50_split <- str_split(prnt50, ",", simplify = TRUE)
  # remove the original inside spaces
  prnt50_split <- str_trim(prnt50_split)
  prnt50 <- str_c(prnt50_split, collapse = ",")
  
  prnt50
}

## Apply to the cleaned data set
DENV_cleaned_comp <- DENV_cleaned_comp %>%
  mutate(PRNT50 = sapply(PRNT50, clean_prnt50))

## Function: check_PRNT50_match###############################################################################
## Function to check whether PRNT50 column matches PRNT50 DENV1-4 columns
## Input: row - A row of "PRNT50", "PRNT50 DENV1", "PRNT50 DENV2", "PRNT50 DENV3", "PRNT50 DENV4" 
## Output (logic): TRUE if PRNT50 matches individual columns, FALSE otherwise
##############################################################################################################
check_PRNT50_match <- function(row) {
  # Extract serotypes in PRNT50, separated by ","
  prnt50_serotypes <- unlist(strsplit(as.character(row$PRNT50), ","))
  
  # An empty list to store actual serotypes based on PRNT50 DENV1-4 columns
  actual_serotypes <- c()
  if (!is.na(row$`PRNT50 DENV1`) && row$`PRNT50 DENV1` != "") actual_serotypes <- c(actual_serotypes, "D1")
  if (!is.na(row$`PRNT50 DENV2`) && row$`PRNT50 DENV2` != "") actual_serotypes <- c(actual_serotypes, "D2")
  if (!is.na(row$`PRNT50 DENV3`) && row$`PRNT50 DENV3` != "") actual_serotypes <- c(actual_serotypes, "D3")
  if (!is.na(row$`PRNT50 DENV4`) && row$`PRNT50 DENV4` != "") actual_serotypes <- c(actual_serotypes, "D4")
  
  # If no PRNT50 values against all serotypes, then PRNT50 column stores "Neg"
  if (!is.na(row$PRNT50) && row$PRNT50 == "Neg" && 
      all(is.na(row$`PRNT50 DENV1`) | row$`PRNT50 DENV1` == "" | row$`PRNT50 DENV1` == "Neg", 
          is.na(row$`PRNT50 DENV2`) | row$`PRNT50 DENV2` == "" | row$`PRNT50 DENV2` == "Neg", 
          is.na(row$`PRNT50 DENV3`) | row$`PRNT50 DENV3` == "" | row$`PRNT50 DENV3` == "Neg", 
          is.na(row$`PRNT50 DENV4`) | row$`PRNT50 DENV4` == "" | row$`PRNT50 DENV4` == "Neg")) {
    return(TRUE)
  }
  # Compare serotypes in the data set vs actual serotypes
  setequal(prnt50_serotypes, actual_serotypes)
}


## Apply the function to each rows of the data and
# create a column to store whether the values in PRNT50 matches the actual serotype
DENV_cleaned_comp <- DENV_cleaned_comp %>%
  rowwise() %>%
  # pick: returns a data frame of the selected columns 
  mutate(PRNT50_check = check_PRNT50_match(pick(c("PRNT50", "PRNT50 DENV1", "PRNT50 DENV2", 
                                                         "PRNT50 DENV3", "PRNT50 DENV4")))) %>%
  ungroup()


## Function: update_wrong_PRNT50###############################################################################
## Function to extract the wrong PRNT50 rows and update them
## Input: row - A row of `PRNT50 DENV1`, `PRNT50 DENV2`, `PRNT50 DENV3`, `PRNT50 DENV4`
## Output (string): prnt50 - Updated PRNT50 characters that matches the PRNT50 DENV1-4 columns
##############################################################################################################
update_wrong_PRNT50 <- function(row) {
  serotypes <- c()
  if (!is.na(row$`PRNT50 DENV1`) && row$`PRNT50 DENV1` != "" && row$`PRNT50 DENV1` != "Neg") 
      serotypes <- c(serotypes, "D1")
  if (!is.na(row$`PRNT50 DENV2`) && row$`PRNT50 DENV2` != "" && row$`PRNT50 DENV2` != "Neg") 
      serotypes <- c(serotypes, "D2")
  if (!is.na(row$`PRNT50 DENV3`) && row$`PRNT50 DENV3` != "" && row$`PRNT50 DENV3` != "Neg") 
      serotypes <- c(serotypes, "D3")
  if (!is.na(row$`PRNT50 DENV4`) && row$`PRNT50 DENV4` != "" && row$`PRNT50 DENV4` != "Neg") 
      serotypes <- c(serotypes, "D4")
  
  paste(serotypes, collapse = ",")
}

## Update the wrong PRNT50 rows
DENV_cleaned_comp <- DENV_cleaned_comp %>%
  rowwise() %>%
  mutate(
    PRNT50 = if_else(
      !PRNT50_check, 
      update_wrong_PRNT50(pick(c(`PRNT50 DENV1`, `PRNT50 DENV2`, `PRNT50 DENV3`, `PRNT50 DENV4`))), 
      PRNT50
    )
  ) %>%
  ungroup()


##############################################################
## Section 2: Update PRNT50 columns based on the thresholds###
##############################################################

##############################################################################################################
## Convert PRNT50 DENVi columns into numeric values for comparing with thresholds
DENV_cleaned_comp$`PRNT50 DENV1` <- as.numeric(DENV_cleaned_comp$`PRNT50 DENV1`)
DENV_cleaned_comp$`PRNT50 DENV2` <- as.numeric(DENV_cleaned_comp$`PRNT50 DENV2`)


DENV_cleaned_comp$`PRNT50 DENV3` <- as.numeric(DENV_cleaned_comp$`PRNT50 DENV3`)
DENV_cleaned_comp$`PRNT50 DENV4` <- as.numeric(DENV_cleaned_comp$`PRNT50 DENV4`)


## Function: replace_with_na##################################################################################
## Function to replace PRNT50 values less than or equal to the lower threshold with NA
## Input: 
##  - df: Cleaned DENV data frame
##  - threshold_low: numeric threshold value for differentiating dengue cases and non-cases, 
##                   values below the lower threshold indicating non-cases (replace with NA)
## Output: 
##  - df: updated data frame by replacing the rows where PRNT50 values smaller than threshold_low by NA
##############################################################################################################
# Replace those values less than or equal to the threshold with NA
replace_with_na <- function(df, threshold_low) {
  prnt50_cols <- c("PRNT50 DENV1", "PRNT50 DENV2", "PRNT50 DENV3", "PRNT50 DENV4")
  
  # replace the values less than the threshold for each column
  df[prnt50_cols] <- lapply(df[prnt50_cols], function(column) {
    column[column <= threshold_low] <- NA
    return(column)
  })
  
  return(df)
}


## Function: update_prnt50_column#############################################################################
## Function to update the PRNT50 column based on individual PRNT50 DENVi columns after replacing the 
## lower-than threshold values with NA
## Input: 
##  - row: Each row of the data frame after applying function replace_with_na
##  - threshold_low: numeric threshold value for differentiating dengue cases and non-cases 
## Output: 
##  - updated_prnt50: updated PRNT50 characters after applying the threshold 
##############################################################################################################
## Update PRNT50 column
update_prnt50_column <- function(row, threshold_low) {
  # Split the PRNT50 column into individual serotypes by ","
  prnt50_serotypes <- unlist(strsplit(row["PRNT50"], ","))
  
  # Initialize a vector to store the updated serotypes
  updated_serotypes <- vector("character")
  
  # Loop through each serotype and check its corresponding PRNT50 value
  for (serotype in prnt50_serotypes) {
    # Get the PRNT50 value for the current serotype
    prnt50_value <- as.numeric(row[paste0("PRNT50 DENV", substr(serotype, 2, 2))])
    
    # Check if the PRNT50 value is greater than the threshold
    if (!is.na(prnt50_value) && prnt50_value > threshold_low) {
      # If it is, add the serotype to the updated serotypes list
      updated_serotypes <- c(updated_serotypes, serotype)
    }
  }
  
  # Update the PRNT50 column by combining the updated serotypes back into a string, separated by ","
  updated_prnt50 <- paste(updated_serotypes, collapse = ",")
  
  # Return the updated PRNT50 column value
  return(updated_prnt50)
}


# Define the thresholds
## Remember to re-run the code from beginning (right after the importation of data) if changing the thresholds
threshold_low <- 40 # lower threshold for differentiating between cases and non-cases
threshold_high <- 160# upper threshold for differentiating between low and high PRNT50 levels

# Apply the function to replace values with NA if not larger than low threshold
DENV_cleaned_comp <- replace_with_na(DENV_cleaned_comp, threshold_low)
# Apply the function to update the PRNT50 column for each row
DENV_cleaned_comp$PRNT50 <- apply(DENV_cleaned_comp, 1, update_prnt50_column, threshold_low)
# Replace empty PRNT50 values with "Neg", indicating negative against all serotypes
DENV_cleaned_comp$PRNT50[DENV_cleaned_comp$PRNT50 == ""] <- "Neg"
##############################################################################################################


#########################################################
## Section 3: Differentiate current and pat infections###
#########################################################

## Function: classify_PRNT50##################################################################################
## Function to classify PRNT50 values to differentiate between current and past infections
## based in PCR result and PRNT50 results
## Input:
##   - prnt50: PRNT50 values for a given individual 
##   - serotype: the PCR result for that individual
## Output:
##   - classification: A character indicating the classification of the infection status as 
##     "None", "1 homo (homologous)", "1 heter (heterlogous)", "2+ homo", or "2+ heter" 
##############################################################################################################
classify_PRNT50 <- function(prnt50, serotype) {
  
  # Split PRNT50 into individual serotypes, separated by ","
  prnt50_serotypes <- strsplit(prnt50, ",")
  
  # Determine the number of unique serotypes (excluding 'Neg' and 'NA') 
  # to tell whether it is infected by 1 or 2+ serotypes
  serotype_count <- sapply(prnt50_serotypes, function(x) length(unique(x[x %in% c("D1", "D2", "D3", "D4")])))
  
  # Determine if the serotype in the PRNT50 column is homologous or heterologous to the PCR result
  is_homo <- sapply(seq_along(serotype), function(i) serotype[i] %in% prnt50_serotypes[[i]])
  homo_heter <- ifelse(is_homo, "homo", "heter")
  
  # Create classification based on number of serotypes and is_homo status
  # Note that is PCR is negative, then the classification is always heterlogous
  classification <- ifelse(serotype_count == 0, "None",
                           ifelse(serotype_count == 1, paste("1", homo_heter, sep=" "),
                                  ifelse(serotype_count >= 2, paste("2+", homo_heter, sep=" "), NA)))
  
  return(classification)
}

# Apply the function to data
DENV_cleaned_comp$PRNT_classification <- mapply(classify_PRNT50, 
                                                DENV_cleaned_comp$PRNT50, DENV_cleaned_comp$Serotype)

DENV_cleaned_comp <- DENV_cleaned_comp %>%
  ## filter out co-infections
  filter(Serotype %in% c("D1", "D2", "D3", "D4", "Neg")) %>%
  ## Create 2 columns Current and Past to store the infections
  mutate(Current = NA_character_,
         Past = NA_character_)


## Categorize serotype to see if the individual is PCR positive or negative
DENV_cleaned_comp <- DENV_cleaned_comp %>%
  mutate(
    PCR_result = if_else(
      str_detect(Serotype, "D1|D2|D3|D4"),
      "positive",
      "negative"
    )
  )
##############################################################################################################


##############################################################################################################
## Basic Data Summary
##############################################################################################################
summary_table <- DENV_cleaned_comp %>%
  group_by(OnsetYear) %>%
  summarise(
    Number_of_Enrollees = n(),
    Male = as.integer(sum(Gender == "Male", na.rm = TRUE)),
    Female = sum(Gender == "Female", na.rm = TRUE),
    Mean_Age = if(all(is.na(Age_clean))) NA else round(mean(Age_clean, na.rm = TRUE)),
    SD_Age = if(all(is.na(Age_clean))) NA else round(sd(Age_clean, na.rm = TRUE)),
    #Min_Age = if(all(is.na(Age_clean))) NA else min(Age_clean, na.rm = TRUE),
    #Max_Age = if(all(is.na(Age_clean))) NA else max(Age_clean, na.rm = TRUE),
    Median_Age = if(all(is.na(Age_clean))) NA else median(Age_clean, na.rm = TRUE),
    Q1_Age = if(all(is.na(Age_clean))) NA else quantile(Age_clean, probs = 0.25, na.rm = TRUE),
    Q3_Age = if(all(is.na(Age_clean))) NA else quantile(Age_clean, probs = 0.75, na.rm = TRUE),
    IgM = count(`MAC-ELISA- Final result` == "Positive"),
    IgG = count(`GAC-ELISA-Final result` == "Positive"),
    PCR_D1 = count(Serotype == "D1"),
    PCR_D2 = count(Serotype == "D2"),
    PCR_D3 = count(Serotype == "D3"),
    PCR_D4 = count(Serotype == "D4"),
    PRNT50_D1 = count(!c(is.na(`PRNT50 DENV1`))),
    PRNT50_D2 = count(!c(is.na(`PRNT50 DENV2`))),
    PRNT50_D3 = count(!c(is.na(`PRNT50 DENV3`))),
    PRNT50_D4 = count(!c(is.na(`PRNT50 DENV4`)))
  )
summary_table_long <- summary_table %>%
  pivot_longer(cols = -OnsetYear, names_to = "Variable", values_to = "Value") %>%
  pivot_wider(names_from = OnsetYear, values_from = Value)

table(DENV_cleaned_comp$`MAC-ELISA- Final result`, DENV_cleaned_comp$`GAC-ELISA-Final result`)



## Group summary ##########################################################################################
## Summarize the number of individuals for each scenarios in differentiating current and past infections
##############################################################################################################
group_summary <- DENV_cleaned_comp %>%
  mutate(Category = case_when(
    PCR_result == "positive" & PRNT_classification %in% c("None", "1 homo") ~ 
      "G1: Positive PCR, No past infections",
    PCR_result == "positive" & PRNT_classification %in% c("1 heter", "2+ heter", "2+ homo") ~ 
      "G2: Positive PCR with past infections",
    PCR_result == "negative" & `MAC-ELISA- Final result` %in% c("Negative", "Equivocal") & 
      `GAC-ELISA-Final result` == "Negative" ~
      "G3: Negative PCR, Negative/Equivocal IgM, Negative IgG, no PRNT50",
    PCR_result == "negative" & `MAC-ELISA- Final result` == "Negative" & 
      `GAC-ELISA-Final result` %in% c("Positive", "Equivocal") ~ 
      "G4: Negative PCR, Negative IgM, Positive/Equivocal IgG",
    PCR_result == "negative" & `MAC-ELISA- Final result` == "Positive" & PRNT_classification == "None" &

      `GAC-ELISA-Final result`  %in% c("Positive", "Equivocal") ~ 
      "G5: Negative PCR, Positive IgM, Positive/Equivocal IgG, no PRNT50 (excluding?)",
    PCR_result == "negative" & `MAC-ELISA- Final result`  %in% c("Positive", "Equivocal") & 
      PRNT_classification == "1 heter" & `GAC-ELISA-Final result` %in% c("Positive", "Equivocal") ~ 
      "G6: Negative PCR, Positive/Equivocal IgM, Positive/Equivocal IgG, 1PRNT (Past infection)",
    PCR_result == "negative" & `MAC-ELISA- Final result` %in% c("Positive", "Equivocal") & 
      PRNT_classification == "None" & `GAC-ELISA-Final result` == "Negative" ~ 
      "G7: Negative PCR, Positive IgM, Negative IgG, no PRNT50 (excluding?)",
    PCR_result == "negative" & `MAC-ELISA- Final result` == "Positive" & PRNT_classification == "1 heter" &
      `GAC-ELISA-Final result` == "Negative" ~ 
      "G8: Negative PCR, Positive IgM, Negative IgG, 1PRNT (Current infection)",
    PCR_result == "negative" & `MAC-ELISA- Final result` %in% c("Positive", "Equivocal") & 
      PRNT_classification == "2+ heter" ~ "G9: Cannot differentiate",
    TRUE ~ "Other"
  )) %>%
  group_by(Category) %>%
  summarise(Count = n())

## Differentiate current and pat infections###################################################################
## Differentiate between current and past infections based on the standard for each group
##############################################################################################################

## For PCR positive cases (G1 and G2), the current infection is the PCR positive serotype (confirmed cases)
## and the rest are categorized as past infection
DENV_cleaned_comp <- DENV_cleaned_comp %>%
  mutate(
    Current = if_else(
      PCR_result == "positive",
      case_when(
        PRNT_classification %in% c("None","1 homo", "2+ homo", "2+ heter", "1 heter") ~ Serotype, 
        TRUE ~ "Neg" 
      ),
      NA_character_  
    ),
    Past = if_else(
      PCR_result == "positive",
      case_when(
        PRNT_classification %in% c("None", "1 homo") ~ "Neg",
        PRNT_classification %in% c("1 heter", "2+ heter") ~ PRNT50, 
        
        # PRNT50 column excludes infected serotype if "2+ homo"
        PRNT_classification == "2+ homo" ~ str_replace_all(
          PRNT50, 
          paste0("(^|,)", Serotype, "(,|$)"), 
          # replace the current infected serotype in PRNT50 column with ","
          function(x) ifelse(nchar(x) > 1, ",", "")  
        ),  
        TRUE ~ "Neg"  
      ),
      NA_character_  
    )
  ) 

## G3: Negative PCR, Negative/Equivocal IgM, Negative IgG, no PRNT50
# Treat IgG and IgM both negative samples as PRNT50 negative 
DENV_cleaned_comp <- DENV_cleaned_comp %>%
  mutate(
    Current = if_else(
      `MAC-ELISA- Final result` %in% c("Negative", "Equivocal") & 
        PCR_result == "negative"  & `GAC-ELISA-Final result` == "Negative",
      "Neg",
      Current
    ),
    Past = if_else(
      `MAC-ELISA- Final result` %in% c("Negative", "Equivocal") & 
        PCR_result == "negative"  & `GAC-ELISA-Final result` %in% c("Negative", "Equivocal"),
      "Neg",
      Past
    )
  )

## G4: Negative PCR, Negative IgM, Positive/Equivocal IgG
DENV_cleaned_comp <- DENV_cleaned_comp %>%
  mutate(
    Current = if_else(
      `MAC-ELISA- Final result` == "Negative" & PCR_result == "negative"  &
        `GAC-ELISA-Final result` %in% c("Positive", "Equivocal"),
      "Neg",
      Current
    ),
    Past = if_else(
      `MAC-ELISA- Final result` == "Negative" & PCR_result == "negative"  & 
        `GAC-ELISA-Final result`  %in% c("Positive", "Equivocal"),
      PRNT50,
      Past
    )
  )

# G5: Negative PCR, Positive IgM, Positive/Equivocal IgG, no PRNT50 (exclude from analysis)
DENV_cleaned_comp <- DENV_cleaned_comp %>%
  filter(!(`MAC-ELISA- Final result` == "Positive" &
             PCR_result == "negative"  &
             `GAC-ELISA-Final result` %in% c("Positive" ,"Equivocal") &
             PRNT_classification == "None"))

# G6: Negative PCR, Positive IgM, Positive/Equivocal IgG, 1PRNT (Past infection)
# G8: Negative PCR, Positive IgM, Negative IgG, 1PRNT (Current infection)
DENV_cleaned_comp <- DENV_cleaned_comp %>% 
  mutate(
    Current = if_else(
      `MAC-ELISA- Final result` == "Positive" & PCR_result == "negative",
      case_when(
        PRNT_classification == "1 heter" & `GAC-ELISA-Final result` == "Negative" ~ PRNT50, # S8
        PRNT_classification == "1 heter" & 
          `GAC-ELISA-Final result` %in% c("Positive" ,"Equivocal") ~ "Neg", #S6
        TRUE ~ NA_character_
      ),
      Current  # Keep existing conditions if not applicable
    ),
    Past = if_else(
      `MAC-ELISA- Final result` == "Positive" & PCR_result == "negative",
      case_when(
        PRNT_classification == "1 heter" & 
          `GAC-ELISA-Final result` %in% c("Positive" ,"Equivocal") ~ PRNT50, #S6
        PRNT_classification == "1 heter" & `GAC-ELISA-Final result` == "Negative" ~ "Neg", #S8
        TRUE ~ NA_character_
      ),
      Past  
    )
  )

# G7: Negative PCR, Positive/Equivocal IgM, Negative IgG, no PRNT50
DENV_cleaned_comp <- DENV_cleaned_comp %>%
  mutate(
    Current = if_else(
      PCR_result == "negative" & `MAC-ELISA- Final result` %in% c("Positive", "Equivocal") &
        `GAC-ELISA-Final result` == "Negative" & PRNT_classification == "None",
      "Neg",
      Current
    ),
    Past = if_else(
      PCR_result == "negative" & `MAC-ELISA- Final result` %in% c("Positive", "Equivocal") &
        `GAC-ELISA-Final result` == "Negative" & PRNT_classification == "None",
      "Neg",
      Past
    )
  )

# Negative PCR, Equivocal IgM, Positive IgG, 1 heter (past)
DENV_cleaned_comp <- DENV_cleaned_comp %>%
  mutate(
    Current = if_else(
      PCR_result == "negative" & `MAC-ELISA- Final result` == "Equivocal" &
        `GAC-ELISA-Final result` == "Positive" & PRNT_classification == "1 heter",
      "Neg",
      Current
    ),
    Past = if_else(
      PCR_result == "negative" & `MAC-ELISA- Final result` == "Equivocal" &
        `GAC-ELISA-Final result` == "Positive" & PRNT_classification == "1 heter",
      PRNT50,
      Past
    )
  )


############################################################
## Section 4: Data Preparation for multinomial regression###
############################################################

## Function: update_factors ##################################################################################
## Function to update and relevel the categorical data for regression
## Input: 
##   data - The dataset containing columns to be updated as factors
## Output:
##   data - The dataset with updated factor levels with "Neg" being the reference category
##############################################################################################################
update_factors <- function(data) {
  
  data$Current <- factor(data$Current, levels = c("Neg", "D1", "D2", "D3", "D4"))
  data$Gender <- factor(data$Gender)
  data$Past_PRNT50_D1 <- factor(data$Past_PRNT50_D1)
  data$Past_PRNT50_D2 <- factor(data$Past_PRNT50_D2)
  data$Past_PRNT50_D3 <- factor(data$Past_PRNT50_D3)
  data$Past_PRNT50_D4 <- factor(data$Past_PRNT50_D4)
  data$OnsetYear <- factor(data$OnsetYear)
  
  data$Past_PRNT50_D1 <- relevel(data$Past_PRNT50_D1, "Neg")
  data$Past_PRNT50_D2 <- relevel(data$Past_PRNT50_D2, "Neg")
  data$Past_PRNT50_D3 <- relevel(data$Past_PRNT50_D3, "Neg")
  data$Past_PRNT50_D4 <- relevel(data$Past_PRNT50_D4, "Neg")
  
  return(data)
}


## Function: process_PRNT_data################################################################################
## Function to categorize PRNT50 levels based on the serotypes in the isolated past columns
## Input: 
##   df (data.frame) - The dataset after differentiating between past and current infections
##  threshold_low (numeric) - The lower threshold value for categorization (cases v.s non-cases)
##   threshold_high (numeric) - The upper threshold value for categorization (low v.s high antibody level)
## Output:
##   df (data.frame) - The dataset with categorized past PRNT50 levels
##############################################################################################################
process_PRNT_data <- function(df, threshold_low, threshold_high) {
  # Convert PRNT50 columns to numeric
  df$`PRNT50 DENV1` <- as.numeric(df$`PRNT50 DENV1`)
  df$`PRNT50 DENV2` <- as.numeric(df$`PRNT50 DENV2`)
  df$`PRNT50 DENV3` <- as.numeric(df$`PRNT50 DENV3`)
  df$`PRNT50 DENV4` <- as.numeric(df$`PRNT50 DENV4`)
  
  # remove all the comma in the Past column
  df$Past <- str_trim(gsub(",", "", df$Past))
  
  # Categorize PRNT50 values into "none", "low" and "high" levels based on the thresholds
  df <- df %>%
    mutate(
      PRNT50_D1 = case_when(
        is.na(`PRNT50 DENV1`) | `PRNT50 DENV1` == " " | `PRNT50 DENV1` <= threshold_low ~ "none",
        (`PRNT50 DENV1` <= threshold_high) & `PRNT50 DENV1` > threshold_low ~ "low",
        `PRNT50 DENV1` > threshold_high ~ "high"
      ),
      PRNT50_D2 = case_when(
        is.na(`PRNT50 DENV2`) | `PRNT50 DENV2` == " " | `PRNT50 DENV2` <= threshold_low ~ "none",
        (`PRNT50 DENV2` <= threshold_high) & `PRNT50 DENV2` > threshold_low ~ "low",
        `PRNT50 DENV2` > threshold_high ~ "high"
      ),
      PRNT50_D3 = case_when(
        is.na(`PRNT50 DENV3`) | `PRNT50 DENV3` == " " | `PRNT50 DENV3` <= threshold_low ~ "none",
        (`PRNT50 DENV3` <= threshold_high) & `PRNT50 DENV3` > threshold_low ~ "low",
        `PRNT50 DENV3` > threshold_high ~ "high"
      ),
      PRNT50_D4 = case_when(
        is.na(`PRNT50 DENV4`) | `PRNT50 DENV4` == " " | `PRNT50 DENV4` <= threshold_low ~ "none",
        (`PRNT50 DENV4` <= threshold_high) & `PRNT50 DENV4` > threshold_low ~ "low",
        `PRNT50 DENV4` > threshold_high ~ "high"
      )
    )
  
  # Assign pre-existing PRNT50 levels if the serotype is in Past column, 
  # "Neg" indicates no antibody level against that serotype
  df <- df %>%
    mutate(
      Past_PRNT50_D1 = case_when(
        str_detect(Past, "D1") ~ PRNT50_D1,
        TRUE ~ "Neg"
      ),
      Past_PRNT50_D2 = case_when(
        str_detect(Past, "D2") ~ PRNT50_D2,
        TRUE ~  "Neg"
      ),
      Past_PRNT50_D3 = case_when(
        str_detect(Past, "D3") ~ PRNT50_D3,
        TRUE ~  "Neg"
      ),
      Past_PRNT50_D4 = case_when(
        str_detect(Past, "D4") ~ PRNT50_D4,
        TRUE ~  "Neg"
      )
    )
  
  return(df)
}
##############################################################################################################


##############################################################################################################
## Analysis 1: PCR+ cases v.s PCR- & IgM
## The AIC is much higher than it resulted from the original data due to the lack of inherent structure
## Remember to re-run the code from beginning (right after the importation of data) if changing the thresholds
##############################################################################################################


##############################################################################################################
## Regression for analysis 1
##############################################################################################################
## Filter out the data
DENV_analysis1 <- DENV_cleaned_comp %>%
  filter(PCR_result == "positive" | (PCR_result == "negative" & `MAC-ELISA- Final result` == "Negative"))


## Model numbering corresponds to S2.2 Method-table: Multinomial logistic regression models explored

## Model 1
## For Model 1 only: we tried to use continuous PRNT50 values instead of categorical
# process_PRNT_continuous_data <- function(df, threshold_low, threshold_high) {
#   # Convert PRNT50 columns to numeric
#   df$`PRNT50 DENV1` <- as.numeric(df$`PRNT50 DENV1`)
#   df$`PRNT50 DENV2` <- as.numeric(df$`PRNT50 DENV2`)
#   df$`PRNT50 DENV3` <- as.numeric(df$`PRNT50 DENV3`)
#   df$`PRNT50 DENV4` <- as.numeric(df$`PRNT50 DENV4`)
#   
#   # remove all the comma in the Past column
#   df$Past <- str_trim(gsub(",", "", df$Past))
#   
#   
#   # Assign pre-existing PRNT50 values on log2 scale if the serotype is in Past column, 
#   # 0 indicates no antibody level against that serotype
#   df <- df %>%
#     mutate(
#       Past_PRNT50_D1_value = case_when(
#         str_detect(Past, "D1") ~ log2(`PRNT50 DENV1`),
#         TRUE ~ 0
#       ),
#       Past_PRNT50_D2_value = case_when(
#         str_detect(Past, "D2") ~ log2(`PRNT50 DENV2`),
#         TRUE ~ 0
#       ),
#       Past_PRNT50_D3_value = case_when(
#         str_detect(Past, "D3") ~ log2(`PRNT50 DENV3`),
#         TRUE ~  0
#       ),
#       Past_PRNT50_D4_value = case_when(
#         str_detect(Past, "D4") ~ log2(`PRNT50 DENV4`),
#         TRUE ~  0
#       )
#     )
#   
#   return(df)
# }
# 
# DENV_analysis1_low_high <- process_PRNT_continuous_data(DENV_analysis1, threshold_low, threshold_high)
# 
# multi_model1_low_high <- multinom(Current ~ Past_PRNT50_D1_value + Past_PRNT50_D2_value + Past_PRNT50_D3_value + 
#                                     Past_PRNT50_D4_value + Age_clean + Gender + OnsetYear , 
#                                   data = DENV_analysis1_low_high, maxit = 2000, Hess = TRUE)

## Data preparation for Model 2-6
## Process the PRNT data and convert factors
## The lower threshold should be the same to the previous threshold when doing data cleaning
DENV_analysis1_low_high <- process_PRNT_data(DENV_analysis1, threshold_low, threshold_high)
## convert the data column to categorical for regression
DENV_analysis1_low_high <- update_factors(DENV_analysis1_low_high)

## Model 2
# multi_model1_low_high <- multinom(Current ~ Past_PRNT50_D1 + Past_PRNT50_D2 + Past_PRNT50_D3 +
#                                            Past_PRNT50_D4  + Gender + Age_clean + OnsetYear,
#                                          data = DENV_analysis1_low_high, maxit = 200, Hess = TRUE)

## Model 3
# multi_model1_low_high <- multinom(Current ~ Past_PRNT50_D1 + Past_PRNT50_D2 + Past_PRNT50_D3 +
#                                            Past_PRNT50_D4  + Age_clean + OnsetYear,
#                                          data = DENV_analysis1_low_high, maxit = 200, Hess = TRUE)

## Model 4
# multi_model1_low_high <- multinom(Current ~ Past_PRNT50_D1 + Past_PRNT50_D2 + Past_PRNT50_D3 +
#                                            Past_PRNT50_D4 + Age_clean,
#                                          data = DENV_analysis1_low_high, maxit = 200, Hess = TRUE)

## Model 5
# multi_model1_low_high <- multinom(Current ~ Past_PRNT50_D1 + Past_PRNT50_D2 + Past_PRNT50_D3 + 
#                                            Past_PRNT50_D4  + Gender + poly(Age_clean, 2) + OnsetYear,
#                                          data = DENV_analysis1_low_high, maxit = 2000, Hess = TRUE)

## Model 6 Best model used in our study
multi_model1_low_high <- multinom(Current ~ Past_PRNT50_D1 + Past_PRNT50_D2 + Past_PRNT50_D3 + 
                                    Past_PRNT50_D4 + poly(Age_clean, 2)  + OnsetYear, 
                                  data = DENV_analysis1_low_high, maxit = 2000, Hess = TRUE)
summary(multi_model1_low_high)
estimates <- coef(multi_model1_low_high) # for bootstrapping


##############################################################################################################
## Data Set 1 excluding PCR positive serotype with a notable antibody titer (corresponding to Condition 1(e)
##############################################################################################################
## Filter out those people
DENV_PCR_Hprnt50 <- DENV_analysis1 %>%
  filter((Current == "D1" & `PRNT50 DENV1` > threshold_high) |
           (Current == "D2" & `PRNT50 DENV2` > threshold_high) |
           (Current == "D3" & `PRNT50 DENV3` > threshold_high) |
           (Current == "D4" & `PRNT50 DENV4` > threshold_high))

# Classify infections
DENV_PCR_Hprnt50 <- DENV_PCR_Hprnt50 %>%
  mutate(`Infected Type` = ifelse(Past == 'Neg', 'Primary', 'Subsequent'))

# Count the number of primary and subsequent infections by serotype
count_pri_sub_data <- DENV_PCR_Hprnt50 %>%
  group_by(Current, `Infected Type`) %>%
  summarise(Count_pri_sub = n()) %>%
  ungroup()

# Plotting the bar chart
ggplot(count_pri_sub_data, aes(x = Current, y = Count_pri_sub, fill = `Infected Type`)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "Serotype", y = "Number of Cases", 
       title = "Number of Primary and Subsequent Infections by Serotype") +
  geom_text(aes(label = Count_pri_sub), position = position_stack(vjust = 0.5), color = "red") +
  theme_minimal() +
  scale_fill_manual(values = c("Primary" = "#fbb4ae", "Subsequent" = "#b3cde3"))
##############################################################################################################
## Exclude those PCR positive but have natable PRNT50 cases
DENV_analysis1e <- DENV_analysis1 %>%
  anti_join(DENV_PCR_Hprnt50, by = "ID-NIHE")

## Process the PRNT data and convert factors
## The lower threshold should be the same to the previous threshold when doing data cleaning
DENV_analysis1e_low_high <- process_PRNT_data(DENV_analysis1e, threshold_low, threshold_high)
## convert the data column to categorical for regression
DENV_analysis1e_low_high <- update_factors(DENV_analysis1e_low_high)

## Model 
multi_model1e_low_high <- multinom(Current ~ Past_PRNT50_D1 + Past_PRNT50_D2 + Past_PRNT50_D3 + 
                                           Past_PRNT50_D4 + poly(Age_clean, 2)  + OnsetYear, 
                                         data = DENV_analysis1e_low_high, maxit = 2000, Hess = TRUE)
summary(multi_model1e_low_high)


##############################################################################################################
## Visualize the predicted age distribution
##############################################################################################################
# Predict probabilities
predicted_probs <- predict(multi_model1_low_high, newdata = DENV_analysis1_low_high, type = "probs")

# Add corresponding predicted probability to the data set
DENV_PCR_pos_low_high <- DENV_analysis1_low_high %>%
  mutate(predicted_D1 = predicted_probs[, "D1"],
         predicted_D2 = predicted_probs[, "D2"],
         predicted_D3 = predicted_probs[, "D3"],
         predicted_D4 = predicted_probs[, "D4"])

# Convert the data into long format for plotting
DENV_PCR_pos_low_high_long <- DENV_PCR_pos_low_high %>%
  select(Age_clean, starts_with("predicted_")) %>%
  pivot_longer(cols = starts_with("predicted_"), names_to = "Serotype", values_to = "Probability") %>%
  mutate(Serotype = case_when(
    str_detect(Serotype, "predicted_D1") ~ "D1",
    str_detect(Serotype, "predicted_D2") ~ "D2",
    str_detect(Serotype, "predicted_D3") ~ "D3",
    str_detect(Serotype, "predicted_D4") ~ "D4"
  ))

# Plot the predicted probabilities by age
ggplot(DENV_PCR_pos_low_high_long, aes(x = Age_clean, y = Probability)) +
  geom_line(aes(color = Serotype)) +
  facet_wrap(~ Serotype, scales = "fixed") +
  labs(title = "Predicted Probabilities of Dengue Serotypes by Age",
       x = "Age",
       y = "Predicted Probability") +
  theme_minimal() +
  theme(legend.position = "none")


##############################################################################################################
## Regression for analysis 2
##############################################################################################################
## Filter out the data
DENV_analysis2 <- DENV_cleaned_comp %>%
  filter(!is.na(Current))
## Process the PRNT data and convert factors
## The lower threshold should be the same to the previous threshold when doing data cleaning
DENV_analysis2_low_high <- process_PRNT_data(DENV_analysis2, threshold_low, threshold_high)
## convert the data column to categorical for regression
DENV_analysis2_low_high <- update_factors(DENV_analysis2_low_high)

## Fit the model
multi_model2_low_high <- multinom(Current ~ Past_PRNT50_D1 + Past_PRNT50_D2 + Past_PRNT50_D3 + 
                                           Past_PRNT50_D4 + poly(Age_clean, 2) + OnsetYear, 
                                         data = DENV_analysis2_low_high, maxit = 200, Hess = TRUE)

summary(multi_model2_low_high)
estimates <- coef(multi_model2_low_high) # for bootstrapping
##############################################################################################################



#####################################################
## Analysis 3: Use all the data with criteria G9.1###
#####################################################

## Prepare data for regression################################################################################
## Define the current infection when PCR- & IgM+ & 2+ serotypes with the one matches the threshold the most
##############################################################################################################
## Find the distribution when PCR+ serotype with a PRNT50 value
high_prnt50_D1 <- DENV_cleaned_comp %>%
  filter((Current == "D1" & `PRNT50 DENV1` > threshold_low)) 

high_prnt50_D2 <- DENV_cleaned_comp %>%
  filter((Current == "D2" & `PRNT50 DENV2` > threshold_low)) 

## Only 1 observation which is too little, so we will use the same threshold as DENV2
## since they are antigenically closer
high_prnt50_D4 <- DENV_cleaned_comp %>%
  filter((Current == "D4" & `PRNT50 DENV4` > threshold_low))

## Visualise the distribution
ggplot(high_prnt50_D1, aes(x = log2(`PRNT50 DENV1`))) +
  geom_histogram(binwidth = 0.5, fill = "blue") +
  labs(title = "Distribution of PRNT50 DENV1", x = "PRNT50 DENV1", y = "Frequency") +
  theme_minimal()

ggplot(high_prnt50_D2, aes(x = log2(`PRNT50 DENV2`))) +
  geom_histogram(binwidth = 0.5, fill = "blue") +
  labs(title = "Distribution of PRNT50 DENV2", x = "PRNT50 DENV2", y = "Frequency") 

##############################################################################################################
## Use the median of the PRNT50 distribution for serotype 1 and 2 as the threshold 
## for isolating current infection for Scenario 9 (S9). 
current_threshold_D1 <- summary(high_prnt50_D1$`PRNT50 DENV1`)["Median"]
current_threshold_D24 <- summary(high_prnt50_D2$`PRNT50 DENV2`)["Median"]

# Due to the samll sample size in the data for DENV4 infection, we used the same thresholds for DENV2 and DENV4,
# since they are antigentically similar. The prevalance of DENV3 infections was too low, so we assumed no 
# individuals were infected by DENV3 in the imputed data.
thresholds <- c(D1 = round(log2(current_threshold_D1)), 
                D2 = round(log2(current_threshold_D24)),
                D4 = round(log2(current_threshold_D24)))

## Function: closest_serotype_func#############################################################################
## Function to determine the closest serotype to the threshold based on PRNT50 values and predefined thresholds
## It calculates the normalized differences between the PRNT50 values and the thresholds 
## and returns the serotype with the smallest difference.
##
## Inputs:
## - d1: PRNT50 value for DENV1
## - d2: PRNT50 value for DENV2
## - d3: PRNT50 value for DENV3
## - d4: PRNT50 value for DENV4
## - year: The onset year of the individual
##
## Output:
## - closest_serotype: The serotype with the closest PRNT50 value to the threshold
###############################################################################################################
closest_serotype_func <- function(d1, d2, d3, d4, year) {
  ## Calculate the difference on a log2 scale 
  differences <- c(abs(round(log2(d1)) - thresholds["D1.Median"]), 
                   abs(round(log2(d2)) - thresholds["D2.Median"]), 
                   Inf, abs(round(log2(d4)) - thresholds["D4.Median"]))
  serotypes <- c("D1", "D2", "D3", "D4")
  
  # Filter out D3 as current prevalance of D3 infection is extremely small
  differences <- differences[c(1, 2, 4)]
  serotypes <- serotypes[c(1, 2, 4)]
  
  # Find the serotype with the smallest difference
  min_diff <- min(differences, na.rm = TRUE)
  min_diff_id <- which(differences == min_diff)
  closest_serotypes <- serotypes[min_diff_id]
  
  # If multiple serotype have the same difference to the threshold, 
  # choose the one based on the prevalence in that year
  if (length(min_diff_id) > 1) {
    if (year %in% c(2020, 2021)) {
      # the order of the prevalance from the original data
      order_preference <- c("D1", "D2", "D4")  
    } else {
      order_preference <- c("D2", "D1", "D4") # 2022
    }
    min_diff_id <- which(serotypes == order_preference[1])
  }
  
  return(serotypes[min_diff_id])
}

DENV_analysis3 <- DENV_cleaned_comp
DENV_analysis3$closest_serotype <- NA

# Apply the function to assign the closest serotype
for (i in 1:nrow(DENV_analysis3)) {
  if (DENV_analysis3$PCR_result[i] == "negative" &
      DENV_analysis3$`MAC-ELISA- Final result`[i] %in% c("Positive", "Equivocal") &
      DENV_analysis3$PRNT_classification[i] == "2+ heter") {
    DENV_analysis3$closest_serotype[i] <- closest_serotype_func(
      DENV_analysis3$`PRNT50 DENV1`[i], DENV_analysis3$`PRNT50 DENV2`[i],
      DENV_analysis3$`PRNT50 DENV3`[i], DENV_analysis3$`PRNT50 DENV4`[i], DENV_analysis3$OnsetYear[i]
    )
  }
}

###############################################################################################################
## G9: Cannot differentiate, choose the serotype with the closest PRNT50 values to the threshold among D1,2,4 
## as current infection
###############################################################################################################
DENV_analysis3 <- DENV_analysis3 %>%
  mutate(
    Current = if_else(
      PCR_result == "negative" & `MAC-ELISA- Final result` %in% c("Positive", "Equivocal") & 
        PRNT_classification == "2+ heter",
      closest_serotype,
      Current
    ),
    Past = if_else(
      PCR_result == "negative" & `MAC-ELISA- Final result` %in% c("Positive", "Equivocal") & 
        PRNT_classification == "2+ heter",
      case_when(
        !is.na(Current)  ~ str_replace_all(
          PRNT50, 
          paste0("(^|,)", Current, "(,|$)"),  
          function(x) ifelse(nchar(x) > 1, ",", "") 
        ),
        TRUE ~ NA_character_
      ),
      Past
    )
  )

###############################################################################################################
## Fit to the model
DENV_analysis3_low_high <- process_PRNT_data(DENV_analysis3, threshold_low, threshold_high)
DENV_analysis3_low_high <- update_factors(DENV_analysis3_low_high)
multi_model3_low_high <- multinom(Current ~ Past_PRNT50_D1 + Past_PRNT50_D2 + Past_PRNT50_D3 + Past_PRNT50_D4 
                                    + poly(Age_clean,2) + OnsetYear, 
                                   data = DENV_analysis3_low_high, maxit = 200, Hess = TRUE)
summary(multi_model3_low_high)
estimates <- coef(multi_model3_low_high) # for bootstrapping
###############################################################################################################

##############################################################################################################
##Bootstrapping to get confidence intervals
##############################################################################################################
# Change model names for each analysis
num_coef <- length(as.vector(coef(multi_model1_low_high))) ## number of coefficients
coef_template <- coef(multi_model1_low_high) # for template
#estimates <- coef(multi_model1_low_high) # for estimates if not store in the analysis section
coef_rownames <- rownames(coef_template)
coef_colnames <- colnames(coef_template)

# Number of bootstrapping samples
num_bootstrap <- 1000

# Initialize a matrix to store bootstrap estimates
bootstrap_estimates <- matrix(NA, nrow = num_bootstrap, ncol = num_coef)
set.seed(123) 

# Bootstrapping
for (i in 1:num_bootstrap) {
  # Sample with replacement from the data, change the model names accordingly
  sample_indices <- sample(1:nrow(DENV_analysis1_low_high), replace = TRUE)
  bootstrap_sample <- DENV_analysis1_low_high[sample_indices, ]
  bootstrap_sample <- update_factors(bootstrap_sample)
  
  # Fit the multinomial regression model on the bootstrap sample
  model <- try(multinom(Current ~ Past_PRNT50_D1 + Past_PRNT50_D2 + Past_PRNT50_D3 + Past_PRNT50_D4 
                        + poly(Age_clean,2) + OnsetYear, 
                        data = bootstrap_sample, maxit = 2000, Hess = TRUE), silent = TRUE)
  # produce error message
  if (inherits(model, "try-error")) {
    cat("Model fitting failed at iteration:", i, "\n")
    next # Skip this iteration if the model didn't fit
  }
  
  
  # Initialize the coefficient matrix for the current bootstrap sample
  current_coefs <- coef_template
  current_coefs[,] <- NA
  
  # Get the coefficients from the fitted model
  fitted_coefs <- coef(model)
  
  # Fill the template with the fitted coefficients
  common_levels <- intersect(rownames(current_coefs), rownames(fitted_coefs))
  current_coefs[common_levels, ] <- fitted_coefs[common_levels, ]
  
  # Store the coefficients
  bootstrap_estimates[i, ] <- as.vector(current_coefs)
}

# Remove rows with NA values 
bootstrap_estimates <- bootstrap_estimates[complete.cases(bootstrap_estimates), ]

# Calculate the confidence intervals for each coefficient
conf_intervals <- apply(bootstrap_estimates, 2, function(x) {
  quantile(x, probs = c(0.025, 0.975))
})

## bootstrapping estimates in the format as the template
bstp_temp_estimates <- matrix(NA, nrow = length(coef_rownames), ncol = length(coef_colnames))
rownames(bstp_temp_estimates) <- coef_rownames
colnames(bstp_temp_estimates) <- coef_colnames

## Structure the estimate matrix into the format: estimates(lower bound, upper bound)
for (i in 1:length(coef_colnames)) {
  lower_upper <- conf_intervals[, (i - 1) * length(coef_rownames) + (1:length(coef_rownames))]
  for (j in 1:length(coef_rownames)) {
    estimate <- estimates[j, i]
    lower <- lower_upper[1, j]
    upper <- lower_upper[2, j]
    bstp_temp_estimates[j, i] <- paste0(round(estimate, 3), " (", round(lower, 3), ", ", round(upper, 3), ")")
  }
}
bstp_temp_estimates <- t(bstp_temp_estimates)

## Function: check_in_CI######################################################################################
## Function to check whether the estimates is within the confidence interval#
## Input: value - each value in the resulted bootstrapping matrix in the format 
## estimates(lower bound, upper bound)
## Output: TRUE if the confidence interval includes the estimates; otherwise, FALSE 
##############################################################################################################
check_in_CI <- function(value) {
  # Regular expression for extracting estimates, lower and upper bound
  matches <- regmatches(value, 
             regexec("\\s*(-?\\d+\\.\\d+)\\s*\\(\\s*(-?\\d+\\.\\d+)\\s*,\\s*(-?\\d+\\.\\d+)\\s*\\)", value))
  if (length(matches[[1]]) == 4) {
    estimate <- as.numeric(matches[[1]][2])  # estimate
    ci_lower <- as.numeric(matches[[1]][3])  # lower CI 
    ci_upper <- as.numeric(matches[[1]][4])  # upper CI 
    return(estimate >= ci_lower & estimate <= ci_upper)
  } else {
    return(NA) 
  }
}

# Apply the function to each value in the resulted bootstrapping estimates matrix
results_in_CI <- matrix(nrow = nrow(bstp_temp_estimates), ncol = ncol(bstp_temp_estimates),
                        dimnames = dimnames(bstp_temp_estimates))
for (i in 1:nrow(bstp_temp_estimates)) {
  for (j in 1:ncol(bstp_temp_estimates)) {
    results_in_CI[i, j] <- check_in_CI(bstp_temp_estimates[i, j])
  }
}
results_in_CI
##############################################################################################################

###################################################################################################
## Analysis 4: Use all the data from S1-S8 with 20% individuals from S10 using EM algorithm     ###
## Note: The generated data may not converge in this analysis since the test results are random ###
## which may lead to fail in EM algorithm.                                                      ###
###################################################################################################

## Function: log_likelihood###################################################################################
## Function to calculate the true likelihood of the EM algorithm
##
## Input:
## - params (vector): Model parameters 
## - data: A dataframe containing columns: Current (infection status), id (observation identifier), 
##         weights (observation weights), and other predictor columns.
##
## Output:
## - total_logllh: Total log likelihood value for the dataset given the parameters.
##############################################################################################################
log_likelihood <- function(params, data) {
  # make params of the hessian function a matrix
  vec_par <- c(params)
  num_predictors <- 12
  num_categories <- 4
  param_matrix <- matrix(vec_par, ncol = num_predictors + 1, nrow = num_categories, byrow = FALSE)
  
  rownames(param_matrix) <- c("D1", "D2", "D3", "D4")
  
  # Assign row names 
  colnames(param_matrix) <- c("(Intercept)", "Past_PRNT50_D1high", "Past_PRNT50_D1low",
                              "Past_PRNT50_D2high", "Past_PRNT50_D2low",
                              "Past_PRNT50_D3high", "Past_PRNT50_D3low", "Past_PRNT50_D4high", 
                              "Past_PRNT50_D4low", "poly(Age_clean, 2)1", "poly(Age_clean, 2)2",
                              "OnsetYear2021", "OnsetYear2022")
  
  # transpose and add 0s
  # Add a row of zeros at the top
  param_matrix <- rbind(rep(0, ncol(param_matrix)), param_matrix)
  # Add a column of zeros at the beginning
  param_matrix <- cbind(rep(0, nrow(param_matrix)), param_matrix)
  
  # vectorise
  new_wts <- c(t(param_matrix))
  model_result$wts <- new_wts
  
  # Extract parameters
  probabilities <- predict(model_result, newdata = data, type = "probs") ## take log
  # define the mapping of 'Current' categories to the columns in 'probabilities'
  level_mapping <- c( "Neg" = 1, "D1" = 2, "D2" = 3, "D3" = 4, "D4" = 5)
  level_indices <- sapply(data$Current, function(x) level_mapping[x])
  #  Extract the specific probability for each row based on the mapped 'Current'
  llhs <- sapply(seq_along(level_indices), function(i) probabilities[i, level_indices[i]])
  
  weights <- data$weights
  
  # Create a data frame with NIHE ID, likelihoods, and weights
  llh_data <- data.frame(ID = data$id,
                         llh = llhs,
                         weights = weights)
  
  # Group by NIHE ID and calculate the sum of individual llhs
  sum_llh <- llh_data %>%
    group_by(ID) %>%
    # small number issue: log the product and use logSumExp()
    summarise(sum_llh = -logSumExp(log(weights) + log(llh))) 
  # summarise(sum_llh = -sum(weights * log(llh)))
  
  # Calculate the total log-likelihood
  total_logllh <- sum(sum_llh$sum_llh)
  
  return(total_logllh)
}


###############################################################################################################
## Prepare for the data
###############################################################################################################
## Filter out complete data
DENV_analysis4_comp <- DENV_cleaned_comp %>%
  filter(!is.na(Current))
## Process the PRNT data and convert factors
DENV_analysis4_comp_low_high <- process_PRNT_data(DENV_analysis4_comp, threshold_low, threshold_high)
## convert the data column to categorical for regression
DENV_analysis4_comp_low_high <- update_factors(DENV_analysis4_comp_low_high)

## Filter out data for EM algorithm
DENV_analysis4_EM <- DENV_cleaned_comp %>%
  filter(is.na(Current))

## Add id to identify individuals for bootstrapping
DENV_analysis4_EM$id <- c(1:nrow(DENV_analysis4_EM))
DENV_analysis4_comp_low_high$id <- c((nrow(DENV_analysis4_EM) + 1):
                              (nrow(DENV_analysis4_EM)  + nrow(DENV_analysis4_comp_low_high)))
# Weights for complete data is 1 for each individual
DENV_analysis4_comp_low_high$weights <- 1


## Function: expand_rows######################################################################################
## Function to expand the data for imputation by splitting serotypes in PRNT50 column into multiple rows
## Assigning each serotype in PRNT50 as current infection in turn and treat the rest as past infection(s)
##
## Input:
## - df (dataframe): Data frame to be used for EM algorithm (DENV_analysis4_EM)
##
## Output:
## - Expanded dataframe where each row represents a possible infection status for an individual individual 
##############################################################################################################
expand_rows <- function(df) {
  df %>%
    # Create a temporary split version of PRNT50
    mutate(PRNT50_split = strsplit(PRNT50, ",")) %>%
    # unnest_longer() turns each element of a list-column into a row
    # Expand the data frame so that each PRNT50 element has its own row
    unnest_longer(col = PRNT50_split, names_repair = "minimal") %>%
    # Create current and past infection columns
    group_by(`ID-NIHE`) %>%
    mutate(Current = PRNT50_split,  # Current infection is the split element
           # past infections are those excluding the current infections
           Past = sapply(seq_along(PRNT50_split), 
                         function(i) paste(setdiff(unlist(PRNT50_split), PRNT50_split[i]), collapse = ","))) %>%
    ungroup() %>%
    # Remove the split column
    select(-PRNT50_split)
}
##############################################################################################################
# Randomly select 20% from the data that requires EM algorithm for btsp to avoid over imputation
set.seed(666)
DENV_analysis4_0.2EM <- DENV_analysis4_EM[sample(seq_len(nrow(DENV_analysis4_EM)), 
                                        size = 0.2 * nrow(DENV_analysis4_EM), replace = FALSE), ]

## convert the data column to categorical for regression
DENV_analysis4_0.2EM_low_high <- expand_rows(DENV_analysis4_0.2EM)
DENV_analysis4_0.2EM_low_high <- process_PRNT_data(DENV_analysis4_0.2EM_low_high, threshold_low, threshold_high)
DENV_analysis4_0.2EM_low_high <- update_factors(DENV_analysis4_0.2EM_low_high)
DENV_analysis4_0.2EM_low_high <- DENV_analysis4_0.2EM_low_high %>%
  mutate(index = row_number()) %>%
  group_by(id) %>%
  mutate(weights = 1/n()) %>%
  ungroup() %>%
  select(-index)
# initial weightes
DENV_analysis4_comp_low_high$weights <- 1
DENV_analysis4_EM_initial <- rbind(DENV_analysis4_comp_low_high, DENV_analysis4_0.2EM_low_high)

##############################################################################################################
## EM algorithm
##############################################################################################################
# To derive the initial parameters
multi_model_EM <- multinom(Current ~ Past_PRNT50_D1 + Past_PRNT50_D2 + Past_PRNT50_D3 + 
                             Past_PRNT50_D4 + poly(Age_clean,2) + OnsetYear, 
                           data = DENV_analysis4_EM_initial, weights = DENV_analysis4_EM_initial$weights, 
                           Hess = TRUE, maxit = 500)
initial_params <- coef(multi_model_EM)


# Number of iterations and tolerance for convergence
max_iter <- 100
tolerance <- 1e-6
converged <- FALSE
iteration <- 1
params_trace <- list()

while (!converged && iteration <= max_iter) {
  # E-step: Predict the probabilities of missing categories using the current model
  probabilities <- predict(multi_model_EM, newdata = DENV_analysis4_0.2EM_low_high, type = "probs")
  
  # Assign probabilities to maximize the expected log-likelihood
  DENV_analysis4_0.2EM_low_high <- DENV_analysis4_0.2EM_low_high %>%
    mutate(weights = probabilities[, 1]) %>%
    group_by(id) %>%
    mutate(weights = weights / sum(weights)) %>%
    ungroup()
  
  # Combine the complete and the imputed data
  updated_data <- rbind(DENV_analysis4_comp_low_high, DENV_analysis4_0.2EM_low_high)
  updated_data <- update_factors(updated_data)
  
  
  multi_model_EM <- multinom(Current ~ Past_PRNT50_D1 + Past_PRNT50_D2 + Past_PRNT50_D3 + 
                               Past_PRNT50_D4 + poly(Age_clean,2) + OnsetYear, 
                             weights = updated_data$weights, data = updated_data, maxit = 2000, Hess = TRUE)
  
  # Check for convergence (change in coefficients)
  new_params <- coef(multi_model_EM)
  params_trace[[iteration]] <- summary(multi_model_EM)
  
  difference <- new_params - initial_params
  if (max(abs(difference)) < tolerance) {
    converged <- TRUE
    print("EM algorithm has converged.")
  } else {
    initial_params <- new_params
    print(max(abs(difference)))
    print(paste("Iteration", iteration, "completed."))
  }
  
  iteration <- iteration + 1
}

summary(multi_model_EM) 
coef_matrix <- coef(multi_model_EM) # for bootstrapping 
##############################################################################################################
## bootstrapping method
##############################################################################################################
n_btsp <- 100 
max_iter <- 100
tolerance <- 1e-6

# Store all the parameter estimates from bootstrapping
param_estimates <- list()

for (i in 1:n_btsp) {
  converged <- FALSE
  iteration <- 1
  
  # Randomly select 20% from the data that requires EM algorithm 
  # for bootstrapping to avoid over imputation
  btsp_data <- DENV_analysis4_0.2EM[sample(nrow(DENV_analysis4_0.2EM), replace = TRUE), ]
  # add an index to represent each individual, since an individual may be select multiple times 
  # and should be treated as different people therefore the NIHE-ID is not enough to differentiate
  btsp_data$id <- c(1:nrow(btsp_data)) 
  btsp_data <- expand_rows(btsp_data)
  btsp_data <- process_PRNT_data(btsp_data, threshold_low, threshold_high)
  btsp_data <- update_factors(btsp_data)
  
  # Assign initial weights to each row, assume an individual is equally susceptible to each serotype
  btsp_data <- btsp_data %>%
    mutate(idx = row_number()) %>%
    group_by(id) %>%
    mutate(weights = 1/n()) %>%
    ungroup() %>%
    select(-idx)
  
  # Store likelihoods and parameters after 95 iterations (in case of not converge)
  llh_after_95 <- vector("list", 5)
  param_after_95 <- vector("list", 5)
  
  while (!converged && iteration <= max_iter) {
    probabilities <- predict(multi_model_EM, newdata = btsp_data, type = "probs")
    
    # Assign probabilities to maximize the expected log-likelihood
    btsp_data <- btsp_data %>%
      mutate(weights = probabilities[, 1]) %>%
      group_by(id) %>%
      mutate(weights = weights / sum(weights)) %>%
      ungroup()
    
     
    # Combine the complete data and the imputed data to form the updated data
    updated_data <- rbind(DENV_analysis4_comp_low_high, btsp_data)
    updated_data <- update_factors(updated_data)
    
    tryCatch({
      multi_model_EM <- multinom(Current ~ Past_PRNT50_D1 + Past_PRNT50_D2 + Past_PRNT50_D3 + 
                                   Past_PRNT50_D4 + poly(Age_clean,2) + OnsetYear, 
                                 data = updated_data, weights = updated_data$weights, 
                                 Hess = TRUE, maxit = 500)
      # Extract parameters
      new_params <- coef(multi_model_EM)
      model_result <- multi_model_EM
      
      # If after 95 iterations, still not converge, then store the likelihoods for comparison 
      if (iteration > 95) {
        # Store likeligoods and parameters for the last 5 iterations
        current_llh <- log_likelihood(new_params, updated_data)
        idx <- iteration - 95
        if (idx <= 5) {
          llh_after_95[[idx]] <- current_llh
          param_after_95[[idx]] <- new_params
        }
      }
      
      # Stopping criteria
      difference <- max(abs(new_params - initial_params))
      if (difference < tolerance) {
        converged <- TRUE
        print(paste("EM algorithm has converged for bootstrp", i))
      } else {
        initial_params <- new_params
      }
    }, error = function(e) {
      print(paste("Error during fit: ", e$message))
    })
    
    iteration <- iteration + 1
  }
  
  # If not converge, choose the set of parameters with the highest true likelihood
  if (!converged && any(vapply(llh_after_95, is.numeric, logical(1))) ) {
    # Find the idx of the highest likelihood
    print(unlist(llh_after_95))
    max_llh_idx <- which.max(unlist(llh_after_95))
    param_estimates[[i]] <- param_after_95[[max_llh_idx]]
  } else {
    param_estimates[[i]] <- if (converged) new_params else NA
  }
}

###############################################################################################################
## Get Confidence interval from bootstrapping
###############################################################################################################
# Calculate confidence intervals for each parameter
param_matrix <- do.call(rbind, param_estimates)
param_name <- colnames(param_matrix)  
# Initialize a matrix to store CI values
CI <- matrix(nrow = 4, ncol = length(param_name) * 2)  

row_name <- c("D1", "D2", "D3", "D4")  # Outcome categories
col_name <- vector("character", length = length(param_name) * 2)

for (param_idx in 1:ncol(param_matrix)) {
  for (outcome_idx in 1:4) {  
    indices <- seq(outcome_idx, nrow(param_matrix), by = 4)
    values <- param_matrix[indices, param_idx]
    ci_lower <- quantile(values, probs = 0.25, na.rm = TRUE)
    ci_upper <- quantile(values, probs = 0.75, na.rm = TRUE)
    
    # Fill the matrix with CI values
    CI[outcome_idx, (param_idx * 2 - 1)] <- ci_lower
    CI[outcome_idx, (param_idx * 2)] <- ci_upper
    
    # Set column names for CIs
    col_name[(param_idx * 2 - 1)] <- paste(param_name[param_idx], "lower", sep = "_")
    col_name[(param_idx * 2)] <- paste(param_name[param_idx], "upper", sep = "_")
  }
}

# Assign names to rows and columns of the matrix
rownames(CI) <- row_name
colnames(CI) <- col_name

CI_lower <- CI[, seq(1, by = 2, length.out = 13)]
CI_upper <- CI[, seq(2, by = 2, length.out = 13)]

bstp_temp_estimates <- matrix(nrow = 4, ncol = 13)
rownames(bstp_temp_estimates) <- rownames(coef_matrix)
colnames(bstp_temp_estimates) <- param_name

# Form the structured bootstrapping estimates with CI
for (j in 1:ncol(bstp_temp_estimates)) {
  bstp_temp_estimates[, j] <- paste0(
    format(coef_matrix[, j], digits = 4),  # estimates
    " (", 
    format(CI_lower[, j], digits = 4),  # 25% CI
    ", ", 
    format(CI_upper[, j], digits = 4),  # 75% CI
    ")"
  )
}

# transpose
bstp_temp_estimates <- t(bstp_temp_estimates)
# Apply the function to each cell in the results matrix
results_in_CI <- matrix(nrow = nrow(bstp_temp_estimates), ncol = ncol(bstp_temp_estimates),
                        dimnames = dimnames(bstp_temp_estimates))

# check_in_CI is defined in the previous bootstrapping method for Analysis 1-3
for (i in 1:nrow(bstp_temp_estimates)) {
  for (j in 1:ncol(bstp_temp_estimates)) {
    results_in_CI[i, j] <- check_in_CI(bstp_temp_estimates[i, j])
  }
}
results_in_CI
###############################################################################################################


###################################
## Section 5: Data Visualization###
###################################

## Basic Plot 1a###############################################################################################
## Plot 1a (bar): visualizes age distribution by serotype
## Remember to change the data set names accordingly
###############################################################################################################
# Assign age groups
DENV_age_serotype <- DENV_analysis1_low_high %>%
  mutate(Age_Group = cut(Age_clean,
                         breaks = c(-Inf, 10, 20, 30, 40, 50, 60, 70, Inf),
                         labels = c("<10", "10-20", "20-30", "30-40", "40-50", "50-60", "60-70", "70+"),
                         include.lowest = TRUE))

# Calculate the total population in each age group and the proportion for each serotype within each age group
DENV_age_serotype <- DENV_age_serotype %>%
  group_by(Age_Group) %>%
  summarise(total = n()) %>%
  left_join(
    DENV_age_serotype %>%
      group_by(Age_Group, Serotype) %>%
      summarise(count = n()),
    by = "Age_Group"
  )

color_id <- c("D1" = "#fbb4ae", "D2" = "#b3cde3", "D3" = "#ccebc5", "D4" = "#decbe4", "Neg" = "#fed9a6")

DENV_age_serotype <- DENV_age_serotype %>%
  filter(! Serotype == "Neg")

# Age Distribution by Serotype
ggplot(DENV_age_serotype, aes(x = Age_Group, y = count, fill = Serotype)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ Serotype) +
  scale_fill_manual(values = color_id) +
  labs(title = "Age Distribution by Serotype",
       x = "Age Group",
       y = "Count") +
  theme_minimal() 

## Basic Plot 1b###############################################################################################
## Plot 1b (bar): visualizes gender distribution by serotype
###############################################################################################################
# Calculate the total population in each gender group and the proportion 
# for each serotype within each gender group
DENV_gender_serotype <- DENV_analysis1_low_high %>%
  group_by(Gender) %>%
  summarise(total = n()) %>%
  left_join(
    DENV_analysis1_low_high %>%
      group_by(Gender, Serotype) %>%
      summarise(count = n()),
    by = "Gender"
  )%>%
  filter(! Serotype == "Neg")

# Gender Distribution by Serotype
ggplot(DENV_gender_serotype, aes(x = Gender, y = count, fill = Serotype)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ Serotype) +
  scale_fill_manual(values = color_id) +
  labs(title = "Gender Distribution by Serotype",
       x = "Gender",
       y = "Count") +
  theme_minimal()


## Basic Plot 2################################################################################################
## Plot 2 (line + bar): visualizes yearly dengue infection counts and percentages by serotype
###############################################################################################################
## Count of infection by serotype and year
yearly_infection_count <- DENV_analysis1_low_high %>% 
  filter(Current %in% c("D1", "D2", "D3", "D4")) %>%
  group_by(OnsetYear, Current) %>%
  summarise(Count = n())

## Percentage of infection by serotype and year
yearly_infection_perc <- DENV_analysis1_low_high %>%
  group_by(OnsetYear, Current) %>%
  summarise(count = n()) %>%
  group_by(OnsetYear) %>%
  mutate(total = sum(count),
         percentage = (count / total) * 100) %>%
  filter(!Current %in% c("Neg", "D3"))


# Define the maximums for axis
max_count <- 1500
max_percentage <- 50 

color_id <- c("D1" = "#fbb4ae", "D2" = "#b3cde3", "D3" = "#ccebc5", "D4" = "#decbe4", "Unknown" = "#fed9a6")

## Plotting
ggplot() +
  geom_bar(data = yearly_infection_count, aes(x = factor(OnsetYear), y = Count, fill = factor(Current)), 
           stat = "identity", position = "stack") +
  geom_line(data = yearly_infection_perc, aes(x = factor(OnsetYear), 
            y = (percentage / max_percentage) * max_count, group = Current, color = Current), 
            linewidth = 1) +
  scale_y_continuous(
    name = "Count",
    limits = c(0, max_count), 
    # Adjust the primary y-axis (count) values to the scale of the secondary y-axis (percentage)
    sec.axis = sec_axis(~ . / max_count * max_percentage, name = "Percentage", breaks = seq(0, 50, 10))
  ) +
  scale_fill_manual(values = color_id) +
  scale_color_manual(values = c("lightcoral", "cyan2", "mediumpurple1")) + 
  labs(
    title = "Current Infection by Year",
    x = "Year",
    fill = "Serotype",
    color = "Serotype"
  ) +
  theme_minimal() 

## Basic Plot 3################################################################################################
## Plot 3 (violin): visualizes the distribution of PRNT50 values by serotype
###############################################################################################################
# Get unique PRNT values
PRNT_class <- unique(na.omit(c(DENV_analysis1_low_high$`PRNT50 DENV1`, 
                               DENV_analysis1_low_high$`PRNT50 DENV2`,
                               DENV_analysis1_low_high$`PRNT50 DENV3`,
                               DENV_analysis1_low_high$`PRNT50 DENV4`)))

# Data to produce violin plots
# Replace NA with 0 
violin_data1 <- DENV_analysis1_low_high %>%
  #select(-`PRNT50 DENV1-VN`) %>%  # exist in the orginial data, not in the fake one
  filter(Current %in% c("D1", "D2", "D3", "D4", "Neg")) %>% 
  # avoid the conflict of the function with the same name
  dplyr::count(Current, `PRNT50 DENV1` = `PRNT50 DENV1`) %>% 
  mutate(`PRNT50 DENV1` = replace_na(`PRNT50 DENV1`, 0)) %>%
  mutate(DENV_type = "anti-DENV1")

violin_data2 <- DENV_analysis1_low_high %>%
  #select(-`PRNT50 DENV1-VN`) %>%
  filter(Current %in% c("D1", "D2", "D3", "D4", "Neg")) %>% 
  dplyr::count(Current, `PRNT50 DENV2` = `PRNT50 DENV2`) %>%
  mutate(`PRNT50 DENV2` = replace_na(`PRNT50 DENV2`, 0)) %>%
  mutate(DENV_type = "anti-DENV2")

violin_data3 <- DENV_analysis1_low_high %>%
  #select(-`PRNT50 DENV1-VN`) %>%
  filter(Current %in% c("D1", "D2", "D3", "D4", "Neg")) %>% 
  dplyr::count(Current, `PRNT50 DENV3` = `PRNT50 DENV3`) %>%
  mutate(`PRNT50 DENV3` = replace_na(`PRNT50 DENV3`, 0)) %>%
  mutate(DENV_type = "anti-DENV3")

violin_data4 <- DENV_analysis1_low_high %>%
  #select(-`PRNT50 DENV1-VN`) %>%
  filter(Current %in% c("D1", "D2", "D3", "D4", "Neg")) %>% 
  dplyr::count(Current, `PRNT50 DENV4` = `PRNT50 DENV4`) %>%
  mutate(`PRNT50 DENV4` = replace_na(`PRNT50 DENV4`, 0)) %>%
  mutate(DENV_type = "anti-DENV4")

## Combine the violin data 1-4 for plotting
violin_data <- bind_rows(
  violin_data1 %>% rename(PRNT50 = `PRNT50 DENV1`),
  violin_data2 %>% rename(PRNT50 = `PRNT50 DENV2`),
  violin_data3 %>% rename(PRNT50 = `PRNT50 DENV3`),
  violin_data4 %>% rename(PRNT50 = `PRNT50 DENV4`)
) 

## Filter out the D3 infection which only occurred in 2020 (not enough for plotting)
violin_data <- violin_data %>%
  filter(!is.na(Current)) %>%
  filter(!Current == "D3") %>%
  ## Convert 0 values to log2 scale: log2(1) = 0
  mutate(PRNT50 = as.numeric(as.character(PRNT50)),
         PRNT50 = ifelse(PRNT50 == 0, 1, PRNT50)) 

color_id <- c("D1" = "#fbb4ae", "D2" = "#b3cde3", "D3" = "#ccebc5", "D4" = "#decbe4", "Neg" = "#fed9a6")

## Plotting
violin_data %>%
  uncount(weights = n, .id = "id") %>%
  ggplot(aes(x = Current, y = log2(PRNT50), fill = Current)) +
  geom_violin() +
  facet_wrap(~DENV_type, scales = "fixed") +
  labs(title = "Distribution of PRNT50 Levels by Infection Type",
       x = "Dengue Serotype",
       y = "Log2(PRNT50)") +
  theme_minimal() +
  scale_fill_manual(values = color_id) 
###############################################################################################################


## Basic Plot 4################################################################################################
## Plot 4 (forest): visualizes the odds ratios with confidence intervals
###############################################################################################################
# Remove the comma in Past column
forest <- DENV_analysis1_low_high %>%
  mutate(Past = strsplit(as.character(Past), ",")) 
 
# Extract the PRNT50 values for past infections 
forest <- forest %>%
  mutate(
    Past_PRNT50_D1_values = case_when(
      str_detect(Past, "D1") ~ `PRNT50 DENV1`,
      TRUE ~ 0
    ),
    Past_PRNT50_D2_values = case_when(
      str_detect(Past, "D2") ~ `PRNT50 DENV2`,
      TRUE ~ 0
    ),
    Past_PRNT50_D3_values = case_when(
      str_detect(Past, "D3") ~ `PRNT50 DENV3`,
      TRUE ~ 0
    ),
    Past_PRNT50_D4_values = case_when(
      str_detect(Past, "D4") ~ `PRNT50 DENV4`,
      TRUE ~ 0
    )
  )

# An empty list to store forest data
forest_data <- list()

serotypes <- c("D1", "D2", "D3", "D4")

## Count the number of cases and non-cases based on thresholds for those non-immune to "current" serotype
for (current in serotypes) {
  temp_data <- list()
  for (other in serotypes) {
    if (current != other) {
      # Filter to include only non-immunes to 'current'
      non_immunes <- forest %>%
        filter(!!sym(paste0("Past_PRNT50_", current, "_values")) == 0)
      
      ## Calculate counts for cases and non-cases
      # less than lower threshold
      cases_le_low <- sum(non_immunes$Current == current & 
                            non_immunes[[paste0("Past_PRNT50_", other, "_values")]] <= threshold_low)
      non_cases_le_low <- sum(non_immunes$Current == "Neg" & 
                                non_immunes[[paste0("Past_PRNT50_", other, "_values")]] <= threshold_low)
      cases_gt_low_le_high <- sum(non_immunes$Current == current & 
                                    non_immunes[[paste0("Past_PRNT50_", other, "_values")]] > threshold_low & 
                                    non_immunes[[paste0("Past_PRNT50_", other, "_values")]] <= threshold_high)
      # greater than lower and less than upper threshold
      non_cases_gt_low_le_high <- sum(non_immunes$Current == "Neg" & 
                                    non_immunes[[paste0("Past_PRNT50_", other, "_values")]] > threshold_low & 
                                    non_immunes[[paste0("Past_PRNT50_", other, "_values")]] <= threshold_high)
      # greater than upper threshold
      cases_gt_high <- sum(non_immunes$Current == current & 
                           non_immunes[[paste0("Past_PRNT50_", other, "_values")]] > threshold_high)
      non_cases_gt_low <- sum(non_immunes$Current == "Neg" & 
                              non_immunes[[paste0("Past_PRNT50_", other, "_values")]] > threshold_high)
      
      # Store the counts in a list
      temp_data[[other]] <- c(cases_le_low, non_cases_le_low, cases_gt_low_le_high, 
                              non_cases_gt_low_le_high, cases_gt_high, non_cases_gt_low)
    }
  }
  forest_data[[current]] <- temp_data
}

## Initialize the data frame template
forest_df <- data.frame(
  Current = character(),
  Anti_DENV = character(),
  PRNT50_le_low_Cases = integer(),
  PRNT50_le_low_NonCases = integer(),
  PRNT50_gt_low_le_high_Cases = integer(),
  PRNT50_gt_low_le_high_NonCases = integer(),
  PRNT50_gt_high_Cases = integer(),
  PRNT50_gt_high_NonCases = integer(),
  stringsAsFactors = FALSE
)

# Convert list into data frame
for (current in serotypes) {
  for (other in serotypes) {
    if (current != other) {
      temp_data <- unlist(forest_data[[current]][[other]])
      new_row <- data.frame(
        Current = current,
        Anti_DENV = other,
        PRNT50_le_low_Cases = temp_data[1],
        PRNT50_le_low_NonCases = temp_data[2],
        PRNT50_gt_low_le_high_Cases = temp_data[3],
        PRNT50_gt_low_le_high_NonCases = temp_data[4],
        PRNT50_gt_high_Cases = temp_data[5],
        PRNT50_gt_high_NonCases = temp_data[6]
      )
      forest_df <- rbind(forest_df, new_row)
    }
  }
}

## Calculate Odds Ratios and Confidence Intervals
forest_df <- forest_df %>%
  mutate(
    # Calculate Odds Ratio (reference level is PRNT50 values less than threshold, i.e "None")
    Comparison = paste("anti-", Anti_DENV, "titers in non immunes to", Current),
    
    ## Odds ratio between a low antibody level and reference level
    OR_low = (PRNT50_gt_low_le_high_Cases / PRNT50_gt_low_le_high_NonCases) /
      (PRNT50_le_low_Cases / PRNT50_le_low_NonCases),
    OR_low = ifelse(PRNT50_le_low_Cases == 0 | PRNT50_le_low_NonCases == 0 | 
                    PRNT50_gt_low_le_high_Cases == 0 | PRNT50_gt_low_le_high_NonCases == 0,
                    NA,  
                    OR_low),
    ## Odds ratio between a high antibody level and reference level
    OR_high = (PRNT50_gt_high_Cases / PRNT50_gt_high_NonCases) / 
      (PRNT50_le_low_Cases / PRNT50_le_low_NonCases),
    OR_high = ifelse(PRNT50_le_low_Cases == 0 | PRNT50_le_low_NonCases == 0 | 
                     PRNT50_gt_high_Cases == 0 | PRNT50_gt_high_NonCases == 0,
                     NA,  
                     OR_high),
    
    # Calculate lower bound of 95% CI for low antibody level
    CI_lower_low = exp(log(OR_low) - 1.96 * sqrt(1/PRNT50_le_low_Cases + 1/PRNT50_le_low_NonCases + 
                                            1/PRNT50_gt_low_le_high_Cases + 1/PRNT50_gt_low_le_high_NonCases)),
    # Calculate upper bound of 95% CI for low antibody level
    CI_upper_low = exp(log(OR_low) + 1.96 * sqrt(1/PRNT50_le_low_Cases + 1/PRNT50_le_low_NonCases + 
                                            1/PRNT50_gt_low_le_high_Cases + 1/PRNT50_gt_low_le_high_NonCases)),
    # Calculate lower bound of 95% CI for high antibody level
    CI_lower_high = exp(log(OR_high) - 1.96 * sqrt(1/PRNT50_le_low_Cases + 1/PRNT50_le_low_NonCases + 
                                               1/PRNT50_gt_high_Cases + 1/PRNT50_gt_high_NonCases)),
    # Calculate upper bound of 95% CI for high antibody level
    CI_upper_high = exp(log(OR_high) + 1.96 * sqrt(1/PRNT50_le_low_Cases + 1/PRNT50_le_low_NonCases + 
                                               1/PRNT50_gt_high_Cases + 1/PRNT50_gt_high_NonCases))
  ) 


# convert to factor
forest_df$Comparison <- factor(forest_df$Comparison, levels = unique(forest_df$Comparison))


## Plotting
ggplot(forest_df) +
  # Add points for Odds Ratios (lower threshold)
  geom_point(aes(x = log(OR_low), y = Comparison), color = "blue", size = 2) +  
  # Horizontal error bars (low threshold)
  geom_errorbarh(aes(xmin = log(CI_lower_low), xmax = log(CI_upper_low), y = Comparison), 
                 color = "blue", height = 0.2) + 
  # Add points for Odds Ratios (upper threshold)
  geom_point(aes(x = log(OR_high), y = Comparison), color = "red", size = 2) + 
  # Horizontal error bars (upper threshold)
  geom_errorbarh(aes(xmin = log(CI_lower_high), xmax = log(CI_upper_high), y = Comparison), 
                 color = "red", height = 0.2) + 
  # Vertical line where log odds = 0
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  scale_x_continuous(limits = c(-3, 5)) +
  labs(
    x = "Log Odds Ratio",
    y = "Comparison",
    #title = "Forest Plot of Dengue Infection Odds Ratios",
    color = "Threshold"
  ) +
  theme_minimal() 


## Basic Plot 5################################################################################################
## Plot 5 (bar): Distribution between examination date and onset date
###############################################################################################################
## Re-format the date columns
onset_exam_date_data <- DENV_cleaned_comp %>%
  mutate(
    # Clean 'Date of onset'
    onset_date = case_when(
      str_detect(`Date of onset`, "-") ~ as.Date(`Date of onset`, format = "%Y-%m-%d"),
      str_detect(`Date of onset`, "/") ~ as.Date(`Date of onset`, format = "%d/%m/%Y"),
      str_detect(`Date of onset`, 
                 "\\d{1,2} [A-Za-z]{3} \\d{2}") ~ as.Date(`Date of onset`, format = "%d %b %y"),
      TRUE ~ NA  
    ),
    # Clean 'Examination date'
    exam_date = case_when(
      str_detect(`Examination date`, "-") ~ as.Date(`Examination date`, format = "%Y-%m-%d"),
      str_detect(`Examination date`, "/") ~ as.Date(`Examination date`, format = "%d/%m/%Y"),
      str_detect(`Examination date`, 
                 "\\d{1,2} [A-Za-z]{3} \\d{2}") ~ as.Date(`Examination date`, format = "%d %b %y"),
      TRUE ~ NA   
    )
  )


## Calculate the date difference
onset_exam_date_data <- onset_exam_date_data %>% 
  mutate(
    date_difference = as.integer(exam_date - onset_date)
  )

###############################################################################################################
# Filter out the targeted data set
DENV_PCR_pos <- onset_exam_date_data %>% filter(PCR_result == "positive")
DENV_non_case <- onset_exam_date_data %>% filter(Current == "Neg")
PCR_neg_IgM_pos <- onset_exam_date_data %>% filter(PCR_result == "negative" 
                                                   & `MAC-ELISA- Final result` == "Positive")

# Create frequency tables for each data sets
date_diff_table_all <- table(onset_exam_date_data$date_difference)
date_diff_table_pcr_pos <- table(DENV_PCR_pos$date_difference)
date_diff_table_non_case <- table(DENV_non_case$date_difference)
date_diff_table_pcrn_igmp <- table(PCR_neg_IgM_pos$date_difference)

# For all individuals
filtered_date_diff_table_all <- date_diff_table_all[as.numeric(names(date_diff_table_all)) >= 0 
                                                    & as.numeric(names(date_diff_table_all)) <= 10]
percent_date_diff_all <- (filtered_date_diff_table_all / sum(filtered_date_diff_table_all)) * 100

# For PCR positive cases
filtered_date_diff_table_pcr_pos <- date_diff_table_pcr_pos[as.numeric(names(date_diff_table_pcr_pos)) >= 0
                                                        & as.numeric(names(date_diff_table_pcr_pos)) <= 10]
percent_date_diff_pcr_pos <- (filtered_date_diff_table_pcr_pos / sum(filtered_date_diff_table_pcr_pos)) * 100

# For Non-cases
filtered_date_diff_table_non_case <- date_diff_table_non_case[as.numeric(names(date_diff_table_non_case)) >= 0 
                                                      & as.numeric(names(date_diff_table_non_case)) <= 10]
percent_date_diff_non_case <- (filtered_date_diff_table_non_case / sum(filtered_date_diff_table_non_case)) * 100

# For PCR negative but IgM positive individuals
filtered_date_diff_table_pcrn_igmp <- date_diff_table_pcrn_igmp[as.numeric(names(date_diff_table_pcrn_igmp)) >= 0 
                                                         & as.numeric(names(date_diff_table_pcrn_igmp)) <= 10]
percent_date_diff_pcrn_igmp <- (filtered_date_diff_table_pcrn_igmp / sum(filtered_date_diff_table_pcrn_igmp)) * 100


# Bar plot
barplot(rbind(percent_date_diff_all, percent_date_diff_pcr_pos, 
              percent_date_diff_non_case, percent_date_diff_pcrn_igmp),
        beside = TRUE,
        col = c("#fbb4ae", "#b3cde3", "#ccebc5", "#decbe4"),
        names.arg = names(filtered_date_diff_table_all),
        xlab = "Date Difference",
        ylab = "Percentage",
        main = "Distribution of Date Differences",
        legend.text = c("All the individuals", "PCR positive cases", "Non-cases", "PCR- but IgM+"),
        args.legend = list(x = "topright", bty = "n"))

## Basic Plot 6################################################################################################
## Plot 6 (line): Cumulative distribution plots for PRNT50 values 
###############################################################################################################
# Function to generate the cummulative distribution of the PRNT50 plots for each serotype given
generate_cumulative_prnt50_plot <- function(DENV_cleaned_comp, serotype) {
  # Define the column name for the selected serotype
  serotype_col <- paste0("PRNT50 DENV", serotype)
  
  # Filter out targeted data sets 
  PCRn_IgMp <- DENV_cleaned_comp %>% filter(PCR_result == "negative" & `MAC-ELISA- Final result` == "Positive")
  PCRp <- DENV_cleaned_comp %>% filter(PCR_result == "positive")
  
  # Calculate frequency tables for log2(PRNT50) values
  table_PCRn_IgMp <- table(round(log2(PCRn_IgMp[[serotype_col]])))
  table_PCRp <- table(round(log2(PCRp[[serotype_col]])))
  
  # Include counts for NA values (indicating "0" level)
  table_PCRn_IgMp <- c("0" = sum(is.na(PCRn_IgMp[[serotype_col]])), table_PCRn_IgMp)
  table_PCRp <- c("0" = sum(is.na(PCRp[[serotype_col]])), table_PCRp)
  
  # Calculate percentages
  table_PCRn_IgMp_perc <- table_PCRn_IgMp / sum(table_PCRn_IgMp) * 100
  table_PCRp_perc <- table_PCRp / sum(table_PCRp) * 100
  
  # Calculate cumulative percentages
  table_PCRn_IgMp_cumperc <- cumsum(table_PCRn_IgMp_perc)
  table_PCRp_cumperc <- cumsum(table_PCRp_perc)
  
  # Create data stes for plotting
  df_PCRn_IgMp <- data.frame(Value = as.numeric(names(table_PCRn_IgMp_cumperc)),
                     cum_perc = table_PCRn_IgMp_cumperc, data_set = "PCR- IgM+")
  df_PCRp <- data.frame(Value = as.numeric(names(table_PCRp_cumperc)),
                     cum_perc = table_PCRp_cumperc, data_set = "PCR+")
  
  
  # Combine data sets
  combined_df <- rbind(df_PCRn_IgMp, df_PCRp)
  
  # Plot cumulative distribution
  ggplot(combined_df, aes(x = Value, y = cum_perc, color = data_set)) +
    geom_step() +
    labs(title = paste("Cumulative Distribution of PRNT50 DENV", serotype, "Values"),
         x = paste("log2(PRNT50 DENV", serotype, ")"),
         y = "Cumulative Percentage",
         color = "Dataset") +
    ylim(c(0, 100)) +
    theme_minimal()
}

# Loop throgh each serotype and apply the function
cum_dist <- list()
for (i in 1:4) {
  cum_dist[[i]] <- generate_cumulative_prnt50_plot(DENV_cleaned_comp, i)
}

# Display plots for each serotype on the same scene
gridExtra::grid.arrange(grobs = cum_dist, ncol = 2)


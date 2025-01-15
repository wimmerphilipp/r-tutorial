# always load your packages at the beginning of your script. You only need to
# activate the package once per sessions. Also, if you are working in a team
# your colleagues will get a better overview. 
library(openxlsx)
library(dplyr)
library(ggplot2)
library(doBy)
library(convey)
library(survey)

# or with 

# function to install and load packages
install_and_load <- function(package_names) {
  for (package in package_names) {
    if (!require(package, character.only = TRUE, quietly = TRUE)) {
      install.packages(package, dependencies = TRUE)
    }
    library(package, character.only = TRUE, quietly = TRUE)
  }
}


# list of packages to install and load
packages <- c(
  "webshot", "ineq", "dineq","broom",
  "sf", 
  "convey", "descr",
  "EnvStats", 
  "forcats", "ggstance",
  "Hmisc", 
  "ggtext", "ggthemes", 
  "grid", "gridExtra", 
  "knitr", "lemon",
  "matrixStats", "mitools", 
  "plotly", "readxl", 
  "spatstat", "stringr",
  "stargazer", "survey", 
  "tidyverse", "tinytex", "vtable", 
  "xtable",
  "viridis", 
  "bibtex", "binsreg", "writexl", "knitcitations", "highcharter","tigris","rmapshaper",
  "rbenchmark","tmap",
  "leaflet","terra", "pdftools","stringr","openxlsx", "wesanderson"
)


# install and load packages
install_and_load(packages)

## Read Data ##
# if you always save your data in the same path as the project, you can easily
# access the data, without specifying the full path. and, again this is good
# for sharing work, as you only have to copy the project's folder:
df <- read.xlsx("example.xlsx")
# obviously you could still use the full path. Note, that you need to change \ to /
df_1 <- read.xlsx("C:/Users/diexl/Documents/Unijob/WU Department/Tutor WiPol/R Tutorium/R-Tutorium-SoSe23/example.xlsx")
# if you need a specific sheet of your excel-file use: (this is only an example)
df_2 <- read.xlsx("data/example.xlsx", sheet = "Sheet2")
# getwd gets you your current working directory so you can check where your files
# are saved. use setwd() to change your current working directory.
getwd()
# within the brackets you should specify the link to your directory.
setwd()
# read CSV-files:
# this will be our practice file: it contains data on daily covid numbers in
# Austria. The data starts on 26.02.2020 and ends on 03.05.2021. We will 
# have a closer look on what data is contained within the dataset, manipulate
# it and do some descriptives. 
df<- read.csv("ts_covid_sortiert.csv", sep=";")
str(df$SiebenTageInzidenzFaelle)
table(df$SiebenTageInzidenzFaelle)
df$SiebenTageInzidenzFaelle <- gsub(",",".", df$SiebenTageInzidenzFaelle)
df$SiebenTageInzidenzFaelle <- as.numeric(df$SiebenTageInzidenzFaelle)
# alternatively you can also use the dplyr package
df <- df %>%
  mutate(SiebenTageInzidenzFaelle = as.numeric(gsub(",", ".", SiebenTageInzidenzFaelle)))
str(df$SiebenTageInzidenzFaelle)

## Data Cleaning ##
### have a first look at the data ###
# use View() to open a new tab with your data table, or click on it in the 
# environment
View(df)
# get column names
colnames(df)
# have a look at the structure of your data. Note that this is the same as the
# information in your environment
str(df)
# get a summary for each column
summary(df)
# have a look at the first few rows
head(df)

table(df$Bezirk)

# scatter plot: why is there such a gap in between the points?
plot(df$AnzahlFaelleSum)

summary(df$AnzahlFaelleSum)

dftest <- df %>% filter(AnzahlFaelleSum > 50000)|>
  select(AnzahlFaelleSum, Bezirk, Time)

# now we want to get rid of the time dimension. Therefore, we extract the day
wien <- subset(df, Bezirk=="Wien")
plot(as.ts(wien$AnzahlFaelleSum), lwd=2, col = "red")

# frequency distribution of numerical vectors
# you can also check via histogram
hist(df$AnzahlFaelle)
hist(df$AnzahlFaelle, breaks = "sturges")
hist(df$AnzahlFaelle, breaks = "fd")
hist(df$AnzahlFaelle, breaks = 3)
hist(df$AnzahlFaelle, breaks = 100)
# now we want to get rid of the time dimension. Therefore, we extract the day
# when the covid-incidence was the highest
# first we need to find out which day this is. max() gives you the highest 
# covid-incidence contained in the data.
max_value <- max(df$SiebenTageInzidenzFaelle)
# now we use the returned value to subset our dataset. [,] this kind of brackets
# accesses the dataframe. Entries before the comma refer to the rows of your
# dataframe, whereas after the comma refer to the column. Think of this code 
# line as: subset the data at the row where "SiebenTageInzidenzFaelle" equals
# our max()-value and also return the corresponding column entries for Time and
# district.
max_date <- df[df$SiebenTageInzidenzFaelle == max_value, 
               c("Time", "Bezirk", "SiebenTageInzidenzFaelle")]
max_date
# now we know that on the 12.11.2020 Rohrbach had the highest covid incidence 
# in our observation period. therefore, we now get rid of the time dimension
# by extracting this day.
data_1 <- subset(df, Time == "2020-11-12")
head(data_1)
# however, we no longer care about the Time-column and we also do not need the
# number of people who have recovered. (look at the changes in your
# envrionment). Also, note that instead of using the subset()-function we
# directly adressed the index of the columns of interest.
data_2 <- data_1[,-c(1, 12)]
head(data_2)
# however, we could have done this more efficiently
data <- subset(df, Time == "2020-11-12", select = c(2:10))
# or equivalently using dplyr
data <- df %>%
  filter(Time == "2020-11-12") %>%
  select(2:10)
head(data)
# You might also see a different type of "coding": the dplyr-version. "dplyr" is
# a package designed for data manipulation which uses tubes (%>%). In most cases
# it does the same thing. However, the subset()-function sometimes struggles 
# with large datafiles - more than a million entries. The dplyr-package seems to
# be able to handle this better.
?dplyr
# now follow the instructions of the description: enter 
# "browseVignettes(package = "dplyr")" into your console. Your browser will 
# open and display a help page. Choose "introduction to dplyr" and open the 
# R-code version. Here you get usefull advice on how to handle the package.

# when looking at the data, we can also see that the incidence uses decimal 
# values with a comma instead of a point. we need to change this in order to 
# use the variable
# you can also combine tubes
# 
rm(list = ls())

# Load the data and fix encoding issues
df <- read.csv("ts_covid_sortiert.csv", sep = ";")  # Load data from CSV file with semicolon as separator
df$Bezirk <- iconv(df$Bezirk, from = "latin1", to = "UTF-8", sub = "byte")  # Convert 'Bezirk' column encoding to UTF-8

# Prepare the data for plotting
data <- df %>% 
  filter(Time == "2020-11-12") %>%  # Filter rows for a specific date: 12th November 2020
  select(Bezirk, GKZ, AnzEinwohner, AnzahlFaelle, AnzahlFaelleSum, SiebenTageInzidenzFaelle) %>%  # Select relevant columns
  rename(Inzidenz = SiebenTageInzidenzFaelle) %>%  # Rename 'SiebenTageInzidenzFaelle' to 'Inzidenz' for simplicity
  mutate(Inzidenz = as.numeric(gsub(",", ".", Inzidenz)))  # Convert Inzidenz to numeric, replacing ',' with '.'

# Quick barplot to check the data
barplot(data$Inzidenz)  # Basic barplot of incidence values

# Improved barplot with axis labels
barplot(data$Inzidenz, 
        names.arg = data$Bezirk,  # Use district names ('Bezirk') as labels
        ylab = "Number of positive cases on the reference day")  # Y-axis label

# Define colors for barplot (9 unique colors for differentiation)
colors <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", 
            "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22")

# Create a color indicator based on the first digit of 'GKZ' (district code)
country_indicator <- as.integer(substr(data$GKZ, 1, 1))  # Extract the first character of 'GKZ' and convert to integer
bar_colors <- colors[(country_indicator %% length(colors)) + 1]  # Map colors cyclically using modulo operation

# Final barplot with enhancements
par(las = 2, mar = c(12, 5, 2, 2))  # Rotate axis labels (las = 2) and adjust margins
barplot(data$Inzidenz, 
        names.arg = data$Bezirk,  # District names as labels
        ylab = "Number of positive cases on the reference day",  # Y-axis label
        col = bar_colors,  # Apply colors to bars
        border = NA,  # Remove borders around bars
        space = 0.5,  # Adjust spacing between bars
        cex.names = 0.7)  # Reduce font size of axis labels for readability

# Scatterplot: Number of inhabitants vs. number of positive cases
scatter.smooth(data$AnzahlFaelle, data$AnzEinwohner,  # Smooth scatterplot of positive cases vs inhabitants
               xlab = "Number of people tested positive",  # X-axis label
               ylab = "Number of inhabitants")  # Y-axis label

# Identify and exclude Vienna as an outlier (assuming Vienna is at row 94)
data_novie <- data[-94,]  # Exclude row 94 (Vienna) from the dataset

# Scatterplot without Vienna (outlier removed)
scatter.smooth(data_novie$AnzahlFaelle, data_novie$AnzEinwohner,  # Smooth scatterplot without Vienna
               family = "gaussian",  # Use Gaussian smoothing
               xlab = "Number of people tested positive",  # X-axis label
               ylab = "Number of inhabitants",  # Y-axis label
               lpars = list(col = "red", lwd = 3, lty = 1))  # Customize line: red color, thicker line, solid style

################################################################################
# 1. Clean up the environment
rm(list = ls())  # Remove all objects from the R environment to start fresh.

# https://www.gesis.org/en/missy/metadata/EU-SILC/2013/Cross-sectional/original

# 2. Introduction to EU-SILC data
# EU-SILC contains data across four different sheets:
# - 'p': Individual-level data (e.g., income, demographics)
# - 'r': Register data (technical details about individuals)
# - 'd': Household-level data (economic data at the household level)
# - 'h': Technical household data

# 3. Load the four datasets
load("AT_CD_2013_small.RData")  # Load household-level data ('d')
load("AT_CH_2013_small.RData")  # Load technical household-level data ('h')
load("AT_CP_2013_small.RData")  # Load individual-level data ('p')
load("AT_CR_2013_small.Rdata")  # Load register data ('r')

# Note: These files are sample data for Austria (2013) and are not representative.

# 4. Create unique ID variables for merging datasets later
# Create a unique household ID ('id_h') for each file
aut_h$id_h <- paste(aut_h$HB020, aut_h$HB030, sep="")  # Combine Country-Code + Household ID from 'h'
aut_d$id_h <- paste(aut_d$DB020, aut_d$DB030, sep="")  # Combine Country-Code + Household ID from 'd'
aut_p$id_h <- paste(aut_p$PB020, aut_p$PX030, sep="")  # Combine Country-Code + Household ID from 'p'
aut_r$id_h <- paste(aut_r$RB020, aut_r$RX030, sep="")  # Combine Country-Code + Household ID from 'r'

# Create a unique personal ID ('id_p') for each file
aut_p$id_p <- paste(aut_p$PB020, aut_p$PB030, sep="")  # Combine Country-Code + Personal ID from 'p'
aut_r$id_p <- paste(aut_r$RB020, aut_r$RB030, sep="")  # Combine Country-Code + Personal ID from 'r'

# 5. Compare and merge household data
# Check for households that are in 'd' but not in 'h'
testjoin_hh <- anti_join(aut_d, aut_h, by="id_h")  # Find rows in 'aut_d' with no match in 'aut_h'
testjoin_hh  # Display the unmatched rows

# Merge the household datasets: 'd' (economic) and 'h' (technical)
aut_hh <- right_join(aut_d, aut_h, by="id_h")  # Keep all rows from 'h', add matching rows from 'd'

# 6. Compare and merge individual-level data
# Check for individuals that are in 'r' but not in 'p'
testjoin_p <- anti_join(aut_r, aut_p, by = c("id_p" = "id_p", "id_h" = "id_h"))  
# `anti_join` finds rows in 'aut_r' with no match in 'aut_p'

# The result shows that some individuals appear in the register ('r') but do not have full data in 'p'.
# This happens when interviews were incomplete, so only technical data exists.

# Merge the individual datasets: 'r' (register) and 'p' (data)
aut_indiv <- right_join(aut_r, aut_p, by = c("id_p" = "id_p", "id_h" = "id_h"))  
# Keep all rows from 'p', add matching rows from 'r'

# 7. Merge household data with individual data
# Combine the merged household dataset (aut_hh) with the merged individual dataset (aut_indiv)
data <- right_join(aut_hh, aut_indiv, by = "id_h")  
# Keep all rows from 'aut_indiv', add matching household data from 'aut_hh'

# 8. Clean up unnecessary objects
rm(aut_d, aut_h, aut_p, aut_r, aut_indiv, testjoin_hh, testjoin_p)  # Remove temporary objects to save memory

# 9. View the final combined dataset
View(data)  # Open the merged dataset in the R Viewer

# 1. Subset the dataset and rename variables for clarity
# Select relevant variables and rename them for readability
data.temp <- data %>% 
  select(id_p, id_h, RB080, RB090, PB210, PE040, PL031, HY020, HS120, HX040, 
         hh_weight = DB090, p_weight = RB050, DB040, RB050, DB020, PY010G, 
         PL040, PL060, DB020, HX090, HH050) %>%
  rename(birthyear = RB080,          # Year of birth
         gender = RB090,             # Gender: 1 - male, 2 - female
         birthcountry = PB210,       # Birth country
         educ = PE040,               # Education level
         economicstatus = PL031,     # Economic activity status
         hh_income = HY020,          # Household income
         makeendsmeet = HS120,       # Household financial situation
         hh_size = HX040,            # Household size
         region = DB040,             # Region
         country = DB020,            # Country
         p_income = PY010G,          # Personal income
         emp_status = PL040,         # Employment status
         h_worked_week = PL060,      # Hours worked per week
         arp = HX090,                # At-risk-of-poverty indicator
         warm = HH050)               # Ability to keep home warm
View(data.temp)

# 2. Create an age variable
# Calculate age by subtracting the year of birth from 2013
data.temp <- data.temp %>%
  mutate(age = 2013 - birthyear)

# 3. Inspect the dataset
str(data.temp)       # Structure of the dataset
glimpse(data.temp)   # Compact summary of variables
head(data.temp)      # First few rows of the dataset
dim(data.temp)       # Dimensions: number of rows and columns
is.na(data.temp)     # Check for missing values

# 4. Remove rows with missing values
data.temp <- na.omit(data.temp)  # Remove rows with NAs
View(data.temp)
sum(is.na(data.temp))  # Verify no missing values remain

# 5. Population counts using survey weights
# Total population represented by individual weights
n.pop <- sum(data.temp$p_weight)  # Sum of individual weights
# Total population using household weights
n.pop.h <- sum(data.temp$hh_weight)
n.pop; n.pop.h  # Compare both values

# 6. Gender distribution
# Count males (1) and females (2) in the sample
table(data.temp$gender)

# Use weights to calculate population counts for males and females
sum((data.temp$gender == 1) * data.temp$p_weight)  # Male population
sum((data.temp$gender == 2) * data.temp$p_weight)  # Female population

# Alternative with dplyr
males <- data.temp %>% 
  filter(gender == 1) %>%
  summarise(males_in_pop = sum(p_weight))
males

# Calculate population for both genders using group_by
pop.by.gender <- data.temp %>%
  group_by(factor(gender)) %>%
  summarise(population = sum(p_weight))|>
  ungroup()# Sum weights by gender
pop.by.gender

# Calculate the share of women in the population
n.fem.pop <- sum((data.temp$gender == 2) * data.temp$p_weight)
share.fem.pop <- n.fem.pop / n.pop  # Share of women
share.fem.pop

# Share of women in the sample (not weighted)
n <- nrow(data.temp)                       # Total sample size
n.fem.sample <- sum(data.temp$gender == 2) # Count females
share.fem.sampl <- n.fem.sample / n        # Share of women in sample
share.fem.sampl

# 7. Education summary
summary(data.temp$educ)  # Summary of education variable
table(data.temp$educ)    # Frequency table of education levels

# Population share by education level
out <- data.temp %>%
  group_by(educ) %>%
  summarise(pop = sum(p_weight))|>
  ungroup()# Sum weights for each education level
out$share <- out$pop / n.pop * 100  # Calculate share as a percentage
sum(out$share)  # Verify shares sum to 100%

# 8. Barplot for education levels
barplot(out$share, 
        names.arg = c("Pre-Primary", "Low. Second", "Second", "Post Second", "Tertiary"), 
        ylab = "% of Population", ylim = c(0, 60), 
        main = "Education", col = "firebrick")

# Add colors with RColorBrewer
coul <- RColorBrewer::brewer.pal(5, "Set2")
barplot(out$share, 
        names.arg = c("Pre-Primary", "Low. Second", "Second", "Post Second", "Tertiary and Higher"), 
        ylab = "% of Population", ylim = c(0, 60), 
        main = "Education", col = coul)

# 9. Barplot for gender shares
out.gender <- data.temp %>%
  group_by(gender) %>%
  summarise(pop = sum(p_weight))|> # Population by gender
  ungroup ()
out.gender$share <- out.gender$pop / n.pop * 100  # Share as percentage

barplot(out.gender$share, 
        names.arg = c("Males", "Females"), 
        ylab = "% Percent", main = "Gender", col = "firebrick",ylim = c(0, 60))

# Using ggplot2 for gender shares
values <- c("#999999", "#E69F00")
lab <- c("Males", "Females")

ggplot(data = out.gender, aes(x = factor(gender), y = share, fill = factor(gender))) +
  geom_bar(stat = "identity") +
  scale_x_discrete("Gender") +
  scale_y_continuous("Percent") +
  scale_fill_manual("Pop by Gender", values = values, labels = lab)

# 10. Create survey designs for weighted analysis
# Filter rows where personal income is positive
data.temp <- data.temp %>% filter(p_income > 0)  

# Individual-level survey design
data.pd.svy <- svydesign(
  ids = ~id_p,           # Primary Sampling Unit (PSU) is the individual ID
  strata = ~region,      # Stratification variable: 'region'
  weights = ~p_weight,   # Use personal survey weights to ensure representation
  data = data.temp       # Input dataset
) %>% 
  convey_prep()      # Prepare for inequality analysis (e.g., Lorenz curve)

# Household-level survey design
data.hd.svy <- svydesign(
  ids = ~id_h,           # PSU is the household ID
  strata = ~region,      # Stratification variable: 'region'
  weights = ~hh_weight,  # Use household weights
  data = data.temp       # Input dataset
) %>% 
  convey_prep()

# 11. Lorenz Curve for household income
# Calculate and plot the Lorenz curve to show income inequality
lorenz <- svylorenz(
  ~hh_income,                                      # Variable: household income
  design = subset(data.hd.svy, !is.na(hh_income)), # Exclude rows with missing income
  quantiles = seq(0, 1, 0.1),                      # Deciles for cumulative distribution
  alpha = 0.05,                                    # Confidence interval level
  curve.col = "red"                                # Color for the curve
)

# 12. Histograms for household income
# Weighted histogram for household income
svyhist(
  ~hh_income, 
  design = data.hd.svy, 
  main = "Household Income", 
  col = "deeppink4", 
  breaks = 100           # Number of bins
)

# Subset for incomes under 75,000 and plot histogram
svyhist(
  ~hh_income, 
  design = subset(data.hd.svy, hh_income < 75000), 
  main = "Household Income (< 75k)", 
  col = "deeppink4", 
  breaks = 75
)

# Add a smoothed density line to visualize distribution
lines(
  svysmooth(~hh_income, design = subset(data.hd.svy, hh_income > 0 & hh_income < 75000)), 
  lwd = 2, 
  col = "gray"           # Line color
)

# 13. Wage regressions
# Prepare the data for wage regression analysis
data.reg <- data.temp %>% 
  filter(emp_status == 3, h_worked_week > 0) %>%     # Only employed individuals working > 0 hours
  mutate(hwages = p_income / (h_worked_week * 52)) %>% # Calculate hourly wage
  mutate(gender = factor(gender, labels = c("Male", "Female"))) # Convert gender to factor

# OLS regression without survey design (unweighted)
nonsvyols <- lm(hwages ~ age + gender, data = data.reg)

# OLS regression with survey design (weighted)
reg.pd.svy <- svydesign(
  ids = ~id_p, strata = ~region, weights = ~p_weight, data = data.reg
)
svyols <- svyglm(hwages ~ age + gender, design = reg.pd.svy)

# Compare regression results
summary(svyols)       # Weighted regression results
summary(nonsvyols)    # Unweighted regression results

# Compare residual distributions
par(mfrow = c(1, 2))  # Set layout for two plots
hist(svyols$residuals, main = "Residuals (Weighted)", col = "lightblue")
hist(nonsvyols$residuals, main = "Residuals (Unweighted)", col = "lightgreen")



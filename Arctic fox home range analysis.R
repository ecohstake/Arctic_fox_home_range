# load packages
library(adehabitatHR)
library(dplyr)
library(sp)
library(sf)
library(lubridate)
library(grid)
library(tibble)
library(glm2)
library(tidyverse)
library(car)
library(ctmm)
library(lme4)

# set working directory

# load raw data
raw_data <- read.csv("raw fox data 2018-2025.csv")

# load filtered data
fox_data <- read.csv("Anonymized Arctic fox data.csv")
fox_data$ID <- as.factor(fox_data$ID)
fox_data$timestamp <- as.POSIXct(fox_data$timestamp, format = "%Y-%m-%d %H:%M:%S")
fox_data <- fox_data[,-1]

# load AKDE 90% home range data
data <- read.csv("wAKDE90 monthly (until Feb 2025).csv")

# load overlap data
bc <- read.csv("monthly BC, AKDE, and size diff.csv")

###########################################
#### HOME RANGE AND OVERLAP ESTIMATION ####
###########################################
#### Monthly home range estimation (AKDE) ####

# Example:
month_data <- fox_data %>%
  filter(month(timestamp) %in% c(1))   # 1 = January. Repeat for 1-12 for all twelve months.
month_data <- as.telemetry(month_data)

FITS <- list()
for(i in 1:length(month_data)) {
  GUESS <- ctmm.guess(month_data[[i]], interactive = FALSE)
  FITS[[i]] <- ctmm.fit(month_data[[i]], GUESS, method = "pHREML", trace = 2)
}

AKDES <- list()
AKDES <- akde(month_data, FITS, weights = TRUE, trace = 2)

spdf_list <- list()
for(i in 1:length(AKDES)) {
  # level.UD = 0.90 estimates the 90% isopleth. Change to 0.50 for the 50% isopleth.
  spdf_list[[i]] <- SpatialPolygonsDataFrame.UD(AKDES[[i]],level.UD=0.90,level=0.95)
}



#### Monthly overlap (BC) ####
# Example: overlap between May-June
# Repeat for all pairs of successive months

month_1 <- fox_data %>%
  filter(month(timestamp) %in% c(5)) %>%
  mutate(ID = paste0(ID, " 1"))

month_2 <- fox_data %>%
  filter(month(timestamp) %in% c(6)) %>%
  mutate(ID = paste0(ID, " 2"))

fox_months <- rbind(month_1, month_2)
fox_months <- as.telemetry(fox_months)

FITS <- list()
for(i in 1:length(fox_months)) {
  GUESS <- ctmm.guess(fox_months[[i]], interactive = FALSE)
  FITS[[i]] <- ctmm.fit(fox_months[[i]], GUESS, method = "pHREML", trace = 2)
}

AKDES <- list()
AKDES <- akde(fox_months, FITS, weights = TRUE, trace = 2)
OVERLAP <- overlap(AKDES, method = "Bhattacharyya", level = .95)


#### Annual overlap (BC) ####

fox_data$timestamp <- as.POSIXct(fox_data$timestamp, format = "%Y-%m-%d %H:%M:%S")

# Fox 1
subset_data_1 <- subset(fox_data, ID %in% c("Fox 1 2018", "Fox 1 2019")) %>%  
  # Select period to compare between seasons
  filter(
    (timestamp >= as.POSIXct("2018-07-11") & timestamp <= as.POSIXct("2018-08-20")) |
      (timestamp >= as.POSIXct("2019-07-11") & timestamp <= as.POSIXct("2019-08-20"))
  )

subset_data_1 <- as.telemetry(subset_data_1)

GUESS_1 <- ctmm.guess(subset_data_1[[1]], interactive = FALSE)
GUESS_2 <- ctmm.guess(subset_data_1[[2]], interactive = FALSE)
#Fit and select movement models for each animal
FITS_1 <- ctmm.select(subset_data_1[[1]], GUESS_1)
FITS_2 <- ctmm.select(subset_data_1[[2]], GUESS_2)
# Return a summary of the selected models
summary(FITS_1)
summary(FITS_2)
# Estimate home ranges (utilization distributions)
FITS <- list(FITS_1, FITS_2)
HR_UDS_1 <- akde(subset_data_1, FITS, weights = TRUE)

# Calculate home range overlap (BC)
overlap(HR_UDS_1, method = "Bhattacharyya", level = .95)



# Fox 5
subset_data_5 <- subset(data_all, ID %in% c("Fox 5 2019", "Fox 5 2020")) %>%
  filter(
    (timestamp >= as.POSIXct("2019-07-22") & timestamp <= as.POSIXct("2019-09-15")) |
      (timestamp >= as.POSIXct("2020-07-22") & timestamp <= as.POSIXct("2020-09-15"))
  )

subset_data_5 <- as.telemetry(subset_data_5)

GUESS_1 <- ctmm.guess(subset_data_5[[1]], interactive = FALSE)
GUESS_2 <- ctmm.guess(subset_data_5[[2]], interactive = FALSE)
#Fit and select movement models for each animal
FITS_1 <- ctmm.select(subset_data_5[[1]], GUESS_1)
FITS_2 <- ctmm.select(subset_data_5[[2]], GUESS_2)
# Return a summary of the selected model
summary(FITS_1)
summary(FITS_2)

FITS <- list(FITS_1, FITS_2)
HR_UDS_5 <- akde(subset_data_5, FITS, weights = TRUE)


overlap(HR_UDS_5, method = "Bhattacharyya", level = .95)


# Fox 7
subset_data_7 <- subset(data_all, ID %in% c("Fox 7 2019", "Fox 7 2020")) %>%
  filter(
    (timestamp >= as.POSIXct("2019-07-25") & timestamp <= as.POSIXct("2019-10-19")) |
      (timestamp >= as.POSIXct("2020-07-25") & timestamp <= as.POSIXct("2020-10-19"))
  )

subset_data_7 <- as.telemetry(subset_data_7)

GUESS_1 <- ctmm.guess(subset_data_7[[1]], interactive = FALSE)
GUESS_2 <- ctmm.guess(subset_data_7[[2]], interactive = FALSE)
#Fit and select movement models for each animal
FITS_1 <- ctmm.select(subset_data_7[[1]], GUESS_1)
FITS_2 <- ctmm.select(subset_data_7[[2]], GUESS_2)
# Return a summary of the selected model
summary(FITS_1)
summary(FITS_2)

FITS <- list(FITS_1, FITS_2)
HR_UDS_7 <- akde(subset_data_7, FITS, weights = TRUE)

overlap(HR_UDS_7, method = "Bhattacharyya", level = .95)

# Fox 8
subset_data_8 <- subset(data_all, ID %in% c("Fox 8 2019", "Fox 8 2020")) %>%  
  filter(
    (timestamp >= as.POSIXct("2019-07-28") & timestamp <= as.POSIXct("2019-12-20")) |
      (timestamp >= as.POSIXct("2020-07-28") & timestamp <= as.POSIXct("2020-12-20"))
  )

subset_data_8 <- as.telemetry(subset_data_8)

GUESS_1 <- ctmm.guess(subset_data_8[[1]], interactive = FALSE)
GUESS_2 <- ctmm.guess(subset_data_8[[2]], interactive = FALSE)
#Fit and select movement models for each animal
FITS_1 <- ctmm.select(subset_data_8[[1]], GUESS_1)
FITS_2 <- ctmm.select(subset_data_8[[2]], GUESS_2)
# Return a summary of the selected model
summary(FITS_1)
summary(FITS_2)

FITS <- list(FITS_1, FITS_2)
HR_UDS_8 <- akde(subset_data_8, FITS, weights = TRUE)

overlap(HR_UDS_8, method = "Bhattacharyya", level = .95)


# Fox 16
subset_data_16 <- subset(data_all, ID %in% c("Fox 16 2023", "Fox 16 2024")) %>%
  filter(
    (timestamp >= as.POSIXct("2023-07-26") & timestamp <= as.POSIXct("2024-02-28")) |
      (timestamp >= as.POSIXct("2024-07-26") & timestamp <= as.POSIXct("2025-02-28"))
  )

subset_data_16 <- as.telemetry(subset_data_16)

GUESS_1 <- ctmm.guess(subset_data_16[[1]], interactive = FALSE)
GUESS_2 <- ctmm.guess(subset_data_16[[2]], interactive = FALSE)
#Fit and select movement models for each animal
FITS_1 <- ctmm.select(subset_data_16[[1]], GUESS_1)
FITS_2 <- ctmm.select(subset_data_16[[2]], GUESS_2)
# Return a summary of the selected model
summary(FITS_1)
summary(FITS_2)

FITS <- list(FITS_1, FITS_2)
HR_UDS_16 <- akde(subset_data_16, FITS, weights = TRUE)

overlap(HR_UDS_16, method = "Bhattacharyya", level = .95)

# Fox 17
subset_data_17 <- subset(data_all, ID %in% c("Fox 17 2023", "Fox 17 2024")) %>%
  filter(
    (timestamp >= as.POSIXct("2023-07-28") & timestamp <= as.POSIXct("2023-11-18")) |
      (timestamp >= as.POSIXct("2024-07-28") & timestamp <= as.POSIXct("2024-11-18"))
  )

subset_data_17 <- as.telemetry(subset_data_17)

GUESS_1 <- ctmm.guess(subset_data_17[[1]], interactive = FALSE)
GUESS_2 <- ctmm.guess(subset_data_17[[2]], interactive = FALSE)
#Fit and select movement models for each animal
FITS_1 <- ctmm.select(subset_data_17[[1]], GUESS_1)
FITS_2 <- ctmm.select(subset_data_17[[2]], GUESS_2)
# Return a summary of the selected model
summary(FITS_1)
summary(FITS_2)

FITS <- list(FITS_1, FITS_2)
HR_UDS_17 <- akde(subset_data_17, FITS, weights = TRUE)

overlap(HR_UDS_17, method = "Bhattacharyya", level = .95)


#### Meta-analysis ####
#remove non-residents
raw_data <- raw_data %>%
  filter(!(ID %in% "Fox 14") & 
           !(ID %in% "Fox 11") & 
           !(ID %in% "Fox 12") &
           !(ID %in% "Fox 18"))

# Load the tracking dataset:
foxes <- as.telemetry(raw_data)

# fit movement models
FITS <- list()
for(i in 1:length(foxes)) {
  GUESS <- ctmm.guess(foxes[[i]], interactive = FALSE)
  FITS[[i]] <- ctmm.fit(foxes[[i]], GUESS, trace = 2)
}

# calculate AKDEs on a consistent grid
AKDES <- list()
AKDES <- akde(foxes, FITS, trace = 2)

# cluster-analysis of foxes
cluster(AKDES, sort = TRUE)

# meta-analysis of fox home-range areas
## calculate mean home range area of an average individual:

meta(AKDES,
     col = c(COL,"black"), 
     verbose = TRUE, # verbose output with CIs
     sort = TRUE) 

# Population models / population range:

MEAN.FITS <- mean(FITS)
summary(MEAN.FITS) 

## What is the population range?
MEAN <- mean(AKDES) # distribution of the sample
plot(foxes, MEAN)

PKDE <- pkde(foxes, AKDES) # distribution of the population
plot(foxes, PKDE)

EXT <- extent(list(MEAN, PKDE))
COL <- c("red", "black", "blue")

####################
#### STATISTICS ####
####################
#### ANOVA: AKDE ~ month ####


# check assumptions of normality
hist(sqrt(data$AKDE90_km2))
qqnorm(sqrt(data$AKDE90_km2))
qqline(sqrt(data$AKDE90_km2))

# Set ID and month as factors
data$month <- as.factor(data$month)

# one-way ANOVA
model1 <- lm(sqrt(data$AKDE90_km2) ~ data$month)
anova(model1)

summary(data$AKDE90_km2)
sd(data$AKDE90_km2)

#### Correlation AKDE 90% and 50% ####
data_90 <- read.csv("wAKDE90 monthly (until Feb 2025).csv")
data_50 <- read.csv("wAKDE50 monthly (until Feb 2025).csv")
data_50 <- data_50 %>%
  mutate(ID = ifelse(ID == "Fox 15 2023", "Fox 15", ID))

combined_data <- left_join(data_50, data_90, by = c("ID", "month"))

cor.test(sqrt(combined_data$AKDE90_km2), sqrt(combined_data$AKDE50_km2))

#### ANOVA: monthly AKDE ~ month * sex ####

# remove foxes of unknown sex
data2 <- data[!(data$ID %in% c("Fox 17 2023", "Fox 17 2024")), ]
# Add a new column "sex" to the existing data frame
data2$sex <- ifelse(data2$ID %in% c("Fox 1 2018", "Fox 1 2019", "Fox 2", "Fox 3", "Fox 4", "Fox 5 2019", "Fox 5 2020", "Fox 6", "Fox 8 2019", "Fox 8 2020",
                                    "Fox 9", "Fox 13", "Fox 16 2023", "Fox 16 2024", "Fox 16 2025"), "Female", "Male")

data2$sex <- as.factor(data2$sex)
data2$month <- as.factor(data2$month)

hist(sqrt(data2$AKDE90_km2))

# non-parametric two-way ANOVA
data3 <- data2%>%
  filter(month %in% c("5", "6", "7", "8", "9", "10", "11"))

hist(sqrt(data3$AKDE90_km2))

lm2 <- lm(sqrt(AKDE90_km2) ~ month * sex, data = data3)
anova(lm2)


library(ARTool)
art_model <- art(AKDE90_km2 ~ sex * month, data = data3)
anova(art_model)


#### t-test: monthly AKDE ~ sex ####

# remove foxes of unknown sex
data2 <- data[!(data$ID %in% c("Fox 17 2023", "Fox 17 2024")), ]
# Add a new column "sex" to the existing data frame
data2$sex <- ifelse(data2$ID %in% c("Fox 1 2018", "Fox 1 2019", "Fox 2", "Fox 3", "Fox 4", "Fox 5 2019", "Fox 5 2020", "Fox 6", "Fox 8 2019", "Fox 8 2020",
                                    "Fox 9", "Fox 13", "Fox 16 2023", "Fox 16 2024", "Fox 16 2025"), "Female", "Male")

# only include months where both sexes are represented
data3 <- data2%>%
  filter(month %in% c("5", "6", "7", "8", "9", "10", "11"))

data3$sex <- as.factor(data3$sex)
data3$month <- as.factor(data3$month)

hist(sqrt(data3$AKDE90_km2))

# non-parametric test
wilcox.test(AKDE90_km2 ~ sex, data = data3)




#### Monthly overlap ####

# mean monthly bc
mean(bc$overlap)


# linear regression: Monthly overlap ~ hr size
hist(bc$overlap)
lm_mean <- lm(overlap ~ AKDE90_mean, data = bc)
summary(lm_mean)
plot(lm_mean)

# linear regression: Monthly overlap ~ hr size change
bc$AKDE90_diff <- abs(bc$AKDE90_diff)   # make  diff values positive
lm_diff <- lm(overlap ~ AKDE90_diff, data = bc)
summary(lm_diff)

# Extract residuals
residuals <- residuals(lm_mean)
qqnorm(residuals)  # Create the Q-Q plot
qqline(residuals, col = "red")  # Add a reference line
# Histogram of residuals
hist(residuals, breaks = 20, main = "Histogram of Residuals", 
     xlab = "Residuals", col = "lightblue", border = "black")


# overlap ~ months
hist(bc$overlap)
kruskal.test(overlap ~ months, data = bc)



# overlap ~ sex
# remove foxes of unknown sex
bc2 <- bc[!(bc$ID %in% c("Fox 17 2023", "Fox 17 2024")), ]
# Add a new column "sex" to the existing data frame
bc2$sex <- ifelse(bc2$ID %in% c("Fox 1 2018", "Fox 1 2019", "Fox 2", "Fox 3", "Fox 4", "Fox 5 2019", "Fox 5 2020", "Fox 6", "Fox 8 2019", "Fox 8 2020",
                                "Fox 9", "Fox 13", "Fox 16 2023", "Fox 16 2024", "Fox 16 2025"), "Female", "Male")
# only include months where both sexes are represented
bc2 <- bc2%>%
  filter(month %in% c("5", "6", "7", "8", "9", "10", "11"))
hist(bc2$overlap)
wilcox.test(overlap ~ sex, data = bc2)
bc2 %>%
  group_by(sex) %>%
  summarise(mean_overlap = mean(overlap, na.rm = TRUE),
            median_overlap = median(overlap, na.rm = TRUE),
            n = n())



#### Rodent density ####
#remove Råsto foxes
data_rodent <- data %>% filter(!(ID %in% c("Fox 16 2023", "Fox 16 2024", "Fox 16 2025", "Fox 17 2023", "Fox 17 2024")))

data_rodent$ID <- as.factor(data_rodent$ID)

# add year column
data_rodent$ID <- sub("Fox 2$", "Fox 2 2018", data_rodent$ID)
data_rodent$ID <- sub("Fox 3$", "Fox 3 2018", data_rodent$ID)
data_rodent$ID <- sub("Fox 4$", "Fox 4 2018", data_rodent$ID)
data_rodent$ID <- sub("Fox 6$", "Fox 6 2019", data_rodent$ID)
data_rodent$ID <- sub("Fox 9$", "Fox 9 2019", data_rodent$ID)
data_rodent$ID <- sub("Fox 12$", "Fox 12 2019", data_rodent$ID)
data_rodent$ID <- sub("Fox 13$", "Fox 13 2023", data_rodent$ID)
data_rodent$ID <- sub("Fox 15$", "Fox 15 2023", data_rodent$ID)
data_rodent <- data_rodent %>%
  mutate(year = as.numeric(sub(".* (\\d{4})$", "\\1", ID)),  # Extract the year
         ID = sub(" \\d{4}$", "", ID))  # Remove the year from ID

# add "rodents" column
data_rodent$rodents <- ifelse(data_rodent$year == 2018, "high",
                          ifelse(data_rodent$year == 2019, "high",
                                 ifelse(data_rodent$year == 2020, "low",
                                        ifelse(data_rodent$year == 2023, "low", NA))))

data_rodent$rodents <- as.factor(data_rodent$rodents)

# LMM
lmm <- lmer(AKDE90_km2 ~ rodents + (1 | ID) + (1 | month), data = data_rodent)
summary(lmm)





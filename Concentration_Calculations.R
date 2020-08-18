library(ggpmisc)
library(tidyverse)
options(scipen = 999)

# Linear fitting + prediction ---------------------------------------------
x <- 1:100
y <- 2 + 8 * x + rnorm(100, mean = 0, sd = 20)
plot(x, y)  
fit <- lm(y ~ x)   
abline(fit)        
predict(fit, data.frame(x = 10, interval = "confidence")) 

# Example from real data ---------------------------------------------

# This is an example of reading in a file from your CURRENT directory.
# Notice you don't need the tilde or the whole file path!

# Replace those values that act as characters in the dataset- they mess up the calculations!
replace_nonvalues <- function(x) (gsub("#N/A", NA, x))
replace_dashes <- function(x) (gsub("-", ".", x))
concentrations <- c("0uM", "0.5uM", "1uM", "2.5uM")
vial_types <- c("Durapore", "GF75", "Omnipore", "Vial")

Filter_Test <- read.csv("191121_BSA-Filter-Test_All.csv", stringsAsFactors = FALSE) %>% 
  select(Replicate.Name, Precursor.Ion.Name, Area) %>%
  mutate_at(c("Area"), replace_nonvalues) %>% # Replace the non-standard #N/A values. Can take multiple columns.
  mutate_at(c("Replicate.Name"), replace_dashes) %>%
  mutate(Area = as.numeric(Area)) %>%
  filter(!str_detect(Precursor.Ion.Name, "IS"))

# Isolate Durapores as an example
chooseFilter <- function(df, filter.type) {
  df.w.filter <- df %>%
    mutate(runtype = ifelse(str_detect(Replicate.Name, "Blk"), "Blank", "Sample")) %>%
    filter(str_detect(Replicate.Name, filter.type))
}

# Subtract the blanks from the areas
subtractBlanks <- function(df.w.filter) {
  df.blank.subtracted <- df.w.filter %>%
    group_by(Precursor.Ion.Name) %>%
    mutate(blankArea = Area[which(runtype == "Blank")]) %>%
    mutate(Area_noBlank = Area - blankArea)
  return(df.blank.subtracted)
}

# Rename the concentrations
createConcentration <- function(df, myconcentrations) {
  df.w.concentration <- df %>%
    mutate(Concentration = ifelse(str_detect(Replicate.Name,
                                             myconcentrations[1]), myconcentrations[1],
                                  ifelse(str_detect(Replicate.Name,
                                                    myconcentrations[2]), myconcentrations[2],
                                         ifelse(str_detect(Replicate.Name,
                                                           myconcentrations[3]), myconcentrations[3],
                                                ifelse(str_detect(Replicate.Name,
                                                                  myconcentrations[4]), myconcentrations[4], NA))))) %>%
    mutate(Concentration = substr(Concentration, 1, nchar(Concentration)-2),
           Concentration = as.numeric(Concentration)) %>%
    filter(runtype != "Blank") %>%
    na.omit() 
  return(df.w.concentration)
}


filters.chosen <- chooseFilter(Filter_Test, filter.type = "Durapore")
blanks.subtracted <- subtractBlanks(filters.chosen)
concentrations.added <- createConcentration(blanks.subtracted, myconcentrations = c("0uM", "0.5uM", "1uM", "2.5uM"))

# plots
ggplot(concentrations.added, aes(x=Concentration, y=Area_noBlank, group = Precursor.Ion.Name)) +
  facet_wrap(~Precursor.Ion.Name) +
  geom_point() + 
  geom_smooth(method=lm, se=TRUE, fullrange=TRUE) +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5)) +
  # stat_poly_eq(parse=T, aes(label = ..eq.label..), formula=y~x) +
  stat_poly_eq(aes(label=..eq.label..),
               geom="label", alpha=0.33, formula=(y ~ x),
               label.y = 0.9 * max(concentrations.added$Area_noBlank), 
               label.x = 0.5 * max(concentrations.added$Concentration),
               size = 2.5, parse=TRUE) +
  theme(text=element_text(size=10)) 


# Extract slope values
Slope.Values <- concentrations.added %>% 
  group_by(Precursor.Ion.Name) %>% 
  do({
    mod = lm(Area_noBlank ~ Concentration, data = .)
    data.frame(Intercept = coef(mod)[1],
               Slope = coef(mod)[2])
  })

All.Info <- concentrations.added %>%
  left_join(Slope.Values) %>%
  select(Replicate.Name, Precursor.Ion.Name, Area_noBlank, Concentration, Intercept, Slope)

Final.Concentrations <- All.Info %>%
  group_by(Precursor.Ion.Name) %>%
  mutate(Calculated.Concentration = (Area_noBlank - Intercept) / Slope)


# Use one set of concentrations to calculate another ------------------------------

New.Filters.Chosen <- chooseFilter(Filter_Test, filter.type = "GlassVial|GlassBlk_Vial")
New.Blanks.Subtracted <- subtractBlanks(New.Filters.Chosen)

New.Final.Concentration <- New.Blanks.Subtracted %>%
  left_join(Final.Concentrations %>% select(Precursor.Ion.Name, Intercept, Slope) %>% unique()) %>%
  group_by(Precursor.Ion.Name) %>%
  mutate(Calculated.Concentration = (Area_noBlank - Intercept) / Slope)




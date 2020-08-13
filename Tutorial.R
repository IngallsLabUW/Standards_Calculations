## This is the source editor, where R scripts are written.

library(tidyverse)
options(scipen = 999)

# Uploading data ----------------------------------------------------------

# This is an example of reading in a file from anywhere on your computer.
# The tilde "~" indicates your home, or starting, directory. The filepath proceeds from there.
mydata <- read.csv("~/work/Ingalls_Standards/Ingalls_Lab_Standards_NEW.csv",
                   stringsAsFactors = FALSE)

# Subsetting data
mydata_subset <- mydata %>%
  select(Compound.Type, Compound.Name, m.z) %>% 
  arrange(Compound.Name)

# Plotting data ----------------------------------------------------------

# Cool, but hard to see the text. 
mydata_plot <- ggplot(data = mydata_subset, aes(x = Compound.Name, y = m.z)) +
  geom_point()
mydata_plot

# Text is rotated 90 degrees, and color is added to differentiate the Compound.Type
mydata_plot <- ggplot(data = mydata_subset, 
                      aes(x = Compound.Name, y = m.z, color = Compound.Type)) +
  geom_point() +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5))
mydata_plot

# Specific data plotting ----------------------------------------------------------

# Let's focus on a few specific types of data to make this clearer
mydata_subset.AminoAcids <- mydata_subset %>%
  select(Compound.Type, Compound.Name, m.z) %>%
  arrange(Compound.Name) %>%
  filter(Compound.Type == "Amino Acid") %>%
  mutate(m.z = as.numeric(m.z))

mydata_plot2.AminoAcids <- ggplot(data = mydata_subset.AminoAcids, 
                      aes(reorder(Compound.Name, -m.z), y = m.z)) +
  geom_point() +
  theme(axis.text.x = element_text(angle = 90))
mydata_plot2.AminoAcids



# Facet wrapping ----------------------------------------------------------
my.data <- read.csv("191121_BSA-Filter-Test_All.csv") %>%
  select(Precursor.Ion.Name, Replicate.Name, Area) 

ggplot(data = my.data, aes(x = Replicate.Name, y = Area)) +
  geom_bar(stat = "identity") +
  facet_wrap(~Precursor.Ion.Name) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank())


# Taking averages of triplicates ---------------------------------------------
Replicates_Split <- my.data %>%
  filter(str_detect(Replicate.Name, "Smp")) %>%
  # Option 1: Split by underscore and drop everything but the SampID.
  separate(Replicate.Name,
           into = c("Date", "Type", "SampID", "Replicate"), sep = "_") %>%
  # Option 2: Drop last two characters in column.
  # mutate(SampID = substr(Replicate.Name, 1, nchar(Replicate.Name)-2)) %>%
  group_by(Precursor.Ion.Name, SampID) %>%
  mutate(Average = mean(Area, na.rm = TRUE))

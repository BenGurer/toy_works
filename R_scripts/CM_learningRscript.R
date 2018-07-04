# Cortical Magnification

# load tidyverse library
# install.packages("tidyverse")
library(tidyverse)

# set working directory to correct folder
setwd("C:/UoN/OneDrive - The University of Nottingham/data/testData")

# Importing the dataset
dataset = read.csv('11108_006_CM.csv')

ggplot(data = dataset) + 
  geom_point(mapping = aes(x = CorticalDistance, y = Frequency))

ggplot(data = dataset, mapping = aes(x = CorticalDistance, y = Frequency, color = Analysis)) + 
  geom_point() +
  geom_smooth() +
  facet_wrap(~ ROI, nrow = 2)+
  theme_minimal()

ggplot(data = dataset, mapping = aes(x = Frequency, y = r2, color = Analysis)) + 
  geom_point() +
  geom_smooth() +
  facet_wrap(~ ROI, nrow = 2)

ggplot(data = dataset, mapping = aes(x = Frequency, y = TuningWidth, color = Analysis)) + 
  geom_point() +
  geom_smooth() +
  facet_wrap(~ ROI, nrow = 2) +
  theme_minimal()




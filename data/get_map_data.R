install.packages("malariaAtlas")
library(malariaAtlas)
africadata <- getVecOcc(continent = "Africa")
write.csv(africadata, "data/map_data.csv"; row.names = FALSE)

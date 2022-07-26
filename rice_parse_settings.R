SPECIES.MAP <- list("A"="OSJ", 
                    "P"="OSI", 
                    "N"="OSI", 
                    "D"="OSI", 
                    "E"="OG", 
                    "B"="OG", 
                    "C"="OG", 
                    "S"="OG",
                    "W"="OSJ",
                    "T"="OSI",
                    "K"="OSI",
                    "L"="OSJ")

CONDITIONS <<- list(
  "W" = "Water",
  "S" = "NaCl"
)

TIMES <<- list(
  "1" = 0,
  "2" = 15,
  "3" = 45,
  "4" = 60,
  "5" = 90,
  "6" = 180,
  "7" = 240,
  "8" = 7200,
  "9" = 30,
  "15" = 15,
  "60" = 60,
  "180" = 180,
  "5D" = 7200
)

GENOTYPES <<- list(
  "A" = "Azucena",
  "P" = "Pokkali",
  "N" = "Nona Bokra",
  "D" = "Dhagad deshi",
  "E" = "Ex Wurno",
  "B" = "BG71",
  "C" = "Cuisson 2",
  "S" = "Salifore",
  "W" = "Pandan wangi",
  "T" = "Tadukan",
  "K" = "Kinandang puti",
  "L" = "Palawan"
)

parse.sample.mz2022 <- function(snames) {
  
  mdata <- data.frame(row.names=snames)

  mdata$Genotype <- vapply(snames, function(x) {substr(x, 1, 1)}, character(1))
  mdata$Genotype_Name <- unlist(GENOTYPES[mdata$Genotype])
  
  mdata$Replicate <- as.integer(substring(snames, 2, 2))
  mdata$Condition <- unlist(CONDITIONS[substring(snames, 3, 3)])
  mdata$Time <- unlist(TIMES[substring(snames, 4, 4)])
  
  return(mdata)
}

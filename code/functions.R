# functions for chapter 3 

# list of all functions used in chapter 3 analyses
# call from main.R

# Extract and organize SALMOD output data ####

SALMODextract <- function() {
  
  # A function to 1) extract annual summary data on the freshwater stages 
  # of spring run chinook salmon in Butte Creek for each mgmt scenario, 
  # climate model, and emissions scenario  
  
  # note that summary.out file, the year corresponds to the year - 1, 
  # except for the final year. This goofiness must have to do with the 
  # difference between model year and calendar year. There are two 2099 years. 
  # The first is 2098, and the second 2099]
  
  # Things to calculate
  # fraction of females (males) surviving to spawn each year
  # relationship between number of females surviving to spawn and 
  
  # Definitions
  # Adult Female & Adult Male Entrants - each year 7568 F and 7578 M are seeded
  # of those Adults
  # there is a background mortality rate - can the input value be calculated based on annual data?
  # I think that likely it won't, since all mortality acts at weekly time steps
  
  # separate functions for cumulative mortality stats for :  
  # Spawning Males
  # Spawning Females
  # Egg 
  # Fry
  # library(reshape2)
  # relative path to simulations results: "./data/simulation_results")
  
  sim_path <- file.path(".", "data", "simulation_results")
  year <- 2010:2099
  
  mgmtScenarios <- c("BAU", "ColdWater", "ForeCast", "ForecastError", "NoDiversion", 
                     "RaisePhilbrook", "RaisePhilbrook1710", "ReservoirCover", "Shade")
  climateScenarios <- c("A2","B1")
  gcms <- c("cnrmcm3", "gfdlcm21", "miroc32med", "mpiecham5", "ncarccsm3", "ncarpcm1")
  gcmScenarios <- numeric(length(climateScenarios)*length(gcms))
  for (i in 1:length(climateScenarios)) {
    for (j in 1:length(gcms)) {
      if (i == 1) gcmScenarios[j] <- paste(climateScenarios[i], "_", gcms[j], sep = "")
      if (i == 2) gcmScenarios[6+j] <- paste(climateScenarios[i], "_", gcms[j], sep = "")
    }
  }
  
  for (i in 1:length(mgmtScenarios)) {
    for (j in 1:length(gcmScenarios)) {
      path2sumout <- file.path(sim_path, mgmtScenarios[i], gcmScenarios[j], "summary.out")
      temp <- readLines(path2sumout)
      # temp <- readLines(paste(mgmtScenarios[i], "/", gcmScenarios[j],"/summary.out", sep = ""))
      # (paste(mgmtScenarios[i], "/", gcmScenarios[j],"/summary.out", sep = ""))
      # extract data from temp
      #From Mass Balance 'for Chinook Salmon' table
      
      # Extract lines with the number of Adult females from summary.out and get 
      # numbers of Adult females and Adult female deaths
      ind <- grep("^\\sAdult\\sFemales\\s+[0-9]+\\s.+$", temp)
      
      AF <- as.numeric(gsub("(^\\sAdult\\sFemales\\s+)([0-9]+)(\\s.+$)","\\2" , temp[ind]))
      
      AFdeath <- as.numeric(gsub("(^\\sAdult\\sFemales\\s+)([0-9]+)(\\s+)([0])(\\s+)([0])(\\s+)([0-9]+)(.+$)","\\8" , temp[ind]))
      
      # Extract lines with the number of Adult males from summary.out and get 
      # numbers of Adult males and Adult male deaths
      ind <- grep("^\\sAdult\\sMales\\s+[0-9]+\\s.+$", temp)
      
      AM <- as.numeric(gsub("(^\\sAdult\\sMales\\s+)([0-9]+)(\\s.+$)","\\2" , temp[ind]))
      
      AMdeath <- as.numeric(gsub("(^\\sAdult\\sMales\\s+)([0-9]+)(\\s+)([0])(\\s+)([0])(\\s+)([0-9]+)(.+$)","\\8" , temp[ind]))
      
      # Extract lines with the number of Eggs/Embryos from summary.out and get
      # numbers of entrants and deaths
      ind <- grep("^\\sEggs/Embryos\\s+[0-9]+\\s.+$", temp)
      EE <- as.numeric(gsub("(^\\sEggs/Embryos\\s+)([0-9]+)(\\s.+$)","\\2" , temp[ind]))
      
      EEdeath <- as.numeric(gsub("(^\\sEggs/Embryos\\s+)([0-9]+)(\\s+)([0])(\\s+)([0])(\\s+)([0-9]+)(.+$)","\\8" , temp[ind]))
      
      # Extract lines with the number of Fry from summary.out and get
      # numbers of Graduates   Exiters    Deaths  InStream
      
      ind <- grep("^\\sFry\\s+[0-9]+\\s.+$", temp)
      FryGrad <- as.numeric(gsub("(^\\sFry\\s+)([0]\\s+)([0-9]+)(\\s.+$)","\\3" , temp[ind]))
      FryExit <- as.numeric(gsub("(^\\sFry\\s+)([0]\\s+)([0-9]+)(\\s+)([0-9]+)(.+$)","\\5" , temp[ind]))
      FryDeath <- as.numeric(gsub("(^\\sFry\\s+)([0]\\s+)([0-9]+)(\\s+)([0-9]+)(\\s+)([0-9]+)(.+$)","\\7" , temp[ind]))
      FryInStr <- as.numeric(gsub("(^\\sFry\\s+)([0]\\s+)([0-9]+)(\\s+)([0-9]+)(\\s+)([0-9]+)(\\s+)([0-9]+)(.+$)","\\9" , temp[ind]))
      
      # Put it all together in a dataframe
      MassBal <- data.frame(Year = year, 
                            AFem = AF, 
                            AFemDeath = AFdeath, 
                            AMal = AM,
                            AMalDeath = AMdeath,
                            Eggs = EE,
                            EggsDeath = EEdeath,
                            FryGrad = FryGrad,
                            FryExit = FryExit,
                            FryDeath = FryDeath,
                            FryInStr = FryInStr)
      
      # Cumulative Mortality Statistics
      
      # Find mortality statistics for Adult females
      ind <- grep("\"     Adult Females             \"", temp)  # finds the row with total mortalities
      # the grep indexes all rows where Adult Female mortality
      # data are tabulated in summary.out
      allMortAF <- as.numeric(gsub("(\"          All                  \")(\\s+)([0-9]+)(.+$)", "\\3", temp[ind + 1])) # adding 1 gets total mortality
      baseMortAF <- as.numeric(gsub("(\"          Base                 \")(\\s+)([0-9]+)(.+$)", "\\3", temp[ind + 2])) # adding 2 gets base mortality 
      tempMortAF <- as.numeric(gsub("(\"          Temperature          \")(\\s+)([0-9]+)(.+$)", "\\3", temp[ind + 3])) # adding 3 gets temperature mortality
      
      # Find mortality statistics for Adult males
      ind <- grep("\"     Adult Males               \"", temp)  # finds the row with total mortalities
      # the grep indexes all rows where Adult Male mortality
      # data are tabulated in summary.out
      allMortAM <- as.numeric(gsub("(\"          All                  \")(\\s+)([0-9]+)(.+$)", "\\3", temp[ind + 1])) # adding 1 gets total mortality
      baseMortAM <- as.numeric(gsub("(\"          Base                 \")(\\s+)([0-9]+)(.+$)", "\\3", temp[ind + 2])) # adding 2 gets base mortality 
      tempMortAM <- as.numeric(gsub("(\"          Temperature          \")(\\s+)([0-9]+)(.+$)", "\\3", temp[ind + 3])) # adding 3 gets temperature mortality
      
      # Find mortality statistics for InVivo Eggs
      ind <- grep("\"     InVivo Eggs               \"", temp)  # finds the row with total mortalities
      InVivo <- as.numeric(gsub("(\"     InVivo Eggs               \")(\\s+)([0-9]+)(.+)","\\3" , temp[ind]))
      
      cumMort <- data.frame(allMortAF = allMortAF,
                            baseMortAF = baseMortAF,
                            tempMortAF = tempMortAF,
                            allMortAM = allMortAM,
                            baseMortAM = baseMortAM,
                            tempMortAM = tempMortAM,
                            InVivo = InVivo)
      
      # combine data from temp
      
      sumOutdf <- cbind(MassBal, cumMort)
      sumOutdf[is.na(sumOutdf)] <- 0
      
      # save it
      
      write.table(sumOutdf, file = file.path(sim_path, 
                                             mgmtScenarios[i], 
                                             gcmScenarios[j], 
                                             "SALMODsumOut1.txt"),
                  sep = ",",
                  row.names = FALSE)
    }
  }   
} # end curly bracket for SALMODextract()

SALMODsumSpawnJuv <- function () {
  
  sim_path <- file.path(".", "data", "simulation_results")
  year <- 2010:2099
  
  mgmtScenarios <- c("BAU", "ColdWater", "ForeCast", "ForecastError", "NoDiversion", 
                     "RaisePhilbrook", "RaisePhilbrook1710", "ReservoirCover", "Shade")
  climateScenarios <- c("A2","B1")
  gcms <- c("cnrmcm3", "gfdlcm21", "miroc32med", "mpiecham5", "ncarccsm3", "ncarpcm1")
  gcmScenarios <- numeric(length(climateScenarios)*length(gcms))
  for (i in 1:length(climateScenarios)) {
    for (j in 1:length(gcms)) {
      if (i == 1) gcmScenarios[j] <- paste(climateScenarios[i], "_", gcms[j], sep = "")
      if (i == 2) gcmScenarios[6+j] <- paste(climateScenarios[i], "_", gcms[j], sep = "")
    }
  }
  
  for (i in 1:length(mgmtScenarios)) {
    for (j in 1:length(gcmScenarios)) {
      temp <- readLines(file.path(sim_path, 
                                  mgmtScenarios[i],
                                  gcmScenarios[j],
                                  "summary.out"))
      # print(paste(mgmtScenarios[i], gcmScenarios[j]), sep = " ")
      # Find mortality statistics for Spawning males
      ind <- grep("\"     Spawning Males            \"", temp)  # finds the row with total mortalities
      allMortSM <- as.numeric(gsub("(\"          All                  \")(\\s+)([0-9]+)(.+$)", "\\3", temp[ind + 1]))# adding 1 gets total mortality
      # determine the correct years
      allMortSMyrs <- (as.numeric(gsub("(\"Cumulative Mortality Statistics\")(\\s+)([0-9]{4})", "\\3", temp[ind - 28]  )) - 1)
      
      # Since Salmod summary.out files have 2 2099 and the years in summary.out are offset ahead by one year,
      # and the different egg, juevenile, and spanwing classes are not present in all years, 
      # this means that care must be taken to assign years correctly
      # The preceding code accounts for the offset by one year, but not the goofiness by having summary.out for two 
      # 2099 (one which is 2098)
      # If there are data in the final two years, then need to check out how this will work
      # A) The first case when there's spawning, egg and juvenile data in both of the final two years (the heat doesn't kill fish before they spawn)
      # That just requires reassignment of the latter year to 2099 in the '...yrs' object
      # B) The next case is when there is spawning, egg and juvenile data in one of the final two years, but not the other
      # This means that we have to test if the year from the summary.out file = the year from the previous year in the summary.out file
      # If they ARE EQUAL then assign the last element of '...yrs' = 2099
      # If they ARE NOT EQUAL then assign the last element of '...yrs' = 2098
      # This step requires referencing the summary.out file directly
      # There should not be a problem if there are no data for the final years
      
      yrs <- allMortSMyrs; rm(allMortSMyrs) 
      caseInd <- ind-28
      next2last <- as.numeric(gsub("(\"Cumulative Mortality Statistics\")(\\s+)([0-9]{4})", "\\3",  temp[caseInd[length(caseInd)]-57] ) )
      # print(paste("Spawn males", next2last))
      if (yrs[length(yrs)] == 2098 & yrs[(length(yrs) -1) ] == 2098 ) { # here the last two
        yrs[(length(yrs)) ] <- 2099
      } 
      
      if (!is.na(next2last)) {
        if (next2last == 2099) yrs[length(yrs)] <- 2099 else yrs[length(yrs)] <- 2098
      } else { } # print(paste("fry final year before 2098 or 2099 for ", mgmtScenarios[i], gcms[j])) }
      rm(next2last)
      SMmort <- data.frame(Year = yrs, variable = "allMortSM", value = allMortSM)
      
      # Find mortality statistics for Spawning females
      ind <- grep("\"     Spawning Females          \"", temp)  # finds the row with total mortalities
      allMortSF <- as.numeric(gsub("(\"          All                  \")(\\s+)([0-9]+)(.+$)", "\\3", temp[ind + 1]))# adding 1 gets total mortality
      allMortSFyrs <- as.numeric(gsub("(\"Cumulative Mortality Statistics\")(\\s+)([0-9]{4})", "\\3", temp[ind - 11]  ))-1
      
      yrs <- allMortSFyrs; rm(allMortSFyrs) 
      caseInd <- ind-11
      next2last <- as.numeric(gsub("(\"Cumulative Mortality Statistics\")(\\s+)([0-9]{4})", "\\3",  temp[caseInd[length(caseInd)]-57] ) )
      # print(paste("Spawn females", next2last))
      if (yrs[length(yrs)] == 2098 & yrs[(length(yrs) -1) ] == 2098 ) { # here the last two
        yrs[(length(yrs)) ] <- 2099
      } 
      
      if (!is.na(next2last)) {
        if (next2last == 2099) yrs[length(yrs)] <- 2099 else yrs[length(yrs)] <- 2098
      } else {} # print(paste("fry final year before 2098 or 2099 for ", mgmtScenarios[i], gcms[j])) }
      
      rm(next2last)
      SFmort <- data.frame(Year = yrs, variable = "allMortSF", value = allMortSF )
      
      
      # Find mortality statistics for Eggs/Embryos
      ind <- grep("\"     Eggs/Embryos              \"", temp)  # finds the row with total mortalities
      indm <- ind - 38
      eggMortYr <- (as.numeric(gsub("(\"Cumulative Mortality Statistics\")(\\s+)([0-9]{4})", "\\3", temp[indm]  )) -1)
      if(any(is.na(eggMortYr))) {
        indm[which(is.na(eggMortYr))] <- indm[which(is.na(eggMortYr))]+9
        eggMortYr <- as.numeric(gsub("(\"Cumulative Mortality Statistics\")(\\s+)([0-9]{4})", "\\3", temp[indm]  ))-1
      }
      
      yrs <- eggMortYr; rm(eggMortYr) 
      caseInd <- ind-38
      next2last <- as.numeric(gsub("(\"Cumulative Mortality Statistics\")(\\s+)([0-9]{4})", "\\3",  temp[caseInd[length(caseInd)]-57] ) )
      # print(paste("Eggs", next2last))
      if (yrs[length(yrs)] == 2098 & yrs[(length(yrs) -1) ] == 2098 ) { # here the last two
        yrs[(length(yrs)) ] <- 2099
      } 
      
      if (!is.na(next2last)) {
        if (next2last == 2099) yrs[length(yrs)] <- 2099 else yrs[length(yrs)] <- 2098
      } else {} # print(paste("fry final year before 2098 or 2099 for ", mgmtScenarios[i], gcms[j])) }
      rm(next2last)
      
      allMortEgg <- as.numeric(gsub("(\"          All                  \")(\\s+)([0-9]+)(.+$)", "\\3", temp[ind + 1])) # adding 1 gets total mortality
      baseMortEgg <- as.numeric(gsub("(\"          Base                 \")(\\s+)([0-9]+)(.+$)", "\\3", temp[ind + 2])) # adding 2 gets base mortality 
      tempMortEgg <- as.numeric(gsub("(\"          Temperature          \")(\\s+)([0-9]+)(.+$)", "\\3", temp[ind + 3])) # adding 3 gets temperature mortality
      denMortEgg <- as.numeric(gsub("(\"          Density Dependent    \")(\\s+)([0-9]+)(.+$)", "\\3", temp[ind + 4])) # adding 4 gets density mortality
      incMortEgg <- as.numeric(gsub("(\"          Incubation Loss      \")(\\s+)([0-9]+)(.+$)", "\\3", temp[ind + 5])) # adding 5 gets incubation loss mortality
      superMortEgg <- as.numeric(gsub("(\"          Superimposition      \")(\\s+)([0-9]+)(.+$)", "\\3", temp[ind + 6])) # adding 6 gets superimposition mortality
      lostMortEgg <- as.numeric(gsub("(\"          Lost Eggs            \")(\\s+)([0-9]+)(.+$)", "\\3", temp[ind + 7])) # adding 7 gets lost eggs mortality
      catMortEgg <- as.numeric(gsub("(\"          Catastrophic         \")(\\s+)([0-9]+)(.+$)", "\\3", temp[ind + 8])) # adding 8 gets catastrophic mortality
      
      eggDF <- data.frame(Year = yrs,
                          allMortEgg = allMortEgg,
                          baseMortEgg = baseMortEgg,
                          tempMortEgg = tempMortEgg,
                          denMortEgg = denMortEgg,
                          incMortEgg = incMortEgg,
                          superMortEgg = superMortEgg,
                          lostMortEgg = lostMortEgg,
                          catMortEgg = catMortEgg)
      
      # Find mortality statistics for Fry
      ind <- grep("\"     Fry                       \"", temp)  # finds the row with total mortalities
      indm <- ind - 47
      fryMortYr <- as.numeric(gsub("(\"Cumulative Mortality Statistics\")(\\s+)([0-9]{4})", "\\3", temp[indm]  ))-1
      if(any(is.na(fryMortYr))) {
        indm[which(is.na(fryMortYr))] <- indm[which(is.na(fryMortYr))]+9
        fryMortYr <- as.numeric(gsub("(\"Cumulative Mortality Statistics\")(\\s+)([0-9]{4})", "\\3", temp[indm]  ))-1
      }
      yrs <- fryMortYr; rm(fryMortYr) 
      caseInd <- ind-47
      next2last <- as.numeric(gsub("(\"Cumulative Mortality Statistics\")(\\s+)([0-9]{4})", "\\3",  temp[caseInd[length(caseInd)]-57] ) )
      # print(paste("Fry", next2last))
      if (yrs[length(yrs)] == 2098 & yrs[(length(yrs) -1) ] == 2098 ) { # here the last two
        yrs[(length(yrs)) ] <- 2099
      } 
      
      if (!is.na(next2last)) {
        if (next2last == 2099) yrs[length(yrs)] <- 2099 else yrs[length(yrs)] <- 2098
      } else { } #print(paste("fry final year before 2098 or 2099 for ", mgmtScenarios[i], gcms[j])) }
      
      rm(next2last)
      
      allMortFry <- as.numeric(gsub("(\"          All                  \")(\\s+)([0-9]+)(.+$)", "\\3", temp[ind + 1])) # adding 1 gets total mortality
      baseMortFry <- as.numeric(gsub("(\"          Base                 \")(\\s+)([0-9]+)(.+$)", "\\3", temp[ind + 2])) # adding 2 gets base mortality 
      tempMortFry <- as.numeric(gsub("(\"          Temperature          \")(\\s+)([0-9]+)(.+$)", "\\3", temp[ind + 3])) # adding 3 gets temperature mortality
      denMortFry <- as.numeric(gsub("(\"          Density Dependent    \")(\\s+)([0-9]+)(.+$)", "\\3", temp[ind + 4])) # adding 4 gets density mortality
      habMortFry <- as.numeric(gsub("(\"          Habitat Movement     \")(\\s+)([0-9]+)(.+$)", "\\3", temp[ind + 5])) # adding 3 gets incubation loss mortality
      freshMortFry <- as.numeric(gsub("(\"          Freshet Movement     \")(\\s+)([0-9]+)(.+$)", "\\3", temp[ind + 6])) # adding 3 gets superimposition mortality
      seasMortFry <- as.numeric(gsub("(\"          Seasonal Movement    \")(\\s+)([0-9]+)(.+$)", "\\3", temp[ind + 7])) # adding 3 gets lost eggs mortality
      catMortFry <- as.numeric(gsub("(\"          Catastrophic         \")(\\s+)([0-9]+)(.+$)", "\\3", temp[ind + 8])) # adding 3 gets catastrophic mortality
      
      fryDF <- data.frame(Year = yrs,
                          allMortFry = allMortFry,
                          baseMortFry = baseMortFry,
                          tempMortFry = tempMortFry,
                          denMortFry = denMortFry,
                          habMortFry = habMortFry,
                          freshMortFry = freshMortFry,
                          seasMortFry = seasMortFry,
                          catMortFry = catMortFry) 
      
      # Melt dfs and use a dummy so you can fill in the missing years in 2010:2099
      
      eggDFm <- melt(eggDF, id = c("Year"))
      fryDFm <- melt(fryDF, id = c("Year"))               
      
      # the dummy data 
      fakeout <- data.frame(Year = as.numeric(2010:2099), variable = factor("placeholder"), value = rnorm(1) )
      
      tricky <- rbind(eggDFm, fryDFm, SFmort, SMmort, fakeout)
      
      sumSpawnJuv <- dcast(tricky, Year ~ variable)
      
      sumSpawnJuv$placeholder <- NULL 
      
      sumSpawnJuv[is.na(sumSpawnJuv)] <- 0
      
      # save it
      
      write.table(sumSpawnJuv, file = file.path(sim_path, 
                                                mgmtScenarios[i], 
                                                gcmScenarios[j],
                                                "SALMODsumOut2.txt"),
                  sep = ",",
                  row.names = FALSE)
      
    }
  }
  
  
} # end curly bracket for SALMODsumSpawnJuv


mergeSALMODextract <- function() {
  # Merge SALMODsumOut1.txt with SALMODsumOut2.txt
  sim_path <- file.path(".", "data", "simulation_results")
  # code to run through the files.
  mgmtScenarios <- c("BAU", "ColdWater", "ForeCast", "ForecastError", "NoDiversion", 
                     "RaisePhilbrook", "RaisePhilbrook1710", "ReservoirCover", "Shade")
  climateScenarios <- c("A2","B1")
  gcms <- c("cnrmcm3", "gfdlcm21", "miroc32med", "mpiecham5", "ncarccsm3", "ncarpcm1")
  gcmScenarios <- numeric(length(climateScenarios)*length(gcms))
  
  for (i in 1:length(climateScenarios)) {
    for (j in 1:length(gcms)) {
      if (i == 1) gcmScenarios[j] <- paste(climateScenarios[i], "_", gcms[j], sep = "")
      if (i == 2) gcmScenarios[6+j] <- paste(climateScenarios[i], "_", gcms[j], sep = "")
    }
  }
  
  for (i in 1:length(mgmtScenarios)) {
    for (j in 1:length(gcmScenarios)) {
      temp1 <- read.delim(file.path(sim_path, mgmtScenarios[i], gcmScenarios[j], "SALMODsumOut1.txt"), sep = ",")
      temp2 <- read.delim(file.path(sim_path, mgmtScenarios[i], gcmScenarios[j], "SALMODsumOut2.txt"), sep = ",")
      
      temp <- cbind(temp1, temp2)
      rm(temp1, temp2)
      write.table(temp, file = file.path(sim_path, mgmtScenarios[i], gcmScenarios[j], "SALMODsumOutMerge.txt"),
                  sep = ",",
                  row.names = FALSE)
    }
  }
  
} # end mergeSALMODextract

spawnFemSurvSALMOD <- function() {
  sim_path <- file.path(".", "data", "simulation_results")
  # code to run through the files.
  mgmtScenarios <- c("BAU", "ColdWater", "ForeCast", "ForecastError", "NoDiversion", 
                     "RaisePhilbrook", "RaisePhilbrook1710", "ReservoirCover", "Shade")
  climateScenarios <- c("A2","B1")
  gcms <- c("cnrmcm3", "gfdlcm21", "miroc32med", "mpiecham5", "ncarccsm3", "ncarpcm1")
  gcmScenarios <- numeric(length(climateScenarios)*length(gcms))
  
  # Calc annual predicted survival of spawning females for each water management scenario, climate model and emissions scenario
  for (i in 1:length(mgmtScenarios)) {
    A2 <- B1 <- as.data.frame(matrix(NA, nrow = 90, ncol = length(gcms)))
    for (j in 1:length(climateScenarios)) {
      for (k in 1:length(gcms)) {
        # For each management & emissions scenario combination, store annual predicted survival for each climate model in an 
        # appropriately named object  
        names(A2) <- names(B1) <- gcms
        merged <- read.delim(file.path(sim_path, 
                                       mgmtScenarios[i], 
                                       paste0(climateScenarios[j],"_", gcms[k]), 
                                       "SALMODsumOutMerge.txt"), 
                             sep = ",")[,c("Year", "AFem", "allMortSF")]
        if (j == 1) {
          A2[,k] <- with(merged, allMortSF/AFem)
        } else { 
          B1[,k] <- with(merged, allMortSF/AFem)
        } # end if else conditional statement
      }
    }
    
    # For each management & emissions scenario combination, calculate mean and sd of annual predicted survival across models    
    A2$meanSFsurv <- rowMeans(A2)
    A2$sdSFsurv <- apply(A2, 1, sd)
    B1$meanSFsurv <- rowMeans(B1)
    B1$sdSFsurv <- apply(B1, 1, sd)
    write.table(A2, file = file.path(sim_path, 
                                     mgmtScenarios[i], 
                                     "SALMODspawnFemSurvA2.txt"),
                sep = ",",
                row.names = FALSE)
    write.table(B1, file = file.path(sim_path, 
                                     mgmtScenarios[i], 
                                     "SALMODspawnFemSurvB1.txt"),
                sep = ",",
                row.names = FALSE)
  }
} # end curly bracket spawnFemSurvSALMOD()

bootSRsurvSALMOD <- function (nboot = 10000, seed = 27) {
  # Resample, with replacement, the annual estiamtes from each of the six models, to create N alternate scenarios for 
  # each management & emissions scenario
  sim_path <- file.path(".", "data", "simulation_results")
  set.seed(seed)
  # code to run through the files.
  mgmtScenarios <- c("BAU", "ColdWater", "ForeCast", "ForecastError", "NoDiversion", 
                     "RaisePhilbrook", "RaisePhilbrook1710", "ReservoirCover", "Shade")
  climateScenarios <- c("A2","B1")
  #   gcms <- c("cnrmcm3", "gfdlcm21", "miroc32med", "mpiecham5", "ncarccsm3", "ncarpcm1")
  #   gcmScenarios <- numeric(length(climateScenarios)*length(gcms))
  bootA2 <- bootB1 <- matrix(NA, nrow = 90, ncol = nboot)
  for (i in 1:length(mgmtScenarios)) {
    path2data <- file.path(sim_path, mgmtScenarios[i])
    A2 <- read.delim(file.path(path2data, "SALMODspawnFemSurvA2.txt"), sep = ",") 
    B1 <- read.delim(file.path(path2data, "SALMODspawnFemSurvB1.txt"), sep = ",") 
    for (Nrow in 1:nrow(A2)) {
      bootA2[Nrow,] <- as.numeric(sample(A2[Nrow,1:6], size = nboot, replace = TRUE))
      bootB1[Nrow,] <- as.numeric(sample(B1[Nrow,1:6], size = nboot, replace = TRUE))
    }
    file.create(file.path(path2data, "bootFemSurvA2.txt"))
    file.create(file.path(path2data, "bootFemSurvB1.txt"))
    write.table(bootA2, file = file.path(path2data, "bootFemSurvA2.txt"), append = TRUE)
    write.table(bootB1, file = file.path(path2data, "bootFemSurvB1.txt"), append = TRUE)
  }
} # end curly bracket for bootSRsurvSALMOD()

# op <- par(mfrow = c(2,1));
# matplot(2010:2099, bootA2, type = "l", pch = 19, col = "slateblue");
# lines(2010:2099, rowMeans(bootA2), col = "black", lwd = 3)
# matplot(2010:2099, bootB1, type = "l", pch = 19, col = "slateblue");
# lines(2010:2099, rowMeans(bootB1), col = "black", lwd = 3)
# par(op)

egg2fry2outSurv <- function() {
  sim_path <- file.path(".", "data", "simulation_results")
  # code to run through the files.
  mgmtScenarios <- c("BAU", "ColdWater", "ForeCast", "ForecastError", "NoDiversion", 
                     "RaisePhilbrook", "RaisePhilbrook1710", "ReservoirCover", "Shade")
  climateScenarios <- c("A2","B1")
  gcms <- c("cnrmcm3", "gfdlcm21", "miroc32med", "mpiecham5", "ncarccsm3", "ncarpcm1")
  gcmScenarios <- numeric(length(climateScenarios)*length(gcms))
  
  # Calc annual predicted survival of spawning females for each water management scenario, climate model and emissions scenario
  for (i in 1:length(mgmtScenarios)) {
    A2egg2fry <- A2fryout <- B1egg2fry <- B1fryout <- as.data.frame(matrix(NA, nrow = 90, ncol = length(gcms)))
    for (j in 1:length(climateScenarios)) {
      for (k in 1:length(gcms)) {
        # For each management & emissions scenario combination, store annual predicted survival for each climate model in an 
        # appropriately named object  
        names(A2egg2fry) <- names(A2fryout) <- names(B1egg2fry) <- names(B1fryout) <- gcms
        path_full <- file.path(mgmtScenarios[i], paste0(climateScenarios[j],"_", gcms[k]), "SALMODsumOutMerge.txt")
        merged <- read.delim(file.path(sim_path, path_full), sep = ",")
        if (j == 1) {
          A2egg2fry[,k] <- with(merged, FryGrad/Eggs)
          A2fryout[,k] <- with(merged, FryExit/FryGrad)
        } else { 
          B1egg2fry[,k] <- with(merged, FryGrad/Eggs)
          B1fryout[,k] <- with(merged, FryExit/FryGrad)
        } # end if else conditional statement
      }
    }
    # Write file For each management & emissions scenario combination, calculate mean and sd of annual predicted survival across models    
    write.table(A2egg2fry, file = file.path(sim_path, mgmtScenarios[i], "A2egg2fry.txt"),
                sep = ",",
                row.names = FALSE)
    write.table(B1egg2fry, file = file.path(sim_path, mgmtScenarios[i], "B1egg2fry.txt"),
                sep = ",",
                row.names = FALSE) 
    write.table(A2fryout, file = file.path(sim_path, mgmtScenarios[i], "A2fryout.txt"),
                sep = ",",
                row.names = FALSE)
    write.table(B1fryout, file = file.path(sim_path, mgmtScenarios[i], "B1fryout.txt"),
                sep = ",",
                row.names = FALSE)  
  }
}

SALMODfecundity <- function() {
  sim_path <- file.path(".", "data", "simulation_results")
  # code to run through the files.
  mgmtScenarios <- c("BAU", "ColdWater", "ForeCast", "ForecastError", "NoDiversion", 
                     "RaisePhilbrook", "RaisePhilbrook1710", "ReservoirCover", "Shade")
  climateScenarios <- c("A2","B1")
  gcms <- c("cnrmcm3", "gfdlcm21", "miroc32med", "mpiecham5", "ncarccsm3", "ncarpcm1")
  gcmScenarios <- numeric(length(climateScenarios)*length(gcms))
  
  for (i in 1:length(mgmtScenarios)) {
    A2fecundity <- B1fecundity <- as.data.frame(matrix(NA, nrow = 90, ncol = length(gcms)))
    for (j in 1:length(climateScenarios)) {
      for (k in 1:length(gcms)) {
        # For each management & emissions scenario combination, store annual predicted survival for each climate model in an 
        # appropriately named object  
        names(A2fecundity)  <- names(B1fecundity) <- gcms
        full_path <- file.path(mgmtScenarios[i], paste0(climateScenarios[j],"_", gcms[k]), "SALMODsumOutMerge.txt")
        merged <- read.delim(file.path(sim_path, full_path), sep = ",")
        if (j == 1) {
          A2fecundity[,k] <- with(merged, Eggs/allMortSF)
          
        } else { 
          B1fecundity[,k] <- with(merged, Eggs/allMortSF)
        } # end if else conditional statement
      }
    }
    # Write file For each management & emissions scenario combination, calculate mean and sd of annual predicted survival across models    
    
    write.table(A2fecundity, file = file.path(sim_path, mgmtScenarios[i], "A2fecundity.txt"),
                sep = ",",
                row.names = FALSE)
    write.table(B1fecundity, file = file.path(sim_path, mgmtScenarios[i], "B1fecundity.txt"),
                sep = ",",
                row.names = FALSE)  
  }
} # end courly bracker for SALMODfecundity()

SALMODrelMort <- function() {
  sim_path <- file.path(".", "data", "simulation_results")
  # code to run through the files.
  mgmtScenarios <- c("BAU", "ColdWater", "ForeCast", "ForecastError", "NoDiversion", 
                     "RaisePhilbrook", "RaisePhilbrook1710", "ReservoirCover", "Shade")
  climateScenarios <- c("A2","B1")
  gcms <- c("cnrmcm3", "gfdlcm21", "miroc32med", "mpiecham5", "ncarccsm3", "ncarpcm1")
  gcmScenarios <- numeric(length(climateScenarios)*length(gcms))
  
  for (i in 1:length(mgmtScenarios)) {
    A2relEggMortDD <- B1relEggMortDD <- A2relFryMortDD <- B1relFryMortDD <- as.data.frame(matrix(NA, nrow = 90, ncol = length(gcms)))
    A2relEggMortTemp <- B1relEggMortTemp <- A2relFryMortTemp <- B1relFryMortTemp <- as.data.frame(matrix(NA, nrow = 90, ncol = length(gcms)))
    for (j in 1:length(climateScenarios)) {
      for (k in 1:length(gcms)) {
        # For each management & emissions scenario combination, store annual predicted survival for each climate model in an 
        # appropriately named object  
        names(A2relEggMortDD)  <- names(B1relEggMortDD) <- names(A2relFryMortDD)  <- names(B1relFryMortDD) <- gcms
        names(A2relEggMortTemp)  <- names(B1relEggMortTemp) <- names(A2relFryMortTemp)  <- names(B1relFryMortTemp) <- gcms
        full_path <- file.path(mgmtScenarios[i], paste0(climateScenarios[j],"_", gcms[k]), "SALMODsumOutMerge.txt")
        merged <- read.delim(file.path(sim_path, full_path), sep = ",")
        if (j == 1) {
          A2relEggMortDD[,k] <- with(merged, denMortEgg/allMortEgg)
          A2relFryMortDD[,k] <- with(merged, denMortFry/allMortFry)
          A2relEggMortTemp[,k] <- with(merged, tempMortEgg/allMortEgg)
          A2relFryMortTemp[,k] <- with(merged, tempMortFry/allMortFry)
        } else { 
          B1relEggMortDD[,k] <- with(merged, denMortEgg/allMortEgg)
          B1relFryMortDD[,k] <- with(merged, denMortFry/allMortFry)
          B1relEggMortTemp[,k] <- with(merged, tempMortEgg/allMortEgg)
          B1relFryMortTemp[,k] <- with(merged, tempMortFry/allMortFry)
        } # end if else conditional statement
      }
    }
    # Write file For each management & emissions scenario combination, calculate mean and sd of annual predicted survival across models    
    write.table(A2relEggMortDD, file = file.path(sim_path, mgmtScenarios[i], "A2EggMortDD.txt"),
                sep = ",",
                row.names = FALSE)
    write.table(A2relEggMortTemp, file = file.path(sim_path, mgmtScenarios[i], "A2EggMortTemp.txt"),
                sep = ",",
                row.names = FALSE)
    write.table(A2relFryMortDD, file = file.path(sim_path, mgmtScenarios[i], "A2FryMortDD.txt"),
                sep = ",",
                row.names = FALSE)
    write.table(A2relFryMortTemp, file = file.path(sim_path, mgmtScenarios[i], "A2FryMortTemp.txt"),
                sep = ",",
                row.names = FALSE)
    write.table(B1relEggMortDD, file = file.path(sim_path, mgmtScenarios[i], "B1EggMortDD.txt"),
                sep = ",",
                row.names = FALSE) 
    write.table(B1relEggMortTemp, file = file.path(sim_path, mgmtScenarios[i], "B1EggMortTemp.txt"),
                sep = ",",
                row.names = FALSE)  
    write.table(B1relFryMortDD, file = file.path(sim_path, mgmtScenarios[i], "B1FryMortDD.txt"),
                sep = ",",
                row.names = FALSE)  
    write.table(B1relFryMortTemp, file = file.path(sim_path, mgmtScenarios[i], "B1FryMortTemp.txt"),
                sep = ",",
                row.names = FALSE)  
  }
  
}  # end curly bracket for SALMODrelMort()

# Plotting SALMOD output ####

plotEggsPerSpawner <- function() {
  sim_path <- file.path(".", "data", "simulation_results")
  # code to run through the files.
  mgmtScenarios <- c("BAU", "ColdWater", "ForeCast", "ForecastError", "NoDiversion", 
                     "RaisePhilbrook", "RaisePhilbrook1710", "ReservoirCover", "Shade")
  climateScenarios <- c("A2","B1")
  gcms <- c("cnrmcm3", "gfdlcm21", "miroc32med", "mpiecham5", "ncarccsm3", "ncarpcm1")
  gcmScenarios <- numeric(length(climateScenarios)*length(gcms))
  
  for (i in 1:length(mgmtScenarios)) {
    for (j in 1:length(climateScenarios)) {
      temp <- read.delim(file.path(sim_path, mgmtScenarios[i], paste0(climateScenarios[j], "fecundity.txt")), 
                         sep = ",") 
      if (!dir.exists(file.path("output", mgmtScenarios[i]))) { 
        dir.create(file.path("output", mgmtScenarios[i])) }
      
      pdf(file.path("output", mgmtScenarios[i], paste0(climateScenarios[j],"annFecundity.pdf")), 
          width=16/2.54, 
          height=10/2.54, 
          pointsize=10)
      par(mgp=c(2.2,0.45,0), tcl=-0.4, mar=c(3.3,3.6,1.1,1.1))
      matplot(2010:2099,temp, 
              type = "p",
              pch = 19,
              col = "black",
              main = paste("Annual Fecundity for", mgmtScenarios[i], climateScenarios[j], sep = " "),
              ylab = "Eggs per Female",
              xlab = "Year")
      lines(2010:2099, rowMeans(temp, na.rm = TRUE), col = "red", lwd = 2)
      dev.off()
    }
  }
  
} # end curly bracker for plotEggsPerSpawner()

plotMortEggsFryDD <- function() {
  sim_path <- file.path(".", "data", "simulation_results")
  # code to run through the files.
  mgmtScenarios <- c("BAU", "ColdWater", "ForeCast", "ForecastError", "NoDiversion", 
                     "RaisePhilbrook", "RaisePhilbrook1710", "ReservoirCover", "Shade")
  climateScenarios <- c("A2","B1")
  gcms <- c("cnrmcm3", "gfdlcm21", "miroc32med", "mpiecham5", "ncarccsm3", "ncarpcm1")
  gcmScenarios <- numeric(length(climateScenarios)*length(gcms))
  
  for (i in 1:length(mgmtScenarios)) {
    for (j in 1:length(climateScenarios)) {
      densE <- read.delim(file.path(sim_path, mgmtScenarios[i], paste0(climateScenarios[j], "EggMortDD.txt")), sep = ",") 
      densF <- read.delim(file.path(sim_path, mgmtScenarios[i], paste0(climateScenarios[j], "FryMortDD.txt")), sep = ",") 
      tempE <- read.delim(file.path(sim_path, mgmtScenarios[i], paste0(climateScenarios[j], "EggMortTemp.txt")), sep = ",") 
      tempF <- read.delim(file.path(sim_path, mgmtScenarios[i], paste0(climateScenarios[j], "FryMortTemp.txt")), sep = ",") 
      # plot density dependent mortality
      # Eggs
      if (!dir.exists(file.path("output", mgmtScenarios[i]))) { 
        dir.create(file.path("output", mgmtScenarios[i])) }
      
      pdf(file.path("output", mgmtScenarios[i], paste0(climateScenarios[j],"EggMortDensDep.pdf")), 
          width=16/2.54, 
          height=10/2.54, 
          pointsize=10)

      par(mgp=c(2.2,0.45,0), tcl=-0.4, mar=c(3.3,3.6,1.1,1.1))
      matplot(2010:2099,densE, 
              type = "p",
              pch = 19,
              col = "black",
              main = paste("Density dependent related egg mortality", mgmtScenarios[i], climateScenarios[j], sep = " "),
              ylab = "Fraction of Egg Mortality from Density Dependence",
              xlab = "Year")
      lines(2010:2099, rowMeans(densE, na.rm = TRUE), col = "red", lwd = 2)
      dev.off()
      # Fry
      pdf(file.path("output", mgmtScenarios[i], paste0(climateScenarios[j],"FryMortDensDep.pdf")), 
          width=16/2.54, 
          height=10/2.54, 
          pointsize=10)
      par(mgp=c(2.2,0.45,0), tcl=-0.4, mar=c(3.3,3.6,1.1,1.1))
      matplot(2010:2099,densF, 
              type = "p",
              pch = 19,
              col = "black",
              main = paste("Density dependent related fry mortality", mgmtScenarios[i], climateScenarios[j], sep = " "),
              ylab = "Fraction of Fry Mortality from Density Dependence",
              xlab = "Year")
      lines(2010:2099, rowMeans(densF, na.rm = TRUE), col = "red", lwd = 2)
      dev.off()
      
      # Plot temperature related mortality
      # Eggs
      pdf(file.path("output", mgmtScenarios[i], paste0(climateScenarios[j],"EggMortTemperature.pdf")), 
          width=16/2.54, 
          height=10/2.54, 
          pointsize=10)
      par(mgp=c(2.2,0.45,0), tcl=-0.4, mar=c(3.3,3.6,1.1,1.1))
      matplot(2010:2099,tempE, 
              type = "p",
              pch = 19,
              col = "black",
              main = paste("Temperature related egg mortality\n", mgmtScenarios[i], climateScenarios[j], sep = " "),
              ylab = "Fraction of Egg Mortality from Temperature",
              xlab = "Year")
      lines(2010:2099, rowMeans(tempE, na.rm = TRUE), col = "red", lwd = 2)
      dev.off()
      # Fry
      pdf(file.path("output", mgmtScenarios[i], paste0(climateScenarios[j],"FryMortTemperature.pdf")), 
          width=16/2.54, 
          height=10/2.54, 
          pointsize=10)
      par(mgp=c(2.2,0.45,0), tcl=-0.4, mar=c(3.3,3.6,1.1,1.1))
      matplot(2010:2099,tempF, 
              type = "p",
              pch = 19,
              col = "black",
              main = paste("Temperature related fry mortality\n", mgmtScenarios[i], climateScenarios[j], sep = " "),
              ylab = "Fraction of Fry Mortality from Temperature",
              xlab = "Year")
      lines(2010:2099, rowMeans(tempF, na.rm = TRUE), col = "red", lwd = 2)
      dev.off()
    }
  }
  
} # end curly bracket for plotMortEggsFryDD

plotMortEggsFryTemp <- function() {
  sim_path <- file.path(".", "data", "simulation_results")
  # code to run through the files.
  mgmtScenarios <- c("BAU", "ColdWater", "ForeCast", "ForecastError", "NoDiversion") 
  #"RaisePhilbrook", "RaisePhilbrook1710", "ReservoirCover", "Shade")
  climateScenarios <- c("A2","B1")
  gcms <- c("cnrmcm3", "gfdlcm21", "miroc32med", "mpiecham5", "ncarccsm3", "ncarpcm1")
  gcmScenarios <- numeric(length(climateScenarios)*length(gcms))
  
  for (i in 1:length(mgmtScenarios)) {
    for (j in 1:length(climateScenarios)) {
      temp <- read.delim(file.path(sim_path, mgmtScenarios[i], paste0(climateScenarios[j], "fecundity.txt")), sep = ",") 
    
      pdf(file.path("output", mgmtScenarios[i], paste0(climateScenarios[j], "annFecundity.pdf")), 
          width=16/2.54, 
          height=10/2.54, 
          pointsize=10)
      par(mgp=c(2.2,0.45,0), tcl=-0.4, mar=c(3.3,3.6,1.1,1.1))
      matplot(2010:2099,temp, 
              type = "p",
              pch = 19,
              col = "black",
              main = paste("Annual Fecundity for", mgmtScenarios[i], climateScenarios[j], sep = " "),
              ylab = "Eggs per Female",
              xlab = "Year")
      lines(2010:2099, rowMeans(temp, na.rm = TRUE), col = "red", lwd = 2)
      dev.off()
    }
  }
  
} # end curly bracket for plotMortEggFryTemp


plotSpawnEggFate <- function() {
  
  egg <- read.delim(file.path(mgmtScenarios[i], paste0(climateScenarios[j],"_", gcms[k]), "SALMODsumOutMerge.txt"), sep = ",")#[,c("Eggs", "EggsDeath", "InVivo", "allMortEgg", "baseMortEgg", "tempMortEgg", "incMortEgg", "superMortEgg", "lostMortEgg", "FryGrad")]
  
  junk <- with(egg, data.frame(allMortSF=allMortSF, Eggs=Eggs, InVivo=InVivo, EpF = Eggs/allMortSF, EIpF = (Eggs + InVivo)/allMortSF) )
  
} # end curly bracket for plotSpawnEggFate()


plotPrespawnMort <- function() {
  
  sim_path <- file.path(".", "data", "simulation_results")
  mgmtScenarios <- c("BAU", "NoDiversion") #"ColdWater", "ForeCast", "ForecastError", "NoDiversion", 
  #"RaisePhilbrook", "RaisePhilbrook1710", "ReservoirCover", "Shade")
  climateScenarios <- c("A2","B1")
  gcms <- c("cnrmcm3", "gfdlcm21", "miroc32med", "mpiecham5", "ncarccsm3", "ncarpcm1")
  gcmScenarios <- numeric(length(climateScenarios)*length(gcms))
  for (i in 1:length(mgmtScenarios)) {
    for (j in 1:length(climateScenarios)) {
      mortFW <- 
        for (k in 1:length(gcms)) {
          salmod_out <- read.delim(file.path(sim_path, mgmtScenarios[i], paste0(climateScenarios[j],"_", gcms[k]), "SALMODsumOutMerge.txt"), 
                                   sep = ",")
          mortFW <- subset(salmod_out, select = c("AFem", "allMortSF", "Eggs", "FryGrad" ,"FryExit"))
          mortFW$gcms <- gcms[k]
          
          if(k == 1) mortFWgcms <- mortFW else mortFWgcms <- rbind(mortFWgcms, mortFW)
        } # end of k loop
      
      mortFWgcms$year <- 2010:2099
      mortFWgcms <- transform(mortFWgcms, PreSpawn = allMortSF/AFem)
      mortFWgcms <- transform(mortFWgcms, Egg = FryGrad/Eggs)
      mortFWgcms <- transform(mortFWgcms, Fry = FryExit/FryGrad)
      
      ann_meanFWps <- ddply(mortFWgcms, .(year), summarise, meanPreSpawn = mean(PreSpawn, na.rm = TRUE))
      ann_meanFWegg <- ddply(mortFWgcms, .(year), summarise, meanEgg = mean(Egg, na.rm = TRUE))
      ann_meanFWfry <- ddply(mortFWgcms, .(year), summarise, meanFry = mean(Fry, na.rm = TRUE))
      
      p1 <- ggplot(data = mortFWgcms, aes(x = year, y = PreSpawn)) +
        geom_point() +
        geom_line(data = ann_meanFWps, aes(x = year, y = meanPreSpawn), colour = "red") +
        theme_bw() + 
        ggtitle(paste("PreSpawn survival for ", mgmtScenarios[i], " water mgmt scenario & \n", climateScenarios[j], " emissions scenario", sep = ""))
      
      p2 <- ggplot(data = mortFWgcms, aes(x = year, y = Egg)) +
        geom_point() +
        geom_line(data = ann_meanFWegg, aes(x = year, y = meanEgg), colour = "red") +
        theme_bw() + 
        ggtitle(paste("Egg survival for ", mgmtScenarios[i], " water mgmt scenario & \n", climateScenarios[j], " emissions scenario", sep = ""))
      
      p3 <- ggplot(data = mortFWgcms, aes(x = year, y = Fry)) +
        geom_point() +
        geom_line(data = ann_meanFWfry, aes(x = year, y = meanFry), colour = "red") +
        theme_bw() + 
        ggtitle(paste("Fry survival for ", mgmtScenarios[i], " water mgmt scenario & \n", climateScenarios[j], " emissions scenario", sep = ""))

      pdf(file = file.path("output", mgmtScenarios[i], paste0(climateScenarios[j],"_PSsurv.pdf")),
          width = 10, height = 6) 
      print(p1)
      dev.off()
      
      pdf(file = file.path("output", mgmtScenarios[i], paste0(climateScenarios[j],"_Eggsurv.pdf")),
          width = 10, height = 6) 
      print(p2)
      dev.off()
      
      pdf(file = file.path("output", mgmtScenarios[i], paste0(climateScenarios[j],"_Frysurv.pdf")),
          width = 10, height = 6)
      print(p3)
      dev.off()
      
    }
  }
  
}

# Population vulnerability modelling ####

# Determine the delta_e (fraction of spawning females spawning at age-3)
# returns a plot of N3/(N3+N4) vs delta_e for a given age-3 survival
# Also returns the de for the the wanted fraction (N3/(N3+N4)) and age-3 survival
calc_de <- function(surv3, wanted_frac, plot_wanted = "no") {
  # For Butte Creek chinook salmon this function calculates the delta_e
  # value that returns the "wanted fraction", which is the fraction of 
  # spawners that return at age three assuming equilibrium conditions
  # for a given/assumed value of ocean survival from age-3 to age-4 (surv3)
  # The number of age-3 spawners at EQ is: surv1*surv2*delta_e
  # The number of age-3 spawners at EQ is: surv1*surv2*(1-delta_e)*surv3
  # 
  # Limited data from coded wire tag returns of spawning fish suggest that 
  # BC spring-run chinook spawn at ages 3 and 4, with sometimes as few as 20% 
  # spawning at age-3 in a given year to as much as 80% spawning at age-3.
  #
  #
  # It also returns a plot of the fraction of age-3 and age-4 to total spawnwers
  # [age-3 + age-4 spawners] vs. delta_e
  
  n <- 100
  r3_S <- r4_S <- numeric(n)
  de <- (1:n)/n
  
  for (i in 1:n) {
    r3_S[i] <- de[i]/(de[i] + surv3*(1-de[i]))
    r4_S[i] <- (surv3*(1-de[i]))/(de[i] + surv3*(1-de[i]))
  }
  
  if (plot_wanted == "yes") {
    plot(de, r3_S, type = "l", col = "slateblue", lwd = 3,
         xlab = "Fraction spawning at age-3 (delta_e)",
         ylab = "Number of Age-a Spawners/Total Number of Spawners")
    lines(de, r4_S, type = "l", col = "red", lwd = 3)
    abline(h=0.5, lty = 2)
    abline(v=0.5, lty = 2)
    title(paste("Age-3 ocean survival = ", surv3, sep = ""))
    legend("right", c("a=3","a=4"), # puts text in the legend 
           lty=c(1,1), # gives the legend appropriate symbols (lines)
           lwd=c(2.5,2.5),
           col=c("slateblue","red"),
           cex = .7,
           bty = "n") # gives the legend lines the correct color and width
  }
  
  desired_de <- wanted_frac*surv3/(wanted_frac*(surv3 - 1) + 1)
  #print("The delta_e that you want is:")
  #print(desired_de)
  return(desired_de)
} # end bracket of calc_de()

# calculate spawners per recruit (1/spr = replacement line, baseline/floor for alpha)
SPR_srcs <- function(surv1, surv2, surv3, delta_e, survPS = 1) {
  SPR <- surv1*surv2*delta_e*survPS + surv1*surv2*surv3*(1-delta_e)*survPS
  return(SPR)
} # end bracket for SPR_scrs

# plot how SPR changes with early ocean survival and fraction of age-3 spawners
plotSPR_de_surv1 <- function (sprDF) {
  ggplot(sprDF, aes(x = delta_e, y = SPR, colour = factor(surv1))) +
    geom_point(size = 3) +
    scale_colour_manual(values = c("black","grey20", "grey40", "grey70"), name = "Early Ocean \nSurvival") +
    ylab("Spawners Per Recruit") +
    xlab("Fraction of spawners returning at age-3") +
    xlim(c(0, 1)) +
    ylim(c(0, 0.08)) +
    theme_bw() +
    theme(panel.grid.major = theme_line(colour = NA), panel.grid.minor = theme_line(colour = NA)) 
}


# function to calculate beta of Beverton-Holt (where alpha/beta = max recruits)
# based on a given alpha, equilibrium spawning stock size (EQsp) and SPR level

calc_beta <-function(alpha, EQsp, spr) {
  # For a given SPR, alpha and EQsp, we can calculate beta, for Butte Creek 
  # spring run chinook salmon
  # To calculate beta for multiple alpha values, just pass a vecter of alphas
  # to the alpha argument
  #if (length(EQsp) == 1 & all(EQsp == "infinity"))  beta <- 0 else beta <- (alpha * spr - 1)/EQsp
  beta <- (alpha * spr - 1)/EQsp

  if (any(beta < 0)) warning("beta is negative, that's not good - population parameters not sustainable?.?.?")
  return(beta)
}

# function to plot SR relationships investigated in population model

plotStockRecruit <- function(alpha, beta, EQsp, delta_e, spr, make_plot) {
  
  spawners <- seq(0, max(EQsp)*2, len = 200)
  
  if (length(alpha) > 1 & length(EQsp) > 1) stop("Only alpha or beta (EQsp) can be varied in a function call")
  if (make_plot == TRUE) {
    if (length(alpha) > 5) warning("only 1 to 5 values of alpha may be plotted at a time")
    if (length(EQsp) > 5) warning("only 1 to 5 values of EQsp may be plotted at a time")
  }
  #if (is.empty(alpha) | is.empty(beta)) stop("need to enter alpha and beta")
  
  # Set plotting colors
  if (length(alpha)  == 2 | length(EQsp) == 2) cols <- c("blue", "red")
  if (length(alpha)  == 3 | length(EQsp) == 3) cols <- c("blue", "black", "red")
  if (length(alpha)  == 4 | length(EQsp) == 4) cols <- c("blue", "dodgerblue", "plum", "red")
  if (length(alpha)  == 5 | length(EQsp) == 5) cols <- c("blue", "dodgerblue", "black", "plum",  "red")
  if (length(alpha) == 1 & length(EQsp) == 1) cols <- "black" 
  #   if (length(alpha) > 5) cols <- rainbow(n)
  #   if (length(EQsp) > 5) cols <- rainbow(n)
  # makes a plot of underlying SR relationship for population projections
  # if only one alpha and one beta
  
  if (length(alpha) > 1 & length(EQsp) == 1) {
    recruits <- matrix(0, nrow = length(spawners), ncol = length(alpha))
    alpha_names <- list()
    for (i in 1:length(spawners)) {
      for (j in 1:length(alpha)) {
        recruits[i,j] <- alpha[j]*spawners[i]/(1 + beta[j]*spawners[i])
        alpha_names[j] <- paste("alpha_",toString(signif(alpha[j], 2)), sep = "")
      }
    }
    
    # make and melt df
    sr <- data.frame(spawners, recruits)
    names(sr)[2:(length(alpha)+1)] <- alpha_names
    srm <- melt(sr, id.vars = "spawners")
    names(srm)[2:3] <- c("Alpha", "recruits")
    spawn4spr <- (alpha*spr - 1)/beta
    
    p <- ggplot(data = srm, aes(spawners, recruits, colour = Alpha) ) +
      geom_line() +
      scale_colour_manual(values = rev(cols)) +
      theme_bw() +
      xlab("Spawners") +
      ylab("Recruits") +
      theme(panel.grid.major = element_line(colour = NA), panel.grid.minor = element_line(colour = NA)) + 
      theme(legend.position = "right") +
      geom_abline(intercept = 0, slope = 1/spr, colour = "black", linetype = 2) +
      geom_point(data = data.frame(x = spawn4spr, y = alpha*spawn4spr/(1+beta*spawn4spr)), aes(x = x, y = y), colour = "black", shape = 22, size = 4) +
      scale_x_continuous(expand=c(0,0)) + 
      ylim(c(0, max(recruits)*1.1))
  }
  
  if (length(EQsp) > 1 & length(alpha) == 1) {
    recruits <- matrix(0, nrow = length(spawners), ncol = length(EQsp))
    beta_names <- list()
    for (i in 1:length(spawners)) {
      for (j in 1:length(beta)) {
        recruits[i,j] <- alpha*spawners[i]/(1 + beta[j]*spawners[i])
        beta_names[j] <- paste("EQsp_",toString(EQsp[j]), sep = "")
      }
    }
    
    # make and melt df
    sr <- data.frame(spawners, recruits)
    names(sr)[2:(length(beta)+1)] <- beta_names
    srm <- melt(sr, id.vars = "spawners")
    names(srm)[2:3] <- c("Beta", "recruits")
    spawn4spr <- (alpha*spr - 1)/beta
    
    p <- ggplot(data = srm, aes(spawners, recruits, colour = Beta) ) +
      geom_line() +
      scale_colour_manual(values = rev(cols)) +
      theme_bw() +
      xlab("Spawners") +
      ylab("Recruits") +
      theme(panel.grid.major = element_line(colour = NA), panel.grid.minor = element_line(colour = NA)) + 
      theme(legend.position = "right") +
      geom_abline(intercept = 0, slope = 1/spr, colour = "black", linetype = 2) +
      geom_point(data = data.frame(x = spawn4spr, y = alpha*spawn4spr/(1+beta*spawn4spr)), aes(x = x, y = y), colour = "black", shape = 22, size = 4) +
      scale_x_continuous(expand=c(0,0)) + 
      ylim(c(0, max(recruits)*1.1))
    
  }
  
  if (length(beta) == 1 & length(alpha) == 1) {
    recruits <- alpha*spawners/(1 + beta*spawners)
    # make df for ggplot call
    sr <- data.frame(spawners, recruits)
    
    # Spawning stock for given SPR
    spawn4spr <- (alpha*spr - 1)/beta
    #spawn4spr == EQsp
    #all.equal(EQsp, spawn4spr) # nearly equal, there's a slight rounding error
    p <- ggplot(data = sr, aes(spawners, recruits) ) +
      geom_line(colour = "slateblue") +
      theme_bw() +
      xlab("Spawners") +
      ylab("Recruits") +
      theme(panel.grid.major = element_line(colour = NA), panel.grid.minor = element_line(colour = NA)) + 
      theme(legend.position = "right") +
      geom_abline(intercept = 0, slope = 1/spr, colour = "red", line_type = 2) +
      geom_point(data = data.frame(x = spawn4spr, y = alpha*spawn4spr/(1+beta*spawn4spr)), aes(x = x, y = y), colour = "black", shape = 22, size = 4) +
      scale_x_continuous(expand=c(0,0)) + 
      ylim(c(0, max(sr$recruits)*1.05)) 
  } 
  if (make_plot == TRUE) {
    print(p)
  }
} # end of plotStockRecruit()

plotStockSpEsc <- function(alpha, beta, EQsp, delta_e, surv1, surv2, surv3, spr, make_plot) {
  
  spawners <- seq(0, max(EQsp)*2, len = 200)
  
  if (length(alpha) > 1 & length(EQsp) > 1) stop("Only alpha or EQsp can be varied in a function call")
  if (make_plot == TRUE) {
    if (length(alpha) > 5) stop("only 1 to 5 values of alpha may be plotted at a time")
    if (length(EQsp) > 5) stop("only 1 to 5 values of EQsp may be plotted at a time")
  }
  #if (is.empty(alpha) | is.empty(beta)) stop("need to enter alpha and beta")
  
  # Set plotting colors
  if (length(alpha)  == 2 | length(EQsp) == 2) cols <- c("blue", "red")
  if (length(alpha)  == 3 | length(EQsp) == 3) cols <- c("blue", "black", "red")
  if (length(alpha)  == 4 | length(EQsp) == 4) cols <- c("blue", "dodgerblue", "plum", "red")
  if (length(alpha)  == 5 | length(EQsp) == 5) cols <- c("blue", "dodgerblue", "black", "plum",  "red")
  if (length(alpha) == 1 & length(EQsp) == 1) cols <- "black"  
  
  # makes a plot of underlying SR relationship for population projections
  # if only one alpha and one beta
  
  if (length(alpha) > 1 & length(EQsp) == 1) {
    recruits <- matrix(0, nrow = length(spawners), ncol = length(alpha))
    alpha_names <- list()
    for (i in 1:length(spawners)) {
      for (j in 1:length(alpha)) {
        recruits[i,j] <- alpha[j]*spawners[i]/(1 + beta[j]*spawners[i])
        alpha_names[j] <- paste("alpha_",toString(signif(alpha[j], 2)), sep = "")
      }
    }
    # make and melt df
    sr <- data.frame(spawners, recruits)
    names(sr)[2:(length(alpha)+1)] <- alpha_names
    srm <- melt(sr, id.vars = "spawners")
    srm$spEsc <- with(srm, value*surv1*surv2*delta_e + value*surv1*surv2*surv3*(1-delta_e))
    names(srm)[2:3] <- c("Alpha", "recruits")
    
    p <- ggplot(data = srm, aes(spawners, spEsc, colour = Alpha) ) +
      geom_line() +
      scale_colour_manual(values = rev(cols)) +
      theme_bw() +
      xlab("Spawners") +
      ylab("Recruits: Returning Spawners") +
      theme(panel.grid.major = element_line(colour = NA), panel.grid.minor = element_line(colour = NA)) + 
      theme(legend.position = "right") +
      geom_abline(intercept = 0, slope = 1, colour = "black", linetype = 2) +
      #geom_point(aes(x = EQsp, y = EQsp), colour = "black", shape = 22, size = 4) +
      scale_x_continuous(expand=c(0,0)) + 
      ylim(c(0, max(srm$spEsc)*1.1))
  }
  
  if (length(EQsp) > 1 & length(alpha) == 1) {
    recruits <- matrix(0, nrow = length(spawners), ncol = length(EQsp))
    beta_names <- list()
    for (i in 1:length(spawners)) {
      for (j in 1:length(beta)) {
        recruits[i,j] <- alpha*spawners[i]/(1 + beta[j]*spawners[i])
        beta_names[j] <- paste("EQsp_",toString(EQsp[j]), sep = "")
      }
    }
    
    # make and melt df
    sr <- data.frame(spawners, recruits)
    names(sr)[2:(length(beta)+1)] <- beta_names
    srm <- melt(sr, id.vars = "spawners")
    srm$spEsc <- with(srm, value*surv1*surv2*delta_e + value*surv1*surv2*surv3*(1-delta_e))
    names(srm)[2:3] <- c("Beta", "recruits")
    
    p <- ggplot(data = srm, aes(spawners, spEsc, colour = Beta) ) +
      geom_line() +
      scale_colour_manual(values = rev(cols)) +
      theme_bw() +
      xlab("Spawners") +
      ylab("Recruits: Returning Spawners") +
      theme(panel.grid.major = element_line(colour = NA), panel.grid.minor = element_line(colour = NA)) + 
      theme(legend.position = "right") +
      geom_abline(intercept = 0, slope = 1, colour = "black", linetype = 2) +
      #geom_point(aes(x = EQsp, y = EQsp), colour = "black", shape = 22, size = 4) +
      xlim(c(0, max(srm$spEsc)*1.1)) + 
      ylim(c(0, max(srm$spEsc)*1.1))
    
  }
  
  if (length(beta) == 1 & length(alpha) == 1) {
    recruits <- alpha*spawners/(1 + beta*spawners)
    # make df for ggplot call
    sr <- data.frame(spawners, recruits)
    sr$spEsc <- with(sr, recruits*surv1*surv2*delta_e + recruits*surv1*surv2*surv3*(1-delta_e))
    # Spawning stock for given SPR
    spawn4spr <- (alpha*spr - 1)/beta
    #spawn4spr == EQsp
    #all.equal(EQsp, spawn4spr) # nearly equal, there's a slight rounding error
    p <- ggplot(data = sr, aes(spawners, spEsc) ) + 
      geom_line(colour = "slateblue") +
      theme_bw() +
      xlab("Spawners") +
      ylab("Recruits: Returning Spawners") +
      theme(panel.grid.major = element_line(colour = NA), panel.grid.minor = element_line(colour = NA)) + 
      theme(legend.position = "right") +
      geom_abline(intercept = 0, slope = 1, colour = "red", line_type = 2) +
      #geom_point(aes(x = EQsp, y = EQsp), colour = "black", shape = 22, size = 4) +
      scale_x_continuous(expand=c(0,0)) + 
      ylim(c(0, max(sr$spEsc)*1.05)) 
  } 
  if (make_plot == TRUE) {
    print(p)
  }
} # end of plotStockSpEsc()


# age-structured model that holds early ocean survival constant [low, mid, high, very high]
# only time varying parameter is pre-spawning survival
# Vary: alpha or beta f
# Constant: EQsp, wanted_frac, water management scenario, climate scenario, and initial age-structure

# surv1 <- 0.05
# surv2 <- 0.8
# surv3 <- 0.8
# survPS <- "vary"
# wanted_frac <- 0.5 
# EQsp <- 7500; alpha_vec <- c(1.5, 3, 5, 10)
# # EQsp <- c(2500, 7500, 15000); alpha_vec <- 2
# # EQsp <- 7500; alpha_vec <- 2
# wtrMgmt <- "NoD" # "NoD" or "BAU" ""
# climate <- "A2" # "B1" or "A2" ""
# variableSR <- FALSE
# yrs <- 2010:2099

popProjConstMarine <- function (surv1, 
                                surv2 = 0.8, 
                                surv3 = 0.8,
                                survPS = 1,
                                wanted_frac,
                                alpha_vec,
                                EQsp,
                                wtrMgmt = "", # can be "" for test case of survPS = 1
                                climate = "", # can be "" for test case of survPS = 1
                                yrs = 2010:2099,
                                make_plot,
                                variableSR = FALSE) {
  sim_path <- file.path(".", "data", "simulation_results")
  # some checks on model parameterization
  if (length(surv1) != 1) stop("This function assumes constant early ocean survival, please enter a single value")
  
  if (is.numeric(survPS)) {
    print("survPS = 1 and should yield constant number of spawners over time (EQsp), except with beta is 0")
  } 
  
  if (survPS == "vary") {
    print(paste("Investigating", wtrMgmt, "and", climate, "prespawn mortality scenarios"))
  }
  
  if (any(alpha_vec <= 1) | any(alpha_vec > 10)) stop("1/spr multiplier must be greater than 1 or less than or equal to 10") 
  
  
  # start parameterizing model
  # determine the value for delta_e based on the fraction of age-3 spawner to total spawners of interest
  # The delta_e calcs are based on longterm/equilibrium assumptions.
  # CADFG data are based on recovery BC SRCS CWT recoveries. Might not be definitive, but provide 
  # a range to test
  delta_e <- calc_de(surv3, wanted_frac)  
  
  # Calculate SPR, use survPS = 1 (no oversummer prespawning mortality) to get maximum SPR for assumptions
  # about surv1 and delta_e (assume surv2 and surv3 fixed)
  spr <- SPR_srcs(surv1, surv2, surv3, delta_e = delta_e, survPS = 1)
  
  # alpha_vec contains multipliers that are different assumptions 
  # about the productivity of the population at
  # low abundances relatiove to the slope of the replacement
  # line: alpha_vec = 2 indicates --> alpha = 2*(1/SPR)
  if (exists("alpha")) rm(alpha)
  alpha <- 1/spr*alpha_vec
  
  # For a given SPR, alpha and EQsp, we can calculate beta, for Butte Creek 
  # spring run chinook salmon
  # To calculate beta for multiple alpha values, just pass a vecter of alphas
  # to the alpha argument
  
  beta <- calc_beta(alpha = alpha, EQsp = EQsp, spr = spr)
  
  # Plot underlying stock recruitment functions 
  # Spawners vs. Age-1 
  
  plotStockRecruit(alpha = alpha, beta = beta, EQsp = EQsp, delta_e = delta_e, spr = spr, make_plot = make_plot)
  #print("works to plotStockRecruit") 
  # Spawners vs. returning spawners (3 and 4 years later from resultant cohort)
  plotStockSpEsc(alpha = alpha, beta = beta, EQsp = EQsp, delta_e = delta_e, surv1 = surv1, surv2 = surv2, surv3 = surv3, spr = spr, make_plot = make_plot)
  #print("works to plotStockSpEsc") 
  # Select PreSpawn survival scenarios
  
  if (wtrMgmt == "NoD" & climate == "A2") {
    survPS <- as.matrix(read.delim(file.path(sim_path, "NoDiversion", "SALMODspawnFemSurvA2.txt"), sep = ","  )[,1:6])
  } 
  
  if (wtrMgmt == "NoD" & climate == "B1") {
    survPS <- as.matrix(read.delim(file.path(sim_path, "NoDiversion", "SALMODspawnFemSurvB1.txt"), sep = ","  )[,1:6])
  } 
  
  if (wtrMgmt == "BAU" & climate == "A2") {
    survPS <- as.matrix(read.delim(file.path(sim_path, "BAU", "SALMODspawnFemSurvA2.txt"), sep = ","  )[,1:6])     
  }
  
  if (wtrMgmt == "BAU" & climate == "B1") {
    survPS <- as.matrix(read.delim(file.path(sim_path, "BAU", "SALMODspawnFemSurvB1.txt"), sep = ","  )[,1:6])     
  }
  
  if (wtrMgmt == "" & climate == "") survPS <- 1
  
  gcm_names <- dimnames(survPS)[[2]]
  
  # Set up population projection model ####
  
  # Set up an initial age vector
  # start from equilibrium age structure in ocean just before returning to spawn
  
  n_alpha <- length(alpha)
  n_EQsp <- length(EQsp)
  
  if (n_EQsp > 1 & n_alpha > 1) stop("n_alpha and n_EQsp cannot both be greater than 1")
  
  if (n_EQsp > 1 & n_alpha == 1)  n0 <- matrix(0, nrow = 4, ncol = n_EQsp)
  if (n_EQsp == 1 & n_alpha >= 1)  n0 <- matrix(0, nrow = 4, ncol = 1)
  
  # When using multiple EQ spawning levels, need to adjust the intial age distribution
  # based on higher or lower equilibrium abundances
  # There will be an n0 for each EQsp
  # For a single EQsp, there is only a single 
  
  if (n_EQsp > 1) {
    for (i in 1:ncol(n0)) {
      n0[1,i] <- (EQsp[i]*(1-wanted_frac))
      n0[2,i] <- n0[1,i]/(surv3*(1-delta_e))
      n0[3,i] <- n0[2,i]/surv2
      n0[4,i] <- n0[3,i]/surv1
    } 
  } else {
    n0[1] <- (EQsp*(1-wanted_frac))
    n0[2] <- n0[1]/(surv3*(1-delta_e))
    n0[3] <- n0[2]/surv2
    n0[4] <- n0[3]/surv1
  }
  
  n0 <- apply(n0,2,rev) # This flips the matix, so that age-1 is on the first row rather than the 4th row 
  
  # Storage array for population simulations
  # 1-d for each climate model (6)
  # 1-d for each productivity hypothesis (1-5) 
  # 1-d for time (90)
  # 1-d for each age (4)
  # so need a 6 X 1-5 X 90 array
  
  # Create arrays for storing model output
  if (n_EQsp > 1 & n_alpha == 1) {
    pop <- array(0, c(4, length(yrs), 6, n_EQsp))
    spawners <-  array(0, c(length(yrs), 6, n_EQsp))
  }
  
  if (n_alpha >= 1 & n_EQsp == 1) {
    pop <- array(0, c(4, length(yrs), 6, n_alpha))
    spawners <-  array(0, c(length(yrs), 6, n_alpha))
  }
  
  if (n_EQsp == 1 & n_alpha == 1) {
    pop <- array(0, c(4, length(yrs), 6, 1))
    spawners <-  array(0, c(length(yrs), 6, 1))
  }
  
  if (n_EQsp == 1 & n_alpha == 1) {prodH0 <- "oneSR"; j_len <- 1}
  if (n_EQsp > 1 & n_alpha == 1)  {prodH0 <- "EQspVary"; j_len <- n_EQsp}
  if (n_EQsp == 1 & n_alpha > 1)  {prodH0 <- "alphaVary"; j_len <- n_alpha}
  
  for (j in 1:j_len) {  
    for (c in 1:6) {
      for (t in 1:length(yrs)) {
        
        if (t == 1) { 
          # pop[a, t, c, j] # a(ge), t(years), c(limate model), j(hypothesis about EQsp or alpha)
          if (n_EQsp == 1) {
            pop[1, t, c, j] <- n0[1]
            pop[2, t, c, j] <- n0[2]
            pop[3, t, c, j] <- n0[3]
            pop[4, t, c, j] <- n0[4]
          } 
          
          if (n_EQsp > 1) {
            pop[1, t, c, j] <- n0[1,j]
            pop[2, t, c, j] <- n0[2,j]
            pop[3, t, c, j] <- n0[3,j]
            pop[4, t, c, j] <- n0[4,j]
          }
          
        } else { 
          if (length(survPS) > 1 ) {
            spawners[t-1, c, j] <- pop[3, t-1, c, j]*delta_e*survPS[t-1, c] + pop[4,t-1, c, j]*survPS[t-1, c] 
          } else {
            spawners[t-1, c, j] <- pop[3, t-1, c, j]*delta_e + pop[4,t-1, c, j]
          }
          
          if (variableSR == FALSE) {
            if (n_alpha > 1 & n_EQsp == 1) pop[1,t, c, j] <- (spawners[t-1, c, j]*alpha[j])/(1+spawners[t-1, c, j]*beta[j])
            if (n_alpha == 1 & n_EQsp > 1) pop[1,t, c, j] <- (spawners[t-1, c, j]*alpha)/(1+spawners[t-1, c, j]*beta[j])
            if (n_alpha == 1 & n_EQsp == 1) pop[1,t, c, j] <- (spawners[t-1, c, j]*alpha)/(1+spawners[t-1, c, j]*beta)
          } else {
            if (n_alpha > 1 & n_EQsp == 1)  pop[1,t, c, j] <- (spawners[t-1, c, j]*alpha[j])/(1+spawners[t-1, c, j]*beta[j])*(exp(sig_r*rnorm(1,mean=0, sd=1)))
            if (n_alpha == 1 & n_EQsp > 1)  pop[1,t, c, j] <- (spawners[t-1, c, j]*alpha)/(1+spawners[t-1, c, j]*beta[j])*(exp(sig_r*rnorm(1,mean=0, sd=1)))
            if (n_alpha == 1 & n_EQsp == 1)  pop[1,t, c, j] <- (spawners[t-1, c, j]*alpha)/(1+spawners[t-1, c, j]*beta)*(exp(sig_r*rnorm(1,mean=0, sd=1)))
          }
          pop[2,t, c, j] <- pop[1,t-1, c, j]*surv1
          pop[3,t, c, j] <- pop[2, t-1, c, j]*surv2
          pop[4,t, c, j] <- pop[3, t-1, c, j]*surv3*(1-delta_e)
        }
        
        if (t == length(yrs) & length(survPS) > 1 )  spawners[t] <- pop[3, t, c, j]*delta_e*survPS[t, c] + pop[4, t, c, j]*survPS[t, c] 
        if (t == length(yrs) & length(survPS) == 1 )  spawners[t] <- pop[3,t, c, j]*delta_e + pop[4, t, c, j]
      }
    }
  }
  
  if (n_alpha > 1) {
    dimnames(pop)[4] <- list(alpha)
    dimnames(pop)[3] <- list(gcm_names)
    dimnames(spawners)[3] <- list(alpha)
    dimnames(spawners)[2] <- list(gcm_names)
    popm <- melt(pop, varnames = c("age", "year", "gcm", "alpha"))
    spm <- melt(spawners, varnames = c("year", "gcm", "alpha"))
    spm_sc <- ddply(spm, .(gcm, alpha), transform, value = value/(max(value, na.rm = TRUE)))
  }
  if (n_EQsp > 1) {
    dimnames(pop)[4] <- list(EQsp)
    dimnames(pop)[3] <- list(gcm_names)
    dimnames(spawners)[3] <- list(EQsp)
    dimnames(spawners)[2] <- list(gcm_names)
    popm <- melt(pop, varnames = c("age", "year", "gcm", "EQsp"))
    
    #popm_sc <- ddply(spm, .(gcm, EQsp), transform, scale_abund = value/(max(value, na.rm = TRUE)))
    spm <- melt(spawners, varnames = c("year", "gcm", "EQsp"))
    spm_sc <- ddply(spm, .(gcm, EQsp), transform, value = value/(max(value, na.rm = TRUE)))
  }
  
  survPS <- data.frame(year = 1:90, survPS)
  survPSm <- melt(survPS, id.var = "year")
  names(survPSm) <- c("year", "gcm", "value")
  
  if (n_alpha > 1)   survPSm$alpha = 1 # 'ref'
  if (n_EQsp > 1)   survPSm$EQsp = 1 # 'ref'
  
  survPS <- survPSm[,c(1,2,4,3)]
  spm_sc <- rbind(spm_sc, survPSm)
  return(list(popm = popm, spm = spm, spm_sc = spm_sc, 
              delta_e = delta_e, spr = spr, alpha=alpha, 
              beta=beta))
  
} # end of popProjConstMarine() 


plotPopProjConst <- function(population, spawners, spawners_sc, vary) {
  if (vary == "alpha") {
    p_pop <- ggplot(population, aes(x = year, y = value, colour = factor(gcm))) + 
      geom_line() + 
      #geom_point() + 
      facet_grid(age ~ alpha, scales = "free_y") +
      ylab("Abundance") + 
      xlab("Year") +
      theme_few() +
      scale_colour_few()
    
    
    print(p_pop)
    
    p_sp <- ggplot(spawners, aes(x = year, y = value, colour = factor(gcm))) + 
      geom_line() + 
      facet_grid(alpha ~ ., scales = "free_y") + 
      ylab("Spawners") + 
      xlab("Year") +
      theme_few() +
      scale_colour_few()
    
    print(p_sp)
    
    p_spm_sc <- ggplot(spawners_sc, aes(x = year, y = value, colour = factor(gcm))) + 
      geom_line() + 
      facet_grid(alpha ~ ., scales = "free_y") + 
      ylab("Scaled Spawners") + 
      xlab("Year") +
      theme_few() +
      scale_colour_few()
    
    print(p_spm_sc)
  }
  
  if (vary == "EQsp") {
    p_pop <- ggplot(population, aes(x = year, y = value, colour = factor(gcm))) + 
      geom_line() + 
      #geom_point() + 
      facet_grid(age ~ EQsp, scales = "free_y") +
      ylab("Abundance") + 
      xlab("Year") +
      theme_few() +
      scale_colour_few()
    
    print(p_pop)
    
    p_sp <- ggplot(spawners, aes(x = year, y = value, colour = factor(gcm))) + 
      geom_line() + 
      facet_grid(EQsp ~ ., scales = "free_y") + 
      ylab("Spawners") + 
      xlab("Year") +
      theme_few() +
      scale_colour_few()
    
    print(p_sp)
    
    p_spm_sc <- ggplot(spawners_sc, aes(x = year, y = value, colour = factor(gcm))) + 
      geom_line() + 
      facet_grid(EQsp ~ ., scales = "free_y") + 
      ylab("Scaled Spawners") + 
      xlab("Year") +
      theme_few() +
      scale_colour_few()
    
    print(p_spm_sc)
  }
  
  if (vary == "oneEach") {
    p_pop <- ggplot(popm, aes(x = year, y = value, colour = factor(gcm))) + 
      geom_line() + 
      #geom_point() + 
      ylab("Abundance") + 
      xlab("Year") +
      theme_few() +
      scale_colour_few()
    
    print(p_pop)
    
    p_sp <- ggplot(spm, aes(x = year, y = value, colour = factor(gcm))) + 
      geom_line() + 
      ylab("Spawners") + 
      xlab("Year")
    print(p_sp) +
      theme_few() +
      scale_colour_few()
  }
  
} # end plotPopProjConst()


# this function counts the first four consecutive years of zero (survival or abundance)
# and then records how many years to the first or last of those consecutive year and the year
# climate <- "B1"
# mgmt <- "NoDiversion"
# 
# data <- case1A2[[3]]
calcExtTimeConst <- function(data, climate, mgmt) {
  
  library(ggthemes)
  if (names(data)[3] == "EQsp") {
    sr <- "EQsp" }  else if (names(data[3]) == "alpha") {
      sr <- "alpha" } else {
        print("names data[3] must be EQsp or alpha")
      }
  names(data)[3] <- "SR"
  datac <- dcast(data, year ~ SR + gcm) 
  extInd <- vector("integer", (ncol(datac)-1) )
  for (j in 2:ncol(datac)) {
    zeroInd <- which(datac[,j] == 0) 
    #zeroInd <- c(1,0,0,1,0,2,3,4,3,4,1)
    extRuns <- vector("integer", 0)
    k <- 1
    for (i in 1:(length(zeroInd)-3)) {
      if (length(zeroInd) >= 4) {
        if (zeroInd[i] == (zeroInd[i+1] - 1) & zeroInd[i] == (zeroInd[i+2] - 2) & 
              zeroInd[i] == (zeroInd[i+3] - 3)  ) {
          extRuns[k] <- (zeroInd[i] + 3)
        } else {
          extRuns[k] <- 90
        }
      } else {
        extRuns[k] <- 90
      }
      k <- k + 1 
    }
    if (length(extRuns > 0)) extInd[j-1] <- min(extRuns) else extInd[j-1] <- NA
  }
  
  yrs <- 2010:2099
  extYr <- vector("integer", length(extInd))
  for (i in 1:length(extInd)) {
    extYr[i] <- yrs[extInd[i]]
  }
  
  if (sr == "EQsp") {
    SRinfo <- trunc(as.numeric(gsub("([0-9]{1,5})(_)([a-z0-9]+)", "\\1" , names(datac)[-1])))
    gcm <- gsub("([0-9]{1,5})(_)([a-z0-9]+)", "\\3" , names(datac)[-1])
  } else if (sr == "alpha") {
    SRinfo <- trunc(as.numeric(gsub("([0-9\\.]{1,7})(_)([a-z0-9]+)", "\\1" , names(datac)[-1])))
    gcm <- gsub("([0-9\\.]{1,7})(_)([a-z0-9]+)", "\\3" , names(datac)[-1])
  }
  outData <- data.frame(SRinfo = SRinfo, gcm = gcm, yrs2ext = extInd, extYr = extYr)
  
  #   time2extGCM <- ggplot(data = outData, aes(x = yrs2ext, y = gcm, colour = gcm)) + 
  #     # geom_point(size = 4, alpha = .7, colour = "black") + # position = position_jitter(w = 0, h = 0.5), 
  #     geom_point(size = 3, alpha = .7, position = position_jitter(w = 0, h = 0.5)) + # 
  #     scale_color_brewer(palette = "Dark2") +
  #     theme_bw() +
  #     xlim(c(0,90)) +
  #     xlab("Years to extinction") +
  #     ggtitle(title = paste("Mean time to extinction for", climate, "&", mgmt, "scenario\nvarying", sr, sep = " "))
  #   print(time2extGCM)
  
  xTicks <- unique(outData$SRinfo)
  
  time2extSR <- ggplot(data = outData, aes(x = yrs2ext, y = SRinfo, colour = gcm)) + 
    geom_point(size = 4, colour = "black") +
    geom_point(size = 3) + 
    scale_color_brewer(palette = "Dark2") +
    theme_bw() +
    xlim(c(0,90)) +
    ylab(paste(sr)) +
    xlab("Years to extinction") +
    ggtitle(paste("Mean time to extinction for", climate, "&", mgmt, "scenario\nvarying", sr, sep = " ")) +
    scale_y_continuous(breaks= xTicks)  
  print(time2extSR)
  
  return(outData)
}


# Lisa says that extinction is a QE of 4 years sans 20 fish

# # Vary alpha, A2 BAU
# data = case1A2#  
# mgmt = "BAU"
# climate = "A2"
# nSpawners = 7500
# working = YES
# getSALMODnums(case1A2) # Works



# Vary beta, A2 BAU
# data <- case2A2
# mgmt = "BAU"
# climate = "A2"
# calcQEtimeConst(case2A2, "A2", "BAU", QEthr = 100)
# working = YES (but should check)
# getSALMODnums(case2A2) #  WORKS

# # Vary alpha, B1 No Diversion
# data = case1B1
# mgmt = "NoDiversion"
# climate = "B1"
# calcQEtimeConst(case1B1, "B1", "No Diversion", QEthr = 100)
# working = NO 
# getSALMODnums(case1B1) # WORKS


# # Vary alpha, B1 No Diversion
# data = case2B1
# mgmt = "NoDiversion"
# climate = "B1"
# calcQEtimeConst(case2B1, "B1", "No Diversion", QEthr = 100)
# working = Yes
# getSALMODnums(case2B1) # WORKS

# Pick a QE threshold
# QEthr <- 2 # 100, 20, 0

getSALMODnums <- function(Sdata) {
  # Function to get the number of females surviving to spawn in the fall calculated by SALMOD
  # to compare SALMOD time to QE with those from the deterministic population projection model
  # Just get the SALMOD prespawn survival output (stored in the [[3]] list from deterministic population projection model)
  names(Sdata[[3]])[3] <- "SR"
  SALMOD_surv  <-  dcast(subset(Sdata[[3]], SR == 1)[,c(1,2,4)], year ~ gcm)
  SALMOD_num <- cbind(SALMOD_surv[,1], SALMOD_surv[,2:ncol(SALMOD_surv)] * 7568) # 7568 is number of spawners used by LCT and CM
  return(SALMOD_num)
} # end getSALMODnums()


customFR2ts <- function(N, # number of time steps
                        reps,
                        r_seed = 1,
                        amp) {
  # LWB notes on generating time series with a specific 
  # So, if you want to create white noise (equal variance at all frequencies) 
  # you should generate a sine wave of frequency .1, .2, .3,....up to the maximum frequency you want, 
  # then give each one a phase picked from a uniform distribution between 0 and 2 pi, 
  # then add them together and divide by whatever it takes to get the variance you want.
  # 
  # If you instead want to create a series to which salmon would be sensitive do not make 
  # all of the sine waves the same amplitude, rather make the sine waves near 1/(generation time) 
  # larger than the others somehow.  One way to do that would be to make them have the a Gaussian shape.
  # 
  # Make sense?
  
  t <- 1:N # time index
  f <- (1:(N/2))/N # frequencies from 1/N to 0.5
  tsReps <- matrix(NA, nrow = N, ncol = reps)
  set.seed(r_seed)
  for (h in 1:reps) {
    
    rtheta <- runif(length(f), 0, 2*pi) # the random phases to apply to each frequency
    
    ts <- matrix(NA, nrow = N, ncol = length(f) ) 
    
    for (i in 1:length(f)) {
      if (length(amp) == 1) {
        ts[,i] <- amp * cos(2*pi*f[i]*t - rtheta[i])
      } else {
        ts[,i] <- amp[i] * cos(2*pi*f[i]*t - rtheta[i])         
      }
    }
    
    noise <- rowSums(ts) # add up the curves
    noise <- (noise - mean(noise, na.rm = TRUE) )/ sd(noise, na.rm = TRUE) # mean = 0, sd of 1
    tsReps[,h] <- noise
  }
  return(tsReps)
  
} # end of customFR2ts

# a wrapper for customFR2ts ####
noiseWrapper <- function(N, reps, r_seed = 1) {
  # create "reps" number of white noise time series of length N using random sine wave approach
  # with selected frequency contents
  white_n <- customFR2ts(N = N, # number of time steps
                         reps = reps,
                         r_seed = r_seed,
                         amp = mk_white(N)) # mean = 0, variance/sd = 1
  
  # Bandpass period 3-4 (frequencies 1/4 to 1/3) in white noise signals.
  
  rsin_34_n <- customFR2ts(N = N, # number of time steps
                           reps = reps,
                           r_seed = r_seed,
                           amp = mk_rsin(N, highF=1/3, lowF=1/4)) # mean = 0, variance/sd = 1
  
  # Bandreject period 3-4 (frequencies 1/4 to 1/3) in white noise signals.
  
  rsin_34_reject_n <- customFR2ts(N = N,
                                  reps = reps,
                                  r_seed = r_seed,
                                  amp = mk_rsin_reject(N, lowF=1/4, highF=1/3) )
  # Band-pass greater than period 4 (frequencies lower than 0.25)
  rsin_gt4_n <- customFR2ts(N = N,
                            reps = reps,
                            r_seed = r_seed,
                            amp = mk_rsin(N, lowF=0, highF=1/4) )
  
  # Band-pass less than period 3 (frequencies higher than 0.33)
  rsin_lt3_n <- customFR2ts(N = N,
                            reps = reps,
                            r_seed = r_seed,
                            amp = mk_rsin(N, lowF=1/3, highF=1/2) )
  
  # Band-pass greater than period 3 (frequencies lower than 0.33)
  rsin_gt3_n <- customFR2ts(N = N,
                            reps = reps,
                            r_seed = r_seed,
                            amp = mk_rsin(N, lowF=0, highF=1/3) )
  
  # Band-pass less than period 4 (frequencies higher than 0.25)
  rsin_lt4_n <- customFR2ts(N = N,
                            reps = reps,
                            r_seed = r_seed,
                            amp = mk_rsin(N, lowF=1/4, highF=1/2) )
  
  # Band-pass greater than period 10 (frequencies lower than 0.1)
  rsin_gt10_n <- customFR2ts(N = N,
                             reps = reps,
                             r_seed = r_seed,
                             amp = mk_rsin(N, lowF=0, highF=1/10) )  
  
  # Band-pass greater than period 10 and period 3-4 (frequencies lower than 0.1 plus period 3-4)
  # By adding the period 10 and greater noise to the period 3-4 only noise (divide by 2) 
  rsin_34_gt10_n <- matrix(NA, nrow = N, ncol = reps)
  for (i in 1:reps) {
    rsin_34_gt10_n[,i] <- (rsin_gt10_n[,i] + rsin_34_n[,i])/2
  }
  rsin_34_gt10_n <- apply(rsin_34_gt10_n, 2, scale)
  
  # Store simulated noise sets in a list for use in simualtions (easier indexing)
  noiseList <- list(noise_white = white_n,
                    noise_34 = rsin_34_n,
                    nose_34_reject = rsin_34_reject_n,
                    noise_gt4 = rsin_gt4_n,
                    noise_lt3 = rsin_lt3_n,
                    noise_gt3 = rsin_gt3_n,
                    noise_lt4 = rsin_lt4_n,
                    noise_gt10 = rsin_gt10_n,
                    noise_34gt10 = rsin_34_gt10_n)
  
  # rm individual sets of noise
  rm(white_n, rsin_34_n, rsin_34_reject_n, rsin_gt4_n, rsin_lt3_n, rsin_gt3_n, rsin_lt4_n, rsin_gt10_n, rsin_34_gt10_n)
  return(noiseList)
} # end of noiseWrapper <- function(N, reps, r_seed = 1)

# plot the mean frequency response of multiple spectra generated by the same noise ####

plotMeanFreqR <- function(dataMat, N) {
  # spectral frequency 
  if (trunc(sqrt(N)) %% 2 == 0) {
    m <- trunc(sqrt(N)) + 1 
  } else {
    m <- trunc(sqrt(N))
  }
  
  spcMean <- matrix(NA, nrow = N/2, ncol = ncol(dataMat))
  
  freq <- (1:(N/2))/N
  for (i in 1:ncol(dataMat)) {
    ifelse(all(dataMat[,i] == 0), 
            spcMean[,i] <- rep(NA, times = N/2),
           spcMean[,i] <- spec.pgram(scale(dataMat[,i]), plot = F, c(m,m))$spec )
   }
  mean_spc <- rowMeans(spcMean, na.rm = TRUE)
  
  spc10 <- apply(spcMean, 1, quantile, probs = c(0.1), na.rm = TRUE)
  spc25 <- apply(spcMean, 1, quantile, probs = c(0.25), na.rm = TRUE)
  spc75 <- apply(spcMean, 1, quantile, probs = c(0.75), na.rm = TRUE)
  spc90 <- apply(spcMean, 1, quantile, probs = c(0.9), na.rm = TRUE)
  
  plot(freq, mean_spc, type = "n", lwd = 3.7, lty = 1, col = "black", xlab = "Frequency", ylab = "Relative Magnitude", las = 1, 
       ylim = c(0, ceiling(max(spc90, na.rm = TRUE))))
  
  lines(freq, spc90, lwd = 2)
  lines(freq, spc10, lwd = 2)
  polygon(x = c(freq, rev(freq)), y = c(spc10, rev(spc90)), col = "lightgrey")
  
  lines(freq, spc75, lwd = 2)
  lines(freq, spc25, lwd = 2)
  polygon(x = c(freq, rev(freq)), y = c(spc25, rev(spc75)), col = "darkgrey")
  
  lines(freq, mean_spc, type = "l", lwd = 3.7, lty = 1, col = "black")
  lines(freq, mean_spc, type = "l", lwd = 2, col = "white")
  lines(freq, mean_spc, lwd = 2, lty = 2)
  
  box(lwd=2)
}

# plot mean Frequency response for data.tables

plotMeanFreqR_DT <- function(dataTable, N, surv, scale = "CV", yaxis_lim = c(0, ceiling(max(spc90, na.rm = TRUE))) ) {
  # spectral frequency 
  if (trunc(sqrt(N)) %% 2 == 0) {
    m <- trunc(sqrt(N)) + 1 
  } else {
    m <- trunc(sqrt(N))
  }
  
  spcMean <- matrix(NA, nrow = N/2, ncol = length(unique(dataTable$reps_c)))
  
  freq <- (1:(N/2))/N
  for (i in 1:ncol(spcMean)) {
    ifelse(all(dataTable[i = reps_c == i & meanPS_c >= surv-0.001 & meanPS_c <= surv+0.001,
                         j = white] == 0), 
           spcMean[,i] <- rep(NA, times = N/2),
           ifelse(scale == "CV", 
                  spcMean[,i] <- periodogram(dataTable[i = reps_c == i & meanPS_c >= surv-0.001 & meanPS_c <= surv+0.001,
                                                             j = white]/mean(dataTable[i = reps_c == i & meanPS_c >= surv-0.001 & meanPS_c <= surv+0.001,
                                                                                       j = white]), plot = F )$spec, 
                  spcMean[,i] <- periodogram(scale(dataTable[i = reps_c == i & meanPS_c >= surv-0.001 & meanPS_c <= surv+0.001,
                                                             j = white]), plot = F )$spec ) )
#            spcMean[,i] <- spec.pgram(scale(dataTable[i = reps_c == i & meanPS_c == surv,
#                                                      j = white]), plot = F, c(m,m))$spec )
  }
  mean_spc <- rowMeans(spcMean, na.rm = TRUE)
  
  spc10 <- apply(spcMean, 1, quantile, probs = c(0.1), na.rm = TRUE)
  spc25 <- apply(spcMean, 1, quantile, probs = c(0.25), na.rm = TRUE)
  spc75 <- apply(spcMean, 1, quantile, probs = c(0.75), na.rm = TRUE)
  spc90 <- apply(spcMean, 1, quantile, probs = c(0.9), na.rm = TRUE)
  
  #if (!exists("yaxis_lim")) yaxis_lim <- c(0, ceiling(max(spc90, na.rm = TRUE)))
  
  plot(freq, mean_spc, type = "n", lwd = 3.7, lty = 1, col = "black", xlab = "Frequency", ylab = "Relative Magnitude", las = 1, 
       ylim = yaxis_lim)
  
  lines(freq, spc90, lwd = 2)
  lines(freq, spc10, lwd = 2)
  polygon(x = c(freq, rev(freq)), y = c(spc10, rev(spc90)), col = "lightgrey")
  
  lines(freq, spc75, lwd = 2)
  lines(freq, spc25, lwd = 2)
  polygon(x = c(freq, rev(freq)), y = c(spc25, rev(spc75)), col = "darkgrey")
  
  lines(freq, mean_spc, type = "l", lwd = 3.7, lty = 1, col = "black")
  lines(freq, mean_spc, type = "l", lwd = 2, col = "white")
  lines(freq, mean_spc, lwd = 2, lty = 2)
  
  box(lwd=2)
}

# Plot only the mean Frequency Response, and allow adding addition Frequency Response plots (lines)


plotMeanFR_DTmanyAR <- function(dataTable, N, surv, scale = "CV", AR_col, yaxis_lim) {
  # spectral frequency 
  if (trunc(sqrt(N)) %% 2 == 0) {
    m <- trunc(sqrt(N)) + 1 
  } else {
    m <- trunc(sqrt(N))
  }
  
  spcMean <- matrix(NA, nrow = N/2, ncol = length(unique(dataTable$reps_c)))
  
  freq <- (1:(N/2))/N
  for (i in 1:ncol(spcMean)) {

    ts <- as.ts(droplevels(subset(dataTable, reps_c == i & meanPS_c >= surv-0.001 & meanPS_c <= surv+0.001, select = AR_col) ))
  
    ifelse(all(ts == 0), 
           spcMean[,i] <- rep(NA, times = N/2),
           ifelse(scale == "CV", 
                  spcMean[,i] <- periodogram(ts/mean(ts), plot = F )$spec, 
                  spcMean[,i] <- periodogram(scale(ts), plot = F )$spec ) )
  }
  mean_spc <- rowMeans(spcMean, na.rm = TRUE)
  
  spc10 <- apply(spcMean, 1, quantile, probs = c(0.1), na.rm = TRUE)
  spc25 <- apply(spcMean, 1, quantile, probs = c(0.25), na.rm = TRUE)
  spc75 <- apply(spcMean, 1, quantile, probs = c(0.75), na.rm = TRUE)
  spc90 <- apply(spcMean, 1, quantile, probs = c(0.9), na.rm = TRUE)

  if ( !exists("yaxis_lim") ) ifelse(scale == "CV", yaxis_lim <- c(0,1.2*max(mean_spc[which(freq > 0.2)])), yaxis_lim <- c(0, 20))
  
  plot(freq, mean_spc, type = "n", xlab = "Frequency", ylab = "Relative Magnitude", las = 1, ylim = yaxis_lim)
  lines(freq, mean_spc, type = "l", lwd = 3.7, lty = 1, col = "black")
  lines(freq, mean_spc, type = "l", lwd = 2, col = "white")
  lines(freq, mean_spc, lwd = 2, lty = 2)

  box(lwd=2)
}

plotMeanFR_DTmany <- function(dataTable, N, surv, scale = "CV", yaxis_lim) {
  # spectral frequency 
  if (trunc(sqrt(N)) %% 2 == 0) {
    m <- trunc(sqrt(N)) + 1 
  } else {
    m <- trunc(sqrt(N))
  }
  
  spcMean <- matrix(NA, nrow = N/2, ncol = length(unique(dataTable$reps_c)))
  
  freq <- (1:(N/2))/N
  for (i in 1:ncol(spcMean)) {
    ifelse(all(dataTable[i = reps_c == i & meanPS_c >= surv-0.001 & meanPS_c <= surv+0.001,
                         j = white] == 0), 
           spcMean[,i] <- rep(NA, times = N/2),
           ifelse(scale == "CV", 
                  spcMean[,i] <- periodogram(dataTable[i = reps_c == i & meanPS_c >= surv-0.001 & meanPS_c <= surv+0.001,
                                                       j = white]/mean(dataTable[i = reps_c == i & meanPS_c >= surv-0.001 & meanPS_c <= surv+0.001,
                                                                                 j = white]), plot = F )$spec, 
                  spcMean[,i] <- periodogram(scale(dataTable[i = reps_c == i & meanPS_c >= surv-0.001 & meanPS_c <= surv+0.001,
                                                             j = white]), plot = F )$spec ) )
    
  }
  mean_spc <- rowMeans(spcMean, na.rm = TRUE)

  spc10 <- apply(spcMean, 1, quantile, probs = c(0.1), na.rm = TRUE)
  spc25 <- apply(spcMean, 1, quantile, probs = c(0.25), na.rm = TRUE)
  spc75 <- apply(spcMean, 1, quantile, probs = c(0.75), na.rm = TRUE)
  spc90 <- apply(spcMean, 1, quantile, probs = c(0.9), na.rm = TRUE)
 
  if ( !exists("yaxis_lim") ) ifelse(scale == "CV", yaxis_lim <- c(0,1.2*max(mean_spc[which(freq > 0.2)])), yaxis_lim <- c(0, 20))

  plot(freq, mean_spc, type = "n", xlab = "Frequency", ylab = "Relative Magnitude", las = 1, ylim = yaxis_lim)
  lines(freq, mean_spc, type = "l", lwd = 3.7, lty = 1, col = "black")
  lines(freq, mean_spc, type = "l", lwd = 2, col = "white")
  lines(freq, mean_spc, lwd = 2, lty = 2)
  
  box(lwd=2)
}

linesMeanFR_DTmanyAR <- function(dataTable, N, surv, scale = "CV", line_color, AR_col) {
  # spectral frequency 
  if (trunc(sqrt(N)) %% 2 == 0) {
    m <- trunc(sqrt(N)) + 1 
  } else {
    m <- trunc(sqrt(N))
  }
  
  spcMean <- matrix(NA, nrow = N/2, ncol = length(unique(dataTable$reps_c)))
  
  freq <- (1:(N/2))/N
  for (i in 1:ncol(spcMean)) {
    
    ts <- as.ts(droplevels(subset(dataTable, reps_c == i & meanPS_c >= surv-0.001 & meanPS_c <= surv+0.001, select = AR_col) ))
        
    ifelse(all(ts == 0), 
           spcMean[,i] <- rep(NA, times = N/2),
           ifelse(scale == "CV", 
                  spcMean[,i] <- periodogram(ts/mean(ts), plot = F )$spec, 
                  spcMean[,i] <- periodogram(scale(ts), plot = F )$spec ) )
  }
  mean_spc <- rowMeans(spcMean, na.rm = TRUE)
  
  spc10 <- apply(spcMean, 1, quantile, probs = c(0.1), na.rm = TRUE)
  spc25 <- apply(spcMean, 1, quantile, probs = c(0.25), na.rm = TRUE)
  spc75 <- apply(spcMean, 1, quantile, probs = c(0.75), na.rm = TRUE)
  spc90 <- apply(spcMean, 1, quantile, probs = c(0.9), na.rm = TRUE)
  
  lines(freq, mean_spc, type = "l", lwd = 3.7, lty = 1, col = line_color)
  
  box(lwd=2)
}


linesMeanFR_DTmany <- function(dataTable, N, surv, scale = "CV", line_color) {
  # spectral frequency 
  if (trunc(sqrt(N)) %% 2 == 0) {
    m <- trunc(sqrt(N)) + 1 
  } else {
    m <- trunc(sqrt(N))
  }
  
  spcMean <- matrix(NA, nrow = N/2, ncol = length(unique(dataTable$reps_c)))
  
  freq <- (1:(N/2))/N
  for (i in 1:ncol(spcMean)) {
    ifelse(all(dataTable[i = reps_c == i & meanPS_c >= surv-0.001 & meanPS_c <= surv+0.001,
                         j = white] == 0), 
           spcMean[,i] <- rep(NA, times = N/2),
           ifelse(scale == "CV", 
                  spcMean[,i] <- periodogram(dataTable[i = reps_c == i & meanPS_c >= surv-0.001 & meanPS_c <= surv+0.001,
                                                             j = white]/mean(dataTable[i = reps_c == i & meanPS_c >= surv-0.001 & meanPS_c <= surv+0.001,
                                                                                       j = white]), plot = F )$spec, 
                  spcMean[,i] <- periodogram(scale(dataTable[i = reps_c == i & meanPS_c >= surv-0.001 & meanPS_c <= surv+0.001,
                                                             j = white]), plot = F )$spec ) )
    
    #            spcMean[,i] <- spec.pgram(scale(dataTable[i = reps_c == i & meanPS_c >= surv-0.001 & meanPS_c <= surv+0.001,
    #                                                      j = white]), plot = F, c(m,m))$spec )
  }
  mean_spc <- rowMeans(spcMean, na.rm = TRUE)
  
  spc10 <- apply(spcMean, 1, quantile, probs = c(0.1), na.rm = TRUE)
  spc25 <- apply(spcMean, 1, quantile, probs = c(0.25), na.rm = TRUE)
  spc75 <- apply(spcMean, 1, quantile, probs = c(0.75), na.rm = TRUE)
  spc90 <- apply(spcMean, 1, quantile, probs = c(0.9), na.rm = TRUE)
  
  lines(freq, mean_spc, type = "l", lwd = 3.7, lty = 1, col = line_color)
  
  box(lwd=2)
}


calcQEtimeConstPers <- function(data, climate, mgmt, QEthr) {
  
  library(ggthemes)
  
  sr <- "alpha" 
  names(data[[2]])[3] <- "SR"
  datac <- dcast(data[[2]], year ~ SR + gcm) # numbers from deterministic population model 
  salmod_nums <- getSALMODnums(data) # SALMOD Numbers
  
  # Calculate time to Quasi-extinction for the deterministic population model
  extInd <- vector("integer", (ncol(datac)-1) )
  for (j in 2:ncol(datac)) {
    qeInd <- which(datac[,j] <= QEthr) 
    extRuns <- vector("integer", 0)
    k <- 1
    if (length(qeInd) >= 4) {
      for (i in 1:(length(qeInd)-3)) {
        if (qeInd[i] == (qeInd[i+1] - 1) & qeInd[i] == (qeInd[i+2] - 2) & 
              qeInd[i] == (qeInd[i+3] - 3)  ) {
          extRuns[k] <- (qeInd[i] + 3)
          k <- k + 1
        }  
      }
    }else {
      extRuns[k] <- 90
      k <- k + 1
    }
    if (length(extRuns > 0)) extInd[j-1] <- min(extRuns) else extInd[j-1] <- 90
  }
  # Calculate time to Quasi-extinction for the SALMOD data
  SALMODextInd <- vector("integer", (ncol(salmod_nums)-1) )
  for (j in 2:ncol(salmod_nums)) {
    qeInd <- which(salmod_nums[,j] <= QEthr) 
    SALMODextRuns <- vector("integer", 0)
    k <- 1
    if (length(qeInd) >= 4) {
      for (i in 1:(length(qeInd)-3)) {
        if (qeInd[i] == (qeInd[i+1] - 1) & qeInd[i] == (qeInd[i+2] - 2) & 
              qeInd[i] == (qeInd[i+3] - 3)  ) {
          SALMODextRuns[k] <- (qeInd[i] + 3)
          k <- k + 1
        }  
      } 
    } else {
      SALMODextRuns[k] <- 90
      k <- k + 1
    }
    if (length(SALMODextRuns > 0)) SALMODextInd[j-1] <- min(SALMODextRuns) else SALMODextInd[j-1] <- 90
  }
  
  yrs <- 2010:2099
  # Vector of extinction year for each climate model/stock recruit (pop dy) scenario
  extYr <- vector("integer", length(extInd))
  for (i in 1:length(extInd)) {
    extYr[i] <- yrs[extInd[i]]
  }
  # Vector of extinction year for each climate model/stock recruit (pop dy) scenario
  SALMODextYr <- vector("integer", length(SALMODextInd))
  for (i in 1:length(SALMODextInd)) {
    SALMODextYr[i] <- yrs[SALMODextInd[i]]
  }
  
  SRinfo <- trunc(as.numeric(gsub("([0-9\\.]{1,7})(_)([a-z0-9]+)", "\\1" , names(datac)[-1])))
  gcm <- gsub("([0-9\\.]{1,9})(_)([a-z0-9]+)", "\\3" , names(datac)[-1])

  outData <- data.frame(SRinfo = SRinfo*data[[5]], gcm = gcm, yrs2ext = extInd, extYr = extYr)
  SALMOD <- data.frame(gcm = unique(gcm), yrs2ext = SALMODextInd, extYr = SALMODextYr)
  xTicks <- 1:10#seq(0,1, by = 0.1) # unique(outData$SRinfo) 

  time2extSR <- ggplot(data = outData, aes(y = yrs2ext, x = SRinfo, colour = gcm)) + 
    geom_line(linetype = 1, size = 1) +
    geom_hline(aes(yintercept=90), size = 2) +
    scale_color_brewer(palette = "Dark2", name = "") +
    theme_classic() +
    xlab(expression(paste(alpha,  " x SPR"))) +
    ylab("Years to QE") +
    labs(title = paste(climate, " & ", mgmt, ", QE = ", QEthr, sep = "")) +
    scale_x_continuous(breaks = xTicks) +
    ylim(c(0,90)) + 
    geom_hline(data = SALMOD, aes(yintercept = yrs2ext, colour = gcm), linetype = 2,
               size = 1, alpha = 0.8) 
         
  print(time2extSR)
  
  
  extYrSR <- ggplot(data = outData, aes(y = extYr, x = SRinfo, colour = gcm)) + 
    geom_line(linetype = 1, size = 1) +
    scale_color_brewer(palette = "Dark2", name = "") +
    theme_classic() +
    xlab(expression(paste(alpha,  " x SPR"))) +
    ylab("Year of QE") +
    labs(title = paste(climate, " & ", mgmt, ", QE = ", QEthr, sep = "")) +
    scale_x_continuous(breaks= xTicks) +
    ylim(c(2010,2100)) +
    geom_hline(data = SALMOD, aes(yintercept = extYr, colour = gcm), linetype = 2,
               size = 1, alpha = 0.8) 
  print(extYrSR)
  
  return(outData) 
  
} # end calcQEtimeConst 

calcQEtimeConstEQvQE <- function(data, climate, mgmt) {
    
  # data = spawner output for many runs for each climate and mgmt e.g.,  EQsp <- seq(500, 50000, by = 500)
  # for each climate and mgmt loop over each EQsp and QEthr 
    
  QEthr <- seq(0, 1000, by = 10)
  sr <- "EQsp" 
  names(data[[2]])[3] <- "SR"
  datac <- dcast(data[[2]], year ~ SR + gcm) # numbers from deterministic population model 
    
  # Calculate time to Quasi-extinction for the deterministic population model
  extInd <- extYr <- matrix(NA, nrow = length(QEthr), ncol = (ncol(datac)-1))
  
  for (q in 1:length(QEthr)) {
    for (j in 2:ncol(datac)) {
      qeInd <- which(datac[,j] <= QEthr[q]) 
      extRuns <- vector("integer", 0)
      k <- 1
      if (length(qeInd) >= 4) {
        for (i in 1:(length(qeInd)-3)) {
          if (qeInd[i] == (qeInd[i+1] - 1) & qeInd[i] == (qeInd[i+2] - 2) & 
                qeInd[i] == (qeInd[i+3] - 3)  ) {
            extRuns[k] <- (qeInd[i] + 3)
            k <- k + 1
          }  
        }
      }else {
        extRuns[k] <- 90
        k <- k + 1
      }
      if (length(extRuns > 0)) extInd[q, j-1] <- min(extRuns) else extInd[q, j-1] <- 90
    }
  } # end of q loop
  
  yrs <- 2010:2099
  # Matrix of extinction year for each climate model/stock recruit (pop dy) scenario
  for (i in 1:nrow(extInd)) {
    for (j in 1:ncol(extInd)) {
      extYr[i,j] <- yrs[extInd[i,j]]
    }
  }
 
  # DATA should look like this:
  #   GCM  |  SRinfo  | QEthr | extInd  | extYr
  #   a    |   1000   | 10    |   50    |  ...
  #   a    |   1000   | 20    |   50    |  ...
  #   a    |   1000   | 50    |   50    |  ...
  #   a    |   2000   | 10    |   50    |  ...
  #   a    |   2000   | 20    |   50    |  ...
  #   a    |   2000   | 50    |   50    |  ...
  #   b    |   1000   | 10    |   50    |  ...
  #   b    |   1000   | 20    |   50    |  ...
  #   b    |   1000   | 50    |   50    |  ...
  #   b    |   2000   | 10    |   50    |  ...
  #   b    |   2000   | 20    |   50    |  ...
  #   b    |   2000   | 50    |   50    |  ...
  # and have GCM * SRinfo * QEthr rows
  
  # wrangle extInd data into above format
  extIndDF <- as.data.frame(extInd)
  extYrDF <- as.data.frame(extYr)
  names(extIndDF) <- names(extYrDF) <- names(datac)[-1]
  extIndDF <- cbind(QEthr = QEthr, extIndDF)
  extYrDF <- cbind(QEthr = QEthr, extYrDF)
  extIndDF_m <- melt(extIndDF, id.vars = "QEthr")
  extYrDF_m <- melt(extYrDF, id.vars = "QEthr")
  
  SRinfo <- trunc(as.numeric(gsub("([0-9]{1,16})(_)([a-z0-9]+)", "\\1" ,  extIndDF_m$variable)))#names(datac)[-1])))
  gcm <- gsub("([0-9.]{1,16})(_)([a-z0-9]+)", "\\3" , extIndDF_m$variable)#names(datac)[-1])
  
  outData <- data.frame(gcm = gcm, SRinfo = as.numeric(SRinfo), QEthr = extIndDF_m$QEthr, yrs2ext = extIndDF_m$value, extYr = extYrDF_m$value)

  time2extSR <- ggplot(data = outData, aes(x = QEthr, y = SRinfo, z = yrs2ext)) +
    geom_tile(aes(fill = yrs2ext)) +
    scale_fill_gradient2(low = "red", mid = "white", high = "black", midpoint = 45) +
    #stat_contour(aes(colour = ..level..), binwidth = 10) +
    facet_grid(gcm~.) +
    xlab("Quasi-extinction threshold (numbers)") +
    ylab("Equilibrium abundance of spawners") +
    theme_classic() +
    scale_x_continuous(expand=c(0,0)) + 
    scale_y_continuous(expand=c(0,0)) 
  
  print(time2extSR)
  
  extYrSR <- ggplot(data = outData, aes(x = QEthr, y = SRinfo, z = extYr)) + 
    geom_tile(aes(fill = extYr)) +
    scale_fill_gradient2(low="red", mid = "white", high = "black", midpoint = 2055) +
    facet_grid(gcm~.) +
    xlab("Quasi-extinction threshold (numbers)") +
    ylab("Equilibrium abundance of spawners") +
    theme_classic() +
    scale_x_continuous(expand=c(0,0)) + 
    scale_y_continuous(expand=c(0,0)) 
  print(extYrSR)
  
  return(outData)
  
  
} # end calcQEtimeConstEQvQE



QE_counter <- function(data, qeLev = 100, run_length = 4) {
  
  extInd <- vector("integer", ncol(data))
  for (h in 1:ncol(data)) {  
    qeInd <- which(data[,h] <= qeLev)
    extRuns <- vector("integer", 0)
    k <- 1
    if (length(qeInd) >= run_length) {
      for (i in 1:(length(qeInd)-3)) {
        if (qeInd[i] == (qeInd[i+1] - 1) & qeInd[i] == (qeInd[i+2] - 2) & 
              qeInd[i] == (qeInd[i+3] - 3)  ) {
          extRuns[k] <- (qeInd[i] + 3)
          k <- k + 1
        }  
      }
    } else {
      extRuns[k] <- (nrow(data) + 1) 
      k <- k+1
    }
    if(length(extRuns > 0)) extInd[h] <- min(extRuns) else extInd[h] <- NA
  }
  return(extInd)
} # end QE_count

# A better QE counter ####
# courtesy Jaime Ashander

JA_consec <- function(z, qeLev = 100, run_length = 4) {
  i.low <- which(z <= qeLev)
  
  if (length(i.low) <= run_length) return(as.integer(NA))
  
  i.last <- i.low[1]
  counter <- 1
  for( i in i.low[-1]){

    if(counter == run_length)  return(as.integer(i.last)) #  need the double to ensure that data.table gets what it expects
    
    if(i-1 == i.last)  counter <- counter + 1 else counter <- 1
    
    if(i == i.low[length(i.low)] & counter < run_length) i.last <- NA else i.last <- i
  }
  if (is.null(class(i.last))) print("error here")
  return(as.integer(i.last))
}
#if (length(which(z < qeLev)) < run_length)  return(length(z)+1) else i.low <- which(z < 100)


# Stochastic population projection models ####

# create random time series
# white noise
mk_white <- function(N) {
  return(rep(1/N, N/2))
}

# create "band-pass" style time series
mk_rsin <- function(N, lowF, highF) {
  temp <- (1:(N/2))/(N)
  ns <- length(which(temp <= highF & temp >= lowF))
  # print(ns)
  goodInd <- which(temp <= highF & temp >= lowF)
  temp[-goodInd] <- 0
  temp[goodInd] <- 1/ns
  #print(temp)
  return(temp)  }

# create "band-reject" style time series

mk_rsin_reject <- function(N, lowF, highF) {
  temp <- (1:(N/2))/(N)
  ns <- length(which(temp >= highF & temp <= lowF))
  # print(ns)
  goodInd <- which(temp >= highF & temp <= lowF)
  temp[-goodInd] <- 0
  temp[goodInd] <- 1/ns
  #print(temp)
  return(temp)  
}

# Wavelet filtering of time series to do "band-pass" style filtering
# do band pass filter on the wavelet decomposition of a white noise series.
# This reconstructs the signal using only the pd 3 to 4 components of the
# wavelet power spectrum (T&C 1998 equation 11)
waveFlt34 <- function(white_n) {
  # Wavelet filtering pd 3-4 signal from white noise signal above
  # White noise filtering out pd 3-4 variability from the wavelet power spectrum 
  # vs generating period 3-4 noise using the random sine wave method
  
  # Create a "time" vector
  N <- nrow(white_n)
  x <- 1:N
  
  # I copied/translated Flora's matlab code (Will, Matt, Lauren, T&C) for J1
  J1 <- trunc((log(32/(2 * 1))/log(2))/0.01) # number of scales minus - 1
  # Torrence and Compo 1998 equation 11 constants for reconstruction of time series for
  # Morlet mother wavelets using a delta function 
  # Morlet wavelet constants
  C_delta <- 0.776 #  The factor C_delta comes from the reconstruction of a
  # delta function from its wavelet transform using the function Psi_0(eta)
  Psi_0 <- pi^(-0.25) # removes the energy scaling
  dj <- 0.01
  dt <- 1
  
  reCon <- (dj*dt^(0.5))/(C_delta*Psi_0) # Constant for reconstruction of wavelets
  
  library(biwavelet) # R package for continuous wavelet analysis based on T & C 1998
  
  wave34_n <- matrix(NA, nrow = nrow(white_n), ncol = ncol(white_n))
  
  for (i in 1:ncol(white_n)) {
    white.wt <- wt(cbind(x, white_n[,i]), dj = 0.01, J1 = J1, max.scale = 32, mother = "morlet", sig.test = 0, sig.level = 0.95)
    pd34ind <- which(white.wt$period >= 3 & white.wt$period <= 4)
    # new_y <- as.numeric(scale(colMeans(white.wt$power.corr[pd34ind,])))
    # new.wt <- wt(cbind(x,new_y), dj = 0.01, J1 = J1, max.scale = 32, mother = "morlet", sig.test = 0, sig.level = 0.95)
    
    sumWave <- matrix(NA, nrow = nrow(white.wt$wave),
                      ncol = ncol(white.wt$wave))
    
    for (j in 1:nrow(white.wt$wave)) {
      sumWave[j,] <- Re(white.wt$wave[j,])/(white.wt$scale[j]^(0.5))
      # scaling by s^(1/2) converts the wavelet transform to an energy density
    }
    
    # reCon_y <- reCon*colSums(sumWave) # reconstruct the original white noise signal
    sumWave34 <- sumWave[pd34ind,]
    temp <- reCon*colSums(sumWave34)
    wave34_n[,i] <- (temp - mean(temp))/sd(temp)
  }
  return(wave34_n)
} # end of waveFlt34()


popSimPSvary <- function(rand_surv, surv1, surv2, surv3, EQsp, wanted_frac, alpha_scale) {
  # popSimPSvary() is a wrapper for PopProjPSvar().
  # PopProjPSvar() runs a single stochasitic run of the population model with time varying survPS
  # for this specified model parameters: surv1, surv2, surv3, EQsp, wanted_frac, alpha_scale
  # popSimPSvary() just creates output containers (outPop: contains abundance at age for each year;
  # and outSpawn: contains the number of age_3 and age_4 critters surviving to spawn) to collect the output
  # of the replicate runs [for each columon of the rand_surv input]
  outPop <- array(NA, c(reps, 4, nrow(rand_surv)))
  outSpawn <- matrix(NA, ncol = reps, nrow = nrow(rand_surv))
  
  for (i in 1:ncol(rand_surv)) {
    outList <- PopProjPSvar(survPS = rand_surv[,i], 
                            surv1 = surv1, 
                            surv2 = surv2, 
                            surv3 = surv3,
                            EQsp = EQsp, 
                            wanted_frac = wanted_frac,
                            alpha_scale = alpha_scale)
    outPop[i,,] <- outList[[1]]
    outSpawn[,i] <- outList[[2]]
  }
  return(list(outPop, outSpawn))
} # end of popSimPSvary()

PopProjPSvar <- function(survPS, surv1, surv2, surv3, EQsp, wanted_frac, alpha_scale) {
  # This function projects model BC spring run chinook salmon (females) populations with either constant 
  # or variable survival for a !single run!
  
  # USe popSimPSvary() to run multiple simulations and collect output for analysis
  
  # PopProjPSvar() runs a single stochasitic run of the population model with time varying survPS
  # for this specified model parameters: surv1, surv2, surv3, EQsp, wanted_frac, alpha_scale
  
  # survPS = either constant or variable 
  # surv1, surv2 and surv3 = constant. surv1 ~ early ocean survival (0.01-0.05). surv2 == surv3 (0.8 high)
  # EQsp = the number of age-3 and age-4 females returning to spawn in Butte Creek, susceptible to prespawning 
  # mortality
  # wanted_frac = the ratio of age-3 spawners to total spawners (age-3 + age-4), used to determine the fraction 
  # of age-3 fish in the ocean returning to spawn each year
  # alpha_scale = a multiplier of the 1/spr to get the slope of the Beverton-Holt stock recruitment curve
  
  # May 20: check that function arguments are of length == 1
  if ( length(surv1) != 1 | length(surv2) != 1 | length(surv3) != 1 | 
         length(EQsp) != 1 | length(wanted_frac) != 1 | length(alpha_scale) != 1) stop("Function arguments survPS, surv1, surv2, 
                                                                                       surv3, EQsp, wanted_frac, 
                                                                                       be length 1! Doofus...")
  
  delta_e <- calc_de(wanted_frac, surv3)
  # assuming survPS = 1 produces maximum spr, 
  # hence least steep replacement line
  # use mean(survPS, na.rm = ...)
  if (length(survPS) > 1) { 
    spr <- SPR_srcs(surv1, surv2, surv3, delta_e, survPS=1) 
  } else {
    spr <- SPR_srcs(surv1, surv2, surv3, delta_e, survPS=survPS) 
  } 
  alpha <- 1/spr*alpha_scale 
  beta <- calc_beta(alpha, EQsp, spr)  
  
  # Set up an initial age vector
  # start from equilibrium age structure in ocean just before returning to spawn
  
  n0 <- numeric(4) #matrix(0, nrow = 4, ncol = 1)
  
  # There will be an n0 for a given EQsp, 
  # but have to calculate from oldest age
  # to youngest age
  
#   Test code to get initial abundance at age vector from reduced EQ based on reduced meanPS level
#   ---
#   spr_loweq <- SPR_srcs(surv1, surv2, surv3, delta_e, survPS = PSmean) 
#   
#   init_eq <- (alpha*spr_loweq - 1)/beta # equilibrium number of spawners
#   n0[1] <- (init_eq*(1-wanted_frac)) # age-4
  # ---
  n0[1] <- (EQsp*(1-wanted_frac)) # age-4
  n0[2] <- n0[1]/(surv3*(1-delta_e))
  n0[3] <- n0[2]/surv2
  n0[4] <- n0[3]/surv1
  
  n0 <- rev(n0) # This flips the n0 vector, so that age-1 is now in the 
  # first position rather than the 4th position 
  # Storage array for population simulations
  # 1-d for each age (4)
  # 1-d for time (nyrs)
  # so need a 4 X nyrs dimension MATRIX for pop and 
  # a VECTOR nyrs long for spawners
  
  # Create arrays for storing model output
  #burn_in <- trunc(length(envF)/3)
  nyrs <- length(survPS) #burn_in + length(envF)
  #envF <- c(sample(envF, size = trunc(length(envF)/3), replace = TRUE), envF)
  pop <- matrix(0, nrow = length(n0), ncol = nyrs) #array(0, c(4, n,  reps))
  sp <-  numeric(nyrs) # array(0, c(n, reps))
  
  for (t in 1:nyrs) {
    if (t == 1) {   
      pop[1, t] <- n0[1]
      pop[2, t] <- n0[2]
      pop[3, t] <- n0[3]
      pop[4, t] <- n0[4]
    } else {
      sp[t-1] <- pop[3, t-1]*delta_e*(survPS[t-1]) + pop[4,t-1]*(survPS[t-1]) 
      pop[1,t] <- (sp[t-1]*alpha)/(1+sp[t-1]*beta)  
      pop[2,t] <- pop[1, t-1]*surv1
      pop[3,t] <- pop[2, t-1]*surv2
      pop[4,t] <- pop[3, t-1]*surv3*(1-delta_e)
    }
    sp[t] <- pop[3, t]*delta_e*(survPS[t]) + pop[4, t]*(survPS[t])   
  }  
  
  return(list(pop = pop, sp = sp))
  
} # end PopProjPSvar()

##### beta distribution approach #####

range01 <- function(x){(x-min(x))/(max(x)-min(x))}

# a fucntion for generation correlated RV from two normal distributions
# http://www.sitmo.com/article/generating-correlated-random-numbers/
# http://r.789695.n4.nabble.com/generate-two-sets-of-random-numbers-that-are-correlated-td3736161.html
corrNorm <- function(n, rho, seed) { 
  
  if ( rho < 0 | rho > 1 ) stop("rho must be between 0 and 1")
  set.seed(seed)
  X1 <- rnorm(n, 0, 1)
  X2 <- rnorm(n, 0, 1) 
  Z <- rho*X1+sqrt(1-rho^2)*X2
  return(list(X1, X2, Z) )
  
} 

# function for ar-1 noise translated from Matt Holland's Matlab code

ar1rand <- function(n, mu, sigma, phi, seed) {
  # R version of Matt Hollands Matlab function for making
  # AR(1) noise
  # calculate AR(1) model parameters
  c <- mu*(1-phi)
  sigma_xi <- sqrt(sigma^2*(1-phi^2))
  set.seed(seed)
  xi <- sigma_xi * rnorm(n, 0, 1)
  # initialize y and define first element
  y <- numeric(n)
  y[1] <- c + phi*mu + xi[1]
  # iterate over the remaining samples of xi, assigning subsequent
  # values of y
  for (i in 2:n) {
    y[i] <- c + phi * y[i-1] + xi[i]
  }
  return(y)
}

ar1_redden <- function(ts, mu, sigma, phi) {
  # Variant of R version of Matt Hollands Matlab function for making
  # AR(1) noise
  # This one reddens an existing time series
  
  # calculate AR(1) model parameters
  n <- length(ts)
  n_phi <- length(phi)
  c <- mu*(1-phi)
  sigma_xi <- sqrt(sigma^2*(1-phi^2))
  
  # initialize y and define first element
  y <- xi <- matrix(0, nrow = n, ncol = n_phi)
  
  for (k in 1:n_phi) {
    xi[,k] <- sigma_xi[k] * ts#rnorm(n, 0, 1)
    y[1,k] <- c[k] + phi[k]*mu + xi[1,k]
    # iterate over the remaining samples of xi, assigning subsequent
    # values of y
    for (i in 2:n) {
      y[i,k] <- c[k] + phi[k] * y[i-1,k] + xi[i,k]
    }  
  }
  colnames(y) <- paste("R=", phi, sep="")
  return(y)
}  

# calculate spawners per recruit (1/spr = replacement line, baseline/floor for alpha)
SPR_srcs <- function(surv1, surv2, surv3, delta_e, survPS = 1) {
  SPR <- surv1*surv2*delta_e*survPS + surv1*surv2*surv3*(1-delta_e)*survPS
  return(SPR)
} # end bracket for SPR_scrs

# plot how SPR changes with early ocean survival and fraction of age-3 spawners
plotSPR_de_surv1 <- function (sprDF) {
  ggplot(sprDF, aes(x = delta_e, y = SPR, colour = factor(surv1))) +
    geom_point(size = 3) +
    scale_colour_manual(values = c("black","grey20", "grey40", "grey70"), name = "Early Ocean \nSurvival") +
    ylab("Spawners Per Recruit") +
    xlab("Fraction of spawners returning at age-3") +
    xlim(c(0, 1)) +
    ylim(c(0, 0.08)) +
    theme_bw() +
    theme(panel.grid.major = theme_line(colour = NA), panel.grid.minor = theme_line(colour = NA)) 
}

Zuur_acm_univar <- function(vector, path4file) {
  
  df <- data.frame(years = 2010:2099, surv = vector, period=gl(n = 3, k = 30, labels = c("Early", "Middle", "Late")))
  
  #Makes 9 diagnostic plots of SALMOD output data 
  pdf(file = path4file)
  
  # Outliers
  # 1a Boxplot 
  
  df_sd <- aggregate(df$surv, by = list(df$period), sd)
  df_mn <- aggregate(df$surv, by = list(df$period), mean)
  boxplot <- ggplot(data = df, aes(x = period, y = surv)) +
    geom_jitter(colour = "slateblue", alpha = .8, size = 2.5) +
    geom_boxplot(outlier.colour = "red", alpha=.8, fill = "grey", outlier.size = 3, outlier.shape = 3 ) + 
    theme_classic() +
    ylab("Survival") +
    xlab("Time period\n(2010-2039, 2040-2069, 2070-2099)") +
    annotate("text", x = df_sd$Group.1, y =0.9, label = paste("sd =",round(df_sd$x, 2)), size = 6,
             colour = "red") +
    annotate("text", x = df_mn$Group.1, y =0.1, label = paste("mean =",round(df_mn$x, 2)), size = 6,
             colour = "red")
  
  print(boxplot)
  
  
  # 1b Clevelend dotplot 
  
  dotchart(df$surv)
  # conditioned on "period"  
  dotchart(df$surv, groups = df$period)
  
  # Homegeneity
  # 2. conditional boxplot
  
  # See above
  
  # Normality
  # 3a. Histogram geom_dotplot
  
  histogram <- ggplot(data = df, aes(x = surv)) +
    geom_histogram(aes(y = ..density..), fill = "grey30", binwidth = 0.1) +
    geom_density(size = 1.5, fill = "red", alpha = .5) +
    theme_few() +
    ylab("Frequency") +
    xlab("Survival")
  
  print(histogram)
  
  # 3b. QQ-plot
  
  qqnorm(vector)
  
  # Zero-trouble
  # 4. frequency plot
  
  dotplot <-  ggplot(data = df, aes(x = surv)) + 
    geom_dotplot(binwidth = 0.05) +
    theme_bw() +
    ylab("Frequency") +
    xlab("Survival")
  print(dotplot)    
  
  cond_dotplot <-  ggplot(data = df, aes(x = surv, fill = period)) + 
    geom_dotplot(binwidth = 0.05, stackgroups = TRUE, method = "histodot") +
    scale_fill_tableau("colorblind10") + 
    theme_bw() +
    ylab("Frequency") +
    xlab("Survival") 
  print(cond_dotplot)
  
  # Interactions
  
  coplots <- ggplot(data = df, aes(x = years, y = surv, colour = period)) +
    geom_line(size = 0.4, linetype = 2) +
    scale_colour_tableau("colorblind10") +
    geom_point(size = 3) +
    theme_few() +
    ylab("Survival") +
    xlab("Year")
  
  print(coplots)
  
  # Independence
  
  acf(df$surv)
  
  dev.off()  
} # end of Zuur_acm_univar()

salmod_surv_dp <- function() {
  sim_path <- file.path(".", "data", "simulation_results")
  mgmtScenarios <- c("BAU", "NoDiversion") 
  #"ColdWater", "ForeCast", "ForecastError", "NoDiversion", 
  #"RaisePhilbrook", "RaisePhilbrook1710", "ReservoirCover", "Shade")
  climateScenarios <- c("A2","B1")
  gcms <- c("cnrmcm3", "gfdlcm21", "miroc32med", "mpiecham5", "ncarccsm3", "ncarpcm1")
  gcmScenarios <- numeric(length(climateScenarios)*length(gcms))
  stage <- c("prespawn", "egg", "fry")
  
  for (m in 1:length(mgmtScenarios))  {
    for (c in 1:length(climateScenarios)) {
      
      for (k in 1:length(gcms)) {
        salmod_out <- read.delim(file.path(sim_path, mgmtScenarios[m], paste0(climateScenarios[c],"_", gcms[k]), "SALMODsumOutMerge.txt"), sep = ",")
        mortFW <- subset(salmod_out, select = c("AFem", "allMortSF", "Eggs", "FryGrad" ,"FryExit"))
        rm(salmod_out)
        mortFW$gcms <- gcms[k]
        
        if(k == 1) mortFWgcms <- mortFW else mortFWgcms <- rbind(mortFWgcms, mortFW)
      } # end of k loop
      
      prespSurv <-  with(mortFWgcms, data.frame(year = rep(2010:2099, time = 6), gcms = gcms, PreSpawn = allMortSF/AFem))
      eggSurv <- with(mortFWgcms, data.frame(year = rep(2010:2099, time = 6), gcms = gcms, Egg = FryGrad/Eggs))
      frySurv <- with(mortFWgcms, data.frame(year = rep(2010:2099, time = 6), gcms = gcms, Fry = FryExit/FryGrad))
      prespSurvC <- dcast(prespSurv, year ~ gcms, value.var = "PreSpawn")
      eggSurv$Egg[which(is.na(eggSurv$Egg))] <- 0
      eggSurvC <- dcast(eggSurv, year ~ gcms, value.var = "Egg")
      frySurv$Fry[which(is.na(frySurv$Fry))] <- 0
      frySurvC <- dcast(frySurv, year ~ gcms, value.var = "Fry")
      
      for ( g in 2:7 ) {
        for ( s in 1:3 ) {    
          if (!dir.exists(file.path(".", "output", "diagnostic_plots", mgmtScenarios[m], climateScenarios[c]))) { 
            dir.create(file.path(".", "output", "diagnostic_plots", mgmtScenarios[m], climateScenarios[c]), recursive = TRUE) }
          dp_path <- file.path(".", "output", "diagnostic_plots", mgmtScenarios[m], climateScenarios[c])
          path4file <- file.path(dp_path, paste0(gcms[g-1], "_", stage[s], ".pdf")) 
          ifelse(stage[s] == "prespawn", vector <- prespSurvC[,g],
                 ifelse(stage[s] == "egg", vector <- eggSurvC[,g],
                        ifelse(stage[s] == "fry", vector <- frySurvC[,g])) )
          
          Zuur_acm_univar(vector, path4file)
          
        }
      }
    }
  }
  
} # end of salmod_surv_dp


wavelet_SALMOD_series <- function(management = c("BAU", "NoDiversion"),
                                  emissions = c("A2","B1"),
                                  global_mods = c("cnrmcm3", "gfdlcm21", "miroc32med", "mpiecham5", "ncarccsm3", "ncarpcm1") 
) {
  
  # Function that extracts a specific survival time series (Prespawn, egg and fry) 
  # from SALMOD output files for each climate and water management scenario.
  # Then it prints a pdf with plots of the ts, WPS, global WPS and scale-averaged WPS
  
  sim_path <- file.path(".", "data", "simulation_results")
  
  mgmtScenarios <- management #"ColdWater", "ForeCast", "ForecastError", "NoDiversion", 
  #"RaisePhilbrook", "RaisePhilbrook1710", "ReservoirCover", "Shade")
  climateScenarios <- emissions
  gcms <- global_mods 
  gcmScenarios <- numeric(length(climateScenarios)*length(gcms))
  
  J1 <- trunc((log(32/(2 * 1))/log(2))/0.01) # number of scales minus - 1; I translated from your matlab code
  
  for (m in 1:length(mgmtScenarios))  {
    for (c in 1:length(climateScenarios)) {
      
      for (k in 1:length(gcms)) {
        salmod_out <- read.delim(file.path(sim_path, mgmtScenarios[m], 
                                           paste0(climateScenarios[c],"_", gcms[k]), "SALMODsumOutMerge.txt"), 
                                 sep = ",")
        mortFW <- subset(salmod_out, select = c("AFem", "allMortSF", "Eggs", "FryGrad" ,"FryExit"))
        mortFW$gcms <- gcms[k]
        
        if(k == 1) mortFWgcms <- mortFW else mortFWgcms <- rbind(mortFWgcms, mortFW)
      } # end of k loop
      
      prespSurv <-  with(mortFWgcms, data.frame(year = rep(2010:2099, time = 6), gcms = gcms, PreSpawn = allMortSF/AFem))
      eggSurv <- with(mortFWgcms, data.frame(year = rep(2010:2099, time = 6), gcms = gcms, Egg = FryGrad/Eggs))
      frySurv <- with(mortFWgcms, data.frame(year = rep(2010:2099, time = 6), gcms = gcms, Fry = FryExit/FryGrad))
      
      prespSurvC <- dcast(prespSurv, year ~ gcms, value.var = "PreSpawn")
      eggSurv$Egg[which(is.na(eggSurv$Egg))] <- 0
      eggSurvC <- dcast(eggSurv, year ~ gcms, value.var = "Egg")
      frySurv$Fry[which(is.na(frySurv$Fry))] <- 0
      frySurvC <- dcast(frySurv, year ~ gcms, value.var = "Fry")
      
      # For each cliamte scenario and water management alternative, plot the
      # the time series and wavelet spectra for each gcm for
      # prespawn survival, egg survival and fry survival
      # Calculate wavelet spectra
      # Visualize spectra
      dims <- dim(prespSurvC)
      out_wavelet <- file.path(".", "output", "wavelet_plots")
      if (!dir.exists(out_wavelet)) { dir.create(out_wavelet, recursive = TRUE) }
    
      path4file <- file.path(out_wavelet, paste0(mgmtScenarios[m], "_", climateScenarios[c], ".pdf"))
      pdf(file = path4file)
      
      for (i in 2:dims[2]) {
        #i <- 2
        # Time series
        par(fig= c(0,0.65,0.6,1))
        old <- par(mar = c(3, 3, 3, 1) )
        plot(2010:2099, prespSurvC[,i], type = "l", lwd = 2, col = "darkgrey",
             axes = FALSE, ylab = "Survival Rate", xlab = "Year", xlim = c(2010, 2100), xaxs = "i")
        axis(1, at = c(2010, 2020, 2030, 2040, 2050, 2060, 2070, 2080, 2090, 2100), 
             lab = c(2010, 2020, 2030, 2040, 2050, 2060, 2070, 2080, 2090, 2100), 
             las = 1, cex.axis = 1, tck = 0.02)
        axis(2, at = c(0, 0.2, 0.4, 0.6, 0.8, 1.0), 
             lab =  c(0, 0.2, 0.4, 0.6, 0.8, 1.0), 
             las = 2, cex.axis = 1, tck = 0.02)
        box(lwd = 2)
        title(paste("GCM: ", gcms[i-1], " Prespawn survival emissions scenario: ", climateScenarios[c], " \n Water managament alternative: ", mgmtScenarios[m], sep = ""))
        par(old)
        
        # Wavelet Power Spectrum
        par(fig= c(0,0.65,0.3,0.6), new = TRUE)
        old <- par(mar = c(3, 3, 1, 1) )
        test_cwt <- wt(cbind(2010:2099, detrend(as.numeric(scale(prespSurvC[,i])))), dj = 0.01, J1 = J1, max.scale = 32, mother = "morlet", sig.test = 0, sig.level = 0.95)
        plot(test_cwt, type = "power.corr.norm", xaxt = 'n', xlim = c(2010, 2100))
        axis(1, at = c(2010, 2020, 2030, 2040, 2050, 2060, 2070, 2080, 2090, 2100), 
             lab = c(2010, 2020, 2030, 2040, 2050, 2060, 2070, 2080, 2090, 2100), 
             las = 1, cex.axis = 1, tck = 0.02)
        par(old)
        # Scale-averaged Power Spectrum - periods of greatest variability
        par(fig= c(0,0.65,0,0.3), new = TRUE)
        old <- par(mar = c(3, 3, 1, 1) )
        plot(2010:2099, colSums(test_cwt$power.corr), type = "l", lwd = 2, xlab = "Years", ylab = "Power")
        par(old)
        # Global Wavelet Power Spectrum
        par(fig= c(.65,1,0.3,.6), new = TRUE)
        old <- par(mar = c(3, 3, 1, 1) )
        yrange <- NULL 
        y_ticks <- 2^(floor(log2(min(test_cwt$period, yrange))):(floor(log2(max(test_cwt$period, yrange)))+1))
        FR <- log(rowSums(test_cwt$power.corr))
        plot(FR, log2(test_cwt$period), type = "l", lwd = 2, ylim = rev(range(log2(test_cwt$period))), axes = FALSE)
        axis(1)
        axis(2,  log2(y_ticks[length(y_ticks):1]), y_ticks[length(y_ticks):1])
        box()
        par(old)
      }
      for (i in 2:dims[2]) {
        #i <- 2
        par(fig= c(0,0.65,0.6,1))
        old <- par(mar = c(3, 3, 3, 1) )
        # plot egg surv ts
        plot(2010:2099, eggSurvC[,i], type = "l", lwd = 2, col = "darkgrey",
             axes = FALSE, ylab = "Survival Rate", xlab = "Year", xlim = c(2010, 2100), xaxs = "i")
        axis(1, at = c(2010, 2020, 2030, 2040, 2050, 2060, 2070, 2080, 2090, 2100), 
             lab = c(2010, 2020, 2030, 2040, 2050, 2060, 2070, 2080, 2090, 2100), 
             las = 1, cex.axis = 1, tck = 0.02)
        axis(2, at = c(0, 0.2, 0.4, 0.6, 0.8, 1.0), 
             lab =  c(0, 0.2, 0.4, 0.6, 0.8, 1.0), 
             las = 2, cex.axis = 1, tck = 0.02)
        box(lwd = 2)
        title(paste("GCM: ", gcms[i-1], " Eggs survival emissions scenario: ", climateScenarios[c], " \n Water managament alternative: ", mgmtScenarios[m], sep = ""))
        par(old)
        # plot wavelet power spectrum of egg survival
        par(fig= c(0,0.65,0.3,0.6), new = TRUE)
        old <- par(mar = c(3, 3, 1, 1) )
        test_cwt <- wt(cbind(2010:2099, detrend(as.numeric(scale(eggSurvC[,i])))), dj = 0.01, J1 = J1, max.scale = 32, mother = "morlet", sig.test = 0, sig.level = 0.95)
        plot(test_cwt, type = "power.corr.norm", xaxt = 'n', xlim = c(2010, 2100))
        axis(1, at = c(2010, 2020, 2030, 2040, 2050, 2060, 2070, 2080, 2090, 2100), 
             lab = c(2010, 2020, 2030, 2040, 2050, 2060, 2070, 2080, 2090, 2100), 
             las = 1, cex.axis = 1, tck = 0.02)
        par(old)
        # Scale-averaged Power Spectrum - periods of greatest variability
        par(fig= c(0,0.65,0,0.3), new = TRUE)
        old <- par(mar = c(3, 3, 1, 1) )
        plot(2010:2099, colSums(test_cwt$power.corr), type = "l", lwd = 2, xlab = "Years", ylab = "Power")
        par(old)
        # Global Wavelet Power Spectrum
        par(fig= c(.65,1,0.3,0.6), new = TRUE)
        old <- par(mar = c(3, 3, 1, 1) )
        yrange <- NULL 
        y_ticks <- 2^(floor(log2(min(test_cwt$period, yrange))):(floor(log2(max(test_cwt$period, yrange)))+1))
        FR <- log(rowSums(test_cwt$power.corr))
        plot(FR, log2(test_cwt$period), type = "l", lwd = 2, ylim = rev(range(log2(test_cwt$period))), axes = FALSE)
        axis(1)
        axis(2,  log2(y_ticks[length(y_ticks):1]), y_ticks[length(y_ticks):1])
        box()
        par(old)
      }
      
      for (i in 2:dims[2]) {
        #i <- 2
        # plot fry survival time series 
        par(fig= c(0,0.65,0.6,1))
        old <- par(mar = c(3, 3, 3, 1) )
        plot(2010:2099, frySurvC[,i], type = "l", lwd = 2, col = "darkgrey",
             axes = FALSE, ylab = "Survival Rate", xlab = "Year", xlim = c(2010, 2100), xaxs = "i")
        axis(1, at = c(2010, 2020, 2030, 2040, 2050, 2060, 2070, 2080, 2090, 2100), 
             lab = c(2010, 2020, 2030, 2040, 2050, 2060, 2070, 2080, 2090, 2100), 
             las = 1, cex.axis = 1, tck = 0.02)
        axis(2, at = c(0, 0.2, 0.4, 0.6, 0.8, 1.0), 
             lab =  c(0, 0.2, 0.4, 0.6, 0.8, 1.0), 
             las = 2, cex.axis = 1, tck = 0.02)
        box(lwd = 2)
        title(paste("GCM: ", gcms[i-1], " Fry survival emissions scenario: ", climateScenarios[c], " \n Water managament alternative: ", mgmtScenarios[m], sep = ""))        
        par(old)
        # plot Wavelet power spectrum
        par(fig= c(0,0.65,0.3,0.6), new = TRUE)
        old <- par(mar = c(3, 3, 1, 1) )
        test_cwt <- wt(cbind(2010:2099, detrend(as.numeric(scale(frySurvC[,i])))), dj = 0.01, J1 = J1, max.scale = 32, mother = "morlet", sig.test = 0, sig.level = 0.95)
        plot(test_cwt, type = "power.corr.norm", xaxt = 'n', xlim = c(2010, 2100))
        axis(1, at = c(2010, 2020, 2030, 2040, 2050, 2060, 2070, 2080, 2090, 2100), 
             lab = c(2010, 2020, 2030, 2040, 2050, 2060, 2070, 2080, 2090, 2100), 
             las = 1, cex.axis = 1, tck = 0.02)
        par(old)
        # Scale-averaged Power Spectrum - periods of greatest variability
        par(fig= c(0,0.65,0,0.3), new = TRUE)
        old <- par(mar = c(3, 3, 1, 1) )
        plot(2010:2099, colSums(test_cwt$power.corr), type = "l", lwd = 2, xlab = "Years", ylab = "Power")
        par(old)
        # Global Wavelet Power Spectrum
        par(fig= c(.65,1,0.3,.6), new = TRUE)
        old <- par(mar = c(3, 3, 1, 1) )
        yrange <- NULL 
        y_ticks <- 2^(floor(log2(min(test_cwt$period, yrange))):(floor(log2(max(test_cwt$period, yrange)))+1))
        FR <- log(rowSums(test_cwt$power.corr))
        plot(FR, log2(test_cwt$period), type = "l", lwd = 2, ylim = rev(range(log2(test_cwt$period))), axes = FALSE)
        axis(1)
        axis(2,  log2(y_ticks[length(y_ticks):1]), y_ticks[length(y_ticks):1])
        box()
        par(old)
      }
      
      
      dev.off()
      
    } # end of c climate emissions loop
  } # end of m water management loop
  
} # end of wavelet_SALMOD_series()

PDO_wvlt_test <- function() {
  # sim_path <- file.path(".", "data", "simulation_results")
  # setwd("~/NEP_salmon/chap_3/data/PDO/")
  # 
  # Testing the biwavelet function to show Loo 
  # Use PDO time series 
  pdo <- read.table(url("http://research.jisao.washington.edu/pdo/PDO.latest.txt"))
                 # row.names = FALSE,
                 # stringsAsFactors = FALSE,
                 # sep = "\t",
                 # skip = 34)
  pdo <- readLines(url("http://research.jisao.washington.edu/pdo/PDO.latest.txt"))
  lpdo <- length(pdo)
  pdo <- (pdo[-c(1:34, (lpdo-20):lpdo)])
  columns <- c("YEAR", "JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", "SEP", "OCT", "NOV", "DEC")
  pdo_df <- data.frame(matrix(ncol = length(columns), nrow = 0))
  
  for (i in 1:length(pdo) ) {
    df <- read.table(textConnection(pdo[[i]]))
    pdo_df <- rbind(pdo_df, df)
  }
  colnames(pdo_df) <- columns
  
  # pdo <- read.table("PDO_2012SEP.txt", header = TRUE, sep = "\t", fill = TRUE)
  rowMeans(pdo[,2:13], na.rm = TRUE)
  # wavelet analysis of the PDO annual signal
  ann_pdo <- scale(rowMeans(pdo[,2:13], na.rm = TRUE)) # mean of 0, sd of 1, drop year 2012 - since incomplete
  years <- pdo[,1] # drop 2012 since incomplete
  J1 <- trunc((log(32/(2 * 1))/log(2))/0.01) # number of scales minus - 1; I translated from your matlab code
  # Calculate the wavelet, using parameters as in your matlab code
  pdo.wt <- wt(cbind(years, ann_pdo), dj = 0.01, J1 = J1, max.scale = 32, 
               mother = "morlet", sig.test = 0, sig.level = 0.95) # acf(ann_pdo)$acf[2]
  
  # Compare to Flora's matlab: April 24
  # outputs from wavelet function call on PDO through 2012. Inputs in R and Matlab are the same
  # Wavelet coefficients: look the same
  # Period: same
  # Scale: same
  # Wavelet Power: different. matlab code rescales wavelet power spectrum to have units of variance.
  # R {biwavelet} also calculates bias corrected power
  # COI: same
  # signif: different. to get signif R divides wavelet power spectrum by product of var(t.s.) and signif (Pk), 
  # while matlab only divides by signif.
  # Plot time series and the wavelet power spectrum
  old <- par(mfrow = c(2,1), mar = c(3,4,0.5,0.5))
  plot(years, ann_pdo, type = "l", lwd = 3.7, lty = 1, col = "black", xlab = "Year", ylab = "PDO")
  lines(years, ann_pdo, type = "l", lwd = 2, col = "white")
  lines(years, ann_pdo, lwd = 2, lty = 2)
  box(lwd = 2)
  plot(pdo.wt, type = "power.corr.norm", tol = 0.95)
  box(lwd = 2)
  par(old)
  
} # end of PDO_wvlt_test()

# Plot global wavelt power spectrum

wt_plot_GWPS <- function(biwv_obj, labs = TRUE, small = TRUE) {
  library(biwavelet)
  if (small == TRUE) vold <- par(mar = c(1,2,1,1), cex = .7) else vold <- par(mar = c(4,4,0.5,0.5))
  yrange <- NULL 
  y_ticks <- 2^(floor(log2(min(biwv_obj$period, yrange))):(floor(log2(max(biwv_obj$period, yrange)))+1))
  FR <- log(rowMeans(biwv_obj$power.corr))
  if (labs == TRUE) {
    xlabel <- "Relative Magnitude"
    ylabel <- "Period"
  } else {
    xlabel <- ""
    ylabel <- ""
  }
  plot(FR, log2(biwv_obj$period), type = "n", ylim = rev(range(log2(biwv_obj$period))), 
       axes = FALSE, ylab = ylabel, xlab = ylabel)
  lines(FR, log2(biwv_obj$period), lty = 1, lwd = 4)
  lines(FR, log2(biwv_obj$period), lty = 1, lwd = 3, col = "white")
  lines(FR, log2(biwv_obj$period), lty = 2, lwd = 2.5)
  axis(1)
  axis(2,  log2(y_ticks[length(y_ticks):1]), y_ticks[length(y_ticks):1])
  box(lwd = 2)
  par(vold)
}

summaryCh3tsPlot <- function(noise = noiseList, fish = spawners, n = 1) {
  
  library(biwavelet)
  
  ylim = c(0,2)
  lwd_ts <- 1.5
  J1 <- trunc((log(32/(2 * 1))/log(2))/0.01)
  
  # pdf()
  old <- par(mar = c(1,2,1,1), cex = .7)
  
  # white noise
  # Generating spectrum
  par(fig=c(0, 0.2, 0.9, 1))
  plot(x = seq(0,0.5, length = 50), y = runif(50, min = 0, max = 2), type = "n",
       xlim = c(0,0.5), ylim = ylim, xlab = "Frequency", ylab = "",
       xaxs = "i", yaxs = "i")
  rect(xleft = 0, xright  = 0.5, ybottom = 0,  ytop = 1, col="gray")
  
  # time series plot
  par(fig=c(0.20, 0.40, 0.9, 1), new = TRUE)
  plot(1:nrow(noiseList[[1]]), noiseList[[1]][,n], type = "l", col = "slateblue", lwd = lwd_ts, xlab = "Time (years)", ylab = "")
  
  # plot wavelet power spectrum
  par(fig=c(0.40, 0.6, 0.9, 1), new = TRUE)
  white.wt <- wt(cbind(1:nrow(noiseList[[1]]), noiseList[[1]][,n]), dj = 0.01, J1 = J1, max.scale = 32, mother = "morlet", sig.test = 0, sig.level = 0.95)
  plot(white.wt)
  
  # plot spawner abundance time series
  par(fig=c(0.6, 0.8, 0.9, 1), new = TRUE)
  plot(1:nrow(spawners), spawners[,white], type = "l", col = "slateblue", lwd = lwd_ts, xlab = "Time (years)", ylab = "")
  
  # plot wavelet power spectrum of spawner abundance
  par(fig=c(0.8, 1.0, 0.9, 1), new = TRUE)
  white.wt <- wt(cbind(1:nrow(spawners), spawners[,white]), dj = 0.01, J1 = J1, max.scale = 32, mother = "morlet", sig.test = 0, sig.level = 0.95)
  plot(white.wt)
  
  # Bandpass period 3-4
  # Generating spectrum
  par(fig=c(0, 0.2, 0.8, 0.9), new = TRUE)
  plot(x = seq(0,0.5, length = 50), y = runif(50, min = 0, max = 2), type = "n",
       xlim = c(0,0.5), ylim = ylim, xlab = "Frequency", ylab = "",
       xaxs = "i", yaxs = "i")
  rect(xleft = 0.25, xright  = 0.33, ybottom = 0,  ytop = 1, col="gray")
  
  # time series plot
  par(fig=c(0.2, 0.4, 0.8, 0.9), new = TRUE)
  plot(1:nrow(noiseList[[2]]), noiseList[[2]][,n], type = "l", col = "slateblue", lwd = lwd_ts, xlab = "Time (years)", ylab = "")
  
  # plot wavelet power spectrum
  par(fig=c(0.4, 0.6, 0.8, 0.9), new = TRUE)
  p34.wt <- wt(cbind(1:nrow(noiseList[[2]]), noiseList[[2]][,n]), dj = 0.01, J1 = J1, max.scale = 32, mother = "morlet", sig.test = 0, sig.level = 0.95)
  plot(p34.wt)
  
  # plot spawner abundance time series
  par(fig=c(0.6, 0.8, 0.8, 0.9), new = TRUE)
  plot(1:nrow(spawners), spawners[,p34], type = "l", col = "slateblue", lwd = lwd_ts, xlab = "Time (years)", ylab = "")
  
  # plot wavelet power spectrum of spawner abundance
  par(fig=c(0.8, 1.0, 0.8, 0.9), new = TRUE)
  p34.wt <- wt(cbind(1:nrow(spawners), spawners[,p34]), dj = 0.01, J1 = J1, max.scale = 32, mother = "morlet", sig.test = 0, sig.level = 0.95)
  plot(p34.wt)
  
  # Band reject period 3-4
  # Generating spectrum
  par(fig=c(0, 0.2, 0.7, 0.8), new = TRUE)
  plot(x = seq(0,0.5, length = 50), y = runif(50, min = 0, max = 2), type = "n",
       xlim = c(0,0.5), ylim = ylim, xlab = "Frequency", ylab = "",
       xaxs = "i", yaxs = "i")
  rect(xleft = 0, xright  = 0.25, ybottom = 0,  ytop = 1, col="gray")
  rect(xleft = 0.33, xright  = 0.5, ybottom = 0,  ytop = 1, col="gray")
  
  # time series plot
  par(fig=c(0.2, 0.4, 0.7, 0.8), new = TRUE)
  plot(1:nrow(noiseList[[3]]), noiseList[[3]][,n], type = "l", col = "slateblue", lwd = lwd_ts, xlab = "Time (years)", ylab = "")
  
  # plot wavelet power spectrum
  par(fig=c(0.4, 0.6, 0.7, 0.8), new = TRUE)
  pNo34.wt <- wt(cbind(1:nrow(noiseList[[3]]), noiseList[[3]][,n]), dj = 0.01, J1 = J1, max.scale = 32, mother = "morlet", sig.test = 0, sig.level = 0.95)
  plot(pNo34.wt)
  
  # plot spawner abundance time series
  par(fig=c(0.6, 0.8, 0.7, 0.8), new = TRUE)
  plot(1:nrow(spawners), spawners[,pNo34], type = "l", col = "slateblue", lwd = lwd_ts, xlab = "Time (years)", ylab = "")
  
  # plot wavelet power spectrum of spawner abundance
  par(fig=c(0.8, 1.0, 0.7, 0.8), new = TRUE)
  pNo34.wt <- wt(cbind(1:nrow(spawners), spawners[,pNo34]), dj = 0.01, J1 = J1, max.scale = 32, mother = "morlet", sig.test = 0, sig.level = 0.95)
  plot(pNo34.wt)
  
  
  # Bandpass greater than period 4
  # Generating spectrum
  par(fig=c(0, 0.2, 0.6, 0.7), new = TRUE)
  plot(x = seq(0,0.5, length = 50), y = runif(50, min = 0, max = 2), type = "n",
       xlim = c(0,0.5), ylim = ylim, xlab = "Frequency", ylab = "",
       xaxs = "i", yaxs = "i")
  rect(xleft = 0, xright  = 0.25, ybottom = 0,  ytop = 1, col="gray")
  
  # time series plot
  par(fig=c(0.2, 0.4, 0.6, 0.7), new = TRUE)
  plot(1:nrow(noiseList[[4]]), noiseList[[4]][,n], type = "l", col = "slateblue", lwd = lwd_ts, xlab = "Time (years)", ylab = "")
  
  # plot wavelet power spectrum
  par(fig=c(0.4, 0.6, 0.6, 0.7), new = TRUE)
  pgt4.wt <- wt(cbind(1:nrow(noiseList[[4]]), noiseList[[4]][,n]), dj = 0.01, J1 = J1, max.scale = 32, mother = "morlet", sig.test = 0, sig.level = 0.95)
  plot(pgt4.wt)
  
  # plot spawner abundance time series
  par(fig=c(0.6, 0.8, 0.6, 0.7), new = TRUE)
  plot(1:nrow(spawners), spawners[,pgt3], type = "l", col = "slateblue", lwd = lwd_ts, xlab = "Time (years)", ylab = "")
  
  # plot wavelet power spectrum of spawner abundance
  par(fig=c(0.8, 1.0, 0.6, 0.7), new = TRUE)
  pgt4.wt <- wt(cbind(1:nrow(spawners), spawners[,pgt3]), dj = 0.01, J1 = J1, max.scale = 32, mother = "morlet", sig.test = 0, sig.level = 0.95)
  plot(pgt4.wt)
  
  
  # Bandpass less than period 3
  # Generating spectrum
  par(fig=c(0, 0.2, 0.5, 0.6), new = TRUE)
  plot(x = seq(0,0.5, length = 50), y = runif(50, min = 0, max = 2), type = "n",
       xlim = c(0,0.5), ylim = ylim, xlab = "Frequency", ylab = "",
       xaxs = "i", yaxs = "i")
  rect(xleft = 0.33, xright  = 0.5, ybottom = 0,  ytop = 1, col="gray")
  
  # time series plot
  par(fig=c(0.2, 0.4, 0.5, 0.6), new = TRUE)
  plot(1:nrow(noiseList[[5]]), noiseList[[5]][,n], type = "l", col = "slateblue", lwd = lwd_ts, xlab = "Time (years)", ylab = "")
  
  # plot wavelet power spectrum
  par(fig=c(0.4, 0.6, 0.5, 0.6), new = TRUE)
  plt4.wt <- wt(cbind(1:nrow(noiseList[[5]]), noiseList[[5]][,n]), dj = 0.01, J1 = J1, max.scale = 32, mother = "morlet", sig.test = 0, sig.level = 0.95)
  plot(plt4.wt)
  
  # plot spawner abundance time series
  par(fig=c(0.6, 0.8, 0.5, 0.6), new = TRUE)
  plot(1:nrow(spawners), spawners[,plt4], type = "l", col = "slateblue", lwd = lwd_ts, xlab = "Time (years)", ylab = "")
  
  # plot wavelet power spectrum of spawner abundance
  par(fig=c(0.8, 1.0, 0.5, 0.6), new = TRUE)
  plt4.wt <- wt(cbind(1:nrow(spawners), spawners[,plt4]), dj = 0.01, J1 = J1, max.scale = 32, mother = "morlet", sig.test = 0, sig.level = 0.95)
  plot(plt4.wt)
  
  # Bandpass greater than period 3
  # Generating spectrum
  par(fig=c(0, 0.2, 0.4, 0.5), new = TRUE)
  plot(x = seq(0,0.5, length = 50), y = runif(50, min = 0, max = 2), type = "n",
       xlim = c(0,0.5), ylim = ylim, xlab = "Frequency", ylab = "",
       xaxs = "i", yaxs = "i")
  rect(xleft = 0, xright  = 0.33, ybottom = 0,  ytop = 1, col="gray")
  
  # time series plot
  par(fig=c(0.2, 0.4, 0.4, 0.5), new = TRUE)
  plot(1:nrow(noiseList[[6]]), noiseList[[6]][,n], type = "l", col = "slateblue", lwd = lwd_ts, xlab = "Time (years)", ylab = "")
  
  # plot wavelet power spectrum
  par(fig=c(0.4, 0.6, 0.4, 0.5), new = TRUE)
  pgt4.wt <- wt(cbind(1:nrow(noiseList[[6]]), noiseList[[6]][,n]), dj = 0.01, J1 = J1, max.scale = 32, mother = "morlet", sig.test = 0, sig.level = 0.95)
  plot(pgt4.wt)
  
  # plot spawner abundance time series
  par(fig=c(0.6, 0.8, 0.4, 0.5), new = TRUE)
  plot(1:nrow(spawners), spawners[,pgt4], type = "l", col = "slateblue", lwd = lwd_ts, xlab = "Time (years)", ylab = "")
  
  # plot wavelet power spectrum of spawner abundance
  par(fig=c(0.8, 1.0, 0.4, 0.5), new = TRUE)
  plgt4.wt <- wt(cbind(1:nrow(spawners), spawners[,pgt4]), dj = 0.01, J1 = J1, max.scale = 32, mother = "morlet", sig.test = 0, sig.level = 0.95)
  plot(pgt4.wt)
  
  
  # Bandpass less than period 4
  # Generating spectrum
  par(fig=c(0, 0.2, 0.3, 0.4), new = TRUE)
  plot(x = seq(0,0.5, length = 50), y = runif(50, min = 0, max = 2), type = "n",
       xlim = c(0,0.5), ylim = ylim, xlab = "Frequency", ylab = "",
       xaxs = "i", yaxs = "i")
  rect(xleft = 0.25, xright  = 0.5, ybottom = 0,  ytop = 1, col="gray")
  
  # time series plot
  par(fig=c(0.2, 0.4, 0.3, 0.4), new = TRUE)
  plot(1:nrow(noiseList[[7]]), noiseList[[7]][,n], type = "l", col = "slateblue", lwd = lwd_ts, xlab = "Time (years)", ylab = "")
  
  # plot wavelet power spectrum
  par(fig=c(0.4, 0.6, 0.3, 0.4), new = TRUE)
  plt3.wt <- wt(cbind(1:nrow(noiseList[[7]]), noiseList[[7]][,n]), dj = 0.01, J1 = J1, max.scale = 32, mother = "morlet", sig.test = 0, sig.level = 0.95)
  plot(plt3.wt)
  
  # plot spawner abundance time series
  par(fig=c(0.6, 0.8, 0.3, 0.4), new = TRUE)
  plot(1:nrow(spawners), spawners[,plt3], type = "l", col = "slateblue", lwd = lwd_ts, xlab = "Time (years)", ylab = "")
  
  # plot wavelet power spectrum of spawner abundance
  par(fig=c(0.8, 1.0, 0.3, 0.4), new = TRUE)
  plt3.wt <- wt(cbind(1:nrow(spawners), spawners[,plt3]), dj = 0.01, J1 = J1, max.scale = 32, mother = "morlet", sig.test = 0, sig.level = 0.95)
  plot(plt3.wt)
  
  
  # Bandpass greater than period 10
  # Generating spectrum
  par(fig=c(0, 0.2, 0.2, 0.3), new = TRUE)
  plot(x = seq(0,0.5, length = 50), y = runif(50, min = 0, max = 2), type = "n",
       xlim = c(0,0.5), ylim = ylim, xlab = "Frequency", ylab = "",
       xaxs = "i", yaxs = "i")
  rect(xleft = 0, xright  = 0.1, ybottom = 0,  ytop = 1, col="gray")
  
  # time series plot
  par(fig=c(0.2, 0.4, 0.2, 0.3), new = TRUE)
  plot(1:nrow(noiseList[[8]]), noiseList[[8]][,n], type = "l", col = "slateblue", lwd = lwd_ts, xlab = "Time (years)", ylab = "")
  
  # plot wavelet power spectrum
  par(fig=c(0.4, 0.6, 0.2, 0.3), new = TRUE)
  pgt10.wt <- wt(cbind(1:nrow(noiseList[[8]]), noiseList[[8]][,n]), dj = 0.01, J1 = J1, max.scale = 32, mother = "morlet", sig.test = 0, sig.level = 0.95)
  plot(pgt10.wt)
  
  # plot spawner abundance time series
  par(fig=c(0.6, 0.8, 0.2, 0.3), new = TRUE)
  plot(1:nrow(spawners), spawners[,pgt10], type = "l", col = "slateblue", lwd = lwd_ts, xlab = "Time (years)", ylab = "")
  
  # plot wavelet power spectrum of spawner abundance
  par(fig=c(0.8, 1.0, 0.2, 0.3), new = TRUE)
  pgt10.wt <- wt(cbind(1:nrow(spawners), spawners[,pgt10]), dj = 0.01, J1 = J1, max.scale = 32, mother = "morlet", sig.test = 0, sig.level = 0.95)
  plot(pgt10.wt)
  
  # Bandpass greater than period 10 and period 3-4
  # Generating spectrum
  par(fig=c(0, 0.2, 0.1, 0.2), new = TRUE)
  plot(x = seq(0,0.5, length = 50), y = runif(50, min = 0, max = 2), type = "n",
       xlim = c(0,0.5), ylim = ylim, xlab = "Frequency", ylab = "",
       xaxs = "i", yaxs = "i")
  rect(xleft = 0, xright  = 0.1, ybottom = 0,  ytop = 1, col="gray")
  rect(xleft = 0.25, xright  = 0.33, ybottom = 0,  ytop = 1, col="gray")
  
  # time series plot
  par(fig=c(0.2, 0.4, 0.1, 0.2), new = TRUE)
  plot(1:nrow(noiseList[[9]]), noiseList[[9]][,n], type = "l", col = "slateblue", lwd = lwd_ts, xlab = "Time (years)", ylab = "")
  
  # plot wavelet power spectrum
  par(fig=c(0.4, 0.6, 0.1, 0.2), new = TRUE)
  p34gt10.wt <- wt(cbind(1:nrow(noiseList[[9]]), noiseList[[9]][,n]), dj = 0.01, J1 = J1, max.scale = 32, mother = "morlet", sig.test = 0, sig.level = 0.95)
  plot(p34gt10.wt)
  
  # plot spawner abundance time series
  par(fig=c(0.6, 0.8, 0.1, 0.2), new = TRUE)
  plot(1:nrow(spawners), spawners[,p34gt10], type = "l", col = "slateblue", lwd = lwd_ts, xlab = "Time (years)", ylab = "")
  
  # plot wavelet power spectrum of spawner abundance
  par(fig=c(0.8, 1.0, 0.1, 0.2), new = TRUE)
  p34gt10.wt <- wt(cbind(1:nrow(spawners), spawners[,p34gt10]), dj = 0.01, J1 = J1, max.scale = 32, mother = "morlet", sig.test = 0, sig.level = 0.95)
  plot(p34gt10.wt)
  
  par(old)
  
} # end summaryCh3tsPlot

summaryCh3tsPlotNoise <- function(noise = noiseList, n = 1, J1 = trunc((log(32/(2 * 1))/log(2))/0.01)) {
  
  library(biwavelet)
  
  ylim = c(0,2)
  lwd_ts <- 1.5
  
  old <- par(mar = c(1,2,1,1), cex = .7)
  
  # white noise
  # Generating spectrum
  par(fig=c(0, 0.33, 0.9, 1))
  plot(x = seq(0,0.5, length = 50), y = runif(50, min = 0, max = 2), type = "n",
       xlim = c(0,0.5), ylim = ylim, xlab = "Frequency", ylab = "",
       xaxs = "i", yaxs = "i")
  rect(xleft = 0, xright  = 0.5, ybottom = 0,  ytop = 1, col="gray")
  
  # time series plot
  par(fig=c(0.33, 0.66, 0.9, 1), new = TRUE)
  plot(1:nrow(noiseList[[1]]), noiseList[[1]][,n], type = "l", col = "slateblue", lwd = lwd_ts, xlab = "Time (years)", ylab = "")
  
  # plot wavelet power spectrum
  par(fig=c(0.66, 1, 0.9, 1), new = TRUE)
  white.wt <- wt(cbind(1:nrow(noiseList[[1]]), noiseList[[1]][,n]), dj = 0.01, J1 = J1, max.scale = 32, mother = "morlet", sig.test = 0, sig.level = 0.95)
  plot(white.wt)
  
  # Bandpass period 3-4
  # Generating spectrum
  par(fig=c(0, 0.33, 0.8, 0.9), new = TRUE)
  plot(x = seq(0,0.5, length = 50), y = runif(50, min = 0, max = 2), type = "n",
       xlim = c(0,0.5), ylim = ylim, xlab = "Frequency", ylab = "",
       xaxs = "i", yaxs = "i")
  rect(xleft = 0.25, xright  = 0.33, ybottom = 0,  ytop = 1, col="gray")
  
  # time series plot
  par(fig=c(0.33, 0.66, 0.8, 0.9), new = TRUE)
  plot(1:nrow(noiseList[[2]]), noiseList[[2]][,n], type = "l", col = "slateblue", lwd = lwd_ts, xlab = "Time (years)", ylab = "")
  
  # plot wavelet power spectrum
  par(fig=c(0.66, 1, 0.8, 0.9), new = TRUE)
  p34.wt <- wt(cbind(1:nrow(noiseList[[2]]), noiseList[[2]][,n]), dj = 0.01, J1 = J1, max.scale = 32, mother = "morlet", sig.test = 0, sig.level = 0.95)
  plot(p34.wt)
  
  
  # Band reject period 3-4
  # Generating spectrum
  par(fig=c(0, 0.33, 0.7, 0.8), new = TRUE)
  plot(x = seq(0,0.5, length = 50), y = runif(50, min = 0, max = 2), type = "n",
       xlim = c(0,0.5), ylim = ylim, xlab = "Frequency", ylab = "",
       xaxs = "i", yaxs = "i")
  rect(xleft = 0, xright  = 0.25, ybottom = 0,  ytop = 1, col="gray")
  rect(xleft = 0.33, xright  = 0.5, ybottom = 0,  ytop = 1, col="gray")
  
  # time series plot
  par(fig=c(0.33, 0.66, 0.7, 0.8), new = TRUE)
  plot(1:nrow(noiseList[[3]]), noiseList[[3]][,n], type = "l", col = "slateblue", lwd = lwd_ts, xlab = "Time (years)", ylab = "")
  
  # plot wavelet power spectrum
  par(fig=c(0.66, 1, 0.7, 0.8), new = TRUE)
  pNo34.wt <- wt(cbind(1:nrow(noiseList[[3]]), noiseList[[3]][,n]), dj = 0.01, J1 = J1, max.scale = 32, mother = "morlet", sig.test = 0, sig.level = 0.95)
  plot(pNo34.wt)
  
  
  # Bandpass greater than period 4
  # Generating spectrum
  par(fig=c(0, 0.33, 0.6, 0.7), new = TRUE)
  plot(x = seq(0,0.5, length = 50), y = runif(50, min = 0, max = 2), type = "n",
       xlim = c(0,0.5), ylim = ylim, xlab = "Frequency", ylab = "",
       xaxs = "i", yaxs = "i")
  rect(xleft = 0, xright  = 0.25, ybottom = 0,  ytop = 1, col="gray")
  
  # time series plot
  par(fig=c(0.33, 0.66, 0.6, 0.7), new = TRUE)
  plot(1:nrow(noiseList[[4]]), noiseList[[4]][,n], type = "l", col = "slateblue", lwd = lwd_ts, xlab = "Time (years)", ylab = "")
  
  # plot wavelet power spectrum
  par(fig=c(0.66, 1, 0.6, 0.7), new = TRUE)
  pgt4.wt <- wt(cbind(1:nrow(noiseList[[4]]), noiseList[[4]][,n]), dj = 0.01, J1 = J1, max.scale = 32, mother = "morlet", sig.test = 0, sig.level = 0.95)
  plot(pgt4.wt)
  
  
  # Bandpass less than period 3
  # Generating spectrum
  par(fig=c(0, 0.33, 0.5, 0.6), new = TRUE)
  plot(x = seq(0,0.5, length = 50), y = runif(50, min = 0, max = 2), type = "n",
       xlim = c(0,0.5), ylim = ylim, xlab = "Frequency", ylab = "",
       xaxs = "i", yaxs = "i")
  rect(xleft = 0.33, xright  = 0.5, ybottom = 0,  ytop = 1, col="gray")
  
  # time series plot
  par(fig=c(0.33, 0.66, 0.5, 0.6), new = TRUE)
  plot(1:nrow(noiseList[[5]]), noiseList[[5]][,n], type = "l", col = "slateblue", lwd = lwd_ts, xlab = "Time (years)", ylab = "")
  
  # plot wavelet power spectrum
  par(fig=c(0.66, 1, 0.5, 0.6), new = TRUE)
  plt4.wt <- wt(cbind(1:nrow(noiseList[[5]]), noiseList[[5]][,n]), dj = 0.01, J1 = J1, max.scale = 32, mother = "morlet", sig.test = 0, sig.level = 0.95)
  plot(plt4.wt)
  
  
  # Bandpass greater than period 3
  # Generating spectrum
  par(fig=c(0, 0.33, 0.4, 0.5), new = TRUE)
  plot(x = seq(0,0.5, length = 50), y = runif(50, min = 0, max = 2), type = "n",
       xlim = c(0,0.5), ylim = ylim, xlab = "Frequency", ylab = "",
       xaxs = "i", yaxs = "i")
  rect(xleft = 0, xright  = 0.33, ybottom = 0,  ytop = 1, col="gray")
  
  # time series plot
  par(fig=c(0.33, 0.66, 0.4, 0.5), new = TRUE)
  plot(1:nrow(noiseList[[6]]), noiseList[[6]][,n], type = "l", col = "slateblue", lwd = lwd_ts, xlab = "Time (years)", ylab = "")
  
  # plot wavelet power spectrum
  par(fig=c(0.66, 1, 0.4, 0.5), new = TRUE)
  pgt4.wt <- wt(cbind(1:nrow(noiseList[[6]]), noiseList[[6]][,n]), dj = 0.01, J1 = J1, max.scale = 32, mother = "morlet", sig.test = 0, sig.level = 0.95)
  plot(pgt4.wt)
  
  
  # Bandpass less than period 4
  # Generating spectrum
  par(fig=c(0, 0.33, 0.3, 0.4), new = TRUE)
  plot(x = seq(0,0.5, length = 50), y = runif(50, min = 0, max = 2), type = "n",
       xlim = c(0,0.5), ylim = ylim, xlab = "Frequency", ylab = "",
       xaxs = "i", yaxs = "i")
  rect(xleft = 0.25, xright  = 0.5, ybottom = 0,  ytop = 1, col="gray")
  
  # time series plot
  par(fig=c(0.33, 0.66, 0.3, 0.4), new = TRUE)
  plot(1:nrow(noiseList[[7]]), noiseList[[7]][,n], type = "l", col = "slateblue", lwd = lwd_ts, xlab = "Time (years)", ylab = "")
  
  # plot wavelet power spectrum
  par(fig=c(0.66, 1, 0.3, 0.4), new = TRUE)
  plt3.wt <- wt(cbind(1:nrow(noiseList[[7]]), noiseList[[7]][,n]), dj = 0.01, J1 = J1, max.scale = 32, mother = "morlet", sig.test = 0, sig.level = 0.95)
  plot(plt3.wt)
  
  
  # Bandpass greater than period 10
  # Generating spectrum
  par(fig=c(0, 0.33, 0.2, 0.3), new = TRUE)
  plot(x = seq(0,0.5, length = 50), y = runif(50, min = 0, max = 2), type = "n",
       xlim = c(0,0.5), ylim = ylim, xlab = "Frequency", ylab = "",
       xaxs = "i", yaxs = "i")
  rect(xleft = 0, xright  = 0.1, ybottom = 0,  ytop = 1, col="gray")
  
  # time series plot
  par(fig=c(0.33, 0.66, 0.2, 0.3), new = TRUE)
  plot(1:nrow(noiseList[[8]]), noiseList[[8]][,n], type = "l", col = "slateblue", lwd = lwd_ts, xlab = "Time (years)", ylab = "")
  
  # plot wavelet power spectrum
  par(fig=c(0.66, 1, 0.2, 0.3), new = TRUE)
  pgt10.wt <- wt(cbind(1:nrow(noiseList[[8]]), noiseList[[8]][,n]), dj = 0.01, J1 = J1, max.scale = 32, mother = "morlet", sig.test = 0, sig.level = 0.95)
  plot(pgt10.wt)
  
  
  # Bandpass greater than period 10 and period 3-4
  # Generating spectrum
  par(fig=c(0, 0.33, 0.1, 0.2), new = TRUE)
  plot(x = seq(0,0.5, length = 50), y = runif(50, min = 0, max = 2), type = "n",
       xlim = c(0,0.5), ylim = ylim, xlab = "Frequency", ylab = "",
       xaxs = "i", yaxs = "i")
  rect(xleft = 0, xright  = 0.1, ybottom = 0,  ytop = 1, col="gray")
  rect(xleft = 0.25, xright  = 0.33, ybottom = 0,  ytop = 1, col="gray")
  
  # time series plot
  par(fig=c(0.33, 0.66, 0.1, 0.2), new = TRUE)
  plot(1:nrow(noiseList[[9]]), noiseList[[9]][,n], type = "l", col = "slateblue", lwd = lwd_ts, xlab = "Time (years)", ylab = "")
  
  # plot wavelet power spectrum
  par(fig=c(0.66, 1, 0.1, 0.2), new = TRUE)
  p34gt10.wt <- wt(cbind(1:nrow(noiseList[[9]]), noiseList[[9]][,n]), dj = 0.01, J1 = J1, max.scale = 32, mother = "morlet", sig.test = 0, sig.level = 0.95)
  plot(p34gt10.wt)
  
  par(old)
  
} # end summaryCh3tsPlotNoise


summaryCh3tsPlotSpawners <- function(spawners = storage, 
                                     meanSurv = 0.35,
                                     sigma = 0.2, 
                                     n = 1, 
                                     J1 = trunc((log(32/(2 * 1))/log(2))/0.01)) {
  
  spPlot <- copy(spawners[ i = N > 400 & reps_c == 1 & sigPSmult_c == sigma & meanPS_c == meanSurv ])
  
  library(biwavelet)
  
  ylim = c(0,2)
  lwd_ts <- 1.5
  
  old <- par(mar = c(1,2,1,1), cex = .7)
  
  # white noise
  # Generating spectrum
  par(fig=c(0, 0.33, 0.9, 1))
  plot(x = seq(0,0.5, length = 50), y = runif(50, min = 0, max = 2), type = "n",
       xlim = c(0,0.5), ylim = ylim, xlab = "Frequency", ylab = "",
       xaxs = "i", yaxs = "i")
  rect(xleft = 0, xright  = 0.5, ybottom = 0,  ytop = 1, col="gray")
  
  # time series plot
  par(fig=c(0.33, 0.66, 0.9, 1), new = TRUE)
  plot(1:length(spPlot[ j = white]), spPlot[ j = white], type = "l", col = "slateblue", lwd = lwd_ts, xlab = "Time (years)", ylab = "")
  
  # plot wavelet power spectrum
  par(fig=c(0.66, 1, 0.9, 1), new = TRUE)
  white.wt <- wt(cbind(1:length(spPlot[ j = white]), spPlot[ j = white]), dj = 0.01, J1 = J1, max.scale = 32, mother = "morlet", sig.test = 0, sig.level = 0.95)
  plot(white.wt)
  
  # Bandpass period 3-4
  # Generating spectrum
  par(fig=c(0, 0.33, 0.8, 0.9), new = TRUE)
  plot(x = seq(0,0.5, length = 50), y = runif(50, min = 0, max = 2), type = "n",
       xlim = c(0,0.5), ylim = ylim, xlab = "Frequency", ylab = "",
       xaxs = "i", yaxs = "i")
  rect(xleft = 0.25, xright  = 0.33, ybottom = 0,  ytop = 1, col="gray")
  
  # time series plot
  par(fig=c(0.33, 0.66, 0.8, 0.9), new = TRUE)
  plot(1:length(spPlot[ j = p34]), spPlot[ j = p34], type = "l", col = "slateblue", lwd = lwd_ts, xlab = "Time (years)", ylab = "")
  
  # plot wavelet power spectrum
  par(fig=c(0.66, 1, 0.8, 0.9), new = TRUE)
  p34.wt <- wt(cbind(1:length(spPlot[ j = p34]), spPlot[ j = p34]), dj = 0.01, J1 = J1, max.scale = 32, mother = "morlet", sig.test = 0, sig.level = 0.95)
  plot(p34.wt)
  
  
  # Band reject period 3-4
  # Generating spectrum
  par(fig=c(0, 0.33, 0.7, 0.8), new = TRUE)
  plot(x = seq(0,0.5, length = 50), y = runif(50, min = 0, max = 2), type = "n",
       xlim = c(0,0.5), ylim = ylim, xlab = "Frequency", ylab = "",
       xaxs = "i", yaxs = "i")
  rect(xleft = 0, xright  = 0.25, ybottom = 0,  ytop = 1, col="gray")
  rect(xleft = 0.33, xright  = 0.5, ybottom = 0,  ytop = 1, col="gray")
  
  # time series plot
  par(fig=c(0.33, 0.66, 0.7, 0.8), new = TRUE)
  plot(1:length(spPlot[ j = pNo34]), spPlot[ j = pNo34], type = "l", col = "slateblue", lwd = lwd_ts, xlab = "Time (years)", ylab = "")
  
  # plot wavelet power spectrum
  par(fig=c(0.66, 1, 0.7, 0.8), new = TRUE)
  pNo34.wt <- wt(cbind(1:length(spPlot[ j = pNo34]), spPlot[ j = pNo34]), dj = 0.01, J1 = J1, max.scale = 32, mother = "morlet", sig.test = 0, sig.level = 0.95)
  plot(pNo34.wt)
  
  
  # Bandpass greater than period 4
  # Generating spectrum
  par(fig=c(0, 0.33, 0.6, 0.7), new = TRUE)
  plot(x = seq(0,0.5, length = 50), y = runif(50, min = 0, max = 2), type = "n",
       xlim = c(0,0.5), ylim = ylim, xlab = "Frequency", ylab = "",
       xaxs = "i", yaxs = "i")
  rect(xleft = 0, xright  = 0.25, ybottom = 0,  ytop = 1, col="gray")
  
  # time series plot
  par(fig=c(0.33, 0.66, 0.6, 0.7), new = TRUE)
  plot(1:length(spPlot[ j = pgt4]), spPlot[ j = pgt4], type = "l", col = "slateblue", lwd = lwd_ts, xlab = "Time (years)", ylab = "")
  
  # plot wavelet power spectrum
  par(fig=c(0.66, 1, 0.6, 0.7), new = TRUE)
  pgt4.wt <- wt(cbind(1:length(spPlot[ j = pgt4]), spPlot[ j = pgt4]), dj = 0.01, J1 = J1, max.scale = 32, mother = "morlet", sig.test = 0, sig.level = 0.95)
  plot(pgt4.wt)
  
  
  # Bandpass less than period 3
  # Generating spectrum
  par(fig=c(0, 0.33, 0.5, 0.6), new = TRUE)
  plot(x = seq(0,0.5, length = 50), y = runif(50, min = 0, max = 2), type = "n",
       xlim = c(0,0.5), ylim = ylim, xlab = "Frequency", ylab = "",
       xaxs = "i", yaxs = "i")
  rect(xleft = 0.33, xright  = 0.5, ybottom = 0,  ytop = 1, col="gray")
  
  # time series plot
  par(fig=c(0.33, 0.66, 0.5, 0.6), new = TRUE)
  plot(1:length(spPlot[ j = plt3]), spPlot[ j = plt3], type = "l", col = "slateblue", lwd = lwd_ts, xlab = "Time (years)", ylab = "")
  
  # plot wavelet power spectrum
  par(fig=c(0.66, 1, 0.5, 0.6), new = TRUE)
  plt4.wt <- wt(cbind(1:length(spPlot[ j = plt3]), spPlot[ j = plt3]), dj = 0.01, J1 = J1, max.scale = 32, mother = "morlet", sig.test = 0, sig.level = 0.95)
  plot(plt4.wt)
  
  
  # Bandpass greater than period 3
  # Generating spectrum
  par(fig=c(0, 0.33, 0.4, 0.5), new = TRUE)
  plot(x = seq(0,0.5, length = 50), y = runif(50, min = 0, max = 2), type = "n",
       xlim = c(0,0.5), ylim = ylim, xlab = "Frequency", ylab = "",
       xaxs = "i", yaxs = "i")
  rect(xleft = 0, xright  = 0.33, ybottom = 0,  ytop = 1, col="gray")
  
  # time series plot
  par(fig=c(0.33, 0.66, 0.4, 0.5), new = TRUE)
  plot(1:length(spPlot[ j = pgt3]), spPlot[ j = pgt3], type = "l", col = "slateblue", lwd = lwd_ts, xlab = "Time (years)", ylab = "")
  
  # plot wavelet power spectrum
  par(fig=c(0.66, 1, 0.4, 0.5), new = TRUE)
  pgt4.wt <- wt(cbind(1:length(spPlot[ j = pgt3]), spPlot[ j = pgt3]), dj = 0.01, J1 = J1, max.scale = 32, mother = "morlet", sig.test = 0, sig.level = 0.95)
  plot(pgt4.wt)
  
  
  # Bandpass less than period 4
  # Generating spectrum
  par(fig=c(0, 0.33, 0.3, 0.4), new = TRUE)
  plot(x = seq(0,0.5, length = 50), y = runif(50, min = 0, max = 2), type = "n",
       xlim = c(0,0.5), ylim = ylim, xlab = "Frequency", ylab = "",
       xaxs = "i", yaxs = "i")
  rect(xleft = 0.25, xright  = 0.5, ybottom = 0,  ytop = 1, col="gray")
  
  # time series plot
  par(fig=c(0.33, 0.66, 0.3, 0.4), new = TRUE)
  plot(1:length(spPlot[ j = plt4]), spPlot[ j = plt4], type = "l", col = "slateblue", lwd = lwd_ts, xlab = "Time (years)", ylab = "")
  
  # plot wavelet power spectrum
  par(fig=c(0.66, 1, 0.3, 0.4), new = TRUE)
  plt3.wt <- wt(cbind(1:length(spPlot[ j = plt4]), spPlot[ j = plt4]), dj = 0.01, J1 = J1, max.scale = 32, mother = "morlet", sig.test = 0, sig.level = 0.95)
  plot(plt3.wt)
  
  
  # Bandpass greater than period 10
  # Generating spectrum
  par(fig=c(0, 0.33, 0.2, 0.3), new = TRUE)
  plot(x = seq(0,0.5, length = 50), y = runif(50, min = 0, max = 2), type = "n",
       xlim = c(0,0.5), ylim = ylim, xlab = "Frequency", ylab = "",
       xaxs = "i", yaxs = "i")
  rect(xleft = 0, xright  = 0.1, ybottom = 0,  ytop = 1, col="gray")
  
  # time series plot
  par(fig=c(0.33, 0.66, 0.2, 0.3), new = TRUE)
  plot(1:length(spPlot[ j = pgt10]), spPlot[ j = pgt10], type = "l", col = "slateblue", lwd = lwd_ts, xlab = "Time (years)", ylab = "")
  
  # plot wavelet power spectrum
  par(fig=c(0.66, 1, 0.2, 0.3), new = TRUE)
  pgt10.wt <- wt(cbind(1:length(spPlot[ j = pgt10]), spPlot[ j = pgt10]), dj = 0.01, J1 = J1, max.scale = 32, mother = "morlet", sig.test = 0, sig.level = 0.95)
  plot(pgt10.wt)
  
  
  # Bandpass greater than period 10 and period 3-4
  # Generating spectrum
  par(fig=c(0, 0.33, 0.1, 0.2), new = TRUE)
  plot(x = seq(0,0.5, length = 50), y = runif(50, min = 0, max = 2), type = "n",
       xlim = c(0,0.5), ylim = ylim, xlab = "Frequency", ylab = "",
       xaxs = "i", yaxs = "i")
  rect(xleft = 0, xright  = 0.1, ybottom = 0,  ytop = 1, col="gray")
  rect(xleft = 0.25, xright  = 0.33, ybottom = 0,  ytop = 1, col="gray")
  
  # time series plot
  par(fig=c(0.33, 0.66, 0.1, 0.2), new = TRUE)
  plot(1:length(spPlot[ j = p34gt10]), spPlot[ j = p34gt10], type = "l", col = "slateblue", lwd = lwd_ts, xlab = "Time (years)", ylab = "")
  
  # plot wavelet power spectrum
  par(fig=c(0.66, 1, 0.1, 0.2), new = TRUE)
  p34gt10.wt <- wt(cbind(1:length(spPlot[ j = p34gt10]), spPlot[ j = p34gt10]), dj = 0.01, J1 = J1, max.scale = 32, mother = "morlet", sig.test = 0, sig.level = 0.95)
  plot(p34gt10.wt)
  
  par(old)
  
} # end summaryCh3tsPlotSpawners

corrCoefDist <- function(white, wave34, rsw34, hi  = 1, lo = -(hi)) {
  
  w <- cor(white)
  wtCorrs <- w[col(w) < row(w)]
  print(range(wtCorrs))
  wt <- cor(wave34)
  wv34Corrs <- wt[col(wt) < row(wt)]
  rsw <- cor(rsw34)
  rsw34Corrs <- rsw[col(rsw) < row(rsw)]
  
  brk <- seq(lo, hi, by = 0.05)
  
  empCIwt  <- quantile(wtCorrs, c(0.05, 0.95))
  empCIwv34  <- quantile(wv34Corrs, c(0.05, 0.95))
  empCIrsw34  <- quantile(rsw34Corrs, c(0.05, 0.95))
  
  hist(wtCorrs, breaks = brk, col = rgb(0,0,0,.8), xlim = c(lo, hi), main = "", 
       xlab = expression(paste(rho)))
  hist(wv34Corrs, breaks = brk, col = rgb(1,0,0,.6), add = TRUE)
  hist(rsw34Corrs, breaks = brk, col = rgb(0,0,1,.4), add = TRUE)
  
  abline(v = empCIwt[1], col = rgb(0,0,0,1), lty = 2, lwd = 2)
  abline(v = empCIwt[2], col = rgb(0,0,0,1), lty = 2, lwd = 2)
  abline(v = empCIwv34[1], col = rgb(1,0,0,1), lty = 2, lwd = 2)
  abline(v = empCIwv34[2], col = rgb(1,0,0,1), lty = 2, lwd = 2)
  abline(v = empCIrsw34[1], col = rgb(0,0,1,1), lty = 2, lwd = 2)
  abline(v = empCIrsw34[2], col = rgb(0,0,1,1), lty = 2, lwd = 2)
  legend("topright", c("white", "wavelet34", "rsw34"), col = c("black", "red", "blue"), lty = 1, lwd = 3, cex = .8)
}


# age-structured models that allows time-varying ocean survival 
# vary: alpha, EQsp and delta_e
# 
# variability about [mean: low, mid, high, very high: var: high and low for mean]
# 
# 





# # simple deterministic Leslie matrix ####
# baseDetLM <- function (mgmtScenarios = "BAU",
#                        climateScenarios = "B1",
#                        n0 = c(10000, 5000, 5000, 5000),
#                        f3 = 5530/2,
#                        f4 = 5530/2,
#                        s_eo = 0.02,
#                        s2 = 0.5,
#                        s3 = 0.8,
#                        d3 = 0.4) {
#   setwd("~/NEP_salmon/chap_3/data/simulation_results/")
#   temp <- 
#     # Use mean value from all climate models over first five years
#     # Mean prespawn survival rate (spawners/total number of returning females)   
#     FWps <- mean(as.matrix(read.delim(paste( mgmtScenarios, "/SALMODspawnFemSurv", climateScenarios,".txt", sep = ""), sep = ","  )[1:5,]), na.rm = TRUE) 
#   # Mean egg to fry survival
#   FWef <- mean(as.matrix(read.delim(paste( mgmtScenarios, "/", climateScenarios, "egg2fry.txt", sep = ""), sep = ","  )[1:5,]), na.rm = TRUE) 
#   # Mean survival rate of fry exiting system
#   FWf <- mean(as.matrix(read.delim(paste( mgmtScenarios, "/", climateScenarios, "fryout.txt", sep = ""), sep = ","  )[1:5,]), na.rm = TRUE) 
#   
#   A <-  matrix(c(0, 0, FWps*f3*FWef*FWf*d3, FWps*f3*FWef*FWf,
#                  s_eo, 0, 0, 0,
#                  0, s2, 0, 0,
#                  0, 0, s3*(1-d3), 0), 
#                nrow = 4,
#                byrow = TRUE)
#   Aeigs <- eigen(A)
#   nyr <- 89
#   Nage <- matrix(NA, nrow = nyr, ncol = length(n0))
#   
#   for (i in 1:nyr) {
#     if (i == 1) Nage[i,] <- A%*%n0 else Nage[i,] <- A%*%Nage[i-1,]
#   }
#   
#   pdf(paste("~/NEP_salmon/chap_3/output/", "/prelim_matrix_model/",mgmtScenarios, "_", climateScenarios, "_detMatrixSim.pdf", sep = ""), 
#       width=16/2.54, 
#       height=10/2.54, 
#       pointsize=10)
#   par(mgp=c(2.2,0.45,0), tcl=-0.4, mar=c(3.3,3.6,1.1,1.1))
#   matplot(2010:2099,rbind(n0,Nage)[,2:4], type = "l", col = 1:3,
#           xlab = "Years", ylab = "Numbers")
#   ditle(paste("s_eo=", s_eo, ", d3=", d3, ", fec=", f3,", EV1=", round(as.real(Aeigs$values[1]), 2),sep = ""), cex.main = 0.7)
#   legend( c("topleft"), c("Age-2", "Age-3", "Age-4"), col = 1:3, lty = 1, cex = .5)
#   dev.off()
# }


# # Age structured model analysis ####
# 
# # Leslie Matrix
# L <- matrix(c(0,0.576,1.152, 1.152,
#               0.8,0,0,0,
#               0,0.8,0,0,
#               0,0,0.8,0),
#             nrow = 4,
#             byrow = T)
# s0 <- 0.576
# 
# # do lifetable calculations for leslie matrix Chap x from Quinn and Deriso 1999
# lesMatCalcs <- function(L, s0) {
#   print("Values in first row are products of fecundity at age and age-0 survival (s0).
#         You should enter s0 as an arguement in the function for the function to make
#         correct calculations")
#   if (!is.matrix(L))  stop("must use a matrix for calculations")
#   a <- 1:ncol(L)
#   f <- L[1,]/s0
#   s <- diag(L[-1,-ncol(L)])
#   l <- numeric(nrow(L))
#   for (i in 1:nrow(L)) {
#     if (i == 1) l[i] <- 1 else l[i] <- prod(s[1:(i-1)])/l[1]
#   }
#   
#   L_eigs <- eigen(L)
#   
#   Evals <- L_eigs$values
#   lam1 <- as.real(L_eigs$values[1])
#   Rev1 <- as.real(L_eigs$vectors[,1]/L_eigs$vectors[1,1])
#   Lev1 <- as.real(solve(L_eigs$vectors)[1,])/as.real(solve(L_eigs$vectors)[1,1])
#   
#   faladivlam <- fala <- numeric(nrow(L))
#   for (i in 1:nrow(L)) {
#     faladivlam[i] <- f[i]*l[i]/(lam1^i)  
#   }
#   
#   for (i in 1:nrow(L)) {
#     fala[i] <- Reduce(`+`, f[i:length(f)]*l[i:length(l)]) # my first use of R's higher order funcitons!
#     
#   }
#   
#   R_age <- numeric(nrow(L)+1)
#   for (i in 1:length(R_age)) {
#     if (i == 1) {
#       R_age[i] <- s0*Reduce(`+`, f[i:length(f)]*l[i:length(l)])
#     } else  {
#       R_age[i] <- Reduce(`+`, f[(i-1):length(f)]*l[(i-1):length(l)])/l[i-1]
#     }
#   }
#   R0 <- R_age[1]
#   
#   afaRev1 <- a*f*Rev1
#   faRev1 <- f*Rev1
#   
# }
# 
# # Fitting stock recruitment models ####
# 
# linBHfit2 <- function(data) {
#   # Fits a lnearized for me of Beverton-Holt
#   # remove missing data
#   # Quinn and Deriso 1999, chapter 3 "Stock and Recruitment"
#   
#   # Two formats of BH
#   # 1. rec = alpha*stock/(1 + beta*stock)
#   # 2. rec = stock/(alpha_star + beta_star*stock)
#   
#   # This function uses form 2.
#   # relationships of parameters
#   # alpha_star = 1/alpha (inverse of slope of SR at origin)
#   # beta_star = beta/alpha
#   # Asympotic recruitment (maximum recruitment)
#   # rec_max = 1/beta_star = alpha/beta  
#   
#   # Only use complete cases with S & R data (no missing values)
#   data <- data[complete.cases(data),] 
#   data$inRec <- 1/(data$rec)
#   data$inStock <- 1/(data$stock)
#   
#   # fit the linear model  
#   fit <- lm(inRec ~ inStock, data = data) #
#   alpha_star <- coef(fit)[2] 
#   beta_star <- coef(fit)[1]
#   
#   return(list(fit = fit, alpha_star = alpha_star, beta_star = beta_star, data = data))
# }  
# 
# plotBH2 <- function(data, lmfit) {
#   # plot BH curve 
#   plot(data$stock, data$rec, pch = 19,
#        ylab = "Recruits",
#        xlab = "Stock",
#   )
#   xrange <- seq(0,max(data$stock), length.out = 100)
#   lines(xrange, 
#         xrange/(lmfit$alpha_star + lmfit$beta_star*xrange),
#         col = "slateblue",
#         lwd = 2)
#   #lines(data$stock, ((1/alpha_star)*data$stock)/ (1 + (beta_star/alpha_star)*data$stock) , col = "blue", lty = 2)  
#   
#   # output slope of SR at origin
#   text(x = 0.7*max(xrange),
#        y = 0.2*max(data$rec),
#        labels = paste("Slope SR at origin = ", round(1/(lmfit$alpha_star), digits = 3)),
#        cex = 0.8)
#   # output max recruitment
#   text(x = 0.7*max(xrange),
#        y = 0.1*max(data$rec),
#        labels = paste("Max recruit = ", round(1/(lmfit$beta_star), digits = 1)),
#        cex = 0.8)
# }
# 
# linBHfitdiag <- function(data, lmfit) {
#   
#   # plot fit of linearized BH model version 2 in Q&D
#   # Check: does the line fit the data?
#   plot(data$inStock, data$inRec)
#   lines(data$inStock, lmfit$beta_star + lmfit$alpha_star*data$inStock)
#   
#   # plot.lm() diagnostic plots
#   op <- par(mfrow = c(2,2))
#   plot(fit) 
#   par(op) # see notes below
#   
#   # Plot leverages
#   fit_mat <- model.matrix(fit)
#   lev <- hat(fit_mat)
#   p <- 2 # num params [slope and intercept]
#   # sum(lev)
#   plot(lev,ylab="Leverages",main="Index plot of Leverages", ylim = c(0,.5))
#   abline(h = 2*p/length(lev))
#   
#   
#   plot(density(residuals(fit),na.rm=TRUE))
#   plot(sort(residuals(fit)),pch= 19, col = "slateblue")
#   
#   # notes on the diagnostic plot.lm() output
#   # http://www.stat.berkeley.edu/classes/s133/Lr0.html
#   # 1. resid v fitted: randomly distributed above/below evenly, var constant over time
#   # 2. Normal QQ: departures from straight line potential violations of normality assumption
#   # 3. sqrt(std resid) v. fitted (scale-location): no dicernable pattern; sim to res v fitted
#   # 4. cook's distance: identify points which have more influence than other points 
#   # (points that are distant from other points in the data, either for the 
#   # dependent variable or one or more independent variables)
#   # large values might require additional investigation
#   # 5. resid v leverage: Labeled points on this plot represent cases we may want 
#   # to investigate as possibly having undue influence on the regression 
#   # relationship
#   # 6. cook's distance v. leverage:  (wikipedia)
#   # Cook's distance measures the effect of deleting a given observation.
#   # a commonly used estimate of the influence of a data point
#   # to indicate data points that are particularly worth checking for validity; 
#   # to indicate regions of the design space where it would be good to be able 
#   # to obtain more data points.
#   # versus
#   # those observations, if any, made at extreme or outlying values of the 
#   # independent variables such that the lack of neighboring observations means 
#   # that the fitted regression model will pass close to that particular observation.
#   
#   
# }
# 
# nlsBH_FSA <- function(data) {
#   # nonlinear BH from FSA package
#   library(FSA)
#   library(nlstools)
#   data$logR <-log(data$rec)
#   
#   bh1s <- srStarts(rec ~ stock, data = data, type = "BevertonHolt", param = 1)
#   
#   bh1 <- logR ~ log((a*stock)/(1+b*stock))
#   
#   bh1nls <- nls(bh1, data = data, start = list(a = 1, b = 0.0006), algorithm = "port", lower = c(0,0)) 
#   
#   return(list(start = bh1s, model = bh1, fittedmodel = bh1nls))
#   
# }  
# 
# 
# 
# plotBHnls <- function(data, nlsList) {
#   
#   plot(rec ~ stock, data = data, 
#        pch = 19,
#        col = "slateblue",
#        xlab = "Stock",
#        ylab = "Recruits 3 years later",
#        main = "Snorkel survey data")
#   curve((coef(nlsList$fittedmodel)[1]*x)/(1+coef(nlsList$fittedmodel)[2]*x), 
#         from = 0, 
#         to = 20000,
#         col = "red4",
#         lwd = 3,
#         add = TRUE)
# }
# 
# # diagnostic plots of nls BH fit
# diagBHnls <- function(nlsList, nBoot = 1000) {
#   
#   print(summary(nlsList$fittedmodel))
#   
#   bootbh <- nlsBoot(nlsList$fittedmodel, niter = nBoot)
#   plot(1:length(bootbh$coefboot[,2]), 
#        bootbh$coefboot[,1]/bootbh$coefboot[,2],
#        type = "l",
#        ylab = "Max Recruits (bootstrapped a/b)",
#        xlab = "",
#        main = "Bootstrapped Snorkel Survey\nRecruits 3 years later\na/b")
#   maxR <- bootbh$coefboot[,1]/bootbh$coefboot[,2]
#   lines(1:length(bootbh$coefboot[,2]),
#         rep(mean(maxR),
#             times = length(bootbh$coefboot[,2])),
#         col = "red",
#         lwd = 3)
#   lines(1:length(bootbh$coefboot[,2]),
#         rep(median(maxR),
#             times = length(bootbh$coefboot[,2])),
#         col = "blue",
#         lwd = 3)
#   plot(bootbh)
#   confint(bootbh, plot = TRUE)
#   print(confint(bootbh, plot = FALSE))
#   par(mfrow = c(1,1)) # confint divides the plotting window, so 'undo' that
#   
# }
# 
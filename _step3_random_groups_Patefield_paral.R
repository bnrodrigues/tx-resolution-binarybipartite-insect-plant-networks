###############################################################################
# load work space do "several-webs_metrics_win"
#setwd("<your_work_diretory_path>")

webType <- "antagonism" #"mutualism"

# carrega a work space da etapa 1 (step1)
padrao <- paste0("ws_", webType)
arq <- list.files("txtoutput/", pattern = padrao)
arq <- arq[length(arq)] # arquivo por ultimo (ordem alfabetica)
cat( paste0(arq,"\n") )
load(paste0("txtoutput/",arq))

# limpa as var 
rm(list = setdiff(ls(), c("listMtMdOrigSemProblem","lstGps", 
                          "arqs","gps","webType")))

padrao <- paste0("charact_mtMetrics_", webType)
arq <- list.files("../output/", pattern = padrao)
arq <- arq[length(arq)] # arquivo por ultimo (ordem alfabetica)
info <- read.csv( paste0("../output/", arq) )

numMt <- sum(info$GroupType == "Original")

vrsn <- "_v01"

subconj <- 1:47 # 1:47 (47/2 = 24)
# subconj <- 25:47
###############################################################################
# bibliotecas calculo das metricas
#install.packages("bipartite")
library(bipartite)
# biblioteca para a paralelizacao
library(parallel)

### var relacionadas a paralelizacao
# fala qts cores fisicos existem
#P <- detectCores(logical = F)
# deixa 1/4 dos cores para outros usuarios
#P <- as.integer(P-ceiling(P/4) )
P <- as.integer(7)

# fala para construir as unidades de processamento (cl)
# com base na quantidade de cores fisicos (P)
cl <- makeCluster(P)
# cl eh uma lista, cada elemento é um nó. É uma lista
# com dados da conexao e rank

# para cada core entender as bibliotecas, precisa falar para 
# ele entender a biblioteca
clusterEvalQ(cl, library(bipartite))


# var necessarias
# listMtMdOrigSemProblem
# arqs
namesInfo <- c('Tipo' , 'GroupType', 'ArqName', 'RandNum', 'numAnimals', 
               'numPlants', 'numSpecies', 'numInterac')
metric <- c("connectance","links per species","nestedness", "NODF",
            "number of compartments","web asymmetry")
metricL <- c("robustness","niche overlap")

# quantidade de repeticoes do processo e agrupamento aleaotrio
qtRandMt <- 100 

###############################################################################
###############################################################################
# funcao
geraMtRand <- function(z,temp,arq,gp) {
    # temp <- listMtMdOrigSemProblem[[j]] # mto gd. Ocupa mta RAM. Pssr so 1
  
    #--------------------------------------------------------------------------
    # matriz aleatoria (Patefield's 1981)
    temp <- (r2dtable(1, 
                   rowSums(temp), # ttl marginal linhas
                   colSums(temp) )[[1]] # ttl marginal col
                   >0)*1 # tornar binaria
    #--------------------------------------------------------------------------
    
    # salva a matriz aleatoria
    write.table(temp, paste0("../output/mt_simulation/", webType, vrsn,
                                "_infoRandPatefield_",
                                "_", 
                                as.character(gp),"_",
                                # as.character(gps$ID[w]),"_",
                                sub("....$","",arq),"-",
                                # sub("....$","",arqs[j]),"-",
                                formatC(z, width=3, flag="0"),".txt"),
                sep = "\t", row.names = F, col.names = F)

    ###########################################################
    # calculo das metricas
    # list para armazenar info das matrizes aleatorias
    rwRandMtInfo <- vector("list", length= length(namesInfo))
    names(rwRandMtInfo) <- namesInfo
    # armazena calculos metricas das matrizes aleatorias
    rwRandLstMetrics <- vector("list", length= length(metric)) 
    names(rwRandLstMetrics) <- metric
    rwRandLstMetricsL <- vector("list", length= length(metricL))
	    
    rwRandMtInfo$Tipo <- webType
    rwRandMtInfo$RandNum <- z
    rwRandMtInfo$ArqName <- arq
    rwRandMtInfo$GroupType <- as.character(gp)
    rwRandMtInfo$numPlants <- length(temp[,1])
    rwRandMtInfo$numAnimals <- length(temp[1,])
    rwRandMtInfo$numSpecies <- rwRandMtInfo$numPlants[length(rwRandMtInfo$numPlants)] + 
		               rwRandMtInfo$numAnimals[length(rwRandMtInfo$numAnimals)]
    rwRandMtInfo$numInterac <- sum(temp)
    for(i in 1:length(metric)){ # calculo das metricas para cd mt aleatoria
      rwRandLstMetrics[c(metric[i])] <- networklevel(temp, index = metric[i])
    }
    for(i in 1:length(metricL)){ # calculo das metricas para cd mt aleatoria
      rwRandLstMetricsL[[i]] <- networklevel(temp, index = metricL[i])
    }
    # fim calculo das metricas
    ###########################################################
	    # lista com uma lista de listas...
	    return(list("rwRandMtInfo" = rwRandMtInfo, 
		        "rwRandLstMetrics" = rwRandLstMetrics,
				    "rwRandLstMetricsL" = rwRandLstMetricsL))
}
# Fim da funcao
##############################################################################################################################################################
##############################################################################################################################################################


# inicializa algumas var
wichGroup <- 1
contNumRandMt <- 0 # quatidade de redes aleatorias feitas (qd a mt agrupada eh menor que a original)
contNumRandNA <- 0

###########################################################
# inicializa o arquivo que armazena as informacoes
# list para armazenar info das matrizes aleatorias
randMtInfo <- vector("list", length= length(namesInfo))
names(randMtInfo) <- namesInfo
# armazena calculos metricas das matrizes aleatorias
RandLstMetrics <- vector("list", length= length(metric)) 
names(RandLstMetrics) <- metric
RandLstMetricsL <- vector("list", length= length(metricL))
###########################################################
randMtInfo$Tipo <- append(randMtInfo$Tipo, NA)
randMtInfo$RandNum <- append(randMtInfo$RandNum, NA)
randMtInfo$ArqName <- append(randMtInfo$ArqName, NA)
randMtInfo$GroupType <- append(randMtInfo$GroupType,NA)
randMtInfo$numPlants <- append(randMtInfo$numPlants,NA)
randMtInfo$numAnimals <- append(randMtInfo$numAnimals,NA)
randMtInfo$numSpecies <- append(randMtInfo$numSpecies,NA)
randMtInfo$numInterac <- append(randMtInfo$numInterac,NA)
for(i in 1:length(metric)){ # calculo das metricas para cd mt aleatoria
  RandLstMetrics[[i]] <- append(RandLstMetrics[[i]], NA)
}
for(i in 1:length(metricL)){ # calculo das metricas para cd mt aleatoria
  RandLstMetricsL[[i]] <- rbind(RandLstMetricsL[[i]], c(NA,NA))
}
RandLstMetricsL <- as.data.frame(RandLstMetricsL)
colnames(RandLstMetricsL) <- paste( rep(metricL, each= 2), 
                                    rep(c(".HL",".LL")), sep="" )
infoRand <- cbind(as.data.frame(randMtInfo),
                  as.data.frame(RandLstMetrics),
                  RandLstMetricsL)
write.table(infoRand, paste0("../output/", webType, vrsn,
                             "_infoRandPatefield_", 
                             sub(":","-",paste(data.frame(subconj))),
                             ".csv"), sep = ",")
                             
                             
                             
###########################################################
## passa as var do ambiente para a func executada nos nos
# var mto grande, passar so a parte necessaria como 
# listMtMdOrigSemProblem[[j]]
clusterExport(cl,c("webType","metric", "metricL","namesInfo","gps","vrsn"))





###########################################################

## w ## tipo de agrupamento: 4 tipos
# sp x sp (1), sp x ge (2), ge x sp (5) e ge x ge (6)
for(w in c(1,2,5,6)){
  cat(paste0("\n",as.character(gps$ID[w]))) # nome do agrupamento
  ## j ## qual matriz vai ser simulada
  for(j in subconj){
    
    ###########################################################
    # list para armazenar info das matrizes aleatorias
    randMtInfo <- vector("list", length= length(namesInfo))
    names(randMtInfo) <- namesInfo
    # armazena calculos metricas das matrizes aleatorias
    RandLstMetrics <- vector("list", length= length(metric)) 
    names(RandLstMetrics) <- metric
    RandLstMetricsL <- vector("list", length= length(metricL))
    ###########################################################
    
    # exibe no "log"
    cat(paste("\n\t",as.character(info$ArqName[j]))) # nome do arquivo
    
    temp <- lstGps[[w]][[j]]

    # se tem diferencas entre a rede original e a agrupada
    # loops1 <- info$numAnimals[j]-info$numAnimals[j+numMt*w] # colunas
    # loops2 <- info$numPlants[j]-info$numPlants[j+numMt*w] # linhas
    # if(loops1+loops2 > 0){
      # quantidade de interacoes
      tam <- sum(temp)
      ## z ou qtRandMt ## quantidade de redes aleatorias
        # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        # gera as redes e calcula as metricas dentro de um func q eh 
        # paralelizada. Depois só concatena os resultados das mt aleatorias
        # funcao recebe temp,arqs[j],gps$ID[w]
        tempRand <- clusterApply(cl, 1:qtRandMt, fun = geraMtRand, 
                               temp = temp,
                               arq = arqs[j],
                               gp = gps$ID[w])

        # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++	
        contNumRandMt <- contNumRandMt + qtRandMt  
        # concatena os resultados das mt aleatorias
      for(z in 1:qtRandMt){
        randMtInfo$Tipo <- append(randMtInfo$Tipo, tempRand[[z]]$rwRandMtInfo$Tipo)
        randMtInfo$RandNum <- append(randMtInfo$RandNum, tempRand[[z]]$rwRandMtInfo$RandNum)
        randMtInfo$ArqName <- append(randMtInfo$ArqName, tempRand[[z]]$rwRandMtInfo$ArqName)
        randMtInfo$GroupType <- append(randMtInfo$GroupType, tempRand[[z]]$rwRandMtInfo$GroupType)
        randMtInfo$numPlants <- append(randMtInfo$numPlants, tempRand[[z]]$rwRandMtInfo$numPlants)
        randMtInfo$numAnimals <- append(randMtInfo$numAnimals, tempRand[[z]]$rwRandMtInfo$numAnimals)
        randMtInfo$numSpecies <- append(randMtInfo$numSpecies,
                                        randMtInfo$numPlants[length(randMtInfo$numPlants)] + 
                                          randMtInfo$numAnimals[length(randMtInfo$numAnimals)])
        randMtInfo$numInterac <- append(randMtInfo$numInterac, tempRand[[z]]$rwRandMtInfo$numInterac)
        for(i in 1:length(metric)){ # calculo das metricas para cd mt aleatoria
          RandLstMetrics[[i]] <- append(RandLstMetrics[[i]], tempRand[[z]]$rwRandLstMetrics[[i]])
        }
        for(i in 1:length(metricL)){ # calculo das metricas para cd mt aleatoria
          RandLstMetricsL[[i]] <- rbind(RandLstMetricsL[[i]], tempRand[[z]]$rwRandLstMetricsL[[i]])
        }        
      }
	  
      # se nao tem diferencas entre a rede original e a agrupada
    # } else {
    #   contNumRandNA <- contNumRandNA + 1
    #   cat(" - ",contNumRandNA," no change")
    #   
    #   ###########################################################
    #   # matriz nao calcula as metricas, mas registra NA
    #   randMtInfo$Tipo <- append(randMtInfo$Tipo, webType)
    #   randMtInfo$RandNum <- append(randMtInfo$RandNum, NA)
    #   randMtInfo$ArqName <- append(randMtInfo$ArqName, arqs[j])
    #   randMtInfo$GroupType <- append(randMtInfo$GroupType,as.character(gps$ID[w]))
    #   randMtInfo$numPlants <- append(randMtInfo$numPlants,NA)
    #   randMtInfo$numAnimals <- append(randMtInfo$numAnimals,NA)
    #   randMtInfo$numSpecies <- append(randMtInfo$numSpecies,NA)
    #   randMtInfo$numInterac <- append(randMtInfo$numInterac,NA)
    #   for(i in 1:length(metric)){ # calculo das metricas para cd mt aleatoria
    #     RandLstMetrics[[i]] <- append(RandLstMetrics[[i]], NA)
    #   }
    #   for(i in 1:length(metricL)){ # calculo das metricas para cd mt aleatoria
    #     RandLstMetricsL[[i]] <- rbind(RandLstMetricsL[[i]], c(NA,NA))
    #   }
    #   # fim "matriz nao calcula as metricas, mas registra NA"
    #   ###########################################################
    # }
    ###########################################################
    # grava as info da matriz
    infoRand <- cbind(as.data.frame(randMtInfo),
                      as.data.frame(RandLstMetrics),
                      as.data.frame(RandLstMetricsL))
    write.table(infoRand, paste0("../output/", webType, vrsn,
                                 "_infoRandPatefield_", 
                                 sub(":","-",paste(data.frame(subconj))),
                                 ".csv"), sep = ",", col.names = F, append = T)
    ###########################################################
  }
  wichGroup <- wichGroup + 1
}

# parar de usar CPUs separados
stopCluster(cl) # limpar o que ja tem

numMtRand <- (wichGroup-1)*numMt*qtRandMt
cat(numMtRand == contNumRandMt + qtRandMt*contNumRandNA)

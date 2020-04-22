#####################################################################
# grafico de deltas de conectancia, NODF e modularidade             #
#####################################################################
#setwd("<your_work_diretory_path>")


#####################################################################
## funcao calculo dos deltas
dltCalc <- function(mtor, mttx, mtpr, mtpt, nMetric){
  temp <- data.frame(grp = as.character(mttx$GroupType), delta = vector("numeric", length = nrow(mttx)))
  for(i in 1:nrow(mttx)){
    temp$delta[i] <- mtor[mtor$ArqName == as.character(mttx$ArqName[i]), nMetric] - 
      mttx[i, nMetric]
  }
  dltMean <- tapply(na.exclude(temp$delta), temp$grp[!is.na(temp$delta)], mean)
  dltSD <- tapply(na.exclude(temp$delta), temp$grp[!is.na(temp$delta)], sd)
  # par a par
  temp <- data.frame(grp = as.character(mtpr$GroupType), delta = vector("numeric", length = nrow(mtpr)))
  for(i in 1:nrow(mtpr)){
    temp$delta[i] <- mtor[mtor$ArqName == as.character(mtpr$ArqName[i]), nMetric] - 
      mtpr[i, nMetric]
  }
  dltMean <- cbind.data.frame(dltMean, tapply(na.exclude(temp$delta), temp$grp[!is.na(temp$delta)], mean))
  dltSD <- cbind.data.frame(dltSD, tapply(na.exclude(temp$delta), temp$grp[!is.na(temp$delta)], sd))
  colnames(dltMean)[1] <- "tax"
  colnames(dltSD)[1] <- "tax"
  colnames(dltMean)[2] <- "par"
  colnames(dltSD)[2] <- "par"
  # Patefield
  temp <- data.frame(grp = as.character(mtpt$GroupType), delta = vector("numeric", length = nrow(mtpt)))
  for(i in 1:nrow(mtpt)){
    temp$delta[i] <- mtor[mtor$ArqName == as.character(mtpt$ArqName[i]), nMetric] - 
      mtpt[i, nMetric]
  }
  dltMean <- cbind.data.frame(dltMean, tapply(na.exclude(temp$delta), temp$grp[!is.na(temp$delta)], mean))
  dltSD <- cbind.data.frame(dltSD, tapply(na.exclude(temp$delta), temp$grp[!is.na(temp$delta)], sd))
  colnames(dltMean)[3] <- "pat"
  colnames(dltSD)[3] <- "pat"
  
  
  rownames(dltMean) <- unique(as.character(mttx$GroupType))
  rownames(dltMean) <- sub("A_", "", rownames(dltMean))
  rownames(dltMean) <- sub("\\.x\\.", "\nX\n", rownames(dltMean))
  rownames(dltMean) <- sub("P_", "", rownames(dltMean))
  rownames(dltMean) <- gsub("sp", "species", rownames(dltMean))
  rownames(dltMean) <- gsub("ge", "genus", rownames(dltMean))
  rownames(dltMean) <- gsub("fa", "family", rownames(dltMean))
  rownames(dltMean) <- gsub("or", "order", rownames(dltMean))
  
  rownames(dltSD) <- rownames(dltMean)
  
  # inverte a ordem das linhas
  dltMean <- dltMean[nrow(dltMean):1,]
  dltSD <- dltSD[nrow(dltSD):1,]
  
  # # "A_sp.x.P_sp" "A_sp.x.P_ge" "A_ge.x.P_sp" "A_ge.x.P_ge"
  # rownames(dltMean) <- c("species\nX\nspecies","specie\nX\ngenus",
  #                        "genus\nX\nspecie","genus\nX\ngenus")
  
  dlt <- list(Mean = t(dltMean), SD = t(dltSD))
  return(dlt)
}

#####################################################################
## funcao calculo dos deltas
dltPrctCalc <- function(mtor, mttx, nMetric){
  temp <- data.frame(grp = as.character(mttx$GroupType), 
                     delta = vector("numeric", length = nrow(mttx)))
  for(i in 1:nrow(mttx)){
    temp$delta[i] <- ((mtor[mtor$ArqName == 
                              as.character(mttx$ArqName[i]), 
                            nMetric] - mttx[i, 
                                            nMetric])/mtor[mtor$ArqName == 
                                                  as.character(mttx$ArqName[i]), 
                                                  nMetric])*100
  }
  dltMean <- tapply(na.exclude(temp$delta), temp$grp[!is.na(temp$delta)], mean)
  dltSD <- tapply(na.exclude(temp$delta), temp$grp[!is.na(temp$delta)], sd)
  
  #colnames(dltMean)[1] <- "tax"
  #colnames(dltSD)[1] <- "tax"
  
  names(dltMean) <- unique(as.character(mttx$GroupType))
  names(dltMean) <- sub("A_", "", rownames(dltMean))
  names(dltMean) <- sub("\\.x\\.", "\nX\n", rownames(dltMean))
  names(dltMean) <- sub("P_", "", rownames(dltMean))
  names(dltMean) <- gsub("sp", "species", rownames(dltMean))
  names(dltMean) <- gsub("ge", "genus", rownames(dltMean))
  names(dltMean) <- gsub("fa", "family", rownames(dltMean))
  names(dltMean) <- gsub("or", "order", rownames(dltMean))
  
  names(dltSD) <- names(dltMean)
  
  # inverte a ordem das linhas
  dltMean <- dltMean[length(dltMean):1]
  dltSD <- dltSD[length(dltSD):1]
  
  # # "A_sp.x.P_sp" "A_sp.x.P_ge" "A_ge.x.P_sp" "A_ge.x.P_ge"
  # rownames(dltMean) <- c("specie\nX\nspecie","specie\nX\ngenus",
  #                        "genus\nX\nspecie","genus\nX\ngenus")
  
  dlt <- list(Mean = t(dltMean), SD = t(dltSD))
  return(dlt)
}


## Funcao para grafico do delta
pltGrfDlt <- function(dltLst, metricN){
  # limites superior e inforior do grafico
  mini <- floor(min(dltLst$Mean - dltLst$SD)*100)/100
  maxi <- ceiling(max(dltLst$Mean + dltLst$SD)*100)/100
  # mar=c(5.1, 4.1, 4.1, 2.1)
  par(mar=c(2.5, 4.1, 1.5, 2.1))
  # barCenters recebe as posicoes das colunas
  barCenters <- barplot(as.matrix(dltLst$Mean), beside = T,
                        col = c("#434343","#707070","#c8c8c8"),
                        ylim= c(mini,maxi), border = NA, xaxt="n")
  title(ylab=paste0("Δ of ", metricN))
  axis(1,at=barCenters[2,], labels=colnames(dltLst$Mean), tick=F ,pos=mini-((maxi-mini)/6))
  segments(barCenters, dltLst$Mean - dltLst$SD, barCenters,
           dltLst$Mean + dltLst$SD, lwd = 1.0)
  
  arrows(barCenters, dltLst$Mean - dltLst$SD, barCenters,
         dltLst$Mean + dltLst$SD, lwd = 1.5, angle = 90,
         code = 3, length = 0.05)
  # legend(1,2*(mini/3), # posiciona no 2/3 do lim inferior
  #        c("Taxonomic","Paired random","Patefield's algorithm"), 
  #        fill = c("#434343","#707070","#c8c8c8"), 
  #        bty = "n", border = "white")
  # points(barCenters, dltLst$Mean, pch= 18) # pto na marcando a media
}

## Funcao para grafico do delta so tx
pltDltTx <- function(dltLst, metricN){
  # limites superior e inforior do grafico
  mini <- floor(min(dltLst$Mean[1,] - dltLst$SD[1,])*100)/100
  maxi <- ceiling(max(dltLst$Mean[1,] + dltLst$SD[1,])*100)/100
  # mar=c(5.1, 4.1, 4.1, 2.1)
  par(mar=c(2.5, 4.1, 1.5, 2.1))
  # barCenters recebe as posicoes das colunas
  barCenters <- barplot(as.matrix(dltLst$Mean[1,]), beside = T,
                        col = c("#434343","#707070","#c8c8c8"),
                        ylim= c(mini,maxi), border = NA, xaxt="n")
  title(ylab=paste0("Δ of ", metricN))
  axis(1,at=barCenters, labels=colnames(dltLst$Mean), tick=F ,pos=mini-((maxi-mini)/6))
  segments(barCenters, dltLst$Mean[1,] - dltLst$SD[1,], barCenters,
           dltLst$Mean[1,] + dltLst$SD[1,], lwd = 1.0)
  
  arrows(barCenters, dltLst$Mean[1,] - dltLst$SD[1,], barCenters,
         dltLst$Mean[1,] + dltLst$SD[1,], lwd = 1.5, angle = 90,
         code = 3, length = 0.05)
  # legend(1,2*(mini/3), # posiciona no 2/3 do lim inferior
  #        c("Taxonomic","Paired random","Patefield's algorithm"), 
  #        fill = c("#434343","#707070","#c8c8c8"), 
  #        bty = "n", border = "white")
  # points(barCenters, dltLst$Mean, pch= 18) # pto na marcando a media
}

# Fim das funcoes
#####################################################################



#####################################################################
# Leitura das informacoes das redes
## leitura dos dados de redes mutualisticas
mutu_tax <- read.csv("../output/mutualism_infoOriTxSpGe+modul_2019-10-d01.csv")[,-1] # so tem 
mutu_par <- read.csv("../output/mutualism_v01_infoRandParPar+modul_1-47_2019-10-d01.csv", row.names = NULL)[,-1]
mutu_pat <- read.csv("../output/mutualism_v02_infoRandPatefield+modul_1-47_2019-10-d14.csv", row.names = NULL)[,-1]

# separa original das agrupadas por taxonomia
mutu_ori <- mutu_tax[mutu_tax$GroupType == "Original",] # so info das originais
mutu_tax <- mutu_tax[mutu_tax$GroupType != "Original",] # so dos agrup taxo

# tira a coluna "RandNum" para as colunas ficarem equivalentes
#mutu_pat <- mutu_pat[,-3]
mutu_par <- mutu_par[,-3]


## leitura dos dados de redes antagonicas
anta_tax <- read.csv("../output/antagonism_infoOriTxSpGe+modul_2019-10-d01.csv")[,-1] # so tem 
anta_par <- read.csv("../output/antagonism_v01_infoRandParPar+modul_1-47_2019-10-d01.csv", row.names = NULL)[,-1]
anta_pat <- read.csv("../output/antagonism_v02_infoRandPatefield+modul_1-47_2019-10-d14.csv", row.names = NULL)[,-1]
# remove as redes do HOSTS
anta_tax <- anta_tax[-grep("herb-hosts",anta_tax[,1]), ]
anta_par <- anta_par[-grep("herb-hosts",anta_par[,1]), ]
anta_pat <- anta_pat[-grep("herb-hosts",anta_pat[,1]), ]

# separa original das agrupadas por taxonomia
anta_ori <- anta_tax[anta_tax$GroupType == "Original",] # so info das originais
anta_tax <- anta_tax[anta_tax$GroupType != "Original",] # so dos agrup taxo

# tira a coluna "RandNum" para as colunas ficarem equivalentes
#anta_pat <- anta_pat[,-3]
anta_par <- anta_par[,-3]


#####################################################################
# calculo dos deltas 
# uso da func dltCalc. Entrada da func: mtor, mttx, mtpr, mtpt, nMetric

## Mutualistica
# conectance
dltMCon <- dltCalc(mutu_ori, mutu_tax, mutu_par, mutu_pat, 8)
# NDOF
dltMNODF <- dltCalc(mutu_ori, mutu_tax, mutu_par, mutu_pat, 11)
# Modularidade
dltMMod <- dltCalc(mutu_ori, mutu_tax, mutu_par, mutu_pat, 19)

## Antagonica
# conectance
dltACon <- dltCalc(anta_ori, anta_tax, anta_par, anta_pat, 8)
# NDOF
dltANODF <- dltCalc(anta_ori, anta_tax, anta_par, anta_pat, 11)
# Modularidade
dltAMod <- dltCalc(anta_ori, anta_tax, anta_par, anta_pat, 19)


# save.image(paste0("txtoutput/grafDelta-ws_","v01_",
#                   format(Sys.time(), "%Y-%m-d%d"),
#                   ".RData"))

arq <- list.files("txtoutput/", pattern = "grafDelta-ws_")
arq <- arq[length(arq)] # arquivo por ultimo (ordem alfabetica)
load(paste0("txtoutput/",arq))
#####################################################################
# Grafico dos deltas
# uso da funcao pltGrfDlt. Entradas da func: dlt, metricN

# para uma figura grande de width = 174 e height = 174
svg(filename = paste0("../output/mutu-anta_delta_con-NODF-mod_", format(Sys.time(), "%Y-%m-d%d"),".svg"),
    width = 7, height = 7, # tamanho equivalente a 177x177 mm 
    pointsize = 14) # 18 equivale a uma fonte de 12 de Arial
par(mfcol= c(4,2), oma = c(0, 0, 2, 0))
pltGrfDlt(dltMCon, "Connectance")
pltGrfDlt(dltMNODF, "NODF")
pltGrfDlt(dltMMod, "Modularity")

plot.new() # plot vazio, sem info

pltGrfDlt(dltACon, "Connectance")
pltGrfDlt(dltANODF, "NODF")
pltGrfDlt(dltAMod, "Modularity")

plot.new()

legend("topright",
       c("Taxonomic","Paired random","Patefield's algorithm"),
       fill = c("#434343","#707070","#c8c8c8"),
       bty = "n", border = "white")
mtext("Mutualism   Antagonism", outer = TRUE, cex = 0.7)

dev.off()

par(mfcol= c(1,1))

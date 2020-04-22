###############################################################################
# Tabela com numero de interacoes                                             #
###############################################################################

#setwd("<your_work_diretory_path>")



###############################################################################
# MUTUALISM
webType <- "mutualism"

# carrega a work space da etapa 1 (step1)
padrao <- paste0("ws_", webType)
arq <- list.files("txtoutput/", pattern = padrao)
arq <- arq[length(arq)] # arquivo por ultimo (ordem alfabetica)
cat( paste0(arq,"\n") )
load(paste0("txtoutput/",arq))
rm(list = setdiff(ls(),
                  # deixa so as var de interesse
                  c("webType","arqs","gps","listMtOrig","listMtMdOrigSemProblem","lstGps")) )

numNos <- cbind.data.frame(arqs, 
                           Source = sapply(listMtOrig, nrow) + sapply(listMtOrig, ncol),
                           Final = sapply(listMtMdOrigSemProblem, nrow) + sapply(listMtMdOrigSemProblem, ncol))

for(i in 1:16) {
  numNos <- cbind.data.frame( numNos, sapply(lstGps[[i]], nrow) + sapply(lstGps[[i]], ncol) )
  colnames(numNos)[ncol(numNos)] <- as.character(gps[i,1])
}

perNos1 <- numNos
# qt tirou
for(i in 2:ncol(perNos1)) {
  perNos1[,i] <- ((numNos[,2]-perNos1[,i])/numNos[,2])*100
}

perNos2 <- numNos
# qt manteve
for(i in 2:ncol(perNos2)) {
  perNos2[,i] <- (perNos2[,i]/numNos[,2])*100
}


###############################################################################
# ordena colunas com base na abrangencia
# por combinacao
# ord <- c("A_sp.x.P_sp", 
#          "A_sp.x.P_ge", "A_ge.x.P_sp", "A_ge.x.P_ge", 
#          "A_sp.x.P_fa", "A_fa.x.P_sp", "A_ge.x.P_fa", "A_fa.x.P_ge", "A_fa.x.P_fa", 
#          "A_sp.x.P_or", "A_or.x.P_sp", "A_ge.x.P_or", "A_or.x.P_ge", "A_fa.x.P_or", "A_or.x.P_fa", "A_or.x.P_or")

# por valores de abrangencia
# ord <- c("A_sp.x.P_sp", "A_ge.x.P_sp", "A_sp.x.P_ge", 
#              "A_fa.x.P_sp", "A_sp.x.P_fa", 
#          "A_ge.x.P_ge", "A_fa.x.P_ge", "A_ge.x.P_fa", 
#              "A_or.x.P_sp", "A_sp.x.P_or", 
#              "A_or.x.P_ge", "A_ge.x.P_or", 
#          "A_fa.x.P_fa", "A_or.x.P_fa", "A_fa.x.P_or", 
#          "A_or.x.P_or")

ord <- c("A_sp.x.P_sp", "A_sp.x.P_ge", "A_sp.x.P_fa", "A_sp.x.P_or",
         "A_ge.x.P_sp", "A_ge.x.P_ge", "A_ge.x.P_fa", "A_ge.x.P_or",
         "A_fa.x.P_sp", "A_fa.x.P_ge", "A_fa.x.P_fa", "A_fa.x.P_or",
         "A_or.x.P_sp", "A_or.x.P_ge", "A_or.x.P_fa", "A_or.x.P_or")

numNos <- numNos[, c("arqs","Source", "Final", ord)]
perNos1 <- perNos1[, c("arqs","Source", "Final", ord)]
perNos2 <- perNos2[, c("arqs","Source", "Final", ord)]


###############################################################################
# mean and sd
numNos.m <- list(Mean = apply(numNos[,-1], 2, mean), SD = apply(numNos[,-1], 2, sd) )
perNos1.m <- list(Mean = apply(perNos1[,-1], 2, mean), SD = apply(perNos1[,-1], 2, sd) )
perNos2.m <- list(Mean = apply(perNos2[,-1], 2, mean), SD = apply(perNos2[,-1], 2, sd) )

# write.csv(cbind.data.frame(Mean = numNos.m$Mean, 
#                            How.much.represents = perNos2.m$Mean, 
#                            How.much.differs= perNos1.m$Mean), 
#           paste0("../output/", webType, "_mean_numNodes.csv") )
# 
# write.csv(numNos, 
#           paste0("../output/", webType, "_numNodes_aqs.csv") )
# write.csv(perNos1, 
#           paste0("../output/", webType, "_numNodesPer1_aqs.csv") )
# write.csv(perNos2, 
#           paste0("../output/", webType, "_numNodesPer2_aqs.csv") )

mutu_numNos.m <- numNos.m
mutu_perNos1.m <- perNos1.m
mutu_perNos2.m <- perNos2.m
rm(list = setdiff(ls(),
                  # deixa so as var de interesse
                  c("mutu_numNos.m","mutu_perNos1.m","mutu_perNos2.m")) )


###############################################################################
# ANTAGONISM
webType <- "antagonism" #"mutualism"

# carrega a work space da etapa 1 (step1)
padrao <- paste0("ws_", webType)
arq <- list.files("txtoutput/", pattern = padrao)
arq <- arq[length(arq)] # arquivo por ultimo (ordem alfabetica)
cat( paste0(arq,"\n") )
load(paste0("txtoutput/",arq))
rm(list = setdiff(ls(),
                  # deixa so as var de interesse
                  c("webType","arqs","gps","listMtOrig","listMtMdOrigSemProblem","lstGps",
                    "mutu_numNos.m","mutu_perNos1.m","mutu_perNos2.m")) )
# remove os dados do HOSTS
smHSTS <- (1:47)[-grep("herb-hosts",arqs)]
arqs <- arqs[smHSTS]
listMtOrig <- listMtOrig[smHSTS]
listMtMdOrigSemProblem <- listMtMdOrigSemProblem[smHSTS]
for(i in 1:length(lstGps)){
  lstGps[[i]] <- lstGps[[i]][smHSTS]
}

numNos <- cbind.data.frame(arqs, 
                           Source = sapply(listMtOrig, nrow) + sapply(listMtOrig, ncol),
                           Final = sapply(listMtMdOrigSemProblem, nrow) + sapply(listMtMdOrigSemProblem, ncol))

for(i in 1:16) {
  numNos <- cbind.data.frame( numNos, sapply(lstGps[[i]], nrow) + sapply(lstGps[[i]], ncol) )
  colnames(numNos)[ncol(numNos)] <- as.character(gps[i,1])
}

perNos1 <- numNos
# qt tirou
for(i in 2:ncol(perNos1)) {
  perNos1[,i] <- ((numNos[,2]-perNos1[,i])/numNos[,2])*100
}

perNos2 <- numNos
# qt manteve
for(i in 2:ncol(perNos2)) {
  perNos2[,i] <- (perNos2[,i]/numNos[,2])*100
}


###############################################################################
# ordena colunas com base na abrangencia
# por combinacao
# ord <- c("A_sp.x.P_sp", 
#          "A_sp.x.P_ge", "A_ge.x.P_sp", "A_ge.x.P_ge", 
#          "A_sp.x.P_fa", "A_fa.x.P_sp", "A_ge.x.P_fa", "A_fa.x.P_ge", "A_fa.x.P_fa", 
#          "A_sp.x.P_or", "A_or.x.P_sp", "A_ge.x.P_or", "A_or.x.P_ge", "A_fa.x.P_or", "A_or.x.P_fa", "A_or.x.P_or")
## col = c("gray10","#c8c8c8","#c8c8c8","#c8c8c8","#c8c8c8","#c8c8c8","#c8c8c8",
##         "gray10","#c8c8c8","#c8c8c8","#c8c8c8","#c8c8c8",
##         "gray10","#c8c8c8","#c8c8c8",
##         "gray10", # sp X sp
##         "white", "white", # com NA
##         "#434343","#707070") # originais (do artigo e padronizado)

# por valores de abrangencia
# ord <- c("A_sp.x.P_sp", "A_ge.x.P_sp", "A_sp.x.P_ge", 
#              "A_fa.x.P_sp", "A_sp.x.P_fa", 
#          "A_ge.x.P_ge", "A_fa.x.P_ge", "A_ge.x.P_fa", 
#              "A_or.x.P_sp", "A_sp.x.P_or", 
#              "A_or.x.P_ge", "A_ge.x.P_or", 
#          "A_fa.x.P_fa", "A_or.x.P_fa", "A_fa.x.P_or", 
#          "A_or.x.P_or")

ord <- c("A_sp.x.P_sp", "A_sp.x.P_ge", "A_sp.x.P_fa", "A_sp.x.P_or",
         "A_ge.x.P_sp", "A_ge.x.P_ge", "A_ge.x.P_fa", "A_ge.x.P_or",
         "A_fa.x.P_sp", "A_fa.x.P_ge", "A_fa.x.P_fa", "A_fa.x.P_or",
         "A_or.x.P_sp", "A_or.x.P_ge", "A_or.x.P_fa", "A_or.x.P_or")

numNos <- numNos[, c("arqs","Source", "Final", ord)]
perNos1 <- perNos1[, c("arqs","Source", "Final", ord)]
perNos2 <- perNos2[, c("arqs","Source", "Final", ord)]


###############################################################################
# mean and sd
numNos.m <- list(Mean = apply(numNos[,-1], 2, mean), SD = apply(numNos[,-1], 2, sd) )
perNos1.m <- list(Mean = apply(perNos1[,-1], 2, mean), SD = apply(perNos1[,-1], 2, sd) )
perNos2.m <- list(Mean = apply(perNos2[,-1], 2, mean), SD = apply(perNos2[,-1], 2, sd) )

# write.csv(cbind.data.frame(Mean = numNos.m$Mean, 
#                            How.much.represents = perNos2.m$Mean, 
#                            How.much.differs= perNos1.m$Mean), 
#           paste0("../output/", webType, "_mean_numNodes.csv") )
# 
# write.csv(numNos, 
#           paste0("../output/", webType, "_numNodes_aqs.csv") )
# write.csv(perNos1, 
#           paste0("../output/", webType, "_numNodesPer1_aqs.csv") )
# write.csv(perNos2, 
#           paste0("../output/", webType, "_numNodesPer2_aqs.csv") )

anta_numNos.m <- numNos.m
anta_perNos1.m <- perNos1.m
anta_perNos2.m <- perNos2.m



















###############################################################################
## Funcao para grafico da quantidade de nos com sd
pltGrfNodesSD <- function(dltLst, label){
  # dltLst <- numNos.m
  # limites superior e inforior do grafico
  mini <- floor(min(dltLst$Mean - dltLst$SD, na.rm = T)*100)/100
  maxi <- ceiling(max(dltLst$Mean + dltLst$SD, na.rm = T)*100)/100
  # renomeia as colunas
  names(dltLst$Mean)[1:2] <- c("source format", "original (ungrouped)") 
  temp <- sub("A_", "", names(dltLst$Mean)[3:length(dltLst$Mean)] )
  temp <- sub("\\.x\\.", " X ", temp)
  temp <- sub("P_", "", temp)
  temp <- gsub("sp", "species", temp)
  temp <- gsub("ge", "genus", temp)
  temp <- gsub("fa", "family", temp)
  temp <- gsub("or", "order", temp)
  names(dltLst$Mean)[3:length(dltLst$Mean)] <- temp
  # acrescenta duas colunas com NA
  dltLst$Mean <- c(dltLst$Mean[1:2], NA, NA, dltLst$Mean[3:length(dltLst$Mean)])
  dltLst$SD <- c(dltLst$SD[1:2], NA, NA, dltLst$SD[3:length(dltLst$SD)])
  # mar=c(5.1, 4.1, 4.1, 2.1)
  par(mar=c(9, 3, 0, 2.1))
  # barCenters recebe as posicoes das colunas
  barCenters <- barplot(dltLst$Mean, las = 2, #
                        col = c("#707070","#434343",# originais (do artigo e padronizado)
                                "white", "white", # com NA
                                "gray10", # sp X sp
                                "#c8c8c8","#c8c8c8","#c8c8c8","#c8c8c8","gray10",
                                "#c8c8c8","#c8c8c8","#c8c8c8","#c8c8c8","gray10",
                                "#c8c8c8","#c8c8c8","#c8c8c8","#c8c8c8","gray10"), 
                        ylim= c(mini,maxi), border = NA, xaxt="n", las = 3)
  title(ylab=label)
  axis(1,at=barCenters, labels=names(dltLst$Mean), tick=F ,pos=mini, las = 2)
  segments(barCenters, dltLst$Mean - dltLst$SD, barCenters,
           dltLst$Mean + dltLst$SD, lwd = 1.0)
  
  arrows(barCenters, dltLst$Mean - dltLst$SD, barCenters,
         dltLst$Mean + dltLst$SD, lwd = 1.5, angle = 90,
         code = 3, length = 0.05)
}

###############################################################################
## Funcao para grafico da quantidade de nos sem SD
pltGrfNodes <- function(dltLst, label){
  # dltLst <- numNos.m
  # limites superior e inforior do grafico
  mini <- floor(min(dltLst$Mean, na.rm = T)*100)/100
  maxi <- ceiling(max(dltLst$Mean, na.rm = T)*100)/100
  # renomeia as colunas
  names(dltLst$Mean)[1:2] <- c("source format", "original (ungrouped)") 
  temp <- sub("A_", "", names(dltLst$Mean)[3:length(dltLst$Mean)] )
  temp <- sub("\\.x\\.", " X ", temp)
  temp <- sub("P_", "", temp)
  temp <- gsub("sp", "species", temp)
  temp <- gsub("ge", "genus", temp)
  temp <- gsub("fa", "family", temp)
  temp <- gsub("or", "order", temp)
  names(dltLst$Mean)[3:length(dltLst$Mean)] <- temp
  # acrescenta duas colunas com NA
  dltLst$Mean <- c(dltLst$Mean[1:2], NA, NA, dltLst$Mean[3:length(dltLst$Mean)])
  # inverte a ordem do vetor (para grafico de barras " horiz = T "
  dltLst$Mean <- dltLst$Mean[length(dltLst$Mean):1]
  # mar=c(5.1, 4.1, 4.1, 2.1)
  par(mar=c(3, 9, 1, 2.1))
  # barCenters recebe as posicoes das colunas
  barCenters <- barplot(dltLst$Mean, las = 1, horiz = T,
                        col = c("gray10","#c8c8c8","#c8c8c8","#c8c8c8","#c8c8c8",
                                "gray10","#c8c8c8","#c8c8c8","#c8c8c8","#c8c8c8",
                                "gray10","#c8c8c8","#c8c8c8","#c8c8c8","#c8c8c8",
                                "gray10", # sp X sp
                                "white", "white", # com NA
                                "#434343","#707070"), # originais (do artigo e padronizado)
                        xlim= c(mini,maxi), border = NA, yaxt = "n")
  axis(2,at=barCenters, labels=names(dltLst$Mean), tick=F, las = 1, pos = 5)
  title(xlab=label,line = 2)
}

###############################################################################
## Funcao para grafico perdentagem de nos
pltGrfPercent <- function(dltLst, label){
  # dltLst <- numNos.m
  # limites superior e inforior do grafico
  mini <- 0
  maxi <- 100
  # renomeia as colunas
  names(dltLst$Mean)[1:2] <- c("source format", "original (ungrouped)") 
  temp <- sub("A_", "", names(dltLst$Mean)[3:length(dltLst$Mean)] )
  temp <- sub("\\.x\\.", " X ", temp)
  temp <- sub("P_", "", temp)
  temp <- gsub("sp", "species", temp)
  temp <- gsub("ge", "genus", temp)
  temp <- gsub("fa", "family", temp)
  temp <- gsub("or", "order", temp)
  names(dltLst$Mean)[3:length(dltLst$Mean)] <- temp
  # acrescenta duas colunas com NA
  dltLst$Mean <- c(dltLst$Mean[1:2], NA, NA, dltLst$Mean[3:length(dltLst$Mean)])
  dltLst$SD <- c(dltLst$SD[1:2], NA, NA, dltLst$SD[3:length(dltLst$SD)])
  # inverte a ordem do vetor (para grafico de barras " horiz = T "
  dltLst$Mean <- dltLst$Mean[length(dltLst$Mean):1]
  dltLst$SD <- dltLst$SD[length(dltLst$SD):1]
  # mar=c(5.1, 4.1, 4.1, 2.1)
  par(mar=c(3, 7.8, 1, 2.1))
  # barCenters recebe as posicoes das colunas
  barCenters <- barplot(dltLst$Mean, las = 1, horiz = T,
                        col = c("gray10","#c8c8c8","#c8c8c8","#c8c8c8","#c8c8c8",
                                "gray10","#c8c8c8","#c8c8c8","#c8c8c8","#c8c8c8",
                                "gray10","#c8c8c8","#c8c8c8","#c8c8c8","#c8c8c8",
                                "gray10", # sp X sp
                                "white", "white", # com NA
                                "#434343","#707070"), # originais (do artigo e padronizado)
                        xlim= c(mini,maxi), border = NA, yaxt = "n")
  axis(2,at=barCenters, labels=names(dltLst$Mean), tick=F, las = 1, pos = 5)
  title(xlab=label,line = 2)
}
# fim das funcoes #############################################################
###############################################################################



###############################################################################
# graficos
# para uma figura grande de width = 174 e height = 174
svg(filename = paste0("../output/mutu-anta-numNos-per_", format(Sys.time(), "%Y-%m-d%d"),".svg"),
    width = 7, height = 10, # tamanho equivalente a 177x177 mm 
    pointsize = 14) # 18 equivale a uma fonte de 12 de Arial
par(mfcol= c(3,2), oma = c(0, 0, 2, 0))

# mutualism
pltGrfNodes(mutu_numNos.m, "num. of nodes")
pltGrfPercent(mutu_perNos2.m, "% of nodes")

plot.new()

# antagonism
pltGrfNodes(anta_numNos.m, "num. of nodes")
pltGrfPercent(anta_perNos2.m, "% of nodes")

# legend
par(mar=c(0.5, 0.5, 1, 1))
plot.new()

legend("topright",
       c("Nodes of the original source",
         "Nodes remaining from original source\nafter standardization of nomenclature\n",
         "Both insect and plants grouped under\nthe same taxonomic level",
         "Others taxonomic grouping types"),
       fill = c("#707070","#434343","gray10","#c8c8c8"),
       bty = "n", border = "white")
mtext("Mutualism           Antagonism", outer = TRUE, cex = 0.7, line = 1)

dev.off()

par(mfcol= c(1,1))






svg(filename = paste0("../output/mutu-anta-per_", format(Sys.time(), "%Y-%m-d%d"),".svg"),
    width = 7, height = 8, # tamanho equivalente a 177x177 mm 
    pointsize = 14) # 18 equivale a uma fonte de 12 de Arial
par(mfcol= c(2,2), oma = c(0, 0, 2, 0))

# mutualism
pltGrfPercent(mutu_perNos2.m, "% of nodes")

plot.new()

# antagonism
pltGrfPercent(anta_perNos2.m, "% of nodes")

# legend
par(mar=c(0.5, 0.5, 1, 1))
plot.new()

legend("topright",
       c("Nodes of the original source",
         "Nodes remaining from original source\nafter standardization of nomenclature\n",
         "Both insect and plants grouped under\nthe same taxonomic level",
         "Others taxonomic grouping types"),
       fill = c("#707070","#434343","gray10","#c8c8c8"),
       bty = "n", border = "white")
mtext("Mutualism           Antagonism", outer = TRUE, cex = 0.7, line = 1)

dev.off()

par(mfcol= c(1,1))







pltGrfPercent(perNos1.m, "% of nodes")



pltGrfNodesSD(mutu_perNos2.m, "% of nodes")


svg(filename = paste0("../output/mutu-anta-numNosSD_", format(Sys.time(), "%Y-%m-d%d"),".svg"),
    width = 7, height = 8, # tamanho equivalente a 177x177 mm 
    pointsize = 14) # 18 equivale a uma fonte de 12 de Arial
par(mfcol= c(2,2), oma = c(0, 0, 2, 0))
pltGrfNodesSD(mutu_numNos.m, "num. of nodes")
plot.new()
pltGrfNodesSD(anta_numNos.m, "num. of nodes")


par(mar=c(0.5, 0.5, 1, 1))
plot.new()

legend("topright",
       c("Nodes of the original source",
         "Nodes remaining from original source\nafter standardization of nomenclature\n",
         "Both insect and plants grouped under\nthe same taxonomic level",
         "Others taxonomic grouping types"),
       fill = c("#707070","#434343","gray10","#c8c8c8"),
       bty = "n", border = "white")
mtext("Mutualism           Antagonism", outer = TRUE, cex = 0.7, line = 1)

dev.off()

par(mfcol= c(1,1))


svg(filename = paste0("../output/mutu-anta-perNosSD_", format(Sys.time(), "%Y-%m-d%d"),".svg"),
    width = 7, height = 8, # tamanho equivalente a 177x177 mm 
    pointsize = 14) # 18 equivale a uma fonte de 12 de Arial
par(mfcol= c(2,2), oma = c(0, 0, 2, 0))
pltGrfNodesSD(mutu_perNos2.m, "% of nodes")
plot.new()
pltGrfNodesSD(anta_perNos2.m, "% of nodes")


par(mar=c(0.5, 0.5, 1, 1))
plot.new()

legend("topright",
       c("Nodes of the original source",
         "Nodes remaining from original source\nafter standardization of nomenclature\n",
         "Both insect and plants grouped under\nthe same taxonomic level",
         "Others taxonomic grouping types"),
       fill = c("#707070","#434343","gray10","#c8c8c8"),
       bty = "n", border = "white")
mtext("Mutualism           Antagonism", outer = TRUE, cex = 0.7, line = 1)

dev.off()

par(mfcol= c(1,1))

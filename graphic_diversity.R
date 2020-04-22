###############################################################################
# change the work diretory
#setwd("<your_work_diretory_path>")


###############################################################################
# MUTUALISM
# load the image (all vars)
load("txtoutput/ws_mutualism_2019-07-d03_v01.RData")
# retirar var que nao serao utilizadas
rm(list = setdiff( ls(), c("listMtMdOrig","arqs","plotdiversity")) )
# importa a classificacao taxonomica
origNvplanTx <- read.csv("../output/mutualism_gNvplanTx_2019-07-d03.csv")[,-1]
origNvanimTx <- read.csv("../output/mutualism_gNvanimTx_2019-07-d03.csv")[,-1]
## listas de especies dos arquivos sem HOSTS
# vetor que vao armazenar nomes de todas sp
AnimNm <- vector(mode= "character")
PlanNm <- vector(mode= "character")
for(i in 1:47){ # mutualismo
  # armazena as sp
  # para plantas nas linhas
  PlanNm <- append(PlanNm, as.vector(rownames(listMtMdOrig[[i]])))
  # para animas nas colunas
  AnimNm <- append(AnimNm, as.vector(colnames(listMtMdOrig[[i]])))
}
# select all classification of these sp from the gNvplanTx and gNvanimTx
mutuPlt <- origNvplanTx[is.element(origNvplanTx[,1],PlanNm),]
mutuAnm <- origNvanimTx[is.element(origNvanimTx[,1],AnimNm),]


# ANTAGONISM
# load the image (all vars)
load("txtoutput/ws_antagonism_2019-07-d25_v01.RData")
# retirar var que nao serao utilizadas
rm(list = setdiff( ls(), c("listMtMdOrig","arqs", "plotdiversity", 
                           "mutuPlt", "mutuAnm")) )
# importa a classificacao taxonomica
origNvplanTx <- read.csv("../output/antagonism_gNvplanTx_2019-07-d25.csv")[,-1]
origNvanimTx <- read.csv("../output/antagonism_gNvanimTx_2019-07-d25.csv")[,-1]
## listas de especies dos arquivos sem HOSTS
# vetor que vao armazenar nomes de todas sp
AnimNm <- vector(mode= "character")
PlanNm <- vector(mode= "character")
for(i in (1:47)[-grep("herb-hosts",arqs)]){
  # armazena as sp
  # para plantas nas linhas
  PlanNm <- append(PlanNm, as.vector(rownames(listMtMdOrig[[i]])))
  # para animas nas colunas
  AnimNm <- append(AnimNm, as.vector(colnames(listMtMdOrig[[i]])))
}
# select all classification of these sp from the gNvplanTx and gNvanimTx
antaPlt <- origNvplanTx[is.element(origNvplanTx[,1],PlanNm),]
antaAnm <- origNvanimTx[is.element(origNvanimTx[,1],AnimNm),]

# junta as sp das mt antagonista e mutualistica
ttlPlt <- unique(rbind(antaPlt,mutuPlt))
ttlAnm <- unique(rbind(antaAnm,mutuAnm))

###############################################################################
# tabela do artigo
tbl <- data.frame(Organism = rep(c("Animals","Plants"), each = 3),
                  Data.set = rep(c("Mutualistic","Antagonistic","Total"), time = 2),
                  Genus = vector("numeric", length = 6),
                  Family = vector("numeric", length = 6),
                  Order = vector("numeric", length = 6))
# Animal
tbl$Genus[1] <- length(unique(mutuAnm$genus)) # mutu
tbl$Genus[2] <- length(unique(antaAnm$genus)) # anta
tbl$Genus[3] <- length(unique(ttlAnm$genus)) # ttl
tbl$Family[1] <- length(unique(mutuAnm$family)) # mutu
tbl$Family[2] <- length(unique(antaAnm$family)) # anta
tbl$Family[3] <- length(unique(ttlAnm$family)) # ttl
tbl$Order[1] <- length(unique(mutuAnm$order)) # mutu
tbl$Order[2] <- length(unique(antaAnm$order)) # anta
tbl$Order[3] <- length(unique(ttlAnm$order)) # ttl
# planta
tbl$Genus[4] <- length(unique(mutuPlt$genus)) # mutu
tbl$Genus[5] <- length(unique(antaPlt$genus)) # anta
tbl$Genus[6] <- length(unique(ttlPlt$genus)) # ttl
tbl$Family[4] <- length(unique(mutuPlt$family)) # mutu
tbl$Family[5] <- length(unique(antaPlt$family)) # anta
tbl$Family[6] <- length(unique(ttlPlt$family)) # ttl
tbl$Order[4] <- length(unique(mutuPlt$order)) # mutu
tbl$Order[5] <- length(unique(antaPlt$order)) # anta
tbl$Order[6] <- length(unique(ttlPlt$order)) # ttl
# salva csv
write.csv(tbl, file = "../output/table_quantClassTax.csv")





# Funcao para plotar grafico
###############################################################################  
# graficos apresentando a diversidade taxonomica das redes
## grafico de barras de diversidade
###############################################################################
plotdiversity <- function(gNvplanTx,gNvanimTx,webType){
  tiff(paste0("../output/", webType, "_numSpColor.tif"),
       width = 174,height = 233,units = "mm",res = 300)
  par(mfcol= c(3,2), mar= c(2,8,1.8,1)+0.1, mgp= c(1,1,0))
  # plantas
  for (i in c(9,7,5)){
    temp <- unique(gNvplanTx[ ,i])
    soma <- vector(mode = "numeric", length = length(temp))
    for (j in 1:length(temp)){
      soma[j] <- sum(gNvplanTx[ ,i] == temp[j])
    }
    soma <- data.frame(name = temp[order(soma, decreasing =T)], 
                       sum = soma[order(soma, decreasing =T)])
    if(i == 9){
      barplot(soma[25:1,2], names.arg = soma[25:1,1], las = 1, horiz = T,
              beside = T, col = "#008033", border = NA,
              main = paste0("Number of plant species by\n ", colnames(gNvplanTx)[i],
                            " - total of ",length(soma[,1]) ))
    }else{
      barplot(soma[25:1,2], names.arg = soma[25:1,1], las = 1, horiz = T,
              beside = T, col = "#008033", border = NA,
              main = paste0(colnames(gNvplanTx)[i],
                            " - total of ",length(soma[,1]) ))
    }
  }
  # animais
  for (i in c(9,7,5)){
    temp <- unique(gNvanimTx[ ,i])
    soma <- vector(mode = "numeric", length = length(temp))
    for (j in 1:length(temp)){
      soma[j] <- sum(gNvanimTx[ ,i] == temp[j])
    }
    soma <- data.frame(name = temp[order(soma, decreasing =T)], 
                       sum = soma[order(soma, decreasing =T)])
    if(i == 9){
      barplot(soma[25:1,2], names.arg = soma[25:1,1], las = 1, horiz = T,
              beside = T, col = "#ff2a2a", border = NA,
              main = paste0("Number of animal species by\n ", colnames(gNvplanTx)[i],
                            " - total of ",length(soma[,1]) ))
    }else{
      barplot(soma[25:1,2], names.arg = soma[25:1,1], las = 1, horiz = T,
              beside = T, col = "#ff2a2a", border = NA,
              main = paste0(colnames(gNvplanTx)[i],
                            " - total of ",length(soma[,1]) ))
    }
  }
  par(mfcol= c(1,1), mar= c(5,4,4,2)+0.1, mgp= c(3,1,0))
  dev.off()
}
###############################################################################
# Fim da funcao


###############################################################################
# graficos de quantidade de especies e demais classificacoes
plotdiversity(antaPlt,antaAnm,"antagonism")
plotdiversity(mutuPlt,mutuAnm,"mutualism")

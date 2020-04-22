###############################################################################
# Calculo da correlacao entre as metricas                                     #
###############################################################################

# define diretorio de trabalho
#setwd("<your_work_diretory_path>")



#####################################################################
# Leitura das informacoes das redes
## leitura dos dados de redes mutualisticas
mutu_tax <- read.csv("../output/mutualism_allMetricsAllGps.csv")[,-1] # so tem 

# colunas para analise dos agrupamentos
mutu_tax <- mutu_tax[, -3] # tira a col tipo (3)

# separa original das agrupadas por taxonomia
mutu_ori <- mutu_tax[mutu_tax$GroupType == "Original",] # so info das originais
#mutu_tax <- mutu_tax[mutu_tax$GroupType != "Original",] # so dos agrup taxo
rownames(mutu_ori) <- mutu_ori$ArqName

# metricas presentes
metrics <- colnames(mutu_ori)[5:20]


## leitura dos dados de redes antagonicas
anta_tax <- read.csv("../output/antagonism_allMetricsAllGps.csv")[,-1] # so tem 
# remove as redes do HOSTS
anta_tax <- anta_tax[-grep("herb-hosts",anta_tax[,1]), ]
# atualiza factors
temp <- as.character(anta_tax$ArqName)
anta_tax$ArqName <- factor(temp, levels = unique(temp))

# colunas para analise dos agrupamentos
anta_tax <- anta_tax[, -3] # tira a col tipo (3)

# separa original das agrupadas por taxonomia
anta_ori <- anta_tax[anta_tax$GroupType == "Original",] # so info das originais
#anta_tax <- anta_tax[anta_tax$GroupType != "Original",] # so dos agrup taxo
rownames(anta_ori) <- anta_ori$ArqName


###############################################################################
# matriz so com os valores das metricas

# mutu
mutu_mt <- as.matrix(mutu_tax[,7:18])

tiff(paste0("../output/_scatterplot_", "mutualism","_newMetrics_",
            format(Sys.time(), "_%Y-%m-d%d"),"-",
            gsub(":", "", format(Sys.time(), "%X") ), ".tif"),
     width = 87,height = 87,units = "mm",res = 300)
pairs(~connectance+NODF+Modularity, data = mutu_mt)
dev.off()

pairs(~nestedness+links.per.species+web.asymmetry, data = mutu_mt)
pairs(~niche.overlap.HL+niche.overlap.LL+robustness.HL+robustness.LL, data = mutu_mt)

tiff(paste0("../output/_scatterplot_", "mutualism", "_oldMetrics_",
            format(Sys.time(), "_%Y-%m-d%d"),"-",
            gsub(":", "", format(Sys.time(), "%X") ), ".tif"),
     width = 174,height = 174,units = "mm",res = 300)
pairs(~connectance+nestedness+links.per.species+web.asymmetry+
        niche.overlap.HL+niche.overlap.LL+robustness.HL+robustness.LL, data = mutu_mt)
dev.off()


mutu_cor <- cor(mutu_mt, use = "complete.obs", method = "spearman") 
#  "then missing values are handled by casewise deletion (and if there are no 
#complete cases, that gives an error)"
write.csv(mutu_cor, paste0("../output/_corSpearman_", "mutualism", "_allMetrics_",
                           format(Sys.time(), "_%Y-%m-d%d"),"-",
                           gsub(":", "", format(Sys.time(), "%X") ), ".csv"))

# anta
anta_mt <- as.matrix(anta_tax[,7:18])

tiff(paste0("../output/_scatterplot_", "antagonis","_newMetrics_",
            format(Sys.time(), "_%Y-%m-d%d"),"-",
            gsub(":", "", format(Sys.time(), "%X") ), ".tif"),
     width = 87,height = 87,units = "mm",res = 300)
pairs(~connectance+NODF+Modularity, data = anta_mt)
dev.off()

pairs(~nestedness+links.per.species+web.asymmetry, data = anta_mt)
pairs(~niche.overlap.HL+niche.overlap.LL+robustness.HL+robustness.LL, data = anta_mt)

tiff(paste0("../output/_scatterplot_", "antagonis", "_oldMetrics_",
            format(Sys.time(), "_%Y-%m-d%d"),"-",
            gsub(":", "", format(Sys.time(), "%X") ), ".tif"),
     width = 174,height = 174,units = "mm",res = 300)
pairs(~connectance+nestedness+links.per.species+web.asymmetry+
        niche.overlap.HL+niche.overlap.LL+robustness.HL+robustness.LL, data = anta_mt)
dev.off()


anta_cor <- cor(anta_mt, use = "complete.obs", method = "spearman")
write.csv(anta_cor, paste0("../output/_corSpearman_", "antagonis", "_allMetrics_",
                           format(Sys.time(), "_%Y-%m-d%d"),"-",
                           gsub(":", "", format(Sys.time(), "%X") ), ".csv"))


###############################################################################
# grafico linear do valor de correlacao por metrica X todas as outras metricas

# tipo de agrupamento em ordem
ord <- c("original",
         "A_sp.x.P_sp", "A_sp.x.P_ge", "A_sp.x.P_fa", "A_sp.x.P_or",
         "A_ge.x.P_sp", "A_ge.x.P_ge", "A_ge.x.P_fa", "A_ge.x.P_or",
         "A_fa.x.P_sp", "A_fa.x.P_ge", "A_fa.x.P_fa", "A_fa.x.P_or",
         "A_or.x.P_sp", "A_or.x.P_ge", "A_or.x.P_fa", "A_or.x.P_or")
# tipo de agrupamento em ordem abreviado
ordAbrv <- sub("A_", "", ord)
ordAbrv <- sub("\\.x\\.", "×", ordAbrv)
ordAbrv <- sub("P_", "", ordAbrv)
ordAbrv[1] <- "Ori"

# nomes das metricas 3X4
mts <- c("connectance", "Modularity", "NODF", "nestedness", 
         "links.per.species", "web.asymmetry","number.of.compartments", "N.modules",
         "robustness.HL", "robustness.LL", 
         "niche.overlap.HL", "niche.overlap.LL")
         
# cores
crMts <- c("#404040", "#737373", "#A6A6A6", "#ff6600",
           "#aa8800", "#00aa00","#cc00ff", "#ff00e6",
           "#ff0000", "#aa0000", 
           "#0000ff", "#367ef1")
## ANTA
# cria uma lista para armazenar a correlacao para:
# cada um dos 16 tipos de agrup + original
lstCorrAnta <- vector("list", length = 17)
lstCorrAnta[[1]] <- cor(as.matrix(anta_ori[,7:18]), 
                        use = "complete.obs", method = "spearman")
for(i in 2:17){
  lstCorrAnta[[i]] <- cor(as.matrix(anta_tax[anta_tax$GroupType == ord[i],
                                               7:18]),
                            use = "complete.obs", method = "spearman")
}
names(lstCorrAnta) <- ord

# transforma em um data frame
dfCorrAnt <- as.data.frame(lstCorrAnta[[1]])
dfCorrAnt <- cbind.data.frame(metrics = rownames(dfCorrAnt), dfCorrAnt)
dfCorrAnt <- cbind.data.frame(group = names(lstCorrAnta)[1], dfCorrAnt)
rownames(dfCorrAnt) <- 1:12
for(i in 2:17){
  temp <- as.data.frame(lstCorrAnta[[i]])
  temp <- cbind.data.frame(metrics = rownames(temp), temp)
  temp <- cbind.data.frame(group = names(lstCorrAnta)[i], temp)
  rownames(temp) <- (1:12)+12*(i-1)
  dfCorrAnt <- rbind.data.frame(dfCorrAnt, temp)
}
# ordena os grupos
dfCorrAnt$group <- factor(dfCorrAnt$group, levels = ord)
dfCorrAnt <- dfCorrAnt[order(dfCorrAnt$group),]





## MUTU
# cria uma lista para armazenar a correlacao para:
# cada um dos 16 tipos de agrup + original
lstCorrMut <- vector("list", length = 17)
lstCorrMut[[1]] <- cor(as.matrix(mutu_ori[,7:18]), 
                       use = "complete.obs", method = "spearman")
for(i in 2:17){
  lstCorrMut[[i]] <- cor(as.matrix(mutu_tax[mutu_tax$GroupType == ord[i],
                                            7:18]),
                         use = "complete.obs", method = "spearman")
}
names(lstCorrMut) <- ord

# transforma em um data frame
dfCorrMut <- as.data.frame(lstCorrMut[[1]])
dfCorrMut <- cbind.data.frame(metrics = rownames(dfCorrMut), dfCorrMut)
dfCorrMut <- cbind.data.frame(group = names(lstCorrMut)[1], dfCorrMut)
rownames(dfCorrMut) <- 1:12
for(i in 2:17){
  temp <- as.data.frame(lstCorrMut[[i]])
  temp <- cbind.data.frame(metrics = rownames(temp), temp)
  temp <- cbind.data.frame(group = names(lstCorrMut)[i], temp)
  rownames(temp) <- (1:12)+12*(i-1)
  dfCorrMut <- rbind.data.frame(dfCorrMut, temp)
}
# ordena os grupos
dfCorrMut$group <- factor(dfCorrMut$group, levels = ord)
dfCorrMut <- dfCorrMut[order(dfCorrMut$group),]




###############################################################################
# grafico

svg(paste0("../output/corrLinePlotsAnta",
           format(Sys.time(), "_%Y-%m-d%d"),
           gsub(":", "", format(Sys.time(), "%X") ), ".svg"),
    width = 8, height = 9, pointsize = 14) # 18 equivale a fonte arial 12

  layout(matrix(c(1,2,3,
                  4,5,6,
                  7,8,9,
                  10,11,12,
                  13,14,15), # legenda
                nrow = 5, ncol = 3, byrow = TRUE))

  # 5.1 4.1 4.1 2.1
  par(mar=c(4.5,3.5,0.5,0.5)) # margin around the plot
  # graficos
  # ant
  for(j in 1:12){
    plot(rep(0,17),1:17,
         ylim= c(-1,1), 
         xlim= c(0,17), xaxt = "n",
         type = "n", bty= "l", ylab = "", xlab = "")
    title(ylab = mts[j], mgp=c(2.2,0,0)) # posicao
    title(xlab = "grouping type", mgp=c(3.3,0,0)) # posicao
    # eixo com os agrupamentos
    axis(1,at= seq(1, 17, 1), # posicao de inicio, passos, fim
         labels= ordAbrv, las = 2) # qtdade de ticks
    # 70% de correlacao
    abline(0.7,0, lty = 2); abline(-0.7,0, lty = 2)
    
    # cada a corr dessa metrica com cada metrica
    for(i in 4:12){
      temp <- dfCorrAnt[dfCorrAnt$metrics == mts[i], mts[j]]
      temp <- cbind(1:17,temp)
      lines(x = temp[,1], y = temp[,2], col = crMts[i],lwd = 1.5)
    }
    for(i in 1:3){
      temp <- dfCorrAnt[dfCorrAnt$metrics == mts[i], mts[j]]
      temp <- cbind(1:17,temp)
      lines(x = temp[,1], y = temp[,2], col = crMts[i],lwd = 2.6)
    }
  }
  plot.new()
  legend("topleft", bty= "n", # posicao e sem borda
         lwd= c(2.6,2.6,2.6, 1.5), seg.len= 2,  # espessura e comprimento
         lty= 1, # tipo do traco das linhas
         legend= mts[1:4],
         col= crMts[1:4])
  plot.new()
  legend("topleft", bty= "n", # posicao e sem borda
         lwd= 1.5, seg.len= 2,  # espessura e comprimento
         lty= 1, # tipo do traco das linhas
         legend= mts[5:8],
         col= crMts[5:8])
  plot.new()
  legend("topleft", bty= "n", # posicao e sem borda
         lwd= 1.5, seg.len= 2,  # espessura e comprimento
         lty= 1, # tipo do traco das linhas
         legend= mts[9:12],
         col= crMts[9:12])
  dev.off()
  
  
  
  
  #########################################
  ## MUTU
  svg(paste0("../output/corrLinePlotsMutu",
             format(Sys.time(), "_%Y-%m-d%d"),
             gsub(":", "", format(Sys.time(), "%X") ), ".svg"),
      width = 8, height = 9, pointsize = 14) # 18 equivale a fonte arial 12
  
  layout(matrix(c(1,2,3,
                  4,5,6,
                  7,8,9,
                  10,11,12,
                  13,14,15), # legenda
                nrow = 5, ncol = 3, byrow = TRUE))
  
  # 5.1 4.1 4.1 2.1
  par(mar=c(4.5,3.5,0.5,0.5)) # margin around the plot
  # graficos
  #mutu
  for(j in 1:12){
    plot(rep(0,17),1:17,
         ylim= c(-1,1), 
         xlim= c(0,17), xaxt = "n",
         type = "n", bty= "l", ylab = "", xlab = "")
    title(ylab = mts[j], mgp=c(2.2,0,0)) # posicao
    title(xlab = "grouping type", mgp=c(3.3,0,0)) # posicao
    #title(ylab= mts[j]); title(ylab= mts[j])
    # eixo com os agrupamentos
    axis(1,at= seq(1, 17, 1), # posicao de inicio, passos, fim
         labels= ordAbrv, las = 2) # qtdade de ticks
    # 70% de correlacao
    abline(0.7,0, lty = 2); abline(-0.7,0, lty = 2)
    
    # cada a corr dessa metrica com cada metrica
    for(i in 4:12){
      temp <- dfCorrMut[dfCorrMut$metrics == mts[i], mts[j]]
      temp <- cbind(1:17,temp)
      lines(x = temp[,1], y = temp[,2], col = crMts[i],lwd = 1.5)
    }
    for(i in 1:3){
      temp <- dfCorrMut[dfCorrMut$metrics == mts[i], mts[j]]
      temp <- cbind(1:17,temp)
      lines(x = temp[,1], y = temp[,2], col = crMts[i],lwd = 2.6)
    }
  }
  plot.new()
  legend("topleft", bty= "n", # posicao e sem borda
         lwd= c(2.6,2.6,2.6, 1.5), seg.len= 2,  # espessura e comprimento
         lty= 1, # tipo do traco das linhas
         legend= mts[1:4],
         col= crMts[1:4])
  plot.new()
  legend("topleft", bty= "n", # posicao e sem borda
         lwd= 1.5, seg.len= 2,  # espessura e comprimento
         lty= 1, # tipo do traco das linhas
         legend= mts[5:8],
         col= crMts[5:8])
  plot.new()
  legend("topleft", bty= "n", # posicao e sem borda
         lwd= 1.5, seg.len= 2,  # espessura e comprimento
         lty= 1, # tipo do traco das linhas
         legend= mts[9:12],
         col= crMts[9:12])
dev.off()
  


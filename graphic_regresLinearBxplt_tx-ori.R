###############################################################################
# Regressao linear                                                            #
# 1 - como os dados agrupados representam os originais                        #
# metricas originais ~ metricas agrupadas                                     #
###############################################################################

# define diretorio de trabalho
#setwd("<your_work_diretory_path>")


###############################################################################
###############################################################################
# colunas para analise dos agrupamentos                                       #
###############################################################################
adColAgrup <- function(dtMetrics){
  pltGrp <- sub(".........","", dtMetrics$GroupType)
  pltGrp <- sub("ginal","", pltGrp)
  anmGrp <- sub("A_","", dtMetrics$GroupType)
  anmGrp <- sub(".x.P_..","", anmGrp)
  anmGrp <- sub("ginal","", anmGrp)
  dtMetrics <- cbind(anmGrp,pltGrp,dtMetrics)
  txNm <- c('Ori','sp','ge','fa','or')
  dtMetrics$pltGrp <- factor(as.character(dtMetrics$pltGrp),
                             levels = txNm)
  dtMetrics$anmGrp <- factor(as.character(dtMetrics$anmGrp),
                             levels = txNm)
  
  return(dtMetrics)
}

###############################################################################
# Funcao calculo da Regressao linear                                          #
###############################################################################
# valores das metricas nas redes originais X redes agrupadas
# precisa dos valores para as matrizes originais e das agrupadas
lstReg <- function(orMt,txGrp){
  # lista de data frames para armazenar os valores da regressao
  # cada item da lista
  temp <- as.data.frame(unique(txGrp[,c(1:2,4)])) # tipos de agrupamento (ex.: 16)
  tam <- length(temp[,1]) # qtdade de agrupamentos (ex.: 16)
  j <- 4
  # construi o data.frame a ser preenchido com as info da regressao
  for(i in (1:ncol(txGrp))[colnames(txGrp)==
                           "connectance"]:ncol(txGrp)){
    temp <- cbind(temp, vector('numeric', length = tam))
    colnames(temp)[j] <- colnames(orMt)[i]
    j <- j + 1
  }
  ### regressao linear
  # parametros da regressao para analise + modelo de regressao
  metrics.lm <- vector('list', length = 7) 
  # sigma: Residual standard error
  # adj.r.squared: r quadrado ajustado
  # a$df[2]: graus de liberdade
  # coefficients[2,4]: significancia
  # coefficients[1,1]: intersept (b)
  # coefficients[2,1]: angulo da reta (a)
  names(metrics.lm) <- c('sigma','adj.r.squared','f.degree',
                         'signif','intersept','angle','lm')
  # data.frame para os parametros das regressoes
  for(i in 1:6){
    metrics.lm[[i]] <- temp # temp = data.frame vazio (so com info de agrupamento)
  }
  # lista de regressoes por agrupamento...
  metrics.lm[[7]] <- vector('list', length = tam)
  names(metrics.lm[[7]]) <- as.character(temp$GroupType)
  # ... e por metrica
  for( i in 1: tam){ #length(8:ncol(txGrp)) connectancia ate modularidade
    metrics.lm[[7]][[i]] <- vector('list', length = length((1:ncol(txGrp))[colnames(txGrp)==
                                                                             "connectance"]:
                                                             ncol(txGrp)) )
    # nomes das metricas
    names(metrics.lm[[7]][[i]]) <- colnames(metrics.lm[[1]])[ 4: ncol(metrics.lm[[1]]) ]
  }
  # calculo dos modelos
  for(i in 4: ncol(metrics.lm[[1]])){ # cd metrica
    for(j in 1:tam){ # cd agrupamento 
      # regresssao lm: dados originais ~ agrupado[j]
      metrics.lm[[7]][[j]][[i-3]] <- lm(orMt[, i+5]~ # valores ori
                                   txGrp[txGrp$GroupType == # valores por agrupamento
                                           as.character(metrics.lm[[1]]$GroupType[j]), i+5],
                                   na.action = na.omit)
      # dados da regressao
      a <- summary(metrics.lm[[7]][[j]][[i-3]])
      metrics.lm$sigma[j,i] <- a$sigma # Residual standard error
      metrics.lm$adj.r.squared[j,i] <- a$adj.r.squared # r quadrado ajustado
      metrics.lm$f.degree[j,i] <- a$df[2] # graus de liberdade
      metrics.lm$signif[j,i] <- a$coefficients[2,4] # significancia
      metrics.lm$intersept[j,i] <- a$coefficients[1,1] # intersept (b)
      metrics.lm$angle[j,i] <- a$coefficients[2,1] # angulo da reta (a)
    }
  }
  return(metrics.lm)
}

###############################################################################
# Cria e salva um grafico em svg para cada metrica                            #
###############################################################################
grfLm <- function(lmMetrics, orMt,txGrp, metric,
                  txNm, webType, fourColors, agrp = "ori-tx", frmto = "tela"){
  # parametros para o grafico
  cores <- fourColors # 4 cores para representar cd agrup de planta
  simbo <- c(15,2,3,4) # simbolos para os ptos do grafico
  
  ### Grafico
  if (frmto == "tiff") {
    tiff(paste0("../output/", webType, "_regressionGraph_",metric,"_",agrp,
                format(Sys.time(), "_%Y-%m-d%d-"),
                gsub(":", "", format(Sys.time(), "%X") ), ".tif"),
         width = 174,height = 55,units = "mm",res = 300)
  }
  
  if (frmto == "svg") {
    svg(paste0("../output/", webType, "_regressionGraph_",metric,"_",agrp,
               format(Sys.time(), "_%Y-%m-d%d"),
               gsub(":", "", format(Sys.time(), "%X") ), ".svg"),
        width = 7, height = 3, pointsize = 14) # 18 equivale a fonte arial 12
  }

  # from https://www.statmethods.net/advgraphs/layout.html
  layout(matrix(c(1,2,3,4, # fig 1, 2, 3 e 4
                  5,5,6,6), # fig 5 e 6 (cd uma ocupando "2 espacos")
                nrow = 2, ncol = 4, byrow = TRUE))
  
  for(j in 2:length(txNm) ){ # agrupamento de animais
    # teste para saber se plota ou nao (ate q class vai sp->ge ou sp->or)
    tmpAnm <- txGrp[txGrp$anmGrp == txNm[j],]
    i <- (1:ncol(orMt))[colnames(orMt) == metric] # metrica foco
    par(mar=c(2.5,3.5,1.5,0.5), mgp= c(1.5,0.25,0))
    plot(orMt[,i]~tmpAnm[tmpAnm$pltGrp == txNm[j],i],
         ylim= c(0,max(c(txGrp[,i],orMt[,i]), na.rm = T)), 
         xlim= c(0,max(c(txGrp[,i],orMt[,i]), na.rm = T)),
         type = "n", bty= "l", ylab = "", xlab = "Grouped metric values")
    abline(0,1, lty = 2)
    # tipo do agrupamento de animais
    title(paste0("Insects grouped by ", txNm[j]), 
          adj=0.8, # posicao
          cex.main = 1, # tamanho em rel ao cex
          font.main = 1) # 1: normal; 2: bold; 3: italic; 4: bold italic 
    # nome da metrica no eixo y e que sao os valores do agrupa original (y)
    if(j == 2){
      title(ylab = paste0(metric,"\nOriginal metric values") ) # nome da metrica
    }
    # pontos e linhas (valore de y [original] em funcao de x [agrupado])
    for(w in 2:length(txNm)){ # para sp, ge, fa e or do grupamento de plantas
      # tipo do agrupamento (spXsp, spXge, geXsp...)
      a <- as.character(tmpAnm$GroupType[tmpAnm$pltGrp == txNm[w]])[1]
      # se a regressao eh significafiva (p < 0.05)
      if (lmMetrics$signif[lmMetrics$signif$GroupType == a,i-5] < 0.05){
        temp <- !is.na(tmpAnm[tmpAnm$pltGrp == txNm[w],i])
        points(orMt[,i][temp]~tmpAnm[tmpAnm$pltGrp == txNm[w],i][temp], 
               pch = simbo[w-1], col = cores[w-1])
        abline(lmMetrics$intersept[lmMetrics$intersept$GroupType == a,i-5], 
               lmMetrics$angle[lmMetrics$angle$GroupType == a,i-5],
               col = cores[w-1]) # cores das linhas
      }else{
        # se a regressao NAO eh significafiva (p > 0.05)
        temp <- !is.na(tmpAnm[tmpAnm$pltGrp == txNm[w],i])
        points(orMt[,i][temp]~tmpAnm[tmpAnm$pltGrp == txNm[w],i][temp], 
               pch = simbo[w-1], col = cores[w-1])
        abline(lmMetrics$intersept[lmMetrics$intersept$GroupType == a,i-5], 
               lmMetrics$angle[lmMetrics$angle$GroupType == a,i-5],
               col = cores[w-1], lty = 3) # cores das linha e tipo da linha
      }
    }
    i <- length(txNm)
  } 
  # se sao menos que ori + 4 taxonomias (menos que 5)
  if(i < 5){
    for(j in 2:i){
      plot.new()
    }
  }
  plot.new()
  legend("topleft", bty= "n", # posicao e sem borda
         lwd= 1, seg.len= 2,  # espessura e comprimento
         lty= c(2,1,3), # tipos do traco das linhas
         title = "Line type",
         legend= c("original metric values",
                   "p-value < 0.05",
                   "p-value > 0.05"), 
         col=   c("black",
                  "gray42",
                  "gray42"))
  plot.new()
  legend("topleft", bty= "n", # posicao e sem borda
         lwd= 1, seg.len= 2,  # espessura e comprimento
         lty= rep(1,5), # tipo do traco das linhas
         pch= simbo,
         title = "Plants grouped by",
         legend= c("specie", 
                   "genus",
                   "family", 
                   "order"), 
         col= c(cores[1], 
                cores[2],
                cores[3], 
                cores[4]))
  par(mfcol= c(1,1), mar= c(5,4,4,2)+0.1, mgp= c(3,1,0))
  # se foi criado um arquivo, finalizar o arquivo
  if(frmto == "svg" | frmto == "tiff"){
    dev.off()
  }
}



###############################################################################
# Grafico com boxplot com dados da regressao (tons de cinza) e com outline    #
# 1 - representacao dos dados agrupados e dos dados originais                 #
###############################################################################
## 1 - graficos originais e agrupados
grfBxPlt <- function(metrics.lm, orMt, gpMt, metric, title.metric = F, 
                     axisX.name = F){
  txNmPr <- levels(gpMt$anmGrp)
  txNmPr <- c(txNmPr[1], txNmPr[is.element(txNmPr, as.character(unique(gpMt$anmGrp) ) )])
    # ordena os dados das redes originais
  orInf <- orMt[order(orMt$pltGrp),]
  orInf <- orInf[order(orInf$anmGrp),]
  # ordena os dados das redes agrupadas
  gpInf <- gpMt[order(gpMt$pltGrp),]
  gpInf <- gpInf[order(gpInf$anmGrp),]
  # dados de p-value e r2
  pVlue <- metrics.lm$signif[order(metrics.lm$signif$pltGrp),]
  pVlue <- pVlue[order(pVlue$anmGrp),]
  adjR2 <- metrics.lm$adj.r.squared[order(metrics.lm$adj.r.squared$pltGrp),]
  adjR2 <- adjR2[order(adjR2$anmGrp),]
  

  # margens
  par(mar= c(3,3,1.5,0.5)+0.1, mgp= c(2,0.5,0))
  
  w <- 1
  i <- (1:ncol(orMt))[colnames(orMt) == metric] # metrica foco
  for(j in 2:length(txNmPr) ){
    # seleciona animais agrupamento ori txNm[1] e o agrupamento "j" (txNm[j])
    temp <- rbind.data.frame( orInf, gpInf[gpInf$anmGrp == txNmPr[j],] )
    # cores representado dados da regressao
    a <- rep("#ffffff", 4)
    a[ pVlue[(1:4)+4*(j-2), i-5] < 0.05 ] <- "#7f7e7e"
    a[ adjR2[(1:4)+4*(j-2), i-5] > 0.7 ] <- "#1f1f1f"
    a <- c("#000000",a)
    # boxplot
    boxplot(temp[,i]~temp$pltGrp,col= a,  
            las= 1, border="#404040", outline = T)
    # tipo do agrupamento de plantas (cada boxplot - eixo x)
    title(xlab = paste0("Type of plant grouping"))
    # tipo do agrupamento de animais (cada plot)
    title(paste0("Insects grouped by ", txNmPr[j] ), 
          adj=0.8, # posicao
          cex.main = 1, # tamanho em rel ao cex
          font.main = 1) # 1: normal; 2: bold; 3: italic; 4: bold italic 
    if(j == 2 & title.metric){
      title(ylab = metric)
    }
  }
  if(axisX.name){
    title(xlab = paste0("Type of grouping"), outer = T)
  }
}
# fim das funcoes
###############################################################################
###############################################################################








# -----------------------------------------------------------------------------









#####################################################################
# Leitura das informacoes das redes
## leitura dos dados de redes mutualisticas
mutu_tax <- read.csv("../output/mutualism_allMetricsAllGps.csv")[,-1] # so tem 

# colunas para analise dos agrupamentos
mutu_tax <- mutu_tax[, -3] # tira a col tipo (3)
mutu_tax <- adColAgrup(mutu_tax)

# separa original das agrupadas por taxonomia
mutu_ori <- mutu_tax[mutu_tax$GroupType == "Original",] # so info das originais
mutu_tax <- mutu_tax[mutu_tax$GroupType != "Original",] # so dos agrup taxo
rownames(mutu_ori) <- mutu_ori$ArqName



## leitura dos dados de redes antagonicas
anta_tax <- read.csv("../output/antagonism_allMetricsAllGps.csv")[,-1] # so tem 
anta_tax <- anta_tax[-grep("herb-hosts",anta_tax[,1]), ]
# colunas para analise dos agrupamentos
anta_tax <- anta_tax[, -3] # tira a col tipo (3)
anta_tax <- adColAgrup(anta_tax)

# separa original das agrupadas por taxonomia
anta_ori <- anta_tax[anta_tax$GroupType == "Original",] # so info das originais
anta_tax <- anta_tax[anta_tax$GroupType != "Original",] # so dos agrup taxo
rownames(anta_ori) <- anta_ori$ArqName

# metricas presentes
metrics <- colnames(anta_ori)[5:20]


###############################################################################
# calcula a regressao
mutu_lstTx.lm <- lstReg(orMt = mutu_ori, txGrp = mutu_tax)

anta_lstTx.lm <- lstReg(orMt = anta_ori, txGrp = anta_tax)

# -----------------------------------------------------------------------------




###############################################################################
# Grafico de boxplot dos dados originais x agrupados                          #
###############################################################################
## metricas foco
if (frmto == "tiff") {
  tiff(paste0("../output/_bxplt-oriXtx_", "mutu_regressionGraph_",
              format(Sys.time(), "_%Y-%m-d%d"),"-",
              gsub(":", "", format(Sys.time(), "%X") ), ".tif"),
       width = 174,height = 87,units = "mm",res = 300)
}

if (frmto == "svg") {
  svg(paste0("../output/_bxplt-oriXtx_", "mutu_regressionGraph_",
             format(Sys.time(), "_%Y-%m-d%d"),"-",
             gsub(":", "", format(Sys.time(), "%X") ), ".svg"),
      width = 7, height = 4.5, pointsize = 14) # 18 equivale a fonte arial 12
}
par(mfrow = c(3,4), oma = c(0, 0, 1, 0))

## mutualistic
### connectance (5)
grfBxPlt(mutu_lstTx.lm, mutu_ori, mutu_tax, metrics[5], title.metric = T)
### NODF (8)
grfBxPlt(mutu_lstTx.lm, mutu_ori, mutu_tax, metrics[8], title.metric = T)
### Modularity (16)
grfBxPlt(mutu_lstTx.lm, mutu_ori, mutu_tax, metrics[16], title.metric = T)
title(main = "Mutualistic networks", outer = T)

dev.off()



if (frmto == "tiff") {
  tiff(paste0("../output/_bxplt-oriXtx_", "anta_regressionGraph_",
              format(Sys.time(), "_%Y-%m-d%d"),"-",
              gsub(":", "", format(Sys.time(), "%X") ), ".tif"),
       width = 174,height = 87,units = "mm",res = 300)
}

if (frmto == "svg") {
  svg(paste0("../output/_bxplt-oriXtx_", "anta_regressionGraph_",
             format(Sys.time(), "_%Y-%m-d%d"),"-",
             gsub(":", "", format(Sys.time(), "%X") ), ".svg"),
      width = 7, height = 4.5, pointsize = 14) # 18 equivale a fonte arial 12
}
par(mfrow = c(3,4), oma = c(0, 0, 1, 0))
## Antagonistic
### connectance (5)
grfBxPlt(anta_lstTx.lm, anta_ori, anta_tax, metrics[5], title.metric = T, 
         axisX.name = T)
### NODF (8)
grfBxPlt(anta_lstTx.lm, anta_ori, anta_tax, metrics[8], title.metric = T, 
         axisX.name = T)
### Modularity (16)
grfBxPlt(anta_lstTx.lm, anta_ori, anta_tax, metrics[16], title.metric = T, 
         axisX.name = T)
title(main = "Antagonistic networks", outer = T)

dev.off()


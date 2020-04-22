###############################################################################
# Regressao linear                                                            #
# 1 - como os dados agrupados representam os originais                        #
# metricas originais ~ metricas agrupadas                                     #
###############################################################################

# define diretorio de trabalho
#setwd("<your_work_diretory_path>")


###############################################################################
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
# grava em arquivo os dados da regressao                                      #
###############################################################################
# salva os parametros da regressao linear em um arquivo csv
# precisa:
#  - da lista retornada pela funcao lstReg
#  - se e mutualism ou antagonism
#  - tipo do agrupamento
svReg <- function(lmMetrics, webType, agrp = "ori-tx"){
  rgrs <- lmMetrics[[1]]
  rgrs <- cbind(rep(names(lmMetrics)[1], length(rgrs[,1])), rgrs)
  colnames(rgrs)[1] <- 'X'
  for (i in 2:6){
    a <- lmMetrics[[i]]
    a <- cbind(rep(names(lmMetrics)[i], length(a[,1])), a)
    colnames(a)[1] <- 'X'
    rgrs <- rbind(rgrs, a)
  }
  write.csv(rgrs, paste0("../output/", webType, "Regression_",
                         agrp, format(Sys.time(), "_%Y-%m-d%d"),".csv"))
}
###############################################################################
# Cria e salva um grafico em svg para cada metrica                            #
###############################################################################
grfLm <- function(lmMetrics, orMt,txGrp, metric,
                  txNm, webType, fourColors, agrp = "ori-tx", frmto = "tela",
                  mostrarTitle = F, mostrarLegend = F){
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
  if(frmto != "tela"){
    # from https://www.statmethods.net/advgraphs/layout.html
    layout(matrix(c(1,2,3,4, # fig 1, 2, 3 e 4
                    5,5,6,6), # fig 5 e 6 (cd uma ocupando "2 espacos")
                  nrow = 2, ncol = 4, byrow = TRUE))
  }
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
    if(mostrarTitle == T){
      title(paste0("Animals grouped by ", txNm[j]), 
            adj=0.8, # posicao
            cex.main = 1, # tamanho em rel ao cex
            font.main = 1) # 1: normal; 2: bold; 3: italic; 4: bold italic 
      # nome da metrica no eixo y e que sao os valores do agrupa original (y)
    }
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
  
  if(mostrarLegend == T){
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
  } else {
    if (frmto == "svg") {
      plot.new()
      plot.new()
    }
  }
  if (frmto == "svg") {
    par(mfcol= c(1,1), mar= c(5,4,4,2)+0.1, mgp= c(3,1,0))
  }
  # se foi criado um arquivo, finalizar o arquivo
  if(frmto == "svg" | frmto == "tiff"){
    dev.off()
  }
}



###############################################################################
# Grafico de barras com a correspondencia                                     #
###############################################################################
# ate 4 metricas por grafico?                                                 #
# Regressao linear - r2 (r quadrado)  ajustado
grfBrR2 <- function(metrics.lm, metrics, webType, fourColors,
                    agrp = "Ori-GpTax", frmto = "tela"){
  # data.frame base para os graficos
  temp <- metrics.lm$adj.r.squared[ order(metrics.lm$adj.r.squared$pltGrp),]
  temp <- temp[ order(temp$anmGrp),-c(1:2) ]
  a <- metrics.lm$signif[ order(metrics.lm$adj.r.squared$pltGrp),]
  a <- a[ order(a$anmGrp),-c(1:2) ]
  # p-value < 0.005 como significativo (zerei os r2 com p-value maior
  # para nao plotar no grafico)
  temp[a > 0.05] <- 0
  rownames(temp) <- temp$GroupType
  # coloca os agrupamentos nas colunas e as metricas na linhas
  temp <- t(temp[, -1]) # tira a coluna GroupType
  
  ### Formato do Grafico
  if (frmto == "tiff") {
    tiff(paste0("../output/", webType, "_brpltR2Graph_",
                paste0(substr(gsub("\\.", "", metrics),1,4), collapse = "-"), # nome das metricas
                "_",agrp,
                format(Sys.time(), "_%Y-%m-d%d-"),
                gsub(":", "", format(Sys.time(), "%X") ), ".tif"),
         width = 174,height = 55,units = "mm",res = 300)
  }
  
  if (frmto == "svg") {
    svg(paste0("../output/", webType, "_brpltR2Graph_",
               paste0(substr(gsub("\\.", "", metrics),1,4), collapse = "-"), # nome das metricas
               "_",agrp,
               format(Sys.time(), "_%Y-%m-d%d-"),
               gsub(":", "", format(Sys.time(), "%X") ), ".svg"),
        width = 7, height = 3, pointsize = 14) # 18 equivale a fonte arial 12
  }
  
  
  # from https://www.statmethods.net/advgraphs/layout.html
  layout(matrix(c(1,2,3,4, # fig 1, 2, 3 e 4
                  5,5,6,6), # fig 5 e 6 (cd uma ocupando "2 espacos")
                nrow = 2, ncol = 4, byrow = TRUE))
  par(mar= c(3,3,1.5,0.5)+0.1, mgp= c(2,0.5,0))
  # regressao: r quadrado
  for(i in 0:(sqrt(ncol(temp))-1) ){ # para cada agrupamento de animais
    a <- temp[metrics, # metricas a serem plotadas
              (1:sqrt(ncol(temp)))+sqrt(ncol(temp))*i] # as 4 combinacoes por animal spxsp spxge spxfa spxor
    mini <- min(a)
    mini <- ifelse(mini > 0, 0, mini)
    # nome das colunas de "a" (plantas) 
    colnames(a) <- sub(".........","", colnames(a) )
    lines(x = barplot(a, beside = T, 
                      col = rep(c(fourColors, "#ef0059"),3)[1:length(metrics)],
                      ylim= c(mini,1), border = NA), 
          y = rep(0.7, sqrt(ncol(temp))*length(metrics)),lty="dashed")
    axis(1,at= seq(0.5, sqrt(ncol(temp))*(length(metrics)+1) +0.5, (length(metrics)+1) ), # posicao de inicio, passos, fim
         labels= rep("", (sqrt(ncol(temp))+1) )) # qtdade de ticks
    lines(x = c(0,sqrt(ncol(temp))*6), y = c(0,0) )
    title(xlab = paste0("Type of plant grouping"))
    if(i == 0){
      title(ylab = "adj. r squared")
    }
    # tipo do agrupamento de animais
    title(paste0("Animals grouped by ", colnames(a)[i+1] ), 
          adj=0.8, # posicao
          cex.main = 1, # tamanho em rel ao cex
          font.main = 1) # 1: normal; 2: bold; 3: italic; 4: bold italic 
  }
  # se sao menos que as 4 classificacoes taxonomicas
  i <- sqrt(ncol(temp))
  if(i < 4 ){
    for(j in i:3){
      plot.new()
    }
  }
  par( mar= c(1,1,1,0.5) )
  plot.new()
  legend("topright",
         bty = "n", fill= c(fourColors, "#ef0059"),
         legend= metrics, border = NA)
  plot.new()
  par(mfcol= c(1,1), mar= c(5,4,4,2)+0.1, mgp= c(3,1,0))
  # se foi criado um arquivo, finalizar o arquivo
  if(frmto == "svg" | frmto == "tiff"){
    dev.off()
  }
}


###############################################################################
# Tira a media das 100 redes aleatorias equivalentes a cada agrupamento       #
###############################################################################
meanRand <- function(rndall, metrics){
  # rndall <- patGrp
  # cria o data.frame que vai receber as medias com base no tapply
  temp <- tapply(rndall[, metrics[1] ], list(rndall$ArqName, rndall$GroupType), mean)
  rndMean <- cbind.data.frame(rep(rownames(temp),each = ncol(temp)), # "ArqName" 
                              rep(colnames(temp),time = nrow(temp))) # "GroupType"
  colnames(rndMean) <- colnames(rndall)[1:2] # cols "ArqName" e "GroupType"
  temp2 <- temp[1,]
  for ( j in 2:nrow(temp) ){
    temp2 <- c(temp2, temp[j,])
  }
  rndMean <- cbind.data.frame(rndMean, temp2)
  
  # preenche o data.frame
  for(i in 2:length(metrics)){
    temp <- tapply(rndall[, metrics[i] ], list(rndall$ArqName, rndall$GroupType), 
                   function(x) mean(x, na.rm = T) ) # tira a media desconsiderando NA's
    temp2 <- temp[1,]
    for ( j in 2:nrow(temp) ){
      temp2 <- c(temp2, temp[j,])
    }
    rndMean <- cbind.data.frame(rndMean, temp2)
  }
  
  colnames(rndMean)[3:ncol(rndMean)] <- metrics # 1 metrica por coluna
  # ordena
  rndMean <- rndMean[order(rndMean$ArqName),]
  # colunas para analise dos agrupamentos
  pltGrp <- sub(".........","", rndMean$GroupType)
  pltGrp <- sub("ginal","", pltGrp)
  anmGrp <- sub("A_","", rndMean$GroupType)
  anmGrp <- sub(".x.P_..","", anmGrp)
  anmGrp <- sub("ginal","", anmGrp)
  rndMean <- cbind(anmGrp,pltGrp,rndMean)
  txNm <- c('Ori','sp','ge','fa','or')
  rndMean$pltGrp <- factor(as.character(rndMean$pltGrp),
                           levels = txNm)
  rndMean$anmGrp <- factor(as.character(rndMean$anmGrp),
                           levels = txNm)
  return(rndMean)
}


###############################################################################
# Deixa a matriz rand equivalente a matriz de agrupamento taxonomico          #
###############################################################################
# n ori = n rnd por agrupamento (so e necessario para o par a par)
padroRndGrp <- function(RndGrp, txData){
  # so os agrupamentos que foram aleatorizados
  txData <- txData[is.element(as.character(txData$GroupType),
                              as.character(unique(RndGrp$GroupType))) , ]
  # colunas que estao faltando
  clFaltantes <- txData[ , c(3:8)]
  # juntar
  RndGrp <- merge(RndGrp,clFaltantes)
  # ordena as colunas
  RndGrp <- RndGrp[, colnames(txData)]
  # se falta linhas (q sao = ao grupamento tx), acrescentar
  setRnd <- paste(RndGrp$ArqName, RndGrp$GroupType)
  setAll <- paste(txData$ArqName, txData$GroupType)
  arqFalta <- setdiff(setAll, setRnd)
  if (length(arqFalta) > 0) {
    RndGrp <- rbind.data.frame(RndGrp, txData[is.element(setAll, arqFalta),])
  }
  # ordena pelos arquivos
  RndGrp <- RndGrp[order(RndGrp$ArqName),]
  return(RndGrp)
}
# fim das funcoes
###############################################################################
###############################################################################








# -----------------------------------------------------------------------------











###############################################################################
# variaveis a serem utilizadas
webType <- "antagonism" #"mutualism"
fourColors <- c("#131313","#ff6f00","#0018ce","#59ab00")
# c("#ef0059","#ff6f00","#00d1c1","#7ff300")
# c("#000000","#606060","#959595","#d0d0d0")



###############################################################################
# Leitura dos dados
# carrega matrizes com o calculo das metricas
txGrp <- read.csv(paste0('../output/', webType, '_allMetricsAllGps.csv'))[,-1]
numMt <- sum(txGrp$GroupType == "Original")
# salva numa var o nome dos arquivos das matrizes
# ordena a matriz
txGrp <- txGrp[order(txGrp$ArqName),]

# carrega o calculo das metricas para redes aleatorias Patefield
patGrp <- read.csv(paste0('../output/', webType, 
                          '_v02_infoRandPatefield+modul_1-47_2019-10-d14.csv'))[,-1]
# remove os dados do HOSTS
txGrp <- txGrp[-grep("herb-hosts",txGrp[,1]), ]
patGrp <- patGrp[-grep("herb-hosts",patGrp[,1]), ]
# atualiza factors
temp <- as.character(txGrp$ArqName)
txGrp$ArqName <- factor(temp, levels = unique(temp))
temp <- as.character(patGrp$ArqName)
patGrp$ArqName <- factor(temp, levels = unique(temp))
# salva o nome, em ordem
nArqs <- txGrp$ArqName

###############################################################################
# colunas para analise dos agrupamentos
pltGrp <- sub(".........","", txGrp$GroupType)
pltGrp <- sub("ginal","", pltGrp)
anmGrp <- sub("A_","", txGrp$GroupType)
anmGrp <- sub(".x.P_..","", anmGrp)
anmGrp <- sub("ginal","", anmGrp)
txGrp <- cbind(anmGrp,pltGrp,txGrp[,-3]) # tira a col tipo (3)
txNm <- c('Ori','sp','ge','fa','or')
txGrp$pltGrp <- factor(as.character(txGrp$pltGrp),
                       levels = txNm)
txGrp$anmGrp <- factor(as.character(txGrp$anmGrp),
                       levels = txNm)


###############################################################################
# nomes das metricas
metrics <- colnames(txGrp)[(1:ncol(txGrp))[colnames(txGrp)=="connectance"]:ncol(txGrp)]




# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

###############################################################################
# Permite comparar como o original difere do sp X sp)
# separa a matriz original das agrupadas taxonomicamente
bckup <- txGrp[txGrp$GroupType != "Original",] # TIRA OS ORIGINAIS

orMt <- txGrp[txGrp$GroupType == "Original",] # so info das originais
#orMt <- bckup[bckup$GroupType == "A_sp.x.P_sp",] # so info das "originais"
#txGrp[bckup$GroupType == "A_sp.x.P_sp", 9:ncol(bckup) ] <- NA
txGrp <- bckup

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++





###############################################################################
# arruma matriz de dados aleatorios
# meanRand <- function(rndall, metrics){
patGrp.OK <- meanRand(patGrp, metrics)
patGrp.OK <- padroRndGrp(patGrp.OK, txGrp)


###############################################################################
# calcula a regressao
lstTx.lm <- lstReg(orMt = orMt, txGrp = txGrp)
lstPat.lm <- lstReg(orMt = orMt, txGrp = patGrp.OK)


###############################################################################
# grava em arquivo os dados da regressao
# svReg <- function(lmMetrics, webType, agrp = "ori-tx")
# svReg(lstTx.lm,webType, agrp = "ori-tx")
# svReg(lstPat.lm,webType, agrp = "ori-Patefield")


# -----------------------------------------------------------------------------

###############################################################################
# Grafico da Regressao linear                                                 #
###############################################################################
# plota os pontos e retas da regressao
# grfLm <- function(lmMetrics, orMt,txGrp, metric,
#                   txNm, webType, fourColors, agrp = "ori-tx"){
# # Original X Taxonomico
# 3 metricas
tiff(paste0("../output/", webType, "_regressionGraph_3metrics_", "Ori-GpTax",
            format(Sys.time(), "_%Y-%m-d%d-"),
            gsub(":", "", format(Sys.time(), "%X") ), ".tif"),
     width = 174,height = 75,units = "mm",res = 300)
par(mfrow= c(3,4), mar= c(3,3,1.5,0.5)+0.1, mgp= c(2,0.5,0),oma = c(0, 0, 1, 0))
mostrarTitle = T
for(i in c(1, 4, 12) ){
  mostrarTitle <- ifelse(i == 1, T, F)
  grfLm(lstTx.lm,orMt,txGrp, metrics[i],txNm, webType, fourColors,
                agrp = "Ori-GpTax", mostrarTitle = mostrarTitle)
  title(main = paste0(sub("ism","istic",webType)," networks"), outer = T)
}
dev.off()
# globais
tiff(paste0("../output/", webType, "_regressionGraph_globais_", "Ori-GpTax",
            format(Sys.time(), "_%Y-%m-d%d-"),
            gsub(":", "", format(Sys.time(), "%X") ), ".tif"),
     width = 174,height = 125,units = "mm",res = 300)
par(mfrow= c(5,4), mar= c(3,3,1.5,0.5)+0.1, mgp= c(2,0.5,0),oma = c(0, 0, 1, 0))
mostrarTitle = T
for(i in c(3, 11, 5, 2, 6) ){
  mostrarTitle <- ifelse(i == 3, T, F)
  grfLm(lstTx.lm,orMt,txGrp, metrics[i],txNm, webType, fourColors,
        agrp = "Ori-GpTax", mostrarTitle = mostrarTitle)
  title(main = paste0(sub("ism","istic",webType)," networks"), outer = T)
}
dev.off()
# level
tiff(paste0("../output/", webType, "_regressionGraph_level_", "Ori-GpTax",
            format(Sys.time(), "_%Y-%m-d%d-"),
            gsub(":", "", format(Sys.time(), "%X") ), ".tif"),
     width = 174,height = 100,units = "mm",res = 300)
par(mfrow= c(4,4), mar= c(3,3,1.5,0.5)+0.1, mgp= c(2,0.5,0),oma = c(0, 0, 1, 0))
mostrarTitle = T
for(i in c(7, 8, 9, 10 ) ){
  mostrarTitle <- ifelse(i == 7, T, F)
  grfLm(lstTx.lm,orMt,txGrp, metrics[i],txNm, webType, fourColors,
        agrp = "Ori-GpTax", mostrarTitle = mostrarTitle)
  title(main = paste0(sub("ism","istic",webType)," networks"), outer = T)
}
dev.off()


for(i in 1:12 ){
  grfLm(lstTx.lm,orMt,txGrp, metrics[i],txNm, webType, fourColors,
        agrp = "Ori-GpTax", frmto = "tiff")
}
grfLm(lstTx.lm,orMt,txGrp, metrics[i],txNm, webType, fourColors,
      agrp = "Ori-GpTax")

# Original X Patefield
for(i in 1:length(metrics) ){
  grfLm(lstPat.lm,orMt,patGrp.OK, metrics[i],txNm[1:3], webType, fourColors, 
        agrp = "Ori-GpPatefield", frmto = "tiff")
  
}



###############################################################################
# Grafico de barras com a correspondencia                                     #
# ate 4 metricas por grafico                                                  #
# Regressao linear - r2 (r quadrado)  ajustado                                #
###############################################################################
# grfBrR2 <- function(metrics.lm, metrics, webType, fourColors,
#                     agrp = "Ori-GpTax", frmto = "tela"){

## DA UM WARNING, MAS TA TUDO BEM!!!!!!
# Warning message:
# In Ops.factor(left, right) : '>' not meaningful for factors

### Original X Taxonomico
# "connectance" "NODF" "Modularity"
grfBrR2(lstTx.lm, metrics[c(1,4,12)], webType, fourColors,
        agrp = "Ori-GpTax")
grfBrR2(lstTx.lm, metrics[c(1,4,12)], webType, fourColors, 
        agrp = "Ori-GpTax", frmto = "tiff")


### Original X Patefield
# "connectance" "NODF" "Modularity"
grfBrR2(lstPat.lm, metrics[c(1,4,12)], webType, fourColors,
        agrp = "Ori-GpPatefield")
grfBrR2(lstPat.lm, metrics[c(1,4,12)], webType, fourColors, 
        agrp = "Ori-GpPatefield", frmto = "tiff")

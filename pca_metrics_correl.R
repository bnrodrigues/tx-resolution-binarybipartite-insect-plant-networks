###############################################################################
# Calculo do PCA para a analise do agrupamento                                #
###############################################################################

# define diretorio de trabalho
#setwd("<your_work_diretory_path>")


###############################################################################
# Funcao calculo da Regressao linear                                          #
###############################################################################
# valores das metricas nas redes originais X redes agrupadas
# precisa dos valores para as matrizes originais e das agrupadas
lstRegPC1 <- function(or.pca.pc1, dfGrp){
  # lista de data frames para armazenar os valores da regressao
  # cada item da lista
  temp <- as.data.frame(unique(dfGrp[,c(1:2,4)])) # tipos de agrupamento (ex.: 16)
  tam <- length(temp[,1]) # qtdade de agrupamentos (ex.: 16)
  # sera preenchido com os dados de regressao do PC1
  temp <- cbind.data.frame(temp, vector('numeric', length = tam))
  
  ### regressao linear
  # parametros da regressao para analise + modelo de regressao + PCA
  metrics.lm <- vector('list', length = 8) 
  # 1- sigma: Residual standard error
  # 2- adj.r.squared: r quadrado ajustado
  # 3- a$df[2]: graus de liberdade
  # 4- coefficients[2,4]: significancia
  # 5- coefficients[1,1]: intersept (b)
  # 6- coefficients[2,1]: angulo da reta (a)
  # 7- lm
  # 8- pca
  names(metrics.lm) <- c('sigma','adj.r.squared','f.degree',
                         'signif','intersept','angle','lm','PCAg')
  # data.frame para os parametros das regressoes
  for(i in 1:6){
    metrics.lm[[i]] <- temp # temp = data.frame vazio (so com info de agrupamento)
  }
  # lista de regressoes por agrupamento
  metrics.lm[[7]] <- vector('list', length = tam)
  names(metrics.lm[[7]]) <- as.character(temp$GroupType)
  # lista de PCA por agrupamento
  metrics.lm[[8]] <- vector('list', length = tam)
  names(metrics.lm[[8]]) <- as.character(temp$GroupType)
  # calculo dos modelos
  for(j in 1:tam){ # cd agrupamento
    # PCA do agrupado
    sb.dfGrp <- dfGrp[dfGrp$GroupType == as.character(metrics.lm[[1]]$GroupType[j]),]
    metrics.lm[[8]][[j]] <- prcomp(~ numSpecies + connectance + NODF + Modularity + 
                                     links.per.species + nestedness + web.asymmetry + 
                                     robustness.HL + robustness.LL +
                                     niche.overlap.HL + niche.overlap.LL + N.modules, 
                                   data = sb.dfGrp,
                                   scale = TRUE)
    ### regresssao lm: PC1 originais ~ PC1 agrupado[j]
    ## iguala o comprimento dos vetores
    temp.na <- vector("numeric", length = length(or.pca.pc1) )
    # se uma das colunas tem NA o pca nao e calculado para a coluna e nao aparece no pc1
    fco.na <- is.na( apply(sb.dfGrp[,c("connectance","NODF","Modularity",
                                       "links.per.species","nestedness",
                                       "web.asymmetry","robustness.HL",
                                       "robustness.LL","niche.overlap.HL",
                                       "niche.overlap.LL","N.modules")],1,sum) )
    temp.na[ (1:length(or.pca.pc1))[fco.na] ] <- NA
    temp.na[ (1:length(or.pca.pc1))[!fco.na] ] <- metrics.lm[[8]][[j]]$x[,1]
    ## faz a regressao
    metrics.lm[[7]][[j]] <- lm(or.pca.pc1 ~ # valores ori
                                 temp.na,# valores por agrupamento
                               na.action = na.omit)
    # dados da regressao
    a <- summary(metrics.lm[[7]][[j]])
    metrics.lm$sigma[j,4] <- a$sigma # Residual standard error
    metrics.lm$adj.r.squared[j,4] <- a$adj.r.squared # r quadrado ajustado
    metrics.lm$f.degree[j,4] <- a$df[2] # graus de liberdade
    metrics.lm$signif[j,4] <- a$coefficients[2,4] # significancia
    metrics.lm$intersept[j,4] <- a$coefficients[1,1] # intersept (b)
    metrics.lm$angle[j,4] <- a$coefficients[2,1] # angulo da reta (a)
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
svReg <- function(lmMetrics, webType, agrp){#agrp = "ori-tx"){
  rgrs <- lmMetrics[[1]]
  rgrs <- cbind(rep(names(lmMetrics)[1], length(rgrs[,1])), rgrs)
  colnames(rgrs)[1] <- 'X'
  for (i in 2:6){
    a <- lmMetrics[[i]]
    a <- cbind(rep(names(lmMetrics)[i], length(a[,1])), a)
    colnames(a)[1] <- 'X'
    rgrs <- rbind(rgrs, a)
  }
  write.csv(rgrs, paste0("../output/", webType, "Regression_pc1nNodes_AllMetrics_",
                         agrp, format(Sys.time(), "_%Y-%m-d%d"),".csv"))
}
###############################################################################
# Cria e salva um grafico em svg para cada metrica                            #
###############################################################################
grfLmPC1 <- function(gp.pca, or.pca, dfGrp, cores1,
                     webType, agrp = "ori-tx"){
  txNm <- as.character(unique(sort(gp.pca[[1]]$anmGrp)))

  for(j in 1:length(txNm) ){ # agrupamento de animais
    # teste para saber se plota ou nao (ate q class vai sp->ge ou sp->or)
    tmpAnm <- gp.pca[[8]][ gp.pca[[1]]$anmGrp == txNm[j] ] # sel pca
    tmpAnm.gp <- gp.pca[[1]][ gp.pca[[1]]$anmGrp == txNm[j], ] # sel sigma (so para ter o nome dos grupos)
    par(mar=c(3,3.5,2,0.5), mgp= c(1.5,0.3,0))
    # lim do grafico
    mini <- min(c(tmpAnm[tmpAnm.gp$pltGrp == "sp"]$x[,1],
                   tmpAnm[tmpAnm.gp$pltGrp == "ge"]$x[,1],
                   or.pca$x[,1]), na.rm = T)
    maxi <- max(c(tmpAnm[tmpAnm.gp$pltGrp == "sp"]$x[,1],
                  tmpAnm[tmpAnm.gp$pltGrp == "ge"]$x[,1],
                  or.pca$x[,1]), na.rm = T)
    
    
    # pontos e linhas (valore de y [original] em funcao de x [agrupado])
    for(w in 1:length(txNm)){ # para sp, ge, fa e or do grupamento de plantas
      # tipo do agrupamento (spXsp, spXge, geXsp...)
      a <- as.character(tmpAnm.gp$GroupType[tmpAnm.gp$pltGrp == txNm[w]])[1]
      # separa os arquivos que tem pca
      arqs <- as.character(dfGrp[rownames(tmpAnm[a][[1]]$x), "ArqName"])
      # plot "vazio" (so com a linha dos valores originais)
      plot(or.pca$x[arqs,1] ~ tmpAnm[tmpAnm.gp$pltGrp == txNm[w]][[1]]$x[,1],
           ylim= c(mini,maxi), 
           xlim= c(mini,maxi),
           type = "n", bty= "l", ylab = "", xlab = "")
      abline(0,1, lty = 2)
      # tipo do agrupamento de animais
      title(paste0("Insects grouped by ", txNm[j], " and\nplants grouped by ", txNm[w]), 
            adj=0.8, # posicao horizontal
            cex.main = 1, # tamanho em rel ao cex
            font.main = 1) # 1: normal; 2: bold; 3: italic; 4: bold italic 
      # nome da metrica no eixo y e que sao os valores do agrupa original (y)
      if(j == 1 & w == 1){
        title(ylab = paste0("PC1: ",# % de var explicada pelo pc1
                            round(or.pca$sdev^2/sum(or.pca$sdev^2)*100, 1)[1],
                            "%","\nOriginal PC1 values"), line = 1.3  ) 
      }

      # se a regressao eh significafiva (p < 0.05)
      if (gp.pca$signif[gp.pca$signif$GroupType == a,4] < 0.05){
        # se uma das colunas tem NA o pca nao e calculado para a coluna e nao aparece no pc1
        cat(paste0("- ",a,"\n"))
        temp <- !is.na( apply(dfGrp[dfGrp$GroupType == a,
                                    c("connectance","NODF","Modularity",
                                      "links.per.species","nestedness",
                                      "web.asymmetry","robustness.HL",
                                      "robustness.LL","niche.overlap.HL",
                                      "niche.overlap.LL","N.modules")],
                              1,sum) )
        pc1.gp <- tmpAnm[tmpAnm.gp$pltGrp == txNm[w]][[1]]
        points(or.pca$x[arqs,1][temp] ~ pc1.gp$x[,1][temp],
               pch = 16, col = cores1)
        abline(gp.pca[[7]][names(gp.pca[[7]]) == a][[1]],
               col = cores1) # cores das linhas
        title(xlab = paste0("Grouped  PC1 values\n", "PC1: ",# % de var explicada pelo pc1
                            round(pc1.gp$sdev^2/sum(pc1.gp$sdev^2)*100, 1)[1], "%"),
              line = 2  )
      }else{
        # se a regressao NAO eh significafiva (p > 0.05)
        cat(paste0("- no sig ",a,"\n"))
        temp <- !is.na( apply(dfGrp[dfGrp$GroupType == a,
                                    c("connectance","NODF","Modularity",
                                      "links.per.species","nestedness",
                                      "web.asymmetry","robustness.HL",
                                      "robustness.LL","niche.overlap.HL",
                                      "niche.overlap.LL","N.modules")],
                              1,sum) )
        pc1.gp <- tmpAnm[tmpAnm.gp$pltGrp == txNm[w]][[1]]
        points(or.pca$x[arqs,1][temp] ~ pc1.gp$x[,1][temp],
               pch = 16, col = cores1)
        abline(gp.pca[[7]][[names(gp.pca[[7]]) == a]],
               col = cores1, lty= "dotted") # cores das linha e tipo da linha
        title(xlab = paste0("Grouped  PC1 values\n", "PC1: ",# % de var explicada pelo pc1
                            round(pc1.gp$sdev^2/sum(pc1.gp$sdev^2)*100, 1)[1], "%"),
              line = 2  )
      }
    }
    i <- 1+length(txNm)
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
    temp <- tapply(rndall[, metrics[i] ], list(rndall$ArqName, rndall$GroupType), mean)
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
  rndMean <- adColAgrup(rndMean)
  return(rndMean)
}


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
# Deixa a matriz rand equivalente a matriz de agrupamento taxonomico          #
###############################################################################
# n ori = n rnd por agrupamento (so e necessario para o par a par)
padroRndGrp <- function(RndGrp, txData){
  # so os agrupamentos que foram aleatorizados
  txData <- txData[is.element(as.character(txData$GroupType),
                              as.character(unique(RndGrp$GroupType))) , ]
    # ordena as colunas
  RndGrp <- RndGrp[, colnames(txData)]
  RndGrp <- RndGrp[!is.na(RndGrp$numSpecies), ]
  # se falta linhas (q sao = ao grupamento tx), acrescentar
  # setRnd <- paste(RndGrp$ArqName, RndGrp$GroupType)
  # setAll <- paste(txData$ArqName, txData$GroupType)
  # arqFalta <- setdiff(setAll, setRnd)
  # if (length(arqFalta) > 0) {
  #   temp <- txData[is.element(setAll, arqFalta),]
  #   temp[,9:20] <- NA
  #   RndGrp <- rbind.data.frame(RndGrp, temp)
  # }
  
  # se falta linhas (q sao = ao grupamento tx), acrescentar
  setRnd <- paste(RndGrp$ArqName, RndGrp$GroupType)
  setAll <- paste(txData$ArqName, txData$GroupType)
  arqFalta <- setdiff(setAll, setRnd)
  if (length(arqFalta) > 0) {
    RndGrp <- rbind.data.frame(RndGrp, txData[is.element(setAll, arqFalta),])
  }
  
  # ordena pelos arquivos
  # ordena
  RndGrp$ArqName <- factor(as.character(RndGrp$ArqName), 
                           levels = levels(txData$ArqName) )
  RndGrp <- RndGrp[order(RndGrp$ArqName),]
  return(RndGrp)
}

# fim das funcoes
###############################################################################
###############################################################################








# -----------------------------------------------------------------------------









#####################################################################
# Leitura das informacoes das redes
## leitura dos dados de redes mutualisticas
mutu_tax <- read.csv("../output/mutualism_infoOriTxSpGe+modul_2019-10-d01.csv")[,-1] # so tem 
mutu_par <- read.csv("../output/mutualism_v01_infoRandParPar+modul_1-47_2019-10-d01.csv", row.names = NULL)[,-1]
mutu_pat <- read.csv("../output/mutualism_v02_infoRandPatefield+modul_1-47_2019-10-d14.csv", row.names = NULL)[,-1]

# colunas para analise dos agrupamentos
mutu_tax <- mutu_tax[, -3] # tira a col tipo (3)
mutu_tax <- adColAgrup(mutu_tax)

# separa original das agrupadas por taxonomia
mutu_ori <- mutu_tax[mutu_tax$GroupType == "Original",] # so info das originais
mutu_tax <- mutu_tax[mutu_tax$GroupType != "Original",] # so dos agrup taxo
rownames(mutu_ori) <- mutu_ori$ArqName

# metricas presentes
metrics <- colnames(mutu_ori)[5:20]

# tira a media dos calculos das matrizes randomicas
mutu_pat.ok <- meanRand(mutu_pat, metrics)
mutu_pat.ok <- padroRndGrp(mutu_pat.ok, mutu_tax)
mutu_par.ok <- meanRand(mutu_par, metrics)
mutu_par.ok <- padroRndGrp(mutu_par.ok, mutu_tax)



## leitura dos dados de redes antagonicas
anta_tax <- read.csv("../output/antagonism_infoOriTxSpGe+modul_2019-10-d01.csv")[,-1] # so tem 
anta_par <- read.csv("../output/antagonism_v01_infoRandParPar+modul_1-47_2019-10-d01.csv", row.names = NULL)[,-1]
anta_pat <- read.csv("../output/antagonism_v02_infoRandPatefield+modul_1-47_2019-10-d14.csv", row.names = NULL)[,-1]
# remove as redes do HOSTS
anta_tax <- anta_tax[-grep("herb-hosts",anta_tax[,1]), ]
anta_par <- anta_par[-grep("herb-hosts",anta_par[,1]), ]
anta_pat <- anta_pat[-grep("herb-hosts",anta_pat[,1]), ]
# atualiza factors
temp <- as.character(anta_tax$ArqName)
anta_tax$ArqName <- factor(temp, levels = unique(temp))
temp <- as.character(anta_par$ArqName)
anta_par$ArqName <- factor(temp, levels = unique(temp))
temp <- as.character(anta_pat$ArqName)
anta_pat$ArqName <- factor(temp, levels = unique(temp))

# colunas para analise dos agrupamentos
anta_tax <- anta_tax[, -3] # tira a col tipo (3)
anta_tax <- adColAgrup(anta_tax)

# separa original das agrupadas por taxonomia
anta_ori <- anta_tax[anta_tax$GroupType == "Original",] # so info das originais
anta_tax <- anta_tax[anta_tax$GroupType != "Original",] # so dos agrup taxo
rownames(anta_ori) <- anta_ori$ArqName

# metricas presentes
metrics <- colnames(anta_ori)[5:20]

# tira a media dos calculos das matrizes randomicas
anta_pat.ok <- meanRand(anta_pat, metrics)
anta_pat.ok <- padroRndGrp(anta_pat.ok, anta_tax)
anta_par.ok <- meanRand(anta_par, metrics)
anta_par.ok <- padroRndGrp(anta_par.ok, anta_tax)


###############################################################################
# calculo do PCA PC1 para cada agrupamento
mutu_or.pca <- prcomp(~ numSpecies + connectance + NODF + Modularity + 
                        links.per.species + nestedness + web.asymmetry + 
                        robustness.HL + robustness.LL +
                        niche.overlap.HL + niche.overlap.LL + N.modules, 
                      data = mutu_ori,
                      scale = TRUE)
anta_or.pca <- prcomp(~ numSpecies + connectance + NODF + Modularity + 
                        links.per.species + nestedness + web.asymmetry + 
                        robustness.HL + robustness.LL +
                        niche.overlap.HL + niche.overlap.LL + N.modules, 
                      data = anta_ori,
                      scale = TRUE)

###############################################################################
# calculo do PCA PC1 para cada agrupamento
# lstRegPC1 <- function(or.pca.pc1, dfGrp)
mutu_tx.pca <- lstRegPC1(or.pca.pc1 = mutu_or.pca$x[,1], 
                         dfGrp = mutu_tax)
mutu_pt.pca <- lstRegPC1(or.pca.pc1 = mutu_or.pca$x[,1], 
                         dfGrp = mutu_pat.ok)
mutu_pr.pca <- lstRegPC1(or.pca.pc1 = mutu_or.pca$x[,1], 
                         dfGrp = mutu_par.ok)

anta_tx.pca <- lstRegPC1(or.pca.pc1 = anta_or.pca$x[,1], 
                         dfGrp = anta_tax)
anta_pt.pca <- lstRegPC1(or.pca.pc1 = anta_or.pca$x[,1], 
                         dfGrp = anta_pat.ok)
anta_pr.pca <- lstRegPC1(or.pca.pc1 = anta_or.pca$x[,1], 
                         dfGrp = anta_par.ok)


###############################################################################
# Salva os dados da regressao
# function(lmMetrics, webType, agrp)
svReg(mutu_tx.pca, "mutualism", "ori-tx")
svReg(mutu_pr.pca, "mutualism", "ori-pr")
svReg(mutu_pt.pca, "mutualism", "ori-pt")

svReg(anta_tx.pca, "antagonism", "ori-tx")
svReg(anta_pr.pca, "antagonism", "ori-pr")
svReg(anta_pt.pca, "antagonism", "ori-pt")


###############################################################################
# grafico 

if (frmto == "tiff") {
  tiff(paste0("../output/_pc1nNode-AllMetricss_", "_regressionGraph_",
              format(Sys.time(), "_%Y-%m-d%d"),"-",
              gsub(":", "", format(Sys.time(), "%X") ), ".tif"),
       width = 174,height = 200,units = "mm",res = 300)
}

if (frmto == "svg") {
  svg(paste0("../output/_pc1nNodes-AllMetrics_", "_regressionGraph_",
             format(Sys.time(), "_%Y-%m-d%d"),"-",
             gsub(":", "", format(Sys.time(), "%X") ), ".svg"),
      width = 7, height = 9, pointsize = 14) # 18 equivale a fonte arial 12
}

# from https://www.statmethods.net/advgraphs/layout.html
layout(matrix(c(1,2,3,4, # fig 1, 2, 3 e 4 - mutualismo
                5,6,7,8,
                9,10,11,12,
              #13,13,13,13, # legenda
              13,14,15,16, # antagonismo
              17,18,19,20,
              21,22,23,24,
              25,25,25,25
              ), # fig 5 e 6 (cd uma ocupando "2 espacos")
              nrow = 7, ncol = 4, byrow = TRUE))




# fourColors <- c("#131313","#ff6f00","#0018ce","#59ab00")
# grey(seq(0.25,0.75, 0.6/3)) "#404040" "#737373" "#A6A6A6"
grfLmPC1(mutu_tx.pca, mutu_or.pca, mutu_tax, "#404040",
         "mutualism", agrp = "ori-tx")
grfLmPC1(mutu_pt.pca, mutu_or.pca, mutu_pat.ok, "#737373",
         "mutualism", agrp = "ori-pt")
grfLmPC1(mutu_pr.pca, mutu_or.pca, mutu_par.ok, "#A6A6A6",
         "mutualism", agrp = "ori-pr")

# par(mar=c(3,3.5,0.5,0.5))
# plot.new()
# legend("topleft", bty= "n", # posicao e sem borda
#        lwd= 1, seg.len= 2,  # espessura e comprimento
#        lty= rep(1,5), # tipo do traco das linhas
#        pch= 16,
#        title = "Mutualistic",
#        legend= c("Taxonomic aggregation",
#                  "Patefield's model",
#                  "Randomic-paired aggregation"),
#        col= c("#4a7fb0",
#               "#5b9bd5",
#               "#abc6e5"))


grfLmPC1(anta_tx.pca, anta_or.pca, anta_tax, "#404040",
         "antagonism", agrp = "ori-tx")
grfLmPC1(anta_pt.pca, anta_or.pca, anta_pat.ok, "#737373",
         "antagonism", agrp = "ori-pt")
grfLmPC1(anta_pr.pca, anta_or.pca, anta_par.ok, "#A6A6A6",
         "antagonism", agrp = "ori-pr")

par(mar=c(3,3.5,0.5,0.5))
plot.new()
legend("topleft", bty= "n", # posicao e sem borda
       lwd= 1, seg.len= 2,  # espessura e comprimento
       lty= rep(1,5), # tipo do traco das linhas
       pch= 16,
       # title = "Antagonistic",
       legend= c("Taxonomic aggregation",
                 "Patefield's model",
                 "Randomic-paired aggregation"),
       col= c("#404040",
              "#737373",
              "#A6A6A6"))

dev.off()

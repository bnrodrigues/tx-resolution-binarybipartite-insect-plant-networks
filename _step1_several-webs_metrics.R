###############################################################################
################ Redes de interacao planta e visitante floral #################
# Pega todos os arquivos csv (separados por ",") presentes numa pasta e le    #
# para uma lista. Os arquivos devem estar com a plantas nas linhas e animais  #
# nas colunas.                                                                #
# Procura a grafia correta e classificacao taxonomica de cada nome.           #
# Desconsidera nomes que nao sao encontrados em bases de dados confiaveis.    #
# Cria matrizes com os nomes corretos e matrizes por genero e familia.        #
###############################################################################

webType <- "antagonism" #"mutualism"
###############################################################################
# entra no diretorio de trabalho ##############################################
#setwd("<your_work_diretory_path>")

###############################################################################
# carrega as bibliotecas ######################################################

# biblioteca para trabalhar com os nomes cientificos e taxonomia ##############
#Scott Chamberlain and Eduard Szocs (2013). taxize - taxonomic search and 
##retrieval in R. F1000Research, 2:191. 
##URL: http://f1000research.com/articles/2-191/v2. 
#Scott Chamberlain, Eduard Szocs, Carl Boettiger, Karthik Ram, Ignasi Bartomeus, 
##John Baumgartner, Zachary Foster, James O’Donnell, and Jari Oksanen (2017). 
##taxize: Taxonomic information from around the web. R package version 0.9.0. 
##https://github.com/ropensci/taxize
#install.packages("taxize")
library(taxize)

# biblioteca para redes bipartidas e metricas de redes
# To cite bipartite:
#       Dormann, C.F., Gruber B. & Fruend, J. (2008). Introducing the
#       bipartite Package: Analysing Ecological Networks. R news Vol 8/2,
#       8 - 11.
# For network-level analyses:
#       Dormann, C.F., Fruend, J., Bluethgen, N. & Gruber B. 2009.
#       Indices, graphs and null models: analyzing bipartite ecological
#       networks. The Open Ecology Journal, 2, 7-24.
library(bipartite)


###############################################################################
# funcoes a serem utilizadas ##################################################

# tira os nomes que nao tem significado para alguma identificacao
##por exemplo: unidentified, morphotype (no inicio do nome)
#temp[966:972, ] # unidentified
tiraSemSent <- function(vetorSp){
  vetorSpSemSet <- vetorSp
  tam <- length(vetorSp[grep("^[A-z]nidentified", vetorSp)])
  if(tam > 0){
    vetorSpSemSet <- vetorSp[-grep("^[A-z]nidentified", vetorSp)]
    write(c(paste0("It was found ", tam, " unidentified sp"),
            vetorSp[grep("^[A-z]nidentified", vetorSp)]), 
          file = paste0("unidentified", tam,".txt"),
          sep = "\n")
  }
  tam <- length(vetorSp[grep("^[A-z]orphotype", vetorSp)])
  if(tam > 0){
    vetorSpSemSet <- vetorSp[-grep("^[A-z]orphotype", vetorSp)]
    write(c(paste0("It was found ", tam, " mophotypes"),
            vetorSp[grep("^[A-z]orphotype", vetorSp)]), 
          file = paste0("morphotype", tam,".txt"),
          sep = "\n")
  }                  
  return(vetorSpSemSet)
}

# retira todos os nomes que nao tem resolução até especie no vetor com os nomes
#originais
selecSp <- function(nomesOriginais){
  a <- strsplit(nomesOriginais, " ")
  temp1 <- vector()
  temp2 <- vector()
  for(i in 1: length(nomesOriginais)){
    if(length(a[[i]]) < 2){ # se nao tem pelo menos 2 palavras no nome
      temp1 <- append(temp1, a[[i]])
    } else { # se TEM pelo menos 2 palavras no nome
      temp2 <- append(temp2, nomesOriginais[i])
    }
  }
  tam <- length(temp1)
  write(c(paste0("It was found ", tam, " non-species names"),
          temp1), 
        file = paste0("non-speciesNames", tam,".txt"),
        sep = "\n")
  return(temp2)
}

# normalizacao simples dos nomes das especies 
simplNorm <- function(nomesOriginais){
  # tira ponto do final
  simpl <- sub("\\.", "",  nomesOriginais)
  # tira espacos do final
  simpl <- sub(" $", "",  simpl)
  # relaciona o nome corrigido com o original
  simpl <- cbind(nomesOriginais,  simpl)
  # ordena o vetor (pela coluna com os nomes originais)
  simpl <- simpl[order(simpl[, 1]), ]
  return(simpl)
}

# grafia correta dos nomes segundo a funcao gnr_resolve 
#(qt mais nomes, mais demora)
# quando o vetor de nomes esta muito grande a coneccao eh interrompida
# eh necessario separar em conjuntos menores e depois juntar
## (de 250 em 250)
## de 500 em 500 excede o tempo limite do servidor e nao conclui a busca
nomeCorreto <- function(nomeOriginal){
  # bases de dados que retornam a classificacao taxonomica tbm
  #databases <- gnr_datasources()[c(1,4,11),] # cat of life, NCBI, GBIF
  databases <- gnr_datasources()[c(1,11),] # cat of life, GBIF
  # dependendo da base de dados e do nome, as vezes pode nao ter o mesmo num
  #de colunas, nem na mesma ordem 
  tam <- length(nomeOriginal)
  qtpg <- 100
  i <- 1
  temp <- data.frame()
  while ((i+qtpg) < tam){
    gnrTtl <- gnr_resolve(names = nomeOriginal[i:(i+qtpg-1)], 
                          best_match_only = T, 
                          fields = "all", # busca todas as info possiveis
                          data_source_ids = databases[, 1])
    temp <- rbind(temp,
                  data.frame(user_supplied_name = gnrTtl$user_supplied_name,
                             matched_name = gnrTtl$matched_name,
                             classification_path = gnrTtl$classification_path,
                             classification_path_ranks = gnrTtl$classification_path_ranks,
                             data_source_title = gnrTtl$data_source_title,
                             taxon_id = gnrTtl$taxon_id,
                             score = gnrTtl$score)
    )
    cat(paste(" ",i," |"))
    i <- i+qtpg
  }
  cat("done")
  gnrTtl <- gnr_resolve(names = nomeOriginal[i:tam], 
                        best_match_only = T, 
                        fields = "all", # busca todas as info possiveis
                        data_source_ids = databases[, 1])
  temp <- rbind(temp,
                data.frame(user_supplied_name = gnrTtl$user_supplied_name,
                           matched_name = gnrTtl$matched_name,
                           classification_path = gnrTtl$classification_path,
                           classification_path_ranks = gnrTtl$classification_path_ranks,
                           data_source_title = gnrTtl$data_source_title,
                           taxon_id = gnrTtl$taxon_id,
                           score = gnrTtl$score)
  )
  return(temp)
}

# Arruma os nomes de sp retornados pelo gnr_resolve
## deixa no formato genero + epiteto especifico ou so um (nome do gen, fam...)
# normaliza o retorno das bases de dados
# site https://regexr.com/ ajudou a pensar na expressao regular
#a <- normSpAni[c(66,46,32,33,71,86)] # para testar
transNome <- function(nomesOriginais){
  # tira caracteres entre parenteses
  a1 <-  sub(" \\([A-Za-z0-9]+\\)", "", nomesOriginais)
  # posicao de onde esta o nome da especie (genero + epiteto sp)
  pscao <- regexpr("^([A-Za-z]+ ?[a-z]*)", a1)
  a2 <- vector("character", length = length(a1))
  # separo so o nome da sp (genero + epiteto sp)
  for (i in 1:length(a1)){
    # assumindo q o nome comeca sempre com o nome do genero
    a2[i] <- substring(a1[i], 1, attr(pscao,"match.length")[i])
  }
  a2 <- sub(" $",  "", a2)
  return(a2)
}

# constroi uma matriz com as classificacoes taxonomicas disponiveis
# utiliza a saida da funcao gnr_resolve com o parametro fields= "all"
classifTx <- function(gnrOutput_alterado){
  tam <- length(gnrOutput_alterado[,1])
  nomeCol <- c("original_name", "matched_name", "species", "subgenus", 
               "genus", "subfamily", "family", "suborder", "order",
               "class", "phylum", "kingdom")
  numCol <- length(nomeCol)
  # retirar# nao precisa da funcao transNome
  # retirar#aTx <- matrix(c(b, gnrOutput_alterado$matched_name, 
  # retirar### (nome igual ao presente na base)
  # retirar##(por vezes desenecessariamente extenso, com o nome de quem decobriu e data)
  # precisa da funcao transNome
  aTx <- matrix(c(as.vector(gnrOutput_alterado$user_supplied_name),# era factor 
                  transNome(gnrOutput_alterado$matched_name), 
                  rep("",tam*(numCol-2))), # preenche com as info que ja tem 
                nrow = tam, ncol = numCol)
  colnames(aTx) <- nomeCol
  for(i in 1:tam){
    vtClass <- strsplit(as.vector(gnrOutput_alterado$classification_path[i]),
                        "\\|")[[1]]
    vtRanks <- strsplit(as.vector(gnrOutput_alterado$classification_path_ranks[i]),
                        "\\|")[[1]]
    cat(paste0(i, " | "))
    for(j in 1:numCol){
      if(length(vtClass[vtRanks == colnames(aTx)[j]]) > 0){
        aTx[i,j] <- vtClass[vtRanks == colnames(aTx)[j]]
      }
    }
  }
  return(aTx)
}

# funcao para lidar com as linhas vazias na col. "species"
# recebe a matriz e retorna ela mesma, mas com as lin. dessa col. preenchidas
# sao os nomes corretos e atuais q vao montar as redes
setSpecies <- function(matrizClass){
  matrizClassNova <- matrizClass
  for(i in 1:length(matrizClass[, 1])){ # para percorrer a matriz
    if(matrizClass[i, 3] == ""){ # onde esta vazio
      # Se o nome eh unico
      testVF <- strsplit(matrizClass[i, 2], " ")
      if(!(length(testVF[[1]]) > 1)){
        #como as col. 1 e 2 sao preenchidas com os nomes mais proximos ao 
        # fornecido, a proxima col. com informacao eh a que vai ser pega
        # pelo indice 3
        # ex. se o prox. nome eh genus a col sp sera preenc. com genus
        # ex. se o prox. nome eh fam a col sp sera preenc. com fam
        matrizClassNova[i, 3] <- 
          as.vector(matrizClass[i, matrizClass[i ,] != ""])[3]
      }else{
        # Se eh binomial (nome de especie)
        # Se o 1o nome do binomio corresponde ao genero da col. genus
        if(testVF[[1]][1] == matrizClass[i, 5]){
          matrizClassNova[i, 3] <- matrizClass[i, 2]
        }else{
          # Se o 1o nome do binomio nao corresponde ao genero da col. genus
          # substitui o 1o nome do binomio pelo nome presente na col. genus 
          matrizClassNova[i, 3] <- paste(matrizClass[i, 5], testVF[[1]][2])
        }
      }
    }
  }
  return(matrizClassNova)
}

# retirar as linhas da matriz de classificacao q apresentam alguma 
#inconsistencia:
# sem genero e sem familia
# sem classificacoes (alem da coluna de especie q ja foi preenchida 
#na outra func, a "setSpecies")
tiraLnhRuins <- function(mtTx){
  # todas as linhas q nao estao com genero & sem familia
  novaMtTx <- mtTx[!(mtTx[, 5] !=  "" & mtTx[, 7] ==  ""), ]
  # todas as linhas q tem classificacao de filo
  novaMtTx <- mtTx[mtTx[, 11] !=  "", ]
  return(novaMtTx)
}

# junta a matriz com os nomes dos especimes nos arquivos 
#(nome original|nome simplificado) com a matriz com as classificacoes 
#taxonomicas
juntaMt <- function(mtArtigo, mtTx){
  numcol <- length(mtTx[1, ])
  numrow <- length(mtArtigo[, 1])
  nvMt <- matrix(rep("", numcol*numrow), nrow = numrow, ncol = numcol)
  colnames(nvMt) <- colnames(mtTx)
  rownames(nvMt) <- rownames(mtArtigo)
  nvMt[, 1:2] <- mtArtigo[, 1:2] # teoricamente td a matriz
  for(i in 1:numrow){
    linha <- mtTx[mtTx[, 1] == mtArtigo[i, 2], 3:numcol]
    if(length(linha) > 0) {
      nvMt[i, 3:numcol] <- mtTx[mtTx[, 1] == mtArtigo[i, 2], 3:numcol]
    }
  }
  nvMt <- nvMt[!(nvMt[, 3] == ""), ] # tira nomes sem taxonomia
  return(nvMt)
}


# colocar uma identificacao diferentes para especimes identificados ate genero
## ex: se ID "Alloscirtetica" e portanto o genero eh "Alloscirtetica" os 2 
## especimes do mesmo genero serao 
## "Alloscirtetica sp.1" e "Alloscirtetica sp.2"
difSpGen <- function(mtClassUni){
  # ordena pela col de interesse
  mtClassUni <- mtClassUni[order(mtClassUni[, 3]), ] 
  ## primeiro soh coloca "sp."
  # quando so esta o genero
  for (i in 1:length(mtClassUni[, 1])) {
    if (mtClassUni[i, 3] == mtClassUni[i, 5]) { # species = genus
      mtClassUni[i, 3] <- paste0(mtClassUni[i, 3]," sp.")
    }
  }
  # quando so esta o subgenero
  for (i in 1:length(mtClassUni[, 1])) {
    if (mtClassUni[i, 3] == mtClassUni[i, 4]) {
      # nao vou trab com resolucao ate sub genero (linha logo abaixo)
      # mtClassUni[i, 3] <- paste0(mtClassUni[i, 5], "(", mtClassUni[i, 4],") sp.")
      mtClassUni[i, 3] <- paste0(mtClassUni[i, 5]," sp.")
    }
  }
  # depois verifica se ha mais de uma ocorrencia do genero sp. (para colocar num)
  i <- 1
  while (i < length(mtClassUni[, 1])) {
    # verifica se tem sp. no nome
    if (length(grep("sp.", mtClassUni[i, 3])) > 0) {
      j <- 1
      # verifica se ha iguais e armazena (!!precisa estar em ordem alfabetica!!)
      # para nao ficar para sempre aqui: "i < length(mtClassUni[, 1]) no 
      ## while de cima"
      while (mtClassUni[i, 3] == mtClassUni[i+1, 3] & i < length(mtClassUni[, 3])) { 
        i <- i + 1
        j <- j + 1
      }
      # se tem nomes iguais, passa a ter uma numeracao diferente
      if (j > 1) {
        # como esta em i, e tem qts sao iguais pelo valor de j, entao
        numRep <- j
        while (j > 0){
          mtClassUni[i-(numRep-j), 3] <- paste0(mtClassUni[i-(numRep-j), 3],j)
          j <- j - 1
        }
      }
    }
    i <- i + 1
  }
  return(mtClassUni)
}

# tira linhas com nomes repetidos
# somando os valores presentes nas colunas com linhas com o mesmo nome
linhSemRep <- function(mtLnhNm){
  idUnico <- unique(rownames(mtLnhNm))
  tam <- length(idUnico)
  tamCol <- length(mtLnhNm[1, ])
  mttemp <- matrix(data= rep(0, tam*tamCol), ncol = tamCol, nrow = tam)
  # ve se tem mais de uma coluna
  if (tamCol > 1){
    # cria uma matriz de 0 sem linhas repetidas (matriz final "vazia")
    colnames(mttemp) <- colnames(mtLnhNm)
    rownames(mttemp) <- idUnico
    i <- 1
    for (i in 1:tam){
      test <- mtLnhNm[idUnico[i] == rownames(mtLnhNm), ]
      # se test eh um vecto (so tem 1 linha), se nao eh matriz
      if (is.vector(test)){
        mttemp[i,] <- test # mantem a linha 
        i <- i + 1
      } else {
        mttemp[i,] <- apply(test, 2, sum) # tira repeticao
        i <- i + length(test[, 1])
      }
      #cat(paste0(i, " | "))
    }
  }else{
    for (i in 1:tam){
      mttemp[i,] <- sum(mtLnhNm[idUnico[i] == rownames(mtLnhNm),1])
    }
  }
  return(mttemp)
}

# tira inconsitencias nas matrizes de interacoes
# tira linhas e colunas sem interacao
# soma linhas e colunas com nomes repetidos (usa func linhSemRep)
# torna a matriz binaria (usa func linhSemRep)
ajeitaListMt <- function(listaMt){
  tam <- length(listaMt)
  submtx <- vector("list", length = tam)
  cat(paste0("mt without inconsistencies\n "))
  for(i in 1:tam){
    #cat(paste0("|",i, "|"))
    submtx[[i]] <- listaMt[[i]]
    # sem linhas rep 
    submtx[[i]] <- linhSemRep(submtx[[i]])
    # sem colunas repetidas
    submtx[[i]] <- t(linhSemRep(t(submtx[[i]])))
    # faz virar matriz binaria (1: interage, 0: nao interage)
    # F ou T para >0 e vezes 1 para virar numero
    submtx[[i]] <- (submtx[[i]] > 0) * 1 
  }
  return(submtx)
}

# substitui os nomes das linhas e colunas pelos nomes presentes na matriz de
# classificacao taxonomica da coluna desejada (numero): 
#       col 3: "species" (nomes com a grafia correta)
#       col >3: qualquer taxonomias presentes na matriz
# retira linhas/colunas sem nomes corretos correspondentes
# precisa da lista de matrizes, matriz de 2 colunas com nome original e 
#simplificado para os nomes nas colunas e nas linhas, bem como as matrizes
#de classificacao taxonomica respectivas
# combinar diferentes agrupamentos de animais e de plantas
# ex: planta agrupada por genero e animais por sp
# 
# para matrizes que deixaram de ser matriz (só uma entrada com nome 
#correspondente nas linhas ou colunas: 
# - mantem a matriz original
# - ultimo item da lista eh um vetor com o indice das matrizes problematicas
# - retirar essas matrizes, pois nessa etapa eh necessario que se tenha ao 
#   menos 2 colunas ou linhas para ter a possibilidade do efeito do agrupamento
subClass <- function(listMtMdOri,
                     mtClassAni, # matriz c as classificacoes dos animais
                     mtClassPlt, # matriz c as classificacoes dos plantas
                     nColA = 3, # qual col da mt sera utilizada p agr. animais
                     nColP = 3){# qual col da mt sera utilizada p agr. plantas
  # inicia variaveis
  tamList <- length(listMtMdOri)
  ltModif <- vector("list", length = tamList)
  nmNaoEncontradosCol <- vector("character")
  nmNaoEncontradosRow <- vector("character")
  mtSemLinhas <- vector("numeric") # nenhum dos nomes foi encontrado mas bases de dados
  for(j in 1:tamList){
    submtx <- listMtMdOri[[j]]
    # colunas (com animais)
    temp <- colnames(submtx)
    nmCol <- vector("character") # sera preenchido com os nomes desejados
    colDesj <- vector("numeric") # onde foi encontrada correspondencia
    tam <- length(temp)
    for(i in 1:tam){
      nmDesej <- mtClassAni[temp[i] == mtClassAni[, 1], nColA]
      if(length(nmDesej)>0){
        nmCol <- c(nmCol, nmDesej)
        colDesj <- c(colDesj, i)
      }else{
        nmNaoEncontradosCol <- c(nmNaoEncontradosCol, temp[i])
      }
      # cat(paste0(i, " | "))
    }
	submtx <- submtx[ ,colDesj]
	if (length(nmCol) > 0 & is.matrix(submtx)){
	   colnames(submtx) <- nmCol
	} else {
	   mtSemLinhas <- c(mtSemLinhas, j)
	}
    
    # linhas (com plantas)
	if(is.matrix(submtx)){ # so faz para linha se ainda eh uma matriz
		temp <- rownames(submtx)
		nmRow <- vector("character") # sera preenchido com os nomes desejados
		rowDesj <- vector("numeric")
		tam <- length(temp)
		for(i in 1:tam){
		  nmDesej <- mtClassPlt[temp[i] == mtClassPlt[, 1], nColP]
		  if(length(nmDesej)>0){
			nmRow <- c(nmRow, nmDesej)
			rowDesj <- c(rowDesj, i)
		  }else{
			nmNaoEncontradosRow <- c(nmNaoEncontradosRow, temp[i])
		  }
		  # cat(paste0(i, " | "))
		}
		submtx <- submtx[rowDesj, ]
		if (length(nmRow) > 0 & is.matrix(submtx)){
		   rownames(submtx) <- nmRow
		} else {
		   mtSemLinhas <- c(mtSemLinhas, j)
		}
	}else{
		submtx <- listMtMdOri[[j]] # mantem a original
	}
	
	# retira linhas ou colunas sem interacoes
	if (is.matrix(submtx)){
		# seleciona as linhas (plantas) que tem alguma interação
		submtx <- submtx[apply(submtx,1,sum)>0,]
	    if(is.matrix(submtx)){
			# grep com 1 para selecionar so as colunas (animais) com interacao 
			#(TRUE p >0)
			submtx <- submtx[,grep(1,(apply(submtx,2,sum)>0)*1)]
		}
	}
	
	# se a matriz sem linhas ou colunas de soma 0 virou um vetor = matriz probl
	# se a matriz virou uma matriz de dimensao 0 = matriz probl
    if (!is.matrix(submtx) | length(submtx) < 1){
		mtSemLinhas <- c(mtSemLinhas, j)
	}
	
    # passa para a lista a matriz alterada
    ltModif[[j]] <- submtx
    cat(paste0(j, " | "))
  }
  
  # acrescenta matrizes que precisam ser retiradas (nao foram alteradas)
  temp <- vector(mode='list', length= 1)
  temp[[1]] <- unique(mtSemLinhas)
  ltModif <- c(ltModif, temp)
  return(ltModif)
}
# fim das funcoes #############################################################
###############################################################################



###############################################################################
# leitura das redes e armazenamento em lista e sp em vetores ##################
# criacao de 2 tabelas
#       tabela (inicialmente um vetor) com todas as sp de animais
#       tabela (inicialmente um vetor) com todas as sp de plantas

# nomes dos arquivos
caminhoArqs <- paste0("../adjMatrices/",webType,"-all/")
arqs <- dir(caminhoArqs)

# lista que vai armazenar os arquivos
listArqs <- vector("list", length = length(arqs))

# vetor que vao armazenar nomes de todas sp
AnimNm <- vector(mode= "character")
PlanNm <- vector(mode= "character")

# preenche as listas 
for(i in 1:length(arqs)){
  # importa cada artigo
  listArqs[[i]] <- read.csv(paste0(caminhoArqs,arqs[i]), header = F)
  # armazena as sp
  # para plantas nas linhas
  PlanNm <- append(PlanNm, as.vector(listArqs[[i]][-1, 1]))
  # para animas nas colunas
  AnimNm <- append(AnimNm, as.vector(t(listArqs[[i]][1, -1])))
}
# mudar o identificador de cada uma das redes presente na lista
names(listArqs) <- arqs


# tira repeticoes
PlanNm <- unique(PlanNm)
AnimNm <- unique(AnimNm)


###############################################################################
# tira os nomes que nao tem significado algum para identificacao ##############
PlanNm <- tiraSemSent(PlanNm)
AnimNm <- tiraSemSent(AnimNm)


###############################################################################
# tira os nomes que nao sao de especie (nao sao nomes compostos) ##############
PlanNm <- selecSp(PlanNm)
AnimNm <- selecSp(AnimNm)


###############################################################################
# normalizacao dos nomes das especies #########################################
# mantem os nomes da forma que estavam no arquivo original 
#(vector -> matrix com 2 col)
#(para facilitar montar a tabela com os nomes corrigidos)
sPlanNm <- simplNorm(PlanNm)
sAnimNm <- simplNorm(AnimNm)


###############################################################################
# procurar a grafia correta para cada nome ####################################
# uso da func gnr_resolve do pacote taxize (demora para terminar de rodar)
# <http = "post"> esse parametro para queries grandes nao funcionou...
gnrSpPln <- nomeCorreto(unique(sPlanNm[, 2])) 
gnrSpAni <- nomeCorreto(unique(sAnimNm[, 2]))



###############################################################################
# constroi uma matriz com as classificacoes taxonomicas disponiveis ###########
tplanTx <- classifTx(gnrSpPln)
tanimTx <- classifTx(gnrSpAni)


###############################################################################
# lidar com nomes que nao sao mais aceitos (ex. genero nao existe mais)########
# na coluna "species" algumas linhas estao em branco
# preenche com a informacao mais correta possivel
# (pois e a partir de onde serao construidas as redes)
tplanTx <- setSpecies(tplanTx)
tanimTx <- setSpecies(tanimTx)


###############################################################################
# tirar inconsistencias na matriz de classificacao ############################
# com genero, mas sem familia
# sem classificacoes (alem da coluna de especie q ja foi preenchida 
#na outra func, a "setSpecies")
tplanTx <- tiraLnhRuins(tplanTx)
tanimTx <- tiraLnhRuins(tanimTx)

###############################################################################
# selecao de grupos de interesse ##############################################
# nomes de plantas q produzem flores (Magnoliopsida)
# nomes de animais q correspondem a insetos (Insecta)
#planTx <- tplanTx[tplanTx[,10] == "Magnoliopsida",] 
planTx <- tplanTx[tplanTx[,10] == "Magnoliopsida" |
                    tplanTx[,10] == "Liliopsida",]
animTx <- tanimTx[tanimTx[,10] == "Insecta",]

###############################################################################
# registra a quantidade de nomes nos arquivos #################################
write(c(paste0("Animals:\n There are ", 
               length(sAnimNm[,1])," names in all matrices, \n",
               length(tanimTx[,1])," names have a correct classification and ",
               length(animTx[,1])," names were selected."),
        paste0("Plants:\n There are ", 
               length(sPlanNm[,1])," names in all matrices, \n",
               length(tplanTx[,1])," names have a correct classification and ",
               length(planTx[,1])," names were selected.")), 
      file = paste0("txtoutput/namesNumberOf", length(listArqs),
                    "matrices.txt"), sep = "\n")

###############################################################################
# arrumar os nomes da matriz de cada classificacao para:
#   - sem as classificacoes mais abrangentes (> genero)

# na lista de classificacoes taxonomicas
## 1 - todo especime com nome diferente
# 1o passo
# matriz com a juncao dos nomes como aparecem no arquivos e a taxonomia
nvplanTx <- juntaMt(sPlanNm, planTx)
nvanimTx <- juntaMt(sAnimNm, animTx)

# 2o passo
# colocar uma identificacao diferentes para especimes identificados ate genero
# coloca sp. e um numero para diferenciar
nvplanTx <- difSpGen(nvplanTx)
nvanimTx <- difSpGen(nvanimTx)
# 3o passo
# sem as classificacoes mais abrangentes (> genero)
# linhas em q a coluna genus (col. 5) nao esta em branco
gNvplanTx <- nvplanTx[!(nvplanTx[, 5] == ""), ]
gNvanimTx <- nvanimTx[!(nvanimTx[, 5] == ""), ]

# 4o passo
# tem algumas linhas que estao sem a classif por familia para animais
a <- gNvanimTx[gNvanimTx[, 7] == "", ]
# busca em cada uma das linhas q a familia esta faltando
for(i in 1:length(a[, 1])){
  temp <- classification(a[i,5], db = "gbif", 
                         rows = 1) # considera apenas a primeira opcao 
                          #(qd tem mais de um match)
  if (length(temp[[1]]$name[temp[[1]]$rank == "family"]) > 0){
    a[i, 7] <- temp[[1]]$name[temp[[1]]$rank == "family"]
  }else{
    a[i, 7] <- paste0(a[i, 9], " (awaiting allocation)")
  }
}
gNvanimTx[gNvanimTx[, 7] == "", ] <- a 
# algumas familias ainda permaceram em branco pois sao como o exemplo abaixo:
# # no site: https://www.gbif-uat.org/species/107204798
# Family  Hymenoptera (awaiting allocation)

write.csv(gNvplanTx, file = paste0('../output/',webType,'_gNvplanTx_',
                                   format(Sys.time(), "%Y-%m-d%d"),'.csv'))
write.csv(gNvanimTx, file = paste0('../output/',webType,'_gNvanimTx_',
                                   format(Sys.time(), "%Y-%m-d%d"),'.csv'))

###############################################################################
# criar as diferentes listas de matrizes com base 
## na lista de arquivos e na matriz class

# cria lista com as matrizes com os nomes originais
#(presentes nos arquivos lidos como data.frames, nomes dados pelos autores das 
#matrizes)
tam <- length(listArqs)
listMtOrig <- vector("list", length = tam)
for(i in 1:tam){
  listMtOrig[[i]] <- as.matrix(listArqs[[i]][-1,-1])
  listMtOrig[[i]] <- matrix(as.numeric(listMtOrig[[i]]),
                            ncol = length(listMtOrig[[i]][1,]),
                            nrow = length(listMtOrig[[i]][,1]))
  colnames(listMtOrig[[i]]) <- t(listArqs[[i]][1, -1])
  rownames(listMtOrig[[i]]) <- listArqs[[i]][-1, 1]
  cat(i, " | ")
}

###############################################################################
# torna a matriz binaria
listMtMdOrig <- vector("list", length = tam)
for(i in 1:tam){
  listMtMdOrig[[i]] <- (listMtOrig[[i]] > 0)*1
}

###############################################################################
# data.frame com os tipo de agrupamento taxonomicos ###########################
nmTxInteres <- c("species","genus","family","order")
nuTxInteres <- match(nmTxInteres, colnames(gNvplanTx))
tam <- length(nmTxInteres)
# codigo p identificar o tipo de agrupamento de maneira abreviada
ID <- c("sp","ge","fa","or")
ID <- as.matrix(data.frame(A = "A_",
                              AGpBy =  rep(ID, times= tam),
                              X = ".x.",
                              P = "P_",
                              PGpBy =  rep(ID, each = tam)))
for(i in 1:length(ID[,1])){
  ID[i,1] <- paste(ID[i,], collapse= "")
}
ID <- ID[,1]
gps <- data.frame(ID = ID,
                  AnimalnumGp =  rep(nuTxInteres, times= tam),
                  PlantnumGp =  rep(nuTxInteres, each= tam),
                  AnimalGroupedBy =  rep(nmTxInteres, times= tam),
                  PlantGroupedBy =  rep(nmTxInteres, each = tam))

###############################################################################
# lista de lista de matrizes e calculo das metricas p cd matriz ###############
# cada conjunto de matrizes dos arquivos em um item da lista
# lista com as listas de matrizes + calculo das metricas
# cada item da lista de matrizes (lstGps) tem um lista com as matrizes
lstGps <- vector("list", length = tam^2)
# armazena vetor com os indices da matrizes com problemas
# matrizes com problemas: matrizes q não foi encontrado nomes de plantas ou 
# animais nas bases de dados. Resultam em matrizes com 0 linhas ou colunas.
mtProblem <- vector('numeric')
numMt <- length(listMtMdOrig) # qtdd de matrizes
for(i in 1:tam^2){
  # cria uma lista de matrizes com base no agrupamento do tipo "i"
  lstGps[[i]] <- subClass(listMtMdOri = listMtMdOrig,
                       mtClassAni = gNvanimTx,
                       mtClassPlt = gNvplanTx,
                       # ql coluna vai pegar da matriz de taxonomia
                       nColA = gps$AnimalnumGp[i], 
                       nColP = gps$PlantnumGp[i])
 

  # pega o ultimo item da lista (vetor com os indices da matrizes com 
  #problemas)
  mtProblem <- c(mtProblem, lstGps[[i]][[numMt+1]])
  lstGps[[i]] <- lstGps[[i]][-(numMt+1)] # tira ultimo elemento da lista
}

# matriz com linhas e colunas validas (id de sp), 
#mas sem agrupar nomes iguais (sem passar pela func ajeitaListMt)
semMod <- listMtMdOrig
listMtMdOrig <- lstGps[[1]]

# se tem matriz com problema...
if(length(mtProblem) > 0){
  # todas as matrizes com problemas
  mtProblem <- unique(mtProblem)
  semModVal <- semMod[-mtProblem]
  for(i in 1:tam^2){
  	# tira matrizes com problemas
  	lstGps[[i]] <- lstGps[[i]][-mtProblem]
  	# tira as entradas repetidas nas linhas e colunas da lista de matrizes
  	lstGps[[i]] <- ajeitaListMt(lstGps[[i]])
  }
} else { # se nao...
  semModVal <- semMod
  for(i in 1:tam^2){
    # tira as entradas repetidas nas linhas e colunas da lista de matrizes
    lstGps[[i]] <- ajeitaListMt(lstGps[[i]])
  }
  mtProblem <- -c(1:numMt)
}
arqsComProbl <- arqs[mtProblem]
arqs <- arqs[-mtProblem]
write.csv(arqsComProbl, file = '../output/mtProblem.csv')

# sumario das redes com agrupamento
info <- data.frame(Tipo= rep(webType, numMt*tam^2),
                   GroupType=  vector('character', length = numMt*tam^2),
                   ArqName= vector('character', length = numMt*tam^2),
                   numAnimals= vector('numeric', length = numMt*tam^2),
                   numPlants= vector('numeric', length = numMt*tam^2),
                   numSpecies= vector('numeric', length = numMt*tam^2),
                   numInterac= vector('numeric', length = numMt*tam^2))
info$ArqName <- as.character(info$ArqName)
info$GroupType <- as.character(info$GroupType)
for(j in 1:tam^2){
  for(i in 1:numMt){
    info$ArqName[(j-1)*numMt+i] <- arqs[i]
    info$GroupType[(j-1)*numMt+i] <- as.character(gps$ID[j])
    info$numPlants[(j-1)*numMt+i] <- length(lstGps[[j]][[i]][,1])
    info$numAnimals[(j-1)*numMt+i] <- length(lstGps[[j]][[i]][1,])
    info$numSpecies[(j-1)*numMt+i] <- info$numPlants[(j-1)*numMt+i] + 
      info$numAnimals[(j-1)*numMt+i]
    info$numInterac[(j-1)*numMt+i] <- sum(lstGps[[j]][[i]])
  }
}


# sumario das redes sem problemas
listMtMdOrigSemProblem <- listMtMdOrig[-mtProblem]
numMt <- length(listMtMdOrigSemProblem) # qtdd de matrizes
info1 <- data.frame(Tipo= rep(webType, numMt),
                    GroupType=  "Original",
                    ArqName= vector('character', length = numMt),
                    numAnimals= vector('numeric', length = numMt),
                    numPlants= vector('numeric', length = numMt),
                    numSpecies= vector('numeric', length = numMt),
                    numInterac= vector('numeric', length = numMt))
info1$ArqName <- as.character(info1$ArqName)
for(i in 1:numMt){
  info1$ArqName[i] <- arqs[i]
  info1$numPlants[i] <- length(listMtMdOrigSemProblem[[i]][,1])
  info1$numAnimals[i] <- length(listMtMdOrigSemProblem[[i]][1,])
  info1$numSpecies[i] <- info1$numPlants[i] + info1$numAnimals[i]
  info1$numInterac[i] <- sum(listMtMdOrigSemProblem[[i]])
}
# junta info das redes originais com as agrupadas
info <- rbind(info1,info)
write.csv(info, file = paste0('../output/charact_mtInfo',
                              format(Sys.time(), "%Y-%m-d%d"),'.csv'))

###############################################################################
# matriz com as diferentes formas de agrupar ##################################
# matrizes para o calculo das metricas
# nao calcula para as matrizes com problemas (zero linha ou colunas)

#metricas gerais
#metric <- c("connectance","links per species")
metric <- c("connectance","links per species","nestedness", "NODF",
            "number of compartments","web asymmetry")

# calculo das metricas
for(j in 1:length(metric)){
  cat("\n",metric[j],": ")
  vtAux <- vector("numeric")
  # matrizes sem alteracoes
  for(i in 1:numMt){
    vtAux[i] <- networklevel(listMtMdOrigSemProblem[[i]],
                             index = metric[j]) # nome na metrica
  }
  # matrizes agrupadas
  for(w in 1:tam^2){
    for(i in 1:numMt){
      vtAux[(w-1)*numMt+i+numMt] <- networklevel(lstGps[[w]][[i]],
                                                 index = metric[j])
      cat("=")
    }
    cat("\n")
  }
  info <- cbind(info,vtAux)
  colnames(info)[length(info[1,])] <- metric[j]
}

#metricas por level
metricL <- c("robustness","niche overlap")
vtAux <- vector("list", length =  length(metricL))
for(j in 1:length(metricL)){
  cat("\n",metricL[j],": \n")
  # matrizes sem alteracoes
  for(i in 1:numMt){
    vtAux[[j]] <- rbind(vtAux[[j]],
                        networklevel(listMtMdOrigSemProblem[[i]], 
                                     index = metricL[j])) # nome da metrica
    cat("=")
  }
  cat("\n")
  # matrizes agrupadas
  for(w in 1:tam^2){
    cat("w")
    for(i in 1:numMt){
      vtAux[[j]] <- rbind(vtAux[[j]],
                          networklevel(lstGps[[w]][[i]], 
                                       index = metricL[j]))
      cat("=")
    }
    cat("\n")
  }
}
info <- cbind(info,as.data.frame(vtAux))

write.csv(info, file = paste0('../output/charact_mtMetrics_',webType,"_",
                              format(Sys.time(), "%Y-%m-d%d"), '.csv'))

# fim da parte sem supervisionamento ##########################################
###############################################################################
save.image(paste0("txtoutput/ws_",webType,"_", format(Sys.time(), "%Y-%m-d%d"),
                  "_v01.RData"))

#_____________________________________________________________________________#

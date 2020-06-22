##### USE THIS ONE!


# fixed yibar2 issue in gamma.sex and betas .sex
## using glimmer - it's a named integer from Fstat.R 
#### it is part of whats causing the problems (I think), so it's part of
#### troubleshooting... to find it, need to go back to Fstat.R... i think...
#### it represents portoricensis males - the 7 males, only 1 of which tested positive



setwd("/Users/annarosesjodin/Documents/Dissertation/Herpes")

raw.data <- read.csv("OSMat91_No_NoclepSteruf.csv", header=TRUE)

# get rid of extra rows
raw.data <- raw.data[-42,]
raw.data <- raw.data[-42,]
raw.data <- raw.data[,-1]
raw.data <- raw.data[,-1088]
raw.data <- raw.data[,-1087]


# rename the rows
rownames(raw.data) <- c("OTU 01", "OTU 10", "OTU 11", "OTU 12", "OTU 13", "OTU 14", "OTU 15",
                        "OTU 16", "OTU 17", "OTU 18", "OTU 19", "OTU 02", "OTU 20", "OTU 21",
                        "OTU 22", "OTU 23", "OTU 24", "OTU 25", "OTU 26", "OTU 27", "OTU 28",
                        "OTU 29", "OTU 03", "OTU 30", "OTU 31", "OTU 32", "OTU 33", "OTU 34",
                        "OTU 35", "OTU 36", "OTU 37", "OTU 38", "OTU 39", "OTU 04", "OTU 40",
                        "OTU 41", "OTU 05", "OTU 06", "OTU 07", "OTU 08", "OTU 09")

# make sure only bats from Culebrones are included
raw.data.Cul <- raw.data[, -grep("Artibeus", colnames(raw.data))]
raw.data.Cul <- raw.data.Cul[, -grep("Eptesicus", colnames(raw.data.Cul))]
raw.data.Cul <- raw.data.Cul[which(rowSums(raw.data.Cul) !=0), ]

# or Larva
raw.data.Lar <- raw.data[, -grep("Pteronotus", colnames(raw.data))]
raw.data.Lar <- raw.data.Lar[, -grep("Monophyllus", colnames(raw.data.Lar))]
raw.data.Lar <- raw.data.Lar[, -grep("Mormoops", colnames(raw.data.Lar))]
raw.data.Lar <- raw.data.Lar[, -grep("Erophylla", colnames(raw.data.Lar))]
raw.data.Lar <- raw.data.Lar[, -grep("Brachyphylla", colnames(raw.data.Lar))]
raw.data.Lar <- raw.data.Lar[which(rowSums(raw.data.Lar) !=0), ]

# add full word for sex so that we don't have the M being Monophyllus or Mormoops issue
females <- grep("F", colnames(raw.data.Cul))
colnames(raw.data.Cul)[-females] <- paste(colnames(raw.data.Cul)[-females], "ale", sep="")
colnames(raw.data.Cul)[females] <- paste(colnames(raw.data.Cul)[females], "emale", sep="")

females <- grep("F", colnames(raw.data))
colnames(raw.data)[-females] <- paste(colnames(raw.data)[-females], "ale", sep="")
colnames(raw.data)[females] <- paste(colnames(raw.data)[females], "emale", sep="")

females <- grep("F", colnames(raw.data.Lar))
colnames(raw.data.Lar)[-females] <- paste(colnames(raw.data.Lar)[-females], "ale", sep="")
colnames(raw.data.Lar)[females] <- paste(colnames(raw.data.Lar)[females], "emale", sep="")


# create a function that separates out groups (e.g. species, sex, treatment)

find_group7 <- function(mat, var1, var2=NULL){
  if(class(mat) == "matrix" | class(mat) == "data.frame"){
    group <- mat[ , grep(var1 , colnames(mat))]
    if(! is.null(var2)){
      group <- group[ , grep(var2, colnames(group))]
    }
  } else {
    if(class(mat) == "vector" | class(mat) == "integer"){
      group <- mat[grep(var1 , names(mat))]
      if(! is.null(var2)){
        group <- group[grep(var2, names(group))]
      }
    }
  }
  
  return(group)
}

freud7 <- find_group7(raw.data.Cul, var1 = "quadridens")
frodo <- find_group7(raw.data.Cul, var1 = "quadridens", var2 = "Male")
fifi <- find_group7(raw.data.Cul, var1 = "Male")

### create a function that calculates alpha

remove_garbage <- function(x){
  if(length(x)==0){
    x <- 0
  } else if(is.na(x)){
    x <- 0
  } else if(length(x)!=1){
    x <- 0
  } else{x <- x}
  
  return(x)}

# this alpha is calculated from groups, as output using find_group5
calc_alpha7 <- function(mat){
    alpha <- mean(colSums(mat))
alpha <- remove_garbage(alpha)
    return(alpha)
}

calc_alpha7(freud7) #pte qua

### create a function that calculates gamma

calc_gamma7 <- function(mat){
    gamm <- (nrow(mat[which(rowSums(mat) !=0),]))
    gamm <- remove_garbage(gamm)
  return(gamm)
}

calc_gamma7(freud7) # pte qua - should be 15

### create a function that calclualtes beta multiplicative

calc_beta_mult7 <- function(mat){
    betamult <- (nrow(mat[which(rowSums(mat) !=0),])) / (mean(colSums(mat))) 
    betamult <- remove_garbage(betamult)
  return(betamult)
}

calc_beta_mult7(freud7)

### create a function that calculates beta additive

calc_beta_add7 <- function(mat){
    betaadd <- (nrow(mat[which(rowSums(mat) !=0),])) - (mean(colSums(mat)))
    betaadd <- remove_garbage(betaadd)
  return(betaadd)
}

calc_beta_add7(freud7)

# calculate stats for each group - get means

calc_alpha7(mat = find_group5(mat = raw.data, var1 = "Male"))
calc_alpha7(mat = find_group5(mat = raw.data, var1 = "Female"))
calc_alpha7(mat = find_group5(mat = raw.data, var1 = "quadridens"))
calc_alpha7(mat = find_group5(mat = raw.data, var1 = "portoricensis"))
calc_alpha7(mat = find_group5(mat = raw.data, var1 = "Mormoops"))
calc_alpha7(mat = find_group5(mat = raw.data, var1 = "Monophyllus"))
calc_alpha7(mat = find_group5(mat = raw.data, var1 = "Erophylla"))
calc_alpha7(mat = find_group5(mat = raw.data, var1 = "Brachyphylla"))
calc_alpha7(mat = find_group5(mat = raw.data, var1 = "Artibeus"))
calc_alpha7(mat = find_group5(mat = raw.data, var1 = "quadridens", var2 = "Male"))
calc_alpha7(mat = find_group5(mat = raw.data, var1 = "portoricensis", var2 = "Male"))
calc_alpha7(mat = find_group5(mat = raw.data, var1 = "Mormoops", var2 = "Male"))
calc_alpha7(mat = find_group5(mat = raw.data, var1 = "Monophyllus", var2 = "Male"))
calc_alpha7(mat = find_group5(mat = raw.data, var1 = "Erophylla", var2 = "Male"))
calc_alpha7(mat = find_group5(mat = raw.data, var1 = "Brachyphylla", var2 = "Male"))
calc_alpha7(mat = find_group5(mat = raw.data, var1 = "Artibeus", var2 = "Male"))
calc_alpha7(mat = find_group5(mat = raw.data, var1 = "quadridens", var2 = "Female"))
calc_alpha7(mat = find_group5(mat = raw.data, var1 = "portoricensis", var2 = "Female"))
calc_alpha7(mat = find_group5(mat = raw.data, var1 = "Mormoops", var2 = "Female"))
calc_alpha7(mat = find_group5(mat = raw.data, var1 = "Monophyllus", var2 = "Female"))
calc_alpha7(mat = find_group5(mat = raw.data, var1 = "Erophylla", var2 = "Female"))
calc_alpha7(mat = find_group5(mat = raw.data, var1 = "Brachyphylla", var2 = "Female"))
calc_alpha7(mat = find_group5(mat = raw.data, var1 = "Artibeus", var2 = "Female"))

calc_gamma7(mat = find_group5(mat = raw.data, var1 = "Male"))
calc_gamma7(mat = find_group5(mat = raw.data, var1 = "Female"))
calc_gamma7(mat = find_group5(mat = raw.data, var1 = "quadridens"))
calc_gamma7(mat = find_group5(mat = raw.data, var1 = "portoricensis"))
calc_gamma7(mat = find_group5(mat = raw.data, var1 = "Mormoops"))
calc_gamma7(mat = find_group5(mat = raw.data, var1 = "Monophyllus"))
calc_gamma7(mat = find_group5(mat = raw.data, var1 = "Erophylla"))
calc_gamma7(mat = find_group5(mat = raw.data, var1 = "Brachyphylla"))
calc_gamma7(mat = find_group5(mat = raw.data, var1 = "Artibeus"))
calc_gamma7(mat = find_group5(mat = raw.data, var1 = "quadridens", var2 = "Male"))
calc_gamma7(mat = find_group5(mat = raw.data, var1 = "portoricensis", var2 = "Male"))
calc_gamma7(mat = find_group5(mat = raw.data, var1 = "Mormoops", var2 = "Male"))
calc_gamma7(mat = find_group5(mat = raw.data, var1 = "Monophyllus", var2 = "Male"))
calc_gamma7(mat = find_group5(mat = raw.data, var1 = "Erophylla", var2 = "Male"))
calc_gamma7(mat = find_group5(mat = raw.data, var1 = "Brachyphylla", var2 = "Male"))
calc_gamma7(mat = find_group5(mat = raw.data, var1 = "Artibeus", var2 = "Male"))
calc_gamma7(mat = find_group5(mat = raw.data, var1 = "quadridens", var2 = "Female"))
calc_gamma7(mat = find_group5(mat = raw.data, var1 = "portoricensis", var2 = "Female"))
calc_gamma7(mat = find_group5(mat = raw.data, var1 = "Mormoops", var2 = "Female"))
calc_gamma7(mat = find_group5(mat = raw.data, var1 = "Monophyllus", var2 = "Female"))
calc_gamma7(mat = find_group5(mat = raw.data, var1 = "Erophylla", var2 = "Female"))
calc_gamma7(mat = find_group5(mat = raw.data, var1 = "Brachyphylla", var2 = "Female"))
calc_gamma7(mat = find_group5(mat = raw.data, var1 = "Artibeus", var2 = "Female"))

calc_beta_mult7(mat = find_group5(mat = raw.data, var1 = "Male"))
calc_beta_mult7(mat = find_group5(mat = raw.data, var1 = "Female"))
calc_beta_mult7(mat = find_group5(mat = raw.data, var1 = "quadridens"))
calc_beta_mult7(mat = find_group5(mat = raw.data, var1 = "portoricensis"))
calc_beta_mult7(mat = find_group5(mat = raw.data, var1 = "Mormoops"))
calc_beta_mult7(mat = find_group5(mat = raw.data, var1 = "Monophyllus"))
calc_beta_mult7(mat = find_group5(mat = raw.data, var1 = "Erophylla"))
calc_beta_mult7(mat = find_group5(mat = raw.data, var1 = "Brachyphylla"))
calc_beta_mult7(mat = find_group5(mat = raw.data, var1 = "Artibeus"))
calc_beta_mult7(mat = find_group5(mat = raw.data, var1 = "quadridens", var2 = "Male"))
calc_beta_mult7(mat = find_group5(mat = raw.data, var1 = "portoricensis", var2 = "Male"))
calc_beta_mult7(mat = find_group5(mat = raw.data, var1 = "Mormoops", var2 = "Male"))
calc_beta_mult7(mat = find_group5(mat = raw.data, var1 = "Monophyllus", var2 = "Male"))
calc_beta_mult7(mat = find_group5(mat = raw.data, var1 = "Erophylla", var2 = "Male"))
calc_beta_mult7(mat = find_group5(mat = raw.data, var1 = "Brachyphylla", var2 = "Male"))
calc_beta_mult7(mat = find_group5(mat = raw.data, var1 = "Artibeus", var2 = "Male"))
calc_beta_mult7(mat = find_group5(mat = raw.data, var1 = "quadridens", var2 = "Female"))
calc_beta_mult7(mat = find_group5(mat = raw.data, var1 = "portoricensis", var2 = "Female"))
calc_beta_mult7(mat = find_group5(mat = raw.data, var1 = "Mormoops", var2 = "Female"))
calc_beta_mult7(mat = find_group5(mat = raw.data, var1 = "Monophyllus", var2 = "Female"))
calc_beta_mult7(mat = find_group5(mat = raw.data, var1 = "Erophylla", var2 = "Female"))
calc_beta_mult7(mat = find_group5(mat = raw.data, var1 = "Brachyphylla", var2 = "Female"))
calc_beta_mult7(mat = find_group5(mat = raw.data, var1 = "Artibeus", var2 = "Female"))

calc_beta_add7(mat = find_group5(mat = raw.data, var1 = "Male"))
calc_beta_add7(mat = find_group5(mat = raw.data, var1 = "Female"))
calc_beta_add7(mat = find_group5(mat = raw.data, var1 = "quadridens"))
calc_beta_add7(mat = find_group5(mat = raw.data, var1 = "portoricensis"))
calc_beta_add7(mat = find_group5(mat = raw.data, var1 = "Mormoops"))
calc_beta_add7(mat = find_group5(mat = raw.data, var1 = "Monophyllus"))
calc_beta_add7(mat = find_group5(mat = raw.data, var1 = "Erophylla"))
calc_beta_add7(mat = find_group5(mat = raw.data, var1 = "Brachyphylla"))
calc_beta_add7(mat = find_group5(mat = raw.data, var1 = "Artibeus"))
calc_beta_add7(mat = find_group5(mat = raw.data, var1 = "quadridens", var2 = "Male"))
calc_beta_add7(mat = find_group5(mat = raw.data, var1 = "portoricensis", var2 = "Male"))
calc_beta_add7(mat = find_group5(mat = raw.data, var1 = "Mormoops", var2 = "Male"))
calc_beta_add7(mat = find_group5(mat = raw.data, var1 = "Monophyllus", var2 = "Male"))
calc_beta_add7(mat = find_group5(mat = raw.data, var1 = "Erophylla", var2 = "Male"))
calc_beta_add7(mat = find_group5(mat = raw.data, var1 = "Brachyphylla", var2 = "Male"))
calc_beta_add7(mat = find_group5(mat = raw.data, var1 = "Artibeus", var2 = "Male"))
calc_beta_add7(mat = find_group5(mat = raw.data, var1 = "quadridens", var2 = "Female"))
calc_beta_add7(mat = find_group5(mat = raw.data, var1 = "portoricensis", var2 = "Female"))
calc_beta_add7(mat = find_group5(mat = raw.data, var1 = "Mormoops", var2 = "Female"))
calc_beta_add7(mat = find_group5(mat = raw.data, var1 = "Monophyllus", var2 = "Female"))
calc_beta_add7(mat = find_group5(mat = raw.data, var1 = "Erophylla", var2 = "Female"))
calc_beta_add7(mat = find_group5(mat = raw.data, var1 = "Brachyphylla", var2 = "Female"))
calc_beta_add7(mat = find_group5(mat = raw.data, var1 = "Artibeus", var2 = "Female"))


### create a function that creates a dataframe of: 
#### groups, alpha, beta.add, beta.mult, gamma

##### requires an input list defining treatments

species7 <- list(
  c("quadridens"), c("portoricensis"), c("Mormoops"),
  c("Monophyllus"), c("Erophylla"), c("Brachyphylla"),
  c("Artibeus")
)

treatments8 <- list(
  c("quadridens", "Male"), c("quadridens", "Female"),
  c("portoricensis", "Male"), c("portoricensis", "Female"),
  c("Mormoops", "Male"), c("Mormoops", "Female"),
  c("Monophyllus", "Male"), c("Monophyllus", "Female"),
  c("Erophylla", "Male"), c("Erophylla", "Female"),
  c("Brachyphylla", "Male"), c("Brachyphylla", "Female"),
  c("Artibeus", "Male"), c("Artibeus", "Female")
)

sex <- list(
  c("Male"), c("Female")
)


# create a dataframe of empirical values of a,b,g for each group, plus n
abg_summ7 <- function(groups, mat){
  grpnames <- c()
  all_alpha <-c()
  all_beta_add <- c()
  all_beta_mult <- c()
  all_gamma <- c()
  ni <- c()
  for(i in 1:length(groups)){
    current_group <- groups[[i]]
    if(length(current_group) > 1){ # look for those that have species and sex
      x <- find_group7(mat=mat, var1=current_group[1], var2=current_group[2])
    } else {
      x <- find_group7(mat=mat, var1=current_group[1])
    }
    grpnames[i] <- paste(current_group, collapse="_")
    all_alpha[i] <- calc_alpha7(x)
    all_beta_add[i] <- calc_beta_add7(x)
    all_beta_mult[i] <- calc_beta_mult7(x)
    all_gamma[i] <- calc_gamma7(x)
    ni[i] <- ncol(x)
    if(i==length(groups)){
      mydat <- cbind(grpnames, all_alpha, all_beta_add, all_beta_mult, all_gamma, ni)
    }
  }
  return(mydat)
}


### create a function that calculates F statistic

F.alpha8.sex <- function(mat, groups){ # groups must == list(c("Male"), c("Female"))
  # define constants
  K <- 2
  N <- ncol(mat)
  # i <- 2
  yij <- colSums(mat)
  sex1 <- find_group7(mat = mat, var1 = groups[1])
  sex2 <- find_group7(mat = mat, var1 = groups[2])
  yij.males <- colSums(sex1)
  yij.females <- colSums(sex2)
  top <- c()
  yibar <- c() # will have length = 2 
  abg_df <- abg_summ7(mat = mat, groups = groups)
  ybar <- mean(yij)
  males <- c()
  females <- c()
  # calculate numerator
  for(i in 1:K){
    ni <- as.numeric(abg_df[i, 6])
    yibar[i] <- as.numeric(abg_df[i, 2])
    top[i] <- ((ni*((yibar[i] - ybar)^2))/(K-1))
  }
  # calculate denom
  for(j in 1:length(yij.males)){
    males[j] <- (((yij.males[j]-yibar[1])^2)/(N-K))
  }
  for(k in 1:length(yij.females)){
    females[k] <- (((yij.females[k]-yibar[2])^2)/(N-K))
  }
  num <- sum(top)
  denom <- sum(males) + sum(females)
  F_stat <- num/denom
  return(F_stat)
}

F.alpha8.species <- function(mat, groups){ # groups must == species7
  # define constants
  K <- 7
  N <- ncol(mat)
  # i <- 7
  yij <- colSums(mat)
  spp1 <- find_group7(mat = mat, var1 = groups[1])
  spp2 <- find_group7(mat = mat, var1 = groups[2])
  spp3 <- find_group7(mat = mat, var1 = groups[3])
  spp4 <- find_group7(mat = mat, var1 = groups[4])
  spp5 <- find_group7(mat = mat, var1 = groups[5])
  spp6 <- find_group7(mat = mat, var1 = groups[6])
  spp7 <- find_group7(mat = mat, var1 = groups[7])
  yij.spp1 <- colSums(spp1)
  yij.spp2 <- colSums(spp2)
  yij.spp3 <- colSums(spp3)
  yij.spp4 <- colSums(spp4)
  yij.spp5 <- colSums(spp5)
  yij.spp6 <- colSums(spp6)
  yij.spp7 <- colSums(spp7)
  top <- c()
  yibar <- c() # will have length = 8
  abg_df <- abg_summ7(mat = mat, groups = groups)
  ybar <- mean(yij)
  spp.1 <- c()
  spp.2 <- c()
  spp.3 <- c()
  spp.4 <- c()
  spp.5 <- c()
  spp.6 <- c()
  spp.7 <- c()
  # calculate numerator
  for(i in 1:K){
    ni <- as.numeric(abg_df[i, 6])
    yibar <- as.numeric(abg_df[, 2])
    top[i] <- ((ni*((yibar[i] - ybar)^2))/(K-1))
  }
  # calculate denom
  for(j in 1:length(yij.spp1)){
    spp.1[j] <- (((yij.spp1[j]-yibar[1])^2)/(N-K))
  }
  for(k in 1:length(yij.spp2)){
    spp.2[k] <- (((yij.spp2[k]-yibar[2])^2)/(N-K))
  }
  for(l in 1:length(yij.spp3)){
    spp.3[l] <- (((yij.spp3[l]-yibar[3])^2)/(N-K))
  }
  for(m in 1:length(yij.spp4)){
    spp.4[m] <- (((yij.spp4[m]-yibar[4])^2)/(N-K))
  }
  for(n in 1:length(yij.spp5)){
    spp.5[n] <- (((yij.spp5[n]-yibar[5])^2)/(N-K))
  }
  for(o in 1:length(yij.spp6)){
    spp.6[o] <- (((yij.spp6[o]-yibar[6])^2)/(N-K))
  }
  for(p in 1:length(yij.spp7)){
    spp.7[p] <- (((yij.spp7[p]-yibar[7])^2)/(N-K))
  }
  num <- sum(top)
  denom <- sum(spp.1) + sum(spp.2) + sum(spp.3) + sum(spp.4) + sum(spp.5) + sum(spp.6) + sum(spp.7)
  F_stat <- num/denom
  return(F_stat)
}

F.alpha8.trt <- function(mat, groups){ # groups must == treatments8
  # define constants
  K <- 14
  N <- ncol(mat)
  # i <- 14
  yij <- colSums(mat)
  trt1 <- find_group7(mat = mat, var1 = groups[[1]][1], var2 = groups[[1]][2])
  trt2 <- find_group7(mat = mat, var1 = groups[[2]][1], var2 = groups[[2]][2])
  trt3 <- find_group7(mat = mat, var1 = groups[[3]][1], var2 = groups[[3]][2])
  trt4 <- find_group7(mat = mat, var1 = groups[[4]][1], var2 = groups[[4]][2])
  trt5 <- find_group7(mat = mat, var1 = groups[[5]][1], var2 = groups[[5]][2])
  trt6 <- find_group7(mat = mat, var1 = groups[[6]][1], var2 = groups[[6]][2])
  trt7 <- find_group7(mat = mat, var1 = groups[[7]][1], var2 = groups[[7]][2])
  trt8 <- find_group7(mat = mat, var1 = groups[[8]][1], var2 = groups[[8]][2])
  trt9 <- find_group7(mat = mat, var1 = groups[[9]][1], var2 = groups[[9]][2])
  trt10 <- find_group7(mat = mat, var1 = groups[[10]][1], var2 = groups[[10]][2])
  trt11 <- find_group7(mat = mat, var1 = groups[[11]][1], var2 = groups[[11]][2])
  trt12 <- find_group7(mat = mat, var1 = groups[[12]][1], var2 = groups[[12]][2])
  trt13 <- find_group7(mat = mat, var1 = groups[[13]][1], var2 = groups[[13]][2])
  trt14 <- find_group7(mat = mat, var1 = groups[[14]][1], var2 = groups[[14]][2])
  yij.trt1 <- colSums(trt1)
  yij.trt2 <- colSums(trt2)
  yij.trt3 <- colSums(trt3)
  yij.trt4 <- colSums(trt4)
  yij.trt5 <- colSums(trt5)
  yij.trt6 <- colSums(trt6)
  yij.trt7 <- colSums(trt7)
  yij.trt8 <- colSums(trt8)
  yij.trt9 <- colSums(trt9)
  yij.trt10 <- colSums(trt10)
  yij.trt11 <- colSums(trt11)
  yij.trt12 <- colSums(trt12)
  yij.trt13 <- colSums(trt13)
  yij.trt14 <- colSums(trt14)
  top <- c()
  yibar <- c() # will have length = 14
  abg_df <- abg_summ7(mat = mat, groups = groups)
  ybar <- mean(yij)
  trt.1 <- c()
  trt.2 <- c()
  trt.3 <- c()
  trt.4 <- c()
  trt.5 <- c()
  trt.6 <- c()
  trt.7 <- c()
  trt.8 <- c()
  trt.9 <- c()
  trt.10 <- c()
  trt.11 <- c()
  trt.12 <- c()
  trt.13 <- c()
  trt.14 <- c()
  # calculate numerator
  for(i in 1:K){
    ni <- as.numeric(abg_df[i, 6])
    yibar <- as.numeric(abg_df[, 2])
    top[i] <- ((ni*((yibar[i] - ybar)^2))/(K-1))
  }
  # calculate denom
  for(j in 1:length(yij.trt1)){
    trt.1[j] <- (((yij.trt1[j]-yibar[1])^2)/(N-K))
  }
  for(k in 1:length(yij.trt2)){
    trt.2[k] <- (((yij.trt2[k]-yibar[2])^2)/(N-K))
  }
  for(l in 1:length(yij.trt3)){
    trt.3[l] <- (((yij.trt3[l]-yibar[3])^2)/(N-K))
  }
  for(m in 1:length(yij.trt4)){
    trt.4[m] <- (((yij.trt4[m]-yibar[4])^2)/(N-K))
  }
  for(n in 1:length(yij.trt5)){
    trt.5[n] <- (((yij.trt5[n]-yibar[5])^2)/(N-K))
  }
  for(o in 1:length(yij.trt6)){
    trt.6[o] <- (((yij.trt6[o]-yibar[6])^2)/(N-K))
  }
  for(p in 1:length(yij.trt7)){
    trt.7[p] <- (((yij.trt7[p]-yibar[7])^2)/(N-K))
  }
  for(q in 1:length(yij.trt8)){
    trt.8[q] <- (((yij.trt8[q]-yibar[8])^2)/(N-K))
  }
  for(r in 1:length(yij.trt9)){
    trt.9[r] <- (((yij.trt9[r]-yibar[9])^2)/(N-K))
  }
  for(s in 1:length(yij.trt10)){
    trt.10[s] <- (((yij.trt10[s]-yibar[10])^2)/(N-K))
  }
  for(t in 1:length(yij.trt11)){
    trt.11[t] <- (((yij.trt11[t]-yibar[11])^2)/(N-K))
  }
  for(u in 1:length(yij.trt12)){
    trt.12[u] <- (((yij.trt12[u]-yibar[12])^2)/(N-K))
  }
  for(v in 1:length(yij.trt13)){
    trt.13[v] <- (((yij.trt13[v]-yibar[13])^2)/(N-K))
  }
  for(w in 1:length(yij.trt14)){
    trt.14[w] <- (((yij.trt14[w]-yibar[14])^2)/(N-K))
  }
  num <- sum(top)
  denom <- sum(trt.1) + sum(trt.2) + sum(trt.3) + sum(trt.4) + 
    sum(trt.5) + sum(trt.6) + sum(trt.7) + sum(trt.8) +
    sum(trt.9) + sum(trt.10) + sum(trt.11) + sum(trt.12) +
    sum(trt.13) + sum(trt.14)
  F_stat <- num/denom
  return(F_stat)
}

F.alpha8.trt(mat = raw.data, groups = treatments8)

## create function to calculate f.gamma.sex

F.gamma.sex8 <- function(mat, groups){ # groups must == treatments8
  # define constants
  K <- 2
  N <- 14
  # i = 2 
  abg_df <- abg_summ7(mat = mat, groups = groups)
  yij <- as.numeric(abg_df[,5]) # treatment gammas - should be 14 values
  # calculate numerator
  top <- c()
  bottom <- c()
  yibar <- c() # will be two values
  ybar <- mean(yij) # mean gamma for all 14 treatments, only a single value
  jdat <- matrix(nrow = (N/K), ncol = K)
  jdat[,1] <- as.numeric(abg_df[c(1,3,5,7,9,11,13),5])
  jdat[,2] <- as.numeric(abg_df[c(2,4,6,8,10,12,14),5])
  jdat2 <- matrix(nrow = (N/K), ncol = K)
  for(i in 1:K){
    ni <- 7
    yibar[i] <- as.numeric(mean(jdat[,i]))
    top[i] <- as.numeric(((ni*((yibar[i] - ybar)^2))/(K-1)))
    # calculate denom
    for(j in 1:nrow(jdat2)){
      jdat2[j, i] <- as.numeric((((jdat[j,i]-yibar[i])^2)/(N-K)))
    }
  }
  bottom <- as.numeric(colSums(jdat)) # sum across all j (7 values)
  num <- sum(top) # sum across all K (2 values)
  denom <- sum(bottom) # sum across all K (2 values)
  F_stat <- num/denom
  return(F_stat)
}

F.gamma.sex8(mat = raw.data, groups = treatments8)

# f stat for beta mult sex

F.betam.sex8 <- function(mat, groups){ # groups must == treatments8
  # define constants
  K <- 2
  N <- 14
  # i = 2
  abg_df <- abg_summ7(mat = mat, groups = groups)
  yij <- as.numeric(abg_df[,4]) # treatment betams - should be 14 values
  # calculate numerator
  top <- c()
  bottom <- c()
  yibar <- c() # will be two values
  ybar <- mean(yij) # mean gamma for all 14 treatments, only a single value
  jdat <- matrix(nrow = (N/K), ncol = K)
  jdat[,1] <- as.numeric(abg_df[c(1,3,5,7,9,11,13),4])
  jdat[,2] <- as.numeric(abg_df[c(2,4,6,8,10,12,14),4])
  jdat2 <- matrix(nrow = (N/K), ncol = K)
  for(i in 1:K){
    ni <- 7
    yibar[i] <- as.numeric(mean(jdat[,i]))
    top[i] <- as.numeric(((ni*((yibar[i] - ybar)^2))/(K-1)))
    # calculate denom
    for(j in 1:nrow(jdat2)){
      jdat2[j, i] <- as.numeric((((jdat[j,i]-yibar[i])^2)/(N-K)))
    }
  }
  bottom <- as.numeric(colSums(jdat)) # sum across all j (7 values)
  num <- sum(top) # sum across all K (2 values)
  denom <- sum(bottom) # sum across all K (2 values)
  F_stat <- num/denom
  return(F_stat)
}

F.betam.sex8(mat = raw.data, groups = treatments8)

# f stat for beta add sex

F.betaa.sex8 <- function(mat, groups){ # groups must == treatments8
  # define constants
  K <- 2
  N <- 14
  # i = 2 
  abg_df <- abg_summ7(mat = mat, groups = groups)
  yij <- as.numeric(abg_df[,3]) # treatment betams - should be 14 values
  # calculate numerator
  top <- c()
  bottom <- c()
  yibar <- c() # will be two values
  ybar <- mean(yij) # mean gamma for all 14 treatments, only a single value
  jdat <- matrix(nrow = (N/K), ncol = K)
  jdat[,1] <- as.numeric(abg_df[c(1,3,5,7,9,11,13),3])
  jdat[,2] <- as.numeric(abg_df[c(2,4,6,8,10,12,14),3])
  jdat2 <- matrix(nrow = (N/K), ncol = K)
  for(i in 1:K){
    ni <- 7
    yibar[i] <- as.numeric(mean(jdat[,i]))
    top[i] <- as.numeric(((ni*((yibar[i] - ybar)^2))/(K-1)))
    # calculate denom
    for(j in 1:nrow(jdat2)){
      jdat2[j, i] <- as.numeric((((jdat[j,i]-yibar[i])^2)/(N-K)))
    }
  }
  bottom <- as.numeric(colSums(jdat)) # sum across all j (8 values)
  num <- sum(top) # sum across all K (2 values)
  denom <- sum(bottom) # sum across all K (2 values)
  F_stat <- num/denom
  return(F_stat)
}

F.betaa.sex8(mat = raw.data, groups = treatments8)

## create function to calculate f.gamma.species

F.gamma.spp8 <- function(mat, groups){ # groups must == treatments8
  # define constants
  K <- 7
  N <- 14
  # i = 7
  abg_df <- abg_summ7(mat = mat, groups = groups)
  yij <- as.numeric(abg_df[,5]) # treatment gammas - should be 14 values
  # calculate numerator
  top <- c()
  bottom <- c()
  yibar <- c() # will be two values
  ybar <- mean(yij) # mean gamma for all 14 treatments, only a single value
  jdat <- matrix(nrow = (N/K), ncol = K)
  jdat[,1] <- as.numeric(abg_df[c(1,2),5])
  jdat[,2] <- as.numeric(abg_df[c(3,4),5])
  jdat[,3] <- as.numeric(abg_df[c(5,6),5])
  jdat[,4] <- as.numeric(abg_df[c(7,8),5])
  jdat[,5] <- as.numeric(abg_df[c(9,10),5])
  jdat[,6] <- as.numeric(abg_df[c(11,12),5])
  jdat[,7] <- as.numeric(abg_df[c(13,14),5])
  jdat2 <- matrix(nrow = (N/K), ncol = K)
  for(i in 1:K){
    ni <- 2
    yibar[i] <- as.numeric(mean(jdat[,i]))
    top[i] <- as.numeric(((ni*((yibar[i] - ybar)^2))/(K-1)))
    # calculate denom
    for(j in 1:nrow(jdat2)){
      jdat2[j, i] <- as.numeric((((jdat[j,i]-yibar[i])^2)/(N-K)))
    }
  }
  bottom <- as.numeric(colSums(jdat)) # sum across all j (2 values)
  num <- sum(top) # sum across all K (7 values)
  denom <- sum(bottom) # sum across all K (7 values)
  F_stat <- num/denom
  return(F_stat)
}

F.gamma.spp8(mat = raw.data, groups = treatments8)

## create function to calculate f.betam.species

F.betam.spp8 <- function(mat, groups){ # groups must == treatments8
  # define constants
  K <- 7
  N <- 14
  # i = 7 
  abg_df <- abg_summ7(mat = mat, groups = groups)
  yij <- as.numeric(abg_df[,4]) # treatment gammas - should be 14 values
  # calculate numerator
  top <- c()
  bottom <- c()
  yibar <- c() # will be two values
  ybar <- mean(yij) # mean gamma for all 14 treatments, only a single value
  jdat <- matrix(nrow = (N/K), ncol = K)
  jdat[,1] <- as.numeric(abg_df[c(1,2),4])
  jdat[,2] <- as.numeric(abg_df[c(3,4),4])
  jdat[,3] <- as.numeric(abg_df[c(5,6),4])
  jdat[,4] <- as.numeric(abg_df[c(7,8),4])
  jdat[,5] <- as.numeric(abg_df[c(9,10),4])
  jdat[,6] <- as.numeric(abg_df[c(11,12),4])
  jdat[,7] <- as.numeric(abg_df[c(13,14),4])
  jdat2 <- matrix(nrow = (N/K), ncol = K)
  for(i in 1:K){
    ni <- 2
    yibar[i] <- as.numeric(mean(jdat[,i]))
    top[i] <- as.numeric(((ni*((yibar[i] - ybar)^2))/(K-1)))
    # calculate denom
    for(j in 1:nrow(jdat2)){
      jdat2[j, i] <- as.numeric((((jdat[j,i]-yibar[i])^2)/(N-K)))
    }
  }
  bottom <- as.numeric(colSums(jdat)) # sum across all j (2 values)
  num <- sum(top) # sum across all K (7 values)
  denom <- sum(bottom) # sum across all K (7 values)
  F_stat <- num/denom
  return(F_stat)
}


F.betam.spp8(mat = raw.data, groups = treatments8)

## create function to calculate f.betaa.species

F.betaa.spp8 <- function(mat, groups){ # groups must == treatments8
  # define constants
  K <- 7
  N <- 14
  # i = 7 
  abg_df <- abg_summ7(mat = mat, groups = groups)
  yij <- as.numeric(abg_df[,3]) # treatment gammas - should be 14 values
  # calculate numerator
  top <- c()
  bottom <- c()
  yibar <- c() # will be two values
  ybar <- mean(yij) # mean gamma for all 14 treatments, only a single value
  jdat <- matrix(nrow = (N/K), ncol = K)
  jdat[,1] <- as.numeric(abg_df[c(1,2),3])
  jdat[,2] <- as.numeric(abg_df[c(3,4),3])
  jdat[,3] <- as.numeric(abg_df[c(5,6),3])
  jdat[,4] <- as.numeric(abg_df[c(7,8),3])
  jdat[,5] <- as.numeric(abg_df[c(9,10),3])
  jdat[,6] <- as.numeric(abg_df[c(11,12),3])
  jdat[,7] <- as.numeric(abg_df[c(13,14),3])
  jdat2 <- matrix(nrow = (N/K), ncol = K)
  for(i in 1:K){
    ni <- 2
    yibar[i] <- as.numeric(mean(jdat[,i]))
    top[i] <- as.numeric(((ni*((yibar[i] - ybar)^2))/(K-1)))
    # calculate denom
    for(j in 1:nrow(jdat2)){
      jdat2[j, i] <- as.numeric((((jdat[j,i]-yibar[i])^2)/(N-K)))
    }
  }
  bottom <- as.numeric(colSums(jdat)) # sum across all j (2 values)
  num <- sum(top) # sum across all K (8 values)
  denom <- sum(bottom) # sum across all K (8 values)
  F_stat <- num/denom
  return(F_stat)
}

F.betaa.spp8(mat = raw.data, groups = treatments8)



#############################################################################
#---------------------------------------------------------------------------#
#############################################################################

## create null distribution of f statistics
### sampling WITHOUT replacement

## simulate the empirical matrix 1000 times

# make an array to store 1000 simulations
sim.array8 <- array(data = NA, dim = c(nrow(raw.data),ncol(raw.data), 1000))
colnames(sim.array8) <- colnames(raw.data)

# fill array with 1000 simulations of re-shuffled columns
for(i in 1:dim(sim.array8)[3]){
  sim.array8[,,i] <- as.matrix(raw.data[,sample(ncol(raw.data), dim(raw.data)[2], replace = F)])
}


## make a df that will be filled in with F statistics for simulations
### needs to be 1000 rows and enough columns for all f.stat measures

sim.fstats8 <- data.frame(matrix(NA, nrow = 1000, ncol = 9))
colnames(sim.fstats8) <- c("a.sex", "a.spp", "a.treat", "bm.sex", "bm.spp",
                          "ba.sex", "ba.spp", "g.sex", "g.spp")

for(i in 1:dim(sim.array8)[3]){
  sim.fstats8[i,1] <- F.alpha8.sex(mat = sim.array8[,,i], groups = sex) #a.sex
  sim.fstats8[i,2] <- F.alpha8.species(mat = sim.array8[,,i], groups = species7) #a.spp
  sim.fstats8[i,3] <- F.alpha8.trt(mat = sim.array8[,,i], groups = treatments8) #a.treat
  sim.fstats8[i,4] <- F.betam.sex8(mat = sim.array8[,,i], groups = treatments8) #bm.sex
  sim.fstats8[i,5] <- F.betam.spp8(mat = sim.array8[,,i], groups = treatments8) #bm.spp
  sim.fstats8[i,6] <- F.betaa.sex8(mat = sim.array8[,,i], groups = treatments8) #ba.sex
  sim.fstats8[i,7] <- F.betaa.spp8(mat = sim.array8[,,i], groups = treatments8) #ba.spp
  sim.fstats8[i,8] <- F.gamma.sex8(mat = sim.array8[,,i], groups = treatments8) #g.sex
  sim.fstats8[i,9] <- F.gamma.spp8(mat = sim.array8[,,i], groups = treatments8) #g.spp
}

### calculate empirical stuff

emp.fstats8 <- data.frame(matrix(NA, nrow = 1, ncol = 9))
colnames(emp.fstats8) <- c("a.sex", "a.spp", "a.treat", "bm.sex", "bm.spp",
                           "ba.sex", "ba.spp", "g.sex", "g.spp")
emp.fstats8[1,1] <- F.alpha8.sex(mat = raw.data, groups = sex)
emp.fstats8[1,2] <- F.alpha8.species(mat = raw.data, groups = species7)
emp.fstats8[1,3] <- F.alpha8.trt(mat = raw.data, groups = treatments8)
emp.fstats8[1,4] <- F.betam.sex8(mat = raw.data, groups = treatments8)
emp.fstats8[1,5] <- F.betam.spp8(mat = raw.data, groups = treatments8)
emp.fstats8[1,6] <- F.betaa.sex8(mat = raw.data, groups = treatments8)
emp.fstats8[1,7] <- F.betaa.spp8(mat = raw.data, groups = treatments8)
emp.fstats8[1,8] <- F.gamma.sex8(mat = raw.data, groups = treatments8)
emp.fstats8[1,9] <- F.gamma.spp8(mat = raw.data, groups = treatments8)


hist(sim.fstats8$a.sex)
abline(v = emp.fstats8[1,1], col = "red", lwd = 3)
hist(sim.fstats8$a.spp, breaks = 20, ylim = c(0,150))
abline(v = emp.fstats8[1,2], col = "red", lwd = 3)
hist(sim.fstats8$a.treat, breaks = 20, xlim = c(0,4), ylim = c(0,200))
abline(v = emp.fstats8[1,3], col = "red", lwd = 3)
hist(sim.fstats8$bm.sex, breaks = 20, xlim = c(4,4.8), ylim = c(0,250))
abline(v = emp.fstats8[1,4], col = "red", lwd = 3)
hist(sim.fstats8$bm.spp, breaks = 20)
abline(v = emp.fstats8[1,5], col = "red", lwd = 3)
hist(sim.fstats8$ba.sex, breaks = 20)
abline(v = emp.fstats8[1,6], col = "red", lwd = 3)
hist(sim.fstats8$ba.spp, breaks = 20)
abline(v = emp.fstats8[1,7], col = "red", lwd = 3)
hist(sim.fstats8$g.sex, breaks = 20)
abline(v = emp.fstats8[1,8], col = "red", lwd = 3)
hist(sim.fstats8$g.spp, breaks = 20)
abline(v = emp.fstats8[1,9], col = "red", lwd = 3)

# how many f stats are greater than empirical
# proportion greater = p-value

(length(which(sim.fstats8$a.sex > emp.fstats8[1,1])))/1000
(length(which(sim.fstats8$a.spp > emp.fstats8[1,2])))/1000
(length(which(sim.fstats8$a.treat > emp.fstats8[1,3])))/1000
(length(which(sim.fstats8$bm.sex > emp.fstats8[1,4])))/1000
(length(which(sim.fstats8$bm.spp > emp.fstats8[1,5])))/1000
(length(which(sim.fstats8$ba.sex > emp.fstats8[1,6])))/1000
(length(which(sim.fstats8$ba.spp > emp.fstats8[1,7])))/1000
(length(which(sim.fstats8$g.sex > emp.fstats8[1,8])))/1000
(length(which(sim.fstats8$g.spp > emp.fstats8[1,9])))/1000


#############################################################################
# do all the simulations with replace = T

# make an array to store 1000 simulations
sim.array6.replaced <- array(data = NA, dim = c(nrow(raw.data),ncol(raw.data), 1000))
colnames(sim.array6.replaced) <- colnames(raw.data)

# fill array with 1000 simulations of re-shuffled columns
for(i in 1:dim(sim.array6.replaced)[3]){
  sim.array6.replaced[,,i] <- as.matrix(raw.data[,sample(ncol(raw.data), dim(raw.data)[2], replace = T)])
}


## make a df that will be filled in with F statistics for simulations
### needs to be 1000 rows and enough columns for all f.stat measures

sim.fstats6.replaced <- data.frame(matrix(NA, nrow = 1000, ncol = 9))
colnames(sim.fstats6.replaced) <- c("a.sex", "a.spp", "a.treat", "bm.sex", "bm.spp",
                                   "ba.sex", "ba.spp", "g.sex", "g.spp")

for(i in 1:dim(sim.array6.replaced)[3]){
  sim.fstats6.replaced[i,1] <- F.alpha8.sex(mat = sim.array6.replaced[,,i], groups = sex) #a.sex
  sim.fstats6.replaced[i,2] <- F.alpha8.species(mat = sim.array6.replaced[,,i], groups = species7) #a.spp
  sim.fstats6.replaced[i,3] <- F.alpha8.trt(mat = sim.array6.replaced[,,i], groups = treatments8) #a.treat
  sim.fstats6.replaced[i,4] <- F.betam.sex8(mat = sim.array6.replaced[,,i], groups = treatments8) #bm.sex
  sim.fstats6.replaced[i,5] <- F.betam.spp8(mat = sim.array6.replaced[,,i], groups = treatments8) #bm.spp
  sim.fstats6.replaced[i,6] <- F.betaa.sex8(mat = sim.array6.replaced[,,i], groups = treatments8) #ba.sex
  sim.fstats6.replaced[i,7] <- F.betaa.spp8(mat = sim.array6.replaced[,,i], groups = treatments8) #ba.spp
  sim.fstats6.replaced[i,8] <- F.gamma.sex8(mat = sim.array6.replaced[,,i], groups = treatments8) #g.sex
  sim.fstats6.replaced[i,9] <- F.gamma.spp8(mat = sim.array6.replaced[,,i], groups = treatments8) #g.spp
}

(length(which(sim.fstats6.replaced$a.sex > emp.fstats8[1,1])))/1000
(length(which(sim.fstats6.replaced$a.spp > emp.fstats8[1,2])))/1000
(length(which(sim.fstats6.replaced$a.treat > emp.fstats8[1,3])))/1000
(length(which(sim.fstats6.replaced$bm.sex > emp.fstats8[1,4])))/1000
(length(which(sim.fstats6.replaced$bm.spp > emp.fstats8[1,5])))/1000
(length(which(sim.fstats6.replaced$ba.sex > emp.fstats8[1,6])))/1000
(length(which(sim.fstats6.replaced$ba.spp > emp.fstats8[1,7])))/1000
(length(which(sim.fstats6.replaced$g.sex > emp.fstats8[1,8])))/1000
(length(which(sim.fstats6.replaced$g.spp > emp.fstats8[1,9])))/1000

hist(sim.fstats6.replaced$a.sex, breaks = 40)
abline(v = emp.fstats8[1,1], col = "red", lwd = 3)
hist(sim.fstats6.replaced$a.spp, breaks = 20)
abline(v = emp.fstats8[1,2], col = "red", lwd = 3)
hist(sim.fstats6.replaced$a.treat, breaks = 20)
abline(v = emp.fstats8[1,3], col = "red", lwd = 3)
hist(sim.fstats6.replaced$bm.sex, breaks = 20)
abline(v = emp.fstats8[1,4], col = "red", lwd = 3)
hist(sim.fstats6.replaced$bm.spp, breaks = 20)
abline(v = emp.fstats8[1,5], col = "red", lwd = 3)
hist(sim.fstats6.replaced$ba.sex, breaks = 20)
abline(v = emp.fstats8[1,6], col = "red", lwd = 3)
hist(sim.fstats6.replaced$ba.spp, breaks = 20)
abline(v = emp.fstats8[1,7], col = "red", lwd = 3)
hist(sim.fstats6.replaced$g.sex, breaks = 20)
abline(v = emp.fstats8[1,8], col = "red", lwd = 3)
hist(sim.fstats6.replaced$g.spp, breaks = 20)
abline(v = emp.fstats8[1,9], col = "red", lwd = 3)

########## for two-step simulation

trt.name6 <- c("quadridens_Male", "quadridens_Female", 
              "portoricensis_Male", "portoricensis_Female",
              "blainvillii_Male", "blainvillii_Female",
              "redmani_Male", "redmani_Female",
              "sezekorni_Male", "sezekorni_Female",
              "cavernarum_Male", "cavernarum_Female",
              "jamaicensis_Male", "jamaicensis_Female"
              )

step1 <- sample(trt.name6, 1)
step2 <- which(stringr::str_detect(colnames(raw.data), step1))
sample(step2, 1) # this will give an index value

which(stringr::str_detect(colnames(raw.data), "jamaicensis_Female"))

# make an array to store 1000 simulations
sim.array6.2step <- array(data = NA, dim = c(nrow(raw.data),ncol(raw.data), 1000))

for(i in 1:dim(sim.array6.2step)[3]){
  for(j in 1:dim(sim.array6.2step)[2]){
    step1 <- sample(trt.name6, 1)
    step2 <- which(stringr::str_detect(colnames(raw.data), step1))
    sim.array6.2step[,j,i] <- raw.data[,sample(step2,1)]
  }
}


# run this at the end so that the colnames arent replaced with the names
## of the randomly drawn columns
colnames(sim.array6.2step) <- colnames(raw.data)

# now do the f-stats
sim.fstats6.2step <- data.frame(matrix(NA, nrow = 1000, ncol = 9))
colnames(sim.fstats6.2step) <- c("a.sex", "a.spp", "a.treat", "bm.sex", "bm.spp",
                                "ba.sex", "ba.spp", "g.sex", "g.spp")

for(i in 1:dim(sim.array6.2step)[3]){
  sim.fstats6.2step[i,1] <- F.alpha8.sex(mat = sim.array6.2step[,,i], groups = sex) #a.sex
  sim.fstats6.2step[i,2] <- F.alpha8.species(mat = sim.array6.2step[,,i], groups = species7) #a.spp
  sim.fstats6.2step[i,3] <- F.alpha8.trt(mat = sim.array6.2step[,,i], groups = treatments8) #a.treat
  sim.fstats6.2step[i,4] <- F.betam.sex8(mat = sim.array6.2step[,,i], groups = treatments8) #bm.sex
  sim.fstats6.2step[i,5] <- F.betam.spp8(mat = sim.array6.2step[,,i], groups = treatments8) #bm.spp
  sim.fstats6.2step[i,6] <- F.betaa.sex8(mat = sim.array6.2step[,,i], groups = treatments8) #ba.sex
  sim.fstats6.2step[i,7] <- F.betaa.spp8(mat = sim.array6.2step[,,i], groups = treatments8) #ba.spp
  sim.fstats6.2step[i,8] <- F.gamma.sex8(mat = sim.array6.2step[,,i], groups = treatments8) #g.sex
  sim.fstats6.2step[i,9] <- F.gamma.spp8(mat = sim.array6.2step[,,i], groups = treatments8) #g.spp
}


(length(which(sim.fstats6.2step$a.sex > emp.fstats8[1,1])))/1000
(length(which(sim.fstats6.2step$a.spp > emp.fstats8[1,2])))/1000
(length(which(sim.fstats6.2step$a.treat > emp.fstats8[1,3])))/1000
(length(which(sim.fstats6.2step$bm.sex > emp.fstats8[1,4])))/1000
(length(which(sim.fstats6.2step$bm.spp > emp.fstats8[1,5])))/1000
(length(which(sim.fstats6.2step$ba.sex > emp.fstats8[1,6])))/1000
(length(which(sim.fstats6.2step$ba.spp > emp.fstats8[1,7])))/1000
(length(which(sim.fstats6.2step$g.sex > emp.fstats8[1,8])))/1000
(length(which(sim.fstats6.2step$g.spp > emp.fstats8[1,9])))/1000

hist(sim.fstats6.2step$a.sex)
abline(v = emp.fstats8[1,1], col = "red", lwd = 3)
hist(sim.fstats6.2step$a.spp, breaks = 20)
abline(v = emp.fstats8[1,2], col = "red", lwd = 3)
hist(sim.fstats6.2step$a.treat, breaks = 20)
abline(v = emp.fstats8[1,3], col = "red", lwd = 3)
hist(sim.fstats6.2step$bm.sex, breaks = 20)
abline(v = emp.fstats8[1,4], col = "red", lwd = 3)
hist(sim.fstats6.2step$bm.spp, breaks = 20)
abline(v = emp.fstats8[1,5], col = "red", lwd = 3)
hist(sim.fstats6.2step$ba.sex, breaks = 20)
abline(v = emp.fstats8[1,6], col = "red", lwd = 3)
hist(sim.fstats6.2step$ba.spp, breaks = 20)
abline(v = emp.fstats8[1,7], col = "red", lwd = 3)
hist(sim.fstats6.2step$g.sex, breaks = 20)
abline(v = emp.fstats8[1,8], col = "red", lwd = 3)
hist(sim.fstats6.2step$g.spp, breaks = 20)
abline(v = emp.fstats8[1,9], col = "red", lwd = 3)


calc_alpha7(find_group7(raw.data, var1 = "Male"))
calc_alpha7(find_group7(raw.data, var1 = "Female"))
calc_gamma7(find_group7(raw.data, var1 = "Male"))
calc_gamma7(find_group7(raw.data, var1 = "Female"))
calc_beta_add7((find_group7(raw.data, var1 = "Male")))
calc_beta_add7((find_group7(raw.data, var1 = "Female")))
calc_beta_mult7(find_group7(raw.data, var1 = "Male"))
calc_beta_mult7(find_group7(raw.data, var1 = "Female"))

calc_alpha7(mat = raw.data.Lar)
calc_alpha7(mat = raw.data.Cul)
calc_gamma7(raw.data.Lar)
calc_gamma7(raw.data.Cul)
calc_beta_add7(raw.data.Lar)
calc_beta_add7(raw.data.Cul)
calc_beta_mult7(raw.data.Lar)
calc_beta_mult7(raw.data.Cul)

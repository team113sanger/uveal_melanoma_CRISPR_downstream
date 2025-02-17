ccr.custom_impactOnPhenotype <- function(
    MO_uncorrectedFile,
    MO_correctedFile,
    sigFDR = 0.05,
    expName = 'expName',
    display = TRUE
) {
  
  beta_item <- paste0(expName, '|beta')
  fdr_item <- paste0(expName, '|fdr')

  pre <- read_tsv(MO_uncorrectedFile)
  pre <- pre[order(rownames(pre)), ]
  
  post <- read_tsv(MO_correctedFile)
  post <- post[order(rownames(post)), ]
  
  preD <- which(pre[[beta_item]] < 0 & pre[[fdr_item]] < sigFDR)
  preE <- which(pre[[beta_item]] >= 0 & pre[[fdr_item]] < sigFDR)
  preNULL <- setdiff(seq_len(nrow(pre)), c(preD, preE))

  postD <- which(post[[beta_item]] < 0 & post[[fdr_item]] < sigFDR)
  postE <- which(post[[beta_item]] >= 0 & post[[fdr_item]] < sigFDR)
  postNULL <- setdiff(seq_len(nrow(post)), c(postD, postE))
  
  aDD <- length(intersect(preD, postD))
  aDN <- length(intersect(preD, postNULL))
  aDE <- length(intersect(preD, postE))
  
  aND <- length(intersect(preNULL, postD))
  aNN <- length(intersect(preNULL, postNULL))
  aNE <- length(intersect(preNULL, postE))
  
  aED <- length(intersect(preE, postD))
  aEN <- length(intersect(preE, postNULL))
  aEE <- length(intersect(preE, postE))
  
  cm <- matrix(
    c(aDD, aDN, aDE, aND, aNN, aNE, aED, aEN, aEE),
    3,
    3,
    dimnames = list(c("cD", "cN", "cE"), c("uD", "uN", "uE"))
  )
  cm[is.na(cm)] <- 0
  
  IMPACTEDg <- 100 * sum(triu(cm, 1) + tril(cm, -1)) / sum(c(cm))
  IMPACTED_phenGenes <- 100 * (
    cm[2, 1] + cm[2, 3] + cm[3, 1] + cm[3, 2]) / sum(c(cm[, c(1, 3)])
    )
  
  IMPACTED_Depletions <- 100 * (cm[2, 1] + cm[2, 3]) / sum(cm[, 1])
  IMPACTED_Enrichments <- 100 * (cm[2, 3] + cm[1, 3]) / sum(cm[, 3])
  
  DISTORTEDg <- 100 * (cm[1, 3] + cm[3, 1]) / sum(c(cm))
  DISTORTED_phenGenes <- 100 * (cm[1, 3] + cm[3, 1]) / sum(c(cm[, c(1, 3)]))
  
  DISTORTED_Depletions <- 100 * cm[3, 1] / sum(cm[, 1])
  DISTORTED_Enrichments <- 100 * cm[1, 3] / sum(cm[, 3])
  
  geneCounts <- cm
  
  colnames(cm) <- paste(
    colSums(cm),
    c("loss of fitness", "no phenotype", "gain of fitness"),
    sep = "\n"
  )
  cm <- cm / t(matrix(rep(colSums(cm), nrow(cm)), 3, 3))
  
  if (display) {
    withr::with_par(
      list(mar = c(5, 4, 4, 10), xpd = TRUE), {
        barplot(
          100 * cm,
          col = c("red", "gray", "blue"),
          border = FALSE,
          main = expName,
          ylab = "%",
          xlab = "original counts"
        )
        legend(
          "right",
          c("loss of fitness", "no phenotype", "gain of fitness"),
          inset = c(-.5, 0),
          title = "Corrected counts",
          fill = c("red", "gray", "blue"),
          border = NA
        )
      },
      no.readonly = TRUE
    )
    
    withr::with_par(
      list(mfrow = c(2, 2), mar = c(0, 0, 2, 0), xpd = TRUE), {
        pie(
          c(IMPACTEDg, 100 - IMPACTEDg),
          col = c("blue", "white"),
          border = "gray",
          labels = c(paste0(format(IMPACTEDg, digits = 4), "%"), ""),
          main = "Overall impact"
        )
        pie(
          c(DISTORTEDg, 100 - DISTORTEDg),
          col = c("blue", "white"),
          border = "gray",
          labels = c(paste0(format(DISTORTEDg, digits = 4), "%"), ""),
          main = "Overall distortion"
        )
        pie(
          c(IMPACTED_phenGenes, 100 - IMPACTED_phenGenes),
          col = c("darkgreen", "white"),
          border = "gray",
          labels = c(paste0(
            format(IMPACTED_phenGenes, digits = 4), "%"
          ), ""),
          main = "Impact (G/L fitness genes)"
        )
        pie(
          c(DISTORTED_phenGenes, 100 - DISTORTED_phenGenes),
          col = c("darkgreen", "white"),
          border = "gray",
          labels = c(paste0(
            format(DISTORTED_phenGenes, digits = 4), "%"
          ), ""),
          main = "Distortion (G/L fitness genes)"
        )
      },
      no.readonly = TRUE
    )
  }
  
  dimnames(geneCounts) <- list(
    `corrected counts` = c("dep.", "null", "enr."),
    `original counts` = c("dep.", "null", "enr.")
  )
  
  id <- intersect(preD, postE)
  to_bind <- cbind(
    pre[id, c(beta_item, fdr_item)],
    post[id, c(beta_item, fdr_item)]
  )
  
  id <- intersect(preE, postD)
  to_bind <- rbind(
    to_bind,
    cbind(
      pre[id, c(beta_item, fdr_item)],
      post[id, c(beta_item, fdr_item)]
    )
  )
  colnames(to_bind) <- paste0(
    c("", "", "ccr.", "ccr."),
    colnames(to_bind)
  )
  
  id <- intersect(preD, postD)
  to_bind_c <- cbind(
    pre[id, c(beta_item, fdr_item)],
    post[id, c(beta_item, fdr_item)]
  )
  
  id <- intersect(preE, postE)
  to_bind_c <- rbind(
    to_bind_c,
    cbind(
      pre[id, c(beta_item, fdr_item)],
      post[id, c(beta_item, fdr_item)]
    )
  )
  colnames(to_bind_c) <- paste0(
    c("", "", "ccr.", "ccr."),
    colnames(to_bind_c)
  )
  
  id <- intersect(preD, postNULL)
  to_bind_A <- cbind(
    pre[id, c(beta_item, fdr_item)],
    post[id, c(beta_item, fdr_item)]
  )
  
  id <- intersect(preE, postNULL)
  to_bind_A <- rbind(
    to_bind_A,
    cbind(
      pre[id, c(beta_item, fdr_item)],
      post[id, c(beta_item, fdr_item)]
    )
  )
  colnames(to_bind) <- paste0(
    c("", "", "ccr.", "ccr."),
    colnames(to_bind)
  )
  
  return(list(
    `GW_impact %` = IMPACTEDg,
    `Phenotype_G_impact %` = IMPACTED_phenGenes,
    `Depleted_G_impact %` = IMPACTED_Depletions,
    `Enriched_G_impact %` = IMPACTED_Enrichments,
    `GW_distortion %` = DISTORTEDg,
    `Phenotype_G_distortion %` = DISTORTED_phenGenes,
    `Depleted_G_distortion %` = DISTORTED_Depletions,
    `Enriched_G_distortion %` = DISTORTED_Enrichments,
    geneCounts = geneCounts,
    distortion = to_bind,
    null = to_bind_c,
    impact = to_bind_A
  ))
}


ccr.fixedFDRthreshold <- function(
  FCsprofile,
  TruePositives,
  TrueNegatives,
  th
) {
  presentGenes <- intersect(c(TruePositives, TrueNegatives), names(FCsprofile))
  predictions <- FCsprofile[presentGenes]
  observations <- is.element(presentGenes, TruePositives) + 0
  names(observations) <- presentGenes
  RES <- roc(observations, predictions, direction = ">", quiet = TRUE)
  COORS <- coords(RES, "all", ret = c("threshold", "ppv"), transpose = TRUE)
  FDRpercTh <- max(COORS["threshold", which(COORS["ppv", ] >= (1 - th))])
  FDRpercRANK <- max(which(sort(FCsprofile) <= FDRpercTh))
  return(list(FCth = FDRpercTh, RANK = FDRpercRANK))
}


ccr.custom_VisDepAndSig <- function(
  FCsprofile,
  SIGNATURES,
  pIs = NULL,
  nIs = NULL,
  th = 0.05
) {

  sigNames <- names(SIGNATURES)

  nsig <- length(SIGNATURES)

  if (length(pIs) > 0 & length(nIs) > 0) {
    RES <- ccr.fixedFDRthreshold(
      FCsprofile = FCsprofile,
      TruePositives = SIGNATURES[[pIs]],
      TrueNegatives = SIGNATURES[[nIs]],
      th = th
    )
    FDR5percRANK <- RES[["RANK"]]
  } else {
    FDR5percRANK <- NULL
  }

  nelements <- length(FCsprofile)
  TPR <- vector()

  for (i in seq_along(SIGNATURES)) {
    hitPositions <- match(SIGNATURES[[i]], names(sort(FCsprofile)))
    hitPositions <- hitPositions[!is.na(hitPositions)]

    TPR[i] <- (
      length(which(hitPositions <= FDR5percRANK)) /
      length(hitPositions)
    )
  }
  names(TPR) <- sigNames

  return(TPR)
}
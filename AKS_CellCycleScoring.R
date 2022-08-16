AKS_CellCycle<-function (object, g1s.features, s.features, g2m.features, m.features, mg1.features, ctrl = NULL, set.ident = FALSE, 
                         ...) 
{
  name <- "AKS_CellCycle"
  
  features <- list(G1S.Score=g1s.features, S.Score = s.features, G2M.Score = g2m.features, M.Score = m.features,MG1.Score=mg1.features)
  if (is.null(x = ctrl)) {
    ctrl <- min(vapply(X = features, FUN = length, FUN.VALUE = numeric(length = 1)))
  }
  object.cc <- AddModuleScore(slot="raw.data",object = object, features = features, 
                              name = name, ctrl = ctrl, ...)
  cc.columns <- grep(pattern = name, x = colnames(x = object.cc[[]]), value = TRUE)
  cc.scores <- object.cc[[cc.columns]]
  rm(object.cc)
  CheckGC()
  assignments <- apply(X = cc.scores, MARGIN = 1, FUN = function(scores,
                                                                 first = "G1S", 
                                                                 second = "S", 
                                                                 third = "G2M",
                                                                 fourth = "M",
                                                                 fifth = "MG1",
                                                                 null = "G0") {
    if (all(scores < 1)) {
      return(null)
    }
    else {
      if (length(which(x = scores == max(scores))) > 1) {
        return("Undecided")
      }
      else {
        return(c(first,second,third,fourth,fifth)[which(x = scores == max(scores))])
      }
    }
  })
  cc.scores <- merge(x = cc.scores, y = data.frame(assignments), by = 0)
  colnames(x = cc.scores) <- c("rownames","G1S.Score","S.Score", "G2M.Score","M.Score","MG1.Score", "Phase")
  rownames(x = cc.scores) <- cc.scores$rownames
  cc.scores <- cc.scores[, c("G1S.Score","S.Score","G2M.Score","M.Score","MG1.Score", "Phase")] #"MG1.Score",
  object[[colnames(x = cc.scores)]] <- cc.scores
  if (set.ident) {
    object[["old.ident"]] <- Idents(object = object)
    Idents(object = object) <- "Phase"
  }
  return(object)
}

CheckGC <- function() {
  if (getOption(x = "Seurat.memsafe")) {
    gc(verbose = FALSE)
  }
}

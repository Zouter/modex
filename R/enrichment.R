# estimating jaccard using matrix multiplication
# see http://stats.stackexchange.com/questions/49453/calculating-jaccard-or-other-association-coefficient-for-binary-data-using-matrix
jaccard <- function(m) {
  A <- Matrix::crossprod(m)
  im <- Matrix::which(A > 0, arr.ind=TRUE)
  b <- Matrix::colSums(m)

  # only calculate for non zero overlaps
  Aim <- A[im]

  J <- Matrix::sparseMatrix(
    i = im[,1],
    j = im[,2],
    x = Aim / (b[im[,1]] + b[im[,2]] - Aim),
    dims = dim(A)
  )
  return(J)
}

#' @export
filterGsets <- function(gsets, background, minsize=5, maxsize=500, maxoverlap=0.5, verbose=F) {
  gsets <- lapply(gsets, function(gset) unique(intersect(gset, background)))
  gsets <- Filter(function(x) between(length(x), minsize, maxsize), gsets)

  genemap <- setNames(1:length(background), background) # map genes to column id

  gsetorder <- order(unlist(lapply(gsets, length)))
  gsets <- gsets[gsetorder]
  i <- unlist(lapply(gsets, function(gset) as.numeric(genemap[gset])))
  j <- unlist(lapply(1:length(gsets), function(i) rep(i, length(gsets[[i]]))))
  membership <- Matrix::sparseMatrix(i, j, x=1) # not boolean, for matrix multiplication

  jac <- jaccard(membership)
  retained <- c(1)

  if(verbose) pb<-utils::txtProgressBar(max=ncol(jac))

  for (i in 2:ncol(jac)) {
    if(max(jac[i, retained]) < maxoverlap) {
      retained <- c(retained, i)
    }
    if((i %% 100 == 0) && verbose) utils::setTxtProgressBar(pb, i)
  }

  if(verbose) close(pb)

  gsetids <- names(gsets)[retained]
  gsets[gsetids]
}


#' Simple test for enrichment in a given set of genes
#'
#' Assumes gsets is already filtered for background, will not check this!
#'
#' @param module character vector containing the genes tested for enrichment
#' @param gsets a list containing lists of genes
#' @param background character vector
#' @return Dataframe containing for every gene set which passed the filtering: \itemize{
#'  \item p-value of enrichment
#'  \item q-value of enrichment (p-value corrected for multiple testing)
#'  \item estimated odds score of enrichment, how much more likely is a given gene set to be part of
#' }
#' @export
testEnrichment = function(module, gsets, background) {
  scores = do.call(rbind.data.frame, lapply(names(gsets), function(gsetid){
    gset = gsets[[gsetid]]

    # much faster than using table
    tp = length(intersect(module, gset))
    fn = length(module) - tp
    fp = length(gset) - tp
    tn = length(background) - tp - fp - fn

    pval = phyper(tp, length(module), length(background) - length(module), length(gset), lower.tail=F)

    #if(tp == 0) return(list(pval=1, odds=1, found=0, gsetid=gsetid))

    contingency_table = matrix(c(tp, fp, fn, tn), 2, 2)

    fisher_result = fisher.test(contingency_table, alternative="greater", conf.int = F)
    list(pval=fisher_result$p.value, odds=fisher_result$estimate, found=contingency_table[1,1], gsetid=gsetid)
    #list(pval=pval, odds=(tp*tn)/(fn*fp), found=tp, gsetid=gsetid)
  }))
  if (nrow(scores) > 0) {
    scores$qval = p.adjust(scores$pval, method="bonferroni")
  }
  scores
}

#' @export
getAucodds <- function(modules, gsets_filtered, background, qvalcutoff=0.05, oddscutoffs = 10^seq(0, 2, length.out=100)) {
  if(length(modules) == 0) {
    return(list(aucodds=0))
  }
  scores = bind_rows(lapply(1:length(modules), function(moduleid) {
    scores = testEnrichment(modules[[moduleid]], gsets_filtered, background)
    scores$moduleid = moduleid
    scores
  }))
  pvals = aacast(scores, moduleid~gsetid, "pval")
  qvals = aacast(scores, moduleid~gsetid, "qval")
  #qvals = matrix(p.adjust(pvals, method="fdr"), ncol=length(modules), byrow = F, dimnames = dimnames(pvals))
  odds = aacast(scores, moduleid~gsetid, "odds")

  newodds = odds
  newodds[qvals > qvalcutoff] = 0

  bestodds = apply(newodds, 1, max)
  stillenriched = unlist(lapply(oddscutoffs, function(cutoff) mean(bestodds > cutoff)))
  aucodds = (stillenriched[[1]] + 2*sum(stillenriched[2:length(stillenriched)-1]) + stillenriched[[length(stillenriched)]])/(2*length(stillenriched))

  bestodds = apply(newodds, 2, max)
  stillenriched = unlist(lapply(oddscutoffs, function(cutoff) mean(bestodds > cutoff)))
  aucodds2 = (stillenriched[[1]] + 2*sum(stillenriched[2:length(stillenriched)-1]) + stillenriched[[length(stillenriched)]])/(2*length(stillenriched))

  return(list(aucodds=aucodds*aucodds2, stillenriched=stillenriched, newodds=newodds, scores=scores))
}

#' @export
testGridResults <- function(results, gsets_filtered, background, parallel="qsub", qsub.conf=PRISM::qsub.configuration()) {
  if (is.null(parallel)) {
    mc.cores = 1
    parallel = "multicore"
  }

  loop_func = function(i) {
    modules = results$modules[[i]]
    params = results$params[i,,drop=F]

    scores = modex::getAucodds(modules, gsets_filtered, background)

    data.frame(aucodds=scores$aucodds,  params, stringsAsFactors = F)
  }

  if(parallel == "multicore") {
    return(dplyr::bind_rows(mclapply(1:length(results$modules), loop_func, mc.cores = getOption("mc.cores", default=1))))
  } else if(parallel == "qsub") {
    return(dplyr::bind_rows(PRISM::qsub.lapply(1:length(results$modules), loop_func, qsub.config=qsub.conf)))
  } else {
    stop("wrong parallel parameter!")
  }
}

#' @export
testGrid <- function(E, gsets, method=agglom, filterGsets=T, ..., verbose=F, background=rownames(E)) {
  if (filterGsets) {
    if (verbose) message("Filtering gene sets")
    gsets = filterGsets(gsets, background, verbose = T)
    if (verbose) message("Gene sets filtered")
  }
  results = do.call(method, c(list(E=E), list(...)))
  if (verbose) message("Clustered")
  if (verbose) message("Calculating enrichment...")
  scores= testGridResults(results, gsets, background)
  if (verbose) message("Enrichment calculated")
  list(scores=scores, results=results)
}

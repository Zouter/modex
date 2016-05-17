agglom <- function(E, k=c(1,2,3,4,5,6,7,8,9,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50,52,54,56,58,60,62,64,66,68,70,72,74,76), agglom.method=c("ward.D2", "single", "complete", "average")) {
  # although this can look quite cumbersome, this kind of setup is easy to expand to multiple parameters and multiple commands which have to be run progressively
  treeparams = expand_grid(agglom.method=agglom.method)

  distances = as.dist(1-cor(t(E)))

  trees = mcdapply(treeparams, function(params) {
    stats::hclust(distances, method = params$agglom.method)
  })

  cutparams = expand_grid(k=k, agglom.method=agglom.method)
  modules = mcdapply(cutparams, function(params) {
    tree = trees[[matchparams(treeparams, params)]]
    labels2modules(stats::cutree(tree, params$k))
  })
  list(params=cutparams, modules=modules)
}

affinity <- function(E, preference_fraction=c(-2, -1.5, -1, 0.5, 0, 0.2, 0.4, 0.6)) {
  similarities = cor(t(E))

  simmax = max(similarities)
  simmin = min(similarities)

  params = expand_grid(preference_fraction=preference_fraction)
  modules = mcdapply(params, function(params) {
    preference = simmin + (simmax - simmin)*params$preference_fraction
    labels2modules(apcluster::labels(apcluster::apcluster(similarities, p=preference), type="enum"), rownames(E))
  })
  list(params=params, modules=modules)
}

ica_fdr <- function(E, k=c(1,2,3,4,5,6,7,8,9,10), qvalcutoff=c(0.5, 0.2, 0.1, 0.05, 0.02, 0.01), icarepeat=3) {
  library(fastICA)
  library(fdrtool)

  icaparams = expand_grid(k=k, icarepeat=seq(icarepeat))

  sources = mcdapply(icaparams, function(params) {
    result = fastICA::fastICA(E, params$k, row.norm=T)
    result$S
  })

  cutoffparams = expand_grid(k=k, qvalcutoff=qvalcutoff, icarepeat=seq(icarepeat))
  modules = mcdapply(cutoffparams, function(params) {
    S = sources[[matchparams(icaparams, params)]]
    modules = as.list(apply(S, 2, function(s) {
      qvals = fdrtool::fdrtool(s, plot=F, verbose=F, cutoff.method="fndr")$qval
      names(qvals)[(qvals < params$qvalcutoff) & (s>0)]
    }))
    modules = c(modules, as.list(apply(S, 2, function(s) {
      qvals = fdrtool::fdrtool(s, plot=F, verbose=F, cutoff.method="fndr")$qval
      names(qvals)[(qvals < params$qvalcutoff) & (s<0)]
    })))
    modules = Filter(function(module) length(module) > 0, modules)
    modules
  })
  list(params=cutoffparams, modules=modules)
}

wgcna <- function(E, power=c(1,2,3,4,5,6,7,8,9,10), mergeCutHeight=c(0, 0.1, 0.2, 0.3, 0.4, 0.5)) {
  library(WGCNA)

  params = expand_grid(power=power, mergeCutHeight=mergeCutHeight)
  modules = mcdapply(params, function(params) {
    result = WGCNA::blockwiseModules(t(E), maxBlockSize = 10000, networkType="signed", power=params$power, mergeCutHeight = params$mergeCutHeight)
    labels2modules(result$colors, rownames(E), ignore_label="grey")
  })
  list(params=params, modules=modules)
}

kmedoids <- function(E, k=c(2,3,4,5,6,7,8,9,10,11,12,13,14,15), kmedoidsrepeat=1) {
  distances = (1-cor(t(E)))

  params = expand_grid(k=k, kmedoidsrepeat=seq(kmedoidsrepeat))

  modules = mcdapply(params, function(params) {
    result = cluster::pam(distances, diss=T, k=params$k)
    labels2modules(result$clustering, rownames(E))
  })
  list(params=params, modules=modules)
}

process_time <- proc.time()[3]

initial_dir <- getwd()

source("otherScripts/R/test.phylosignal.R")
source("otherScripts/R/runIQtree.R")
source("otherScripts/R/clean.gene.R")
source("otherScripts/R/get.model.R")
source("otherScripts/R/sample.quartet.R")
source("testStatistics/get.phylosignal.metrics.R")
setwd("otherScripts/phybase/")
for (i in dir()) source(i)
setwd(initial_dir)
print("Functions required were loaded successfully")

selected_stats <- unlist(input$testStats)
outs <- list()
locilengths <- vector()

if (input$testType == "tree") {
  first_line <- readLines(as.character(input$sptreePath[1, 4]), n = 1)

  if (grepl("[#]NEXUS|[#]nexus", first_line)) {
    trees_format <- "nexus"
    sptreraw <- read.nexus(as.character(input$sptreePath[1, 4]))
  } else if (grepl("[(]", first_line)) {
    trees_format <- "newick"
    sptreraw <- read.tree(as.character(input$sptreePath[1, 4]))
  }

  if (is.rooted(sptreraw)) sptre <- unroot(sptreraw) else sptre <- sptreraw
  # List of edges, each with with a list of n_sims quartets
  allsampquarts <- lapply(
    which(!sptre$edge[, 2] %in% 1:Ntip(sptre)),
    function(x) lapply(1:input$n_sims, function(y) sample.quartet(sptre, x))
  )
  # Matrices where columns are edges and rows are replicates
  medmat <- matrix(NA, ncol = length(allsampquarts), nrow = input$n_sims * nrow(input$dataPath))
  pvalmat <- matrix(NA, ncol = length(allsampquarts), nrow = input$n_sims * nrow(input$dataPath))
  presabs <- matrix(NA, ncol = length(allsampquarts), nrow = input$n_sims * nrow(input$dataPath))
  colnames(medmat) <- colnames(pvalmat) <- colnames(presabs) <- names(allsampquarts) <- which(!sptre$edge[, 2] %in% 1:Ntip(sptre))
  rownames(medmat) <- rownames(pvalmat) <- rownames(presabs) <- as.character(
    sapply(
      seq_along(nrow(input$dataPath)),
      function(x) paste0("rep.", 1:input$n_sims, ".locus.", x)
    )
  )
}

for (j in seq_along(nrow(input$dataPath))) {
  first_line <- readLines(as.character(input$dataPath[j, 4]), n = 1)
  if (grepl("[>]", first_line)) {
    data_format <- "fasta"
  } else if (grepl("[#]NEXUS|[#]nexus", first_line)) {
    data_format <- "nexus"
  } else if (grepl("[(]", first_line)) {
    data_format <- "newick"
  } else {
    data_format <- "phylip"
  }

  model <- "GTR+G"

  print("Model to be assessed was identified")

  if (input$testType %in% c("locus", "tree")) {
    analysisdata <- clean.gene(
      sdata = as.character(input$dataPath[j, 4]),
      format = data_format,
      aadata = input$data_type,
      clean = FALSE
    )
  } else {
    analysisdata <- as.character(input$dataPath[j, 4])
  }

  print("Locus was cleaned successfully")

  setwd(input$outputFolder)

  print("Output folder was identified successfully")

  what_to_output <- unlist(input$what_to_output)

  if (input$testType == "locus") {
    if (!input$overwrite && file.exists(paste0(as.character(input$dataPath[j, 1]), ".phylomad.phylosig"))) {
      stop("Exisitng files will not be overwritten. Aborting.")
    } else {
      system(
        paste0(
          "mkdir ", as.character(input$dataPath[j, 1]), ".phylomad.phylosig"
        )
      )
      setwd(
        paste0(
          as.character(input$dataPath[j, 1]), ".phylomad.phylosig"
        )
      )
    }

    gene_results <- try(
      test.phylosignal(
        sdata = analysisdata,
        format = if (input$testType == "locus") "bin" else data_format,
        testType = input$testType,
        aadata = input$data_type,
        model = model,
        iqtreePath = iqtreePath,
        astralPath = astralPath,
        n_sims = input$n_sims,
        testStats = selected_stats,
        returnSimulations = "simdat" %in% what_to_output
      )
    )
    if (class(gene_results) == "try-error") {
      setwd("..")
      system(paste0("rm -r ", as.character(input$dataPath[j, 1]), ".phylomad.phylosig"))
      print(paste0("Assessment of ", as.character(input$dataPath[j, 1]), " failed"))
      next
    }

    rownames(gene_results[[1]]) <- gene_results[[1]][, 1]
    colnames(gene_results[[1]]) <- gsub("Length", "br.length", colnames(gene_results[[1]]))
    colnames(gene_results[[1]]) <- gsub("Label", "br.support", colnames(gene_results[[1]]))
    for (y in 2:ncol(gene_results[[1]])) {
      gene_results[[1]][, y] <- round(as.numeric(gene_results[[1]][, y]), 3)
    }
    gene_results[[1]] <- rbind(
      gene_results[[1]],
      round(colMeans(gene_results[[1]], na.rm = TRUE), 3)
    )
    rownames(gene_results[[1]])[nrow(gene_results[[1]])] <- "mean"

    ## Annotate empirical tree with bootstrap support ranges and p-valuea
    if (input$testType == "locus" && "phyloempres" %in% what_to_output) {
      nexhead <- c("#NEXUS", "begin trees;")
      antre <- gene_results[[2]]
      antre$node.label <- vector()
      for (y in seq_len(nrow(gene_results[[1]] - 1))) {
        annot <- ""
        for (z in seq_along(length(selected_stats))) {
          if (selected_stats[z] == "CF") {
            fullstat <- as.numeric(gene_results[[1]][y, paste0("sCF.sim.", 1:input$n_sims)])
          } else {
            fullstat <- as.numeric(
              gene_results[[1]][
                y,
                paste0("sim.", 1:input$n_sims, ".", selected_stats[z])
              ]
            )
          }

          statmed <- round(median(fullstat, na.rm = TRUE), 2)
          statran <- round(quantile(fullstat, probs = c(0.025, 0.975), na.rm = TRUE), 2)
          annot <- c(
            annot,
            paste0(
              selected_stats[z], '=\"',
              statmed,
              " [", statran[1],
              ", ", statran[2],
              "] ovcon.p=", gene_results[[1]][
                y, paste0(selected_stats[z], ".p.value")
              ],
              '\"'
            )
          )
        }
        annot <- paste0("<&!", paste0(annot, collapse = ","), ">")
        antre$node.label[
          antre$edge[
            which(gene_results[[3]]$edge.length == gene_results[[1]][y, 1]), 2
          ] - Ntip(antre)
        ] <- annot
      }
      antre <- gsub("[_]", " ", write.tree(antre))
      antre <- gsub("[-]", ",", antre)
      antre <- gsub("[[]", "(", antre)
      antre <- gsub("[]]", ")", antre)
      antre <- gsub("[<]", "[", antre)
      antre <- gsub("[>]", "]", antre)
      antre <- c(nexhead, paste0("tree tree_1 = [&R] ", antre), "end;")
      writeLines(antre, con = "estimated.tree.supports.tre")
    }

    if ("allqp" %in% what_to_output) {
      write.csv(t(gene_results[[1]]), file = "full.results.csv")
    }

    if ("pvals" %in% what_to_output) {
      restab <- gene_results[[1]][, -grep("ID|sDF|sN|gDF|gN", colnames(gene_results[[1]]))]
      write.csv(t(restab), file = "results.summary.csv")
    }

    if ("testPlots" %in% what_to_output) {
      histplotdat <- as.numeric(gene_results[[1]]["mean", ])
      names(histplotdat) <- colnames(gene_results[[1]])
      histplotdat <- histplotdat[-grep("ID|sDF|sN|gDF|gN|value|sdpd", names(histplotdat))]
      statlabels <- vector()
      if ("dnet" %in% selected_stats) statlabels <- c(statlabels, "Distance to network")
      if ("dtree" %in% selected_stats) statlabels <- c(statlabels, "Distance to tree")
      if ("entrop" %in% selected_stats) statlabels <- c(statlabels, "Entropy")
      if ("icert" %in% selected_stats) statlabels <- c(statlabels, "Internode certainty")
      if ("CF" %in% selected_stats) statlabels <- c(statlabels, "Concordance factor")
      if ("binp" %in% selected_stats) {
        statlabels <- c(statlabels, "Binomial P of concordance factor")
      }
      if ("dstat" %in% selected_stats) statlabels <- c(statlabels, "D-statistic")
      if ("kcstat" %in% selected_stats) statlabels <- c(statlabels, "Chifman-Kubatko statistic")

      if (length(selected_stats) == 0) {
        print("Test plots cannot be returned because no test statistics were calculated.")
      } else {
        gene_results[[1]] <- as.matrix(gene_results[[1]])
        pdf("tests.density.plots.pdf", useDingbats = FALSE, height = 4, width = 15)
        par(mfrow = c(1, 3))
        for (i in seq_along(length(selected_stats))) {
          statdat <- histplotdat[grep(selected_stats[i], names(histplotdat))]
          if (is.na(statdat[1]) || is.nan(statdat[1]) || any(is.infinite(statdat)) || all(is.na(statdat)) || all(is.nan(statdat))) {
            print(paste0("Plots of ", statlabels[i], " cannot be created."))
            next
          }

          # Plot of all simulated branch values vs. all empirical branch values
          statsim <- gene_results[[1]][
            -nrow(gene_results[[1]]),
            grep(selected_stats[i], colnames(gene_results[[1]]))
          ]
          statsim <- as.numeric(statsim[, 2:(input$n_sims + 1)])

          statemp <- as.numeric(
            gene_results[[1]][
              -nrow(gene_results[[1]]),
              if (selected_stats[i] == "CF") 2 else paste0("emp.", selected_stats[i])
            ]
          )

          statsimdens <- try(
            density(
              statsim,
              from = min(statsim, na.rm = TRUE),
              to = max(statsim, na.rm = TRUE),
              na.rm = TRUE
            )
          )
          statempdens <- try(
            density(
              statemp,
              from = min(statemp, na.rm = TRUE),
              to = max(statemp, na.rm = TRUE),
              na.rm = TRUE
            )
          )

          if (class(statsimdens) != "try-error" && class(statempdens) != "try-error") {
            plot(
              statsimdens,
              xlim = range(c(statsim, statemp), na.rm = TRUE),
              ylim = c(0, max(c(statsimdens$y, statempdens$y), na.rm = TRUE)),
              xlab = statlabels[i],
              main = ""
            )
            lines(statempdens, col = "red")
            legend(
              "topright",
              legend = c("Empirical branches", "Simlulated branches"),
              col = c("red", "black"),
              lty = 1
            )
          } else {
            frame()
          }

          # Plot of branch-wise p-values
          statpvals <- as.numeric(
            gene_results[[1]][
              -nrow(gene_results[[1]]),
              paste0(selected_stats[i], ".p.value")
            ]
          )
          pdens <- try(
            density(
              statpvals,
              from = min(statpvals, na.rm = TRUE),
              to = max(statpvals, na.rm = TRUE),
              na.rm = TRUE
            )
          )
          if (class(pdens) != "try-error") {
            plot(pdens, main = "", xlab = paste0("Branch P-values\n(", statlabels[i], ")"))
          } else {
            frame()
          }

          # Plot of test of mean across branches
          statmeandens <- try(
            density(
              statdat[-1],
              from = min(statdat[-1], na.rm = TRUE),
              to = max(statdat[-1], na.rm = TRUE),
              na.rm = TRUE
            )
          )
          if (class(statmeandens) != "try-error") {
            plot(
              statmeandens,
              main = "",
              xlim = range(statdat, na.rm = TRUE),
              xlab = paste0("Mean ", statlabels[i], " across branches")
            )
            abline(v = statdat[1], lwd = 2, col = "red")
            legend(
              "topright",
              legend = c("Empirical", "Simulated"),
              lty = 1,
              col = c("red", "black")
            )
          } else {
            frame()
          }
        }
        dev.off()

        pdf(
          "tests.summary.tree.pdf",
          useDingbats = FALSE,
          height = if (length(gene_results[[2]]$edge.length) < 50) 5 else 10,
          width = if (length(gene_results[[2]]$edge.length) < 50) 10 else 20
        )
        for (i in seq_along(length(selected_stats))) {
          par(mfrow = c(1, 2), mar = c(5.1, 4.1, 4.1, 2.1))
          tr <- gene_results[[2]]
          tr$edge.length <- rep(1, length(gene_results[[2]]$edge.length))
          brpvalue <- gene_results[[1]][, paste0(selected_stats[i], ".p.value")]
          brpvalue[which(brpvalue > 0.01)] <- 1
          brpvalue[which(brpvalue <= 0.01)] <- 2
          brsdpd <- round(gene_results[[1]][, paste0(selected_stats[i], ".sdpd")], 1)
          names(brpvalue) <- names(brsdpd) <- gene_results[[1]][, 1]

          brpvalue <- brpvalue[as.character(gene_results[[3]]$edge.length)]
          brsdpd <- brsdpd[as.character(gene_results[[3]]$edge.length)]
          brpvalue[is.na(brpvalue)] <- 1
          brsdpd[is.na(brsdpd)] <- mean(brsdpd, na.rm = TRUE)
          if (all(is.na(brpvalue)) || all(is.nan(brpvalue)) || any(is.infinite(brpvalue)) || all(is.na(brsdpd)) || all(is.nan(brsdpd)) || any(is.infinite(brsdpd))) {
            print(paste0("Tree depicting ", statlabels[i], " statistic cannot be ploted."))
            next
          }

          plot(
            gene_results[[2]],
            edge.color = brpvalue,
            type = "unrooted",
            main = paste0(
              statlabels[i],
              ". In red P-values < 0.01\n(includes branch lengths)"
            )
          )
          plot(tr, edge.color = brpvalue, type = "unrooted", main = "\n(excludes branch lengths)")
          brpvalcols <- brpvalue
          brpvalcols[brpvalcols == 1] <- "white"
          brpvalcols[brpvalcols == "2"] <- "red"
          edgelabels(names(brpvalue), frame = "circle", bg = brpvalcols, cex = 0.7)

          plotBranchbyTrait(
            gene_results[[2]],
            brsdpd,
            mode = "edges",
            palette = "rainbow",
            type = "unrooted",
            legend = FALSE
          )
          plotBranchbyTrait(
            tr,
            brsdpd,
            mode = "edges",
            palette = "rainbow",
            type = "unrooted",
            title = paste0(statlabels[i], "\nz-score\n")
          )
        }
        dev.off()

        # finish up and add custom legend?
        # Add two ternary plots one with all the simulated branches,
        # and one of the mean cfs for each data set, put both in one page with grid.arrange
        pdf("ternary.plots.pdf", height = 4, width = 8)
        # gene_results[[1]] <- as.data.frame(gene_results[[1]])
        if (input$testType == "locus") ttyp <- "s" else ttyp <- "g"
        terndat <- gene_results[[1]][-nrow(gene_results[[1]]), paste0(ttyp, c("CF", "DF1", "DF2"))]
        ternmeandat <- colMeans(terndat)
        nbranch <- nrow(terndat)
        nsimbranch <- nbranch * input$n_sims
        for (i in 1:input$n_sims) {
          terndat <- rbind(
            gene_results[[1]][
              -nrow(gene_results[[1]]),
              c(
                paste0(ttyp, "CF.sim.", i),
                paste0(ttyp, "DF1.sim.", i),
                paste0(ttyp, "DF2.sim.", i)
              )
            ],
            terndat
          )
          ternmeandat <- rbind(
            gene_results[[1]][
              nrow(gene_results[[1]]),
              c(
                paste0(ttyp, "CF.sim.", i),
                paste0(ttyp, "DF1.sim.", i),
                paste0(ttyp, "DF2.sim.", i)
              )
            ],
            ternmeandat
          )
        }
        colnames(terndat) <- colnames(ternmeandat) <- c("CF", "DF1", "DF2")

        terndat <- as.data.frame(terndat)
        # terndat$cols <- c(rep(2, nsimbranch), rep(1, nbranch))
        terndat$cols <- as.factor(rep(1:nbranch, input$n_sims + 1))
        ternallbrs <- ggtern(
          data = terndat,
          aes(x = DF1, y = CF, z = DF2),
          aes(x, y, z)
        ) +
          geom_point(
            aes(fill = cols),
            alpha = c(rep(0.25, nsimbranch), rep(1, nbranch)),
            stroke = 0, size = 2,
            shape = c(rep(21, nsimbranch), rep(23, nbranch))
          ) +
          theme(legend.position = "none") +
          ggtitle("All simulated and empirical\nbranches")

        ternmeandat <- as.data.frame(ternmeandat)
        ternmeandat$cols <- c(rep(2, input$n_sims), 1)
        ternmean <- ggtern(
          data = ternmeandat, aes(x = DF1, y = CF, z = DF2), aes(x, y, z)
        ) +
          geom_point(
            aes(fill = cols),
            alpha = c(rep(0.5, input$n_sims), 1),
            stroke = 0, size = 2,
            shape = c(rep(21, input$n_sims), 23)
          ) +
          theme(legend.position = "none") +
          ggtitle("Mean across each simulation\nand empirical tree")

        grid.arrange(ternallbrs, ternmean, nrow = 1)
        dev.off()
      }
    }

    setwd("..")
  } else if (input$testType == "tree") {
    if (!input$overwrite && file.exists(paste0(as.character(input$sptreePath[1, 1]), ".phylomad.supports.tre"))) {
      stop("Exisitng files will not be overwritten. Aborting.")
    }

    if (input$Ncores == 1) {
      for (x in seq_along(length(allsampquarts))) {
        allqtax <- sapply(
          allsampquarts[[x]],
          function(y) all(y$tip.label %in% rownames(analysisdata))
        )
        if (!all(allqtax)) {
          print("Quartet missing from tree")
          next
        }
        for (z in which(allqtax)) {
          quartloc <- analysisdata[allsampquarts[[x]][[z]]$tip.label, ]
          gene_result <- try(
            test.phylosignal(
              sdata = quartloc,
              format = "bin",
              testType = "locus",
              aadata = input$data_type,
              model = model,
              iqtreePath = iqtreePath,
              astralPath = astralPath,
              n_sims = input$n_sims,
              testStats = "CF",
              returnSimulations = FALSE
            )
          )
          fullstat <- as.numeric(gene_result[[1]][1, paste0("sCF.sim.", 1:input$n_sims)])
          medmat[z + ((j - 1) * input$n_sims), x] <- round(median(fullstat, na.rm = TRUE), 3)
          pvalmat[z + ((j - 1) * input$n_sims), x] <- round(gene_result[[1]][1, "CF.p.value"], 3)
          presabs[z + ((j - 1) * input$n_sims), x] <- (RF.dist(allsampquarts[[x]][[z]], gene_result[[2]]) != 0)
        }
        pvalmat[which(allqtax) + ((j - 1) * input$n_sims), x] <- round(
          p.adjust(
            pvalmat[which(allqtax) + ((j - 1) * input$n_sims), x],
            method = "fdr"
          ),
          3
        )
      }
    } else {
      ### START PARALLEL COMPUTING
      print("Parallel computing started")
      system(paste0("mkdir processing.locus.", j))
      setwd(paste0("processing.locus.", j))
      require(foreach)
      require(doParallel)
      data_type <- input$data_type
      n_sims <- input$n_sims
      run_sim <- function(x) {
        system(paste0("mkdir branchfolder.", x))
        setwd(paste0("branchfolder.", x))
        allqtax <- sapply(
          allsampquarts[[x]],
          function(y) all(y$tip.label %in% rownames(analysisdata))
        )
        if (!all(allqtax)) {
          print("Quartet missing from tree")
          next
        }
        medmatdat <- matrix(NA, ncol = 1, nrow = length(allsampquarts[[x]]))
        pvalmatdat <- matrix(NA, ncol = 1, nrow = length(allsampquarts[[x]]))
        presabsdat <- matrix(NA, ncol = 1, nrow = length(allsampquarts[[x]]))
        for (z in which(allqtax)) {
          quartloc <- analysisdata[allsampquarts[[x]][[z]]$tip.label, ]
          gene_result <- try(
            test.phylosignal(
              sdata = quartloc,
              format = "bin",
              testType = "locus",
              aadata = data_type,
              model = model,
              iqtreePath = iqtreePath,
              astralPath = astralPath,
              n_sims = n_sims,
              testStats = "CF",
              returnSimulations = FALSE
            )
          )
          fullstat <- as.numeric(gene_result[[1]][1, paste0("sCF.sim.", 1:n_sims)])
          medmatdat[z, 1] <- round(median(fullstat, na.rm = TRUE), 3)
          pvalmatdat[z, 1] <- round(gene_result[[1]][1, "CF.p.value"], 3)
          presabsdat[z, 1] <- (RF.dist(allsampquarts[[x]][[z]], gene_result[[2]]) != 0)
        }
        pvalmatdat[which(allqtax), 1] <- round(
          p.adjust(
            pvalmatdat[which(allqtax), 1],
            method = "fdr"
          ),
          3
        )
        setwd("..")
        system(paste0("rm -r branchfolder.", x))
        res <- list(medmatdat, pvalmatdat, presabsdat)
        return(res)
      }

      cl <- makeCluster(input$Ncores)
      registerDoParallel(cl)
      sim_reps <- foreach(
        x = seq_along(length(allsampquarts)),
        .packages = c("phangorn", "ape"),
        .export = c("test.phylosignal", "runIQtree")
      ) %dopar% run_sim(x)

      for (x in seq_along(length(allsampquarts))) {
        medmat[(1:input$n_sims) + ((j - 1) * input$n_sims), x] <- sim_reps[[x]][[1]]
        pvalmat[(1:input$n_sims) + ((j - 1) * input$n_sims), x] <- sim_reps[[x]][[2]]
        presabs[(1:input$n_sims) + ((j - 1) * input$n_sims), x] <- sim_reps[[x]][[3]]
      }

      stopCluster(cl)
      setwd("..")
      system(paste0("rm -r processing.locus.", j))
      print("Parallel computing ended successfully")
      ### END PARALLEL COMPUTING
    }


    if ("pvals" %in% what_to_output) {
      write.csv(
        medmat,
        file = paste0(
          as.character(input$sptreePath[1, 1]),
          ".phylomad.median.CFs.csv"
        )
      )
      write.csv(
        pvalmat,
        file = paste0(
          as.character(input$sptreePath[1, 1]),
          ".phylomad.median.pvals.csv"
        )
      )
      write.csv(
        presabs,
        file = paste0(
          as.character(input$sptreePath[1, 1]),
          ".phylomad.median.presence.csv"
        )
      )
    }
    if ("phyloempres" %in% what_to_output) {
      nexhead <- c("#NEXUS", "begin trees;")
      annot <- sapply(seq_along(length(allsampquarts)), function(x) {
        medmed <- round(median(medmat[, x], na.rm = TRUE), 2)
        ranmed <- round(quantile(medmat[, x], probs = c(0.025, 0.975), na.rm = TRUE), 2)
        medpval <- round(median(pvalmat[, x], na.rm = TRUE), 2)
        ranpval <- round(quantile(pvalmat[, x], probs = c(0.025, 0.975), na.rm = TRUE), 2)
        medabs <- round(mean(presabs[, x], na.rm = TRUE), 2)
        ranabs <- round(quantile(presabs[, x], probs = c(0.025, 0.975), na.rm = TRUE), 2)
        ann <- paste0(
          '<&!,CF=\"', medmed,
          " [", ranmed[1],
          ", ", ranmed[2],
          ']\",OverconP=\"', medpval,
          " [", ranpval[1],
          ", ", ranpval[2],
          ']\",Absence=\"', medabs,
          " [", ranabs[1],
          ", ", ranabs[2],
          ']\">'
        )
        return(ann)
      })
      antre <- sptre
      antre$node.label <- c(NA, annot)
      antre <- gsub("[_]", " ", write.tree(antre))
      antre <- gsub("[-]", ",", antre)
      antre <- gsub("[[]", "(", antre)
      antre <- gsub("[]]", ")", antre)
      antre <- gsub("[<]", "[", antre)
      antre <- gsub("[>]", "]", antre)
      antre <- c(nexhead, paste0("tree tree_1 = [&R] ", antre), "end;")
      writeLines(antre, con = paste0(as.character(input$sptreePath[1, 1]), ".phylomad.support.tre"))
    }
  }
}

setwd(initial_dir)

print(
  paste0(
    "Assessment completed successfully in ",
    round(proc.time()[3] - process_time, 3),
    " seconds."
  )
)

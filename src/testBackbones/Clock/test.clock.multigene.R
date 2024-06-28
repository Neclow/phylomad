process_time <- proc.time()[3]
source("otherScripts/R/run.gene.clock.R")
source("otherScripts/R/get.test.statistics.R")
source("otherScripts/R/runIQtree.R")
source("otherScripts/R/getRatogs.R")
source("testStatistics/get.df.R")
source("testStatistics/stemmy.R")

initial_dir <- getwd()

print("Functions required were loaded successfully")

trees <- read.nexus(as.character(input$trees_path[1, 4]))

selected_stats <- unlist(input$testStats)

outs <- list()

for (j in seq_along(nrow(input$dataPath))) {
  parallelise <- input$Ncores > 1

  setwd(input$outputFolder)

  print("Output folder was identified")

  if (input$overwrite == FALSE && file.exists(paste0(as.character(input$dataPath[j, 1]), ".phylomad.clock"))) {
    stop("Existing file will not be overwritten")
  } else {
    system(paste0("mkdir ", as.character(input$dataPath[j, 1]), ".phylomad.clock"))
    setwd(paste0(as.character(input$dataPath[j, 1]), ".phylomad.clock"))
  }

  if (nrow(input$trees_path) == 1) {
    trees_path <- as.character(input$trees_path[1, 4])
  } else {
    trees_path <- as.character(input$trees_path[j, 4])
  }

  if (nrow(input$posteriorPath) == 1) {
    post_path <- as.character(input$posteriorPath[1, 4])
  } else {
    post_path <- as.character(input$posteriorPath[j, 4])
  }

  first_line <- readLines(as.character(input$dataPath[j, 4]), n = 1)
  if (grepl("[>]", first_line)) {
    data_format <- "fasta"
  } else if (grepl("[#]NEXUS|[#]nexus", first_line)) {
    data_format <- "nexus"
  } else {
    data_format <- "phylip"
  }

  gene_results <- try(
    run.gene.clock(
      sdata = as.character(input$dataPath[j, 4]),
      format = data_format, treesFile = trees_path,
      logFile = post_path, burninpercentage = input$burnin,
      iqtreePath = iqtreePath,
      Nsims = input$Nsims,
      para = parallelise,
      ncore = input$Ncores,
      testStats = selected_stats,
      returnSimPhylo = TRUE,
      returnSimDat = TRUE
    )
  )

  if (class(gene_results) == "try-error") {
    setwd("..")
    system(paste0("rm -r ", as.character(input$dataPath[j, 1]), ".phylomad.clock"))
    print(paste("Analysis of locus", as.character(input$dataPath[j, 1]), "failed"))
    next
  }

  if ("pvals" %in% unlist(input$whatToOutput) || "simple" %in% unlist(input$whatToOutput)) {
    gene_res_mat <- matrix(NA, nrow = 3, ncol = length(selected_stats))
    for (i in seq_along(length(selected_stats))) {
      gene_res_mat[1, i] <- gene_results[[paste0(selected_stats[i], ".tailp")]]
      gene_res_mat[2, i] <- gene_results[[paste0("emp.", selected_stats[i])]]
      gene_res_mat[3, i] <- gene_results[[paste0(selected_stats[i], ".sdpd")]]
    }
    outs[[j]] <- gene_res_mat
    if (length(outs[[j]]) == 0) {
      print("P-values cannot be returned because no test statistics were calculated.")
    } else {
      colnames(outs[[j]]) <- unlist(input$testStats)
      rownames(outs[[j]]) <- c(
        "Tail area probability",
        "Empirical test statistic",
        "Standard deviations from simulated distribution"
      )
      write.csv(t(outs[[j]]), file = "output.pvals.PhyloMAd.csv")
    }
    resvector <- matrix(as.vector(t(outs[[j]])), nrow = 1)
    rownames(resvector) <- as.character(input$dataPath[j, 1])
    colnames(resvector) <- c(
      paste0(colnames(outs[[j]]), ".tail.area.p"),
      paste0(colnames(outs[[j]]), ".empirical.statistic"),
      paste0(colnames(outs[[j]]), ".stdev.from.pred.dist")
    )
    outs[[j]] <- resvector
  }

  if ("phyloempres" %in% unlist(input$whatToOutput)) {
    if ("aindex" %in% unlist(input$testStats)) {
      branchwise_assessment <- rbind(gene_results$aindex.allp, gene_results$aindex.sds)
      rownames(branchwise_assessment) <- c("Branch-wise P-values", "Branch-wise SDPD")
      write.csv(branchwise_assessment, file = "Branch-wise.assessment.results.csv")
    }
    write.tree(gene_results$empirical.tree, file = "estimate.empirical.data.clock.free.tre")
  }

  if ("simdat" %in% unlist(input$whatToOutput)) {
    if (input$outputFormat == "phylip") {
      writer <- write.dna
      ext <- ".phy"
      format <- "interleaved"
    } else if (input$outputFormat == "fasta") {
      writer <- write.dna
      ext <- ".fasta"
      format <- "fasta"
    } else if (input$outputFormat == "nexus") {
      writer <- write.nexus.data
      ext <- ".nex"
      format <- "dna"
    }

    for (i in 1:input$Nsims) {
      writer(
        gene_results$simDat[[i]]$alignment,
        file = paste0("predictive.data.", i, ext), format = format
      )
    }
  }

  if ("phylosimres" %in% unlist(input$whatToOutput)) {
    write.tree(gene_results$simPhylos, file = "estimate.predictive.data.tre")
  }

  if ("testPlots" %in% unlist(input$whatToOutput)) {
    empstats <- unlist(gene_results[grep("emp[.]", names(gene_results))])
    simstats <- do.call(cbind, gene_results[grep("sim[.]", names(gene_results))])
    statsmat <- rbind(empstats, simstats)
    statlabels <- vector()
    if ("imbal" %in% selected_stats) statlabels <- c(statlabels, "Imbalance (Colless index)")
    if ("stemmystat" %in% selected_stats) statlabels <- c(statlabels, "Stemminess")
    if ("df" %in% selected_stats) statlabels <- c(statlabels, "Df neutrality statistic")
    if ("aindex" %in% selected_stats) statlabels <- c(statlabels, "A index")
    if ("trlen" %in% selected_stats) statlabels <- c(statlabels, "Tree length")
    if ("maha" %in% selected_stats) statlabels <- c(statlabels, "Squared Mahalanobis distance")

    if (length(empstats) == 0) {
      print("Test plots cannot be returned because no test statistics were calculated.")
    } else {
      pdf("adequacy.tests.plots.pdf", useDingbats = FALSE, height = 5)
      for (i in seq_along(length(empstats))) {
        sdstat <- sd(statsmat[2:nrow(statsmat), i])
        hist(
          statsmat[2:nrow(statsmat), i],
          xlim = c(
            min(statsmat[, i]) - sdstat,
            max(statsmat[, i]) + sdstat
          ),
          xlab = statlabels[i],
          ylab = "Frequency of predictive simulations",
          main = ""
        )
        abline(v = empstats[i], col = "red", lwd = 3)
      }
      dev.off()

      if ("aindex" %in% unlist(input$testStats)) {
        pdf(
          "adequacy.branch-wise.tests.plots.pdf",
          useDingbats = FALSE,
          height = 5 + (Ntip(gene_results$empirical.tree) * 0.2)
        )

        treetoprint <- gene_results$empirical.tree
        treetoprint$edge.length <- NULL
        palette <- colorRampPalette(c("blue", "green", "yellow", "red"))
        sdcol <- palette(5)[as.numeric(cut(gene_results$aindex.sds, breaks = 5))]

        plot(
          treetoprint,
          col = palette(5)[as.numeric(cut(gene_results$aindex.allp, breaks = 5))],
          main = "Branches labelled by P-value"
        )
        edgelabels(round(gene_results$aindex.allp, 2))

        plot(treetoprint, col = sdcol, main = "Branches labelled by SDPD")
        edgelabels(round(gene_results$aindex.sds, 2))
        dev.off()
      }
    }
  }

  setwd("..")

  if ("simple" %in% unlist(input$whatToOutput)) {
    system(paste0("rm -r ", as.character(input$dataPath[j, 1]), ".phylomad.clock"))
  }
}

if (nrow(input$dataPath) > 1 && "pvals" %in% unlist(input$whatToOutput) || "simple" %in% unlist(input$whatToOutput)) {
  all_output <- do.call("rbind", outs)
  write.csv(t(all_output), file = "output.all.loci.clock.PhyloMAd.csv")
}

setwd(initial_dir)

print(
  paste0(
    "Assessment completed successfully in ",
    round(proc.time()[3] - process_time, 3),
    " seconds."
  )
)

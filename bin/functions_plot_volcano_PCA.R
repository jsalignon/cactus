
# note: these functions function were adapted from sleuth::plot_volcano, sleuth::plot_loadings, sleuth::plot_pc_variance and sleuth::plot_pca (using the script when all parameters are left to their default values)

theme_set(theme_bw() + theme(plot.title = element_text(hjust = 0.5, size = 11)))


##### volcano plots

plot_volcano_custom <- function(res, sig_level, label_column, title, sig_color = 'red', point_alpha = 0.4, label_alpha = 0.9, label_size = 2, label_nb = 15, epsilon = 10^-320){
  res = res[!is.na(res$padj), ]
  res = dplyr::mutate(res, significant = res$padj < sig_level)
  res$min_log10_FDR = -log10(res$padj + epsilon)
  p1 = ggplot(res, aes(L2FC, min_log10_FDR)) + geom_point(aes(colour = significant), alpha = point_alpha) + scale_colour_manual(values = c('black', sig_color))
  # xlab = paste0('log2(', cond1, ') - log(', cond2, ')')
  xlab = paste0('Log2 Fold Change')
  p1 = p1 + xlab(xlab) + ylab('-log10(FDR)') + geom_vline(xintercept = 0, colour = 'black', linetype = 'longdash')
  
  # adding the highlight for the top N up and down
  top_N = seq_len(label_nb)
  top_reg_up = order(res$min_log10_FDR * -sign(res$L2FC), decreasing = F)[top_N]
  top_reg_down = order(res$min_log10_FDR * sign(res$L2FC), decreasing = F)[top_N]
  top_reg = c(top_reg_up, top_reg_down)
  res1 = res[top_reg, ]
  
  nudge_y = max(res$min_log10_FDR) / 60
  p1 = p1 + ggtitle(title) + geom_text(data = res1, aes_string(label = label_column), nudge_y = nudge_y, alpha = label_alpha, size = label_size)
  
  p1
}

# res = res_volcano ; sig_level = FDR_threshold ; label_column = 'gene_name' ; title = paste(COMP, 'mRNA') ; sig_color = 'red' ; point_alpha = 0.4 ; label_alpha = 0.9 ; label_size = 2 ; label_nb = 15 ; epsilon = 10^-320



##### PCA plots


plot_pc_variance_custom <- function(prcomp1, title, max_pc = 5){
  eigenvalues <- (prcomp1$sdev)^2
  var_explained <- eigenvalues * 100/sum(eigenvalues)
  vec = 1:max_pc
  pc_df <- data.frame(PC_count = vec, var = var_explained[vec])
  
  p1 <- ggplot(pc_df, aes(x = PC_count, y = var)) + geom_bar(stat = 'identity')  + scale_x_continuous(breaks = 1:length(eigenvalues)) + ylab('% of variance') + xlab('principal components') + ylim(0, 100) + ggtitle(title) # + ggtitle(paste('PC Variance:', id_comp))
  return(p1)
}


plot_loadings_custom <- function(prcomp1, pc_input = 1, pc_count = 20){
  
  loadings <- prcomp1$x[, pc_input]
  
  loadings_sorted <- sort(abs(loadings), decreasing = TRUE)
  
  dat <- data.frame(entry = names(loadings_sorted), loadings = loadings_sorted)[1:pc_count, ]
  dat$sign = factor(sign(loadings)[dat$entry])
  dat$entry <- factor(dat$entry, levels = unique(dat$entry))
  
  p1 <- ggplot(dat, aes(x = entry, y = loadings, fill = sign)) + geom_bar(stat = 'identity') + ylab(paste0('contribution scores PC', pc_input)) + xlab('entries') + ggtitle('') # + ggtitle(paste0('PC ', pc_input))
  p1 <- p1 + theme(axis.text.x = element_text(hjust = 1, angle = 90, vjust = 0.5)) + ggtitle('') # + ggtitle(paste0('Loadings PC', pc_input, ': ', id_comp))
  
  return(p1)
  
}


plot_pca_custom <- function(prcomp1, pc_x, pc_y){
  pc_x = paste0('PC', pc_x)
  pc_y = paste0('PC', pc_y)
  pcs <- as.data.frame(prcomp1$rotation[, c(pc_x, pc_y)], stringsAsFactors = F)
  pcs$sample = rownames(pcs)
  pcs$condition = sapply(pcs$sample, function(x) strsplit(x, '_')[[1]][1])
  p1 <- ggplot(pcs, aes_string(x = pc_x, y = pc_y, colour = 'condition')) + geom_text(aes(label = sample)) + ggtitle('') # + ggtitle(paste('PCA:', id_comp))
  return(p1)
}


get_lp <- function(prcomp1, pc_x, pc_y, title){
  p1 = plot_pc_variance_custom(prcomp1, title)
  p2 = plot_pca_custom(prcomp1, pc_x, pc_y)
  # p2 = dba.plotPCA(dbo, components = c(pc_x, pc_y))
  p3 = plot_loadings_custom(prcomp1, pc_input = pc_x, pc_count = 20)
  p4 = plot_loadings_custom(prcomp1, pc_input = pc_y, pc_count = 20)
  lp_3_4 = list(p1, p2, p3, p4)
}


make_4_plots <- function(lp){
  pushViewport(viewport(layout = grid.layout(2, 2)))
  vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
  print(lp[[1]], vp = vplayout(1, 1))
  print(lp[[2]], vp = vplayout(1, 2))
  print(lp[[3]], vp = vplayout(2, 1))
  print(lp[[4]], vp = vplayout(2, 2))
}






DiffBind__pv_pcmask__custom <- 
function (pv, numSites, mask, sites, removeComps, cor = FALSE,
    bLog = TRUE)
{
    if (missing(numSites))
        numSites <- nrow(pv$binding)
    if (is.null(numSites))
        numSites <- nrow(pv$binding)
    numSites <- min(numSites, nrow(pv$binding))
    if (missing(sites))
        sites <- 1:numSites
    if (is.null(sites))
        sites <- 1:numSites
    if (missing(mask))
        mask <- rep(TRUE, ncol(pv$class))
    for (i in which(mask)) {
        if (nrow(pv$peaks[[i]]) == 0) {
            mask[i] <- FALSE
        }
    }
    if (sum(mask) < 2) {
        stop("Need at least two samples for PCA.", call. = FALSE)
    }
    res <- NULL
    res$class <- pv$class
    pv$values <- pv$binding[sites, c(FALSE, FALSE, FALSE, mask)]
    active <- apply(pv$values, 1, function(x) DiffBind:::pv.activefun(x))
    numSites <- min(numSites, sum(active))
    pv$values <- pv$values[active, ][1:numSites, ]
    if (!missing(removeComps)) {
        pv$values <- DiffBind:::pv.removeComp(pv$values, numRemove = removeComps)
    }
    if (bLog) {
        if (max(pv$values) > 1) {
            pv$values[pv$values <= 1] <- 1
            pv$values <- log2(pv$values)
        }
    }
    if (nrow(pv$values) >= sum(mask)) {
        res$pc <- prcomp((pv$values))
    }
    res$mask <- mask
    return(res)
}


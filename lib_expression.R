 
##################################################
####
#### functions related to RNAseq analysis
####
##################################################

# == title
# read count from htseq-count 
read.count = function(sample_id, gencode, template) {
	
	gene_id = gencode$getGeneID()
	count = matrix(nrow = length(gene_id), ncol = length(sample_id))
	rownames(count) = gene_id
	colnames(count) = sample_id
	
	for(sid in sample_id) {
		qqcat("[@{sid}] reading count data.\n")
		d = read.table(qq(template), sep = "\t", row.names = 1)
		nr = nrow(d)
		d = d[0:4 - nr, , drop = FALSE]
		count[rownames(d), sid] = d[[1]]
	}
	return(count)
}

# useful links:
# http://bib.oxfordjournals.org/content/14/6/671.full.pdf+html
# https://www.biostars.org/p/56919/
# http://haroldpimentel.wordpress.com/2014/05/08/what-the-fpkm-a-review-rna-seq-expression-units/
normalize.count = function(count, gencode, method = "rpkm", param) {
	method = tolower(method)[1]
	
	gene_length = gencode$geneLength(type = "exon")
	gene_length = gene_length[rownames(count)]

	qqcat("normalizing raw count by @{method}\n")
	
	if(method == "rpkm") {
		all_count = colSums(count)
		
		rpkm = matrix(0, nrow = nrow(count), ncol = ncol(count))
		rownames(rpkm) = rownames(count)
		colnames(rpkm) = colnames(count)
		for(i in seq_len(nrow(count))) {
			rpkm[i, ] = count[i, ] / gene_length[i]
		}
			
		for(j in seq_len(ncol(count))) {
			rpkm[, j] = rpkm[, j] / all_count[j]
		}
		
		rpkm = 10^9 * rpkm
		return(rpkm)
	} else if(method == "voom") {
		require(limma)
		expr = voom(count, normalize.method = param$normalize.method)$E
		if(param$normalizeTOGeneLength) {
			for(i in seq_len(nrow(count))) {
				expr[i, ] = expr[i, ] / gene_length[i] * 1e3
			}
		}
		return(expr)
	} else if(method == "tpm") {
		gene_length = gencode$geneLength(type = "exon")
		gene_length = gene_length[rownames(count)]
	
		tpm = matrix(0, nrow = nrow(count), ncol = ncol(count))
		rownames(tpm) = rownames(count)
		colnames(tpm) = colnames(count)
		for(i in seq_len(nrow(count))) {
			tpm[i, ] = count[i, ] / gene_length[i]
		}
		
		sum_ratio = colSums(tpm)
		for(j in seq_len(ncol(count))) {
			tpm[, j] = tpm[, j] / sum_ratio[j]
		}
		if(param$normalizeTOGeneLength) {
			for(i in seq_len(nrow(count))) {
				tpm[i, ] = tpm[i, ] / gene_length[i] * 1e3
			}
		}
		tpm = tpm * 1e6
		return(tpm)
	} else if(method == "tc") {
		all_count = colSums(count)
		mean_all_count = mean(all_count)
		expr = matrix(0, nrow = nrow(count), ncol = ncol(count))
		rownames(expr) = rownames(count)
		colnames(expr) = colnames(count)
		for(j in seq_len(ncol(count))) {
			expr[, j] = count[, j] / all_count[j] * mean_all_count
		}
		if(param$normalizeTOGeneLength) {
			for(i in seq_len(nrow(count))) {
				expr[i, ] = expr[i, ] / gene_length[i] * 1e3
			}
		}
		return(expr)
	} else if(method == "med") {
		median_count = apply(count, 2, function(x) median(x[x > 0]))
		mean_median_count = median(median_count)
		expr = matrix(0, nrow = nrow(count), ncol = ncol(count))
		rownames(expr) = rownames(count)
		colnames(expr) = colnames(count)
		for(j in seq_len(ncol(count))) {
			expr[, j] = count[, j] / median_count[j] * mean_median_count
		}
		if(param$normalizeTOGeneLength) {
			for(i in seq_len(nrow(count))) {
				expr[i, ] = expr[i, ] / gene_length[i] * 1e3
			}
		}
		return(expr)
	} else if(method == "deseq2") {
		require("DESeq2")
		df = data.frame(sth = rep("foo", ncol(count)))
		rownames(df) = colnames(count)
		cds = DESeqDataSetFromMatrix(countData = count, colData = df, design = ~ 1)
		cds = estimateSizeFactors(cds)
		cts = counts(cds, normalized=TRUE)
		if(param$varianceStabilize) {
			cds = estimateDispersions(cds)
			vsd = getVarianceStabilizedData(cds)
			expr = as.data.frame(vsd)
		} else {
			expr = as.data.frame(cts)
		}
		if(param$normalizeTOGeneLength) {
			for(i in seq_len(nrow(count))) {
				expr[i, ] = expr[i, ] / gene_length[i] * 1e3
			}
		}
		return(expr)
	} else if(method == "tmm") {
		require(edgeR)
		sizes = calcNormFactors(DGEList(counts = count))
		expr = matrix(0, nrow = nrow(count), ncol = ncol(count))
		rownames(expr) = rownames(count)
		colnames(expr) = colnames(count)
		for(j in seq_len(ncol(count))) {
			expr[, j] = count[, j] / sizes$samples$norm.factors[j]
		}
		if(param$normalizeTOGeneLength) {
			for(i in seq_len(nrow(count))) {
				expr[i, ] = expr[i, ] / gene_length[i] * 1e3
			}
		}
		return(expr)
	} else {
		stop(qq("@{method} is not supported.\n"))
	}
}

read_rnaseqqc = function(sample_id, template) {
	matqc = matrix(nrow = length(sample_id), ncol = 23)
	colnames(matqc) = c("Total", "Unique", "Duplicates", "Duplication Rate", "Estimated Library Size",
		"Mapped", "Mapping Rate", "Mapped Unique", "Mapped Unique Rate", "rRNA", "rRNA rate",
		"Intragenic Rate", "Exonic Rate", "Intronic Rate", "Intergenic Rate", "Expression Profiling Efficiency", "Expressed Transcripts",
		"End 1 Sense", "End 1 Antisense", "End 2 Sense", "End 2 Antisense", "End 1 % Sense", "End 2 % Sense")
	rownames(matqc) = sample_id
	for(sid in sample_id) {
		lines = readLines(qq(template))
		txt = strsplit(lines, "\t")
		matqc[sid, txt[[1]]] = as.numeric(txt[[2]])
		matqc[sid, txt[[3]]] = as.numeric(txt[[4]])
		matqc[sid, txt[[5]]] = as.numeric(txt[[6]])
		matqc[sid, txt[[7]]] = as.numeric(txt[[8]])
	}
	return(matqc)
}

# == title
# RNASEQ QC plot
#
# == param
# -count count matrix
# -expr normalized expressed matrix
# -main title
# -is.logged is the normalized expression already log-transformed
#
plot_count_distribution = function(count, expr, sample, main = NULL, is.logged = FALSE) { 
	
	sample_id = colnames(count)

	count = log2(count + 1)
	q90 = quantile(as.matrix(count), 0.9)
	count = apply(count, 2, function(x) {
		x[x > q90] = q90
		x
	})

	layout(cbind(1:2, 3:4))

	par(mar = c(12, 4, 4 ,2), xpd = NA)
	density_heatplot(count, draw.quantiles = TRUE)
	aa = seq(0, 50, by = 2)
	aa = aa[aa <= max(count)]
	axis(side = 2, at = aa, labels = 2^aa-1)
	par(las = 3)
	axis(side = 1, at = seq_along(sample_id), labels = sample_id)
	mtext("log2(count + 1)", side = 2, line = 2) 
	title(qq("density heatmap for count (< q90) distribution\n@{main}"))

	par(xpd = FALSE)
	den_x = matrix(nrow = 512, ncol = dim(count)[2])
	den_y = matrix(nrow = 512, ncol = dim(count)[2])
	for(i in seq_len(dim(count)[2])) {
		den  = density(count[, i])
		den_x[, i] = den$x
		den_y[, i] = den$y
	}
	matplot(den_x, den_y, type = "l", xlab = "log2(count+1)", ylab = "density",main = "density distribution of counts")

	sample_id = colnames(expr)

	if(!is.logged) {
		expr = log2(expr + 1)
	}
	q90 = quantile(as.matrix(expr), 0.9)
	expr = apply(expr, 2, function(x) {
		x[x > q90] = q90
		x
	})

	par(mar = c(12, 4, 4 ,2), xpd = NA)
	density_heatplot(expr, draw.quantiles = TRUE)
	aa = seq(0, 50, by = 2)
	aa = aa[aa <= max(expr)]
	axis(side = 2, at = aa, labels = 2^aa-1)
	par(las = 3)
	axis(side = 1, at = seq_along(sample_id), labels = sample_id)
	mtext("log2(expr + 1)", side = 2, line = 2) 
	title(qq("density heatmap for normalized expression (< q90)\n@{main}"))

	par(xpd = FALSE)
	den_x = matrix(nrow = 512, ncol = dim(expr)[2])
	den_y = matrix(nrow = 512, ncol = dim(expr)[2])
	for(i in seq_len(dim(expr)[2])) {
		den = density(expr[, i])
		den_x[, i] = den$x
		den_y[, i] = den$y
	}
	matplot(den_x, den_y, type = "l", xlab = "log2-expression", ylab = "density",main = "density distribution of normalized expression")

	layout(matrix(1))
}



# == title
# use colors to represent density
#
# == param
# -x              a matrix or a list. If it is a matrix, density will be calculated by columns
# -col            colors theme, can be a vector or a list
# -draw.quantiles whether to draw quantile
# -align          align at the two ends?
# -...            pass to `plot.default`
#
density_heatplot = function(x, col=rev(brewer.pal(10, "Spectral")), draw.quantiles = FALSE, align = TRUE, ...) {
    if(is.vector(x) && class(x) != "list") {
        x = as.matrix(x)
    }
    if(is.matrix(x)) {
        n = dim(x)[2]
        # if different styles of colors are used, it should be formatted as a list
        # because the number of sections of colors may be different

        dx = apply(x, 2, function(x) density(x)$x)  # data value
        dy = apply(x, 2, function(x) density(x)$y)  # density value
        quantile.values = apply(x, 2, quantile)
    }
    if(is.list(x)) {
        n = length(x)

        dx = sapply(x, function(x) density(x)$x)  # data value
        dy = sapply(x, function(x) density(x)$y)  # density value
        quantile.values = sapply(x, quantile)
    }

    if(!is.list(col)) {
        col = rep(list(col), n)
    }
    if(is.list(col) && length(col) != n) {
        stop("Since 'col' is specified as a list, it should has the same length as numbers of columns in 'x'.")
    }
    if(!all(sapply(col, length) > 1)) {
        stop("Length of any color vector should contain at least two colors.")
    }

    min.density = min(as.vector(dy))
    max.density = max(as.vector(dy))
    range.density = max.density - min.density
    dy = (dy - min.density) / range.density

    min.value = min(as.vector(dx))
    max.value = max(as.vector(dx))
    range.value = max.value - min.value

    plot(c(0, n+1), c(min.value, max.value), type = "n", axes=FALSE, ann=FALSE, ...)
    for(j in 1:n) {

    	col_fun = colorRamp2(breaks = seq(0, 1, length.out=length(col[[j]])), col = col[[j]])

        for(i in 2:length(dy[, j])) {
            rect(j-0.5, dx[i-1, j], j+0.5, dx[i, j], col=col_fun(dy[i, j]), border=col_fun(dy[i, j]))
        }

        if(align) {
            rect(j-0.5, min(dx[, j]), j+0.5, min.value, col=col_fun(min.density), border=col_fun(min.density))
            rect(j-0.5, max(dx[, j]), j+0.5, max.value, col=col_fun(min.density), border=col_fun(min.density))
        }
    }
    #axis(side = 2)
    if(draw.quantiles) {
        for(i in 1:dim(quantile.values)[1]) {
            lines(1:n, quantile.values[i, ], col="black", lwd=1)
        }
        text(rep(n+1, dim(quantile.values)[1]), quantile.values[, n], rownames(quantile.values), cex=0.8, adj=c(0, 0.5))
    }
}


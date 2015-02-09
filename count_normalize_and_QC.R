library(GetoptLong)

normalize = "rpkm"
param = NULL

GetoptLong.options("startingMsg" = "Normalize raw counts and make some QC plots\n")
GetoptLong.options("endingMsg" = qq("
parameters for different normalization methods:

voom    normalization.method = none|scale|quantile|cyclicloess
DESeq2   varianceStabilize = 1|0
for all  normalizeTOGeneLength = 1|0

Examples:

Rscript @{get_scriptname()} \
        --input /icgc/dkfzlsdf/analysis/hipo/hipo_016/rnaseq_low_grade_glioma \
        --input_filter 'AK\\\\d+$' \
        --output /icgc/dkfzlsdf/analysis/hipo/hipo_016/rnaseq_low_grade_glioma/hipo_016_ \
        --normalize rpkm voom deseq2
"))

# input = "/icgc/dkfzlsdf/analysis/hipo_035/data_types/RNAseq/year6_gencode19"
# output = "~/lina_year6_"

input_filter = NULL
GetoptLong(c("input=s", "input dir which contains folders for samples",
             "input_filter=s", "a regular expression to filter sample id (which is the dirname under `input`",
	         "output=s", "output prefix, should be a absolute path prefix",
	         "normalize=s{1,}", "normalization method. It can be more than one values. Options are rpkm, voom, tpm, tc, med, deseq2, tmm.",
	         "param=s%", "additional arguments for `normalize` method"))

default_param = list(
	"normalization.method" = "none",
	"varianceStabilize" = 1,
	"normalizeTOGeneLength" = 0
)
if(!is.null(param)) {
	default_param[names(param)] = param
	default_param$varianceStabilize = as.numeric(default_param$varianceStabilize)
	default_param$normalizeTOGeneLength = as.numeric(default_param$normalizeTOGeneLength)
}
param = default_param
 
setwd(dirname(get_scriptname()))


source("lib_expression.R")


library(pheatmap)
library(RColorBrewer)
library(circlize)
library(GTF)


qq.options(cat_prefix = function(x) xterm256::style(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), fg = "green") )

# configuration for RNASeq 
counting_file_template       = qq("`input`/@{sid}/@{sid}.mkdup.exon.count", code.patter = "`CODE`")
nodup_counting_file_template = qq("`input`/@{sid}/@{sid}.mkdup.exon.nodup.count", code.patter = "`CODE`")

print(input_filter)

sample_id = dir(input)
if(!is.null(input_filter)) {
	sample_id = grep(input_filter, sample_id, value = TRUE)
}

n_sample = length(sample_id)
if(n_sample == 0) {
	stop("There is no sample under your settings.\n")
}

qqcat("There @{ifelse(n_sample == 1, 'is', 'are')} @{length(sample_id)} sample@{ifelse(n_sample == 1, '', 's')}.\n")

GENE_NUMBER = c("gencode_v17" = 57281,
				"gencode_v17_lincRNA" = 63802,
				"gencode_v21_lincRNA" = 68150,
	            "gencode_v18" = 57445,
	            "gencode_v19" = 57820,
	            "gencode_v20" = 58688,
	            "gencode_vM1" = 37310,
	            "gencode_vM2" = 38924,
	            "gencode_vM3" = 41128,
	            "refGene" = 25763,
	            "knownGene" = 65723
	            )

sid = sample_id[1]
d = read.table(qq(counting_file_template), sep = "\t", row.names = 1)
ng = nrow(d) - 5

if(sum(GENE_NUMBER == ng)) {
	annotation = names(GENE_NUMBER)[which(GENE_NUMBER == ng)]
	qqcat("I guess your annotation is @{annotation}.\n")
} else {
	stop("Cannot find propor gene annotation file.\n")
}



qqcat("loading /icgc/dkfzlsdf/analysis/B080/guz/gencode/@{annotation}_transcript_merged.RData\n")
obj = load(qq("/icgc/dkfzlsdf/analysis/B080/guz/gencode/@{annotation}_transcript_merged.RData"))


gencode = get(obj, envir = .GlobalEnv)

expression = vector("list", length = length(normalize) + 1)
names(expression) = c("count", normalize)
expression_nodup = vector("list", length = length(normalize) + 1)
names(expression_nodup) = c("count", normalize)

count1 = read.count(sample_id, gencode, template = counting_file_template)
expression[["count"]] = count1
count2 = read.count(sample_id, gencode, template = nodup_counting_file_template)
expression_nodup[["count"]] = count2

expression$normalize = normalize
expression$param = param
expression_nodup$normalize = normalize
expression_nodup$param = param


for(nmlz in normalize) {
	
	pdf(qq("@{output}count_qc_normalize_by_@{nmlz}.pdf"), width = 16, height = 16)

	#############################################
	## count and expr distribution 
	#############################################
	qqcat("count with keeping duplicates.\n")

	expr1 = normalize.count(count1, gencode, method = nmlz, param)
	expression[[nmlz]] = expr1

	if(nmlz %in%  c("voom")) {
		plot_count_distribution(count1, expr1, main = qq("keep duplicate, normalize = @{nmlz}"), is.logged = TRUE)	
	} else {
		plot_count_distribution(count1, expr1, main = qq("keep duplicate, normalize = @{nmlz}"), is.logged = FALSE)
	}

	# three heatmaps for 10% lowest, median, highest expressed genes
	qqcat("generating heatmaps for 10% lowest, median, highest expressed genes\n")
	l = apply(count1, 1, function(x) any(x > 0)) & apply(expr1, 1, sd) > 0
	mat = expr1[l, , drop = FALSE]
	avg = apply(mat, 1, mean)
	l = avg < quantile(avg, 0.1)
	if(sum(l) > 1) {
		pheatmap(mat[l, , drop = FALSE], main = qq("10% low expressed genes, @{nmlz}"), show_rownames = FALSE)
	}
	if(nmlz %in%  c("voom")) {
		l = avg < quantile(avg, 0.55) & avg > quantile(avg, 0.45)
		pheatmap(mat[l, , drop = FALSE], main = qq("10% middle expressed genes, @{nmlz}"), show_rownames = FALSE)
		l = avg > quantile(avg, 0.9)
		pheatmap(mat[l, , drop = FALSE], main = qq("10% high expressed genes, @{nmlz}"), show_rownames = FALSE)
	} else {
		l = avg < quantile(avg, 0.55) & avg > quantile(avg, 0.45)
		pheatmap(log2(mat[l, , drop = FALSE] + 1), main = qq("10% middle expressed genes, log2(@{nmlz} + 1)"), show_rownames = FALSE)
		l = avg > quantile(avg, 0.9)
		pheatmap(log2(mat[l, , drop = FALSE] + 1), main = qq("10% high expressed genes, log2(@{nmlz} + 1)"), show_rownames = FALSE)
	}


	#############################################
	## count and expr distribution, remove duplicates 
	#############################################
	qqcat("count with removing duplicates.\n")

	expr2 = normalize.count(count2, gencode, method = nmlz, param)
	expression_nodup[[nmlz]] = expr2

	if(nmlz %in%  c("voom")) {
		plot_count_distribution(count2, expr2, main = qq("remove duplicate, normalize = @{nmlz}"), is.logged = TRUE)	
	} else {
		plot_count_distribution(count2, expr2, main = qq("remove duplicate, normalize = @{nmlz}"), is.logged = FALSE)
	}

	# three heatmaps for 10% lowest, median, highest expressed genes
	qqcat("generating heatmaps for 10% lowest, median, highest expressed genes\n")
	l = apply(count2, 1, function(x) any(x > 0)) & apply(expr2, 1, sd) > 0
	mat = expr2[l, , drop = FALSE]
	avg = apply(mat, 1, mean)
	l = avg < quantile(avg, 0.1)
	if(sum(l)) {
		pheatmap(mat[l, , drop = FALSE], main = qq("10% low expressed genes, @{nmlz}"), show_rownames = FALSE)
	}
	if(nmlz %in%  c("voom")) {
		l = avg < quantile(avg, 0.55) & avg > quantile(avg, 0.45)
		pheatmap(mat[l, , drop = FALSE], main = qq("10% middle expressed genes, @{nmlz}"), show_rownames = FALSE)
		l = avg > quantile(avg, 0.9)
		pheatmap(mat[l, , drop = FALSE], main = qq("10% high expressed genes, @{nmlz}"), show_rownames = FALSE)
	} else {
		l = avg < quantile(avg, 0.55) & avg > quantile(avg, 0.45)
		pheatmap(log2(mat[l, , drop = FALSE] + 1), main = qq("10% middle expressed genes, log2(@{nmlz} + 1)"), show_rownames = FALSE)
		l = avg > quantile(avg, 0.9)
		pheatmap(log2(mat[l, , drop = FALSE] + 1), main = qq("10% high expressed genes, log2(@{nmlz} + 1)"), show_rownames = FALSE)
	}


	################################################
	### comparison between remove duplicates and not
	################################################
	count2 = count2[rownames(count1), , drop = FALSE]
	expr2 = expr2[rownames(expr1), , drop = FALSE]

	gi = gencode$getGeneID()
	name = gencode$getValueByGeneID(gi, type = "name")
	col = "#00000080"
	pch = rep(".", length(gi))
	pch[grep("^mt(-|_)", name, ignore.case = TRUE)] = "M"
	par(mfrow = c(3, 3))
	qqcat("plot scatter plot for all count\n")
	for(i in seq_len(ncol(count1))) {
		if(i > 9) {
			break
		}
		xmax = max(c(log2(count1[, i] + 1), log2(count2[, i] + 1)))
		plot(NULL, xlim = c(0, xmax), ylim = c(0, xmax), xlab = "keep duplicate", ylab = "remove duplicate", main = qq("log2(count+1): @{colnames(count1)[i]}"))
		abline(a = 0, b = 1, col = "grey")
		points(log2(count1[, i] + 1), log2(count2[, i] + 1), pch = pch, cex = 0.7, col = col)
	}
	par(mfrow = c(1, 1))

	qqcat("plot scatter plot for all normalized expression\n")
	par(mfrow = c(3, 3))
	for(i in seq_len(ncol(count1))) {
		if(i > 9) {
			break
		}
		if(nmlz %in%  c("voom")) {
			xmax = max(c(expr1[, i], expr2[, i]))
			xmin = min(c(expr1[, i], expr2[, i]))
			plot(NULL, xlim = c(xmin, xmax), ylim = c(xmin, xmax), xlab = "keep duplicate", ylab = "remove duplicate", main = qq("@{nmlz}: @{colnames(count1)[i]}"))
			abline(a = 0, b = 1, col = "grey")
			points((expr1[, i] + 1), (expr2[, i] + 1), pch = pch, cex = 0.7, col = col)
		} else {
			xmax = max(c(log2(expr1[, i] + 1), log2(expr2[, i] + 1)))
			xmin = min(c(log2(expr1[, i] + 1), log2(expr2[, i] + 1)))
			plot(NULL, xlim = c(xmin, xmax), ylim = c(xmin, xmax), xlab = "keep duplicate", ylab = "remove duplicate", main = qq("log2(@{nmlz}+1): @{colnames(count1)[i]}"))
			abline(a = 0, b = 1, col = "grey")
			points(log2(expr1[, i] + 1), log2(expr2[, i] + 1), pch = pch, cex = 0.7, col = col)
		}
	}
	par(mfrow = c(1, 1))
	dev.off()
}


save(expression, file = qq("@{output}expression.RData"))
save(expression_nodup, file = qq("@{output}expression_nodup.RData"))


matqc = read_rnaseqqc(sample_id, qq("`input`/@{sid}/rnaseqqc/@{sid}/@{sid}.metrics.txt", code.patter = "`CODE`"))
write.table(matqc, file = qq("@{output}rnaseqqc.txt"), sep = "\t", quote = FALSE)

pdf(qq("@{output}rnaseqqc.pdf"), width = length(sample_id) * 1, height = 8)

par(mar = c(8, 4, 4, 4))
max_number = max(matqc[, 1])
plot(NULL, xlim = c(0.5, length(sample_id) + 0.5), ylim = c(0, max_number*1.2), axes = FALSE, ann = FALSE)
for(i in seq_along(sample_id)) {
	rect(i - 0.4, 0, i, matqc[i, 2], col = "blue", border = NA)
	rect(i - 0.4, matqc[i, 2], i, matqc[i, 2] + matqc[i, 3], col = "green", border = NA)
	rect(i, 0, i + 0.4, matqc[i, 8], col = "red", border = NA)
	rect(i, matqc[i, 8], i + 0.4, matqc[i, 6], col = "pink", border = NA)
}
par(las = 3)
axis(side = 1, at = seq_along(sample_id), labels = sample_id)
axis(side = 2)
title(ylab = "Number of reads", main = "Reads statistics")
legend("topleft", pch = 15, col = c("blue", "green", "red", "pink"), legend = c("Unique", "Duplicates", "Mapped", "Mapped unique"))
box()

par(new = TRUE)
plot(seq_along(sample_id) - 0.2, matqc[, 4], xlim = c(0.5, length(sample_id) + 0.5), ylim = c(0, 1.2), type = "b", pch = "D", axes = FALSE, ann = FALSE)
lines(seq_along(sample_id) + 0.2, matqc[, 7], type = "b", pch = "M") 
axis(side = 4, at = seq(0, 1, by = 0.2), labels = seq(0, 1, by = 0.2))
mtext("Percent", side = 4, line = 2)
legend("topright", pch = c("D", "M"), legend = c("Duplication rate", "Mapped unique rate"))

par(new = FALSE)

barplot(t(matqc[, 13:15]), col = c("blue", "green", "grey"), ylim = c(0, 1.2), main = "exonic/intronic/intergenic rate")
legend("topright", col = c("blue", "green", "grey"), pch = 15, legend = c("Exonic", "Intronic", "Intergenic"))
barplot(t(matqc[, 22:23]/100), col = c("blue", "green"), ylim = c(0, 1.2), main = "strand specificity")
legend("topright", col = c("blue", "green"), pch = 15, legend = c("End 1 % Sense", "End 2 % Sense"))

dev.off()


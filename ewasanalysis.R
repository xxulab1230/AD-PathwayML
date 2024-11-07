library(readxl)
setwd("#####/project_directory/EWAS")  # Set working directory to project folder
BiocManager::install("GEM")
library(GEM)
GEM_GUI()

DATADIR <- system.file('extdata', package='GEM')
RESULTDIR <- getwd()

env_file_name <- paste(DATADIR, "env.txt", sep = .Platform$file.sep)
covariates_file_name <- paste(DATADIR, "cov.txt", sep = .Platform$file.sep)
covariates_file_name_gxe <- paste(DATADIR, "gxe.txt", sep = .Platform$file.sep)
methylation_file_name <- paste(DATADIR, "methylation.txt", sep = .Platform$file.sep)
snp_file_name <- paste(DATADIR, "snp.txt", sep = .Platform$file.sep)

Emodel_pv <- 1
Gmodel_pv <- 1e-04
GxEmodel_pv <- 1

Emodel_result_file_name <- paste(RESULTDIR, "Result_Emodel.txt", sep = .Platform$file.sep)
Gmodel_result_file_name <- paste(RESULTDIR, "Result_Gmodel.txt", sep = .Platform$file.sep)
GxEmodel_result_file_name <- paste(RESULTDIR, "Result_GxEmodel.txt", sep = .Platform$file.sep)

Emodel_qqplot_file_name <- paste(RESULTDIR, "QQplot_Emodel.jpg", sep = .Platform$file.sep)

GEM_Emodel(env_file_name, covariates_file_name, methylation_file_name, Emodel_pv, Emodel_result_file_name, Emodel_qqplot_file_name)
GEM_Gmodel(snp_file_name, covariates_file_name, methylation_file_name, Gmodel_pv, Gmodel_result_file_name)
GEM_GxEmodel(snp_file_name, covariates_file_name_gxe, methylation_file_name, GxEmodel_pv, GxEmodel_result_file_name)

library(meffil)
options(mc.cores=1)
setwd("#####/methy_analysis_directory")  # Updated directory for methylation data

samplesheet <- meffil.create.samplesheet("#####/IDAT_files", recursive=TRUE)
meffil.list.cell.type.references()
qc.objects <- meffil.qc(samplesheet, cell.type.reference="dlpfc_reference", verbose=TRUE)
save(qc.objects, file="#####/output/rush.qc.objects.Robj")
save(samplesheet, file="#####/output/rush.samplesheet.Robj")

length(qc.objects)
names(qc.objects)
names(qc.objects[[1]])
qc.objects[[1]]$sample.name

writeLines(meffil.snp.names(), con="#####/output/snp-names.txt")

qc.parameters <- meffil.qc.parameters(
  beadnum.samples.threshold = 0.1,
  detectionp.samples.threshold = 0.1,
  detectionp.cpgs.threshold = 0.1, 
  beadnum.cpgs.threshold = 0.1,
  sex.outlier.sd = 5,
  snp.concordance.threshold = 0.95,
  sample.genotype.concordance.threshold = 0.8
)

qc.summary <- meffil.qc.summary(qc.objects)
save(qc.summary, file="#####/output/qcsummary.Robj")
meffil.qc.report(qc.summary, output.file="#####/output/qcreport.html")

qc.objects <- meffil.remove.samples(qc.objects, qc.summary$bad.samples$sample.name)
save(qc.objects, file="#####/output/qc.objects.clean.Robj")

qc.summary <- meffil.qc.summary(qc.objects)
save(qc.summary, file="#####/output/qcsummary.clean.Robj")
meffil.qc.report(qc.summary, output.file="#####/output/qc-report.clean.html")

load("#####/output/qc.objects.clean.Robj")
length(qc.objects)
load("#####/output/qcsummary.clean.Robj")

y <- meffil.plot.pc.fit(qc.objects)
ggsave(y$plot, filename="#####/output/pc_fit_plot.pdf", height=6, width=6)

pc <- 20
norm.objects <- meffil.normalize.quantiles(qc.objects, number.pcs=pc)
save(norm.objects, file="#####/output/norm.obj.Robj")

norm.beta <- meffil.normalize.samples(norm.objects, cpglist.remove=qc.summary$bad.cpgs$name)
save(norm.beta, file="#####/output/norm.beta.Robj")

norm.summary <- meffil.normalization.summary(norm.objects)
meffil.normalization.report(norm.summary, output.file="#####/output/normalization_report.html")

cc <- t(sapply(qc.objects, function(obj) obj$cell.counts$counts))
write.table(cc, "#####/output/cell_counts.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

samplephone <- read.csv("#####/SamplePhenos.csv", header=TRUE)
colnames(cc)[1] <- 'barcodes'
phe <- merge(samplephone, cc, by.x='barcodes', by.y=0)
rownames(phe) <- phe$barcodes

ewas.ret <- meffil.ewas(norm.beta, variable=as.factor(phe$DX), covariates=phe[,c('covariate1', 'covariate2')], isva=TRUE, most.variable=20000)
save(ewas.ret, file="#####/output/ewas_results.Robj")

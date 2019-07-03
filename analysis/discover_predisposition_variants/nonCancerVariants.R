##### nonCancerVariants.R #####
# Kuan-lin Huang 2019

source("../global_aes_out.R")

# read in related files
fn = "../../data/charged.PCA.r1.TCGAbarcode.merge.exon.ALL.vcf.samples.expanded.AFcorrected.lowAF.sele.tsv"
variants = read.table(sep="\t",header=T,file=fn, stringsAsFactors=FALSE, quote = "",fill=TRUE)

# cancer predisposition genes for filtering out
gene_fn = "../../data/20160713_Rahman_KJ_KH_152_gene_table_list.txt"
predisposition_genes = as.vector(t(read.table(sep="\t",header=F,file=gene_fn, stringsAsFactors=FALSE, quote = "")))

# subset patients to the ones that passed the final QC
clin_f = "../../data/PanCan_ClinicalData_V4_wAIM_filtered10389.txt"
clin = read.table(header=T, quote = "", sep="\t", fill =T, file = clin_f, stringsAsFactors=FALSE)

variants$bcr_patient_barcode = substr(variants$Sample,1,12)
variants_pca = variants[variants$bcr_patient_barcode %in% clin$bcr_patient_barcode,]

##### Add readcounts and readcount filter #####
# use readcounts (obtained using bamreadcount) for the ref/alt allele in both the normal/tumor samples to ensure the variant is bona fide germline variant
all_rc_fn = "../../data/all.rc_anno"
all_rc = read.table(sep="\t", header=T, file=all_rc_fn , stringsAsFactors=FALSE, quote = "",fill=T)
all_rc[all_rc==-1]=NA # swap out -1 for NA to avoid confusion
all_rc$altBase = gsub("\\-","del",all_rc$altBase)
all_rc$altBase = gsub("\\+","ins",all_rc$altBase)# excel friendly format
variants_pca = merge(variants_pca,all_rc,all.x = T, by=c("Sample","HGVSg"))
variants_pca$normalVAF = variants_pca$normalAltCnt/variants_pca$normalDepth
variants_pca$tumorVAF = variants_pca$tumorAltCnt/variants_pca$tumorDepth

# plots to check
p = ggplot(data=variants_pca)
p = p + geom_point(aes(x=normalVAF,y=tumorVAF), alpha = 0.02,stroke=0)
p = p + theme_bw()
p = p + theme(axis.text.x = element_text(colour="black", size=8,angle=90,vjust=0.5))#,
p
fn = "out/normal_tumor_VAF_prefilter.pdf"
ggsave(file=fn, h=5,w=5,useDingbats=FALSE)

p = ggplot(data=variants_pca)
p = p + geom_point(aes(x=normalAltCnt,y=tumorAltCnt), alpha = 0.02,stroke=0)
p = p + theme_bw() + xlim(0,50) + ylim(0,50)
p = p + theme(axis.text.x = element_text(colour="black", size=8,angle=90,vjust=0.5))#,
p
fn = "out/normal_tumor_altcnt_prefilter.pdf"
ggsave(file=fn, h=5,w=5,useDingbats=FALSE)

variants_pca = variants_pca[is.na(variants_pca$normalAltCnt) | (variants_pca$normalAltCnt > 4),]
cat("Applying normal Alt readcount >= 5, # of variants left: ",nrow(variants_pca),"\n")
variants_pca = variants_pca[is.na(variants_pca$tumorAltCnt) | (variants_pca$tumorAltCnt > 4),]
cat("Applying tumor Alt readcount >= 5, # of variants left: ",nrow(variants_pca),"\n")

variants_pca = variants_pca[is.na(variants_pca$normalVAF) | (variants_pca$normalVAF >= 0.2),]
cat("Applying normal Alt readcount >= 0.2, # of variants left: ",nrow(variants_pca),"\n")
variants_pca = variants_pca[is.na(variants_pca$tumorVAF) | (variants_pca$tumorVAF >= 0.2),]
cat("Applying tumor Alt readcount >= 0.2, # of variants left: ",nrow(variants_pca),"\n")

# plots to check
p = ggplot(data=variants_pca)
p = p + geom_point(aes(x=normalAltCnt,y=tumorAltCnt), alpha = 0.02,stroke=0)
p = p + theme_bw() + xlim(0,50) + ylim(0,50)
p = p + theme(axis.text.x = element_text(colour="black", size=8,angle=90,vjust=0.5))#,
p
fn = "out/normal_tumor_altcnt_postfilter.pdf"
ggsave(file=fn, h=5,w=5,useDingbats=FALSE)

p = ggplot(data=variants_pca)
p = p + geom_point(aes(x=normalVAF,y=tumorVAF), alpha = 0.02,stroke=0)
p = p + theme_bw()
p = p + theme(axis.text.x = element_text(colour="black", size=8,angle=90,vjust=0.5))#,
p
fn = "out/normal_tumor_VAF_postfilter.pdf"
ggsave(file=fn, h=5,w=5,useDingbats=FALSE)

##### classify whether a variant is cancer-relevant #####
cancer_terms = c("tumor","cancer","neoplasia")

variants_pca$predisposition_gene = F
variants_pca$predisposition_gene[variants_pca$HUGO_Symbol %in% predisposition_genes] = TRUE
variants_pca$cancer_term_trait = FALSE
for (term in cancer_terms){
  variants_pca$cancer_term_trait[grep(term,tolower(variants_pca$ClinVar_Traits))] = TRUE
}
variants_pca$cancer_term_trait[grep("oma$",tolower(variants_pca$ClinVar_Traits))] = TRUE

cat("In 152 cancer predisposition gene vs cancer related traits\n")
table(variants_pca$predisposition_gene,variants_pca$cancer_term_trait)
variants_pca$cancer_related = F
variants_pca$cancer_related[variants_pca$predisposition_gene | variants_pca$cancer_term_trait] = T

# # quickly check
# table(variants_pca$ClinVar_Traits[variants_pca$cancer_term_trait])[table(variants_pca$ClinVar_Traits[variants_pca$cancer_term_trait])>3]

# variant frequency annotation
var_freq = data.frame(table(variants_pca$HGVSg))
colnames(var_freq) = c("HGVSg","Cohort_AC")
variants_pca_frq = merge(variants_pca,var_freq,by="HGVSg")
variants_pca_frq$Cohort_AF = variants_pca_frq$Cohort_AC/10389*2
cat("Summary of cohort allelic frequency\n")
summary(variants_pca_frq$Cohort_AF)

# only include those with ACMG classification pass LP or P
variants_pca_frq_strict = variants_pca_frq[variants_pca_frq$ACMG_Classification != "Uncertain Significance",]
cat("ACMG pathogenic or likely pathogenic, # of sample-variants left: ",nrow(variants_pca_frq_strict),"\n")

variants_pca_frq_rare = variants_pca_frq_strict[variants_pca_frq_strict$Cohort_AF < 0.01,]
cat("Applying cohort AF < 1%, # of sample-variants left: ",nrow(variants_pca_frq_rare),"\n")

variants_pca_frq_rare = variants_pca_frq_rare[(is.na(variants_pca_frq_rare$GMAF) | variants_pca_frq_rare$GMAF < 0.01) & (is.na(variants_pca_frq_rare$ExAC_adj_AF) | variants_pca_frq_rare$ExAC_adj_AF < 0.0005),]
cat("Applying 1000G  < 1% & ExAC AF < 0.05%, # of sample-variants left: ",nrow(variants_pca_frq_rare),"\n")

# output all these variants
tn = "out/variants_pca_frq_rare_all.txt"
write.table(variants_pca_frq_rare, quote=F, sep="\t", file = tn, row.names = F)

# anonymized unique variants
variants_pca_frq_anony = variants_pca_frq_rare[,-which(colnames(variants_pca_frq_rare) %in% c("Sample","bcr_patient_barcode")),]
variants_pca_frq_anony_uniq = variants_pca_frq_anony[!duplicated(paste(variants_pca_frq_anony$HUGO_Symbol,variants_pca_frq_anony$HGVSc)),]
cat("# of unique variants: ",nrow(variants_pca_frq_anony_uniq),"\n")
tn = "out/variants_pca_frq_rare_all_uniqVarOnly.txt"
write.table(variants_pca_frq_anony_uniq, quote=F, sep="\t", file = tn, row.names = F)

# non-cancer related variants
variants_pca_frq_rare_nonCancer = variants_pca_frq_rare[!variants_pca_frq_rare$cancer_related,]
cat("None-cancer related, # of variants left: ",nrow(variants_pca_frq_rare_nonCancer),"\n")
tn = "out/variants_pca_frq_rare_nonCancer.txt"
write.table(variants_pca_frq_rare_nonCancer, quote=F, sep="\t", file = tn, row.names = F)

variants_pca_frq_anony_uniq_nonCancer = variants_pca_frq_anony_uniq[!variants_pca_frq_anony_uniq$cancer_related,]
cat("# of unique non-cancer variants: ",nrow(variants_pca_frq_anony_uniq_nonCancer),"\n")
tn = "out/variants_pca_frq_rare_nonCancer_uniqVarOnly.txt"
write.table(variants_pca_frq_anony_uniq_nonCancer, quote=F, sep="\t", file = tn, row.names = F)

top_traits = names(table(variants_pca_frq_rare_nonCancer$ClinVar_Traits)[table(variants_pca_frq_rare_nonCancer$ClinVar_Traits)>10])
top_traits = top_traits[-grep("not provided",top_traits)]
variants_pca_frq_rare_nonCancer_top_traits = variants_pca_frq_rare_nonCancer[variants_pca_frq_rare_nonCancer$ClinVar_Traits %in% top_traits,]

var_counts = data.frame(table(variants_pca_frq_rare_nonCancer_top_traits$HUGO_Symbol,variants_pca_frq_rare_nonCancer_top_traits$ClinVar_Traits))
colnames(var_counts) = c("Gene","ClinVar_Traits","Count")
var_counts$Freq = var_counts$Count/10389
var_counts$ClinVar_Traits_Abbrev = substr(var_counts$ClinVar_Traits,1,20)
var_counts$Count[var_counts$Count ==0] = NA
var_counts$Count_plot = var_counts$Count
#var_counts$Count_plot[var_counts$Count_plot>15] = 15

getPalette = colorRampPalette(c("#FFFFFF","#fed976","#e31a1c"))
p = ggplot(data=var_counts)
p = p + geom_tile(aes(x=ClinVar_Traits_Abbrev, y=Gene, fill=Count), linetype="blank") + scale_fill_gradientn(name= "Count", colours=getPalette(100), na.value=NA, limit=c(0,NA))
p = p + geom_text(aes(x=ClinVar_Traits_Abbrev, y=Gene, label = Count), color="black", size=3)
p = p + labs(title = "Pathogenic or likely pathogenic variant counts in 10,389 TCGA cases")
p = p  + theme_bw() + theme_nogrid() +
  theme(axis.title = element_text(size=16), axis.text.x = element_text(colour="black", size=12, angle=90, vjust = 0.5), axis.text.y = element_text(colour="black", size=12),axis.ticks = element_blank())#element_text(colour="black", size=14))
p
fn = "out/top_clinvar_trait_variants_pancan.pdf"
ggsave(file=fn, h=7,w=7,useDingbats=FALSE)


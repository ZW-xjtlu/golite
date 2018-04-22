devtools::use_data(METTL3_TREW,overwrite = TRUE)

fl <- "http://geneontology.org/ontology/subsets/goslim_generic.obo"

slim_generic <- getOBOCollection(fl)

devtools::use_data(slim_generic, internal = TRUE, overwrite = TRUE)

library(GO.db)

GO_term_indx <- Term(GOTERM)

devtools::use_data(GO_term_indx, internal = TRUE, overwrite = TRUE)


getwd()

#Code for test
library(TxDb.Hsapiens.UCSC.hg19.knownGene)


#An index table that match geneID to GOslim

txdb = TxDb.Hsapiens.UCSC.hg19.knownGene
all_eids_hg19 <- names(genes(txdb))

set.seed(1)
eids_bg <- sample(all_eids_hg19,3500)
eids_sets <- lapply(1:10,function(x)sample(eids_bg,300))

gene_set = eids_sets
back_ground = eids_bg

Gene_ID = eids_bg

Gene_symbol = AnnotationDbi::select(OrgDB, keys = as.character(Gene_ID), columns = "SYMBOL", keytype = "ENTREZID")

Gene_ID = Gene_symbol$SYMBOL

gene2goslim(Gene_ID,"SYMBOL",OrgDB)

library(org.Hs.eg.db)
orgDb = org.Hs.eg.db
OrgDB = org.Hs.eg.db

min_bg_count = 1
max_bg_count = 1000

library(golite)
GOEA(gene_set = eids_sets,
     orgDb = org.Hs.eg.db,
     back_ground = eids_bg
     )


map_index <- gene2go( eids_bg ,OrgDB = org.Hs.eg.db, Category = "BP")

GOEA()

Exact_terms = map_index$GO

exact2slim(NA)



#GO slim mapping for each individual genes

#vi get indx map: BP
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(golite)
txdb = TxDb.Hsapiens.UCSC.hg19.knownGene
all_eids_hg19 <- names(exonsBy(txdb,by = "gene"))
map_index <- gene2go( all_eids_hg19 ,OrgDB = org.Hs.eg.db, Category = "BP")
GO_gene_lst <- split(map_index$GO,map_index$ENTREZID)
GOslim_gene_lst <- lapply(GO_gene_lst,exact2slim,Category = "BP")
saveRDS(GOslim_gene_lst,"GOslim_gene_lst_BP.rds")

#nohup Rscript goslim_BP.R > goslim_BP.out &
#nohup Rscript goslim_CC.R > goslim_CC.out &
#nohup Rscript goslim_MF.R > goslim_MF.out &

#scp zhen@10.7.6.53:/home/zhen/golite_data/All_GOslim.tar.gz /Users/zhenwei/Documents/GitHub

GOslim_gene_lst_BP <- readRDS("GOslim_gene_lst_BP.rds")
GOslim_gene_lst_CC <- readRDS("GOslim_gene_lst_CC.rds")
GOslim_gene_lst_MF <- readRDS("GOslim_gene_lst_MF.rds")

Trans_df <- function(x){
data.frame(
  ENTREZID = rep(names(x),elementNROWS(x)),
  GO_slim = unlist(x)
           )
}

GOslim_gene_lst_BP <- Trans_df(GOslim_gene_lst_BP)
GOslim_gene_lst_CC <- Trans_df(GOslim_gene_lst_CC)
GOslim_gene_lst_MF <- Trans_df(GOslim_gene_lst_MF)

devtools::use_data(GOslim_gene_lst_BP, internal = TRUE, overwrite = TRUE)
devtools::use_data(GOslim_gene_lst_CC, internal = TRUE, overwrite = TRUE)
devtools::use_data(GOslim_gene_lst_MF, internal = TRUE, overwrite = TRUE)


mean(sapply(GOslim_gene_lst_BP,anyNA)) #27% genes are NAs.

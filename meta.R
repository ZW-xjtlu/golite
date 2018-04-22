
library(GSEABase)

fl <- "http://geneontology.org/ontology/subsets/goslim_generic.obo"

slim_generic <- getOBOCollection(fl)

library(GO.db)

GO_term_indx <- Term(GOTERM)

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

devtools::use_data(GO_term_indx,
                     slim_generic,
                     GOslim_gene_lst_BP,
                     GOslim_gene_lst_CC,
                     GOslim_gene_lst_MF, internal = TRUE, overwrite = TRUE)

#GO slim mapping for each individual genes

#vi get indx map: BP
#library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#library(org.Hs.eg.db)
#library(golite)
#txdb = TxDb.Hsapiens.UCSC.hg19.knownGene
#all_eids_hg19 <- names(exonsBy(txdb,by = "gene"))
#map_index <- gene2go( all_eids_hg19 ,OrgDB = org.Hs.eg.db, Category = "BP")
#GO_gene_lst <- split(map_index$GO,map_index$ENTREZID)
#GOslim_gene_lst <- lapply(GO_gene_lst,exact2slim,Category = "BP")
#saveRDS(GOslim_gene_lst,"GOslim_gene_lst_BP.rds")
#nohup Rscript goslim_BP.R > goslim_BP.out &
#nohup Rscript goslim_CC.R > goslim_CC.out &
#nohup Rscript goslim_MF.R > goslim_MF.out &

#scp zhen@10.7.6.53:/home/zhen/golite_data/All_GOslim.tar.gz /Users/zhenwei/Documents/GitHub

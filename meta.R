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

library(org.Hs.eg.db)
orgDb = org.Hs.eg.db

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

map_index_cc <- gene2go( eids_bg ,OrgDB = org.Hs.eg.db, Category = "CC")
map_index_mf <- gene2go( eids_bg ,OrgDB = org.Hs.eg.db, Category = "MF")


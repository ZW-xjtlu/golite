devtools::use_data(METTL3_TREW,overwrite = TRUE)

fl <- "http://geneontology.org/ontology/subsets/goslim_generic.obo"

slim_generic <- getOBOCollection(fl)

devtools::use_data(slim_generic, internal = TRUE, overwrite = TRUE)

getwd()

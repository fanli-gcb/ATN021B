
# store chemical annotations
chem_annot <- read.table("ChemicalAnnotation.txt", header=T, as.is=T, sep="\t", comment.char="", quote="")
rownames(chem_annot) <- chem_annot$CHEM_ID
metabolon_map <- chem_annot[, c("CHEM_ID", "CHEMICAL_NAME", "SUB_PATHWAY", "SUPER_PATHWAY", "COMP_ID", "PLATFORM", "HMDB", "KEGG", "PUBCHEM")]
colnames(metabolon_map) <- c("CHEM.ID", "CHEMICAL.NAME", "SUB.PATHWAY", "SUPER.PATHWAY", "COMP_ID", "PLATFORM", "HMDB", "KEGG", "PUBCHEM")
metabolon_map$CHEMICAL.NAME <- gsub("\"", "", metabolon_map$CHEMICAL.NAME)
rownames(metabolon_map) <- metabolon_map$CHEMICAL.NAME


# load RData object
load("data_for_mixOmics.RData")

# remove drugs from df.multiomics[["BIOCHEMICAL"]]
metabolon <- df.multiomics[["BIOCHEMICAL"]]
colnames(metabolon) <- gsub("prime", "'", colnames(metabolon))

tmp <- metabolon_map[colnames(metabolon),]
to_remove <- rownames(tmp)[which(tmp[,"SUB.PATHWAY"] %in% c("Drug - Antibiotic", "Drug - Antiviral"))]
metabolon <- metabolon[, setdiff(colnames(metabolon), to_remove)]

# replace
df.multiomics[["BIOCHEMICAL"]] <- metabolon

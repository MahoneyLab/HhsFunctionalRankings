#This function gets GO and phenotype information from biomaRt
#gene.id is a vector of gene.ids of one of the types id.type
#lib is the biomart library.

gene.pheno.table <- function(gene.id, lib, 
id.type = c("external_gene_name", "ensembl_gene_id", "entrezgene_id")){

    id.type = id.type[1]
    sig.info <- getBM(c("external_gene_name", "phenotype_description", "goslim_goa_description"),
    id.type, gene.id, lib)

    get.info <- function(gene.name, info.table){
        gen.terms <- c("cytoplasm", "nucleus", "biological_process", "molecular_function",
        "cellular_component", "cell", "intracellular", "organell", "biosynthetic process")
        gene.locale <- which(info.table[,1] == gene.name)
        all.pheno <- paste(unique(info.table[gene.locale,2]), collapse = ", ")
        all.go <- unique(info.table[gene.locale,3])
        trimmed.go <- paste(setdiff(all.go, gen.terms), collapse = ", ")
        table.row <- c(gene.name, all.pheno, trimmed.go)
        return(table.row)
    }

    u_gene <- unique(sig.info[,1])
    gene.info <- t(sapply(u_gene, function(x) get.info(x, sig.info)))
    rownames(gene.info) <- NULL
    colnames(gene.info) <- c("gene.name", "phenotypes", "GO Slim")
    return(gene.info)

}
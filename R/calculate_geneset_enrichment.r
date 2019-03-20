#' Use MAGMA to test enrichment in a geneset
#'
#' Assumes that you have already run map.snps.to.genes()
#'
#' @param geneset Genes which are to be tested (as HGNC / MGI symbols)
#' @param gwas_sumstats_path Filepath of the summary statistics file
#' @param analysis_name Used in filenames which area created
#' @param upstream_kb How many kb upstream of the gene should SNPs be included?
#' @param downstream_kb How many kb downstream of the gene should SNPs be included?
#' @param genome_ref_path Path to the folder containing the 1000 genomes .bed files (which can be downloaded from https://ctg.cncr.nl/software/MAGMA/ref_data/g1000_eur.zip)
#' @param geneset_species Species name relevant to the genes in the geneset, i.e. "mouse" or "human"
#'
#' @return Filepath for the genes.out file
#'
#' @export
calculate_geneset_enrichment <- function(geneset,gwas_sumstats_path,analysis_name,upstream_kb=10,downstream_kb=1.5,genome_ref_path,geneset_species="mouse"){
    gwas_sumstats_path = path.expand(gwas_sumstats_path)
    magmaPaths = get.magma.paths(gwas_sumstats_path,upstream_kb,downstream_kb)
    
    # First, check that the genes are HGNC/MGI IDs
    if(geneset_species=="human"){
        if(sum(geneset %in% all_hgnc_wtEntrez$hgnc_symbol)<0.5){stop("Less than 50% of the geneset are recognised HGNC symbols. Have you entered them in the wrong format? Or wrong species?")}
        geneset_entrez = all_hgnc_wtEntrez[all_hgnc_wtEntrez$hgnc_symbol %in% geneset,]$entrezgene
    }else if(geneset_species=="mouse"){
        if(sum(geneset %in% One2One::ortholog_data_Mouse_Human$orthologs_one2one$mouse.symbol)<0.25){stop("Less than 25% of the geneset are recognised MGI symbols with 1:1 orthologs. Have you entered them in the wrong format? Or wrong species?")}
        geneset_m2h = One2One::ortholog_data_Mouse_Human$orthologs_one2one[One2One::ortholog_data_Mouse_Human$orthologs_one2one$mouse.symbol %in% geneset,]$human.symbol
        geneset_entrez = all_hgnc_wtEntrez[all_hgnc_wtEntrez$hgnc_symbol %in% geneset_m2h,]$entrezgene
    }
    
    ctRows = paste(c(analysis_name,geneset_entrez),collapse=" ")
    
    # Write genes covar file to disk
    geneCovarFile=tempfile()
    #write.table(quantDat2,file=geneCovarFile,quote=FALSE,row.names=FALSE,sep="\t")
    write.table(ctRows,file=geneCovarFile,quote=FALSE,row.names=FALSE,sep="\t",col.names=FALSE)
    
    magma_cmd = sprintf("magma --gene-results '%s.genes.raw' --set-annot '%s' --out '%s.%s'",magmaPaths$filePathPrefix,geneCovarFile,magmaPaths$filePathPrefix,analysis_name)
     
    print(magma_cmd)
    system(magma_cmd)

    #path = sprintf("%s.%s.sets.out",magmaPaths$filePathPrefix,analysis_name)
    path = sprintf("%s.%s.gsa.out",magmaPaths$filePathPrefix,analysis_name)
    res = read.table(path,stringsAsFactors = FALSE)
    colnames(res) = as.character(res[1,])
    res=res[-1,]

    return(res)
}

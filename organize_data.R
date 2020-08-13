library(dplyr)
library(biomaRt)
library(tidyr)
library(DESeq2)
library(IsoformSwitchAnalyzeR)
library(rentrez)
library(purrr)

# DTE ---------------------------------------------------------------------
# Load results from diff expression
load("../results/diff_exp/gene_rin_ph_diff.rda")
df_res_gene <- df_res
load("../results/diff_exp/tx_rin_ph_diff.rda")
df_res_tx <- df_res
load("../results/diff_exp/diff_tx_corrected.rda")
load("../results/diff_exp/canonical_df.rda")

# Create gene names dictionary
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
dict <- getBM(attributes = c("hgnc_symbol", "ensembl_gene_id", "entrezgene_id"), 
              mart = ensembl)

# Filter genes that did converge 
rem <- rownames(dds)[!mcols(dds)$betaConv]
df_res_gene <- df_res_gene[!df_res_gene$gene %in% rem,]
df_res_padj <- df_res_padj[df_res_padj$geneID %in% df_res_gene$gene,]
canonical_df$group <- gsub("_CTRL_", "_", canonical_df$group)
all(df_res_padj$geneID %in% canonical_df$gene) 

# If FALSE, there are tx belonging to non-coding genes in tx diff exp,
# Select only the protein coding genes
df_res_padj <- inner_join(df_res_padj, canonical_df, by = c("geneID" = "gene", "txID" = "tx", "group"))

# Remove transcript versions
df_res_tx$tx <- gsub("\\.\\d+", "", df_res_tx$tx)

# Organize dataframe ------------------------------------------------------
plot_table <- df_res_gene[df_res_gene$gene %in% df_res_padj$geneID, c("gene", "logFC", "group")]
plot_table <- unique(inner_join(plot_table, df_res_padj, by = c("gene" = "geneID", "group")))
colnames(plot_table) <- c("gene", "logFC_gene", "group", "tx", "padj_gene", "padj_tx", "canonical", "prop")
plot_table <- inner_join(plot_table, df_res_tx[, c("logFC", "tx", "group")], by = c("tx", "group"))
colnames(plot_table)[length(plot_table)] <- "logFC_tx"
plot_table <- left_join(plot_table, dict, by = c("gene" = "ensembl_gene_id"))
plot_table <- separate(plot_table, col = group, into = c("region", "gender"), sep = "_")
plot_table$signif <- ifelse(plot_table$padj_tx <= 0.05, "S", "NS")

plot_table$hgnc_symbol <- ifelse(plot_table$hgnc_symbol == "", plot_table$gene, plot_table$hgnc_symbol)

plot_table <- plot_table %>% 
  group_by(gene, region, gender) %>% 
  filter(!all(signif == "NS")) %>% 
  ungroup()


# DTU ---------------------------------------------------------------------
files <- list.files("../data/ISA/obj/")
isa_df <- do.call(rbind, lapply(files, function(x) {
  load(paste0("../data/ISA/obj/", x))
  extractTopSwitches(
    SwitchList_2, filterForConsequences = T, n = Inf, inEachComparison = T
  )
}))

male <- grepl("_CTRL_male", isa_df$condition_1) & grepl("_MDD_male", isa_df$condition_2)
female <- grepl("_CTRL_female", isa_df$condition_1) & grepl("_MDD_female", isa_df$condition_2)

isa_df <- isa_df[male | female,]
isa_df$gene_id <- gsub("\\.\\d+", "", isa_df$gene_id)
isa_df$group <- gsub("CTRL_", "", isa_df$condition_1)

isa_df %>% 
  #dplyr::filter(isoform_switch_q_value <= 0.05) %>% 
  dplyr::select(gene_id, group) %>% 
  unique() %>% 
  dplyr::count(group)

isa_df %>% 
  mutate(id = paste(gene_name, gene_id, sep = "_"),
         entrezgene_id = dict$entrezgene_id[match(isa_df$gene_id, dict$ensembl_gene_id)]) %>% 
  dplyr::select(gene_id, gene_name, entrezgene_id, group, id) -> isa_df

# Get descriptions for all genes
desc_list <- imap_chr(unique(c(isa_df$gene_name, plot_table$hgnc_symbol)), function(x, y){ 
  print(c(x, y))
  
  res <- tryCatch(entrez_summary("gene", dict$entrezgene_id[match(x, dict$hgnc_symbol)])$summary,
                  error = function(e) "No summary available for this gene in RefSeq.")
  if(res[1] == "") {
    res <- "No summary available for this gene in RefSeq."
  }
  return(res)
  
})
names(desc_list) <- unique(c(isa_df$gene_name, plot_table$hgnc_symbol))

# Save data to be used in the app
if(!dir.exists("data/")) {
  dir.create("data/")
}
save(plot_table, isa_df, desc_list, file = "data/plot_tables.rda")



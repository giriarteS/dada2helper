#!/usr/bin/env Rscript

library(optparse)



description = "Write the taxa table and sequences to file."

option_list = list(
  make_option(c("-i", "--input"), type = "character", default = NULL,
              help = "The file path for the sequence table in RDS format.",
              metavar = "input_filepath"),
  make_option(c("-o", "--output_dir"), type = "character", default = "./",
              help = paste0("Directory where the sequence table and sequences",
                            "(txt and fasta formats) will be saved ",
                            "(default = ./)."),
              metavar = "output_directory")
)

opt_parser = OptionParser(option_list = option_list, description = description)
opt = parse_args(opt_parser)

if (is.null(opt$input)) {
  print_help(opt_parser)
  stop("No arguments supplied.\n", call. = FALSE)
}

export_taxa_table_and_seqs = function(seqtab, file_seqtab, file_seqs) {
  seqtab.t = as.data.frame(t(seqtab))
  seqs = row.names(seqtab.t)
  row.names(seqtab.t) = paste0("Taxon_", 1:nrow(seqtab.t))
  outlist = list(data_loaded = seqtab.t)
  mctoolsr::export_taxa_table(outlist, file_seqtab)
  seqs = as.list(seqs)
  seqinr::write.fasta(seqs, row.names(seqtab.t), file_seqs)
}

# output paths
outdir = ifelse(
  substr(opt$output_dir, nchar(opt$output_dir), nchar(opt$output_dir)) == '/',
  opt$output_dir,
  paste0(opt$output_dir, '/')
)
if (!dir.exists(outdir)) dir.create(outdir)

seqtab = readRDS(opt$input)

export_taxa_table_and_seqs(seqtab,
                           paste0(outdir, "taxa_table.txt"),
                           paste0(outdir, "taxa_seqs.fa"))


#!/usr/bin/env Rscript

library(optparse)



description = "Remove chimeras and produce taxa table."

option_list = list(
  make_option(c("-i", "--input"), type = "character", default = NULL,
              help = "The file path for the sequence table in RDS format.",
              metavar = "input_filepath"),
  make_option(c("-o", "--output_dir"), type = "character", default = "./no_chim",
              help = paste0("Directory where the no chimera sequence table ",
                            "file will be saved (default = ./no_chim)."),
              metavar = "output_directory"),
  make_option(c("-v", "--verbose"), default = FALSE, action = 'store_true',
              help = "Print out extra info?")
)

opt_parser = OptionParser(option_list = option_list, description = description)
opt = parse_args(opt_parser)

if (is.null(opt$input)) {
  print_help(opt_parser)
  stop("No arguments supplied.\n", call. = FALSE)
}


seqtab = readRDS(opt$input)

seqtab.nochim = dada2::removeBimeraDenovo(seqtab, verbose = opt$verbose)

# print summary stats
message("ORIGINAL TABLE:\n")
message(paste0("Dimensions: ", paste(dim(seqtab), collapse = "; "), "\n"))
message("TABLE WITH CHIMERAS REMOVED:\n")
message(paste0("Dimensions: ", paste(dim(seqtab.nochim), collapse = "; "), "\n"))
message(paste0(
  "Proportion seqs that were not chimeras: ",
  sum(seqtab.nochim) / sum(seqtab),
  "\n"
))


# output as r data objects
outdir = ifelse(
  substr(opt$output_dir, nchar(opt$output_dir), nchar(opt$output_dir)) == '/',
  opt$output_dir,
  paste0(opt$output_dir, '/')
)
if (!dir.exists(outdir)) dir.create(outdir)
saveRDS(seqtab.nochim, file = paste0(outdir, "seqtab_nochim.RDS"))


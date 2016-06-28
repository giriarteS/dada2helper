#!/usr/bin/env Rscript

library(optparse)



description = "Merge paired end sequences and produce taxa table."

option_list = list(
  make_option(c("-i", "--input_dir_dada"), type = "character", default = NULL,
              help = "Directory with dada RDS objects.",
              metavar = "input_directory"),
  make_option(c("-j", "--input_dir_derep"), type = "character", default = NULL,
              help = "Directory with dereplicated RDS objects.",
              metavar = "input_directory"),
  make_option(c("-o", "--output_dir"), type = "character", default = "./merged",
              help = paste0("Directory where the merged seqs R files will be ",
                            "saved (default = ./merged)."),
              metavar = "output_directory"),
  make_option(c("-v", "--verbose"), type = "logical", default = FALSE,
              help = "Print out extra info?")
)

opt_parser = OptionParser(option_list = option_list, description = description)
opt = parse_args(opt_parser)

if (is.null(opt$input_dir_dada)) {
  print_help(opt_parser)
  stop("No arguments supplied.\n", call. = FALSE)
}


indir_dada = opt$input_dir_dada
indir_derep = opt$input_dir_derep
path_dada = ifelse(substr(indir_dada, nchar(indir_dada),
                          nchar(indir_dada)) == '/',
                   indir_dada,
                   paste0(indir_dada, '/'))
path_derep = ifelse(substr(indir_derep, nchar(indir_derep),
                          nchar(indir_derep)) == '/',
                   indir_derep,
                   paste0(indir_derep, '/'))
dada_files = list.files(path_dada)
for (i in seq_along(dada_files)) {
  if (i == 1)
    dadas1 = readRDS(paste0(path_dada, dada_files[i]))
  else if (i == 2)
    dadas2 = readRDS(paste0(path_dada, dada_files[i]))
}
derep_files = list.files(path_derep)
for (i in seq_along(derep_files)) {
  if (i == 1)
    dereps1 = readRDS(paste0(path_derep, derep_files[i]))
  else if (i == 2)
    dereps2 = readRDS(paste0(path_derep, derep_files[i]))
}

merged = dada2::mergePairs(dadas1, dereps1, dadas2, dereps2,
                           verbose = opt$verbose)

# output as r data objects
outdir = ifelse(
  substr(opt$output_dir, nchar(opt$output_dir), nchar(opt$output_dir)) == '/',
  opt$output_dir,
  paste0(opt$output_dir, '/')
)
if (!dir.exists(outdir)) dir.create(outdir)
saveRDS(merged, file = paste0(outdir, "merged.RDS"))
seqtab = dada2::makeSequenceTable(merged)
saveRDS(seqtab, file = paste0(outdir, "seqtab.RDS"))


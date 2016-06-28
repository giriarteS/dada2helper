#!/usr/bin/env Rscript

library(optparse)


## Quality filter sequences in fastq files

description = "Dereplicate sequences using dada2"

option_list = list(
  make_option(c("-i", "--input_dir"), type = "character", default = NULL,
              help = "Directory with filtered sequences in fastq format",
              metavar = "input_directory"),
  make_option(c("-o", "--output_dir"), type = "character", default = "./derep",
              help = paste0("Directory where the dereplicated files will be ",
                            "saved (default = ./derep)."),
              metavar = "output_directory"),
  make_option(c("-v", "--verbose"), default = FALSE, action = 'store_true',
              help = "Print out extra info?")
)

opt_parser = OptionParser(option_list = option_list, description = description)
opt = parse_args(opt_parser)

if (is.null(opt$input_dir)) {
  print_help(opt_parser)
  stop("No arguments supplied.\n", call. = FALSE)
}


get_sample_names_and_fps2 = function(indir) {
  path = ifelse(substr(indir, nchar(indir), nchar(indir)) == '/',
                indir, paste0(indir, '/'))
  fns = list.files(path)
  fastqs = fns[grepl(".fq.gz$|.fastq.gz$", fns)]
  # print(fastqs)
  fastqs = sort(fastqs) # Sort ensures forward/reverse reads are in same order
  fnFs = fastqs[grepl("__filtR1", fastqs)] # Just the forward read files
  fnRs = fastqs[grepl("__filtR2", fastqs)] # Just the reverse read files
  # Get sample names from the first part of the forward read filenames
  # ensure forward and reverse reads in same order
  if (!identical(sapply(strsplit(fnFs, "__filtR1"), `[`, 1),
                 sapply(strsplit(fnRs, "__filtR2"), `[`, 1))) {
    stop('Forward and reverse reads not sorted in same order. Try simpler ",
         "sample names.')
  }
  sample.names = sapply(strsplit(fnFs, "__filtR1"), `[`, 1)
  # ensure all sample names are unique
  if (length(unique(sample.names)) != length(sample.names)) {
    stop('Make sure all sample IDs are unique.')
  }
  # Fully specify the path for the fnFs and fnRs
  fnFs = paste0(path, fnFs)
  fnRs = paste0(path, fnRs)
  # print(sample.names)
  list(sample.names = sample.names, fnFs = fnFs, fnRs = fnRs)
}

names_fps = get_sample_names_and_fps2(opt$input_dir)

message("Dereplicating forward reads.")
derepFs = dada2::derepFastq(names_fps$fnFs, verbose = opt$verbose)
message("Dereplicating reverse reads.")
derepRs = dada2::derepFastq(names_fps$fnRs, verbose = opt$verbose)
# Name the derep-class objects by the sample names
names(derepFs) = names_fps$sample.names
names(derepRs) = names_fps$sample.names

# output as r data objects
outdir = ifelse(
  substr(opt$output_dir, nchar(opt$output_dir), nchar(opt$output_dir)) == '/',
  opt$output_dir,
  paste0(opt$output_dir, '/')
)
if (!dir.exists(outdir)) dir.create(outdir)
saveRDS(derepFs, file = paste0(outdir, "derepFs.RDS"))
saveRDS(derepRs, file = paste0(outdir, "derepRs.RDS"))

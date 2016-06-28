#!/usr/bin/env Rscript

library(optparse)


## Quality filter sequences in fastq files

description = paste0("Quality filter sequences in fastq files using dada2 ",
                     "paired filtering.")
option_list = list(
  make_option(c("-i", "--input_dir"), type = "character", default = NULL,
              help = "Directory with demultiplexed sequences in fastq format",
              metavar = "input_directory"),
  make_option(c("-f", "--fwd_sfx"), type = "character", default = "_R1",
              help = "Forward sequence files suffix (default = '_R1')",
              metavar = "_R1"),
  make_option(c("-r", "--rev_sfx"), type = "character", default = "_R2",
              help = "Reverse sequence files suffix (default = '_R2')",
              metavar = "_R2"),
  make_option(c("-t", "--fwd_trunc"), type = "integer", default = 0,
              help = paste0("Forward truncate length. How long do you want ",
                            "the filtered sequences to be? Default = no ",
                            "truncation")),
  make_option(c("-y", "--rev_trunc"), type = "integer", default = 0,
              help = paste0("Reverse truncate length. How long do you want ",
                            "the filtered sequences to be? Default = no ",
                            "truncation")),
  make_option(c("-e", "--maxEE"), type = "double", default = 2,
              help = paste0("The maximum expected errors per sequence ",
                            "allowed (dDefault = 2).")),
  make_option(c("-o", "--output_dir"), type = "character", default = NULL,
              help = paste0("Directory where the filtered sequences will be ",
                            "saved (default = <input_dir>_filtered)."),
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


get_sample_names_and_fps = function(indir, fwd_sfx, rev_sfx) {
  path = ifelse(substr(indir, nchar(indir), nchar(indir)) == '/',
                indir, paste0(indir, '/'))
  fns = list.files(path)
  fastqs = fns[grepl(".fq$|.fastq$|.fq.gz$|.fastq.gz$", fns)]
  # print(fastqs)
  fastqs = sort(fastqs) # Sort ensures forward/reverse reads are in same order
  fnFs = fastqs[grepl(fwd_sfx, fastqs)] # Just the forward read files
  fnRs = fastqs[grepl(rev_sfx, fastqs)] # Just the reverse read files
  # Get sample names from the first part of the forward read filenames
  sample.names = sapply(strsplit(fnFs, fwd_sfx), `[`, 1)
  # ensure all sample names are unique
  if(length(unique(sample.names)) != length(sample.names)) {
    stop('Make sure all sample IDs are unique.')
  }
  # Fully specify the path for the fnFs and fnRs
  fnFs = paste0(path, fnFs)
  fnRs = paste0(path, fnRs)
  # print(sample.names)
  list(sample.names = sample.names, fnFs = fnFs, fnRs = fnRs)
}

names_fps = get_sample_names_and_fps(opt$input_dir, opt$fwd_sfx, opt$rev_sfx)



# Make filenames for the filtered fastq files
if (is.null(opt$output_dir)) {
  outpref = ifelse(
    substr(opt$input_dir, nchar(opt$input_dir), nchar(opt$input_dir)) == '/',
    substr(opt$input_dir, 1, nchar(opt$input_dir) - 1),
    opt$input_dir
  )
  outdir = paste0(outpref, "_filtered")
} else {
  outdir = ifelse(
    substr(opt$output_dir, nchar(opt$output_dir), nchar(opt$output_dir)) == '/',
    opt$output_dir,
    paste0(opt$output_dir, '/')
  )
}
if (!dir.exists(outdir)) dir.create(outdir)
filtFs = paste0(outdir, "/", names_fps$sample.names, "__filtR1.fastq.gz")
filtRs = paste0(outdir, "/", names_fps$sample.names, "__filtR2.fastq.gz")
# Filters
for(i in seq_along(names_fps$fnFs)) {
  print(paste0(i, " of ", length(names_fps$fnFs)))
  dada2::fastqPairedFilter(
    c(names_fps$fnFs[i], names_fps$fnRs[i]),
    c(filtFs[i], filtRs[i]),
    truncLen = c(opt$fwd_trunc, opt$rev_trunc),
    maxN = 0,
    maxEE = opt$maxEE,
    truncQ = 2,
    compress = TRUE,
    verbose = opt$verbose
  )
}


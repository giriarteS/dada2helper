#!/usr/bin/env Rscript

library(optparse)


## generate fastq sequence quality plots

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
  make_option(c("-n", "--num_samples"), type = "integer", default = 2,
              help = paste0("the number of samples for which to create ",
                            "quality profile plots. Both forward and reverse ",
                            "read plots will be created for each sample.")),
  make_option(c("-o", "--output_dir"), type = "character", default = "./",
              help = "Directory where the plots will be saved",
              metavar = "output_directory")
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

if (is.null(opt$input_dir)) {
  print_help(opt_parser)
  stop("No arguments supplied.\n", call. = FALSE)
}

get_sample_names_and_fps = function(indir, fwd_sfx, rev_sfx) {
  path = ifelse(substr(indir, nchar(indir), nchar(indir)) == '/',
                indir, paste0(indir, '/'))
  fns = list.files(path)
  fastqs = fns[grepl(".fq$|.fastq$", fns)]
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

sample_idxs = sample(length(names_fps$sample.names), opt$num_samples)

outdir = ifelse(
  substr(opt$output_dir, nchar(opt$output_dir), nchar(opt$output_dir)) == '/',
  opt$output_dir,
  paste0(opt$output_dir, '/')
)
if (!dir.exists(outdir)) dir.create(outdir)
for(i in sample_idxs) {
  pref = paste0(outdir, "qualplot_", names_fps$sample.names[i])
  # print(paste0(pref, "_F.pdf"))
  pdf(file = paste0(pref, "_F.pdf"), width = 5, height = 5)
  print(dada2::plotQualityProfile(names_fps$fnFs[[i]]))
  dev.off()
  pdf(file = paste0(pref, "_R.pdf"), width = 5, height = 5)
  print(dada2::plotQualityProfile(names_fps$fnRs[[i]]))
  dev.off()
}

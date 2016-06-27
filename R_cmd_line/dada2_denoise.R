#!/usr/bin/env Rscript

library(optparse)



description = "Denoise sequences using dada2."

option_list = list(
  make_option(c("-i", "--input_dir"), type = "character", default = NULL,
              help = "Directory with dereplicated RData objects.",
              metavar = "input_directory"),
  make_option(c("-o", "--output_dir"), type = "character", default = "./dada",
              help = paste0("Directory where the dada files will be ",
                            "saved (default = ./dada)."),
              metavar = "output_directory"),
  make_option(c("-s", "--subset"), type = "logical", default = FALSE,
              help = paste0("Whether to parameterize initial error rates on a ",
                            "random 10% subset of the samples (default = ",
                            "FALSE).")),
  make_option(c("-t", "--threads"), type = "integer", default = 1,
              help = "The number of threads (cores) to use (default = 1).")
)

opt_parser = OptionParser(option_list = option_list, description = description)
opt = parse_args(opt_parser)

if (is.null(opt$input_dir)) {
  print_help(opt_parser)
  stop("No arguments supplied.\n", call. = FALSE)
}

indir = opt$input_dir
path = ifelse(substr(indir, nchar(indir), nchar(indir)) == '/',
              indir, paste0(indir, '/'))
derepFs = readRDS(paste0(path, "derepFs.RDS"))
derepRs = readRDS(paste0(path, "derepRs.RDS"))

if (!opt$subset) {
  message("Denoising forward reads.")
  dadaFs = dada2::dada(
    derepFs,
    err = NULL,
    selfConsist = TRUE,
    multithread = opt$threads
  )
  message("Denoising reverse reads.")
  dadaRs = dada2::dada(
    derepRs,
    err = NULL,
    selfConsist = TRUE,
    multithread = opt$threads
  )
} else {
  if (floor(0.1 * length(derepFs)) == 0) {
    stop("Too few samples to use '--subset'")
  } else {
    if (floor(0.1 * length(derepFs)) <= 5) {
      warning("Few samples, so '--subset' not recommended.")
    }
  }
  derepFs.sub = derepFs[sample(1:length(derepFs), floor(0.1 * length(derepFs)))]
  derepRs.sub = derepRs[sample(1:length(derepRs), floor(0.1 * length(derepRs)))]
  message("Paramaterizing error rates for forward reads.")
  dadaFs.sub = dada2::dada(
    derepFs.sub,
    err = NULL,
    selfConsist = TRUE,
    multithread = opt$threads
  )
  message("Paramaterizing error rates for reverse reads.")
  dadaRs.sub = dada2::dada(
    derepRs.sub,
    err = NULL,
    selfConsist = TRUE,
    multithread = opt$threads
  )
  message("Denoising forward reads with existing error rates.")
  dadaFs = dada2::dada(
    derepFs,
    err = dadaFs.sub[[1]]$err_out,
    selfConsist = TRUE,
    multithread = opt$threads
  )
  message("Denoising reverse reads with existing error rates.")
  dadaRs = dada2::dada(
    derepRs,
    err = dadaRs.sub[[1]]$err_out,
    selfConsist = TRUE,
    multithread = opt$threads
  )
}



# output as r data objects
outdir = ifelse(
  substr(opt$output_dir, nchar(opt$output_dir), nchar(opt$output_dir)) == '/',
  opt$output_dir,
  paste0(opt$output_dir, '/')
)
if (!dir.exists(outdir)) dir.create(outdir)
saveRDS(dadaFs, file = paste0(outdir, "dadaFs.RDS"))
saveRDS(dadaRs, file = paste0(outdir, "dadaRs.RDS"))

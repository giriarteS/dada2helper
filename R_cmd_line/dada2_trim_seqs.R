#!/usr/bin/env Rscript

library(optparse)



description = "Trim primer/adapter seqs from dada objects."

option_list = list(
  make_option(c("-i", "--input_dir"), type = "character", default = NULL,
              help = "Directory with dada RData objects.",
              metavar = "input_directory"),
  make_option(c("-o", "--output_dir"), type = "character",
              default = "./dada_trim",
              help = paste0("Directory where the trimmed dada files will be ",
                            "saved (default = ./dada_trim)."),
              metavar = "output_directory"),
  make_option(c("-f", "--fwd_trim_seq"), type = "character", default = NULL,
              help = paste0("The sequence to trim off of the right end of the ",
                            "forward sequence.")),
  make_option(c("-r", "--rev_trim_seq"), type = "character", default = NULL,
              help = paste0("The sequence to trim off of the right end of the ",
                            "reverse sequence."))
)

opt_parser = OptionParser(option_list = option_list, description = description)
opt = parse_args(opt_parser)

if (is.null(opt$input_dir)) {
  print_help(opt_parser)
  stop("No arguments supplied.\n", call. = FALSE)
}


trim_dada_seqs = function(dada, R_trim_pattern) {
  if (methods::is(dada, "dada")) {
    dada = list(dada)
  }
  out_dada = dada
  for (i in seq_along(dada)) {
    # note that 'N's are necessary to match the adapter in the middle of the
    # sequence
    trimmed = Biostrings::trimLRPatterns(
      Lpattern = "",
      Rpattern = paste0(
        R_trim_pattern,
        paste0("NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN",
               "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN",
               "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN",
               "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN")
      ),
      Biostrings::DNAStringSet(dada[[i]]$sequence),
      Rfixed = FALSE
    )
    message(paste0(names(dada)[i], ' median trimmed length: ',
                   mean(Biostrings::width(trimmed))))
    # Check if duplicate sequences after trimming, meaning that dada probably
    # split biological sequence variants into multiple seqs because of junk
    # bases after primer/adapters). Not necessarily a problem since
    # 'makeSequenceTable()' will merge duplicate sequences.
    if (any(duplicated(trimmed))) {
      warning(paste0(names(dada)[i],
                     ' has duplicate seqs. Results might be wrong.'))
    }
    out_dada[[i]]$sequence = as.character(trimmed)
    out_dada[[i]]$clustering$sequence = as.character(trimmed)
    names(out_dada[[i]]$denoised) = as.character(trimmed)
  }
  out_dada
}

indir = opt$input_dir
path = ifelse(substr(indir, nchar(indir), nchar(indir)) == '/',
              indir, paste0(indir, '/'))
dadaFs = readRDS(paste0(path, "dadaFs.RDS"))
dadaRs = readRDS(paste0(path, "dadaRs.RDS"))
message("Trimming forward reads")
dadaFs.trim = trim_dada_seqs(dadaFs, opt$fwd_trim_seq)
message("Trimming reverse reads")
dadaRs.trim = trim_dada_seqs(dadaRs, opt$rev_trim_seq)


# output as r data objects
outdir = ifelse(
  substr(opt$output_dir, nchar(opt$output_dir), nchar(opt$output_dir)) == '/',
  opt$output_dir,
  paste0(opt$output_dir, '/')
)
if (!dir.exists(outdir)) dir.create(outdir)
saveRDS(dadaFs.trim, file = paste0(outdir, "dadaFs_trim.RDS"))
saveRDS(dadaRs.trim, file = paste0(outdir, "dadaRs_trim.RDS"))

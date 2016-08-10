#!/usr/bin/env Rscript

library(optparse)


description = "Combine multiple dada2 taxa tables."

option_list = list(
  make_option(c("-i", "--input_tables"), type = "character", default = NULL,
              help = paste0("File paths for taxa tables to be combined ",
                            "separated by spaces")),
  make_option(c("-f", "--input_fasta"), type = "character", default = NULL,
              help = paste0("File paths for fasta files for taxa sequences ",
                            "separated by spaces. Must be in same order as ",
                            "Taxa tables.")),
  make_option(c("-o", "--output_dir"), type = "character", default = NULL,
              help = paste0("Directory where the combined taxa table and ",
                            "sequences will be written."),
              metavar = "output_directory")
)

opt_parser = OptionParser(option_list = option_list, description = description)
opt = parse_args(opt_parser)

if (is.null(opt$input_tables)) {
  print_help(opt_parser)
  stop("No arguments supplied.\n", call. = FALSE)
}


tablefiles = opt$input_tables
seqfiles = opt$input_fasta
tablefiles_sep = strsplit(tablefiles, " ")[[1]]
seqfiles_sep = strsplit(seqfiles, " ")[[1]]



combine_dada2_output = function(table_files, sequence_files){
  if (!identical(length(table_files), length(sequence_files))) {
    stop("Different number of taxa table and sequence input files.")
  }

  # get consensus taxa
  keys = vector("list", length(sequence_files))
  for (i in seq_along(sequence_files)){
    tmpseqs = Biostrings::readDNAStringSet(sequence_files[i])
    # if adding seq
    if (i > 1) {
      # give repeat seqs previous names and new seqs new names
      repseqs = tmpseqs[as.character(tmpseqs) %in% as.character(all_seqs)]
      new_names = names(all_seqs)[match(as.character(repseqs), as.character(all_seqs))]
      repseqs_names = data.frame(old_names = names(repseqs), new_names)
      newseqs = tmpseqs[!as.character(tmpseqs) %in% as.character(all_seqs)]
      # get last taxon ID
      maxnum = max(as.numeric(sapply(strsplit(
        names(all_seqs), 'Taxon_'
      ), "[", 2)))
      newnames =
        paste0("Taxon_", seq(from = maxnum + 1, to = maxnum + length(newseqs)))
      # make key to translate old names to new names
      keys[[i]] = rbind(repseqs_names,
                        data.frame(old_names = names(newseqs), 
                                   new_names = newnames))
      names(newseqs) = newnames
      # append new seqs to existing seqs
      all_seqs = Biostrings::append(all_seqs, newseqs)
      if (length(unique(names(all_seqs))) != length(names(all_seqs))) {
        stop("A problem occurred.")
      }
    } else {# if first seq set
      all_seqs = tmpseqs
      keys[[i]] = data.frame(old_names = names(tmpseqs), new_names = names(tmpseqs))
    }
  }

  # rename taxa in taxa tables
  tables = vector("list", length(table_files))
  for (i in seq_along(table_files)) {
    tmptab = read.delim(table_files[i], skip = 1, check.names = FALSE)
    tmptab$`#OTU ID` =
      keys[[i]]$new_names[match(tmptab$`#OTU ID`, keys[[i]]$old_names)]
    tables[[i]] = tmptab
  }

  # merge tables
  ## check that sample IDs don't overlap
  sampleIDs = unlist(sapply(tables, colnames))[unlist(sapply(tables, colnames)) != "#OTU ID"]
  if (length(unique(sampleIDs)) != length(sampleIDs)) {
    stop("Duplicate sample IDs exist. Make sure each sample ID is unique",
         "across taxa tables.")
  }
  merged_tables = Reduce(function(...) merge(..., all = TRUE, sort = FALSE), tables)
  merged_tables[is.na(merged_tables)] = 0

  list(combined_taxa_table = merged_tables, combined_seqs = all_seqs)
}



combined_output = combine_dada2_output(tablefiles_sep, seqfiles_sep)

# write output
outdir = ifelse(
  substr(opt$output_dir, nchar(opt$output_dir), nchar(opt$output_dir)) == '/',
  opt$output_dir,
  paste0(opt$output_dir, '/')
)
if (!dir.exists(outdir)) dir.create(outdir)

Biostrings::writeXStringSet(combined_output$combined_seqs,
                            paste0(outdir, "taxa_seqs.fa"))


write("# Taxa table", paste0(outdir, "taxa_table.txt"))
suppressWarnings(
  write.table(
    combined_output$combined_taxa_table,
    paste0(outdir, "taxa_table.txt"),
    row.names = FALSE,
    sep = '\t',
    append = TRUE
  )
)



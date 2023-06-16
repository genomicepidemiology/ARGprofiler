#!/usr/bin/env Rscript
# MapStatFilters - A script for filtering a mapstat file for unwanted hits / references that do not meet certain conditions
# Also adds a combination of derived statistics 
# Author: Patrick Munk
# Edited: March 28, 2023
# Example use : Rscript.exe mapstatFilters.R -i example_data\SRR1027651.mapstat -o thisIsTestData/filtered.flt.out.mapstat -r .\example_data\ResFinder_20200125.refdata

# TODO
# Ensure safety with empty mapstat files, refdata files etc
# Add filter numbers as proper options

# Generates python-like options used when the Rscript from a terminal. Requires the R-package 'optparse'
library("optparse")
option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="An input .mapstat file from read alignment that needs filtering", metavar="character"),
  make_option(c("-o", "--outputr"), type="character", default=NULL,
              help="An output .mapstat file filtered", metavar="character"),
  make_option(c("-r", "--refdata"), type="character", default=NULL,
              help="A refdata file or similar where 1st and 2nd columns are reference names and lengths", metavar="character"),
  make_option(c("-d", "--depth"), type="double", default=6,
              help="%default", metavar="Mean depth coverage of gene"));

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# Set user input - for testing
#opt$input = "data/mapstatfiles/SRR1027651.mapstat"
#opt$output = "mapstatfilters/SRR1027651.flt.mapstat"
#opt$refdata = "data/ResFinder_20200125.refdata"

# More variables
uservars = list()
uservars$outputDir = dirname(opt$output)
uservars$outFileName = basename(opt$output)

################# PROGRAM START #################

# Create output directory if it does not exist
if (dir.exists(uservars$outputDir) == F) {
  print("Creating the output dir")
  dir.create(uservars$outputDir, recursive = T)
}

# Import mapstat file and split into header, colnames and data
mapstat = readLines(opt$input, encoding = "UTF-8")

ms_header = mapstat[grep("##", mapstat)]
ms_colnames = unlist(strsplit(mapstat[length(ms_header)+1], "\t"))
ms_data = read.delim(textConnection(mapstat), comment.char = "#", h = F, 
                     col.names = ms_colnames, check.names = F)

# Read in the refdata / onject with genes and gene lengths
refdata_width = ncol(read.delim(opt$refdata, nrows =  1, comment.char = "#", sep = "\t"))
refdata = read.delim(opt$refdata, comment.char = "#", sep = "\t", h = F,
                     colClasses = c("character", "integer", 
                                    rep("NULL", refdata_width - 2)))
colnames(refdata) = c("# refSequence", "length")
      
# Merge data               
ms_data2 = merge(ms_data, refdata, all.x = T)

# Detect if we have any mapstat genes unmatched by refdata length
# Stop program with error? Warning?
num_unmatched_genes = sum(is.na(ms_data2$length))
if (num_unmatched_genes > 0) {
  print(paste("WARNING:", num_unmatched_genes, "gene lengths not found in refdata!"))
  # Remove NAs as we cant evaluate those
  ms_data2 = na.omit(ms_data2)
}



ms_newcols = with(ms_data2, data.frame(meanDepthCovered = bpTotal / length,
                                       spuriosCovRatio = nucHighDepthVariance / refCoveredPositions,
                                       FragReadRatio = fragmentCountAln / readCountAln,
                                       meanMapScore = mapScoreSum / bpTotal,
                                       readConsensRefIdentity = refConsensusSum /refCoveredPositions,
                                       readRefIdentity = (bpTotal - snpSum - insertSum - deletionSum) / bpTotal,
                                       depthCV = depthVariance,
                                       propCovered = refCoveredPositions / length))

ms_data2 = cbind(ms_data2, ms_newcols)


FilterMapstat = function(ms_dat, prop_cov, read_id, cons_id, spur_ratio, mean_depth) {
  ms_data_flt = ms_dat[ms_dat$propCovered > prop_cov &
                         ms_dat$readRefIdentity > read_id &
                         ms_dat$readConsensRefIdentity > cons_id &
                         ms_dat$spuriosCovRatio < spur_ratio &
                         ms_dat$meanDepthCovered > mean_depth,]
}

ms_data_flt = FilterMapstat(ms_data2, 0.9, 0.9, 0.9, 0.90, opt$depth)

WriteMapstatFile = function(ms_header, ms_data, file_path) {
  # Open a file for writing
  outFile = file(file_path, "w")
  # Add the mapstat header
  writeLines(ms_header, outFile)
  # Add the main body to file and close
  #writeLines(ms_data, outFile)
  #write.table(ms_data, outFile, append = T, sep = "\t")
  #write.table(ms_data, "test.outfile.txt", append = T, sep = "\t")
  suppressWarnings(write.table(ms_data_flt, outFile, sep="\t", 
              append=TRUE, 
              quote = F, row.names = F))
  close(outFile)
}

# Export filtered mapstat file
WriteMapstatFile(ms_header = ms_header, 
                 ms_data = ms_data_flt, 
                 file_path = opt$output)

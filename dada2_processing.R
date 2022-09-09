# dada2_processing.R

library(dada2.pipeline)
library(data.table)
library(magrittr)
library(stringr)

###### Assign/change these variables as needed ######
raw.seq.dirs <- c(
  Run1 = ___,
  Run2 = ___
) # you won't necessarily have more than one directory, but I put this here as a example of how you can process multiple runs together
inDir <- "Input" # directory for all input files (not including raw sequence files)
                 # this should include the files you will assign to `sample.ids.filename` and
                 # `metadata.filename`
sample.ids.filename <-  ___ # for properly grabbing sequence files and assigning the sample
                            # names to be used going forward.
                            # this *may* be the same file as `metadata.filename` below
                            # if you have sample across multiple sequencing runs,
                            # make sure you have a column named `Run` in this table with ids                                 # equivalent to those you use starting up on line 9
# A multi-run sample ids table should look something like this:
# | Sample | FileID   | RunID |
# |--------|----------|-------|
# | s001   | ACAGCGTA | Run1  |
# | s002   | CAGAGTAC | Run1  |
#    ...       ...      ...
# | s100   | CACAGAGA | Run2  |
#    ...       ...      ...

file.ids.col <- ___  # name of column containing the file ID patterns to grab proper files
                     # these IDs *can* be identical to the samples name if that is how the
                     # raw sequence files are labelled
                     # it may also be, e.g., the barcode sequences for each sample
run.ids.col <- ___  # name of column containing the run IDs to grab proper files
                    # if processing a single run, set to NULL
metadata.filename <- ___ # sample metadata to be used in endpoint phyloseq object
                         # this *may be the same file as `sample.ids.filename` above
# If this file will be used as the sample ids file,
# then the above types of columns should be include here as well.
# A metadata table should look something like this:
# | Sample | Covariate1 | Covariate2 | ...
# |--------|------------|------------| ...
# | s001   | wildtype   |      89.30 | ...
# | s002   | knockout   |      32.23 | ...
#    ...        ...             ...

smpl.ids.col <- "Sample" # name of column containing names of samples in
                     # `sample.ids.filename` and `metadata.filename` files
taxa.db.path <- ___ # path to database to be used for assigning taxonomy
                    # this should be a gzipped FASTA file
fasttree.path <- ___ # path to the FastTree executable
cleanup <- TRUE
cores <- 50
split <- "--" # the string used to split the sample name from the read id in the symlink
              # file names: e.g. sample1--R1.fastq.gz

###### Begin script ######
if (!dir.exists(inDir)) { dir.create(inDir) }
sample.ids.dt <- read.file(file.path(inDir, sample.ids.filename))
if (!("data.table" %in% sample.ids.dt)) {
  sample.ids.dt <- as.data.table(sample.ids.dt) %>%
    setkeyv(smpl.ids.col)
}
orig.dir <- getwd()
fastq.dir <- file.path(inDir, "FastQs")
setDTthreads(threads = cores)
if (!dir.exists(fastq.dir)) { dir.create(fastq.dir) }

###### Begin processsing ######
initiate.pipeline()
setwd(fastq.dir)

###### Symlink FASTQs ######
symlink.fastqs(
  seq.dirs = raw.seq.dirs,
  ids.tbl = sample.ids.dt,
  smpl.id.col = smpl.ids.col,
  file.id.col = file.ids.col,
  run.id.col = run.ids.col,
  split.pattern = split,
  quiet = FALSE
)

setwd(orig.dir)

###### Run dada2 up to quality plots ######
dada2.upto.qualPlots(
  fastq.path = fastq.dir,
  file.split.pattern = split,
  maxCores = cores,
  random.seed = 42
)

### If working on a cluster you may need to move PDFs to a local machine, which may be faster if they're in your home directory
### If so, uncomment below
# system(paste0("cp -v dada2*", Sys.Date(), "_output/*.pdf ~/temp"))

### Look at quality PDFs
### Quality of reads looks ___

### Proceeding with ___ reads # (forward only or both)

###### Finish dada2 processing ######
dada2.finish(
  fastq.path = fastq.dir,
  truncFwd.len = ___,
  truncRev.len = ___,
  taxa.db = taxa.db.path,
  metadata.file = file.path(inDir, metadata.filename),
  paired = TRUE,
  maxCores = cores,
  build.tree = TRUE,
  fasttree.path = fasttree.path
)

###### Move files to proper directories ######
file.copy(
  from = file.path(run.env$output.path, "phyloseq.rds"),
  to = file.path(inDir, "phyloseq.rds"),
  overwrite = TRUE
)

###### Compress and archive output, remove unneeded files ######
if (cleanup) {
  system(paste0("tar zvcf ", run.env$output.path, ".tgz ", run.env$output.path))
  system(paste("rm -r", run.env$output.path))
  system(paste("rm -r", fastq.dir))
}

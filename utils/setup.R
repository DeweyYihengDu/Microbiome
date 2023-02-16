# ##############################################################################
#
## Setup the repository
#
# ##############################################################################

library("tidyverse")
library("rjson")
library("tools")
library("filesstrings")
library("here")

dir.create(here('temp'))

# ##############################################################################
# download CRC meta dataset for prevalence calculation
crc.meta <- fromJSON(file=here('files', 'zenodo_crc.json')) %>% 
  pluck('files') %>% 
  bind_rows() %>% 
  filter(key %in% c('meta_all.tsv', 'species_profiles_g2_l75_motus2.0.0.tsv'))


for (i in seq_len(nrow(crc.meta))){
  options(timeout = max(700, getOption("timeout")))
  x <- crc.meta[i,]
  download.file(x$links[[1]], destfile = here('temp', x$key))
  md5ref <- str_remove(x$checksum, 'md5:')
  md5file <- md5sum(here('temp', x$key))
  if (md5ref != md5file){
    stop('MD5 sums did not agree for: ', x$key)
  } else {
    file.move(here('temp', x$key), here('data'))
  }
}


# ##############################################################################
# download data for this actual project
data.files <- fromJSON(file=here('files', 'zenodo_main.json')) %>% 
  pluck('files') %>% 
  bind_rows()

for (i in seq_len(nrow(data.files))){
  dir.create(here('data', 'annotated_CBs'))
  x <- data.files[i,]
  download.file(x$links[[1]], destfile = here('temp', x$key))
  md5ref <- str_remove(x$checksum, 'md5:')
  md5file <- md5sum(here('temp', x$key))
  if (md5ref != md5file){
    stop('MD5 sums did not agree for: ', x$key)
  } else {
    if (str_detect(x$key, '^P[0-9]_')){
      file.move(here('temp', x$key), here('data', 'annotated_CBs/'))
    } else {
      file.move(here('temp', x$key), here('data'))
    }
  }
}


# ##############################################################################
# clean-up
unlink(here('temp'))

message('Script finished without errors!')

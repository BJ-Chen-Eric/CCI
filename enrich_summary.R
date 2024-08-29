suppressWarnings({
  suppressPackageStartupMessages(library(argparse))
  suppressPackageStartupMessages(library(data.table))
  suppressPackageStartupMessages(library(dplyr))
  suppressPackageStartupMessages(library(stringr))
})



parser <- ArgumentParser()
parser$add_argument("-p", "--pathway_summary", type='character', default=TRUE,
                    help="Requied input pathway_summary")

parser$add_argument("-g", "--gene_list", type='character', default=TRUE,
                    help="Requied input pathway_summary")

args <- parser$parse_args()

# engo <- fread('/home/eric/Analysis/spatial/plot/old_10x/deg300/wi_en_summary', header = F, sep = '') %>% as.character()
engo <- fread(args$pathway_summary, header = F, sep = '') %>% as.character()

engo <- sub(engo, replacement = '', pattern='The subnetwork shows the following associations: .*: ')

engo <- strsplit(engo, split = '\\. ')[[1]] %>% as.list()

p <- lapply(engo, function(x){str_extract(string = x, pattern = 'belong to the .*')})
names(engo) <- lapply(p, function(x){sub(x, replacement = '', pattern='belong to the biological process ')}) %>% unlist()
engo <- lapply(engo, function(x){sub(x, replacement = '', pattern = ' belong to the .*')})
engo <- lapply(engo, function(x){sub(x, replacement = '', pattern = 'and ')})
engo <- lapply(engo, function(x){strsplit(x, split = ', ') %>% unlist()})

test <- engo %>% unlist() %>% as.data.frame()
type <- fread(args$gene_list, header = F) %>% as.data.frame()

for(i in 1:nrow(type))  {
  row <- test[,1] %in% type[i, 1] %>% which()
  if(identical(row, integer(0)) == F)  {
    test[row, 'type'] <- type[i, 2]
  }
}

r <- table(complete.cases(test[,2])) %>% as.data.frame() %>% nrow()
if(r!=1) {
  stop("Wrong pathway to gene profile")
}

path <- paste(args$pathway_summary, '_edit.csv', sep = '')
write.table(test, path, quote = F, sep='\t')

cat(paste0("save edited summary", args$out, "\n"))


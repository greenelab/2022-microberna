library(dplyr)
library(tidyr)

read_cdhit_clstr <- function(filepath) {
  clstr <- read.table(filepath, sep = "\t", comment.char = "", quote = "", fill = T, stringsAsFactors = F, col.names = c("col1", "col2")) %>%
    separate(col1, into = c("sequence_num_in_cluster", "cluster"), sep = " ", fill = "right") %>%
    fill(cluster) %>%
    filter(!grepl(">", sequence_num_in_cluster)) %>%
    separate(col2, into = c("sequence_length", "col2"), sep = "nt, >") %>%
    extract(col2, into = c("sequence_name", "representative_sequence", "col2"), regex = "(.*?)[.]{3} ([*]|at) ?(.*)") %>%
    mutate(representative_sequence = representative_sequence == "*", col2 = ifelse(representative_sequence, "100%", col2)) %>%
    group_by(cluster) %>%
    mutate(representative = sequence_name[which(representative_sequence)]) %>%
    separate(col2, into = c("strand", "p_identity"), sep = "/", fill = "left", convert = T) %>%
    mutate(p_identity = sub("%", "", p_identity) %>% as.numeric)
  return(clstr)
}

read_gff <- function(filepath){
  # adapted from https://gist.githubusercontent.com/acvill/03343034392cff158d2369483ed8935f/raw/6ba9df8f2a27d5a4789475668c353d339fe45ca0/gtf2tibble.R
  require(tidyverse)
  
  cnames <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")
  
  # for files that contain fasta sequences at the bottom, determine where to stop reading
  n_max_num <- grep("##FASTA", read_lines(filepath))
  # tidyverse doesn't count comment characters toward line limit.
  # determine correction by counting number of commented lines at beginning of file
  n_max_num <- n_max_num - length(grep("^#", read_lines(filepath)))

  # read in raw gff as tsv and remove comment rows
  if(length(n_max_num) == 0){
    gff_raw <- read_tsv(filepath, col_types =  "cccddcccc", col_names = cnames, comment = "#") 
  } else {
    gff_raw <- read_tsv(filepath, col_types =  "cccddcccc", col_names = cnames, comment = "#", n_max = n_max_num) 
  }
  
  # get the unique attribute types.
  # this assumes there are no spaces in the attribute names.
  # it first parses on the semicolon, and then removes the unique attributes
  # using the = character
  attribute_names <- gff_raw %>%
    select(attribute) %>%
    apply(., MARGIN = 1, FUN = str_split, pattern = ';') %>%
    unlist() %>% 
    trimws() %>% 
    trimws(whitespace = ";") %>%
    sub("=.*$", "", .) %>% 
    unique()
  
  # catch empty strings
  attribute_names <- attribute_names[attribute_names != ""]
  # catch NA strings
  attribute_names <- attribute_names[!is.na(attribute_names)]
  
  # for each attribute type, create column
  # apply over gff to fill in rows where attribute type is found
  for (attribute in attribute_names) {
    #attribute = "ID"
    collect_attributes <- apply(gff_raw, MARGIN = 1, function(x) {
      variable <- str_extract(string = x[9], pattern = paste0(attribute, "=.*?;")) %>% 
        trimws(whitespace = '[";]+', which = 'left') %>% 
        str_extract('(?<==)[^;]+')
    })
    gff_raw <- gff_raw %>% add_column("{attribute}" := collect_attributes)
  }

  # remove original attribute column
  gff_parsed <- gff_raw %>% select(-c(attribute))
  return(gff_parsed)
}

read_gather <- function(path){
  gather <- read_csv(path, col_types = "dddddlllcccddddcccd")
}

# add an R2 and eq to a linear regression ---------------------------------

# From https://gist.github.com/kdauria/524eade46135f6348140
# and  https://gist.github.com/dsaiztc/e615740359fdb8532aca485b5ee990aa
stat_smooth_func <- function(mapping = NULL, data = NULL,
                             geom = "smooth", position = "identity",
                             ...,
                             method = "auto",
                             formula = y ~ x,
                             se = TRUE,
                             n = 80,
                             span = 0.75,
                             fullrange = FALSE,
                             level = 0.95,
                             method.args = list(),
                             na.rm = FALSE,
                             show.legend = NA,
                             inherit.aes = TRUE,
                             xpos = NULL,
                             ypos = NULL) {
  layer(
    data = data,
    mapping = mapping,
    stat = StatSmoothFunc,
    geom = geom,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      method = method,
      formula = formula,
      se = se,
      n = n,
      fullrange = fullrange,
      level = level,
      na.rm = na.rm,
      method.args = method.args,
      span = span,
      xpos = xpos,
      ypos = ypos,
      ...
    )
  )
}


StatSmoothFunc <- ggproto("StatSmooth", Stat,
                          setup_params = function(data, params) {
                            # Figure out what type of smoothing to do: loess for small datasets,
                            # gam with a cubic regression basis for large data
                            # This is based on the size of the _largest_ group.
                            if (identical(params$method, "auto")) {
                              max_group <- max(table(data$group))
                              
                              if (max_group < 1000) {
                                params$method <- "loess"
                              } else {
                                params$method <- "gam"
                                params$formula <- y ~ s(x, bs = "cs")
                              }
                            }
                            if (identical(params$method, "gam")) {
                              params$method <- mgcv::gam
                            }
                            
                            params
                          },
                          
                          compute_group = function(data, scales, method = "auto", formula = y~x,
                                                   se = TRUE, n = 80, span = 0.75, fullrange = FALSE,
                                                   xseq = NULL, level = 0.95, method.args = list(),
                                                   na.rm = FALSE, xpos=NULL, ypos=NULL) {
                            if (length(unique(data$x)) < 2) {
                              # Not enough data to perform fit
                              return(data.frame())
                            }
                            
                            if (is.null(data$weight)) data$weight <- 1
                            
                            if (is.null(xseq)) {
                              if (is.integer(data$x)) {
                                if (fullrange) {
                                  xseq <- scales$x$dimension()
                                } else {
                                  xseq <- sort(unique(data$x))
                                }
                              } else {
                                if (fullrange) {
                                  range <- scales$x$dimension()
                                } else {
                                  range <- range(data$x, na.rm = TRUE)
                                }
                                xseq <- seq(range[1], range[2], length.out = n)
                              }
                            }
                            # Special case span because it's the most commonly used model argument
                            if (identical(method, "loess")) {
                              method.args$span <- span
                            }
                            
                            if (is.character(method)) method <- match.fun(method)
                            
                            base.args <- list(quote(formula), data = quote(data), weights = quote(weight))
                            model <- do.call(method, c(base.args, method.args))
                            
                            m = model
                            eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                                             list(a = format(coef(m)[1], digits = 2), 
                                                  b = format(coef(m)[2], digits = 2), 
                                                  r2 = format(summary(m)$r.squared, digits = 2)))
                            func_string = as.character(as.expression(eq))
                            
                            if(is.null(xpos)) xpos = min(data$x)*0.95
                            if(is.null(ypos)) ypos = max(data$y)*0.95
                            data.frame(x=xpos, y=ypos, label=func_string)
                            
                          },
                          
                          required_aes = c("x", "y")
)

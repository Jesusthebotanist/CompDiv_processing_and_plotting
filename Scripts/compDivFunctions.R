# Dependencies 
library(dplyr)
library(readr)
library(coda)
library(tibble) 
library(ggtree)
library(ggplot2)
library(reshape2) # melt()
library(cowplot) # plot_grid()

# PreProcessing Scripts

# "stolen" with permission from RevGadgets
getMAP <- function(var) {
  d <- density(var)
  f <- approxfun(d$x, d$y)
  op <- stats::optim(
    par = mean(var),
    fn = f,
    method = "SANN",
    control = list(fnscale = -1)
  )
  return(op$par)
}

# "stolen" with permission from RevGadgets
matchNodes <- function(phy) {
  # get some useful info
  num_tips <- length(phy$tip.label)
  num_nodes <- phy$Nnode
  tip_indexes <- 1:num_tips
  node_indexes <- num_tips + num_nodes:1
  
  node_map <-
    data.frame(R = 1:(num_tips + num_nodes),
               Rev = NA,
               visits = 0)
  current_node <- phy$Nnode + 2
  k <- 1
  t <- 1
  
  while (TRUE) {
    if (current_node <= num_tips) {
      node_map$Rev[node_map$R == current_node] <- t
      current_node <- phy$edge[phy$edge[, 2] == current_node, 1]
      t <- t + 1
    } else {
      if (node_map$visits[node_map$R == current_node] == 0) {
        node_map$Rev[node_map$R == current_node] <- node_indexes[k]
        k <- k + 1
      }
      node_map$visits[node_map$R == current_node] <-
        node_map$visits[node_map$R == current_node] + 1
      
      if (node_map$visits[node_map$R == current_node] == 1) {
        # go right
        current_node <- phy$edge[phy$edge[, 1] == current_node, 2][2]
      } else if (node_map$visits[node_map$R == current_node] == 2) {
        # go left
        current_node <- phy$edge[phy$edge[, 1] == current_node, 2][1]
      } else if (node_map$visits[node_map$R == current_node] == 3) {
        # go down
        if (current_node == num_tips + 1) {
          break
        } else {
          current_node <- phy$edge[phy$edge[, 2] == current_node, 1]
        }
      }
    }
  }
  
  return(node_map[, 1:2])
  
}


# MTBD/LSBDS

LSBDSMTBD_combine_burin <- function(mcmc_chain_directory,burnin=.05, MTBDorLSBDS = NULL){
  # mcmc_chain_directory - the path to a directory that contains only the rates chains.
  # burnin - default is 0.05
  # Note we use coda::mcmc.list() to combine chains. This requires chins to have
  # the same number of generations. If they do not this code will find the chian
  # with the lowest number of generations and trim all other chains to match. 
  
  if(dir.exists(mcmc_chain_directory)==F){
    stop("Directory path is not a directory or there is a mistake in path")
  }
  
  if(missing(burnin)) {
    burnin=0.05
    print("Burnin not specified setting default")
  } 
  
  if(missing(MTBDorLSBDS)){
    stop("Please specify MTBDorLSBDS argument, if you are extracting LSBDS or MTBD")
  }
  
  mcmc_chains <- list.files(mcmc_chain_directory, full.names= "T") 
  mcmc_list <- list()
  n_generations <- character()
  count = 0 
  for(i in mcmc_chains){
    count = count + 1
    
    
    
    if(MTBDorLSBDS == "LSBDS"){
    chain_dw <- read.table(i,
                               header=TRUE)  %>% 
        as.data.frame() 
    burnin_c <- burnin*max(chain_dw[,1])
    chain <- chain_dw[-(1:burnin_c), , drop = FALSE] %>%
      as.mcmc()
    
    mcmc_list[[count]] <- chain
    chain_generation <- summary(chain)[["end"]]
    n_generations <- c(n_generations,chain_generation)
    }
    
    if(MTBDorLSBDS == "MTBD"){
    chain_dw <- read.csv(i,header=TRUE)  %>% 
        as.data.frame() 
    burnin_amount <- round((dim(chain_dw )[1]-1)*burnin, digits=0) +1
    chain  <- chain_dw[-(1:burnin_amount),(2:dim(chain_dw )[2]), drop = FALSE] %>%
      as.mcmc()
    mcmc_list[[count]] <- chain
    chain_generation <- summary(chain)[["end"]]
    n_generations <- c(n_generations,chain_generation)
    }
    
    
    
  }
  
  if(length(unique(n_generations))==1) {
    combined_chains <- mcmc.list(mcmc_list) %>% 
      as.matrix() %>%
      as.data.frame()
  } else{
    
    min_generations <- min(n_generations)
    mcmc_list <- lapply(mcmc_list,
                        function(x) as.mcmc(x[1:min_generations,]) )
    
    combined_chains <- mcmc.list(mcmc_list)%>% 
      as.matrix() %>%
      as.data.frame()
    
  }
  
  if(MTBDorLSBDS == "LSBDS"){
    lambds_order_pattern <- "avg_lambda.([0-9]+)"
    mu_order_pattern <- "avg_mu.([0-9]+)"
    net_div_pattern <- "avg_netDiv.\\1"
  }
  
  if(MTBDorLSBDS == "MTBD"){
    lambds_order_pattern <- "lambda_([0-9]+)"
    mu_order_pattern <- "mu_([0-9]+)"
    net_div_pattern <-  "netDiv_\\1"
  }
  
  
# Calculate net Div
  ## Subset out Lambda
  lamba_columns <- grep(pattern = "lambda",
                        x=colnames(combined_chains),
                        value=T)
  lamba_columns_order <- gsub(pattern = lambds_order_pattern,
                              replacement="\\1",
                              x = lamba_columns) %>% as.numeric()
  combined_chains_lambda <- combined_chains[, colnames(combined_chains) %in% lamba_columns]
  
  ## Subset out mu values
  mu_columns <- grep(pattern = "mu",
                     x=colnames(combined_chains),
                     value=T)
  mu_columns_order <- gsub(pattern = mu_order_pattern,
                           replacement="\\1",
                           x = mu_columns) %>% as.numeric()
  combined_chains_mu <- combined_chains[, colnames(combined_chains) %in% mu_columns]
  
  # Calculate net div
  combined_chains_netDiv <- combined_chains_lambda - combined_chains_mu
  colnames(combined_chains_netDiv) <- gsub(pattern = lambds_order_pattern,
                                           replacement = net_div_pattern,
                                           x = colnames(combined_chains_netDiv))
  combined_chains <- cbind(combined_chains_lambda,combined_chains_mu,combined_chains_netDiv)
  return(combined_chains)
  
}


claDS2_combine_burin <- function(claDS2Rdata_path,burnin=.25){
  
  if(dir.exists(claDS2Rdata_path)==F){
    stop("Directory path is not a directory or there is a mistake in path")
  }
  
  if(missing(burnin)) {
    burnin=0.25
    print("Burnin not specified setting default")
  } 
  
  # Load RPANDA output
  claDS2Rdata_path_files <- list.files(claDS2Rdata_path, full.names= "T") 
  if(length(claDS2Rdata_path_files) >1){
    stop("CladS2 directory has more than 1 file. It should only contain ClaDS2 data file")
  }
  
  load(claDS2Rdata_path_files)
  
  n <- length(CladsOutput$chains[[1]][[1]]) # number of iterations in each chain
  id <- length(CladsOutput$tree$edge[,2])  # number of edges in phylogeny
  
  # Initialize dataframe to store extracted values
  clads2_posterior <- data.frame("remove_column!" =
                                   seq(4:(length((length(1:floor(n/4)):n))*3)))
  clads2_posterior_netDiv <- data.frame("remove_column!" =
                                          seq(4:(length((length(1:floor(n/4)):n))*3)))
  
  # Extrac epsilon estimates
  epsilon_posterior = sapply(1:3,function(i){
    CladsOutput$chains[[i]][[3]][-(1:floor(n*burnin))]}) %>% c() %>% as.data.frame()
  
  # Iterate over lambda paramters for each edge
  for (o in 1:id) {
    lambda_edge_posterior = sapply(1:3,function(i){
      CladsOutput$chains[[i]][[4+o]][-(1:floor(n*burnin))]}) %>%
      # Takes burnin for a particular paramters, code from Odile?
      c() %>% 
      # concatinate posterior per chain into a single vector
      as.data.frame()
    # Make dataframe
    
    # Calculate mu per edge
    mu_edge_posterior = epsilon_posterior*lambda_edge_posterior
    
    # Rename Parametesr - they are in the same order  tree$edge[,2]
    colnames(mu_edge_posterior) <- paste0("mu_",
                                          CladsOutput$tree$edge[o,2])
    colnames(lambda_edge_posterior) <- paste0("lambda_",
                                              CladsOutput$tree$edge[o,2])
    
    # Combine
    clads2_posterior <- cbind(clads2_posterior,
                              lambda_edge_posterior,
                              mu_edge_posterior)
    
    net_div_edge_posterior <- lambda_edge_posterior-mu_edge_posterior
    
    colnames(net_div_edge_posterior)  <- paste0("netDiv_",
                                                CladsOutput$tree$edge[o,2])
    clads2_posterior_netDiv <- cbind(clads2_posterior_netDiv,
                                     net_div_edge_posterior)
  }
  
  
  # Remove first column 
  clads2_posterior$remove_column. <- NULL 
  clads2_posterior_netDiv$remove_column. <- NULL 
  
  # Write code to combine clads2_posterior and clads2_posterior_netDiv
  clads2_posterior_combined <- cbind(clads2_posterior,clads2_posterior_netDiv)
  return(clads2_posterior_combined)
}

# MTBD
# Th
.parseTreeString <- function(string) {
  text <- sub("[^(]*", "", string)
  # stats <- treeio:::read.stats_beast_internal( "", text )
  stats <- .read.stats_revbayes_internal( "", text )
  tree <- ape::read.tree(text = text)
  obj <- treeio:::BEAST("", text, stats, tree )
  return(obj)
}

add_pseudo_nodelabel <- function(phylo) {
  if(is.null(phylo$node.label)) {
    nnode <- phylo$Nnode
    phylo$node.label <- paste("X", 1:nnode, sep="")
    ## for (i in 1:nnode) {
    ##     treetext <- sub("\\)([:;])",
    ##                     paste0("\\)", nlab[i], "\\1"),
    ##                     treetext)
    ## }
  }
  ## if tip.label contains () which will broken node name extraction
  phylo$tip.label <- gsub("[\\(\\)]", "_", phylo$tip.label)
  
  treetext <- treeio::write.tree(phylo)
  return(treetext)
}

.read.stats_revbayes_internal <- function(beast, tree) {
  
  phylo <- ape::read.tree(text = tree) # read the tree
  tree2 <- add_pseudo_nodelabel(phylo) # add nodelabels (if there aren't already any)
  
  ## node name corresponding to stats
  nn <- unlist(strsplit(unlist(strsplit(tree2, split=",")), "\\)"))
  nn <- gsub("\\(*", "", nn)
  nn <- gsub("[:;].*", "", nn)
  nn <- gsub(" ", "", nn)
  nn <- gsub("'", "", nn)
  nn <- gsub('"', "", nn)
  
  # nn <- strsplit(tree2, split=",") %>% unlist %>%
  #   strsplit(., split="\\)") %>% unlist %>%
  #   gsub("\\(*", "", .) %>%
  #   gsub("[:;].*", "", .) %>%
  #   gsub(" ", "", .) %>%
  #   gsub("'", "", .) %>%
  #   gsub('"', "", .)
  # get the label of each node (internal and external) in the order
  # they appear in the newick string
  
  phylo <- ape::read.tree(text = tree2)
  root <- tidytree:::rootnode(phylo)
  nnode <- phylo$Nnode
  
  tree_label <- c(phylo$tip.label, phylo$node.label)
  ii <- match(nn, tree_label)
  
  if ( any(grepl("TRANSLATE", beast, ignore.case = TRUE)) ) {
    label2 <- c(phylo$tip.label, root:treeio:::getNodeNum(phylo))
  } else {
    label2 <- as.character(1:treeio:::getNodeNum(phylo))
  }
  node <- label2[match(nn, tree_label)]
  
  ## stats <- unlist(strsplit(tree, "\\["))[-1]
  ## stats <- sub(":.+$", "", stats
  
  ## BEAST1 edge stat fix
  tree <- gsub("\\]:\\[&(.+?\\])", ",\\1:", tree)
  tree <- gsub(":(\\[.+?\\])", "\\1:", tree)
  
  if (grepl("\\]:[0-9\\.eE+\\-]*\\[", tree) || grepl("\\]\\[", tree)) {
    ## MrBayes output
    # stats <- strsplit(tree, "\\]:[0-9\\.eE+\\-]*\\[") %>% unlist
    stats <- unlist(strsplit(tree, "\\]:[0-9\\.eE+\\-]*\\["))
    lstats <- lapply(stats, function(x) {
      unlist(strsplit(x, split="\\][,\\)]"))
    })
    
    for (i in seq_along(stats)) {
      n <- length(lstats[[i]])
      if (i == length(stats)) {
        stats[i] <- lstats[[i]][n]
      } else {
        stats[i] <- paste0(lstats[[i]][n],
                           sub("&", ",", lstats[[i+1]][1])
        )
      }
    }
    stats <- gsub("\\]\\[&", ",", stats)
  } else {
    ## BEAST output
    stats <- unlist(strsplit(tree, ":"))
  }
  
  names(stats) <- node
  
  stats <- stats[grep("\\[", stats)]
  stats <- sub("[^\\[]*\\[", "", stats)
  
  stats <- sub("^&", "", stats)
  stats <- sub("];*$", "", stats)
  stats <- gsub("\"", "", stats)
  
  #this is what is breaking readTrees for the OU output
  stats2 <- lapply(seq_along(stats), function(i) {
    
    x <- stats[[i]]
    y <- unlist(strsplit(x, ","))
    sidx <- grep("=\\{", y)
    eidx <- grep("\\}$", y)
    
    flag <- FALSE
    if (length(sidx) > 0) {
      flag <- TRUE
      # SETS <- lapply(seq_along(sidx), function(k) {
      #   p <- y[sidx[k]:eidx[k]]
      #   gsub(".*=\\{", "", p) %>% gsub("\\}$", "", .)
      # })
      SETS <- lapply(seq_along(sidx), function(k) {
        p <- y[sidx[k]:eidx[k]]
        gsub("\\}$", "",gsub(".*=\\{", "", p))
      })
      names(SETS) <- gsub("=.*", "", y[sidx])
      
      kk <- unlist(lapply(seq_along(sidx), function(k) {
        sidx[k]:eidx[k]
      }))
      y <- y[-kk]
    }
    
    if (length(y) == 0)
      return(SETS)
    
    name <- gsub("=.*", "", y)
    val <- gsub(".*=", "", y)
    val <- gsub("^\\{", "", val)
    val <- gsub("\\}$", "", val)
    
    # %>%
    #   gsub("^\\{", "", .) %>%
    #   gsub("\\}$", "", .)
    
    if (flag) {
      nn <- c(name, names(SETS))
    } else {
      nn <- name
    }
    
    res <- rep(NA, length(nn))
    names(res) <- nn
    
    for (i in seq_along(name)) {
      # res[i] <- if(treeio:::is.numeric(val[i])) as.numeric(val[i]) else val[i]
      res[i] <- if(is.numeric(val[i])) as.numeric(val[i]) else val[i]
    }
    if (flag) {
      j <- i
      for (i in seq_along(SETS)) {
        if(is.numeric(SETS[[i]])) {
          res[i+j] <- list(as.numeric(SETS[[i]]))
        } else {
          res[i+j] <- SETS[i]
        }
      }
    }
    
    return(res)
  })
  
  nn <- sort(unique(unlist(lapply(stats2, names))))
  
  # nn <- lapply(stats2, names) %>% unlist %>%
  #   unique %>% sort
  
  stats2 <- lapply(stats2, function(x) {
    y <- x[nn]
    names(y) <- nn
    y[vapply(y, is.null, logical(1))] <- NA
    y
  })
  
  stats3 <- do.call(rbind, stats2)
  stats3 <- tibble::as_tibble(stats3)
  
  ## no need to extract sd from prob+-sd
  ## as the sd is stored in prob_stddev
  ##
  ## "prob_stddev"   "prob(percent)" "prob+-sd"
  ##
  ##
  ##
  ## idx <- grep("\\+-", colnames(stats3))
  ## if (length(idx)) {
  ##     for (i in idx) {
  ##         stats3[,i] <- as.numeric(gsub("\\d+\\+-", "", stats3[,i]))
  ##     }
  ## }
  
  cn <- gsub("(\\d+)%", "0.\\1", colnames(stats3))
  cn <- gsub("\\(([^\\)]+)\\)", "_\\1", cn)
  ## cn <- gsub("\\+-", "_", cn)
  
  colnames(stats3) <- cn
  stats3$node <- names(stats)
  
  i <- vapply(stats3,
              function(x) max(vapply(x, length, numeric(1))),
              numeric(1))
  
  for (j in which(i==1)) {
    stats3[,j] <- unlist(stats3[,j])
  }
  stats3$node <- as.integer(stats3$node)
  return(stats3)
}

# fake newick string that mimics our data
#tree_string <- "(A[&mu=1,lambda=2]:1.0,(B[&mu=.5,lambda=1]:0.5,C[&mu=5,lambda=0.1]:0.5)[&mu=.75,lambda=1]:0.5)[&mu=.25,lambda=.1];"
#test <- .parseTreeString(tree_string)

# get data by node 
#test@data



MTBD_extract <- function(MTBD_rates_treeFile){
  
  if(file.exists(MTBD_rates_treeFile)==F){
    stop("File does not exisits or there is a mistake in path")
  }
  # Step 0: read Phylo trace

file <- MTBD_rates_treeFile

# Step 1: Identify the text line right before the line where the tree samples begin
con <- file(description=MTBD_rates_treeFile, "r")
line <- readLines(con, 1)
trees_start = 0 
while(!grepl("^tree", line)) {
  trees_start = trees_start + 1
  line <- readLines(con, 1)
}
close(con)

# Step 2: Fine how many lines their are
com <- paste("wc -l ", MTBD_rates_treeFile, " | awk '{ print $1 }'", sep="")
trees_end <- as.numeric(system(command=com, intern=TRUE))-1

# Step 3: Initalize dataframe to store values using first tree as guide
# 1. This codes skips the lines up to tree start and read the next line which is the first tree. 
# 2. Afterwards is extracts just the annotated newick
# 3. Reads it into custome script which extracts data 
tree <- sub(pattern = "tree.* \\(",
            replacement = "\\(", 
            x = read_lines(file, skip=trees_start, n_max = 1)) %>%
  .parseTreeString()
# 4. The following generates an 2*n by m+1 matrix where n is the number of nodes/ tips in the 
# phylogeny (multiplied by two for lambda an mu at each node) and m is the number of MCMC 
# in the matrix (plus 1 which is the node ordering)
BEAST <- c(tree@data$node,tree@data$node) %>% as.data.frame()
rownames(BEAST) <- c(sub(pattern = "^", replacement = "lambda_", x = tree@data$node),
                     sub(pattern = "^", replacement = "mu_", x = tree@data$node))
MSBD_Posterior_node <- t(BEAST) %>% data.frame()
rownames(MSBD_Posterior_node) <- "node"


# Step 5: Using code from step 4 loop over all the phylogenies and get the data! 
count = 1 # first row stores node numbers
percentage <- round(quantile(trees_start:(trees_end-1),
                             +          c(.1,.2,.3,.4,.5,.6,.7,.8,.9,1)))
for(i in trees_start:(trees_end-1)) {
  # Reads 
  tree<- sub(pattern = "tree.* \\(",
             replacement = "\\(", 
             x = read_lines(file, skip=i, n_max = 1)) %>%
    .parseTreeString()
  
  metadata_temp <- c(tree@data$lambda,tree@data$mu) %>% as.data.frame()
  rownames(metadata_temp) <- c(sub(pattern = "^", replacement = "lambda_", x = tree@data$node),
                               sub(pattern = "^", replacement = "mu_", x = tree@data$node))
  metadata_temp <- t(metadata_temp) %>% data.frame()
  
  MSBD_Posterior_node <- rbind(MSBD_Posterior_node, metadata_temp)
  count = count + 1
  rownames(MSBD_Posterior_node)[count] <- count
  
  if(i %in% percentage ){
    quant <- names(percentage)[percentage %in% i]
    print(paste0(quant," Done. Working on tree ",i, " of ", (trees_end-1)))
  }
}

MSBD_Posterior_node[] <- lapply(MSBD_Posterior_node, function(x) as.numeric(as.character(x)))
return(MSBD_Posterior_node)

}

BAMM_combine <- function(BAMM_rates_path){
  
  if(dir.exists(BAMM_rates_path)==F){
    stop("File path is not a File or there is a mistake in path")
  }
  
bamm_posterior <- read_csv(paste0(BAMM_rates_path,"/BAMMrates.csv")) %>%
  as.data.frame()

bamm_mu <- bamm_posterior[ ,grepl( "mu_" , colnames( bamm_posterior ) )   ]
mu_columns_order <- gsub(pattern = "mu_branch_matrix.([0-9]+)",
                         replacement="\\1",
                         x = colnames(bamm_mu)) %>% as.numeric()


bamm_lambda <- bamm_posterior[ ,grepl( "lambda_" , colnames( bamm_posterior ) )  ]
lamba_columns_order <- gsub(pattern = "lambda_branch_matrix.([0-9]+)",
                            replacement="\\1",
                            x = colnames(bamm_lambda)) %>% as.numeric()

if(all(lamba_columns_order==mu_columns_order)) {
  bamm_netDiv <- bamm_lambda - bamm_mu
  
  colnames(bamm_netDiv) <- gsub(pattern ="lambda_branch_matrix.([0-9]+)",
                                replacement="netDiv_branch_matrix.\\1",
                                x=colnames(bamm_netDiv))
  }

  bamm_combined <- cbind(bamm_lambda,bamm_mu,bamm_netDiv)
  return(bamm_combined)


}

#### Summarize Posterior

######### Functions calculates HPD and formats output. Used in summarizePosterior() 
calcHPD <- function(mcmc){
  coda_HPDinterval <- as.list(mcmc) %>%
    lapply(X=.,FUN=as.mcmc) %>%
    lapply(X=.,FUN=coda::HPDinterval,prob= .95) %>%
    as.data.frame()  
  coda_HPDinterval_lower <- coda_HPDinterval[,grep(x = colnames(coda_HPDinterval),
                                                   pattern = "lower")] %>%
    t() %>% 
    as.data.frame()
  rownames(coda_HPDinterval_lower) <- gsub(x = rownames(coda_HPDinterval_lower),
                                           pattern = "(.*).lower",
                                           replacement = "\\1")
  
  coda_HPDinterval_upper <- coda_HPDinterval[,grep(x = colnames(coda_HPDinterval),
                                                   pattern = "upper")] %>%
    t()  %>% 
    as.data.frame()
  rownames(coda_HPDinterval_upper) <- gsub(x = rownames(coda_HPDinterval_upper),
                                           pattern = "(.*).upper",
                                           replacement = "\\1")
  HPD_intervals <- cbind(coda_HPDinterval_upper,coda_HPDinterval_lower)
  return(HPD_intervals)
  
}

# test
LSBDS_posterior <- "/Users/jesusmartinez-gomez/Desktop/A2017/output/combined_posteriors/LSBDS_combinedPosterior.csv"
phylo_path<- "/Users/jesusmartinez-gomez/Desktop/A2017/phylo_A2017.tr.pr_.txt"

LSBDS_summarizePosterior <- function(LSBDS_posterior_path = NULL,
                                     Phylogeny_path = NULL){
  
  if(file.exists(LSBDS_posterior_path)==F){
    stop("LSBDS posterior file not found")
  }
  
  if(file.exists(Phylogeny_path)==F){
    stop("Phylogeny path not found.")
  }
  
  phylo <- read.tree(phylo_path)
  nodeMatchKey <- matchNodes(phylo) 
  root_edge <- setdiff(1:max(phylo$edge[,2]),phylo$edge[,2])
                               
  lsbds_posterior <- read.table(LSBDS_posterior_path,header=T,sep = ",") %>% 
    as.data.frame()
  lsbdsBurnin<-na.omit(lsbds_posterior) # Remove MCMC samples with NA - see notes L2012 had an issues
  
  # Extract Speciation and index
  rev_lambda <- lsbdsBurnin[ , grepl( "avg_lambda" , colnames( lsbdsBurnin ) ) ] # extract lambda values
  rev_summary_lambda <- colMeans(rev_lambda) %>% # calculate average
    as.data.frame()
  colnames(rev_summary_lambda)[1] <- "LSBDS_avg_lambda"
  ## Calculate HPD
  tmp <- calcHPD(rev_lambda)
  colnames(tmp) <- c("LSBDS_HPD0.05_lambda","LSBDS_HPD0.95_lambda")
  rev_summary_lambda <- cbind(rev_summary_lambda, tmp)
  ## Calculate Quantile 
  rev_summary_lambda$LSBDS_Quantile0.05_lambda <- mapply(FUN=quantile,rev_lambda,c(0.05)) 
  rev_summary_lambda$LSBDS_Quantile0.95_lambda <- mapply(FUN=quantile,rev_lambda,c(0.95))
  ## Calculate MAP
  rev_summary_lambda$LSBDS_MAP_lambda <- apply(rev_lambda,2,getMAP) # calculate  MAP
  ## Calculate posterior median
  rev_summary_lambda$LSBDS_median_lambda <- apply(rev_lambda,2,median) 
  ## Index Nodes
  rev_summary_lambda$RevBranches <- grep("([0-9])",rownames(rev_summary_lambda)) 
  rev_summary_lambda$Index <- match(rev_summary_lambda$RevBranches, nodeMatchKey$Rev) 
  ## Removes the root node
  rev_summary_lambda <- rev_summary_lambda[-match(root_edge,rev_summary_lambda$Index),] 
  
  # Extract Extinction and index
  rev_mu <- lsbdsBurnin[ , grepl( "avg_mu" , colnames( lsbdsBurnin ) ) ]
  ## Calculate Posterior Mean
  rev_summary_mu <- colMeans(rev_mu) %>% 
    as.data.frame()
  colnames(rev_summary_mu)[1] <- "LSBDS_avg_mu"
  ## Calculate HPD
  tmp <- calcHPD(rev_mu)
  colnames(tmp) <- c("LSBDS_HPD0.05_mu","LSBDS_HPD0.95_mu")
  rev_summary_mu <- cbind(rev_summary_mu, tmp)
  ## Calculate Quantile
  rev_summary_mu$LSBDS_Quantile0.05_mu <- mapply(FUN=quantile,rev_mu,c(0.05))
  rev_summary_mu$LSBDS_Quantile0.95_mu <- mapply(FUN=quantile,rev_mu,c(0.95))
  ## Calculate MAP
  rev_summary_mu$LSBDS_MAP_mu <- apply(rev_mu,2,getMAP) 
  ## Calculate Posterior Median
  rev_summary_mu$LSBDS_median_mu<- apply(rev_mu,2,median) 
  ## Index Nodes
  rev_summary_mu$RevBranches <- grep("([0-9])",rownames(rev_summary_mu)) 
  rev_summary_mu$Index <- match(rev_summary_mu$RevBranches, nodeMatchKey$Rev)
  ## Removes the root node
  rev_summary_mu <- rev_summary_mu[-match(root_edge,rev_summary_mu$Index),] 
  
  # Extract Net Diversification and index
  rev_netDiv <- lsbds_posterior[ , grepl( "avg_netDiv" , colnames( lsbds_posterior ) ) ] # extract net Div values
  ## Calculate Posterior Mean
  rev_summary_netDiv <- colMeans(rev_netDiv) %>% 
    as.data.frame()
  colnames(rev_summary_netDiv)[1] <- "LSBDS_avg_div"
  # Calculate HPD
  tmp <- calcHPD(rev_netDiv)
  colnames(tmp) <- c("LSBDS_HPD0.05_div","LSBDS_HPD0.95_div")
  rev_summary_netDiv <- cbind(rev_summary_netDiv, tmp)
  ## Calculate Quantile
  rev_summary_netDiv$LSBDS_Quantile0.05_div <- mapply(FUN=quantile,rev_netDiv,c(0.05)) # calculate 0.05 quantile
  rev_summary_netDiv$LSBDS_Quantile0.95_div <- mapply(FUN=quantile,rev_netDiv,c(0.95)) # calculate 0.95 quantile
  ## Calculate MAP
  rev_summary_netDiv$LSBDS_MAP_div <- apply(rev_netDiv,2,getMAP) 
  ## Calculate Posterior Median
  rev_summary_netDiv$LSBDS_median_div <- apply(rev_netDiv,2,median) 
  ## Index Nodes
  rev_summary_netDiv$RevBranches <- grep("([0-9])",rownames(rev_summary_netDiv)) 
  rev_summary_netDiv$Index <- match(rev_summary_netDiv$RevBranches, nodeMatchKey$Rev) 
  ## Remove the root node
  rev_summary_netDiv <- rev_summary_netDiv[-match(root_edge,rev_summary_netDiv$Index),] # Remove the root node 
  
  LSBDS_summaryStats <- merge(rev_summary_lambda,rev_summary_mu, by="Index") %>%
    merge(rev_summary_netDiv,by="Index")
  LSBDS_summaryStats$RevBranches.x <- NULL
  LSBDS_summaryStats$RevBranches.y <- NULL
  return(LSBDS_summaryStats)
  
}

# test files
phylo_path<- "/Users/jesusmartinez-gomez/Desktop/A2017/phylo_A2017.tr.pr_.txt"
bamm_posterior_path <- "/Users/jesusmartinez-gomez/Desktop/A2017/output/combined_posteriors/BAMM_combinedPosterior.csv"

BAMM_summarizePosterior <- function(bamm_posterior_path = NULL,
                                     Phylogeny_path = NULL){
  
  if(file.exists(bamm_posterior_path)==F){
    stop("BAMM posterior file not found")
  }
  
  if(file.exists(Phylogeny_path)==F){
    stop("Phylogeny path not found.")
  }
  # Phylogeny
  phylo <- read.tree(phylo_path)
  nodeMatchKey <- matchNodes(phylo) 
  root_edge <- setdiff(1:max(phylo$edge[,2]),phylo$edge[,2])
  
  
  bamm_posterior <- read_csv(bamm_posterior_path) %>%
    t() %>% 
    as.data.frame()
  
  # Extract speciation and Index
  bamm_lambda <- bamm_posterior[ grepl( "lambda_" , rownames( bamm_posterior ) ) ,  ]
  ## Calculate Posterior Mean
  bamm_summary_lambda <- colMeans(bamm_lambda) %>% 
    as.data.frame()
  colnames(bamm_summary_lambda)[1] <- "BAMM_avg_lambda"
  ## Calculate HPD
  tmp <- calcHPD(bamm_lambda)
  colnames(tmp) <- c("BAMM_HPD0.05_lambda","BAMM_HPD0.95_lambda")
  bamm_summary_lambda <- cbind(bamm_summary_lambda, tmp)
  ## Calculate quantile
  bamm_summary_lambda$BAMM_Quantile0.05_lambda <- mapply(FUN=quantile,bamm_lambda,c(0.05)) # calculate 0.05 quantile
  bamm_summary_lambda$BAMM_Quantile0.95_lambda <- mapply(FUN=quantile,bamm_lambda,c(0.95)) # calculate 0.95 quantile
  ## Calculate MAP
  bamm_summary_lambda$BAMM_MAP_lambda <- apply(bamm_lambda,2,getMAP) 
  ## Calculate Posterior Median
  bamm_summary_lambda$BAMM_median_lambda<- apply(bamm_lambda,2,median) 
  ## Index Nodes
  bamm_summary_lambda$BammBranches <- grep("([0-9])",rownames(bamm_summary_lambda)) 
  bamm_summary_lambda$Index <- phylo$edge[,2]
  
  # Extract Extinction and Index
  bamm_mu <- bamm_posterior[ grepl( "mu_" , rownames( bamm_posterior ) ) ,  ]
  ## Calculate Posterior Mean
  bamm_summary_mu <- colMeans(bamm_mu) %>% 
    as.data.frame()
  colnames(bamm_summary_mu)[1] <- "BAMM_avg_mu"
  ## Calcualte HPD
  tmp <- calcHPD(bamm_mu)
  colnames(tmp) <- c("BAMM_HPD0.05_mu","BAMM_HPD0.95_mu")
  bamm_summary_mu <- cbind(bamm_summary_mu, tmp)
  ## Calculate quantile
  bamm_summary_mu$BAMM_Quantile0.05_mu <- mapply(FUN=quantile,bamm_mu,c(0.05)) # calculate 0.05 quantile
  bamm_summary_mu$BAMM_Quantile0.95_mu <- mapply(FUN=quantile,bamm_mu,c(0.95)) # calculate 0.95 quantile
  ## Calculate MAP
  bamm_summary_mu$BAMM_MAP_mu <- apply(bamm_mu,2,getMAP)
  ## Calculate Posterior Median
  bamm_summary_mu$BAMM_median_mu<- apply(bamm_mu,2,median) 
  ## Index nodes
  bamm_summary_mu$BammBranches <- grep("([0-9])", rownames(bamm_summary_mu)) 
  bamm_summary_mu$Index <- phylo$edge[,2]
  
  
  # Extract Net diverisification and Index
  bamm_netDiv <- bamm_posterior[ grepl( "netDiv_" , rownames( bamm_posterior ) ) ,  ]
  ## Calculate Posterior Mean
  bamm_summary_netDiv <- colMeans(bamm_netDiv) %>% 
    as.data.frame()
  colnames(bamm_summary_netDiv)[1] <- "BAMM_avg_div"
  ## Calcualte HPD
  tmp <- calcHPD(bamm_netDiv)
  colnames(tmp) <- c("BAMM_HPD0.05_div","BAMM_HPD0.95_div")
  bamm_summary_netDiv <- cbind(bamm_summary_netDiv, tmp)
  #Calculate Quantile
  bamm_summary_netDiv$BAMM_Quantile0.05_div <- mapply(FUN=quantile,bamm_netDiv,c(0.05)) # calculate 0.05 quantile
  bamm_summary_netDiv$BAMM_Quantile0.95_div <- mapply(FUN=quantile,bamm_netDiv,c(0.95)) # calculate 0.95 quantile
  ## Calculate MAP
  bamm_summary_netDiv$BAMM_MAP_div <- apply(bamm_netDiv,2,getMAP) 
  ## Calculate Posterior Median
  bamm_summary_netDiv$BAMM_median_div <- apply(bamm_netDiv,2,median) 
  ## Index nodes
  bamm_summary_netDiv$BammBranches <- grep("([0-9])", rownames(bamm_summary_netDiv)) 
  bamm_summary_netDiv$Index <- phylo$edge[,2]
  
  BAMM_summaryStats <- merge(bamm_summary_lambda,bamm_summary_mu, by="Index") %>%
    merge(bamm_summary_netDiv,by="Index")
  BAMM_summaryStats$BammBranches.x <- NULL
  BAMM_summaryStats$BammBranches.y <- NULL
  return(BAMM_summaryStats)
  
}


MTBD_posterior_path <- "/Users/jesusmartinez-gomez/Desktop/A2017/output/combined_posteriors/MTBD_combinedPosterior.csv"
phylo_path<- "/Users/jesusmartinez-gomez/Desktop/A2017/phylo_A2017.tr.pr_.txt"
MTBD_summarizePosterior <- function(MTBD_posterior_path = NULL,
                                    Phylogeny_path = NULL){
  
  if(file.exists(MTBD_posterior_path)==F){
    stop("MTBD posterior file not found")
  }
  
  if(file.exists(Phylogeny_path)==F){
    stop("Phylogeny path not found.")
  }
  
  # Phylogeny
  phylo <- read.tree(phylo_path)
  nodeMatchKey <- matchNodes(phylo) 
  root_edge <- setdiff(1:max(phylo$edge[,2]),phylo$edge[,2])
  
  # Read in Posterior
  msbd_posterior <- read_csv(MTBD_posterior_path) %>%
    as.data.frame()
  
  # Speciation
  msbd_lambda <- msbd_posterior[ , grepl( "lambda" , colnames( msbd_posterior ) ) ] # extract lambda values
  ## Calculate Posterior Mean
  msbd_summary_lambda <- colMeans(msbd_lambda) %>% 
    as.data.frame()
  colnames(msbd_summary_lambda)[1] <- "MTBD_avg_lambda"
  # Calculate HPD
  tmp <- calcHPD(msbd_lambda)
  colnames(tmp) <- c("MTBD_HPD0.05_lambda","MTBD_HPD0.95_lambda")
  msbd_summary_lambda <- cbind(msbd_summary_lambda, tmp)
  ## Calculate quantile
  msbd_summary_lambda$MTBD_Quantile0.05_lambda <- mapply(FUN=quantile,msbd_lambda,c(0.05)) # calculate 0.05 quantile
  msbd_summary_lambda$MTBD_Quantile0.95_lambda <- mapply(FUN=quantile,msbd_lambda,c(0.95)) # calculate 0.95 quantile
  ## Calculate MAP Median
  msbd_summary_lambda$MTBD_MAP_lambda <- apply(msbd_lambda,2,getMAP) 
  ## Calculate Posterior Median 
  msbd_summary_lambda$MTBD_median_lambda<- apply(msbd_lambda,2,median) 
  ## Index node
  msbd_summary_lambda$Index <- gsub(pattern = "lambda_([0-9])",
                                    x = rownames(msbd_summary_lambda),
                                    replacement = "\\1") %>% as.numeric()
  ## Remove Root
  msbd_summary_lambda <- msbd_summary_lambda[-match(root_edge,msbd_summary_lambda$Index),] 
  
  # Extinction
  msbd_mu<- msbd_posterior[ , grepl( "mu" , colnames( msbd_posterior ) ) ] # extract lambda values
  ## Calculate posterior mean
  msbd_summary_mu <- colMeans(msbd_mu) %>% 
    as.data.frame()
  colnames(msbd_summary_mu)[1] <- "MTBD_avg_mu"
  ## Calculate HPD
  tmp <- calcHPD(msbd_mu)
  colnames(tmp) <- c("MTBD_HPD0.05_mu","MTBD_HPD0.95_mu")
  msbd_summary_mu <- cbind(msbd_summary_mu, tmp)
  ## Calculate quantile
  msbd_summary_mu$MTBD_Quantile0.05_mu<- mapply(FUN=quantile,msbd_mu,c(0.05)) # calculate 0.05 quantile
  msbd_summary_mu$MTBD_Quantile0.95_mu <- mapply(FUN=quantile,msbd_mu,c(0.95)) # calculate 0.95 quantile
  ## Calculate MAP
  msbd_summary_mu$MTBD_MAP_mu <- apply(msbd_mu,2,getMAP) 
  ## Calculate Posterior Median
  msbd_summary_mu$MTBD_median_mu<- apply(msbd_mu,2,median) 
  ## Index node
  msbd_summary_mu$Index <- gsub(pattern = "mu_([0-9])",
                                x = rownames(msbd_summary_mu),
                                replacement = "\\1") %>% as.numeric()
  ## Remove the root
  msbd_summary_mu <- msbd_summary_mu[-match(root_edge,msbd_summary_mu$Index),] 
  
  # Net div
  msbd_netDiv <- msbd_posterior[ , grepl( "netDiv" , colnames( msbd_posterior ) ) ] # extract lambda values
  ## Calculate Posterior Mean
  msbd_summary_netDiv <- colMeans(msbd_netDiv) %>% # calculate average
    as.data.frame()
  colnames(msbd_summary_netDiv)[1] <- "MTBD_avg_div"
  ## Calculate HPD
  tmp <- calcHPD(msbd_netDiv)
  colnames(tmp) <- c("MTBD_HPD0.05_div","MTBD_HPD0.95_div")
  msbd_summary_netDiv <- cbind(msbd_summary_netDiv, tmp)
  ## Calculate Quantile
  msbd_summary_netDiv$MTBD_Quantile0.05_div <- mapply(FUN=quantile,msbd_netDiv,c(0.05)) # calculate 0.05 quantile
  msbd_summary_netDiv$MTBD_Quantile0.95_div <- mapply(FUN=quantile,msbd_netDiv,c(0.95)) # calculate 0.95 quantile
  ## Calculate  MAP
  msbd_summary_netDiv$MTBD_MAP_div <- apply(msbd_netDiv,2,getMAP) 
  ## Calculate  Posterior Median
  msbd_summary_netDiv$MTBD_median_div <- apply(msbd_netDiv,2,median) 
  ## Index Nodes
  msbd_summary_netDiv$Index <- gsub(pattern = "netDiv_([0-9])",
                                    x = rownames(msbd_summary_netDiv),
                                    replacement = "\\1") %>% as.numeric()
  ## Remove the nodes
  msbd_summary_netDiv <- msbd_summary_netDiv[-match(root_edge,msbd_summary_netDiv$Index),] 
  
  # Combine Everything
  MTBD_summaryStats <- merge(msbd_summary_lambda,msbd_summary_mu, by="Index") %>%
    merge(msbd_summary_netDiv,by="Index")
  MTBD_summaryStats$BammBranches.x <- NULL
  MTBD_summaryStats$BammBranches.y <- NULL
  return(MTBD_summaryStats)
  
}

clads2_summarizePosterior <- function(clads2_posterior_path = NULL,
                                    Phylogeny_path = NULL){
  
  if(file.exists(clads2_posterior_path)==F){
    stop("Clads2 posterior file not found")
  }
  
  if(file.exists(Phylogeny_path)==F){
    stop("Phylogeny path not found.")
  }
  # Phylogeny
  phylo <- read.tree(phylo_path)
  nodeMatchKey <- matchNodes(phylo) 
  root_edge <- setdiff(1:max(phylo$edge[,2]),phylo$edge[,2])
  
  clads2_posterior <- read_csv(clads2_posterior_path) %>%
    as.data.frame()
  
  # Speciation and Index
  clads2_lambda <- clads2_posterior[ , grepl( "lambda" , colnames( clads2_posterior ) ) ] # extract lambda values
  ## Calculate Posterior Mean
  clads2_summary_lambda <- colMeans(clads2_lambda) %>%
    as.data.frame()
  colnames(clads2_summary_lambda)[1] <- "ClaDS2_avg_lambda"
  # Calculate HPD
  tmp <- calcHPD(clads2_lambda)
  colnames(tmp) <- c("ClaDS2_HPD0.05_lambda","ClaDS2_HPD0.95_lambda")
  clads2_summary_lambda <- cbind(clads2_summary_lambda, tmp)
  ## Calculate Quantile
  clads2_summary_lambda$ClaDS2_Quantile0.05_lambda <- mapply(FUN=quantile,clads2_lambda,c(0.05)) # calculate 0.05 quantile
  clads2_summary_lambda$ClaDS2_Quantile0.95_lambda <- mapply(FUN=quantile,clads2_lambda,c(0.95)) # calculate 0.95 quantile
  ## Calculate MAP
  clads2_summary_lambda$ClaDS2_MAP_lambda <- apply(clads2_lambda,2,getMAP) 
  ## Calculate Posterior Median
  clads2_summary_lambda$ClaDS2_median_lambda <- apply(clads2_lambda,2,median) 
  ## Index nodes
  clads2_summary_lambda$Index <- gsub(pattern = "lambda_([0-9])",
                                      x = rownames(clads2_summary_lambda),
                                      replacement = "\\1") %>% as.numeric()
  
  # Extinction and Index
  clads2_mu<- clads2_posterior[ , grepl( "mu" , colnames( clads2_posterior ) ) ] # extract lambda values
  ## Calculate Posterior Mean
  clads2_summary_mu <- colMeans(clads2_mu) %>% 
    as.data.frame()
  colnames(clads2_summary_mu)[1] <- "ClaDS2_avg_mu"
  # Calculate HPD
  tmp <- calcHPD(clads2_mu)
  colnames(tmp) <- c("ClaDS2_HPD0.05_mu","ClaDS2_HPD0.95_mu")
  clads2_summary_mu <- cbind(clads2_summary_mu, tmp)
  ## Calculate Quantile
  clads2_summary_mu$ClaDS2_Quantile0.05_mu <- mapply(FUN=quantile,clads2_mu,c(0.05)) # calculate 0.05 quantile
  clads2_summary_mu$ClaDS2_Quantile0.95_mu <- mapply(FUN=quantile,clads2_mu,c(0.95)) # calculate 0.95 quantile
  ## Calculate MAP
  clads2_summary_mu$ClaDS2_MAP_mu <- apply(clads2_mu,2,getMAP) 
  ## Calculate Posterior Median
  clads2_summary_mu$ClaDS2_median_mu <- apply(clads2_mu,2,median)
  ## Index
  clads2_summary_mu$Index <- gsub(pattern = "mu_([0-9])",
                                  x = rownames(clads2_summary_mu),
                                  replacement = "\\1") %>% as.numeric()
  
  ## Extract Net Diversification
  clads2_netDiv <- clads2_posterior[ , grepl( "netDiv" , colnames( clads2_posterior ) ) ] # extract lambda values
  ## Calculate Posterior Mean
  clads2_summary_netDiv <- colMeans(clads2_netDiv) %>% 
    as.data.frame()
  colnames(clads2_summary_netDiv)[1] <- "ClaDS2_avg_div"
  # #Calculate HPD
  tmp <- calcHPD(clads2_netDiv)
  colnames(tmp) <- c("ClaDS2_HPD0.05_div","ClaDS2_HPD0.95_div")
  clads2_summary_netDiv <- cbind(clads2_summary_netDiv, tmp)
  ## Calculate Quantile
  clads2_summary_netDiv$ClaDS2_Quantile0.05_div <- mapply(FUN=quantile,clads2_netDiv,c(0.05)) # calculate 0.05 quantile
  clads2_summary_netDiv$ClaDS2_Quantile0.95_div <- mapply(FUN=quantile,clads2_netDiv,c(0.95)) # calculate 0.95 quantile
  ## Calculate MAP
  clads2_summary_netDiv$ClaDS2_MAP_div <- apply(clads2_netDiv,2,getMAP) 
  ## Calculate Posterior Median
  clads2_summary_netDiv$ClaDS2_median_div <- apply(clads2_netDiv,2,median)
  ## Index and Node
  clads2_summary_netDiv$Index <- gsub(pattern = "netDiv_([0-9])",
                                      x = rownames(clads2_summary_netDiv),
                                      replacement = "\\1") %>% as.numeric()
  
  clads2_summaryStats <- merge(clads2_summary_lambda,clads2_summary_mu, by="Index") %>%
    merge(clads2_summary_netDiv,by="Index")
  return(clads2_summaryStats)
  
}


Pesto_summarizePosterior <- function(pesto_path,
                                     Phylogeny_path){
  if(file.exists(pesto_path)==F){
    stop("PESTO posterior file not found")
  }
  
  if(file.exists(Phylogeny_path)==F){
    stop("Phylogeny path not found.")
  }
  # Phylogeny
  phylo <- read.tree(phylo_path)
  nodeMatchKey <- matchNodes(phylo) 
  root_edge <- setdiff(1:max(phylo$edge[,2]),phylo$edge[,2])
  
  # Read in PESTO
  PESTO <- read.table(pesto_path,
                          sep = "\t",
                          header = T) %>% 
    t() %>%
    as.data.frame()
  
  # Speciation      
  PESTO$V2 <- rownames(PESTO)
  PESTO_lambda <- subset(PESTO, grepl(pattern="avg_lambda", x=PESTO$V2))
  PESTO_lambda$RevBranches <- gsub(pattern = "avg_lambda.([0-9]+).",
                                       replacement = "\\1",
                                       x = PESTO_lambda$V2) %>%
    as.numeric()
  ## Index nodes
  PESTO_lambda$Index <- match(PESTO_lambda$RevBranches, nodeMatchKey$Rev) 
  ## Remove the Root Node
  PESTO_lambda <- PESTO_lambda[-match(root_edge,PESTO_lambda$Index),] 
  PESTO_lambda$PESTO_lambda <- PESTO_lambda$V1
  PESTO_lambda$V2 <- NULL
  PESTO_lambda$V1<- NULL
  
  # Extinction
  PESTO$V2 <- rownames(PESTO)
  PESTO_mu <- subset(PESTO, grepl(pattern="avg_mu", x=PESTO$V2))
  PESTO_mu$RevBranches <- gsub(pattern = "avg_mu.([0-9]+).",
                                   replacement = "\\1",
                                   x = PESTO_mu$V2) %>%
    as.numeric()
  
  PESTO_mu$Index <- match(PESTO_mu$RevBranches, nodeMatchKey$Rev) # Index Rev estimates 
  PESTO_mu <- PESTO_mu[-match(root_edge,PESTO_mu$Index),] # Remove the root node
  PESTO_mu$PESTO_mu <- PESTO_mu$V1
  PESTO_mu$V2 <- NULL
  PESTO_mu$V1 <- NULL
  
  # Calculate net div and combine
  PESTO_combined <- merge(PESTO_lambda,PESTO_mu, by = "Index")
  PESTO_combined$PESTO_div <- PESTO_combined$PESTO_lambda - PESTO_combined$PESTO_mu
  return(PESTO_combined)
  
}

# Specify Some Function for Plotting 
plot_rate_phylo <- function(data) {

  ggtree_plot <- ggtree(data, aes(color=value),ladderize = T, continuous = TRUE) + 
    scale_color_gradientn(colours = c("#2c7bb6","#abd9e9","#ffffbf","#fdae61","#d7191c")) +
    # c("blue","white","red"),
    # Some error here - attempted to have each legend have same number of breaks. I think?
    # breaks=seq(min(node_data),max(node_data),
    #            (max(node_data)-min(node_data))/5)) +
    
    theme(legend.position = "right",
          legend.text = element_text(size=3.5),
          legend.title = element_text (size=7.25)) +
    labs(color="rate") 
  return(ggtree_plot)
  
}

plot_rate_phylo_scaled <- function(data,
                                   gradient_max,
                                   gradient_min) {
  ggtree_plot <- ggtree(data, aes(color=value),ladderize = T, continuous = TRUE) + 
    scale_color_gradientn(colours = c("#2c7bb6","#abd9e9","#ffffbf","#fdae61","#d7191c"), #c("blue","red"),#c("blue","white","red"),
                          limits = c(gradient_min, gradient_max)) +
    # scale_color_viridis(option = "D",
    #                     limits = c(gradient_min, gradient_max)) +
    theme(legend.position = "right",
          legend.text = element_text(size=3.5),
          legend.title = element_text (size=7.25)) +
    labs(color="rate") 
  labs(color="rate") 
  return(ggtree_plot)
}

comp_branch_plots <- function(summaryStatisticsData,
                              parameter,
                              Branch_index, 
                              phylo_path,
                              print_scale){
  # summaryStatisticsData - data
  # paramter - "SummaryStatistics is not correctly specified. Shoulds be either, \"avg\" for posterior mean, \"median\" for posterior median or \"MAP\" for maximum a posteriori"
  # index - a vector or that shows phylogeny index 
  
  if(missing(summaryStatisticsData)) {
    stop("SummaryStatistics path not found")
  } 
  
  if(file.exists(phylo_path)==F){
    stop("Phylopgeny not found, please check path.")
  }
  
  if(parameter %in% c("avg","median","map") == F) {
    stop("SummaryStatistics is not correctly specified. Shoulds be either, \"avg\" for posterior mean, \"median\" for posterior median or \"MAP\" for maximum a posteriori")
  } 
  
  if(missing(print_scale)) {
    print_scale = F
    stop("scaled argument not specified. Default is FALSE. If true phylogenies will share a legend.")
  } 
  
  combined <- summaryStatisticsData
  par <- parameter

  
  phylo_tib <- read.tree(phylo_path) %>%
    as.tibble(phylo)
  
  # 1. Check which methods are present in summary statistics
  div_methods <- data.frame(method = c("BAMM","LSBDS","PESTO","MTBD","ClaDS2"))
  
  PESTO_present <- grepl(x=colnames(combined),
                         pattern = "PESTO_lambda") %>%
    any()
  MTBD_present <- grepl(x=colnames(combined),
                        pattern = "MTBD_median_lambda")%>%
    any()
  LSBDS_present <- grepl(x=colnames(combined),
                         pattern = "LSBDS_median_lambda")%>%
    any()
  BAMM_present <- grepl(x=colnames(combined),
                        pattern = "BAMM_median_lambda")%>%
    any()
  ClaDS2_present <- grepl(x=colnames(combined),
                          pattern = "ClaDS2_median_lambda")%>%
    any()
  
  methods_present <- c(BAMM_present,LSBDS_present,PESTO_present,MTBD_present,ClaDS2_present)
  div_methods <- data.frame(method = div_methods$method[methods_present])
  ncolumns <- dim(div_methods)[1]
  
  # Generate labels
  LETTERS_labels <- LETTERS[1:(ncolumns*3)]
  first_row <- paste(LETTERS_labels[1:ncolumns],div_methods$method)
  plot_labels <- data.frame(row1 = first_row,
                            row2 = LETTERS_labels[(ncolumns+1):(ncolumns*2)],
                            row3 = LETTERS_labels[((ncolumns*2)+1):(ncolumns*3)])
  
  rate_labels<- data.frame(div = c("lambda","mu","div"),
                           symbol = c("Speciation","Extinction","Net Diversification"))
  
  
  count_par = 0 
  row_plots <- vector('list', ncolumns) # store row plots
  row_plots_scaled <- vector('list', ncolumns)
  for(b in rate_labels$div) {
    count_par = count_par + 1
    
    # 1. For scaled phylogeny we need to get min and max
    if("PESTO" %in% div_methods$method){
      
      others <- div_methods$method[div_methods$method != "PESTO"]
      combined_pesto <- combined[,paste0("PESTO_",b)] %>% 
        as.data.frame()
      colnames(combined_pesto) <- paste0("PESTO_",b)
      combined_others <- combined[,paste0(others,"_",par,"_",b)]
      # combined_melt <- cbind(combined_pesto,combined_others) %>%
      #   melt()
      combined_scale <- cbind(combined_pesto,combined_others) %>%
        scale() %>% 
        as.data.frame() 
      
    } else{
      
      # combined_melt <- combined[,paste0(div_methods$method,"_",par,"_",b)] %>%
      #   melt()
      combined_scale <- scale(combined[,paste0(div_methods$method,"_",par,"_",b)]) %>% 
        as.data.frame() 
    }
    
    # The scaled_version needs
    combined_scale_melt <- melt(combined_scale)
    gradient_max <- max(combined_scale)
    gradient_min <- min(combined_scale)
    
    
    combined_scale$Index <- Branch_index
    combined$Index <- Branch_index
    
    
    ggplot_list <- list()
    test_list <- list()
    count_method = 0 
    individual_plots <- vector('list', ncolumns)
    individual_plots_scaled <- vector('list', ncolumns)
    
    for(t in div_methods$method){
      count_method = count_method + 1
      
      # non scaled
      if(t == "PESTO") {
        combined_Plot <- combined[, c(paste0(t,"_",b),"Index") ]
      } else{
        combined_Plot <- combined[, c(paste0(t,"_",par,"_",b),"Index") ]
      }
      names(combined_Plot) <- c("value","node") 
      phylo_tib_data <- dplyr::full_join(phylo_tib, combined_Plot, by = 'node') %>% treeio::as.treedata()
      individual_plots[[count_method]] <- local({plot_rate_phylo(phylo_tib_data)})
      #BAMM_rate_Phy_noLeg <- plot_rate_phylo_remove_legend_andScale(BAMM_rate_Phy_mu)
      #names(individual_plots[[count_method]]) <-t
      
      if(t == "PESTO") {
        combined_scale_plot <- combined_scale[ , c(paste0(t,"_",b),"Index")]
      } else{
        combined_scale_plot <- combined_scale[ , c(paste0(t,"_",par,"_",b),"Index")]
      }
      
      names(combined_scale_plot) <- c("value","node") 
      phylo_tib_data <- dplyr::full_join(phylo_tib, combined_scale_plot, by = 'node')%>% treeio::as.treedata()
      
      individual_plots_scaled[[count_method]] <- local({ (plot_rate_phylo_scaled(phylo_tib_data,
                                                                                 gradient_max,
                                                                                 gradient_min))})

    }
    
    # Individual Legends legends phylogeny plots
    row_plots[[count_par]] <- local({ plot_grid(plotlist =individual_plots,
                                                labels = as.character(as.matrix(plot_labels[count_par])),
                                                label_size = 8,
                                                nrow = 1) })
    
    # Scaled Phylogenies Plots
    legend_scaled <- cowplot::get_legend( individual_plots_scaled[[count_method]] ) # Shared legend  
    
    individual_plots_scaled_removed <- lapply(individual_plots_scaled,
                                              FUN=function(x) x + theme(legend.position="none") )
    scale <-plot_grid(plotlist=individual_plots_scaled_removed,
                      labels = as.character(as.matrix(plot_labels[count_par])), 
                      label_size = 8,
                      nrow = 1,
                      label_x = 0, # Makes all labels left justified, I think. 
                      hjust = 0)
    
    # Generates titles for row
    title <- ggdraw() + 
      draw_label( label = rate_labels$symbol[grep(x=rate_labels$div,pattern = b)],
                  fontface = 'bold', x = 0, hjust = .55, size=7.5,
                  angle = 90) +
      theme(
        # add margin on the left of the drawing canvas,
        # so title is aligned with left edge of first plot
        plot.margin = margin(0, 0, 0, 7)
      )
    
    row_plots_scaled[[count_par]] <- local({ 
      plot_grid(title, scale, legend_scaled, rel_widths = c(.1, 3, .4), nrow=1) 
    })
    
  }
  
  individual_legends <- plot_grid(row_plots[[1]],
                                  row_plots[[2]],
                                  row_plots[[3]],
                                  ncol = 1) 
  
  scaled_legends <- plot_grid(row_plots_scaled[[1]],
                              row_plots_scaled[[2]],
                              row_plots_scaled[[3]],
                              ncol = 1) 
  if (print_scale == T){
    return(scaled_legends)
  }
  if(print_scale == F){
    return(individual_legends)
  }
}


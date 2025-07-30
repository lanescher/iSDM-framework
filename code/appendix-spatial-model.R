
# Load functions
# source('functions/FXN-MVPv1.1.R')

# Set common variables
# sp.code <- 'RACA'

# Define functions

# Adapted from plot_pars function
pull_pars <- function(out) {
  
  cov.labs <- read.csv("data/covariate-labels.csv")
  
  dat <- out$process.coef %>%
    mutate(x = covariate)
  xlab <- "Covariate"
  
  # rename interaction reference level
  ints <- grep("_x_", dat$x)
  if (length(ints) == 2) { # interaction with quadratic term
    cov <- gsub("_x_.*", "", dat$x[ints])
    newname <- grep(paste0(cov, collapse = "|"), dat$x)[1:2]
    newname1 <- paste0(dat$x[newname], "_x_reference")
    dat$x[newname] <- newname1
    
  } else if (length(ints) == 1) { # interaction with linear term
    cov <- gsub("_x_.*", "", dat$x[ints])
    newname <- grep(paste0(cov, collapse = "|"), dat$x)[1]
    newname1 <- paste0(dat$x[newname], "_x_reference")
    dat$x[newname] <- newname1
  }
  
  dat <- dat %>%
    mutate(cov1 = gsub("2", "", covariate),
           tmp = gsub("_x_.*", "", covariate),
           quad = case_when(substr(tmp, nchar(tmp), nchar(tmp)) == 2 ~ "^2",
                            T ~ "")) %>%
    left_join(cov.labs, by = c("cov1" = "covariate")) %>%
    mutate(x = paste0(Label, quad)) %>%
    select(!tmp)
  
  return(dat)
  
}

# Adapted from plot_chains function
plot_chains <- function(samples,
                        data,
                        B = T,
                        tau = F) {
  
  chaincols <-  c("1" = "hotpink1", "2" = "olivedrab3", "3" = "deepskyblue3")
  
  do <- c()
  if (B == T) do <- c(do, "B")
  if (tau == T) do <- c(do, "tau")
  
  
  for (d in do) {
    
    if (d == "B") {
      bnames <- data.frame(name = colnames(data$Xz),
                           param = paste0("B[", 1:ncol(data$Xz), "]"))
      bind <- grep("B", colnames(samples[[1]]))
      rm <- grep("XB", colnames(samples[[1]]))
      bind <- bind[-which(bind %in% rm)]
      
    } else if (d == "tau") {
      bind <- grep("tau", colnames(samples[[1]]))
      bnames <- data.frame(name = "tau",
                           param = paste0("tau[", 1:length(bind), "]"))
    }
    
    
    
    # select which chains to use
    chains <- 1:length(samples)
    
    if (length(chains) == 0) chains <- 1:3
    
    # calculate rhat
    all <- list()
    for (i in 1:length(bind)) {
      df <- as.matrix(data.frame(chain1 = as.vector(samples[[1]][,bind[i]]),
                                 chain2 = as.vector(samples[[2]][,bind[i]]),
                                 chain3 = as.vector(samples[[3]][,bind[i]])))
      # df <- df[cutoff:nrow(df),chains]
      all[[i]] <- df
    }
    
    
    rhat <- unlist(lapply(all, rstan::Rhat))
    bnames$rhat <- round(rhat, 2)
    
    
    # plot chains
    call <- c()
    for (c in chains) {
      c1 <- as.data.frame(samples[[c]][,bind]) %>%
        mutate(n = 1:nrow(.)) %>%
        pivot_longer(!n) %>%
        mutate(chain = c)
      call <- bind_rows(call, c1)
    }
    
    # If there's only one, need to rename it
    if (length(bind) == 1) {
      call$name <- paste0(d, "[", 1, "]")
    }
    
    call <- call %>%
      full_join(bnames, by = c("name" = "param"))
    
    if (d == "B") {
      cov.labs <- read.csv("data/covariate-labels.csv")
      
      call <- call %>%
        mutate(cov1 = gsub("2", "", name.y),
               tmp = gsub("_x_.*", "", name.y),
               quad = case_when(substr(tmp, nchar(tmp), nchar(tmp)) == 2 ~ "^2",
                                T ~ "")) %>%
        left_join(cov.labs, by = c("cov1" = "covariate")) %>%
        mutate(name.y = paste0(Label, quad)) %>%
        select(!tmp)
    }
    
    
    call <- call %>%
      mutate(lab = ifelse(name.y != "tau", paste0(name.y, "\nRhat = ", rhat), paste0("Rhat = ", rhat)))
    
    
    pl <- ggplot(call) +
      geom_line(aes(x = n, y = value, color = as.factor(chain)), linewidth = 0.1) +
      # geom_vline(xintercept = cutoff) +
      facet_wrap(~lab, scales = "free_y") +
      labs(x = "Iteration", y = "Parameter value", color = "Chain") +
      scale_color_manual(values = chaincols) +
      theme_bw()
    
    sub <- paste0(chains, collapse = "")
    
    return(list(plot = pl, rhat = bnames))
    
  }
  
}

# Adapted from map_species_data function
pull_data <- function(plot, out, region) {
  # put all values together
  ind <- grep(plot, names(out))
  intensity <- c()
  for (i in 1:length(ind)) {
    
    tmp <- out[[ind[i]]] %>%
      mutate(scenario = gsub(plot, "", names(out)[[ind[i]]]))
    intensity <- bind_rows(intensity, tmp)
    
  }
  
  # get correct column to plot
  intensity$plot.val <- intensity$mean
  
  # adjust values
  if (plot != "spat") {
    q99 <- quantile(intensity$plot.val, 0.99, na.rm = T)
    intensity$plot.val[which(intensity$plot.val > q99)] <- q99
  }
  
  
  intensity <- full_join(region$sp.grid, intensity, by = "conus.grid.id") %>%
    select(conus.grid.id, plot.val)
  
  return(intensity)
}

# Adapted from map_species_data function
map_species_data <- function(intensity,
                             region,
                             plot = "lambda") {
  
  
  # 1. load map features
  na <- rnaturalearth::ne_countries(continent = "North America", 
                                    returnclass = "sf", 
                                    scale = 10) %>%
    st_transform(crs = st_crs(region$range))
  
  wa <- rnaturalearth::ne_download(type = "lakes", 
                                   category = "physical", load = T,
                                   returnclass = "sf",
                                   scale = 10) %>%
    st_transform(crs = st_crs(region$range))
  
  st <- rnaturalearth::ne_states(country = c("Canada", "Mexico", "United States of America"),
                                 returnclass = "sf") %>%
    st_transform(crs = st_crs(region$range))
  
  
  # 2. make base map
  base <- ggplot() +
    geom_sf(data = na, fill = "gray90") +
    geom_sf(data = st, fill = "gray90", color= "gray40") +
    geom_sf(data = wa, fill = "lightsteelblue") +
    theme(panel.background = element_rect(fill = "lightsteelblue"),
          panel.grid = element_blank(),
          axis.text = element_text(size = 10),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 11))
  
  
  # __c. fill color
  
  if (plot == "lambda") {
    fill <- "Relative abundance"
    digits <- 7
    tit.vjust <- 1
  } else if (plot == "spat") {
    fill <- "Spatial effect"
    digits <- 2
    tit.vjust <- 1
  } 
  unc <- ""
  

  # set color scale
  if (plot != "spat") {
    
    # set color scheme
    base <- base +
      viridis::scale_fill_viridis(guide = guide_colorbar(), 
                                  option = "magma",
                                  na.value = NA,
                                  direction = -1) +
      viridis::scale_color_viridis(guide = "none",
                                   option = "magma",
                                   na.value = NA,
                                   direction = -1)
    
    
    base <- base + 
      geom_sf(data = intensity, aes(fill = plot.val, color = plot.val)) +
      facet_wrap(~ tau, nrow = 1, labeller = label_bquote(tau == .(tau))) +
      guides(fill = guide_colorbar(theme = theme(
        legend.key.height = unit(0.75, "lines"),
        legend.key.width = unit(20, "lines")), title.hjust = 0.9, title.vjust = tit.vjust)) +
      labs(fill = paste0(fill, unc), color = paste0(fill, unc))
    
  } else { # if plot == 'spat'
    base <- base +
      scale_fill_gradient2(guide = guide_colorbar(), high = "darkblue", low = "darkred", mid = "white", na.value = NA) +
      scale_color_gradient2(guide = "none", high = "darkblue", low = "darkred", mid = "white", na.value = NA)
    
    intensity <- filter(intensity, is.na(plot.val) == F)
    
    base <- base + 
      geom_sf(data = intensity, aes(fill = plot.val, color = plot.val)) +
      facet_wrap(~ tau, nrow = 1, labeller = label_bquote(tau == .(tau))) +
      guides(fill = guide_colorbar(theme = theme(
        legend.key.height = unit(0.75, "lines"),
        legend.key.width = unit(20, "lines")), title.hjust = 1, title.vjust = tit.vjust)) +
      labs(fill = paste0(fill, unc), color = paste0(fill, unc))
  }
  
  
  # put state lines and water on top of fill color
  base <- base +
    geom_sf(data = st, fill = NA) +
    geom_sf(data = wa, fill = "lightsteelblue")
  
  
  
  # 4. add limits
  bb <- st_bbox(region$region)
  xlim <- c(bb[1], bb[3])
  ylim <- c(bb[2], bb[4])
  
  base <- base +
    coord_sf(xlim = xlim, ylim = ylim)
  
  
  # 5. change theme
  base <- base +
    theme(legend.position = "bottom",
          axis.text.x = element_text(angle = 45,
                                     hjust = 1,
                                     vjust = 1),
          strip.text = element_text(size = 11))
  
  base
  
} 

path = 'outputs/proj08-model-framework/figures/'

# Reference model -----
test <- '31_RACA_test'

load(paste0("~/GitHub/species-futures/outputs/03-species-models/MVPv1/RACA-spatial-model-tests/",test,"/datafull-info.rdata"))
load(paste0("~/GitHub/species-futures/outputs/03-species-models/MVPv1/RACA-spatial-model-tests/",test,"/datafull.rdata"))
load(paste0("~/GitHub/species-futures/outputs/03-species-models/MVPv1/RACA-spatial-model-tests/",test,"/region.rdata"))
samples <- readRDS(paste0("~/GitHub/species-futures/outputs/03-species-models/MVPv1/RACA-spatial-model-tests/",test,"/samples_none.rds"))


# Lambda
ref.lam <- pull_data(plot = 'lambda', out = out, region = region) %>% mutate(tau = '1')

# Spatial effect
ref.spat <- pull_data(plot = 'spat', out = out, region = region) %>% mutate(tau = '1')

# Chains
ref.chains <- plot_chains(samples = samples, data = data, B = T, tau = F)

# Parameters
ref.dat <- pull_pars(out) %>% 
  mutate(tau = '1', zeromean = T) %>%
  select(-rhat) %>%
  left_join(., select(ref.chains$rhat, -param), by = c('covariate' = 'name'))





# FIXED VALUES ; ZERO_MEAN=T -----

### Tau = 0.1 -----
test <- '33_RACA_tau0.1'

load(paste0("~/GitHub/species-futures/outputs/03-species-models/MVPv1/RACA-spatial-model-tests/",test,"/datafull-info.rdata"))
load(paste0("~/GitHub/species-futures/outputs/03-species-models/MVPv1/RACA-spatial-model-tests/",test,"/datafull.rdata"))
load(paste0("~/GitHub/species-futures/outputs/03-species-models/MVPv1/RACA-spatial-model-tests/",test,"/region.rdata"))
samples <- readRDS(paste0("~/GitHub/species-futures/outputs/03-species-models/MVPv1/RACA-spatial-model-tests/",test,"/samples_none.rds"))

# Lambda
tau0.1.lam <- pull_data(plot = 'lambda', out = out, region = region) %>% mutate(tau = '0.1')

# Spatial effect
tau0.1.spat <- pull_data(plot = 'spat', out = out, region = region) %>% mutate(tau = '0.1')

# Chains
tau0.1.chains <- plot_chains(samples = samples, data = data, B = T, tau = F)

# Parameters
pull_pars(out) %>% 
  mutate(tau = '0.1', zeromean = T) %>%
  select(-rhat) %>%
  left_join(., select(tau0.1.chains$rhat, -param), by = c('covariate' = 'name')) %>%
  bind_rows(ref.dat,.) -> dat


### Tau = 10 -----
test <- '34_RACA_tau10'

load(paste0("~/GitHub/species-futures/outputs/03-species-models/MVPv1/RACA-spatial-model-tests/",test,"/datafull-info.rdata"))
load(paste0("~/GitHub/species-futures/outputs/03-species-models/MVPv1/RACA-spatial-model-tests/",test,"/datafull.rdata"))
load(paste0("~/GitHub/species-futures/outputs/03-species-models/MVPv1/RACA-spatial-model-tests/",test,"/region.rdata"))
samples <- readRDS(paste0("~/GitHub/species-futures/outputs/03-species-models/MVPv1/RACA-spatial-model-tests/",test,"/samples_none.rds"))


# Lambda
tau10.lam <- pull_data(plot = 'lambda', out = out, region = region) %>% mutate(tau = '10')

# Spatial effect
tau10.spat <- pull_data(plot = 'spat', out = out, region = region) %>% mutate(tau = '10')

# Chains
tau10.chains <- plot_chains(samples = samples, data = data, B = T, tau = F)

# Parameters
pull_pars(out) %>% 
  mutate(tau = '10', zeromean = T) %>%
  select(-rhat) %>%
  left_join(., select(tau10.chains$rhat, -param), by = c('covariate' = 'name')) %>%
  bind_rows(dat,.) -> dat





# Create plots -----

### Spatial effect ----
spat <- bind_rows(ref.spat, tau0.1.spat, tau10.spat)

map_species_data(spat, region = region, plot = 'spat')

ggsave('FigA1S1-spat.jpg',
       path = path,
       dpi = 600,
       width = 6.5,
       height = 3.5,
       units = 'in')


### Lambda ----
int <- bind_rows(ref.lam, tau0.1.lam, tau10.lam)

map_species_data(int, region = region, plot = 'lambda')

ggsave('FigA1S2-lambda.jpg',
       path = path,
       dpi = 600,
       width = 6.5,
       height = 3.5,
       units = 'in')


### Chains
cowplot::plot_grid(tau0.1.chains$plot + theme(legend.position = 'none'),
                   ref.chains$plot + theme(legend.position = 'none'), 
                   tau10.chains$plot  + theme(legend.position = 'none'), 
                   nrow = 3,labels = 'AUTO')

ggsave('FigA1S3-chains.jpg',
       path = path,
       dpi = 600,
       width = 8,
       height = 11,
       units = 'in')

# Parameter effects
dat %>%
  mutate(converge = ifelse(rhat < 1.1, 'Yes','No')) %>%
  ggplot() +
    geom_hline(yintercept = 0) +
    # scale_x_discrete(labels = scales::label_wrap(15)) +
    geom_point(aes(x = x, y = mean, col = tau, group = tau, shape = converge),
               position = position_dodge(0.8), size = 2.5) +
    geom_linerange(aes(x = x, ymin = lo, ymax = hi, col = tau, group = tau), 
                  position = position_dodge(0.8)) +
    scale_shape_manual('Rhat < 1.1', values = c(1, 19)) +
    labs(y = 'Estimate', x = '', color = expression("Value of "*tau)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45,
                                     hjust = 1,
                                     vjust = 1))

ggsave('FigA1S4-coef.jpg',
       path = path,
       dpi = 600,
       width = 8,
       height = 6,
       units = 'in')




# FIXED VALUES ; ZERO_MEAN=F -----

### Tau = 1 ----
test <- '44_RACA_tau1zmF'

load(paste0("~/GitHub/species-futures/outputs/03-species-models/MVPv1/RACA-spatial-model-tests/",test,"/datafull-info.rdata"))
load(paste0("~/GitHub/species-futures/outputs/03-species-models/MVPv1/RACA-spatial-model-tests/",test,"/datafull.rdata"))
load(paste0("~/GitHub/species-futures/outputs/03-species-models/MVPv1/RACA-spatial-model-tests/",test,"/region.rdata"))
samples <- readRDS(paste0("~/GitHub/species-futures/outputs/03-species-models/MVPv1/RACA-spatial-model-tests/",test,"/samples_none.rds"))


# Lambda
tau1zmf.lam <- pull_data(plot = 'lambda', out = out, region = region) %>% mutate(tau = '1')

# Spatial effect
tau1zmf.spat <- pull_data(plot = 'spat', out = out, region = region) %>% mutate(tau = '1')

# Chains
tau1zmf.chains <- plot_chains(samples = samples, data = data, B = T, tau = F)

# Parameters
pull_pars(out) %>% 
  mutate(tau = '1', zeromean = F) %>%
  select(-rhat) %>%
  left_join(., select(tau1zmf.chains$rhat, -param), by = c('covariate' = 'name')) %>%
  bind_rows(dat,.) -> dat



### Tau = 0.1 ----
test <- '45_RACA_tau0.1zmF'

load(paste0("~/GitHub/species-futures/outputs/03-species-models/MVPv1/RACA-spatial-model-tests/",test,"/datafull-info.rdata"))
load(paste0("~/GitHub/species-futures/outputs/03-species-models/MVPv1/RACA-spatial-model-tests/",test,"/datafull.rdata"))
load(paste0("~/GitHub/species-futures/outputs/03-species-models/MVPv1/RACA-spatial-model-tests/",test,"/region.rdata"))
samples <- readRDS(paste0("~/GitHub/species-futures/outputs/03-species-models/MVPv1/RACA-spatial-model-tests/",test,"/samples_none.rds"))

# Lambda
tau0.1zmf.lam <- pull_data(plot = 'lambda', out = out, region = region) %>% mutate(tau = '0.1')

# Spatial effect
tau0.1zmf.spat <- pull_data(plot = 'spat', out = out, region = region) %>% mutate(tau = '0.1')

# Chains
tau0.1zmf.chains <- plot_chains(samples = samples, data = data, B = T, tau = F)

# Pull out parameters
pull_pars(out) %>% 
  mutate(tau = '0.1', zeromean = F) %>%
  select(-rhat) %>%
  left_join(., select(tau0.1zmf.chains$rhat, -param), by = c('covariate' = 'name')) %>%
  bind_rows(dat,.) -> dat




### Tau = 10 ----
test <- '46_RACA_tau10zmF'

load(paste0("~/GitHub/species-futures/outputs/03-species-models/MVPv1/RACA-spatial-model-tests/",test,"/datafull-info.rdata"))
load(paste0("~/GitHub/species-futures/outputs/03-species-models/MVPv1/RACA-spatial-model-tests/",test,"/datafull.rdata"))
load(paste0("~/GitHub/species-futures/outputs/03-species-models/MVPv1/RACA-spatial-model-tests/",test,"/region.rdata"))
samples <- readRDS(paste0("~/GitHub/species-futures/outputs/03-species-models/MVPv1/RACA-spatial-model-tests/",test,"/samples_none.rds"))


# Lambda
tau10zmf.lam <- pull_data(plot = 'lambda', out = out, region = region) %>% mutate(tau = '10')

# Spatial effect
tau10zmf.spat <- pull_data(plot = 'spat', out = out, region = region) %>% mutate(tau = '10')

# Chains
tau10zmf.chains <- plot_chains(samples = samples, data = data, B = T, tau = F)

# Pull out parameters
pull_pars(out) %>% 
  mutate(tau = '10', zeromean = F) %>%
  select(-rhat) %>%
  left_join(., select(tau10zmf.chains$rhat, -param), by = c('covariate' = 'name')) %>%
  bind_rows(dat,.) -> dat



# Create plots -----

### Spatial effect ----
spat <- bind_rows(tau1zmf.spat, tau0.1zmf.spat, tau10zmf.spat)

map_species_data(spat, region = region, plot = 'spat')

ggsave('FigA1S5-spat.jpg',
       path = path,
       dpi = 600,
       width = 6.5,
       height = 3.5,
       units = 'in')


### Lambda ----
int <- bind_rows(tau1zmf.lam, tau0.1zmf.lam, tau10zmf.lam)

map_species_data(int, region = region, plot = 'lambda')

ggsave('FigA1S6-lambda.jpg',
       path = path,
       dpi = 600,
       width = 6.5,
       height = 3.5,
       units = 'in')


### Chains ----
cowplot::plot_grid(tau0.1zmf.chains$plot + theme(legend.position = 'none'),
                   tau1zmf.chains$plot + theme(legend.position = 'none'),  
                   tau10zmf.chains$plot  + theme(legend.position = 'none'), 
                   nrow = 3, labels = "AUTO")

ggsave('FigA1S7-chains.jpg',
       path = path,
       dpi = 600,
       width = 8,
       height = 11,
       units = 'in')

### Parameter effects ----
dat %>%
  mutate(converge = ifelse(rhat < 1.1, 'Yes','No')) %>%
  ggplot() +
    geom_hline(yintercept = 0) +
    # scale_x_discrete(labels = scales::label_wrap(15)) +
    geom_point(aes(x = x, y = mean, col = tau, group = tau, shape = converge),
               position = position_dodge(0.9)) +
    geom_linerange(aes(x = x, ymin = lo, ymax = hi, col = tau, group=tau), 
                   position = position_dodge(0.9)) +
    scale_shape_manual('Rhat < 1.1', values = c(1, 19)) +
    facet_wrap(~zeromean, nrow = 1) +
    labs(y = 'Estimate', x = '', color = expression("Value of "*tau)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45,
                                     hjust = 1,
                                     vjust = 1))

ggsave('FigA1S8-coef.jpg',
       path = path,
       dpi = 600,
       width = 8,
       height = 6,
       units = 'in')




# PRIOR VALUES; ZERO_MEAN=T -----

### ~ dgamma(5,5) ----
test <- '37_RACA_tauGam5-5'

load(paste0("~/GitHub/species-futures/outputs/03-species-models/MVPv1/RACA-spatial-model-tests/",test,"/datafull-info.rdata"))
load(paste0("~/GitHub/species-futures/outputs/03-species-models/MVPv1/RACA-spatial-model-tests/",test,"/datafull.rdata"))
load(paste0("~/GitHub/species-futures/outputs/03-species-models/MVPv1/RACA-spatial-model-tests/",test,"/region.rdata"))
samples <- readRDS(paste0("~/GitHub/species-futures/outputs/03-species-models/MVPv1/RACA-spatial-model-tests/",test,"/samples_none.rds"))

# Lambda
tauGam5.lam <- pull_data(plot = 'lambda', out = out, region = region) %>% mutate(tau = 'dgamma(5,5)')

# Spatial effect
tauGam5.spat <- pull_data(plot = 'spat', out = out, region = region) %>% mutate(tau = 'dgamma(5,5)')

# Chains
tauGam5.chains <- plot_chains(samples = samples, data = data, B = T, tau = F)

tauGam5.tau <- plot_chains(samples = samples, data = data, B = F, tau = T)

# Parameters
pull_pars(out) %>% 
  mutate(tau = 'dgamma(5,5)', zeromean = T) %>%
  select(-rhat) %>%
  left_join(., select(tauGam5.chains$rhat, -param), by = c('covariate' = 'name')) %>%
  bind_rows(dat,.) -> dat





### ~ dnorm(1,0.1) ----
test <- '38_RACA_tauNorm1-0.1'

load(paste0("~/GitHub/species-futures/outputs/03-species-models/MVPv1/RACA-spatial-model-tests/",test,"/datafull-info.rdata"))
load(paste0("~/GitHub/species-futures/outputs/03-species-models/MVPv1/RACA-spatial-model-tests/",test,"/datafull.rdata"))
load(paste0("~/GitHub/species-futures/outputs/03-species-models/MVPv1/RACA-spatial-model-tests/",test,"/region.rdata"))
samples <- readRDS(paste0("~/GitHub/species-futures/outputs/03-species-models/MVPv1/RACA-spatial-model-tests/",test,"/samples_none.rds"))


# Lambda
tauNorm1.lam <- pull_data(plot = 'lambda', out = out, region = region) %>% mutate(tau = 'dnorm(1,0.1)')

# Spatial effect
tauNorm1.spat <- pull_data(plot = 'spat', out = out, region = region) %>% mutate(tau = 'dnorm(1,0.1)')

# Chains
tauNorm1.chains <- plot_chains(samples = samples, data = data, B = T, tau = F)

tauNorm1.tau <- plot_chains(samples = samples, data = data, B = F, tau = T)

# Parameters
pull_pars(out) %>% 
  mutate(tau = 'dnorm(1,0.1)', zeromean = T) %>%
  select(-rhat) %>%
  left_join(., select(tauNorm1.chains$rhat, -param), by = c('covariate' = 'name')) %>%
  bind_rows(dat,.) -> dat





### ~ dgamma(0.01, 0.01) ----
test <- '39_RACA_tauGam0.01-0.01'

load(paste0("~/GitHub/species-futures/outputs/03-species-models/MVPv1/RACA-spatial-model-tests/",test,"/datafull-info.rdata"))
load(paste0("~/GitHub/species-futures/outputs/03-species-models/MVPv1/RACA-spatial-model-tests/",test,"/datafull.rdata"))
load(paste0("~/GitHub/species-futures/outputs/03-species-models/MVPv1/RACA-spatial-model-tests/",test,"/region.rdata"))
samples <- readRDS(paste0("~/GitHub/species-futures/outputs/03-species-models/MVPv1/RACA-spatial-model-tests/",test,"/samples_none.rds"))


# Lambda
tauGam0.01.lam <- pull_data(plot = 'lambda', out = out, region = region) %>% mutate(tau = 'dgamma(0.01,0.01)')

# Spatial effect
tauGam0.01.spat <- pull_data(plot = 'spat', out = out, region = region) %>% mutate(tau = 'dgamma(0.01,0.01)')

# Chains
tauGam0.01.chains <- plot_chains(samples = samples, data = data, B = T, tau = F)

tauGam0.01.tau <- plot_chains(samples = samples, data = data, B = F, tau = T)

# Parameters
pull_pars(out) %>% 
  mutate(tau = 'dgamma(0.01,0.01)', zeromean = T) %>%
  select(-rhat) %>%
  left_join(., select(tauGam0.01.chains$rhat, -param), by = c('covariate' = 'name')) %>%
  bind_rows(dat,.) -> dat





# Create plots -----

### Spatial effect ----
spat <- bind_rows(tauGam5.spat, tauNorm1.spat, tauGam0.01.spat)

map_species_data(spat, region = region, plot = 'spat')

ggsave('FigA1S9-spat.jpg',
       path = path,
       dpi = 600,
       width = 6.5,
       height = 3.5,
       units = 'in')


### Lambda ----
int <- bind_rows(tauGam5.lam, tauNorm1.lam, tauGam0.01.lam)

map_species_data(int, region = region, plot = 'lambda')

ggsave('FigA1S10-lambda.jpg',
       path = path,
       dpi = 600,
       width = 6.5,
       height = 3.5,
       units = 'in')


### Chains ----
cowplot::plot_grid(tauGam5.chains$plot + theme(legend.position = 'none'),
                   tauNorm1.chains$plot + theme(legend.position = 'none'), 
                   tauGam0.01.chains$plot + theme(legend.position = 'none'), 
                   nrow = 3,labels = 'AUTO')

ggsave('FigA1S11-chains.jpg',
       path = path,
       dpi = 600,
       width = 8,
       height = 11,
       units = 'in')

### Tau ----
cowplot::plot_grid(tauGam5.tau$plot + theme(legend.position = 'none'),
                   tauNorm1.tau$plot + theme(legend.position = 'none'), 
                   tauGam0.01.tau$plot + theme(legend.position = 'none'), 
                   nrow = 1,labels = 'AUTO')

ggsave('FigA1S12-tau.jpg',
       path = path,
       dpi = 600,
       width = 6.5,
       height = 2,
       units = 'in')

### Parameter effects ----
dat %>%
  mutate(converge = ifelse(rhat < 1.1, 'Yes','No')) %>%
  ggplot() +
    geom_hline(yintercept = 0) +
    # scale_x_discrete(labels = scales::label_wrap(15)) +
    geom_point(aes(x = x, y = mean, col = tau, group = tau, shape = converge),
               position = position_dodge(0.9)) +
    geom_linerange(aes(x = x, ymin = lo, ymax = hi, col = tau, group=tau), 
                   position = position_dodge(0.9)) +
    scale_shape_manual('Rhat < 1.1', values = c(1, 19)) +
    facet_wrap(~zeromean, nrow = 1) +
    labs(y = 'Estimate', x = '', color = expression("Value of "*tau)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45,
                                     hjust = 1,
                                     vjust = 1))

ggsave('FigA1S13-coef.jpg',
       path = path,
       dpi = 600,
       width = 8,
       height = 6,
       units = 'in')



# PRIOR VALUES; ZERO_MEAN=F -----

### ~ dgamma(5,5) ----
test <- '40_RACA_tau37zmF'

load(paste0("~/GitHub/species-futures/outputs/03-species-models/MVPv1/RACA-spatial-model-tests/",test,"/datafull-info.rdata"))
load(paste0("~/GitHub/species-futures/outputs/03-species-models/MVPv1/RACA-spatial-model-tests/",test,"/datafull.rdata"))
load(paste0("~/GitHub/species-futures/outputs/03-species-models/MVPv1/RACA-spatial-model-tests/",test,"/region.rdata"))
samples <- readRDS(paste0("~/GitHub/species-futures/outputs/03-species-models/MVPv1/RACA-spatial-model-tests/",test,"/samples_none.rds"))

# Lambda
tauGam5zmf.lam <- pull_data(plot = 'lambda', out = out, region = region) %>% mutate(tau = 'dgamma(5,5)')

# Spatial effect
tauGam5zmf.spat <- pull_data(plot = 'spat', out = out, region = region) %>% mutate(tau = 'dgamma(5,5)')

# Chains
tauGam5zmf.chains <- plot_chains(samples = samples, data = data, B = T, tau = F)

tauGam5zmf.tau <- plot_chains(samples = samples, data = data, B = F, tau = T)

# Parameters
pull_pars(out) %>% 
  mutate(tau = 'dgamma(5,5)', zeromean = F) %>%
  select(-rhat) %>%
  left_join(., select(tauGam5zmf.chains$rhat, -param), by = c('covariate' = 'name')) %>%
  bind_rows(dat,.) -> dat



### ~ dnorm(1,0.1) ----
test <- '41_RACA_tau38zmF'

load(paste0("~/GitHub/species-futures/outputs/03-species-models/MVPv1/RACA-spatial-model-tests/",test,"/datafull-info.rdata"))
load(paste0("~/GitHub/species-futures/outputs/03-species-models/MVPv1/RACA-spatial-model-tests/",test,"/datafull.rdata"))
load(paste0("~/GitHub/species-futures/outputs/03-species-models/MVPv1/RACA-spatial-model-tests/",test,"/region.rdata"))
samples <- readRDS(paste0("~/GitHub/species-futures/outputs/03-species-models/MVPv1/RACA-spatial-model-tests/",test,"/samples_none.rds"))

# Lambda
tauNorm1zmf.lam <- pull_data(plot = 'lambda', out = out, region = region) %>% mutate(tau = 'dnorm(1,0.1)')

# Spatial effect
tauNorm1zmf.spat <- pull_data(plot = 'spat', out = out, region = region) %>% mutate(tau = 'dnorm(1,0.1)')

# Chains
tauNorm1zmf.chains <- plot_chains(samples = samples, data = data, B = T, tau = F)

tauNorm1zmf.tau <- plot_chains(samples = samples, data = data, B = F, tau = T)

# Parameters
pull_pars(out) %>% 
  mutate(tau = 'dnorm(1,0.1)', zeromean = F) %>%
  select(-rhat) %>%
  left_join(., select(tauNorm1zmf.chains$rhat, -param), by = c('covariate' = 'name')) %>%
  bind_rows(dat,.) -> dat



### ~ dgamma(0.01, 0.01) ----
test <- '42_RACA_tau39zmF'

load(paste0("~/GitHub/species-futures/outputs/03-species-models/MVPv1/RACA-spatial-model-tests/",test,"/datafull-info.rdata"))
load(paste0("~/GitHub/species-futures/outputs/03-species-models/MVPv1/RACA-spatial-model-tests/",test,"/datafull.rdata"))
load(paste0("~/GitHub/species-futures/outputs/03-species-models/MVPv1/RACA-spatial-model-tests/",test,"/region.rdata"))
samples <- readRDS(paste0("~/GitHub/species-futures/outputs/03-species-models/MVPv1/RACA-spatial-model-tests/",test,"/samples_none.rds"))


# Lambda
tauGam0.01zmf.lam <- pull_data(plot = 'lambda', out = out, region = region) %>% mutate(tau = 'dgamma(0.01,0.01)')

# Spatial effect
tauGam0.01zmf.spat <- pull_data(plot = 'spat', out = out, region = region) %>% mutate(tau = 'dgamma(0.01,0.01)')

# Chains
tauGam0.01zmf.chains <- plot_chains(samples = samples, data = data, B = T, tau = F)

tauGam0.01zmf.tau <- plot_chains(samples = samples, data = data, B = F, tau = T)

# Parameters
pull_pars(out) %>% 
  mutate(tau = 'dgamma(0.01,0.01)', zeromean = F) %>%
  select(-rhat) %>%
  left_join(., select(tauGam0.01zmf.chains$rhat, -param), by = c('covariate' = 'name')) %>%
  bind_rows(dat,.) -> dat



# Create plots -----

### Spatial effect ---- 
spat <- bind_rows(tauGam5zmf.spat, tauNorm1zmf.spat, tauGam0.01zmf.spat)

map_species_data(spat, region = region, plot = 'spat')

ggsave('FigA1S14-spat.jpg',
       path = path,
       dpi = 600,
       width = 6.5,
       height = 3.5,
       units = 'in')


### Lambda ----
int <- bind_rows(tauGam5zmf.lam, tauNorm1zmf.lam, tauGam0.01zmf.lam)

map_species_data(int, region = region, plot = 'lambda')

ggsave('FigA1S15-lambda.jpg',
       path = path,
       dpi = 600,
       width = 6.5,
       height = 3.5,
       units = 'in')


### Chains ----
cowplot::plot_grid(tauGam5zmf.chains$plot + theme(legend.position = 'none'),
                   tauNorm1zmf.chains$plot + theme(legend.position = 'none'), 
                   tauGam0.01zmf.chains$plot  + theme(legend.position = 'none'), 
                   nrow = 3,labels = 'AUTO')

ggsave('FigA1S16-chains.jpg',
       path = path,
       dpi = 600,
       width = 8,
       height = 11,
       units = 'in')


### Tau ----
cowplot::plot_grid(tauGam5zmf.tau$plot + theme(legend.position = 'none'),
                   tauNorm1zmf.tau$plot + theme(legend.position = 'none'), 
                   tauGam0.01zmf.tau$plot  + theme(legend.position = 'none'), 
                   nrow = 1,labels = 'AUTO')

ggsave('FigA1S17-tau.jpg',
       path = path,
       dpi = 600,
       width = 6.5,
       height = 2,
       units = 'in')



### Parameter effects ----
dat %>%
  mutate(converge = ifelse(rhat < 1.1, 'Yes','No')) %>%
  ggplot() +
  geom_hline(yintercept = 0) +
  # scale_x_discrete(labels = scales::label_wrap(15)) +
  geom_point(aes(x = x, y = mean, col = tau, group = tau, shape = converge),
             position = position_dodge(0.9)) +
  geom_linerange(aes(x = x, ymin = lo, ymax = hi, col = tau, group=tau), 
                 position = position_dodge(0.9)) +
  scale_shape_manual('Rhat < 1.1', values = c(1, 19)) +
  facet_wrap(~zeromean, nrow = 1) +
  labs(y = 'Estimate', x = '', color = expression("Value of "*tau)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1,
                                   vjust = 1))

ggsave('FigA1S18-coef.jpg',
       path = path,
       dpi = 600,
       width = 8,
       height = 6,
       units = 'in')









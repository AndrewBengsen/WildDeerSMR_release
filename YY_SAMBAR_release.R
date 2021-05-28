## Site YY, Species SAMBAR

## Libraries and functions 
source("Scripts/SMR_functions.R")

inDir <- "Input/"
site <- "YY"

# Geographic data: state space and camera locations
dsn <- "Geographic"
cam.region <- readOGR(dsn, "S_YY_GDA94") %>%
  spTransform(CRS("+proj=utm +zone=55 +south +ellps=GRS80 +units=km")) %>%
  as.owin()

cams <- read.csv(paste0(inDir, "camlocs_MW.csv"), header=T) %>%
  filter(substr(as.character(Station), 1, 2) == site) %>%
  mutate(x = X, y = Y) %>%
  select(Station, x, y) %>%
  filter(Station != "YY_A4")

camlocs <- cams %>% 
  dplyr::select(x,y)

camlocs <- camlocs/1000
camlocs <- as.ppp(camlocs, W=cam.region)
locs <- coords(camlocs)

# Yk = 3D array of known ID's
IDarray <- readRDS(paste0(inDir, "SDH_count_SAMBAR_MW"))
YYCams <- grep(site, colnames(IDarray))
YYDeer <- grep(site, row.names(IDarray))
Yk <- IDarray[YYDeer, YYCams, ]

# Unmarked detection history matrix
Yu <- read.csv(paste0(inDir, "Yu_SAMBAR_MW.csv"), header=T) %>%
  filter(substr(station, 1, 2) == site) 

# Change order of Yu rows to match order of cams and Yk
cams_order <- data.frame(station = cams$Station,
                         order = seq(1,nrow(cams), 1))

Yu <- left_join(Yu, cams_order) %>%
  arrange(order) %>%
  select(-station, -order)

# Set model foundations
M <-200
mmax <- 50
delta<-c(0.25, 0.0015, 4) 
sigma.prior <- list("uniform", 0, 10)
lam0.prior <- list("uniform", 0, 5)
xlim <- range(locs[,1])
ylim <- range(locs[,2])
niters <- 150000
nb <- 5000
nc <- 6
set.seed(783)

# Run model
startDate <- Sys.Date()
logfile <- paste0("Logfiles/mylog_YY_SAMBAR_", startDate, ".txt")
cl <- makePSOCKcluster(nc, outfile=logfile)
clusterSetRNGStream(cl)
clusterEvalQ(cl, {library(spatstat)})
clusterExport(cl, c("SMR_dens", "calc.dist", "prior.density", "sample.prior"))

inits <- function(){list(S = cbind(runif(M+mmax, xlim[1], xlim[2]), 
                                   runif(M+mmax, ylim[1], ylim[2])), 
                        lam0=runif(1, 0.01, 0.5),
                        sigma=runif(1, 0.4, 2), psi=runif(1, 0.4, 0.6),
                        psim=runif(1, 0.4, 0.6))}

mod <- clusterCall(cl, SMR_dens, n = Yu, X = locs, y = Yk, M = M, mmax = mmax, 
                  obsmod = "pois", niters = niters, region = cam.region, 
                  sigma.prior = sigma.prior, lam0.prio = lam0.prior, 
                  inits = inits(), delta = delta)
stopCluster(cl)

# drop burnin
for(i in 1:nc) {
  mod[[i]] <- mod[[i]][-c(1:nb),] 
  mod[[i]] <- as.mcmc(mod[[i]])
}
mod <- mcmc.list(mod)

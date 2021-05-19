## Site BL, Species SAMBAR

## Libraries and functions 
source("Scripts/SMR_functions.R")

inDir <- "Input/"
Site <- "BL"
crs <- CRS("+proj=lcc +lat_1=-30.75 +lat_2=-35.75 +lat_0=-33.25 +lon_0=147 +x_0=9300000 +y_0=4500000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
crsKM <- CRS("+proj=lcc +lat_1=-30.75 +lat_2=-35.75 +lat_0=-33.25 +lon_0=147 +x_0=9300000 +y_0=4500000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=km +no_defs")

# Geographic data: 1) state space and 2) camera locations
dsn <- "Geographic"
cam.regionS <- readOGR(dsn, "S_BL_3km_3308")
# Specify CRS to output in km
cam.region <- spTransform(cam.regionS, crsKM) 
cam.region <- as.owin(cam.region)

cams <- read.csv(paste0(inDir, "camlocs_kosi.csv"), header=T) %>%
  mutate(x = X, y = Y) %>%
  filter(substr(Station, 1, 1) == tolower(substr(Site, 1, 1))) %>%
  select(Station, x, y, lat, lon) 
camlocs <- cams %>% dplyr::select(x,y)
camlocs <- camlocs/1000 # km
camlocs<- as.ppp(camlocs, W=cam.region)
locs<- coords(camlocs)

# Yk = 3D array of known ID's, filter out other sites
IDarray <- readRDS(paste0(inDir, "SDH_count_SAMBAR_KOSI"))
BLCams <- grep("b-", colnames(IDarray))
BLSams <- grep("_B", row.names(IDarray))
Yk <- IDarray[BLSams, BLCams, ]

# Unmarked detection history matrix
Yu <- read.csv(paste0(inDir, "Yu_SAMBAR_KOSI.csv"), header=T) %>%
  filter(tolower(substr(X, 1, 1)) == tolower(substr(Site, 1, 1))) %>%
  select(-detection_history.o90, -X)

# Set model foundations
M <- 150
mmax <- 50 
delta <- c(0.12, 0.01, 2) 
xlim <- range(locs[,1])
ylim <- range(locs[,2])
ni <- 350000
nb <- 25000
nc <- 4

# sigma prior from sambar posterior at site GG
sigma.prior <- list("gamma", 8.65, 30)  
lam0.prior <- list("uniform", 0, 5) 

seed <- 90
set.seed(seed)

# Run model
startDate <- Sys.Date()
ncores <- nc
logfile <- paste0("Logfiles/mylog_BL_Sam_", startDate, ".txt")
cl<- makePSOCKcluster(ncores, outfile=logfile)
clusterSetRNGStream(cl)
clusterEvalQ(cl, {library(spatstat)})
clusterExport(cl, c("SMR_dens", "calc.dist", "prior.density", "sample.prior"))

inits <- function(){list(S=cbind(runif(M+mmax, xlim[1], xlim[2]), 
                                runif(M+mmax, ylim[1], ylim[2])), 
                        lam0=runif(1, 0.01, 0.5),
                        sigma=runif(1, 0.4, 2), psi=runif(1, 0.4, 0.6),
                        psim=runif(1, 0.4, 0.6))}

mod <- clusterCall(cl, SMR_dens, n = Yu, X = locs, y = Yk, M = M, mmax = mmax, 
                  obsmod = "pois", niters = ni, region = cam.region, 
                  sigma.prior = sigma.prior, lam0.prio = lam0.prior, 
                  inits = inits(), delta = delta)
stopCluster(cl)

# drop burnin
for(i in 1:nc) {
  mod[[i]] <- mod[[i]][-c(1:nb),] 
  mod[[i]] <- as.mcmc(mod[[i]])
}
mod <- mcmc.list(mod)
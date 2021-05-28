## Site YP, Species RUSA

## Libraries and functions 
source("Scripts/SMR_functions.R")

inDir <- "Input/"
dsn <- "Geographic"

# Geographic data: state space and camera locations
# Projected coordinate system EPSG:28355 = GDA94 MGA55S
crs <- CRS("+proj=utm +zone=55 +south +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
crsKM <- CRS("+proj=utm +zone=55 +south +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=km +no_defs")
wgs84 <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

# Read and centre state space and camera stations
cam.region <- readOGR(dsn, "S_YP_GDA94")

centroidS <- rgeos::gCentroid(cam.region)

cam.regionS <- cam.region %>%
  broom::tidy() %>%
  mutate(long = long -centroidS$x,
         lat = lat - centroidS$y) 

cam.regionS <- Polygon(coords = as.matrix(cbind(cam.regionS$long, cam.regionS$lat)), hole=F)
cam.regionS <- Polygons(srl = list(cam.regionS), ID=1)
cam.regionS <- SpatialPolygons(Srl = list(cam.regionS), proj4string = crs)
cam.regionS <- SpatialPolygonsDataFrame(Sr = cam.regionS, data=as.data.frame(1), match.ID=F) %>%
  spTransform(crsKM) %>%
  as.owin()

camlocs <- read.csv(paste0(inDir, "camlocs_YP.csv")) %>%
  mutate(x = X-centroidS$x, y = Y-centroidS$y) %>% 
  dplyr::select(x,y)
camlocs <- camlocs/1000
camlocs <- as.ppp(camlocs, W=cam.regionS)
locs <- coords(camlocs)

# Load marked (Yk) and unmarked (Yu) detections
Yk <- readRDS(paste0(inDir, "SDH_count_RUSA_YP_subset"))
Yu <- read.csv(paste0(inDir, "Yu_RUSA_YP_subset.csv"), header=T) %>%
  dplyr::select(-X)

# Set model foundations
M <- 300
mmax <- 60
delta <- c(0.08, 0.16, 2)
sigma.prior <- list("uniform",0, 5)
lam0.prior <- list("uniform",0, 5)
xlim <- range(locs[,1])
ylim <- range(locs[,2])
ni <- 200000
nb <- 5000
nc <- 6
set.seed(53)

# Run model
startDate <- Sys.Date()
logfile <- paste0("mylog_YP_RUSA_", startDate, ".txt")
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
                   obsmod = "pois", niters = ni, region = cam.regionS, 
                   sigma.prior = sigma.prior, lam0.prio = lam0.prior, 
                   inits = inits(), delta = delta)
stopCluster(cl)

# drop burnin
for(i in 1:nc) {
  mod[[i]] <- mod[[i]][-c(1:nb),] 
  mod[[i]] <- as.mcmc(mod[[i]])
}
mod <- mcmc.list(mod)
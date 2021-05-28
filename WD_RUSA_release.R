## Site WD, Species RUSA

## Libraries and functions 
source("Scripts/SMR_functions.R")

inDir <- "Input/"
dsn <- "Geographic"

# Geographic data: state space and camera locations
# Projected coordinate system EPSG:28355 = GDA94 MGA56S
crs <- CRS("+proj=utm +zone=55 +south +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
crsKM <- CRS("+proj=utm +zone=55 +south +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=km +no_defs")
wgs84 <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

# Read and centre state space and camera stations
cam.region <- readOGR(dsn, "S_WD_GDA94")

centroidS <- rgeos::gCentroid(cam.region)

cam.region <- cam.region %>%
  broom::tidy() %>%
  mutate(long = long -centroidS$x,
         lat = lat - centroidS$y) 

cam.region <- Polygon(coords = as.matrix(cbind(cam.region$long, cam.region$lat)), hole=F)
cam.region <- Polygons(srl = list(cam.region), ID=1)
cam.region <- SpatialPolygons(Srl = list(cam.region), proj4string = crs)
cam.region <- SpatialPolygonsDataFrame(Sr = cam.region, data=as.data.frame(1), match.ID=F) %>%
  spTransform(crsKM) %>%
  as.owin()

camlocs <- read.csv(paste0(inDir, "camlocs_WD.csv")) %>%
  mutate(x = X-centroidS$x, y = Y-centroidS$y) %>% 
  dplyr::select(x,y)
camlocs <- camlocs/1000
camlocs <- as.ppp(camlocs, W=cam.region)
locs <- coords(camlocs)

# Load marked (Yk) and unmarked (Yu) detections
Yk <- readRDS(paste0(inDir, "SDH_count_RUSA_WD"))
Yu <- read.csv(paste0(inDir, "Yu_RUSA_WD.csv"), header=T) %>%
  dplyr::select(-X)

# Set model foundations
M <- 100
mmax <- 10
delta <- c(0.04, 0.007, 2)
sigma.prior <- list("uniform",0, 5)
lam0.prior <- list("uniform",0, 5)
xlim <- range(locs[,1])
ylim <- range(locs[,2])
ni<- 1000*750
nb<- 5000  
nc<- 6 
set.seed(53)

# Run model
startDate <- Sys.Date()
logfile <- paste0("mylog_WD_RUSA_", startDate, ".txt")
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

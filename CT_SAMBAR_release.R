## Site CT, Species SAMBAR

## Libraries and functions 
source("Scripts/SMR_functions.R")

inDir <- "Input/"

# Geographic data: 1) state space and 2) camera locations
# Projected coordinate system EPSG:3308
crs <- CRS("+proj=lcc +lat_1=-30.75 +lat_2=-35.75 +lat_0=-33.25 +lon_0=147 +x_0=9300000 +y_0=4500000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
crsKM <- CRS("+proj=lcc +lat_1=-30.75 +lat_2=-35.75 +lat_0=-33.25 +lon_0=147 +x_0=9300000 +y_0=4500000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=km +no_defs")
wgs84 <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

dsn <- "Geographic"

camlocs <- read.csv("Input/camlocs_CT.csv",
                    header=T, na.strings="na", skipNul=T) 

camlocsCRS <- WGS84toCRS(camlocs$lon, camlocs$lat)
camlocs$X <- camlocsCRS$X
camlocs$Y <- camlocsCRS$Y

camlocsSPDF <- SpatialPointsDataFrame(coords = data.frame(X = camlocs$X, Y = camlocs$Y), 
                                      data=camlocs, proj4string = crs, match.ID=F)

# Define the state space as a buffer around the outer cameras, trimmed of agricultural, urban and sea
cam.region1 <- readOGR(dsn, "S_CT_3km_3308")

# Find the sample space centroid
centroidS <- rgeos::gCentroid(cam.region1)

# Centre the sample space off the centroid
cam.regionS <- cam.region1 %>%
  broom::tidy() %>%
  mutate(long = long -centroidS$x,
         lat = lat - centroidS$y) 

cam.regionS <- Polygon(coords = as.matrix(cbind(cam.regionS$long, cam.regionS$lat)), hole=F)
cam.regionS <- Polygons(srl = list(cam.regionS), ID=1)
cam.regionS <- SpatialPolygons(Srl = list(cam.regionS), proj4string = crs)
cam.regionS <- SpatialPolygonsDataFrame(Sr = cam.regionS, data=as.data.frame(1), match.ID=F) %>%
  spTransform(crsKM) 

# Centre camlocs off the sample space centroid
camlocs <- camlocs %>%
  mutate(X = X-centroidS$x, Y = Y-centroidS$y)

cam.region <- as.owin(cam.regionS)

# 3D array of known ID's
Yk <- readRDS(paste0(inDir, "SDH_count_SAM_CT"))
Cams <- colnames(Yk)
Deer <- row.names(Yk)

# Unmarked detection history matrix
Yu <- read.csv(paste0(inDir, "Yu_SAMBAR_CT.csv"), header=T) %>%
  select(-X) # drop column of row (camera site) names

# Specify camera locations as spatial point process, in km scale
camlocs <- camlocs %>% dplyr::select(X,Y)
camlocs <- camlocs/1000 
camlocs <- as.ppp(camlocs, W=cam.region)
locs <- coords(camlocs)

# Set model foundations
M <- 150
mmax <- 50
lam0.prior <- list("uniform", 0, 5)
sigma.prior <- list("gamma", 8.6, 30)
delta <- c(0.17, 0.02, 2) 
xlim <- range(locs[,1])
ylim <- range(locs[,2])
ni <- 120000
nb <- 10000
nc <- 6
seed <-90
set.seed(seed)

# Run model
startDate <- Sys.Date()
logfile <- paste0("mylog_CT_SAMBAR_", startDate, ".txt")
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
## spatial climates
library(data.table)
library(sf)
library(RPostgreSQL)
library(dplyr)
library(foreach)
library(rmapshaper)
library(tictoc)
library(rasterVis)
library(raster)
library(ccissdev)
library(RPostgreSQL)
library(sf)
library(pool)

##some setup
con <- dbPool(
  drv = RPostgres::Postgres(),
  dbname = Sys.getenv("BCGOV_DB"),
  host = Sys.getenv("BCGOV_HOST"),
  port = 5432, 
  user = Sys.getenv("BCGOV_USR"),
  password = Sys.getenv("BCGOV_PWD")
)
sppDb <- dbPool(
  drv = RPostgres::Postgres(),
  dbname = "spp_feas",
  host = Sys.getenv("BCGOV_HOST"),
  port = 5432,
  user = Sys.getenv("BCGOV_USR"),
  password = Sys.getenv("BCGOV_PWD")
)

X <- raster("BC_Raster.tif")
X <- raster::setValues(X,NA)
outline <- st_read(con,query = "select * from bc_outline")
S1 <- setDT(dbGetQuery(sppDb,"select bgc,ss_nospace,spp,newfeas from feasorig"))
setnames(S1,c("BGC","SS_NoSpace","Spp","Feasible"))

##adapted feasibility function
ccissMap <- function(SSPred,suit,spp_select){
  ### generate raw feasibility ratios
  
  suit <- suit[Spp == spp_select,.(BGC,SS_NoSpace,Spp,Feasible)]
  suit <- unique(suit)
  suit <- na.omit(suit)
  SSPred <- SSPred[,.(SiteRef,FuturePeriod,BGC,SS_NoSpace,SS.pred,SSprob)]
  Site_BGC <- unique(SSPred[,.(SiteRef,BGC)])
  SSPred <- na.omit(SSPred)
  setkey(SSPred,SS.pred)
  setkey(suit,SS_NoSpace)
  suitMerge <- suit[SSPred, allow.cartesian = T]
  suitMerge <- na.omit(suitMerge)
  setnames(suitMerge, old = c("SS_NoSpace", "i.SS_NoSpace"), new = c("SS.pred", "SS_NoSpace"))
  suitVotes <- data.table::dcast(suitMerge, SiteRef + Spp + FuturePeriod + SS_NoSpace ~ Feasible, 
                                 value.var = "SSprob", fun.aggregate = sum)
  # Fill with 0 if columns does not exist, encountered the error at SiteRef 3104856 
  set(suitVotes, j = as.character(1:5)[!as.character(1:5) %in% names(suitVotes)], value = 0)
  suitVotes[,VoteSum := `1`+`2`+`3`+`4`+`5`]
  suitVotes[,X := 1 - VoteSum]
  suitVotes[,VoteSum := NULL]
  suitVotes[,X := X + `5` + `4`]
  suitVotes[,`:=`(`5` = NULL, `4` = NULL)]
  setkey(suitVotes, SS_NoSpace, Spp)
  setkey(suit, SS_NoSpace, Spp)
  suitVotes[suit, Curr := i.Feasible]
  suitVotes[is.na(Curr), Curr := 5]
  setorder(suitVotes,SiteRef,SS_NoSpace,Spp,FuturePeriod)
  suitVotes[Curr > 3.5, Curr := 4]
  colNms <- c("1","2","3","X")
  suitVotes <- suitVotes[,lapply(.SD, sum),.SDcols = colNms, 
                         by = .(SiteRef,FuturePeriod, SS_NoSpace,Spp,Curr)]
  suitVotes[,NewSuit := `1`+(`2`*2)+(`3`*3)+(X*5)]
  suitRes <- suitVotes[,.(Curr = mean(Curr),NewSuit = mean(NewSuit)), by = .(SiteRef)]
  return(suitRes)
}


################### straight predicted feasibility maps #####################
feasCols <- data.table(Feas = c(1,2,3,4,5),Col = c("limegreen", "deepskyblue", "gold", "grey","grey"))
area <- st_read("~/Downloads/ReburnBC_StudySite1")
area <- st_zm(area)
X <- raster(area, resolution = 400)
values(X) <- 1:ncell(X)
hexPts <- st_read("~/BC_HexGrid/BC_HexPoints400m.gpkg")
hexPts <- st_crop(hexPts, st_bbox(X))
ids <- raster::extract(X, hexPts)
cw_table <- data.table(SiteNo = hexPts$siteno,RastID = ids)
cw_table <- unique(cw_table, by = "RastID")

##gcm and rcp weight
gcm_weight <- data.table(gcm = c("ACCESS-ESM1-5", "BCC-CSM2-MR", "CanESM5", "CNRM-ESM2-1", "EC-Earth3",
                                 "GFDL-ESM4", "GISS-E2-1-G", "INM-CM5-0", "IPSL-CM6A-LR", "MIROC6",
                                 "MPI-ESM1-2-HR", "MRI-ESM2-0", "UKESM1-0-LL"),
                         weight = c(0,0,0,1,1,0,0,0,0,0,0,1,0))
#weight = c(1,1,0,0,1,1,1,0,1,1,1,1,0))
rcp_weight <- data.table(rcp = c("ssp126","ssp245","ssp370","ssp585"),
                         weight = c(0,1,1,0))

all_weight <- as.data.table(expand.grid(gcm = gcm_weight$gcm,rcp = rcp_weight$rcp))
all_weight[gcm_weight,wgcm := i.weight, on = "gcm"]
all_weight[rcp_weight,wrcp := i.weight, on = "rcp"]
all_weight[,weight := wgcm*wrcp]
modWeights <- all_weight

dat <- dbGetCCISS(con, cw_table$SiteNo, avg = F, modWeights = all_weight)
dat[,SiteRef := as.integer(SiteRef)]
timeperiods <- c("1961","2041","2081")
edaPos <- c("B2", "C4", "D6")
species <- c("Fd", "Sx", "Pl", "Bl", "Py", "Lw", "At")

for(tp in timeperiods){
  datTP <- dat[FuturePeriod == tp,]
  for(eda in edaPos){
    edaTemp <- data.table::copy(E1)
    edaTemp <- edaTemp[is.na(SpecialCode),]
    
    edaTemp[,HasPos := if(any(Edatopic == eda)) T else F, by = .(SS_NoSpace)]
    edaZonal <- edaTemp[(HasPos),]
    edaZonal[,HasPos := NULL]
    ##edatopic overlap
    SSPreds <- edatopicOverlap(datTP,edaZonal,E1_Phase,onlyRegular = TRUE) ##takes about 30 seconds
    ##loop through species
    for(spp in species){
      sppFeas <- ccissMap(SSPreds,S1,spp) ##~ 15 seconds
      sppFeas <- unique(sppFeas,by = "SiteRef")
      sppFeas[,SiteRef := as.integer(SiteRef)]
      sppFeas[cw_table, RastID := i.RastID, on = c(SiteRef = "SiteNo")]
      sppFeas <- sppFeas[,.(RastID,NewSuit)]
      X <- raster::setValues(X,NA)
      X[sppFeas$RastID] <- sppFeas$NewSuit
      writeRaster(X, filename = paste("./ReburnBC_Maps/CCISSFeas",tp,eda,spp,".tif",sep = "_"),format = "GTiff", overwrite = T)
      # X2 <- ratify(X)
      # rat <- as.data.table(levels(X2)[[1]])
      # rat[feasCols,`:=`(col = i.Col), on = c(ID = "Feas")]
      # 
      # pdf(file=paste("./FeasibilityMaps/Feasibility",timeperiods,spp,".pdf",sep = "_"), width=6.5, height=7, pointsize=10)
      # plot(X2,col = rat$col,legend = FALSE,axes = FALSE, box = FALSE, main = paste0(spp," (",timeperiods,")"))
      # plot(outline, col = NA, add = T)
      # dev.off()
    }
  }
}


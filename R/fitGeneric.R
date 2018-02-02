# fit data for nlsGFP
source("frapr-pkg.R")
library(gdata)


#' Calculate and plot the variance in F.init and C.init initialisation
#' 
#' @param p params list used in the simulation, with calculate F.init and C.init
plotInitVariance = function(p){
    total = sum(p$F.init[p$shapeMat] + p$C.init[p$shapeMat])
    
    ss = steadyState(total, p$shapeMat, p$k.on, p$k.off)

    F.prop = p$F.init / ss$F
    C.prop = p$C.init / ss$C
    
    F.prop[!p$shapeMat] = 0
    C.prop[!p$shapeMat] = 0
    if(any(is.nan(C.prop)))
        C.prop[is.nan(C.prop)] = 0
    
    F.depth = mean(F.prop)
    C.depth = mean(C.prop)
    
    C.exists = any(ss$C!=0)
    
    if(C.exists){    
        #par(mfrow=c(1,2))
        plotFL(F.prop, bleachMask = p$bleachMask)
        plotFL(C.prop, bleachMask = p$bleachMask)
    } else {
        plotFL(F.prop, bleachMask = p$bleachMask)
    }
    
    list(F.prop=F.prop, C.prop=C.prop, 
         F.depth=F.depth, C.depth=C.depth)
}

#' Calculate the distribution of fluorescence from centre and edge of bleaching area
#' Takes the output of plotInitVariance()
#' @param fl.var output of plotInitVariance()
#' @param p parameters
calculateFSpread = function(fl.var, p){
    f = fl.var$F.prop
    s = p$shapeMat

    dx = p$L.x / p$grid.x
    dy = p$L.y / p$grid.y
        
    # calculate centre
    points = which(s, arr.ind=TRUE)
    bleach = which(p$bleachMask, arr.ind=TRUE)
    bleach[,1] = bleach[,1] * dx
    bleach[,2] = bleach[,2] * dy
    centroid = colSums(bleach) / nrow(bleach)
    
    # calculate the value relative to centroid
    d.centroid = matrix(ncol=4, nrow=nrow(points))
    colnames(d.centroid) = c("dist", "val", "x", "y")
    for(i in 1:nrow(points)){
        p1 = points[i,1] * dx
        p2 = points[i,2] * dy
        d.centroid[i, "dist"] = sqrt(sum((c(p1,p2)-centroid)^2))
        d.centroid[i, "val"] = f[points[i,1], points[i,2]]
        d.centroid[i, 3:4] = points[i,]
    }
    
    list(centroid=centroid, d.centroid=d.centroid)
}

#' A simple wrapper to fit multiple replicates
#' @param replicates a vector of replicates numbers
fitGenericReplicates = function(replicates, ...){
    lapply(replicates, fitGeneric, ...)
}

#' A generic function to fit a FRAP curve
#' 
#' @param replicate the number of the replicate, e.g. 1
#' @param experiment the directory of the experiment
#' @param frame.pat file name pattern used to get frame
#' @param bleach.frames frames where bleaching occurs
#' @param time.inc base time increment
#' @param roi.fg names of ROIs that are foreground, i.e. bleaching spots
#' @param real.len.x Actual length in meters over the x axis
#' @param real.len.y Actual length in meters over the y axis
#' @param roi.exclude Name of the ROI that should be used to exclude an area, e.g. a nucleolus
#' @param roi.total Name of the ROI that contains the total area of the nucleus
#' @param time.inc2.start If not NULL, the start frame index of the second time increment
#' @param time.inc2.value The second time increment in seconds
#' @param base.path Path where the experiments are (normally automatically formed by 'experiment' and 'replicate')
#' @param algorithm Algorithm to use in fitting
#' @param scale.by scale-down factor, higher values indicate larger scaling down
#' @param D.range the range of Diffusion values to consider
#' @param Free.range the range of Free molecules to consider
#' @param res.time.range the range of residence times to consider
#' @param save.to.file if TRUE, will save to a file in data/, if FALSE, will return the fitted values in a list
#' @param use.old.init.rates if to use the old initRatesOld function that initialises the rates using the post-bleach image
fitGeneric = function(replicate = NULL, 
                      experiment = NULL, 
                      frame.pat = "img_t%F.tif", 
                      bleach.frames = 9, 
                      time.inc = 0.098,  
                      roi.fg=c("B", "B"),
                      real.len.x = 64.58 * 1e-6,
                      real.len.y = 37.79 * 1e-6,
                      roi.exclude = "",
                      roi.total = "T", 
                      roi.bg.pat = "BG",
                      roi.file="rois.zip",
                      time.inc2.start=310,
                      time.inc2.value=1,
                      base.path = paste("../images/", experiment, "/", replicate, sep=""),
                      algorithm=c("full", "diffusion", "mcnally-full", "mcnally-diffusion"),
                      scale.by=10,
                      D.range=seq(0.5, 20, 0.5) * 1e-12,
                      Free.range=c(0.05, 0.95, 0.05),
                      res.time.range=c(0.5, 1, 2, 4, 6, 8, 10, 15, 20, 25, 30, 40, 50, 60, 70, 80),
                      res.time2.range=NULL,
                      save.to.file=NA,
                      use.old.init.rates=FALSE,
                      post.avg=3,
                      roi.total.file=NULL,
                      first.ROI.only=TRUE,
                      bleach.norm.sd=1.7*1e-6,
                      bleaching.init.version=1,
                      use.bleaching.init=!is.null(bleaching.init.version),
                      bleach.norm.sd.range=bleach.norm.sd,
                      res.time.c2.range=NULL,
                      prop.c2.range=NULL,
                      use.uniform.init=FALSE
                      ){
    
    cat("Using scaling factor: ", scale.by, "\n")
    
    algorithm = match.arg(algorithm)
    
    # load the images
    img = makeImageProps(experiment, replicate, base.path, frame.pat, bleach.frames, time.inc, real.len.x, real.len.y, 
                         roi.total=roi.total, roi.bg.pat=roi.bg.pat, roi.file=roi.file, roi.exclude = roi.exclude,
                         post.avg=post.avg, roi.total.file = roi.total.file, bleach.norm.sd=bleach.norm.sd,
                         bleach.norm.sd.range = bleach.norm.sd.range)
    
    ####################################################################
    # modify the time stamp of the last 100 frames for longer times
    if(!is.null(time.inc2.start)){
        img$t[time.inc2.start:length(img$t)] = img$t[time.inc2.start-1] + time.inc2.value*(1:(length(img$t)-(time.inc2.start-1)))
    }
    
    # make sure nothing is outside the total area		
    #img$rois$B = img$rois$B & img$rois$T
    #img$rois$UB = img$rois$UB & img$rois$T
    
    img = imageStats(img, roi.fg)
    
    ### Fitting with mcNally algorithm!!!
    if(algorithm == "mcnally-full" | algorithm == "mcnally-diffusion"){
        p1 = initFromImage(img, scale.by=1)
        sim = list()
        sim$t = img$t
        sim$params = p1
        sim$curve = img$stats$rec[,1]
        sim$params$bleach.depth = img$stats$rec[img$bleach.frames+1,1]
        sim$params$bleach.time = img$t[img$bleach.frames]
        
        k.off.all = 1 / res.time.range
        
        if(algorithm == "mcnally-full"){
            frap = fitMcNally2008_full(sim, Df.all = D.range,
                                       Free.all=Free.range, k.off.all=k.off.all)
        } else {
            frap = fitMcNally2008_diffusionOnly(sim, Df.all=D.range)
        }
        
        if(!is.na(save.to.file)){
            save(img, p1, frap, sim, file=paste("data/", save.to.file, "-", algorithm, "-replicate_", replicate, ".RData", sep=""))
        }
        return(list(p=p1, img=img, frap=frap, sim=sim))
    } else {
        # initialise the simulation parameters
        if(algorithm == "diffusion"){
            p = initFromImage(img, scale.by=scale.by, "diffusion")
        } else {
            p = initFromImage(img, scale.by=scale.by)
        }
        
        # evaluate the parameter range
        if(!is.null(res.time.c2.range)){
            e = evalParamsC2(img, p, res.time.range, res.time.c2.range, prop.c2.range, 
                             Free.range, D.range, first.ROI.only=first.ROI.only,
                             use.bleaching.init=use.bleaching.init,
                             bleaching.init.version=bleaching.init.version,
                             use.uniform.init=use.uniform.init)
        } else if(is.null(res.time2.range)){
            e = evalParams(img, p, res.time.range, Free.range, D.range, use.old.init.rates=use.old.init.rates, 
                           first.ROI.only=first.ROI.only, use.bleaching.init=use.bleaching.init,
                           bleaching.init.version=bleaching.init.version,
                           use.uniform.init=use.uniform.init)
        } else {
            e = evalParams2Res(img, p, res.time.range, res.time2.range, names(img$rois)[1], 
                               Free.range, D.range, use.old.init.rates=use.old.init.rates,
                               use.bleaching.init=use.bleaching.init,
                               bleaching.init.version=bleaching.init.version)
        }
        
        if(!is.na(save.to.file)){
            save(img, p, e, file=paste("data/", save.to.file, "-", algorithm, "-replicate_", replicate, ".RData", sep=""))
        }
        return(list(img=img, p=p, e=e))
    }
}
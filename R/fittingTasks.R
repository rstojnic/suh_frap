# fitting tasks

source("fitGeneric.R")

D.range.std = c(0.1, 0.5, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20) * 1e-12
res.time.range.std = c(0.1, 0.25, 0.5, 0.75, 1, 
                       2, 4, 6, 8, 10, 20, 
                       30, 40)

res.time.range.long = c(0.1, 0.25, 0.5, 0.75, 1, 2, 4, 6, 8, 10, 15, 20, 25, 30, 
                        40, 50, 60, 70, 80, 100, 120, 140, 160, 180, 200)

processParamsStd = function(p){
    if(is.null(p$D.range))
        p$D.range = D.range.std
    if(is.null(p$res.time.range))
        p$res.time.range = res.time.range.long
    if(is.null(p$Free.range))
        p$Free.range = seq(0.05, 0.95, 0.05)
    if(is.null(p$use.old.init.rates))
        p$use.old.init.rates = TRUE
    if(is.null(p$scale.by))
        p$scale.by = 10
    if(is.null(p$roi.file))
        p$roi.file = "rois.zip"
    if(is.null(p$roi.fg))
        p$roi.fg = c("B", "B")
    if(is.null(p$roi.total.file))
        p$roi.total.file = NULL
    if(is.null(p$first.ROI.only))
        p$first.ROI.only = TRUE
    if(is.null(p$bleaching.init.version))
        p$bleaching.init.version = NULL
    if(is.null(p$bleach.norm.sd.range))
        p$bleach.norm.sd.range = seq(1.2, 2.2, 0.1)*1e-6
    if(is.null(p$res.time.c2.range))
        p$res.time.c2.range = NULL
    if(is.null(p$prop.c2.range))
        p$prop.c2.range = NULL
    if(is.null(p$roi.exclude))
        p$roi.exclude = ""
    if(is.null(p$use.uniform.init))
        p$use.uniform.init = FALSE
    
    p
}

processParams2Res = function(p){
    if(is.null(p$D.range))
        p$D.range = D.range.std
    if(is.null(p$Free.range))
        p$Free.range = seq(0.05, 0.95, 0.05)
    if(is.null(p$use.old.init.rates))
        p$use.old.init.rates = TRUE
    if(is.null(p$scale.by))
        p$scale.by = 10
    if(is.null(p$roi.file))
        p$roi.file = "rois.zip"
    if(is.null(p$roi.fg))
        p$roi.fg = c("B", "B")
    if(is.null(p$roi.total.file))
        p$roi.total.file = NULL
    if(is.null(p$first.ROI.only))
        p$first.ROI.only = TRUE
    if(is.null(p$bleaching.init.version))
        p$bleaching.init.version = NULL
    if(is.null(p$bleach.norm.sd.range))
        p$bleach.norm.sd.range = seq(1.2, 2.2, 0.1)*1e-6
    if(is.null(p$res.time.c2.range))
        p$res.time.c2.range = NULL
    if(is.null(p$prop.c2.range))
        p$prop.c2.range = NULL
    if(is.null(p$roi.exclude))
        p$roi.exclude = ""
    if(is.null(p$use.uniform.init))
        p$use.uniform.init = FALSE

    if(is.null(p$res.time.range))
        p$res.time.range = 10
    if(is.null(p$res.time2.range))
        p$res.time2.range = res.time.range.long
    
    p
}

fittingBaseGenFit = function(r,
                             p=list(), dataset.name, dataset.path){
    p = processParamsStd(p)
    if(is.null(p$save.base))
        p$save.base = dataset.name

    fitGeneric(replicate=r,
               experiment = dataset.path,
               time.inc = 0.115,
               time.inc2.start = 310,
               time.inc2.value = 2,
               real.len.x = 64.58 * 1e-6,
               real.len.y = 40.30 * 1e-6,
               algorithm= p$algorithm,
               D.range = p$D.range,
               res.time.range=p$res.time.range,
               Free.range=p$Free.range,
               use.old.init.rates = p$use.old.init.rates,
               scale.by = p$scale.by,
               save.to.file=p$save.base,
               roi.file=p$roi.file,
               roi.fg=p$roi.fg,
               roi.total.file=p$roi.total.file,
               first.ROI.only=p$first.ROI.only,
               bleaching.init.version=p$bleaching.init.version,
               bleach.norm.sd.range=p$bleach.norm.sd.range,
               res.time.c2.range=p$res.time.c2.range, 
               prop.c2.range=p$prop.c2.range,
               roi.exclude = p$roi.exclude,
               use.uniform.init=p$use.uniform.init
    )
}

fittingBaseGenFit2017 = function(r,
                             p=list(), dataset.name, dataset.path){
    p = processParamsStd(p)
    if(is.null(p$save.base))
        p$save.base = dataset.name
    
    fitGeneric(replicate=r,
               experiment = dataset.path,
               time.inc = 0.112,
               time.inc2.start = 310,
               time.inc2.value = 2,
               real.len.x = 64.58 * 1e-6,
               real.len.y = 37.79 * 1e-6,
               algorithm= p$algorithm,
               D.range = p$D.range,
               res.time.range=p$res.time.range,
               Free.range=p$Free.range,
               use.old.init.rates = p$use.old.init.rates,
               scale.by = p$scale.by,
               save.to.file=p$save.base,
               roi.file=p$roi.file,
               roi.fg=p$roi.fg,
               roi.total.file=p$roi.total.file,
               first.ROI.only=p$first.ROI.only,
               bleaching.init.version=p$bleaching.init.version,
               bleach.norm.sd.range=p$bleach.norm.sd.range,
               res.time.c2.range=p$res.time.c2.range, 
               prop.c2.range=p$prop.c2.range,
               roi.exclude = p$roi.exclude,
               use.uniform.init=p$use.uniform.init
    )
}

fittingBaseGenFit2Res = function(r,
                             p=list(), dataset.name, dataset.path){
    p = processParams2Res(p)
    if(is.null(p$save.base))
        p$save.base = dataset.name
    
    fitGeneric(replicate=r,
               experiment = dataset.path,
               time.inc = 0.115,
               time.inc2.start = 310,
               time.inc2.value = 2,
               real.len.x = 64.58 * 1e-6,
               real.len.y = 40.30 * 1e-6,
               algorithm= p$algorithm,
               D.range = p$D.range,
               res.time.range=p$res.time.range,
               res.time2.range=p$res.time2.range,
               Free.range=p$Free.range,
               use.old.init.rates = p$use.old.init.rates,
               scale.by = p$scale.by,
               save.to.file=p$save.base,
               roi.file=p$roi.file,
               roi.fg=p$roi.fg,
               roi.total.file=p$roi.total.file,
               first.ROI.only=p$first.ROI.only,
               bleaching.init.version=p$bleaching.init.version,
               bleach.norm.sd.range=p$bleach.norm.sd.range,
               res.time.c2.range=p$res.time.c2.range, 
               prop.c2.range=p$prop.c2.range,
               roi.exclude = p$roi.exclude,
               use.uniform.init=p$use.uniform.init
    )
}

fittingBase = list(
    "GFP"=list(
        "file"="20160912_pointFRAP GFP",
        "fit"=function(r, p){
            if(is.null(p$D.range))
                p$D.range = seq(2, 40, 2) * 1e-12
            if(is.null(p$save.base))
                p$save.base = "GFP-sep2016"
            if(is.null(p$algorithm))
                p$algorithm = "diffusion"
            if(is.null(p$use.old.init.rates))
                p$use.old.init.rates = TRUE
            if(is.null(p$scale.by))
                p$scale.by = 10
            if(is.null(p$roi.file))
                p$roi.file = "rois.zip"
            if(is.null(p$res.time.range))
                p$res.time.range = 0
            if(is.null(p$Free.range))
                p$Free.range = 1
            
            fitGeneric(replicate=r,
                       experiment = "20160912_pointFRAP GFP",
                       time.inc = 0.115,
                       time.inc2.start = 310,
                       time.inc2.value = 2,
                       real.len.x = 64.58 * 1e-6,
                       real.len.y = 40.30 * 1e-6,
                       algorithm=p$algorithm,
                       D.range =p$D.range,
                       res.time.range = p$res.time.range,
                       Free.range = p$Free.range,
                       use.old.init.rates = p$use.old.init.rates,
                       scale.by = p$scale.by,
                       save.to.file=p$save.base
            )
        },
        "reps"=1:21
    ),
    "SuH_DNAmutant"=list(
        "file"="20150818_PointFRAP R2ggH 1copy_ROI",
        "fit"=function(r,
                       p=list()
                       ){
            p = processParamsStd(p)
            if(is.null(p$save.base))
                p$save.base = "SuH_DNAmutant"
                        
            fitGeneric(replicate=r,
                       experiment = "20150818_PointFRAP R2ggH 1copy_ROI",
                       time.inc = 0.098,
                       time.inc2.start = 310,
                       time.inc2.value = 1,
                       real.len.x = 64.58 * 1e-6,
                       real.len.y = 40.30 * 1e-6,
                       algorithm = p$algorithm,
                       D.range = p$D.range,
                       res.time.range = p$res.time.range,
                       Free.range= p$Free.range,
                       use.old.init.rates = p$use.old.init.rates,
                       scale.by = p$scale.by,
                       save.to.file= p$save.base,
                       roi.file=p$roi.file,
                       roi.fg=p$roi.fg,
                       roi.total.file=p$roi.total.file,
                       first.ROI.only=p$first.ROI.only,
                       bleaching.init.version=p$bleaching.init.version,
                       bleach.norm.sd.range=p$bleach.norm.sd.range,
                       res.time.c2.range=p$res.time.c2.range, 
                       prop.c2.range=p$prop.c2.range,
                       roi.exclude = p$roi.exclude,
                       use.uniform.init=p$use.uniform.init
            )
        },
        "reps"=1:17
    ),
    
    "SuH_WT"=list(
        "file"="20150922_pointFRAP_Su(H)_WT",
        "fit"=function(r,
                       p = list()
                       ){
            p = processParamsStd(p)
            if(is.null(p$save.base))
                p$save.base = "SuH_WT"
            
            fitGeneric(replicate=r,
                       experiment = "20150922_pointFRAP_Su(H)_WT",
                       time.inc = 0.098,
                       time.inc2.start = 310,
                       time.inc2.value = 1,
                       real.len.x = 48.44 * 1e-6,
                       real.len.y = 30.23 * 1e-6,
                       algorithm= p$algorithm,
                       D.range = p$D.range,
                       res.time.range= p$res.time.range,
                       Free.range= p$Free.range,
                       use.old.init.rates = p$use.old.init.rates,
                       scale.by = p$scale.by,
                       save.to.file= p$save.base,
                       roi.file=p$roi.file,
                       roi.fg=p$roi.fg,
                       roi.total.file=p$roi.total.file,
                       first.ROI.only=p$first.ROI.only,
                       bleaching.init.version=p$bleaching.init.version,
                       bleach.norm.sd.range=p$bleach.norm.sd.range,
                       res.time.c2.range=p$res.time.c2.range, 
                       prop.c2.range=p$prop.c2.range,
                       roi.exclude = p$roi.exclude,
                       use.uniform.init=p$use.uniform.init
            )
        },
        "reps"=1:19
    ),
    "Smr"=list(
        "file"="20150916_PointFRAP_Smr-YFP",
        "fit"=function(r,
                       p=list()){
            p = processParamsStd(p)
            if(is.null(p$save.base))
                p$save.base = "Smr"
            
            fitGeneric(replicate=r,
                       experiment = "20150916_PointFRAP_Smr-YFP",
                       time.inc = 0.098,
                       time.inc2.start = 310,
                       time.inc2.value = 1,
                       real.len.x = 48.44 * 1e-6,
                       real.len.y = 30.23 * 1e-6,
                       algorithm= p$algorithm,
                       D.range = p$D.range,
                       res.time.range=p$res.time.range,
                       Free.range=p$Free.range,
                       use.old.init.rates = p$use.old.init.rates,
                       scale.by = p$scale.by,
                       save.to.file=p$save.base,
                       roi.file=p$roi.file,
                       roi.fg=p$roi.fg,
                       roi.total.file=p$roi.total.file,
                       first.ROI.only=p$first.ROI.only,
                       bleaching.init.version=p$bleaching.init.version,
                       bleach.norm.sd.range=p$bleach.norm.sd.range,
                       res.time.c2.range=p$res.time.c2.range, 
                       prop.c2.range=p$prop.c2.range,
                       roi.exclude = p$roi.exclude,
                       use.uniform.init=p$use.uniform.init
            )
        },
        "reps"=2:17
    ),
    "H"=list(
        "file"="20150924_pointFRAP_H-GFP",
        "fit"=function(r,
                       p=list()){
            p = processParamsStd(p)
            if(is.null(p$save.base))
                p$save.base = "H"
            
            fitGeneric(replicate=r,
                       experiment = "20150924_pointFRAP_H-GFP",
                       time.inc = 0.098,
                       time.inc2.start = 310,
                       time.inc2.value = 1,
                       real.len.x = 48.44 * 1e-6,
                       real.len.y = 30.23 * 1e-6,
                       algorithm= p$algorithm,
                       D.range = p$D.range,
                       res.time.range=p$res.time.range,
                       Free.range=p$Free.range,
                       use.old.init.rates = p$use.old.init.rates,
                       scale.by = p$scale.by,
                       save.to.file=p$save.base,
                       roi.file=p$roi.file,
                       roi.fg=p$roi.fg,
                       roi.total.file=p$roi.total.file,
                       first.ROI.only=p$first.ROI.only,
                       bleaching.init.version=p$bleaching.init.version,
                       bleach.norm.sd.range=p$bleach.norm.sd.range,
                       res.time.c2.range=p$res.time.c2.range, 
                       prop.c2.range=p$prop.c2.range,
                       roi.exclude = p$roi.exclude,
                       use.uniform.init=p$use.uniform.init
            )
        },
        "reps"=1:21
    ),
    "SuH_locus_necd"=list(
        "file"="20160329_Band FRAP locus tagging NDECD ROIs",
        "fit"=function(r,
                       p=list()){
            p = processParamsStd(p)
            if(is.null(p$save.base))
                p$save.base = "SuH_locus_ndecd"
            
            fitGeneric(replicate=r,
                       experiment = "20160329_Band FRAP locus tagging NDECD ROIs",
                       time.inc = 0.115,
                       time.inc2.start = 310,
                       time.inc2.value = 2,
                       real.len.x = 64.58 * 1e-6,
                       real.len.y = 40.30 * 1e-6,
                       algorithm= p$algorithm,
                       D.range = p$D.range,
                       res.time.range=p$res.time.range,
                       Free.range=p$Free.range,
                       use.old.init.rates = p$use.old.init.rates,
                       scale.by = p$scale.by,
                       save.to.file=p$save.base,
                       roi.file=p$roi.file,
                       roi.fg=p$roi.fg,
                       roi.total.file=p$roi.total.file,
                       first.ROI.only=p$first.ROI.only,
                       bleaching.init.version=p$bleaching.init.version,
                       bleach.norm.sd.range=p$bleach.norm.sd.range,
                       res.time.c2.range=p$res.time.c2.range, 
                       prop.c2.range=p$prop.c2.range,
                       roi.exclude = p$roi.exclude,
                       use.uniform.init=p$use.uniform.init
            )
        },
        "reps"=1:7
    ),
    "SuH_locus_necd_scale6"=list(
        "file"="20160329_Band FRAP locus tagging NDECD ROIs",
        "fit"=function(r,
                       p=list()){
            if(is.null(p$res.time.range))
                p$res.time.range = res.time.range.std
            p = processParamsStd(p)
            if(is.null(p$save.base))
                p$save.base = "SuH_locus_ndecd_scale6"
            
            p$scale.by = 6
            
            fitGeneric(replicate=r,
                       experiment = "20160329_Band FRAP locus tagging NDECD ROIs",
                       time.inc = 0.115,
                       time.inc2.start = 310,
                       time.inc2.value = 2,
                       real.len.x = 64.58 * 1e-6,
                       real.len.y = 40.30 * 1e-6,
                       algorithm= p$algorithm,
                       D.range = p$D.range,
                       res.time.range=p$res.time.range,
                       Free.range=p$Free.range,
                       use.old.init.rates = p$use.old.init.rates,
                       scale.by = p$scale.by,
                       save.to.file=p$save.base,
                       roi.file=p$roi.file,
                       roi.fg=p$roi.fg,
                       roi.total.file=p$roi.total.file,
                       first.ROI.only=p$first.ROI.only,
                       bleaching.init.version=p$bleaching.init.version,
                       bleach.norm.sd.range=p$bleach.norm.sd.range,
                       res.time.c2.range=p$res.time.c2.range, 
                       prop.c2.range=p$prop.c2.range,
                       roi.exclude = p$roi.exclude,
                       use.uniform.init=p$use.uniform.init
            )
        },
        "reps"=1:7
    ),
    "SuH_locus_necd_2res"=list(
        "file"="20160329_Band FRAP locus tagging NDECD ROIs",
        "fit"=function(r,
                       p=list()){
            p = processParams2Res(p)
            if(is.null(p$save.base))
                p$save.base = "SuH_locus_ndecd_2res"
            
            fitGeneric(replicate=r,
                       experiment = "20160329_Band FRAP locus tagging NDECD ROIs",
                       time.inc = 0.115,
                       time.inc2.start = 310,
                       time.inc2.value = 2,
                       real.len.x = 64.58 * 1e-6,
                       real.len.y = 40.30 * 1e-6,
                       algorithm= p$algorithm,
                       D.range = p$D.range,
                       res.time.range=p$res.time.range,
                       res.time2.range=p$res.time2.range,
                       Free.range=p$Free.range,
                       use.old.init.rates = p$use.old.init.rates,
                       scale.by = p$scale.by,
                       save.to.file=p$save.base,
                       roi.file=p$roi.file,
                       roi.fg=p$roi.fg,
                       roi.total.file=p$roi.total.file,
                       first.ROI.only=p$first.ROI.only,
                       bleaching.init.version=p$bleaching.init.version,
                       bleach.norm.sd.range=p$bleach.norm.sd.range,
                       res.time.c2.range=p$res.time.c2.range, 
                       prop.c2.range=p$prop.c2.range,
                       roi.exclude = p$roi.exclude,
                       use.uniform.init=p$use.uniform.init
            )
        },
        "reps"=1:7
    ),
    "SuH_locus_control"=list(
        "file"="20160407_Band FRAP locus tagging control ROIs",
        "fit"=function(r,
                       p=list()){
            p = processParamsStd(p)
            if(is.null(p$save.base))
                p$save.base = "SuH_locus_control"
            
            fitGeneric(replicate=r,
                       experiment = "20160407_Band FRAP locus tagging control ROIs",
                       time.inc = 0.115,
                       time.inc2.start = 310,
                       time.inc2.value = 2,
                       real.len.x = 64.58 * 1e-6,
                       real.len.y = 40.30 * 1e-6,
                       algorithm= p$algorithm,
                       D.range = p$D.range,
                       res.time.range=p$res.time.range,
                       Free.range=p$Free.range,
                       use.old.init.rates = p$use.old.init.rates,
                       scale.by = p$scale.by,
                       save.to.file=p$save.base,
                       roi.file=p$roi.file,
                       roi.fg=p$roi.fg,
                       roi.total.file=p$roi.total.file,
                       first.ROI.only=p$first.ROI.only,
                       bleaching.init.version=p$bleaching.init.version,
                       bleach.norm.sd.range=p$bleach.norm.sd.range,
                       res.time.c2.range=p$res.time.c2.range, 
                       prop.c2.range=p$prop.c2.range,
                       roi.exclude = p$roi.exclude,
                       use.uniform.init=p$use.uniform.init
            )
        },
        "reps"=1:23,
        "clean.reps"=c(1, 2, 6:9, 13, 14, 18, 20:23)
    ),
    "gro1"=list(
        "file"="20151027_FRAP Gro",
        "fit"=function(r,
                       p=list()){
            p = processParamsStd(p)
            if(is.null(p$save.base))
                p$save.base = "gro1"
            
            fitGeneric(replicate=r,
                       experiment = "20151027_FRAP Gro",
                       time.inc = 0.098,
                       time.inc2.start = 310,
                       time.inc2.value = 2,
                       real.len.x = 64.58 * 1e-6,
                       real.len.y = 40.30 * 1e-6,
                       algorithm= p$algorithm,
                       D.range = p$D.range,
                       res.time.range=p$res.time.range,
                       Free.range=p$Free.range,
                       use.old.init.rates = p$use.old.init.rates,
                       scale.by = p$scale.by,
                       save.to.file=p$save.base,
                       frame.pat = "img_%F.tiff",
                       bleach.frames = 10,
                       roi.file=p$roi.file,
                       roi.fg=p$roi.fg,
                       roi.total.file=p$roi.total.file,
                       first.ROI.only=p$first.ROI.only,
                       bleaching.init.version=p$bleaching.init.version,
                       bleach.norm.sd.range=p$bleach.norm.sd.range,
                       res.time.c2.range=p$res.time.c2.range, 
                       prop.c2.range=p$prop.c2.range,
                       roi.exclude = p$roi.exclude,
                       use.uniform.init=p$use.uniform.init
            )
        },
        "reps"=1:15
    ),
    "gro2"=list(
        "file"="20151104_FRAP Gro",
        "fit"=function(r,
                       p=list()){
            p = processParamsStd(p)
            if(is.null(p$save.base))
                p$save.base = "gro2"
            
            fitGeneric(replicate=r,
                       experiment = "20151104_FRAP Gro",
                       time.inc = 0.098,
                       time.inc2.start = 310,
                       time.inc2.value = 2,
                       real.len.x = 64.58 * 1e-6,
                       real.len.y = 40.30 * 1e-6,
                       algorithm= p$algorithm,
                       D.range = p$D.range,
                       res.time.range=p$res.time.range,
                       Free.range=p$Free.range,
                       use.old.init.rates = p$use.old.init.rates,
                       scale.by = p$scale.by,
                       save.to.file=p$save.base,
                       frame.pat = "img_%F.tiff",
                       bleach.frames = 10,
                       roi.file=p$roi.file,
                       roi.fg=p$roi.fg,
                       roi.total.file=p$roi.total.file,
                       first.ROI.only=p$first.ROI.only,
                       bleaching.init.version=p$bleaching.init.version,
                       bleach.norm.sd.range=p$bleach.norm.sd.range,
                       res.time.c2.range=p$res.time.c2.range, 
                       prop.c2.range=p$prop.c2.range,
                       roi.exclude = p$roi.exclude,
                       use.uniform.init=p$use.uniform.init
            )
        },
        "reps"=1:15
    ),
    "fkh1"=list(
        "file"="20151028_FRAP_Fkh",
        "fit"=function(r,
                       p=list()){
            p = processParamsStd(p)
            if(is.null(p$save.base))
                p$save.base = "fkh1"
            
            fitGeneric(replicate=r,
                       experiment = "20151028_FRAP_Fkh",
                       time.inc = 0.098,
                       time.inc2.start = 310,
                       time.inc2.value = 2,
                       real.len.x = 64.58 * 1e-6,
                       real.len.y = 40.30 * 1e-6,
                       algorithm= p$algorithm,
                       D.range = p$D.range,
                       res.time.range=p$res.time.range,
                       Free.range=p$Free.range,
                       use.old.init.rates = p$use.old.init.rates,
                       scale.by = p$scale.by,
                       save.to.file=p$save.base,
                       frame.pat = "img_%F.tiff",
                       bleach.frames = 10,
                       roi.file=p$roi.file,
                       roi.fg=p$roi.fg,
                       roi.total.file=p$roi.total.file,
                       first.ROI.only=p$first.ROI.only,
                       bleaching.init.version=p$bleaching.init.version,
                       bleach.norm.sd.range=p$bleach.norm.sd.range,
                       res.time.c2.range=p$res.time.c2.range, 
                       prop.c2.range=p$prop.c2.range,
                       roi.exclude = p$roi.exclude,
                       use.uniform.init=p$use.uniform.init
            )
        },
        "reps"=1:16
    ),
    "fkh2"=list(
        "file"="20151103_FRAP Fkh",
        "fit"=function(r,
                       p=list()){
            p = processParamsStd(p)
            if(is.null(p$save.base))
                p$save.base = "fkh2"
            
            fitGeneric(replicate=r,
                       experiment = "20151103_FRAP Fkh",
                       time.inc = 0.098,
                       time.inc2.start = 310,
                       time.inc2.value = 2,
                       real.len.x = 64.58 * 1e-6,
                       real.len.y = 40.30 * 1e-6,
                       algorithm= p$algorithm,
                       D.range = p$D.range,
                       res.time.range=p$res.time.range,
                       Free.range=p$Free.range,
                       use.old.init.rates = p$use.old.init.rates,
                       scale.by = p$scale.by,
                       save.to.file=p$save.base,
                       frame.pat = "img_%F.tiff",
                       bleach.frames = 10,
                       roi.file=p$roi.file,
                       roi.fg=p$roi.fg,
                       roi.total.file=p$roi.total.file,
                       first.ROI.only=p$first.ROI.only,
                       bleaching.init.version=p$bleaching.init.version,
                       bleach.norm.sd.range=p$bleach.norm.sd.range,
                       res.time.c2.range=p$res.time.c2.range, 
                       prop.c2.range=p$prop.c2.range,
                       roi.exclude = p$roi.exclude,
                       use.uniform.init=p$use.uniform.init
            )
        },
        "reps"=1:15
    ),
    "SuH_necd_c1" = list(
        "file" = "20160901_SuH WT FRAP_ROIs",
        "fit" = function(r, p=list()){
            fittingBaseGenFit(r, p, "SuH_necd_c1", "20160901_SuH WT FRAP_ROIs")
        },
        "reps"=1:6
    ),
    "SuH_necd_c2" = list(
        "file" = "20160902_FRAP band SuH NDECD_ROIs",
        "fit" = function(r, p=list()){
            fittingBaseGenFit(r, p, "SuH_necd_c2", "20160902_FRAP band SuH NDECD_ROIs")
        },
        "reps"=1:8
    ),
    "SuH_necd_c1_2res" = list(
        "file" = "20160901_SuH WT FRAP_ROIs",
        "fit" = function(r, p=list()){
            fittingBaseGenFit2Res(r, p, "SuH_necd_c1_2res", "20160901_SuH WT FRAP_ROIs")
        },
        "reps"=1:6
    ),
    "SuH_necd_c2_2res" = list(
        "file" = "20160902_FRAP band SuH NDECD_ROIs",
        "fit" = function(r, p=list()){
            fittingBaseGenFit2Res(r, p, "SuH_necd_c2_2res", "20160902_FRAP band SuH NDECD_ROIs")
        },
        "reps"=1:8
    ),
    "SuH_nicdmut_necd1" = list(
        "file"="20160824_SuHallC band FRAP_ROIs",
        "fit" = function(r, p=list()){
            fittingBaseGenFit(r, p, "SuH_nicdmut_necd1", "20160824_SuHallC band FRAP_ROIs")
        },
        "reps"=1:6
    ),
    "SuH_nicdmut_necd1_2res" = list(
        "file"="20160824_SuHallC band FRAP_ROIs",
        "fit" = function(r, p=list()){
            fittingBaseGenFit2Res(r, p, "SuH_nicdmut_necd1_2res", "20160824_SuHallC band FRAP_ROIs")
        },
        "reps"=1:6
    ),
    "SuH_nicdmut_necd2" = list(
        "file"="20160825_SuHallC band FRAP_ROIs",
        "fit" = function(r, p=list()){
            fittingBaseGenFit(r, p, "SuH_nicdmut_necd2", "20160825_SuHallC band FRAP_ROIs")
        },
        "reps"=1:16
    ),
    "SuH_nicdmut_necd2_2res" = list(
        "file"="20160825_SuHallC band FRAP_ROIs",
        "fit" = function(r, p=list()){
            fittingBaseGenFit2Res(r, p, "SuH_nicdmut_necd2_2res", "20160825_SuHallC band FRAP_ROIs")
        },
        "reps"=1:16
    ),
    "SuH_nicdmut_necd3" = list(
        "file"="20160908_FRAP band SuHallC NDECD_ROIs",
        "fit" = function(r, p=list()){
            fittingBaseGenFit(r, p, "SuH_nicdmut_necd3", "20160908_FRAP band SuHallC NDECD_ROIs")
        },
        "reps"=1:14
    ),
    "SuH_nicdmut_necd3_2res" = list(
        "file"="20160908_FRAP band SuHallC NDECD_ROIs",
        "fit" = function(r, p=list()){
            fittingBaseGenFit2Res(r, p, "SuH_nicdmut_necd3_2res", "20160908_FRAP band SuHallC NDECD_ROIs")
        },
        "reps"=1:14
    ),
    "SuH_nicdmut_necd4" = list(
        "file"="20160913_FRAP Band SuHallC_ROIs",
        "fit" = function(r, p=list()){
            fittingBaseGenFit(r, p, "SuH_nicdmut_necd4", "20160913_FRAP Band SuHallC_ROIs")
        },
        "reps"=1:13
    ),
    "SuH_nicdmut_necd4_2res" = list(
        "file"="20160913_FRAP Band SuHallC_ROIs",
        "fit" = function(r, p=list()){
            fittingBaseGenFit2Res(r, p, "SuH_nicdmut_necd4_2res", "20160913_FRAP Band SuHallC_ROIs")
        },
        "reps"=1:13
    ),
    "SuH_mamdn_necd1" = list(
        "file"="20160822_SuHWT Mam NDECD_ROIs",
        "fit" = function(r, p=list()){
            fittingBaseGenFit(r, p, "SuH_mamdn_necd1", "20160822_SuHWT Mam NDECD_ROIs")
        },
        "reps"=1:23
    ),
    "SuH_mamdn_necd1_2res" = list(
        "file"="20160822_SuHWT Mam NDECD_ROIs",
        "fit" = function(r, p=list()){
            fittingBaseGenFit2Res(r, p, "SuH_mamdn_necd1_2res", "20160822_SuHWT Mam NDECD_ROIs")
        },
        "reps"=1:23
    ),
    "SuH_mamdn_necd2" = list(
        "file"="20160901_SuH WT mam NDECD_ROIs",
        "fit" = function(r, p=list()){
            fittingBaseGenFit(r, p, "SuH_mamdn_necd2", "20160901_SuH WT mam NDECD_ROIs")
        },
        "reps"=1:10
    ),
    "SuH_mamdn_necd2_2res" = list(
        "file"="20160901_SuH WT mam NDECD_ROIs",
        "fit" = function(r, p=list()){
            fittingBaseGenFit2Res(r, p, "SuH_mamdn_necd2_2res", "20160901_SuH WT mam NDECD_ROIs")
        },
        "reps"=1:10
    ),
    "SuH_peripodial"=list(
        "file"="20160506_FRAP SuH peripodial cells_ROIs",
        "fit"=function(r,
                       p=list()){
            p = processParamsStd(p)
            if(is.null(p$D.range))
                p$D.range = D.range.std
            if(is.null(p$save.base))
                p$save.base = "SuH_peripodial"
            
            fitGeneric(replicate=r,
                       experiment = "20160506_FRAP SuH peripodial cells_ROIs",
                       time.inc = 0.119,
                       time.inc2.start = NULL,
                       real.len.x = 30.75 * 1e-6,
                       real.len.y = 17.99 * 1e-6,
                       algorithm= p$algorithm,
                       D.range = p$D.range,
                       res.time.range=p$res.time.range,
                       Free.range=p$Free.range,
                       use.old.init.rates = p$use.old.init.rates,
                       scale.by = p$scale.by,
                       save.to.file=p$save.base,
                       roi.file=p$roi.file,
                       roi.fg=p$roi.fg,
                       roi.total.file=p$roi.total.file,
                       first.ROI.only=p$first.ROI.only,
                       bleaching.init.version=p$bleaching.init.version,
                       bleach.norm.sd.range=p$bleach.norm.sd.range,
                       res.time.c2.range=p$res.time.c2.range, 
                       prop.c2.range=p$prop.c2.range,
                       roi.exclude = p$roi.exclude,
                       use.uniform.init=p$use.uniform.init
            )
        },
        "reps"=1:10
    ),
    "SuH_WT_slow"=list(
        "file"="SuHWT",
        "fit"=function(r,
                       p = list()
        ){
            if(is.null(p$D.range))
                p$D.range = D.range.std
            if(is.null(p$res.time.range))
                p$res.time.range = res.time.range.std
            if(is.null(p$Free.range))
                p$Free.range = seq(0.05, 0.95, 0.05)
            if(is.null(p$save.base))
                p$save.base = "SuH_WT_slow"
            if(is.null(p$algorithm))
                p$algorithm = "full"
            if(is.null(p$use.old.init.rates))
                p$use.old.init.rates = TRUE
            if(is.null(p$post.avg))
                p$post.avg = 1
            if(is.null(p$scale.by))
                p$scale.by = 10
            
            fitGeneric(replicate=r,
                       experiment = "SuHWT",
                       time.inc = 0.26,
                       real.len.x = 64.58 * 1e-6,
                       real.len.y = 37.79 * 1e-6,
                       time.inc2.start=NULL,
                       algorithm= p$algorithm,
                       D.range = p$D.range,
                       res.time.range= p$res.time.range,
                       Free.range= p$Free.range,
                       use.old.init.rates = p$use.old.init.rates,
                       scale.by = p$scale.by,
                       save.to.file= p$save.base,
                       post.avg=p$post.avg
            )
        },
        "reps"=1:12
    ),
    "SuH_DNAmutant_slow"=list(
        "file"="R266H",
        "fit"=function(r,
                       p = list()
        ){
            if(is.null(p$D.range))
                p$D.range = D.range.std
            if(is.null(p$res.time.range))
                p$res.time.range = res.time.range.std
            if(is.null(p$Free.range))
                p$Free.range = seq(0.05, 0.95, 0.05)
            if(is.null(p$save.base))
                p$save.base = "SuH_WT_slow"
            if(is.null(p$algorithm))
                p$algorithm = "full"
            if(is.null(p$use.old.init.rates))
                p$use.old.init.rates = TRUE
            if(is.null(p$post.avg))
                p$post.avg = 1
            if(is.null(p$scale.by))
                p$scale.by = 10
            
            fitGeneric(replicate=r,
                       experiment = "R266H",
                       time.inc = 0.26,
                       real.len.x = 64.58 * 1e-6,
                       real.len.y = 37.79 * 1e-6,
                       time.inc2.start=NULL,
                       algorithm= p$algorithm,
                       D.range = p$D.range,
                       res.time.range= p$res.time.range,
                       Free.range= p$Free.range,
                       use.old.init.rates = p$use.old.init.rates,
                       scale.by = p$scale.by,
                       save.to.file= p$save.base,
                       post.avg=p$post.avg
            )
        },
        "reps"=1:12
    ),
    "SuH_locus_control_20170807" = list(
        "file"="20170807_FRAP locus tag ctrl_ROI",
        "fit" = function(r, p=list()){
            fittingBaseGenFit2017(r, p, "SuH_locus_control_20170807", "20170807_FRAP locus tag ctrl_ROI")
        },
        "reps"=1:15
    ),
    "SuH_locus_notch_20170807" = list(
        "file"="20170807_FRAP locus tag Notch",
        "fit" = function(r, p=list()){
            fittingBaseGenFit2017(r, p, "SuH_locus_notch_20170807", "20170807_FRAP locus tag Notch")
        },
        "reps"=1:13
    ),
    "SuH_locus_notch_20170818" = list(
        "file"="20170818_FRAP locus tag",
        "fit" = function(r, p=list()){
            fittingBaseGenFit2017(r, p, "SuH_locus_notch_20170818", "20170818_FRAP locus tag")
        },
        "reps"=1:9
    ),
    "SuH_locus_notch_20170818_5" = list(
        "file"="20170818_FRAP locus tag",
        "fit" = function(r, p=list()){
            p$scale.by = 5
            fittingBaseGenFit2017(r, p, "SuH_locus_notch_20170818", "20170818_FRAP locus tag")
        },
        "reps"=5
    ),
    "SuH_locus_notch_20170818_6" = list(
        "file"="20170818_FRAP locus tag",
        "fit" = function(r, p=list()){
            p$scale.by = 5
            fittingBaseGenFit2017(r, p, "SuH_locus_notch_20170818", "20170818_FRAP locus tag")
        },
        "reps"=6
    ),
    "SuH_mamdn_20170808" = list(
        "file"="20170808_FRAP mamDN",
        "fit" = function(r, p=list()){
            fittingBaseGenFit2017(r, p, "SuH_mamdn_20170808", "20170808_FRAP mamDN")
        },
        "reps"=1:8
    ),
    "SuH_mamdn_20170905" = list(
        "file"="20170905_FRAP mamDN",
        "fit" = function(r, p=list()){
            fittingBaseGenFit2017(r, p, "SuH_mamdn_20170905", "20170905_FRAP mamDN")
        },
        "reps"=1:8
    ),
    "SuH_nbm_20170814" = list(
        "file"="20170814_FRAP SuH NBM",
        "fit" = function(r, p=list()){
            fittingBaseGenFit2017(r, p, "SuH_nbm_20170814", "20170814_FRAP SuH NBM")
        },
        "reps"=1:9
    ),
    "SuH_nbm_20170822" = list(
        "file"="20170822_FRAP SuH NBM",
        "fit" = function(r, p=list()){
            fittingBaseGenFit2017(r, p, "SuH_nbm_20170822", "20170822_FRAP SuH NBM")
        },
        "reps"=1:16
    ),
    "SuH_mamri_20170815" = list(
        "file"="20170815_FRAP mamRi",
        "fit" = function(r, p=list()){
            fittingBaseGenFit2017(r, p, "SuH_mamri_20170815", "20170815_FRAP mamRi")
        },
        "reps"=1:13
    ),
    "SuH_mamri_20170901" = list(
        "file"="20170901_FRAP mamRi",
        "fit" = function(r, p=list()){
            fittingBaseGenFit2017(r, p, "SuH_mamri_20170901", "20170901_FRAP mamRi")
        },
        "reps"=1:15
    ),
    "SuH_sf8_a9_20170809" = list(
        "file"="20170809_FRAP SuH SF8 A9",
        "fit" = function(r, p=list()){
            fittingBaseGenFit2017(r, p, "SuH_sf8_a9_20170809", "20170809_FRAP SuH SF8 A9")
        },
        "reps"=1:10
    ),
    "SuH_sf8_a9_20170810" = list(
        "file"="20170810_FRAP SuH SF8 A9",
        "fit" = function(r, p=list()){
            fittingBaseGenFit2017(r, p, "SuH_sf8_a9_20170810", "20170810_FRAP SuH SF8 A9")
        },
        "reps"=1:10
    ),
    "SuH_sf8_a9_20170810_4" = list(
        "file"="20170810_FRAP SuH SF8 A9",
        "fit" = function(r, p=list()){
            p$scale.by = 5
            fittingBaseGenFit2017(r, p, "SuH_sf8_a9_20170810", "20170810_FRAP SuH SF8 A9")
        },
        "reps"=4
    ),
    "SuH_sf8_a9_20170810_5" = list(
        "file"="20170810_FRAP SuH SF8 A9",
        "fit" = function(r, p=list()){
            p$scale.by = 5
            fittingBaseGenFit2017(r, p, "SuH_sf8_a9_20170810", "20170810_FRAP SuH SF8 A9")
        },
        "reps"=5
    ),
    "SuH_WT_20170809" = list(
        "file"="20170809_FRAP SuH WT",
        "fit" = function(r, p=list()){
            fittingBaseGenFit2017(r, p, "SuH_WT_20170809", "20170809_FRAP SuH WT")
        },
        "reps"=1:10
    ),
    "SuH_WT_20170810" = list(
        "file"="20170810_FRAP SuH WT",
        "fit" = function(r, p=list()){
            fittingBaseGenFit2017(r, p, "SuH_WT_20170810", "20170810_FRAP SuH WT")
        },
        "reps"=1:8
    )
)


fittingFunctions = list(
    "GFP"=list(
        "fit"=function(r){
            fitGeneric(replicate=r,
                       experiment = "20160912_pointFRAP GFP",
                       time.inc = 0.115,
                       time.inc2.start = 310,
                       time.inc2.value = 2,
                       real.len.x = 64.58 * 1e-6,
                       real.len.y = 40.30 * 1e-6,
                       algorithm="diffusion",
                       D.range =seq(10, 40, 1) * 1e-12,
                       res.time.range = 0,
                       Free.range = 1,
                       use.old.init.rates = TRUE,
                       scale.by = 10,
                       save.to.file="GFP"
            )
        },
        "reps"=1:21
    ),
    "SuH_DNAmutant_full"=list(
        "fit"=function(r, 
                       p){
            fittingBase$SuH_DNAmutant$fit(r, p)
        },
        "reps"=fittingBase$SuH_DNAmutant$reps
    )
)


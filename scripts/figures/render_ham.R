# Figure 5 publication layout V2. Reads only validated new map assets; no model fitting.
# Usage: Rscript Figure5_HAM_Layout_V2.R <new-map-output-directory>
suppressPackageStartupMessages({library(grid);library(png);library(jsonlite);library(digest)})
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=1L) stop('Supply the map-data output directory.')
out <- normalizePath(args[1],winslash='/',mustWork=TRUE)
manifest <- fromJSON(file.path(out,'Figure5_Map_Provenance_V1.json'),simplifyVector=FALSE)
stopifnot(manifest$Status=='COMPLETE_MAP_DATA',manifest$MapCount==11L,manifest$CurveRows==1100L)
for(n in names(manifest$OutputsSHA256)) stopifnot(digest(file=file.path(out,n),algo='sha256')==manifest$OutputsSHA256[[n]])
tab <- read.csv(file.path(out,'Figure5_Map_Summary_V1.csv'),check.names=FALSE)
stopifnot(nrow(tab)==11L,identical(tab$Day,c(57L,86L,96L,113L,130L,28L,47L,68L,88L,96L,111L)))
targets <- file.path(out,paste0('Figure5_HAM_Within_Image_V2.',c('png','jpg','pdf')))
if(any(file.exists(targets))) stop('Figure export exists; use a new numbered output directory.')
white <- gpar(col='white',fontfamily='sans')
text_at <- function(label,x,y,size=12,rot=0,bold=FALSE) grid.text(label,x,y,rot=rot,gp=gpar(col='white',fontfamily='sans',fontsize=size,fontface=if(bold) 'bold' else 'plain'))
place_image <- function(path,x,y,w,h) {
  im <- readPNG(file.path(out,path)); ar <- dim(im)[2]/dim(im)[1]
  # Device aspect matters because normalized x/y units have different physical lengths.
  dims <- dev.size('in'); iw <- min(w*dims[1],h*dims[2]*ar); ih <- iw/ar
  grid.raster(im,x=x,y=y,width=unit(iw,'in'),height=unit(ih,'in'),interpolate=FALSE)
}
colorbar <- function(y,h,label) {
  grid.raster(readPNG(file.path(out,'Map_Assets/Colorbar_Jet_V1.png')),x=.913,y=y,width=.018,height=h,interpolate=FALSE)
  for(z in seq(0,1,.2)) text_at(sprintf('%.1f',z),.944,y-h/2+z*h,8)
  text_at(label,.982,y,11,270)
}
draw <- function() {
 grid.newpage(); grid.rect(gp=gpar(fill='black',col=NA))
 # Original two-panel black-background arrangement, with RGB above HAM in each panel.
 text_at('A',.025,.965,22,bold=TRUE)
 xs <- seq(.145,.815,length.out=5)
 for(i in 1:5) {
   text_at(tab$Day[i],xs[i],.965,15)
   place_image(tab$RGBFile[i],xs[i],.882,.130,.142)
   place_image(tab$HAMFile[i],xs[i],.720,.130,.142)
 }
 text_at('RGB',.055,.882,13,90); text_at('HAM',.055,.720,13,90)
 colorbar(.809,.30,expression(H^2~'(PSRI)'))
 text_at('B',.025,.612,22,bold=TRUE)
 xs <- seq(.145,.845,length.out=6)
 for(k in 1:6) {
   i <- k+5
   text_at(tab$Day[i],xs[k],.612,15)
   place_image(tab$RGBFile[i],xs[k],.463,.108,.267)
   place_image(tab$HAMFile[i],xs[k],.170,.108,.267)
 }
 text_at('RGB',.055,.463,13,90); text_at('HAM',.055,.170,13,90)
 colorbar(.317,.55,expression(H^2~'(ExGR)'))
}
png(targets[1],width=4500,height=5220,res=600,type='cairo',bg='black');draw();dev.off()
jpeg(targets[2],width=4500,height=5220,res=600,quality=98,type='cairo',bg='black');draw();dev.off()
cairo_pdf(targets[3],width=7.5,height=8.7,bg='black',onefile=TRUE);draw();dev.off()
script <- sub('^--file=','',commandArgs()[grep('^--file=',commandArgs())])
receipt <- list(Status='COMPLETE_FIGURE_EXPORT',DimensionsInches=c(7.5,8.7),RasterPixels=c(4500,5220),DPI=600,
                Layout='Two species panels; RGB then HAM; maize rotated 90 degrees counterclockwise without resampling numerical data',
                RendererSHA256=digest(file=script,algo='sha256'),InputManifestSHA256=digest(file=file.path(out,'Figure5_Map_Provenance_V1.json'),algo='sha256'),
                OutputSHA256=setNames(lapply(targets,function(p)digest(file=p,algo='sha256')),basename(targets)))
write_json(receipt,file.path(out,'Figure5_Layout_Provenance_V2.json'),auto_unbox=TRUE,pretty=TRUE)
writeLines(capture.output(sessionInfo()),file.path(out,'Figure5_Layout_SessionInfo_V2.txt'))
cat('COMPLETE_FIGURE_EXPORT\n')

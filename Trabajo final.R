rm(list = ls())

library(rspatial)
library(rgdal)
library(RColorBrewer)
library(sf)
library(ggplot2)
library(spdep)
library(maptools)
library(tidyverse)
library(rmarkdown)
library(ggmap)
library(readxl)
library(adespatial)
library(adegenet)
library(spdep)
library(scales)
library(akima)
library(mapview)
library(leaflet)
library(sp)
library(data.table)
library(lattice)

carpeta_raiz = "C:\\Users\\gperetti\\Google Drive\\Tesis\\"

# Cargo el shapefile de las fracciones censales de Córdoba
cba_shp <- readOGR(dsn = paste(carpeta_raiz, "Mapas", sep = ""), 
                   layer = "capa provincial - cordoba - division censal fracciones - 2010 - shp")

# Filtro solo las fracciones censales de la Ciudad de Córdoba (Departamento 014)
cba_shp <- cba_shp[cba_shp$DEPTO == "014",]
cba_shp@data <- filter(cba_shp@data, cba_shp@data$DEPTO == "014")

# Creo una nueva columna con el ID de la fracción
cba_shp@data <- cbind(CODIGO = paste(cba_shp@data$PROV, cba_shp@data$DEPTO, cba_shp@data$FRAC,
                                     sep = ""), cba_shp@data)

# Cargo el dataframe con los datos del censo
df_fracciones <- fread(paste(carpeta_raiz, "df_fracciones_revisado.csv", sep = ""), dec=",")
colnames(df_fracciones)[1] <- "CODIGO"

# Junto la tabla de atributos con el dataframe
cba_shp@data <- merge(cba_shp@data, df_fracciones, by = "CODIGO")

# Elimino columnas que no me interesan
cba_shp@data[c("R14014_", "R14014_ID", "PROV", "CODLOC", "NOMMUNI", "NOMLOC", 
               "FRAC", "RADIO", "FRAC14_", "FRAC14_ID", "fraccion")] <- NULL

# Grafico algunas cosas solo para chequear que esté el mapa completo
centroids <- getSpPPolygonsLabptSlots(cba_shp)
plot(centroids)
mapview(cba_shp)
plot(cba_shp$satisfactoria)

# PCA
colnames(cba_shp@data)

cba_shp_df <- data.frame(cba_shp)[, c(9:45)]
xy = coordinates(cba_shp)
frac_names <- data.frame(cba_shp)[, 1]
col_frac <- colors()[c(149, 254, 468, 552, 26)]

pca <- dudi.pca(cba_shp_df, scannf = FALSE, nf = 3)
scatter(pca, cex = 0.1)
summary(pca)

pca$eig
pca$eig[1]/sum(pca$eig)*100

pca$c1

sum(pca$eig/sum(pca$eig) * 100)

s.corcircle(pca$co, clabel = 0.6)

# Spatial partition
library(ade4)

bet <- bca(pca, factor(frac_names), scannf = FALSE, nf = 2)

bet$eig/sum(bet$eig) * 100

plot(bet)

# PCAIV
poly.xy <- poly(xy, degree = 2)

pcaiv.xy <- pcaiv(pca, poly.xy, scannf = FALSE, nf = 2) 

sum(pcaiv.xy$eig)/sum(pca$eig) * 100

pcaiv.xy$eig/sum(pcaiv.xy$eig) * 100

# PCAIV-MEM

#png(file = "figs/fig-fig7.png", width = 10, height = 4, units = "in", res = 72) 
mem <- scores.listw(lw) 
par(mfrow = c(2, 5), mar = rep(0.1, 4)) 
for (i in 1:10) { 
  plot(cba_shp, col = "grey95", border = "grey") 
  s.value(xy, mem$vectors[, i], add.plot = TRUE, clegend = 0) 
  text(270000, 2600000, bquote(paste("MC=", .(round(mem$values[i], 3)))), cex = 1.5) } 
#dev.off()

pcaiv.mem <- pcaiv(pca, mem$vectors, scannf = FALSE, nf = 2)

pcaiv.mem$eig/sum(pcaiv.mem$eig) * 100


# Multispati

nb <- poly2nb(cba_shp) 
lw <- nb2listw(nb, style = "W") 

ms <- multispati(pca, lw, scannf = FALSE)

plot(ms)

sum.ms <- summary(ms)

sum.ms

s.arrow(ms$c1, clabel = 0.8)

s.match(ms$li, ms$ls, clabel = 0, pch = 15) 
s.match(ms$li[c(10, 41, 27), ], ms$ls[c(10, 41, 27), ], 
        label = frac_names[c(10, 41, 27)], clabel = 0.8, add.plot = TRUE, pch = 15)


# sPCA ------------------------------------------------------------------------------------------------

cba_vecinos_list <- poly2nb(cba_shp)

mySpca <- spca.data.frame(x = scale(cba_shp@data[c(5:41)]), xy = coordinates(cba_shp), 
                                              cn = cba_vecinos_list, scannf = FALSE)

barplot(mySpca$eig, main="Valores propios del sPCA",
        col=spectral(length(mySpca$eig)))
legend("topright", fill=spectral(2),
       leg=c("Estructuras globales", "Estructuras locales"))
abline(h=0,col="grey")

plot(mySpca)

# Test global para contrastar la presencia de estructuras globales
myGtest <- global.rtest(cba_shp@data[c(5:41)], mySpca$lw, nperm=99)
myGtest
plot(myGtest)

# Test local para contrastar la presencia de estructuras locales
myLtest <- local.rtest(cba_shp@data[c(5:41)], mySpca$lw, nperm=99)
myLtest
plot(myLtest)

# Hago diferentes gráficos
# Gráfico 1
colorplot(mySpca, cex=3, main="colorplot of mySpca, first global score")

# Gráfico 2
x <- coordinates(cba_shp)[,1]
y <- coordinates(cba_shp)[,2]
temp <- interp(x, y, mySpca$ls[,1]) 
image(temp, col=azur(100)) 
points(x,y)

# Gráfico 3
interpX <- seq(min(x),max(x),le=200) 
interpY <- seq(min(y),max(y),le=200) 
temp <- interp(x, y, mySpca$ls[,1], xo=interpX, yo=interpY) 
image(temp, col=azur(100)) 
points(x,y)

# Gráfico 4
myPal <- colorRampPalette(c("firebrick2", "white", "lightslateblue")) 
annot <- function(){ 
  title("sPCA - interpolated map of individual scores") 
  points(x,y) 
} 
filled.contour(temp, color.pal=myPal, nlev=50, 
               key.title=title("lagged \nscore 1"), plot.title=annot())

# Gráfico 5
myLoadings <- mySpca$c1[,1]^2 
names(myLoadings) <- rownames(mySpca$c1) 
loadingplot(myLoadings, xlab="Alleles", 
            ylab="Weight of the alleles", 
            main="Contribution of alleles \n to the first sPCA axis")
myLoadings

# -----------------------------------------------------------------------------------------------------

# Regresión espacial
grps <- 10
brks <- quantile(cba_shp$tenencia_inquilino, 0:(grps-1)/(grps-1), na.rm=TRUE)
p <- spplot(cba_shp, "tenencia_inquilino", 
            at=brks, col.regions=rev(brewer.pal(grps, "RdBu")), col="transparent")
p

f1 <- gasto_pc ~ tenencia_propietario + construccion_basica + tenencia_inquilino
m1 <- lm(f1, data=cba_shp)
summary(m1)

cba_shp$residuals <- residuals(m1)
brks <- quantile(cba_shp$residuals, 0:(grps-1)/(grps-1), na.rm=TRUE)
spplot(cba_shp, "residuals", at=brks, col.regions=rev(brewer.pal(grps, "RdBu")), col="black")

plot(cba_shp)
plot(cba_vecinos_list, coordinates(cba_shp), col='red', lwd=2, add=TRUE)

resnb <- sapply(cba_vecinos_list, function(x) mean(cba_shp$residuals[x]))
cor(cba_shp$residuals, resnb)

plot(cba_shp$residuals, resnb, xlab='Residuals', ylab='Mean adjacent residuals')

lw <- nb2listw(cba_vecinos_list)
moran.mc(cba_shp$residuals, lw, 999)

# Hago los test de diagnóstico para ver cual modelo usar
cba_vecinos_list
lmLMtests <- lm.LMtests(m1, lw, test=c("LMerr", "LMlag", "RLMerr", "RLMlag", "SARMA"))
lmLMtests

# Spatial lag model
library(spatialreg)
m1s = lagsarlm(f1, data=cba_shp, lw)
summary(m1s)

impacts(m1s, listw = lw)

# Evaluamos los p-value
summary(impacts(m1s, listw=lw, R=500),zstats=TRUE)

cba_shp$residuals <- residuals(m1s)
moran.mc(cba_shp$residuals, lw, 999)

brks <- quantile(cba_shp$residuals, 0:(grps-1)/(grps-1), na.rm=TRUE)
p <- spplot(cba_shp, "residuals", at=brks, col.regions=rev(brewer.pal(grps, "RdBu")), col="transparent")
p

# Spatial error model
m1e = errorsarlm(f1, data=cba_shp, lw, tol.solve=1.0e-30)
summary(m1e)

cba_shp$residuals <- residuals(m1e)
moran.mc(cba_shp$residuals, lw, 999)

brks <- quantile(cba_shp$residuals, 0:(grps-1)/(grps-1), na.rm=TRUE)
p <- spplot(cba_shp, "residuals", at=brks, col.regions=rev(brewer.pal(grps, "RdBu")), col="transparent")
p


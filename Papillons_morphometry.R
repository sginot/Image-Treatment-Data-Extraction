#Data extraction after segmentation of images and quite a lot of cleaning of  the dataset#

library(Momocs)
library(MASS)
library(scales)
source("/home/samuel/Documents/Rfunctions1.txt") #Text file containing code of function from Claude (2008). Morphometrics with R.

#################################
### Make Momocs-suited object ###
#################################

ar.shapes <- array(NA, dim=c(200,2,length(extr$shapes))) #Create empty array
bugged <- which(unlist(lapply(extr$shapes, is.null))) #Find index of NULL elements in shape list
for (i in 1:length(extr$shapes)){
tryCatch({ar.shapes[,,i] <- extr$shapes[[i]]}, error=function(e){})} #Fill the array, NULL elements will leave the corresponding matrices in the array as NA

ar.shapes <- ar.shapes[,,-bugged] #Remove NULL/NA from array

Out.butterflies <- Out(ar.shapes) #Make "Out" object (for Momocs)

bad<-c(267,258,386,408,859,860,1048,1165,1207,1288,1289,1309,1426,1542) #Vector of 'messed-up' shapes

############################
### Alignment strategies ###
############################

#Baseline alignment (Bookstein style) based on extreme points along x (supposedly tips of wings)

al.ar.shapes <- array(NA, dim=dim(ar.shapes)) #Empty array

for (i in 1:dim(ar.shapes)[3]) {
shp <- ar.shapes[,,i]
maxx <- which.max(shp[,1]) #Index of point with max x value
minx <- which.min(shp[,1]) #Index of point with min x value
al.ar.shapes[,,i] <- booksteinM(shp, maxx, minx) #Baseline alignment function
} #Problem : some butterflies have wider hindwings, mdpor max x value is not of 'tip' of wing

#Same but with slight correction
al.ar.shapes2 <- array(NA, dim=dim(ar.shapes)) #Empty array

for (i in 1:dim(ar.shapes)[3]) {
shp <- ar.shapes[,,i]
meany <- mean(shp[,2]) #Mean y value for points
no <- which(shp[,2] > meany)
nomax <- which.max(shp[no,1]) #Index of point with max x value, and y > mean
nomin <- which.min(shp[no,1]) #Index of point with min x value, and y > mean
maxx <- no[nomax]
minx <- no[nomin]

al.ar.shapes2[,,i] <- booksteinM(shp, maxx, minx) #Baseline alignment function
} #Improvement but biased if maximum x and minimum x points are not horizontally aligned

#Approach with maximum mean x,y coordinate = top right point
shp <- ar.shapes[,,1]
which.max(apply(shp,1,mean)) #134
p1 <- 134
y.line <- shp[134,2]
meanx <- mean(shp[,1]) #Mean x value for points
opoint <- c(meanx-(shp[134,1]-meanx), y.line)
dm <- dist(rbind(opoint,shp))
which.min(dm[1:dim(shp)[1]]) #200
p2 <- 200

#Approach with maximum mean x,y coordinate, and centering and absolute values
al.ar.shapes3 <- array(NA, dim=dim(ar.shapes)) #Empty array

for (i in 1:dim(ar.shapes)[3]) {
shp <- ar.shapes[,,i]
p1 <- which.max(apply(shp,1,mean))

shp.cent <- shp
shp.cent[,1] <- shp[,1]-mean(shp[,1])
left <- shp.cent[which(shp.cent[,1] < 0),]

p2.left <- which.max(apply(abs(left),1,mean))
p2 <- which(shp.cent[,1]==left[p2.left,1] & shp.cent[,2]==left[p2.left,2])

al.ar.shapes3[,,i] <- booksteinM(shp.cent, p1, p2) #Baseline alignment function
} #Aligment looks fairly good BUT will bias the amount of variation is the anterior wing.


#Superimposition on two points (tips of BODY) 

al.ar.shapes4 <- array(NA, dim=dim(ar.shapes))

for (i in 1:dim(ar.shapes)[3]){

shp <- ar.shapes[,,i]

maxx <- max(shp[,1]) #Max x value
minx <- min(shp[,1]) #Min x value
midx <- mean(c(maxx,minx))
meany <- mean(shp[,2])
lower.h <- shp[which(shp[,2] < meany),]
upper.h <- shp[which(shp[,2] > meany),]

lower.tip <- which.min(abs(lower.h[,1]-midx))
upper.tip <- which.min(abs(upper.h[,1]-midx))

p1 <- which(shp[,1] == lower.h[lower.tip,1] & shp[,2] == lower.h[lower.tip,2])
p2 <- which(shp[,1] == upper.h[upper.tip,1] & shp[,2] == upper.h[upper.tip,2])

al.ar.shapes4[,,i] <- booksteinM(shp, p1, p2)
}
#Some outliers to remove
bad2 <- c(201,  208,  687, 1294, 1535)
Out.butt.bookstein <- Out(al.ar.shapes4[,,-bad2])
pca.butt.bookstein <- PCA(efourier(Out.butt.bookstein, norm=F))
#Weird results


#Combine both alignments to get 4 points? 

#Resistant-fit approach, avoids pinocchio effect

grf.ar <- grf2(al.ar.shapes3[,,-bad]) #Does not converge.

#Procrustes

p.ar <- pgpa(al.ar.shapes3[,,-bad]) #Partial Procrustes a la Claude
o.ar <- orp(p.ar$rotated)

fp.ar <- fgProcrustes(al.ar.shapes3[,,-bad]) #Full Procrustes, Momocs version (based on fgpa2)

plot(p.ar$rotated[,1,],p.ar$rotated[,2,], type="l", asp=1) #Probably not the best idea to Procrustes on non-homologous landmarks. Results in 'over-rotated' shapes, because points are not 'in order'.

###### Probably the best approach is the one that follows #######

#Procrustes on 4 landmarks: tips of anterior wings, anterior and posterior tips of body along symmetry axis. 
#In Momocs, Out objects can have a $ldk : a list with index of landmarks that are to be superimposed
ar.shapes <- ar.shapes[,,-bad]

ldk.ls <- vector("list", length=dim(ar.shapes)[3])

for (i in 1:dim(ar.shapes)[3]){

shp <- ar.shapes[,,i]
p1 <- which.max(apply(shp,1,mean)) #Tip of right anterior wing. Mean is useful, because hindwings can sometimes have a larger x value. But combining x and y make it sure that we catch the correct point

shp.cent <- shp
shp.cent[,1] <- shp[,1]-mean(shp[,1])
left <- shp.cent[which(shp.cent[,1] < 0),]

p2.left <- which.max(apply(abs(left),1,mean)) #Equivalent to p1, but abs() mirrors the shape
p2 <- which(shp.cent[,1]==left[p2.left,1] & shp.cent[,2]==left[p2.left,2])

maxx <- max(shp[,1]) #Max x value
minx <- min(shp[,1]) #Min x value
midx <- mean(c(maxx,minx))
meany <- mean(shp[,2])
lower.h <- shp[which(shp[,2] < meany),] #Lower half of shape
upper.h <- shp[which(shp[,2] > meany),] #Upper half of shape

lower.tip <- which.min(abs(lower.h[,1]-midx)) #Difference between x coordinate and mean x allows to find the point closest to the symmetry axis in terms of x.
upper.tip <- which.min(abs(upper.h[,1]-midx))

p3 <- which(shp[,1] == lower.h[lower.tip,1] & shp[,2] == lower.h[lower.tip,2])
p4 <- which(shp[,1] == upper.h[upper.tip,1] & shp[,2] == upper.h[upper.tip,2])

ldk.ls[[i]] <- c(p1,p2,p3,p4)
}
###Illustrating the previous process
shp_center <- cbind(ar.shapes[,1,i]-mean(ar.shapes[,1,i]), ar.shapes[,2,i]-mean(ar.shapes[,2,i]))
plot(shp_center, asp=1, type="b", pch=20,lty=1, cex=0.7, xlab="X centered coordinates", ylab="Y centered cooredinates")
lines(shp_center, pch=20,lty=1)
points(0,0, pch=3)
points(shp_center[which.max(apply(shp_center,1,mean)),1],shp_center[which.max(apply(shp_center,1,mean)),2], col="firebrick", pch=20, cex=1.5)
text(shp_center[which.max(apply(shp_center,1,mean)),1],shp_center[which.max(apply(shp_center,1,mean)),2], labels="p1", col="firebrick", pos=3, cex=1.5)
left <- shp_center[which(shp_center[,1] < 0),]
p2.left <- which.max(apply(abs(left),1,mean)) #Equivalent to p1, but abs() mirrors the shape
p2 <- which(shp_center[,1]==left[p2.left,1] & shp_center[,2]==left[p2.left,2])
points(shp_center[p2,1], shp_center[p2,2], col="firebrick", pch=20, cex=1.5)
text(shp_center[p2,1], shp_center[p2,2],labels="p2", col="firebrick", pos=3, cex=1.5)

maxx <- max(shp_center[,1]) #Max x value
minx <- min(shp_center[,1]) #Min x value
midx <- mean(c(maxx,minx))
meany <- mean(shp_center[,2])
lower.h <- shp_center[which(shp_center[,2] < meany),] #Lower half of shape
upper.h <- shp_center[which(shp_center[,2] > meany),] #Upper half of shape
lower.tip <- which.min(abs(lower.h[,1]-midx)) #Difference between x coordinate and mean x allows to find the point closest to the symmetry axis in terms of x.
upper.tip <- which.min(abs(upper.h[,1]-midx))
p3 <- which(shp_center[,1] == lower.h[lower.tip,1] & shp_center[,2] == lower.h[lower.tip,2])
p4 <- which(shp_center[,1] == upper.h[upper.tip,1] & shp_center[,2] == upper.h[upper.tip,2])
points(shp_center[c(p3,p4),1], shp_center[c(p3,p4),2], col="orange", pch=20, cex=1.5)
text(shp_center[c(p3,p4),1], shp_center[c(p3,p4),2],labels=c("p3","p4"), col="orange", pos=c(1,3), cex=1.5)


Out.butterflies <- Out(ar.shapes)

Out.butterflies$ldk <- ldk.ls

p.butterflies <- fgProcrustes(Out.butterflies)

stack(p.butterflies) #Fair results...
fou.butterflies <- efourier(p.butterflies, norm=F) #NEED TO rm(efourier) and reload Momocs
#Still a problem of order of points. The outline should start from top-left corner of left wing.

#sfou.butterflies <- efourier(p.butterflies, nb.h=5, norm=F)
#pca.sfou <- PCA(sfou.butterflies) #Testing with less harmonics

#Can be solved by half shape superimposition and averaging + reordering of points (or just reordering)

### Reordering of points, starting from the anterior-right tip of the right wing ###
po.butterflies <- array(NA, dim=dim(ar.shapes))

for (i in 1:dim(po.butterflies)[3]) {
start <- ldk.ls[[i]][1] #This is p1 from previous loop
m <- p.butterflies[i]
mo <- m[[1]][c(start:200, 1:(start-1)),]
po.butterflies[,,i] <- mo
}

po.butterflies <- Out(po.butterflies) #RE-convert to momocs object

stack(po.butterflies)

#Checking how fourier deals with butterflies (4 first harmonics)
pdf(file=paste("butterfly", 1427))
calibrate_reconstructions_efourier(x=po.butterflies, range=4, id=1427)
dev.off()


fou.po <- efourier(po.butterflies, norm=F) #Fourier analysis

pca.fou<-PCA(fou.po)

plot(pca.fou, morphospace=T, density=T,lev.density=15, contour=T, lev.contour=15, cex=0.2, pos.shp="range") #momocs style plot
plot(pca.fou, xax=1, yax=2, loadings=T, zoom=0.5, morphospace=F, cex=0.2) #Show loadings
### New plot, more customized ###
hull.pca <- chull(pca.fou$x[,1],pca.fou$x[,2]) #Convex hull points
hull.o <- sort(hull.pca) # Reorder

f <- kde2d(pca.fou$x[,1], pca.fou$x[,2], n = 50) #Compute 2D kernel (function from MASS)

image(f, xlim=c(-0.4,0.25), ylim=c(-0.3,0.3), zlim = c(0.2, 60), xlab="PC1 (41.6%)", ylab="PC2 (20.9%)") #Plot background "landscape" with density of points
abline(h=0, col="gray80")
abline(v=0, col="gray80")
points(pca.fou$x[,1],pca.fou$x[,2], cex=0.2, pch=20, col="gray50")
points(pca.fou$x[hull.pca,1],pca.fou$x[hull.pca,2], cex=0.6, pch=20, col="gray50")
grid()

#Compute "NPC" coordinates for subplots to add shapes of convex hull points to the current plot.
library(grid)
coo.plot<-array(NA, dim=c(2,2,13))
for (i in 1:13) {
xy1 <- grid.locator(unit="npc")
coo.plot[1,1,i] <- as.numeric(substr(xy1$x, 1, 5))
coo.plot[1,2,i] <- as.numeric(substr(xy1$y, 1, 5))
coo.plot[2,,i] <- coo.plot[1,,i] + 0.08
}

#Auto-add shape sub-plots
for (i in 1:13) {
par(new=T, fig=c(coo.plot[,,i]), mar=c(0,0,0,0), xaxt="n", yaxt="n")
n <- hull.pca[i]
m <- po.butterflies[n]
m <- m[[1]]
plot(m[,1], m[,2], bty="n", type="l")
polygon(m[,1], m[,2], col="grey70")
}

#Function to add N arbitrary shape sub-plots
#N=number of points
#l=length of side of subplot in npc coordinates
#A=array of shapes (momocs style)
subplots <- function(x, y, N=20, l=0.08, A=po.butterflies, col="black", opacity=0.5) { 

coo.plots <- array(NA, dim=c(2,2,N))
n <- rep(NA, N)

   for (i in 1:N) { #First loop defines points and positions
   n[i] <- identify(x, y, n=1, plot=F)
   xy <- grid.locator(unit="npc")
   coo.plots[1,1,i] <- as.numeric(substr(xy$x, 1, 5))
   coo.plots[1,2,i] <- as.numeric(substr(xy$y, 1, 5))
   coo.plots[2,,i] <- coo.plots[1,,i] + l
   }
   for (i in 1:N) { #Second loop plots
   par(new=T, fig=c(coo.plots[,,i]), mar=c(0,0,0,0), xaxt="n", yaxt="n")
   m <- A[n[i]]
   m <- m[[1]]
   plot(m[,1], m[,2], bty="n", type="l")
   polygon(m[,1], m[,2], col=alpha(col, opacity))
   }
}

subplots(x=pca.fou$x[,1], y=pca.fou$x[,2], N=30, l=0.05, A=po.butterflies, col="black", opacity=0.2) #Use of function

###Checking for difference when using half-shapes

left.po <- vector("list")

for (i in 1:length(po.butterflies$coo)) {left.po[[i]] <- po.butterflies$coo[[i]][which(po.butterflies$coo[[i]][,1] >= 0),]}

left.po <- Out(left.po)

stack(left.po)

fou.left <- efourier(left.po, norm=F) #10 harmonics kept only

pca.left <- PCA(fou.left)

plot(pca.left, morphospace=T, density=T,lev.density=15, contour=T, lev.contour=15, cex=0.2, pos.shp="range")
cor.test(pca.fou$x[,1], pca.left$x[,1]) #r = -0.9, p = 0
cor.test(pca.fou$x[,2], pca.left$x[,2]) #r = -0.74, p = 0
cor.test(pca.fou$x[,3], pca.left$x[,3]) #r = -0.08, p = 0.001

right.po <- vector("list")

for (i in 1:length(po.butterflies$coo)) {right.po[[i]] <- po.butterflies$coo[[i]][which(po.butterflies$coo[[i]][,1] <= 0),]
right.po[[i]][,1] <- abs(right.po[[i]][,1])
}

right.po <- Out(right.po)

stack(right.po)

fou.right <- efourier(right.po, norm=F) #10 harmonics kept only

pca.right <- PCA(fou.right)


plot(pca.right, morphospace=T, density=T,lev.density=15, contour=T, lev.contour=15, cex=0.2, pos.shp="range")
cor.test(pca.fou$x[,1], pca.right$x[,1]) #r = 0.86, p = 0
cor.test(pca.fou$x[,2], pca.right$x[,2]) #r = 0.83, p = 0
cor.test(pca.fou$x[,3], pca.right$x[,3]) #r = -0.69, p = 0

cor.test(pca.left$x[,1], pca.right$x[,1]) #r = -0.88, p = 0
cor.test(pca.left$x[,2], pca.right$x[,2]) #r = -0.52, p = 0
cor.test(pca.left$x[,3], pca.right$x[,3]) #r = 0.46, p = 0


#Maybe good to try the "Protest" procedure cited in Renaud et al. 2020 to compare the morphospaces of left, right and combined
protest(pca.left$x, pca.fou$x[,1:40], scale=T)

#Call:
#protest(X = pca.left$x, Y = pca.fou$x[, 1:40], scale = T) 
#Procrustes Sum of Squares (m12 squared):        0.3885 
#Correlation in a symmetric Procrustes rotation: 0.782 
#Significance:  0.001 
#Permutation: free
#Number of permutations: 999

protest(pca.right$x, pca.fou$x[,1:40], scale=T)

#Call:
#protest(X = pca.right$x, Y = pca.fou$x[, 1:40], scale = T) 
#Procrustes Sum of Squares (m12 squared):        0.4543 
#Correlation in a symmetric Procrustes rotation: 0.7387 
#Significance:  0.001 
#Permutation: free
#Number of permutations: 999

protest(pca.right$x, pca.left$x, scale=T)

#Call:
#protest(X = pca.right$x, Y = pca.left$x, scale = T) 
#Procrustes Sum of Squares (m12 squared):        0.6899 
#Correlation in a symmetric Procrustes rotation: 0.5569 
#Significance:  0.001 
#Permutation: free
#Number of permutations: 999

###############################################
### Finding ways to solve symmetry problems ###
###############################################

###Possibility 1 : do fourier + PCA on all half shapes (with left mirrored) then use mean coordinates of both halves to get coordinate of individual. If two specific halves are very far, check and delete if necessary
###Possibility 2: do Fourier on halves separately (mirror left half), then average Fourier coefficients by individual. Then do PCA. => It actually doesn't seem that fourier coef for left and right half are strongly correlated, so bad idea.
###Possibility 3: If doing PCA on both right and left halves, ANOVA can discriminate between within individual (i.e. asymmetry) and between individual variation.

left.po <- vector("list")
right.po <- vector("list")

for (i in 1:length(po.butterflies$coo)) {
right.po[[i]] <- po.butterflies$coo[[i]][which(po.butterflies$coo[[i]][,1] >= 0),]


first <- which.max(apply(right.po[[i]], 1, mean))

right.po[[i]] <- right.po[[i]][c(first:dim(right.po[[i]])[1], 1:(first-1)),]

}

for (i in 1:length(po.butterflies$coo)) {
left.po[[i]] <- po.butterflies$coo[[i]][which(po.butterflies$coo[[i]][,1] <= 0),]

left.po[[i]] <- left.po[[i]][nrow(left.po[[i]]):1,]

left.po[[i]][,1] <- abs(left.po[[i]][,1])

first <- which.max(apply(left.po[[i]], 1, mean))

left.po[[i]] <- left.po[[i]][c(first:dim(left.po[[i]])[1], 1:(first-1)),]
}

all.halves <- c(right.po, left.po)

dfac<- data.frame(
       side=c(rep("right", 5791), rep("left", 5791)),
       museum=c(rep("MHNT",1715), rep("NHM", 1597), rep("NHM", 2479), rep("MHNT",1715), rep("NHM", 1597), rep("NHM", 2479)),
       family = c(as.character(fam.MHNT), as.character(fam.NHM), as.character(fam.BW), as.character(fam.MHNT), as.character(fam.NHM), as.character(fam.BW)),
       genus = c(as.character(gen.MHNT), as.character(gen.NHM), as.character(gen.BW), as.character(gen.MHNT), as.character(gen.NHM), as.character(gen.BW)),
       species = c(as.character(sp.MHNT), paste(gen.NHM, sp.NHM), as.character(sp.BW), as.character(sp.MHNT), paste(gen.NHM, sp.NHM), as.character(sp.BW)),
	dataset = c(rep("MHNT",1715), rep("Types", 1597), rep("Birdwing", 2479), rep("MHNT",1715), rep("Types", 1597), rep("Birdwing", 2479))
       )

out.all.halves <- Out(all.halves, fac=dfac)

fou.all.halves <- efourier(out.all.halves, norm=F)
pca.all.halves <- PCA(fou.all.halves) #Good because we can take in account symmetry problems with Anova for example.
side <- dfac$side
plot(pca.all.halves, col=c("red","blue")[side]) #Looks good after adding reordering of points, full overlap between right and left distributions.

oo<-manova(lm(pca.all.halves$x[,1:2] ~ dfac$side * dfac$museum + dfac$family + dfac$side:dfac$family)) #No interaction between museum and family because of singularity problems when using etasq
summary(oo) #WARNING: this is type I SS, so order of factors matters...
#                         Df  Pillai approx F num Df den Df    Pr(>F)    
#dfac$side                 1 0.00577     33.5      2  11527 3.207e-15 ***
#dfac$museum               1 0.37367   3438.6      2  11527 < 2.2e-16 ***
#dfac$family              23 0.66524    249.8     46  23056 < 2.2e-16 ***
#dfac$side:dfac$museum     1 0.00139      8.0      2  11527 0.0003295 ***
#dfac$side:dfac$family    23 0.02315      5.9     46  23056 < 2.2e-16 ***
#Residuals             11528                                             

etasq(oo)
#                            eta^2
#dfac$side             0.005773734
#dfac$museum           0.014168833
#dfac$family           0.332620080 #Taxonomic group is clearly the main effect here
#dfac$side:dfac$museum 0.004802565
#dfac$side:dfac$family 0.011575172

library(car)           

Anova(lm(pca.all.halves$x[,1] ~ dfac$family * dfac$museum * dfac$side))

#Response: pca.all.halves$x[, 1]
#                                  Sum Sq    Df  F value    Pr(>F)    
#dfac$family                       31.862    23 219.1722 < 2.2e-16 ***
#dfac$museum                        0.191     1  30.1841 4.013e-08 ***
#dfac$side                          0.386     1  61.0331 6.099e-15 ***
#dfac$family:dfac$museum            1.685     5  53.3223 < 2.2e-16 ***
#dfac$family:dfac$side              0.769    23   5.2866 2.676e-15 ***
#dfac$museum:dfac$side              0.194     1  30.7514 2.998e-08 ***
#dfac$family:dfac$museum:dfac$side  0.061     5   1.9199   0.08749 .  
#Residuals                         72.801 11518                       
 
Manova(lm(pca.all.halves$x[,1:2] ~ dfac$family + dfac$museum + dfac$side)) #No interaction due to singularity problems. Type II SS does not cause problems with order of factors.
sum(pca.all.halves$sdev[1:8]^2)/sum(pca.all.halves$sdev^2)*100 #91% of variance

Manova(lm(pca.all.halves$x[,1:8] ~ dfac$family + dfac$museum + dfac$side)) 
#Type II MANOVA Tests: Pillai test statistic
#            Df test stat approx F num Df den Df    Pr(>F)    
#dfac$family 23   1.61616  127.154    184  92416 < 2.2e-16 ***
#dfac$museum  1   0.14742  249.531      8  11545 < 2.2e-16 ***
#dfac$side    1   0.01078   15.726      8  11545 < 2.2e-16 ***

etasq(Manova(lm(pca.all.halves$x[,1:8] ~ dfac$family + dfac$museum + dfac$side)))
#                 eta^2
#dfac$family 0.20201990 #Taxonomic group is still the main effect, but it appears that adding PC axes tends to balance it with the museum effect
#dfac$museum 0.14741996
#dfac$side   0.01077982

### In any case, it appears clear that the effect of side is always minor, so let's average both sides and run another round of analyses.

fou.right <- efourier(Out(right.po), norm=F)
fou.left <- efourier(Out(left.po), norm=F)

fou.av <- (fou.right$coe + fou.left$coe)/2

dfac.av <- dfac[1:5791,]

pca.av.halves <- prcomp(fou.av) #Strongly correlated with right and left half PCAs.

protest(pca.right$x, pca.av.halves$x[,1:40], scale=T)
#Procrustes Sum of Squares (m12 squared):        0.5218 
#Correlation in a symmetric Procrustes rotation: 0.6915 
#Significance:  0.001 
#Permutation: free
#Number of permutations: 999

protest(pca.left$x, pca.av.halves$x[,1:40], scale=T)
#Procrustes Sum of Squares (m12 squared):        0.1532 
#Correlation in a symmetric Procrustes rotation: 0.9202 
#Significance:  0.001 
#Permutation: free
#Number of permutations: 999


all.gen <- dfac.av$genus
all.sp <- dfac.av$species

av.gen.fou <- matrix(NA, ncol=40, nrow=length(levels(all.gen)))
av.sp.fou <- matrix(NA, ncol=40, nrow=length(levels(all.sp)))

for (i in 1:40) {
av.gen.fou[,i] <- tapply(fou.av[,i], INDEX=all.gen, FUN=mean)
av.sp.fou[,i] <- tapply(fou.av[,i], INDEX=all.sp, FUN=mean)
}

rownames(av.gen.fou) <- levels(all.gen)
rownames(av.sp.fou) <- levels(all.sp)
colnames(av.gen.fou) <- colnames(av.sp.fou) <- colnames(fou.av)


###Do a nice plot for average half PCA.
pca.av.halves$sdev^2/sum(pca.av.halves$sdev^2)*100 #Variance explained by PCs
library(MASS)
hull.pca <- chull(pca.av.halves$x[,1],pca.av.halves$x[,2])
hull.o <- sort(hull.pca)
f1 <- kde2d(pca.av.halves$x[,1], pca.av.halves$x[,2], n = 50) 

image(f1, xlim=c(-0.2,0.40), ylim=c(-0.25,0.4), zlim = c(0.5, 73), xlab="PC1 (50.2%)", ylab="PC2 (18.5%)")
abline(h=0, col="gray80")
abline(v=0, col="gray80")
points(pca.av.halves$x[,1], pca.av.halves$x[,2], cex=0.2, pch=20, col="gray50")
points(pca.av.halves$x[hull.pca,1], pca.av.halves$x[hull.pca,2], cex=0.6, pch=20, col="gray50")
text(pca.av.halves$x[hull.pca,1], pca.av.halves$x[hull.pca,2],labels=hull.pca, pos=3, cex=0.5)
grid()
###WARNING SOMETHING IS WRONG IN THE MEAN POINTS (at least for genera)
#levels(gen.ind.clean)[gen.converg] <-check
#points(sp.meanPC1, sp.meanPC2, col=rainbow(24)[unique.sp.fam], pch=20, cex=0.8)
###
#points(gen.meanPC1[which(unique.gen.fam=="Sphingidae")], gen.meanPC2[which(unique.gen.fam=="Sphingidae")], col="darkblue", pch=20, cex=0.8)
#points(gen.meanPC1[which(unique.gen.fam=="Nymphalidae")], gen.meanPC2[which(unique.gen.fam=="Nymphalidae")], col="forestgreen", pch=20, cex=0.8)
#points(gen.meanPC1[which(unique.gen.fam=="Papilionidae")], gen.meanPC2[which(unique.gen.fam=="Papilionidae")], col="firebrick", pch=20, cex=0.8)
#points(gen.meanPC1[which(unique.gen.fam=="Pieridae")], gen.meanPC2[which(unique.gen.fam=="Pieridae")], col="darkorange", pch=20, cex=0.8)
#legend("topleft", legend=c("Sphingidae", "Nymphalidae", "Papilionidae", "Pieridae"), col=c("darkblue","forestgreen", "firebrick","darkorange"), pch=20)

#points(gen.meanPC1, gen.meanPC2, col=rainbow(24)[unique.gen.fam], pch=20, cex=0.9)
#legend("topright", legend=levels(unique.gen.fam), col=rainbow(24), pch=20, cex=0.7, bty="n")

#library(grid)
#coo.plot<-array(NA, dim=c(2,2,length(hull.pca)))
#for (i in 1:length(hull.pca)) {
#xy1 <- grid.locator(unit="npc")
#coo.plot[1,1,i] <- as.numeric(substr(xy1$x, 1, 5))
#coo.plot[1,2,i] <- as.numeric(substr(xy1$y, 1, 5))
#coo.plot[2,,i] <- coo.plot[1,,i] + 0.08 #0.08 is arbitrary, can be increased or decreased to change area #of subplots
#}

#Auto-add shape sub-plots
#for (i in 1:length(hull.pca)) {
#par(new=T, fig=c(coo.plot[,,i]), mar=c(0,0,0,0), xaxt="n", yaxt="n")
#n <- hull.pca[i]
#m <- po.butterflies[n]
#m <- m[[1]]
#plot(m[,1], m[,2], bty="n", type="l")
#polygon(m[,1], m[,2], col="grey70")
#}

#mean.nymph <- apply(pca.av.halves$x[which(fam.ind.clean == "Nymphalidae"),], 2, mean)
#mean.sphing <- apply(pca.av.halves$x[which(fam.ind.clean == "Sphingidae"),], 2, mean)
#mean.papi <- apply(pca.av.halves$x[which(fam.ind.clean == "Papilionidae"),], 2, mean)
#mean.ereb <- apply(pca.av.halves$x[which(fam.ind.clean == "Lymantriidae" | fam.ind.clean == "Arctiidae" | fam.ind.clean == "Erebidae"),], 2, mean)

#library(car)
#datell.Nymph <- dataEllipse(pca.av.halves$x[which(fam.ind.clean == "Nymphalidae"),1],pca.av.halves$x[which(fam.ind.clean == "Nymphalidae"),2])
#datell.Sphing <- dataEllipse(pca.av.halves$x[which(fam.ind.clean == "Sphingidae"),1],pca.av.halves$x[which(fam.ind.clean == "Sphingidae"),2])
#datell.Papi <- dataEllipse(pca.av.halves$x[which(fam.ind.clean == "Papilionidae"),1],pca.av.halves$x[which(fam.ind.clean == "Papilionidae"),2])
#datell.Ereb <- dataEllipse(pca.av.halves$x[which(fam.ind.clean == "Lymantriidae" | fam.ind.clean == "Arctiidae" | fam.ind.clean == "Erebidae"),1],pca.av.halves$x[which(fam.ind.clean == "Lymantriidae" | fam.ind.clean == "Arctiidae" | fam.ind.clean == "Erebidae"),2])

#lines(datell.Nymph[['0.95']][,1:2], col="forestgreen", lwd=1.5, lty=2)
#lines(datell.Nymph[['0.5']][,1:2], col="forestgreen", lwd=2)
#points(mean.nymph[1],mean.nymph[2], col="forestgreen", pch=20, cex=2)
#lines(datell.Sphing[['0.95']][,1:2], col="darkblue", lwd=1.5, lty=2)
#lines(datell.Sphing[['0.5']][,1:2], col="darkblue", lwd=2)
#points(mean.sphing[1],mean.sphing[2], col="darkblue", pch=20, cex=2)
#lines(datell.Papi[['0.95']][,1:2], col="firebrick", lwd=1.5, lty=2)
#lines(datell.Papi[['0.5']][,1:2], col="firebrick", lwd=2)
#points(mean.papi[1],mean.papi[2], col="firebrick", pch=20, cex=2)
#lines(datell.Ereb[['0.95']][,1:2], col="darkorange", lwd=1.5, lty=2)
#lines(datell.Ereb[['0.5']][,1:2], col="darkorange", lwd=2)
#points(mean.ereb[1],mean.ereb[2], col="darkorange", pch=20, cex=2)
#legend("topleft", legend=c("Nymphalidae","Sphingidae", "Papilionidae","Erebidae"), col=c("forestgreen", "darkblue", "firebrick", "darkorange"), lty=1, lwd=2)

#nymph.gen <- unique(gen.ind.clean[which(fam.ind.clean=="Nymphalidae")])
#sphing.gen <- unique(gen.ind.clean[which(fam.ind.clean=="Sphingidae")])
#papi.gen <- unique(gen.ind.clean[which(fam.ind.clean=="Papilionidae")])
#ereb.gen <- unique(gen.ind.clean[which(fam.ind.clean == "Lymantriidae" | fam.ind.clean == "Arctiidae" | fam.ind.clean == "Erebidae")])

#nymph.gen.meanPC1 <- gen.meanPC1[match(nymph.gen,names(gen.meanPC1))]
#nymph.gen.meanPC2 <- gen.meanPC2[match(nymph.gen,names(gen.meanPC2))]

#ereb.gen.meanPC1 <- gen.meanPC1[match(ereb.gen,names(gen.meanPC1))]
#ereb.gen.meanPC2 <- gen.meanPC2[match(ereb.gen,names(gen.meanPC2))]

#sphing.gen.meanPC1 <- gen.meanPC1[match(sphing.gen,names(gen.meanPC1))]
#sphing.gen.meanPC2 <- gen.meanPC2[match(sphing.gen,names(gen.meanPC2))]

#papi.gen.meanPC1 <- gen.meanPC1[match(papi.gen,names(gen.meanPC1))]
#papi.gen.meanPC2 <- gen.meanPC2[match(papi.gen,names(gen.meanPC2))]

#points(sphing.gen.meanPC1, sphing.gen.meanPC2, col="black", bg="darkblue", pch=21, cex=0.8)
#points(nymph.gen.meanPC1, nymph.gen.meanPC2, col="black",bg="forestgreen", pch=21, cex=0.8)
#points(papi.gen.meanPC1, papi.gen.meanPC2, col="black", bg="firebrick", pch=21, cex=0.8)
#points(ereb.gen.meanPC1, ereb.gen.meanPC2, col="black", bg="darkorange", pch=21, cex=0.8)


#ereb.gen[identify(ereb.gen.meanPC1, ereb.gen.meanPC2,n=3)]
#nymph.gen[c(24, 38, 51, 52, 63, 64, 68, 70)] #"convergent" nymphalidae


###Re-do the plotting with all shapes (MHNT+NHM)
library(MASS)
pca.av.halves <- prcomp(fou.av)
pca.av.halves$sdev^2/sum(pca.av.halves$sdev^2)*100


hull.pca <- chull(pca.av.halves$x[,1],pca.av.halves$x[,2])
hull.o <- sort(hull.pca)
f1 <- kde2d(pca.av.halves$x[,1], pca.av.halves$x[,2], n = 50)

image(f1, xlim=c(-0.3,0.3), ylim=c(-0.2,0.4), zlim = c(0.4, 50), xlab="PC1 (59.4%)", ylab="PC2 (15.6%)")
abline(h=0, col="gray80")
abline(v=0, col="gray80")
grid()
#points(pca.av.halves$x[,1], pca.av.halves$x[,2], cex=0.6, pch=c(17,16)[dfac$side], col=c("forestgreen","darkblue")[dfac$museum])
points(pca.av.halves$x[,1], pca.av.halves$x[,2], cex=0.2, pch=20, col="gray50")
points(pca.av.halves$x[hull.pca,1], pca.av.halves$x[hull.pca,2], cex=0.6, pch=20, col="gray50")
text(pca.av.halves$x[hull.pca,1], pca.av.halves$x[hull.pca,2],labels=hull.pca, pos=3, cex=0.5)

#New plot with average halves x average genus
pca.av.gen <- prcomp(av.gen.fou)
pca.av.gen$sdev^2/sum(pca.av.gen$sdev^2)*100
fam.gen <- unlist(strsplit(sort(unique(paste(dfac$genus, dfac$family))), " "))[seq(from=2, to=714, by=2)]
fam.gen[258] <- "Nymphalidae"
hull.pca <- chull(pca.av.gen$x[,1],pca.av.gen$x[,2])
hull.o <- sort(hull.pca)
f1 <- kde2d(pca.av.gen$x[,1], pca.av.gen$x[,2], n = 50)

image(f1, xlim=c(-0.25,0.3), ylim=c(-0.2,0.3), zlim = c(2, 35), xlab="PC1 (55.4%)", ylab="PC2 (20.1%)")
abline(h=0, col="gray80")
abline(v=0, col="gray80")
grid()
points(pca.av.gen$x[,1], pca.av.gen$x[,2], cex=0.5, pch=20, col="gray50")

big6 <- names(sort(table(fam.gen), decreasing=T))[1:6]

for (i in 1:6) {
points(pca.av.gen$x[which(fam.gen==big6[i]),1], pca.av.gen$x[which(fam.gen==big6[i]),2], cex=1, pch=20, col=i)
chul <- chull(pca.av.gen$x[which(fam.gen==big6[i]),1], pca.av.gen$x[which(fam.gen==big6[i]),2])
lines(pca.av.gen$x[which(fam.gen==big6[i])[c(chul,chul[1])],1], pca.av.gen$x[which(fam.gen==big6[i])[c(chul,chul[1])],2], col=i)
}


library(car)
datell.Nymph <- dataEllipse(pca.av.gen$x[which(fam.gen == "Nymphalidae"),1],pca.av.gen$x[which(fam.gen == "Nymphalidae"),2])

datell.Sphing <- dataEllipse(pca.av.gen$x[which(fam.gen == "Sphingidae"),1],pca.av.gen$x[which(fam.gen == "Sphingidae"),2])

datell.Lyca <- dataEllipse(pca.av.gen$x[which(fam.gen == "Lycaenidae"),1],pca.av.gen$x[which(fam.gen == "Lycaenidae"),2])

datell.Papi <- dataEllipse(pca.av.gen$x[which(fam.gen == "Papilionidae"),1],pca.av.gen$x[which(fam.gen == "Papilionidae"),2], draw=F)

datell.Pier <- dataEllipse(pca.av.gen$x[which(fam.gen == "Pieridae"),1],pca.av.gen$x[which(fam.gen == "Pieridae"),2], draw=F)

datell.Hesp <- dataEllipse(pca.av.gen$x[which(fam.gen == "Hesperiidae"),1],pca.av.gen$x[which(fam.gen == "Hesperiidae"),2], draw=F)

datell.Ereb <- dataEllipse(pca.av.gen$x[which(fam.gen == "Lymantriidae" | fam.gen == "Arctiidae" | fam.gen == "Erebidae"),1], pca.av.gen$x[which(fam.gen == "Lymantriidae" | fam.gen == "Arctiidae" | fam.gen == "Erebidae"),2], draw=F)

points(pca.av.gen$x[which(fam.gen == "Nymphalidae"),1],pca.av.gen$x[which(fam.gen == "Nymphalidae"),2], col="forestgreen", pch=20, cex=1)
lines(datell.Nymph[['0.95']][,1:2], col="forestgreen", lwd=1.5, lty=2)
lines(datell.Nymph[['0.5']][,1:2], col="forestgreen", lwd=2)

points(pca.av.gen$x[which(fam.gen == "Sphingidae"),1],pca.av.gen$x[which(fam.gen == "Sphingidae"),2], col="darkblue", pch=20, cex=1)
lines(datell.Sphing[['0.95']][,1:2], col="darkblue", lwd=1.5, lty=2)
lines(datell.Sphing[['0.5']][,1:2], col="darkblue", lwd=2)

points(pca.av.gen$x[which(fam.gen == "Lycaenidae"),1],pca.av.gen$x[which(fam.gen == "Lycaenidae"),2], col="firebrick", pch=20, cex=1)
lines(datell.Lyca[['0.95']][,1:2], col="firebrick", lwd=1.5, lty=2)
lines(datell.Lyca[['0.5']][,1:2], col="firebrick", lwd=2)

legend("topleft", legend=c("Nymphalidae (N=126)", "Sphingidae (N=98)", "Lycaenidae (N=42)"), pch=20, bty="n", col=c("forestgreen","darkblue", "firebrick"))

points(pca.av.gen$x[which(fam.gen == "Hesperiidae"),1],pca.av.gen$x[which(fam.gen == "Hesperiidae"),2], col="black", pch=20, cex=1.5)
lines(datell.Hesp[['0.5']][,1:2], col="black", lwd=2)

points(pca.av.gen$x[which(fam.gen == "Lymantriidae" | fam.gen == "Arctiidae" | fam.gen == "Erebidae"),1], pca.av.gen$x[which(fam.gen == "Lymantriidae" | fam.gen == "Arctiidae" | fam.gen == "Erebidae"),2], col='red', pch=20)
lines(datell.Ereb[['0.5']][,1:2], col="red", lwd=2)

points(pca.av.gen$x[which(fam.gen == "Papilionidae"),1],pca.av.gen$x[which(fam.gen == "Papilionidae"),2], col="blue", pch=20, cex=1.5)
lines(datell.Papi[['0.5']][,1:2], col="blue", lwd=2)

points(pca.av.gen$x[which(fam.gen == "Pieridae"),1],pca.av.gen$x[which(fam.gen == "Pieridae"),2], col="darkmagenta", pch=20, cex=1.5)
lines(datell.Pier[['0.5']][,1:2], col="darkmagenta", lwd=2)

legend("topleft", legend=c("Nymphalidae (N=126)", "Sphingidae (N=98)", "Lycaenidae (N=42)", "Hesperiidae", "Erebidae","Papilionidae","Pieridae"), pch=20, bty="n", col=c("forestgreen","darkblue", "firebrick", "black", "red","blue","darkmagenta"))

#New plot with average halves x average species
pca.av.sp <- prcomp(av.sp.fou)
pca.av.sp$sdev^2/sum(pca.av.sp$sdev^2)*100
fam.sp <- unlist(strsplit(sort(unique(paste(dfac$species, dfac$family))), " "))[seq(from=3, to=2553, by=3)]
fam.sp[612] <- "Nymphalidae"
hull.pca <- chull(pca.av.sp$x[,1],pca.av.sp$x[,2])
hull.o <- sort(hull.pca)
f1 <- kde2d(pca.av.sp$x[,1], pca.av.sp$x[,2], n = 50)

image(f1, xlim=c(-0.25,0.3), ylim=c(-0.2,0.3), zlim = c(2, 59), xlab="PC1 (57.4%)", ylab="PC2 (17.0%)")
abline(h=0, col="gray80")
abline(v=0, col="gray80")
grid()
points(pca.av.sp$x[,1], pca.av.sp$x[,2], cex=0.5, pch=20, col="gray50")

big6 <- names(sort(table(fam.sp), decreasing=T))[1:6]

for (i in 1:6) {
points(pca.av.sp$x[which(fam.sp==big6[i]),1], pca.av.sp$x[which(fam.sp==big6[i]),2], cex=1, pch=20, col=i)
chul <- chull(pca.av.sp$x[which(fam.sp==big6[i]),1], pca.av.sp$x[which(fam.sp==big6[i]),2])
lines(pca.av.sp$x[which(fam.sp==big6[i])[c(chul,chul[1])],1], pca.av.sp$x[which(fam.sp==big6[i])[c(chul,chul[1])],2], col=i)
}


library(car)
datell.Nymph <- dataEllipse(pca.av.sp$x[which(fam.sp == "Nymphalidae"),1],pca.av.sp$x[which(fam.sp == "Nymphalidae"),2])

datell.Sphing <- dataEllipse(pca.av.sp$x[which(fam.sp == "Sphingidae"),1],pca.av.sp$x[which(fam.sp == "Sphingidae"),2])

datell.Lyca <- dataEllipse(pca.av.sp$x[which(fam.sp == "Lycaenidae"),1],pca.av.sp$x[which(fam.sp == "Lycaenidae"),2])

datell.Papi <- dataEllipse(pca.av.sp$x[which(fam.sp == "Papilionidae"),1],pca.av.sp$x[which(fam.sp == "Papilionidae"),2], draw=F)

datell.Pier <- dataEllipse(pca.av.sp$x[which(fam.sp == "Pieridae"),1],pca.av.sp$x[which(fam.sp == "Pieridae"),2], draw=F)

datell.Hesp <- dataEllipse(pca.av.sp$x[which(fam.sp == "Hesperiidae"),1],pca.av.sp$x[which(fam.sp == "Hesperiidae"),2], draw=F)

datell.Ereb <- dataEllipse(pca.av.sp$x[which(fam.sp == "Lymantriidae" | fam.sp == "Arctiidae" | fam.sp == "Erebidae"),1], pca.av.sp$x[which(fam.sp == "Lymantriidae" | fam.sp == "Arctiidae" | fam.sp == "Erebidae"),2], draw=F)

points(pca.av.sp$x[which(fam.sp == "Nymphalidae"),1],pca.av.sp$x[which(fam.sp == "Nymphalidae"),2], col="forestgreen", pch=20, cex=1)
lines(datell.Nymph[['0.5']][,1:2], col="forestgreen", lwd=2)

points(pca.av.sp$x[which(fam.sp == "Sphingidae"),1],pca.av.sp$x[which(fam.sp == "Sphingidae"),2], col="darkblue", pch=20, cex=1)
lines(datell.Sphing[['0.5']][,1:2], col="darkblue", lwd=2)

points(pca.av.sp$x[which(fam.sp == "Lycaenidae"),1],pca.av.sp$x[which(fam.sp == "Lycaenidae"),2], col="firebrick", pch=20, cex=1)
lines(datell.Lyca[['0.5']][,1:2], col="firebrick", lwd=2)

points(pca.av.sp$x[which(fam.sp == "Hesperiidae"),1],pca.av.sp$x[which(fam.sp == "Hesperiidae"),2], col="black", pch=20, cex=1)
lines(datell.Hesp[['0.5']][,1:2], col="black", lwd=2)

points(pca.av.sp$x[which(fam.sp == "Lymantriidae" | fam.sp == "Arctiidae" | fam.sp == "Erebidae"),1], pca.av.sp$x[which(fam.sp == "Lymantriidae" | fam.sp == "Arctiidae" | fam.sp == "Erebidae"),2], col='red', pch=20)
lines(datell.Ereb[['0.5']][,1:2], col="red", lwd=2)

points(pca.av.sp$x[which(fam.sp == "Papilionidae"),1],pca.av.sp$x[which(fam.sp == "Papilionidae"),2], col="blue", pch=20, cex=1)
lines(datell.Papi[['0.5']][,1:2], col="blue", lwd=2)

points(pca.av.sp$x[which(fam.sp == "Pieridae"),1],pca.av.sp$x[which(fam.sp == "Pieridae"),2], col="darkmagenta", pch=20, cex=1)
lines(datell.Pier[['0.5']][,1:2], col="darkmagenta", lwd=2)

legend("topleft", legend=c("Nymphalidae (N=255)", "Sphingidae (N=391)", "Lycaenidae (N=59)", "Hesperiidae (N=16)", "Erebidae (N=14)","Papilionidae (N=53)","Pieridae (N=15)"), pch=20, bty="n", col=c("forestgreen","darkblue", "firebrick", "black", "red","blue","darkmagenta"))

#Plot with insets

pdf(file="ACP_insets_ellipses.pdf", width=10, height=10)

layout(matrix(c(1,2,3,5, 4,4,4, 6, 4,4,4, 7, 4,4,4,8), nrow=4, ncol=4, byrow=T))
par(mar=c(1,1,1,1), xaxt="n", yaxt="n")

image(f1, xlim=c(-0.2,0.3), ylim=c(-0.2,0.3), zlim = c(2, 59), xlab="",  ylab="")
abline(h=0, col="gray80")
abline(v=0, col="gray80")
grid()
points(pca.av.sp$x[,1], pca.av.sp$x[,2], cex=0.5, pch=20, col="gray50")

points(pca.av.sp$x[which(fam.sp == "Nymphalidae"),1],pca.av.sp$x[which(fam.sp == "Nymphalidae"),2], col="forestgreen", pch=20, cex=1)
lines(datell.Nymph[['0.5']][,1:2], col="forestgreen", lwd=2)
lines(datell.Nymph[['0.95']][,1:2], col="forestgreen", lwd=2, lty=2)
legend("topleft", legend="Nymphalidae", bty="n")

image(f1, xlim=c(-0.2,0.3), ylim=c(-0.2,0.3), zlim = c(2, 59), xlab="",  ylab="")
abline(h=0, col="gray80")
abline(v=0, col="gray80")
grid()
points(pca.av.sp$x[,1], pca.av.sp$x[,2], cex=0.5, pch=20, col="gray50")

points(pca.av.sp$x[which(fam.sp == "Sphingidae"),1],pca.av.sp$x[which(fam.sp == "Sphingidae"),2], col="darkblue", pch=20, cex=1)
lines(datell.Sphing[['0.5']][,1:2], col="darkblue", lwd=2)
lines(datell.Sphing[['0.95']][,1:2], col="darkblue", lwd=2, lty=2)
legend("topleft", legend="Sphingidae", bty="n")

image(f1, xlim=c(-0.2,0.3), ylim=c(-0.2,0.3), zlim = c(2, 59), xlab="",  ylab="")
abline(h=0, col="gray80")
abline(v=0, col="gray80")
grid()
points(pca.av.sp$x[,1], pca.av.sp$x[,2], cex=0.5, pch=20, col="gray50")

points(pca.av.sp$x[which(fam.sp == "Lycaenidae"),1],pca.av.sp$x[which(fam.sp == "Lycaenidae"),2], col="firebrick", pch=20, cex=1)
lines(datell.Lyca[['0.5']][,1:2], col="firebrick", lwd=2)
lines(datell.Lyca[['0.95']][,1:2], col="firebrick", lwd=2, lty=2)
legend("topleft", legend="Lycaenidae", bty="n")

par(mar=c(5,5,2,2), xaxt="s", yaxt="s")
image(f1, xlim=c(-0.2,0.3), ylim=c(-0.2,0.3), zlim = c(2, 59), xlab="PC1 (57.4%)", ylab="PC2 (17.0%)", cex.lab=1.5)
abline(h=0, col="gray80")
abline(v=0, col="gray80")
grid()
points(pca.av.sp$x[,1], pca.av.sp$x[,2], cex=1, pch=20, col="black")

par(mar=c(1,1,1,1), xaxt="n", yaxt="n")
image(f1, xlim=c(-0.2,0.3), ylim=c(-0.2,0.3), zlim = c(2, 59), xlab="",  ylab="")
abline(h=0, col="gray80")
abline(v=0, col="gray80")
grid()
points(pca.av.sp$x[,1], pca.av.sp$x[,2], cex=0.5, pch=20, col="gray50")

points(pca.av.sp$x[which(fam.sp == "Hesperiidae"),1],pca.av.sp$x[which(fam.sp == "Hesperiidae"),2], col="black", pch=20, cex=1)
lines(datell.Hesp[['0.5']][,1:2], col="black", lwd=2)
lines(datell.Hesp[['0.95']][,1:2], col="black", lwd=2, lty=2)
legend("topleft", legend="Hesperiidae", bty="n")

image(f1, xlim=c(-0.2,0.3), ylim=c(-0.2,0.3), zlim = c(2, 59), xlab="",  ylab="")
abline(h=0, col="gray80")
abline(v=0, col="gray80")
grid()
points(pca.av.sp$x[,1], pca.av.sp$x[,2], cex=0.5, pch=20, col="gray50")

points(pca.av.sp$x[which(fam.sp == "Lymantriidae" | fam.sp == "Arctiidae" | fam.sp == "Erebidae"),1], pca.av.sp$x[which(fam.sp == "Lymantriidae" | fam.sp == "Arctiidae" | fam.sp == "Erebidae"),2], col='red', pch=20)
lines(datell.Ereb[['0.5']][,1:2], col="red", lwd=2)
lines(datell.Ereb[['0.95']][,1:2], col="red", lwd=2, lty=2)
legend("topleft", legend="Erebidae", bty="n")

image(f1, xlim=c(-0.2,0.3), ylim=c(-0.2,0.3), zlim = c(2, 59), xlab="",  ylab="")
abline(h=0, col="gray80")
abline(v=0, col="gray80")
grid()
points(pca.av.sp$x[,1], pca.av.sp$x[,2], cex=0.5, pch=20, col="gray50")

points(pca.av.sp$x[which(fam.sp == "Papilionidae"),1],pca.av.sp$x[which(fam.sp == "Papilionidae"),2], col="blue", pch=20, cex=1)
lines(datell.Papi[['0.5']][,1:2], col="blue", lwd=2)
lines(datell.Papi[['0.95']][,1:2], col="blue", lwd=2, lty=2)
legend("topleft", legend="Papilionidae", bty="n")

image(f1, xlim=c(-0.2,0.3), ylim=c(-0.2,0.3), zlim = c(2, 59), xlab="", ylab="")
abline(h=0, col="gray80")
abline(v=0, col="gray80")
grid()
points(pca.av.sp$x[,1], pca.av.sp$x[,2], cex=0.5, pch=20, col="gray50")

points(pca.av.sp$x[which(fam.sp == "Pieridae"),1],pca.av.sp$x[which(fam.sp == "Pieridae"),2], col="darkmagenta", pch=20, cex=1)
lines(datell.Pier[['0.5']][,1:2], col="darkmagenta", lwd=2)
lines(datell.Pier[['0.95']][,1:2], col="darkmagenta", lwd=2, lty=2)
legend("topleft", legend="Pieridae", bty="n")

dev.off()

#Extreme species to illustrate morphospace variation
image(f1, xlim=c(-0.23,0.36), ylim=c(-0.2,0.3), zlim = c(2, 59), xlab="",  ylab="")
abline(h=0, col="gray80")
abline(v=0, col="gray80")
grid()
points(pca.av.sp$x[which(!butterfly.sp),1], pca.av.sp$x[which(!butterfly.sp),2], cex=0.7, pch=16, col=alpha("black",0.75))
points(pca.av.sp$x[which(butterfly.sp),1], pca.av.sp$x[which(butterfly.sp),2], cex=0.7, pch=17, col=alpha("gray50",0.75))
legend("topleft", pch=c(17,16), cex=1, legend=c("'Butterflies'", "'Moths'"), col=c("gray50", "black"))
points(pca.av.sp$x[chull(pca.av.sp$x[,1:2])[-c(3,6)],1], pca.av.sp$x[chull(pca.av.sp$x[,1:2])[-c(3,6)],2], cex=1.5, pch=0, col="black")
text(pca.av.sp$x[chull(pca.av.sp$x[,1:2])[-c(3,6)],1], pca.av.sp$x[chull(pca.av.sp$x[,1:2])[-c(3,6)],2], cex=0.7, labels=dfac.sp$species[chull(pca.av.sp$x[,1:2])[-c(3,6)]], col="black", pos=1)

dfac.sp$species[chull(pca.av.sp$x[,1:2])]
#Marpesia orsilochus 991      Eunica amelia 578        
#Symmachia miron 3020        Zela steineri 2392           
#Oleria ilerdina 1156          Cautethia spuria 346       
#Zygaena transalpina 1714       Protambulyx goeldii  1313    
#Pterophorus pentadactyla 1335  Stichophthalma devyatkini 1825
#Urania fulgens 1538

coo.plot<-array(NA, dim=c(2,2,9))
for (i in 1:9) {
xy1 <- grid.locator(unit="npc")
coo.plot[1,1,i] <- as.numeric(substr(xy1$x, 1, 5))
coo.plot[1,2,i] <- as.numeric(substr(xy1$y, 1, 5))
coo.plot[2,,i] <- coo.plot[1,,i] + 0.08
}

hull.pca <- c(991,578,3020,2392,1156,346,1714,1313,1335,1825,1538)
#Auto-add shape sub-plots
for (i in 1:11) {
par(new=T, fig=c(coo.plot[,,i]), mar=c(0,0,0,0), xaxt="n", yaxt="n")
n <- hull.pca[i]
m <- po.butterflies[n]
plot(m[,1], m[,2], bty="n", type="l", asp=1)
polygon(m[,1], m[,2], col="grey70")
}
hull.pca <- c(991,578,2392,1156,1714,1313,1335,1825,1538)
coo.plot <- coo.plot[,,-c(3,6)]
for (i in 1:9) {
par(new=T, fig=c(coo.plot[,,i]), mar=c(0,0,0,0), xaxt="n", yaxt="n")
n <- hull.pca[i]
m <- po.butterflies[n]
plot(m[,1], m[,2], bty="n", type="l", asp=1)
polygon(m[,1], m[,2], col="grey70")
}

#Test of effects
dfac.sp <- unique(dfac[,c("family","species")])
dfac.sp <- dfac.sp[order(dfac.sp$species),]
ooo <- unique(dfac[,c("museum","species")]) #WARNING: A handful of species are present in both museums. For simplicity the first museum appearing in dfac.av was chosen.
dfac.sp <- cbind(dfac.sp, ooo$museum[match(dfac.sp$species, ooo$species)])
colnames(dfac.sp)[3] <- "museum"

sum(pca.av.sp$sdev[1:6]^2)/sum(pca.av.sp$sdev^2)*100 #91.7% of variance

Manova(lm(pca.av.sp$x[,1:6] ~ dfac.sp$family + dfac.sp$museum)) #No interaction due to singularity
#Type II MANOVA Tests: Pillai test statistic
#               Df test stat approx F num Df den Df    Pr(>F)    
#dfac.sp$family 23   1.68297   13.984    138   4950 < 2.2e-16 ***
#dfac.sp$museum  1   0.08739   13.087      6    820 3.547e-14 ***

etasq(Manova(lm(pca.av.sp$x[,1:6] ~ dfac.sp$family + dfac.sp$museum)))
#                    eta^2
#dfac.sp$family 0.28049560
#dfac.sp$museum 0.08738859
#Taxonomic signal is 3 times more important than museum bias

#Difference on PC1
fam.sp2 <- fam.sp
fam.sp2[which(fam.sp == "Lymantriidae" | fam.sp == "Arctiidae" | fam.sp == "Erebidae")] <- "Erebidae"
big7 <- c(big6, "Erebidae")

par(mar=c(7,5,1,1))
boxplot(pca.av.sp$x[which(fam.sp %in% big6),1] ~ fam.sp[which(fam.sp %in% big6)], las=2, ylab="PC1 coordinate", xlab="", widht=50, col=c("gray20", "firebrick", "forestgreen", "blue", "darkmagenta", "darkblue"))

TukeyHSD(aov(pca.av.sp$x[which(fam.sp2 %in% big7),1] ~ fam.sp2[which(fam.sp2 %in% big7)])) #Differences in mean position in PC1 are significant for many families.

#Disparity calculation

SOV <- SOR <- rep(NA,6)

for (i in 1:6) {
SOV[i] <- sum(diag(var(pca.av.gen$x[which(fam.gen == big6[i]),]))) #A verifier
SOR[i] <- sum(diff(apply(pca.av.gen$x[which(fam.gen == big6[i]),], 2, range)))
}
names(SOV) <- names(SOR) <- big6

SOVb <- SORb <- matrix(NA, nrow=1000, ncol=6)
for (j in 1:1000) {

   for (i in 1:6) {
   SOVb[j,i] <- sum(diag(var(pca.av.gen$x[which(fam.gen == big6[i])[sample.int(n=11, replace=T)],])))
   SORb[j,i] <- sum(diff(apply(pca.av.gen$x[which(fam.gen == big6[i])[sample.int(n=11, replace=T)],], 2, range)))
   }

}

par(mar=c(7,3,3,2))
layout(matrix(c(1:2), nrow=1, ncol=2))
boxplot(SOVb, names=big6, las=2, main="Sum of Variance")
points(1:6, SOV, col="red", pch=20, cex=3)
boxplot(SODb, names=big6, las=2, main="Sum of Ranges")
points(1:6, SOD, col="red", pch=20, cex=3)

###For species
SOV <- SOR <- rep(NA,7)

for (i in 1:7) {
SOV[i] <- sum(diag(var(pca.av.sp$x[which(fam.sp == big7[i]),])))
SOR[i] <- sum(apply(pca.av.sp$x[which(fam.sp == big7[i]),], 2, max)-apply(pca.av.sp$x[which(fam.sp == big7[i]),], 2, min))
}
names(SOV) <- names(SOR) <- big7

SOVb <- SORb <- matrix(NA, nrow=1000, ncol=7)
for (j in 1:1000) {

   for (i in 1:7) {
   SOVb[j,i] <- sum(diag(var(pca.av.sp$x[which(fam.sp2 == big7[i]),][sample.int(n=14, replace=T),])))
   o <- sample.int(n=14, replace=T)
   SORb[j,i] <- sum(apply(pca.av.sp$x[which(fam.sp2 == big7[i]),][o,], 2, max)-apply(pca.av.sp$x[which(fam.sp2 == big7[i]),][o,], 2, min))
   }

}

par(mar=c(7,3,3,2))
layout(matrix(c(1:2), nrow=1, ncol=2))
boxplot(SOVb, names=big7, las=2, main="Sum of Variance")
points(1:7, SOV, col="red", pch=20, cex=3)
boxplot(SORb, names=big7, las=2, main="Sum of Ranges")
points(1:7, SOR, col="red", pch=20, cex=3)

pdf(file="violin_plots.pdf", width=10, height=3.5)
par(mar=c(7,3,3,1))
layout(matrix(c(1:3), nrow=1, ncol=3))
vioplot(pca.av.sp$x[which(fam.sp2 %in% big7),1] ~ fam.sp2[which(fam.sp2 %in% big7)], las=2, main="Position on PC1", ylab="", xlab="", widht=50, col=c("red", "gray20", "firebrick", "forestgreen", "blue", "darkmagenta", "darkblue"))
vioplot(c(SOVb) ~ rep(big7, each=1000), las=2, main="Sum of Variances", col=c("red", "gray20", "firebrick", "forestgreen", "blue", "darkmagenta", "darkblue"), xlab="", ylab="")
vioplot(c(SORb) ~ rep(big7, each=1000), las=2, main="Sum of Ranges", col=c("red", "gray20", "firebrick", "forestgreen", "blue", "darkmagenta", "darkblue"), xlab="", ylab="")
dev.off()

TukeyHSD(aov(c(SOVb) ~ rep(big7, each=1000))) #Test shows that all bootstraped Sum of variance are different except Papilionidae vs Nymphalidae.
TukeyHSD(aov(c(SORb) ~ rep(big7, each=1000))) #Test shows that all bootstraped Sum of ranges are different except Papilionidae vs Lycaenidae.

###Other thing

#Use it before superimposition
#which.min(extr$shapes[[1]][,1]) #1
#which.max(extr$shapes[[1]][,1]) #130
#mean(extr$shapes[[1]][c(1,130),1]) #1946
#left <- extr$shapes[[1]][which(extr$shapes[[1]][,1] <= 1946),]
#right <- extr$shapes[[1]][which(extr$shapes[[1]][,1] > 1946),]

#rot.left <- cbind(-1*(left[,1]-1946), left[,2])
#trans.right <- cbind(right[,1]-1946, right[,2])
#rli <- coo_interpolate(rot.left, 100)
#ri <- coo_interpolate(trans.right, 100)

#plot(rli, asp=1)
#points(ri, col="red")
#which.max(rli[,1]) #1
#which.min(dist(rbind(rli[1,], ri))[1:101]) #52
#ri.ord <- ri[c(52:1,100:53),]
#aproc<-pPsup(rli, ri.ord)
#ar<-array(data=c(aproc$Mp1,aproc$Mp2), dim=c(100,2,2))
#mshape(ar)


###New thing to test for convergence of Nymphalidae towards Sphingidae###

#Bootstrap estimate of mean along PC1/PC2, to test if some butterflies are out of the intervals

#shapiro.test(pca.av.sp$x[,1])

#	Shapiro-Wilk normality test
#data:  pca.av.sp$x[, 1]
#W = 0.90392, p-value < 2.2e-16

#shapiro.test(pca.av.sp$x[,2])

#	Shapiro-Wilk normality test
#data:  pca.av.sp$x[, 2]
#W = 0.98166, p-value = 7.53e-09

#Non-normal distro so let's use non parametric test => bootstrap

#bmeanPC1 <- bmeanPC2 <- rep(NA, 1000)

#for (i in 1:1000) {
#bmeanPC1[i] <- min(sample(pca.av.sp$x[which(fam.sp == "Nymphalidae"),1], replace=T))
#bmeanPC2[i] <- min(sample(pca.av.sp$x[which(fam.sp == "Nymphalidae"),2], replace=T))
#}

#meanPC1 <- mean(pca.av.sp$x[which(fam.sp == "Nymphalidae"),1])
#meanPC2 <- mean(pca.av.sp$x[which(fam.sp == "Nymphalidae"),2])


#sphin_bmeanPC1 <- sphin_bmeanPC2 <- rep(NA, 1000)

#for (i in 1:1000) {
#sphin_bmeanPC1[i] <- min(sample(pca.av.sp$x[which(fam.sp == "Sphingidae"),1], replace=T))
#sphin_bmeanPC2[i] <- min(sample(pca.av.sp$x[which(fam.sp == "Sphingidae"),2], replace=T))
#}

#sphin_meanPC1 <- mean(pca.av.sp$x[which(fam.sp == "Sphingidae"),1])
#sphin_meanPC2 <- mean(pca.av.sp$x[which(fam.sp == "Sphingidae"),2])

#CI_max_PC1 <- sort(bmeanPC1)[c(25,975)]
#CI_max_PC2 <- sort(bmeanPC2)[c(25,975)]

#image(f1, xlim=c(-0.2,0.3), ylim=c(-0.2,0.3), zlim = c(2, 59), xlab="",  ylab="")
#abline(h=0, col="gray80")
#abline(v=0, col="gray80")
#grid()
#points(pca.av.sp$x[,1], pca.av.sp$x[,2], cex=0.5, pch=20, col="gray50")

#points(pca.av.sp$x[which(fam.sp == "Nymphalidae"),1],pca.av.sp$x[which(fam.sp == "Nymphalidae"),2], col="forestgreen", pch=20, cex=1)

#points(pca.av.sp$x[which(fam.sp == "Sphingidae"),1],pca.av.sp$x[which(fam.sp == "Sphingidae"),2], col="darkblue", pch=20, cex=1)
#points(meanPC1,meanPC2, cex=4, pch=21)
#points(meanPC1,meanPC2, cex=4, pch=22)
#points(meanPC1,meanPC2, cex=4, pch=18)
#points(meanPC1,meanPC2, cex=4, pch=21, col="forestgreen")
#points(meanPC1,meanPC2, cex=4, pch=20, col="forestgreen")
#points(sphin_meanPC1,sphin_meanPC2, cex=4, pch=20, col="darkblue")

###Discriminant analysis to check if some Nymphalids are misidentified as Sphingids (for PC1 and PC2)
two.groups <- fam.sp[which(fam.sp=="Nymphalidae" | fam.sp=="Sphingidae")]
dis <- lda(pca.av.sp$x[which(fam.sp=="Nymphalidae" | fam.sp=="Sphingidae"),1:9], grouping=two.groups)
pred <- predict(dis, pca.av.sp$x[which(fam.sp=="Nymphalidae" | fam.sp=="Sphingidae"),1:9])
dis_CV <- lda(pca.av.sp$x[which(fam.sp=="Nymphalidae" | fam.sp=="Sphingidae"),1:9], grouping=two.groups, CV=T)
two.groups[which(dis_CV$class!=two.groups)] #Wrong identificatins are only Nymphalidae for PC1-PC2. For PC1-PC9 also, but misidentifications are less numerous

discri <- lda(pca.av.sp$x[,1:2], grouping=fam.sp2)
predis <- predict(discri, pca.av.sp$x[,1:2])

discri_CV <- lda(pca.av.sp$x[,1:2], grouping=fam.sp2, CV=T)
length(discri_CV$class)
#851
length(which(discri_CV$class==fam.sp2))
#603
#0.71 of correct identification with cross validation

discri_CV$class[which(fam.sp2=="Sphingidae")]=="Sphingidae"
#All Sphingidae are correctly identified

droplevels(discri_CV$class[which(fam.sp2=="Nymphalidae")])
#Some Nymphalidae are misidentified as Lycaenidae or Sphingidae

table(droplevels(discri_CV$class[which(fam.sp2=="Nymphalidae")]))
# Lycaenidae Nymphalidae  Sphingidae 
#         27         177          51 
#      0.106       0.694         0.2

table(droplevels(discri_CV$class[which(fam.sp2=="Lycaenidae")]))
#  Endromidae   Lycaenidae  Nymphalidae Papilionidae 
#           1           33           21            4 

table(droplevels(discri_CV$class[which(fam.sp2=="Erebidae")]))
#Nymphalidae Papilionidae   Sphingidae 
#          6            1            7 

nymph_as_sphin <- which(discri_CV$class[which(fam.sp2=="Nymphalidae")]=='Sphingidae')

points(pca.av.sp$x[which(fam.sp2 == "Nymphalidae")[nymph_as_sphin],1],pca.av.sp$x[which(fam.sp2 == "Nymphalidae")[nymph_as_sphin],2], col="forestgreen", pch=20, cex=1) #To add those convergent points to the PCA graph

discri$posterior[which(fam.sp2 == "Nymphalidae")[nymph_as_sphin], c("Nymphalidae","Sphingidae")]
#To find which are the taxa concerned, and the posterior probabilities assigned.


###Defining groups of 'butterflies' and 'moths'
butterfly.fam <- c("Papilionidae", "Hesperidae", "Pieridae", "Lycaenidae", "Riodinidae", "Nymphalidae")
butterfly.dfac <- rep(F,dim(dfac)[1])
butterfly.dfac[which(dfac$family %in% butterfly.fam)] <- T

butterfly.av <- rep(F,dim(dfac.av)[1])
butterfly.av[which(dfac.av$family %in% butterfly.fam)] <- T

butterfly.sp <- rep(F,dim(dfac.sp)[1])
butterfly.sp[which(dfac.sp$family %in% butterfly.fam)] <- T

#Butterfly vs Moth PCA
image(f1, xlim=c(-0.23,0.36), ylim=c(-0.2,0.3), zlim = c(2, 59), xlab="",  ylab="")
abline(h=0, col="gray80")
abline(v=0, col="gray80")
grid()
points(pca.av.sp$x[which(!butterfly.sp),1], pca.av.sp$x[which(!butterfly.sp),2], cex=0.7, pch=16, col=alpha("black",0.75))
points(pca.av.sp$x[which(butterfly.sp),1], pca.av.sp$x[which(butterfly.sp),2], cex=0.7, pch=17, col=alpha("gray50",0.75))
legend("topleft", pch=c(17,16), cex=1, legend=c("'Butterflies'", "'Moths'"), col=c("gray50", "black"))

which(pca.av.sp$x[,1] > 0 & !butterfly.sp) #Moths that are on positive side of PC1 -> convergent towards butterflies.

table(dfac.sp[which(pca.av.sp$x[,1] > 0 & !butterfly.sp),1]) #They are part of  Castniidae (2), Endromidae (1), Erebidae(3), Geometridae(5), Lacturidae(1), Lasiocampidae(9), Lymantriidae(7), Notodontidae(1), Pterophoridae(1), Saturnidae(5), Sphingidae(3), Uraniidae(3)

dfac.sp$species[c(1,35,53,106, 195, 300, 321, 328, 343, 421, 440, 441, 456, 462, 522, 684, 711, 717, 718, 751, 778, 790, 791, 795)] #The most butterfly like
# [1] Abraxas grossulariata    Aglia tau                Aloeides thyra          
# [4] Arctornis l-nigrum       Chrysiridia rhipheus     Euproctis chrysorrhoea  
#[7] Ginungagapus awarreni    Graellsia isabellae      Hasora hurama           
#[10] Lasiocampa quercus       Lymantria dispar         Lymantria monacha       
#[13] Macrothylacia rubi       Malacosoma alpicola      Menophra abruptaria     
#[16] Pseudopanthera macularia Rothschildia aurota      Saturnia pavonia        
#[19] Saturnia pyri            Telchin licus            Thyatira batis          
#[22] Urania fulgens           Urania leilus            Xanthocastnia evalthe   

buttormoth <- as.numeric(butterfly.sp)
discri_moth <- lda(pca.av.sp$x[,1:6], grouping=buttormoth, CV=T)
sum(discri_moth$class==buttormoth) / 851 *100 #89% of correct assignments for PC1-PC2, 92.8% for PC1-PC6
dfac.sp[!discri_moth$class==buttormoth,]
sort(table(dfac.sp[!discri_moth$class==buttormoth,1])) #Bad assignements :
#Drepanidae    Lacturidae    Lycaenidae     Noctuidae    Castniidae 
#         1             1             1             1             2 
#  Erebidae  Papilionidae      Pieridae     Uraniidae Lasiocampidae 
#         2             2             3             3             4 
#Geometridae   Saturniidae  Lymantriidae   Nymphalidae 
#         5             5             6            24 
#24+2+1+3=30 Half are butterflies, the other half are moths

#calculate body length / wingspan

BL <- WS <- rep(NA, length(p.butterflies))
for (i in 1:length(p.butterflies)) {
p1 <- p.butterflies$ldk[[i]][1]
p2 <-  p.butterflies$ldk[[i]][2]
p3 <-  p.butterflies$ldk[[i]][3]
p4 <-  p.butterflies$ldk[[i]][4]

WS[i] <- dist(p.butterflies[i][[1]][c(p1,p2),])
BL[i] <- dist(p.butterflies[i][[1]][c(p3,p4),])
}

ratioWSBL <- WS/BL

##################################################################
####### New run with data from Morpho collection of V Debat ######
##################################################################

output_shapes <- auto.LM.cont(x=list.files(), NLM=200, tol=0.2)
bugged <- which(unlist(lapply(output_shapes$shapes, is.null)))
length(bugged)/2545 #Lost 15% of pictures already
ar.shapes <- array(NA, dim=c(200,2,length(output_shapes$shapes)))
for (i in 1:length(output_shapes$shapes)){
tryCatch({ar.shapes[,,i] <- output_shapes$shapes[[i]]}, error=function(e){})}
ar.shapes <- ar.shapes[,,-bugged] #Remove NULL/NA from array
Out.butterflies <- Out(ar.shapes)
save.image("morpho_debat.RData")

pdf("bigpanel.pdf", width=15, height=15)
panel(Out.butterflies, names=T)
dev.off()

list_badshape <- list(NULL)

for (i in 1:dim(ar.shapes)[3]) {

plot(ar.shapes[,,i], asp=1, type="l")
list_badshape[[i]] <- locator()

}

lapply(list_badshape, is.null)
goodshape<-which(unlist(lapply(list_badshape, is.null)))
ar.goodshapes <- array(NA, dim=c(200,2,length(goodshape)))

for (i in 1:length(goodshape)) {ar.goodshapes[,,i] <- ar.shapes[,, goodshape[i]]}
#Missed a few, have to remove them

ar.goodshapes <- ar.goodshapes[,,-c(101,256, 540,565, 834, 1076)]
Out.goodshapes <- Out(ar.goodshapes[,,-c(101,256, 540,565, 834, 1076)])
panel(Out.goodshapes)

a <- list.files("jpgs/")[-bugged]
aa <- a[goodshape]
aaa <- gsub("PB13.", "", aa)
aaaa<-substr(aaa, 1,4)
DV <- as.factor(substr(aaa, 6,6))[-c(101,256, 540,565, 834, 1076)]
indi <- as.factor(aaaa[-c(101,256, 540,565, 834, 1076)])


ldk.ls <- vector("list", length=dim(ar.goodshapes)[3])

for (i in 1:dim(ar.goodshapes)[3]){

shp <- ar.goodshapes[,,i]
p1 <- which.max(apply(shp,1,mean)) #Tip of right anterior wing. Mean is useful, because hindwings can sometimes have a larger x value. But combining x and y make it sure that we catch the correct point

shp.cent <- shp
shp.cent[,1] <- shp[,1]-mean(shp[,1])
left <- shp.cent[which(shp.cent[,1] < 0),]

p2.left <- which.max(apply(abs(left),1,mean)) #Equivalent to p1, but abs() mirrors the shape
p2 <- which(shp.cent[,1]==left[p2.left,1] & shp.cent[,2]==left[p2.left,2])

maxx <- max(shp[,1]) #Max x value
minx <- min(shp[,1]) #Min x value
midx <- mean(c(maxx,minx))
meany <- mean(shp[,2])
lower.h <- shp[which(shp[,2] < meany),] #Lower half of shape
upper.h <- shp[which(shp[,2] > meany),] #Upper half of shape

lower.tip <- which.min(abs(lower.h[,1]-midx)) #Difference between x coordinate and mean x allows to find the point closest to the symmetry axis in terms of x.
upper.tip <- which.min(abs(upper.h[,1]-midx))

p3 <- which(shp[,1] == lower.h[lower.tip,1] & shp[,2] == lower.h[lower.tip,2])
p4 <- which(shp[,1] == upper.h[upper.tip,1] & shp[,2] == upper.h[upper.tip,2])

ldk.ls[[i]] <- c(p1,p2,p3,p4)
}


Out.goodshapes$ldk <- ldk.ls

p.butterflies <- fgProcrustes(Out.goodshapes)

stack(p.butterflies) #Fair results...

### Reordering of points, starting from the anterior-right tip of the right wing ###
po.butterflies <- array(NA, dim=dim(ar.goodshapes))

for (i in 1:dim(po.butterflies)[3]) {
start <- ldk.ls[[i]][1] #This is p1 from previous loop
m <- p.butterflies[i]
mo <- m[[1]][c(start:200, 1:(start-1)),]
po.butterflies[,,i] <- mo
}

po.butterflies <- Out(po.butterflies[,,-c(7,8,17,1, 216, 367, 401, 418, 424, 434, 509, 826, 829)]) #RE-convert to momocs object

stack(po.butterflies)

fou.po <- efourier(po.butterflies, norm=F) #Fourier analysis

pca.fou <- PCA(fou.po)
plot(pca.fou)
plot(pca.fou$x[,1:2], col=c(1,2)[dv])
#Probleme de symmétrie...?

ind <- droplevels(indi[-c(7,8,17,1, 216, 367, 401, 418, 424, 434, 509, 826, 829)])
dv <- DV[-c(7,8,17,1, 216, 367, 401, 418, 424, 434, 509, 826, 829)]

av.ind <- matrix(NA, length(levels(ind)), dim(fou.po$coe)[2])

for (i in 1:dim(av.ind)[2]) {av.ind[,i] <- tapply(fou.po$coe[,i], ind, mean)}

pca.av.ind <- prcomp(av.ind)
plot(pca.av.ind$x[,1:2])
# Not enough, let's try to mirror ventral/dorsal?


stack(Out(ar.goodshapes[,,which(DV=="v")]))
x11()
stack(Out(ar.goodshapes[,,which(DV=="d")]))

po.butterflies.rev <- po.butterflies$coo

for (i in 1:length(po.butterflies.rev)) {
 if (dv[i] == "d") {
 po.butterflies.rev[[i]][,1] <- -po.butterflies.rev[[i]][,1]
 }
}

Out.rev <- Out(po.butterflies.rev, fac=data.frame(dv=dv, ind=ind))
fou.rev <- efourier(Out.rev, norm=F) #Fourier analysis

po.v <- Out(po.butterflies$coo[which(dv=="v")])
po.d <- Out(po.butterflies$coo[which(dv=="d")])
fou.v <- efourier(po.v, norm=F)
fou.d <- efourier(po.d, norm=F)

plot(PCA(fou.v))
plot(PCA(fou.d))
plot(PCA(fou.rev))

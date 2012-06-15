#Scaling
scale(x, center = TRUE, scale = TRUE)

#Diagramme en bâton

plot(table(dataframe[,"TotalDeMots"]), ylab="Number of individuals",xlab="number of words", main="Barplot for the number of words", yaxt="n")

#Mise en place des unités sur l'axe y
axis(2, at =  1:3)


#Fonction de la distribution empirique discrète:
     
plot(ecdf(dataframe[,"TotalDeMots"]), main="Empirical Discrete Distribution for the Number of words",xlab="Number of words")


###########
#Classification ascendante hiérarchique, cf. Stats avec R, p. 220 et sqq.
library(cluster)
CAHbrut1 = agnes(imported,method="ward") # pour Euclid. Pour Manhattan : metric="manhattan",method="ward")
plot(CAHbrut1,xlab="Individus")


#CAH : graphe des hauteurs
CAHbrut1bis = as.hclust(CAHbrut1)
plot(rev(CAHbrut1bis$height), type="h",ylab="hauteurs")
#CAH sur données réduites
CAHreduit = agnes(scale(dff),method="ward")
#Découpages en classes
CAHbrut1classes = cutree(CAHbrut1,k=6)#6 classes
#Insérer le nom des classes au tableau initial
ImportRenverseClasses = cbind.data.frame(ImportRenverse,as.factor(CAHbrut1classes))
#Avoir les éléments caractérisants des classes
library(FactoMineR)
catdes(dff)
###########


#Régression linéaire

RegressionRichesse = lm(Nombre.de.formes.differentes~Total.de.mots, data=Richesse);
summary(RegressionRichesse);
abline(RegressionRichesse)

##### Regression non linéaire ####
## Calcul des paramètres non linéaires Vm et K
Vm = 1/coef(RegressionRichesse)[2]
K = Vm*coef(RegressionRichesse)[1]
## Nonlinear regression model
RegNonLinRichesse = nls(NombreDeFormesDifferentes~Vm*TotalDeMots/(K+TotalDeMots),data=Richesse, start=list(K=55.07962, Vm=2.391919))

x = seq(0, 700, length=100)
y2 <- (coef(RegNonLinRichesse)["Vm"]*x)/(coef(RegNonLinRichesse)["K"]+x)

lines(x, y2)



###Kohonen



## With the Yasomi package

library("yasomi")

#
sg = somgrid(xdim=10,ydim=10,topo="hexagonal")

somtuning <- som.tune(NouvFreqCamps,TEST)

     som.tune(data, somgrid, control = som.tunecontrol(somgrid),
              weights, verbose = FALSE, internalVerbose = FALSE)

som <- somtuning$best.som


plot(x, y, mode=c("prototype","data"),
     type=c("parallel","stars","barplot"), 
     with.grid=FALSE, barplot.align=c("bottom","center"), global.scale=FALSE, ...)

plot(SOMTEST,mode="data",type="stars",with.grid=TRUE)

## Obtenir les distances entre les seuls prototypes correspondant à des observations sur les grandes cartes

est.non.vide = table(factor(som$classif, levels=1:nrow(som$prototypes)))!=0
dp = as.dist(as.matrix(prototype.distances(som))[est.non.vide,est.non.vide])

est.non.vide = table(factor(SOMBestNouvFreqCamps$classif, levels=1:nrow(SOMBestNouvFreqCamps$prototypes)))
SOMBestNouvFreqCamps.distancesnonvide = as.dist(as.matrix(prototype.distances(SOMBestNouvFreqCamps))[est.non.vide,est.non.vide])


#Puis calcul d'une CAH

#som.tune

Parameter tuning for Self-Organising Maps

Description:

     This function tunes some parameters of a Self-Organising Map by
     optimising a specified error measure. The prior structure is not
     optimised by this function.

Usage:

     som.tune(data, somgrid, control = som.tunecontrol(somgrid),
              weights, verbose = FALSE, internalVerbose = FALSE)
     ## S3 method for class 'somtune'
     print(x,...)
     
Arguments:

    data: the data to which the SOM will be fitted. This can be, e.g.,
          a matrix or data frame of observations (which should be
          scaled), or a distance matrix

Value:

     An object of class ‘"somtune"’, a list with components

best.som: the best SOM according to the chosen error criterion

  errors: the error for each configuration

quantisation: the quantisation error for each configuration

 isquant: ‘TRUE’ if the best SOM was chosen according to the
          quantisation error

 control: the control object used by this call

dimensions: a list of strings with the names of the parameters that
          were varied by the function

init,assignement,radii,annealing,kernel: 5 vectors containing the
          parameters used by each tested configuration

best.index: the index of the configuration used by the best SOM




##som.tune

Paramètres sur lesquels jouer:
assignment = single / heskes
kernel = gaussian / linear

test = som.tune(NouvFreqCamps,sg,som.tunecontrol(sg,init="random"),verbose=TRUE)

Description:

     Creates a list of parameters for the ‘som.tune’ function.

Usage:

     som.tunecontrol(somgrid, init = "pca", ninit = 1, assignment = "single",
                     radii = c(2, 2/3 * somgrid$diam), nradii = 10,
                     innernradii = 30, maxiter = 75, annealing = "power",
                     kernel = "gaussian", criterion = error.quantisation)
     
Arguments:

 somgrid: an object of class ‘"somgrid"’

    init: prototypes initialization method. Valid values are ‘"pca"’
          and ‘"random"’. The former corresponds to principal component
          based initialization (see ‘sominit.pca’), while the latter
          uses randomly selected observations as initial values for the
 prototypes (see ‘sominit.random’)

   ninit: number of initial prototype values to test (only relevant for
          ‘init="random"’)

assignment: assignment method with valid values ‘"single"’ and
          ‘"heskes"’ (see ‘batchsom’)

   radii: the range of radii to explore, i.e., a vector of length two
          containing a minimal and a maximal value of radii. The
          default minimum radius is 2 (almost purely local k-means like
          optimization) while the maximum is equal to two third of the
          diameter of the prior struture

  nradii: number of radii to generate from the range specified in
          ‘radii’

innernradii: number of radii to use in the annealing scheme during the
          SOM fitting (see ‘batchsom’)

 maxiter: maximal number of iteration for each radius during fitting

(see ‘batchsom’)

annealing: annealing scheme with valid values ‘"power"’ (exponential
          like annealing) and ‘"linear"’ (linear scheme)

  kernel: kernel chosen between ‘"gaussian"’ and ‘"linear"’

criterion: an error criterion, i.e., a function that evaluate the
          quality of a fitted som on a dataset


batchsom                package:yasomi                 R Documentation

Generic Self-Organising Map fitting function

Description:

     Generic function for fitting a Self-Organising Map to some data,
     using a batch algorithm.

Usage:

     batchsom(data, somgrid, init=c("pca","random"), prototypes, weights,
              mode = c("continuous","stepwise"), min.radius, max.radius,
              steps, decrease = c("power", "linear"), max.iter,
              kernel = c("gaussian", "linear"), normalised,
              assignment = c("single", "heskes"), cut = 1e-07,
              verbose = FALSE, keepdata = TRUE, ...)
     
Arguments:

    data: the data to which the SOM will be fitted. Acceptable data
          type depend on the available methods, see details

somgrid: an object of class ‘"somgrid"’ that specifies the prior
          structure of the Self-Organising Map: see ‘somgrid’

    init: the initialisation method (defaults to ‘"pca"’, see details)

prototypes: Initial values for the prototypes (the exact representation
          of the prototypes depends on the data type). If missing,
          initial prototypes are chosen via the method specified by the
          ‘init’ parameter (see details)

 weights: optional weights for the data points

    mode: annealing mode:

          ‘"continuous"’ (default) this is the standard annealing
              strategy for SOM: the influence of neighbours changes at
              each epoch of the algorithm, from ‘max.radius’ to
              ‘min.radius’ in exactly ‘step’ steps.

          ‘"stepwise"’ in this strategy, the algorithm performs several
              epochs (a maximum of ‘max.iter’) for each of the ‘step’
radii (from ‘max.radius’ to ‘min.radius’). The algorithm
              changes the neighbours influence only when the
              classification remains stable from one epoch to another.
              The ‘max.iter’ parameter provides a safeguard against
              cycling behaviours.

min.radius: the minimum neighbourhood influence radius. If missing, the
          value depends on the one of ‘kernel’ but ensures in practice
          a local learning only (see details)

max.radius: the maximal neighbourhood influence radius. If missing two
          third of the prior structure diameter plus one

   steps: the number of radii to use during annealing

decrease: the radii generating formula (‘"power"’ or ‘"linear"’), i.e.,
          the way the ‘steps’ radii are generated from the extremal
          values given by ‘min.radius’ and ‘max.radius’
max.iter: maximal number of epochs for one radius in the ‘"stepwise"’
          annealing mode (defaults to 75)

  kernel: the kernel used to transform distances in the prior structure
          into influence coefficients

normalised: switch for normalising the neighbouring interactions. Has
          no influence with the ‘"single"’ assignment method

assignment: the assignment method used to compute the best matching
          unit (BMU) of an observation during training:

          ‘"single"’ (default) this is the standard BMU calculation
              approach in which the best unit for an observation is the
              one of the closest prototype of this observation

          ‘"heskes"’ Tom Heskes' variant for the BMU in which a
              weighted fit of all the prototypes to an observation is
              used to compute the best unit. The rationale is that the
              BMU's prototype and its neighbouring units' prototypes
              must be close to the observation.

     cut: minimal value below wich neighbouring interactions are not  take into account

     verbose: switch for tracing the fitting process

keepdata: if ‘TRUE’, the original data are returned as part of the
          result object

     ...: additional arguments to be passed to methods




###MDSMultidimensional Scaling
#R provides functions for both classical and nonmetric multidimensional scaling. Assume that we have N objects measured on p numeric variables. We want to represent the distances among the objects in a parsimonious (and visual) way (i.e., a lower k-dimensional space).
#Classical MDS

#You can perform a classical MDS using the cmdscale( ) function.

# Classical MDS
# N rows (objects) x p columns (variables)
# each row identified by a unique row name

d <- dist(mydata) # euclidean distances between the rows
fit <- cmdscale(d,eig=TRUE, k=2) # k is the number of dim
fit # view results

# plot solution
x <- fit$points[,1]
y <- fit$points[,2]
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2",
  main="Metric MDS", type="n")
text(x, y, labels = row.names(NouvFreqCamps), cex=.7) 

##Nonmetric MDS

##Nonmetric MDS is performed using the isoMDS( ) function in the MASS package.

# Nonmetric MDS
# N rows (objects) x p columns (variables)
# each row identified by a unique row name

library(MASS)
d <- dist(mydata) # euclidean distances between the rows
fit <- isoMDS(d, k=2) # k is the number of dim
fit # view results

# plot solution
x <- fit$points[,1]
y <- fit$points[,2]
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2",
  main="Nonmetric MDS", type="n")
text(x, y, labels = row.names(ComplFreq100Camps), cex=.7) 



#####Fuzzy C-Means Clustering
##Description

The fuzzy version of the known kmeans clustering algorithm as well as its online update (Unsupervised Fuzzy Competitive learning).
Usage

cmeans (dataframe, centers, iter.max=100, verbose=FALSE, dist="euclidean",method="cmeans", m=2, rate.par = NULL)
        
        
        
Arguments
x 	The data matrix where columns correspond to variables and rows to observations
centers 	Number of clusters or initial values for cluster centers
iter.max 	Maximum number of iterations
verbose 	If TRUE, make some output during learning
dist 	Must be one of the following: If "euclidean", the mean square error, if "manhattan", the mean absolute error is computed. Abbreviations are also accepted.
method 	If "cmeans", then we have the cmeans fuzzy clustering method, if "ufcl" we have the On-line Update. Abbreviations in the method names are also accepted.
m 	The degree of fuzzification. It is defined for values greater than 1
rate.par 	The parameter of the learning rate


##Graphe du même


plot(dataframe, main="C-means: 2-way Fuzzy Membership", type="n", xlab="Variable 1", ylab="Variable 2")
points(cc$centers, col = c("red", "blue"), pch = 8, cex=2)

mfuzz.plot(dataframe, cl = POUAITE, mfrow = c(4, 4), time.labels = seq(0, 160, 10))

fanny(dataframe, 9, diss = FALSE, memb.exp = 2,metric = "euclidean",stand = FALSE, iniMem.p = NULL, cluster.only = FALSE,keep.diss = !diss && !cluster.only && n < 100, keep.data = !diss && !cluster.only, maxit = 500, tol = 1e-15, trace.lev = 0)




####
# Distribution
#on imprimer le graphique, sans l'échelle pour y, car on ne veut que des entiers, et pas l'échelle par défaut qui intègre les O.5
plot(table(dataframe[,"NombreDeFormesDifferentes"]), ylab="Number of individuals",xlab="number of word-forms", main="Barplot for the number of word-forms", yaxt="n")
#ajout de l'échelle pour y
axis(2, at =  1:3)

## Fonction de la distribution empirique discrète
plot(ecdf(dataframe[,"NombreDeFormesUniques"]), main="Empirical Discrete Distribution for the Number of hapax",xlab="Number of hapax")

## « Équation de centrage » (centrage de la matrice)

scale(X, center=TRUE, scale=TRUE) # centre et réduit les colonnes d'une matrice

### Spécifier la taille des labels sur un plot
cex=0.75

### Transformer data.frame en valeurs numériques
data.matrix(frame, rownames.force = NA)


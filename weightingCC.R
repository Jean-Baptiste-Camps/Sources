#Transforme un tableau d'occurrences dont les colonnes doivent être les textes et les lignes les mots. Les transforme selon une méthode maison pour éviter que les classifications ne soient influencées par la longueur des textes. Ainsi, le nombre brut d'occurrences est transformé en un taux, selon la formule Nombre d'occurrences dans le texte / nombre de mots dans le texte.

"weightingCC" = function(x) {
	X = x ;
	M = ncol(x);
	for(i in 1:M){
		X[,i] = X[,i]/sum(X[,i])
	}
	return(X)
}



limma.test.paired=function(tab,Conditions){
  require(limma)
  Subject=rep(x = 1:(length(Conditions)/2),2)
  design <- model.matrix(~Subject+Conditions)
  print(design)
  fit <- lmFit(tab, design)
  print("fit")
  print(fit)
  resfit <- eBayes(fit,robust=TRUE)
  print("resfit")
  print(resfit)
  return(resfit)
}



##############################################################
#Fonction réalisant un test de contrastes entre conditions à l'aide du package LIMMA.
#
#En entrée:
#tab: tableau de données sans valeurs manquantes avec chaque réplicat en colonne
#Conditions: indique le numéro de la condition biologique auquel appartient chaque réplicat
#RepBio: indique le numéro du réplicat biologique auquel appartient chaque réplicat
#RepTech: indique le numéro du réplicat technique auquel appartient chaque réplicat
#
#Contrast: indique si l'on souhaite tester chaque condition biologique contre chacune (Contrast=1; par exemple H0:"C1=C2" vs H1:"C1!=C2", etc.) 
#ou chaque condition contre toutes les autres (Contrast=2; par exemple H0:"C1=(C2+C3)/2" vs H1:"C1!=(C2+C3)/2", etc. si on a 3 conditions ).
#
#En sortie :
#Objet fit renvoyé par la fonction eBayes de LIMMA.
##############################################################

limma.test=function(tab,Conditions,RepBio,RepTech,Contrast=1){
  require(limma)
  #####################################################################
  #Fonctions permettant de créer des design et contrastes
  #######################
  #création de la matrice du design hiérarchique : cas de 2 niveaux Bio-Tech ou Bio-Analytique ou Tech-Analytique .
  make.design.2=function(tab,Condition,RepBio){
    #Renomme les niveaux des facteurs
    levels(Condition)=c(1:length(levels(Condition)))
    levels(RepBio)=c(1:length(levels(RepBio)))
    #Matrice de design initiale
    design=model.matrix(tab[1,]~0+Condition:RepBio)
    #Suppression des colonnes nulles de la matrice de design
    design=design[,(apply(design,2,sum)>0)]
    #Suppression des colonnes égales de la matrice de design
    coldel=-1
    for (i in 1:(length(design[1,])-1)){
      d2=as.matrix(design[,(i+1):length(design[1,])]);
      for (j in 1:length(d2[1,])){
        d2[,j]=d2[,j]-design[,i];
      }
      e=as.matrix(rnorm(length(design[,1]),10,1));
      sd2=t(e)%*%d2
      liste=which(sd2==0)
      coldel=c(coldel,liste+i)
    }
    design=design[,(1:length(design[1,]))!=coldel]
    colnames(design)=make.names(colnames(design))
    return(design)
  }
  #######################
  #création de la matrice du design hiérarchique : cas de 3 niveaux Bio, Tech et Analytique.
  make.design.3=function(tab,Condition,RepBio,RepTech){
    #Renomme les niveaux des facteurs
    levels(Condition)=c(1:length(levels(Condition)))
    levels(RepBio)=c(1:length(levels(RepBio)))
    levels(RepTech)=c(1:length(levels(RepTech)))
    #Matrice de design initiale
    design=model.matrix(tab[1,]~0+Condition:RepBio:RepTech)
    #Suppression des colonnes nulles de la matrice de design
    design=design[,(apply(design,2,sum)>0)]
    #Suppression des colonnes égales de la matrice de design
    coldel=-1
    for (i in 1:(length(design[1,])-1)){
      d2=as.matrix(design[,(i+1):length(design[1,])]);
      for (j in 1:length(d2[1,])){
        d2[,j]=d2[,j]-design[,i];
      }
      e=as.matrix(rnorm(length(design[,1]),10,1));
      sd2=t(e)%*%d2
      liste=which(sd2==0)
      coldel=c(coldel,liste+i)
    }
    design=design[,(1:length(design[1,]))!=coldel]
    colnames(design)=make.names(colnames(design))
    return(design)
  }
  #######################
  #aggrégation des colonnes correspondant à une même condition pour créer la matrice de contraste.
  aggreg.column.design=function(design,Condition){
    nb.cond=length(levels(Condition))
    name.col=colnames(design)
    name.cond=NULL
    nb.col=NULL
    for (i in 1:nb.cond){
      col.select=NULL
      col.name.begin=paste("Condition",i, sep = "")
      nc=nchar(col.name.begin)
      for (j in 1:length(design[1,])){
        if (substr(name.col[j], 1, nc)==col.name.begin){col.select=c(col.select,j)}
      }
      name.aggreg=NULL
      for (j in 1:length(col.select)){
        name.aggreg=paste(name.aggreg,name.col[col.select[j]],sep="+")
      }
      name.aggreg=substr(name.aggreg, 2, nchar(name.aggreg))
      name.cond=c(name.cond,name.aggreg)
      nb.col=c(nb.col,length(col.select))
    }
    print(name.cond,nb.col)
    return(list(name.cond,nb.col))
  }
  #######################
  #Matrice de contraste : cas une condition contre une condition.
  make.contraste.1.1=function(design,Condition){
    nb.cond=length(levels(Condition))
    r=aggreg.column.design(design,Condition)
    label.agg=r[[1]]
    nb.agg=r[[2]]
    k=1
    contra=rep(0,sum(1:(nb.cond-1)))
    for (i in 1:(nb.cond-1)){
      for (j in (i+1):nb.cond){
        contra[k]=c(paste("(",label.agg[i],")/",nb.agg[i],"-(",label.agg[j],")/",nb.agg[j]))
        k=k+1
      }
    }
    return(contra)
  }
  #######################
  #Matrice de contraste : cas une condition contre toutes les autres conditions.
  make.contraste.2.1=function(design,Condition){
    nb.cond=length(levels(Condition))
    r=aggreg.column.design(design,Condition)
    label.agg=r[[1]]
    nb.agg=r[[2]]
    k=1
    contra=rep(0,sum(1:(nb.cond-1)))
    for (i in 1:(nb.cond)){
      contra[k]=c(paste("(",label.agg[i],")/",nb.agg[i]))
      nb=sum(nb.agg[(1:nb.cond)[(1:nb.cond)!=i]])
      for (j in (1:nb.cond)[(1:nb.cond)!=i]){
        contra[k]=c(paste(contra[k],"-(",label.agg[j],")/",nb))
      }
      k=k+1
    }
    return(contra)
  }
  #####################################################################
  #Fonction générale
  lt=length(tab[1,])
  lc=length(Conditions)
  lb=length(RepBio)
  lte=length(RepTech)
  #vérification de la bonne longueur des vecteurs
  if (length(levels(Conditions))!=lc){
    if(lt==lc){
      if (lb==lc){
        if (lte==lc){  
          #vérification que le nbre de facteurs dans Conditions est inférieur à celui de RepBio qui est lui-même inférieur à celui de RepTech
          if (length(levels(Conditions))<=length(levels(RepBio))){
            if (length(levels(RepBio))<=length(levels(RepTech))){
              #détermination du nombre de niveaux hiérarchiques
              if (length(levels(RepBio))==lt){niveau=1;}
              else{ if (length(levels(RepTech))==lt){niveau=2;}
                else{niveau=3;}
              }
              #Cas non hiérarchique : on a que des réplicats bio ou tech ou anal.
              if (niveau==1){
                nb_cond=length(levels(Conditions))
                #Calcul du nombre de réplicats par conditions
                nb_Rep=rep(0,nb_cond)
                for (i in 1:nb_cond){
                  nb_Rep[i]=sum((Conditions==levels(Conditions)[i]))
                }
                #####################
                #Calcul de la matrice de design
                design=matrix(0,lt,nb_cond)
                n0=1
                coln=NULL
                for (j in 1:nb_cond){
                  coln=c(coln,paste("Condition",j,collapse=NULL,sep=""))
                  design[(n0:(n0+nb_Rep[j]-1)),j]=rep(1,length((n0:(n0+nb_Rep[j]-1))))
                  n0=n0+nb_Rep[j]
                }
                colnames(design)=coln
                #####################
                #Création des contrastes
                if (Contrast==1){contra=make.contraste.1.1(design,Condition=Conditions);}
                if (Contrast==2){contra=make.contraste.2.1(design,Condition=Conditions);}
                cmtx=makeContrasts(contrasts=contra,levels=make.names(colnames(design)))
                fit=eBayes(contrasts.fit(lmFit(tab, design), cmtx),robust=TRUE)
                return(fit)  
              }
              if (niveau==2){
                #####################
                #Calcul de la matrice de design
                design=make.design.2(tab,Condition=Conditions,RepBio=RepBio)
                #####################
                #Calcul de la matrice de contrastes
                if (Contrast==1){contra=make.contraste.1.1(design,Condition=Conditions);}
                if (Contrast==2){contra=make.contraste.2.1(design,Condition=Conditions);}
                cmtx=makeContrasts(contrasts=contra,levels=make.names(colnames(design)))
                fit <- eBayes(contrasts.fit(lmFit(tab, design), cmtx),robust=TRUE)
                return(fit)
              }
              if (niveau==3){
                #####################
                #Calcul de la matrice de design
                design=make.design.3(tab,Condition=Conditions,RepBio=RepBio,RepTech=RepTech)
                #####################
                #Calcul de la matrice de contrastes
                if (Contrast==1){contra=make.contraste.1.1(design,Condition=Conditions);}
                if (Contrast==2){contra=make.contraste.2.1(design,Condition=Conditions);}
                cmtx=makeContrasts(contrasts=contra,levels=make.names(colnames(design)))
                fit <- eBayes(contrasts.fit(lmFit(tab, design), cmtx),robust=TRUE)
                return(fit)
              }
            }
            else{cat("Problème : le nombre de facteurs dans Réplicats Bio. doit être inférieur ou égale au nombre de facteurs dans Réplicats Tech.!\n")}
            
          }
          else{cat("Problème : le nombre de facteurs dans Conditions doit être inférieur ou égale au nombre de facteurs dans Réplicats Bio.!\n")}
          
        }
        else{cat("Problème : le vecteur-facteur Réplicats Tech. doit être de longueur égale au vecteur-facteur Réplicats Bio.!\n")}
        
      }
      else{cat("Problème : le vecteur-facteur Réplicats Bio. doit être de longueur égale au vecteur-facteur Conditions!\n")}
      
    }
    else{cat("Problème : le vecteur-facteur Conditions doit être de même longueur que le nombre de colonnes de la matrice d'entrée !\n")}
    
  }
  else{cat("Problème : le vecteur-facteur Conditions doit contenir plusieurs réplicats dans chaque condition !\n")}
}


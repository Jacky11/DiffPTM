
# Venn Diagram


makeVenn=function(datpeptC1,nintersects=NA){
  require(VennDiagram)
  require(venn)
  require(UpSetR)
  require(ggplotify)
  require(grid)
  #require(ggimage)
  
  datpeptC1[datpeptC1>0]=1;
  
  ldpc=ncol(datpeptC1);
  print("oo");
  print(ldpc)
  #on limite a 15 le nombre max d'echantillons par condition
  venn.plot=NULL;
  if (ldpc<70){
    #calcul des similitudes entre replicats
    ones=rep(1,nrow(datpeptC1));
    vecto=as.matrix(1:nrow(datpeptC1));
    print("xxx")
    peptid=apply(datpeptC1,2,function(x) x*vecto);
    print("oooooo")
    #calcul des intersections entre deux
    if (ldpc<=5){
    nn=matrix(0,ldpc,ldpc);
    for (i in 1:(ldpc)){
      for (j in (i):ldpc){
        inter=intersect(peptid[,i],peptid[,j]);
        nn[i,j]=length(inter[inter!=0]);
      }
    }
    
    if (ldpc>2){
      #calcul des intersections entre trois
      nnn=array(0,dim=c(ldpc,ldpc,ldpc))
      for (i in 1:(ldpc)){
        for (j in (i):ldpc){
          for (k in (j):ldpc){
            inter=intersect(peptid[,i],intersect(peptid[,j],peptid[,k]));
            nnn[i,j,k]=length(inter[inter!=0]);
          }
        }
      }
    }
    
    if (ldpc>3){
      #calcul des intersections entre quatres
      nnnn=array(0,dim=c(ldpc,ldpc,ldpc,ldpc))
      for (i in 1:(ldpc)){
        for (j in (i):ldpc){
          for (k in (j):ldpc){
            for (l in (k):ldpc){
              inter=intersect(peptid[,i],intersect(peptid[,j],intersect(peptid[,k],peptid[,l])));
              nnnn[i,j,k,l]=length(inter[inter!=0]);
            }
          }
        }
      }
    }
    if (ldpc>4){
      #calcul des intersections entre cinqs
      nnnnn=array(0,dim=c(ldpc,ldpc,ldpc,ldpc,ldpc))
      for (i in 1:(ldpc)){
        for (j in (i):ldpc){
          for (k in (j):ldpc){
            for (l in (k):ldpc){
              for (m in (l):ldpc){
                inter=intersect(peptid[,i],intersect(peptid[,j],intersect(peptid[,k],intersect(peptid[,l],peptid[,m]))));
                nnnnn[i,j,k,l,m]=length(inter[inter!=0]);
              }
            }
          }
        }
      }
    }
    
    if (ldpc==2){
      venn.plot <- draw.pairwise.venn(
        area1 = nn[1,1],
        area2 = nn[2,2],
        cross.area = nn[1,2],
        category = colnames(datpeptC1),
        fill = c("dodgerblue", "goldenrod1"),
        cat.col = c("dodgerblue", "goldenrod4"),
        cat.cex = 1.5,
        cat.dist=0.1,
        margin = 0.1,
        cex = c(1.5, 2, 1.5),
        ind = TRUE
      );
      grid.newpage();
      p1 <- grobTree(venn.plot);
      print(class(p1))
      venn.plot <- as.ggplot(p1);
    }  
    if (ldpc==3){
      venn.plot <- draw.triple.venn(
        area1 = nn[1,1],
        area2 = nn[2,2],
        area3 = nn[3,3],
        n12 = nn[1,2],
        n13 = nn[1,3],
        n23 = nn[2,3],
        n123 = nnn[1,2,3],
        category = colnames(datpeptC1),
        fill = c("dodgerblue", "goldenrod1", "seagreen3"),
        cat.col = c("dodgerblue", "goldenrod4", "seagreen4"),
        cat.cex = 1.5,
        cat.dist=0.1,
        margin = 0.1,
        cex = c(1.5, 1.5, 1.5, 1.5, 2, 1.5, 1.5),
        ind = TRUE
      );
      grid.newpage();
      p1 <- grobTree(venn.plot);
      print(class(p1))
      venn.plot <- as.ggplot(p1);
    }    
    if (ldpc==4){
      venn.plot <- draw.quad.venn(
        area1 = nn[1,1],
        area2 = nn[2,2],
        area3 = nn[3,3],
        area4 = nn[4,4],
        n12 = nn[1,2],
        n13 = nn[1,3],
        n14 = nn[1,4],
        n23 = nn[2,3],
        n24 = nn[2,4],
        n34 = nn[3,4],
        n123 = nnn[1,2,3],
        n124 = nnn[1,2,4],
        n134 = nnn[1,3,4],
        n234 = nnn[2,3,4],
        n1234 = nnnn[1,2,3,4],
        category = colnames(datpeptC1),
        fill = c("dodgerblue", "plum1", "darkorange1", "seagreen3"),
        cat.col = c("dodgerblue", "plum4", "darkorange2", "seagreen4"),
        cat.cex = 1.5,
        cat.dist=0.22,
        margin = 0.1,
        cex = c(1.5, 1.5, 1.5, 1.5, 1.5, 2, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5),
        ind = TRUE
      );
      grid.newpage();
      p1 <- grobTree(venn.plot);
      venn.plot <- as.ggplot(p1);
    }
    if (ldpc==5){
      # Reference five-set diagram
      venn.plot <- draw.quintuple.venn(
        area1 = nn[1,1],
        area2 = nn[2,2],
        area3 = nn[3,3],
        area4 = nn[4,4],
        area5 = nn[5,5],
        n12 = nn[1,2],
        n13 = nn[1,3],
        n14 = nn[1,4],
        n15 = nn[1,5],
        n23 = nn[2,3],
        n24 = nn[2,4],
        n25 = nn[2,5],
        n34 = nn[3,4],
        n35 = nn[3,5],
        n45 = nn[4,5],
        n123 = nnn[1,2,3],
        n124 = nnn[1,2,4],
        n125 = nnn[1,2,5],
        n134 = nnn[1,3,4],
        n135 = nnn[1,3,5],
        n145 = nnn[1,4,5],
        n234 = nnn[2,3,4],
        n235 = nnn[2,3,5],
        n245 = nnn[2,4,5],
        n345 = nnn[3,4,5],
        n1234 = nnnn[1,2,3,4],
        n1235 = nnnn[1,2,3,5],
        n1245 = nnnn[1,2,4,5],
        n1345 = nnnn[1,3,4,5],
        n2345 = nnnn[2,3,4,5],
        n12345 = nnnnn[1,2,3,4,5],
        category = colnames(datpeptC1),
        fill = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
        cat.col = c("dodgerblue", "goldenrod4", "darkorange2", "seagreen4", "orchid4"),
        cat.cex = 1.5,
        cat.dist=0.26,
        margin = 0.05,
        cex = c(1.5, 1.5, 1.5, 1.5, 1.5, 1, 0.8, 1, 0.8, 1, 0.8, 1, 0.8, 1, 0.8,
                1, 0.75, 1, 0.75, 1, 0.75, 1, 0.75, 1, 0.75, 1, 1, 1, 1, 1, 2),
        ind = TRUE
      );
      grid.newpage();
      p1 <- grobTree(venn.plot);
      venn.plot <- as.ggplot(p1);
    }
    }
    # if (ldpc==6){
    #   if (length(which(is.na(datpeptC1[,1])))!=0){
    #     datpeptC1=datpeptC1[-which(is.na(datpeptC1[,1])),];
    #   }
    #   snb=apply(datpeptC1,1,sum)
    #   if (length(which(snb==0))!=0){
    #     datpeptC1=datpeptC1[-which(snb==0),];
    #   }
    #   datpeptC1=as.data.frame(datpeptC1);
    #   venn.plot<-venn::venn(x=datpeptC1,zcolor="style",cexil=0.75);
    # }
    # if (ldpc==7){
    #   if (length(which(is.na(datpeptC1[,1])))!=0){
    #     datpeptC1=datpeptC1[-which(is.na(datpeptC1[,1])),];
    #   }
    #   snb=apply(as.matrix(datpeptC1),1,sum)
    #   if (length(which(snb==0))!=0){
    #     datpeptC1=datpeptC1[-which(snb==0),];
    #   }
    #   datpeptC1=as.data.frame(datpeptC1);
    #   venn.plot<-venn::venn(x=datpeptC1,zcolor="style",cexil=0.75);
    # }
    if (ldpc>5){
      print("prout");
      lp2=as.data.frame(datpeptC1);
      print("oo2");
      names(lp2)=colnames(datpeptC1);
      print(head(lp2));
      
      
      venn.plot <- upset(lp2, 
                         order.by = "freq",number.angles=315,
                         nsets=ldpc,
                         nintersects=25,
                         mb.ratio=c(0.5,0.5),
                         main.bar.color=4,
                         sets.x.label="Nb of proteins",
                         text.scale = c(1, 1, 0.5, 1, 1, 1));
      venn.plot <- as.ggplot(venn.plot);
      venn.plot =  venn.plot + coord_flip();
    }
    #    queries = list(list(query = intersects, params = as.list(names(lp2)), color= "red")));
  }else{
    df=data.frame(x=c(0),y=c(0))
    rownames(df)=c("Too much samples")
    venn.plot <- ggplot(df, aes(x, y, label = rownames(df))) + geom_point()
    venn.plot <- venn.plot + geom_text()
  }
  return(venn.plot)
}

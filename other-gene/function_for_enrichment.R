
require(DBI)
require(AnnotationDbi)
require(Biobase)
require(org.Hs.eg.db)
require(Category)
require(graph)
require(GOstats)
require(GO.db)
require(KEGG.db)

#ontology="BP","CC","MF"
enrichmentTest=function(genes,ontology,outFile)
{
  #genes <- scan(i,quote="\"",sep="\n",what="character");
  goAnn <- get("org.Hs.egGO")
  universe <- Lkeys(goAnn)
  entrezIDs <- mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
  entrezIDs <- as.character(entrezIDs)
  params <- new("GOHyperGParams",
                geneIds=entrezIDs,
                universeGeneIds=universe,
                annotation="org.Hs.eg.db",
                ontology="BP",
                pvalueCutoff=0.05,
                conditional=FALSE,
                testDirection="over")
  over <- hyperGTest(params)
  
  glist <- geneIdsByCategory(over)
  glist <- sapply(glist, function(.ids) {
    .sym <- mget(.ids, envir=org.Hs.egSYMBOL, ifnotfound=NA)
    .sym[is.na(.sym)] <- .ids[is.na(.sym)]
    paste(.sym, collapse=";")
  })
  bp <- summary(over)
  bp$Symbols <- glist[as.character(bp$GOBPID)]
  result<-summary(over,categorySize=0)
  write.table(bp, file=outFile,sep="\t",quote = FALSE,col.names = NA)
  
  bp
  
}





enrichmentTest_yeast=function(genes,ontology,outFile)
{
  
  require(GOstats)
  require(org.Sc.sgd.db)
  require(annotate)
  
  orf2go <- getAnnMap('GO', 'org.Sc.sgd')
  go2orf <- getAnnMap('GO2ORF', 'org.Sc.sgd')
  all.orfs <- keys(orf2go)
  

  selected=c()
  for(i in c(1:length(genes)))
  {
    tmp=orf2go[[genes[i]]]
    #print(tmp)
    if(!is.null(tmp))
    {
      if(!is.na(tmp))
      {
        selected=c(selected,i)
      }
      
    }
  }
  genes=genes[selected]
  
  print(length(genes))
  
  if(length(genes)>0)
  {
    ## Do GOstats test to get significant GO ontologies
    params <- new("GOHyperGParams", geneIds=genes,
                  universeGeneIds=all.orfs, annotation='org.Sc.sgd',
                  ontology=ontology, pvalueCutoff=0.05,
                  conditional=FALSE, testDirection='over')
    
    
    
    
    over <- hyperGTest(params)
    
    glist <- geneIdsByCategory(over)
    bp <- summary(over)
    bp$Symbols <- glist[as.character(bp$GOBPID)]
    result<-summary(over,categorySize=0)
    write.table(as.matrix(bp), file=outFile,sep="\t",quote = FALSE,col.names = NA)
    
  }else{
    bp=NA
  }

  
  bp
}


doEnrichment_yeast=function(geneName,clusterColor,output_dir,prefix)
{
  GO_terms_mat=NA
  geneNameList=as.character(geneName)
  #for each class, do Go enrichment
  classIndex=unique(clusterColor)
  for(i in classIndex)
  {
    A=which(clusterColor==i,arr.ind = T)
    genesInClass=geneNameList[A]
    #convert to gene symbol
    #geneSymbolsInClass=convertID(genesInClass,fromID,toID)
    print(i)
    bp=NA
    bp=enrichmentTest_yeast(genesInClass,ontology = "BP",outFile=paste(output_dir,prefix,i,"BP.txt",sep="_"))
    if(!is.null(dim(bp)))
    {
      if(dim(bp)[1]>0)
      {
        pvalues=bp[,2]
        GO_terms=bp[,1]
        mat1=bp[,1:2]
        if(is.na(GO_terms_mat))
        {
          GO_terms_mat=mat1
        }else{
          GO_terms_mat=rbind(GO_terms_mat,mat1)
        }
      }
      
    }
    
  }
  GO_terms_mat
}

doEnrichment=function(geneName,clusterColor,output_dir,prefix)
{
  GO_terms_mat=NA
  geneNameList=as.character(geneName)
  #for each class, do Go enrichment
  classIndex=unique(clusterColor)
  for(i in classIndex)
  {
    A=which(clusterColor==i,arr.ind = T)
    genesInClass=geneNameList[A]
    #convert to gene symbol
    #geneSymbolsInClass=convertID(genesInClass,fromID,toID)
    print(i)
    bp=enrichmentTest(genesInClass,ontology = "BP",outFile=paste(output_dir,prefix,i,"BP.txt",sep="_"))
    pvalues=bp[,2]
    GO_terms=bp[,1]
    mat1=bp[,1:2]
    if(is.na(GO_terms_mat))
    {
      GO_terms_mat=mat1
    }else{
      GO_terms_mat=rbind(GO_terms_mat,mat1)
    }
  }
  GO_terms_mat
}

#library(Category)
KEGGTest=function(genes,pvalueCutoff=0.05,outFile)
{
  entrezIDs <- mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
  entrezIDs <- as.character(entrezIDs)
  
  keggAnn <- get("org.Hs.egPATH")
  universe <- Lkeys(keggAnn)
  params <- new("KEGGHyperGParams", 
                geneIds=entrezIDs, 
                universeGeneIds=universe, 
                annotation="org.Hs.eg.db", 
                categoryName="KEGG", 
                pvalueCutoff=pvalueCutoff,
                testDirection="over")
  over <- hyperGTest(params)
  kegg <- summary(over)
  
  glist <- geneIdsByCategory(over)
  glist <- sapply(glist, function(.ids) {
    .sym <- mget(.ids, envir=org.Hs.egSYMBOL, ifnotfound=NA)
    .sym[is.na(.sym)] <- .ids[is.na(.sym)]
    paste(.sym, collapse=";")
  })
  kegg$Symbols <- glist[as.character(kegg$KEGGID)]
  
  write.table(kegg, file=outFile,sep="\t",quote = FALSE, col.names = NA)
  kegg
}
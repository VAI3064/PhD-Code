function(workingdirectory,targetfile,designvector,groupsvector,study,patternsvector,colnamesvector,numberofgroups,filename,gene.name.file)
{
    setwd(workingdirectory)
    
    suppressMessages(library(limma))
    suppressMessages(library(Biobase))
    suppressMessages(library(marray))
    suppressMessages(library(stats))
    suppressMessages(library(convert))
    suppressMessages(library(methods))
    suppressMessages(library(lattice))
    suppressMessages(library(coda))
    suppressMessages(library(gaga))
    suppressMessages(library(gplots))
    suppressMessages(library(topGO))
    suppressMessages(library(affycoretools))

    
    targets<- readTargets(targetfile)
    gene_names_file<- read.table(paste(workingdirectory,gene.name.file,sep=""), header = TRUE, fill = TRUE, sep = "\t")
    gene_names_file<-as.matrix(gene_names_file)
    red_green_data<- read.maimages(targets, source = "agilent.median", annotation = c("Row", "Col", "FeatureNum", "ControlType", "ProbeName", "SystematicName")) 
    red<- as.matrix(red_green_data$R) 
    green<- as.matrix(red_green_data$G)
    
    #Extracts red and green median values, puts in a matrix with columns being each array.
    source(paste(workingdirectory,"/Functions/RGdyeswap.R",sep=""))
    design<- designvector
    new_RG<- RGdyeswap(red, green, design)
    new_RG$R<-new_RG[,1:length(groupsvector)]
    new_RG$G<-new_RG[,(length(groupsvector)+1):(length(groupsvector)+length(groupsvector))]
    red_green_data$R<- new_RG$R
    red_green_data$G<- new_RG$G
    positions<-grep("DETOX",red_green_data$genes$ProbeName,value=FALSE)

    source(paste(workingdirectory,"/Functions/correct_gene_names3.R",sep=""))
    new_new_red_green<- correct_gene_names3(red_green_data, gene_names_file,positions)
    new_red_green<- new_new_red_green[which(new_new_red_green$genes$ControlType == 0),]
    source(paste(workingdirectory,"/Functions/MyNormalizeWithinArrays.R",sep=""))
    source(paste(workingdirectory,"/Functions/myMA.RG.R",sep=""))
    source(paste(workingdirectory,"/Functions/myMA2.RG.R",sep=""))
    source(paste(workingdirectory,"/Functions/MyNormalizeBetweenArrays.R",sep=""))
    convert_to_MA<- mynormalizeWithinArrays(new_red_green, method = "loess")
    within_norm_RG<- RG.MA(convert_to_MA)
    within_norm_RG<-backgroundCorrect(within_norm_RG, method="normexp",normexp.method='mle', offset=50,verbose=TRUE)
    between_norm<- mynormalizeBetweenArrays(within_norm_RG, method = "Aquantile")
    between_norm2<- between_norm
    new_offset<- min(between_norm2$M)
    between_norm2$M <- between_norm2$M + abs(new_offset)
    between_norm2$A <- between_norm2$A + abs(new_offset)
    out_RG<-RG.MA(between_norm)
    
    dir.create(paste(workingdirectory,"/",filename,sep=""))
    
    for(i in 1:length(groupsvector))
        {
            output<-between_norm[,i]
            output2<-out_RG[,i]
            nameoffile<-output$targets$FileName
            nameoffile2<- paste(workingdirectory,"/",filename,"/",nameoffile,sep="")
            
            if(file.exists(nameoffile2)==FALSE)
            {
                output2<-cbind(output$genes$ProbeName,output$genes$SystematicName,output$M[,1],output$A[,1],output2$R[,1],output2$G[,1])
                colnames(output2)<-c("Probe","ID_REF","M","A","R","G")
                write.table(output2, file = paste(workingdirectory,"/",filename,"/",nameoffile,sep=""), row.names = FALSE)
            }
        }
    
    #Output a matrix of normalised values for each file. If file exists, skips.
	match_names<- colnames(between_norm2$M)
    row.names(between_norm2$targets)<- match_names
    expression<- as(between_norm2, "ExpressionSet")
    groups<-groupsvector
    expression$groups<- groups
    expression$groups <- as.factor(expression$group)
    
    if(length(study)>0)
    {
        expression$study<-study
        expression$study<-as.factor(expression$study)
    }
    
    par(mfrow=c(2,1))
    
    x11()
    
    x2<- exprs(expression)
    row.names(x2)<- new_red_green$genes$ProbeName
    texts<-1:length(groupsvector)
    texts<-as.character(texts)
    
    if(length(study)!=0)
    {
        plotPCA(x2,groups=study,main="PCA plot Before",groupnames = NULL, addtext = texts, x.coord = NULL,y.coord = NULL, screeplot = FALSE, squarepca = FALSE, pch = NULL, col = NULL,legend=FALSE)
        x2<-c()
        design<-model.matrix(~ -1 + as.factor(expression$study) + as.factor(expression$group))
        lm1<-lmFit(expression,design)
        p<-1:nlevels(expression$study)
        linpred<-coef(lm1)[,p] %*% t(design[,p])
        exprs(expression) <- exprs(expression) - linpred + rowMeans(exprs(expression),na.rm=FALSE)
        if (min(exprs(expression)) < 0)
        {
            exprs(expression) <- exprs(expression) + abs(min(exprs(expression)))
        }
    }
    
    x<- exprs(expression)
    row.names(x)<- new_red_green$genes$ProbeName
    
    #Inserts the Probe name as row names
    plotPCA(x,groups=study,main="PCA plot After",groupnames = NULL, addtext = texts, x.coord = NULL,y.coord = NULL, screeplot = FALSE, squarepca = FALSE, pch = NULL, col = NULL,legend=FALSE)
    par(mfrow=c(1,1))
    
    patterns<- matrix(patternsvector, ncol = length(colnamesvector), byrow = TRUE)
    colnames(patterns)<- colnamesvector
    gaga_fit<- fitGG(x = x, groups = expression$groups, patterns = patterns, nclust = 1, method = "EM", trace = FALSE)
    gaga_model<- parest(gaga_fit, x = x, groups = expression$groups, alpha = 0.05)
    print(gaga_model)
    x11()
    par( mfrow = c( 2, 2 ) )
    checkfit(gaga_model, x = x, groups = expression$groups, type = "data", main = "")
    checkfit(gaga_model, x = x, groups = expression$groups, type = "shape", main = "")
    checkfit(gaga_model, x = x, groups = expression$groups, type = "mean", main = "") 
    checkfit(gaga_model, x = x, groups = expression$groups, type = "shapemean", main = "", xlab = "Mean", ylab = "1/sqrt(CV)")
    
    namesvec<-rep(NA,length(colnamesvector))
    j<-rep(0,length(colnamesvector))
    j2<-rep(0,length(colnamesvector))
    
    for(i in 1:length(colnamesvector))
    {
        if(i==1)
        {
            j[i]<-length(which(groups==i))
            namesvec[i]<-paste("pred",i,sep="")
        }
        else
        {
            j[i]<-j[i-1]+length(which(groups==i))
            namesvec[i]<-paste("pred",i,sep="")
        }
    }  
    
    j2<-(j*-1)
    
    for(i in 1:length(j))
    {
        assign(namesvec[i],classpred(gaga_model, xnew = x[, j[i]], x = x[,j2], groups = expression$groups[j2], ngene = 5000))
        pred<-get(namesvec[i])
        print(pred)
    }
    
    d1<- findgenes(gaga_model, x, expression$groups, fdrmax = 0.05, parametric = TRUE, B = 1500)
    print(table(d1$d))
    
    all_names<-new_red_green$gene$SystematicName
    
    namesDE<-rep(NA,numberofgroups)
    namegenes<-rep(NA,numberofgroups)

    for(i in 1:length(namesDE))
    {
        namesDE[i]<-paste("DE_genes",(i-1),sep="")
        namegenes[i]<-paste("genes",(i-1),sep="")
    }
    
    for(i in 1:numberofgroups)
    {
        assign(namesDE[i],which(d1$d == i))
        gene<-get(namesDE[i])
        assign(namegenes[i],new_red_green$gene$SystematicName[gene])
        head(paste("genes_",i,sep=""))
    }
    
    source(paste(workingdirectory,"/Functions/gene_means.R",sep=""))
    source(paste(workingdirectory,"/Functions/FCs.R",sep=""))
    gene_means<- gene_means(x)
    gene_means<-t(gene_means)
    x2<- cbind(x,gene_means) 
    colnames(x2)<-c(match_names,"means")
    row.names(x2)<- new_red_green$gene$SystematicName
    source(paste(workingdirectory,"/Functions/compare_lists.R",sep=""))
    last<- length(x2[1,])

    genes_of_interest<-c()
    
    for(i in 1:length(namegenes))
    {
        genes_vec<-get(namegenes[i])
        genes_of_interest<-c(genes_of_interest,genes_vec)
        pattern<- as.data.frame(genes_vec)
        write.table(pattern, file = paste(workingdirectory,"/",filename,"/pattern",i-1,".txt",sep=""), row.names = FALSE)
    }
    
    genes_of_interest_2<- unique(genes_of_interest)
    genes_of_interest_2<- genes_of_interest_2[which(genes_of_interest_2 != "")]
    end123<-length(genes_of_interest_2)

    genes_of_interest_ALL<- as.data.frame(genes_of_interest_2)

    write.table(genes_of_interest_ALL, file = paste(workingdirectory,"/",filename,"/genesofinterest1.txt",sep=""), row.names = FALSE)
    FC<-FCs(x,gaga_model,expression$groups,all_names,genes_of_interest_ALL,workingdirectory,filename)
    
    #Ortholog Searching

    source(paste(workingdirectory,"/Functions/fold_changes2.R",sep=""))
    
    changes<- fold_changes2(x2, genes_of_interest_2,j)
    changes<- changes[-1,]
    changes2<- as.numeric(changes)
    changes3<- matrix(changes2, ncol = (max(j)+1), byrow = FALSE)
    changes<- changes3 + new_offset
    genes_of_interest_2<- cbind(genes_of_interest_2, changes)
    genes_of_interest_2<-genes_of_interest_2[1:end123,]
    
    anopheles_description_file<- read.table(paste(workingdirectory,"/Ortholog_Searches/anopheles_description_file.txt",sep=""), header = TRUE, na.strings = "NA", fill = TRUE, sep = "\t")
    source(paste(workingdirectory,"/Functions/orthologs_and_descriptions_vectorbase/anopheles_descriptions2.R",sep=""))
    source(paste(workingdirectory,"/Functions/orthologs_and_descriptions_vectorbase/orthologs_with_anopheles2.R",sep=""))
    anopheles_description_file<-as.matrix(anopheles_description_file)
    anoph<- suppressWarnings(anopheles_descriptions2(anopheles_description_file, genes_of_interest_2,j))
    if(length(anoph)>1)
    {
        anoph<- anoph[-1,]
    }
    else
    {
        anoph<-matrix("NONE",2,30)
    }
    
    anopheles_aedes<- read.table(paste(workingdirectory,"/Ortholog_Searches/anopheles_aedes_orthologs.txt",sep=""), header = TRUE, na.strings = "NA", fill = TRUE)
    aedes_descriptions<- read.delim(paste(workingdirectory,"/Ortholog_Searches/aedes_description.txt",sep=""), header = TRUE)
    aedes_anopheles<- anopheles_aedes[,1:4]
    aedes_anopheles<- as.matrix(aedes_anopheles)
    aedes_descriptions<- as.matrix(aedes_descriptions)
    aedes_ortholog_list<- suppressWarnings(orthologs_with_anopheles2(aedes_anopheles, genes_of_interest_2,j))
    
    if(length(aedes_ortholog_list)>1)
    {
        aedes_ortholog_list<-aedes_ortholog_list[-1,]
        source(paste(workingdirectory,"/Functions/orthologs_and_descriptions_vectorbase/orthologs2.R",sep=""))
        aedes_ortholog_descriptions<- suppressWarnings(orthologs2(aedes_ortholog_list, aedes_descriptions,j))
    }
    else
    {
        aedes_ortholog_descriptions<-matrix("NONE",2,30)
    }
    
    if(length(aedes_ortholog_descriptions)>1)
    {
        aedes_ortholog_descriptions<- aedes_ortholog_descriptions[-1,]
        interesting_aedes_orthologs<- as.data.frame(aedes_ortholog_descriptions)
        write.table(interesting_aedes_orthologs, file = paste(workingdirectory,"/",filename,"/aedes_orthologs.txt",sep=""),row.names = FALSE)
    }
    else
    {
        aedes_ortholog_descriptions<-matrix("NONE",2,30)
    }
    
    anopheles_culex<- read.table(paste(workingdirectory,"/Ortholog_Searches/culex_orthologs.txt",sep=""), header = TRUE, na.strings = "NA", fill = TRUE)
    culex_descriptions<- read.delim(paste(workingdirectory,"/Ortholog_Searches/culex_description.txt",sep=""), header = TRUE)
    culex_anopheles<- anopheles_culex[,1:2]
    culex_anopheles<- as.matrix(culex_anopheles)
    culex_descriptions<- as.matrix(culex_descriptions)
    source(paste(workingdirectory,"/Functions/orthologs_and_descriptions_vectorbase/split_for_drosophila.R",sep=""))
    source(paste(workingdirectory,"/Functions/orthologs_and_descriptions_vectorbase/culex_with_anopheles2.R",sep=""))
    source(paste(workingdirectory,"/Functions/orthologs_and_descriptions_vectorbase/match_drosophila2.R",sep=""))
    new_culex<- suppressWarnings(match_drosophila2(culex_anopheles, culex_descriptions))
    
    new_GOI<- split_for_drosophila(genes_of_interest_2)
    new_GOI<-new_GOI[1:end123,]
    write.table(new_GOI, file = paste(workingdirectory,"/",filename,"/GenesListFile.txt",sep=""), row.names = FALSE)
    new_culex<-new_culex[-1,]
    culex<- suppressWarnings(culex_with_anopheles2(new_culex, new_GOI,j))
    
    if(length(culex)>1)
    {
        interesting_culex_orthologs<- as.data.frame(culex)
        write.table(interesting_culex_orthologs, file = paste(workingdirectory,"/",filename,"/culex_orthologs.txt",sep=""), row.names = FALSE)
    }
    else
    {
        culex<-matrix("NONE",2,30)
    }

    
    anopheles_drosophila<- read.table(paste(workingdirectory,"/Ortholog_Searches/drosophila_orthologs.txt",sep=""), header = TRUE, na.strings = "NA", fill = TRUE)
    drosophila_description<- read.delim(paste(workingdirectory,"/Ortholog_Searches/drosophila_description.txt",sep=""), header = TRUE)
    drosophila<- read.table(paste(workingdirectory,"/Ortholog_Searches/drosophila_orthologs.txt",sep=""), header = TRUE, fill = TRUE)
    descriptions<- read.delim(paste(workingdirectory,"/Ortholog_Searches/drosophila_description.txt",sep=""), header = TRUE)
    new_drosophila<- suppressWarnings(match_drosophila2(drosophila, descriptions))
    
    if(length(new_drosophila)>1)
    {
        new_drosophila<- new_drosophila[-1,]
        drosophila<- suppressWarnings(culex_with_anopheles2(new_drosophila, new_GOI,j))
        drosophila_ortholog_descriptions<- drosophila[-1,]   
        interesting_drosophila_orthologs<- as.data.frame(drosophila_ortholog_descriptions)
        write.table(interesting_drosophila_orthologs, file = paste(workingdirectory,"/",filename,"/drosophila_orthologs.txt",sep=""), row.names = FALSE)
    }
    else
    {
        drosophila_ortholog_descriptions<-matrix("NONE",2,30)
    }
    
    source(paste(workingdirectory,"/Functions/orthologs_and_descriptions_vectorbase/full_descriptions2.R",sep=""))
    all_descriptions<- full_descriptions2(genes_of_interest_2,new_GOI, anoph, aedes_ortholog_descriptions, culex, drosophila_ortholog_descriptions,j)
    all_orthologs<- as.data.frame(all_descriptions)
    names(all_orthologs)<-c("An. Gambiae Identifier","An. Gambiae Putative Funciton","Aedes Identifier","Relation of Aedes-Anopheles","Aedes Putative Function","Culex Identifier","Relation of Culex-Anopheles","Culex Putative Function","Drosophila Melanogaster Identifier","Relation of Drosophila-Anopheles","Drosophila Putative Function",rep("Microarray Data",(max(j)+1)))
    write.table(all_orthologs, file = paste(workingdirectory,"/",filename,"/GOI_ortholog_descriptions.txt",sep=""), row.names = FALSE)
    groups_length<-length(groupsvector)
    
    for(i in 2:length(namegenes))
    {
        genes_vec<-get(namegenes[i])
        source(paste(workingdirectory,"/Functions/orthologs_and_descriptions_vectorbase/colouring.R",sep=""))
        source(paste(workingdirectory,"/Functions/orthologs_and_descriptions_vectorbase/colourmatrix.R",sep=""))
        source(paste(workingdirectory,"/Functions/orthologs_and_descriptions_vectorbase/heatmaps.R",sep=""))
        col.map<-colouring(expression$groups)
        colours2<-unique(col.map)
        number<-1:length(colours2)
        col.label<-colourmatrix(colours2,number)
        col.label<-paste(col.label,sep=" ",collapse=" ")
        heatmapinfo<-heatmaps(genes_vec,all_orthologs)
        source(paste(workingdirectory,"/Functions/orthologs_and_descriptions_vectorbase/heatmap_production.R",sep=""))
        namesvec_up_down<-heatmap_production(heatmapinfo,i,groupsvector,patternsvector,col.label,col.map,filename,workingdirectory)
        
        if(length(namesvec_up_down$up)==0)
        {
            up_reg<-1
        }
        else
        {
            up_reg<-namesvec_up_down$up
        }
        
        if(length(namesvec_up_down$down)==0)
        {
            down_reg<-1
        }
        else
        {
            down_reg<-namesvec_up_down$down
        }
        
        if(length(namesvec_up_down$other)==0)
        {
            other_reg<-1
        }
        
        else
        {
            other_reg<-namesvec_up_down$other
        }
        
        namesvec_up_down<-c()
        source(paste(workingdirectory,"/Functions/orthologs_and_descriptions_vectorbase/pathways2.R",sep=""))
        source(paste(workingdirectory,"/Functions/orthologs_and_descriptions_vectorbase/GO.R",sep=""))
        source(paste(workingdirectory,"/Functions/orthologs_and_descriptions_vectorbase/reactome.R",sep=""))
        sig_genes_pathways<-pathways2(genes_vec,workingdirectory,all_orthologs)
        reactome_genes<-reactome(genes_vec,workingdirectory,all_orthologs)
        source(paste(workingdirectory,"/Functions/orthologs_and_descriptions_vectorbase/GO_test2.R",sep=""))
        source(paste(workingdirectory,"/Functions/orthologs_and_descriptions_vectorbase/reactome_test.R",sep=""))
        source(paste(workingdirectory,"/Functions/orthologs_and_descriptions_vectorbase/pathway_test.R",sep=""))
        source(paste(workingdirectory,"/Functions/orthologs_and_descriptions_vectorbase/all_orthologs_search.R",sep=""))
        all_orthologs_search(workingdirectory,genes_vec,filename,i)
        try(pathway_test(sig_genes_pathways,groups_length,col.label,col.map,i,filename,workingdirectory,up_reg,down_reg,other_reg))
        try(reactome_test(reactome_genes,groups_length,col.label,col.map,i,filename,workingdirectory,up_reg,down_reg,other_reg))
        try(GO_test2(genes_vec,new_new_red_green$genes$SystematicName,workingdirectory,i,all_orthologs,groups_length,col.map,col.label,up_reg,down_reg,other_reg,filename),silent=T)
    }
}

 # install.packages("xlsx")
 # install.packages("dplyr")
 # install.packages("tidyr")
library(dplyr)
library(xlsx)
library(tidyr)

##   输入要处理的压缩包   ###
name<- "UP-1.5-1-N.zip"
nchar.tocut<- (nchar(name)-4)
pure.name<- substring(name,1,nchar.tocut)
unzip(name,files = "metascape_result.xlsx")
metascape <- read.xlsx("metascape_result.xlsx",sheetIndex = 2)
metascape <- metascape[,c(-8,-9)]
metascape <- separate(metascape, col = "InTerm_InList", into = c("counts","ALL_counts"), sep = "/", remove = TRUE, convert = FALSE, extra = "warn", fill = "warn")
metascape$counts<- as.numeric(metascape$counts)
metascape$ALL_counts<- as.numeric(metascape$ALL_counts)
metascape$enrich_factor <- (metascape$counts)/(metascape$ALL_counts)
metascape$Qvalue <- 10^(metascape$Log.q.value.)
metascape<- metascape[order(metascape$Log.q.value., decreasing= F), ]
# df[,-which(names(metascape)%in%c("1_Summary","2_Summary","3_Summary","4_Summary","5_Summary")]
# metascape <- subset(metascape,select = -c("1_Summary","2_Summary","3_Summary","4_Summary","5_Summary") )
# metascape1 <- subset(metascape,which(metascape,GroupId==c("1_Summary","2_Summary","3_Summary","4_Summary","5_Summary")))

#清除Summary
if (T) {
  #清除Summary
metascape <- metascape[!grepl("Summary", metascape$GroupID),]
}

if (T) {
  #清除多余重复Term
  metascape <- metascape[!grepl("hsa", metascape$Term),]
  metascape <- metascape[!grepl("HSA", metascape$Term),]
  metascape <- metascape[!grepl("ko", metascape$Term),]
}
metascape.GO_BP <- na.omit(subset(metascape,Category =="GO Biological Processes",)[c(1:20),])
metascape.GO_MF <-  na.omit(subset(metascape,Category =="GO Molecular Functions")[c(1:20),])
metascape.GO_CC <-  na.omit(subset(metascape,Category =="GO Cellular Components")[c(1:20),])
metascape.KEGG <-  na.omit(subset(metascape,Category =="KEGG Pathway")[c(1:20),])
metascape.REACTOME <-  na.omit(subset(metascape,Category =="Reactome Gene Sets")[c(1:20),])
metascape.CORUM <-  na.omit(subset(metascape,Category =="CORUM")[c(1:20),])
metascape.IS <-  na.omit(subset(metascape,Category =="Immunologic Signatures")[c(1:20),])
metascape.CGP <-  na.omit(subset(metascape,Category =="Chemical And Genetic Perturbations")[c(1:20),])
metascape.KEGG.SC <-  na.omit(subset(metascape,Category =="KEGG Structural Complexes")[c(1:20),])
metascape.BGS <-  na.omit(subset(metascape,Category =="BioCarta Gene Sets")[c(1:20),])
metascape.HGS <-  na.omit(subset(metascape,Category =="Hallmark Gene Sets")[c(1:20),])


#三个画在一起
if(F){
  ##加载ggplot2包
  library(ggplot2)
  ##读取刚才保存的富集分析结果文件
  # goinput <- read.table("table.txt",header = T,sep = "\t")##header = T 第一行是表头，sep = "\t"表示以制表符（tab）分割
  ##利用levels来设置factor中的顺序，保证最后出图时按照我们之前排好的顺序排列
  metascape.GO <- rbind (metascape.GO_BP,metascape.GO_MF)
  metascape.GO <- rbind (metascape.GO,metascape.GO_CC)
  x=metascape.GO$enrich_factor
  y=factor(metascape.GO$Description,levels = metascape.GO$Description)
  ##先出一个框架
  p = ggplot(metascape,aes(x,y))
  ##数据特征包括bubble大小来源为匹配到这个term上的基因数，颜色为qvalue，颜色变化为从低到高:"SpringGreen"到"DeepPink"而类型则是class，
  p1 = p + geom_point(aes(size=Qvalue,color=-1*log(counts),shape=Category,))+
    scale_color_gradient(low = "SpringGreen", high = "DeepPink")
  ##设置横纵坐标名字，标题，legend名字
  p2 = p1 + labs(color=expression(-log[10](Qvalue)),
                 size="Gene Number",
                 x="enrich_factor",
                 y="Go_term",
                 title="")
  ##搞个主题，把边框画出来也可以通过ggplot Theme assistant修改
  p3 = p2 +theme_bw() 
  p3
}

#三个不画在一起

##########
#   BP   #
##########

if(T){
  ##加载ggplot2包
  library(ggplot2)
  ##读取刚才保存的富集分析结果文件
  # goinput <- read.table("table.txt",header = T,sep = "\t")##header = T 第一行是表头，sep = "\t"表示以制表符（tab）分割
  ##利用levels来设置factor中的顺序，保证最后出图时按照我们之前排好的顺序排列
  # library(dplyr)
  # goinput.bp <- subset(goinput,Class == "GO Biological Processes")
  x=metascape.GO_BP$enrich_factor
  y=factor(metascape.GO_BP$Description,levels = metascape.GO_BP$Description)
  ##先出一个框架
  p = ggplot(metascape.GO_BP,aes(x,y))
  ##数据特征包括bubble大小来源为匹配到这个term上的基因数，颜色为qvalue，颜色变化为从低到高:"SpringGreen"到"DeepPink"而类型则是class，
  p1 = p + geom_point(aes(size=counts,color=Qvalue,))+
    scale_color_gradient(low = "#96BCFF", high = "#0B2D69")
  ##设置横纵坐标名字，标题，legend名字
  p2 = p1 + labs(color=expression(-log[10](Qvalue)),
                 size="Gene Number",
                 x="Enrich Factor",
                 y=" ",
                 title=" ")
  ##搞个主题，把边框画出来也可以通过ggplot Theme assistant修改
  p3 = p2 +theme_bw() 
  p3
  ggsave(p3,filename = paste(pure.name,"GO_BP.pdf",sep = "-"),width = 7.52,height = 6.52)
}

##########
#   CC   #
##########

if(T){
  ##加载ggplot2包
  library(ggplot2)
  ##读取刚才保存的富集分析结果文件
  # goinput <- read.table("table.txt",header = T,sep = "\t")##header = T 第一行是表头，sep = "\t"表示以制表符（tab）分割
  # ##利用levels来设置factor中的顺序，保证最后出图时按照我们之前排好的顺序排列
  library(dplyr)
  # metascape.GO_CC <- subset(metascape.GO_CC,Class == "GO Cellular Components")
  x=metascape.GO_CC$enrich_factor
  y=factor(metascape.GO_CC$Description,levels = metascape.GO_CC$Description)
  ##先出一个框架
  p = ggplot(metascape.GO_CC,aes(x,y))
  ##数据特征包括bubble大小来源为匹配到这个term上的基因数，颜色为qvalue，颜色变化为从低到高:"SpringGreen"到"DeepPink"而类型则是class，
  p1 = p + geom_point(aes(size=counts,color=Qvalue,))+
    scale_color_gradient(low = "#F4DBFE", high = "#510E6F")
  ##设置横纵坐标名字，标题，legend名字
  p2 = p1 + labs(color=expression(-log[10](Qvalue)),
                 size="Gene Number",
                 x="Enrich Factor",
                 y=" ",
                 title=" ")
  ##搞个主题，把边框画出来也可以通过ggplot Theme assistant修改
  p3 = p2 +theme_bw() 
  p3
  ggsave(p3,filename = paste(pure.name,"GO_CC.pdf",sep = "-"),width = 6.52,height = 6.52)
}

##########
#   MF   #
##########

if(T){
  ##加载ggplot2包
  library(ggplot2)
  ##读取刚才保存的富集分析结果文件
  # goinput <- read.table("table.txt",header = T,sep = "\t")##header = T 第一行是表头，sep = "\t"表示以制表符（tab）分割
  ##利用levels来设置factor中的顺序，保证最后出图时按照我们之前排好的顺序排列
  library(dplyr)
  # metascape.GO_MF<- subset(goinput,Class == "GO Molecular Functions")
  x=metascape.GO_MF$enrich_factor
  y=factor(metascape.GO_MF$Description,levels = metascape.GO_MF$Description)
  ##先出一个框架
  p = ggplot(metascape.GO_MF,aes(x,y))
  ##数据特征包括bubble大小来源为匹配到这个term上的基因数，颜色为qvalue，颜色变化为从低到高:"SpringGreen"到"DeepPink"而类型则是class，
  p1 = p + geom_point(aes(size=counts,color=Qvalue,))+
    scale_color_gradient(low = "#AEDFB2", high = "#0E4E13")
  ##设置横纵坐标名字，标题，legend名字
  p2 = p1 + labs(color=expression(-log[10](Qvalue)),
                 size="Gene Number",
                 x="Enrich Factor",
                 y=" ",
                 title=" ")
  ##搞个主题，把边框画出来也可以通过ggplot Theme assistant修改
  p3 = p2 +theme_bw() 
  p3
  ggsave(p3,filename = paste(pure.name,"GO_MF.pdf",sep = "-"),width = 6.52,height = 6.52)
}

##########
#  KEGG  #
##########

if(T){
  ##加载ggplot2包
  library(ggplot2)
  ##读取刚才保存的富集分析结果文件
  # goinput <- read.table("table.txt",header = T,sep = "\t")##header = T 第一行是表头，sep = "\t"表示以制表符（tab）分割
  # ##利用levels来设置factor中的顺序，保证最后出图时按照我们之前排好的顺序排列
  library(dplyr)
  # metascape.KEGG<- subset(goinput,Class == "KEGG Pathway")
  x=metascape.KEGG$enrich_factor
  y=factor(metascape.KEGG$Description,levels = metascape.KEGG$Description)
  ##先出一个框架
  p = ggplot(metascape.KEGG,aes(x,y))
  ##数据特征包括bubble大小来源为匹配到这个term上的基因数，颜色为qvalue，颜色变化为从低到高:"SpringGreen"到"DeepPink"而类型则是class，
  p1 = p + geom_point(aes(size=counts,color=Qvalue,))+
    scale_color_gradient(low = "#E2BABA", high = "#FF0606")
  ##设置横纵坐标名字，标题，legend名字
  p2 = p1 + labs(color=expression(-log[10](Qvalue)),
                 size="Gene Number",
                 x="Enrich Factor",
                 y=" ",
                 title=" ")
  ##搞个主题，把边框画出来也可以通过ggplot Theme assistant修改
  p3 = p2 +theme_bw() 
  p3
  ggsave(p3,filename = paste(pure.name,"KEGG.pdf",sep = "-"),width = 6,height = 5.52)
}

##############
#  Reactome  #
##############

if(T){
  ##加载ggplot2包
  library(ggplot2)
  ##读取刚才保存的富集分析结果文件
  # goinput <- read.table("table.txt",header = T,sep = "\t")##header = T 第一行是表头，sep = "\t"表示以制表符（tab）分割 
  ##利用levels来设置factor中的顺序，保证最后出图时按照我们之前排好的顺序排列
  library(dplyr)
  # metascape.REACTOME<- subset(goinput,Class == "Reactome Gene Sets")
  x=metascape.REACTOME$enrich_factor
  y=factor(metascape.REACTOME$Description,levels = metascape.REACTOME$Description)
  ##先出一个框架
  p = ggplot(metascape.REACTOME,aes(x,y))
  ##数据特征包括bubble大小来源为匹配到这个term上的基因数，颜色为qvalue，颜色变化为从低到高:"SpringGreen"到"DeepPink"而类型则是class，
  p1 = p + geom_point(aes(size=counts,color=Qvalue,))+
    scale_color_gradient(low = "#CBCAF1", high = "#1710F6")
  ##设置横纵坐标名字，标题，legend名字
  p2 = p1 + labs(color=expression(-log[10](Qvalue)),
                 size="Gene Number",
                 x="Enrich Factor",
                 y=" ",
                 title=" ")
  ##搞个主题，把边框画出来也可以通过ggplot Theme assistant修改
  p3 = p2 +theme_bw() 
  p3
  ggsave(p3,filename = paste(pure.name,"Reactome.pdf",sep = "-"),width = 7.52,height = 6.52)
}
############################
#  Immunologic Signatures  #
############################


if(F){
  ##加载ggplot2包
  library(ggplot2)
  ##读取刚才保存的富集分析结果文件
  # goinput <- read.table("table.txt",header = T,sep = "\t")##header = T 第一行是表头，sep = "\t"表示以制表符（tab）分割 
  ##利用levels来设置factor中的顺序，保证最后出图时按照我们之前排好的顺序排列
  library(dplyr)
  # metascape.REACTOME<- subset(goinput,Class == "Immunologic Signatures")
  x=metascape.IS$enrich_factor
  y=factor(metascape.IS$Description,levels = metascape.IS$Description)
  ##先出一个框架
  p = ggplot(metascape.IS,aes(x,y))
  ##数据特征包括bubble大小来源为匹配到这个term上的基因数，颜色为qvalue，颜色变化为从低到高:"SpringGreen"到"DeepPink"而类型则是class，
  p1 = p + geom_point(aes(size=counts,color=Qvalue,))+
    scale_color_gradient(low = "#CBCAF1", high = "#1710F6")
  ##设置横纵坐标名字，标题，legend名字
  p2 = p1 + labs(color=expression(-log[10](Qvalue)),
                 size="Gene Number",
                 x="Enrich Factor",
                 y=" ",
                 title=" ")
  ##搞个主题，把边框画出来也可以通过ggplot Theme assistant修改
  p3 = p2 +theme_bw() 
  p3
  ggsave(p3,filename = paste(pure.name,"Immunologic Signatures.pdf",sep = "-"),width = 7.52,height = 6.52)
}

########################################
#  Chemical And Genetic Perturbations  #
########################################


if(F){
  ##加载ggplot2包
  library(ggplot2)
  ##读取刚才保存的富集分析结果文件
  # goinput <- read.table("table.txt",header = T,sep = "\t")##header = T 第一行是表头，sep = "\t"表示以制表符（tab）分割 
  ##利用levels来设置factor中的顺序，保证最后出图时按照我们之前排好的顺序排列
  library(dplyr)
  # metascape.REACTOME<- subset(goinput,Class == "Immunologic Signatures")
  x=metascape.CGP$enrich_factor
  y=factor(metascape.CGP$Description,levels = metascape.CGP$Description)
  ##先出一个框架
  p = ggplot(metascape.CGP,aes(x,y))
  ##数据特征包括bubble大小来源为匹配到这个term上的基因数，颜色为qvalue，颜色变化为从低到高:"SpringGreen"到"DeepPink"而类型则是class，
  p1 = p + geom_point(aes(size=counts,color=Qvalue,))+
    scale_color_gradient(low = "#CBCAF1", high = "#1710F6")
  ##设置横纵坐标名字，标题，legend名字
  p2 = p1 + labs(color=expression(-log[10](Qvalue)),
                 size="Gene Number",
                 x="Enrich Factor",
                 y=" ",
                 title=" ")
  ##搞个主题，把边框画出来也可以通过ggplot Theme assistant修改
  p3 = p2 +theme_bw() 
  p3
  ggsave(p3,filename = paste(pure.name,"Chemical And Genetic Perturbations.pdf",sep = "-"),width = 7.52,height = 6.52)
}


########################################
#  KEGG Structural Complexes  #
########################################

if(F){
  ##加载ggplot2包
  library(ggplot2)
  ##读取刚才保存的富集分析结果文件
  # goinput <- read.table("table.txt",header = T,sep = "\t")##header = T 第一行是表头，sep = "\t"表示以制表符（tab）分割 
  ##利用levels来设置factor中的顺序，保证最后出图时按照我们之前排好的顺序排列
  library(dplyr)
  # metascape.REACTOME<- subset(goinput,Class == "Immunologic Signatures")
  x=metascape.KEGG.SC$enrich_factor
  y=factor(metascape.KEGG.SC$Description,levels = metascape.KEGG.SC$Description)
  ##先出一个框架
  p = ggplot(metascape.KEGG.SC,aes(x,y))
  ##数据特征包括bubble大小来源为匹配到这个term上的基因数，颜色为qvalue，颜色变化为从低到高:"SpringGreen"到"DeepPink"而类型则是class，
  p1 = p + geom_point(aes(size=counts,color=Qvalue,))+
    scale_color_gradient(low = "#CBCAF1", high = "#1710F6")
  ##设置横纵坐标名字，标题，legend名字
  p2 = p1 + labs(color=expression(-log[10](Qvalue)),
                 size="Gene Number",
                 x="Enrich Factor",
                 y=" ",
                 title=" ")
  ##搞个主题，把边框画出来也可以通过ggplot Theme assistant修改
  p3 = p2 +theme_bw() 
  p3
  ggsave(p3,filename = paste(pure.name,"KEGG Structural Complexes.pdf",sep = "-"),width = 6.52,height = 6.52)
}

########################
#  BioCarta Gene Sets  #
########################

if(F){
  ##加载ggplot2包
  library(ggplot2)
  ##读取刚才保存的富集分析结果文件
  # goinput <- read.table("table.txt",header = T,sep = "\t")##header = T 第一行是表头，sep = "\t"表示以制表符（tab）分割 
  ##利用levels来设置factor中的顺序，保证最后出图时按照我们之前排好的顺序排列
  library(dplyr)
  # metascape.REACTOME<- subset(goinput,Class == "Immunologic Signatures")
  x=metascape.BGS$enrich_factor
  y=factor(metascape.BGS$Description,levels = metascape.BGS$Description)
  ##先出一个框架
  p = ggplot(metascape.BGS,aes(x,y))
  ##数据特征包括bubble大小来源为匹配到这个term上的基因数，颜色为qvalue，颜色变化为从低到高:"SpringGreen"到"DeepPink"而类型则是class，
  p1 = p + geom_point(aes(size=counts,color=Qvalue,))+
    scale_color_gradient(low = "#CBCAF1", high = "#1710F6")
  ##设置横纵坐标名字，标题，legend名字
  p2 = p1 + labs(color=expression(-log[10](Qvalue)),
                 size="Gene Number",
                 x="Enrich Factor",
                 y=" ",
                 title=" ")
  ##搞个主题，把边框画出来也可以通过ggplot Theme assistant修改
  p3 = p2 +theme_bw() 
  p3
  ggsave(p3,filename = paste(pure.name,"BioCarta Gene Sets.pdf",sep = "-"),width = 6.52,height = 6.52)
}

########################
#  BioCarta Gene Sets  #
########################

if(F){
  ##加载ggplot2包
  library(ggplot2)
  ##读取刚才保存的富集分析结果文件
  # goinput <- read.table("table.txt",header = T,sep = "\t")##header = T 第一行是表头，sep = "\t"表示以制表符（tab）分割 
  ##利用levels来设置factor中的顺序，保证最后出图时按照我们之前排好的顺序排列
  library(dplyr)
  # metascape.REACTOME<- subset(goinput,Class == "Immunologic Signatures")
  x=metascape.HGS$enrich_factor
  y=factor(metascape.HGS$Description,levels = metascape.HGS$Description)
  ##先出一个框架
  p = ggplot(metascape.HGS,aes(x,y))
  ##数据特征包括bubble大小来源为匹配到这个term上的基因数，颜色为qvalue，颜色变化为从低到高:"SpringGreen"到"DeepPink"而类型则是class，
  p1 = p + geom_point(aes(size=counts,color=Qvalue,))+
    scale_color_gradient(low = "#CBCAF1", high = "#1710F6")
  ##设置横纵坐标名字，标题，legend名字
  p2 = p1 + labs(color=expression(-log[10](Qvalue)),
                 size="Gene Number",
                 x="Enrich Factor",
                 y=" ",
                 title=" ")
  ##搞个主题，把边框画出来也可以通过ggplot Theme assistant修改
  p3 = p2 +theme_bw() 
  p3
  ggsave(p3,filename = paste(pure.name,"Hallmark Gene Sets.pdf",sep = "-"),width = 6.52,height = 6.52)
}
















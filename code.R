

####gene####

####使用ggpubr包绘制带组内差异比较的柱形图或箱线图
rm(list=ls())
setwd("E:/X-14952B/gene/")


####水杨酸途径24 h####

SA1 <- read.delim("SA_24h.txt", sep = '\t', stringsAsFactors = FALSE)
SA1

####组内涉及多组比较时的示例
#####大组是处理类型，小组是处理梯度

##使用 ggpubr 包绘制统计图
library(ggpubr)

#使用 ggpubr 包的方法绘制柱形图
p2_bar <- ggbarplot(SA1, x = 'gene', y = 'value', fill = 'treatment', add = 'mean_sd', 
                    color = 'gray10', palette= "lancet", position = position_dodge(0.8), width = 0.6, size = 1, legend = 'right') +
  labs(x = '24 h', y = 'Relative expression',expand = FALSE, fill = 'Treatment')

p2_bar

##使用 rstatix 包进行两组比较的统计分析
library(rstatix)

#如果期望使用参数检验，以 t 检验为例
stat_t <- t_test(group_by(SA1, gene), value~treatment)  #获取组内两两比较的 p 值
stat_t <- adjust_pvalue(stat_t, method = 'fdr')  #多组比较时，可能需要添加 p 值校正
stat_t <- add_significance(stat_t, 'p.adj')  #根据 p 校正值添加显著性标记 * 符号

#根据 t 检验的结果，在柱形图中添加 p 校正值或者显著性标记 * 符号
stat_t.test <-  add_xy_position(stat_t, x = 'gene', dodge = 0.8)
p1 <- p2_bar + stat_pvalue_manual(stat_t.test, label = 'p.adj', tip.length = 0.02)  #显示为 p 校正值
p11 <- p2_bar + stat_pvalue_manual(stat_t.test, label = 'p.adj.signif', tip.length = 0.05)  #显示为 * 符号
#p 校正值和 * 符号一起带着
p111 <- p2_bar + stat_pvalue_manual(stat_t.test, label = '{p.adj} {p.adj.signif}', tip.length = 0.05) 
#输出图像。
ggsave("SA_24h.pdf", p1, width = 10, height = 8, units="cm")
ggsave("SA_24h.tiff",p1, width = 10, height = 8, units = "cm")
ggsave("SA_24h.png",p1, width = 10, height = 8, units = "cm")



####水杨酸途径72 h####
SA2 <- read.delim("SA_72h.txt", sep = '\t', stringsAsFactors = FALSE)
SA2

####组内涉及多组比较时的示例
#####大组是处理类型，小组是处理梯度

##使用 ggpubr 包绘制统计图
library(ggpubr)

#使用 ggpubr 包的方法绘制柱形图
p2_bar <- ggbarplot(SA2, x = 'gene', y = 'value', fill = 'treatment', add = 'mean_sd', 
                    color = 'gray10', palette= "lancet", position = position_dodge(0.8), width = 0.6, size = 1, legend = 'right') +
  labs(x = '72 h', y = 'Relative expression', fill = 'Treatment')

p2_bar

##使用 rstatix 包进行两组比较的统计分析
library(rstatix)

#如果期望使用参数检验，以 t 检验为例
stat_t <- t_test(group_by(SA2, gene), value~treatment)  #获取组内两两比较的 p 值
stat_t <- adjust_pvalue(stat_t, method = 'fdr')  #多组比较时，可能需要添加 p 值校正
stat_t <- add_significance(stat_t, 'p.adj')  #根据 p 校正值添加显著性标记 * 符号

#根据 t 检验的结果，在柱形图中添加 p 校正值或者显著性标记 * 符号
stat_t.test <-  add_xy_position(stat_t, x = 'gene', dodge = 0.8)
p2 <- p2_bar + stat_pvalue_manual(stat_t.test, label = 'p.adj', tip.length = 0.02)  #显示为 p 校正值
p22 <- p2_bar + stat_pvalue_manual(stat_t.test, label = 'p.adj.signif', tip.length = 0.05)  #显示为 * 符号
#p 校正值和 * 符号一起带着
p222 <- p2_bar + stat_pvalue_manual(stat_t.test, label = '{p.adj} {p.adj.signif}', tip.length = 0.05) 
#输出图像。
ggsave("SA_72h.pdf",p2, width = 10, height = 8, units="cm")
ggsave("SA_72h.tiff",p2, width = 10, height = 8, units = "cm")
ggsave("SA_72h.png",p2, width = 10, height = 8, units = "cm")



####茉莉酸途径24 h####
JA1 <- read.delim("JA_24h.txt", sep = '\t', stringsAsFactors = FALSE)
JA1

####组内涉及多组比较时的示例
#####大组是处理类型，小组是处理梯度

##使用 ggpubr 包绘制统计图
library(ggpubr)

#使用 ggpubr 包的方法绘制柱形图
p2_bar <- ggbarplot(JA1, x = 'gene', y = 'value', fill = 'treatment', add = 'mean_sd', 
                    color = 'gray10', palette= "lancet", position = position_dodge(0.8), width = 0.6, size = 1, legend = 'right') +
  labs(x = '24 h', y = 'Relative expression', fill = 'Treatment')

p2_bar

##使用 rstatix 包进行两组比较的统计分析
library(rstatix)

#如果期望使用参数检验，以 t 检验为例
stat_t <- t_test(group_by(JA1, gene), value~treatment)  #获取组内两两比较的 p 值
stat_t <- adjust_pvalue(stat_t, method = 'fdr')  #多组比较时，可能需要添加 p 值校正
stat_t <- add_significance(stat_t, 'p.adj')  #根据 p 校正值添加显著性标记 * 符号

#根据 t 检验的结果，在柱形图中添加 p 校正值或者显著性标记 * 符号
stat_t.test <-  add_xy_position(stat_t, x = 'gene', dodge = 0.8)
p3 <- p2_bar + stat_pvalue_manual(stat_t.test, label = 'p.adj', tip.length = 0.02)  #显示为 p 校正值
p33 <- p2_bar + stat_pvalue_manual(stat_t.test, label = 'p.adj.signif', tip.length = 0.05)  #显示为 * 符号
#p 校正值和 * 符号一起带着
p333 <- p2_bar + stat_pvalue_manual(stat_t.test, label = '{p.adj} {p.adj.signif}', tip.length = 0.05) 
#输出图像。
ggsave("JA_24h.pdf",p3, width = 10, height = 8, units="cm")
ggsave("JA_24h.tiff",p3, width = 10, height = 8, units = "cm")
ggsave("JA_24h.png",p3, width = 10, height = 8, units = "cm")


####茉莉酸途径72 h####
JA2 <- read.delim("JA_72h.txt", sep = '\t', stringsAsFactors = FALSE)
JA2

####组内涉及多组比较时的示例
#####大组是处理类型，小组是处理梯度

##使用 ggpubr 包绘制统计图
library(ggpubr)

#使用 ggpubr 包的方法绘制柱形图
p2_bar <- ggbarplot(JA2, x = 'gene', y = 'value', fill = 'treatment', add = 'mean_sd', 
                    color = 'gray10', palette= "lancet", position = position_dodge(0.8), width = 0.6, size = 1, legend = 'right') +
  labs(x = '72 h', y = 'Relative expression', fill = 'Treatment')

p2_bar

##使用 rstatix 包进行两组比较的统计分析
library(rstatix)

#如果期望使用参数检验，以 t 检验为例
stat_t <- t_test(group_by(JA2, gene), value~treatment)  #获取组内两两比较的 p 值
stat_t <- adjust_pvalue(stat_t, method = 'fdr')  #多组比较时，可能需要添加 p 值校正
stat_t <- add_significance(stat_t, 'p.adj')  #根据 p 校正值添加显著性标记 * 符号

#根据 t 检验的结果，在柱形图中添加 p 校正值或者显著性标记 * 符号
stat_t.test <-  add_xy_position(stat_t, x = 'gene', dodge = 0.8)
p4 <- p2_bar + stat_pvalue_manual(stat_t.test, label = 'p.adj', tip.length = 0.02)  #显示为 p 校正值
p44 <- p2_bar + stat_pvalue_manual(stat_t.test, label = 'p.adj.signif', tip.length = 0.05)  #显示为 * 符号
#p 校正值和 * 符号一起带着
p444 <- p2_bar + stat_pvalue_manual(stat_t.test, label = '{p.adj} {p.adj.signif}', tip.length = 0.05) 
#输出图像。
ggsave("JA_72h.pdf",p4, width = 10, height = 8, units="cm")
ggsave("JA_72h.tiff",p4, width = 10, height = 8, units = "cm")
ggsave("JA_72h.png",p4, width = 10, height = 8, units = "cm")











####微生物组分析####
rm(list=ls())
####α-diversity####

setwd("E:/X-14952B/otus/α-diversity/")
#首先导入数据文件


#1.微生物群落α-多样性指数的计算
#1.1R包的安装及加载
#install.packages("vegan")
#install.packages("picante")
library(vegan)
library(picante)
#注：picante 包加载时默认同时加载 vegan，如果加载了它，可省略“library(vegan)”这一步

#1.2读入物种数据
otu <- read.delim('./OTU_Tables.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
otu
otu <- t(otu)
otu

alpha <- function(x, tree = NULL, base = exp(1)) {
  est <- estimateR(x)
  Richness <- est[1, ]
  Chao1 <- est[2, ]
  ACE <- est[4, ]
  Shannon <- diversity(x, index = 'shannon', base = base)
  Simpson <- diversity(x, index = 'simpson')    #Gini-Simpson 指数
  Pielou <- Shannon / log(Richness, base)
  goods_coverage <- 1 - rowSums(x == 1) / rowSums(x)
  
  result <- data.frame(Richness, Shannon, Simpson, Pielou, Chao1, ACE, goods_coverage)
  if (!is.null(tree)) {
    PD_whole_tree <- pd(x, tree, include.root = FALSE)[1]
    names(PD_whole_tree) <- 'PD_whole_tree'
    result <- cbind(result, PD_whole_tree)
  }
  result
}

#1.4生成阿尔法多样性指数表格并保存文件
alpha_all <- alpha(otu, base = 2)
write.table(alpha_all,file="alpha-indexsummary.txt",sep="\t",col.names = NA)
dat <- read.delim('./alpha-indexsummary.txt')#读入α多样性数据
head(dat, n = 3)
design <- read.delim('./mapping.txt')#读入试验设计文件
head(design, n = 3)

dat <- design %>%
  dplyr::select(OTU.ID,group) %>%
  inner_join(dat,by='OTU.ID')#按照OTU.ID内连接文件
head(dat, n = 3)
write.table(dat,file="alpha_all.xls",sep="\t",col.names = NA)


#2.α-多样性指数统计学分析
#表格格式转换
library(tidyverse)
library(ggprism)

dat <- read.delim('alpha_all.xls', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE) #读入α多样性数据
dat
dat <- gather(dat,alpha,v,-(OTU.ID:group))#宽数据变长数据
head(dat, n = 3)

#2.2α多样性指数统计学分析

kruskalmc_compare1 <- function(data,group,compare,value){
  library(multcompView)
  library(pgirmess)##多组两辆比较函数用到的包
  library(multcomp)
  a <- data.frame(stringsAsFactors = F)
  type <- unique(data[,group])
  for (i in type)
  {
    g1=compare
    sub_dat <- data[data[,group]==i,]
    names(sub_dat)[names(sub_dat)==compare] <- 'g1'
    names(sub_dat)[names(sub_dat)==value] <- 'value'
    sub_dat$g1 <- factor(sub_dat$g1)
    options(warn = -1)
    
    k <- kruskalmc(value ~ g1, data=sub_dat, probs=0.05)
    dif <- k$dif.com[['difference']]
    names(dif) <- rownames(k$dif.com)
    difL <- multcompLetters(dif)
    label <- data.frame(difL['Letters'], stringsAsFactors = FALSE)
    label$compare = rownames(label)
    label$type <- i
    
    mean_sd <- merge(aggregate(sub_dat[['value']],by=list(sub_dat[,'g1']),FUN=sd),
                     aggregate(sub_dat[['value']],by=list(sub_dat[,'g1']),FUN=mean),by="Group.1"
    )
    names(mean_sd) <- c('compare','std','mean')
    a <- rbind(a,merge(mean_sd,label,by='compare'))
  }
  names(a) <- c(compare,'std','mean','Letters',group)
  return(a)
}

df1 <- kruskalmc_compare1(dat,'alpha','group','v')
df1######字母是反着标的(a<b<c)



#导出α-多样性统计学检验结果
write.table(df1,file="alpha-kruskaltest.xls",sep="\t",col.names = NA)

#α-多样性出图
#所有α多样性指数统一出图
p1 = ggplot(dat)+geom_boxplot(aes(x=group,y=v,fill=group))+
  geom_text(data=df1,aes(x=group,y=mean+5.0*std,label=Letters),size=4,position=position_dodge(0.2))+
  facet_wrap(~alpha,scales = "free_y")+ labs(x='',y='AlphaDiv')+
  ggprism::theme_prism()+theme(axis.text.x = element_text(angle = 45))
p1

#图形整体调整
mytheme<-theme_bw()+theme(axis.title=element_text(size=12),axis.text = element_text(size = 12),
                          panel.grid.major = element_line(color = "white"),
                          panel.grid.minor = element_line(color = "white"),
                          axis.text.x = element_text(size = 12, angle=30,vjust = 0.5,hjust = 0.5,color = "black"),
                          axis.text.y = element_text(size=12,color = "black"),
                          legend.position = "none",
                          plot.margin = unit(c(0.5, 0.5, 0.5, 0.5),"cm"))
p11 <- p1+mytheme
p11

ggsave("./alpha.pdf", p11, width = 350, height = 200, units = "mm")
ggsave("./alpha.png", p11, width = 350, height = 200, units = "mm")
ggsave("./alpha.tiff", p11, width = 350, height = 200, units = "mm")


#3.α多样性指数挑选后的单独出图
dat2 <- read.delim('alpha_RSS.txt', row.names = 1, sep = '\t') #读入α多样性数据
dat2
dat2 <- gather(dat2,alpha,v,-(OTU.ID:group))#宽数据变长数据
head(dat2, n = 3)


kruskalmc_compare1 <- function(data,group,compare,value){
  library(multcompView)
  library(pgirmess)##多组两辆比较函数用到的包
  library(multcomp)
  a <- data.frame(stringsAsFactors = F)
  type <- unique(data[,group])
  for (i in type)
  {
    g1=compare
    sub_dat <- data[data[,group]==i,]
    names(sub_dat)[names(sub_dat)==compare] <- 'g1'
    names(sub_dat)[names(sub_dat)==value] <- 'value'
    sub_dat$g1 <- factor(sub_dat$g1)
    options(warn = -1)
    
    k <- kruskalmc(value ~ g1, data=sub_dat, probs=0.05)
    dif <- k$dif.com[['difference']]
    names(dif) <- rownames(k$dif.com)
    difL <- multcompLetters(dif)
    label <- data.frame(difL['Letters'], stringsAsFactors = FALSE)
    label$compare = rownames(label)
    label$type <- i
    
    mean_sd <- merge(aggregate(sub_dat[['value']],by=list(sub_dat[,'g1']),FUN=sd),
                     aggregate(sub_dat[['value']],by=list(sub_dat[,'g1']),FUN=mean),by="Group.1"
    )
    names(mean_sd) <- c('compare','std','mean')
    a <- rbind(a,merge(mean_sd,label,by='compare'))
  }
  names(a) <- c(compare,'std','mean','Letters',group)
  return(a)
}

df2 <- kruskalmc_compare1(dat2,'alpha','group','v')
df2######字母是反着标的(a<b<c)

#导出α-多样性统计学检验结果
write.table(df2,file="alpha_RSS_kruskaltest.xls",sep="\t",col.names = NA)

#α-多样性出图
#所有α多样性指数统一出图
p2 = ggplot(dat2)+geom_boxplot(aes(x=group,y=v,fill=group), size=0.5, width=0.4)+
  geom_text(data=df2,aes(x=group,y=mean+2.0*std,label=Letters),size=4,color = "black",fontface = "bold",position=position_dodge(0.2))+
  facet_wrap(~alpha,nrow=1,scales = "free_y")+ labs(x='',y='AlphaDiv')+
  labs(x='',y='AlphaDiv')+
  scale_fill_manual(values =c('#CD5B45', '#228B22', '#00688B',"orange"))+
  ggprism::theme_prism()+theme(axis.text.x = element_text(angle = 45))
p2

#图形整体调整
mytheme<-theme_bw()+theme(axis.title=element_text(size=14),axis.text = element_text(size = 14),
                          panel.grid.major = element_line(color = "white"),
                          panel.grid.minor = element_line(color = "white"),
                          axis.ticks = element_line(color='black'),
                          axis.line = element_line(colour = "black"),
                          axis.text.x = element_text(size = 14, angle=45,vjust = 0.5,hjust = 0.5,color = "black"),
                          axis.text.y = element_text(size=14,color = "black"),
                          legend.position = "",
                          plot.margin = unit(c(0.5, 0.5, 0.5, 0.5),"cm"))


p3<-p2+mytheme
p3
ggsave("./alpha_RSS.pdf", p3, width = 120, height = 80, units = "mm")
ggsave("./alpha_RSS.png", p3, width = 120, height = 80, units = "mm")
ggsave("./alpha_RSS.tiff", p3, width = 120, height = 80, units = "mm")

#其他α多样性指数
dat3 <- read.delim('alpha_others.txt', row.names = 1, sep = '\t') #读入α多样性数据
dat3
dat3 <- gather(dat3,alpha,v,-(OTU.ID:group))#宽数据变长数据
head(dat3, n = 3)


kruskalmc_compare1 <- function(data,group,compare,value){
  library(multcompView)
  library(pgirmess)##多组两辆比较函数用到的包
  library(multcomp)
  a <- data.frame(stringsAsFactors = F)
  type <- unique(data[,group])
  for (i in type)
  {
    g1=compare
    sub_dat <- data[data[,group]==i,]
    names(sub_dat)[names(sub_dat)==compare] <- 'g1'
    names(sub_dat)[names(sub_dat)==value] <- 'value'
    sub_dat$g1 <- factor(sub_dat$g1)
    options(warn = -1)
    
    k <- kruskalmc(value ~ g1, data=sub_dat, probs=0.05)
    dif <- k$dif.com[['difference']]
    names(dif) <- rownames(k$dif.com)
    difL <- multcompLetters(dif)
    label <- data.frame(difL['Letters'], stringsAsFactors = FALSE)
    label$compare = rownames(label)
    label$type <- i
    
    mean_sd <- merge(aggregate(sub_dat[['value']],by=list(sub_dat[,'g1']),FUN=sd),
                     aggregate(sub_dat[['value']],by=list(sub_dat[,'g1']),FUN=mean),by="Group.1"
    )
    names(mean_sd) <- c('compare','std','mean')
    a <- rbind(a,merge(mean_sd,label,by='compare'))
  }
  names(a) <- c(compare,'std','mean','Letters',group)
  return(a)
}

df2 <- kruskalmc_compare1(dat3,'alpha','group','v')
df2######字母是反着标的(a<b<c)

#导出α-多样性统计学检验结果
write.table(df2,file="alpha_others_kruskaltest.xls",sep="\t",col.names = NA)

#α-多样性出图
#所有α多样性指数统一出图
p2 = ggplot(dat3)+geom_boxplot(aes(x=group,y=v,fill=group), size=0.5, width=0.4)+
  geom_text(data=df2,aes(x=group,y=mean+2.0*std,label=Letters),size=4,color = "black",fontface = "bold",position=position_dodge(0.2))+
  facet_wrap(~alpha,nrow=1,scales = "free_y")+ labs(x='',y='AlphaDiv')+
  labs(x='',y='AlphaDiv')+
  ggprism::theme_prism()+theme(axis.text.x = element_text(angle = 45))
p2

#图形整体调整
mytheme<-theme_bw()+theme(axis.title=element_text(size=14),axis.text = element_text(size = 14),
                          panel.grid.major = element_line(color = "white"),
                          panel.grid.minor = element_line(color = "white"),
                          axis.ticks = element_line(color='black'),
                          axis.line = element_line(colour = "black"),
                          axis.text.x = element_text(size = 14, angle=45,vjust = 0.5,hjust = 0.5,color = "black"),
                          axis.text.y = element_text(size=14,color = "black"),
                          legend.position = "",
                          plot.margin = unit(c(0.5, 0.5, 0.5, 0.5),"cm"))


p3<-p2+mytheme
p3
ggsave("./alpha_others.pdf", p3, width = 150, height = 80, units = "mm")
ggsave("./alpha_others.png", p3, width = 150, height = 80, units = "mm")
ggsave("./alpha_others.tiff", p3, width = 150, height = 80, units = "mm")



#4.α多样性指数第二次挑选后的单独出图
dat4 <- read.delim('shannon_simpson.txt', row.names = 1, sep = '\t') #读入α多样性数据
dat4
dat4 <- gather(dat4,alpha,v,-(OTU.ID:group))#宽数据变长数据
head(dat4, n = 3)


kruskalmc_compare1 <- function(data,group,compare,value){
  library(multcompView)
  library(pgirmess)##多组两辆比较函数用到的包
  library(multcomp)
  a <- data.frame(stringsAsFactors = F)
  type <- unique(data[,group])
  for (i in type)
  {
    g1=compare
    sub_dat <- data[data[,group]==i,]
    names(sub_dat)[names(sub_dat)==compare] <- 'g1'
    names(sub_dat)[names(sub_dat)==value] <- 'value'
    sub_dat$g1 <- factor(sub_dat$g1)
    options(warn = -1)
    
    k <- kruskalmc(value ~ g1, data=sub_dat, probs=0.05)
    dif <- k$dif.com[['difference']]
    names(dif) <- rownames(k$dif.com)
    difL <- multcompLetters(dif)
    label <- data.frame(difL['Letters'], stringsAsFactors = FALSE)
    label$compare = rownames(label)
    label$type <- i
    
    mean_sd <- merge(aggregate(sub_dat[['value']],by=list(sub_dat[,'g1']),FUN=sd),
                     aggregate(sub_dat[['value']],by=list(sub_dat[,'g1']),FUN=mean),by="Group.1"
    )
    names(mean_sd) <- c('compare','std','mean')
    a <- rbind(a,merge(mean_sd,label,by='compare'))
  }
  names(a) <- c(compare,'std','mean','Letters',group)
  return(a)
}

df4 <- kruskalmc_compare1(dat4,'alpha','group','v')
df4######字母是反着标的(a<b<c)

#导出α-多样性统计学检验结果
write.table(df4,file="shannon_simpson_kruskaltest.xls",sep="\t",col.names = NA)

#α-多样性出图
#所有α多样性指数统一出图
p2 = ggplot(dat4)+geom_boxplot(aes(x=group,y=v,fill=group), size=0.5, width=0.4)+
  geom_text(data=df4,aes(x=group,y=mean+2.0*std,label=Letters),size=4,color = "black",fontface = "bold",position=position_dodge(0.2))+
  facet_wrap(~alpha,nrow=1,scales = "free_y")+ labs(x='',y='AlphaDiv')+
  labs(x='',y='AlphaDiv')+
  ggprism::theme_prism()+theme(axis.text.x = element_text(angle = 45))
p2

#图形整体调整
mytheme<-theme_bw()+theme(axis.title=element_text(size=14),axis.text = element_text(size = 14),
                          panel.grid.major = element_line(color = "white"),
                          panel.grid.minor = element_line(color = "white"),
                          axis.ticks = element_line(color='black'),
                          axis.line = element_line(colour = "black"),
                          axis.text.x = element_text(size = 14, angle=90,vjust = 0.5,hjust = 0.5,color = "black"),
                          axis.text.y = element_text(size=14,color = "black"),
                          legend.position = "",
                          plot.margin = unit(c(0.5, 0.5, 0.5, 0.5),"cm"))


p3<-p2+mytheme
p3
ggsave("./shannon_simpson.pdf", p3, width = 80, height = 80, units = "mm")
ggsave("./shannon_simpson.png", p3, width = 80, height = 80, units = "mm")
ggsave("./shannon_simpson.tiff", p3, width = 80, height = 80, units = "mm")

#其他α多样性指数
dat5 <- read.delim('alpha_others2.txt', row.names = 1, sep = '\t') #读入α多样性数据
dat5
dat5 <- gather(dat5,alpha,v,-(OTU.ID:group))#宽数据变长数据
head(dat5, n = 3)


kruskalmc_compare1 <- function(data,group,compare,value){
  library(multcompView)
  library(pgirmess)##多组两辆比较函数用到的包
  library(multcomp)
  a <- data.frame(stringsAsFactors = F)
  type <- unique(data[,group])
  for (i in type)
  {
    g1=compare
    sub_dat <- data[data[,group]==i,]
    names(sub_dat)[names(sub_dat)==compare] <- 'g1'
    names(sub_dat)[names(sub_dat)==value] <- 'value'
    sub_dat$g1 <- factor(sub_dat$g1)
    options(warn = -1)
    
    k <- kruskalmc(value ~ g1, data=sub_dat, probs=0.05)
    dif <- k$dif.com[['difference']]
    names(dif) <- rownames(k$dif.com)
    difL <- multcompLetters(dif)
    label <- data.frame(difL['Letters'], stringsAsFactors = FALSE)
    label$compare = rownames(label)
    label$type <- i
    
    mean_sd <- merge(aggregate(sub_dat[['value']],by=list(sub_dat[,'g1']),FUN=sd),
                     aggregate(sub_dat[['value']],by=list(sub_dat[,'g1']),FUN=mean),by="Group.1"
    )
    names(mean_sd) <- c('compare','std','mean')
    a <- rbind(a,merge(mean_sd,label,by='compare'))
  }
  names(a) <- c(compare,'std','mean','Letters',group)
  return(a)
}

df5 <- kruskalmc_compare1(dat5,'alpha','group','v')
df5######字母是反着标的(a<b<c)

#导出α-多样性统计学检验结果
write.table(df5,file="alpha_others2_kruskaltest.xls",sep="\t",col.names = NA)

#α-多样性出图
#所有α多样性指数统一出图
p2 = ggplot(dat5)+geom_boxplot(aes(x=group,y=v,fill=group), size=0.5, width=0.4)+
  geom_text(data=df5,aes(x=group,y=mean+2.0*std,label=Letters),size=4,color = "black",fontface = "bold",position=position_dodge(0.2))+
  facet_wrap(~alpha,nrow=1,scales = "free_y")+ labs(x='',y='AlphaDiv')+
  labs(x='',y='AlphaDiv')+
  scale_fill_manual(values =c('#CD5B45', '#228B22', '#00688B',"orange"))
  ggprism::theme_prism()+theme(axis.text.x = element_text(angle = 45))
p2

#图形整体调整
mytheme<-theme_bw()+theme(axis.title=element_text(size=14),axis.text = element_text(size = 14),
                          panel.grid.major = element_line(color = "white"),
                          panel.grid.minor = element_line(color = "white"),
                          axis.ticks = element_line(color='black'),
                          axis.line = element_line(colour = "black"),
                          axis.text.x = element_text(size = 14, angle=90,vjust = 0.5,hjust = 0.5,color = "black"),
                          axis.text.y = element_text(size=14,color = "black"),
                          legend.position = "",
                          plot.margin = unit(c(0.5, 0.5, 0.5, 0.5),"cm"))


p3<-p2+mytheme
p3
ggsave("./others2.pdf", p3, width = 180, height = 80, units = "mm")
ggsave("./others2.png", p3, width = 180, height = 80, units = "mm")
ggsave("./others2.tiff", p3, width = 180, height = 80, units = "mm")










####β-diversity####

library("phyloseq")
library("ggplot2")
library("vegan")
library("ggpubr")
#加载必要的包

# 加载github包安装工具
library(devtools)

# 检测amplicon包是否安装，没有从源码安装
if (!requireNamespace("amplicon", quietly=TRUE))
  install_github("microbiota/amplicon")
# 提示升级，选择3 None不升级；升级会容易出现报错
# library加载包，suppress不显示消息和警告信息
suppressWarnings(suppressMessages(library(amplicon)))
library(amplicon)

setwd("E:/博士数据2/小麦微生物组数据/customer_backup/otus/β-diversity/")
#首先导入数据文件

####β多样性距离的计算（跟着生信小白鱼学习）####
##vegan 包 vegdist() 计算群落相似性/距离，详情 ?vegdist
#以下相似度（S）和距离测度（D）的转换公式，均使用 S = 1-D
library(vegan)

#读入物种数据
otu <- read.delim('OTU_Tables.txt', row.names = 1, sep = '\t')
otu <- t(otu)

#jaccard 相异度
jacc_dis <- vegdist(otu, method = 'jaccard', binary = TRUE)
#jaccard 相似度
jacc_sim <- 1 - jacc_dis

#欧几里得距离
eucl_dis <- vegdist(otu, method = 'euclidean')
#这种距离转换为相似性时，需要首先作个标准化，例如
eucl_dis_norm <- eucl_dis/max(eucl_dis)
eucl_sim <- 1 - eucl_dis_norm

#弦距离
otu_chord <- decostand(otu, method = 'normalize')
chord_dis <- vegdist(otu_chord, method = 'euclidean')

#Hellinger 距离
otu_hell <- decostand(otu, method = 'hellinger')
hell_dis <- vegdist(otu_hell, method = 'euclidean')

#卡方距离
#先做卡方标准化，再计算欧几里得距离，即可得
otu_chi <- decostand(otu, method = 'chi.square')
otu_chi_dis <- vegdist(otu_chi, 'euclidean')

#Bray-curtis 距离
bray_dis <- vegdist(otu, method = 'bray')
#Bray-curtis 相似度
bray_sim <- 1 - bray_dis

#更多示例不再展示了
#以上结果均为 dis 类型，若要输出在本地，需要首先转化为 matrix 类型
#以 Bray-curtis 距离为例
bray <- as.matrix(bray_dis)
write.table(bray, 'bray-curtis_distance.txt', col.names = NA, sep = '\t', quote = FALSE)

#距离矩阵热图
heatmap(bray, scale = 'none', RowSideColors = rep(c('red', 'blue'), each = 18), ColSideColors = rep(c('red', 'blue'), each = 18))

##phyloseq 包 distance() 计算群落相似性/距离，详情查看 ?distance
library(phyloseq)

#读入物种数据
otu <- read.delim('OTU_Tables.txt', row.names = 1, sep = '\t')
#读入进化树文件（这里为有根树）
tree <- read_tree('rep_set.fasta')
#不过师兄说无根树也可以，只是每次计算的 UUF 和 WUF 会不一样，他会随机指定根

#首先构建 phyloseq 对象
physeq <- phyloseq(otu_table(otu, taxa_are_rows = TRUE), phy_tree(tree))

#计算加权 Unifrac 距离
wei_unif_dis <- distance(physeq, method = 'wunifrac')
#计算非加权 Unifrac 距离
unwei_unif_dis <- distance(physeq, method = 'unifrac')

#也可用其计算其它类型距离，如 Bray-curtis 距离
bray_dis <- distance(physeq, method = 'bray')

#更多示例不再展示了
#注：vegan 包 vegdist() 中的 method 参数同样适用 phyloseq 包 distance() 中的 method 参数

#以上结果均为 dis 类型，若要输出在本地，需要首先转化为 matrix 或类型
#以加权 Unifrac 距离为例
wei_unif <- as.matrix(wei_unif_dis)
write.table(wei_unif, 'weighted-unifrac_distance.txt', col.names = NA, sep = '\t', quote = FALSE)

##GUniFrac 包 GUniFrac() 计算微生物群落的 Unifrac 距离，详情查看 ?GUniFrac
library(GUniFrac)

#读入物种数据
otu <- read.delim('otu_table.txt', row.names = 1, sep = '\t')
otu <- data.frame(t(otu))
#读入进化树文件（这个只能使用有根树）
tree <- read.tree('otu_tree.tre')

#计算加权 Unifrac 距离
unifrac <- GUniFrac(otu, tree)$unifracs
wei_unif_dis <- as.dist(unifrac[, , 'd_1'])	
#计算非加权 Unifrac 距离
unwei_unif_dis <- as.dist(unifrac[, , 'd_UW'])

#以上结果均为 dis 类型，若要输出在本地，需要首先转化为 matrix 或类型
#以非加权 Unifrac 距离为例
unwei_unif <- as.matrix(unwei_unif_dis)
write.table(unwei_unif, 'unweighted-unifrac_distance.txt', col.names = NA, sep = '\t', quote = FALSE)





####β多样性排序分析####
##非限制性排序PCA/PCoA/NMDS
####方法一：NMDS分析####
#方法一：使用R语言vegan包进行NMDS分析，输入文件为样本的OTU丰度表格。
#载入分析包
library(vegan)
#载入分析数据
otu <- read.delim('./OTU_Tables.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
otu
otu <- t(otu)
otu
#计算与 Beta 多样性有关的群落相异指数，例如使用 vegan 包计算 Bray-curtis 距离，详情加载 vegan 包后 ?vegdist
dis <- vegan::vegdist(otu, method = 'bray')

#以矩阵形式输出
dis <- as.matrix(dis)
write.table(dis, 'Bray-curtis.txt', sep = '\t', col.names = NA, quote = FALSE)


#对分析数据进行转置、去除缺失值
#otu <- t(otu)
#otu[is.na(otu)] <- 0
#otu <- data.frame(t(otu))
#otu
#进行NMDS分析
nmds <- metaMDS(otu)
#保存stress结果
capture.output(nmds,file = "Stress.txt")
#保存scores结果
nmds_scores=scores(nmds, choices = c(1,2))
write.table(nmds_scores,file="NMDS_scores.txt")


#结果图绘制
#与PCoA一样，结果图的绘制还需要一个分组文件。
#首先制作用于绘图的数据文件。

scores <-  read.table("NMDS_scores.txt",header=T)
scores
groups <- read.table("mapping.txt", head=T)
groups
show_point_shapes() + theme_classic()
pich=c(21:24)
Palette <- c("#CC6666","#999999","#B2182B","#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7","#9999CC","#66CC99","#ADD1E5")
#cbbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442","#808000")
Palette1 <- c("#000000","#000000","#000000","#000000","#000000","#000000","#000000","#000000","#000000","#000000","#000000","#000000")
plotdata <- data.frame(rownames(scores),scores$NMDS1,scores$NMDS2,groups$Group)
colnames(plotdata)=c("sample","MDS1","MDS2","Group")
plotdata$sample <- factor(plotdata$sample)
plotdata$MDS1 <- as.numeric(as.vector(plotdata$MDS1))
plotdata$MDS2 <- as.numeric(as.vector(plotdata$MDS2))

library(ggplot2)
#mi = c("#FFF5EB" ,"#FEE6CE" ,"#FDD0A2", "#FDAE6B", "#FD8D3C", "#F16913", "#D94801", "#A63603", "#7F2704","black")
p<-ggplot(plotdata, aes(MDS1, MDS2)) +
  geom_point(aes(colour=Group,shape=Group,fill=Group),alpha=.7, size=4)+
  scale_shape_manual(values=pich)+
  scale_colour_manual(values=Palette1)+
  scale_fill_manual(values=Palette)+
  labs(title="") + xlab(paste("NMDS1")) + ylab(paste("NMDS2"))+
  theme(text=element_text(size=10))+
  #scale_x_continuous(breaks = seq(-2,2,0.5)) +
  geom_hline(aes(yintercept=0), colour="black", linetype=1) +
  geom_vline(aes(xintercept=max(plotdata$MDS2/2)), colour="black", linetype="dashed")+
  #geom_vline(aes(xintercept = 0),linetype="dotted")+
  #geom_hline(aes(yintercept = 0),linetype="dotted")+
  geom_text(aes(x=max(MDS1),y=max(MDS2)),hjust=3,vjust=1,size=5,label=paste("Stress = ",round(nmds$stress,3),sep = ""),colour="black")+
  #hjust为正数，向左移动，数字越大，越靠右；vjust为正数向下移动，数值越大，越靠下。
  theme(panel.background = element_rect(fill='white', colour='black'), panel.grid=element_blank(), 
        axis.title = element_text(color='black',size=18),axis.ticks.length = unit(0.4,"lines"), 
        axis.ticks = element_line(color='black'),axis.line = element_line(colour = "black"), 
        axis.title.x=element_text(colour='black', size=18),axis.title.y=element_text(colour='black', size=18),axis.text=element_text(colour='black',size=18),
        legend.title=element_blank(),legend.text=element_text(size=18),legend.key=element_blank())+
  theme(plot.title = element_text(size=20,colour = "black",face = "bold"))
p

#使用vegan包中的adonis函数进行PERMANOVA分析。
otu.adonis=adonis(otu~Group,data = groups,permutations = 10000,distance = "bray")
otu.adonis
#p1<-p+geom_text(aes(x = 0.10,y = 0.05,label = paste('PERMANOVA:\n  Rseq = 0.880\n  p= 0.001',sep = "")),hjust=2,vjust=1,size=5)
#p1

ggsave(plot = p,"nmds_of_bary_curtis.pdf",dpi = 300,width = 7,height = 4)
ggsave(plot = p,"nmds_of_bary_curtis.png",dpi = 300,width = 7,height = 4)
ggsave(plot = p,"nmds_of_bary_curtis.tiff",dpi = 300,width = 7,height = 4)
#保存图片


####方法二：PCoA分析####
#方法二：使用R语言vegan包进行PCoA分析，输入文件为样本的OTU丰度表格。
#载入分析包
library(vegan)
#载入分析数据
otu <- read.delim('./OTU_Tables.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
otu

#对分析数据进行转置、去除缺失值
#otu <- t(otu)
otu[is.na(otu)] <- 0
otu <- data.frame(t(otu))

#计算Bray-Curtis距离
otu <- vegdist(otu, method = "bray")
#进行PCoA分析
pcoa<- pcoa(otu, correction = "none", rn = NULL)

#载入分组文件
groups <- read.table("mapping.txt", head=T,sep = "\t",colClasses = c("character"))
#将分组文件转换为列表
groups <- as.list(groups)
#制作绘图文件
PCOA1 = pcoa$vectors[,1]
PCOA2 = pcoa$vectors[,2]
plotdata2 <- data.frame(rownames(pcoa$vectors),PCOA1,PCOA2,groups$Group)
colnames(plotdata2) <-c("sample","PCOA1","PCOA2","Group")

#用于填充样本点的颜色
Palette <- c("#CC6666","#999999","#B2182B","#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7","#9999CC","#66CC99","#ADD1E5")
#样本点的边框颜色
Palette1 <- c("#000000","#000000","#000000","#000000","#000000","#000000","#000000","#000000","#000000","#000000","#000000","#000000")

#用于绘制横纵坐标label的文本，以显示解释比例
PCOA1 <-floor(pcoa$values$Relative_eig[1]*100)
PCOA2 <-floor(pcoa$values$Relative_eig[2]*100)

#载入绘图包
library(ggplot2)
#匹配横纵坐标
pp<-ggplot(plotdata2, aes(PCOA1, PCOA2)) +
  #绘制样本点，根据分组匹配颜色和形状，size调整点的大小
  geom_point(aes(colour=Group,shape=Group,fill=Group),size=5)+
  #匹配形状、边框和填充的图例
  scale_shape_manual(values=pich)+
  scale_colour_manual(values=Palette)+
  scale_fill_manual(values=Palette)+
  #设置标题和横纵坐标label文字
  labs(title="") + 
  xlab(paste("PCOA1 ( ",PCOA1,"%"," )",sep="")) + 
  ylab(paste("PCOA2 ( ",PCOA2,"%"," )",sep=""))+
  theme(text=element_text(size=18))+
  #添加横纵两条虚线
  scale_y_continuous(breaks = seq(-0.2,0.5,0.1)) +
  geom_vline(aes(xintercept = 0),linetype="dotted")+
  geom_hline(aes(yintercept = 0),linetype="dotted")+
  #geom_vline(aes(xintercept=0), colour="black", linetype=1) +
  #geom_hline(aes(yintercept=max(plotdata$MDS2/2)), colour="black", linetype=1)+
  #调整背景、坐标轴、图例的格式
  theme(panel.background = element_rect(fill='white', colour='black'),
        panel.grid=element_blank(), 
        axis.title = element_text(color='black',size=18),
        axis.ticks.length = unit(0.4,"lines"), axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"), 
        axis.title.x=element_text(colour='black', size=18),
        axis.title.y=element_text(colour='black', size=18),
        axis.text=element_text(colour='black',size=18),
        legend.title=element_blank(),
        legend.text=element_text(size=15),
        legend.key=element_blank(),
        #legend.background = element_rect(colour = "black"),
        legend.key.height=unit(0.8,"cm"))+
  #legend.title=element_blank(),legend.text=element_text(size=18),legend.key=element_blank())+
  #设置标题的格式
  theme(plot.title = element_text(size=20,colour = "black",hjust = 0.5,face = "bold"))

pp

#使用vegan包中的adonis函数进行PERMANOVA分析。
otu.adonis=adonis(otu~Group,data = groups,permutations = 1000,distance = "bray")
pcoa1<-p+geom_text(aes(label = paste(x=0,y=1.5,'Rseq = 0.880 \n p= 0.001',sep = "")),hjust=2,vjust=1,size=5)
pcoa1

ggsave(plot = pp,"PCoA_of_bary_curtis.pdf",dpi = 300,width = 7,height = 6)
ggsave(plot = pp,"PCoA_of_bary_curtis.png",dpi = 300,width = 7,height = 6)
ggsave(plot = pp,"PCoA_of_bary_curtis.tiff",dpi = 300,width = 7,height = 6)
#保存图片


####限制性排序####
####cpcoa####
#方法一：限制性主坐标排序分析
# 读取元数据
metadata=read.table("mapping.txt", header=T, row.names=1, sep="\t", comment.char="", stringsAsFactors=F)
metadata
# 设置实验分组列名
group="Group"
group
# 图片16：10，即半版89 mm x 56 mm; 183 mm x 114 mm
w=89
h=59

##限制性PCoA(Constrained PCoA)，要求至少分组数量 >= 3时才可用，否则请跳过此步分析。

#方法1. 基于特征OTU/ASV表使用vegan包计算矩阵矩阵并进行限制性主坐标分析。

library(amplicon)
# 读取特征ASV表，多样性分析输入抽平标准化的表
otutab=read.table("OTU_Tables.txt", header=T, row.names=1, sep="\t", comment.char="")
# 基于特征表、原数据、距离类型、分组列名、是否添加置信椭圆，是否添加样本标签
(p = beta_cpcoa(otutab, metadata, dis="bray", groupID="Group", ellipse=T, label=F))
ggsave("p1.cpcoa_otutab.pdf", p, width=w, height=h, units="mm")
ggsave("p1.cpcoa_otutab.png", p, width=w, height=h, units="mm", dpi=300)

##基于weighted_unifrac距离
distance_mat <- read.table("weighted_unifrac.txt", header=T, row.names=1, sep="\t", comment.char="")
distance_mat[1:3, 1:3]
(p=beta_cpcoa_dis(distance_mat, metadata, groupID=group))
ggsave("p2.cpcoa_weighted_unifrac_distance.pdf", p, width=w, height=h, units="mm")
ggsave("p2.cpcoa_weighted_unifrac_distance.png", p, width=w, height=h, units="mm", dpi=300)


##基于unweighted_unifrac距离
distance_mat2 <- read.table("unweighted_unifrac.txt", header=T, row.names=1, sep="\t", comment.char="")
distance_mat2[1:3, 1:3]
(p=beta_cpcoa_dis(distance_mat2, metadata, groupID=group))
ggsave("p2.cpcoa_unweighted_unifrac_distance.pdf", p, width=w, height=h, units="mm")
ggsave("p2.cpcoa_unweighted_unifrac_distance.png", p, width=w, height=h, units="mm", dpi=300)

##基于binary_jaccard距离
distance_mat3 <- read.table("binary_jaccard.txt", header=T, row.names=1, sep="\t", comment.char="")
distance_mat3[1:3, 1:3]
(p=beta_cpcoa_dis(distance_mat3, metadata, groupID=group))
ggsave("p2.cpcoa_binary_jaccard_distance.pdf", p, width=w, height=h, units="mm")
ggsave("p2.cpcoa_binary_jaccard_distance.png", p, width=w, height=h, units="mm", dpi=300)


distance_mat4 <- read.table("bray_curtis2.txt", header=T, row.names=1, sep="\t", comment.char="")
distance_mat4[1:3, 1:3]
(p=beta_cpcoa_dis(distance_mat4, metadata, groupID=group))



##通过相似性或相异指数的数值分布比较群落β多样性高低
#使用距离矩阵进行计算和作图

#读取 Bray-curtis 距离矩阵
dis <- read.delim('Bray-curtis.txt', row.names = 1)

#读取样本分组信息
group <- read.delim('mapping.txt', stringsAsFactors = FALSE)

##例如，比较 Ggt、H、SA、S_135 四组之间，群落的 Beta 多样性差异
#根据分组获得组内距离矩阵
Ggt <- subset(group, Group == 'Ggt')$SampleID
dis_Ggt <- dis[Ggt,Ggt]

H <- subset(group, Group == 'H')$SampleID
dis_H <- dis[H,H]

SA <- subset(group, Group == 'SA')$SampleID
dis_SA <- dis[SA,SA]

S_135 <- subset(group, Group == 'S_135')$SampleID
dis_S_135 <- dis[S_135,S_135]


#将矩阵转化为向量，以便用于作图和统计
dis_Ggt <- as.vector(as.dist(dis_Ggt))
dis_Ggt
dis_H <- as.vector(as.dist(dis_H))
dis_H
dis_SA <- as.vector(as.dist(dis_SA))
dis_SA
dis_S_135 <- as.vector(as.dist(dis_S_135))
dis_S_135

#构建作图数据集
dat2 <- data.frame(
  dis = c(dis_Ggt, dis_H, dis_SA, dis_S_135),
  group = factor(c(
    rep('Ggt', length(dis_Ggt)), 
    rep('H', length(dis_H)), 
    rep('SA', length(dis_SA)),
    rep('S_135', length(dis_S_135))
  ), levels = c('Ggt', 'H', 'SA', 'S_135'))
)

#使用 ggplot2 绘制各组内 Bray-curtis 距离指数分布的箱线图
library(ggplot2)

p <- ggplot(dat2, aes(group, dis)) +
  geom_boxplot(aes(fill = group), width = 0.3) +
  scale_fill_manual(values =c('#CD5B45', '#228B22', '#00688B',"orange")) +
  theme(panel.grid = element_blank(), panel.background = element_blank(), 
        axis.line = element_line(colour = 'black'), legend.position = 'none') +
  labs(x = NULL, y = 'Bray-Curtis dissimilarity\n')

p

#三组的整体差异分析，使用 Kruskal-Wallis Test 执行，详情 ?kruskal.test
kruskal.test(dis~group, data = dat2)

#如果整体显著再进行两两分组的比较，使用 Wilcoxon 秩和检验执行双侧检验，详情 ?wilcox.test
wilcox.test(dis_Ggt, dis_H, alternative = 'two.sided')
wilcox.test(dis_Ggt, dis_SA, alternative = 'two.sided')
wilcox.test(dis_Ggt, dis_S_135, alternative = 'two.sided')
wilcox.test(dis_H, dis_SA, alternative = 'two.sided')
wilcox.test(dis_H, dis_S_135, alternative = 'two.sided')
wilcox.test(dis_SA, dis_S_135, alternative = 'two.sided')

#考虑到 Wilcoxon 秩和检验体现了中位数的差异，因此计算三组数据的中位数以评估 Beta 多样性的高低水平
median(dis_Ggt)
median(dis_H)
median(dis_SA)
median(dis_S_135)

#基于上述统计结果，判断好组间差异后，将差异分析结果添加到箱线图中
p +
  annotate('text', label = 'Kruskal-Wallis Test', x = 1, y = 0.56, size = 3) +
  annotate('text', label = sprintf('italic(P) = %.3f', 0.05991), x = 1, y = 0.53, size = 3, parse = TRUE)

ggsave("Bray-curtis dissimilarity.pdf", width = 6, height = 7, units="cm")
ggsave("Bray-curtis dissimilarity.tiff", width = 6, height = 7, units = "cm")
ggsave("Bray-curtis dissimilarity.png", width = 6, height = 7, units = "cm")


#p +
annotate('text', label = 'Kruskal-Wallis Test', x = 1, y = 0.56, size = 3) +
  annotate('text', label = sprintf('italic(P) = %.3f', 0.05991), x = 1, y = 0.53, size = 3, parse = TRUE) +
  annotate('text', label = 'c', x = 1, y = max(dis_env1)+0.05, size = 3) +
  annotate('text', label = 'a', x = 2, y = max(dis_env2)+0.05, size = 3) +
  annotate('text', label = 'b', x = 3, y = max(dis_env3)+0.05, size = 3) +
  annotate('text', label = 'b', x = 3, y = max(dis_env3)+0.05, size = 3)


##通过置换多元离散度（PERMDISP）分析群落β多样性差异
##示例
library(vegan)
#24行样本和44列物种丰度的示例数据集，详情加载vegan包后？varespec
data(varespec)
varespec[1:6,1:6]

#24个样本的分组信息，放牧（grazed）和非放牧（ungrazed）
groups <- factor(c(rep(1,16),rep(2,8)),labels=c("grazed","ungrazed"))
groups

#首先根据丰度表计算群落距离，这里以常见的bray-curtis相异指数为例
dis <-vegdist(varespec,method = "bray")
dis
#分析两组群落距离的多元离散度（MDISP），type中指定计算质心的方式，详情？betadisper
mod <- betadisper(d=dis,group=groups,type="centroid")
mod

##继续检验两组群落距离的多元离散度（MDISP）的P值
#参数方法（数据满足多元正态性时），F检验的原理
#由于群落数据通常拒绝正态性，因此群落分析中不是很常用它
anova(mod)

#非参数的方法（数据拒绝多元正态性），基于伪F值的置换检验原理，以下通过置换999次获取P值的估计
#群落数据分析中经常用到置换检验的思维，这个方法就是本篇提到的置换多元离散度（PERMDISP）分析
set(123)
permutest(mod,pairwise = TRUE,permutations = 999)

#配合PCOA观测群落beta多样性的差异
plot(mod,ellipse=TRUE,hull= FALSE,conf=0.95)

#群落组内相异度的箱线图
boxplot(mod)



##实战
rm(list=ls())
library(vegan)

#载入分析数据
otu <- read.delim('./OTU_Tables.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
otu
otu <- t(otu)
otu[1:12,1:6]
#计算与 Beta 多样性有关的群落相异指数，例如使用 vegan 包计算 Bray-curtis 距离，详情加载 vegan 包后 ?vegdist
dis <- vegan::vegdist(otu, method = 'bray')
dis

#载入分组文件
groups <- factor(c(rep(1,3),rep(2,3),rep(3,3),rep(4,3)),labels=c("Ggt","S_135","SA","H"))
groups


#分析两组群落距离的多元离散度（MDISP），type中指定计算质心的方式，详情？betadisper
mod <- betadisper(d=dis,group=groups,type="centroid")
mod

##继续检验两组群落距离的多元离散度（MDISP）的P值
#参数方法（数据满足多元正态性时），F检验的原理
#由于群落数据通常拒绝正态性，因此群落分析中不是很常用它
anova(mod)

#非参数的方法（数据拒绝多元正态性），基于伪F值的置换检验原理，以下通过置换999次获取P值的估计
#群落数据分析中经常用到置换检验的思维，这个方法就是本篇提到的置换多元离散度（PERMDISP）分析
set(123)
permutest(mod,pairwise = TRUE,permutations = 999)

#配合PCOA观测群落beta多样性的差异
plot(mod,ellipse=TRUE,hull= FALSE,conf=0.95)

#群落组内相异度的箱线图
pdf(paste0("beta_dispersion-bray_curtis.pdf"),width=3.5,height=6)

boxplot(mod,col = c("#CC6666","#9999CC","#66CC99","#999999"))

dev.off()

####重大发现：上述计算过程中我们并不知道箱线图中四个数值的具体数值，但是经过下面的代码，就可以呈现出具体数值
p <- boxplot(mod,col = c("#CC6666","#9999CC","#66CC99","#999999"))
p

boxplot(mod)




dis2 <- read.table("bray_curtis2.txt", header=T, row.names=1, sep="\t", comment.char="")
dis2

groups2 <- read.table("mapping2.txt", header=T, row.names=1, sep="\t", comment.char="")
groups2

#分析两组群落距离的多元离散度（MDISP），type中指定计算质心的方式，详情？betadisper
mod2 <- betadisper(d=dis2,group=groups2$Group,type="centroid")
mod2





####物种组成分析####
#堆叠柱状图、堆叠状图中间连线、聚类+堆叠柱状图、弦图、气泡图/树状图、热图

####堆叠柱状图中间连线####

setwd("E:/X-14952B/otus/microbial composition/")
#首先导入数据文件

##带物种注释的相对丰度表的准备，如果有则不需要这一步
#如果是count的绝对丰度表格，需要转换为相对丰度表格，转换方法是使得某一个样本中所有otu的和为1，每个otu在该样本中的百分比即为相对丰度  
#将绝对丰度转化为百分比形式的相对丰度
phylum_per <- as.data.frame(lapply(OTUID, function(x) x / sum(x)))
row.names(OTUID_per) <- row.names(OTUID) 
#加一下行名
#不带物种注释的OTU相对丰度表格
otu <- read.delim('./OTU_Tables_relative.txt', row.names = 1)
otu
#物种注释表格
tax <- read.delim('./taxonomy2.txt',row.names = 1)
head(tax, n = 3)
#将
dat <- merge(x=otu,y=tax,by='row.names')
head(dat, n = 3)
dat =dplyr::rename(dat,OTUID = Row.names)
head(dat, n = 3)
write.table(dat, 'otu_AAT_tax.txt', row.names = FALSE, sep = '\t', quote = FALSE)



##直接读入带物种注释的OTU相对丰度表格
dat <- read.delim('OTU_Tables_relative.txt', row.names = 1)
dat

##按分类水平分组汇总(根据自己需求更改要展示的物种水平)
phylum <-aggregate(dat[,1:12],by=list(dat$Phylum),FUN=sum)
phylum
write.table(phylum, 'phylum.txt', row.names = FALSE, sep = '\t', quote = FALSE)

class <-aggregate(dat[,1:12],by=list(dat$Class),FUN=sum)
class
write.table(class, 'class.txt', row.names = FALSE, sep = '\t', quote = FALSE)

order <-aggregate(dat[,1:12],by=list(dat$Order),FUN=sum)
order
write.table(order, 'order.txt', row.names = FALSE, sep = '\t', quote = FALSE)


family <-aggregate(dat[,1:12],by=list(dat$Family),FUN=sum)
family
write.table(family, 'family.txt', row.names = FALSE, sep = '\t', quote = FALSE)


genus <-aggregate(dat[,1:12],by=list(dat$Genus),FUN=sum)
genus
write.table(genus, 'genus.txt', row.names = FALSE, sep = '\t', quote = FALSE)

species <-aggregate(dat[,1:12],by=list(dat$Species),FUN=sum)
species
write.table(species, 'species.txt', row.names = FALSE, sep = '\t', quote = FALSE)


####门水平top10####
#色卡，建议收藏留用，

Palette <- c("#B2182B","#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7","#CC6666","#9999CC","#66CC99","#999999","#ADD1E5")

#数据表格读取
phylum <- read.table ("phylum.txt", header = TRUE, sep ="\t")
phylum
#首先倒入绘图数据，因为是用门水平做例子，因此导入的是L2水平的结果文件，L2-L6分别代表门、纲、目、科和属。
#as.data.frame(),as.array()等函数，可以对已有的常用的类型进行转换。as.data.frame()检查对象是否是数据帧，或者是否可以强制执行
phylum <- as.data.frame(phylum)
?as.data.frame
#删除OTU相对丰度表中，有关p_的字样
taxon <- gsub(".*(p__)","",phylum$Taxon)
taxon
#gsub()可以用于字段的删减、增补、替换和切割，可以处理一个字段也可以处理由字段组成的向量。
?gsub
#删除Taxon中有关Bacteria;的字样
#taxon <- gsub(".*(Bacteria;)","",taxon)
#删除taxon中的Other字样
#taxon <- gsub(".*(Other)","",taxon)

taxon <- gsub("(\\[)","",taxon)
taxon <- gsub("(\\])","",taxon)
#将新产生的表头与原来的表头进行替换
phylum$Taxon <- taxon
phylum
phylum[phylum==""] <- NA
phylum <- na.omit(phylum)
phylum <- t(phylum)
colnames(phylum) <- phylum[1,]
phylum <- phylum[-1,]
phylum
#使用gsub从分类学名称一列中将所需要的物种名提取出来，??????第二行命令，门水平的数据为p，其它分类学水平依次为c、o、f和g。

#接下来需要对物种根据其在所有样本中的平均丰度进行排序，保留丰度排在前10的物种，对其它物种的丰度进行求和，将其命名为Others。


phylum <- as.matrix(phylum)
phylum1 <- matrix(as.numeric(phylum),nrow = nrow(phylum))
rownames(phylum1) <- rownames(phylum)
colnames(phylum1) <- colnames(phylum)
phylum1 <- t(phylum1)
sum <- apply(phylum1,1,sum) 
phylum1 <- cbind(phylum1,sum)
phylum1 <- as.data.frame(phylum1)
phylum1 <- phylum1[order(phylum1[,"sum"],decreasing = T),]
write.table(phylum1, file = "Phylumorder.txt",sep = '\t',)

phylum1 <- subset(phylum1, select = -sum)
phylum1 <- t(phylum1)
phylum1 <- subset(phylum1, select = -Unassigned)
phylum2 <- t(phylum1)
phylum2
#如果需要增减物种的数目，请调整此处的数值
phylum3 <- phylum2[1:10,]
phylum3
phylum3 <- t(phylum3)
sum <- apply(phylum3,1,sum) 
Others <- 1-sum
phylum10 <- cbind(phylum3,Others)
phylum10
write.table(phylum10, file = "PhylumTop10.txt",sep = '\t',)
#接下来进行新的样本名称的替换，并保证图像中样本的绘制顺序与预期一致。

sample.name <- read.table("mapping.txt",head=T,colClasses=c("character","character"),sep = "\t")
sample.name
phylum10 <- as.data.frame(phylum10)
phylum10$SampleID <- rownames(phylum10)
phylum10 <- merge(sample.name,phylum10,by = "SampleID",sort = F)
rownames(phylum10) <- phylum10$SampleID
phylum10 <- subset(phylum10, select = -Group)
phylum10 <- subset(phylum10, select = -SampleID)
#phylum10 <- subset(phylum10, select = -V3)
phylum10 <- t(phylum10)
phylum10
#图像绘制

#绘制的是冲积条形图所使用的绘图包是“ggplot2+ggalluvial”。
library(reshape2)
taxon <- melt(phylum10)
taxon
colnames(taxon) <- c("Taxon","variable","value")
library(ggplot2)
#install.packages("ggalluvial")
library(ggalluvial)
library("ggsci")
library("ggplot2")
library("gridExtra")

####堆叠柱状图####

p1 <- ggplot(data = taxon,aes(variable, value, fill = Taxon, color = Taxon))+
  geom_bar(stat = "identity", position = "stack")
#+
#  scale_fill_brewer(palette = "RdYlGn")+
#  scale_color_brewer(palette = "RdYlGn")
p1



####冲积条形图####

p <- ggplot(data = taxon,aes(x = variable, y = value, alluvium = Taxon, stratum = Taxon))
p
p1 <- p + geom_alluvium(aes(fill = Taxon),alpha = .5,width = 0.6) + geom_stratum(aes(fill = Taxon),width = 0.6)
p1

#修改横纵坐标轴标题。

p2 <- p1 + ylab(label = "Relative Abundance") + xlab(label = "")
p2
#修改图像配色。
p3 <- p2 + scale_fill_manual(values = Palette)
#p3 <- p2 + scale_fill_brewer(palette = "RdYlGn")
p3
#调整绘图区主题。

p4 <- p3 + theme_bw()+ theme(panel.grid=element_blank()) + theme(panel.border = element_blank()) + theme(panel.background=element_rect(fill='transparent',color='black'),plot.margin = unit(c(3,5,1,1),"mm"))

#调整坐标轴文字。

p5 <- p4 + theme(axis.text.x=element_text(colour="black",size=12,angle = 270,vjust = 0.5,hjust = 0)) + 
  theme(axis.text.y=element_text(colour = "black",size = 12)) + 
  theme(axis.line = element_line(colour = "black"))+ 
  theme(axis.title.y = element_text(size = 14,
                                    margin = unit(c(0,1,0,1),"lines"))) 
p5

#去除横纵坐标与边框之间的空白。

p6 <- p5 + scale_y_continuous(limits = c(0,1),expand = c(0,0))
p6
#调整图例。

p7 <- p6 + theme(legend.text = element_text(color="black",size = 12)) + theme(legend.title = element_text(size = 18,colour = "black"))
p7
#修改字体。

p8 <- p7 + theme(text = element_text(family = "Times"))
p8


#输出图像。
ggsave("phylumTop104.pdf", width = 18, height = 10, units="cm")
ggsave("phylumTop104.tiff", width = 18, height = 10, units = "cm")
ggsave("phylumTop104.png", width = 18, height = 10, units = "cm")



##纲水平top10
#色卡，建议收藏留用，

Palette <- c("#B2182B","#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7","#CC6666","#9999CC","#66CC99","#999999","#ADD1E5")

#数据表格读取
class <- read.table ("class.txt", header = TRUE, sep ="\t")
class
#首先倒入绘图数据，因为是用门水平做例子，因此导入的是L2水平的结果文件，L2-L6分别代表门、纲、目、科和属。
#as.data.frame(),as.array()等函数，可以对已有的常用的类型进行转换。as.data.frame()检查对象是否是数据帧，或者是否可以强制执行
class <- as.data.frame(class)
?as.data.frame
#删除OTU相对丰度表中，有关p_的字样
taxon <- gsub(".*(c__)","",class$Taxon)
#gsub()可以用于字段的删减、增补、替换和切割，可以处理一个字段也可以处理由字段组成的向量。
?gsub
#删除Taxon中有关Bacteria;的字样
#taxon <- gsub(".*(Bacteria;)","",taxon)
#删除taxon中的Other字样
#taxon <- gsub(".*(Other)","",taxon)

taxon <- gsub("(\\[)","",taxon)
taxon <- gsub("(\\])","",taxon)
#将新产生的表头与原来的表头进行替换
class$Taxon <- taxon
class
class[class==""] <- NA
class <- na.omit(class)
class <- t(class)
colnames(class) <- class[1,]
class <- class[-1,]

#使用gsub从分类学名称一列中将所需要的物种名提取出来，??????第二行命令，门水平的数据为p，其它分类学水平依次为c、o、f和g。

#接下来需要对物种根据其在所有样本中的平均丰度进行排序，保留丰度排在前10的物种，对其它物种的丰度进行求和，将其命名为Others。


class <- as.matrix(class)
class1 <- matrix(as.numeric(class),nrow = nrow(class))
rownames(class1) <- rownames(class)
colnames(class1) <- colnames(class)
class1 <- t(class1)
sum <- apply(class1,1,sum) 
class1 <- cbind(class1,sum)
class1 <- as.data.frame(class1)
class1 <- class1[order(class1[,"sum"],decreasing = T),]
write.table(class1, file = "Classorder.txt",sep = '\t',)

class1 <- subset(class1, select = -sum)
class1 <- t(class1)
class2 <- subset(class1, select = -Unassigned)
class2 <- t(class2)
#如果需要增减物种的数目，请调整此处的数值
class2 <- class2[1:10,]
class2 <- t(class2)
sum <- apply(class2,1,sum) 
Others <- 1-sum
class10 <- cbind(class2,Others)
class10
write.table(class10, file = "ClassTop10.txt",sep = '\t',)
#接下来进行新的样本名称的替换，并保证图像中样本的绘制顺序与预期一致。

sample.name <- read.table("mapping.txt",head=T,colClasses=c("character","character"),sep = "\t")
sample.name
class10 <- as.data.frame(class10)
class10$SampleID <- rownames(class10)
class10 <- merge(sample.name,class10,by = "SampleID",sort = F)
rownames(class10) <- class10$SampleID
class10 <- subset(class10, select = -Group)
class10 <- subset(class10, select = -SampleID)
#f.abundance10 <- subset(f.abundance10, select = -V3)
class10 <- t(class10)
class10
#图像绘制

#绘制的是冲积条形图所使用的绘图包是“ggplot2+ggalluvial”。
library(reshape2)
taxon <- melt(class10)
taxon
colnames(taxon) <- c("Taxon","variable","value")
library(ggplot2)
#install.packages("ggalluvial")

library(ggalluvial)
#进行绘图坐标的匹配并生成基本的冲积条形图。

p <- ggplot(data = taxon,aes(x = variable, y = value, alluvium = Taxon, stratum = Taxon))
p
p1 <- p + geom_alluvium(aes(fill = Taxon),alpha = .5,width = 0.6) + geom_stratum(aes(fill = Taxon),width = 0.6)
p1

#修改横纵坐标轴标题。

p2 <- p1 + ylab(label = "Relative Abundance") + xlab(label = "")
p2
#修改图像配色。

p3 <- p2 + scale_fill_manual(values = Palette)
p3
#调整绘图区主题。

p4 <- p3 + theme_bw()+ theme(panel.grid=element_blank()) + theme(panel.border = element_blank()) + theme(panel.background=element_rect(fill='transparent',color='black'),plot.margin = unit(c(3,5,1,1),"mm"))

#调整坐标轴文字。

p5 <- p4 + theme(axis.text.x=element_text(colour="black",size=12,angle = 270,vjust = 0.5,hjust = 0)) + 
  theme(axis.text.y=element_text(colour = "black",size = 12)) + 
  theme(axis.line = element_line(colour = "black"))+ 
  theme(axis.title.y = element_text(size = 14,
                                    margin = unit(c(0,1,0,1),"lines"))) 
p5

#去除横纵坐标与边框之间的空白。

p6 <- p5 + scale_y_continuous(limits = c(0,1),expand = c(0,0))
p6
#调整图例。

p7 <- p6 + theme(legend.text = element_text(colour = "black",size = 12)) + theme(legend.title = element_text(size = 18,colour = "black"))
p7
#修改字体。

p8 <- p7 + theme(text = element_text(family = "Times"))
p8

#输出图像。
ggsave("classTop10.pdf", width = 18, height = 10, units="cm")
ggsave("classTop10.tiff", width = 18, height = 10, units = "cm")
ggsave("classTop10.png", width = 18, height = 10, units = "cm")



##目水平top10
#色卡，建议收藏留用，

Palette <- c("#B2182B","#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7","#CC6666","#9999CC","#66CC99","#999999","#ADD1E5")

#数据表格读取
order <- read.table ("order.txt", header = TRUE, sep ="\t")
order
#首先倒入绘图数据，因为是用门水平做例子，因此导入的是L2水平的结果文件，L2-L6分别代表门、纲、目、科和属。
#as.data.frame(),as.array()等函数，可以对已有的常用的类型进行转换。as.data.frame()检查对象是否是数据帧，或者是否可以强制执行
order <- as.data.frame(order)
?as.data.frame
#删除OTU相对丰度表中，有关p_的字样
taxon <- gsub(".*(o__)","",order$Taxon)
#gsub()可以用于字段的删减、增补、替换和切割，可以处理一个字段也可以处理由字段组成的向量。
?gsub
#删除Taxon中有关Bacteria;的字样
#taxon <- gsub(".*(Bacteria;)","",taxon)
#删除taxon中的Other字样
#taxon <- gsub(".*(Other)","",taxon)

taxon <- gsub("(\\[)","",taxon)
taxon <- gsub("(\\])","",taxon)
#将新产生的表头与原来的表头进行替换
order$Taxon <- taxon
order
order[order==""] <- NA
order <- na.omit(order)
order <- t(order)
colnames(order) <- order[1,]
order <- order[-1,]

#使用gsub从分类学名称一列中将所需要的物种名提取出来，??????第二行命令，门水平的数据为p，其它分类学水平依次为c、o、f和g。

#接下来需要对物种根据其在所有样本中的平均丰度进行排序，保留丰度排在前10的物种，对其它物种的丰度进行求和，将其命名为Others。


order <- as.matrix(order)
order1 <- matrix(as.numeric(order),nrow = nrow(order))
rownames(order1) <- rownames(order)
colnames(order1) <- colnames(order)
order1 <- t(order1)
sum <- apply(order1,1,sum) 
order1 <- cbind(order1,sum)
order1 <- as.data.frame(order1)
order1 <- order1[order(order1[,"sum"],decreasing = T),]
write.table(order1, file = "Orderorder.txt",sep = '\t',)

order1 <- subset(order1, select = -sum)
order2 <- t(order1)
order2 <- subset(order2, select = -Unassigned)
order2 <- subset(order2, select = -Unclassified)
order2 <- t(order2)
#如果需要增减物种的数目，请调整此处的数值
order2 <- order2[1:10,]
order2 <- t(order2)
sum <- apply(order2,1,sum) 
Others <- 1-sum
order10 <- cbind(order2,Others)
order10
write.table(order10, file = "OrderTop10.txt",sep = '\t',)
#接下来进行新的样本名称的替换，并保证图像中样本的绘制顺序与预期一致。

sample.name <- read.table("mapping.txt",head=T,colClasses=c("character","character"),sep = "\t")
sample.name
order10 <- as.data.frame(order10)
order10$SampleID <- rownames(order10)
order10 <- merge(sample.name,order10,by = "SampleID",sort = F)
rownames(order10) <- order10$SampleID
order10 <- subset(order10, select = -Group)
order10 <- subset(order10, select = -SampleID)
#order10 <- subset(order10, select = -V3)
order10 <- t(order10)
#图像绘制

#绘制的是冲积条形图所使用的绘图包是“ggplot2+ggalluvial”。
library(reshape2)
taxon <- melt(order10)
taxon
colnames(taxon) <- c("Taxon","variable","value")
library(ggplot2)
#install.packages("ggalluvial")

library(ggalluvial)
#进行绘图坐标的匹配并生成基本的冲积条形图。

p <- ggplot(data = taxon,aes(x = variable, y = value, alluvium = Taxon, stratum = Taxon))
p
p1 <- p + geom_alluvium(aes(fill = Taxon),alpha = .5,width = 0.6) + geom_stratum(aes(fill = Taxon),width = 0.6)
p1

#修改横纵坐标轴标题。

p2 <- p1 + ylab(label = "Relative Abundance") + xlab(label = "")
p2
#修改图像配色。

p3 <- p2 + scale_fill_manual(values = Palette)
p3
#调整绘图区主题。

p4 <- p3 + theme_bw()+ theme(panel.grid=element_blank()) + theme(panel.border = element_blank()) + theme(panel.background=element_rect(fill='transparent',color='black'),plot.margin = unit(c(3,5,1,1),"mm"))
p4
#调整坐标轴文字。

p5 <- p4 + theme(axis.text.x=element_text(colour="black",size=12,angle = 270,vjust = 0.5,hjust = 0)) + 
  theme(axis.text.y=element_text(colour = "black",size = 12)) + 
  theme(axis.line = element_line(colour = "black"))+ 
  theme(axis.title.y = element_text(size = 14,
                                    margin = unit(c(0,1,0,1),"lines"))) 
p5

#去除横纵坐标与边框之间的空白。

p6 <- p5 + scale_y_continuous(limits = c(0,1),expand = c(0,0))
p6
#调整图例。

p7 <- p6 + theme(legend.text = element_text(colour = "black",size = 12)) + theme(legend.title = element_text(size = 18,colour = "black"))
p7
#修改字体。

p8 <- p7 + theme(text = element_text(family = "Times"))
p8

#输出图像。
ggsave("orderTop10.pdf", width = 18, height = 10, units="cm")
ggsave("orderTop10.tiff", width = 18, height = 10, units = "cm")
ggsave("orderTop10.png", width = 18, height = 10, units = "cm")




##科水平top10
#色卡，建议收藏留用，

Palette <- c("#B2182B","#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7","#CC6666","#9999CC","#66CC99","#999999","#ADD1E5")

#数据表格读取
family <- read.table ("family.txt", header = TRUE, sep ="\t")
family
#首先倒入绘图数据，因为是用门水平做例子，因此导入的是L2水平的结果文件，L2-L6分别代表门、纲、目、科和属。
#as.data.frame(),as.array()等函数，可以对已有的常用的类型进行转换。as.data.frame()检查对象是否是数据帧，或者是否可以强制执行
family <- as.data.frame(family)
?as.data.frame
#删除OTU相对丰度表中，有关p_的字样
taxon <- gsub(".*(f__)","",family$Taxon)
#gsub()可以用于字段的删减、增补、替换和切割，可以处理一个字段也可以处理由字段组成的向量。
?gsub
#删除Taxon中有关Bacteria;的字样
#taxon <- gsub(".*(Bacteria;)","",taxon)
#删除taxon中的Other字样
#taxon <- gsub(".*(Other)","",taxon)

taxon <- gsub("(\\[)","",taxon)
taxon <- gsub("(\\])","",taxon)
#将新产生的表头与原来的表头进行替换
family$Taxon <- taxon
family
family[family==""] <- NA
family <- na.omit(family)
family <- t(family)
colnames(family) <- family[1,]
family <- family[-1,]

#使用gsub从分类学名称一列中将所需要的物种名提取出来，??????第二行命令，门水平的数据为p，其它分类学水平依次为c、o、f和g。

#接下来需要对物种根据其在所有样本中的平均丰度进行排序，保留丰度排在前10的物种，对其它物种的丰度进行求和，将其命名为Others。


family <- as.matrix(family)
family1 <- matrix(as.numeric(family),nrow = nrow(family))
rownames(family1) <- rownames(family)
colnames(family1) <- colnames(family)
family1 <- t(family1)
sum <- apply(family1,1,sum) 
family1 <- cbind(family1,sum)
family1 <- as.data.frame(family1)
family1 <- family1[order(family1[,"sum"],decreasing = T),]
write.table(family1, file = "Familyorder.txt",sep = '\t',)

family1 <- subset(family1, select = -sum)
family2 <- t(family1)
family2 <- subset(family2, select = -Unassigned)
family2 <- subset(family2, select = -Unclassified)
family2 <- t(family2)
#如果需要增减物种的数目，请调整此处的数值
family2 <- family2[1:12,]
family2 <- t(family2)
sum <- apply(family2,1,sum) 
Others <- 1-sum
family10 <- cbind(family2,Others)
family10
write.table(family10, file = "FamilyTop10.txt",sep = '\t',)
#接下来进行新的样本名称的替换，并保证图像中样本的绘制顺序与预期一致。

sample.name <- read.table("mapping.txt",head=T,colClasses=c("character","character"),sep = "\t")
sample.name
family10 <- as.data.frame(family10)
family10$SampleID <- rownames(family10)
family10 <- merge(sample.name,family10,by = "SampleID",sort = F)
rownames(family10) <- family10$SampleID
family10 <- subset(family10, select = -Group)
family10 <- subset(family10, select = -SampleID)
#family10 <- subset(family10, select = -V3)
family10 <- t(family10)
family10
#图像绘制

#绘制的是冲积条形图所使用的绘图包是“ggplot2+ggalluvial”。
library(reshape2)
taxon <- melt(family10)
taxon
colnames(taxon) <- c("Taxon","variable","value")
library(ggplot2)
#install.packages("ggalluvial")

library(ggalluvial)
#进行绘图坐标的匹配并生成基本的冲积条形图。

p <- ggplot(data = taxon,aes(x = variable, y = value, alluvium = Taxon, stratum = Taxon))
p
p1 <- p + geom_alluvium(aes(fill = Taxon),alpha = .5,width = 0.6) + geom_stratum(aes(fill = Taxon),width = 0.6)
p1

#修改横纵坐标轴标题。

p2 <- p1 + ylab(label = "Relative Abundance") + xlab(label = "")
p2
#修改图像配色。

p3 <- p2 + scale_fill_manual(values = Palette)
p3
#调整绘图区主题。

p4 <- p3 + theme_bw()+ theme(panel.grid=element_blank()) + theme(panel.border = element_blank()) + theme(panel.background=element_rect(fill='transparent',color='black'),plot.margin = unit(c(3,5,1,1),"mm"))
p4
#调整坐标轴文字。

p5 <- p4 + theme(axis.text.x=element_text(colour="black",size=12,angle = 270,vjust = 0.5,hjust = 0)) + 
  theme(axis.text.y=element_text(colour = "black",size = 12)) + 
  theme(axis.line = element_line(colour = "black"))+ 
  theme(axis.title.y = element_text(size = 14,
                                    margin = unit(c(0,1,0,1),"lines"))) 
p5

#去除横纵坐标与边框之间的空白。

p6 <- p5 + scale_y_continuous(limits = c(0,1),expand = c(0,0))
p6
#调整图例。

p7 <- p6 + theme(legend.text = element_text(colour = "black",size = 12)) + theme(legend.title = element_text(size = 18,colour = "black"))
p7
#修改字体。

p8 <- p7 + theme(text = element_text(family = "Times"))
p8

#输出图像。
ggsave("familyTop10.pdf", width = 18, height = 10, units="cm")
ggsave("familyTop10.tiff", width = 18, height = 10, units = "cm")
ggsave("familyTop10.png", width = 18, height = 10, units = "cm")



##属水平top10
#色卡，建议收藏留用，

Palette <- c("#B2182B","#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7","#CC6666","#9999CC","#66CC99","#999999","#ADD1E5")

#数据表格读取
genus <- read.table ("genus.txt", header = TRUE, sep ="\t")
genus
#首先倒入绘图数据，因为是用门水平做例子，因此导入的是L2水平的结果文件，L2-L6分别代表门、纲、目、科和属。
#as.data.frame(),as.array()等函数，可以对已有的常用的类型进行转换。as.data.frame()检查对象是否是数据帧，或者是否可以强制执行
genus <- as.data.frame(genus)
?as.data.frame
#删除OTU相对丰度表中，有关p_的字样
taxon <- gsub(".*(g__)","",genus$Taxon)
#gsub()可以用于字段的删减、增补、替换和切割，可以处理一个字段也可以处理由字段组成的向量。
?gsub
#删除Taxon中有关Bacteria;的字样
#taxon <- gsub(".*(Bacteria;)","",taxon)
#删除taxon中的Other字样
#taxon <- gsub(".*(Other)","",taxon)

taxon <- gsub("(\\[)","",taxon)
taxon <- gsub("(\\])","",taxon)
#将新产生的表头与原来的表头进行替换
genus$Taxon <- taxon
genus
genus[genus==""] <- NA
genus <- na.omit(genus)
genus <- t(genus)
colnames(genus) <- genus[1,]
genus <- genus[-1,]

#使用gsub从分类学名称一列中将所需要的物种名提取出来，??????第二行命令，门水平的数据为p，其它分类学水平依次为c、o、f和g。

#接下来需要对物种根据其在所有样本中的平均丰度进行排序，保留丰度排在前10的物种，对其它物种的丰度进行求和，将其命名为Others。


genus <- as.matrix(genus)
genus1 <- matrix(as.numeric(genus),nrow = nrow(genus))
rownames(genus1) <- rownames(genus)
colnames(genus1) <- colnames(genus)
genus1 <- t(genus1)
sum <- apply(genus1,1,sum) 
genus1 <- cbind(genus1,sum)
genus1 <- as.data.frame(genus1)
genus1 <- genus1[order(genus1[,"sum"],decreasing = T),]
write.table(genus1, file = "Genusorder.txt",sep = '\t',)


genus1 <- subset(genus1, select = -sum)
genus2 <- t(genus1)
genus2 <- subset(genus2, select = -Unassigned)
genus2 <- t(genus2)
#如果需要增减物种的数目，请调整此处的数值
genus2 <- genus2[1:10,]
genus2 <- t(genus2)
sum <- apply(genus2,1,sum) 
Others <- 1-sum
genus10 <- cbind(genus2,Others)
genus10
write.table(genus10, file = "GenusTop10.txt",sep = '\t',)
#接下来进行新的样本名称的替换，并保证图像中样本的绘制顺序与预期一致。

sample.name <- read.table("mapping.txt",head=T,colClasses=c("character","character"),sep = "\t")
sample.name
genus10 <- as.data.frame(genus10)
genus10$SampleID <- rownames(genus10)
genus10 <- merge(sample.name,genus10,by = "SampleID",sort = F)
rownames(genus10) <- genus10$SampleID
genus10 <- subset(genus10, select = -Group)
genus10 <- subset(genus10, select = -SampleID)
#genus10 <- subset(genus10, select = -V3)
genus10 <- t(genus10)
#图像绘制

#绘制的是冲积条形图所使用的绘图包是“ggplot2+ggalluvial”。
library(reshape2)
taxon <- melt(genus10)
taxon
colnames(taxon) <- c("Taxon","variable","value")
library(ggplot2)
#install.packages("ggalluvial")

library(ggalluvial)
#进行绘图坐标的匹配并生成基本的冲积条形图。

p <- ggplot(data = taxon,aes(x = variable, y = value, alluvium = Taxon, stratum = Taxon))
p
p1 <- p + geom_alluvium(aes(fill = Taxon),alpha = .5,width = 0.6) + geom_stratum(aes(fill = Taxon),width = 0.6)
p1

#修改横纵坐标轴标题。

p2 <- p1 + ylab(label = "Relative Abundance") + xlab(label = "")
p2
#修改图像配色。

p3 <- p2 + scale_fill_manual(values = Palette)
p3
#调整绘图区主题。

p4 <- p3 + theme_bw()+ theme(panel.grid=element_blank()) + theme(panel.border = element_blank()) + theme(panel.background=element_rect(fill='transparent',color='black'),plot.margin = unit(c(3,5,1,1),"mm"))
p4
#调整坐标轴文字。

p5 <- p4 + theme(axis.text.x=element_text(colour="black",size=12,angle = 270,vjust = 0.5,hjust = 0)) + 
  theme(axis.text.y=element_text(colour = "black",size = 12)) + 
  theme(axis.line = element_line(colour = "black"))+ 
  theme(axis.title.y = element_text(size = 14,
                                    margin = unit(c(0,1,0,1),"lines"))) 
p5

#去除横纵坐标与边框之间的空白。

p6 <- p5 + scale_y_continuous(limits = c(0,1),expand = c(0,0))
p6
#调整图例。

p7 <- p6 + theme(legend.text = element_text(colour = "black",size = 12)) + theme(legend.title = element_text(size = 18,colour = "black"))
p7
#修改字体。

p8 <- p7 + theme(text = element_text(family = "Times"))
p8

#输出图像。
ggsave("genusTop10.pdf", width = 18, height = 10, units="cm")
ggsave("genusTop10.tiff", width = 18, height = 10, units = "cm")
ggsave("genusTop10.png", width = 18, height = 10, units = "cm")




##种水平top10
#色卡，建议收藏留用，

Palette <- c("#B2182B","#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7","#CC6666","#9999CC","#66CC99","#999999","#ADD1E5")

#数据表格读取
species <- read.table ("species.txt", header = TRUE, sep ="\t")
species
#首先倒入绘图数据，因为是用门水平做例子，因此导入的是L2水平的结果文件，L2-L6分别代表门、纲、目、科和属。
#as.data.frame(),as.array()等函数，可以对已有的常用的类型进行转换。as.data.frame()检查对象是否是数据帧，或者是否可以强制执行
species <- as.data.frame(species)
?as.data.frame
#删除OTU相对丰度表中，有关p_的字样
taxon <- gsub(".*(s__)","",species$Taxon)
#gsub()可以用于字段的删减、增补、替换和切割，可以处理一个字段也可以处理由字段组成的向量。
?gsub
#删除Taxon中有关Bacteria;的字样
#taxon <- gsub(".*(Bacteria;)","",taxon)
#删除taxon中的Other字样
#taxon <- gsub(".*(Other)","",taxon)

taxon <- gsub("(\\[)","",taxon)
taxon <- gsub("(\\])","",taxon)
#将新产生的表头与原来的表头进行替换
species$Taxon <- taxon
species
species[species==""] <- NA
species <- na.omit(species)
species <- t(species)
colnames(species) <- species[1,]
species <- species[-1,]

#使用gsub从分类学名称一列中将所需要的物种名提取出来，??????第二行命令，门水平的数据为p，其它分类学水平依次为c、o、f和g。

#接下来需要对物种根据其在所有样本中的平均丰度进行排序，保留丰度排在前10的物种，对其它物种的丰度进行求和，将其命名为Others。


species <- as.matrix(species)
species1 <- matrix(as.numeric(species),nrow = nrow(species))
rownames(species1) <- rownames(species)
colnames(species1) <- colnames(species)
species1 <- t(species1)
sum <- apply(species1,1,sum) 
species1 <- cbind(species1,sum)
species1 <- as.data.frame(species1)
species1 <- species1[order(species1[,"sum"],decreasing = T),]
write.table(species1, file = "Speciesorder.txt",sep = '\t',)

species1 <- subset(species1, select = -sum)
species2 <- t(species1)
species2 <- subset(species2, select = -Unassigned)
species2 <- subset(species2, select = -Unclassified)
species2 <- t(species2)
#如果需要增减物种的数目，请调整此处的数值
species2 <- species2[1:10,]
species2 <- t(species2)
sum <- apply(species2,1,sum) 
Others <- 1-sum
species10 <- cbind(species2,Others)
species10
write.table(species10, file = "SpeciesTop10.txt",sep = '\t',)
#接下来进行新的样本名称的替换，并保证图像中样本的绘制顺序与预期一致。

sample.name <- read.table("mapping.txt",head=T,colClasses=c("character","character"),sep = "\t")
sample.name
species10 <- as.data.frame(species10)
species10$SampleID <- rownames(species10)
species10 <- merge(sample.name,species10,by = "SampleID",sort = F)
rownames(species10) <- species10$SampleID
species10 <- subset(species10, select = -Group)
species10 <- subset(species10, select = -SampleID)
#species10 <- subset(species10, select = -V3)
species10 <- t(species10)
#图像绘制

#绘制的是冲积条形图所使用的绘图包是“ggplot2+ggalluvial”。
library(reshape2)
taxon <- melt(species10)
taxon
colnames(taxon) <- c("Taxon","variable","value")
library(ggplot2)
#install.packages("ggalluvial")

library(ggalluvial)
#进行绘图坐标的匹配并生成基本的冲积条形图。

p <- ggplot(data = taxon,aes(x = variable, y = value, alluvium = Taxon, stratum = Taxon))
p
p1 <- p + geom_alluvium(aes(fill = Taxon),alpha = .5,width = 0.6) + geom_stratum(aes(fill = Taxon),width = 0.6)
p1

#修改横纵坐标轴标题。

p2 <- p1 + ylab(label = "Relative Abundance") + xlab(label = "")
p2
#修改图像配色。

p3 <- p2 + scale_fill_manual(values = Palette)
p3
#调整绘图区主题。

p4 <- p3 + theme_bw()+ theme(panel.grid=element_blank()) + theme(panel.border = element_blank()) + theme(panel.background=element_rect(fill='transparent',color='black'),plot.margin = unit(c(3,5,1,1),"mm"))
p4
#调整坐标轴文字。

p5 <- p4 + theme(axis.text.x=element_text(colour="black",size=12,angle = 270,vjust = 0.5,hjust = 0)) + 
  theme(axis.text.y=element_text(colour = "black",size = 12)) + 
  theme(axis.line = element_line(colour = "black"))+ 
  theme(axis.title.y = element_text(size = 14,
                                    margin = unit(c(0,1,0,1),"lines"))) 
p5

#去除横纵坐标与边框之间的空白。

p6 <- p5 + scale_y_continuous(limits = c(0,1),expand = c(0,0))
p6
#调整图例。

p7 <- p6 + theme(legend.text = element_text(colour = "black",size = 12)) + theme(legend.title = element_text(size = 18,colour = "black"))
p7
#修改字体。

p8 <- p7 + theme(text = element_text(family = "Times"))
p8

#输出图像。
ggsave("speciesTop10.pdf", width = 18, height = 10, units="cm")
ggsave("speciesTop10.tiff", width = 18, height = 10, units = "cm")
ggsave("speciesTop10.png", width = 18, height = 10, units = "cm")

####热图####
##TOP40属在不同处理间的丰度用热图展示





####弦图####
##TOP20属在不同处理间的丰度用弦图展示
library(circlize) #使用该包绘制 circos 图
library(reshape2)
library(ComplexHeatmap) #可用此包添加图例
library(grid) #可用此包调整画板

#读取 taxonomy.txt 的内容，获取“OTU/分类”排序，OTU 名称
taxonomy <- read.delim('taxonomy_circlize2.txt', sep = '\t', stringsAsFactors = FALSE)
taxonomy
tax_phylum <- unique(taxonomy$phylum) #分类单位的内容及数量
tax_phylum
taxonomy$phylum <- factor(taxonomy$phylum, levels = tax_phylum)
all_otu <- taxonomy$OTU_ID
all_otu
taxonomy$OTU_ID <- factor(taxonomy$OTU_ID, levels = all_otu)
taxonomy

#读取 group.txt 的内容，获取“样本/分组”排序，样本名称
group <- read.delim('group_circlize.txt', sep = '\t', stringsAsFactors = FALSE)
group
all_group <- unique(group$group_ID) #样本的数量及内容
all_group
group$group_ID <- factor(group$group_ID, levels = all_group)
all_sample <- group$SampleID
all_sample
#读取 otu_table.txt，排序 OTU 和样本
otu_table <- read.delim('otu_table_circlize.txt', sep = '\t')
otu_table
#otu_table <- merge(taxonomy, otu_table, by = 'OTU_ID')
#otu_table
#otu_table <- otu_table[order(otu_table$phylum, otu_table$OTU_ID), ]#将phylum根据genus进行排序
#otu_table
rownames(otu_table) <- otu_table$OTU_ID
otu_table <- otu_table[all_sample]
otu_table

##生成作图数据
#circlize 外圈属性数据
all_ID <- c(all_otu, all_sample)
all_ID
accum_otu <- rowSums(otu_table)#计算每个属在所有样本中的数量
accum_otu
accum_sample <- colSums(otu_table)#计算每个样本中所有otu的数量
accum_sample 
all_ID_xlim <- cbind(rep(0, length(all_ID)),data.frame(c(accum_otu, accum_sample)))
all_ID_xlim

#circlize 内圈连线数据
otu_table$OTU_ID <- all_otu
otu_table
plot_data <- melt(otu_table, id = 'OTU_ID') #此处使用了reshape2包中的melt()命令
plot_data
colnames(plot_data)[2] <- 'sample_ID'
plot_data$OTU_ID <- factor(plot_data$OTU_ID, levels = all_otu)
plot_data$sample_ID <- factor(plot_data$sample_ID, levels = all_sample)
plot_data <- plot_data[order(plot_data$OTU_ID, plot_data$sample_ID), ]
plot_data <- plot_data[c(2, 1, 3, 3)]
plot_data

#颜色设置
color_otu <- c('#8DD3C7', '#FFFFB3', '#BEBADA', '#FB8072', '#80B1D3', '#FDB462', '#B3DE69', '#FCCDE5', '#BC80BD', '#CCEBC5', '#FFED6F', '#E41A1C', '#377EB8', '#4DAF4A', '#984EA3', '#FF7F00', '#FFFF33', '#A65628', '#F781BF', '#66C2A5')
color_sample <- c('#6181BD', '#F34800', '#64A10E', '#FF00FF', '#c7475b', '#049a0b',"#B2182B","#E69F00","#56B4E9","#009E73","#F0E442","#0072B2")
color_phylum <- c('#BEAED4', '#FDC086', '#FFFF99', '#386CB0')
color_group <- c('#4253ff', '#ff4308',"#CC6666","#9999CC")

names(color_otu) <- all_otu
color_otu
names(color_sample) <- all_sample
color_sample

####circlize 绘图
pdf('circlize_plot.pdf', width = 12, height = 8)
circle_size = unit(1, 'snpc')



##整体布局
gap_size <- c(rep(3, length(all_otu) - 1), 6, rep(3, length(all_sample) - 1), 6)
circos.par(cell.padding = c(0, 0, 0, 0), start.degree = 270, gap.degree = gap_size)
circos.initialize(factors = factor(all_ID, levels = all_ID), xlim = all_ID_xlim)



##绘制 OTU 分类、样本分组区块（第一圈）
circos.trackPlotRegion(
  ylim = c(0, 1), track.height = 0.03, bg.border = NA, 
  panel.fun = function(x, y) {
    sector.index = get.cell.meta.data('sector.index')
    xlim = get.cell.meta.data('xlim')
    ylim = get.cell.meta.data('ylim')
  } )

for (i in 1:length(tax_phylum)) {
  tax_OTU <- {subset(taxonomy, phylum == tax_phylum[i])}$OTU_ID
  highlight.sector(tax_OTU, track.index = 1, col = color_phylum[i], text = tax_phylum[i], cex = 0.5, text.col = 'black', niceFacing = FALSE)
}

for (i in 1:length(all_group)) {
  group_sample <- {subset(group, group_ID == all_group[i])}$SampleID
  highlight.sector(group_sample, track.index = 1, col = color_group[i], text = all_group[i], cex = 0.7, text.col = 'black', niceFacing = FALSE)
}



##各 OTU、样本绘制区
#添加百分比注释（第二圈）

circos.trackPlotRegion(
  ylim = c(0, 1), track.height = 0.05, bg.border = NA, 
  panel.fun = function(x, y) {
    sector.index = get.cell.meta.data('sector.index')
    xlim = get.cell.meta.data('xlim')
    ylim = get.cell.meta.data('ylim')
  } )

circos.track(
  track.index = 2, bg.border = NA, 
  panel.fun = function(x, y) {
    xlim = get.cell.meta.data('xlim')
    ylim = get.cell.meta.data('ylim')
    sector.name = get.cell.meta.data('sector.index')
    xplot = get.cell.meta.data('xplot')
    
    by = ifelse(abs(xplot[2] - xplot[1]) > 30, 0.25, 1)
    for (p in c(0, seq(by, 1, by = by))) circos.text(p*(xlim[2] - xlim[1]) + xlim[1], mean(ylim) + 0.4, paste0(p*100, '%'), cex = 0.4, adj = c(0.5, 0), niceFacing = FALSE)
    
    circos.lines(xlim, c(mean(ylim), mean(ylim)), lty = 3)
  } )

#绘制 OTU、样本主区块（第三圈）
circos.trackPlotRegion(
  ylim = c(0, 1), track.height = 0.03, bg.col = c(color_otu, color_sample), bg.border = NA, track.margin = c(0, 0.01),
  panel.fun = function(x, y) {
    xlim = get.cell.meta.data('xlim')
    sector.name = get.cell.meta.data('sector.index')
    circos.axis(h = 'top', labels.cex = 0.4, major.tick.length = 0.4, labels.niceFacing = FALSE)
    circos.text(mean(xlim), 0.2, sector.name, cex = 0.4, niceFacing = FALSE, adj = c(0.5, 0))
  } )



#绘制 OTU、样本副区块（第四圈）
circos.trackPlotRegion(ylim = c(0, 1), track.height = 0.03, track.margin = c(0, 0.01))



##绘制 OTU-样本关联连线（最内圈）
for (i in seq_len(nrow(plot_data))) {
  circos.link(
    plot_data[i,2], c(accum_otu[plot_data[i,2]], accum_otu[plot_data[i,2]] - plot_data[i,4]),
    plot_data[i,1], c(accum_sample[plot_data[i,1]], accum_sample[plot_data[i,1]] - plot_data[i,3]),
    col = paste0(color_otu[plot_data[i,2]], '70'), border = NA )
  
  circos.rect(accum_otu[plot_data[i,2]], 0, accum_otu[plot_data[i,2]] - plot_data[i,4], 1, sector.index = plot_data[i,2], col = color_sample[plot_data[i,1]], border = NA)
  circos.rect(accum_sample[plot_data[i,1]], 0, accum_sample[plot_data[i,1]] - plot_data[i,3], 1, sector.index = plot_data[i,1], col = color_otu[plot_data[i,2]], border = NA)
  
  accum_otu[plot_data[i,2]] = accum_otu[plot_data[i,2]] - plot_data[i,4]
  accum_sample[plot_data[i,1]] = accum_sample[plot_data[i,1]] - plot_data[i,3]
}



##添加图例

otu_legend <- Legend(
  at = all_otu, labels = taxonomy$detail, labels_gp = gpar(fontsize = 8),    
  grid_height = unit(0.5, 'cm'), grid_width = unit(0.5, 'cm'), type = 'points', pch = NA, background = color_otu)

pushViewport(viewport(x = 0.85, y = 0.5))
grid.draw(otu_legend)
upViewport()


##清除 circlize 样式并关闭画板
circos.clear()
dev.off()


####差异分析####
##差异分析，包括丰度差异显著性比较、物种丰度上下调表达（edgeR、DESeq2）、经典的Lefse方法
#有许多可用的工具可用于差异丰度评估和推断。例如，基于广义线性模型(GLM)的方法，包括通过负二项式分布(DESeq2和edgeR)或使用对数计数(百万分之一)和高斯分布(limma)假设的基于微生物组数据或基因表达数据的DESeq2、edgeR和limma模型计数。可以确定在不同组中具有不同丰度的核心微生物。

####ANOVA####

##方法一：进行处理间各物种ANOVA检验,不进行物种丰富度划分的原始数据，往往数据量很大而无法进行完全的展示，故进行数据筛选或者丰富度的划分后结果更好
rm(list=ls())
setwd("E:/博士数据2/小麦微生物组数据/customer_backup/otus/ANOVA/")
#首先导入数据文件

##门水平phylum
#采用scale_test对数据进行对数转换
phylum <- read.delim('./Phylum.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
phylum

##排序，便于筛选TOP物种
phylum <- as.matrix(phylum)
phylum1 <- matrix(as.numeric(phylum),nrow = nrow(phylum))
rownames(phylum1) <- rownames(phylum)
colnames(phylum1) <- colnames(phylum)
sum <- apply(phylum1,1,sum) 
phylum1 <- cbind(phylum1,sum)
phylum1 <- as.data.frame(phylum1)
phylum1 <- phylum1[order(phylum1[,"sum"],decreasing = T),]
phylum1
write.table(phylum1, file = "phylum_order.txt",sep = '\t',)

#如果需要增减物种的数目，请调整此处的数值
#phylum2 <- phylum1[1:20,]
#phylum2 <- t(phylum2)
#sum <- apply(phylum2,1,sum) 
#Others <- 1-sum
#phylum20 <- cbind(phylum2,Others)
#phylum20
#write.table(phylum20, file = "phylum_AT_Top20.txt",sep = '\t',)



##纲水平class
#采用scale_test对数据进行对数转换
class <- read.delim('./class.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
class

##排序，便于筛选TOP物种
class <- as.matrix(class)
class1 <- matrix(as.numeric(class),nrow = nrow(class))
rownames(class1) <- rownames(class)
colnames(class1) <- colnames(class)
sum <- apply(class1,1,sum) 
class1 <- cbind(class1,sum)
class1 <- as.data.frame(class1)
class1 <- class1[order(class1[,"sum"],decreasing = T),]
class1
write.table(class1, file = "class_order.txt",sep = '\t',)

#如果需要增减物种的数目，请调整此处的数值
class2 <- subset(class1,select = -sum)
class2
class40 <- class2[1:40,]
class40 <- t(class40)
class40
write.table(class40, file = "class_Top40.txt",sep = '\t',)



##目水平order
#采用scale_test对数据进行对数转换
order <- read.delim('./order.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
order

##排序，便于筛选TOP物种
order <- as.matrix(order)
order1 <- matrix(as.numeric(order),nrow = nrow(order))
rownames(order1) <- rownames(order)
colnames(order1) <- colnames(order)
sum <- apply(order1,1,sum) 
order1 <- cbind(order1,sum)
order1 <- as.data.frame(order1)
order1 <- order1[order(order1[,"sum"],decreasing = T),]
order1
write.table(order1, file = "order_order.txt",sep = '\t',)

#如果需要增减物种的数目，请调整此处的数值
order2 <- subset(order1,select = -sum)
order2
order40 <- order2[1:40,]
order40 <- t(order40)
order40
write.table(order40, file = "order_Top40.txt",sep = '\t',)


##科水平family
#采用scale_test对数据进行对数转换
family <- read.delim('./family.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
family

##排序，便于筛选TOP物种
family <- as.matrix(family)
family1 <- matrix(as.numeric(family),nrow = nrow(family))
rownames(family1) <- rownames(family)
colnames(family1) <- colnames(family)
sum <- apply(family1,1,sum) 
family1 <- cbind(family1,sum)
family1 <- as.data.frame(family1)
family1 <- family1[order(family1[,"sum"],decreasing = T),]
family1
write.table(family1, file = "family_order.txt",sep = '\t',)

#如果需要增减物种的数目，请调整此处的数值
family2 <- subset(family1,select = -sum)
family2
family40 <- family2[1:40,]
family40 <- t(family40)
family40
write.table(family40, file = "family_Top40.txt",sep = '\t',)



##属水平genus
#采用scale_test对数据进行对数转换
genus <- read.delim('./genus.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
genus

##排序，便于筛选TOP物种
genus <- as.matrix(genus)
genus1 <- matrix(as.numeric(genus),nrow = nrow(genus))
rownames(genus1) <- rownames(genus)
colnames(genus1) <- colnames(genus)
sum <- apply(genus1,1,sum) 
genus1 <- cbind(genus1,sum)
genus1 <- as.data.frame(genus1)
genus1 <- genus1[order(genus1[,"sum"],decreasing = T),]
genus1
write.table(genus1, file = "genus_order.txt",sep = '\t',)

#如果需要增减物种的数目，请调整此处的数值
genus2 <- subset(genus1,select = -sum)
genus2
genus40 <- genus2[1:40,]
genus40 <- t(genus40)
genus40
write.table(genus40, file = "genus_Top40.txt",sep = '\t',)


##种水平species
#采用scale_test对数据进行对数转换
species <- read.delim('./species.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
species

##排序，便于筛选TOP物种
species <- as.matrix(species)
species1 <- matrix(as.numeric(species),nrow = nrow(species))
rownames(species1) <- rownames(species)
colnames(species1) <- colnames(species)
sum <- apply(species1,1,sum) 
species1 <- cbind(species1,sum)
species1 <- as.data.frame(species1)
species1 <- species1[order(species1[,"sum"],decreasing = T),]
species1
write.table(species1, file = "species_order.txt",sep = '\t',)

#如果需要增减物种的数目，请调整此处的数值
species2 <- subset(species1,select = -sum)
species2
species40 <- species2[1:40,]
species40 <- t(species40)
species40
write.table(species40, file = "species_Top40.txt",sep = '\t',)

####对数据进行对数转换和开平方根

##门水平phylum
phylum <- read.delim('./phylum_order.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
phylum
##首先进行对数转换
phylum <- apply(phylum, 2, function(x) log2(x+1))
phylum
##由于仅仅利用对数转换后的数值绘制出的热图颜色区分度不大，故对对数转换后的数值进一步求平方根
phylum <- apply(phylum, 2, function(x) sqrt(x))
phylum

write.table(phylum, 'phylum_scale_test.txt', row.names = TRUE, sep = '\t', quote = FALSE)



##纲水平class
class40 <- read.delim('./class_Top40.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
class40
##首先进行对数转换
class40 <- apply(class40, 2, function(x) log2(x+1))
class40
##由于仅仅利用对数转换后的数值绘制出的热图颜色区分度不大，故对对数转换后的数值进一步求平方根
class40 <- apply(class40, 2, function(x) sqrt(x))
class40

write.table(class40, 'class_scale_test.txt', row.names = TRUE, sep = '\t', quote = FALSE)



##目水平order
order40 <- read.delim('./order_Top40.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
order40
##首先进行对数转换
order40 <- apply(order40, 2, function(x) log2(x+1))
order40
##由于仅仅利用对数转换后的数值绘制出的热图颜色区分度不大，故对对数转换后的数值进一步求平方根
order40 <- apply(order40, 2, function(x) sqrt(x))
order40

write.table(order40, 'order_scale_test.txt', row.names = TRUE, sep = '\t', quote = FALSE)




##科水平family
family40 <- read.delim('./family_Top40.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
family40
##首先进行对数转换
family40 <- apply(family40, 2, function(x) log2(x+1))
family40
##由于仅仅利用对数转换后的数值绘制出的热图颜色区分度不大，故对对数转换后的数值进一步求平方根
family40 <- apply(family40, 2, function(x) sqrt(x))
family40

write.table(family40, 'family_scale_test.txt', row.names = TRUE, sep = '\t', quote = FALSE)


##属水平genus
genus40 <- read.delim('./genus_Top20.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
genus40
##首先进行对数转换
genus40 <- apply(genus40, 2, function(x) log2(x+1))
genus40
##由于仅仅利用对数转换后的数值绘制出的热图颜色区分度不大，故对对数转换后的数值进一步求平方根
genus40 <- apply(genus40, 2, function(x) sqrt(x))
genus40

write.table(genus40, 'genus_scale_test.txt', row.names = TRUE, sep = '\t', quote = FALSE)


##种水平species
species40 <- read.delim('./species_Top20.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
species40
##首先进行对数转换
species40 <- apply(species40, 2, function(x) log2(x+1))
species40
##由于仅仅利用对数转换后的数值绘制出的热图颜色区分度不大，故对对数转换后的数值进一步求平方根
species40 <- apply(species40, 2, function(x) sqrt(x))
species40

write.table(species40, 'species_scale_test.txt', row.names = TRUE, sep = '\t', quote = FALSE)



metadata<- read.delim('./mapping.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
metadata


phylum2 <- read.delim('./phylum_scale_test.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
phylum2
phylum3 <-t(phylum2)
phylum3

data1 <- merge(x=phylum3,y=metadata,by='row.names')
data1
write.table(data1, 'phylum_data.txt', row.names = FALSE, sep = '\t', quote = FALSE)


class2 <- read.delim('./class_scale_test.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
class2

data2 <- merge(x=class2,y=metadata,by='row.names')
data2
write.table(data2, 'class_data.txt', row.names = FALSE, sep = '\t', quote = FALSE)



order2 <- read.delim('./order_scale_test.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
order2

data3 <- merge(x=order2,y=metadata,by='row.names')
data3
write.table(data3, 'order_data.txt', row.names = FALSE, sep = '\t', quote = FALSE)



family2 <- read.delim('./family_scale_test.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
family2

data4 <- merge(x=family2,y=metadata,by='row.names')
data4
write.table(data4, 'family_data.txt', row.names = FALSE, sep = '\t', quote = FALSE)



genus2 <- read.delim('./genus_scale_test.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
genus2

data5 <- merge(x=genus2,y=metadata,by='row.names')
data5
write.table(data5, 'genus_data.txt', row.names = FALSE, sep = '\t', quote = FALSE)


species2 <- read.delim('./species_scale_test.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
species2

data6 <- merge(x=species2,y=metadata,by='row.names')
data6
write.table(data6, 'species_data.txt', row.names = FALSE, sep = '\t', quote = FALSE)


#生信小白鱼方差分析source
#构建函数，统计各组均值、标准差，执行 ANOVA 以及 Tukey HSD 获得显著性
library(reshape2)
library(multcomp)
library(ggplot2)

aov_tukey <- function(data, groups, values, p = 0.05) {
  stat_anova <- NULL
  abc_list <- NULL
  
  #统计检验
  for (i in groups) {
    for (j in values) {
      
      #单因素 ANOVA，整体差异
      dat <- data[c(i, j)]
      names(dat) <- c('class', 'var')
      dat$class <- factor(dat$class)
      fit <- aov(var~class, dat)
      p_value <- summary(fit)[[1]][1,5]
      
      #单因素 ANOVA 结果整理
      if (p_value < 0.001) sig <- '***'
      else if (p_value >= 0.001 & p_value < 0.01) sig <- '**'
      else if (p_value >= 0.01 & p_value < 0.05) sig <- '*'
      else sig <- ''
      stat_anova <- rbind(stat_anova, c(paste(i, j, sep = '/'), p_value, sig))
      
      #Tukey HSD 检验（multcomp 包），多重比较
      tuk <- cld(glht(fit, alternative = 'two.sided', linfct = mcp(class = 'Tukey')), sig = p, decreasing = TRUE)
      
      #cld() 自动得到了显著性“abc”水平，提取显著性标记“abc”
      sig <- data.frame(tuk$mcletters$Letters, stringsAsFactors = FALSE)
      names(sig) <- 'sig'
      sig$class <- rownames(sig)
      sig$var <- j
      
      #均值标准差统计
      abc_ij <- cbind(aggregate(dat$var, by = list(dat$class), FUN = mean), aggregate(dat$var, by = list(dat$class), FUN = sd)[2])
      names(abc_ij) <- c('class', 'mean', 'sd')
      abc_ij$class <- as.character(abc_ij$class)
      
      #合并结果
      abc_ij$group <- i
      abc_ij <- merge(abc_ij, sig, by = 'class')
      abc_ij <- abc_ij[c(4, 6, 1, 2, 3, 5)]
      abc_list <- rbind(abc_list, abc_ij)
    }
  }
  
  #ggplot2 作图，柱状图
  plot_bar <- ggplot(data = abc_list, aes(x = class, y = mean)) +
    geom_col(aes(fill = class), color = NA, show.legend = FALSE) +
    geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.2) +
    facet_grid(var~group, scale = 'free' , space = 'free_x') +
    geom_text(aes(label = sig, y = mean + sd + 0.3*(mean+sd))) +
    labs(x = '', y = '')
  
  #ggplot2 作图，箱线图
  data <- melt(data[c(groups, values)], id = groups)
  data <- melt(data, id = c('variable', 'value'))
  data <- data[c(3, 1, 4, 2)]
  names(data) <- c('group', 'var', 'class', 'value')
  
  value_max <- aggregate(data$value, by = list(data$group, data$var), FUN = max)
  names(value_max) <- c('group', 'var', 'value')
  value_max <- merge(value_max, abc_list[c(-4, -5)], by = c('group', 'var'))
  
  plot_box <- ggplot(data = data, aes(x = class, y = value)) +
    geom_boxplot(aes(fill = class), show.legend = FALSE) +
    facet_grid(var~group, scale = 'free' , space = 'free_x') +
    geom_text(data = value_max, aes(label = sig, y = value + 0.3*value)) +
    labs(x = '', y = '')
  
  #return
  stat_anova <- data.frame(stat_anova, stringsAsFactors = FALSE)
  names(stat_anova) <- c('group', 'anova_pvalue', 'sig')
  stat_anova$anova_pvalue <- as.numeric(stat_anova$anova_pvalue)
  
  list(stat = list(anova = stat_anova, tukey = abc_list), plot = list(plot_bar = plot_bar, plot_box = plot_box))
}


####phulum门水平
##调用函数运行
data1 <- read.table('phylum_data.txt', header = TRUE)
data1

##将微生物的分类名称提出来进行方差分析时使用
value1 <- colnames(data1)
value1 <- t(value1)
value1
value11 <- paste(value1,collapse = "','")
value11

#groups 指定分组列，可指定多列；values 指定变量列，可指定多列；p 指定显著性水平
#这里作为测试，就先假设这些数据满足方差分析的前提条件了
stat_result <- aov_tukey(data = data1, groups = 'Group', values = c('p__Proteobacteria','p__Acidobacteria','p__Bacteroidetes','p__Firmicutes',
                                                                    'p__Chloroflexi','p__Actinobacteria','p__Fusobacteria','p__Gemmatimonadetes',
                                                                    'p__Patescibacteria','p__Cyanobacteria','p__Verrucomicrobia','p__Planctomycetes',
                                                                    'p__Rokubacteria','p__Nitrospirae','p__Armatimonadetes','p__Fibrobacteres',
                                                                    'p__Spirochaetes','p__Latescibacteria','p__Epsilonbacteraeota','p__WPS.2',
                                                                    'p__Euryarchaeota','p__FCPU426','p__WS2','Unassigned','p__Elusimicrobia',
                                                                    'p__Dependentiae','p__Hydrogenedentes'),p=0.05)

#差异分析统计结果
stat_result$stat$anova
stat_result$stat$tukey

write.table(stat_result$stat$anova, 'phylum_anova.txt', row.names = FALSE, sep = '\t', quote = FALSE)
write.table(stat_result$stat$tukey, 'phylum_tukey.txt', row.names = FALSE, sep = '\t', quote = FALSE)



#显著性标记“abc”作图结果
stat_result$plot$plot_bar
stat_result$plot$plot_box



####class纲水平
##调用函数运行
data2 <- read.table('class_data.txt', header = TRUE)
data2

##将微生物的分类名称提出来进行方差分析时使用
value2 <- colnames(data2)
value2 <- t(value2)
value2
value22 <- paste(value2,collapse = "','")
value22

#groups 指定分组列，可指定多列；values 指定变量列，可指定多列；p 指定显著性水平
#这里作为测试，就先假设这些数据满足方差分析的前提条件了
stat_result <- aov_tukey(data = data2, groups = 'Group', values = c('c__Gammaproteobacteria','c__Bacteroidia','c__Alphaproteobacteria','c__Bacilli',
                                                                    'c__Subgroup_6','c__Anaerolineae','c__Fusobacteriia','c__Deltaproteobacteria',
                                                                    'c__Actinobacteria','c__Acidobacteriia','c__Gemmatimonadetes','c__Blastocatellia_Subgroup_4',
                                                                    'c__Acidimicrobiia','c__Oxyphotobacteria','c__Verrucomicrobiae','c__Clostridia',
                                                                    'c__Thermoleophilia','c__Chloroflexia','c__KD4.96','c__Saccharimonadia',
                                                                    'c__Holophagae','c__Thermoanaerobaculia','c__Negativicutes','c__Parcubacteria',
                                                                    'c__Phycisphaerae','c__S0134_terrestrial_group','c__Longimicrobia','c__NC10',
                                                                    'c__Nitrospira','c__MB.A2.108','c__Dehalococcoidia','c__BD2.11_terrestrial_group',
                                                                    'c__TK10','c__Gitt.GS.136','c__0319.7L14','c__Microgenomatia','c__Subgroup_17',
                                                                    'c__JG30.KF.CM66','c__Ignavibacteria','c__Fibrobacteria'),p=0.05)



#差异分析统计结果
stat_result$stat$anova
stat_result$stat$tukey

write.table(stat_result$stat$anova, 'class_anova.txt', row.names = FALSE, sep = '\t', quote = FALSE)
write.table(stat_result$stat$tukey, 'class_tukey.txt', row.names = FALSE, sep = '\t', quote = FALSE)


####order目水平
##调用函数运行
data3 <- read.table('order_data.txt', header = TRUE)
data3

##将微生物的分类名称提出来进行方差分析时使用
value3 <- colnames(data3)
value3 <- t(value3)
value3
value33 <- paste(value3,collapse = "','")
value33

#groups 指定分组列，可指定多列；values 指定变量列，可指定多列；p 指定显著性水平
#这里作为测试，就先假设这些数据满足方差分析的前提条件了
stat_result <- aov_tukey(data = data3, groups = 'Group', values = c('o__Betaproteobacteriales','o__Bacteroidales','o__Lactobacillales','o__Subgroup_6','o__Rhizobiales',
                                                                    'o__Fusobacteriales','o__SBR1031','o__Solibacterales','o__Bacillales','o__Pasteurellales',
                                                                    'o__Gemmatimonadales','o__Sphingomonadales','o__Cytophagales','o__Myxococcales',
                                                                    'o__Micrococcales','o__Pseudomonadales','o__Blastocatellales','o__Chitinophagales',
                                                                    'o__Microtrichales','o__Clostridiales','o__Azospirillales','o__Anaerolineales',
                                                                    'o__Caulobacterales','o__uncultured_bacterium_c_KD4.96','o__Nostocales','o__Enterobacteriales',
                                                                    'o__Saccharimonadales','o__Subgroup_7','o__Pyrinomonadales','o__Thermomicrobiales',
                                                                    'o__Thermoanaerobaculales','o__Selenomonadales','o__Xanthomonadales','o__Propionibacteriales',
                                                                    'o__Desulfuromonadales','o__Gaiellales','o__Pedosphaerales','o__Solirubrobacterales',
                                                                    'o__uncultured_bacterium_c_S0134_terrestrial_group','o__Phycisphaerales'),p=0.05)

#差异分析统计结果
stat_result$stat$anova
stat_result$stat$tukey

write.table(stat_result$stat$anova, 'order_anova.txt', row.names = FALSE, sep = '\t', quote = FALSE)
write.table(stat_result$stat$tukey, 'order_tukey.txt', row.names = FALSE, sep = '\t', quote = FALSE)


#显著性标记“abc”作图结果
stat_result$plot$plot_bar
stat_result$plot$plot_box

####family目水平
##调用函数运行
data4 <- read.table('family_data.txt', header = TRUE)
data4

##将微生物的分类名称提出来进行方差分析时使用
value4 <- colnames(data4)
value4 <- t(value4)
value4
value44 <- paste(value4,collapse = "','")
value44

#groups 指定分组列，可指定多列；values 指定变量列，可指定多列；p 指定显著性水平
#这里作为测试，就先假设这些数据满足方差分析的前提条件了
stat_result <- aov_tukey(data = data4, groups = 'Group', values = c('f__Prevotellaceae','f__Subgroup_6','f__Streptococcaceae','f__Rhodocyclaceae','f__Burkholderiaceae',
                                                                    'f__Solibacteraceae','f__Pasteurellaceae','f__Gemmatimonadaceae','f__Nitrosomonadaceae',
                                                                    'f__Leptotrichiaceae','f__A4b','f__Bacillaceae','f__Sphingomonadaceae','f__Carnobacteriaceae',
                                                                    'f__uncultured_bacterium_o_SBR1031','f__Porphyromonadaceae','f__Microscillaceae','f__Fusobacteriaceae',
                                                                    'f__Blastocatellaceae','f__Pseudomonadaceae','f__Xanthobacteraceae','f__Neisseriaceae','f__Anaerolineaceae',
                                                                    'f__Micrococcaceae','f__Rhizobiaceae','f__Chitinophagaceae','f__uncultured_bacterium_c_KD4.96',
                                                                    'f__Enterobacteriaceae','f__uncultured_bacterium_o_Subgroup_7','f__uncultured_bacterium_o_Saccharimonadales',
                                                                    'f__Ilumatobacteraceae','f__Pyrinomonadaceae','f__Thermoanaerobaculaceae','f__Veillonellaceae','f__Devosiaceae',
                                                                    'f__Clostridiaceae_1','f__TRA3.20','f__Nocardioidaceae','f__Azospirillaceae','f__Geobacteraceae'),p=0.05)

#差异分析统计结果
stat_result$stat$anova
stat_result$stat$tukey

write.table(stat_result$stat$anova, 'family_anova.txt', row.names = FALSE, sep = '\t', quote = FALSE)
write.table(stat_result$stat$tukey, 'family_tukey.txt', row.names = FALSE, sep = '\t', quote = FALSE)


#显著性标记“abc”作图结果
stat_result$plot$plot_bar
stat_result$plot$plot_box


####genus目水平
##调用函数运行
data5 <- read.delim("genus_data.txt", header = TRUE, sep = '\t', stringsAsFactors = FALSE)
data5

##将微生物的分类名称提出来进行方差分析时使用
value5 <- colnames(data5)
value5 <- t(value5)
value5
value55 <- paste(value5,collapse = "','")
value55


#groups 指定分组列，可指定多列；values 指定变量列，可指定多列；p 指定显著性水平
#这里作为测试，就先假设这些数据满足方差分析的前提条件了
stat_result <- aov_tukey(data = data5, groups = 'Group',values = c('g__Subgroup_6','g__Alloprevotella','g__Streptococcus','g__Haemophilus','g__Bryobacter','g__Gemmatimonadaceae','g__A4b','g__Azoarcus','g__Bacillus',
                                                                   'g__Granulicatella','g__Dechloromonas','g__uncultured_bacterium_o_SBR1031','g__Porphyromonas','g__MND1','g__Fusobacterium','g__Streptobacillus',
                                                                   'g__uncultured_bacterium_f_Microscillaceae','g__Sphingomonas','g__Neisseria','g__Leptotrichia','g__uncultured_bacterium_f_Anaerolineaceae',
                                                                   'g__uncultured_bacterium_c_KD4.96','g__uncultured_bacterium_o_Subgroup_7','g__uncultured_bacterium_o_Saccharimonadales','g__Rothia','g__RB41',
                                                                   'g__Aridibacter','g__Pseudomonas','g__Subgroup_10','g__Escherichia.Shigella','g__Veillonella','g__uncultured_bacterium_f_Xanthobacteraceae',
                                                                   'g__Devosia','g__uncultured_bacterium_f_TRA3.20','g__Hydrogenophaga','g__uncultured_bacterium_f_Pseudomonadaceae','g__Ellin6067','g__Ottowia',
                                                                   'g__Geobacter','g__Azospirillum'),p=0.05)

#差异分析统计结果
stat_result$stat$anova
stat_result$stat$tukey

write.table(stat_result$stat$anova, 'genus_anova.txt', row.names = FALSE, sep = '\t', quote = FALSE)
write.table(stat_result$stat$tukey, 'genus_tukey.txt', row.names = FALSE, sep = '\t', quote = FALSE)


#显著性标记“abc”作图结果
stat_result$plot$plot_bar
stat_result$plot$plot_box


####species目水平
##调用函数运行
data6 <- read.delim("species_data.txt", header = TRUE, sep = '\t', stringsAsFactors = FALSE)
data6

##将微生物的分类名称提出来进行方差分析时使用
value6 <- colnames(data6)
value6 <- t(value6)
value6
value66 <- paste(value6,collapse = "','")
value66


#groups 指定分组列，可指定多列；values 指定变量列，可指定多列；p 指定显著性水平
#这里作为测试，就先假设这些数据满足方差分析的前提条件了
stat_result <- aov_tukey(data = data6, groups = 'Group',values = c('s__uncultured_bacterium_c_Subgroup_6','s__uncultured_bacterium_g_Alloprevotella','s__uncultured_bacterium_g_Streptococcus',
                                                                   's__uncultured_bacterium_g_Haemophilus','s__uncultured_bacterium_g_Bryobacter','s__uncultured_bacterium_f_Gemmatimonadaceae',
                                                                   's__uncultured_bacterium_f_A4b','s__uncultured_bacterium_g_Azoarcus','s__uncultured_bacterium_g_Granulicatella',
                                                                   's__uncultured_bacterium_g_Bacillus','s__uncultured_bacterium_o_SBR1031','s__uncultured_bacterium_g_Porphyromonas',
                                                                   's__uncultured_bacterium_g_MND1','s__uncultured_bacterium_g_Streptobacillus','s__uncultured_bacterium_g_Dechloromonas',
                                                                   's__uncultured_bacterium_f_Microscillaceae','s__uncultured_bacterium_g_Fusobacterium','s__uncultured_bacterium_g_Neisseria',
                                                                   's__uncultured_bacterium_g_Sphingomonas','s__uncultured_bacterium_g_Leptotrichia','s__uncultured_bacterium_f_Anaerolineaceae',
                                                                   's__uncultured_bacterium_c_KD4.96','s__uncultured_bacterium_o_Subgroup_7','s__uncultured_bacterium_o_Saccharimonadales',
                                                                   's__uncultured_bacterium_g_Rothia','s__uncultured_bacterium_g_RB41','s__uncultured_bacterium_g_Aridibacter',
                                                                   's__uncultured_bacterium_g_Escherichia.Shigella','s__uncultured_bacterium_g_Subgroup_10',
                                                                   's__uncultured_bacterium_g_Pseudomonas','s__uncultured_bacterium_g_Veillonella','s__uncultured_bacterium_f_Xanthobacteraceae',
                                                                   's__uncultured_bacterium_g_Devosia','s__uncultured_bacterium_f_TRA3.20','s__uncultured_bacterium_g_Hydrogenophaga',
                                                                   's__uncultured_bacterium_f_Pseudomonadaceae','s__uncultured_bacterium_g_Ellin6067','s__uncultured_bacterium_g_Ottowia',
                                                                   's__uncultured_bacterium_g_Geobacter','s__uncultured_bacterium_f_Pedosphaeraceae'),p=0.05)

#差异分析统计结果
stat_result$stat$anova
stat_result$stat$tukey

write.table(stat_result$stat$anova, 'species_anova.txt', row.names = FALSE, sep = '\t', quote = FALSE)
write.table(stat_result$stat$tukey, 'species_tukey.txt', row.names = FALSE, sep = '\t', quote = FALSE)


#显著性标记“abc”作图结果
stat_result$plot$plot_bar
stat_result$plot$plot_box


####根据log(x+1)对数转换后差异分析得到的数据，进行热图的绘制
library(pheatmap)

phylum_heatmap <- read.table('phylum_data.txt', header = TRUE)
phylum_heatmap

phylum_heatmap <- subset(phylum_heatmap,select= -Group)
phylum_heatmap

phylum_heatmap <- as.data.frame(phylum_heatmap)

rownames(phylum_heatmap) <- phylum_heatmap$Row.names
rownames(phylum_heatmap) 
phylum_heatmap

phylum_heatmap <- subset(phylum_heatmap,select= -Row.names)
phylum_heatmap

phylum_heatmap <-t(phylum_heatmap)
phylum_heatmap

pheatmap(phylum_heatmap)
pheatmap(phylum_heatmap,border_color=NA,cellheight=10,cellwidth=16,
         display_numbers=FALSE,
         cluster_col=F, show_rownames=T,show_colnames=T,
         fontsize_number=7,width=10,height=8,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50))

####添加行列注释信息
# 读取实验设计注释列分组
design<- read.delim('./mapping2.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
design
annot_data1 <- data.frame(row.names=colnames(phylum_heatmap), Treatment=design$Group)

#Treatment <- c("#CC6666","#D0DFE6FF","#999999","#0073c2","black")
#names(Treatment) <- c("Ggt","H","SA","S_135","BS")


pheatmap(phylum_heatmap,
         border_color=NA,
         cellheight=10,
         cellwidth=16,
         cluster_col=F, show_rownames=T,show_colnames=F,treeheight_row=15,treeheight_col=15,
         fontsize_number=7,
         width=10,height=8,
         frontsize=10,
         fontsize_row = 10,
         fontsize_col=10,
         annotation_col=annot_data1,
         angle_col = "270",
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50))

pheatmap(phylum_heatmap,
         border_color=NA,
         cellheight=10,
         cellwidth=16,
         cluster_col=F, show_rownames=T,show_colnames=F,treeheight_row=15,
         fontsize_number=7,
         width=10,height=8,
         frontsize=10,
         fontsize_row = 10,
         fontsize_col=10,
         annotation_col=annot_data1,
         angle_col = "270",
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
         filename="phylum_heatmap.pdf")

pheatmap(phylum_heatmap,
         border_color=NA,
         cellheight=10,
         cellwidth=16,
         cluster_col=F, show_rownames=T,show_colnames=F,treeheight_row=15,
         fontsize_number=7,
         width=10,height=8,
         frontsize=10,
         fontsize_row = 10,
         fontsize_col=10,
         annotation_col=annot_data1,
         angle_col = "270",
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
         filename="phylum_heatmap.png")

pheatmap(phylum_heatmap,
         border_color=NA,
         cellheight=10,
         cellwidth=16,
         cluster_col=F, show_rownames=T,show_colnames=F,treeheight_row=15,
         fontsize_number=7,
         width=10,height=8,
         frontsize=10,
         fontsize_row = 10,
         fontsize_col=10,
         annotation_col=annot_data1,
         angle_col = "270",
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
         filename="phylum_heatmap.tiff")


####纲水平class####
rm(list=ls())
setwd("E:/博士数据2/小麦微生物组数据/customer_backup/otus/ANOVA/")
library(pheatmap)


class_heatmap <- read.table('class_data.txt', header = TRUE)
class_heatmap

class_heatmap <- subset(class_heatmap,select= -Group)
class_heatmap

class_heatmap <- as.data.frame(class_heatmap)

rownames(class_heatmap) <- class_heatmap$Row.names
rownames(class_heatmap) 
class_heatmap

class_heatmap <- subset(class_heatmap,select= -Row.names)
class_heatmap

class_heatmap <-t(class_heatmap)
class_heatmap

pheatmap(class_heatmap)
pheatmap(class_heatmap,border_color="white",cellheight=15,cellwidth=15,
         display_numbers=FALSE,
         cluster_col=F, show_rownames=T,show_colnames=T,
         fontsize_number=7,width=10,height=8,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50))

####添加行列注释信息
# 读取实验设计注释列分组
design<- read.delim('./mapping2.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
design
annot_data2 <- data.frame(row.names=colnames(class_heatmap), group=design$Group)

#group <- c("#CC6666","#D0DFE6FF","#999999","#0073c2","black")
#names(group) <- c("Ggt","H","SA","S_135","BS")


pheatmap(class_heatmap,
         border_color=NA,
         cellheight=10,
         cellwidth=16,
         cluster_col=F, show_rownames=T,show_colnames=F,treeheight_row=15,
         fontsize_number=7,
         width=10,height=8,
         frontsize=10,
         fontsize_row = 10,
         fontsize_col=10,
         annotation_col=annot_data2,
         angle_col = "270")

pheatmap(class_heatmap,
         border_color=NA,
         cellheight=10,
         cellwidth=16,
         cluster_col=F, show_rownames=T,show_colnames=F,treeheight_row=15,
         fontsize_number=7,
         width=10,height=8,
         frontsize=10,
         fontsize_row = 10,
         fontsize_col=10,
         annotation_col=annot_data2,
         angle_col = "270",
         filename="class_heatmap.pdf")

pheatmap(class_heatmap,
         border_color=NA,
         cellheight=10,
         cellwidth=16,
         cluster_col=F, show_rownames=T,show_colnames=F,treeheight_row=15,
         fontsize_number=7,
         width=10,height=8,
         frontsize=10,
         fontsize_row = 10,
         fontsize_col=10,
         annotation_col=annot_data2,
         angle_col = "270",
         filename="class_heatmap.png")

pheatmap(class_heatmap,
         border_color=NA,
         cellheight=10,
         cellwidth=16,
         cluster_col=F, show_rownames=T,show_colnames=F,treeheight_row=15,
         fontsize_number=7,
         width=10,height=8,
         frontsize=10,
         fontsize_row = 10,
         fontsize_col=10,
         annotation_col=annot_data2,
         angle_col = "270",
         filename="class_heatmap.tiff")

####目水平order####
order_heatmap <- read.table('order_data.txt', header = TRUE)
order_heatmap

order_heatmap <- subset(order_heatmap,select= -Group)
order_heatmap

order_heatmap <- as.data.frame(order_heatmap)

rownames(order_heatmap) <- order_heatmap$Row.names
rownames(order_heatmap) 
order_heatmap

order_heatmap <- subset(order_heatmap,select= -Row.names)
order_heatmap

order_heatmap <-t(order_heatmap)
order_heatmap

pheatmap(order_heatmap)
pheatmap(order_heatmap,border_color="white",cellheight=15,cellwidth=15,
         display_numbers=FALSE,
         cluster_col=F, show_rownames=T,show_colnames=T,
         fontsize_number=7,width=10,height=8)

####添加行列注释信息
# 读取实验设计注释列分组
design<- read.delim('./mapping2.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
design
annot_data3 <- data.frame(row.names=colnames(order_heatmap), group=design$Group)

#group <- c("#CC6666","#D0DFE6FF","#999999","#0073c2","black")
#names(group) <- c("Ggt","H","SA","S_135","BS")


pheatmap(order_heatmap,
         border_color=NA,
         cellheight=10,
         cellwidth=16,
         cluster_col=F, show_rownames=T,show_colnames=F,treeheight_row=15,
         fontsize_number=7,
         width=10,height=8,
         frontsize=10,
         fontsize_row = 10,
         fontsize_col=10,
         annotation_col=annot_data3,
         angle_col = "270")

pheatmap(order_heatmap,
         border_color=NA,
         cellheight=10,
         cellwidth=16,
         cluster_col=F, show_rownames=T,show_colnames=F,treeheight_row=15,
         fontsize_number=7,
         width=10,height=8,
         frontsize=10,
         fontsize_row = 10,
         fontsize_col=10,
         annotation_col=annot_data3,
         angle_col = "270",
         filename="order_heatmap.pdf")

pheatmap(order_heatmap,
         border_color=NA,
         cellheight=10,
         cellwidth=16,
         cluster_col=F, show_rownames=T,show_colnames=F,treeheight_row=15,
         fontsize_number=7,
         width=10,height=8,
         frontsize=10,
         fontsize_row = 10,
         fontsize_col=10,
         annotation_col=annot_data3,
         angle_col = "270",
         filename="order_heatmap.png")

pheatmap(order_heatmap,
         border_color=NA,
         cellheight=10,
         cellwidth=16,
         cluster_col=F, show_rownames=T,show_colnames=F,treeheight_row=15,
         fontsize_number=7,
         width=10,height=8,
         frontsize=10,
         fontsize_row = 10,
         fontsize_col=10,
         annotation_col=annot_data3,
         angle_col = "270",
         filename="order_heatmap.tiff")


####科水平family####
family_heatmap <- read.table('family_data.txt', header = TRUE)
family_heatmap

family_heatmap <- subset(family_heatmap,select= -Group)
family_heatmap

family_heatmap <- as.data.frame(family_heatmap)

rownames(family_heatmap) <- family_heatmap$Row.names
rownames(family_heatmap) 
family_heatmap

family_heatmap <- subset(family_heatmap,select= -Row.names)
family_heatmap

family_heatmap <-t(family_heatmap)
family_heatmap


pheatmap(family_heatmap)
pheatmap(family_heatmap,border_color="white",cellheight=15,cellwidth=15,
         display_numbers=FALSE,
         cluster_col=F, show_rownames=T,show_colnames=T,
         fontsize_number=7,width=10,height=8)

####添加行列注释信息
# 读取实验设计注释列分组
design<- read.delim('./mapping.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
design
annot_data4 <- data.frame(row.names=colnames(family_heatmap), group=design$Group)

#group <- c("#CC6666","#D0DFE6FF","#999999","#0073c2","black")
#names(group) <- c("Ggt","H","SA","S_135","BS")


pheatmap(family_heatmap,
         border_color=NA,
         cellheight=10,
         cellwidth=16,
         cluster_col=F, show_rownames=T,show_colnames=F,treeheight_row=15,
         fontsize_number=7,
         width=10,height=8,
         frontsize=10,
         fontsize_row = 10,
         fontsize_col=10,
         annotation_col=annot_data4,
         angle_col = "270")

pheatmap(family_heatmap,
         border_color=NA,
         cellheight=10,
         cellwidth=16,
         cluster_col=F, show_rownames=T,show_colnames=F,treeheight_row=15,
         fontsize_number=7,
         width=10,height=8,
         frontsize=10,
         fontsize_row = 10,
         fontsize_col=10,
         annotation_col=annot_data4,
         angle_col = "270",
         filename="family_heatmap.pdf")

pheatmap(family_heatmap,
         border_color=NA,
         cellheight=10,
         cellwidth=16,
         cluster_col=F, show_rownames=T,show_colnames=F,treeheight_row=15,
         fontsize_number=7,
         width=10,height=8,
         frontsize=10,
         fontsize_row = 10,
         fontsize_col=10,
         annotation_col=annot_data4,
         angle_col = "270",
         filename="family_heatmap.png")

pheatmap(family_heatmap,
         border_color=NA,
         cellheight=10,
         cellwidth=16,
         cluster_col=F, show_rownames=T,show_colnames=F,treeheight_row=15,
         fontsize_number=7,
         width=10,height=8,
         frontsize=10,
         fontsize_row = 10,
         fontsize_col=10,
         annotation_col=annot_data4,
         angle_col = "270",
         filename="family_heatmap.tiff")


####属水平genus####
genus_heatmap <- read.table('genus_data.txt', header = TRUE)
genus_heatmap

genus_heatmap <- subset(genus_heatmap,select= -Group)
genus_heatmap

genus_heatmap <- as.data.frame(genus_heatmap)

rownames(genus_heatmap) <- genus_heatmap$Row.names
rownames(genus_heatmap) 
genus_heatmap

genus_heatmap <- subset(genus_heatmap,select= -Row.names)
genus_heatmap

genus_heatmap <-t(genus_heatmap)
genus_heatmap


pheatmap(genus_heatmap)
pheatmap(genus_heatmap,border_color="white",cellheight=15,cellwidth=15,
         display_numbers=FALSE,
         cluster_col=F, show_rownames=T,show_colnames=T,
         fontsize_number=7,width=10,height=8)

####添加行列注释信息
# 读取实验设计注释列分组
design<- read.delim('./mapping2.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
design
annot_data5 <- data.frame(row.names=colnames(genus_heatmap), group=design$Group)

#group <- c("#CC6666","#D0DFE6FF","#999999","#0073c2","black")
#names(group) <- c("Ggt","H","SA","S_135","BS")


pheatmap(genus_heatmap,
         border_color=NA,
         cellheight=10,
         cellwidth=16,
         cluster_col=F, show_rownames=T,show_colnames=F,treeheight_row=15,
         fontsize_number=7,
         width=10,height=8,
         frontsize=10,
         fontsize_row = 10,
         fontsize_col=10,
         annotation_col=annot_data5,
         angle_col = "270")

pheatmap(genus_heatmap,
         border_color=NA,
         cellheight=10,
         cellwidth=16,
         cluster_col=F, show_rownames=T,show_colnames=F,treeheight_row=15,
         fontsize_number=7,
         width=10,height=8,
         frontsize=10,
         fontsize_row = 10,
         fontsize_col=10,
         annotation_col=annot_data5,
         angle_col = "270",
         filename="genus_heatmap.pdf")

pheatmap(genus_heatmap,
         border_color=NA,
         cellheight=10,
         cellwidth=16,
         cluster_col=F, show_rownames=T,show_colnames=F,treeheight_row=15,
         fontsize_number=7,
         width=10,height=8,
         frontsize=10,
         fontsize_row = 10,
         fontsize_col=10,
         annotation_col=annot_data5,
         angle_col = "270",
         filename="genus_heatmap.png")

pheatmap(genus_heatmap,
         border_color=NA,
         cellheight=10,
         cellwidth=16,
         cluster_col=F, show_rownames=T,show_colnames=F,treeheight_row=15,
         fontsize_number=7,
         width=10,height=8,
         frontsize=10,
         fontsize_row = 10,
         fontsize_col=10,
         annotation_col=annot_data5,
         angle_col = "270",
         filename="genus_heatmap.tiff")

####高丰度物种差异分析####

##基于丰富物种的物种差异分析（丰富物种占比大于千分之一，即0.1%）

####丰度物种划分####

setwd("E:/博士数据2/小麦微生物组数据/customer_backup/otus/difference/")
#首先导入数据文件
otu <- read.delim('OTU_Tables.txt', row.names = 1)
otu
#稀有物种（RT），在所有样本中相对丰度都低于0.01%
otu_ART <- otu[apply(otu, 1, function(x) max(x)<0.0001), ]
otu_ART
write.table(otu_ART, 'otu_ART.txt', col.names = NA, sep = '\t', quote = FALSE)

#丰富物种（AT），在所有样本中相对丰度都高于1%
otu_AAT <- otu[apply(otu, 1, function(x) min(x)>0.01), ]
otu_AAT
write.table(otu_AAT, 'otu_AAT.txt', col.names = NA, sep = '\t', quote = FALSE)


#中等丰度物种（MT），在所有样本中相对丰度介于0.01%和1%之间
otu_MT <- otu[apply(otu, 1, function(x) max(x)<0.01 & min(x)>0.0001), ]
otu_MT
write.table(otu_MT, 'otu_MT.txt', col.names = NA, sep = '\t', quote = FALSE)

#条件稀有物种（CRT），仅在一些样本中相对丰度低于0.01%，但从不高于1%
otu_CRT <- otu[apply(otu, 1, function(x) max(x)<=0.01 & min(x)<0.0001), ]
otu_CRT
write.table(otu_CRT, 'otu_CRT.txt', col.names = NA, sep = '\t', quote = FALSE)

#条件丰富物种（CAT），仅在一些样本中相对丰度高于1%，但从不低于0.01%
otu_CAT <- otu[apply(otu, 1, function(x) min(x)>=0.0001 & max(x)>0.01), ]
otu_CAT <- otu_CAT[!rownames(otu_CAT) %in% rownames(otu_AAT),]
write.table(otu_CAT, 'otu_CAT.txt', col.names = NA, sep = '\t', quote = FALSE)

#条件稀有或丰富物种（CRAT），最低丰度低于0.01，最高丰度高于1%
otu_CRAT <- otu[apply(otu, 1, function(x) max(x)>0.01 & min(x)<0.0001), ]
otu_CRT <- otu_CRT[!rownames(otu_CRT) %in% rownames(otu_ART),]
write.table(otu_CRAT, 'otu_CRAT.txt', col.names = NA, sep = '\t', quote = FALSE)

tax <- read.delim('./taxonomy2.txt',row.names = 1)
head(tax, n = 3)
#删除OTU丰度表中，有关p_的字样
phylum <- gsub(".*(p__)","",tax$phylum)
phylum
phylum <- gsub("(\\[)","",phylum)
phylum <- gsub("(\\])","",phylum)

class <- gsub(".*(c__)","",tax$class)
class
class <- gsub("(\\[)","",class)
class <- gsub("(\\])","",class)

order <- gsub(".*(o__)","",tax$order)
order
order <- gsub("(\\[)","",order)
order <- gsub("(\\])","",order)

family <- gsub(".*(f__)","",tax$family)
family
family <- gsub("(\\[)","",family)
family <- gsub("(\\])","",family)

genus <- gsub(".*(g__)","",tax$genus)
genus
genus <- gsub("(\\[)","",genus)
genus <- gsub("(\\])","",genus)

species <- gsub(".*(s__)","",tax$species)
species
species <- gsub("(\\[)","",species)
species <- gsub("(\\])","",species)


#将新产生的表头与原来的表头进行替换
tax$phylum <- phylum
tax$class <- class
tax$order <- order
tax$family <- family
tax$genus <- genus
tax$species <- species
tax


####基于丰富物种丰度的分析AT

dat <- merge(x=otu_AAT,y=tax,by="row.names")
dat
dat =dplyr::rename(dat,OTUID = Row.names)
dat
write.table(dat, 'otu_AAT_tax.txt', row.names = FALSE, sep = '\t', quote = FALSE)


otu <- read.delim('otu_AT.txt', row.names = 1)
otu
##排序，便于筛选TOP物种
otu <- as.matrix(otu)
otu1 <- matrix(as.numeric(otu),nrow = nrow(otu))
rownames(otu1) <- rownames(otu)
colnames(otu1) <- colnames(otu)
sum <- apply(otu1,1,sum) 
otu1 <- cbind(otu1,sum)
otu1 <- as.data.frame(otu1)
otu1 <- otu1[order(otu1[,"sum"],decreasing = T),]
write.table(otu1, file = "otu_AT_order.txt",sep = '\t',)


##按分类水平分组汇总(根据自己需求更改要展示的物种水平)
phylum <-aggregate(dat[,2:13],by=list(dat$phylum),FUN=sum)
phylum
write.table(phylum, 'phylum_AT.txt', row.names = FALSE, sep = '\t', quote = FALSE)

phylum <- read.delim('phylum_AT.txt', row.names = 1)
phylum
##排序，便于筛选TOP物种
phylum <- as.matrix(phylum)
phylum1 <- matrix(as.numeric(phylum),nrow = nrow(phylum))
rownames(phylum1) <- rownames(phylum)
colnames(phylum1) <- colnames(phylum)
sum <- apply(phylum1,1,sum) 
phylum1 <- cbind(phylum1,sum)
phylum1 <- as.data.frame(phylum1)
phylum1 <- phylum1[order(phylum1[,"sum"],decreasing = T),]
write.table(phylum1, file = "phylum_AT_order.txt",sep = '\t',)

#如果需要增减物种的数目，请调整此处的数值
#phylum2 <- phylum1[1:20,]
#phylum2 <- t(phylum2)
#sum <- apply(phylum2,1,sum) 
#Others <- 1-sum
#phylum20 <- cbind(phylum2,Others)
#phylum20
#write.table(phylum20, file = "phylum_AT_Top20.txt",sep = '\t',)



class <-aggregate(dat[,2:13],by=list(dat$class),FUN=sum)
class
write.table(class, 'class_AT.txt', row.names = FALSE, sep = '\t', quote = FALSE)

class <- read.delim('class_AT.txt', row.names = 1)
class
##排序，便于筛选TOP物种
class <- as.matrix(class)
class1 <- matrix(as.numeric(class),nrow = nrow(class))
rownames(class1) <- rownames(class)
colnames(class1) <- colnames(class)
sum <- apply(class1,1,sum) 
class1 <- cbind(class1,sum)
class1 <- as.data.frame(class1)
class1 <- class1[order(class1[,"sum"],decreasing = T),]
write.table(class1, file = "class_AT_order.txt",sep = '\t',)

#如果需要增减物种的数目，请调整此处的数值
class2 <- subset(class1,select = -sum)
class2
class20 <- class2[1:20,]
class20 <- t(class20)
class20
write.table(class20, file = "class_AT_Top20.txt",sep = '\t',)


order <-aggregate(dat[,2:13],by=list(dat$order),FUN=sum)
order
write.table(order, 'order_AT.txt', row.names = FALSE, sep = '\t', quote = FALSE)

order <- read.delim('order_AT.txt', row.names = 1)
order
##排序，便于筛选TOP物种
order <- as.matrix(order)
order1 <- matrix(as.numeric(order),nrow = nrow(order))
rownames(order1) <- rownames(order)
colnames(order1) <- colnames(order)
sum <- apply(order1,1,sum) 
order1 <- cbind(order1,sum)
order1 <- as.data.frame(order1)
order1 <- order1[order(order1[,"sum"],decreasing = T),]
write.table(order1, file = "order_AT_order.txt",sep = '\t',)

#如果需要增减物种的数目，请调整此处的数值
#如果需要增减物种的数目，请调整此处的数值
order2 <- subset(order1,select = -sum)
order2
order20 <- order2[1:20,]
order20 <- t(order20)
order20
write.table(order20, file = "order_AT_Top20.txt",sep = '\t',)



family <-aggregate(dat[,2:13],by=list(dat$family),FUN=sum)
family
write.table(family, 'family_AT.txt', row.names = FALSE, sep = '\t', quote = FALSE)

family <- read.delim('family_AT.txt', row.names = 1)
family
##排序，便于筛选TOP物种
family <- as.matrix(family)
family1 <- matrix(as.numeric(family),nrow = nrow(family))
rownames(family1) <- rownames(family)
colnames(family1) <- colnames(family)
sum <- apply(family1,1,sum) 
family1 <- cbind(family1,sum)
family1 <- as.data.frame(family1)
family1 <- family1[order(family1[,"sum"],decreasing = T),]
write.table(family1, file = "family_AT_order.txt",sep = '\t',)

#如果需要增减物种的数目，请调整此处的数值
family2 <- subset(family1,select = -sum)
family2
family20 <- family2[1:20,]
family20 <- t(family20)
family20
write.table(family20, file = "family_AT_Top20.txt",sep = '\t',)




genus <-aggregate(dat[,2:13],by=list(dat$genus),FUN=sum)
genus
write.table(genus, 'genus_AT.txt', row.names = FALSE, sep = '\t', quote = FALSE)

genus <- read.delim('genus_AT.txt', row.names = 1)
genus
##排序，便于筛选TOP物种
genus <- as.matrix(genus)
genus1 <- matrix(as.numeric(genus),nrow = nrow(genus))
rownames(genus1) <- rownames(genus)
colnames(genus1) <- colnames(genus)
sum <- apply(genus1,1,sum) 
genus1 <- cbind(genus1,sum)
genus1 <- as.data.frame(genus1)
genus1 <- genus1[order(genus1[,"sum"],decreasing = T),]
write.table(genus1, file = "genus_AT_order.txt",sep = '\t',)

#如果需要增减物种的数目，请调整此处的数值
genus2 <- subset(genus1,select = -sum)
genus2
genus20 <-genus2[1:20,]
genus20 <- t(genus20)
genus20
write.table(genus20, file = "genus_AT_Top20.txt",sep = '\t',)



##采用scale_test对数据进行对数转换
phylum_AT2 <- read.delim('./phylum_AT.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
phylum_AT2
##首先进行对数转换
phylum_AT2 <- apply(phylum_AT2, 2, function(x) log2(x+1))
phylum_AT2
##由于仅仅利用对数转换后的数值绘制出的热图颜色区分度不大，故对对数转换后的数值进一步求平方根
phylum_AT2 <- apply(phylum_AT2, 2, function(x) sqrt(x))
phylum_AT2

write.table(phylum_AT2, 'phylum_AT_scale_test.txt', row.names = TRUE, sep = '\t', quote = FALSE)


class_AT2 <- read.delim('./class_AT_Top20.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
class_AT2
##首先进行对数转换
class_AT2 <- apply(class_AT2, 2, function(x) log2(x+1))
class_AT2
##由于仅仅利用对数转换后的数值绘制出的热图颜色区分度不大，故对对数转换后的数值进一步求平方根
class_AT2 <- apply(class_AT2, 2, function(x) sqrt(x))
class_AT2

write.table(class_AT2, 'class_AT_scale_test.txt', row.names = TRUE, sep = '\t', quote = FALSE)


order_AT2 <- read.delim('./order_AT_Top20.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
order_AT2
##首先进行对数转换
order_AT2 <- apply(order_AT2, 2, function(x) log2(x+1))
order_AT2
##由于仅仅利用对数转换后的数值绘制出的热图颜色区分度不大，故对对数转换后的数值进一步求平方根
order_AT2 <- apply(order_AT2, 2, function(x) sqrt(x))
order_AT2

write.table(order_AT2, 'order_AT_scale_test.txt', row.names = TRUE, sep = '\t', quote = FALSE)


family_AT2 <- read.delim('./family_AT_Top20.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
family_AT2
##首先进行对数转换
family_AT2 <- apply(family_AT2, 2, function(x) log2(x+1))
family_AT2
##由于仅仅利用对数转换后的数值绘制出的热图颜色区分度不大，故对对数转换后的数值进一步求平方根
family_AT2 <- apply(family_AT2, 2, function(x) sqrt(x))
family_AT2

write.table(family_AT2, 'family_AT_scale_test.txt', row.names = TRUE, sep = '\t', quote = FALSE)


genus_AT2 <- read.delim('./genus_AT_Top20.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
genus_AT2
##首先进行对数转换
genus_AT2 <- apply(genus_AT2, 2, function(x) log2(x+1))
genus_AT2
##由于仅仅利用对数转换后的数值绘制出的热图颜色区分度不大，故对对数转换后的数值进一步求平方根
genus_AT2 <- apply(genus_AT2, 2, function(x) sqrt(x))
genus_AT2

write.table(genus_AT2, 'genus_AT_scale_test.txt', row.names = TRUE, sep = '\t', quote = FALSE)


otu_AT2 <- read.delim('./otu_AAT.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
otu_AT2
##首先进行对数转换
otu_AT2 <- apply(otu_AT2, 2, function(x) log2(x+1))
otu_AT2
##由于仅仅利用对数转换后的数值绘制出的热图颜色区分度不大，故对对数转换后的数值进一步求平方根
otu_AT2 <- apply(otu_AT2, 2, function(x) sqrt(x))
otu_AT2

write.table(otu_AT2, 'otu_AT_scale_test.txt', row.names = TRUE, sep = '\t', quote = FALSE)



##数据转置，为方差分析做准备
metadata<- read.delim('./mapping.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
metadata


phylum_AT2 <- read.delim('./phylum_AT_scale_test.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
phylum_AT2
phylum_AT3 <-t(phylum_AT2)
phylum_AT3

data1 <- merge(x=phylum_AT3,y=metadata,by='row.names')
data1
write.table(data1, 'phylun_AT_data.txt', row.names = FALSE, sep = '\t', quote = FALSE)


class_AT2 <- read.delim('./class_AT_scale_test.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
class_AT2
data2 <- merge(x=class_AT2,y=metadata,by='row.names')
data2
write.table(data2, 'class_AT_data.txt', row.names = FALSE, sep = '\t', quote = FALSE)



order_AT2 <- read.delim('./order_AT_scale_test.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
order_AT2

data3 <- merge(x=order_AT2,y=metadata,by='row.names')
data3
write.table(data3, 'order_AT_data.txt', row.names = FALSE, sep = '\t', quote = FALSE)



family_AT2 <- read.delim('./family_AT_scale_test.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
family_AT2

data4 <- merge(x=family_AT2,y=metadata,by='row.names')
data4
write.table(data4, 'family_AT_data.txt', row.names = FALSE, sep = '\t', quote = FALSE)


genus_AT2 <- read.delim('./genus_AT_scale_test.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
genus_AT2
data5 <- merge(x=genus_AT2,y=metadata,by='row.names')
data5
write.table(data5, 'genus_AT_data.txt', row.names = FALSE, sep = '\t', quote = FALSE)


otu_AT2 <- read.delim('./otu_AT_scale_test.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
otu_AT2

otu_AT3 <-t(otu_AT2)
otu_AT3
data6 <- merge(x=otu_AT3,y=metadata,by='row.names')
data6
write.table(data6, 'otu_AT_data.txt', row.names = FALSE, sep = '\t', quote = FALSE)






#生信小白鱼方差分析source
#构建函数，统计各组均值、标准差，执行 ANOVA 以及 Tukey HSD 获得显著性
library(reshape2)
library(multcomp)
library(ggplot2)

aov_tukey <- function(data, groups, values, p = 0.05) {
  stat_anova <- NULL
  abc_list <- NULL
  
  #统计检验
  for (i in groups) {
    for (j in values) {
      
      #单因素 ANOVA，整体差异
      dat <- data[c(i, j)]
      names(dat) <- c('class', 'var')
      dat$class <- factor(dat$class)
      fit <- aov(var~class, dat)
      p_value <- summary(fit)[[1]][1,5]
      
      #单因素 ANOVA 结果整理
      if (p_value < 0.001) sig <- '***'
      else if (p_value >= 0.001 & p_value < 0.01) sig <- '**'
      else if (p_value >= 0.01 & p_value < 0.05) sig <- '*'
      else sig <- ''
      stat_anova <- rbind(stat_anova, c(paste(i, j, sep = '/'), p_value, sig))
      
      #Tukey HSD 检验（multcomp 包），多重比较
      tuk <- cld(glht(fit, alternative = 'two.sided', linfct = mcp(class = 'Tukey')), sig = p, decreasing = TRUE)
      
      #cld() 自动得到了显著性“abc”水平，提取显著性标记“abc”
      sig <- data.frame(tuk$mcletters$Letters, stringsAsFactors = FALSE)
      names(sig) <- 'sig'
      sig$class <- rownames(sig)
      sig$var <- j
      
      #均值标准差统计
      abc_ij <- cbind(aggregate(dat$var, by = list(dat$class), FUN = mean), aggregate(dat$var, by = list(dat$class), FUN = sd)[2])
      names(abc_ij) <- c('class', 'mean', 'sd')
      abc_ij$class <- as.character(abc_ij$class)
      
      #合并结果
      abc_ij$group <- i
      abc_ij <- merge(abc_ij, sig, by = 'class')
      abc_ij <- abc_ij[c(4, 6, 1, 2, 3, 5)]
      abc_list <- rbind(abc_list, abc_ij)
    }
  }
  
  #ggplot2 作图，柱状图
  plot_bar <- ggplot(data = abc_list, aes(x = class, y = mean)) +
    geom_col(aes(fill = class), color = NA, show.legend = FALSE) +
    geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.2) +
    facet_grid(var~group, scale = 'free' , space = 'free_x') +
    geom_text(aes(label = sig, y = mean + sd + 0.3*(mean+sd))) +
    labs(x = '', y = '')
  
  #ggplot2 作图，箱线图
  data <- melt(data[c(groups, values)], id = groups)
  data <- melt(data, id = c('variable', 'value'))
  data <- data[c(3, 1, 4, 2)]
  names(data) <- c('group', 'var', 'class', 'value')
  
  value_max <- aggregate(data$value, by = list(data$group, data$var), FUN = max)
  names(value_max) <- c('group', 'var', 'value')
  value_max <- merge(value_max, abc_list[c(-4, -5)], by = c('group', 'var'))
  
  plot_box <- ggplot(data = data, aes(x = class, y = value)) +
    geom_boxplot(aes(fill = class), show.legend = FALSE) +
    facet_grid(var~group, scale = 'free' , space = 'free_x') +
    geom_text(data = value_max, aes(label = sig, y = value + 0.3*value)) +
    labs(x = '', y = '')
  
  #return
  stat_anova <- data.frame(stat_anova, stringsAsFactors = FALSE)
  names(stat_anova) <- c('group', 'anova_pvalue', 'sig')
  stat_anova$anova_pvalue <- as.numeric(stat_anova$anova_pvalue)
  
  list(stat = list(anova = stat_anova, tukey = abc_list), plot = list(plot_bar = plot_bar, plot_box = plot_box))
}

####phulum门水平####
##调用函数运行
data1 <- read.table('phylun_AT_data.txt', header = TRUE)
data1

##将微生物的分类名称提出来进行方差分析时使用
value1 <- colnames(data1)
value1 <- t(value1)
value1
value11 <- paste(value1,collapse = "','")
value11

#groups 指定分组列，可指定多列；values 指定变量列，可指定多列；p 指定显著性水平
#这里作为测试，就先假设这些数据满足方差分析的前提条件了
#复制粘贴后记得换行，不然会报错
stat_result <- aov_tukey(data = data1, groups = 'Group', values = c('Acidobacteria','Actinobacteria','Armatimonadetes','Bacteroidetes','Chloroflexi',
                                                                    'Cyanobacteria','Epsilonbacteraeota','Fibrobacteres','Firmicutes','Fusobacteria','Gemmatimonadetes','Nitrospirae','Patescibacteria','Planctomycetes','Proteobacteria','Rokubacteria','Verrucomicrobia'), p = 0.05)

#差异分析统计结果
stat_result$stat$anova
stat_result$stat$tukey

write.table(stat_result$stat$anova, 'phylum_AT_anova.txt', row.names = FALSE, sep = '\t', quote = FALSE)
write.table(stat_result$stat$tukey, 'phylum_AT_tukey.txt', row.names = FALSE, sep = '\t', quote = FALSE)



#显著性标记“abc”作图结果
stat_result$plot$plot_bar
stat_result$plot$plot_box


####class纲水平####
##调用函数运行
data2 <- read.table('class_AT_data.txt', header = TRUE)
data2

##将微生物的分类名称提出来进行方差分析时使用
value2 <- colnames(data2)
value2 <- t(value2)
value2
value22 <- paste(value2,collapse = "','")
value22

#groups 指定分组列，可指定多列；values 指定变量列，可指定多列；p 指定显著性水平
#这里作为测试，就先假设这些数据满足方差分析的前提条件了
stat_result <- aov_tukey(data = data2, groups = 'Group', values = c('Bacteroidia','Bacilli','Gammaproteobacteria','Alphaproteobacteria','Subgroup_6','Anaerolineae','Acidobacteriia','Blastocatellia_Subgroup_4',
                                                                    'Deltaproteobacteria','Acidimicrobiia','Actinobacteria','Fusobacteriia','Gemmatimonadetes','Oxyphotobacteria','Negativicutes','KD4.96',
                                                                    'Thermoleophilia','Clostridia','Nitrospira','Holophagae'), p = 0.05)

#差异分析统计结果
stat_result$stat$anova
stat_result$stat$tukey

write.table(stat_result$stat$anova, 'class_AT_anova.txt', row.names = FALSE, sep = '\t', quote = FALSE)
write.table(stat_result$stat$tukey, 'class_AT_tukey.txt', row.names = FALSE, sep = '\t', quote = FALSE)

#显著性标记“abc”作图结果
stat_result$plot$plot_bar
stat_result$plot$plot_box



####order目水平####
##调用函数运行
data3 <- read.table('order_AT_data.txt', header = TRUE)
data3

##将微生物的分类名称提出来进行方差分析时使用
value3 <- colnames(data3)
value3 <- t(value3)
value3
value33 <- paste(value3,collapse = "','")
value33

#groups 指定分组列，可指定多列；values 指定变量列，可指定多列；p 指定显著性水平
#这里作为测试，就先假设这些数据满足方差分析的前提条件了
stat_result <- aov_tukey(data = data3, groups = 'Group', values = c('Bacteroidales','Lactobacillales','Betaproteobacteriales','uncultured_bacterium_c_Subgroup_6',
                                                                    'Rhizobiales','Bacillales','Solibacterales','SBR1031','Sphingomonadales','Blastocatellales',
                                                                    'Pseudomonadales','Fusobacteriales','Enterobacteriales','Gemmatimonadales','Selenomonadales',
                                                                    'Microtrichales','uncultured_bacterium_c_KD4.96','Azospirillales','Cytophagales','Myxococcales'), p = 0.05)

#差异分析统计结果
stat_result$stat$anova
stat_result$stat$tukey

write.table(stat_result$stat$anova, 'order_AT_anova.txt', row.names = FALSE, sep = '\t', quote = FALSE)
write.table(stat_result$stat$tukey, 'order_AT_tukey.txt', row.names = FALSE, sep = '\t', quote = FALSE)

#显著性标记“abc”作图结果
stat_result$plot$plot_bar
stat_result$plot$plot_box




####family目水平####
##调用函数运行
data4 <- read.table('family_AT_data.txt', header = TRUE)
data4

##将微生物的分类名称提出来进行方差分析时使用
value4 <- colnames(data4)
value4 <- t(value4)
value4
value44 <- paste(value4,collapse = "','")
value44


#groups 指定分组列，可指定多列；values 指定变量列，可指定多列；p 指定显著性水平
#这里作为测试，就先假设这些数据满足方差分析的前提条件了
stat_result <- aov_tukey(data = data4, groups = 'Group', values = c('Prevotellaceae','Streptococcaceae','uncultured_bacterium_c_Subgroup_6','Rhodocyclaceae','Burkholderiaceae',
                                                                    'Nitrosomonadaceae','Solibacteraceae_Subgroup_3','Bacillaceae','Sphingomonadaceae','Blastocatellaceae',
                                                                    'Pseudomonadaceae','Family_XI','A4b','Xanthobacteraceae','Fusobacteriaceae','Enterobacteriaceae',
                                                                    'Gemmatimonadaceae','Veillonellaceae','uncultured_bacterium_c_KD4.96','Rhizobiaceae'), p = 0.05)

#差异分析统计结果
stat_result$stat$anova
stat_result$stat$tukey

write.table(stat_result$stat$anova, 'family_AT_anova.txt', row.names = FALSE, sep = '\t', quote = FALSE)
write.table(stat_result$stat$tukey, 'family_AT_tukey.txt', row.names = FALSE, sep = '\t', quote = FALSE)


#显著性标记“abc”作图结果
stat_result$plot$plot_bar
stat_result$plot$plot_box


####genus目水平####
##调用函数运行
data5 <- read.delim("genus_AT_data.txt", header = TRUE, sep = '\t', stringsAsFactors = FALSE)
data5

##将微生物的分类名称提出来进行方差分析时使用
value5 <- colnames(data5)
value5 <- t(value5)
value5
value55 <- paste(value5,collapse = "','")
value55

#groups 指定分组列，可指定多列；values 指定变量列，可指定多列；p 指定显著性水平
#这里作为测试，就先假设这些数据满足方差分析的前提条件了
stat_result <- aov_tukey(data = data5, groups = 'Group',values = c('Alloprevotella','Streptococcus','uncultured_bacterium_c_Subgroup_6','Azoarcus','Bacillus',
                                                                   'Dechloromonas','Bryobacter','MND1','Gemella','uncultured_bacterium_f_A4b','Fusobacterium',
                                                                   'Sphingomonas','uncultured_bacterium_f_Gemmatimonadaceae','Veillonella','Escherichia.Shigella',
                                                                   'uncultured_bacterium_c_KD4.96','Aridibacter','uncultured_bacterium_o_SBR1031','Pseudomonas',
                                                                   'uncultured_bacterium_f_Microscillaceae'), p = 0.05)

#差异分析统计结果
stat_result$stat$anova
stat_result$stat$tukey

write.table(stat_result$stat$anova, 'genus_AT_anova.txt', row.names = FALSE, sep = '\t', quote = FALSE)
write.table(stat_result$stat$tukey, 'genus_AT_tukey.txt', row.names = FALSE, sep = '\t', quote = FALSE)


#显著性标记“abc”作图结果
stat_result$plot$plot_bar
stat_result$plot$plot_box


####otu水平####
##调用函数运行
data6 <- read.delim("otu_AT_data.txt", header = TRUE, sep = '\t', stringsAsFactors = FALSE)
data6

##将微生物的分类名称提出来进行方差分析时使用
value6 <- colnames(data6)
value6 <- t(value6)
value6
value66 <- paste(value6,collapse = "','")
value66

#groups 指定分组列，可指定多列；values 指定变量列，可指定多列；p 指定显著性水平
#这里作为测试，就先假设这些数据满足方差分析的前提条件了
stat_result <- aov_tukey(data = data6, groups = 'Group',values = c('OTU1','OTU10','OTU100','OTU10032','OTU10057','OTU10211','OTU10255','OTU1026','OTU10283','OTU103','OTU10310','OTU104','OTU105','OTU106','OTU1062','OTU1066','OTU107','OTU10746','OTU1078','OTU1088','OTU109','OTU10916','OTU10935','OTU10977','OTU110','OTU1104','OTU111','OTU1113','OTU112','OTU11296','OTU113',
                                                                   'OTU1132','OTU1142','OTU1147','OTU11474','OTU116','OTU1169','OTU117','OTU118','OTU119','OTU12','OTU121','OTU122','OTU123','OTU12334','OTU1238','OTU125','OTU1250','OTU12517','OTU12545','OTU126','OTU127','OTU13','OTU130','OTU132','OTU133','OTU134','OTU1342','OTU135','OTU136','OTU137','OTU1382','OTU1387',
                                                                   'OTU1390','OTU1392','OTU140','OTU141','OTU1427','OTU145','OTU146','OTU1468','OTU1471','OTU1479','OTU148','OTU1481','OTU149','OTU15','OTU151','OTU1512','OTU153','OTU154','OTU155','OTU156','OTU157','OTU158','OTU163','OTU1632','OTU164','OTU165','OTU1658','OTU167','OTU1670','OTU168','OTU17','OTU171','OTU173',
                                                                   'OTU178','OTU1789','OTU179','OTU180','OTU1819','OTU1824','OTU183','OTU1838','OTU184','OTU186','OTU188','OTU19','OTU190','OTU191','OTU192','OTU196','OTU1975','OTU1990','OTU2','OTU200','OTU201','OTU2014','OTU2015','OTU2020','OTU203','OTU204','OTU2044','OTU2045','OTU206','OTU2069','OTU207','OTU2077','OTU2086',
                                                                   'OTU2089','OTU209','OTU21','OTU210','OTU212','OTU213','OTU2135','OTU214','OTU2151','OTU2156','OTU2185','OTU2188','OTU219','OTU2203','OTU222','OTU2225','OTU226','OTU23','OTU230','OTU238','OTU2383','OTU239','OTU241','OTU242','OTU245','OTU2466','OTU247','OTU248','OTU2487','OTU251','OTU254','OTU256','OTU258',
                                                                   'OTU259','OTU2599','OTU26','OTU260','OTU263','OTU2632','OTU264','OTU265','OTU266','OTU2668','OTU267','OTU268','OTU269','OTU270','OTU2729','OTU273','OTU2732','OTU274','OTU275','OTU2771','OTU278','OTU28','OTU2849','OTU285','OTU288','OTU2894','OTU290','OTU2933','OTU295','OTU297','OTU30','OTU300','OTU301',
                                                                   'OTU304','OTU3060','OTU308','OTU3166','OTU318','OTU3180','OTU3195','OTU32','OTU320','OTU323','OTU325','OTU327','OTU329','OTU33','OTU331','OTU333','OTU336','OTU337','OTU338','OTU34','OTU340','OTU342','OTU344','OTU345','OTU3457','OTU349','OTU35','OTU352','OTU354','OTU355','OTU3555','OTU357','OTU36',
                                                                   'OTU361','OTU364','OTU367','OTU37','OTU372','OTU374','OTU375','OTU376','OTU38','OTU380','OTU382','OTU385','OTU386','OTU3872','OTU388','OTU390','OTU395','OTU3953','OTU397','OTU398','OTU399','OTU40','OTU402','OTU405','OTU408','OTU41','OTU4109','OTU4135','OTU415','OTU418','OTU419','OTU42','OTU422',
                                                                   'OTU426','OTU428','OTU43','OTU431','OTU432','OTU434','OTU44','OTU446','OTU45','OTU450','OTU4509','OTU4519','OTU46','OTU460','OTU462','OTU4721','OTU481','OTU486','OTU489','OTU49','OTU494','OTU495','OTU497','OTU501','OTU504','OTU506','OTU508','OTU52','OTU523','OTU527','OTU528','OTU53','OTU5382',
                                                                   'OTU544','OTU5472','OTU55','OTU550','OTU551','OTU553','OTU5539','OTU56','OTU565','OTU567','OTU577','OTU578','OTU58','OTU580','OTU5877','OTU59','OTU593','OTU594','OTU597','OTU60','OTU607','OTU609','OTU61','OTU613','OTU6158','OTU62','OTU622','OTU627','OTU6284','OTU63','OTU631','OTU636','OTU638',
                                                                   'OTU64','OTU642','OTU643','OTU6482','OTU65','OTU6564','OTU658','OTU661','OTU662','OTU6632','OTU6675','OTU668','OTU67','OTU680','OTU682','OTU6884','OTU690','OTU694','OTU6982','OTU70','OTU708','OTU709','OTU71','OTU712','OTU72','OTU724','OTU73','OTU732','OTU733','OTU74','OTU7414','OTU745','OTU7494',
                                                                   'OTU7498','OTU751','OTU754','OTU7548','OTU756','OTU758','OTU7617','OTU7633','OTU7647','OTU770','OTU781','OTU7826','OTU789','OTU793','OTU798','OTU800','OTU808','OTU811','OTU82','OTU830','OTU833','OTU8344','OTU8362','OTU839','OTU84','OTU8425','OTU845','OTU848','OTU85','OTU852','OTU86','OTU866',
                                                                   'OTU87','OTU881','OTU882','OTU8821','OTU9','OTU90','OTU900','OTU9025','OTU908','OTU91','OTU9101','OTU915','OTU92','OTU93','OTU941','OTU9452','OTU950','OTU958','OTU9618','OTU97','OTU972','OTU9727','OTU98','OTU99','OTU992','OTU9999'), p = 0.05)

#差异分析统计结果
stat_result$stat$anova
stat_result$stat$tukey

write.table(stat_result$stat$anova, 'otu_AT_anova.txt', row.names = FALSE, sep = '\t', quote = FALSE)
write.table(stat_result$stat$tukey, 'otu_AT_tukey.txt', row.names = FALSE, sep = '\t', quote = FALSE)


#显著性标记“abc”作图结果
stat_result$plot$plot_bar
stat_result$plot$plot_box



####门水平phylum
####根据log(x+1)对数转换得到的数据，进行热图的绘制
library(pheatmap)

phylum_heatmap_AAT <- read.table('phylun_AT_data.txt', header = TRUE)
phylum_heatmap_AAT
phylum_heatmap_AAT <- subset(phylum_heatmap_AAT,select= -Group)
phylum_heatmap_AAT

phylum_heatmap_AAT <- as.data.frame(phylum_heatmap_AAT)

rownames(phylum_heatmap_AAT) <- phylum_heatmap_AAT$Row.names
rownames(phylum_heatmap_AAT) 
phylum_heatmap_AAT

phylum_heatmap_AAT <- subset(phylum_heatmap_AAT,select= -Row.names)
phylum_heatmap_AAT

phylum_heatmap_AAT <-t(phylum_heatmap_AAT)
phylum_heatmap_AAT

pheatmap(phylum_heatmap_AAT)
pheatmap(phylum_heatmap_AAT,border_color="white",cellheight=15,cellwidth=15,
         display_numbers=FALSE,
         cluster_col=F, show_rownames=T,show_colnames=T,
         fontsize_number=7,width=10,height=8)

####添加行列注释信息
# 读取实验设计注释列分组
design<- read.delim('./mapping.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
design
annot_data1 <- data.frame(row.names=colnames(phylum_heatmap_AAT), group=design$Group)
annot_data1

#group <- c("#CC6666","#D0DFE6FF","#999999","#0073c2","black")
#names(group) <- c("Ggt","H","SA","S_135","BS")


pheatmap(phylum_heatmap_AAT,
         border_color=NA,
         cellheight=10,
         cellwidth=16,
         cluster_col=F, show_rownames=T,show_colnames=F,treeheight_row=15,
         fontsize_number=7,
         width=10,height=8,
         frontsize=10,
         fontsize_row = 10,
         fontsize_col=10,
         annotation_col=annot_data1,
         angle_col = "270")

pheatmap(phylum_heatmap_AAT,
         border_color=NA,
         cellheight=10,
         cellwidth=16,
         cluster_col=F, show_rownames=T,show_colnames=F,treeheight_row=15,
         fontsize_number=7,
         width=10,height=8,
         frontsize=10,
         fontsize_row = 10,
         fontsize_col=10,
         annotation_col=annot_data1,
         angle_col = "270",
         filename="phylum_heatmap_AAT.pdf")

pheatmap(phylum_heatmap_AAT,
         border_color=NA,
         cellheight=10,
         cellwidth=16,
         cluster_col=F, show_rownames=T,show_colnames=F,treeheight_row=15,
         fontsize_number=7,
         width=10,height=8,
         frontsize=10,
         fontsize_row = 10,
         fontsize_col=10,
         annotation_col=annot_data1,
         angle_col = "270",
         filename="phylum_heatmap_AAT.png")

pheatmap(phylum_heatmap_AAT,
         border_color=NA,
         cellheight=10,
         cellwidth=16,
         cluster_col=F, show_rownames=T,show_colnames=F,treeheight_row=15,
         fontsize_number=7,
         width=10,height=8,
         frontsize=10,
         fontsize_row = 10,
         fontsize_col=10,
         annotation_col=annot_data1,
         angle_col = "270",
         filename="phylum_heatmap_AAT.tiff")




####纲水平class
####根据log(x+1)对数转换得到的数据，进行热图的绘制
library(pheatmap)

class_heatmap_AAT <- read.table('class_AT_data.txt', header = TRUE)
class_heatmap_AAT

class_heatmap_AAT <- subset(class_heatmap_AAT,select= -Group)
class_heatmap_AAT

class_heatmap_AAT <- as.data.frame(class_heatmap_AAT)

rownames(class_heatmap_AAT) <- class_heatmap_AAT$Row.names
rownames(class_heatmap_AAT) 
class_heatmap_AAT

class_heatmap_AAT <- subset(class_heatmap_AAT,select= -Row.names)
class_heatmap_AAT

class_heatmap_AAT <-t(class_heatmap_AAT)
class_heatmap_AAT


pheatmap(class_heatmap_AAT)
pheatmap(class_heatmap_AAT,border_color="white",cellheight=15,cellwidth=15,
         display_numbers=FALSE,
         cluster_col=F, show_rownames=T,show_colnames=T,
         fontsize_number=7,width=10,height=8)

####添加行列注释信息
# 读取实验设计注释列分组
design<- read.delim('./mapping.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
design
annot_data2 <- data.frame(row.names=colnames(class_heatmap_AAT), group=design$Group)
annot_data2

#group <- c("#CC6666","#D0DFE6FF","#999999","#0073c2","black")
#names(group) <- c("Ggt","H","SA","S_135","BS")


pheatmap(class_heatmap_AAT,
         border_color=NA,
         cellheight=10,
         cellwidth=16,
         cluster_col=F, show_rownames=T,show_colnames=F,treeheight_row=15,
         fontsize_number=7,
         width=10,height=8,
         frontsize=10,
         fontsize_row = 10,
         fontsize_col=10,
         annotation_col=annot_data2,
         angle_col = "270")

pheatmap(class_heatmap_AAT,
         border_color=NA,
         cellheight=10,
         cellwidth=16,
         cluster_col=F, show_rownames=T,show_colnames=F,treeheight_row=15,
         fontsize_number=7,
         width=10,height=8,
         frontsize=10,
         fontsize_row = 10,
         fontsize_col=10,
         annotation_col=annot_data2,
         angle_col = "270",
         filename="class_heatmap_AAT.pdf")

pheatmap(class_heatmap_AAT,
         border_color=NA,
         cellheight=10,
         cellwidth=16,
         cluster_col=F, show_rownames=T,show_colnames=F,treeheight_row=15,
         fontsize_number=7,
         width=10,height=8,
         frontsize=10,
         fontsize_row = 10,
         fontsize_col=10,
         annotation_col=annot_data2,
         angle_col = "270",
         filename="class_heatmap_AAT.png")

pheatmap(class_heatmap_AAT,
         border_color=NA,
         cellheight=10,
         cellwidth=16,
         cluster_col=F, show_rownames=T,show_colnames=F,treeheight_row=15,
         fontsize_number=7,
         width=10,height=8,
         frontsize=10,
         fontsize_row = 10,
         fontsize_col=10,
         annotation_col=annot_data2,
         angle_col = "270",
         filename="class_heatmap_AAT.tiff")





####目水平order
####根据log(x+1)对数转换得到的数据，进行热图的绘制
library(pheatmap)

order_heatmap_AAT <- read.table('order_AT_data.txt', header = TRUE)
order_heatmap_AAT

order_heatmap_AAT <- subset(order_heatmap_AAT,select= -Group)
order_heatmap_AAT

order_heatmap_AAT <- as.data.frame(order_heatmap_AAT)

rownames(order_heatmap_AAT) <- order_heatmap_AAT$Row.names
rownames(order_heatmap_AAT) 
order_heatmap_AAT

order_heatmap_AAT <- subset(order_heatmap_AAT,select= -Row.names)
order_heatmap_AAT

order_heatmap_AAT <-t(order_heatmap_AAT)
order_heatmap_AAT

pheatmap(order_heatmap_AAT)
pheatmap(order_heatmap_AAT,border_color="white",cellheight=15,cellwidth=15,
         display_numbers=FALSE,
         cluster_col=F, show_rownames=T,show_colnames=T,
         fontsize_number=7,width=10,height=8)

####添加行列注释信息
# 读取实验设计注释列分组
design<- read.delim('./mapping.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
design
annot_data3 <- data.frame(row.names=colnames(order_heatmap_AAT), group=design$Group)
annot_data3

#group <- c("#CC6666","#D0DFE6FF","#999999","#0073c2","black")
#names(group) <- c("Ggt","H","SA","S_135","BS")


pheatmap(order_heatmap_AAT,
         border_color=NA,
         cellheight=10,
         cellwidth=16,
         cluster_col=F, show_rownames=T,show_colnames=F,treeheight_row=15,
         fontsize_number=7,
         width=10,height=8,
         frontsize=10,
         fontsize_row = 10,
         fontsize_col=10,
         annotation_col=annot_data3,
         angle_col = "270")

pheatmap(order_heatmap_AAT,
         border_color=NA,
         cellheight=10,
         cellwidth=16,
         cluster_col=F, show_rownames=T,show_colnames=F,treeheight_row=15,
         fontsize_number=7,
         width=10,height=8,
         frontsize=10,
         fontsize_row = 10,
         fontsize_col=10,
         annotation_col=annot_data3,
         angle_col = "270",
         filename="order_heatmap_AAT.pdf")

pheatmap(order_heatmap_AAT,
         border_color=NA,
         cellheight=10,
         cellwidth=16,
         cluster_col=F, show_rownames=T,show_colnames=F,treeheight_row=15,
         fontsize_number=7,
         width=10,height=8,
         frontsize=10,
         fontsize_row = 10,
         fontsize_col=10,
         annotation_col=annot_data3,
         angle_col = "270",
         filename="order_heatmap_AAT.png")

pheatmap(order_heatmap_AAT,
         border_color=NA,
         cellheight=10,
         cellwidth=16,
         cluster_col=F, show_rownames=T,show_colnames=F,treeheight_row=15,
         fontsize_number=7,
         width=10,height=8,
         frontsize=10,
         fontsize_row = 10,
         fontsize_col=10,
         annotation_col=annot_data3,
         angle_col = "270",
         filename="order_heatmap_AAT.tiff")





####科水平family
####根据log(x+1)对数转换得到的数据，进行热图的绘制
library(pheatmap)

family_heatmap_AAT <- read.table('family_AT_data.txt', header = TRUE)
family_heatmap_AAT

family_heatmap_AAT <- subset(family_heatmap_AAT,select= -Group)
family_heatmap_AAT

family_heatmap_AAT <- as.data.frame(family_heatmap_AAT)

rownames(family_heatmap_AAT) <- family_heatmap_AAT$Row.names
rownames(family_heatmap_AAT) 
family_heatmap_AAT

family_heatmap_AAT <- subset(family_heatmap_AAT,select= -Row.names)
family_heatmap_AAT

family_heatmap_AAT <-t(family_heatmap_AAT)
family_heatmap_AAT

pheatmap(family_heatmap_AAT)
pheatmap(family_heatmap_AAT,border_color="white",cellheight=15,cellwidth=15,
         display_numbers=FALSE,
         cluster_col=F, show_rownames=T,show_colnames=T,
         fontsize_number=7,width=10,height=8)

####添加行列注释信息
# 读取实验设计注释列分组
design<- read.delim('./mapping.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
design
annot_data4 <- data.frame(row.names=colnames(family_heatmap_AAT), group=design$Group)
annot_data4

#group <- c("#CC6666","#D0DFE6FF","#999999","#0073c2","black")
#names(group) <- c("Ggt","H","SA","S_135","BS")


pheatmap(family_heatmap_AAT,
         border_color=NA,
         cellheight=10,
         cellwidth=16,
         cluster_col=F, show_rownames=T,show_colnames=F,treeheight_row=15,
         fontsize_number=7,
         width=10,height=8,
         frontsize=10,
         fontsize_row = 10,
         fontsize_col=10,
         annotation_col=annot_data4,
         angle_col = "270")

pheatmap(family_heatmap_AAT,
         border_color=NA,
         cellheight=10,
         cellwidth=16,
         cluster_col=F, show_rownames=T,show_colnames=F,treeheight_row=15,
         fontsize_number=7,
         width=10,height=8,
         frontsize=10,
         fontsize_row = 10,
         fontsize_col=10,
         annotation_col=annot_data4,
         angle_col = "270",
         filename="family_heatmap_AAT.pdf")

pheatmap(family_heatmap_AAT,
         border_color=NA,
         cellheight=10,
         cellwidth=16,
         cluster_col=F, show_rownames=T,show_colnames=F,treeheight_row=15,
         fontsize_number=7,
         width=10,height=8,
         frontsize=10,
         fontsize_row = 10,
         fontsize_col=10,
         annotation_col=annot_data4,
         angle_col = "270",
         filename="family_heatmap_AAT.png")

pheatmap(family_heatmap_AAT,
         border_color=NA,
         cellheight=10,
         cellwidth=16,
         cluster_col=F, show_rownames=T,show_colnames=F,treeheight_row=15,
         fontsize_number=7,
         width=10,height=8,
         frontsize=10,
         fontsize_row = 10,
         fontsize_col=10,
         annotation_col=annot_data4,
         angle_col = "270",
         filename="family_heatmap_AAT.tiff")



####属水平genus
####根据log(x+1)对数转换得到的数据，进行热图的绘制
library(pheatmap)

genus_heatmap_AAT <- read.table('genus_AT_data.txt', header = TRUE)
genus_heatmap_AAT

genus_heatmap_AAT <- subset(genus_heatmap_AAT,select= -Group)
genus_heatmap_AAT

genus_heatmap_AAT <- as.data.frame(genus_heatmap_AAT)

rownames(genus_heatmap_AAT) <- genus_heatmap_AAT$Row.names
rownames(genus_heatmap_AAT) 
genus_heatmap_AAT

genus_heatmap_AAT <- subset(genus_heatmap_AAT,select= -Row.names)
genus_heatmap_AAT

genus_heatmap_AAT <-t(genus_heatmap_AAT)
genus_heatmap_AAT

pheatmap(genus_heatmap_AAT)
pheatmap(genus_heatmap_AAT,border_color="white",cellheight=15,cellwidth=15,
         display_numbers=FALSE,
         cluster_col=F, show_rownames=T,show_colnames=T,
         fontsize_number=7,width=10,height=8)

####添加行列注释信息
# 读取实验设计注释列分组
design<- read.delim('./mapping.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
design
annot_data5 <- data.frame(row.names=colnames(genus_heatmap_AAT), group=design$Group)
annot_data5

#group <- c("#CC6666","#D0DFE6FF","#999999","#0073c2","black")
#names(group) <- c("Ggt","H","SA","S_135","BS")


pheatmap(genus_heatmap_AAT,
         border_color=NA,
         cellheight=10,
         cellwidth=16,
         cluster_col=F, show_rownames=T,show_colnames=F,treeheight_row=15,
         fontsize_number=7,
         width=10,height=8,
         frontsize=10,
         fontsize_row = 10,
         fontsize_col=10,
         annotation_col=annot_data5,
         angle_col = "270")

pheatmap(genus_heatmap_AAT,
         border_color=NA,
         cellheight=10,
         cellwidth=16,
         cluster_col=F, show_rownames=T,show_colnames=F,treeheight_row=15,
         fontsize_number=7,
         width=10,height=8,
         frontsize=10,
         fontsize_row = 10,
         fontsize_col=10,
         annotation_col=annot_data5,
         angle_col = "270",
         filename="genus_heatmap_AAT.pdf")

pheatmap(family_heatmap_AAT,
         border_color=NA,
         cellheight=10,
         cellwidth=16,
         cluster_col=F, show_rownames=T,show_colnames=F,treeheight_row=15,
         fontsize_number=7,
         width=10,height=8,
         frontsize=10,
         fontsize_row = 10,
         fontsize_col=10,
         annotation_col=annot_data5,
         angle_col = "270",
         filename="genus_heatmap_AAT.png")

pheatmap(genus_heatmap_AAT,
         border_color=NA,
         cellheight=10,
         cellwidth=16,
         cluster_col=F, show_rownames=T,show_colnames=F,treeheight_row=15,
         fontsize_number=7,
         width=10,height=8,
         frontsize=10,
         fontsize_row = 10,
         fontsize_col=10,
         annotation_col=annot_data5,
         angle_col = "270",
         filename="genus_heatmap_AAT.tiff")


####方法二：DESeq2差异分析####
setwd("E:/博士数据2/小麦微生物组数据/customer_backup/otus/DESeq2/")

library(DESeq2)
library(ggplot2)
# 读入实验设计

mapping  =read.table("mapping.txt", header=T, row.names= 1,sep="\t")
mapping

# 读取OTU表,全部otu表没有抽平，基于count的数据，不可用相对丰度数据

otu =read.delim("otu_AAT.txt", row.names= 1,sep="\t",header=T,check.names=F)
otu

#第一步，构建 DESeqDataSet 对象，详见 ?DESeqDataSetFromMatrix
dds <- DESeqDataSetFromMatrix(countData = otu, colData = mapping, design = ~Group)
dds

#第二步，差异分析，详见 ?DESeq 和 ?results
#标准方法
dds <- DESeq(dds, parallel = FALSE)	#parallel = TRUE 将启用多线程模式
suppressMessages(dds)

#第三步，两两之间的比较

####Ggt vs S_135####

#以Ggt为对照进行比较——Ggt vs S_135
res_Ggt_S_135 <- results(dds, contrast = c('Group', 'S_135', 'Ggt'), pAdjustMethod = 'fdr', alpha = 0.05)
res_Ggt_S_135
#或者summary(res_Ggt_S_135)

#可以先按校正和 p 值由小到大排个序，方便查看
Ggt_S_135 <- as.data.frame(res_Ggt_S_135[order(res_Ggt_S_135$padj), ])
Ggt_S_135
Ggt_S_135$OTU.ID <- rownames(Ggt_S_135)
write.table(Ggt_S_135, 'Ggt_S_135.txt', row.names = FALSE, sep = '\t', quote = FALSE)

#例如这里自定义根据 |log2FC| >= 1 和 adj.P.Val < 0.01 标记差异类型
#此处有关Up、Down的定义是根据对照和处理而言的，由于本人前期将对照和
Ggt_S_135[which(Ggt_S_135$padj %in% NA),"sig"] <- "None"
Ggt_S_135[which(Ggt_S_135$padj < 0.05 & Ggt_S_135$log2FoldChange <= -1),"sig"] <- "Down"
Ggt_S_135[which(Ggt_S_135$padj < 0.05 & Ggt_S_135$log2FoldChange >= 1),"sig"] <- "Up"
Ggt_S_135[which(Ggt_S_135$padj >= 0.05 | abs(Ggt_S_135$log2FoldChange) < 1),"sig"] <- "None"
Ggt_S_135

taxonomy <- read.delim("taxonomy2.txt", sep="\t")
taxonomy

Ggt_S_135_1 <- merge(Ggt_S_135, taxonomy, by="OTU.ID")
Ggt_S_135_1
Ggt_S_135_1 <- as.data.frame(Ggt_S_135_1[order(Ggt_S_135_1$padj), ])
Ggt_S_135_1
write.table(Ggt_S_135_1, 'Ggt_S_135_1.txt', row.names = FALSE, sep = '\t', quote = FALSE)


##标注上、下调
up <- subset(Ggt_S_135_1, sig == 'Up')
#up <- up[order(up$padj), ][1:15, ]
down <- subset(Ggt_S_135_1, sig == 'Down')
#down <- down[order(down$padj), ][1:10, ]


##门水平的注释
##对上下调的颜色进行标注
Ggt_S_135_1[which(Ggt_S_135_1$sig == "Up"),"color1"] <- up$phylum
Ggt_S_135_1[which(Ggt_S_135_1$sig == "Down" ),"color1"] <- down$phylum

color1 <- gsub(".*(p__)","",Ggt_S_135_1$color1)
color1
#gsub()可以用于字段的删减、增补、替换和切割，可以处理一个字段也可以处理由字段组成的向量。
?gsub
#删除Taxon中有关Bacteria;的字样
#taxon <- gsub(".*(Bacteria;)","",taxon)
#删除taxon中的Other字样
#taxon <- gsub(".*(Other)","",taxon)

color1 <- gsub("(\\[)","",color1)
color1 <- gsub("(\\])","",color1)
#将新产生的表头与原来的表头进行替换
Ggt_S_135_1$color1 <- color1
Ggt_S_135_1


##纲水平的注释
##对上下调的颜色进行标注
Ggt_S_135_1[which(Ggt_S_135_1$sig == "Up"),"color2"] <- up$class
Ggt_S_135_1[which(Ggt_S_135_1$sig == "Down" ),"color2"] <- down$class

color2 <- gsub(".*(c__)","",Ggt_S_135_1$color2)
color2
#gsub()可以用于字段的删减、增补、替换和切割，可以处理一个字段也可以处理由字段组成的向量。
?gsub
#删除Taxon中有关Bacteria;的字样
#taxon <- gsub(".*(Bacteria;)","",taxon)
#删除taxon中的Other字样
#taxon <- gsub(".*(Other)","",taxon)

color2 <- gsub("(\\[)","",color2)
color2 <- gsub("(\\])","",color2)
#将新产生的表头与原来的表头进行替换
Ggt_S_135_1$color2 <- color2
Ggt_S_135_1


##目水平的注释
##对上下调的颜色进行标注
Ggt_S_135_1[which(Ggt_S_135_1$sig == "Up"),"color3"] <- up$order
Ggt_S_135_1[which(Ggt_S_135_1$sig == "Down" ),"color3"] <- down$order

color3 <- gsub(".*(o__)","",Ggt_S_135_1$color3)
color3
#gsub()可以用于字段的删减、增补、替换和切割，可以处理一个字段也可以处理由字段组成的向量。
?gsub
#删除Taxon中有关Bacteria;的字样
#taxon <- gsub(".*(Bacteria;)","",taxon)
#删除taxon中的Other字样
#taxon <- gsub(".*(Other)","",taxon)

color3 <- gsub("(\\[)","",color3)
color3 <- gsub("(\\])","",color3)
#将新产生的表头与原来的表头进行替换
Ggt_S_135_1$color3 <- color3
Ggt_S_135_1

##科水平的注释
##对上下调的颜色进行标注
Ggt_S_135_1[which(Ggt_S_135_1$sig == "Up"),"color4"] <- up$family
Ggt_S_135_1[which(Ggt_S_135_1$sig == "Down" ),"color4"] <- down$family

color4 <- gsub(".*(f__)","",Ggt_S_135_1$color4)
color4
#gsub()可以用于字段的删减、增补、替换和切割，可以处理一个字段也可以处理由字段组成的向量。
?gsub
#删除Taxon中有关Bacteria;的字样
#taxon <- gsub(".*(Bacteria;)","",taxon)
#删除taxon中的Other字样
#taxon <- gsub(".*(Other)","",taxon)

color4 <- gsub("(\\[)","",color4)
color4 <- gsub("(\\])","",color4)
#将新产生的表头与原来的表头进行替换
Ggt_S_135_1$color4 <- color4
Ggt_S_135_1

##属水平的注释
##对上下调的颜色进行标注
Ggt_S_135_1[which(Ggt_S_135_1$sig == "Up"),"color5"] <- up$genus
Ggt_S_135_1[which(Ggt_S_135_1$sig == "Down" ),"color5"] <- down$genus

color5 <- gsub(".*(g__)","",Ggt_S_135_1$color5)
color5
#gsub()可以用于字段的删减、增补、替换和切割，可以处理一个字段也可以处理由字段组成的向量。
?gsub
#删除Taxon中有关Bacteria;的字样
#taxon <- gsub(".*(Bacteria;)","",taxon)
#删除taxon中的Other字样
#taxon <- gsub(".*(Other)","",taxon)

color5 <- gsub("(\\[)","",color5)
color5 <- gsub("(\\])","",color5)
#将新产生的表头与原来的表头进行替换
Ggt_S_135_1$color5 <- color5
Ggt_S_135_1


write.table(Ggt_S_135_1, 'Ggt_S_135_2.txt', row.names = FALSE, sep = '\t', quote = FALSE)


##第四步，比较结果出图
#门水平的注释

#横轴 log2FoldChange，纵轴 -log10(adj.P.Val)，颜色表示差异
pp1 <- ggplot(Ggt_S_135_1, aes(x = log2FoldChange, y = -log10(padj)),col=mycolor) +
  geom_point(aes(shape = sig, color=color1), alpha = 0.6, size = 2) +
  scale_shape_manual(values = c(17, 19, 15),limits = c("Up","None","Down" ))
theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent')) +
  theme(legend.title = element_blank(), legend.key = element_rect(fill = 'transparent'), legend.background = element_rect(fill = 'transparent'), legend.position = c(0.9, 0.9))

pp1


#调整绘图区主题。
pp2 <- pp1 + theme_bw()+ theme(panel.grid=element_blank()) + 
  theme(panel.border = element_blank()) + 
  xlim(-5, 5)+ylim(0, 6)+
  theme(panel.background=element_rect(fill='transparent',color='black'),plot.margin = unit(c(3,5,1,1),"mm"))
pp2

#调整坐标轴文字。

pp3 <- pp2 + theme(axis.text.x=element_text(colour="black",angle = 0,vjust = 0.5,hjust = 0)) + 
  theme(axis.text.y=element_text(colour = "black")) + 
  theme(axis.line = element_line(colour = "black"))+ 
  theme(axis.title.y = element_text(colour = "black", margin = unit(c(0,1,0,1),"lines"))) 
pp3

pp4<-pp3+geom_vline(xintercept = c(0, 0), color = 'gray', size = 0.25) +
  geom_hline(yintercept = -log(0.05, 10), color = 'gray', size = 0.25) +
  labs(x = 'log2 Fold Change', y = '-log10 adjp-value', color = NA,title ="Ggt vs S_135") 

pp4
#输出图像。
ggsave("Ggt_s_135_phylum.pdf", width = 14, height = 10, units="cm")
ggsave("Ggt_s_135_phylum.tiff", width = 14, height = 10, units = "cm")
ggsave("Ggt_s_135_phylum.png", width = 14, height = 10, units = "cm")



##纲水平的注释

#横轴 log2FoldChange，纵轴 -log10(adj.P.Val)，颜色表示差异
pc1 <- ggplot(Ggt_S_135_1, aes(x = log2FoldChange, y = -log10(padj)),col=mycolor) +
  geom_point(aes(shape = sig, color=color2), alpha = 0.6, size = 2) +
  scale_shape_manual(values = c(17, 19, 15),limits = c("Up","None","Down" ))
theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent')) +
  theme(legend.title = element_blank(), legend.key = element_rect(fill = 'transparent'), legend.background = element_rect(fill = 'transparent'), legend.position = c(0.9, 0.9))

pc1


#调整绘图区主题。
pc2 <- pc1 + theme_bw()+ theme(panel.grid=element_blank()) + 
  theme(panel.border = element_blank()) + 
  xlim(-5, 5)+ylim(0, 6)+
  theme(panel.background=element_rect(fill='transparent',color='black'),plot.margin = unit(c(3,5,1,1),"mm"))
pc2

#调整坐标轴文字。

pc3 <- pc2 + theme(axis.text.x=element_text(colour="black",angle = 0,vjust = 0.5,hjust = 0)) + 
  theme(axis.text.y=element_text(colour = "black")) + 
  theme(axis.line = element_line(colour = "black"))+ 
  theme(axis.title.y = element_text(colour = "black", margin = unit(c(0,1,0,1),"lines"))) 
pc3

pc4<-pc3+geom_vline(xintercept = c(0, 0), color = 'gray', size = 0.25) +
  geom_hline(yintercept = -log(0.05, 10), color = 'gray', size = 0.25) +
  labs(x = 'log2 Fold Change', y = '-log10 adjp-value', color = NA,title ="Ggt vs S_135") 

pc4
#输出图像。
ggsave("Ggt_s_135_class.pdf", width = 14, height = 10, units="cm")
ggsave("Ggt_s_135_class.tiff", width = 14, height = 10, units = "cm")
ggsave("Ggt_s_135_class.png", width = 14, height = 10, units = "cm")



##目水平的注释

#横轴 log2FoldChange，纵轴 -log10(adj.P.Val)，颜色表示差异
po1 <- ggplot(Ggt_S_135_1, aes(x = log2FoldChange, y = -log10(padj)),col=mycolor) +
  geom_point(aes(shape = sig, color=color3), alpha = 0.6, size = 2) +
  scale_shape_manual(values = c(17, 19, 15),limits = c("Up","None","Down" ))
theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent')) +
  theme(legend.title = element_blank(), legend.key = element_rect(fill = 'transparent'), legend.background = element_rect(fill = 'transparent'), legend.position = c(0.9, 0.9))

po1


#调整绘图区主题。
po2 <- po1 + theme_bw()+ theme(panel.grid=element_blank()) + 
  theme(panel.border = element_blank()) + 
  xlim(-5, 5)+ylim(0, 6)+
  theme(panel.background=element_rect(fill='transparent',color='black'),plot.margin = unit(c(3,5,1,1),"mm"))
po2

#调整坐标轴文字。

po3 <- po2 + theme(axis.text.x=element_text(colour="black",angle = 0,vjust = 0.5,hjust = 0)) + 
  theme(axis.text.y=element_text(colour = "black")) + 
  theme(axis.line = element_line(colour = "black"))+ 
  theme(axis.title.y = element_text(colour = "black", margin = unit(c(0,1,0,1),"lines"))) 
po3

po4<-po3+geom_vline(xintercept = c(0, 0), color = 'gray', size = 0.25) +
  geom_hline(yintercept = -log(0.05, 10), color = 'gray', size = 0.25) +
  labs(x = 'log2 Fold Change', y = '-log10 adjp-value', color = NA,title ="Ggt vs S_135") 

po4
#输出图像。
ggsave("Ggt_s_135_order.pdf", width = 14, height = 10, units="cm")
ggsave("Ggt_s_135_order.tiff", width = 14, height = 10, units = "cm")
ggsave("Ggt_s_135_order.png", width = 14, height = 10, units = "cm")


##科水平的注释
#横轴 log2FoldChange，纵轴 -log10(adj.P.Val)，颜色表示差异
pf1 <- ggplot(Ggt_S_135_1, aes(x = log2FoldChange, y = -log10(padj)),col=mycolor) +
  geom_point(aes(shape = sig, color=color4), alpha = 0.6, size = 2) +
  scale_shape_manual(values = c(17, 19, 15),limits = c("Up","None","Down" ))
theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent')) +
  theme(legend.title = element_blank(), legend.key = element_rect(fill = 'transparent'), legend.background = element_rect(fill = 'transparent'), legend.position = c(0.9, 0.9))

pf1


#调整绘图区主题。
pf2 <- pf1 + theme_bw()+ theme(panel.grid=element_blank()) + 
  theme(panel.border = element_blank()) + 
  xlim(-5, 5)+ylim(0, 6)+
  theme(panel.background=element_rect(fill='transparent',color='black'),plot.margin = unit(c(3,5,1,1),"mm"))
pf2

#调整坐标轴文字。

pf3 <- pf2 + theme(axis.text.x=element_text(colour="black",angle = 0,vjust = 0.5,hjust = 0)) + 
  theme(axis.text.y=element_text(colour = "black")) + 
  theme(axis.line = element_line(colour = "black"))+ 
  theme(axis.title.y = element_text(colour = "black", margin = unit(c(0,1,0,1),"lines"))) 
pf3

pf4<-pf3+geom_vline(xintercept = c(0, 0), color = 'gray', size = 0.25) +
  geom_hline(yintercept = -log(0.05, 10), color = 'gray', size = 0.25) +
  labs(x = 'log2 Fold Change', y = '-log10 adjp-value', color = NA,title ="Ggt vs S_135") 

pf4
#输出图像。
ggsave("Ggt_s_135_family.pdf", width = 14, height = 10, units="cm")
ggsave("Ggt_s_135_family.tiff", width = 14, height = 10, units = "cm")
ggsave("Ggt_s_135_family.png", width = 14, height = 10, units = "cm")


##属水平的注释
#横轴 log2FoldChange，纵轴 -log10(adj.P.Val)，颜色表示差异
pg1 <- ggplot(Ggt_S_135_1, aes(x = log2FoldChange, y = -log10(padj)),col=mycolor) +
  geom_point(aes(shape = sig, color=color5), alpha = 0.6, size = 2) +
  scale_shape_manual(values = c(17, 19, 15),limits = c("Up","None","Down" ))
theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent')) +
  theme(legend.title = element_blank(), legend.key = element_rect(fill = 'transparent'), legend.background = element_rect(fill = 'transparent'), legend.position = c(0.9, 0.9))

pg1


#调整绘图区主题。
pg2 <- pg1 + theme_bw()+ theme(panel.grid=element_blank()) + 
  theme(panel.border = element_blank()) + 
  xlim(-5, 5)+ylim(0, 6)+
  theme(panel.background=element_rect(fill='transparent',color='black'),plot.margin = unit(c(3,5,1,1),"mm"))
pg2

#调整坐标轴文字。

pg3 <- pg2 + theme(axis.text.x=element_text(colour="black",angle = 0,vjust = 0.5,hjust = 0)) + 
  theme(axis.text.y=element_text(colour = "black")) + 
  theme(axis.line = element_line(colour = "black"))+ 
  theme(axis.title.y = element_text(colour = "black", margin = unit(c(0,1,0,1),"lines"))) 
pg3

pg4<-pg3+geom_vline(xintercept = c(0, 0), color = 'gray', size = 0.25) +
  geom_hline(yintercept = -log(0.05, 10), color = 'gray', size = 0.25) +
  labs(x = 'log2 Fold Change', y = '-log10 adjp-value', color = NA,title ="Ggt vs S_135") 

pg4

##在图中展示OTU的名称
library(ggrepel)

up <- subset(phylum, sig == 'Up')
up <- up[order(up$padj), ][1:10, ]
down <- subset(phylum, sig == 'Down')
#down <- down[order(down$padj), ][1:10, ]

#一种样式，借助 ggrepel 包中的函数
pg5 <- pg4 + theme(legend.position = 'right') +
  geom_text_repel(data = rbind(up, down), aes(x = log2FoldChange, y = -log10(padj), label = Features, color= Features),
                  size = 2,box.padding = unit(0.5, 'lines'), segment.color = 'black', show.legend = FALSE)

pg5

#输出图像。
ggsave("Ggt_s_135_genus.pdf", width = 14, height = 10, units="cm")
ggsave("Ggt_s_135_genus.tiff", width = 14, height = 10, units = "cm")
ggsave("Ggt_s_135_genus.png", width = 14, height = 10, units = "cm")



#纵轴为基因表达值的 log10
volcano_count <- ggplot(Ggt_S_135_1, aes(log2FoldChange, log10(baseMean))) +
  geom_point(aes(shape = sig, color=color1), alpha = 0.6, size = 2) +
  scale_shape_manual(values = c(17, 19, 15),limits = c("Up","None","Down" )) +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), legend.position = c(0.8, 0.5)) +
  theme(legend.title = element_blank(), legend.key = element_rect(fill = 'transparent'), legend.background = element_rect(fill = 'transparent')) +
  geom_vline(xintercept = c(-1, 1), color = 'gray', size = 0.25) + 
  labs(x = 'log2 Fold Change', y = 'Average log10 baseMean') +
  xlim(-5, 5) +
  coord_flip()
volcano_count 
ggsave('volcano_count.df', volcano_count, width = 7, height = 5)
ggsave('volcano_count.png', volcano_count, width = 7, height = 5)


####Ggt vs SA####

#以Ggt为对照进行比较——Ggt vs SA
res_Ggt_SA <- results(dds, contrast = c('Group', 'SA', 'Ggt'), pAdjustMethod = 'fdr', alpha = 0.05)
res_Ggt_SA
#或者summary(res_Ggt_SA)

#可以先按校正和 p 值由小到大排个序，方便查看
Ggt_SA <- as.data.frame(res_Ggt_SA[order(res_Ggt_SA$padj), ])
Ggt_SA
Ggt_SA$OTU.ID <- rownames(Ggt_SA)
write.table(Ggt_SA, 'Ggt_SA.txt', row.names = FALSE, sep = '\t', quote = FALSE)

#例如这里自定义根据 |log2FC| >= 1 和 adj.P.Val < 0.01 标记差异类型
#此处有关Up、Down的定义是根据对照和处理而言的，由于本人前期将对照和
Ggt_SA[which(Ggt_SA$padj %in% NA),"sig"] <- "None"
Ggt_SA[which(Ggt_SA$padj < 0.05 & Ggt_SA$log2FoldChange <= -1),"sig"] <- "Down"
Ggt_SA[which(Ggt_SA$padj < 0.05 & Ggt_SA$log2FoldChange >= 1),"sig"] <- "Up"
Ggt_SA[which(Ggt_SA$padj >= 0.05 | abs(Ggt_SA$log2FoldChange) < 1),"sig"] <- "None"
Ggt_SA

taxonomy <- read.delim("taxonomy2.txt", sep="\t")
taxonomy

Ggt_SA_1 <- merge(Ggt_SA, taxonomy, by="OTU.ID")
Ggt_SA_1
Ggt_SA_1 <- as.data.frame(Ggt_SA_1[order(Ggt_SA_1$padj), ])
Ggt_SA_1
write.table(Ggt_SA_1, 'Ggt_SA_1.txt', row.names = FALSE, sep = '\t', quote = FALSE)


##标注上、下调
up <- subset(Ggt_SA_1, sig == 'Up')
#up <- up[order(up$padj), ][1:15, ]
down <- subset(Ggt_SA_1, sig == 'Down')
#down <- down[order(down$padj), ][1:10, ]


##门水平的注释
##对上下调的颜色进行标注
Ggt_SA_1[which(Ggt_SA_1$sig == "Up"),"color1"] <- up$phylum
Ggt_SA_1[which(Ggt_SA_1$sig == "Down" ),"color1"] <- down$phylum

color1 <- gsub(".*(p__)","",Ggt_SA_1$color1)
color1
#gsub()可以用于字段的删减、增补、替换和切割，可以处理一个字段也可以处理由字段组成的向量。
?gsub
#删除Taxon中有关Bacteria;的字样
#taxon <- gsub(".*(Bacteria;)","",taxon)
#删除taxon中的Other字样
#taxon <- gsub(".*(Other)","",taxon)

color1 <- gsub("(\\[)","",color1)
color1 <- gsub("(\\])","",color1)
#将新产生的表头与原来的表头进行替换
Ggt_SA_1$color1 <- color1
Ggt_SA_1


##纲水平的注释
##对上下调的颜色进行标注
Ggt_SA_1[which(Ggt_SA_1$sig == "Up"),"color2"] <- up$class
Ggt_SA_1[which(Ggt_SA_1$sig == "Down" ),"color2"] <- down$class

color2 <- gsub(".*(c__)","",Ggt_SA_1$color2)
color2
#gsub()可以用于字段的删减、增补、替换和切割，可以处理一个字段也可以处理由字段组成的向量。
?gsub
#删除Taxon中有关Bacteria;的字样
#taxon <- gsub(".*(Bacteria;)","",taxon)
#删除taxon中的Other字样
#taxon <- gsub(".*(Other)","",taxon)

color2 <- gsub("(\\[)","",color2)
color2 <- gsub("(\\])","",color2)
#将新产生的表头与原来的表头进行替换
Ggt_SA_1$color2 <- color2
Ggt_SA_1


##目水平的注释
##对上下调的颜色进行标注
Ggt_SA_1[which(Ggt_SA_1$sig == "Up"),"color3"] <- up$order
Ggt_SA_1[which(Ggt_SA_1$sig == "Down" ),"color3"] <- down$order

color3 <- gsub(".*(o__)","",Ggt_SA_1$color3)
color3
#gsub()可以用于字段的删减、增补、替换和切割，可以处理一个字段也可以处理由字段组成的向量。
?gsub
#删除Taxon中有关Bacteria;的字样
#taxon <- gsub(".*(Bacteria;)","",taxon)
#删除taxon中的Other字样
#taxon <- gsub(".*(Other)","",taxon)

color3 <- gsub("(\\[)","",color3)
color3 <- gsub("(\\])","",color3)
#将新产生的表头与原来的表头进行替换
Ggt_SA_1$color3 <- color3
Ggt_SA_1

##科水平的注释
##对上下调的颜色进行标注
Ggt_SA_1[which(Ggt_SA_1$sig == "Up"),"color4"] <- up$family
Ggt_SA_1[which(Ggt_SA_1$sig == "Down" ),"color4"] <- down$family

color4 <- gsub(".*(f__)","",Ggt_SA_1$color4)
color4
#gsub()可以用于字段的删减、增补、替换和切割，可以处理一个字段也可以处理由字段组成的向量。
?gsub
#删除Taxon中有关Bacteria;的字样
#taxon <- gsub(".*(Bacteria;)","",taxon)
#删除taxon中的Other字样
#taxon <- gsub(".*(Other)","",taxon)

color4 <- gsub("(\\[)","",color4)
color4 <- gsub("(\\])","",color4)
#将新产生的表头与原来的表头进行替换
Ggt_SA_1$color4 <- color4
Ggt_SA_1

##属水平的注释
##对上下调的颜色进行标注
Ggt_SA_1[which(Ggt_SA_1$sig == "Up"),"color5"] <- up$genus
Ggt_SA_1[which(Ggt_SA_1$sig == "Down" ),"color5"] <- down$genus

color5 <- gsub(".*(g__)","",Ggt_SA_1$color5)
color5
#gsub()可以用于字段的删减、增补、替换和切割，可以处理一个字段也可以处理由字段组成的向量。
?gsub
#删除Taxon中有关Bacteria;的字样
#taxon <- gsub(".*(Bacteria;)","",taxon)
#删除taxon中的Other字样
#taxon <- gsub(".*(Other)","",taxon)

color5 <- gsub("(\\[)","",color5)
color5 <- gsub("(\\])","",color5)
#将新产生的表头与原来的表头进行替换
Ggt_SA_1$color5 <- color5
Ggt_SA_1


write.table(Ggt_SA_1, 'Ggt_SA_2.txt', row.names = FALSE, sep = '\t', quote = FALSE)


##第四步，比较结果出图
#门水平的注释

#横轴 log2FoldChange，纵轴 -log10(adj.P.Val)，颜色表示差异
pp1 <- ggplot(Ggt_SA_1, aes(x = log2FoldChange, y = -log10(padj)),col=mycolor) +
  geom_point(aes(shape = sig, color=color1), alpha = 0.6, size = 2) +
  scale_shape_manual(values = c(17, 19, 15),limits = c("Up","None","Down" ))
theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent')) +
  theme(legend.title = element_blank(), legend.key = element_rect(fill = 'transparent'), legend.background = element_rect(fill = 'transparent'), legend.position = c(0.9, 0.9))

pp1


#调整绘图区主题。
pp2 <- pp1 + theme_bw()+ theme(panel.grid=element_blank()) + 
  theme(panel.border = element_blank()) + 
  xlim(-5, 5)+ylim(0, 6)+
  theme(panel.background=element_rect(fill='transparent',color='black'),plot.margin = unit(c(3,5,1,1),"mm"))
pp2

#调整坐标轴文字。

pp3 <- pp2 + theme(axis.text.x=element_text(colour="black",angle = 0,vjust = 0.5,hjust = 0)) + 
  theme(axis.text.y=element_text(colour = "black")) + 
  theme(axis.line = element_line(colour = "black"))+ 
  theme(axis.title.y = element_text(colour = "black", margin = unit(c(0,1,0,1),"lines"))) 
pp3

pp4<-pp3+geom_vline(xintercept = c(0, 0), color = 'gray', size = 0.25) +
  geom_hline(yintercept = -log(0.05, 10), color = 'gray', size = 0.25) +
  labs(x = 'log2 Fold Change', y = '-log10 adjp-value', color = NA,title ="Ggt vs SA") 

pp4
#输出图像。
ggsave("Ggt_sA_phylum.pdf", width = 14, height = 10, units="cm")
ggsave("Ggt_sA_phylum.tiff", width = 14, height = 10, units = "cm")
ggsave("Ggt_sA_phylum.png", width = 14, height = 10, units = "cm")



####Ggt vs H####

#以Ggt为对照进行比较——Ggt vs H
res_Ggt_H <- results(dds, contrast = c('Group', 'H', 'Ggt'), pAdjustMethod = 'fdr', alpha = 0.05)
res_Ggt_H
#或者summary(res_Ggt_H)

#可以先按校正和 p 值由小到大排个序，方便查看
Ggt_H <- as.data.frame(res_Ggt_H[order(res_Ggt_H$padj), ])
Ggt_H
Ggt_H$OTU.ID <- rownames(Ggt_H)
write.table(Ggt_H, 'Ggt_H.txt', row.names = FALSE, sep = '\t', quote = FALSE)

#例如这里自定义根据 |log2FC| >= 1 和 adj.P.Val < 0.01 标记差异类型
#此处有关Up、Down的定义是根据对照和处理而言的，由于本人前期将对照和
Ggt_H[which(Ggt_H$padj %in% NA),"sig"] <- "None"
Ggt_H[which(Ggt_H$padj < 0.05 & Ggt_H$log2FoldChange <= -1),"sig"] <- "Down"
Ggt_H[which(Ggt_H$padj < 0.05 & Ggt_H$log2FoldChange >= 1),"sig"] <- "Up"
Ggt_H[which(Ggt_H$padj >= 0.05 | abs(Ggt_H$log2FoldChange) < 1),"sig"] <- "None"
Ggt_H

taxonomy <- read.delim("taxonomy2.txt", sep="\t")
taxonomy

Ggt_H_1 <- merge(Ggt_H, taxonomy, by="OTU.ID")
Ggt_H_1
Ggt_H_1 <- as.data.frame(Ggt_H_1[order(Ggt_H_1$padj), ])
Ggt_H_1
write.table(Ggt_H_1, 'Ggt_H_1.txt', row.names = FALSE, sep = '\t', quote = FALSE)


##标注上、下调
up <- subset(Ggt_H_1, sig == 'Up')
#up <- up[order(up$padj), ][1:15, ]
down <- subset(Ggt_H_1, sig == 'Down')
#down <- down[order(down$padj), ][1:10, ]


##门水平的注释
##对上下调的颜色进行标注
Ggt_H_1[which(Ggt_H_1$sig == "Up"),"color1"] <- up$phylum
Ggt_H_1[which(Ggt_H_1$sig == "Down" ),"color1"] <- down$phylum

color1 <- gsub(".*(p__)","",Ggt_H_1$color1)
color1
#gsub()可以用于字段的删减、增补、替换和切割，可以处理一个字段也可以处理由字段组成的向量。
?gsub
#删除Taxon中有关Bacteria;的字样
#taxon <- gsub(".*(Bacteria;)","",taxon)
#删除taxon中的Other字样
#taxon <- gsub(".*(Other)","",taxon)

color1 <- gsub("(\\[)","",color1)
color1 <- gsub("(\\])","",color1)
#将新产生的表头与原来的表头进行替换
Ggt_H_1$color1 <- color1
Ggt_H_1


##纲水平的注释
##对上下调的颜色进行标注
Ggt_H_1[which(Ggt_H_1$sig == "Up"),"color2"] <- up$class
Ggt_H_1[which(Ggt_H_1$sig == "Down" ),"color2"] <- down$class

color2 <- gsub(".*(c__)","",Ggt_H_1$color2)
color2
#gsub()可以用于字段的删减、增补、替换和切割，可以处理一个字段也可以处理由字段组成的向量。
?gsub
#删除Taxon中有关Bacteria;的字样
#taxon <- gsub(".*(Bacteria;)","",taxon)
#删除taxon中的Other字样
#taxon <- gsub(".*(Other)","",taxon)

color2 <- gsub("(\\[)","",color2)
color2 <- gsub("(\\])","",color2)
#将新产生的表头与原来的表头进行替换
Ggt_H_1$color2 <- color2
Ggt_H_1


##目水平的注释
##对上下调的颜色进行标注
Ggt_H_1[which(Ggt_H_1$sig == "Up"),"color3"] <- up$order
Ggt_H_1[which(Ggt_H_1$sig == "Down" ),"color3"] <- down$order

color3 <- gsub(".*(o__)","",Ggt_H_1$color3)
color3
#gsub()可以用于字段的删减、增补、替换和切割，可以处理一个字段也可以处理由字段组成的向量。
?gsub
#删除Taxon中有关Bacteria;的字样
#taxon <- gsub(".*(Bacteria;)","",taxon)
#删除taxon中的Other字样
#taxon <- gsub(".*(Other)","",taxon)

color3 <- gsub("(\\[)","",color3)
color3 <- gsub("(\\])","",color3)
#将新产生的表头与原来的表头进行替换
Ggt_H_1$color3 <- color3
Ggt_H_1

##科水平的注释
##对上下调的颜色进行标注
Ggt_H_1[which(Ggt_H_1$sig == "Up"),"color4"] <- up$family
Ggt_H_1[which(Ggt_H_1$sig == "Down" ),"color4"] <- down$family

color4 <- gsub(".*(f__)","",Ggt_H_1$color4)
color4
#gsub()可以用于字段的删减、增补、替换和切割，可以处理一个字段也可以处理由字段组成的向量。
?gsub
#删除Taxon中有关Bacteria;的字样
#taxon <- gsub(".*(Bacteria;)","",taxon)
#删除taxon中的Other字样
#taxon <- gsub(".*(Other)","",taxon)

color4 <- gsub("(\\[)","",color4)
color4 <- gsub("(\\])","",color4)
#将新产生的表头与原来的表头进行替换
Ggt_H_1$color4 <- color4
Ggt_H_1

##属水平的注释
##对上下调的颜色进行标注
Ggt_H_1[which(Ggt_H_1$sig == "Up"),"color5"] <- up$genus
Ggt_H_1[which(Ggt_H_1$sig == "Down" ),"color5"] <- down$genus

color5 <- gsub(".*(g__)","",Ggt_H_1$color5)
color5
#gsub()可以用于字段的删减、增补、替换和切割，可以处理一个字段也可以处理由字段组成的向量。
?gsub
#删除Taxon中有关Bacteria;的字样
#taxon <- gsub(".*(Bacteria;)","",taxon)
#删除taxon中的Other字样
#taxon <- gsub(".*(Other)","",taxon)

color5 <- gsub("(\\[)","",color5)
color5 <- gsub("(\\])","",color5)
#将新产生的表头与原来的表头进行替换
Ggt_H_1$color5 <- color5
Ggt_H_1


write.table(Ggt_H_1, 'Ggt_H_2.txt', row.names = FALSE, sep = '\t', quote = FALSE)


##第四步，比较结果出图
#门水平的注释

#横轴 log2FoldChange，纵轴 -log10(adj.P.Val)，颜色表示差异
pp1 <- ggplot(Ggt_H_1, aes(x = log2FoldChange, y = -log10(padj)),col=mycolor) +
  geom_point(aes(shape = sig, color=color1), alpha = 0.6, size = 2) +
  scale_shape_manual(values = c(17, 19, 15),limits = c("Up","None","Down" ))
theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent')) +
  theme(legend.title = element_blank(), legend.key = element_rect(fill = 'transparent'), legend.background = element_rect(fill = 'transparent'), legend.position = c(0.9, 0.9))

pp1


#调整绘图区主题。
pp2 <- pp1 + theme_bw()+ theme(panel.grid=element_blank()) + 
  theme(panel.border = element_blank()) + 
  xlim(-5, 5)+ylim(0, 16)+
  theme(panel.background=element_rect(fill='transparent',color='black'),plot.margin = unit(c(3,5,1,1),"mm"))
pp2

#调整坐标轴文字。

pp3 <- pp2 + theme(axis.text.x=element_text(colour="black",angle = 0,vjust = 0.5,hjust = 0)) + 
  theme(axis.text.y=element_text(colour = "black")) + 
  theme(axis.line = element_line(colour = "black"))+ 
  theme(axis.title.y = element_text(colour = "black", margin = unit(c(0,1,0,1),"lines"))) 
pp3

pp4<-pp3+geom_vline(xintercept = c(0, 0), color = 'gray', size = 0.25) +
  geom_hline(yintercept = -log(0.05, 10), color = 'gray', size = 0.25) +
  labs(x = 'log2 Fold Change', y = '-log10 adjp-value', color = NA,title ="Ggt vs H") 

pp4
#输出图像。
ggsave("Ggt_H_phylum.pdf", width = 14, height = 10, units="cm")
ggsave("Ggt_H_phylum.tiff", width = 14, height = 10, units = "cm")
ggsave("Ggt_H_phylum.png", width = 14, height = 10, units = "cm")



#以Ggt为对照进行比较——Ggt vs H
res_Ggt_H <- results(dds, contrast = c('Group', 'H', 'Ggt'), pAdjustMethod = 'fdr', alpha = 0.05)
res_Ggt_H
#或者summary(res_Ggt_H)

#可以先按校正和 p 值由小到大排个序，方便查看
Ggt_H <- as.data.frame(res_Ggt_H[order(res_Ggt_H$padj), ])
Ggt_H
Ggt_H$OTU.ID <- rownames(Ggt_H)
write.table(Ggt_H, 'Ggt_H.txt', row.names = FALSE, sep = '\t', quote = FALSE)

#例如这里自定义根据 |log2FC| >= 1 和 adj.P.Val < 0.01 标记差异类型
#此处有关Up、Down的定义是根据对照和处理而言的，由于本人前期将对照和
Ggt_H[which(Ggt_H$padj %in% NA),"sig"] <- "None"
Ggt_H[which(Ggt_H$padj < 0.05 & Ggt_H$log2FoldChange <= -1),"sig"] <- "Down"
Ggt_H[which(Ggt_H$padj < 0.05 & Ggt_H$log2FoldChange >= 1),"sig"] <- "Up"
Ggt_H[which(Ggt_H$padj >= 0.05 | abs(Ggt_H$log2FoldChange) < 1),"sig"] <- "None"
Ggt_H

taxonomy <- read.delim("taxonomy2.txt", sep="\t")
taxonomy

Ggt_H_1 <- merge(Ggt_H, taxonomy, by="OTU.ID")
Ggt_H_1
Ggt_H_1 <- as.data.frame(Ggt_H_1[order(Ggt_H_1$padj), ])
Ggt_H_1
write.table(Ggt_H_1, 'Ggt_H_1.txt', row.names = FALSE, sep = '\t', quote = FALSE)


##标注上、下调
up <- subset(Ggt_H_1, sig == 'Up')
#up <- up[order(up$padj), ][1:15, ]
down <- subset(Ggt_H_1, sig == 'Down')
#down <- down[order(down$padj), ][1:10, ]


##门水平的注释
##对上下调的颜色进行标注
Ggt_H_1[which(Ggt_H_1$sig == "Up"),"color1"] <- up$phylum
Ggt_H_1[which(Ggt_H_1$sig == "Down" ),"color1"] <- down$phylum

color1 <- gsub(".*(p__)","",Ggt_H_1$color1)
color1
#gsub()可以用于字段的删减、增补、替换和切割，可以处理一个字段也可以处理由字段组成的向量。
?gsub
#删除Taxon中有关Bacteria;的字样
#taxon <- gsub(".*(Bacteria;)","",taxon)
#删除taxon中的Other字样
#taxon <- gsub(".*(Other)","",taxon)

color1 <- gsub("(\\[)","",color1)
color1 <- gsub("(\\])","",color1)
#将新产生的表头与原来的表头进行替换
Ggt_H_1$color1 <- color1
Ggt_H_1


##纲水平的注释
##对上下调的颜色进行标注
Ggt_H_1[which(Ggt_H_1$sig == "Up"),"color2"] <- up$class
Ggt_H_1[which(Ggt_H_1$sig == "Down" ),"color2"] <- down$class

color2 <- gsub(".*(c__)","",Ggt_H_1$color2)
color2
#gsub()可以用于字段的删减、增补、替换和切割，可以处理一个字段也可以处理由字段组成的向量。
?gsub
#删除Taxon中有关Bacteria;的字样
#taxon <- gsub(".*(Bacteria;)","",taxon)
#删除taxon中的Other字样
#taxon <- gsub(".*(Other)","",taxon)

color2 <- gsub("(\\[)","",color2)
color2 <- gsub("(\\])","",color2)
#将新产生的表头与原来的表头进行替换
Ggt_H_1$color2 <- color2
Ggt_H_1


##目水平的注释
##对上下调的颜色进行标注
Ggt_H_1[which(Ggt_H_1$sig == "Up"),"color3"] <- up$order
Ggt_H_1[which(Ggt_H_1$sig == "Down" ),"color3"] <- down$order

color3 <- gsub(".*(o__)","",Ggt_H_1$color3)
color3
#gsub()可以用于字段的删减、增补、替换和切割，可以处理一个字段也可以处理由字段组成的向量。
?gsub
#删除Taxon中有关Bacteria;的字样
#taxon <- gsub(".*(Bacteria;)","",taxon)
#删除taxon中的Other字样
#taxon <- gsub(".*(Other)","",taxon)

color3 <- gsub("(\\[)","",color3)
color3 <- gsub("(\\])","",color3)
#将新产生的表头与原来的表头进行替换
Ggt_H_1$color3 <- color3
Ggt_H_1

##科水平的注释
##对上下调的颜色进行标注
Ggt_H_1[which(Ggt_H_1$sig == "Up"),"color4"] <- up$family
Ggt_H_1[which(Ggt_H_1$sig == "Down" ),"color4"] <- down$family

color4 <- gsub(".*(f__)","",Ggt_H_1$color4)
color4
#gsub()可以用于字段的删减、增补、替换和切割，可以处理一个字段也可以处理由字段组成的向量。
?gsub
#删除Taxon中有关Bacteria;的字样
#taxon <- gsub(".*(Bacteria;)","",taxon)
#删除taxon中的Other字样
#taxon <- gsub(".*(Other)","",taxon)

color4 <- gsub("(\\[)","",color4)
color4 <- gsub("(\\])","",color4)
#将新产生的表头与原来的表头进行替换
Ggt_H_1$color4 <- color4
Ggt_H_1

##属水平的注释
##对上下调的颜色进行标注
Ggt_H_1[which(Ggt_H_1$sig == "Up"),"color5"] <- up$genus
Ggt_H_1[which(Ggt_H_1$sig == "Down" ),"color5"] <- down$genus

color5 <- gsub(".*(g__)","",Ggt_H_1$color5)
color5
#gsub()可以用于字段的删减、增补、替换和切割，可以处理一个字段也可以处理由字段组成的向量。
?gsub
#删除Taxon中有关Bacteria;的字样
#taxon <- gsub(".*(Bacteria;)","",taxon)
#删除taxon中的Other字样
#taxon <- gsub(".*(Other)","",taxon)

color5 <- gsub("(\\[)","",color5)
color5 <- gsub("(\\])","",color5)
#将新产生的表头与原来的表头进行替换
Ggt_H_1$color5 <- color5
Ggt_H_1


write.table(Ggt_H_1, 'Ggt_H_2.txt', row.names = FALSE, sep = '\t', quote = FALSE)


##第四步，比较结果出图
#门水平的注释

#横轴 log2FoldChange，纵轴 -log10(adj.P.Val)，颜色表示差异
pp1 <- ggplot(Ggt_H_1, aes(x = log2FoldChange, y = -log10(padj)),col=mycolor) +
  geom_point(aes(shape = sig, color=color1), alpha = 0.6, size = 2) +
  scale_shape_manual(values = c(17, 19, 15),limits = c("Up","None","Down" ))
theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent')) +
  theme(legend.title = element_blank(), legend.key = element_rect(fill = 'transparent'), legend.background = element_rect(fill = 'transparent'), legend.position = c(0.9, 0.9))

pp1


#调整绘图区主题。
pp2 <- pp1 + theme_bw()+ theme(panel.grid=element_blank()) + 
  theme(panel.border = element_blank()) + 
  xlim(-5, 5)+ylim(0, 16)+
  theme(panel.background=element_rect(fill='transparent',color='black'),plot.margin = unit(c(3,5,1,1),"mm"))
pp2

#调整坐标轴文字。

pp3 <- pp2 + theme(axis.text.x=element_text(colour="black",angle = 0,vjust = 0.5,hjust = 0)) + 
  theme(axis.text.y=element_text(colour = "black")) + 
  theme(axis.line = element_line(colour = "black"))+ 
  theme(axis.title.y = element_text(colour = "black", margin = unit(c(0,1,0,1),"lines"))) 
pp3

pp4<-pp3+geom_vline(xintercept = c(0, 0), color = 'gray', size = 0.25) +
  geom_hline(yintercept = -log(0.05, 10), color = 'gray', size = 0.25) +
  labs(x = 'log2 Fold Change', y = '-log10 adjp-value', color = NA,title ="Ggt vs H") 

pp4
#输出图像。
ggsave("Ggt_H_phylum.pdf", width = 14, height = 10, units="cm")
ggsave("Ggt_H_phylum.tiff", width = 14, height = 10, units = "cm")
ggsave("Ggt_H_phylum.png", width = 14, height = 10, units = "cm")


####H VS Ggt####

#以H为对照进行比较——H vs Ggt
res_H_Ggt<- results(dds, contrast = c('Group', 'Ggt', 'H'), pAdjustMethod = 'fdr', alpha = 0.05)
res_H_Ggt
#或者summary(res_Ggt_H)

#可以先按校正和 p 值由小到大排个序，方便查看
H_Ggt <- as.data.frame(res_H_Ggt[order(res_H_Ggt$padj), ])
H_Ggt
H_Ggt$OTU.ID <- rownames(H_Ggt)
write.table(H_Ggt, 'H_Ggt.txt', row.names = FALSE, sep = '\t', quote = FALSE)

#例如这里自定义根据 |log2FC| >= 1 和 adj.P.Val < 0.01 标记差异类型
#此处有关Up、Down的定义是根据对照和处理而言的，由于本人前期将对照和
H_Ggt[which(H_Ggt$padj %in% NA),"sig"] <- "None"
H_Ggt[which(H_Ggt$padj < 0.05 & H_Ggt$log2FoldChange <= -1),"sig"] <- "Down"
H_Ggt[which(H_Ggt$padj < 0.05 & H_Ggt$log2FoldChange >= 1),"sig"] <- "Up"
H_Ggt[which(H_Ggt$padj >= 0.05 | abs(H_Ggt$log2FoldChange) < 1),"sig"] <- "None"
H_Ggt

taxonomy <- read.delim("taxonomy2.txt", sep="\t")
taxonomy

H_Ggt_1 <- merge(H_Ggt, taxonomy, by="OTU.ID")
H_Ggt_1
H_Ggt_1 <- as.data.frame(H_Ggt_1[order(H_Ggt_1$padj), ])
H_Ggt_1
write.table(H_Ggt_1, 'H_Ggt_1.txt', row.names = FALSE, sep = '\t', quote = FALSE)


##标注上、下调
up <- subset(H_Ggt_1, sig == 'Up')
#up <- up[order(up$padj), ][1:15, ]
down <- subset(H_Ggt_1, sig == 'Down')
#down <- down[order(down$padj), ][1:10, ]


##门水平的注释
##对上下调的颜色进行标注
H_Ggt_1[which(H_Ggt_1$sig == "Up"),"color1"] <- up$phylum
H_Ggt_1[which(H_Ggt_1$sig == "Down" ),"color1"] <- down$phylum

color1 <- gsub(".*(p__)","",H_Ggt_1$color1)
color1
#gsub()可以用于字段的删减、增补、替换和切割，可以处理一个字段也可以处理由字段组成的向量。
?gsub
#删除Taxon中有关Bacteria;的字样
#taxon <- gsub(".*(Bacteria;)","",taxon)
#删除taxon中的Other字样
#taxon <- gsub(".*(Other)","",taxon)

color1 <- gsub("(\\[)","",color1)
color1 <- gsub("(\\])","",color1)
#将新产生的表头与原来的表头进行替换
H_Ggt_1$color1 <- color1
H_Ggt_1


##纲水平的注释
##对上下调的颜色进行标注
H_Ggt_1[which(H_Ggt_1$sig == "Up"),"color2"] <- up$class
H_Ggt_1[which(H_Ggt_1$sig == "Down" ),"color2"] <- down$class

color2 <- gsub(".*(c__)","",H_Ggt_1$color2)
color2
#gsub()可以用于字段的删减、增补、替换和切割，可以处理一个字段也可以处理由字段组成的向量。
?gsub
#删除Taxon中有关Bacteria;的字样
#taxon <- gsub(".*(Bacteria;)","",taxon)
#删除taxon中的Other字样
#taxon <- gsub(".*(Other)","",taxon)

color2 <- gsub("(\\[)","",color2)
color2 <- gsub("(\\])","",color2)
#将新产生的表头与原来的表头进行替换
H_Ggt_1$color2 <- color2
H_Ggt_1


##目水平的注释
##对上下调的颜色进行标注
H_Ggt_1[which(H_Ggt_1$sig == "Up"),"color3"] <- up$order
H_Ggt_1[which(H_Ggt_1$sig == "Down" ),"color3"] <- down$order

color3 <- gsub(".*(o__)","",H_Ggt_1$color3)
color3
#gsub()可以用于字段的删减、增补、替换和切割，可以处理一个字段也可以处理由字段组成的向量。
?gsub
#删除Taxon中有关Bacteria;的字样
#taxon <- gsub(".*(Bacteria;)","",taxon)
#删除taxon中的Other字样
#taxon <- gsub(".*(Other)","",taxon)

color3 <- gsub("(\\[)","",color3)
color3 <- gsub("(\\])","",color3)
#将新产生的表头与原来的表头进行替换
H_Ggt_1$color3 <- color3
H_Ggt_1

##科水平的注释
##对上下调的颜色进行标注
H_Ggt_1[which(H_Ggt_1$sig == "Up"),"color4"] <- up$family
H_Ggt_1[which(H_Ggt_1$sig == "Down" ),"color4"] <- down$family

color4 <- gsub(".*(f__)","",H_Ggt_1$color4)
color4
#gsub()可以用于字段的删减、增补、替换和切割，可以处理一个字段也可以处理由字段组成的向量。
?gsub
#删除Taxon中有关Bacteria;的字样
#taxon <- gsub(".*(Bacteria;)","",taxon)
#删除taxon中的Other字样
#taxon <- gsub(".*(Other)","",taxon)

color4 <- gsub("(\\[)","",color4)
color4 <- gsub("(\\])","",color4)
#将新产生的表头与原来的表头进行替换
H_Ggt_1$color4 <- color4
H_Ggt_1

##属水平的注释
##对上下调的颜色进行标注
H_Ggt_1[which(H_Ggt_1$sig == "Up"),"color5"] <- up$genus
H_Ggt_1[which(H_Ggt_1$sig == "Down" ),"color5"] <- down$genus

color5 <- gsub(".*(g__)","",H_Ggt_1$color5)
color5
#gsub()可以用于字段的删减、增补、替换和切割，可以处理一个字段也可以处理由字段组成的向量。
?gsub
#删除Taxon中有关Bacteria;的字样
#taxon <- gsub(".*(Bacteria;)","",taxon)
#删除taxon中的Other字样
#taxon <- gsub(".*(Other)","",taxon)

color5 <- gsub("(\\[)","",color5)
color5 <- gsub("(\\])","",color5)
#将新产生的表头与原来的表头进行替换
H_Ggt_1$color5 <- color5
H_Ggt_1


write.table(H_Ggt_1, 'H_Ggt_2.txt', row.names = FALSE, sep = '\t', quote = FALSE)


##第四步，比较结果出图
#门水平的注释

#横轴 log2FoldChange，纵轴 -log10(adj.P.Val)，颜色表示差异
pp1 <- ggplot(H_Ggt_1, aes(x = log2FoldChange, y = -log10(padj)),col=mycolor) +
  geom_point(aes(shape = sig, color=color1), alpha = 0.6, size = 2) +
  scale_shape_manual(values = c(17, 19, 15),limits = c("Up","None","Down" ))
theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent')) +
  theme(legend.title = element_blank(), legend.key = element_rect(fill = 'transparent'), legend.background = element_rect(fill = 'transparent'), legend.position = c(0.9, 0.9))

pp1


#调整绘图区主题。
pp2 <- pp1 + theme_bw()+ theme(panel.grid=element_blank()) + 
  theme(panel.border = element_blank()) + 
  xlim(-5, 5)+ylim(0, 16)+
  theme(panel.background=element_rect(fill='transparent',color='black'),plot.margin = unit(c(3,5,1,1),"mm"))
pp2

#调整坐标轴文字。

pp3 <- pp2 + theme(axis.text.x=element_text(colour="black",angle = 0,vjust = 0.5,hjust = 0)) + 
  theme(axis.text.y=element_text(colour = "black")) + 
  theme(axis.line = element_line(colour = "black"))+ 
  theme(axis.title.y = element_text(colour = "black", margin = unit(c(0,1,0,1),"lines"))) 
pp3

pp4<-pp3+geom_vline(xintercept = c(0, 0), color = 'gray', size = 0.25) +
  geom_hline(yintercept = -log(0.05, 10), color = 'gray', size = 0.25) +
  labs(x = 'log2 Fold Change', y = '-log10 adjp-value', color = NA,title ="H vs Ggt") 

pp4
#输出图像。
ggsave("H_Ggt_phylum.pdf", width = 14, height = 10, units="cm")
ggsave("H_Ggt_phylum.tiff", width = 14, height = 10, units = "cm")
ggsave("H_Ggt_phylum.png", width = 14, height = 10, units = "cm")


####H vs SA####

#以H为对照进行比较——H vs SA
res_H_SA<- results(dds, contrast = c('Group', 'SA', 'H'), pAdjustMethod = 'fdr', alpha = 0.05)
res_H_SA
#或者summary(res_Ggt_H)

#可以先按校正和 p 值由小到大排个序，方便查看
H_SA <- as.data.frame(res_H_SA[order(res_H_SA$padj), ])
H_SA
H_SA$OTU.ID <- rownames(H_SA)
write.table(H_SA, 'H_SA.txt', row.names = FALSE, sep = '\t', quote = FALSE)

#例如这里自定义根据 |log2FC| >= 1 和 adj.P.Val < 0.01 标记差异类型
#此处有关Up、Down的定义是根据对照和处理而言的，由于本人前期将对照和
H_SA[which(H_SA$padj %in% NA),"sig"] <- "None"
H_SA[which(H_SA$padj < 0.05 & H_SA$log2FoldChange <= -1),"sig"] <- "Down"
H_SA[which(H_SA$padj < 0.05 & H_SA$log2FoldChange >= 1),"sig"] <- "Up"
H_SA[which(H_SA$padj >= 0.05 | abs(H_SA$log2FoldChange) < 1),"sig"] <- "None"
H_SA

taxonomy <- read.delim("taxonomy2.txt", sep="\t")
taxonomy

H_SA_1 <- merge(H_SA, taxonomy, by="OTU.ID")
H_SA_1
H_SA_1 <- as.data.frame(H_SA_1[order(H_SA_1$padj), ])
H_SA_1
write.table(H_SA_1, 'H_SA_1.txt', row.names = FALSE, sep = '\t', quote = FALSE)


##标注上、下调
up <- subset(H_SA_1, sig == 'Up')
#up <- up[order(up$padj), ][1:15, ]
down <- subset(H_SA_1, sig == 'Down')
#down <- down[order(down$padj), ][1:10, ]


##门水平的注释
##对上下调的颜色进行标注
H_SA_1[which(H_SA_1$sig == "Up"),"color1"] <- up$phylum
H_SA_1[which(H_SA_1$sig == "Down" ),"color1"] <- down$phylum

color1 <- gsub(".*(p__)","",H_SA_1$color1)
color1
#gsub()可以用于字段的删减、增补、替换和切割，可以处理一个字段也可以处理由字段组成的向量。
?gsub
#删除Taxon中有关Bacteria;的字样
#taxon <- gsub(".*(Bacteria;)","",taxon)
#删除taxon中的Other字样
#taxon <- gsub(".*(Other)","",taxon)

color1 <- gsub("(\\[)","",color1)
color1 <- gsub("(\\])","",color1)
#将新产生的表头与原来的表头进行替换
H_SA_1$color1 <- color1
H_SA_1


##纲水平的注释
##对上下调的颜色进行标注
H_SA_1[which(H_SA_1$sig == "Up"),"color2"] <- up$class
H_SA_1[which(H_SA_1$sig == "Down" ),"color2"] <- down$class

color2 <- gsub(".*(c__)","",H_SA_1$color2)
color2
#gsub()可以用于字段的删减、增补、替换和切割，可以处理一个字段也可以处理由字段组成的向量。
?gsub
#删除Taxon中有关Bacteria;的字样
#taxon <- gsub(".*(Bacteria;)","",taxon)
#删除taxon中的Other字样
#taxon <- gsub(".*(Other)","",taxon)

color2 <- gsub("(\\[)","",color2)
color2 <- gsub("(\\])","",color2)
#将新产生的表头与原来的表头进行替换
H_SA_1$color2 <- color2
H_SA_1


##目水平的注释
##对上下调的颜色进行标注
H_SA_1[which(H_SA_1$sig == "Up"),"color3"] <- up$order
H_SA_1[which(H_SA_1$sig == "Down" ),"color3"] <- down$order

color3 <- gsub(".*(o__)","",H_SA_1$color3)
color3
#gsub()可以用于字段的删减、增补、替换和切割，可以处理一个字段也可以处理由字段组成的向量。
?gsub
#删除Taxon中有关Bacteria;的字样
#taxon <- gsub(".*(Bacteria;)","",taxon)
#删除taxon中的Other字样
#taxon <- gsub(".*(Other)","",taxon)

color3 <- gsub("(\\[)","",color3)
color3 <- gsub("(\\])","",color3)
#将新产生的表头与原来的表头进行替换
H_SA_1$color3 <- color3
H_SA_1

##科水平的注释
##对上下调的颜色进行标注
H_SA_1[which(H_SA_1$sig == "Up"),"color4"] <- up$family
H_SA_1[which(H_SA_1$sig == "Down" ),"color4"] <- down$family

color4 <- gsub(".*(f__)","",H_SA_1$color4)
color4
#gsub()可以用于字段的删减、增补、替换和切割，可以处理一个字段也可以处理由字段组成的向量。
?gsub
#删除Taxon中有关Bacteria;的字样
#taxon <- gsub(".*(Bacteria;)","",taxon)
#删除taxon中的Other字样
#taxon <- gsub(".*(Other)","",taxon)

color4 <- gsub("(\\[)","",color4)
color4 <- gsub("(\\])","",color4)
#将新产生的表头与原来的表头进行替换
H_SA_1$color4 <- color4
H_SA_1

##属水平的注释
##对上下调的颜色进行标注
H_SA_1[which(H_SA_1$sig == "Up"),"color5"] <- up$genus
H_SA_1[which(H_SA_1$sig == "Down" ),"color5"] <- down$genus

color5 <- gsub(".*(g__)","",H_SA_1$color5)
color5
#gsub()可以用于字段的删减、增补、替换和切割，可以处理一个字段也可以处理由字段组成的向量。
?gsub
#删除Taxon中有关Bacteria;的字样
#taxon <- gsub(".*(Bacteria;)","",taxon)
#删除taxon中的Other字样
#taxon <- gsub(".*(Other)","",taxon)

color5 <- gsub("(\\[)","",color5)
color5 <- gsub("(\\])","",color5)
#将新产生的表头与原来的表头进行替换
H_SA_1$color5 <- color5
H_SA_1


write.table(H_SA_1, 'H_SA_2.txt', row.names = FALSE, sep = '\t', quote = FALSE)


##第四步，比较结果出图
#门水平的注释

#横轴 log2FoldChange，纵轴 -log10(adj.P.Val)，颜色表示差异
pp1 <- ggplot(H_SA_1, aes(x = log2FoldChange, y = -log10(padj)),col=mycolor) +
  geom_point(aes(shape = sig, color=color1), alpha = 0.6, size = 2) +
  scale_shape_manual(values = c(17, 19, 15),limits = c("Up","None","Down" ))
theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent')) +
  theme(legend.title = element_blank(), legend.key = element_rect(fill = 'transparent'), legend.background = element_rect(fill = 'transparent'), legend.position = c(0.9, 0.9))

pp1


#调整绘图区主题。
pp2 <- pp1 + theme_bw()+ theme(panel.grid=element_blank()) + 
  theme(panel.border = element_blank()) + 
  xlim(-4, 8)+ylim(0, 40)+
  theme(panel.background=element_rect(fill='transparent',color='black'),plot.margin = unit(c(3,5,1,1),"mm"))
pp2

#调整坐标轴文字。

pp3 <- pp2 + theme(axis.text.x=element_text(colour="black",angle = 0,vjust = 0.5,hjust = 0)) + 
  theme(axis.text.y=element_text(colour = "black")) + 
  theme(axis.line = element_line(colour = "black"))+ 
  theme(axis.title.y = element_text(colour = "black", margin = unit(c(0,1,0,1),"lines"))) 
pp3

pp4<-pp3+geom_vline(xintercept = c(0, 0), color = 'gray', size = 0.25) +
  geom_hline(yintercept = -log(0.05, 10), color = 'gray', size = 0.25) +
  labs(x = 'log2 Fold Change', y = '-log10 adjp-value', color = NA,title ="H vs SA") 

pp4
#输出图像。
ggsave("H_SA_phylum.pdf", width = 14, height = 10, units="cm")
ggsave("H_SA_phylum.tiff", width = 14, height = 10, units = "cm")
ggsave("H_SA_phylum.png", width = 14, height = 10, units = "cm")

####H vs S_135####

#以H为对照进行比较——H vs S_135
res_H_S_135<- results(dds, contrast = c('Group', 'S_135', 'H'), pAdjustMethod = 'fdr', alpha = 0.05)
res_H_S_135
#或者summary(res_Ggt_H)

#可以先按校正和 p 值由小到大排个序，方便查看
H_S_135 <- as.data.frame(res_H_S_135[order(res_H_S_135$padj), ])
H_S_135
H_S_135$OTU.ID <- rownames(H_S_135)
write.table(H_S_135, 'H_S_135.txt', row.names = FALSE, sep = '\t', quote = FALSE)

#例如这里自定义根据 |log2FC| >= 1 和 adj.P.Val < 0.01 标记差异类型
#此处有关Up、Down的定义是根据对照和处理而言的，由于本人前期将对照和
H_S_135[which(H_S_135$padj %in% NA),"sig"] <- "None"
H_S_135[which(H_S_135$padj < 0.05 & H_S_135$log2FoldChange <= -1),"sig"] <- "Down"
H_S_135[which(H_S_135$padj < 0.05 & H_S_135$log2FoldChange >= 1),"sig"] <- "Up"
H_S_135[which(H_S_135$padj >= 0.05 | abs(H_S_135$log2FoldChange) < 1),"sig"] <- "None"
H_S_135

taxonomy <- read.delim("taxonomy2.txt", sep="\t")
taxonomy

H_S_135_1 <- merge(H_S_135, taxonomy, by="OTU.ID")
H_S_135_1
H_S_135_1 <- as.data.frame(H_S_135_1[order(H_S_135_1$padj), ])
H_S_135_1
write.table(H_S_135_1, 'H_S_135_1.txt', row.names = FALSE, sep = '\t', quote = FALSE)


##标注上、下调
up <- subset(H_S_135_1, sig == 'Up')
#up <- up[order(up$padj), ][1:15, ]
down <- subset(H_S_135_1, sig == 'Down')
#down <- down[order(down$padj), ][1:10, ]


##门水平的注释
##对上下调的颜色进行标注
H_S_135_1[which(H_S_135_1$sig == "Up"),"color1"] <- up$phylum
H_S_135_1[which(H_S_135_1$sig == "Down" ),"color1"] <- down$phylum

color1 <- gsub(".*(p__)","",H_S_135_1$color1)
color1
#gsub()可以用于字段的删减、增补、替换和切割，可以处理一个字段也可以处理由字段组成的向量。
?gsub
#删除Taxon中有关Bacteria;的字样
#taxon <- gsub(".*(Bacteria;)","",taxon)
#删除taxon中的Other字样
#taxon <- gsub(".*(Other)","",taxon)

color1 <- gsub("(\\[)","",color1)
color1 <- gsub("(\\])","",color1)
#将新产生的表头与原来的表头进行替换
H_S_135_1$color1 <- color1
H_S_135_1


##纲水平的注释
##对上下调的颜色进行标注
H_S_135_1[which(H_S_135_1$sig == "Up"),"color2"] <- up$class
H_S_135_1[which(H_S_135_1$sig == "Down" ),"color2"] <- down$class

color2 <- gsub(".*(c__)","",H_S_135_1$color2)
color2
#gsub()可以用于字段的删减、增补、替换和切割，可以处理一个字段也可以处理由字段组成的向量。
?gsub
#删除Taxon中有关Bacteria;的字样
#taxon <- gsub(".*(Bacteria;)","",taxon)
#删除taxon中的Other字样
#taxon <- gsub(".*(Other)","",taxon)

color2 <- gsub("(\\[)","",color2)
color2 <- gsub("(\\])","",color2)
#将新产生的表头与原来的表头进行替换
H_S_135_1$color2 <- color2
H_S_135_1


##目水平的注释
##对上下调的颜色进行标注
H_S_135_1[which(H_S_135_1$sig == "Up"),"color3"] <- up$order
H_S_135_1[which(H_S_135_1$sig == "Down" ),"color3"] <- down$order

color3 <- gsub(".*(o__)","",H_S_135_1$color3)
color3
#gsub()可以用于字段的删减、增补、替换和切割，可以处理一个字段也可以处理由字段组成的向量。
?gsub
#删除Taxon中有关Bacteria;的字样
#taxon <- gsub(".*(Bacteria;)","",taxon)
#删除taxon中的Other字样
#taxon <- gsub(".*(Other)","",taxon)

color3 <- gsub("(\\[)","",color3)
color3 <- gsub("(\\])","",color3)
#将新产生的表头与原来的表头进行替换
H_S_135_1$color3 <- color3
H_S_135_1

##科水平的注释
##对上下调的颜色进行标注
H_S_135_1[which(H_S_135_1$sig == "Up"),"color4"] <- up$family
H_S_135_1[which(H_S_135_1$sig == "Down" ),"color4"] <- down$family

color4 <- gsub(".*(f__)","",H_S_135_1$color4)
color4
#gsub()可以用于字段的删减、增补、替换和切割，可以处理一个字段也可以处理由字段组成的向量。
?gsub
#删除Taxon中有关Bacteria;的字样
#taxon <- gsub(".*(Bacteria;)","",taxon)
#删除taxon中的Other字样
#taxon <- gsub(".*(Other)","",taxon)

color4 <- gsub("(\\[)","",color4)
color4 <- gsub("(\\])","",color4)
#将新产生的表头与原来的表头进行替换
H_S_135_1$color4 <- color4
H_S_135_1

##属水平的注释
##对上下调的颜色进行标注
H_S_135_1[which(H_S_135_1$sig == "Up"),"color5"] <- up$genus
H_S_135_1[which(H_S_135_1$sig == "Down" ),"color5"] <- down$genus

color5 <- gsub(".*(g__)","",H_S_135_1$color5)
color5
#gsub()可以用于字段的删减、增补、替换和切割，可以处理一个字段也可以处理由字段组成的向量。
?gsub
#删除Taxon中有关Bacteria;的字样
#taxon <- gsub(".*(Bacteria;)","",taxon)
#删除taxon中的Other字样
#taxon <- gsub(".*(Other)","",taxon)

color5 <- gsub("(\\[)","",color5)
color5 <- gsub("(\\])","",color5)
#将新产生的表头与原来的表头进行替换
H_S_135_1$color5 <- color5
H_S_135_1


write.table(H_S_135_1, 'H_S_135_2.txt', row.names = FALSE, sep = '\t', quote = FALSE)


##第四步，比较结果出图
#门水平的注释

#横轴 log2FoldChange，纵轴 -log10(adj.P.Val)，颜色表示差异
pp1 <- ggplot(H_S_135_1, aes(x = log2FoldChange, y = -log10(padj)),col=mycolor) +
  geom_point(aes(shape = sig, color=color1), alpha = 0.6, size = 2) +
  scale_shape_manual(values = c(17, 19, 15),limits = c("Up","None","Down" ))
theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent')) +
  theme(legend.title = element_blank(), legend.key = element_rect(fill = 'transparent'), legend.background = element_rect(fill = 'transparent'), legend.position = c(0.9, 0.9))

pp1


#调整绘图区主题。
pp2 <- pp1 + theme_bw()+ theme(panel.grid=element_blank()) + 
  theme(panel.border = element_blank()) + 
  xlim(-4, 8)+ylim(0, 40)+
  theme(panel.background=element_rect(fill='transparent',color='black'),plot.margin = unit(c(3,5,1,1),"mm"))
pp2

#调整坐标轴文字。

pp3 <- pp2 + theme(axis.text.x=element_text(colour="black",angle = 0,vjust = 0.5,hjust = 0)) + 
  theme(axis.text.y=element_text(colour = "black")) + 
  theme(axis.line = element_line(colour = "black"))+ 
  theme(axis.title.y = element_text(colour = "black", margin = unit(c(0,1,0,1),"lines"))) 
pp3

pp4<-pp3+geom_vline(xintercept = c(0, 0), color = 'gray', size = 0.25) +
  geom_hline(yintercept = -log(0.05, 10), color = 'gray', size = 0.25) +
  labs(x = 'log2 Fold Change', y = '-log10 adjp-value', color = NA,title ="H vs S_135") 

pp4
#输出图像。
ggsave("H_S_135_phylum.pdf", width = 14, height = 10, units="cm")
ggsave("H_S_135_phylum.tiff", width = 14, height = 10, units = "cm")
ggsave("H_S_135_phylum.png", width = 14, height = 10, units = "cm")

####SA VS S_135####

#以H为对照进行比较——SA vs S_135
res_SA_S_135<- results(dds, contrast = c('Group', 'S_135', 'SA'), pAdjustMethod = 'fdr', alpha = 0.05)
res_SA_S_135
#或者summary(res_Ggt_H)

#可以先按校正和 p 值由小到大排个序，方便查看
SA_S_135 <- as.data.frame(res_SA_S_135[order(res_SA_S_135$padj), ])
SA_S_135
SA_S_135$OTU.ID <- rownames(SA_S_135)
write.table(SA_S_135, 'SA_S_135.txt', row.names = FALSE, sep = '\t', quote = FALSE)

#例如这里自定义根据 |log2FC| >= 1 和 adj.P.Val < 0.01 标记差异类型
#此处有关Up、Down的定义是根据对照和处理而言的，由于本人前期将对照和
SA_S_135[which(SA_S_135$padj %in% NA),"sig"] <- "None"
SA_S_135[which(SA_S_135$padj < 0.05 & SA_S_135$log2FoldChange <= -1),"sig"] <- "Down"
SA_S_135[which(SA_S_135$padj < 0.05 & SA_S_135$log2FoldChange >= 1),"sig"] <- "Up"
SA_S_135[which(SA_S_135$padj >= 0.05 | abs(SA_S_135$log2FoldChange) < 1),"sig"] <- "None"
SA_S_135

taxonomy <- read.delim("taxonomy2.txt", sep="\t")
taxonomy

SA_S_135_1 <- merge(SA_S_135, taxonomy, by="OTU.ID")
SA_S_135_1
SA_S_135_1 <- as.data.frame(SA_S_135_1[order(SA_S_135_1$padj), ])
SA_S_135_1
write.table(SA_S_135_1, 'SA_S_135_1.txt', row.names = FALSE, sep = '\t', quote = FALSE)


##标注上、下调
up <- subset(SA_S_135_1, sig == 'Up')
#up <- up[order(up$padj), ][1:15, ]
down <- subset(SA_S_135_1, sig == 'Down')
#down <- down[order(down$padj), ][1:10, ]


##门水平的注释
##对上下调的颜色进行标注
SA_S_135_1[which(SA_S_135_1$sig == "Up"),"color1"] <- up$phylum
SA_S_135_1[which(SA_S_135_1$sig == "Down" ),"color1"] <- down$phylum

color1 <- gsub(".*(p__)","",SA_S_135_1$color1)
color1
#gsub()可以用于字段的删减、增补、替换和切割，可以处理一个字段也可以处理由字段组成的向量。
?gsub
#删除Taxon中有关Bacteria;的字样
#taxon <- gsub(".*(Bacteria;)","",taxon)
#删除taxon中的Other字样
#taxon <- gsub(".*(Other)","",taxon)

color1 <- gsub("(\\[)","",color1)
color1 <- gsub("(\\])","",color1)
#将新产生的表头与原来的表头进行替换
SA_S_135_1$color1 <- color1
SA_S_135_1


##纲水平的注释
##对上下调的颜色进行标注
SA_S_135_1[which(SA_S_135_1$sig == "Up"),"color2"] <- up$class
SA_S_135_1[which(SA_S_135_1$sig == "Down" ),"color2"] <- down$class

color2 <- gsub(".*(c__)","",SA_S_135_1$color2)
color2
#gsub()可以用于字段的删减、增补、替换和切割，可以处理一个字段也可以处理由字段组成的向量。
?gsub
#删除Taxon中有关Bacteria;的字样
#taxon <- gsub(".*(Bacteria;)","",taxon)
#删除taxon中的Other字样
#taxon <- gsub(".*(Other)","",taxon)

color2 <- gsub("(\\[)","",color2)
color2 <- gsub("(\\])","",color2)
#将新产生的表头与原来的表头进行替换
SA_S_135_1$color2 <- color2
SA_S_135_1


##目水平的注释
##对上下调的颜色进行标注
SA_S_135_1[which(SA_S_135_1$sig == "Up"),"color3"] <- up$order
SA_S_135_1[which(SA_S_135_1$sig == "Down" ),"color3"] <- down$order

color3 <- gsub(".*(o__)","",SA_S_135_1$color3)
color3
#gsub()可以用于字段的删减、增补、替换和切割，可以处理一个字段也可以处理由字段组成的向量。
?gsub
#删除Taxon中有关Bacteria;的字样
#taxon <- gsub(".*(Bacteria;)","",taxon)
#删除taxon中的Other字样
#taxon <- gsub(".*(Other)","",taxon)

color3 <- gsub("(\\[)","",color3)
color3 <- gsub("(\\])","",color3)
#将新产生的表头与原来的表头进行替换
SA_S_135_1$color3 <- color3
SA_S_135_1

##科水平的注释
##对上下调的颜色进行标注
SA_S_135_1[which(SA_S_135_1$sig == "Up"),"color4"] <- up$family
SA_S_135_1[which(SA_S_135_1$sig == "Down" ),"color4"] <- down$family

color4 <- gsub(".*(f__)","",SA_S_135_1$color4)
color4
#gsub()可以用于字段的删减、增补、替换和切割，可以处理一个字段也可以处理由字段组成的向量。
?gsub
#删除Taxon中有关Bacteria;的字样
#taxon <- gsub(".*(Bacteria;)","",taxon)
#删除taxon中的Other字样
#taxon <- gsub(".*(Other)","",taxon)

color4 <- gsub("(\\[)","",color4)
color4 <- gsub("(\\])","",color4)
#将新产生的表头与原来的表头进行替换
SA_S_135_1$color4 <- color4
SA_S_135_1

##属水平的注释
##对上下调的颜色进行标注
SA_S_135_1[which(SA_S_135_1$sig == "Up"),"color5"] <- up$genus
SA_S_135_1[which(SA_S_135_1$sig == "Down" ),"color5"] <- down$genus

color5 <- gsub(".*(g__)","",SA_S_135_1$color5)
color5
#gsub()可以用于字段的删减、增补、替换和切割，可以处理一个字段也可以处理由字段组成的向量。
?gsub
#删除Taxon中有关Bacteria;的字样
#taxon <- gsub(".*(Bacteria;)","",taxon)
#删除taxon中的Other字样
#taxon <- gsub(".*(Other)","",taxon)

color5 <- gsub("(\\[)","",color5)
color5 <- gsub("(\\])","",color5)
#将新产生的表头与原来的表头进行替换
SA_S_135_1$color5 <- color5
SA_S_135_1


write.table(SA_S_135_1, 'SA_S_135_2.txt', row.names = FALSE, sep = '\t', quote = FALSE)


##第四步，比较结果出图
#门水平的注释

#横轴 log2FoldChange，纵轴 -log10(adj.P.Val)，颜色表示差异
pp1 <- ggplot(SA_S_135_1, aes(x = log2FoldChange, y = -log10(padj)),col=mycolor) +
  geom_point(aes(shape = sig, color=color1), alpha = 0.6, size = 2) +
  scale_shape_manual(values = c(17, 19, 15),limits = c("Up","None","Down" ))
theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent')) +
  theme(legend.title = element_blank(), legend.key = element_rect(fill = 'transparent'), legend.background = element_rect(fill = 'transparent'), legend.position = c(0.9, 0.9))

pp1


#调整绘图区主题。
pp2 <- pp1 + theme_bw()+ theme(panel.grid=element_blank()) + 
  theme(panel.border = element_blank()) + 
  xlim(-2, 2)+ylim(0, 1.5)+
  theme(panel.background=element_rect(fill='transparent',color='black'),plot.margin = unit(c(3,5,1,1),"mm"))
pp2

#调整坐标轴文字。

pp3 <- pp2 + theme(axis.text.x=element_text(colour="black",angle = 0,vjust = 0.5,hjust = 0)) + 
  theme(axis.text.y=element_text(colour = "black")) + 
  theme(axis.line = element_line(colour = "black"))+ 
  theme(axis.title.y = element_text(colour = "black", margin = unit(c(0,1,0,1),"lines"))) 
pp3

pp4<-pp3+geom_vline(xintercept = c(0, 0), color = 'gray', size = 0.25) +
  geom_hline(yintercept = -log(0.05, 10), color = 'gray', size = 0.25) +
  labs(x = 'log2 Fold Change', y = '-log10 adjp-value', color = NA,title ="SA vs S_135") 

pp4
#输出图像。
ggsave("SA_S_135_phylum.pdf", width = 14, height = 10, units="cm")
ggsave("SA_S_135_phylum.tiff", width = 14, height = 10, units = "cm")
ggsave("SA_S_135_phylum.png", width = 14, height = 10, units = "cm")


####功能分析####

######功能分析function analysis/prediction，包括Tax4Fun、FAPROTAX、PICRUST2三种细菌功能预测方法和FunGuide真菌功能预测
#此外，也可以基于BugBase进行细菌功能预测

setwd("E:/博士数据2/小麦微生物组数据/customer_backup/otus/function/")
library(reshape2)
library(multcomp)
library(ggplot2)
library(pheatmap)


##基于picrust2的功能分析
##采用scale_test对数据进行对数转换
picrust2 <- read.delim('./picrust2_heatmap2.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
picrust2
##首先进行对数转换
picrust2 <- apply(picrust2, 2, function(x) log2(x+1))
picrust2
##由于仅仅利用对数转换后的数值绘制出的热图颜色区分度不大，故对对数转换后的数值进一步求平方根
picrust2 <- apply(picrust2, 2, function(x) sqrt(x))
picrust2

pheatmap(picrust2)
pheatmap(picrust2,border_color="white",cellheight=15,cellwidth=15,
         display_numbers=FALSE,
         cluster_col=T, show_rownames=T,show_colnames=T,
         fontsize_number=7,width=10,height=8)

####添加行列注释信息
# 读取实验设计注释列分组
design<- read.delim('./mapping.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
design
annot_data1 <- data.frame(row.names=colnames(picrust2), group=design$Group)
annot_data1

#group <- c("#CC6666","#D0DFE6FF","#999999","#0073c2","black")
#names(group) <- c("Ggt","H","SA","S_135","BS")


pheatmap(picrust2,
         border_color=NA,
         cellheight=10,
         cellwidth=16,
         cluster_col= F, cluster_rows = T, show_rownames=T,show_colnames=F,treeheight_row=15,
         fontsize_number=7,
         width=10,height=8,
         frontsize=10,
         fontsize_row = 10,
         fontsize_col=10,
         annotation_col=annot_data1,
         angle_col = "270")

pheatmap(picrust2,
         border_color=NA,
         cellheight=10,
         cellwidth=16,
         cluster_col=F, show_rownames=T,show_colnames=F,treeheight_row=15,
         fontsize_number=7,
         width=10,height=8,
         frontsize=10,
         fontsize_row = 10,
         fontsize_col=10,
         annotation_col=annot_data1,
         angle_col = "270",
         filename="picrust2_heatmap.pdf")

pheatmap(picrust2,
         border_color=NA,
         cellheight=10,
         cellwidth=16,
         cluster_col=F, show_rownames=T,show_colnames=F,treeheight_row=15,
         fontsize_number=7,
         width=10,height=8,
         frontsize=10,
         fontsize_row = 10,
         fontsize_col=10,
         annotation_col=annot_data1,
         angle_col = "270",
         filename="picrust2_heatmap.png")

pheatmap(picrust2,
         border_color=NA,
         cellheight=10,
         cellwidth=16,
         cluster_col=F, show_rownames=T,show_colnames=F,treeheight_row=15,
         fontsize_number=7,
         width=10,height=8,
         frontsize=10,
         fontsize_row = 10,
         fontsize_col=10,
         annotation_col=annot_data1,
         angle_col = "270",
         filename="picrust2_heatmap.tiff")








##基于FAPROTAX的功能分析
##采用scale_test对数据进行对数转换
faprotax2 <- read.delim('./functional_table_faprotax2.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
faprotax2
##首先进行对数转换
faprotax2 <- apply(faprotax2, 2, function(x) log2(x+1))
faprotax2
##由于仅仅利用对数转换后的数值绘制出的热图颜色区分度不大，故对对数转换后的数值进一步求平方根
faprotax2 <- apply(faprotax2, 2, function(x) sqrt(x))
faprotax2

pheatmap(faprotax2)
pheatmap(faprotax2,border_color="white",cellheight=15,cellwidth=15,
         display_numbers=FALSE,
         cluster_col=T, show_rownames=T,show_colnames=T,
         fontsize_number=7,width=10,height=8)

####添加行列注释信息
# 读取实验设计注释列分组
design<- read.delim('./mapping.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
design
annot_data1 <- data.frame(row.names=colnames(faprotax2), group=design$Group)
annot_data1

#group <- c("#CC6666","#D0DFE6FF","#999999","#0073c2","black")
#names(group) <- c("Ggt","H","SA","S_135","BS")


pheatmap(faprotax2,
         border_color=NA,
         cellheight=10,
         cellwidth=16,
         cluster_col=F, cluster_rows = T, show_rownames=T,show_colnames=F,treeheight_row=15,
         fontsize_number=7,
         width=10,height=8,
         frontsize=10,
         fontsize_row = 10,
         fontsize_col=10,
         annotation_col=annot_data1,
         angle_col = "270")

pheatmap(faprotax2,
         border_color=NA,
         cellheight=10,
         cellwidth=16,
         cluster_col=F, show_rownames=T,show_colnames=F,treeheight_row=15,
         fontsize_number=7,
         width=10,height=8,
         frontsize=10,
         fontsize_row = 10,
         fontsize_col=10,
         annotation_col=annot_data1,
         angle_col = "270",
         filename="faprotax2_heatmap.pdf")

pheatmap(faprotax2,
         border_color=NA,
         cellheight=10,
         cellwidth=16,
         cluster_col=F, show_rownames=T,show_colnames=F,treeheight_row=15,
         fontsize_number=7,
         width=10,height=8,
         frontsize=10,
         fontsize_row = 10,
         fontsize_col=10,
         annotation_col=annot_data1,
         angle_col = "270",
         filename="faprotax2_heatmap.png")

pheatmap(faprotax2,
         border_color=NA,
         cellheight=10,
         cellwidth=16,
         cluster_col=F, show_rownames=T,show_colnames=F,treeheight_row=15,
         fontsize_number=7,
         width=10,height=8,
         frontsize=10,
         fontsize_row = 10,
         fontsize_col=10,
         annotation_col=annot_data1,
         angle_col = "270",
         filename="faprotax2_heatmap.tiff")



######网络分析，包括CoNet、SPIEC-EASI、分子生态网络分析（MENA）、SparCC

setwd("E:/博士数据2/小麦微生物组数据/customer_backup/otus/networks/")
#首先导入数据文件

#############基于丰度相关性的微生物共发生网络
##计算微生物丰度间的相关系数
library(Hmisc)

#门水平
phylum <- read.delim('phylum.txt', row.name = 1, check.names = FALSE)
phylum

#可选事先过滤一些低丰度或低频的类群
#例如只保留在所有样本中相对丰度总和高于 0.005 的门
phylum <- phylum[which(rowSums(phylum) >= 0.005), ]    
phylum

phylum1 <- phylum
phylum1[phylum1>0] <- 1
phylum <- phylum[which(rowSums(phylum1) >= 4), ]    #例如只保留在 4 个及以上样本中出现的属

#计算两属之间是否存在丰度变化的相关性，以 spearman 相关系数为例
phylum_corr <- rcorr(t(phylum), type = 'spearman')
phylum_corr
#阈值筛选
#将 spearman 相关系数低于 0.7 的关系剔除，即 r>=0.7
r <- phylum_corr$r
r[abs(r) < 0.7] <- 0

#选取显著性 p 值小于 0.05 的相关系数，即 p<0.05
p <- phylum_corr$P
p <- p.adjust(p, method = 'BH')    #可选 p 值校正，这里使用 BH 法校正 p 值
p[p>=0.05] <- -1
p[p<0.05 & p>=0] <- 1
p[p==-1] <- 0

#根据上述筛选的 r 值和 p 值保留数据
z <- r * p
diag(z) <- 0    #将相关矩阵中对角线中的值（代表了自相关）转为 0
head(z)[1:6,1:6]

#如此便得到了邻接矩阵格式的网络文件（微生物属的相关系数矩阵）
write.table(data.frame(z, check.names = FALSE), 'phylum_corr.matrix.txt', col.names = NA, sep = '\t', quote = FALSE)

##获得网络
library(igraph)

#将邻接矩阵转化为 igraph 网络的邻接列表
#构建含权的无向网络，权重代表了微生物属间丰度的 spearman 相关系数 
phylum <- graph.adjacency(z, weighted = TRUE, mode = 'undirected')
phylum

#自相关也可以通过该式去除
phylum <- simplify(phylum)

#孤立节点的删除（删除度为 0 的节点）
phylum <- delete.vertices(phylum, names(degree(phylum)[degree(phylum) == 0]))

#该模式下，边权重代表了相关系数
#由于权重通常为正值，因此最好取个绝对值，相关系数重新复制一列
E(phylum)$correlation <- E(phylum)$weight
E(phylum)$weight <- abs(E(phylum)$weight)

#为节点（微生物属）添加属性信息（界门纲目科属水平注释）
#“genus_taxonomy.txt” 记录了微生物的属性，读入该表后根据已知网络节点匹配对应的行
tax <- read.delim('taxonomy2.txt', row.name = 1, check.names = FALSE, stringsAsFactors = FALSE)
tax <- tax[as.character(V(phylum)$name), ]

#V(phylum)$kingdom <- tax$kingdom
V(phylum)$phylum <- tax$phylum
V(phylum)$class <- tax$class
V(phylum)$order <- tax$order
V(phylum)$family <- tax$family
V(phylum)$genus <- tax$genus

#查看网络图
phylum
plot(phylum)

##网络文件输出，输出特定的网络文件类型，便于后续数据分析需求
#邻接矩阵，出了上述提到的在计算相关系数后，输出筛选后的相关系数矩阵外
#还可以由 igraph 的邻接列表转换
adj_matrix <- as.matrix(get.adjacency(phylum, attr = 'correlation'))
write.table(data.frame(adj_matrix, check.names = FALSE), 'phylum.adj_matrix.txt', col.names = NA, sep = '\t', quote = FALSE)

#边列表
edge <- data.frame(as_edgelist(phylum))    #igraph 的邻接列表转为边列表

edge_list <- data.frame(
  source = edge[[1]],
  target = edge[[2]],
  weight = E(phylum)$weight,
  correlation = E(phylum)$correlation
)
head(edge_list)

write.table(edge_list, 'phylum.edge_list.txt', sep = '\t', row.names = FALSE, quote = FALSE)

#节点属性列表
node_list <- data.frame(
  label = names(V(phylum)),
  phylum = V(phylum)$phylum,
  class = V(phylum)$class,
  order = V(phylum)$order,
  family = V(phylum)$family,
  genus = V(phylum)$genus
)
head(node_list)

write.table(node_list, 'phylum.node_list.txt', sep = '\t', row.names = FALSE, quote = FALSE)

#边列表节点属性列表可以导入至 gephi 或 cytoscape 等网络可视化软件中进行编辑
#此外 igraph 也提供了可以被 gephi 或 cytoscape 等直接识别的格式
#graphml 格式，可使用 gephi 软件打开并进行可视化编辑
write.graph(phylum, 'phylum.graphml', format = 'graphml')

#gml 格式，可使用 cytoscape 软件打开并进行可视化编辑
write.graph(phylum, 'phylum.gml', format = 'gml')


#纲水平
class <- read.delim('class.txt', row.name = 1, check.names = FALSE)
class

#可选事先过滤一些低丰度或低频的类群
#例如只保留在所有样本中相对丰度总和高于 0.005 的纲，即选取了丰度大于0.5%的纲
class <- class[which(rowSums(class) >= 0.005), ]    
class

class1 <- class
class1[class1>0] <- 1
class1 <- class1[which(rowSums(class1) >= 4), ]    #例如只保留在 4 个及以上样本中出现的属

#计算两属之间是否存在丰度变化的相关性，以 spearman 相关系数为例
class_corr <- rcorr(t(class1), type = 'spearman')
class_corr
#阈值筛选
#将 spearman 相关系数低于 0.7 的关系剔除，即 r>=0.7
r <- class_corr$r
r[abs(r) < 0.7] <- 0

#选取显著性 p 值小于 0.05 的相关系数，即 p<0.05
p <- class_corr$P
p <- p.adjust(p, method = 'BH')    #可选 p 值校正，这里使用 BH 法校正 p 值
p[p>=0.05] <- -1
p[p<0.05 & p>=0] <- 1
p[p==-1] <- 0

#根据上述筛选的 r 值和 p 值保留数据
z <- r * p
diag(z) <- 0    #将相关矩阵中对角线中的值（代表了自相关）转为 0
head(z)[1:6,1:6]

#如此便得到了邻接矩阵格式的网络文件（微生物属的相关系数矩阵）
write.table(data.frame(z, check.names = FALSE), 'class_corr.matrix.txt', col.names = NA, sep = '\t', quote = FALSE)

##获得网络
library(igraph)

#将邻接矩阵转化为 igraph 网络的邻接列表
#构建含权的无向网络，权重代表了微生物属间丰度的 spearman 相关系数 
class <- graph.adjacency(z, weighted = TRUE, mode = 'undirected')
class

#自相关也可以通过该式去除
class <- simplify(class)

#孤立节点的删除（删除度为 0 的节点）
class <- delete.vertices(class, names(degree(class)[degree(class) == 0]))

#该模式下，边权重代表了相关系数
#由于权重通常为正值，因此最好取个绝对值，相关系数重新复制一列
E(class)$correlation <- E(class)$weight
E(class)$weight <- abs(E(class)$weight)

#为节点（微生物属）添加属性信息（界门纲目科属水平注释）
#“genus_taxonomy.txt” 记录了微生物的属性，读入该表后根据已知网络节点匹配对应的行
tax <- read.delim('taxonomy2.txt', row.name = 1, check.names = FALSE, stringsAsFactors = FALSE)
tax <- tax[as.character(V(class)$name), ]

#V(phylum)$kingdom <- tax$kingdom
V(class)$phylum <- tax$phylum
V(class)$class <- tax$class
V(class)$order <- tax$order
V(class)$family <- tax$family
V(class)$genus <- tax$genus

#查看网络图
class
plot(class)

##网络文件输出，输出特定的网络文件类型，便于后续数据分析需求
#邻接矩阵，出了上述提到的在计算相关系数后，输出筛选后的相关系数矩阵外
#还可以由 igraph 的邻接列表转换
adj_matrix <- as.matrix(get.adjacency(class, attr = 'correlation'))
write.table(data.frame(adj_matrix, check.names = FALSE), 'class.adj_matrix.txt', col.names = NA, sep = '\t', quote = FALSE)

#边列表
edge <- data.frame(as_edgelist(class))    #igraph 的邻接列表转为边列表

edge_list <- data.frame(
  source = edge[[1]],
  target = edge[[2]],
  weight = E(class)$weight,
  correlation = E(class)$correlation
)
head(edge_list)

write.table(edge_list, 'class.edge_list.txt', sep = '\t', row.names = FALSE, quote = FALSE)

#节点属性列表
node_list <- data.frame(
  label = names(V(class)),
  phylum = V(class)$phylum,
  class = V(class)$class,
  order = V(class)$order,
  family = V(class)$family,
  genus = V(class)$genus
)
head(node_list)

write.table(node_list, 'class.node_list.txt', sep = '\t', row.names = FALSE, quote = FALSE)

#边列表节点属性列表可以导入至 gephi 或 cytoscape 等网络可视化软件中进行编辑
#此外 igraph 也提供了可以被 gephi 或 cytoscape 等直接识别的格式
#graphml 格式，可使用 gephi 软件打开并进行可视化编辑
write.graph(class, 'class.graphml', format = 'graphml')

#gml 格式，可使用 cytoscape 软件打开并进行可视化编辑
write.graph(class, 'class.gml', format = 'gml')



#目水平
order <- read.delim('order.txt', row.name = 1, check.names = FALSE)
order

#可选事先过滤一些低丰度或低频的类群
#例如只保留在所有样本中相对丰度总和高于 0.005 的门
order <- order[which(rowSums(order) >= 0.005), ]    
order

order1 <- order
order1[order1>0] <- 1
order1 <- order1[which(rowSums(order1) >= 4), ]    #例如只保留在 4 个及以上样本中出现的属

#计算两属之间是否存在丰度变化的相关性，以 spearman 相关系数为例
order_corr <- rcorr(t(order1), type = 'spearman')
order_corr
#阈值筛选
#将 spearman 相关系数低于 0.7 的关系剔除，即 r>=0.7
r <- order_corr$r
r[abs(r) < 0.7] <- 0

#选取显著性 p 值小于 0.05 的相关系数，即 p<0.05
p <- order_corr$P
p <- p.adjust(p, method = 'BH')    #可选 p 值校正，这里使用 BH 法校正 p 值
p[p>=0.05] <- -1
p[p<0.05 & p>=0] <- 1
p[p==-1] <- 0

#根据上述筛选的 r 值和 p 值保留数据
z <- r * p
diag(z) <- 0    #将相关矩阵中对角线中的值（代表了自相关）转为 0
head(z)[1:6,1:6]

#如此便得到了邻接矩阵格式的网络文件（微生物属的相关系数矩阵）
write.table(data.frame(z, check.names = FALSE), 'order_corr.matrix.txt', col.names = NA, sep = '\t', quote = FALSE)

##获得网络
library(igraph)

#将邻接矩阵转化为 igraph 网络的邻接列表
#构建含权的无向网络，权重代表了微生物属间丰度的 spearman 相关系数 
order <- graph.adjacency(z, weighted = TRUE, mode = 'undirected')
order

#自相关也可以通过该式去除
order <- simplify(order)

#孤立节点的删除（删除度为 0 的节点）
order <- delete.vertices(order, names(degree(order)[degree(order) == 0]))

#该模式下，边权重代表了相关系数
#由于权重通常为正值，因此最好取个绝对值，相关系数重新复制一列
E(order)$correlation <- E(order)$weight
E(order)$weight <- abs(E(order)$weight)

#为节点（微生物属）添加属性信息（界门纲目科属水平注释）
#“genus_taxonomy.txt” 记录了微生物的属性，读入该表后根据已知网络节点匹配对应的行
tax <- read.delim('taxonomy2.txt', row.name = 1, check.names = FALSE, stringsAsFactors = FALSE)
tax <- tax[as.character(V(order)$name), ]

#V(phylum)$kingdom <- tax$kingdom
V(order)$phylum <- tax$phylum
V(order)$class <- tax$class
V(order)$order <- tax$order
V(order)$family <- tax$family
V(order)$genus <- tax$genus

#查看网络图
order
plot(order)

##网络文件输出，输出特定的网络文件类型，便于后续数据分析需求
#邻接矩阵，出了上述提到的在计算相关系数后，输出筛选后的相关系数矩阵外
#还可以由 igraph 的邻接列表转换
adj_matrix <- as.matrix(get.adjacency(order, attr = 'correlation'))
write.table(data.frame(adj_matrix, check.names = FALSE), 'order.adj_matrix.txt', col.names = NA, sep = '\t', quote = FALSE)

#边列表
edge <- data.frame(as_edgelist(order))    #igraph 的邻接列表转为边列表

edge_list <- data.frame(
  source = edge[[1]],
  target = edge[[2]],
  weight = E(order)$weight,
  correlation = E(order)$correlation
)
head(edge_list)

write.table(edge_list, 'order.edge_list.txt', sep = '\t', row.names = FALSE, quote = FALSE)

#节点属性列表
node_list <- data.frame(
  label = names(V(order)),
  phylum = V(order)$phylum,
  class = V(order)$class,
  order = V(order)$order,
  family = V(order)$family,
  genus = V(order)$genus
)
head(node_list)

write.table(node_list, 'order.node_list.txt', sep = '\t', row.names = FALSE, quote = FALSE)

#边列表节点属性列表可以导入至 gephi 或 cytoscape 等网络可视化软件中进行编辑
#此外 igraph 也提供了可以被 gephi 或 cytoscape 等直接识别的格式
#graphml 格式，可使用 gephi 软件打开并进行可视化编辑
write.graph(order, 'order.graphml', format = 'graphml')

#gml 格式，可使用 cytoscape 软件打开并进行可视化编辑
write.graph(order, 'order.gml', format = 'gml')



#科水平
family <- read.delim('family.txt', row.name = 1, check.names = FALSE)
family

#可选事先过滤一些低丰度或低频的类群
#例如只保留在所有样本中相对丰度总和高于 0.005 的科，即某个科在所有样本中的丰度在0.5%以上
family <- family[which(rowSums(family) >= 0.005), ]    
family

family1 <- family
family1[family1>0] <- 1
family <- family[which(rowSums(family1) >= 4), ]    #例如只保留在 4 个及以上样本中出现的属
family

#计算两属之间是否存在丰度变化的相关性，以 spearman 相关系数为例
family_corr <- rcorr(t(family), type = 'spearman')
family_corr
#阈值筛选
#将 spearman 相关系数低于 0.7 的关系剔除，即 r>=0.7
r <- family_corr$r
r[abs(r) < 0.7] <- 0

#选取显著性 p 值小于 0.05 的相关系数，即 p<0.05
p <- family_corr$P
p <- p.adjust(p, method = 'BH')    #可选 p 值校正，这里使用 BH 法校正 p 值
p[p>=0.05] <- -1
p[p<0.05 & p>=0] <- 1
p[p==-1] <- 0

#根据上述筛选的 r 值和 p 值保留数据
z <- r * p
diag(z) <- 0    #将相关矩阵中对角线中的值（代表了自相关）转为 0
head(z)[1:6,1:6]

#如此便得到了邻接矩阵格式的网络文件（微生物属的相关系数矩阵）
write.table(data.frame(z, check.names = FALSE), 'family_corr.matrix.txt', col.names = NA, sep = '\t', quote = FALSE)

##获得网络
library(igraph)

#将邻接矩阵转化为 igraph 网络的邻接列表
#构建含权的无向网络，权重代表了微生物属间丰度的 spearman 相关系数 
family <- graph.adjacency(z, weighted = TRUE, mode = 'undirected')
family

#自相关也可以通过该式去除
family <- simplify(family)

#孤立节点的删除（删除度为 0 的节点）
family <- delete.vertices(family, names(degree(family)[degree(family) == 0]))

#该模式下，边权重代表了相关系数
#由于权重通常为正值，因此最好取个绝对值，相关系数重新复制一列
E(family)$correlation <- E(family)$weight
E(family)$weight <- abs(E(family)$weight)

#为节点（微生物属）添加属性信息（界门纲目科属水平注释）
#“genus_taxonomy.txt” 记录了微生物的属性，读入该表后根据已知网络节点匹配对应的行
tax <- read.delim('taxonomy2.txt', row.name = 1, check.names = FALSE, stringsAsFactors = FALSE)
tax <- tax[as.character(V(family)$name), ]

#V(phylum)$kingdom <- tax$kingdom
V(family)$phylum <- tax$phylum
V(family)$class <- tax$class
V(family)$order <- tax$order
V(family)$family <- tax$family
V(family)$genus <- tax$genus

#查看网络图
order
plot(family)

##网络文件输出，输出特定的网络文件类型，便于后续数据分析需求
#邻接矩阵，出了上述提到的在计算相关系数后，输出筛选后的相关系数矩阵外
#还可以由 igraph 的邻接列表转换
adj_matrix <- as.matrix(get.adjacency(family, attr = 'correlation'))
write.table(data.frame(adj_matrix, check.names = FALSE), 'family.adj_matrix.txt', col.names = NA, sep = '\t', quote = FALSE)

#边列表
edge <- data.frame(as_edgelist(family))    #igraph 的邻接列表转为边列表

edge_list <- data.frame(
  source = edge[[1]],
  target = edge[[2]],
  weight = E(family)$weight,
  correlation = E(family)$correlation
)
head(edge_list)

write.table(edge_list, 'family.edge_list.txt', sep = '\t', row.names = FALSE, quote = FALSE)

#节点属性列表
node_list <- data.frame(
  label = names(V(family)),
  phylum = V(family)$phylum,
  class = V(family)$class,
  order = V(family)$order,
  family = V(family)$family,
  genus = V(family)$genus
)
head(node_list)

write.table(node_list, 'family.node_list.txt', sep = '\t', row.names = FALSE, quote = FALSE)

#边列表节点属性列表可以导入至 gephi 或 cytoscape 等网络可视化软件中进行编辑
#此外 igraph 也提供了可以被 gephi 或 cytoscape 等直接识别的格式
#graphml 格式，可使用 gephi 软件打开并进行可视化编辑
write.graph(family, 'family.graphml', format = 'graphml')

#gml 格式，可使用 cytoscape 软件打开并进行可视化编辑
write.graph(family, 'family.gml', format = 'gml')





#属水平
genus <- read.delim('genus.txt', row.name = 1, check.names = FALSE)
genus

#可选事先过滤一些低丰度或低频的类群
#例如只保留在所有样本中相对丰度总和高于 0.005 的属
genus <- genus[which(rowSums(genus) >= 0.005), ]    
genus

genus1 <- genus
genus1[genus1>0] <- 1
genus1 <- genus1[which(rowSums(genus1) >= 4), ]    #例如只保留在 4 个及以上样本中出现的属

#计算两属之间是否存在丰度变化的相关性，以 spearman 相关系数为例
genus_corr <- rcorr(t(genus1), type = 'spearman')
genus_corr
#阈值筛选
#将 spearman 相关系数低于 0.7 的关系剔除，即 r>=0.7
r <- genus_corr$r
r[abs(r) < 0.7] <- 0

#选取显著性 p 值小于 0.05 的相关系数，即 p<0.05
p <- genus_corr$P
p <- p.adjust(p, method = 'BH')    #可选 p 值校正，这里使用 BH 法校正 p 值
p[p>=0.05] <- -1
p[p<0.05 & p>=0] <- 1
p[p==-1] <- 0

#根据上述筛选的 r 值和 p 值保留数据
z <- r * p
diag(z) <- 0    #将相关矩阵中对角线中的值（代表了自相关）转为 0
head(z)[1:6,1:6]

#如此便得到了邻接矩阵格式的网络文件（微生物属的相关系数矩阵）
write.table(data.frame(z, check.names = FALSE), 'genus_corr.matrix.txt', col.names = NA, sep = '\t', quote = FALSE)

##获得网络
library(igraph)

#将邻接矩阵转化为 igraph 网络的邻接列表
#构建含权的无向网络，权重代表了微生物属间丰度的 spearman 相关系数 
genus <- graph.adjacency(z, weighted = TRUE, mode = 'undirected')
genus

#自相关也可以通过该式去除
genus <- simplify(genus)

#孤立节点的删除（删除度为 0 的节点）
genus <- delete.vertices(genus, names(degree(genus)[degree(genus) == 0]))

#该模式下，边权重代表了相关系数
#由于权重通常为正值，因此最好取个绝对值，相关系数重新复制一列
E(genus)$correlation <- E(genus)$weight
E(genus)$weight <- abs(E(genus)$weight)

#为节点（微生物属）添加属性信息（界门纲目科属水平注释）
#“genus_taxonomy.txt” 记录了微生物的属性，读入该表后根据已知网络节点匹配对应的行
#tax <- read.delim('taxonomy2.txt', row.name = 1, check.names = FALSE, stringsAsFactors = FALSE)
#tax <- tax[as.character(V(genus)$name), ]

#V(genus)$kingdom <- tax$kingdom
#V(genus)$phylum <- tax$phylum
#V(genus)$class <- tax$class
#V(genus)$order <- tax$order
#V(genus)$family <- tax$family
#V(genus)$genus <- tax$genus

#查看网络图
genus
plot(genus)

##网络文件输出，输出特定的网络文件类型，便于后续数据分析需求
#邻接矩阵，出了上述提到的在计算相关系数后，输出筛选后的相关系数矩阵外
#还可以由 igraph 的邻接列表转换
adj_matrix <- as.matrix(get.adjacency(genus, attr = 'correlation'))
write.table(data.frame(adj_matrix, check.names = FALSE), 'genus.adj_matrix.txt', col.names = NA, sep = '\t', quote = FALSE)

#边列表
edge <- data.frame(as_edgelist(genus))    #igraph 的邻接列表转为边列表

edge_list <- data.frame(
  source = edge[[1]],
  target = edge[[2]],
  weight = E(genus)$weight,
  correlation = E(genus)$correlation
)
head(edge_list)

write.table(edge_list, 'genus.edge_list.txt', sep = '\t', row.names = FALSE, quote = FALSE)

#节点属性列表
node_list <- data.frame(
  label = names(V(genus)),
  kingdom = V(genus)$kingdom,
  phylum = V(genus)$phylum,
  class = V(genus)$class,
  order = V(genus)$order,
  family = V(genus)$family,
  genus = V(genus)$genus
)
head(node_list)

write.table(node_list, 'genus.node_list.txt', sep = '\t', row.names = FALSE, quote = FALSE)

#边列表节点属性列表可以导入至 gephi 或 cytoscape 等网络可视化软件中进行编辑
#此外 igraph 也提供了可以被 gephi 或 cytoscape 等直接识别的格式
#graphml 格式，可使用 gephi 软件打开并进行可视化编辑
write.graph(genus, 'genus.graphml', format = 'graphml')

#gml 格式，可使用 cytoscape 软件打开并进行可视化编辑
write.graph(genus, 'genus.gml', format = 'gml')




rm(list=ls())
#1%高丰度OTU水平,特别值得注意的是要将OTU表格注释到属水平即可
#############基于丰度相关性的微生物共发生网络
##计算微生物丰度间的相关系数
setwd("E:/博士数据2/小麦微生物组数据/customer_backup/otus/networks/")
library(Hmisc)

#以1%高丰度OTU水平丰度为例，“otu_AAT.txt” 是一个1%高丰度OTU水平的微生物丰度表
otu_AAT <- read.delim('otu_AAT.txt', row.name = 1, check.names = FALSE)
otu_AAT

#可选事先过滤一些低丰度或低频的类群
#otu_AAT <- otu_AAT[which(rowSums(otu_AAT) >= 0.005), ]    #例如只保留相对丰度总和高于 0.005 的属

#otu_AAT1 <- otu_AAT
#otu_AAT1[otu_AAT1>0] <- 1
#otu_AAT1 <- otu_AAT1[which(rowSums(otu_AAT1) >= 5), ]    #例如只保留在 5 个及以上样本中出现的属

#计算两属之间是否存在丰度变化的相关性，以 spearman 相关系数为例
otu_AAT_corr <- rcorr(t(otu_AAT), type = 'spearman')

#阈值筛选
#将 spearman 相关系数低于 0.7 的关系剔除，即 r>=0.7
r <- otu_AAT_corr$r
r[abs(r) < 0.7] <- 0

#选取显著性 p 值小于 0.01 的相关系数，即 p<0.01
p <- otu_AAT_corr$P
p <- p.adjust(p, method = 'BH')    #可选 p 值校正，这里使用 BH 法校正 p 值
p[p>=0.01] <- -1
p[p<0.01 & p>=0] <- 1
p[p==-1] <- 0

#根据上述筛选的 r 值和 p 值保留数据
z <- r * p
diag(z) <- 0    #将相关矩阵中对角线中的值（代表了自相关）转为 0
head(z)[1:6,1:6]

#如此便得到了邻接矩阵格式的网络文件（微生物属的相关系数矩阵）
write.table(data.frame(z, check.names = FALSE), 'otu_AAT_corr.matrix.txt', col.names = NA, sep = '\t', quote = FALSE)


##将微生物互作网络转换为0、1矩阵，数值“1”表示微生物 OTU 之间存在互作，“0”表示无互作
otu_AAT_corr.matrix2 <- read.delim('otu_AAT_corr.matrix.txt', row.name = 1, check.names = FALSE)
head(otu_AAT_corr.matrix2)[1:6,1:6]
otu_AAT_corr.matrix2[otu_AAT_corr.matrix2 > 0] <- 1
otu_AAT_corr.matrix2[otu_AAT_corr.matrix2 < 0] <- 1
head(otu_AAT_corr.matrix2)[1:6,1:6]
write.table(data.frame(otu_AAT_corr.matrix2, check.names = FALSE), 'otu_AAT_adjacency_unweight.txt', col.names = NA, sep = '\t', quote = FALSE)


##获得网络
library(igraph)

#将邻接矩阵转化为 igraph 网络的邻接列表
#构建含权的无向网络，权重代表了微生物属间丰度的 spearman 相关系数 
igraph <- graph.adjacency(z, weighted = TRUE, mode = 'undirected')
igraph

#自相关也可以通过该式去除
igraph <- simplify(igraph)

#孤立节点的删除（删除度为 0 的节点）
igraph <- delete.vertices(igraph, names(degree(igraph)[degree(igraph) == 0]))
igraph
#该模式下，边权重代表了相关系数
#由于权重通常为正值，因此最好取个绝对值，相关系数重新复制一列
E(igraph)$corr <- E(igraph)$weight
E(igraph)$weight <- abs(E(igraph)$weight)
E(igraph)$cor <- E(igraph)$corr
E(igraph)$cor[E(igraph)$cor > 0] <- 1
E(igraph)$cor[E(igraph)$cor < 0] <- -1

#为节点（微生物属）添加属性信息（界门纲目科属水平注释）
#“genus_taxonomy.txt” 记录了微生物的属性，读入该表后根据已知网络节点匹配对应的行
tax <- read.delim('taxonomy_AAT.txt', row.name = 1, check.names = FALSE, stringsAsFactors = FALSE)
tax
tax <- tax[as.character(V(igraph)$name), ]
tax


V(igraph)$phylum <- tax$phylum
V(igraph)$class <- tax$class
V(igraph)$order <- tax$order
V(igraph)$family <- tax$family
V(igraph)$genus <- tax$genus
V(igraph)$species <- tax$species

####R语言计算节点和边特征####
##节点特征
#节点数量
length(V(igraph)$name)
#或
vcount(igraph)


#节点度（Degree）
#由于本示例是个无向网络，故无出度和入度之分
V(igraph)$degree <- degree(igraph)
V(igraph)$degree

#查看度分布
#可观察到微生物相关网络通常服从幂律分布，这个下节再讲怎样通过计算验证
degree_dist <- degree.distribution(igraph)[-1]
degree_num <- 1:max(V(igraph)$degree)

par(mfrow = c(1, 2))
hist(V(igraph)$degree, xlab = 'Degree', ylab = 'Frequency', 
     main = 'Degree distribution')
plot(degree_num, degree_dist, log = 'xy', xlab = 'Log-degree', 
     ylab = 'Log-intensity', main = 'Log-log degree distribution')

#查看节点度与其“邻居”的平均度的关系
#微生物网络中高度值的节点更倾向连接在一起，是普遍现象吗？
neighbor_degree <- graph.knn(igraph, V(igraph))$knn
plot(V(igraph)$degree, neighbor_degree, log = 'xy', 
     xlab = 'Log degree', ylab = 'Log average neighbor degree')

#加权度（Weighted degree）
V(igraph)$weight_degree <- strength(igraph)
V(igraph)$weight_degree

#接近中心性（Closeness centrality）
V(igraph)$closeness_centrality <- closeness(igraph)
V(igraph)$closeness_centrality

#介数中心性（Betweenness centrality）
V(igraph)$betweenness_centrality <- betweenness(igraph)
V(igraph)$betweenness_centrality

#特征向量中心性（Eigenvector centrality）
V(igraph)$eigenvector_centrality <- evcent(igraph)$vector
V(igraph)$eigenvector_centrality

#探索三种描述节点中心性的特征的关系
library(car)

scatter3d(V(igraph)$closeness_centrality, V(igraph)$betweenness_centrality, V(igraph)$eigenvector_centrality, 
          xlab =  'Closeness centrality', ylab = 'Betweenness centrality', zlab = 'Eigenvector centrality', 
          surface = FALSE)

#探索节点度和节点中心性的关系，如与特征向量中心性的关系
plot(V(igraph)$degree, V(igraph)$eigenvector_centrality, 
     xlab = 'Degree', ylab = 'Eigenvector centrality')


#模块划分，详情 ?cluster_fast_greedy，有多种模型
set.seed(123)
V(igraph)$modularity <- membership(cluster_fast_greedy(igraph))


#输出列表
node_list <- data.frame(
  node_id = V(igraph)$name, 
  degree = V(igraph)$degree, 
  weight_degree = V(igraph)$weight_degree, 
  closeness_centrality = V(igraph)$closeness_centrality, 
  betweenness_centrality = V(igraph)$betweenness_centrality, 
  eigenvector_centrality = V(igraph)$eigenvector_centrality,
  modularity = V(igraph)$modularity
)

head(node_list)
write.table(node_list, 'otu_AAT_node_list.txt', sep = '\t', row.names = FALSE, quote = FALSE)



##边特征
#边的数量
ecount(igraph)

#权重（Weighted），已在数据读入时转化获得
E(igraph)$weight

#边介数中心性（Edge betweenness centrality）
E(igraph)$betweenness_centrality <- edge.betweenness(igraph)
E(igraph)$betweenness_centrality

#输出列表
edge <- data.frame(as_edgelist(igraph))    #igraph 的邻接列表转为边列表

edge_list <- data.frame(
  source = edge[[1]],
  target = edge[[2]],
  weight = E(igraph)$weight,
  correlation = E(igraph)$corr,
  cor= E(igraph)$cor,
  betweenness_centrality = E(igraph)$betweenness_centrality
)
head(edge_list)

write.table(edge_list, 'otu_AAT_edge_list.txt', sep = '\t', row.names = FALSE, quote = FALSE)



##子图与普查
#所有尺寸的团的普查可以提供一个快照，将显示各尺寸的团的数量
census <- table(sapply(cliques(igraph), length))
census
plot(census)

#k 核
cores <- graph.coreness(igraph)
cores
sna::gplot.target(adjacency_unweight, cores, usearrows = FALSE, vertex.col = cores)

#二元组（dyad）和三元组（triad）
dyad.census(simplify(igraph))
triad.census(simplify(igraph))

#节点数量（number of nodes）和边数量（number of edges）
nodes_num <- length(V(igraph))
nodes_num

edges_num <- length(E(igraph))
edges_num

#平均度（average degree）
average_degree <- mean(degree(igraph))
#或者，2x边数量/节点数量
average_degree <- 2*edges_num/nodes_num



####模块内连通度（Zi）和模块间连通度（Pi）识别关键节点####
##计算模块内连通度（Zi）和模块间连通度（Pi）
source('zi_pi.r')

#上述的邻接矩阵类型的网络文件
adjacency_unweight <- read.delim('otu_AAT_adjacency_unweight_2.txt', row.names = 1, sep = '\t', check.names = FALSE)
head(adjacency_unweight)[1:6,1:6]
#节点属性列表，包含节点所划分的模块
nodes_list <- read.delim('otu_AAT_node_list_2.txt', row.names = 1, sep = '\t', check.names = FALSE)
nodes_list
#两个文件的节点顺序要一致
nodes_list <- nodes_list[rownames(adjacency_unweight), ]
#计算模块内连通度（Zi）和模块间连通度（Pi）
#指定邻接矩阵、节点列表、节点列表中节点度和模块度的列名称
zi_pi <- zi.pi(nodes_list, adjacency_unweight, degree = 'degree', modularity_class = 'modularity')
head(zi_pi)

write.table(zi_pi, 'otu_AAT_zi_pi_result.txt', sep = '\t', row.names = FALSE, quote = FALSE)


##可再根据阈值对节点划分为 4 种类型，并作图展示其分布
library(ggplot2)

zi_pi <- na.omit(zi_pi)   #NA 值最好去掉，不要当 0 处理
zi_pi[which(zi_pi$within_module_connectivities < 2.5 & zi_pi$among_module_connectivities < 0.62),'type'] <- 'Peripherals'
zi_pi[which(zi_pi$within_module_connectivities < 2.5 & zi_pi$among_module_connectivities > 0.62),'type'] <- 'Connectors'
zi_pi[which(zi_pi$within_module_connectivities > 2.5 & zi_pi$among_module_connectivities < 0.62),'type'] <- 'Module hubs'
zi_pi[which(zi_pi$within_module_connectivities > 2.5 & zi_pi$among_module_connectivities > 0.62),'type'] <- 'Network hubs'

ggplot(zi_pi, aes(among_module_connectivities, within_module_connectivities)) +
  geom_point(aes(color = type), alpha = 0.5, size = 2) +
  scale_color_manual(values = c('gray','red','blue','purple'), 
                     limits = c('Peripherals', 'Connectors', 'Module hubs', 'Network hubs'))+
  theme(panel.grid = element_blank(), axis.line = element_line(colour = 'black'), 
        panel.background = element_blank(), legend.key = element_blank()) +
  labs(x = 'Among-module connectivities', y = 'Within-module connectivities', color = '') +
  geom_vline(xintercept = 0.62) +
  geom_hline(yintercept = 2.5)


#查看网络图
igraph
plot(igraph)


#边列表节点属性列表可以导入至 gephi 或 cytoscape 等网络可视化软件中进行编辑
#此外 igraph 也提供了可以被 gephi 或 cytoscape 等直接识别的格式
#graphml 格式，可使用 gephi 软件打开并进行可视化编辑
write.graph(igraph, 'otu_AAT_2.graphml', format = 'graphml')

#gml 格式，可使用 cytoscape 软件打开并进行可视化编辑
write.graph(igraph, 'otu_AAT_2.gml', format = 'gml')



#1%高丰度OTU水平,特别值得注意的是要将OTU表格注释到属水平即可
#############基于丰度相关性的微生物共发生网络
##计算微生物丰度间的相关系数
setwd("E:/博士数据2/小麦微生物组数据/customer_backup/otus/networks/")
library(Hmisc)

#以1%高丰度OTU水平丰度为例，“otu_AAT.txt” 是一个1%高丰度OTU水平的微生物丰度表
genus_AT <- read.delim('genus_AT.txt', row.name = 1, check.names = FALSE)
genus_AT

#可选事先过滤一些低丰度或低频的类群
#otu_AAT <- otu_AAT[which(rowSums(otu_AAT) >= 0.005), ]    #例如只保留相对丰度总和高于 0.005 的属

#otu_AAT1 <- otu_AAT
#otu_AAT1[otu_AAT1>0] <- 1
#otu_AAT1 <- otu_AAT1[which(rowSums(otu_AAT1) >= 5), ]    #例如只保留在 5 个及以上样本中出现的属

#计算两属之间是否存在丰度变化的相关性，以 spearman 相关系数为例
genus_AT_corr <- rcorr(t(genus_AT), type = 'spearman')

#阈值筛选
#将 spearman 相关系数低于 0.7 的关系剔除，即 r>=0.7
r <- genus_AT_corr$r
r[abs(r) < 0.7] <- 0

#选取显著性 p 值小于 0.05 的相关系数，即 p<0.05
p <- genus_AT_corr$P
p <- p.adjust(p, method = 'BH')    #可选 p 值校正，这里使用 BH 法校正 p 值
p[p>=0.05] <- -1
p[p<0.05 & p>=0] <- 1
p[p==-1] <- 0

#根据上述筛选的 r 值和 p 值保留数据
z <- r * p
diag(z) <- 0    #将相关矩阵中对角线中的值（代表了自相关）转为 0
head(z)[1:6,1:6]

#如此便得到了邻接矩阵格式的网络文件（微生物属的相关系数矩阵）
write.table(data.frame(z, check.names = FALSE), 'genus_AT_corr.matrix.txt', col.names = NA, sep = '\t', quote = FALSE)


##将微生物互作网络转换为0、1矩阵，数值“1”表示微生物 OTU 之间存在互作，“0”表示无互作
genus_AT_corr.matrix2 <- read.delim('genus_AT_corr.matrix.txt', row.name = 1, check.names = FALSE)
head(genus_AT_corr.matrix2)[1:6,1:6]
genus_AT_corr.matrix2[genus_AT_corr.matrix2 > 0] <- 1
genus_AT_corr.matrix2[genus_AT_corr.matrix2 < 0] <- 1
head(genus_AT_corr.matrix2)[1:6,1:6]
write.table(data.frame(genus_AT_corr.matrix2, check.names = FALSE), 'genus_AT_adjacency_unweight.txt', col.names = NA, sep = '\t', quote = FALSE)


##获得网络
library(igraph)


#将邻接矩阵转化为 igraph 网络的邻接列表
#构建含权的无向网络，权重代表了微生物属间丰度的 spearman 相关系数 
igraph <- graph.adjacency(z, weighted = TRUE, mode = 'undirected')
igraph
#自相关也可以通过该式去除
igraph <- simplify(igraph)

#孤立节点的删除（删除度为 0 的节点）
igraph <- delete.vertices(igraph, names(degree(igraph)[degree(igraph) == 0]))
igraph
#该模式下，边权重代表了相关系数
#由于权重通常为正值，因此最好取个绝对值，相关系数重新复制一列
E(igraph)$correlation <- E(igraph)$weight
E(igraph)$weight <- abs(E(igraph)$weight)
E(igraph)$cor <- E(igraph)$correlation 
E(igraph)$cor[E(igraph)$cor > 0] <- 1
E(igraph)$cor[E(igraph)$cor < 0] <- -1

#为节点（微生物属）添加属性信息（界门纲目科属水平注释）
#“genus_taxonomy.txt” 记录了微生物的属性，读入该表后根据已知网络节点匹配对应的行
tax <- read.delim('taxonomy_genus_AT.txt', row.name = 1, check.names = FALSE, stringsAsFactors = FALSE)
tax
tax <- tax[as.character(V(igraph)$name), ]
tax


V(igraph)$phylum <- tax$phylum
V(igraph)$class <- tax$class
V(igraph)$order <- tax$order
V(igraph)$family <- tax$family
V(igraph)$genus <- tax$genus
#查看网络图
igraph
plot(igraph)


#边列表节点属性列表可以导入至 gephi 或 cytoscape 等网络可视化软件中进行编辑
#此外 igraph 也提供了可以被 gephi 或 cytoscape 等直接识别的格式
#graphml 格式，可使用 gephi 软件打开并进行可视化编辑
write.graph(igraph, 'genus_AT_2.graphml', format = 'graphml')

#gml 格式，可使用 cytoscape 软件打开并进行可视化编辑
write.graph(igraph, 'genus_AT_2.gml', format = 'gml')




####R语言计算节点和边特征####
##节点特征
#节点数量
length(V(igraph)$name)
#或
vcount(igraph)

#节点度（Degree）
#由于本示例是个无向网络，故无出度和入度之分
V(igraph)$degree <- degree(igraph)
V(igraph)$degree

#查看度分布
#可观察到微生物相关网络通常服从幂律分布，这个下节再讲怎样通过计算验证
degree_dist <- degree.distribution(igraph)[-1]
degree_num <- 1:max(V(igraph)$degree)

par(mfrow = c(1, 2))
hist(V(igraph)$degree, xlab = 'Degree', ylab = 'Frequency', 
     main = 'Degree distribution')
plot(degree_num, degree_dist, log = 'xy', xlab = 'Log-degree', 
     ylab = 'Log-intensity', main = 'Log-log degree distribution')

#查看节点度与其“邻居”的平均度的关系
#微生物网络中高度值的节点更倾向连接在一起，是普遍现象吗？
neighbor_degree <- graph.knn(igraph, V(igraph))$knn
plot(V(igraph)$degree, neighbor_degree, log = 'xy', 
     xlab = 'Log degree', ylab = 'Log average neighbor degree')

#加权度（Weighted degree）
V(igraph)$weight_degree <- strength(igraph)
V(igraph)$weight_degree

#接近中心性（Closeness centrality）
V(igraph)$closeness_centrality <- closeness(igraph)
V(igraph)$closeness_centrality

#介数中心性（Betweenness centrality）
V(igraph)$betweenness_centrality <- betweenness(igraph)
V(igraph)$betweenness_centrality

#特征向量中心性（Eigenvector centrality）
V(igraph)$eigenvector_centrality <- evcent(igraph)$vector
V(igraph)$eigenvector_centrality

#探索三种描述节点中心性的特征的关系
library(car)

scatter3d(V(igraph)$closeness_centrality, V(igraph)$betweenness_centrality, V(igraph)$eigenvector_centrality, 
          xlab =  'Closeness centrality', ylab = 'Betweenness centrality', zlab = 'Eigenvector centrality', 
          surface = FALSE)

#探索节点度和节点中心性的关系，如与特征向量中心性的关系
plot(V(igraph)$degree, V(igraph)$eigenvector_centrality, 
     xlab = 'Degree', ylab = 'Eigenvector centrality')


#模块划分，详情 ?cluster_fast_greedy，有多种模型
set.seed(123)
V(igraph)$modularity <- membership(cluster_fast_greedy(igraph))


#输出列表
node_list <- data.frame(
  node_id = V(igraph)$name, 
  degree = V(igraph)$degree, 
  weight_degree = V(igraph)$weight_degree, 
  closeness_centrality = V(igraph)$closeness_centrality, 
  betweenness_centrality = V(igraph)$betweenness_centrality, 
  eigenvector_centrality = V(igraph)$eigenvector_centrality,
  modularity = V(igraph)$modularity
)

head(node_list)
write.table(node_list, 'node_list.txt', sep = '\t', row.names = FALSE, quote = FALSE)



##边特征
#边的数量
ecount(igraph)

#权重（Weighted），已在数据读入时转化获得
E(igraph)$weight

#边介数中心性（Edge betweenness centrality）
E(igraph)$betweenness_centrality <- edge.betweenness(igraph)
E(igraph)$betweenness_centrality

#输出列表
edge <- data.frame(as_edgelist(igraph))    #igraph 的邻接列表转为边列表

edge_list <- data.frame(
  source = edge[[1]],
  target = edge[[2]],
  weight = E(igraph)$weight,
  correlation = E(igraph)$corr, 
  betweenness_centrality = E(igraph)$betweenness_centrality
)
head(edge_list)

write.table(edge_list, 'edge_list.txt', sep = '\t', row.names = FALSE, quote = FALSE)



##子图与普查
#所有尺寸的团的普查可以提供一个快照，将显示各尺寸的团的数量
census <- table(sapply(cliques(igraph), length))
census
plot(census)

#k 核
cores <- graph.coreness(igraph)
cores
sna::gplot.target(adjacency_unweight, cores, usearrows = FALSE, vertex.col = cores)

#二元组（dyad）和三元组（triad）
dyad.census(simplify(igraph))
triad.census(simplify(igraph))

#节点数量（number of nodes）和边数量（number of edges）
nodes_num <- length(V(igraph))
nodes_num

edges_num <- length(E(igraph))
edges_num

#平均度（average degree）
average_degree <- mean(degree(igraph))
#或者，2x边数量/节点数量
average_degree <- 2*edges_num/nodes_num



####模块内连通度（Zi）和模块间连通度（Pi）识别关键节点####
##计算模块内连通度（Zi）和模块间连通度（Pi）
source('zi_pi.r')

#上述的邻接矩阵类型的网络文件
adjacency_unweight <- read.delim('genus_AT_adjacency_unweight.txt', row.names = 1, sep = '\t', check.names = FALSE)

#节点属性列表，包含节点所划分的模块
nodes_list <- read.delim('nodes_list.txt', row.names = 1, sep = '\t', check.names = FALSE)

#两个文件的节点顺序要一致
nodes_list <- nodes_list[rownames(adjacency_unweight), ]

#计算模块内连通度（Zi）和模块间连通度（Pi）
#指定邻接矩阵、节点列表、节点列表中节点度和模块度的列名称
zi_pi <- zi.pi(nodes_list, adjacency_unweight, degree = 'degree', modularity_class = 'modularity')
head(zi_pi)

write.table(zi_pi, 'zi_pi_result.txt', sep = '\t', row.names = FALSE, quote = FALSE)


##可再根据阈值对节点划分为 4 种类型，并作图展示其分布
library(ggplot2)

zi_pi <- na.omit(zi_pi)   #NA 值最好去掉，不要当 0 处理
zi_pi[which(zi_pi$within_module_connectivities < 2.5 & zi_pi$among_module_connectivities < 0.62),'type'] <- 'Peripherals'
zi_pi[which(zi_pi$within_module_connectivities < 2.5 & zi_pi$among_module_connectivities > 0.62),'type'] <- 'Connectors'
zi_pi[which(zi_pi$within_module_connectivities > 2.5 & zi_pi$among_module_connectivities < 0.62),'type'] <- 'Module hubs'
zi_pi[which(zi_pi$within_module_connectivities > 2.5 & zi_pi$among_module_connectivities > 0.62),'type'] <- 'Network hubs'

ggplot(zi_pi, aes(among_module_connectivities, within_module_connectivities)) +
  geom_point(aes(color = type), alpha = 0.5, size = 2) +
  scale_color_manual(values = c('gray','red','blue','purple'), 
                     limits = c('Peripherals', 'Connectors', 'Module hubs', 'Network hubs'))+
  theme(panel.grid = element_blank(), axis.line = element_line(colour = 'black'), 
        panel.background = element_blank(), legend.key = element_blank()) +
  labs(x = 'Among-module connectivities', y = 'Within-module connectivities', color = '') +
  geom_vline(xintercept = 0.62) +
  geom_hline(yintercept = 2.5)
















#1%高丰度属水平top100,特别值得注意的是要将OTU表格注释到属水平即可
#############基于丰度相关性的微生物共发生网络
##计算微生物丰度间的相关系数
setwd("E:/博士数据2/小麦微生物组数据/customer_backup/otus/networks/")
library(Hmisc)

#以1%高丰度OTU水平丰度为例，“otu_AAT.txt” 是一个1%高丰度OTU水平的微生物丰度表
genus_AT <- read.delim('genus_AT_order.txt', row.name = 1, check.names = FALSE)
genus_AT

#可选事先过滤一些低丰度或低频的类群
#otu_AAT <- otu_AAT[which(rowSums(otu_AAT) >= 0.005), ]    #例如只保留相对丰度总和高于 0.005 的属

#otu_AAT1 <- otu_AAT
#otu_AAT1[otu_AAT1>0] <- 1
#otu_AAT1 <- otu_AAT1[which(rowSums(otu_AAT1) >= 5), ]    #例如只保留在 5 个及以上样本中出现的属

#计算两属之间是否存在丰度变化的相关性，以 spearman 相关系数为例
genus_AT_corr <- rcorr(t(genus_AT), type = 'spearman')

#阈值筛选
#将 spearman 相关系数低于 0.7 的关系剔除，即 r>=0.7
r <- genus_AT_corr$r
r[abs(r) < 0.7] <- 0

#选取显著性 p 值小于 0.05 的相关系数，即 p<0.05
p <- genus_AT_corr$P
p <- p.adjust(p, method = 'BH')    #可选 p 值校正，这里使用 BH 法校正 p 值
p[p>=0.05] <- -1
p[p<0.05 & p>=0] <- 1
p[p==-1] <- 0

#根据上述筛选的 r 值和 p 值保留数据
z <- r * p
diag(z) <- 0    #将相关矩阵中对角线中的值（代表了自相关）转为 0
head(z)[1:6,1:6]

#如此便得到了邻接矩阵格式的网络文件（微生物属的相关系数矩阵）
write.table(data.frame(z, check.names = FALSE), 'genus_AT_corr.matrix.txt', col.names = NA, sep = '\t', quote = FALSE)


##将微生物互作网络转换为0、1矩阵，数值“1”表示微生物 OTU 之间存在互作，“0”表示无互作
genus_AT_corr.matrix2 <- read.delim('genus_AT_corr.matrix.txt', row.name = 1, check.names = FALSE)
head(genus_AT_corr.matrix2)[1:6,1:6]
genus_AT_corr.matrix2[genus_AT_corr.matrix2 > 0] <- 1
genus_AT_corr.matrix2[genus_AT_corr.matrix2 < 0] <- 1
head(genus_AT_corr.matrix2)[1:6,1:6]
write.table(data.frame(genus_AT_corr.matrix2, check.names = FALSE), 'genus_AT_adjacency_unweight.txt', col.names = NA, sep = '\t', quote = FALSE)


##获得网络
library(igraph)


#将邻接矩阵转化为 igraph 网络的邻接列表
#构建含权的无向网络，权重代表了微生物属间丰度的 spearman 相关系数 
igraph <- graph.adjacency(z, weighted = TRUE, mode = 'undirected')
igraph
#自相关也可以通过该式去除
igraph <- simplify(igraph)

#孤立节点的删除（删除度为 0 的节点）
igraph <- delete.vertices(igraph, names(degree(igraph)[degree(igraph) == 0]))
igraph
#该模式下，边权重代表了相关系数
#由于权重通常为正值，因此最好取个绝对值，相关系数重新复制一列
E(igraph)$correlation <- E(igraph)$weight
E(igraph)$weight <- abs(E(igraph)$weight)
E(igraph)$cor <- E(igraph)$correlation 
E(igraph)$cor[E(igraph)$cor > 0] <- 1
E(igraph)$cor[E(igraph)$cor < 0] <- -1

#为节点（微生物属）添加属性信息（界门纲目科属水平注释）
#“genus_taxonomy.txt” 记录了微生物的属性，读入该表后根据已知网络节点匹配对应的行
tax <- read.delim('taxonomy_genus_AT.txt', row.name = 1, check.names = FALSE, stringsAsFactors = FALSE)
tax
tax <- tax[as.character(V(igraph)$name), ]
tax


V(igraph)$phylum <- tax$phylum
V(igraph)$class <- tax$class
V(igraph)$order <- tax$order
V(igraph)$family <- tax$family
V(igraph)$genus <- tax$genus
#查看网络图
igraph
plot(igraph)


#边列表节点属性列表可以导入至 gephi 或 cytoscape 等网络可视化软件中进行编辑
#此外 igraph 也提供了可以被 gephi 或 cytoscape 等直接识别的格式
#graphml 格式，可使用 gephi 软件打开并进行可视化编辑
write.graph(igraph, 'genus_AT_2.graphml', format = 'graphml')

#gml 格式，可使用 cytoscape 软件打开并进行可视化编辑
write.graph(igraph, 'genus_AT_2.gml', format = 'gml')




####对gephi结果进行模块的手动的划分

rm(list=ls())
setwd("E:/博士数据2/小麦微生物组数据/customer_backup/otus/networks2/R_Input1/")
output <- "E:/博士数据2/小麦微生物组数据/customer_backup/otus/networks2/R_Output1/"

##由于在R语言中未计算出节点度的特征，故需要在gephi中计算好后进行数据表格的输入
network_feature <- read.table("all_network_gephi2.txt", row.names=1,sep="\t", header=T, blank.lines.skip=F, check.names=F)
network_feature
##如何将数据框的行名进行修改、重新命名
rownames(network_feature) <- network_feature$v_name
network_feature

##查看模块的划分
unique(network_feature$modularity_class)

module0_gephi <- network_feature[network_feature$modularity_class=="0",]
module0_gephi 
write.table(module0_gephi,paste0(output,"All_module0_gephi.txt"),sep="\t",quote=F)


module1_gephi <- network_feature[network_feature$modularity_class=="1",]
module1_gephi 
write.table(module1_gephi,paste0(output,"All_module1_gephi.txt"),sep="\t",quote=F)


module2_gephi <- network_feature[network_feature$modularity_class=="2",]
module2_gephi 
write.table(module2_gephi,paste0(output,"All_module2_gephi.txt"),sep="\t",quote=F)

module3_gephi <- network_feature[network_feature$modularity_class=="3",]
module3_gephi 
write.table(module3_gephi,paste0(output,"All_module3_gephi.txt"),sep="\t",quote=F)


module4_gephi <- network_feature[network_feature$modularity_class=="4",]
module4_gephi 
write.table(module4_gephi,paste0(output,"All_module4_gephi.txt"),sep="\t",quote=F)


module5_gephi <- network_feature[network_feature$modularity_class=="5",]
module5_gephi 
write.table(module5_gephi,paste0(output,"All_module5_gephi.txt"),sep="\t",quote=F)


module6_gephi <- network_feature[network_feature$modularity_class=="6",]
module6_gephi 
write.table(module6_gephi,paste0(output,"All_module6_gephi.txt"),sep="\t",quote=F)


module7_gephi <- network_feature[network_feature$modularity_class=="7",]
module7_gephi 
write.table(module7_gephi,paste0(output,"All_module7_gephi.txt"),sep="\t",quote=F)

module8_gephi <- network_feature[network_feature$modularity_class=="8",]
module8_gephi 
write.table(module8_gephi,paste0(output,"All_module8_gephi.txt"),sep="\t",quote=F)


module9_gephi <- network_feature[network_feature$modularity_class=="9",]
module9_gephi 
write.table(module9_gephi,paste0(output,"All_module9_gephi.txt"),sep="\t",quote=F)

module10_gephi <- network_feature[network_feature$modularity_class=="10",]
module10_gephi 
write.table(module10_gephi,paste0(output,"All_module10_gephi.txt"),sep="\t",quote=F)


module11_gephi <- network_feature[network_feature$modularity_class=="11",]
module11_gephi 
write.table(module11_gephi,paste0(output,"All_module11_gephi.txt"),sep="\t",quote=F)

#不同模块核心微生物出图
pdf(paste0(output,"all_module0_CD.pdf"), width=7, height=7/6)
par(mfrow=c(1,7), mar=c(2,3,1,0) )
#par(mfrow=c(1,1), mar=c(0,0,0,0))

plot(module0_gephi$closnesscentrality, module0_gephi$degree, 
     col= module0_gephi$v_cols,
     pch = 20,cex =0.5,
     xlab = 'Closeness Centrality', ylab = 'Degree',
     #ylim = c (0,50),
     main="Module 0")

dev.off()


pdf(paste0(output,"all_module1_CD.pdf"), width=7, height=7/6)
par(mfrow=c(1,7), mar=c(2,3,1,0) )
#par(mfrow=c(1,1), mar=c(0,0,0,0))

plot(module1_gephi$closnesscentrality, module1_gephi$degree, 
     col= module1_gephi$v_cols,
     pch = 20,cex =0.5,
     xlab = 'Closeness Centrality', ylab = 'Degree',
     #ylim = c (0,50),
     main="Module 1")

dev.off()


pdf(paste0(output,"all_module2_CD.pdf"), width=7, height=7/6)
par(mfrow=c(1,7), mar=c(2,3,1,0) )

plot(module2_gephi$closnesscentrality, module2_gephi$degree, 
     col= module2_gephi$v_cols,
     pch = 20,cex =0.5,
     xlab = '', ylab = '',
     #ylim = c (0,120),
     main="Module 2")
dev.off()

pdf(paste0(output,"all_module3_CD.pdf"), width=7, height=7/6)
par(mfrow=c(1,7), mar=c(2,3,1,0) )
plot(module3_gephi$closnesscentrality, module3_gephi$degree, 
     col= module3_gephi$v_cols,
     pch = 20,cex =0.5,
     xlab = '', ylab = '',
     #ylim = c (0,60),
     main="Module 3")
dev.off()

pdf(paste0(output,"all_module4_CD.pdf"), width=7, height=7/6)
par(mfrow=c(1,7), mar=c(2,3,1,0) )
plot(module4_gephi$closnesscentrality, module4_gephi$degree, 
     col= module4_gephi$v_cols,
     pch = 20,cex =0.5,
     xlab = '', ylab = '',
     #ylim = c (0,60),
     main="Module 4")
dev.off()


pdf(paste0(output,"all_module5_CD.pdf"), width=7, height=7/6)
par(mfrow=c(1,7), mar=c(2,3,1,0) )
plot(module5_gephi$closnesscentrality, module5_gephi$degree, 
     col= module5_gephi$v_cols,
     pch = 20,cex =0.5,
     xlab = '', ylab = '',
     #ylim = c (0,60),
     main="Module 5")
dev.off()

pdf(paste0(output,"all_module6_CD.pdf"), width=7, height=7/6)
par(mfrow=c(1,7), mar=c(2,3,1,0) )
plot(module6_gephi$closnesscentrality, module6_gephi$degree, 
     col= module6_gephi$v_cols,
     pch = 20,cex =0.5,
     xlab = '', ylab = '',
     #ylim = c (0,60),
     main="Module 6")
dev.off()

pdf(paste0(output,"all_module7_CD.pdf"), width=7, height=7/6)
par(mfrow=c(1,7), mar=c(2,3,1,0) )
plot(module7_gephi$closnesscentrality, module7_gephi$degree, 
     col= module7_gephi$v_cols,
     pch = 20,cex =0.5,
     xlab = '', ylab = '',
     #ylim = c (0,60),
     main="Module 7")
dev.off()


pdf(paste0(output,"all_module8_CD.pdf"), width=7, height=7/6)
par(mfrow=c(1,7), mar=c(2,3,1,0) )
plot(module8_gephi$closnesscentrality, module8_gephi$degree, 
     col= module8_gephi$v_cols,
     pch = 20,cex =0.5,
     xlab = '', ylab = '',
     #ylim = c (0,60),
     main="Module 8")
dev.off()

pdf(paste0(output,"all_module9_CD.pdf"), width=7, height=7/6)
par(mfrow=c(1,7), mar=c(2,3,1,0) )
plot(module9_gephi$closnesscentrality, module9_gephi$degree, 
     col= module9_gephi$v_cols,
     pch = 20,cex =0.5,
     xlab = '', ylab = '',
     #ylim = c (0,60),
     main="Module 9")
dev.off()

pdf(paste0(output,"all_module10_CD.pdf"), width=7, height=7/6)
par(mfrow=c(1,7), mar=c(2,3,1,0) )
plot(module10_gephi$closnesscentrality, module10_gephi$degree, 
     col= module10_gephi$v_cols,
     pch = 20,cex =0.5,
     xlab = '', ylab = '',
     #ylim = c (0,60),
     main="Module 10")
dev.off()


pdf(paste0(output,"all_module11_CD.pdf"), width=7, height=7/6)
par(mfrow=c(1,7), mar=c(2,3,1,0) )
plot(module11_gephi$closnesscentrality, module11_gephi$degree, 
     col= module11_gephi$v_cols,
     pch = 20,cex =0.5,
     xlab = '', ylab = '',
     #ylim = c (0,60),
     main="Module 11")
dev.off()


#总核心微生物出图
pdf(paste0(output,"all_module_all_CD.pdf"), width=7, height=7/6)
par(mfrow=c(1,7), mar=c(2,3,1,0) )
plot(network_feature$closnesscentrality, network_feature$degree, 
     col= network_feature$v_cols,
     pch = 20,cex =0.5,
     xlab = '', ylab = '',
     ylim = c (0,120),
     main="All")

dev.off()



###每个模块中各个处理间OTU的丰度差异，首先得到每个模块中OTU所对应的物种丰度，再做堆叠柱状图
### plotting average module response

## ROOT module 0
module0_rownames <- rownames(module0_gephi)
module0_rownames
module0_abudance <- otu_norm_16s[module0_rownames,]
module0_abudance
write.table(module0_abudance,paste0(output,"all_module0_abudance.txt"),sep="\t",quote=F)

pdf(paste0(output,"negative_module0_abudance.pdf"),width=7,height=7/6)
par(mfrow=c(1,7), mar=c(0.5,3.5,2,0))

CS_cols <- c("dodgerblue4","orange","green","violetred4")
names(CS_cols) <- c("Ggt","S_135", "SA","H")

bargraph.CI(design_filter_16s$Group, colSums(module0_abudance)/1000, 
            las=2, ylab="cumulative relative abundance", cex.lab=.5, cex.axis=.7, cex.names=.7,
            err.width=.025, main="Module 0", col=CS_cols, border=F)
dev.off()



## ROOT module 1
module1_rownames <- rownames(module1_gephi)
module1_rownames
module1_abudance <- otu_norm_16s[module1_rownames,]
module1_abudance
write.table(module1_abudance,paste0(output,"all_module1_abudance.txt"),sep="\t",quote=F)

pdf(paste0(output,"all_module1_abudance.pdf"),width=7,height=7/6)
par(mfrow=c(1,7), mar=c(0.5,3.5,2,0))

CS_cols <- c("dodgerblue4","orange","green","violetred4")
names(CS_cols) <- c("Ggt","S_135", "SA","H")

bargraph.CI(design_filter_16s$Group, colSums(module1_abudance)/1000, 
            las=2, ylab="cumulative relative abundance", cex.lab=.5, cex.axis=.7, cex.names=.7,
            err.width=.025, main="Module 1", col=CS_cols, border=F)
dev.off()




## ROOT module 2
module2_rownames <- rownames(module2_gephi)
module2_rownames
#简单粗暴求两个数据集的交集，首先将所含内容较小的数据集的名称单独赋值提取出来，然后用语法，大数据集[小数据集名称]即可
module2_abudance <- otu_norm_16s[module2_rownames,]
module2_abudance
write.table(module2_abudance,paste0(output,"all_module2_abudance.txt"),sep="\t",quote=F)

pdf(paste0(output,"all_module2_abudance.pdf"),width=7,height=7/6)
par(mfrow=c(1,7), mar=c(0.5,3.5,2,0))

CS_cols <- c("dodgerblue4","orange","green","violetred4")
names(CS_cols) <- c("Ggt","S_135", "SA","H")


bargraph.CI(design_filter_16s$Group, colSums(module2_abudance)/1000, 
            las=2, ylab="", cex.lab=.5, cex.axis=.7, cex.names=.7,
            err.width=.025, main="Module 2", col=CS_cols, border=F)

dev.off()


## ROOT module 3
module3_rownames <- rownames(module3_gephi)
module3_rownames
#简单粗暴求两个数据集的交集，首先将所含内容较小的数据集的名称单独赋值提取出来，然后用语法，大数据集[小数据集名称]即可
module3_abudance <- otu_norm_16s[module3_rownames,]
module3_abudance
write.table(module3_abudance,paste0(output,"all_module3_abudance.txt"),sep="\t",quote=F)


pdf(paste0(output,"all_module3_abudance.pdf"),width=7,height=7/6)
par(mfrow=c(1,7), mar=c(0.5,3.5,2,0))

CS_cols <- c("dodgerblue4","orange","green","violetred4")
names(CS_cols) <- c("Ggt","S_135", "SA","H")


bargraph.CI(design_filter_16s$Group, colSums(module3_abudance)/1000, 
            las=2, ylab="", cex.lab=.5, cex.axis=.7, cex.names=.7,
            err.width=.025, main="Module 3", col=CS_cols, border=F)

dev.off()


## ROOT module 4
module4_rownames <- rownames(module4_gephi)
module4_rownames
#简单粗暴求两个数据集的交集，首先将所含内容较小的数据集的名称单独赋值提取出来，然后用语法，大数据集[小数据集名称]即可
module4_abudance <- otu_norm_16s[module4_rownames,]
module4_abudance
write.table(module4_abudance,paste0(output,"all_module4_abudance.txt"),sep="\t",quote=F)


pdf(paste0(output,"all_module4_abudance.pdf"),width=7,height=7/6)
par(mfrow=c(1,7), mar=c(0.5,3.5,2,0))

CS_cols <- c("dodgerblue4","orange","green","violetred4")
names(CS_cols) <- c("Ggt","S_135", "SA","H")


bargraph.CI(design_filter_16s$Group, colSums(module4_abudance)/1000, 
            las=2, ylab="", cex.lab=.5, cex.axis=.7, cex.names=.7,
            err.width=.025, main="Module 4", col=CS_cols, border=F)

dev.off()




## ROOT module 5
module5_rownames <- rownames(module5_gephi)
module5_rownames
#简单粗暴求两个数据集的交集，首先将所含内容较小的数据集的名称单独赋值提取出来，然后用语法，大数据集[小数据集名称]即可
module5_abudance <- otu_norm_16s[module5_rownames,]
module5_abudance
write.table(module5_abudance,paste0(output,"all_module5_abudance.txt"),sep="\t",quote=F)


pdf(paste0(output,"all_module5_abudance.pdf"),width=7,height=7/6)
par(mfrow=c(1,7), mar=c(0.5,3.5,2,0))

CS_cols <- c("dodgerblue4","orange","green","violetred4")
names(CS_cols) <- c("Ggt","S_135", "SA","H")


bargraph.CI(design_filter_16s$Group, colSums(module5_abudance)/1000, 
            las=2, ylab="", cex.lab=.5, cex.axis=.7, cex.names=.7,
            err.width=.025, main="Module 5", col=CS_cols, border=F)

dev.off()




## ROOT module 6
module6_rownames <- rownames(module6_gephi)
module6_rownames
#简单粗暴求两个数据集的交集，首先将所含内容较小的数据集的名称单独赋值提取出来，然后用语法，大数据集[小数据集名称]即可
module6_abudance <- otu_norm_16s[module6_rownames,]
module6_abudance
write.table(module6_abudance,paste0(output,"all_module6_abudance.txt"),sep="\t",quote=F)


pdf(paste0(output,"all_module6_abudance.pdf"),width=7,height=7/6)
par(mfrow=c(1,7), mar=c(0.5,3.5,2,0))

CS_cols <- c("dodgerblue4","orange","green","violetred4")
names(CS_cols) <- c("Ggt","S_135", "SA","H")


bargraph.CI(design_filter_16s$Group, colSums(module6_abudance)/1000, 
            las=2, ylab="", cex.lab=.5, cex.axis=.7, cex.names=.7,
            err.width=.025, main="Module 6", col=CS_cols, border=F)

dev.off()




## ROOT module 7
module7_rownames <- rownames(module7_gephi)
module7_rownames
#简单粗暴求两个数据集的交集，首先将所含内容较小的数据集的名称单独赋值提取出来，然后用语法，大数据集[小数据集名称]即可
module7_abudance <- otu_norm_16s[module7_rownames,]
module7_abudance
write.table(module7_abudance,paste0(output,"all_module7_abudance.txt"),sep="\t",quote=F)


pdf(paste0(output,"all_module7_abudance.pdf"),width=7,height=7/6)
par(mfrow=c(1,7), mar=c(0.5,3.5,2,0))

CS_cols <- c("dodgerblue4","orange","green","violetred4")
names(CS_cols) <- c("Ggt","S_135", "SA","H")


bargraph.CI(design_filter_16s$Group, colSums(module7_abudance)/1000, 
            las=2, ylab="", cex.lab=.5, cex.axis=.7, cex.names=.7,
            err.width=.025, main="Module 7", col=CS_cols, border=F)

dev.off()




## ROOT module 8
module8_rownames <- rownames(module8_gephi)
module8_rownames
#简单粗暴求两个数据集的交集，首先将所含内容较小的数据集的名称单独赋值提取出来，然后用语法，大数据集[小数据集名称]即可
module8_abudance <- otu_norm_16s[module8_rownames,]
module8_abudance
write.table(module8_abudance,paste0(output,"all_module8_abudance.txt"),sep="\t",quote=F)


pdf(paste0(output,"all_module8_abudance.pdf"),width=7,height=7/6)
par(mfrow=c(1,7), mar=c(0.5,3.5,2,0))

CS_cols <- c("dodgerblue4","orange","green","violetred4")
names(CS_cols) <- c("Ggt","S_135", "SA","H")


bargraph.CI(design_filter_16s$Group, colSums(module8_abudance)/1000, 
            las=2, ylab="", cex.lab=.5, cex.axis=.7, cex.names=.7,
            err.width=.025, main="Module 8", col=CS_cols, border=F)

dev.off()




## ROOT module 9
module9_rownames <- rownames(module9_gephi)
module9_rownames
#简单粗暴求两个数据集的交集，首先将所含内容较小的数据集的名称单独赋值提取出来，然后用语法，大数据集[小数据集名称]即可
module9_abudance <- otu_norm_16s[module9_rownames,]
module9_abudance
write.table(module9_abudance,paste0(output,"all_module9_abudance.txt"),sep="\t",quote=F)


pdf(paste0(output,"all_module9_abudance.pdf"),width=7,height=7/6)
par(mfrow=c(1,7), mar=c(0.5,3.5,2,0))

CS_cols <- c("dodgerblue4","orange","green","violetred4")
names(CS_cols) <- c("Ggt","S_135", "SA","H")


bargraph.CI(design_filter_16s$Group, colSums(module9_abudance)/1000, 
            las=2, ylab="", cex.lab=.5, cex.axis=.7, cex.names=.7,
            err.width=.025, main="Module 9", col=CS_cols, border=F)

dev.off()




## ROOT module 10
module10_rownames <- rownames(module10_gephi)
module10_rownames
#简单粗暴求两个数据集的交集，首先将所含内容较小的数据集的名称单独赋值提取出来，然后用语法，大数据集[小数据集名称]即可
module10_abudance <- otu_norm_16s[module10_rownames,]
module10_abudance
write.table(module10_abudance,paste0(output,"all_module10_abudance.txt"),sep="\t",quote=F)


pdf(paste0(output,"all_module10_abudance.pdf"),width=7,height=7/6)
par(mfrow=c(1,7), mar=c(0.5,3.5,2,0))

CS_cols <- c("dodgerblue4","orange","green","violetred4")
names(CS_cols) <- c("Ggt","S_135", "SA","H")


bargraph.CI(design_filter_16s$Group, colSums(module10_abudance)/1000, 
            las=2, ylab="", cex.lab=.5, cex.axis=.7, cex.names=.7,
            err.width=.025, main="Module 10", col=CS_cols, border=F)

dev.off()




## ROOT module 11
module11_rownames <- rownames(module11_gephi)
module11_rownames
#简单粗暴求两个数据集的交集，首先将所含内容较小的数据集的名称单独赋值提取出来，然后用语法，大数据集[小数据集名称]即可
module11_abudance <- otu_norm_16s[module11_rownames,]
module11_abudance
write.table(module11_abudance,paste0(output,"all_module11_abudance.txt"),sep="\t",quote=F)


pdf(paste0(output,"all_module11_abudance.pdf"),width=7,height=7/6)
par(mfrow=c(1,7), mar=c(0.5,3.5,2,0))

CS_cols <- c("dodgerblue4","orange","green","violetred4")
names(CS_cols) <- c("Ggt","S_135", "SA","H")


bargraph.CI(design_filter_16s$Group, colSums(module11_abudance)/1000, 
            las=2, ylab="", cex.lab=.5, cex.axis=.7, cex.names=.7,
            err.width=.025, main="Module 11", col=CS_cols, border=F)

dev.off()


## All

all_rownames <- rownames(network_feature)
all_rownames
#简单粗暴求两个数据集的交集，首先将所含内容较小的数据集的名称单独赋值提取出来，然后用语法，大数据集[小数据集名称]即可
all_abudance <- otu_norm_16s[all_rownames,]
all_abudance
write.table(all_abudance,paste0(output,"all_abudance.txt"),sep="\t",quote=F)


pdf(paste0(output,"all_abudance.pdf"),width=7,height=7/6)
par(mfrow=c(1,7), mar=c(0.5,3.5,2,0))

CS_cols <- c("dodgerblue4","orange","green","violetred4")
names(CS_cols) <- c("Ggt","S_135", "SA","H")


bargraph.CI(design_filter_16s$Group, colSums(all_abudance)/1000, 
            las=2, ylab="", cex.lab=.5, cex.axis=.7, cex.names=.7,
            err.width=.025, main="All", col=CS_cols, border=F)
dev.off()



## Legend
plot.new()
par(mar=c(0.5,0,2,0))
legend("left",  bty="n", cex=1, #x.intersp=0.1, y.intersp=1,
       legend=names(CS_cols), 
       fill=CS_cols, 
       border=CS_cols , xpd=T)

dev.off()




###不同模块OTU的物种分类单位
#### taxonomies of module OTUs 

## ROOT module 0
## defining bacteria of root modules

module0_rownames <- rownames(module0_gephi)
module0_rownames

### counts of bacteria OTUs
module0_taxonomy <- as.data.frame(table(tax_filter_16s[module0_rownames, "labels"] ) )
module0_taxonomy
colnames(module0_taxonomy) <- c("phylum", "Module 0")
module0_taxonomy



## ROOT module 1
## defining bacteria of root modules

module1_rownames <- rownames(module1_gephi)
module1_rownames

### counts of bacteria OTUs
module1_taxonomy <- as.data.frame(table(tax_filter_16s[module1_rownames, "labels"] ) )
module1_taxonomy
colnames(module1_taxonomy) <- c("phylum", "Module 1")
module1_taxonomy



## ROOT module 2
## defining bacteria of root modules

module2_rownames <- rownames(module2_gephi)
module2_rownames

### counts of bacteria OTUs
module2_taxonomy <- as.data.frame(table(tax_filter_16s[module2_rownames, "labels"] ) )
module2_taxonomy
colnames(module2_taxonomy) <- c("phylum", "Module 2")
module2_taxonomy



## ROOT module 3
## defining bacteria of root modules

module3_rownames <- rownames(module3_gephi)
module3_rownames

### counts of bacteria OTUs
module3_taxonomy <- as.data.frame(table(tax_filter_16s[module3_rownames, "labels"] ) )
module3_taxonomy
colnames(module3_taxonomy) <- c("phylum", "Module 3")
module3_taxonomy



## ROOT module 4
## defining bacteria of root modules

module4_rownames <- rownames(module4_gephi)
module4_rownames

### counts of bacteria OTUs
module4_taxonomy <- as.data.frame(table(tax_filter_16s[module4_rownames, "labels"] ) )
module4_taxonomy
colnames(module4_taxonomy) <- c("phylum", "Module 4")
module4_taxonomy



## ROOT module 5
## defining bacteria of root modules

module5_rownames <- rownames(module5_gephi)
module5_rownames

### counts of bacteria OTUs
module5_taxonomy <- as.data.frame(table(tax_filter_16s[module5_rownames, "labels"] ) )
module5_taxonomy
colnames(module5_taxonomy) <- c("phylum", "Module 5")
module5_taxonomy



## ROOT module 6
## defining bacteria of root modules

module6_rownames <- rownames(module6_gephi)
module6_rownames

### counts of bacteria OTUs
module6_taxonomy <- as.data.frame(table(tax_filter_16s[module6_rownames, "labels"] ) )
module6_taxonomy
colnames(module6_taxonomy) <- c("phylum", "Module 6")
module6_taxonomy



## ROOT module 7
## defining bacteria of root modules

module7_rownames <- rownames(module7_gephi)
module7_rownames

### counts of bacteria OTUs
module7_taxonomy <- as.data.frame(table(tax_filter_16s[module7_rownames, "labels"] ) )
module7_taxonomy
colnames(module7_taxonomy) <- c("phylum", "Module 7")
module7_taxonomy



## ROOT module 8
## defining bacteria of root modules

module8_rownames <- rownames(module8_gephi)
module8_rownames

### counts of bacteria OTUs
module8_taxonomy <- as.data.frame(table(tax_filter_16s[module8_rownames, "labels"] ) )
module8_taxonomy
colnames(module8_taxonomy) <- c("phylum", "Module 8")
module8_taxonomy



## ROOT module 9
## defining bacteria of root modules

module9_rownames <- rownames(module9_gephi)
module9_rownames

### counts of bacteria OTUs
module9_taxonomy <- as.data.frame(table(tax_filter_16s[module9_rownames, "labels"] ) )
module9_taxonomy
colnames(module9_taxonomy) <- c("phylum", "Module 9")
module9_taxonomy



## ROOT module 10
## defining bacteria of root modules

module10_rownames <- rownames(module10_gephi)
module10_rownames

### counts of bacteria OTUs
module10_taxonomy <- as.data.frame(table(tax_filter_16s[module10_rownames, "labels"] ) )
module10_taxonomy
colnames(module10_taxonomy) <- c("phylum", "Module 10")
module10_taxonomy



## ROOT module 11
## defining bacteria of root modules

module11_rownames <- rownames(module11_gephi)
module11_rownames

### counts of bacteria OTUs
module11_taxonomy <- as.data.frame(table(tax_filter_16s[module11_rownames, "labels"] ) )
module11_taxonomy
colnames(module11_taxonomy) <- c("phylum", "Module 11")
module11_taxonomy


module <- merge(module0_taxonomy, module1_taxonomy, all=T, by="phylum") 
module


module <- merge(module, module2_taxonomy, all=T, by="phylum") 
module

module <- merge(module, module3_taxonomy, all=T, by="phylum") 
module

module <- merge(module, module4_taxonomy, all=T, by="phylum") 
module

module <- merge(module, module5_taxonomy, all=T, by="phylum") 
module

module <- merge(module, module6_taxonomy, all=T, by="phylum") 
module

module <- merge(module, module7_taxonomy, all=T, by="phylum") 
module

module <- merge(module, module8_taxonomy, all=T, by="phylum") 
module

module <- merge(module, module9_taxonomy, all=T, by="phylum") 
module


module <- merge(module, module10_taxonomy, all=T, by="phylum") 
module

module <- merge(module, module11_taxonomy, all=T, by="phylum") 
module


all_taxonomy <- as.data.frame(table(tax_filter_16s[, "labels"] ) )
colnames(all_taxonomy) <- c("phylum", "all")

module <- merge(module, all_taxonomy, all=T, by="phylum") 
module

bacteria_modules <- module
bacteria_modules

bacteria_modules_mat <- bacteria_modules[2:14]
rownames(bacteria_modules_mat) <- bacteria_modules$phylum
bacteria_modules_mat[is.na(bacteria_modules_mat)] <- 0
colSums(bacteria_modules_mat)

bacteria_modules_prop <- t(t(bacteria_modules_mat)/colSums(bacteria_modules_mat) ) * 1
bacteria_modules_prop
colSums(bacteria_modules_prop)

##图4c
pdf(paste0(output,"all_module_taxonomy.pdf"), width=7, height=7/6)
par(mfrow=c(1,4), mar=c(2,4,1,0))

### bacteria
# PHYLA_label_cols_16s_legend
table(rownames(bacteria_modules_prop) %in% PHYLA_label_cols_16s_legend$labels) 
bp <- barplot(bacteria_modules_prop[,1:13],
              las=2, border=NA, axes=F, cex.names=.5,
              col=PHYLA_label_cols_16s_legend[rownames(bacteria_modules_prop),]$cols )
text(bp, 1.1, labels=c(colSums(bacteria_modules_mat)[1:13], xpd=T, cex=.6, las=2))

## legend from Fig. S3
plot.new()
legend("left", bty="n", cex=0.5, #x.intersp=0.1, #y.intersp=0.75,
       legend=rev(PHYLA_label_cols_16s_legend$labels), 
       fill=rev(PHYLA_label_cols_16s_legend$cols), 
       border=rev(PHYLA_label_cols_16s_legend$cols) )
dev.off()

#cbind(bacteria_modules_prop[,1:16], bacteria_modules_prop[,4]),



####不同模块核心微生物出图
##barplot绘制堆叠柱状图
barplot(VADeaths, 
        col = NULL, #col: 设置条形底纹或者填充颜色。
        border =par("fg"),#border：设置条形边框颜色。如果设置为NA，则消除了边缘
        main = NULL, sub = NULL, 
        xlab = NULL,ylab = NULL, #xlab和ylab：设置x轴与y轴的lable
        xlim = NULL, ylim = NULL,  #xlim和ylim:设置图形x轴与y轴的范围。
        
        beside = FALSE,#beside:逻辑参数。如果FALSE，那么将绘画堆叠式的条形；如果是TRUE，将绘画并列式条形。
        horiz = FALSE, #horiz：逻辑参数。设置图形是水平或是垂直
        
        width = 1, #width：设置条形的宽度
        space = NULL, #space：设置各个条形间的间隔宽度。相当于各个条形宽度的一部分。
        names.arg = NULL,  #names.arg:设置条形标签（barlabels）。
        
        
        #density 和 angle : 设置柱子用线条填充，density 控制线条的密度， angel 控制线条的角度
        density = NULL,  #density:底纹的密度。默认值为NULL
        angle =45,  #angle：设置底纹的斜率
        
        
        axes = TRUE, #axes:逻辑参数。设置图形是否显示坐标轴。
        las=1,            #设置刻度值的方向,  0表示总是平行于坐标轴；1表示总是水平方向；2表示总是垂直于坐标轴；3表示总是垂直方向。
        
        yaxt= "s", #是否绘制Y坐标轴，s 绘制，n不绘制
        
        axisnames = TRUE, #axisnames：逻辑参数。设置是否显示x轴条形标签
        cex.axis=par("cex.axis"),#cex.axis:设置坐标轴数值的大小。
        cex.names=par("cex.axis"), #cex.names: 设置条形标签（barlabels）的大小。
        
        add = FALSE #add = “TRUE”将barplot加在目前已经有的图上
        
)


####module0
#读取数据
module0 <-read.table("module0_top20_hubtaxa.txt",sep="\t",header = TRUE,check.names=FALSE,comment.char="@")
module0
#数据处理一下，以适合用barplot绘制柱状图
module0 <-module0[order(module0[,2],decreasing =T),]
module01 <-module0[,2:ncol(module0)]
rownames(module01)<-module0[,1]
module01<-as.matrix(module01)
module01
colnames(module01)


##module0核心微生物出图
pdf(paste0(output,"module0_top20_hubtaxa.pdf"), width=7, height=7/6)
par(mfrow=c(1,4), mar=c(2,4,1,0))
mycol<- c("green","red") 

p <- barplot(module01,
             las=2, border=NA, axes=T, cex.names=.5,ylim=c(0,130),main="Module0",
             col=mycol[1:nrow(module01)],ylab="Degree ")

## legend 
plot.new()
legend("top", bty="n", cex=0.5, #x.intersp=0.1, #y.intersp=0.75,
       legend=rev(rownames(module01)), 
       fill=rev(mycol), 
       border=rev(mycol) )
dev.off()



####module1
#读取数据
module1 <-read.table("module1_top20_hubtaxa.txt",sep="\t",header = TRUE,check.names=FALSE,comment.char="@")
module1
#数据处理一下，以适合用barplot绘制柱状图
module1 <-module1[order(module1[,2],decreasing =T),]
module1
module11 <-module1[,2:ncol(module1)]
rownames(module11)<-module1[,1]
module11<-as.matrix(module11)
module11
colnames(module11)


##module1核心微生物出图
pdf(paste0(output,"module1_top20_hubtaxa.pdf"), width=7, height=7/6)
par(mfrow=c(1,4), mar=c(2,4,1,0))
mycol<- c("green","red") 

p1 <- barplot(module11,
              las=2, border=NA, axes=T, cex.names=.5,ylim=c(0,130),main="Module1",
              col=mycol[1:nrow(module11)],ylab="Degree ")

## legend 
plot.new()
legend("top", bty="n", cex=0.5, #x.intersp=0.1, #y.intersp=0.75,
       legend=rev(rownames(module11)), 
       fill=rev(mycol), 
       border=rev(mycol) )
dev.off()





####module2
#读取数据
module2 <-read.table("module2_top20_hubtaxa.txt",sep="\t",header = TRUE,check.names=FALSE,comment.char="@")
module2
#数据处理一下，以适合用barplot绘制柱状图
module2 <-module2[order(module2[,2],decreasing =T),]
module2
module21 <-module2[,2:ncol(module2)]
rownames(module21)<-module2[,1]
module21<-as.matrix(module21)
module21
colnames(module21)


##module2核心微生物出图
pdf(paste0(output,"module2_top20_hubtaxa.pdf"), width=7, height=7/6)
par(mfrow=c(1,4), mar=c(2,4,1,0))
mycol<- c("green","red") 

p2 <- barplot(module21,
              las=2, border=NA, axes=T, cex.names=.5,ylim=c(0,130),main="Module2",
              col=mycol[1:nrow(module21)],ylab="Degree ")

## legend 
plot.new()
legend("top", bty="n", cex=0.5, #x.intersp=0.1, #y.intersp=0.75,
       legend=rev(rownames(module21)), 
       fill=rev(mycol), 
       border=rev(mycol) )
dev.off()





####module3
#读取数据
module3 <-read.table("module3_top20_hubtaxa.txt",sep="\t",header = TRUE,check.names=FALSE,comment.char="@")
module3
#数据处理一下，以适合用barplot绘制柱状图
module31 <-module3[,2:ncol(module3)]
rownames(module31)<-module3[,1]
module31<-as.matrix(module31)
module31
colnames(module31)


##module3核心微生物出图
pdf(paste0(output,"module3_top20_hubtaxa.pdf"), width=7, height=7/6)
par(mfrow=c(1,4), mar=c(2,4,1,0))
mycol<- c("green","red") 

p3 <- barplot(module31,
              las=2, border=NA, axes=T, cex.names=.5,ylim=c(0,130),main="Module3",
              col=mycol[1:nrow(module31)],ylab="Degree ")

## legend 
plot.new()
legend("top", bty="n", cex=0.5, #x.intersp=0.1, #y.intersp=0.75,
       legend=rev(rownames(module31)), 
       fill=rev(mycol), 
       border=rev(mycol) )
dev.off()







####module4
#读取数据
module4 <-read.table("module4_top20_hubtaxa.txt",sep="\t",header = TRUE,check.names=FALSE,comment.char="@")
module4
#数据处理一下，以适合用barplot绘制柱状图
module41 <-module4[,2:ncol(module4)]
rownames(module41)<-module4[,1]
module41<-as.matrix(module41)
module41
colnames(module41)


##module4核心微生物出图
pdf(paste0(output,"module4_top20_hubtaxa.pdf"), width=7, height=7/6)
par(mfrow=c(1,4), mar=c(2,4,1,0))
mycol<- c("green","red") 

p4 <- barplot(module41,
              las=2, border=NA, axes=T, cex.names=.5,ylim=c(0,130),main="Module4",
              col=mycol[1:nrow(module41)],ylab="Degree ")

## legend 
plot.new()
legend("top", bty="n", cex=0.5, #x.intersp=0.1, #y.intersp=0.75,
       legend=rev(rownames(module41)), 
       fill=rev(mycol), 
       border=rev(mycol) )
dev.off()







####module5
#读取数据
module5 <-read.table("module5_top20_hubtaxa.txt",sep="\t",header = TRUE,check.names=FALSE,comment.char="@")
module5
#数据处理一下，以适合用barplot绘制柱状图
module51 <-module5[,2:ncol(module5)]
rownames(module51)<-module5[,1]
module51<-as.matrix(module51)
module51
colnames(module51)


##module5核心微生物出图
pdf(paste0(output,"module5_top20_hubtaxa.pdf"), width=7, height=7/6)
par(mfrow=c(1,4), mar=c(2,4,1,0))
mycol<- c("green","red") 

p5 <- barplot(module51,
              las=2, border=NA, axes=T, cex.names=.5,ylim=c(0,130),main="Module5",
              col=mycol[1:nrow(module51)],ylab="Degree ")

## legend 
plot.new()
legend("top", bty="n", cex=0.5, #x.intersp=0.1, #y.intersp=0.75,
       legend=rev(rownames(module51)), 
       fill=rev(mycol), 
       border=rev(mycol) )
dev.off()







####module6
#读取数据
module6 <-read.table("module6_top20_hubtaxa.txt",sep="\t",header = TRUE,check.names=FALSE,comment.char="@")
module6
#数据处理一下，以适合用barplot绘制柱状图
module61 <-module6[,2:ncol(module6)]
rownames(module61)<-module6[,1]
module61<-as.matrix(module61)
module61
colnames(module61)


##module6核心微生物出图
pdf(paste0(output,"module6_top20_hubtaxa.pdf"), width=7, height=7/6)
par(mfrow=c(1,4), mar=c(2,4,1,0))
mycol<- c("green","red") 

p6 <- barplot(module61,
              las=2, border=NA, axes=T, cex.names=.5,ylim=c(0,130),main="Module6",
              col=mycol[1:nrow(module61)],ylab="Degree ")

## legend 
plot.new()
legend("top", bty="n", cex=0.5, #x.intersp=0.1, #y.intersp=0.75,
       legend=rev(rownames(module61)), 
       fill=rev(mycol), 
       border=rev(mycol) )
dev.off()


####module7没有核心物种




####module8
#读取数据
module8 <-read.table("module8_top20_hubtaxa.txt",sep="\t",header = TRUE,check.names=FALSE,comment.char="@")
module8
#数据处理一下，以适合用barplot绘制柱状图
module81 <-module8[,2:ncol(module8)]
rownames(module81)<-module8[,1]
module81<-as.matrix(module81)
module81
colnames(module81)


##module8核心微生物出图
pdf(paste0(output,"module8_top20_hubtaxa.pdf"), width=7, height=7/6)
par(mfrow=c(1,4), mar=c(2,4,1,0))
mycol<- c("green","red") 

p8 <- barplot(module81,
              las=2, border=NA, axes=T, cex.names=.5,ylim=c(0,130),main="Module8",
              col=mycol[1:nrow(module81)],ylab="Degree ")

## legend 
plot.new()
legend("top", bty="n", cex=0.5, #x.intersp=0.1, #y.intersp=0.75,
       legend=rev(rownames(module81)), 
       fill=rev(mycol), 
       border=rev(mycol) )
dev.off()







####module9
#读取数据
module9 <-read.table("module9_top20_hubtaxa.txt",sep="\t",header = TRUE,check.names=FALSE,comment.char="@")
module9
#数据处理一下，以适合用barplot绘制柱状图
module91 <-module9[,2:ncol(module9)]
rownames(module91)<-module9[,1]
module91<-as.matrix(module91)
module91
colnames(module91)


##module9核心微生物出图
pdf(paste0(output,"module9_top20_hubtaxa.pdf"), width=7, height=7/6)
par(mfrow=c(1,4), mar=c(2,4,1,0))
mycol<- c("green","red") 

p9 <- barplot(module91,
              las=2, border=NA, axes=T, cex.names=.5,ylim=c(0,130),main="Module9",
              col=mycol[1:nrow(module91)],ylab="Degree ")

## legend 
plot.new()
legend("top", bty="n", cex=0.5, #x.intersp=0.1, #y.intersp=0.75,
       legend=rev(rownames(module91)), 
       fill=rev(mycol), 
       border=rev(mycol) )
dev.off()







####module10
#读取数据
module10 <-read.table("module10_top20_hubtaxa.txt",sep="\t",header = TRUE,check.names=FALSE,comment.char="@")
module10
#数据处理一下，以适合用barplot绘制柱状图
module101 <-module10[,2:ncol(module10)]
rownames(module101)<-module10[,1]
module101<-as.matrix(module101)
module101
colnames(module101)


##module10核心微生物出图
pdf(paste0(output,"module10_top20_hubtaxa.pdf"), width=7, height=7/6)
par(mfrow=c(1,4), mar=c(2,4,1,0))
mycol<- c("green","red") 

p10 <- barplot(module101,
               las=2, border=NA, axes=T, cex.names=.5,ylim=c(0,130),main="Module10",
               col=mycol[1:nrow(module101)],ylab="Degree ")

## legend 
plot.new()
legend("top", bty="n", cex=0.5, #x.intersp=0.1, #y.intersp=0.75,
       legend=rev(rownames(module101)), 
       fill=rev(mycol), 
       border=rev(mycol) )
dev.off()







####module11
#读取数据
module11 <-read.table("module11_top20_hubtaxa.txt",sep="\t",header = TRUE,check.names=FALSE,comment.char="@")
module11
#数据处理一下，以适合用barplot绘制柱状图
module111 <-module11[,2:ncol(module11)]
rownames(module111)<-module11[,1]
module111<-as.matrix(module111)
module111
colnames(module111)


##module11核心微生物出图
pdf(paste0(output,"module11_top20_hubtaxa.pdf"), width=7, height=7/6)
par(mfrow=c(1,4), mar=c(2,4,1,0))
mycol<- c("green","red") 

p11 <- barplot(module111,
               las=2, border=NA, axes=T, cex.names=.5,ylim=c(0,130),main="Module11",
               col=mycol[1:nrow(module111)],ylab="Degree ")

## legend 
plot.new()
legend("top", bty="n", cex=0.5, #x.intersp=0.1, #y.intersp=0.75,
       legend=rev(rownames(module111)), 
       fill=rev(mycol), 
       border=rev(mycol) )
dev.off()






####科研白君的土壤世界——加权基因共表达网络分析


setwd("E:/博士数据2/小麦微生物组数据/customer_backup/otus/networks/")
rm(list = ls())
#----------import data and load packages--------------
library(tidyverse)
library(igraph)
library(ggraph)
#install.packages("colormap")
library(colormap)
#install.packages("wesanderson")
library(wesanderson)

###load function
source("calNetwork.R")
source("calNetModule.R")
source("calNetParameters.R")
source("extract_subNet.R")
source("cal_zi_pi.R")
bac <- read.table("otu_AAT_2.txt",header = T,sep= "\t",row.names = 1)
dim(bac)#12 418
head(bac[,1:6])

bac_net <- calNetwork(data = bac,method = "spearman",
                      cutoff_cor = 0.78, cutoff_p = 0.001)#0.6和0.05仅用于演示

dim(bac_net)#107
bac_net[1:6,1:4]
#2------------------node and graph-------------------
detach("package:Hmisc", unload = TRUE)
c(as.character(bac_net$from), as.character(bac_net$to)) %>%
  as_tibble() %>%
  group_by(value) %>%
  summarize(n=n()) -> bac_nodes

colnames(bac_nodes) <- c("name", "degree")

bac_graph <- graph_from_data_frame(bac_net, vertices = bac_nodes, directed = FALSE )
bac_graph

#3----------------network module----------------
bac_nodes2 <- calNetModule(data = bac_graph)
bac_nodes2

#4----------------network parameters----------------
calNetParameters(data = bac_graph)

#5----------------------extract subnetwork----------------------
bac_subNet <- extract_subNet(graph = bac_graph, data = bac, node = bac_nodes)
bac_subNet$sub_networkParameters#subnetwork parameters
bac_subNet$node#nodes
bac_subNet$edge#edges

#6----------------obtain adjacency matrix-----------------
adj2 <- data.frame(get.adjacency(bac_graph,sparse=FALSE))#获得邻接矩阵
adj2_w<-data.frame(get.adjacency(bac_graph,sparse=FALSE,attr="correlation"))#获得相关邻接矩阵

#7----------------------zi and pi-------------------
cal_zi_pi(edge = bac_net,node = bac_nodes)
#reference:
#https://github.com/cwatson/brainGraph/blob/master/R/vertex_roles.R




####R语言NetCoMi包
#install.packages("BiocManager")
#library(BiocManager)
#devtools::install_github("stefpeschel/NetCoMi", dependencies = TRUE,
#                         repos = c("https://cloud.r-project.org/",
#                                   BiocManager::repositories()))
#
#library(NetCoMi)









rm(list=ls())
setwd("E:/X-14952B/otus/networks2/R_Input1/")
output <- "E:/X-14952B/otus/networks2/R_Output1/"
#install.packages("BiodiversityR")
library(BiodiversityR) 
library(ggplot2)
library(vegan)
#install.packages("reshape2")
library(reshape2)
library(Hmisc)
#install.packages("plotrix")
library(plotrix)
#install.packages("phyloseq")
library(phyloseq)
library(MASS)
# library(bioDist) #package ???bioDist??? is not available (for R version 3.2.2)
library(igraph)
library(car)
#install.packages("coin")
library(coin)
#install.packages("edgeR")
library(edgeR)
#install.packages("formatR")
library(formatR)
library(gridExtra)
#install.packages("gplots")
library(gplots)
#install.packages("indicspecies")
library(indicspecies)
#install.packages("sciplot")
library(sciplot)
library(ape)
library(grid)
#install.packages("RVAideMemoire")
library(RVAideMemoire)
#install.packages("gridBase")
library(gridBase)
#install.packages("TukeyC")
library(TukeyC)
#install.packages("corrplot")
library(corrplot)
#install.packages("userfriendlyscience")
#install.packages("devtools")
#library(devtools)
#devtools::install_github("matherion/userfriendlyscience", dependencies=TRUE)
library(userfriendlyscience)
#install.packages("caret")
library(caret)
library(multcompView)
source("vennDia.R")
source("CorrDF.R")
source("plotOTU.R")
source("cor.mtest.R")
source("variance_functions.R")
source("maPalette.R")
source("triangle_shape.R")
source("star_shape.R")
debuggingState(on=FALSE)
options(scipen=10)


##### Import Data #####
#读入otu表格
otu_16s <- read.table("OTU_Tables.txt", row.names=1,sep="\t", header=T, blank.lines.skip=F, check.names=F)
otu_16s
#将OTU表格转为矩阵
otu_16s <- as.matrix(otu_16s)
otu_16s
dim(otu_16s)
#读入实验设计表格,特别需要注意的是，我们的实验处理虽然只有一个因素，但我们的实验设计表格一定得是两维度的，不然后边会报错
design_16s <- read.table("mapping.txt", row.names=1,sep="\t", header=T, blank.lines.skip=F, check.names=F)
design_16s
design_16s$Group <- factor(design_16s$Group ,c("Ggt","S_135","SA","H"))
design_16s$Treatment <- factor(design_16s$Treatment ,c("Ggt","S_135","SA","H"))
str(design_16s)


#读入分类表格,注意此时的物种分类表格已经将界、门、纲、目、科、属、种分开了
tax_16s <- read.table("taxonomy2.txt", row.names=1,sep="\t", header=T, blank.lines.skip=F, check.names=F)
tax_16s
tax_16s[tax_16s==""] <- "Unassigned"
colnames(tax_16s) <- c("Kingdom","Phylum","Class","Order", "Family", "Genus", "Species")
tax_16s

##注意：后续分析不要对变形菌门进行细分
# create separate taxonomy label specifying classes of Proteobacteria
tax_16s$labels <- tax_16s$Phylum
tax_16s[ rownames(tax_16s)[tax_16s$Class=="Alphaproteobacteria" ], ]$labels <- "Alphaproteobacteria"
#tax_16s[ rownames(tax_16s)[tax_16s$Class=="Betaproteobacteria" ], ]$labels <- "Betaproteobacteria"
tax_16s[ rownames(tax_16s)[tax_16s$Class=="Gammaproteobacteria" ], ]$labels <- "Gammaproteobacteria"
tax_16s[ rownames(tax_16s)[tax_16s$Class=="Deltaproteobacteria" ], ]$labels <- "Deltaproteobacteria"
table(tax_16s$labels)


##### Store Cyanobacteria and mitrochrondia sequeces #####
#剔除序列中与蓝藻细菌（Cyanobacteria）和线粒体（mitrochrondia）相关的序列
#unique函数主要用于查看数据框的整体情况
unique(tax_16s$Kingdom)
table(tax_16s$Kingdom)
##发现在Kingdom水平上存在一个unassigned的物种，故需要将其删除
r3 <- rownames(tax_16s[tax_16s$Kingdom=="Unassigned",])

table(tax_16s$Phylum)
r1 <- rownames(tax_16s[tax_16s$Phylum=="Cyanobacteria",])

unique(tax_16s$Family)
r2 <- c(rownames(tax_16s[tax_16s$Family=="mitochondria",]))

##创建需要移除的otu的集合
otus_remove_16s <- c(r1,r2,r3)
otus_remove_16s

## Remove these from otu table, tax table
otu_filter_16s <- otu_16s[-which(rownames(otu_16s) %in% otus_remove_16s),]
tax_filter_16s <- tax_16s[rownames(otu_filter_16s),]
design_filter_16s <- droplevels(design_16s[rownames(design_16s) %in% colnames(otu_filter_16s),])
design_filter_16s <- design_filter_16s[colnames(otu_filter_16s),]

#dim对数据框进行总体统计
dim(otu_filter_16s)
dim(tax_filter_16s)
dim(design_filter_16s)


###### 16S sequence and OTU counts ######
#对每个样本中所含otu进行数量进行统计，论文中描述：细菌群落中含有1,299,436条高质量的序列
sum(colSums(otu_filter_16s))
#按每列otu的总和进行样本排序，序列范围为32990~233359条
sort(colSums(otu_filter_16s))
#求每列样本otu总和 的平均值，平均值为90934
median(colSums(otu_filter_16s))

#16S总物种数量为1621个
nrow(tax_filter_16s)
#其中古菌1个，细菌1620个
table(tax_filter_16s$Kingdom)

## Export filtered OTU table for multiple rarefactions in QIIME
write.table(otu_filter_16s,paste(output,"otu_filter_16s.txt",sep=""),sep="\t",quote=F)

## Order taxonmy file by OTU
#将otu表格与物种分类表格的otu名称进行匹配
otu_order_16s <- match(rownames(otu_filter_16s), rownames(tax_filter_16s))
tax_filter_16s <- tax_filter_16s[otu_order_16s,]
write.table(tax_filter_16s,paste(output,"tax_filter_16s.txt",sep=""),sep="\t",quote=F)

##### Calculate % sequences removed #####
cyano_otu <- otu_16s[r1,]; dim(cyano_otu); length(r1)

sum(colSums(cyano_otu))/sum(colSums(otu_16s))*100

mito_otu <- otu_16s[r2,]; dim(mito_otu); length(r2)

sum(sum(mito_otu))/sum(colSums(otu_16s))*100


##### Plot sequence counts across sample types #####
par(las=2,mar=c(6,6,1,5))
plot((1:4),ylim=c(0,30000),type="n",ann=F,axes=F)

boxplot(colSums(otu_filter_16s[,rownames(design_filter_16s)[which(design_filter_16s$Group=="Ggt")]]),add=T,at=1,ann=F,axes=F)
stripchart(colSums(otu_filter_16s[,rownames(design_filter_16s)[which(design_filter_16s$Group=="Ggt")]]),pch=19,add=T,vertical=T,at=1)

boxplot(colSums(otu_filter_16s[,rownames(design_filter_16s)[which(design_filter_16s$Group=="S_135")]]),add=T,at=2,ann=F,axes=F)
stripchart(colSums(otu_filter_16s[,rownames(design_filter_16s)[which(design_filter_16s$Group=="S_135")]]),pch=19,add=T,vertical=T,at=2)

boxplot(colSums(otu_filter_16s[,rownames(design_filter_16s)[which(design_filter_16s$Group=="SA")]]),add=T,at=3,ann=F,axes=F)
stripchart(colSums(otu_filter_16s[,rownames(design_filter_16s)[which(design_filter_16s$Group=="SA")]]),pch=19,add=T,vertical=T,at=3)

boxplot(colSums(otu_filter_16s[,rownames(design_filter_16s)[which(design_filter_16s$Group=="H")]]),add=T,at=4,ann=F,axes=F)
stripchart(colSums(otu_filter_16s[,rownames(design_filter_16s)[which(design_filter_16s$Group=="H")]]),pch=19,add=T,vertical=T,at=4)


axis(1,at=seq(1,12,1),labels=FALSE)
staxlab(side=1,at=1:4,labels=levels(design_filter_16s$TT),srt=45)
axis(2,at=seq(0,30000,5000))
title(ylab="Number of Sequences",line=5)
title(main="Raw sequence counts")

kruskal.test(colSums(otu_filter_16s)~design_filter_16s$Group)



#################################
##### 2.Analysis and Figures #####
################################

#################################
##### 2.1 microbial component 物种组成分析 #####
################################

###对OTU表格进行标准化
##### TMM normalize 16S/ITS counts for whole community beta diversity analysis #####

## Apply TMM normalization to entire 16s data set and create phyloseq objects for later analysis
group_16s <- design_filter_16s$Group
group_16s
edgeR_16s <- DGEList(counts=otu_filter_16s, 
                     group=design_filter_16s$Group, 
                     genes=tax_filter_16s)

edgeR_16s

edgeR_16s <- calcNormFactors(edgeR_16s)
edgeR_16s

## Extract normalized counts
otu_norm_16s <- cpm(edgeR_16s, normalized.lib.sizes=T, log=F)

## Create phyloseq objects
physeq_16s_norm <- phyloseq(otu_table(otu_norm_16s, taxa_are_rows=T),
                            tax_table(as.matrix(tax_filter_16s)),
                            sample_data(design_filter_16s))


## Create bray-curtis dissimiliartiy matrix
all_dis_16s <- vegdist(t(otu_table(physeq_16s_norm)),method="bray")
all_dis_16s

##### 16s overall PERMANOVA #####
## Perform PERMANVOA testing for sample type and cropping system effects 
paov_all_16s <- adonis(t(otu_table(physeq_16s_norm)) ~Group, data=design_filter_16s, permutations=9999)
paov_all_16s

##### Supplementary Table S2: global PERMANOVA #####
##注意：切记此时要手动将置换多元方差分析结果粘贴到excel表格
paov_all_16s  



##### Supplementary Figure S3: phyla relative abundance plots #####

##### Bacteria

## Express 16S OTU counts as relative abunance percent
otu_16s_RA <- t(t(otu_filter_16s)/colSums(otu_filter_16s)) * 100
otu_16s_RA
colSums(otu_16s_RA)
nrow(otu_16s_RA)

## Get names of bacteria phyla present (use 'labels' as this specifies class within Proteobacteria)
PHYLAnames_16s <- names(sort(table(tax_filter_16s[,"labels"]), decr=T))
PHYLAnames_16s
length(PHYLAnames_16s)
sort(table(tax_filter_16s[,"labels"]), decr=T)

## Preparation of matrix with relative abundance by phylum
y <- NULL
otunames <- rownames(otu_16s_RA)
for (i in PHYLAnames_16s){
  x <- array(colSums(otu_16s_RA[rownames(tax_filter_16s)[which(tax_filter_16s$labels == paste(i))],,drop=FALSE]))
  y <- rbind(y,x)
}

## Create matrix
rownames(y) <- paste(PHYLAnames_16s)
colnames(y) <- paste(colnames(otu_16s_RA))
PHYLUM_mat_16s <- y
PHYLUM_mat_16s[,1:3]
colSums(PHYLUM_mat_16s)
PHYLUM_mat_16s_mean <- sort(apply(PHYLUM_mat_16s,1,mean),decr=T)
PHYLUM_mat_16s <- PHYLUM_mat_16s[names(PHYLUM_mat_16s_mean),]

## Use hierarchical clustering to order samples
norm_clu_16s <- hclust(all_dis_16s, "average")
norm_clu_16s$height
Beta_labs_16s <- norm_clu_16s$labels
Beta_order_16s <- norm_clu_16s$order
Beta_labs_16s <- Beta_labs_16s[Beta_order_16s]
PHYLUM_mat_16s <- PHYLUM_mat_16s[ ,Beta_labs_16s]


### Defining bOTU colors by phylum (using the taxonomy file)
tax_filter_16s$cols <- tax_filter_16s$labels
table(tax_filter_16s$cols)

# Phyla with MEAN abundances lower than 1% relative abundances
table(apply(PHYLUM_mat_16s, 1, mean) < 1)
low_count_phyla_16s <- rownames(PHYLUM_mat_16s)[sort(apply(PHYLUM_mat_16s, 1, mean), decr=T) < 1]
low_count_phyla_16s

# attribute grey color,
##此处特别注意，给稀有物种进行颜色赋值的时候，要注意我们的变量phylum的大小写细节，不然无法给稀有物种进行颜色赋值
for(i in low_count_phyla_16s){
  tax_filter_16s[ rownames(tax_filter_16s)[tax_filter_16s$Phylum==paste(i) ], ]$cols <- "lightgrey"
}
table(tax_filter_16s$cols)

# Phyla with MEAN abundances higher than 1% relative abundances
abundant_phyla_16s <- rownames(PHYLUM_mat_16s)[sort(apply(PHYLUM_mat_16s, 1, mean), decr=T) > 1 ]
abundant_phyla_16s
tax_filter_16s[ rownames(tax_filter_16s)[tax_filter_16s$labels=="Gammaproteobacteria" ], ]$cols <- "firebrick2"
tax_filter_16s[ rownames(tax_filter_16s)[tax_filter_16s$labels=="Acidobacteria" ], ]$cols <- "tan2"
tax_filter_16s[ rownames(tax_filter_16s)[tax_filter_16s$labels=="Bacteroidetes" ], ]$cols <- "palegreen2"
tax_filter_16s[ rownames(tax_filter_16s)[tax_filter_16s$labels=="Alphaproteobacteria" ], ]$cols <- "palegreen4"
tax_filter_16s[ rownames(tax_filter_16s)[tax_filter_16s$labels=="Chloroflexi" ], ]$cols <- "steelblue1"
tax_filter_16s[ rownames(tax_filter_16s)[tax_filter_16s$labels=="Actinobacteria" ], ]$cols <- "tan1"
tax_filter_16s[ rownames(tax_filter_16s)[tax_filter_16s$labels=="Firmicutes" ], ]$cols <- "lightsalmon4"
tax_filter_16s[ rownames(tax_filter_16s)[tax_filter_16s$labels=="Gemmatimonadetes" ], ]$cols <- "firebrick1"
tax_filter_16s[ rownames(tax_filter_16s)[tax_filter_16s$labels=="Fusobacteria" ], ]$cols <- "orchid3"
tax_filter_16s[ rownames(tax_filter_16s)[tax_filter_16s$labels=="Deltaproteobacteria" ], ]$cols <- "palevioletred2"
tax_filter_16s[ rownames(tax_filter_16s)[tax_filter_16s$labels=="Patescibacteria" ], ]$cols <- "peachpuff3"
tax_filter_16s[ rownames(tax_filter_16s)[tax_filter_16s$labels=="Verrucomicrobia" ], ]$cols <- "peachpuff4"

## collaps OTU colors to prepare Phylum level colors
label_cols_16s <- tax_filter_16s[, c("labels", "cols") ]
label_cols_16s
library(plyr)
PHYLA_label_cols_16s <- ddply(label_cols_16s, .variables="cols", .fun=unique)
rownames(PHYLA_label_cols_16s) <- PHYLA_label_cols_16s[,1]
PHYLA_label_cols_16s <- PHYLA_label_cols_16s[c(abundant_phyla_16s, low_count_phyla_16s),]
PHYLA_label_cols_16s

## Legend for Phylum colors
PHYLA_label_cols_16s_legend <- PHYLA_label_cols_16s[1:13,]
PHYLA_label_cols_16s_legend[13,1] <- "other"
rownames(PHYLA_label_cols_16s_legend)[13] <- "other"
PHYLA_label_cols_16s_legend




##### Plot Supplementary Figure S3
pdf(paste0(output,"FigureS3.pdf"), encoding="MacRoman", height=6, width=10, paper="a4r")

layout(matrix(c(1,2),1,2, byrow=F))

par(oma=c(0,0,0,0), mar=c(6,2,1,5), xpd=NA)
phylum_bar_16s <- barplot(as.matrix(PHYLUM_mat_16s), col=PHYLA_label_cols_16s[rownames(PHYLUM_mat_16s),]$cols,
                          ylim=c(0,100), xaxt="n", border=NA, las=2)
axis(1, at=phylum_bar_16s, labels=design_16s[Beta_labs_16s,]$Group, col.axis="black", las=2, cex.axis=0.6)
title(ylab="Relative abundance (%)")
title(main="Bacteria Community")
legend(43, 100, bty="n", cex=0.7, x.intersp=0.1, y.intersp=1,
       legend=rev(PHYLA_label_cols_16s_legend$labels), 
       fill=rev(PHYLA_label_cols_16s_legend$cols), 
       border=rev(PHYLA_label_cols_16s_legend$cols) )
dev.off()


#################################
##### 2.3 Beta diversity β-多样性分析 #####
################################


##### Supplementary Figure S4: Beta diversity (unconstrained PCoA) of bacteria/fungi soil and root communities #####

## Perform unconstrained PCoA of entire bacteria data set
pcoa_norm_16s <- ordinate(physeq_16s_norm,"PCoA","bray")
pcoa_all_16s <- plot_ordination(physeq_16s_norm, pcoa_norm_16s, type="sites", color="Group", shape="Group")
pcoa_all_16s <- pcoa_all_16s+
  geom_point(size=5)+
  scale_shape_manual(name="Group",values=c(15,16,17,18,15,16,17,18,15,16,17,18))+
  scale_color_manual(name="Group",values=c(rep("tan4",4),rep("darkgreen",4),rep("red",4)))+
  xlab(paste("PCo 1", paste("(",round(pcoa_norm_16s$values[1,2]*100,1),"%",")",sep=""),sep=" "))+
  ylab(paste("PCo 2", paste("(",round(pcoa_norm_16s$values[2,2]*100,1),"%",")",sep=""),sep=" "))+
  theme(legend.position="bottom",legend.title=element_blank(),legend.key = element_blank())+
  guides(color=guide_legend(nrow=2,byrow=TRUE))+
  guides(shape=guide_legend(nrow=2,byrow=TRUE))+
  theme(plot.title = element_text(face="bold", hjust = 0.5))+
  ggtitle("Bacteria Community")
pcoa_all_16s


## Plot
plot.new()
pdf(paste0(output,"FigureS4.pdf"), encoding="MacRoman", onefile=F, height=10, width=10)
grid.arrange(pcoa_all_16s,nrow=2, ncol=2)
dev.off()


##### Figure 2 Quantifying cropping system effects in bulk soil bacteria and fungi communities with CAP #####


### Rhizosphere BACTERIA
## Perform CAP analysis constrained to cropping system on Rhizosphere soil bacteria community
cap_16s <- ordinate(physeq_16s_norm,"CAP","bray", ~Group)

## Calculate proporition of variation table
var_tab_16s <- variability_table(cap_16s)

## Set seed for reproducibilty of permutational testing and permute CAP ordination for significance testing
##注意：此处生成的结果要手动复制粘贴到excel，用于后期文章的写作
set.seed(8046)  # zip code Agroscope ZH
permutest(cap_16s, permutations=how(nperm=9999))

cap_16s <- plot_ordination(physeq_16s_norm,cap_16s,color="Group",shape="Group")+
  geom_point(size=5)+
  scale_shape_manual(values=c(15,16,17,18))+
  scale_color_manual(values=c("dodgerblue4","dodgerblue1","firebrick4","firebrick2"))+
  theme(legend.position="bottom",legend.title=element_blank(),legend.key = element_blank())+
  guides(shape=guide_legend(nrow=1,byrow=TRUE))+
  ggtitle(bquote(paste(.(round(var_tab_16s["constrained","proportion"],2)*100),"% of total variation (CI=",
                       .(round(cca_ci(cap_16s)[1],2)*100),"%, ",
                       .(round(cca_ci(cap_16s)[2],2)*100),"%); ", "p=0.001" )))+
  theme(plot.title = element_text(size = 10, face = "bold"))

cap_16s


### Plot Figure: CAP analysis of rhizosphere soil bacteria and fungi 
plot.new()
pdf(paste0(output,"Figure_CAPs.pdf"),encoding="MacRoman",onefile=F,height=10,width=5)
grid.arrange(cap_16s,nrow=2, ncol=1)
dev.off()

##### Supplementary Table S4 (rhizosphere): statistical testing of cropping system effects on rhizosphere soil bacteria and fungi with PERMANOVA & BETADISP #####
##此处所得到的数据需要手动的复制粘贴到excel
## Perform PERMANVOA testing for cropping system effects on rhizosphere soil bacteria community
dis_16s <- vegdist(t(otu_table(physeq_16s_norm)),method="bray")
paov_16s <- adonis(t(otu_table(physeq_16s_norm)) ~ Group, data=design_16s, permutations=9999, method= "bray")
paov_16s

## Perform pairwise comparisions
pairwise.perm.manova(t(otu_table(physeq_16s_norm)), design_16s$Group, nperm=9999)

## Perform BETADISP test of multivariate dispersions
bdi_16s <- betadisper(dis_16s,design_16s$Group,type = "centroid")
permutest(bdi_16s,pairwise=T,permutations=how(nperm=9999))


##### Identifiying cropping system responsive OTUs with indicator species analysis #####

## Identify indicator species in root-associated bacteria communities
indic_root_16s <- as.data.frame(t(otu_norm_16s))
indic_root_groups_16s <- design_filter_16s$Group
length(unique(indic_root_groups_16s))

## Define indicator species for root bacteria community.
## Note: These calculations can be time and processor intensive
set.seed(8046)
indicatorsp_root_16s <- multipatt(indic_root_16s,indic_root_groups_16s,func = "r.g",control=how(nperm=9999))
summary(indicatorsp_root_16s,alpha=0.05,indvalcomp=T)
indicatorsp_root_16s
indic_root_df_16s <- indicatorsp_root_16s$sign
indic_root_df_16s

##indicator阈值
net_16s <- as.matrix(indic_root_df_16s[which(indic_root_df_16s$p.value<0.05),])
net_16s
#获得indicator物种表
indicator_taxa <- subset(tax_filter_16s,rownames(tax_filter_16s) %in% rownames(net_16s))
indicator_taxa
indictor <- cbind(net_16s,indicator_taxa)
indictor
write.table(indictor,paste0(output,"indictor.txt"),sep="\t",quote=F)

## Import data frame of indicator species to save time
indic_root_df_16s <- read.table("indictor.txt",header=T,sep="\t")

Ggt_root_16s <- as.matrix(indic_root_df_16s[which(indic_root_df_16s$s.Ggt == 1 & indic_root_df_16s$p.value < 0.05),])
Ggt_root_16s
S_135_root_16s <- as.matrix(indic_root_df_16s[which(indic_root_df_16s$s.S_135 == 1 & indic_root_df_16s$p.value < 0.05),])
S_135_root_16s
SA_root_16s <- as.matrix(indic_root_df_16s[which(indic_root_df_16s$s.SA == 1 & indic_root_df_16s$p.value < 0.05),])
SA_root_16s
H_root_16s <- as.matrix(indic_root_df_16s[which(indic_root_df_16s$s.H == 1 & indic_root_df_16s$p.value < 0.05),])
H_root_16s

root_r_values_16s <- rbind(Ggt_root_16s, S_135_root_16s, SA_root_16s, H_root_16s)
colnames(root_r_values_16s)[1:4] <-c("Ggt","S_135","SA","H")
root_r_values_16s

## Range of correlation coefficients
range(root_r_values_16s[,"stat"])

## Total number of indicator OTUS
length(unique(rownames(root_r_values_16s)))

## Proportion of root associated bacteria OTUs responding to cropping system
length(unique(rownames(root_r_values_16s)))/nrow(otu_norm_16s)

## Proportion of root-associated bacteria sequences responding to cropping system
root_16s_ra <- t(t(otu_filter_16s)/colSums(otu_filter_16s)) * 100
sum(colSums(root_16s_ra[unique(rownames(root_r_values_16s)),]))/sum(colSums(root_16s_ra))







##### Figure 3: Bipartite networks of OTUs associated with cropping systems #####

### Root associated BACTERIA community
## Construct node table for root-associated bacteria communities from indicator species data
root_bipartite_16s <- data.frame(from= c(rep("Ggt",length(which(root_r_values_16s[,"Ggt"]==1))),
                                         rep("S_135",length(which(root_r_values_16s[,"S_135"]==1))),
                                         rep("SA",length(which(root_r_values_16s[,"SA"]==1))),
                                         rep("H",length(which(root_r_values_16s[,"H"]==1)))),
                                 to= c(rownames(root_r_values_16s)[which(root_r_values_16s[,"Ggt"]==1)],
                                       rownames(root_r_values_16s)[which(root_r_values_16s[,"S_135"]==1)],
                                       rownames(root_r_values_16s)[which(root_r_values_16s[,"SA"]==1)],
                                       rownames(root_r_values_16s)[which(root_r_values_16s[,"H"]==1)]),
                                 r= c(root_r_values_16s[which(root_r_values_16s[,"Ggt"]==1),"stat"],
                                      root_r_values_16s[which(root_r_values_16s[,"S_135"]==1),"stat"],
                                      root_r_values_16s[which(root_r_values_16s[,"SA"]==1),"stat"],
                                      root_r_values_16s[which(root_r_values_16s[,"H"]==1),"stat"]))
##注意，此时计算构建的attribute表格不是用于二分网络的构建，而是
## make node attribute table for each OTU
root_bipartite_attrib_16s <- data.frame(node=unique(rownames(root_r_values_16s)),indicgroup=0)

for (i in as.character(root_bipartite_attrib_16s$node))
{
  root_bipartite_attrib_16s[root_bipartite_attrib_16s$node==i,"indicgroup"] <- paste(colnames(root_r_values_16s)[which(root_r_values_16s[i,1:4]==1)],collapse = "_")
}

root_bipartite_attrib_16s <- cbind(root_bipartite_attrib_16s,tax_filter_16s[as.character(root_bipartite_attrib_16s$node),])

## Create bipartite network with igraph
root_bi_16s <- graph.data.frame(root_bipartite_16s,directed=F)
V(root_bi_16s)$type <- V(root_bi_16s)$name %in% root_bipartite_16s[,1]
root_bi_16s <- simplify(root_bi_16s, remove.multiple=T, remove.loops=T)


## Set node labels
V(root_bi_16s)$label <- V(root_bi_16s)$name
V(root_bi_16s)$label <- gsub("OTU*",NA,V(root_bi_16s)$label)

## Set node sizes,
#注意，此处的（10，4）中的10是节点的大小，4代表自己试验处理的数目
#注意，此处的（2,777）中的2是节点的大小，777代表指示物种的数目，注意length(unique(rownames(root_r_values_16s)))运行后的值
V(root_bi_16s)$size <- c(rep(10,4),rep(2,765))

## Set node shape
V(root_bi_16s)$shape <- c(rep("circle",4), rep("circle",765))

## Define node colors based upon phylum/class taxonomy assignment

V(root_bi_16s)$color <- V(root_bi_16s)$name
V(root_bi_16s)$color[1:4] <- "white"
V(root_bi_16s)$color <- tax_filter_16s[ V(root_bi_16s)$name, ]$cols
V(root_bi_16s)$frame.color <- V(root_bi_16s)$color
V(root_bi_16s)$frame.color

##### Plot Figure 3 : Soil/Root Bipartite Networks
set.seed(8046)
root_layout_16s <- layout_with_fr(root_bi_16s, niter=9999)


#graphml 格式，可使用 gephi 软件打开并进行可视化编辑
write.graph(root_bi_16s, 'Bipartite networks.graphml', format = 'graphml')

#gml 格式，可使用 cytoscape 软件打开并进行可视化编辑
write.graph(root_bi_16s, 'Bipartite networks.gml', format = 'gml')



pdf(paste0(output,"Figure3.pdf"), encoding="MacRoman", height=6, width=6)
layout(matrix(c(1,2,3,4,5,6), nrow=3, byrow=T), c(2,2), c(4,4,3))

par(mar=c(0,0,0,0), oma=c(0,0,0,0))
plot(root_bi_16s, vertex.label.cex=2, layout=root_layout_16s, asp=0)

plot(1, type="n", ann=F, axes=F)
legend("center", pch=19, bty="n", ncol=2, horiz=F, x.intersp=0.4, 
       legend=PHYLA_label_cols_16s_legend$labels, 
       col=PHYLA_label_cols_16s_legend$cols)

dev.off()



##### Cropping system responsive OTUs with edgeR #####
## Get cropping system responsive root-associated bacteria OTUs with likelihood ratio testing in edgeR
model_matroot_16s <- model.matrix(~ Group, data=design_filter_16s)
model_matroot_16s
edgeR_16s_root_fsyst <- DGEList(counts=otu_filter_16s, group=design_filter_16s$Group, genes=tax_filter_16s)
edgeR_16s_root_fsyst <- calcNormFactors(edgeR_16s_root_fsyst)

dge_rootfsyst_16s <- estimateGLMRobustDisp(edgeR_16s_root_fsyst, design=model_matroot_16s)

fit_rootfsyst_16s <- glmFit(dge_rootfsyst_16s, design=model_matroot_16s)
lrt_rootfsyst_16s <- glmLRT(fit_rootfsyst_16s, coef=2:4)
fsyst_root_16s <- topTags(lrt_rootfsyst_16s, n=Inf, p.value=0.05)
fsyst_root_16s <- fsyst_root_16s$table
fsyst_root_16s
write.table(fsyst_root_16s,paste0(output,"edgeR.txt"),sep="\t",quote=F)



##### Figure S6:  management sensitive OTUs validated by both indicator species and edgeR #####
## Management sensitive root-associated bacteria OTus
indic_edge_16s_root <- intersect(rownames(root_r_values_16s),rownames(fsyst_root_16s))
sum(colSums(root_16s_ra[indic_edge_16s_root,]))/sum(colSums(root_16s_ra))*100


###获得具有显著差异的indicator表
indicator <- subset(net_16s,rownames(net_16s) %in% indic_edge_16s_root)
indicator

indicator_taxa <- subset(tax_filter_16s,rownames(tax_filter_16s) %in% rownames(indicator))
indicator_taxa

indicator2 <- cbind(indicator, indicator_taxa)
indicator2

write.table(indicator2,paste0(output,"sensitiv-otu.txt"),sep="\t",quote=F)

## Plot Supplementary sensitive otu
pdf(paste0(output,"sensitiv-otu.pdf"),height=10, width=10)
par(mfrow=c(3,2),xpd=NA)


venndiagram(rownames(root_r_values_16s),rownames(fsyst_root_16s),type=2,
            printsub=F,lcol="black",tcol="black",diacol="black",lines="black",
            labels=c("Indicator Species Analysis","edgeR Analysis"),title="Treatment Sensitive \nRhizosphere Bacteria OTUs")

dev.off()









## Define root samples by management system
Ggt_samples_root <- c("Ggt1","Ggt2","Ggt3")
#rownames(design_filter_16s[design_filter_16s$Group=="Ggt",])
S_135_samples_root <- c("S_135_1","S_135_2","S_135_3")
SA_samples_root <- c("SA1","SA2","SA3")
H_samples_root <- c("H1","H2","H3")

## Root bacteria: Calculate percentage of sequences classified into each phylum
otu_16s_root_mso <- otu_filter_16s[indic_edge_16s_root,]
otu_16s_root_mso_ra <- t(t(otu_16s_root_mso)/colSums(otu_16s_root_mso)) * 100

PHYLAnames_root_mso_16s <- names(sort(table(tax_filter_16s[indic_edge_16s_root,]$Phylum),decr=T))
PHYLAnames_root_mso_16s
length(PHYLAnames_root_mso_16s)

## Preparation of matrix with relative abundance by phylum
y <- NULL
otunames <- rownames(otu_16s_root_mso_ra)
for (i in PHYLAnames_root_mso_16s){
  x <- array(colSums(otu_16s_root_mso_ra[rownames(tax_filter_16s[indic_edge_16s_root,])[which(tax_filter_16s[indic_edge_16s_root,]$Phylum == paste(i))], ,drop=FALSE]))
  y <- rbind(y,x)
}

## Create matrix
rownames(y) <- paste(PHYLAnames_root_mso_16s)
colnames(y) <- paste(colnames(otu_16s_root_mso_ra))
PHYLUM_mat_otu_root_16s_mso <- y
PHYLUM_mat_otu_root_16s_mso_mean <- sort(apply(PHYLUM_mat_otu_root_16s_mso,1,mean),decr=T)
PHYLUM_mat_otu_root_16s_mso <- PHYLUM_mat_otu_root_16s_mso[names(PHYLUM_mat_otu_root_16s_mso_mean),]

## Get sequences abundances by phylum
PHYLUM_mat_otu_root_16s_mso_mean

## Phylum matrix for heatmap of phyla abundances across 
root_16s_mso_ra <- otu_norm_16s[indic_edge_16s_root,]

## Preparation of matrix with relative abundance by phylum
y <- NULL
otunames <- rownames(root_16s_mso_ra)
for (i in PHYLAnames_root_mso_16s){
  x <- array(colSums(root_16s_mso_ra[rownames(tax_filter_16s[indic_edge_16s_root,])[which(tax_filter_16s[indic_edge_16s_root,]$Phylum == paste(i))], ,drop=FALSE]))
  y <- rbind(y,x)
}

## Create matrix
rownames(y) <- paste(PHYLAnames_root_mso_16s)
colnames(y) <- paste(colnames(root_16s_mso_ra))
PHYLUM_mat_root_16s_mso <- y
PHYLUM_mat_root_16s_mso_mean <- sort(apply(PHYLUM_mat_root_16s_mso,1,mean),decr=T)
PHYLUM_mat_root_16s_mso <- PHYLUM_mat_root_16s_mso[names(PHYLUM_mat_root_16s_mso_mean),]

## Make matrix of phyla abundances by cropping system
root_16s_mso_phylum <- as.matrix(cbind(`Ggt`=apply(PHYLUM_mat_root_16s_mso[,Ggt_samples_root],1,mean),
                                       `S_135`=apply(PHYLUM_mat_root_16s_mso[,S_135_samples_root],1,mean),
                                       `SA`=apply(PHYLUM_mat_root_16s_mso[,SA_samples_root],1,mean),
                                       `H`=apply(PHYLUM_mat_root_16s_mso[,H_samples_root],1,mean)))


## Plot Figure S7c
pdf(paste0(output,"FigureS7c.pdf"),height=10, width=10)
par(las=1)
colors <- maPalette(l="lightgray",m="yellow",h="navy", k=60)
heatmap.2(log2(root_16s_mso_phylum+1),col=colors,Rowv=F,Colv=F,margins=c(5,12),scale="none", trace="none"
          ,density.info="none",key.title=NA,srtCol=45,
          labRow=rownames(root_16s_mso_phylum),
          dendrogram="none",key.xlab=expression(paste("Relative OTU abundance"," (log"[2]," CPM)")))
dev.off()



## 16s root-associated management sensitive OTU abundances heatmap 
root_16s_mso <- as.matrix(cbind(`Ggt`=apply(otu_norm_16s[indic_edge_16s_root,Ggt_samples_root],1,mean),
                                `S_135`=apply(otu_norm_16s[indic_edge_16s_root,S_135_samples_root],1,mean),
                                `SA`=apply(otu_norm_16s[indic_edge_16s_root,SA_samples_root],1,mean),
                                `H`=apply(otu_norm_16s[indic_edge_16s_root,H_samples_root],1,mean)))

table(tax_filter_16s[indic_edge_16s_root,]$Phylum)
root_16s_sidecol <- tax_filter_16s[indic_edge_16s_root,]$Phylum
names(root_16s_sidecol) <- root_16s_sidecol

root_16s_sidecol[root_16s_sidecol == "Acidobacteria"] <- "seashell4"
root_16s_sidecol[root_16s_sidecol == "Actinobacteria"] <- "red"
root_16s_sidecol[root_16s_sidecol == "Armatimonadetes"] <- "darkslateblue"
root_16s_sidecol[root_16s_sidecol == "Bacteroidetes"] <- "tan4"
root_16s_sidecol[root_16s_sidecol == "Chloroflexi"] <- "blue"
root_16s_sidecol[root_16s_sidecol == "Epsilonbacteraeota"] <- "violetred4"
root_16s_sidecol[root_16s_sidecol == "Euryarchaeota"] <- "orchid"
root_16s_sidecol[root_16s_sidecol == "Fibrobacteres"] <- "navajowhite4"
root_16s_sidecol[root_16s_sidecol == "Firmicutes"] <- "palegreen"
root_16s_sidecol[root_16s_sidecol == "Fusobacteria"] <- "aquamarine4"
root_16s_sidecol[root_16s_sidecol == "Gemmatimonadetes"] <- "darkslateblue"
root_16s_sidecol[root_16s_sidecol == "Patescibacteria"] <- "tan2"
root_16s_sidecol[root_16s_sidecol == "Planctomycetes"] <- "blue2"
root_16s_sidecol[root_16s_sidecol == "Proteobacteria"] <- "violetred2"
root_16s_sidecol[root_16s_sidecol == "Rokubacteria"] <- "orchid3"
root_16s_sidecol[root_16s_sidecol == "Spirochaetes"] <- "navajowhite2"
root_16s_sidecol[root_16s_sidecol == "Verrucomicrobia"] <- "palegreen4"
root_16s_sidecol[root_16s_sidecol == "WPS-2"] <- "aquamarine1"

## Plot Figure S8c
pdf(paste0(output,"FigureS8c.pdf"),height=10, width=10)
par(las=1)
colors <- maPalette(l="lightgray",m="yellow",h="navy", k=60)
heatmap.2(log2(root_16s_mso+1),col=colors,Rowv=F,Colv=F,margins=c(5,12),scale="none", trace="none"
          ,density.info="none",key.title=NA,srtCol=45,
          labRow=paste(rownames(root_16s_mso),tax_filter_16s[indic_edge_16s_root,]$Family), RowSideColors = root_16s_sidecol,
          dendrogram="none",key.xlab=expression(paste("Relative OTU abundance"," (log"[2]," CPM)")))
dev.off()


## Make vector of OTUs mentioned in text for box plots
root16s_otu <-c("OTU_247","OTU_695","OTU_105","OTU_46","OTU_80","OTU_127","OTU_758","OTU_2632")

## Plot Figure S11
pdf(paste0(output,"FigureS11.pdf"),encoding="MacRoman",height=10,width=10)
par(mfrow=c(3,3))
for(i in root16s_otu)
  plotOTU(sample_table = design_filter_16s, count_matrix = otu_norm_16s,
          OTU = i, sample_factor = "Group",
          main = paste(i, tax_filter_16s[i,"Family"]), ylab = "CPM",
          type = "box")
dev.off()






##### Individual co-occurence network creation and analysis  #####

## Create root-associated bacteria co-occurrence network with igraph package ##

## Define cropping system senstive bulk soil fungi OTUs (cropping system response validated by both statisical methods)
## This is important for highlighting their position in the subsequent figures
indic_edge_16s_root <- intersect(rownames(root_r_values_16s),rownames(fsyst_root_16s))
indic_edge_16s_root
## Perform pairwise Spearman correlations on bulk root bacteria community on TMM normalized CPM counts
root_16s_otu_cor <- rcorr(t(otu_norm_16s),type=c("spearman"))
root_16s_otu_cor



## Create data frame of co-occurring OTUs
root_16s_cor_df <- CorrDF(root_16s_otu_cor$r,root_16s_otu_cor$P)
root_16s_cor_df
root_16s_cor_df$corr <- root_16s_cor_df$cor
root_16s_cor_df$corr[root_16s_cor_df$corr > 0] <- 1
root_16s_cor_df$corr[root_16s_cor_df$corr < 0] <- -1
root_16s_cor_df
root_16s_cor_df$weight <- abs(root_16s_cor_df$cor)
root_16s_cor_df
root_16s_cor_df$padj <- p.adjust(root_16s_cor_df$p, method="none")
root_16s_cor_df
write.table(root_16s_cor_df,paste0(output,"root_16s_cor_df.txt"),sep="\t",quote=F)


## Subset data frame for co-occurring OTUs with Spearman's rho > 0.7 and a p-value < 0.001
#为了能够同时体现正负相关关系，所以此处要对相关系数取绝对值，采用语法abs(correlation)>0.7
root_16s_cor_df_padj <- root_16s_cor_df[which(abs(root_16s_cor_df$cor) > 0.7),]
root_16s_cor_df_padj <- root_16s_cor_df_padj[which(root_16s_cor_df_padj$padj < 0.001),]


write.table(root_16s_cor_df_padj,paste0(output,"root_16s_cor_df_padj_all.txt"),sep="\t",quote=F)







## Make node attribute table
nodeattrib_root_16s <- data.frame(node = union(root_16s_cor_df_padj$from,root_16s_cor_df_padj$to))
nodeattrib_root_16s$indicgroup <- 0

for (i in as.character(nodeattrib_root_16s$node))
{
  if (i %in% indic_edge_16s_root == TRUE)
  {nodeattrib_root_16s[nodeattrib_root_16s$node==i,"indicgroup"] <- paste(colnames(root_r_values_16s)[which(root_r_values_16s[i,1:4]==1)],collapse = "_")}
  else
  { nodeattrib_root_16s[nodeattrib_root_16s$node==i,"indicgroup"]<- "NA"}
}

## Add node taxonomy information to attribute table
nodeattrib_root_16s <- cbind(nodeattrib_root_16s,tax_filter_16s[as.character(nodeattrib_root_16s$node),])
nodeattrib_root_16s
write.table(nodeattrib_root_16s,paste0(output,"nodeattrib_root_16s.txt"),sep="\t",quote=F)

## Create co-occurrence network with igraph
root_net_16s <- graph_from_data_frame(root_16s_cor_df_padj,direct=F, vertices = nodeattrib_root_16s)
root_net_16s

#graphml 格式，可使用 gephi 软件打开并进行可视化编辑
write.graph(root_net_16s, 'all_otu1.graphml', format = 'graphml')

#gml 格式，可使用 cytoscape 软件打开并进行可视化编辑
write.graph(root_net_16s, 'all_otu1.gml', format = 'gml')




## Calculate relative abundance of OTU nodes
root_ra_16s <- apply(otu_norm_16s,1,mean)
root_ra_16s <- root_ra_16s[V(root_net_16s)$name]

write.table(root_ra_16s,paste0(output,"OTU_nodes_relative.txt"),sep="\t",quote=F)

## Network properties ##

## Number of nodes in network
length(V(root_net_16s))

## Number of edges in network
length(E(root_net_16s))

## Calculate network node degrees/max/min
root_16s_deg <- sort(degree(root_net_16s,mode="all"),decr=T)
max(root_16s_deg)
mean(root_16s_deg)

## Number of cropping system sensitive OTUs in the network
length(V(root_net_16s)$name[V(root_net_16s)$name %in% indic_edge_16s_root])

## Degree of cropping system sensitive OTUs
root_16s_deg[names(root_16s_deg) %in% indic_edge_16s_root]

## Relative abundance of cropping system sensitive OTUs
root_ra_16s[names(root_ra_16s) %in% indic_edge_16s_root]

## Taxonomy of cropping system sensitive OTUs
tax_filter_16s[names(sort(root_16s_deg[indic_edge_16s_root],decr=T)),]

## Set node colors based upon sensitivity to management system
cs <- c("Ggt","S_135","SA","H","Ggt_S_135","Ggt_SA","Ggt_H","S_135_SA","S_135_H","SA_H","Ggt_S_135_SA","Ggt_S_135_H","S_135_SA_H")
unique(V(root_net_16s)$indicgroup)
V(root_net_16s)$color <- V(root_net_16s)$indicgroup
V(root_net_16s)$color[!V(root_net_16s)$color %in% cs] <- "gray20"
V(root_net_16s)$color[V(root_net_16s)$color == "Ggt"] <- "dodgerblue4"
V(root_net_16s)$color[V(root_net_16s)$color == "S_135"] <- "orange"
V(root_net_16s)$color[V(root_net_16s)$color == "SA"] <- "green"
V(root_net_16s)$color[V(root_net_16s)$color == "H"] <- "violetred4"
V(root_net_16s)$color[V(root_net_16s)$color == "Ggt_S_135"] <- "dodgerblue1"
V(root_net_16s)$color[V(root_net_16s)$color == "Ggt_SA"] <- "firebrick4"
V(root_net_16s)$color[V(root_net_16s)$color == "Ggt_H"] <- "lightgoldenrod4"
V(root_net_16s)$color[V(root_net_16s)$color == "S_135_SA"] <- "tan4"
V(root_net_16s)$color[V(root_net_16s)$color == "S_135_H"] <- "firebrick2"
V(root_net_16s)$color[V(root_net_16s)$color == "SA_H"] <- "chartreuse"
V(root_net_16s)$color[V(root_net_16s)$color == "Ggt_S_135_SA"] <- "seashell4"
V(root_net_16s)$color[V(root_net_16s)$color == "Ggt_S_135_H"] <- "sienna4"
V(root_net_16s)$color[V(root_net_16s)$color == "S_135_SA_H"] <- "orchid"

V(root_net_16s)$frame.color <- V(root_net_16s)$color
V(root_net_16s)$frame.color
## Set vector of nodes responding to certain management systems
root16s_nodes <- rownames(nodeattrib_root_16s[nodeattrib_root_16s$indicgroup %in% cs,])

## Set node shape
V(root_net_16s)$shape <- "circle"

## Set node size
V(root_net_16s)$size <- V(root_net_16s)$name
V(root_net_16s)$size[!V(root_net_16s)$size %in% root16s_nodes] <- 2
V(root_net_16s)$size[V(root_net_16s)$size %in% root16s_nodes] <- 4
root16s_nodesizes <- as.numeric(V(root_net_16s)$size)
root16s_nodesizes 



#graphml 格式，可使用 gephi 软件打开并进行可视化编辑
write.graph(root_net_16s, 'otu2.graphml', format = 'graphml')

#gml 格式，可使用 cytoscape 软件打开并进行可视化编辑
write.graph(root_net_16s, 'otu2.gml', format = 'gml')



##### Explore community structure of root meta-network by defining modules #####

## Make vectors of network nodes responding to different cropping systems
Ggt_nodes_root <- rownames(nodeattrib_root_16s[nodeattrib_root_16s$indicgroup=="Ggt",])
Ggt_H_nodes_root <- rownames(nodeattrib_root_16s[nodeattrib_root_16s$indicgroup=="Ggt_H",])
Ggt_S_135_nodes_root <- rownames(nodeattrib_root_16s[nodeattrib_root_16s$indicgroup=="Ggt_S_135",])
Ggt_S_135_H_nodes_root <- rownames(nodeattrib_root_16s[nodeattrib_root_16s$indicgroup=="Ggt_S_135_H",])
Ggt_S_135_SA_nodes_root <- rownames(nodeattrib_root_16s[nodeattrib_root_16s$indicgroup=="Ggt_S_135_SA",])
Ggt_SA_nodes_root <- rownames(nodeattrib_root_16s[nodeattrib_root_16s$indicgroup=="Ggt_SA",])
H_nodes_root <- rownames(nodeattrib_root_16s[nodeattrib_root_16s$indicgroup=="H",])
S_135_H_nodes_root <- rownames(nodeattrib_root_16s[nodeattrib_root_16s$indicgroup=="S_135_H",])
S_135_SA_nodes_root <- rownames(nodeattrib_root_16s[nodeattrib_root_16s$indicgroup=="S_135_SA",])
S_135_SA_H_nodes_root <- rownames(nodeattrib_root_16s[nodeattrib_root_16s$indicgroup=="S_135_SA_H",])
SA_nodes_root <- rownames(nodeattrib_root_16s[nodeattrib_root_16s$indicgroup=="SA",])
SA_H_nodes_root <- rownames(nodeattrib_root_16s[nodeattrib_root_16s$indicgroup=="SA_H",])
S_135_nodes_root <- rownames(nodeattrib_root_16s[nodeattrib_root_16s$indicgroup=="S_135",])


cs_nodes_root_all <- c(Ggt_nodes_root, Ggt_H_nodes_root, Ggt_S_135_nodes_root, Ggt_S_135_H_nodes_root, Ggt_S_135_SA_nodes_root, Ggt_SA_nodes_root, H_nodes_root, S_135_H_nodes_root, S_135_SA_nodes_root, S_135_SA_H_nodes_root, SA_nodes_root, SA_H_nodes_root,S_135_nodes_root)
cs_nodes_root_all <- cs_nodes_root_all[grep("OTU*", cs_nodes_root_all)]
cs_nodes_root_all
## Perform cluster analysis using greedy clustering algorithm 
cfg_root <- cluster_fast_greedy(as.undirected(root_net_16s))
cfg_root

## Subset for 20 biggest nodes
root_modules <- sort(table(membership(cfg_root)),decr=T)
root_modules
root_modules_11 <- root_modules[1:11]
sum(root_modules_11)/sum(root_modules)
rm11_plot <- root_modules_11
names(rm11_plot) <- as.factor(1:11)
root_modules_cs <- table(factor(membership(cfg_root)[cs_nodes_root_all],levels=names(root_modules)))
root_modules_cs_11 <- root_modules_cs[names(root_modules_11)]
rmcs11_plot <- root_modules_cs_11
names(rmcs11_plot) <- as.factor(1:11)
rmcs11_plot
## Make vector of nodes in top 20 modules
root_modules_points <- membership(cfg_root)[membership(cfg_root) %in% names(root_modules_11)]
root_modules2 <- as.data.frame(root_modules_points)
root_modules2
write.table(root_modules2,paste0(output,"root_modules2.txt"),sep="\t",quote=F)


root_points <- NULL

for(i in root_modules_points){
  rootx <- which(names(root_modules_11)==i)
  root_points <- c(root_points,rootx)
}

names(root_points) <- names(root_modules_points)

## Set node colors by cropping system sensitivity 
root_all_cols <- sort(root_points)
root_all_cols[!names(root_all_cols) %in% cs] <- "gray20"
root_all_cols[names(root_all_cols) %in% Ggt_nodes_root] <- "dodgerblue4"
root_all_cols[names(root_all_cols) %in% S_135_nodes_root] <- "orange"
root_all_cols[names(root_all_cols) %in% SA_nodes_root] <- "green"
root_all_cols[names(root_all_cols) %in% H_nodes_root] <- "violetred4"
root_all_cols[names(root_all_cols) %in% Ggt_S_135_nodes_root] <- "dodgerblue1"
root_all_cols[names(root_all_cols) %in% Ggt_SA_nodes_root] <- "firebrick4"
root_all_cols[names(root_all_cols) %in% Ggt_H_nodes_root] <- "lightgoldenrod4"
root_all_cols[names(root_all_cols) %in% S_135_SA_nodes_root] <- "tan4"
root_all_cols[names(root_all_cols) %in% S_135_H_nodes_root] <- "firebrick2"
root_all_cols[names(root_all_cols) %in% SA_H_nodes_root] <- "chartreuse"
root_all_cols[names(root_all_cols) %in% Ggt_S_135_SA_nodes_root] <- "seashell4"
root_all_cols[names(root_all_cols) %in% Ggt_S_135_H_nodes_root] <- "sienna4"
root_all_cols[names(root_all_cols) %in% S_135_SA_H_nodes_root] <- "orchid"


root_all_pch <- sort(root_points)
root_all_pch[names(root_all_pch) %in% rownames(otu_norm_16s)] <- 1
root_all_pch[names(root_all_pch) %in% intersect(rownames(otu_norm_16s),cs_nodes_root_all)] <- 16
#root_all_pch[names(root_all_pch) %in% names(root_all_keystone)] <- 8

root_all_cex <- sort(root_points)
root_all_cex[!names(root_all_cex) %in% cs_nodes_root_all] <- 1
root_all_cex[names(root_all_cex) %in% cs_nodes_root_all] <- 2

root_mods_list_cs <- list()

for (i in names(root_modules_cs_11)){
  x1 <- names(membership(cfg_root)[membership(cfg_root) == i])
  x2 <- x1[x1 %in% cs_nodes_root_all]
  root_mods_list_cs[[i]] <- as.numeric(V(root_net_16s)[x2])
}


## Note: the permuations for the layout of the network can be very time consuming and processor intensive
set.seed(8051)
coords_root_16s <- layout_(root_net_16s, with_fr(niter=9999, grid="nogrid"))
write.table(coords_root_16s,paste0(output,"coords_root_16s2.txt"),sep="\t",row.names=F,col.names=F,quote=F)

## Import pre-calculated FR layout coordinates to save time 
coords_soil_16s <- as.matrix(read.table("coords_soil_16s2.txt"))
dimnames(coords_root_16s) <-  NULL

#graphml 格式，可使用 gephi 软件打开并进行可视化编辑
write.graph(root_net_16s, 'all_otu3.graphml', format = 'graphml')

#gml 格式，可使用 cytoscape 软件打开并进行可视化编辑
write.graph(root_net_16s, 'all_otu3.gml', format = 'gml')





## Plot Figure 4
pdf(paste0(output,"all-network.pdf"),width=7,height=3.5)
par(mfrow=c(1,1), mar=c(0,0,0,0))


root_cols <- c("darkseagreen4", "darkseagreen3", "darkseagreen1")
plot(root_net_16s,vertex.label=NA,vertex.size=root16s_nodesizes,layout=coords_root_16s,
     mark.groups=list(root_mods_list_cs$`2`,root_mods_list_cs$`1`,root_mods_list_cs$`3`),
     mark.col=root_cols, mark.border=root_cols)
legend("bottomright",legend=c("Module 2", "Module 1", "Module 3"), col=root_cols,
       bty="n", fill=root_cols, border=root_cols)

cols <- c("dodgerblue4","orange", "green", "violetred4", "dodgerblue1","firebrick4", "lightgoldenrod4", "tan4", "firebrick2","chartreuse", "seashell4", "sienna4","orchid" )
names(cols) <- c("Ggt", "S_135", "SA", "H", "Ggt_S_135", "Ggt_SA", "Ggt_H", "S_135_SA", "S_135_H", "SA_H","Ggt_S_135_SA","Ggt_S_135_H", "S_135_SA_H")

legend2 <- as.data.frame(cols)
legend2$rownames <- rownames(legend2)
legend2
legend("right",legend=legend2$rownames, col=legend2$cols,
       bty="n", fill=legend2$cols, border=legend2$cols)


dev.off()




####子网络提取
####方法一：直接从网络中提取
## Create co-occurrence network with igraph
root_16s_cor_df_padj<-read.delim('root_16s_cor_df_padj_all.txt', row.names = 1, check.names = FALSE)
root_16s_cor_df_padj

nodeattrib_root_16s<-read.delim('nodeattrib_root_16s.txt', row.names = 1, check.names = FALSE)
nodeattrib_root_16s


root_net_16s <- graph_from_data_frame(root_16s_cor_df_padj,direct=F, vertices = nodeattrib_root_16s)
root_net_16s

#graphml 格式，可使用 gephi 软件打开并进行可视化编辑
#write.graph(root_net_16s, 'all_otu1.graphml', format = 'graphml')

#gml 格式，可使用 cytoscape 软件打开并进行可视化编辑
#write.graph(root_net_16s, 'all_otu1.gml', format = 'gml')


#读取每个样品的代表 OTU
sensitiv_otu <- read.delim('sensitiv-otu2.txt', row.names = 1, check.names = FALSE)
sensitiv_otu
#根据样本代表节点，从网络中划分与各样本有关的“子网络”，详情 ?subgraph
#并将所有子网络存储到一个列表（list）里面
sub_graph <- list()
for (i in names(sensitiv_otu)) {
  sample_i <- sensitiv_otu[i]
  select_node <- rownames(sample_i)[which(sample_i != 0)]
  sub_graph[[i]] <- subgraph(root_net_16s, select_node)
}
sub_graph

##简单画图展示全局网络或各子网络的结构，具体的可视化细节如有需要请自行摸索
#library(phyloseq)

#plot_network(g)  #全局网络
#plot_network(sub_graph$Ggt)  #以“sample1”的代表节点划分的“子网络”
#plot_network(sub_graph$S_135)  #以“sample2”的代表节点划分的“子网络”
#plot_network(sub_graph$SA)  #以“sample18”的代表节点划分的“子网络”
#plot_network(sub_graph$H)  #以“sample18”的代表节点划分的“子网络

##Ggt处理子网络(采用gephi软件进行边列表和节点列表的提取)
sub_graph$Ggt

#边列表
edge <- data.frame(as_edgelist(sub_graph$Ggt))    #igraph 的邻接列表转为边列表

edge_list <- data.frame(
  source = edge[[1]],
  target = edge[[2]],
  weight = E(sub_graph$Ggt)$weight,
  cor = E(sub_graph$Ggt)$cor,
  corr = E(sub_graph$Ggt)$corr
)
head(edge_list)
write.table(edge_list, 'Ggt_network.edge_list.txt', sep = '\t', row.names = FALSE, quote = FALSE)


#节点属性列表
#node_list <- data.frame(
#  name = names(V(sub_graph$Ggt)),
#  indicgroup = V(sub_graph$Ggt)$indicgroup,
#  phylum = V(sub_graph$Ggt)$phylum,
# class = V(sub_graph$Ggt)$class,
#  order = V(sub_graph$Ggt)$order,
#  family = V(sub_graph$Ggt)$family,
#  genus = V(sub_graph$Ggt)$genus,
#  species = V(sub_graph$Ggt)$species, 
#  labels = V(sub_graph$Ggt)$labels, 
#  cols = V(sub_graph$Ggt)$cols, 
#  padj = V(sub_graph$Ggt)$padj
#)
#head(node_list)


#graphml 格式，可使用 gephi 软件打开并进行可视化编辑
write.graph(sub_graph$Ggt, 'sub_graph_Ggt.graphml', format = 'graphml')

#gml 格式，可使用 cytoscape 软件打开并进行可视化编辑
write.graph(sub_graph$Ggt, 'sub_graph_Ggt.gml', format = 'gml')


##S_135处理子网络(采用gephi软件进行边列表和节点列表的提取)
sub_graph$S_135

#graphml 格式，可使用 gephi 软件打开并进行可视化编辑
write.graph(sub_graph$S_135, 'sub_graph_S_135.graphml', format = 'graphml')

#gml 格式，可使用 cytoscape 软件打开并进行可视化编辑
write.graph(sub_graph$S_135, 'sub_graph_S_135.gml', format = 'gml')


##SA处理子网络(采用gephi软件进行边列表和节点列表的提取)
sub_graph$SA

#graphml 格式，可使用 gephi 软件打开并进行可视化编辑
write.graph(sub_graph$SA, 'sub_graph_SA.graphml', format = 'graphml')

#gml 格式，可使用 cytoscape 软件打开并进行可视化编辑
write.graph(sub_graph$SA, 'sub_graph_SA.gml', format = 'gml')

##H处理子网络(采用gephi软件进行边列表和节点列表的提取)
sub_graph$H

#graphml 格式，可使用 gephi 软件打开并进行可视化编辑
write.graph(sub_graph$H, 'sub_graph_H.graphml', format = 'graphml')

#gml 格式，可使用 cytoscape 软件打开并进行可视化编辑
write.graph(sub_graph$H, 'sub_graph_H.gml', format = 'gml')



##方法二：利用igraph包对网络文件与邻接矩阵的相互转换功能，可以先生成邻接矩阵列表，再作图

##网络文件输出，输出特定的网络文件类型，便于后续数据分析需求
#邻接矩阵，出了上述提到的在计算相关系数后，输出筛选后的相关系数矩阵外
#还可以由 igraph 的邻接列表转换
#adj_matrix <- as.matrix(get.adjacency(root_net_16s, attr = 'cor'))
#write.table(data.frame(adj_matrix, check.names = FALSE), 'all_network.adj_matrix.txt', col.names = NA, sep = '\t', quote = FALSE)


####核心物种的识别——zi_pi####
####模块内连通度（Zi）和模块间连通度（Pi）识别关键节点####
##计算模块内连通度（Zi）和模块间连通度（Pi）
source('zi_pi.r')

#第一步:从网络文件中提取出邻接矩阵
#注意：root_net_16s为所构建的网路名称，cor是最初构建网络名称时对相关性系数的命名，需要核查原始代码或者边列表来看，否则报错
#大网络邻接矩阵提取meta-network
adj_matrix <- as.matrix(get.adjacency(root_net_16s, attr = 'cor'))
write.table(data.frame(adj_matrix, check.names = FALSE), 'all_network.adj_matrix.txt', col.names = NA, sep = '\t', quote = FALSE)

#Ggt子网络邻接矩阵提取Ggt_network
Ggt_adj_matrix <- as.matrix(get.adjacency(sub_graph$Ggt, attr = 'cor'))
write.table(data.frame(Ggt_adj_matrix, check.names = FALSE), 'Ggt_network.adj_matrix.txt', col.names = NA, sep = '\t', quote = FALSE)


#S_135子网络邻接矩阵提取S_135_network
S_135_adj_matrix <- as.matrix(get.adjacency(sub_graph$S_135, attr = 'cor'))
write.table(data.frame(S_135_adj_matrix, check.names = FALSE), 'S_135_network.adj_matrix.txt', col.names = NA, sep = '\t', quote = FALSE)


#SA子网络邻接矩阵提取SA_network
SA_adj_matrix <- as.matrix(get.adjacency(sub_graph$SA, attr = 'cor'))
write.table(data.frame(SA_adj_matrix, check.names = FALSE), 'SA_network.adj_matrix.txt', col.names = NA, sep = '\t', quote = FALSE)

#H子网络邻接矩阵提取H_network
H_adj_matrix <- as.matrix(get.adjacency(sub_graph$H, attr = 'cor'))
write.table(data.frame(H_adj_matrix, check.names = FALSE), 'H_network.adj_matrix.txt', col.names = NA, sep = '\t', quote = FALSE)



#第二步：将邻接矩阵转换为0-1矩阵
#获取0-1矩阵，1表示节点间存在边，即存在互作，0表示不存在边，即无互作
#大网络邻接矩阵转换为0-1矩阵meta-network
adjacency_matrix<-as.matrix(adj_matrix)
adjacency_matrix[abs(adjacency_matrix)!=0]<-1

write.table(adjacency_matrix,"all_adjacency_unweight.txt",sep = '\t',quote = FALSE,col.names = NA)


#Ggt网络邻接矩阵转换为0-1矩阵
Ggt_adjacency_matrix<-as.matrix(Ggt_adj_matrix)
Ggt_adjacency_matrix[abs(Ggt_adjacency_matrix)!=0]<-1

write.table(Ggt_adjacency_matrix,"Ggt_adjacency_unweight.txt",sep = '\t',quote = FALSE,col.names = NA)


#S_135网络邻接矩阵转换为0-1矩阵
S_135_adjacency_matrix<-as.matrix(S_135_adj_matrix)
S_135_adjacency_matrix[abs(S_135_adjacency_matrix)!=0]<-1

write.table(S_135_adjacency_matrix,"S_135_adjacency_unweight.txt",sep = '\t',quote = FALSE,col.names = NA)


#SA网络邻接矩阵转换为0-1矩阵
SA_adjacency_matrix<-as.matrix(SA_adj_matrix)
SA_adjacency_matrix[abs(SA_adjacency_matrix)!=0]<-1

write.table(SA_adjacency_matrix,"SA_adjacency_unweight.txt",sep = '\t',quote = FALSE,col.names = NA)

#H网络邻接矩阵转换为0-1矩阵
H_adjacency_matrix<-as.matrix(H_adj_matrix)
H_adjacency_matrix[abs(H_adjacency_matrix)!=0]<-1

write.table(H_adjacency_matrix,"H_adjacency_unweight.txt",sep = '\t',quote = FALSE,col.names = NA)

#第三步邻接矩阵与节点矩阵的对齐，为计算zi_pi做准备

####大网络meta####

#读取上述的邻接矩阵类型的网络文件
all_adjacency_unweight <- read.delim('all_adjacency_unweight.txt', row.names = 1, sep = '\t', check.names = FALSE)
head(all_adjacency_unweight)
#读取节点属性列表，包含节点所划分的模块
all_nodes_list <- read.delim('all_network_node_list-4.txt', row.names = 1, sep = '\t', check.names = FALSE)
rownames(all_nodes_list)<-all_nodes_list$v_name
all_nodes_list
#两个文件的节点顺序要一致
all_nodes_list <- all_nodes_list[rownames(all_adjacency_unweight), ]

#计算模块内连通度（Zi）和模块间连通度（Pi）
#指定邻接矩阵、节点列表、节点列表中节点度和模块度的列名称
all_zi_pi <- zi.pi(all_nodes_list, all_adjacency_unweight, degree = 'degree', modularity_class = 'modularity_class')
head(all_zi_pi)

##可再根据阈值对节点划分为 4 种类型，并作图展示其分布
library(ggplot2)

#all_zi_pi <- na.omit(all_zi_pi)   #NA 值最好去掉，不要当 0 处理
all_zi_pi[which(all_zi_pi$within_module_connectivities < 2.5 & all_zi_pi$among_module_connectivities < 0.62),'type'] <- 'Peripherals'
all_zi_pi[which(all_zi_pi$within_module_connectivities < 2.5 & all_zi_pi$among_module_connectivities > 0.62),'type'] <- 'Connectors'
all_zi_pi[which(all_zi_pi$within_module_connectivities > 2.5 & all_zi_pi$among_module_connectivities < 0.62),'type'] <- 'Module hubs'
all_zi_pi[which(all_zi_pi$within_module_connectivities > 2.5 & all_zi_pi$among_module_connectivities > 0.62),'type'] <- 'Network hubs'

write.table(all_zi_pi, 'all_zi_pi_result.txt', sep = '\t', row.names = FALSE, quote = FALSE)


all_zi_pi_p <-ggplot(all_zi_pi, aes(among_module_connectivities, within_module_connectivities)) +
  geom_point(aes(color = type), alpha = 0.5, size = 2) +
  scale_color_manual(values = c('gray','red','blue','purple'), 
                     limits = c('Peripherals', 'Connectors', 'Module hubs', 'Network hubs'))+
  theme(panel.grid = element_blank(), axis.line = element_line(colour = 'black'), 
        panel.background = element_blank(), legend.key = element_blank()) +
  labs(x = 'Among-module connectivities', y = 'Within-module connectivities', color = '') +
  geom_vline(xintercept = 0.62) +
  geom_hline(yintercept = 2.5)

all_zi_pi_p

ggsave("./all_zi_pi_p.pdf", all_zi_pi_p, width = 150, height = 100, units = "mm")
ggsave("./all_zi_pi_p.png", all_zi_pi_p, width = 150, height = 100, units = "mm")
ggsave("./all_zi_pi_p.tiff", all_zi_pi_p, width = 150, height = 100, units = "mm")

ggsave("./all_zi_pi_p2.pdf", all_zi_pi_p, width = 250, height = 200, units = "mm")
ggsave("./all_zi_pi_p2.png", all_zi_pi_p, width = 250, height = 200, units = "mm")
ggsave("./all_zi_pi_p2.tiff", all_zi_pi_p, width = 250, height = 200, units = "mm")




####Ggt网络####

#读取上述的邻接矩阵类型的网络文件
Ggt_adjacency_unweight <- read.delim('Ggt_adjacency_unweight.txt', row.names = 1, sep = '\t', check.names = FALSE)
head(Ggt_adjacency_unweight)
#读取节点属性列表，包含节点所划分的模块
Ggt_nodes_list <- read.delim('Ggt_network_node_list.txt', row.names = 1, sep = '\t', check.names = FALSE)
rownames(Ggt_nodes_list)<-Ggt_nodes_list$v_name
Ggt_nodes_list
#两个文件的节点顺序要一致
Ggt_nodes_list <- Ggt_nodes_list[rownames(Ggt_adjacency_unweight), ]

#计算模块内连通度（Zi）和模块间连通度（Pi）
#指定邻接矩阵、节点列表、节点列表中节点度和模块度的列名称
Ggt_zi_pi <- zi.pi(Ggt_nodes_list, Ggt_adjacency_unweight, degree = 'degree', modularity_class = 'modularity_class')
head(Ggt_zi_pi)

##可再根据阈值对节点划分为 4 种类型，并作图展示其分布
library(ggplot2)

#all_zi_pi <- na.omit(all_zi_pi)   #NA 值最好去掉，不要当 0 处理
Ggt_zi_pi[which(Ggt_zi_pi$within_module_connectivities < 2.5 & Ggt_zi_pi$among_module_connectivities < 0.62),'type'] <- 'Peripherals'
Ggt_zi_pi[which(Ggt_zi_pi$within_module_connectivities < 2.5 & Ggt_zi_pi$among_module_connectivities > 0.62),'type'] <- 'Connectors'
Ggt_zi_pi[which(Ggt_zi_pi$within_module_connectivities > 2.5 & Ggt_zi_pi$among_module_connectivities < 0.62),'type'] <- 'Module hubs'
Ggt_zi_pi[which(Ggt_zi_pi$within_module_connectivities > 2.5 & Ggt_zi_pi$among_module_connectivities > 0.62),'type'] <- 'Network hubs'

write.table(Ggt_zi_pi, 'Ggt_zi_pi_result.txt', sep = '\t', row.names = FALSE, quote = FALSE)


Ggt_zi_pi_p <-ggplot(Ggt_zi_pi, aes(among_module_connectivities, within_module_connectivities)) +
  geom_point(aes(color = type), alpha = 0.5, size = 2) +
  scale_color_manual(values = c('gray','red','blue','purple'), 
                     limits = c('Peripherals', 'Connectors', 'Module hubs', 'Network hubs'))+
  theme(panel.grid = element_blank(), axis.line = element_line(colour = 'black'), 
        panel.background = element_blank(), legend.key = element_blank()) +
  labs(x = 'Among-module connectivities', y = 'Within-module connectivities', color = '') +
  geom_vline(xintercept = 0.62) +
  geom_hline(yintercept = 2.5)

Ggt_zi_pi_p

ggsave("./Ggt_zi_pi_p.pdf", Ggt_zi_pi_p, width = 150, height = 100, units = "mm")
ggsave("./Ggt_zi_pi_p.png", Ggt_zi_pi_p, width = 150, height = 100, units = "mm")
ggsave("./Ggt_zi_pi_p.tiff", Ggt_zi_pi_p, width = 150, height = 100, units = "mm")

ggsave("./Ggt_zi_pi_p2.pdf", Ggt_zi_pi_p, width = 250, height = 200, units = "mm")
ggsave("./Ggt_zi_pi_p2.png", Ggt_zi_pi_p, width = 250, height = 200, units = "mm")
ggsave("./Ggt_zi_pi_p2.tiff", Ggt_zi_pi_p, width = 250, height = 200, units = "mm")


####S_135网络####

#读取上述的邻接矩阵类型的网络文件
S_135_adjacency_unweight <- read.delim('S_135_adjacency_unweight.txt', row.names = 1, sep = '\t', check.names = FALSE)
head(S_135_adjacency_unweight)
#读取节点属性列表，包含节点所划分的模块
S_135_nodes_list <- read.delim('135_network_node_list.txt', row.names = 1, sep = '\t', check.names = FALSE)
rownames(S_135_nodes_list)<-S_135_nodes_list$v_name
S_135_nodes_list
#两个文件的节点顺序要一致
S_135_nodes_list <- S_135_nodes_list[rownames(S_135_adjacency_unweight), ]

#计算模块内连通度（Zi）和模块间连通度（Pi）
#指定邻接矩阵、节点列表、节点列表中节点度和模块度的列名称
S_135_zi_pi <- zi.pi(S_135_nodes_list, S_135_adjacency_unweight, degree = 'degree', modularity_class = 'modularity_class')
head(S_135_zi_pi)

##可再根据阈值对节点划分为 4 种类型，并作图展示其分布
library(ggplot2)

#all_zi_pi <- na.omit(all_zi_pi)   #NA 值最好去掉，不要当 0 处理
S_135_zi_pi[which(S_135_zi_pi$within_module_connectivities < 2.5 & S_135_zi_pi$among_module_connectivities < 0.62),'type'] <- 'Peripherals'
S_135_zi_pi[which(S_135_zi_pi$within_module_connectivities < 2.5 & S_135_zi_pi$among_module_connectivities > 0.62),'type'] <- 'Connectors'
S_135_zi_pi[which(S_135_zi_pi$within_module_connectivities > 2.5 & S_135_zi_pi$among_module_connectivities < 0.62),'type'] <- 'Module hubs'
S_135_zi_pi[which(S_135_zi_pi$within_module_connectivities > 2.5 & S_135_zi_pi$among_module_connectivities > 0.62),'type'] <- 'Network hubs'

write.table(S_135_zi_pi, 'S_135_zi_pi_result.txt', sep = '\t', row.names = FALSE, quote = FALSE)


S_135_zi_pi_p <-ggplot(S_135_zi_pi, aes(among_module_connectivities, within_module_connectivities)) +
  geom_point(aes(color = type), alpha = 0.5, size = 2) +
  scale_color_manual(values = c('gray','red','blue','purple'), 
                     limits = c('Peripherals', 'Connectors', 'Module hubs', 'Network hubs'))+
  theme(panel.grid = element_blank(), axis.line = element_line(colour = 'black'), 
        panel.background = element_blank(), legend.key = element_blank()) +
  labs(x = 'Among-module connectivities', y = 'Within-module connectivities', color = '') +
  geom_vline(xintercept = 0.62) +
  geom_hline(yintercept = 2.5)

S_135_zi_pi_p

ggsave("./S_135_zi_pi_p.pdf", S_135_zi_pi_p, width = 150, height = 100, units = "mm")
ggsave("./S_135_zi_pi_p.png", S_135_zi_pi_p, width = 150, height = 100, units = "mm")
ggsave("./S_135_zi_pi_p.tiff", S_135_zi_pi_p, width = 150, height = 100, units = "mm")

ggsave("./S_135_zi_pi_p2.pdf", S_135_zi_pi_p, width = 250, height = 200, units = "mm")
ggsave("./S_135_zi_pi_p2.png", S_135_zi_pi_p, width = 250, height = 200, units = "mm")
ggsave("./S_135_zi_pi_p2.tiff", S_135_zi_pi_p, width = 250, height = 200, units = "mm")


####SA网络####

#读取上述的邻接矩阵类型的网络文件
SA_adjacency_unweight <- read.delim('SA_adjacency_unweight.txt', row.names = 1, sep = '\t', check.names = FALSE)
head(SA_adjacency_unweight)
#读取节点属性列表，包含节点所划分的模块
SA_nodes_list <- read.delim('SA_network_node_list.txt', row.names = 1, sep = '\t', check.names = FALSE)
rownames(SA_nodes_list)<-SA_nodes_list$v_name
SA_nodes_list
#两个文件的节点顺序要一致
SA_nodes_list <- SA_nodes_list[rownames(SA_adjacency_unweight), ]

#计算模块内连通度（Zi）和模块间连通度（Pi）
#指定邻接矩阵、节点列表、节点列表中节点度和模块度的列名称
SA_zi_pi <- zi.pi(SA_nodes_list, SA_adjacency_unweight, degree = 'degree', modularity_class = 'modularity_class')
head(SA_zi_pi)

##可再根据阈值对节点划分为 4 种类型，并作图展示其分布
library(ggplot2)

#all_zi_pi <- na.omit(all_zi_pi)   #NA 值最好去掉，不要当 0 处理
SA_zi_pi[which(SA_zi_pi$within_module_connectivities < 2.5 & SA_zi_pi$among_module_connectivities < 0.62),'type'] <- 'Peripherals'
SA_zi_pi[which(SA_zi_pi$within_module_connectivities < 2.5 & SA_zi_pi$among_module_connectivities > 0.62),'type'] <- 'Connectors'
SA_zi_pi[which(SA_zi_pi$within_module_connectivities > 2.5 & SA_zi_pi$among_module_connectivities < 0.62),'type'] <- 'Module hubs'
SA_zi_pi[which(SA_zi_pi$within_module_connectivities > 2.5 & SA_zi_pi$among_module_connectivities > 0.62),'type'] <- 'Network hubs'

write.table(SA_zi_pi, 'SA_zi_pi_result.txt', sep = '\t', row.names = FALSE, quote = FALSE)


SA_zi_pi_p <-ggplot(SA_zi_pi, aes(among_module_connectivities, within_module_connectivities)) +
  geom_point(aes(color = type), alpha = 0.5, size = 2) +
  scale_color_manual(values = c('gray','red','blue','purple'), 
                     limits = c('Peripherals', 'Connectors', 'Module hubs', 'Network hubs'))+
  theme(panel.grid = element_blank(), axis.line = element_line(colour = 'black'), 
        panel.background = element_blank(), legend.key = element_blank()) +
  labs(x = 'Among-module connectivities', y = 'Within-module connectivities', color = '') +
  geom_vline(xintercept = 0.62) +
  geom_hline(yintercept = 2.5)

SA_zi_pi_p

ggsave("./SA_zi_pi_p.pdf", SA_zi_pi_p, width = 150, height = 100, units = "mm")
ggsave("./SA_zi_pi_p.png", SA_zi_pi_p, width = 150, height = 100, units = "mm")
ggsave("./SA_zi_pi_p.tiff", SA_zi_pi_p, width = 150, height = 100, units = "mm")

ggsave("./SA_zi_pi_p2.pdf", SA_zi_pi_p, width = 250, height = 200, units = "mm")
ggsave("./SA_zi_pi_p2.png", SA_zi_pi_p, width = 250, height = 200, units = "mm")
ggsave("./SA_zi_pi_p2.tiff", SA_zi_pi_p, width = 250, height = 200, units = "mm")





####H网络####

#读取上述的邻接矩阵类型的网络文件
H_adjacency_unweight <- read.delim('H_adjacency_unweight.txt', row.names = 1, sep = '\t', check.names = FALSE)
head(H_adjacency_unweight)
#读取节点属性列表，包含节点所划分的模块
H_nodes_list <- read.delim('H_network_node_list.txt', row.names = 1, sep = '\t', check.names = FALSE)
rownames(H_nodes_list)<-H_nodes_list$v_name
H_nodes_list
#两个文件的节点顺序要一致
H_nodes_list <- H_nodes_list[rownames(H_adjacency_unweight), ]

#计算模块内连通度（Zi）和模块间连通度（Pi）
#指定邻接矩阵、节点列表、节点列表中节点度和模块度的列名称
H_zi_pi <- zi.pi(H_nodes_list, H_adjacency_unweight, degree = 'degree', modularity_class = 'modularity_class')
head(H_zi_pi)

##可再根据阈值对节点划分为 4 种类型，并作图展示其分布
library(ggplot2)

#all_zi_pi <- na.omit(all_zi_pi)   #NA 值最好去掉，不要当 0 处理
H_zi_pi[which(H_zi_pi$within_module_connectivities < 2.5 & H_zi_pi$among_module_connectivities < 0.62),'type'] <- 'Peripherals'
H_zi_pi[which(H_zi_pi$within_module_connectivities < 2.5 & H_zi_pi$among_module_connectivities > 0.62),'type'] <- 'Connectors'
H_zi_pi[which(H_zi_pi$within_module_connectivities > 2.5 & H_zi_pi$among_module_connectivities < 0.62),'type'] <- 'Module hubs'
H_zi_pi[which(H_zi_pi$within_module_connectivities > 2.5 & H_zi_pi$among_module_connectivities > 0.62),'type'] <- 'Network hubs'

write.table(H_zi_pi, 'H_zi_pi_result.txt', sep = '\t', row.names = FALSE, quote = FALSE)


H_zi_pi_p <-ggplot(H_zi_pi, aes(among_module_connectivities, within_module_connectivities)) +
  geom_point(aes(color = type), alpha = 0.5, size = 2) +
  scale_color_manual(values = c('gray','red','blue','purple'), 
                     limits = c('Peripherals', 'Connectors', 'Module hubs', 'Network hubs'))+
  theme(panel.grid = element_blank(), axis.line = element_line(colour = 'black'), 
        panel.background = element_blank(), legend.key = element_blank()) +
  labs(x = 'Among-module connectivities', y = 'Within-module connectivities', color = '') +
  geom_vline(xintercept = 0.62) +
  geom_hline(yintercept = 2.5)

H_zi_pi_p

ggsave("./H_zi_pi_p.pdf", H_zi_pi_p, width = 150, height = 100, units = "mm")
ggsave("./H_zi_pi_p.png", H_zi_pi_p, width = 150, height = 100, units = "mm")
ggsave("./H_zi_pi_p.tiff", H_zi_pi_p, width = 150, height = 100, units = "mm")

ggsave("./H_zi_pi_p2.pdf", H_zi_pi_p, width = 250, height = 200, units = "mm")
ggsave("./H_zi_pi_p2.png", H_zi_pi_p, width = 250, height = 200, units = "mm")
ggsave("./H_zi_pi_p2.tiff", H_zi_pi_p, width = 250, height = 200, units = "mm")



####核心物种的识别——网络参与系数（Pi）和模块参与度（PM）####
#定义函数计算 Pi
#kiS 为节点 i 与模块 S 中节点的边数
#ki 是节点 i 的总边数
Pi <- function(all_adjacency_unweight, all_nodes_list) {
  Pi <- NULL
  NM <- unique(all_nodes_list$modularity_class) #注意，此处节点列表中的模块名称为modularity_class,而非module
  for (i in rownames(all_nodes_list)) {
    ki <- all_nodes_list[i,'degree']
    kiS_ki = 0
    for (S in NM) {
      S <- subset(all_nodes_list, modularity_class == S)
      i_S <- all_adjacency_unweight[i,rownames(S)]
      i_S[i_S != 0] <- 1
      kiS <- sum(i_S)
      kiS_ki = kiS_ki + (kiS/ki)^2
    }
    Pi <- c(Pi, 1 - kiS_ki)
  }
  Pi <- as.data.frame(Pi)
  rownames(Pi) <- rownames(all_nodes_list)
  Pi
}

#定义函数计算 PM
#kim 是节点 i 与模块 m 中节点的边数
#ki 是节点 i 的总边数
PM <- function(all_adjacency_unweight, all_nodes_list) {
  PM <- NULL
  NM <- unique(all_nodes_list$modularity_class)
  for (i in rownames(all_nodes_list)) {
    ki <- all_nodes_list[i,'degree']
    PM_i <- c()
    for (S in NM) {
      S <- all_nodes_list[which(all_nodes_list$modularity_class == S), ]
      i_S <- all_adjacency_unweight[i,rownames(S)]
      i_S[i_S != 0] <- 1
      kim <- sum(i_S)
      PM_i <- c(PM_i, kim/ki)
    }
    PM <- rbind(PM, PM_i)
  }
  PM <- as.data.frame(PM)
  rownames(PM) <- rownames(all_nodes_list)
  colnames(PM) <- NM
  PM
}

#读取上述的邻接矩阵类型的网络文件
all_adjacency_unweight <- read.delim('all_adjacency_unweight.txt', row.names = 1, sep = '\t', check.names = FALSE)
head(all_adjacency_unweight)
#读取节点属性列表，包含节点所划分的模块
all_nodes_list <- read.delim('all_network_node_list-4.txt', row.names = 1, sep = '\t', check.names = FALSE)
rownames(all_nodes_list)<-all_nodes_list$v_name
all_nodes_list
#两个文件的节点顺序要一致
all_nodes_list <- all_nodes_list[rownames(all_adjacency_unweight), ]


#指定邻接矩阵，节点属性列表计算 Pi 和 PM
net_Pi <- Pi(all_adjacency_unweight, all_nodes_list)
head(net_Pi)  #各节点及其 Pi 值

net_PM <- PM(all_adjacency_unweight, all_nodes_list)
head(net_PM)  #各节点（行）在各模块（列）的 PM 值

#输出
write.table(net_Pi, 'net_Pi.txt', sep = '\t', col.names = NA, quote = FALSE)
write.table(net_PM, 'net_PM.txt', sep = '\t', col.names = NA, quote = FALSE)





####核心物种的识别——SPEC-OCCU plot####
#今天介绍一个用于识别潜在关键物种的方法，叫做特异性-占有率图（Specificity-Occupancy, SPEC-OCCU plot）
#Specificity为一个OTU在分组A中各样品的平均相对丰度与该OTU在所有研究分组中各分组平均相对丰度之和的比值。
#Occupancy为一个OTU在分组A中检出样本的数目与在所有研究样本中检出样本的数目的比值。
#在本文中，将Specificity和Occupancy均高于0.7的OTU归类为一组样本中的潜在关键物种。

#除了需要计算OTU的Specificity和Occupancy之外，还需要得到OTU在组内的平均丰度对点的大小进行赋值，还需要得到OTU对应的门水平分类对点的颜色进行赋值。


#行为OTU/ASV，列为样本，数值为OTU/ASV在样本中的相对丰度，在数据的最后，给出每一个OTU/ASV对应的分类学水平。
#制作此文件十分容易，将我们测序后注释得到的带分类学水平的OTU丰度表输入到R中之后，使用如下代码即可：
##读入带物种注释的原始OTU表格
#library(tidyverse)
#library(dplyr)
## 转化为相对丰度,注意-ncol(OTU)的目的是去除带注释的物种分度表格中物种注释部分，如果不带，直接删除
#一个数据表格中，如何删除某几列数据或者某几行数据
#otu[,-ncol(OTU)] <- t(t(otu[,-ncol(OTU)])/colSums(otu[,-ncol(OTU)])*100)
## 得到分类学水平
#bac3 <- otu %>% separate(taxonomy,c("Domain","Phylum","Class","Order","Family","Genus","Species"),"; ")

##读入OTU表格
otu_filter_16s <- read.delim("otu_filter_16s.txt", row.names=1,sep="\t", header=T, blank.lines.skip=F, check.names=F)
otu_filter_16s

otu_filter_16s_relative <- t(t(otu_filter_16s)/colSums(otu_filter_16s)*100)

##读入tax表格
tax_filter_16s <- read.delim("tax_filter_16s.txt", row.names=1,sep="\t", header=T, blank.lines.skip=F, check.names=F)
tax_filter_16s
tax_filter_16s <- tax_filter_16s[,-ncol(tax_filter_16s)]
tax_filter_16s

#将otu丰度表格与物种注释表格合并
otu_tax <- merge(x=otu_filter_16s_relative,y=tax_filter_16s,by='row.names')
head(otu_tax)
otu_tax =dplyr::rename(otu_tax,OTUID = Row.names)
head(otu_tax, n = 3)
write.table(otu_tax, 'otu_tax.txt', row.names = FALSE, sep = '\t', quote = FALSE)

##读入分组文件
group<- read.delim("mapping.txt", row.names=1,sep="\t", header=T, blank.lines.skip=F, check.names=F)
group


##计算各组样本中，OTU 的特异性（specificity）和占有率（occupancy）

Nindividuals_S <- rep(0, nrow(otu_tax))
for (i in unique(group$Group)) {
  otu_tax_group_i <- otu_tax[ ,rownames(subset(group, Group == i))]
  rownames(otu_tax_group_i) <- otu_tax$OTUID
  Nindividuals_S <- Nindividuals_S + rowMeans(otu_tax_group_i)  #计算 Nindividuals S
}

spec_occu <- NULL
for (i in unique(group$Group)) {
  otu_tax_group_i <- otu_tax[ ,rownames(subset(group, Group == i))]
  Nindividuals_SH <- apply(otu_tax_group_i, 1, mean)  #计算 Nindividuals SH
  Specificity <- Nindividuals_SH / Nindividuals_S  #计算 Specificity
  Nsites_H <- ncol(otu_tax_group_i)  #计算 Nsites H
  Nsites_SH <- apply(otu_tax_group_i, 1, function(x) sum(x>0))  #计算 Nsites SH
  Occupancy <- Nsites_SH / Nsites_H  #计算 Occupancy
  spec_occu_group_i <- data.frame(Group = i, OTU = otu_tax$OTUID, Specificity = Specificity, Occupancy = Occupancy, Abundance_mean = rowMeans(otu_tax_group_i), Taxonomy = otu_tax$Phylum)
  spec_occu <- rbind(spec_occu, spec_occu_group_i)  #合并各组统计
}
head(spec_occu)  #该数据框包含各组中各 OTU 的名称、特异性、占有率、平均丰度、类群信息等

#输出统计表格
write.table(spec_occu, 'spec_occu.txt', row.names = TRUE, sep = '\t', quote = FALSE)

#绘制 SPEC-OCCU 图
library(ggplot2)

p <- ggplot(spec_occu, aes(Occupancy, Specificity)) +
  geom_point(aes(size = log10(Abundance_mean), color = Taxonomy)) +
  #geom_jitter(aes(size = log10(Abundance_mean), color = Taxonomy)) +  #如果觉得普通散点图中的点重叠严重不好看，可以仿照 Gweon et al (2021) 使用抖动点图来展示
  scale_size(breaks = c(-1, -2, -3, -4), labels = c(expression(10^{-1}), expression(10^{-2}), expression(10^{-3}), expression(10^{-4})), range = c(0, 4)) +
  scale_color_manual(values = c('#E7272E', '#F59D1F', '#768CC5', '#9BC648', '#794779', '#A19E9D'), limits = c('Proteobacteria', 'Actinobacteria', 'Acidobacteria', 'Bacteroidetes', 'Firmicutes', 'Others')) +
  theme(panel.grid = element_blank(), panel.background = element_rect(fill = 'white', color = 'gray30'), legend.key = element_blank()) +
  facet_wrap(~Group, ncol = 2) +
  scale_x_continuous(breaks = c(0, 0.5, 1.0), expand = c(0, 0), limit = c(0, 1.2)) +
  scale_y_continuous(breaks = c(0, 0.5, 1.0), expand = c(0, 0), limit = c(0, 1.2)) +
  labs(x = 'Occupancy', y = 'Specificity', size = 'Average relative abundance', color = 'Taxonomy') +
  coord_cartesian(clip = 'off')

p


##识别各组样本中的特化种（specialist species）

#在上述统计表格中直接根据特异性和占有率 ≥0.7 做筛选，保留下的 OTU 即为特化种
spec_occu_specialist <- subset(spec_occu, Specificity >= 0.7 & Occupancy >= 0.7)
head(spec_occu_specialist)

write.table(spec_occu_specialist, 'spec_occu_specialist.txt', row.names = FALSE, sep = '\t', quote = FALSE)


#输出统计表格
#write.table(spec_occu_specialist, 'spec_occu_specialist.txt', row.names = FALSE, sep = '\t', quote = FALSE)

#在上述 SPEC-OCCU 图中添加阈值线
p1<-p + 
  geom_segment(aes(x = 0.7, xend = 1.2, y = 0.7, yend = 0.7), linetype = 2)+
  geom_segment(aes(x = 0.7, xend = 0.7, y = 0.7, yend = 1.2), linetype = 2)

p1

ggsave("./spec_occu.pdf", p1, width = 150, height = 100, units = "mm")
ggsave("./spec_occu.png", p1, width = 150, height = 100, units = "mm")
#ggsave("./spec_occu.tiff", p1, width = 150, height = 100, units = "mm")








####核心物种的识别——接近中心性和度####

####根据接近中心性和度划分不同子网络核心微生物，并出图####
##Meta_network
##由于在R语言中未计算出节点度的特征，故需要在gephi中计算好后,将节点数据表格输入
network_feature <- read.table("all_network_node_list-4.txt", row.names=1,sep="\t", header=T, blank.lines.skip=F, check.names=F)
network_feature
##如何将数据框的行名进行修改、重新命名
rownames(network_feature) <- network_feature$v_name
network_feature

network_feature$v_cols2 <- network_feature$v_cols
network_feature
####特别值得注意的是，前期对门水平颜色进行了自定义，但出图的时候挑选了特定的颜色组合，
#为满足出图颜色一致，需要修改颜色，但是CMYK颜色值数据表格读入会报错，只能修改了
network_feature[ rownames(network_feature)[network_feature$v_labels=="Gammaproteobacteria"], ]$v_cols2 <- "#394A92"
network_feature[ rownames(network_feature)[network_feature$v_labels=="Acidobacteria"], ]$v_cols2 <- "#68AC57"
network_feature[ rownames(network_feature)[network_feature$v_labels=="Bacteroidetes"], ]$v_cols2 <- "#497EB2"
network_feature[ rownames(network_feature)[network_feature$v_labels=="Alphaproteobacteria" ], ]$v_cols2 <- "#394A92"
network_feature[ rownames(network_feature)[network_feature$v_labels=="Chloroflexi"], ]$v_cols2 <- "#8E549E"
network_feature[ rownames(network_feature)[network_feature$v_labels=="Actinobacteria"], ]$v_cols2 <- "#F4C28F"
network_feature[ rownames(network_feature)[network_feature$v_labels=="Firmicutes"], ]$v_cols2 <- "#E600E6"
network_feature[ rownames(network_feature)[network_feature$v_labels=="Gemmatimonadetes" ], ]$v_cols2 <- "#D2352C"
network_feature[ rownames(network_feature)[network_feature$v_labels=="Fusobacteria"], ]$v_cols2 <- "lightgray"
network_feature[ rownames(network_feature)[network_feature$v_labels=="Deltaproteobacteria" ], ]$v_cols2 <- "#394A92"
network_feature[ rownames(network_feature)[network_feature$v_labels=="Patescibacteria"], ]$v_cols2 <- "lightgray"
network_feature[ rownames(network_feature)[network_feature$v_labels=="Verrucomicrobia"], ]$v_cols2 <- "lightgray"

network_feature$v_cols2

network_feature


##all_network##
#pdf(paste0(output,"all_network_CD.pdf"), width=7, height=7/6)
#par(mfrow=c(1,7), mar=c(2,3,1,0) )
#par(mfrow=c(1,1), mar=c(0,0,0,0))
pdf(paste0(output,"all_network_CD.pdf"), width=7, height=7/3)
par(mfrow=c(1,3), mar=c(2,3,2,0.5) )

plot(network_feature$closnesscentrality, network_feature$degree, 
     col= network_feature$v_cols2,
     pch = 20,cex =2.0,
     xlab = 'Closeness Centrality', ylab = 'Degree',
     #ylim = c (0,50),
     main="Meta-network")

dev.off()






##Ggt_network
##由于在R语言中未计算出节点度的特征，故需要在gephi中计算好后,将节点数据表格输入
Ggt_network_feature <- read.table("Ggt_network_node_list.txt", row.names=1,sep="\t", header=T, blank.lines.skip=F, check.names=F)
Ggt_network_feature
##如何将数据框的行名进行修改、重新命名
rownames(Ggt_network_feature) <- Ggt_network_feature$v_name
Ggt_network_feature

unique(Ggt_network_feature$v_labels)

Ggt_network_feature$v_cols2 <- Ggt_network_feature$v_cols
Ggt_network_feature
####特别值得注意的是，前期对门水平颜色进行了自定义，但出图的时候挑选了特定的颜色组合，
#为满足出图颜色一致，需要修改颜色，但是CMYK颜色值数据表格读入会报错，只能修改了
Ggt_network_feature[ rownames(Ggt_network_feature)[Ggt_network_feature$v_labels=="Gammaproteobacteria"], ]$v_cols2 <- "#394A92"
Ggt_network_feature[ rownames(Ggt_network_feature)[Ggt_network_feature$v_labels=="Acidobacteria"], ]$v_cols2 <- "#68AC57"
Ggt_network_feature[ rownames(Ggt_network_feature)[Ggt_network_feature$v_labels=="Bacteroidetes"], ]$v_cols2 <- "#497EB2"
Ggt_network_feature[ rownames(Ggt_network_feature)[Ggt_network_feature$v_labels=="Alphaproteobacteria" ], ]$v_cols2 <- "#394A92"
Ggt_network_feature[ rownames(Ggt_network_feature)[Ggt_network_feature$v_labels=="Chloroflexi"], ]$v_cols2 <- "#8E549E"
Ggt_network_feature[ rownames(Ggt_network_feature)[Ggt_network_feature$v_labels=="Actinobacteria"], ]$v_cols2 <- "#F4C28F"
Ggt_network_feature[ rownames(Ggt_network_feature)[Ggt_network_feature$v_labels=="Firmicutes"], ]$v_cols2 <- "#E600E6"
#Ggt_network_feature[ rownames(Ggt_network_feature)[Ggt_network_feature$v_labels=="Gemmatimonadetes" ], ]$v_cols2 <- "#D2352C"
Ggt_network_feature[ rownames(Ggt_network_feature)[Ggt_network_feature$v_labels=="Fusobacteria"], ]$v_cols2 <- "lightgray"
Ggt_network_feature[ rownames(Ggt_network_feature)[Ggt_network_feature$v_labels=="Deltaproteobacteria" ], ]$v_cols2 <- "#394A92"
Ggt_network_feature[ rownames(Ggt_network_feature)[Ggt_network_feature$v_labels=="Patescibacteria"], ]$v_cols2 <- "lightgray"
#Ggt_network_feature[ rownames(Ggt_network_feature)[Ggt_network_feature$v_labels=="Verrucomicrobia"], ]$v_cols2 <- "lightgray"

Ggt_network_feature$v_cols2

Ggt_network_feature


#pdf(paste0(output,"Ggt_network_CD.pdf"), width=7, height=7/6)
#par(mfrow=c(1,7), mar=c(2,3,1,0) )
#par(mfrow=c(1,1), mar=c(0,0,0,0))
pdf(paste0(output,"Ggt_network_CD.pdf"), width=7, height=7/3)
par(mfrow=c(1,3), mar=c(2,3,2,0.5) )

plot(Ggt_network_feature$closnesscentrality, Ggt_network_feature$degree, 
     col= Ggt_network_feature$v_cols2,
     pch = 20,cex =2.0,
     xlab = 'Closeness Centrality', ylab = 'Degree',
     #ylim = c (0,50),
     main="Ggt-network")

dev.off()




##S135_network
##由于在R语言中未计算出节点度的特征，故需要在gephi中计算好后,将节点数据表格输入
S135_network_feature <- read.table("135_network_node_list.txt", row.names=1,sep="\t", header=T, blank.lines.skip=F, check.names=F)
S135_network_feature
##如何将数据框的行名进行修改、重新命名
rownames(S135_network_feature) <- S135_network_feature$v_name
S135_network_feature



S135_network_feature$v_cols2 <- S135_network_feature$v_cols
S135_network_feature

unique(S135_network_feature$v_labels)

####特别值得注意的是，前期对门水平颜色进行了自定义，但出图的时候挑选了特定的颜色组合，
#为满足出图颜色一致，需要修改颜色，但是CMYK颜色值数据表格读入会报错，只能修改了
S135_network_feature[ rownames(S135_network_feature)[S135_network_feature$v_labels=="Gammaproteobacteria"], ]$v_cols2 <- "#394A92"
S135_network_feature[ rownames(S135_network_feature)[S135_network_feature$v_labels=="Acidobacteria"], ]$v_cols2 <- "#68AC57"
S135_network_feature[ rownames(S135_network_feature)[S135_network_feature$v_labels=="Bacteroidetes"], ]$v_cols2 <- "#497EB2"
S135_network_feature[ rownames(S135_network_feature)[S135_network_feature$v_labels=="Alphaproteobacteria" ], ]$v_cols2 <- "#394A92"
S135_network_feature[ rownames(S135_network_feature)[S135_network_feature$v_labels=="Chloroflexi"], ]$v_cols2 <- "#8E549E"
S135_network_feature[ rownames(S135_network_feature)[S135_network_feature$v_labels=="Actinobacteria"], ]$v_cols2 <- "#F4C28F"
S135_network_feature[ rownames(S135_network_feature)[S135_network_feature$v_labels=="Firmicutes"], ]$v_cols2 <- "#E600E6"
S135_network_feature[ rownames(S135_network_feature)[S135_network_feature$v_labels=="Gemmatimonadetes" ], ]$v_cols2 <- "#D2352C"
#S135_network_feature[ rownames(S135_network_feature)[S135_network_feature$v_labels=="Fusobacteria"], ]$v_cols2 <- "lightgray"
S135_network_feature[ rownames(S135_network_feature)[S135_network_feature$v_labels=="Deltaproteobacteria" ], ]$v_cols2 <- "#394A92"
S135_network_feature[ rownames(S135_network_feature)[S135_network_feature$v_labels=="Patescibacteria"], ]$v_cols2 <- "lightgray"
S135_network_feature[ rownames(S135_network_feature)[S135_network_feature$v_labels=="Verrucomicrobia"], ]$v_cols2 <- "lightgray"

S135_network_feature$v_cols2

S135_network_feature





#pdf(paste0(output,"S135_network_CD.pdf"), width=7, height=7/6)
#par(mfrow=c(1,7), mar=c(2,3,1,0) )
#par(mfrow=c(1,1), mar=c(0,0,0,0))
pdf(paste0(output,"S135_network_CD.pdf"), width=7, height=7/3)
par(mfrow=c(1,3), mar=c(2,3,2,0.5) )


plot(S135_network_feature$closnesscentrality, S135_network_feature$degree, 
     col= S135_network_feature$v_cols2,
     pch = 20,cex =2.0,
     xlab = 'Closeness Centrality', ylab = 'Degree',
     #ylim = c (0,50),
     main="S135-network")

dev.off()




##SA_network
##由于在R语言中未计算出节点度的特征，故需要在gephi中计算好后,将节点数据表格输入
SA_network_feature <- read.table("SA_network_node_list.txt", row.names=1,sep="\t", header=T, blank.lines.skip=F, check.names=F)
SA_network_feature
##如何将数据框的行名进行修改、重新命名
rownames(SA_network_feature) <- SA_network_feature$v_name
SA_network_feature

SA_network_feature$v_cols2 <- SA_network_feature$v_cols
SA_network_feature

unique(SA_network_feature$v_labels)

####特别值得注意的是，前期对门水平颜色进行了自定义，但出图的时候挑选了特定的颜色组合，
#为满足出图颜色一致，需要修改颜色，但是CMYK颜色值数据表格读入会报错，只能修改了
SA_network_feature[ rownames(SA_network_feature)[SA_network_feature$v_labels=="Gammaproteobacteria"], ]$v_cols2 <- "#394A92"
SA_network_feature[ rownames(SA_network_feature)[SA_network_feature$v_labels=="Acidobacteria"], ]$v_cols2 <- "#68AC57"
SA_network_feature[ rownames(SA_network_feature)[SA_network_feature$v_labels=="Bacteroidetes"], ]$v_cols2 <- "#497EB2"
SA_network_feature[ rownames(SA_network_feature)[SA_network_feature$v_labels=="Alphaproteobacteria" ], ]$v_cols2 <- "#394A92"
SA_network_feature[ rownames(SA_network_feature)[SA_network_feature$v_labels=="Chloroflexi"], ]$v_cols2 <- "#8E549E"
SA_network_feature[ rownames(SA_network_feature)[SA_network_feature$v_labels=="Actinobacteria"], ]$v_cols2 <- "#F4C28F"
SA_network_feature[ rownames(SA_network_feature)[SA_network_feature$v_labels=="Firmicutes"], ]$v_cols2 <- "#E600E6"
SA_network_feature[ rownames(SA_network_feature)[SA_network_feature$v_labels=="Gemmatimonadetes" ], ]$v_cols2 <- "#D2352C"
#SA_network_feature[ rownames(SA_network_feature)[SA_network_feature$v_labels=="Fusobacteria"], ]$v_cols2 <- "lightgray"
SA_network_feature[ rownames(SA_network_feature)[SA_network_feature$v_labels=="Deltaproteobacteria" ], ]$v_cols2 <- "#394A92"
SA_network_feature[ rownames(SA_network_feature)[SA_network_feature$v_labels=="Patescibacteria"], ]$v_cols2 <- "lightgray"
SA_network_feature[ rownames(SA_network_feature)[SA_network_feature$v_labels=="Verrucomicrobia"], ]$v_cols2 <- "lightgray"

SA_network_feature$v_cols2

SA_network_feature



#pdf(paste0(output,"SA_network_CD.pdf"), width=7, height=7/6)
#par(mfrow=c(1,7), mar=c(2,3,1,0) )
#par(mfrow=c(1,1), mar=c(0,0,0,0))
pdf(paste0(output,"SA_network_CD.pdf"), width=7, height=7/3)
par(mfrow=c(1,3), mar=c(2,3,2,0.5) )

plot(SA_network_feature$closnesscentrality, SA_network_feature$degree, 
     col= SA_network_feature$v_cols2,
     pch = 20,cex =2.0,
     xlab = 'Closeness Centrality', ylab = 'Degree',
     #ylim = c (0,50),
     main="SA-network")

dev.off()


##H_network
##由于在R语言中未计算出节点度的特征，故需要在gephi中计算好后,将节点数据表格输入
H_network_feature <- read.table("H_network_node_list.txt", row.names=1,sep="\t", header=T, blank.lines.skip=F, check.names=F)
H_network_feature
##如何将数据框的行名进行修改、重新命名
rownames(H_network_feature) <- H_network_feature$v_name
H_network_feature


H_network_feature$v_cols2 <- H_network_feature$v_cols
H_network_feature

unique(H_network_feature$v_labels)

####特别值得注意的是，前期对门水平颜色进行了自定义，但出图的时候挑选了特定的颜色组合，
#为满足出图颜色一致，需要修改颜色，但是CMYK颜色值数据表格读入会报错，只能修改了
H_network_feature[ rownames(H_network_feature)[H_network_feature$v_labels=="Gammaproteobacteria"], ]$v_cols2 <- "#394A92"
H_network_feature[ rownames(H_network_feature)[H_network_feature$v_labels=="Acidobacteria"], ]$v_cols2 <- "#68AC57"
H_network_feature[ rownames(H_network_feature)[H_network_feature$v_labels=="Bacteroidetes"], ]$v_cols2 <- "#497EB2"
H_network_feature[ rownames(H_network_feature)[H_network_feature$v_labels=="Alphaproteobacteria" ], ]$v_cols2 <- "#394A92"
H_network_feature[ rownames(H_network_feature)[H_network_feature$v_labels=="Chloroflexi"], ]$v_cols2 <- "#8E549E"
H_network_feature[ rownames(H_network_feature)[H_network_feature$v_labels=="Actinobacteria"], ]$v_cols2 <- "#F4C28F"
#H_network_feature[ rownames(H_network_feature)[H_network_feature$v_labels=="Firmicutes"], ]$v_cols2 <- "#E600E6"
H_network_feature[ rownames(H_network_feature)[H_network_feature$v_labels=="Gemmatimonadetes" ], ]$v_cols2 <- "#D2352C"
#H_network_feature[ rownames(H_network_feature)[H_network_feature$v_labels=="Fusobacteria"], ]$v_cols2 <- "lightgray"
H_network_feature[ rownames(H_network_feature)[H_network_feature$v_labels=="Deltaproteobacteria" ], ]$v_cols2 <- "#394A92"
H_network_feature[ rownames(H_network_feature)[H_network_feature$v_labels=="Patescibacteria"], ]$v_cols2 <- "lightgray"
H_network_feature[ rownames(H_network_feature)[H_network_feature$v_labels=="Verrucomicrobia"], ]$v_cols2 <- "lightgray"

H_network_feature$v_cols2

H_network_feature





#pdf(paste0(output,"H_network_CD.pdf"), width=7, height=7/6)
#par(mfrow=c(1,7), mar=c(2,3,1,0) )
#par(mfrow=c(1,1), mar=c(0,0,0,0))

pdf(paste0(output,"H_network_CD.pdf"), width=7, height=7/3)
par(mfrow=c(1,3), mar=c(2,3,2,0.5) )


plot(H_network_feature$closnesscentrality, H_network_feature$degree, 
     col= H_network_feature$v_cols2,
     pch = 20,cex =2.0,
     xlab = 'Closeness Centrality', ylab = 'Degree',
     #ylim = c (0,50),
     main="H-network")

dev.off()






####不同子网络核心微生物出图####
##根据度和介数中心性出图可知，核心微生物的划分标准为degree>60,且closeness>0.30

##barplot绘制堆叠柱状图
barplot(VADeaths, 
        col = NULL, #col: 设置条形底纹或者填充颜色。
        border =par("fg"),#border：设置条形边框颜色。如果设置为NA，则消除了边缘
        main = NULL, sub = NULL, 
        xlab = NULL,ylab = NULL, #xlab和ylab：设置x轴与y轴的lable
        xlim = NULL, ylim = NULL,  #xlim和ylim:设置图形x轴与y轴的范围。
        
        beside = FALSE,#beside:逻辑参数。如果FALSE，那么将绘画堆叠式的条形；如果是TRUE，将绘画并列式条形。
        horiz = FALSE, #horiz：逻辑参数。设置图形是水平或是垂直
        
        width = 1, #width：设置条形的宽度
        space = NULL, #space：设置各个条形间的间隔宽度。相当于各个条形宽度的一部分。
        names.arg = NULL,  #names.arg:设置条形标签（barlabels）。
        
        
        #density 和 angle : 设置柱子用线条填充，density 控制线条的密度， angel 控制线条的角度
        density = NULL,  #density:底纹的密度。默认值为NULL
        angle =45,  #angle：设置底纹的斜率
        
        
        axes = TRUE, #axes:逻辑参数。设置图形是否显示坐标轴。
        las=1,            #设置刻度值的方向,  0表示总是平行于坐标轴；1表示总是水平方向；2表示总是垂直于坐标轴；3表示总是垂直方向。
        
        yaxt= "s", #是否绘制Y坐标轴，s 绘制，n不绘制
        
        axisnames = TRUE, #axisnames：逻辑参数。设置是否显示x轴条形标签
        cex.axis=par("cex.axis"),#cex.axis:设置坐标轴数值的大小。
        cex.names=par("cex.axis"), #cex.names: 设置条形标签（barlabels）的大小。
        
        add = FALSE #add = “TRUE”将barplot加在目前已经有的图上
        
)


####Ggt_network####
#提取数据
Ggt_network_hub <- Ggt_network_feature[which(Ggt_network_feature$degree >= 60 & Ggt_network_feature$closnesscentrality >= 0.3),] 
Ggt_network_hub
dim(Ggt_network_hub) #Ggt_network_中没有核心微生物
write.table(Ggt_network_hub,paste0(output,"Ggt_network_hub.txt"),sep="\t",quote=F)

#从大网络边数据列表中提取出核心微生物所对应的边数据
Ggt_network_edge <- read.table("Ggt_network_edge_list.txt", row.names=1,sep="\t", header=T, blank.lines.skip=F, check.names=F)
#读取数据的时候，header=T表示保留原始表格的表头.row.names=1表示保留数据表格的第一列，=2表示去除表格的第一列
Ggt_network_edge

####循环的编写####
# Get unique labels
unique_labels <- unique(Ggt_network_hub$Label)

# Initialize empty matrix
Ggt_network_hub_matrix <- matrix(0, nrow=2, ncol=length(unique_labels), dimnames=list(c("positive","negative"), unique_labels))

# Loop through each label
for (i in unique_labels) {
  
  # Subset edges with the current label
  Ggt_network_hub_degree <- Ggt_network_edge[Ggt_network_edge$Source == i | Ggt_network_edge$Target == i,] 
  
  # Subset edges with positive correlation
  Ggt_network_hub_positive_degree <- Ggt_network_hub_degree[Ggt_network_hub_degree$e_cor >= 0,]
  
  
  # Subset edges with negative correlation
  Ggt_network_hub_negative_degree <- Ggt_network_hub_degree[Ggt_network_hub_degree$e_cor <= 0,]
  
  # Create matrix with positive and negative degrees
  Ggt_network_hub_i <- matrix(c(dim(Ggt_network_hub_positive_degree)[1], dim(Ggt_network_hub_negative_degree)[1]), nrow=2, ncol=1, byrow=TRUE, dimnames=list(c("positive","negative"), i))
  Ggt_network_hub_i
  # Fill in matrix with values from current label
  Ggt_network_hub_matrix[, i] <- Ggt_network_hub_i[, 1]
}

print(Ggt_network_hub_matrix)
colnames(Ggt_network_hub_matrix) <- rownames(Ggt_network_hub)
Ggt_network_hub_matrix
write.table(Ggt_network_hub_matrix,paste0(output,"Ggt_network_hub_matrix.txt"),sep="\t",quote=F)

##特别需要注意的是，由于该模块核心微生物为4个，少于10个，所以人工补全了表格进行出图，完了在AI中删除即可
Ggt_network_hub_matrix3 <- read.table("Ggt_network_hub_matrix2.txt", row.names=1,sep="\t", header=T, blank.lines.skip=F, check.names=F)
Ggt_network_hub_matrix3

#数据处理一下，以适合用barplot绘制柱状图
Ggt_network_hub_matrix <- as.data.frame(t(Ggt_network_hub_matrix))
Ggt_network_hub_matrix
sum <- apply(Ggt_network_hub_matrix,1,sum) 
Ggt_network_hub_matrix1 <- cbind(Ggt_network_hub_matrix,sum)
Ggt_network_hub_matrix1 <- as.data.frame(Ggt_network_hub_matrix1)
Ggt_network_hub_matrix1 <- Ggt_network_hub_matrix1[order(Ggt_network_hub_matrix1[,"sum"],decreasing = T),]
Ggt_network_hub_matrix1

Ggt_network_hub_matrix1 <- subset(Ggt_network_hub_matrix1, select = -sum)
Ggt_network_hub_matrix1 <- t(Ggt_network_hub_matrix1)


#如果需要增减物种的数目，请调整此处的数值
Ggt_network_hub_matrix2 <- Ggt_network_hub_matrix1[,1:10]
Ggt_network_hub_matrix2



##module0核心微生物出图
#pdf(paste0(output,"Ggt_network_top10_hubtaxa.pdf"), width=7, height=7/6)
#par(mfrow=c(1,4), mar=c(2,4,1,0))
#mycol<- c("green","red") 
pdf(paste0(output,"Ggt_network_top10_hubtaxa.pdf"), width=7, height=7/3)
par(mfrow=c(1,3), mar=c(5,3,2,0.5))
mycol<- c("#D55E86","#63D07F") 





p <- barplot(Ggt_network_hub_matrix3,
             las=2, border=NA, axes=T, cex.names=.5,ylim=c(0,130),main="Module0",
             col=mycol[1:nrow(Ggt_network_hub_matrix3)],ylab="Degree ")

## legend 
plot.new()
legend("top", bty="n", cex=0.5, #x.intersp=0.1, #y.intersp=0.75,
       legend=rev(rownames(Ggt_network_hub_matrix3)), 
       fill=rev(mycol), 
       border=rev(mycol) )
dev.off()




####S135_network####
#提取数据
S135_network_hub <- S135_network_feature[which(S135_network_feature$degree >= 60 & S135_network_feature$closnesscentrality >= 0.3),] 
S135_network_hub
dim(S135_network_hub) #S135_network中含有19个核心微生物
write.table(S135_network_hub,paste0(output,"S135_network_hub.txt"),sep="\t",quote=F)

#从大网络边数据列表中提取出核心微生物所对应的边数据
#原始的边列表数据第一列中有重复的名字，所以需要手动插入一列，然后读入
S135_network_edge <- read.table("135_network_edge_list2.txt", row.names=1,sep="\t", header=T, blank.lines.skip=F, check.names=F)
#读取数据的时候，header=T表示保留原始表格的表头.row.names=1表示保留数据表格的第一列，=2表示去除表格的第一列
S135_network_edge

####循环的编写####
# Get unique labels
unique_labels <- unique(S135_network_hub$Label)

# Initialize empty matrix
S135_network_hub_matrix <- matrix(0, nrow=2, ncol=length(unique_labels), dimnames=list(c("positive","negative"), unique_labels))

# Loop through each label
for (i in unique_labels) {
  
  # Subset edges with the current label
  S135_network_hub_degree <- S135_network_edge[S135_network_edge$Source == i | S135_network_edge$Target == i,] 
  
  # Subset edges with positive correlation
  S135_network_hub_positive_degree <- S135_network_hub_degree[S135_network_hub_degree$e_cor >= 0,]
  
  
  # Subset edges with negative correlation
  S135_network_hub_negative_degree <- S135_network_hub_degree[S135_network_hub_degree$e_cor <= 0,]
  
  # Create matrix with positive and negative degrees
  S135_network_hub_i <- matrix(c(dim(S135_network_hub_positive_degree)[1], dim(S135_network_hub_negative_degree)[1]), nrow=2, ncol=1, byrow=TRUE, dimnames=list(c("positive","negative"), i))
  S135_network_hub_i
  # Fill in matrix with values from current label
  S135_network_hub_matrix[, i] <- S135_network_hub_i[, 1]
}

print(S135_network_hub_matrix)
colnames(S135_network_hub_matrix) <- rownames(S135_network_hub)
S135_network_hub_matrix
write.table(S135_network_hub_matrix,paste0(output,"S135_network_hub_matrix.txt"),sep="\t",quote=F)

##特别需要注意的是，由于该模块核心微生物为4个，少于10个，所以人工补全了表格进行出图，完了在AI中删除即可
#S135_network_hub_matrix3 <- read.table("S135_network_hub_matrix2.txt", row.names=1,sep="\t", header=T, blank.lines.skip=F, check.names=F)
#S135_network_hub_matrix3

#数据处理一下，以适合用barplot绘制柱状图
S135_network_hub_matrix <- as.data.frame(t(S135_network_hub_matrix))
S135_network_hub_matrix
sum <- apply(S135_network_hub_matrix,1,sum) 
S135_network_hub_matrix1 <- cbind(S135_network_hub_matrix,sum)
S135_network_hub_matrix1 <- as.data.frame(S135_network_hub_matrix1)
S135_network_hub_matrix1 <- S135_network_hub_matrix1[order(S135_network_hub_matrix1[,"sum"],decreasing = T),]
S135_network_hub_matrix1

S135_network_hub_matrix1 <- subset(S135_network_hub_matrix1, select = -sum)
S135_network_hub_matrix1 <- t(S135_network_hub_matrix1)


#如果需要增减物种的数目，请调整此处的数值
S135_network_hub_matrix2 <- S135_network_hub_matrix1[,1:10]
S135_network_hub_matrix2



##S135_network核心微生物出图
pdf(paste0(output,"S135_network_top10_hubtaxa.pdf"), width=7, height=7/3)
par(mfrow=c(1,3), mar=c(7,3,2,0.5))
mycol<- c("#D55E86","#63D07F") 


p <- barplot(S135_network_hub_matrix2,
             las=2, border=NA, axes=T, cex.names=.5,ylim=c(0,130),main="X_14952B",
             col=mycol[1:nrow(S135_network_hub_matrix2)],ylab="Degree ")

## legend 
plot.new()
legend("top", bty="n", cex=0.5, #x.intersp=0.1, #y.intersp=0.75,
       legend=rev(rownames(S135_network_hub_matrix3)), 
       fill=rev(mycol), 
       border=rev(mycol) )
dev.off()



####SA_network####
#提取数据
SA_network_hub <- SA_network_feature[which(SA_network_feature$degree >= 60 & SA_network_feature$closnesscentrality >= 0.3),] 
SA_network_hub
dim(SA_network_hub) #SA_network中含有40个核心微生物
write.table(SA_network_hub,paste0(output,"SA_network_hub.txt"),sep="\t",quote=F)

#从大网络边数据列表中提取出核心微生物所对应的边数据
#原始的边列表数据第一列中有重复的名字，所以需要手动插入一列，然后读入
SA_network_edge <- read.table("SA_network_edge_list2.txt", row.names=1,sep="\t", header=T, blank.lines.skip=F, check.names=F)
#读取数据的时候，header=T表示保留原始表格的表头.row.names=1表示保留数据表格的第一列，=2表示去除表格的第一列
SA_network_edge

####循环的编写####
# Get unique labels
unique_labels <- unique(SA_network_hub$Label)

# Initialize empty matrix
SA_network_hub_matrix <- matrix(0, nrow=2, ncol=length(unique_labels), dimnames=list(c("positive","negative"), unique_labels))

# Loop through each label
for (i in unique_labels) {
  
  # Subset edges with the current label
  SA_network_hub_degree <- SA_network_edge[SA_network_edge$Source == i | SA_network_edge$Target == i,] 
  
  # Subset edges with positive correlation
  SA_network_hub_positive_degree <- SA_network_hub_degree[SA_network_hub_degree$e_cor >= 0,]
  
  
  # Subset edges with negative correlation
  SA_network_hub_negative_degree <- SA_network_hub_degree[SA_network_hub_degree$e_cor <= 0,]
  
  # Create matrix with positive and negative degrees
  SA_network_hub_i <- matrix(c(dim(SA_network_hub_positive_degree)[1], dim(SA_network_hub_negative_degree)[1]), nrow=2, ncol=1, byrow=TRUE, dimnames=list(c("positive","negative"), i))
  SA_network_hub_i
  # Fill in matrix with values from current label
  SA_network_hub_matrix[, i] <- SA_network_hub_i[, 1]
}

print(SA_network_hub_matrix)
colnames(SA_network_hub_matrix) <- rownames(SA_network_hub)
SA_network_hub_matrix
write.table(SA_network_hub_matrix,paste0(output,"SA_network_hub_matrix.txt"),sep="\t",quote=F)

##特别需要注意的是，由于该模块核心微生物为4个，少于10个，所以人工补全了表格进行出图，完了在AI中删除即可
#SA_network_hub_matrix3 <- read.table("SA_network_hub_matrix2.txt", row.names=1,sep="\t", header=T, blank.lines.skip=F, check.names=F)
#SA_network_hub_matrix3

#数据处理一下，以适合用barplot绘制柱状图
SA_network_hub_matrix <- as.data.frame(t(SA_network_hub_matrix))
SA_network_hub_matrix
sum <- apply(SA_network_hub_matrix,1,sum) 
SA_network_hub_matrix1 <- cbind(SA_network_hub_matrix,sum)
SA_network_hub_matrix1 <- as.data.frame(SA_network_hub_matrix1)
SA_network_hub_matrix1 <- SA_network_hub_matrix1[order(SA_network_hub_matrix1[,"sum"],decreasing = T),]
SA_network_hub_matrix1

SA_network_hub_matrix1 <- subset(SA_network_hub_matrix1, select = -sum)
SA_network_hub_matrix1 <- t(SA_network_hub_matrix1)


#如果需要增减物种的数目，请调整此处的数值
SA_network_hub_matrix2 <- SA_network_hub_matrix1[,1:10]
SA_network_hub_matrix2



##S135_network核心微生物出图
pdf(paste0(output,"SA_network_top10_hubtaxa.pdf"), width=7, height=7/3)
par(mfrow=c(1,3), mar=c(7,3,2,0.5))
mycol<- c("#D55E86","#63D07F") 


p <- barplot(SA_network_hub_matrix2,
             las=2, border=NA, axes=T, cex.names=.5,ylim=c(0,130),main="SA",
             col=mycol[1:nrow(SA_network_hub_matrix2)],ylab="Degree ")

## legend 
plot.new()
legend("top", bty="n", cex=0.5, #x.intersp=0.1, #y.intersp=0.75,
       legend=rev(rownames(SA_network_hub_matrix3)), 
       fill=rev(mycol), 
       border=rev(mycol) )
dev.off()



####H_network####
#提取数据
H_network_hub <- H_network_feature[which(H_network_feature$degree >= 60 & H_network_feature$closnesscentrality >= 0.3),] 
H_network_hub
dim(H_network_hub) #H_network中含有76个核心微生物
write.table(H_network_hub,paste0(output,"H_network_hub.txt"),sep="\t",quote=F)

#从大网络边数据列表中提取出核心微生物所对应的边数据
#原始的边列表数据第一列中有重复的名字，所以需要手动插入一列，然后读入
H_network_edge <- read.table("H_network_edge_list2.txt", row.names=1,sep="\t", header=T, blank.lines.skip=F, check.names=F)
#读取数据的时候，header=T表示保留原始表格的表头.row.names=1表示保留数据表格的第一列，=2表示去除表格的第一列
H_network_edge

####循环的编写####
# Get unique labels
unique_labels <- unique(H_network_hub$Label)

# Initialize empty matrix
H_network_hub_matrix <- matrix(0, nrow=2, ncol=length(unique_labels), dimnames=list(c("positive","negative"), unique_labels))

# Loop through each label
for (i in unique_labels) {
  
  # Subset edges with the current label
  H_network_hub_degree <- H_network_edge[H_network_edge$Source == i | H_network_edge$Target == i,] 
  
  # Subset edges with positive correlation
  H_network_hub_positive_degree <- H_network_hub_degree[H_network_hub_degree$e_cor >= 0,]
  
  
  # Subset edges with negative correlation
  H_network_hub_negative_degree <- H_network_hub_degree[H_network_hub_degree$e_cor <= 0,]
  
  # Create matrix with positive and negative degrees
  H_network_hub_i <- matrix(c(dim(H_network_hub_positive_degree)[1], dim(H_network_hub_negative_degree)[1]), nrow=2, ncol=1, byrow=TRUE, dimnames=list(c("positive","negative"), i))
  H_network_hub_i
  # Fill in matrix with values from current label
  H_network_hub_matrix[, i] <- H_network_hub_i[, 1]
}

print(H_network_hub_matrix)
colnames(H_network_hub_matrix) <- rownames(H_network_hub)
H_network_hub_matrix
write.table(H_network_hub_matrix,paste0(output,"H_network_hub_matrix.txt"),sep="\t",quote=F)

##特别需要注意的是，由于该模块核心微生物为4个，少于10个，所以人工补全了表格进行出图，完了在AI中删除即可
#H_network_hub_matrix3 <- read.table("H_network_hub_matrix2.txt", row.names=1,sep="\t", header=T, blank.lines.skip=F, check.names=F)
#H_network_hub_matrix3

#数据处理一下，以适合用barplot绘制柱状图
H_network_hub_matrix <- as.data.frame(t(H_network_hub_matrix))
H_network_hub_matrix
sum <- apply(H_network_hub_matrix,1,sum) 
H_network_hub_matrix1 <- cbind(H_network_hub_matrix,sum)
H_network_hub_matrix1 <- as.data.frame(H_network_hub_matrix1)
H_network_hub_matrix1 <- H_network_hub_matrix1[order(H_network_hub_matrix1[,"sum"],decreasing = T),]
H_network_hub_matrix1

H_network_hub_matrix1 <- subset(H_network_hub_matrix1, select = -sum)
H_network_hub_matrix1 <- t(H_network_hub_matrix1)


#如果需要增减物种的数目，请调整此处的数值
H_network_hub_matrix2 <- H_network_hub_matrix1[,1:10]
H_network_hub_matrix2



##S135_network核心微生物出图
pdf(paste0(output,"H_network_top10_hubtaxa.pdf"), width=7, height=7/3)
par(mfrow=c(1,3), mar=c(7,3,2,0.5))
mycol<- c("#D55E86","#63D07F") 


p <- barplot(H_network_hub_matrix2,
             las=2, border=NA, axes=T, cex.names=.5,ylim=c(0,130),main="H",
             col=mycol[1:nrow(H_network_hub_matrix2)],ylab="Degree ")

## legend 
plot.new()
legend("top", bty="n", cex=0.5, #x.intersp=0.1, #y.intersp=0.75,
       legend=rev(rownames(H_network_hub_matrix3)), 
       fill=rev(mycol), 
       border=rev(mycol) )
dev.off()







####R语言统计各模块内的节点分类数量并绘制饼图####

#读取节点属性列表
node <- read.delim('all_network_node_list-4.txt', stringsAsFactors = FALSE)
node

####门水平####
#查看网络中一共包含哪些门类群细菌
#并统计各门类群细菌所含节点数量，由大到小排个序
phylum <- table(node$v_phylum)
phylum <- sort(phylum, decreasing = TRUE)

#这个示例中共计 12 类细菌门，因此手动分配 12 种颜色代指
#大家使用自己的数据时，视情况修改颜色数量和赋值
color <- c('#394A92', '#68AC57', '#8E549E', '#F4C28F', '#497EB2', '#D2352C','#E600E6')
#, '#B3DE69', '#FDB462'
names(color) <- c(names(phylum[1]),names(phylum[2]),names(phylum[3]),names(phylum[4]),names(phylum[5]),names(phylum[6]),names(phylum[7]))

#names(color) <- names(c(phylum[1],phylum[2],phylum[3],phylum[4],phylum[5],phylum[6],phylum[7],phylum[8],phylum[9],phylum[10]))

color

#最后统计网络中各模块内，各门水平细菌所含节点数量
#并绘制饼图展示，通过 ggplot2 实现
library(ggplot2)

dir.create('plot', recursive = TRUE)  #创建目录“plot”用于存放图片

for (module in unique(network_feature$modularity_class)) {
  network_feature_module <- subset(network_feature, modularity_class == 5)
  network_feature_phylum <- data.frame(table(network_feature_module$v_phylum))#注意检查node_module中物种水平注释文件的命名
  
  
  p <- ggplot(network_feature_phylum, aes(x = '', y = Freq, fill = Var1)) + 
    geom_bar(stat = 'identity', width = 1) +
    coord_polar(theta = 'y') +
    scale_fill_manual(limits = names(color), values = color) + 
    theme(panel.grid = element_blank(), panel.background = element_blank(), 
          axis.text.x = element_blank()) +
    labs(x = '', y = '')
  ggsave(paste('plot/', module,'.pdf', sep = ''), p)
}


#备注：
#scale_fill_manual() 指定颜色时，limit 和 values 搭配，可以使颜色和微生物一一对应，即使数据中缺失该微生物也不会使顺序错乱
#输出的每个图片，文件名称均以模块名称命名


####纲水平####
#查看网络中一共包含哪些门类群细菌
#并统计各门类群细菌所含节点数量，由大到小排个序
class <- table(node$v_class)
class <- sort(class, decreasing = TRUE)

#这个示例中共计 15 类细菌纲，因此手动分配 12 种颜色代指
#大家使用自己的数据时，视情况修改颜色数量和赋值
color <- c('#FF7F00', '#FFFF33', '#984EA3', '#4DAF4A', '#377EB8', '#E41A1C',
           '#FFED6F', '#CCEBC5', '#BC80BD', '#FCCDE5', '#B3DE69', '#FDB462',"#CC6666","#999999","#808000")

#Palette <- c("#CC6666","#999999","#B2182B","#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7","#9999CC","#66CC99","#ADD1E5")
#cbbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442","#808000")

#names(color) <- c(names(phylum[1]),names(phylum[2]),names(phylum[3]),names(phylum[4]),names(phylum[5]),names(phylum[6]),names(phylum[7]),names(phylum[8]),names(phylum[9]),names(phylum[10]))

names(color) <- names(c(class[1],class[2],class[3],class[4],class[5],class[6],class[7],class[8],class[9],class[10],class[11],class[12],class[13],class[14],class[15]))

color

#最后统计网络中各模块内，各门水平细菌所含节点数量
#并绘制饼图展示，通过 ggplot2 实现
library(ggplot2)

dir.create('class_plot', recursive = TRUE)  #创建目录“plot”用于存放图片

for (module in unique(node$modularity_class)) {
  node_module <- subset(node, modularity_class == module)
  node_class <- data.frame(table(node_module$v_class))#注意检查node_module中物种水平注释文件的命名
  
  p <- ggplot(node_class, aes(x = '', y = Freq, fill = Var1)) + 
    geom_bar(stat = 'identity', width = 1) +
    coord_polar(theta = 'y') +
    scale_fill_manual(limits = names(color), values = color) + 
    theme(panel.grid = element_blank(), panel.background = element_blank(), 
          axis.text.x = element_blank()) +
    labs(x = '', y = '')
  ggsave(paste('class_plot/', module,'.pdf', sep = ''), p)
}

####属水平####
#查看网络中一共包含哪些门类群细菌
#并统计各门类群细菌所含节点数量，由大到小排个序
genus <- table(node$v_genus)
genus <- sort(genus, decreasing = TRUE)
genus
#这个示例中共计 15 类细菌纲，因此手动分配 12 种颜色代指
#大家使用自己的数据时，视情况修改颜色数量和赋值
color <- c('#FF7F00', '#FFFF33', '#984EA3', '#4DAF4A', '#377EB8', '#E41A1C',
           '#FFED6F', '#CCEBC5', '#BC80BD', '#FCCDE5', '#B3DE69', '#FDB462',"#CC6666","#999999","#808000")

#Palette <- c("#CC6666","#999999","#B2182B","#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7","#9999CC","#66CC99","#ADD1E5")
#cbbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442","#808000")

#names(color) <- c(names(phylum[1]),names(phylum[2]),names(phylum[3]),names(phylum[4]),names(phylum[5]),names(phylum[6]),names(phylum[7]),names(phylum[8]),names(phylum[9]),names(phylum[10]))

names(color) <- names(c(genus[1],genus[2],genus[3],genus[4],genus[5],genus[6],genus[7],genus[8],genus[9],genus[10],genus[11],genus[12],genus[13],genus[14],genus[15]))

color

#最后统计网络中各模块内，各门水平细菌所含节点数量
#并绘制饼图展示，通过 ggplot2 实现
library(ggplot2)

dir.create('genus_plot', recursive = TRUE)  #创建目录“plot”用于存放图片

for (module in unique(node$modularity_class)) {
  node_module <- subset(node, modularity_class == module)
  node_class <- data.frame(table(node_module$v_genus))#注意检查node_module中物种水平注释文件的命名
  
  p <- ggplot(node_class, aes(x = '', y = Freq, fill = Var1)) + 
    geom_bar(stat = 'identity', width = 1) +
    coord_polar(theta = 'y') +
    scale_fill_manual(limits = names(color), values = color) + 
    theme(panel.grid = element_blank(), panel.background = element_blank(), 
          axis.text.x = element_blank()) +
    labs(x = '', y = '')
  ggsave(paste('genus_plot/', module,'.pdf', sep = ''), p)
}







####对gephi结果进行模块的手动的划分

rm(list=ls())
setwd("E:/X-14952B/otus/networks2/R_Input1/")
output <- "E:/X-14952B/otus/networks2/R_Output1/"


####模块核心微生物#### 

##由于在R语言中未计算出节点度的特征，故需要在gephi中计算好后,将节点数据表格输入
network_feature <- read.table("all_network_node_list-4.txt", row.names=1,sep="\t", header=T, blank.lines.skip=F, check.names=F)
network_feature
##如何将数据框的行名进行修改、重新命名
rownames(network_feature) <- network_feature$v_name
network_feature
network_feature$v_cols2 <- network_feature$v_cols
network_feature
####特别值得注意的是，前期对门水平颜色进行了自定义，但出图的时候挑选了特定的颜色组合，
#为满足出图颜色一致，需要修改颜色，但是CMYK颜色值数据表格读入会报错，只能修改了
network_feature[ rownames(network_feature)[network_feature$v_labels=="Gammaproteobacteria"], ]$v_cols2 <- "#394A92"
network_feature[ rownames(network_feature)[network_feature$v_labels=="Acidobacteria"], ]$v_cols2 <- "#68AC57"
network_feature[ rownames(network_feature)[network_feature$v_labels=="Bacteroidetes"], ]$v_cols2 <- "#497EB2"
network_feature[ rownames(network_feature)[network_feature$v_labels=="Alphaproteobacteria" ], ]$v_cols2 <- "#394A92"
network_feature[ rownames(network_feature)[network_feature$v_labels=="Chloroflexi"], ]$v_cols2 <- "#8E549E"
network_feature[ rownames(network_feature)[network_feature$v_labels=="Actinobacteria"], ]$v_cols2 <- "#F4C28F"
network_feature[ rownames(network_feature)[network_feature$v_labels=="Firmicutes"], ]$v_cols2 <- "#E600E6"
network_feature[ rownames(network_feature)[network_feature$v_labels=="Gemmatimonadetes" ], ]$v_cols2 <- "#D2352C"
network_feature[ rownames(network_feature)[network_feature$v_labels=="Fusobacteria"], ]$v_cols2 <- "lightgray"
network_feature[ rownames(network_feature)[network_feature$v_labels=="Deltaproteobacteria" ], ]$v_cols2 <- "#394A92"
network_feature[ rownames(network_feature)[network_feature$v_labels=="Patescibacteria"], ]$v_cols2 <- "lightgray"
network_feature[ rownames(network_feature)[network_feature$v_labels=="Verrucomicrobia"], ]$v_cols2 <- "lightgray"

network_feature$v_cols2

network_feature


##查看模块的划分
unique(network_feature$modularity_class)#模块数量为0-10

####模块节点数据表格####

module0_gephi2 <- network_feature[network_feature$modularity_class=="0",]
module0_gephi2 
write.table(module0_gephi2,paste0(output,"All_module0_gephi2.txt"),sep="\t",quote=F)

####module0模块度####
#用于绘制模块度的箱线图
module0_degree <- as.data.frame(module0_gephi2$degree) 
 
write.table(module0_degree,paste0(output,"All_module0_degree.txt"),sep="\t",quote=F)



####module0网络提取####

sub_graph <- list()

select_module0 <- rownames(module0_gephi2)
sub_graph_module0 <- subgraph(root_net_16s, select_module0)

sub_graph_module0

#graphml 格式，可使用 gephi 软件打开并进行可视化编辑
write.graph(sub_graph_module0, 'sub_graph_module0.graphml', format = 'graphml')

#gml 格式，可使用 cytoscape 软件打开并进行可视化编辑
write.graph(sub_graph_module0, 'sub_graph_module0.gml', format = 'gml')



module1_gephi2 <- network_feature[network_feature$modularity_class=="1",]
module1_gephi2 
write.table(module1_gephi2,paste0(output,"All_module1_gephi2.txt"),sep="\t",quote=F)

####module1模块度####
#用于绘制模块度的箱线图
module1_degree <- as.data.frame(module1_gephi2$degree) 

write.table(module1_degree,paste0(output,"All_module1_degree.txt"),sep="\t",quote=F)





####module1网络提取####

sub_graph <- list()

select_module1 <- rownames(module1_gephi2)
sub_graph_module1 <- subgraph(root_net_16s, select_module1)

sub_graph_module1

#graphml 格式，可使用 gephi 软件打开并进行可视化编辑
write.graph(sub_graph_module1, 'sub_graph_module1.graphml', format = 'graphml')

#gml 格式，可使用 cytoscape 软件打开并进行可视化编辑
write.graph(sub_graph_module1, 'sub_graph_module1.gml', format = 'gml')




module2_gephi2 <- network_feature[network_feature$modularity_class=="2",]
module2_gephi2 
write.table(module2_gephi2,paste0(output,"All_module2_gephi2.txt"),sep="\t",quote=F)

####module2模块度####
#用于绘制模块度的箱线图
module2_degree <- as.data.frame(module2_gephi2$degree) 

write.table(module2_degree,paste0(output,"All_module2_degree.txt"),sep="\t",quote=F)


####module2网络提取####

sub_graph <- list()

select_module2 <- rownames(module2_gephi2)
sub_graph_module2 <- subgraph(root_net_16s, select_module2)

sub_graph_module2

#graphml 格式，可使用 gephi 软件打开并进行可视化编辑
write.graph(sub_graph_module2, 'sub_graph_module2.graphml', format = 'graphml')

#gml 格式，可使用 cytoscape 软件打开并进行可视化编辑
write.graph(sub_graph_module2, 'sub_graph_module2.gml', format = 'gml')




module3_gephi2 <- network_feature[network_feature$modularity_class=="3",]
module3_gephi2 
write.table(module3_gephi2,paste0(output,"All_module3_gephi2.txt"),sep="\t",quote=F)

####module3模块度####
#用于绘制模块度的箱线图
module3_degree <- as.data.frame(module3_gephi2$degree) 

write.table(module3_degree,paste0(output,"All_module3_degree.txt"),sep="\t",quote=F)

####module3网络提取####

sub_graph <- list()

select_module3 <- rownames(module3_gephi2)
sub_graph_module3 <- subgraph(root_net_16s, select_module3)

sub_graph_module3

#graphml 格式，可使用 gephi 软件打开并进行可视化编辑
write.graph(sub_graph_module3, 'sub_graph_module3.graphml', format = 'graphml')

#gml 格式，可使用 cytoscape 软件打开并进行可视化编辑
write.graph(sub_graph_module3, 'sub_graph_module3.gml', format = 'gml')





module4_gephi2 <- network_feature[network_feature$modularity_class=="4",]
module4_gephi2 
write.table(module4_gephi2,paste0(output,"All_module4_gephi2.txt"),sep="\t",quote=F)

####module4模块度####
#用于绘制模块度的箱线图
module4_degree <- as.data.frame(module4_gephi2$degree) 

write.table(module4_degree,paste0(output,"All_module4_degree.txt"),sep="\t",quote=F)



####module4网络提取####

sub_graph <- list()

select_module4 <- rownames(module4_gephi2)
sub_graph_module4 <- subgraph(root_net_16s, select_module4)

sub_graph_module4

#graphml 格式，可使用 gephi 软件打开并进行可视化编辑
write.graph(sub_graph_module4, 'sub_graph_module4.graphml', format = 'graphml')

#gml 格式，可使用 cytoscape 软件打开并进行可视化编辑
write.graph(sub_graph_module4, 'sub_graph_module4.gml', format = 'gml')




##提取模块

module5_gephi2 <- network_feature[network_feature$modularity_class=="5",]
module5_gephi2 
write.table(module5_gephi2,paste0(output,"All_module5_gephi2.txt"),sep="\t",quote=F)

####module4模块度####
#用于绘制模块度的箱线图
module5_degree <- as.data.frame(module5_gephi2$degree) 

write.table(module5_degree,paste0(output,"All_module5_degree.txt"),sep="\t",quote=F)



####module5网络提取####

sub_graph <- list()

select_module5 <- rownames(module5_gephi2)
sub_graph_module5 <- subgraph(root_net_16s, select_module5)

sub_graph_module5

#graphml 格式，可使用 gephi 软件打开并进行可视化编辑
write.graph(sub_graph_module5, 'sub_graph_module5.graphml', format = 'graphml')

#gml 格式，可使用 cytoscape 软件打开并进行可视化编辑
write.graph(sub_graph_module5, 'sub_graph_module5.gml', format = 'gml')





module6_gephi2 <- network_feature[network_feature$modularity_class=="6",]
module6_gephi2 
write.table(module6_gephi2,paste0(output,"All_module6_gephi2.txt"),sep="\t",quote=F)

####module4模块度####
#用于绘制模块度的箱线图
module6_degree <- as.data.frame(module6_gephi2$degree) 

write.table(module6_degree,paste0(output,"All_module6_degree.txt"),sep="\t",quote=F)


####module6网络提取####

sub_graph <- list()

select_module6 <- rownames(module6_gephi2)
sub_graph_module6 <- subgraph(root_net_16s, select_module6)

sub_graph_module6

#graphml 格式，可使用 gephi 软件打开并进行可视化编辑
write.graph(sub_graph_module6, 'sub_graph_module6.graphml', format = 'graphml')

#gml 格式，可使用 cytoscape 软件打开并进行可视化编辑
write.graph(sub_graph_module6, 'sub_graph_module6.gml', format = 'gml')




module7_gephi2 <- network_feature[network_feature$modularity_class=="7",]
module7_gephi2 
write.table(module7_gephi2,paste0(output,"All_module7_gephi2.txt"),sep="\t",quote=F)

####module4模块度####
#用于绘制模块度的箱线图
module7_degree <- as.data.frame(module7_gephi2$degree) 

write.table(module7_degree,paste0(output,"All_module7_degree.txt"),sep="\t",quote=F)

####module7网络提取####

sub_graph <- list()

select_module7 <- rownames(module7_gephi2)
sub_graph_module7 <- subgraph(root_net_16s, select_module7)

sub_graph_module7

#graphml 格式，可使用 gephi 软件打开并进行可视化编辑
write.graph(sub_graph_module7, 'sub_graph_module7.graphml', format = 'graphml')

#gml 格式，可使用 cytoscape 软件打开并进行可视化编辑
write.graph(sub_graph_module7, 'sub_graph_module7.gml', format = 'gml')






module8_gephi2 <- network_feature[network_feature$modularity_class=="8",]
module8_gephi2 
write.table(module8_gephi2,paste0(output,"All_module8_gephi2.txt"),sep="\t",quote=F)

####module4模块度####
#用于绘制模块度的箱线图
module8_degree <- as.data.frame(module8_gephi2$degree) 

write.table(module8_degree,paste0(output,"All_module8_degree.txt"),sep="\t",quote=F)


####module8网络提取####

sub_graph <- list()

select_module8 <- rownames(module8_gephi2)
sub_graph_module8 <- subgraph(root_net_16s, select_module8)

sub_graph_module8

#graphml 格式，可使用 gephi 软件打开并进行可视化编辑
write.graph(sub_graph_module8, 'sub_graph_module8.graphml', format = 'graphml')

#gml 格式，可使用 cytoscape 软件打开并进行可视化编辑
write.graph(sub_graph_module8, 'sub_graph_module8.gml', format = 'gml')








module9_gephi2 <- network_feature[network_feature$modularity_class=="9",]
module9_gephi2 
write.table(module9_gephi2,paste0(output,"All_module9_gephi2.txt"),sep="\t",quote=F)

####module4模块度####
#用于绘制模块度的箱线图
module9_degree <- as.data.frame(module9_gephi2$degree) 

write.table(module9_degree,paste0(output,"All_module9_degree.txt"),sep="\t",quote=F)

####module9网络提取####

sub_graph <- list()

select_module9 <- rownames(module9_gephi2)
sub_graph_module9 <- subgraph(root_net_16s, select_module9)

sub_graph_module9

#graphml 格式，可使用 gephi 软件打开并进行可视化编辑
write.graph(sub_graph_module9, 'sub_graph_module9.graphml', format = 'graphml')

#gml 格式，可使用 cytoscape 软件打开并进行可视化编辑
write.graph(sub_graph_module9, 'sub_graph_module9.gml', format = 'gml')






module10_gephi2 <- network_feature[network_feature$modularity_class=="10",]
module10_gephi2 
write.table(module10_gephi2,paste0(output,"All_module10_gephi2.txt"),sep="\t",quote=F)

####module4模块度####
#用于绘制模块度的箱线图
module10_degree <- as.data.frame(module10_gephi2$degree) 

write.table(module10_degree,paste0(output,"All_module10_degree.txt"),sep="\t",quote=F)


####module10网络提取####

sub_graph <- list()

select_module10 <- rownames(module10_gephi2)
sub_graph_module10 <- subgraph(root_net_16s, select_module10)

sub_graph_module10

#graphml 格式，可使用 gephi 软件打开并进行可视化编辑
write.graph(sub_graph_module10, 'sub_graph_module10.graphml', format = 'graphml')

#gml 格式，可使用 cytoscape 软件打开并进行可视化编辑
write.graph(sub_graph_module10, 'sub_graph_module10.gml', format = 'gml')







####根据接近中心性和度划分不同模块核心微生物，并出图####
pdf(paste0(output,"all_module0_CD4.pdf"), width=7, height=7/3)
par(mfrow=c(1,3), mar=c(2,3,2,0.5) )
#par(mfrow=c(1,1), mar=c(0,0,0,0))

plot(module0_gephi2$closnesscentrality, module0_gephi2$degree, 
     col= module0_gephi2$v_cols2,
     pch = 20,cex =2.0,
     xlab = 'Closeness Centrality', ylab = 'Degree',
     #ylim = c (0,50),
     main="Module 0")

dev.off()




pdf(paste0(output,"all_module1_CD4.pdf"), width=7, height=7/3)
par(mfrow=c(1,3), mar=c(2,2.5,2,0) )
#par(mfrow=c(1,1), mar=c(0,0,0,0))

plot(module1_gephi2$closnesscentrality, module1_gephi2$degree, 
     col= module1_gephi2$v_cols2,
     pch = 20,cex =2.0,
     xlab = 'Closeness Centrality', ylab = 'Degree',
     #ylim = c (0,50),
     main="Module 1")

dev.off()



pdf(paste0(output,"all_module2_CD4.pdf"), width=7, height=7/3)
par(mfrow=c(1,3), mar=c(2,3,2,0.5) )

plot(module2_gephi2$closnesscentrality, module2_gephi2$degree, 
     col= module2_gephi2$v_cols2,
     pch = 20,cex =2.0,
     xlab = '', ylab = '',
     #ylim = c (0,120),
     main="Module 2")
dev.off()


pdf(paste0(output,"all_module3_CD4.pdf"), width=7, height=7/3)
par(mfrow=c(1,3), mar=c(2,3,2,0.5) )
plot(module3_gephi2$closnesscentrality, module3_gephi2$degree, 
     col= module3_gephi2$v_cols2,
     pch = 20,cex =2.0,
     xlab = '', ylab = '',
     #ylim = c (0,60),
     main="Module 3")
dev.off()




pdf(paste0(output,"all_module4_CD4.pdf"), width=7, height=7/3)
par(mfrow=c(1,3), mar=c(2,3,2,0.5) )
plot(module4_gephi2$closnesscentrality, module4_gephi2$degree, 
     col= module4_gephi2$v_cols2,
     pch = 20,cex =2.0,
     xlab = '', ylab = '',
     #ylim = c (0,60),
     main="Module 4")
dev.off()



pdf(paste0(output,"all_module5_CD4.pdf"), width=7, height=7/3)
par(mfrow=c(1,3), mar=c(2,3,2,0.5) )
plot(module5_gephi2$closnesscentrality, module5_gephi2$degree, 
     col= module5_gephi2$v_cols2,
     pch = 20,cex =2.0,
     xlab = '', ylab = '',
     #ylim = c (0,60),
     main="Module 5")
dev.off()



pdf(paste0(output,"all_module6_CD4.pdf"), width=7, height=7/3)
par(mfrow=c(1,3), mar=c(2,3,2,0.5) )
plot(module6_gephi2$closnesscentrality, module6_gephi2$degree, 
     col= module6_gephi2$v_cols2,
     pch = 20,cex =2.0,
     xlab = '', ylab = '',
     #ylim = c (0,60),
     main="Module 6")
dev.off()



pdf(paste0(output,"all_module7_CD4.pdf"), width=7, height=7/3)
par(mfrow=c(1,3), mar=c(2,3,2,0.5) )
plot(module7_gephi2$closnesscentrality, module7_gephi2$degree, 
     col= module7_gephi2$v_cols2,
     pch = 20,cex =2.0,
     xlab = '', ylab = '',
     #ylim = c (0,60),
     main="Module 7")
dev.off()


pdf(paste0(output,"all_module8_CD4.pdf"), width=7, height=7/3)
par(mfrow=c(1,3), mar=c(2,3,2,0.5) )
plot(module8_gephi2$closnesscentrality, module8_gephi2$degree, 
     col= module8_gephi2$v_cols2,
     pch = 20,cex =2.0,
     xlab = '', ylab = '',
     #ylim = c (0,60),
     main="Module 8")
dev.off()



pdf(paste0(output,"all_module9_CD4.pdf"), width=7, height=7/3)
par(mfrow=c(1,3), mar=c(2,3,2,0.5) )
plot(module9_gephi2$closnesscentrality, module9_gephi2$degree, 
     col= module9_gephi2$v_cols2,
     pch = 20,cex =2.0,
     xlab = '', ylab = '',
     #ylim = c (0,60),
     main="Module 9")
dev.off()



pdf(paste0(output,"all_module10_CD4.pdf"), width=7, height=7/3)
par(mfrow=c(1,3), mar=c(2,3,2,0.5) )
plot(module10_gephi2$closnesscentrality, module10_gephi2$degree, 
     col= module10_gephi2$v_cols2,
     pch = 20,cex =2.0,
     xlab = '', ylab = '',
     #ylim = c (0,60),
     main="Module 10")
dev.off()



#总核心微生物出图
pdf(paste0(output,"all_CD4.pdf"), width=7, height=7/3)
par(mfrow=c(1,3), mar=c(2,3,2,0.5) )
plot(network_feature$closnesscentrality, network_feature$degree, 
     col= network_feature$v_cols2,
     pch = 20,cex =2.0,
     xlab = '', ylab = '',
     ylim = c (0,120),
     main="All")

dev.off()



####不同模块核心微生物出图####
##根据度和介数中心性出图可知，核心微生物的划分标准为degree>60,且closeness>0.30


##barplot绘制堆叠柱状图
barplot(VADeaths, 
        col = NULL, #col: 设置条形底纹或者填充颜色。
        border =par("fg"),#border：设置条形边框颜色。如果设置为NA，则消除了边缘
        main = NULL, sub = NULL, 
        xlab = NULL,ylab = NULL, #xlab和ylab：设置x轴与y轴的lable
        xlim = NULL, ylim = NULL,  #xlim和ylim:设置图形x轴与y轴的范围。
        
        beside = FALSE,#beside:逻辑参数。如果FALSE，那么将绘画堆叠式的条形；如果是TRUE，将绘画并列式条形。
        horiz = FALSE, #horiz：逻辑参数。设置图形是水平或是垂直
        
        width = 1, #width：设置条形的宽度
        space = NULL, #space：设置各个条形间的间隔宽度。相当于各个条形宽度的一部分。
        names.arg = NULL,  #names.arg:设置条形标签（barlabels）。
        
        
        #density 和 angle : 设置柱子用线条填充，density 控制线条的密度， angel 控制线条的角度
        density = NULL,  #density:底纹的密度。默认值为NULL
        angle =45,  #angle：设置底纹的斜率
        
        
        axes = TRUE, #axes:逻辑参数。设置图形是否显示坐标轴。
        las=1,            #设置刻度值的方向,  0表示总是平行于坐标轴；1表示总是水平方向；2表示总是垂直于坐标轴；3表示总是垂直方向。
        
        yaxt= "s", #是否绘制Y坐标轴，s 绘制，n不绘制
        
        axisnames = TRUE, #axisnames：逻辑参数。设置是否显示x轴条形标签
        cex.axis=par("cex.axis"),#cex.axis:设置坐标轴数值的大小。
        cex.names=par("cex.axis"), #cex.names: 设置条形标签（barlabels）的大小。
        
        add = FALSE #add = “TRUE”将barplot加在目前已经有的图上
        
)


####module0####
#提取数据
module0_hub <- module0_gephi2[which(module0_gephi2$degree >= 60 & module0_gephi2$closnesscentrality >= 0.3),] 
module0_hub
dim(module0_hub) #module0中含有4个核心微生物，分别为OTU186（n347）、OTU208（n408）、OTU324(n623)和OTU389（n707）
write.table(module0_hub,paste0(output,"module0_hub.txt"),sep="\t",quote=F)

#从大网络边数据列表中提取出核心微生物所对应的边数据
all_network_edge <- read.table("all_network_edge_list-4-2.txt", row.names=1,sep="\t", header=T, blank.lines.skip=F, check.names=F)
#读取数据的时候，header=T表示保留原始表格的表头.row.names=1表示保留数据表格的第一列，=2表示去除表格的第一列
all_network_edge

####循环的编写####
# Get unique labels
unique_labels <- unique(module0_hub$Label)

# Initialize empty matrix
module0_hub_matrix <- matrix(0, nrow=2, ncol=length(unique_labels), dimnames=list(c("positive","negative"), unique_labels))

# Loop through each label
for (i in unique_labels) {
  
  # Subset edges with the current label
  module0_hub_degree <- all_network_edge[all_network_edge$Source == i | all_network_edge$Target == i,] 
  
  # Subset edges with positive correlation
  module0_hub_positive_degree <- module0_hub_degree[module0_hub_degree$e_cor >= 0,]
  
  
  # Subset edges with negative correlation
  module0_hub_negative_degree <- module0_hub_degree[module0_hub_degree$e_cor <= 0,]
  
  # Create matrix with positive and negative degrees
  module0_hub_i <- matrix(c(dim(module0_hub_positive_degree)[1], dim(module0_hub_negative_degree)[1]), nrow=2, ncol=1, byrow=TRUE, dimnames=list(c("positive","negative"), i))
  module0_hub_i
  # Fill in matrix with values from current label
  module0_hub_matrix[, i] <- module0_hub_i[, 1]
}

print(module0_hub_matrix)
colnames(module0_hub_matrix) <- rownames(module0_hub)
module0_hub_matrix
write.table(module0_hub_matrix,paste0(output,"module0_hub_matrix.txt"),sep="\t",quote=F)


####循环的拆解####

#module0_hub_n347_degree <- all_network_edge[all_network_edge$Source == "n347" | all_network_edge$Target == "n347",] 
#module0_hub_n347_degree
#dim(module0_hub_n347_degree)#节点n347含有76条边
#正度
#module0_hub_n347_positive_degree <- module0_hub_n347_degree[module0_hub_n347_degree$e_cor >=0,]
#module0_hub_n347_positive_degree 
#dim(module0_hub_n347_positive_degree)#7

#负度
#module0_hub_n347_negative_degree <- module0_hub_n347_degree[module0_hub_n347_degree$e_cor <=0,]
#module0_hub_n347_negative_degree 
#dim(module0_hub_n347_negative_degree)#69


#module0_hub_n408_degree <- all_network_edge[all_network_edge$Source == "n408" | all_network_edge$Target == "n408",] 
#module0_hub_n408_degree
#dim(module0_hub_n408_degree)#75条边

#正度
#module0_hub_n408_positive_degree <- module0_hub_n408_degree[module0_hub_n408_degree$e_cor >=0,]
#module0_hub_n408_positive_degree 
#dim(module0_hub_n408_positive_degree)#60

#负度
#module0_hub_n408_negative_degree <- module0_hub_n408_degree[module0_hub_n408_degree$e_cor <=0,]
#module0_hub_n408_negative_degree 
#dim(module0_hub_n408_negative_degree)#15


#module0_hub_n623_degree <- all_network_edge[all_network_edge$Source == "n623" | all_network_edge$Target == "n623",] 
#module0_hub_n623_degree
#dim(module0_hub_n623_degree)#69条边

#正度
#module0_hub_n623_positive_degree <- module0_hub_n623_degree[module0_hub_n623_degree$e_cor >=0,]
#module0_hub_n623_positive_degree 
#dim(module0_hub_n623_positive_degree)#62

#负度
#module0_hub_n623_negative_degree <- module0_hub_n623_degree[module0_hub_n623_degree$e_cor <=0,]
#module0_hub_n623_negative_degree 
#dim(module0_hub_n623_negative_degree)#7


#module0_hub_n707_degree <- all_network_edge[all_network_edge$Source == "n707" | all_network_edge$Target == "n707",] 
#module0_hub_n707_degree
#dim(module0_hub_n707_degree)#63条边

#正度
#module0_hub_n707_positive_degree <- module0_hub_n707_degree[module0_hub_n707_degree$e_cor >=0,]
#module0_hub_n707_positive_degree 
#dim(module0_hub_n707_positive_degree)#55

#负度
#module0_hub_n707_negative_degree <- module0_hub_n707_degree[module0_hub_n707_degree$e_cor <=0,]
#module0_hub_n707_negative_degree 
#dim(module0_hub_n707_negative_degree)#8






##特别需要注意的是，由于该模块核心微生物为4个，少于10个，所以人工补全了表格进行出图，完了在AI中删除即可
#module0_hub_matrix3 <- read.table("module0_hub_matrix2.txt", row.names=1,sep="\t", header=T, blank.lines.skip=F, check.names=F)
#module0_hub_matrix3


#数据处理一下，以适合用barplot绘制柱状图
module0_hub_matrix <- as.data.frame(t(module0_hub_matrix))
module0_hub_matrix
sum <- apply(module0_hub_matrix,1,sum) 
module0_hub_matrix1 <- cbind(module0_hub_matrix,sum)
module0_hub_matrix1 <- as.data.frame(module0_hub_matrix1)
module0_hub_matrix1 <- module0_hub_matrix1[order(module0_hub_matrix1[,"sum"],decreasing = T),]
module0_hub_matrix1

module0_hub_matrix1 <- subset(module0_hub_matrix1, select = -sum)
module0_hub_matrix1 <- t(module0_hub_matrix1)
module0_hub_matrix1

#如果需要增减物种的数目，请调整此处的数值
module0_hub_matrix2 <- module0_hub_matrix1[,1:10]
module0_hub_matrix2



#数据处理一下，以适合用barplot绘制柱状图
#module0_hub_matrix3 <-module0_hub_matrix3[,1:ncol(module0_hub_matrix3)]
#rownames(module0_hub_matrix3)<-module0_hub_matrix3[,1]
#module0_hub_matrix3<-as.matrix(module0_hub_matrix3)
#module0_hub_matrix3
#colnames(module0_hub_matrix3)


##module0核心微生物出图
pdf(paste0(output,"module0_top10_hubtaxa.pdf"), width=7, height=7/3)
par(mfrow=c(1,3), mar=c(5,3,2,0.5))
mycol<- c("#D55E86","#63D07F") 

p <- barplot(module0_hub_matrix2,
             las=2, border=NA, axes=T, cex.names=.5,ylim=c(0,130),main="Module0",
             col=mycol[1:nrow(module0_hub_matrix2)],ylab="Degree ")

## legend 
plot.new()
legend("top", bty="n", cex=1.0, x.intersp=0.1, y.intersp=0.75,
       legend=rev(rownames(module0_hub_matrix3)), 
       fill=rev(mycol), 
       border=rev(mycol) )
dev.off()


####module0正负边/度####

##方法一：根据从总网络的提取结果
#136
unique_labels <- unique(module0_gephi2$Label)

# Initialize empty matrix
module0_all_matrix <- matrix(0, nrow=2, ncol=length(unique_labels), dimnames=list(c("positive","negative"), unique_labels))

# Loop through each label
for (i in unique_labels) {
  
  # Subset edges with the current label
  module0_all_degree <- all_network_edge[all_network_edge$Source == i | all_network_edge$Target == i,] 
  
  # Subset edges with positive correlation
  module0_all_positive_degree <- module0_all_degree[module0_all_degree$e_cor >= 0,]
  
  
  # Subset edges with negative correlation
  module0_all_negative_degree <- module0_all_degree[module0_all_degree$e_cor <= 0,]
  
  # Create matrix with positive and negative degrees
  module0_all_i <- matrix(c(dim(module0_all_positive_degree)[1], dim(module0_all_negative_degree)[1]), nrow=2, ncol=1, byrow=TRUE, dimnames=list(c("positive","negative"), i))
  module0_all_i
  # Fill in matrix with values from current label
  module0_all_matrix[, i] <- module0_all_i[, 1]
}

print(module0_all_matrix)
colnames(module0_all_matrix) <- rownames(module0_gephi2)
module0_all_matrix

dim(module0_all_matrix)
write.table(module0_all_matrix,paste0(output,"module0_all_matrix.txt"),sep="\t",quote=F)


##方法二：基于网络重新计算得到的数据

module0_network_edge <- read.table("module0_edge_list2.txt", row.names=1,sep="\t", header=T, blank.lines.skip=F, check.names=F)
#读取数据的时候，header=T表示保留原始表格的表头.row.names=1表示保留数据表格的第一列，=2表示去除表格的第一列
module0_network_edge


# Initialize empty matrix
module0_all_matrix2 <- matrix(0, nrow=2, ncol=length(unique_labels), dimnames=list(c("positive","negative"), unique_labels))

# Loop through each label
for (i in unique_labels) {
  
  # Subset edges with the current label
  module0_all_degree <- module0_network_edge[module0_network_edge$Source == i | module0_network_edge$Target == i,] 
  
  # Subset edges with positive correlation
  module0_all_positive_degree <- module0_all_degree[module0_all_degree$e_cor >= 0,]
  
  
  # Subset edges with negative correlation
  module0_all_negative_degree <- module0_all_degree[module0_all_degree$e_cor <= 0,]
  
  # Create matrix with positive and negative degrees
  module0_all_i <- matrix(c(dim(module0_all_positive_degree)[1], dim(module0_all_negative_degree)[1]), nrow=2, ncol=1, byrow=TRUE, dimnames=list(c("positive","negative"), i))
  module0_all_i
  # Fill in matrix with values from current label
  module0_all_matrix2[, i] <- module0_all_i[, 1]
}

print(module0_all_matrix2)
colnames(module0_all_matrix2) <- rownames(module0_gephi2)
module0_all_matrix2

dim(module0_all_matrix2)

write.table(module0_all_matrix2,paste0(output,"module0_all_matrix2.txt"),sep="\t",quote=F)





####module1
#提取数据
module1_hub <- module1_gephi2[which(module1_gephi2$degree >= 60 & module1_gephi2$closnesscentrality >= 0.3),] 
module1_hub
dim(module1_hub) #module1中含有19个核心微生物
write.table(module1_hub,paste0(output,"module1_hub.txt"),sep="\t",quote=F)

#从大网络边数据列表中提取出核心微生物所对应的边数据
all_network_edge <- read.table("all_network_edge_list-4-2.txt", row.names=1,sep="\t", header=T, blank.lines.skip=F, check.names=F)
#读取数据的时候，header=T表示保留原始表格的表头.row.names=1表示保留数据表格的第一列，=2表示去除表格的第一列
all_network_edge

####循环的编写####
# Get unique labels
unique_labels <- unique(module1_hub$Label)
unique_labels
# Initialize empty matrix
module1_hub_matrix <- matrix(0, nrow=2, ncol=length(unique_labels), dimnames=list(c("positive","negative"), unique_labels))

# Loop through each label
for (i in unique_labels) {
  
  # Subset edges with the current label
  module1_hub_degree <- all_network_edge[all_network_edge$Source == i | all_network_edge$Target == i,] 
  
  # Subset edges with positive correlation
  module1_hub_positive_degree <- module1_hub_degree[module1_hub_degree$e_cor >= 0,]
  
  
  # Subset edges with negative correlation
  module1_hub_negative_degree <- module1_hub_degree[module1_hub_degree$e_cor <= 0,]
  
  # Create matrix with positive and negative degrees
  module1_hub_i <- matrix(c(dim(module1_hub_positive_degree)[1], dim(module1_hub_negative_degree)[1]), nrow=2, ncol=1, byrow=TRUE, dimnames=list(c("positive","negative"), i))
  module1_hub_i
  # Fill in matrix with values from current label
  module1_hub_matrix[, i] <- module1_hub_i[, 1]
}

print(module1_hub_matrix)
colnames(module1_hub_matrix) <- rownames(module1_hub)
module1_hub_matrix
write.table(module1_hub_matrix,paste0(output,"module1_hub_matrix.txt"),sep="\t",quote=F)


#数据处理一下，以适合用barplot绘制柱状图
module1_hub_matrix <- as.data.frame(t(module1_hub_matrix))
module1_hub_matrix
sum <- apply(module1_hub_matrix,1,sum) 
module1_hub_matrix1 <- cbind(module1_hub_matrix,sum)
module1_hub_matrix1 <- as.data.frame(module1_hub_matrix1)
module1_hub_matrix1 <- module1_hub_matrix1[order(module1_hub_matrix1[,"sum"],decreasing = T),]
module1_hub_matrix1

module1_hub_matrix1 <- subset(module1_hub_matrix1, select = -sum)
module1_hub_matrix1 <- t(module1_hub_matrix1)


#如果需要增减物种的数目，请调整此处的数值
module1_hub_matrix2 <- module1_hub_matrix1[,1:10]
module1_hub_matrix2


##module1核心微生物出图
pdf(paste0(output,"module1_top10_hubtaxa.pdf"), width=7, height=7/3)
par(mfrow=c(1,3), mar=c(5,3,2,0.5))
mycol<- c("#D55E86","#63D07F") 


p1 <- barplot(module1_hub_matrix2,
              las=2, border=NA, axes=T, cex.names=.5,ylim=c(0,130),main="Module1",
              col=mycol[1:nrow(module1_hub_matrix2)],ylab="Degree ")

## legend 
plot.new()
legend("top", bty="n", cex=0.5, #x.intersp=0.1, #y.intersp=0.75,
       legend=rev(rownames(module1_hub_matrix2)), 
       fill=rev(mycol), 
       border=rev(mycol) )
dev.off()





####module1正负边/度####

##方法一：根据从总网络的提取结果
#162
unique_labels <- unique(module1_gephi2$Label)

# Initialize empty matrix
module1_all_matrix <- matrix(0, nrow=2, ncol=length(unique_labels), dimnames=list(c("positive","negative"), unique_labels))

# Loop through each label
for (i in unique_labels) {
  
  # Subset edges with the current label
  module1_all_degree <- all_network_edge[all_network_edge$Source == i | all_network_edge$Target == i,] 
  
  # Subset edges with positive correlation
  module1_all_positive_degree <- module1_all_degree[module1_all_degree$e_cor >= 0,]
  
  
  # Subset edges with negative correlation
  module1_all_negative_degree <- module1_all_degree[module1_all_degree$e_cor <= 0,]
  
  # Create matrix with positive and negative degrees
  module1_all_i <- matrix(c(dim(module1_all_positive_degree)[1], dim(module1_all_negative_degree)[1]), nrow=2, ncol=1, byrow=TRUE, dimnames=list(c("positive","negative"), i))
  module1_all_i
  # Fill in matrix with values from current label
  module1_all_matrix[, i] <- module1_all_i[, 1]
}

print(module1_all_matrix)
colnames(module1_all_matrix) <- rownames(module1_gephi2)
module1_all_matrix

dim(module1_all_matrix)
write.table(module1_all_matrix,paste0(output,"module1_all_matrix.txt"),sep="\t",quote=F)


##方法二：基于网络重新计算得到的数据

module1_network_edge <- read.table("module1_edge_list2.txt", row.names=1,sep="\t", header=T, blank.lines.skip=F, check.names=F)
#读取数据的时候，header=T表示保留原始表格的表头.row.names=1表示保留数据表格的第一列，=2表示去除表格的第一列
module1_network_edge


# Initialize empty matrix
module1_all_matrix2 <- matrix(0, nrow=2, ncol=length(unique_labels), dimnames=list(c("positive","negative"), unique_labels))

# Loop through each label
for (i in unique_labels) {
  
  # Subset edges with the current label
  module1_all_degree <- module1_network_edge[module1_network_edge$Source == i | module1_network_edge$Target == i,] 
  
  # Subset edges with positive correlation
  module1_all_positive_degree <- module1_all_degree[module1_all_degree$e_cor >= 0,]
  
  
  # Subset edges with negative correlation
  module1_all_negative_degree <- module1_all_degree[module1_all_degree$e_cor <= 0,]
  
  # Create matrix with positive and negative degrees
  module1_all_i <- matrix(c(dim(module1_all_positive_degree)[1], dim(module1_all_negative_degree)[1]), nrow=2, ncol=1, byrow=TRUE, dimnames=list(c("positive","negative"), i))
  module1_all_i
  # Fill in matrix with values from current label
  module1_all_matrix2[, i] <- module1_all_i[, 1]
}

print(module1_all_matrix2)
colnames(module1_all_matrix2) <- rownames(module1_gephi2)
module1_all_matrix2

dim(module1_all_matrix2)

write.table(module1_all_matrix2,paste0(output,"module1_all_matrix2.txt"),sep="\t",quote=F)






####module2
#提取数据
module2_hub <- module2_gephi2[which(module2_gephi2$degree >= 60 & module2_gephi2$closnesscentrality >= 0.3),] 
module2_hub
dim(module2_hub) #module2中没有核心微生物存在，故不运行下列代码
write.table(module2_hub,paste0(output,"module2_hub.txt"),sep="\t",quote=F)

#从大网络边数据列表中提取出核心微生物所对应的边数据
all_network_edge <- read.table("all_network_edge_list-4-2.txt", row.names=1,sep="\t", header=T, blank.lines.skip=F, check.names=F)
#读取数据的时候，header=T表示保留原始表格的表头.row.names=1表示保留数据表格的第一列，=2表示去除表格的第一列
all_network_edge

####循环的编写####
# Get unique labels
unique_labels <- unique(module2_hub$Label)
unique_labels
# Initialize empty matrix
module2_hub_matrix <- matrix(0, nrow=2, ncol=length(unique_labels), dimnames=list(c("positive","negative"), unique_labels))

# Loop through each label
for (i in unique_labels) {
  
  # Subset edges with the current label
  module2_hub_degree <- all_network_edge[all_network_edge$Source == i | all_network_edge$Target == i,] 
  
  # Subset edges with positive correlation
  module2_hub_positive_degree <- module2_hub_degree[module1_hub_degree$e_cor >= 0,]
  
  
  # Subset edges with negative correlation
  module2_hub_negative_degree <- module2_hub_degree[module1_hub_degree$e_cor <= 0,]
  
  # Create matrix with positive and negative degrees
  module2_hub_i <- matrix(c(dim(module2_hub_positive_degree)[1], dim(module2_hub_negative_degree)[1]), nrow=2, ncol=1, byrow=TRUE, dimnames=list(c("positive","negative"), i))
  module2_hub_i
  # Fill in matrix with values from current label
  module2_hub_matrix[, i] <- module2_hub_i[, 1]
}

print(module2_hub_matrix)
colnames(module2_hub_matrix) <- rownames(module2_hub)
module2_hub_matrix
write.table(module2_hub_matrix,paste0(output,"module2_hub_matrix.txt"),sep="\t",quote=F)



#数据处理一下，以适合用barplot绘制柱状图
module2_hub_matrix <- as.data.frame(t(module2_hub_matrix))
module2_hub_matrix
sum <- apply(module2_hub_matrix,1,sum) 
module2_hub_matrix1 <- cbind(module2_hub_matrix,sum)
module2_hub_matrix1 <- as.data.frame(module2_hub_matrix1)
module2_hub_matrix1 <- module2_hub_matrix2[order(module2_hub_matrix1[,"sum"],decreasing = T),]
module2_hub_matrix1

module2_hub_matrix1 <- subset(module2_hub_matrix1, select = -sum)
module2_hub_matrix1 <- t(module2_hub_matrix1)


##特别需要注意的是，由于该模块核心微生物为4个，少于10个，所以人工补全了表格进行出图，完了在AI中删除即可
module2_hub_matrix3 <- read.table("module2_hub_matrix2.txt", row.names=1,sep="\t", header=T, blank.lines.skip=F, check.names=F)
module2_hub_matrix3


#如果需要增减物种的数目，请调整此处的数值
module2_hub_matrix2 <- module2_hub_matrix1[,1:10]
module2_hub_matrix2


##module2核心微生物出图
pdf(paste0(output,"module2_top10_hubtaxa.pdf"), width=7, height=7/3)
par(mfrow=c(1,3), mar=c(5,3,2,0.5))
mycol<- c("red","green") 


p1 <- barplot(module2_hub_matrix2,
              las=2, border=NA, axes=T, cex.names=.5,ylim=c(0,130),main="Module1",
              col=mycol[1:nrow(module2_hub_matrix2)],ylab="Degree ")

## legend 
plot.new()
legend("top", bty="n", cex=0.5, #x.intersp=0.1, #y.intersp=0.75,
       legend=rev(rownames(module2_hub_matrix2)), 
       fill=rev(mycol), 
       border=rev(mycol) )
dev.off()




####module2正负边/度####

##方法一：根据从总网络的提取结果
#74
unique_labels <- unique(module2_gephi2$Label)
unique_labels
# Initialize empty matrix
module2_all_matrix <- matrix(0, nrow=2, ncol=length(unique_labels), dimnames=list(c("positive","negative"), unique_labels))

# Loop through each label
for (i in unique_labels) {
  
  # Subset edges with the current label
  module2_all_degree <- all_network_edge[all_network_edge$Source == i | all_network_edge$Target == i,] 
  
  # Subset edges with positive correlation
  module2_all_positive_degree <- module2_all_degree[module2_all_degree$e_cor >= 0,]
  
  
  # Subset edges with negative correlation
  module2_all_negative_degree <- module2_all_degree[module2_all_degree$e_cor <= 0,]
  
  # Create matrix with positive and negative degrees
  module2_all_i <- matrix(c(dim(module2_all_positive_degree)[1], dim(module2_all_negative_degree)[1]), nrow=2, ncol=1, byrow=TRUE, dimnames=list(c("positive","negative"), i))
  module2_all_i
  # Fill in matrix with values from current label
  module2_all_matrix[, i] <- module2_all_i[, 1]
}

print(module2_all_matrix)
colnames(module2_all_matrix) <- rownames(module2_gephi2)
module2_all_matrix

dim(module2_all_matrix)
write.table(module2_all_matrix,paste0(output,"module2_all_matrix.txt"),sep="\t",quote=F)


##方法二：基于网络重新计算得到的数据

module2_network_edge <- read.table("module2_edge_list2.txt", row.names=1,sep="\t", header=T, blank.lines.skip=F, check.names=F)
#读取数据的时候，header=T表示保留原始表格的表头.row.names=1表示保留数据表格的第一列，=2表示去除表格的第一列
module2_network_edge


# Initialize empty matrix
module2_all_matrix2 <- matrix(0, nrow=2, ncol=length(unique_labels), dimnames=list(c("positive","negative"), unique_labels))

# Loop through each label
for (i in unique_labels) {
  
  # Subset edges with the current label
  module2_all_degree <- module2_network_edge[module2_network_edge$Source == i | module2_network_edge$Target == i,] 
  
  # Subset edges with positive correlation
  module2_all_positive_degree <- module2_all_degree[module2_all_degree$e_cor >= 0,]
  
  
  # Subset edges with negative correlation
  module2_all_negative_degree <- module2_all_degree[module2_all_degree$e_cor <= 0,]
  
  # Create matrix with positive and negative degrees
  module2_all_i <- matrix(c(dim(module2_all_positive_degree)[1], dim(module2_all_negative_degree)[1]), nrow=2, ncol=1, byrow=TRUE, dimnames=list(c("positive","negative"), i))
  module2_all_i
  # Fill in matrix with values from current label
  module2_all_matrix2[, i] <- module2_all_i[, 1]
}

print(module2_all_matrix2)
colnames(module2_all_matrix2) <- rownames(module2_gephi2)
module2_all_matrix2

dim(module2_all_matrix2)

write.table(module2_all_matrix2,paste0(output,"module2_all_matrix2.txt"),sep="\t",quote=F)






####module3
#提取数据
module3_hub <- module3_gephi2[which(module3_gephi2$degree >= 60 & module3_gephi2$closnesscentrality >= 0.3),] 
module3_hub
dim(module3_hub) #module2中含有35个核心微生物存在
write.table(module3_hub,paste0(output,"module3_hub.txt"),sep="\t",quote=F)

#从大网络边数据列表中提取出核心微生物所对应的边数据
all_network_edge <- read.table("all_network_edge_list-4-2.txt", row.names=1,sep="\t", header=T, blank.lines.skip=F, check.names=F)
#读取数据的时候，header=T表示保留原始表格的表头.row.names=1表示保留数据表格的第一列，=2表示去除表格的第一列
all_network_edge

####循环的编写####
# Get unique labels
unique_labels <- unique(module3_hub$Label)
unique_labels
# Initialize empty matrix
module3_hub_matrix <- matrix(0, nrow=2, ncol=length(unique_labels), dimnames=list(c("positive","negative"), unique_labels))

# Loop through each label
for (i in unique_labels) {
  
  # Subset edges with the current label
  module3_hub_degree <- all_network_edge[all_network_edge$Source == i | all_network_edge$Target == i,] 
  
  # Subset edges with positive correlation
  module3_hub_positive_degree <- module3_hub_degree[module3_hub_degree$e_cor >= 0,]
  
  
  # Subset edges with negative correlation
  module3_hub_negative_degree <- module3_hub_degree[module3_hub_degree$e_cor <= 0,]
  
  # Create matrix with positive and negative degrees
  module3_hub_i <- matrix(c(dim(module3_hub_positive_degree)[1], dim(module3_hub_negative_degree)[1]), nrow=2, ncol=1, byrow=TRUE, dimnames=list(c("positive","negative"), i))
  module3_hub_i
  # Fill in matrix with values from current label
  module3_hub_matrix[, i] <- module3_hub_i[, 1]
}

print(module3_hub_matrix)
colnames(module3_hub_matrix) <- rownames(module3_hub)
module3_hub_matrix
write.table(module3_hub_matrix,paste0(output,"module3_hub_matrix.txt"),sep="\t",quote=F)


#数据处理一下，以适合用barplot绘制柱状图
module3_hub_matrix <- as.data.frame(t(module3_hub_matrix))
module3_hub_matrix
sum <- apply(module3_hub_matrix,1,sum) 
module3_hub_matrix1 <- cbind(module3_hub_matrix,sum)
module3_hub_matrix1 <- as.data.frame(module3_hub_matrix1)
module3_hub_matrix1 <- module3_hub_matrix1[order(module3_hub_matrix1[,"sum"],decreasing = T),]
module3_hub_matrix1

module3_hub_matrix1 <- subset(module3_hub_matrix1, select = -sum)
module3_hub_matrix1 <- t(module3_hub_matrix1)
module3_hub_matrix1
#如果需要增减物种的数目，请调整此处的数值
module3_hub_matrix2 <- module3_hub_matrix1[,1:10]
module3_hub_matrix2


##module3核心微生物出图
pdf(paste0(output,"module3_top10_hubtaxa.pdf"), width=7, height=7/3)
par(mfrow=c(1,3), mar=c(7,3,2,0.5))
mycol<- c("#D55E86","#63D07F") 


p1 <- barplot(module3_hub_matrix2,
              las=2, border=NA, axes=T, cex.names=.5,ylim=c(0,130),main="Module3",
              col=mycol[1:nrow(module3_hub_matrix2)],ylab="Degree ")

## legend 
plot.new()
legend("top", bty="n", cex=0.5, #x.intersp=0.1, #y.intersp=0.75,
       legend=rev(rownames(module3_hub_matrix2)), 
       fill=rev(mycol), 
       border=rev(mycol) )
dev.off()






####module3正负边/度####

##方法一：根据从总网络的提取结果
#184
unique_labels <- unique(module3_gephi2$Label)
unique_labels
# Initialize empty matrix
module3_all_matrix <- matrix(0, nrow=2, ncol=length(unique_labels), dimnames=list(c("positive","negative"), unique_labels))

# Loop through each label
for (i in unique_labels) {
  
  # Subset edges with the current label
  module3_all_degree <- all_network_edge[all_network_edge$Source == i | all_network_edge$Target == i,] 
  
  # Subset edges with positive correlation
  module3_all_positive_degree <- module3_all_degree[module3_all_degree$e_cor >= 0,]
  
  
  # Subset edges with negative correlation
  module3_all_negative_degree <- module3_all_degree[module3_all_degree$e_cor <= 0,]
  
  # Create matrix with positive and negative degrees
  module3_all_i <- matrix(c(dim(module3_all_positive_degree)[1], dim(module3_all_negative_degree)[1]), nrow=2, ncol=1, byrow=TRUE, dimnames=list(c("positive","negative"), i))
  module3_all_i
  # Fill in matrix with values from current label
  module3_all_matrix[, i] <- module3_all_i[, 1]
}

print(module3_all_matrix)
colnames(module3_all_matrix) <- rownames(module3_gephi2)
module3_all_matrix

dim(module3_all_matrix)
write.table(module3_all_matrix,paste0(output,"module3_all_matrix.txt"),sep="\t",quote=F)


##方法二：基于网络重新计算得到的数据

module3_network_edge <- read.table("module3_edge_list2.txt", row.names=1,sep="\t", header=T, blank.lines.skip=F, check.names=F)
#读取数据的时候，header=T表示保留原始表格的表头.row.names=1表示保留数据表格的第一列，=2表示去除表格的第一列
module3_network_edge


# Initialize empty matrix
module3_all_matrix2 <- matrix(0, nrow=2, ncol=length(unique_labels), dimnames=list(c("positive","negative"), unique_labels))

# Loop through each label
for (i in unique_labels) {
  
  # Subset edges with the current label
  module3_all_degree <- module3_network_edge[module3_network_edge$Source == i | module3_network_edge$Target == i,] 
  
  # Subset edges with positive correlation
  module3_all_positive_degree <- module3_all_degree[module3_all_degree$e_cor >= 0,]
  
  
  # Subset edges with negative correlation
  module3_all_negative_degree <- module3_all_degree[module3_all_degree$e_cor <= 0,]
  
  # Create matrix with positive and negative degrees
  module3_all_i <- matrix(c(dim(module3_all_positive_degree)[1], dim(module3_all_negative_degree)[1]), nrow=2, ncol=1, byrow=TRUE, dimnames=list(c("positive","negative"), i))
  module3_all_i
  # Fill in matrix with values from current label
  module3_all_matrix2[, i] <- module3_all_i[, 1]
}

print(module3_all_matrix2)
colnames(module3_all_matrix2) <- rownames(module3_gephi2)
module3_all_matrix2

dim(module3_all_matrix2)

write.table(module3_all_matrix2,paste0(output,"module3_all_matrix2.txt"),sep="\t",quote=F)




####module4
#提取数据
module4_hub <- module4_gephi2[which(module4_gephi2$degree >= 60 & module4_gephi2$closnesscentrality >= 0.3),] 
module4_hub
dim(module4_hub) #module4中没有核心微生物存在，以下代码不运行
write.table(module4_hub,paste0(output,"module4_hub.txt"),sep="\t",quote=F)

#从大网络边数据列表中提取出核心微生物所对应的边数据
all_network_edge <- read.table("all_network_edge_list-4-2.txt", row.names=1,sep="\t", header=T, blank.lines.skip=F, check.names=F)
#读取数据的时候，header=T表示保留原始表格的表头.row.names=1表示保留数据表格的第一列，=2表示去除表格的第一列
all_network_edge

####循环的编写####
# Get unique labels
unique_labels <- unique(module4_hub$Label)
unique_labels
# Initialize empty matrix
module4_hub_matrix <- matrix(0, nrow=2, ncol=length(unique_labels), dimnames=list(c("positive","negative"), unique_labels))

# Loop through each label
for (i in unique_labels) {
  
  # Subset edges with the current label
  module4_hub_degree <- all_network_edge[all_network_edge$Source == i | all_network_edge$Target == i,] 
  
  # Subset edges with positive correlation
  module4_hub_positive_degree <- module4_hub_degree[module4_hub_degree$e_cor >= 0,]
  
  
  # Subset edges with negative correlation
  module4_hub_negative_degree <- module4_hub_degree[module4_hub_degree$e_cor <= 0,]
  
  # Create matrix with positive and negative degrees
  module4_hub_i <- matrix(c(dim(module4_hub_positive_degree)[1], dim(module4_hub_negative_degree)[1]), nrow=2, ncol=1, byrow=TRUE, dimnames=list(c("positive","negative"), i))
  module4_hub_i
  # Fill in matrix with values from current label
  module4_hub_matrix[, i] <- module4_hub_i[, 1]
}

print(module4_hub_matrix)
colnames(module4_hub_matrix) <- rownames(module4_hub)
module4_hub_matrix
write.table(module4_hub_matrix,paste0(output,"module4_hub_matrix.txt"),sep="\t",quote=F)


#数据处理一下，以适合用barplot绘制柱状图
module4_hub_matrix <- as.data.frame(t(module4_hub_matrix))
module4_hub_matrix
sum <- apply(module4_hub_matrix,1,sum) 
module4_hub_matrix1 <- cbind(module4_hub_matrix,sum)
module4_hub_matrix1 <- as.data.frame(module4_hub_matrix1)
module4_hub_matrix1 <- module4_hub_matrix1[order(module4_hub_matrix1[,"sum"],decreasing = T),]
module4_hub_matrix1

module4_hub_matrix1 <- subset(module4_hub_matrix1, select = -sum)
module4_hub_matrix1 <- t(module4_hub_matrix1)
module4_hub_matrix1
#如果需要增减物种的数目，请调整此处的数值
module4_hub_matrix2 <- module4_hub_matrix1[,1:10]
module4_hub_matrix2


##module4核心微生物出图
pdf(paste0(output,"module4_top10_hubtaxa.pdf"), width=7, height=7/3)
par(mfrow=c(1,3), mar=c(5,3,2,0.5))
mycol<- c("#D55E86","#63D07F") 


p1 <- barplot(module4_hub_matrix2,
              las=2, border=NA, axes=T, cex.names=.5,ylim=c(0,130),main="Module4",
              col=mycol[1:nrow(module4_hub_matrix2)],ylab="Degree ")

## legend 
plot.new()
legend("top", bty="n", cex=0.5, #x.intersp=0.1, #y.intersp=0.75,
       legend=rev(rownames(module4_hub_matrix2)), 
       fill=rev(mycol), 
       border=rev(mycol) )
dev.off()



####module4正负边/度####

##方法一：根据从总网络的提取结果
#89
unique_labels <- unique(module4_gephi2$Label)
unique_labels
# Initialize empty matrix
module4_all_matrix <- matrix(0, nrow=2, ncol=length(unique_labels), dimnames=list(c("positive","negative"), unique_labels))

# Loop through each label
for (i in unique_labels) {
  
  # Subset edges with the current label
  module4_all_degree <- all_network_edge[all_network_edge$Source == i | all_network_edge$Target == i,] 
  
  # Subset edges with positive correlation
  module4_all_positive_degree <- module4_all_degree[module4_all_degree$e_cor >= 0,]
  
  
  # Subset edges with negative correlation
  module4_all_negative_degree <- module4_all_degree[module4_all_degree$e_cor <= 0,]
  
  # Create matrix with positive and negative degrees
  module4_all_i <- matrix(c(dim(module4_all_positive_degree)[1], dim(module4_all_negative_degree)[1]), nrow=2, ncol=1, byrow=TRUE, dimnames=list(c("positive","negative"), i))
  module4_all_i
  # Fill in matrix with values from current label
  module4_all_matrix[, i] <- module4_all_i[, 1]
}

print(module4_all_matrix)
colnames(module4_all_matrix) <- rownames(module4_gephi2)
module4_all_matrix

dim(module4_all_matrix)
write.table(module4_all_matrix,paste0(output,"module4_all_matrix.txt"),sep="\t",quote=F)


##方法二：基于网络重新计算得到的数据

module4_network_edge <- read.table("module4_edge_list2.txt", row.names=1,sep="\t", header=T, blank.lines.skip=F, check.names=F)
#读取数据的时候，header=T表示保留原始表格的表头.row.names=1表示保留数据表格的第一列，=2表示去除表格的第一列
module4_network_edge


# Initialize empty matrix
module4_all_matrix2 <- matrix(0, nrow=2, ncol=length(unique_labels), dimnames=list(c("positive","negative"), unique_labels))

# Loop through each label
for (i in unique_labels) {
  
  # Subset edges with the current label
  module4_all_degree <- module4_network_edge[module4_network_edge$Source == i | module4_network_edge$Target == i,] 
  
  # Subset edges with positive correlation
  module4_all_positive_degree <- module4_all_degree[module4_all_degree$e_cor >= 0,]
  
  
  # Subset edges with negative correlation
  module4_all_negative_degree <- module4_all_degree[module4_all_degree$e_cor <= 0,]
  
  # Create matrix with positive and negative degrees
  module4_all_i <- matrix(c(dim(module4_all_positive_degree)[1], dim(module4_all_negative_degree)[1]), nrow=2, ncol=1, byrow=TRUE, dimnames=list(c("positive","negative"), i))
  module4_all_i
  # Fill in matrix with values from current label
  module4_all_matrix2[, i] <- module4_all_i[, 1]
}

print(module4_all_matrix2)
colnames(module4_all_matrix2) <- rownames(module4_gephi2)
module4_all_matrix2

dim(module4_all_matrix2)

write.table(module4_all_matrix2,paste0(output,"module4_all_matrix2.txt"),sep="\t",quote=F)






####module5####
#提取数据
module5_hub <- module5_gephi2[which(module5_gephi2$degree >= 60 & module5_gephi2$closnesscentrality >= 0.28),] 
module5_hub
dim(module5_hub) #module5中4个核心微生物
write.table(module5_hub,paste0(output,"module5_hub.txt"),sep="\t",quote=F)

#从大网络边数据列表中提取出核心微生物所对应的边数据
all_network_edge <- read.table("all_network_edge_list-4-2.txt", row.names=1,sep="\t", header=T, blank.lines.skip=F, check.names=F)
#读取数据的时候，header=T表示保留原始表格的表头.row.names=1表示保留数据表格的第一列，=2表示去除表格的第一列
all_network_edge

####循环的编写####
# Get unique labels
unique_labels <- unique(module5_hub$Label)
unique_labels
# Initialize empty matrix
module5_hub_matrix <- matrix(0, nrow=2, ncol=length(unique_labels), dimnames=list(c("positive","negative"), unique_labels))

# Loop through each label
for (i in unique_labels) {
  
  # Subset edges with the current label
  module5_hub_degree <- all_network_edge[all_network_edge$Source == i | all_network_edge$Target == i,] 
  
  # Subset edges with positive correlation
  module5_hub_positive_degree <- module5_hub_degree[module5_hub_degree$e_cor >= 0,]
  
  
  # Subset edges with negative correlation
  module5_hub_negative_degree <- module5_hub_degree[module5_hub_degree$e_cor <= 0,]
  
  # Create matrix with positive and negative degrees
  module5_hub_i <- matrix(c(dim(module5_hub_positive_degree)[1], dim(module5_hub_negative_degree)[1]), nrow=2, ncol=1, byrow=TRUE, dimnames=list(c("positive","negative"), i))
  module5_hub_i
  # Fill in matrix with values from current label
  module5_hub_matrix[, i] <- module5_hub_i[, 1]
}

print(module5_hub_matrix)
colnames(module5_hub_matrix) <- rownames(module5_hub)
module5_hub_matrix
write.table(module5_hub_matrix,paste0(output,"module5_hub_matrix.txt"),sep="\t",quote=F)


#数据处理一下，以适合用barplot绘制柱状图
module5_hub_matrix <- as.data.frame(t(module5_hub_matrix))
module5_hub_matrix
sum <- apply(module5_hub_matrix,1,sum) 
module5_hub_matrix1 <- cbind(module5_hub_matrix,sum)
module5_hub_matrix1 <- as.data.frame(module5_hub_matrix1)
module5_hub_matrix1 <- module5_hub_matrix1[order(module5_hub_matrix1[,"sum"],decreasing = T),]
module5_hub_matrix1

module5_hub_matrix1 <- subset(module5_hub_matrix1, select = -sum)
module5_hub_matrix1 <- t(module5_hub_matrix1)
module5_hub_matrix1
#如果需要增减物种的数目，请调整此处的数值，由于核心微生物的数量为4，故无需进行top10
module5_hub_matrix2 <- module5_hub_matrix1[,1:10]
module5_hub_matrix2


##特别需要注意的是，由于该模块核心微生物为4个，少于10个，所以人工补全了表格进行出图，完了在AI中删除即可
module5_hub_matrix3 <- read.table("module5_hub_matrix2.txt", row.names=1,sep="\t", header=T, blank.lines.skip=F, check.names=F)
module5_hub_matrix3

#数据处理一下，以适合用barplot绘制柱状图
module5_hub_matrix3 <-module5_hub_matrix3[,1:ncol(module5_hub_matrix3)]
rownames(module5_hub_matrix3)<-module5_hub_matrix3[,1]
module5_hub_matrix3<-as.matrix(module5_hub_matrix3)
module5_hub_matrix3
colnames(module5_hub_matrix3)



##module5核心微生物出图
pdf(paste0(output,"module5_top10_hubtaxa_legend.pdf"), width=7, height=7/3)
par(mfrow=c(1,3), mar=c(7,3,2,0.5))
mycol<- c("#D55E86","#63D07F") 


p1 <- barplot(module5_hub_matrix3,
              las=2, border=NA, axes=T, cex.names=.5,ylim=c(0,130),main="Module5",
              col=mycol[1:nrow(module5_hub_matrix3)],ylab="Degree ")

## legend 
plot.new()
legend("top", bty="n", cex=0.5, #x.intersp=0.1, #y.intersp=0.75,
       legend=rev(rownames(module5_hub_matrix3)), 
       fill=rev(mycol), 
       border=rev(mycol) )
dev.off()



####module5正负边/度####

##方法一：根据从总网络的提取结果
#206
unique_labels <- unique(module5_gephi2$Label)
unique_labels
# Initialize empty matrix
module5_all_matrix <- matrix(0, nrow=2, ncol=length(unique_labels), dimnames=list(c("positive","negative"), unique_labels))

# Loop through each label
for (i in unique_labels) {
  
  # Subset edges with the current label
  module5_all_degree <- all_network_edge[all_network_edge$Source == i | all_network_edge$Target == i,] 
  
  # Subset edges with positive correlation
  module5_all_positive_degree <- module5_all_degree[module5_all_degree$e_cor >= 0,]
  
  
  # Subset edges with negative correlation
  module5_all_negative_degree <- module5_all_degree[module5_all_degree$e_cor <= 0,]
  
  # Create matrix with positive and negative degrees
  module5_all_i <- matrix(c(dim(module5_all_positive_degree)[1], dim(module5_all_negative_degree)[1]), nrow=2, ncol=1, byrow=TRUE, dimnames=list(c("positive","negative"), i))
  module5_all_i
  # Fill in matrix with values from current label
  module5_all_matrix[, i] <- module5_all_i[, 1]
}

print(module5_all_matrix)
colnames(module5_all_matrix) <- rownames(module5_gephi2)
module5_all_matrix

dim(module5_all_matrix)
write.table(module5_all_matrix,paste0(output,"module5_all_matrix.txt"),sep="\t",quote=F)


##方法二：基于网络重新计算得到的数据

module5_network_edge <- read.table("module5_edge_list2.txt", row.names=1,sep="\t", header=T, blank.lines.skip=F, check.names=F)
#读取数据的时候，header=T表示保留原始表格的表头.row.names=1表示保留数据表格的第一列，=2表示去除表格的第一列
module5_network_edge


# Initialize empty matrix
module5_all_matrix2 <- matrix(0, nrow=2, ncol=length(unique_labels), dimnames=list(c("positive","negative"), unique_labels))

# Loop through each label
for (i in unique_labels) {
  
  # Subset edges with the current label
  module5_all_degree <- module5_network_edge[module5_network_edge$Source == i | module5_network_edge$Target == i,] 
  
  # Subset edges with positive correlation
  module5_all_positive_degree <- module5_all_degree[module5_all_degree$e_cor >= 0,]
  
  
  # Subset edges with negative correlation
  module5_all_negative_degree <- module5_all_degree[module5_all_degree$e_cor <= 0,]
  
  # Create matrix with positive and negative degrees
  module5_all_i <- matrix(c(dim(module5_all_positive_degree)[1], dim(module5_all_negative_degree)[1]), nrow=2, ncol=1, byrow=TRUE, dimnames=list(c("positive","negative"), i))
  module5_all_i
  # Fill in matrix with values from current label
  module5_all_matrix2[, i] <- module5_all_i[, 1]
}

print(module5_all_matrix2)
colnames(module5_all_matrix2) <- rownames(module5_gephi2)
module5_all_matrix2

dim(module5_all_matrix2)

write.table(module5_all_matrix2,paste0(output,"module5_all_matrix2.txt"),sep="\t",quote=F)




####module6####
#提取数据
module6_hub <- module6_gephi2[which(module6_gephi2$degree >= 60 & module6_gephi2$closnesscentrality >= 0.3),] 
module6_hub
dim(module6_hub) #module6中57个核心微生物
write.table(module6_hub,paste0(output,"module6_hub.txt"),sep="\t",quote=F)

#从大网络边数据列表中提取出核心微生物所对应的边数据
all_network_edge <- read.table("all_network_edge_list-4-2.txt", row.names=1,sep="\t", header=T, blank.lines.skip=F, check.names=F)
#读取数据的时候，header=T表示保留原始表格的表头.row.names=1表示保留数据表格的第一列，=2表示去除表格的第一列
all_network_edge

####循环的编写####
# Get unique labels
unique_labels <- unique(module6_hub$Label)
unique_labels
# Initialize empty matrix
module6_hub_matrix <- matrix(0, nrow=2, ncol=length(unique_labels), dimnames=list(c("positive","negative"), unique_labels))

# Loop through each label
for (i in unique_labels) {
  
  # Subset edges with the current label
  module6_hub_degree <- all_network_edge[all_network_edge$Source == i | all_network_edge$Target == i,] 
  
  # Subset edges with positive correlation
  module6_hub_positive_degree <- module6_hub_degree[module6_hub_degree$e_cor >= 0,]
  
  
  # Subset edges with negative correlation
  module6_hub_negative_degree <- module6_hub_degree[module6_hub_degree$e_cor <= 0,]
  
  # Create matrix with positive and negative degrees
  module6_hub_i <- matrix(c(dim(module6_hub_positive_degree)[1], dim(module6_hub_negative_degree)[1]), nrow=2, ncol=1, byrow=TRUE, dimnames=list(c("positive","negative"), i))
  module6_hub_i
  # Fill in matrix with values from current label
  module6_hub_matrix[, i] <- module6_hub_i[, 1]
}

print(module6_hub_matrix)
colnames(module6_hub_matrix) <- rownames(module6_hub)
module6_hub_matrix
write.table(module6_hub_matrix,paste0(output,"module6_hub_matrix.txt"),sep="\t",quote=F)


#数据处理一下，以适合用barplot绘制柱状图
module6_hub_matrix <- as.data.frame(t(module6_hub_matrix))
module6_hub_matrix
sum <- apply(module6_hub_matrix,1,sum) 
module6_hub_matrix1 <- cbind(module6_hub_matrix,sum)
module6_hub_matrix1 <- as.data.frame(module6_hub_matrix1)
module6_hub_matrix1 <- module6_hub_matrix1[order(module6_hub_matrix1[,"sum"],decreasing = T),]
module6_hub_matrix1

module6_hub_matrix1 <- subset(module6_hub_matrix1, select = -sum)
module6_hub_matrix1 <- t(module6_hub_matrix1)
module6_hub_matrix1
#如果需要增减物种的数目，请调整此处的数值，由于核心微生物的数量为4，故无需进行top10
module6_hub_matrix2 <- module6_hub_matrix1[,1:10]
module6_hub_matrix2


##module6核心微生物出图
pdf(paste0(output,"module6_top10_hubtaxa.pdf"), width=7, height=7/3)
par(mfrow=c(1,3), mar=c(7,3,2,0.5))
mycol<- c("#D55E86","#63D07F") 


p1 <- barplot(module6_hub_matrix2,
              las=2, border=NA, axes=T, cex.names=.5,ylim=c(0,130),main="Module6",
              col=mycol[1:nrow(module6_hub_matrix2)],ylab="Degree ")

## legend 
plot.new()
legend("top", bty="n", cex=0.5, #x.intersp=0.1, #y.intersp=0.75,
       legend=rev(rownames(module6_hub_matrix2)), 
       fill=rev(mycol), 
       border=rev(mycol) )
dev.off()



####module6正负边/度####

##方法一：根据从总网络的提取结果
#196
unique_labels <- unique(module6_gephi2$Label)
unique_labels
# Initialize empty matrix
module6_all_matrix <- matrix(0, nrow=2, ncol=length(unique_labels), dimnames=list(c("positive","negative"), unique_labels))

# Loop through each label
for (i in unique_labels) {
  
  # Subset edges with the current label
  module6_all_degree <- all_network_edge[all_network_edge$Source == i | all_network_edge$Target == i,] 
  
  # Subset edges with positive correlation
  module6_all_positive_degree <- module6_all_degree[module6_all_degree$e_cor >= 0,]
  
  
  # Subset edges with negative correlation
  module6_all_negative_degree <- module6_all_degree[module6_all_degree$e_cor <= 0,]
  
  # Create matrix with positive and negative degrees
  module6_all_i <- matrix(c(dim(module6_all_positive_degree)[1], dim(module6_all_negative_degree)[1]), nrow=2, ncol=1, byrow=TRUE, dimnames=list(c("positive","negative"), i))
  module6_all_i
  # Fill in matrix with values from current label
  module6_all_matrix[, i] <- module6_all_i[, 1]
}

print(module6_all_matrix)
colnames(module6_all_matrix) <- rownames(module6_gephi2)
module6_all_matrix

dim(module6_all_matrix)
write.table(module6_all_matrix,paste0(output,"module6_all_matrix.txt"),sep="\t",quote=F)


##方法二：基于网络重新计算得到的数据

module6_network_edge <- read.table("module6_edge_list2.txt", row.names=1,sep="\t", header=T, blank.lines.skip=F, check.names=F)
#读取数据的时候，header=T表示保留原始表格的表头.row.names=1表示保留数据表格的第一列，=2表示去除表格的第一列
module6_network_edge


# Initialize empty matrix
module6_all_matrix2 <- matrix(0, nrow=2, ncol=length(unique_labels), dimnames=list(c("positive","negative"), unique_labels))

# Loop through each label
for (i in unique_labels) {
  
  # Subset edges with the current label
  module6_all_degree <- module6_network_edge[module6_network_edge$Source == i | module6_network_edge$Target == i,] 
  
  # Subset edges with positive correlation
  module6_all_positive_degree <- module6_all_degree[module6_all_degree$e_cor >= 0,]
  
  
  # Subset edges with negative correlation
  module6_all_negative_degree <- module6_all_degree[module6_all_degree$e_cor <= 0,]
  
  # Create matrix with positive and negative degrees
  module6_all_i <- matrix(c(dim(module6_all_positive_degree)[1], dim(module6_all_negative_degree)[1]), nrow=2, ncol=1, byrow=TRUE, dimnames=list(c("positive","negative"), i))
  module6_all_i
  # Fill in matrix with values from current label
  module6_all_matrix2[, i] <- module6_all_i[, 1]
}

print(module6_all_matrix2)
colnames(module6_all_matrix2) <- rownames(module6_gephi2)
module6_all_matrix2

dim(module6_all_matrix2)

write.table(module6_all_matrix2,paste0(output,"module6_all_matrix2.txt"),sep="\t",quote=F)






####module7####
#提取数据
module7_hub <- module7_gephi2[which(module7_gephi2$degree >= 60 & module7_gephi2$closnesscentrality >= 0.3),] 
module7_hub
dim(module7_hub) #module7中0个核心微生物
write.table(module7_hub,paste0(output,"module7_hub.txt"),sep="\t",quote=F)

#从大网络边数据列表中提取出核心微生物所对应的边数据
all_network_edge <- read.table("all_network_edge_list-4-2.txt", row.names=1,sep="\t", header=T, blank.lines.skip=F, check.names=F)
#读取数据的时候，header=T表示保留原始表格的表头.row.names=1表示保留数据表格的第一列，=2表示去除表格的第一列
all_network_edge

####循环的编写####
# Get unique labels
unique_labels <- unique(module7_hub$Label)
unique_labels
# Initialize empty matrix
module7_hub_matrix <- matrix(0, nrow=2, ncol=length(unique_labels), dimnames=list(c("positive","negative"), unique_labels))

# Loop through each label
for (i in unique_labels) {
  
  # Subset edges with the current label
  module7_hub_degree <- all_network_edge[all_network_edge$Source == i | all_network_edge$Target == i,] 
  
  # Subset edges with positive correlation
  module7_hub_positive_degree <- module7_hub_degree[module7_hub_degree$e_cor >= 0,]
  
  
  # Subset edges with negative correlation
  module7_hub_negative_degree <- module7_hub_degree[module7_hub_degree$e_cor <= 0,]
  
  # Create matrix with positive and negative degrees
  module7_hub_i <- matrix(c(dim(module7_hub_positive_degree)[1], dim(module7_hub_negative_degree)[1]), nrow=2, ncol=1, byrow=TRUE, dimnames=list(c("positive","negative"), i))
  module7_hub_i
  # Fill in matrix with values from current label
  module7_hub_matrix[, i] <- module6_hub_i[, 1]
}

print(module7_hub_matrix)
colnames(module7_hub_matrix) <- rownames(module7_hub)
module7_hub_matrix
write.table(module7_hub_matrix,paste0(output,"module7_hub_matrix.txt"),sep="\t",quote=F)


#数据处理一下，以适合用barplot绘制柱状图
module7_hub_matrix <- as.data.frame(t(module7_hub_matrix))
module7_hub_matrix
sum <- apply(module7_hub_matrix,1,sum) 
module7_hub_matrix1 <- cbind(module7_hub_matrix,sum)
module7_hub_matrix1 <- as.data.frame(module7_hub_matrix1)
module7_hub_matrix1 <- module7_hub_matrix1[order(module7_hub_matrix1[,"sum"],decreasing = T),]
module7_hub_matrix1

module7_hub_matrix1 <- subset(module7_hub_matrix1, select = -sum)
module7_hub_matrix1 <- t(module7_hub_matrix1)
module7_hub_matrix1
#如果需要增减物种的数目，请调整此处的数值，由于核心微生物的数量为4，故无需进行top10
#module5_hub_matrix2 <- module5_hub_matrix1[,1:10]
#module5_hub_matrix2


##module7核心微生物出图
pdf(paste0(output,"module7_top10_hubtaxa.pdf"), width=7, height=7/3)
par(mfrow=c(1,3), mar=c(5,3,2,0.5))
mycol<- c("#D55E86","#63D07F") 


p1 <- barplot(module7_hub_matrix1,
              las=2, border=NA, axes=T, cex.names=.5,ylim=c(0,130),main="Module7",
              col=mycol[1:nrow(module7_hub_matrix1)],ylab="Degree ")

## legend 
plot.new()
legend("top", bty="n", cex=0.5, #x.intersp=0.1, #y.intersp=0.75,
       legend=rev(rownames(module7_hub_matrix1)), 
       fill=rev(mycol), 
       border=rev(mycol) )
dev.off()






####module7正负边/度####

##方法一：根据从总网络的提取结果
#122
unique_labels <- unique(module7_gephi2$Label)
unique_labels
# Initialize empty matrix
module7_all_matrix <- matrix(0, nrow=2, ncol=length(unique_labels), dimnames=list(c("positive","negative"), unique_labels))

# Loop through each label
for (i in unique_labels) {
  
  # Subset edges with the current label
  module7_all_degree <- all_network_edge[all_network_edge$Source == i | all_network_edge$Target == i,] 
  
  # Subset edges with positive correlation
  module7_all_positive_degree <- module7_all_degree[module7_all_degree$e_cor >= 0,]
  
  
  # Subset edges with negative correlation
  module7_all_negative_degree <- module7_all_degree[module7_all_degree$e_cor <= 0,]
  
  # Create matrix with positive and negative degrees
  module7_all_i <- matrix(c(dim(module7_all_positive_degree)[1], dim(module7_all_negative_degree)[1]), nrow=2, ncol=1, byrow=TRUE, dimnames=list(c("positive","negative"), i))
  module7_all_i
  # Fill in matrix with values from current label
  module7_all_matrix[, i] <- module7_all_i[, 1]
}

print(module7_all_matrix)
colnames(module7_all_matrix) <- rownames(module7_gephi2)
module7_all_matrix

dim(module7_all_matrix)
write.table(module7_all_matrix,paste0(output,"module7_all_matrix.txt"),sep="\t",quote=F)


##方法二：基于网络重新计算得到的数据

module7_network_edge <- read.table("module7_edge_list2.txt", row.names=1,sep="\t", header=T, blank.lines.skip=F, check.names=F)
#读取数据的时候，header=T表示保留原始表格的表头.row.names=1表示保留数据表格的第一列，=2表示去除表格的第一列
module7_network_edge


# Initialize empty matrix
module7_all_matrix2 <- matrix(0, nrow=2, ncol=length(unique_labels), dimnames=list(c("positive","negative"), unique_labels))

# Loop through each label
for (i in unique_labels) {
  
  # Subset edges with the current label
  module7_all_degree <- module7_network_edge[module7_network_edge$Source == i | module7_network_edge$Target == i,] 
  
  # Subset edges with positive correlation
  module7_all_positive_degree <- module7_all_degree[module7_all_degree$e_cor >= 0,]
  
  
  # Subset edges with negative correlation
  module7_all_negative_degree <- module7_all_degree[module7_all_degree$e_cor <= 0,]
  
  # Create matrix with positive and negative degrees
  module7_all_i <- matrix(c(dim(module7_all_positive_degree)[1], dim(module7_all_negative_degree)[1]), nrow=2, ncol=1, byrow=TRUE, dimnames=list(c("positive","negative"), i))
  module7_all_i
  # Fill in matrix with values from current label
  module7_all_matrix2[, i] <- module7_all_i[, 1]
}

print(module7_all_matrix2)
colnames(module7_all_matrix2) <- rownames(module7_gephi2)
module7_all_matrix2

dim(module7_all_matrix2)

write.table(module7_all_matrix2,paste0(output,"module7_all_matrix2.txt"),sep="\t",quote=F)






####module8####
#提取数据
module8_hub <- module8_gephi2[which(module8_gephi2$degree >= 60 & module8_gephi2$closnesscentrality >= 0.3),] 
module8_hub
dim(module8_hub) #module8中0个核心微生物
write.table(module8_hub,paste0(output,"module8_hub.txt"),sep="\t",quote=F)

#从大网络边数据列表中提取出核心微生物所对应的边数据
all_network_edge <- read.table("all_network_edge_list-4-2.txt", row.names=1,sep="\t", header=T, blank.lines.skip=F, check.names=F)
#读取数据的时候，header=T表示保留原始表格的表头.row.names=1表示保留数据表格的第一列，=2表示去除表格的第一列
all_network_edge

####循环的编写####
# Get unique labels
unique_labels <- unique(module8_hub$Label)
unique_labels
# Initialize empty matrix
module8_hub_matrix <- matrix(0, nrow=2, ncol=length(unique_labels), dimnames=list(c("positive","negative"), unique_labels))

# Loop through each label
for (i in unique_labels) {
  
  # Subset edges with the current label
  module8_hub_degree <- all_network_edge[all_network_edge$Source == i | all_network_edge$Target == i,] 
  
  # Subset edges with positive correlation
  module8_hub_positive_degree <- module8_hub_degree[module8_hub_degree$e_cor >= 0,]
  
  
  # Subset edges with negative correlation
  module8_hub_negative_degree <- module8_hub_degree[module8_hub_degree$e_cor <= 0,]
  
  # Create matrix with positive and negative degrees
  module8_hub_i <- matrix(c(dim(module8_hub_positive_degree)[1], dim(module8_hub_negative_degree)[1]), nrow=2, ncol=1, byrow=TRUE, dimnames=list(c("positive","negative"), i))
  module8_hub_i
  # Fill in matrix with values from current label
  module8_hub_matrix[, i] <- module8_hub_i[, 1]
}

print(module8_hub_matrix)
colnames(module8_hub_matrix) <- rownames(module8_hub)
module8_hub_matrix
write.table(module8_hub_matrix,paste0(output,"module8_hub_matrix.txt"),sep="\t",quote=F)


#数据处理一下，以适合用barplot绘制柱状图
module8_hub_matrix <- as.data.frame(t(module8_hub_matrix))
module8_hub_matrix
sum <- apply(module8_hub_matrix,1,sum) 
module8_hub_matrix1 <- cbind(module8_hub_matrix,sum)
module8_hub_matrix1 <- as.data.frame(module8_hub_matrix1)
module8_hub_matrix1 <- module8_hub_matrix1[order(module8_hub_matrix1[,"sum"],decreasing = T),]
module8_hub_matrix1

module8_hub_matrix1 <- subset(module8_hub_matrix1, select = -sum)
module8_hub_matrix1 <- t(module8_hub_matrix1)
module8_hub_matrix1
#如果需要增减物种的数目，请调整此处的数值，由于核心微生物的数量为4，故无需进行top10
module8_hub_matrix2 <- module8_hub_matrix1[,1:10]
module8_hub_matrix2


##module8核心微生物出图
pdf(paste0(output,"module8_top10_hubtaxa.pdf"), width=7, height=7/3)
par(mfrow=c(1,3), mar=c(5,3,2,0.5))
mycol<- c("#D55E86","#63D07F")  

p1 <- barplot(module8_hub_matrix1,
              las=2, border=NA, axes=T, cex.names=.5,ylim=c(0,130),main="Module8",
              col=mycol[1:nrow(module8_hub_matrix1)],ylab="Degree ")

## legend 
plot.new()
legend("top", bty="n", cex=0.5, #x.intersp=0.1, #y.intersp=0.75,
       legend=rev(rownames(module8_hub_matrix1)), 
       fill=rev(mycol), 
       border=rev(mycol) )
dev.off()







####module8正负边/度####

##方法一：根据从总网络的提取结果
#78
unique_labels <- unique(module8_gephi2$Label)
unique_labels
# Initialize empty matrix
module8_all_matrix <- matrix(0, nrow=2, ncol=length(unique_labels), dimnames=list(c("positive","negative"), unique_labels))

# Loop through each label
for (i in unique_labels) {
  
  # Subset edges with the current label
  module8_all_degree <- all_network_edge[all_network_edge$Source == i | all_network_edge$Target == i,] 
  
  # Subset edges with positive correlation
  module8_all_positive_degree <- module8_all_degree[module8_all_degree$e_cor >= 0,]
  
  
  # Subset edges with negative correlation
  module8_all_negative_degree <- module8_all_degree[module8_all_degree$e_cor <= 0,]
  
  # Create matrix with positive and negative degrees
  module8_all_i <- matrix(c(dim(module8_all_positive_degree)[1], dim(module8_all_negative_degree)[1]), nrow=2, ncol=1, byrow=TRUE, dimnames=list(c("positive","negative"), i))
  module8_all_i
  # Fill in matrix with values from current label
  module8_all_matrix[, i] <- module8_all_i[, 1]
}

print(module8_all_matrix)
colnames(module8_all_matrix) <- rownames(module8_gephi2)
module8_all_matrix

dim(module8_all_matrix)
write.table(module8_all_matrix,paste0(output,"module8_all_matrix.txt"),sep="\t",quote=F)


##方法二：基于网络重新计算得到的数据

module8_network_edge <- read.table("module8_edge_list2.txt", row.names=1,sep="\t", header=T, blank.lines.skip=F, check.names=F)
#读取数据的时候，header=T表示保留原始表格的表头.row.names=1表示保留数据表格的第一列，=2表示去除表格的第一列
module8_network_edge


# Initialize empty matrix
module8_all_matrix2 <- matrix(0, nrow=2, ncol=length(unique_labels), dimnames=list(c("positive","negative"), unique_labels))

# Loop through each label
for (i in unique_labels) {
  
  # Subset edges with the current label
  module8_all_degree <- module8_network_edge[module8_network_edge$Source == i | module8_network_edge$Target == i,] 
  
  # Subset edges with positive correlation
  module8_all_positive_degree <- module8_all_degree[module8_all_degree$e_cor >= 0,]
  
  
  # Subset edges with negative correlation
  module8_all_negative_degree <- module8_all_degree[module8_all_degree$e_cor <= 0,]
  
  # Create matrix with positive and negative degrees
  module8_all_i <- matrix(c(dim(module8_all_positive_degree)[1], dim(module8_all_negative_degree)[1]), nrow=2, ncol=1, byrow=TRUE, dimnames=list(c("positive","negative"), i))
  module8_all_i
  # Fill in matrix with values from current label
  module8_all_matrix2[, i] <- module8_all_i[, 1]
}

print(module8_all_matrix2)
colnames(module8_all_matrix2) <- rownames(module8_gephi2)
module8_all_matrix2

dim(module8_all_matrix2)

write.table(module8_all_matrix2,paste0(output,"module8_all_matrix2.txt"),sep="\t",quote=F)





####module9####
#提取数据
module9_hub <- module9_gephi2[which(module9_gephi2$degree >= 60 & module9_gephi2$closnesscentrality >= 0.3),] 
module9_hub
dim(module9_hub) #module9中7个核心微生物
write.table(module9_hub,paste0(output,"module9_hub.txt"),sep="\t",quote=F)

#从大网络边数据列表中提取出核心微生物所对应的边数据
all_network_edge <- read.table("all_network_edge_list-4-2.txt", row.names=1,sep="\t", header=T, blank.lines.skip=F, check.names=F)
#读取数据的时候，header=T表示保留原始表格的表头.row.names=1表示保留数据表格的第一列，=2表示去除表格的第一列
all_network_edge

####循环的编写####
# Get unique labels
unique_labels <- unique(module9_hub$Label)
unique_labels
# Initialize empty matrix
module9_hub_matrix <- matrix(0, nrow=2, ncol=length(unique_labels), dimnames=list(c("positive","negative"), unique_labels))

# Loop through each label
for (i in unique_labels) {
  
  # Subset edges with the current label
  module9_hub_degree <- all_network_edge[all_network_edge$Source == i | all_network_edge$Target == i,] 
  
  # Subset edges with positive correlation
  module9_hub_positive_degree <- module9_hub_degree[module9_hub_degree$e_cor >= 0,]
  
  
  # Subset edges with negative correlation
  module9_hub_negative_degree <- module9_hub_degree[module9_hub_degree$e_cor <= 0,]
  
  # Create matrix with positive and negative degrees
  module9_hub_i <- matrix(c(dim(module9_hub_positive_degree)[1], dim(module9_hub_negative_degree)[1]), nrow=2, ncol=1, byrow=TRUE, dimnames=list(c("positive","negative"), i))
  module9_hub_i
  # Fill in matrix with values from current label
  module9_hub_matrix[, i] <- module9_hub_i[, 1]
}

print(module9_hub_matrix)
colnames(module9_hub_matrix) <- rownames(module9_hub)
module9_hub_matrix
write.table(module9_hub_matrix,paste0(output,"module9_hub_matrix.txt"),sep="\t",quote=F)


#数据处理一下，以适合用barplot绘制柱状图
module9_hub_matrix <- as.data.frame(t(module9_hub_matrix))
module9_hub_matrix
sum <- apply(module9_hub_matrix,1,sum) 
module9_hub_matrix1 <- cbind(module9_hub_matrix,sum)
module9_hub_matrix1 <- as.data.frame(module9_hub_matrix1)
module9_hub_matrix1 <- module9_hub_matrix1[order(module9_hub_matrix1[,"sum"],decreasing = T),]
module9_hub_matrix1

module9_hub_matrix1 <- subset(module9_hub_matrix1, select = -sum)
module9_hub_matrix1 <- t(module9_hub_matrix1)
module9_hub_matrix1
#如果需要增减物种的数目，请调整此处的数值，由于核心微生物的数量为4，故无需进行top10
module9_hub_matrix2 <- module9_hub_matrix1[,1:10]
module9_hub_matrix2


##特别需要注意的是，由于该模块核心微生物为9个，少于10个，所以人工补全了表格进行出图，完了在AI中删除即可
module9_hub_matrix3 <- read.table("module9_hub_matrix2.txt", row.names=1,sep="\t", header=T, blank.lines.skip=F, check.names=F)
module9_hub_matrix3

#数据处理一下，以适合用barplot绘制柱状图
module9_hub_matrix3 <-module9_hub_matrix3[,1:ncol(module9_hub_matrix3)]
module9_hub_matrix3<-as.matrix(module9_hub_matrix3)
module9_hub_matrix3
colnames(module9_hub_matrix3)


##module9核心微生物出图
pdf(paste0(output,"module9_top10_hubtaxa.pdf"), width=7, height=7/3)
par(mfrow=c(1,3), mar=c(7,3,2,0.5))
mycol<- c("#D55E86","#63D07F") 

p1 <- barplot(module9_hub_matrix3,
              las=2, border=NA, axes=T, cex.names=.5,ylim=c(0,130),main="Module9",
              col=mycol[1:nrow(module9_hub_matrix3)],ylab="Degree ")

## legend 
plot.new()
legend("top", bty="n", cex=0.5, #x.intersp=0.1, #y.intersp=0.75,
       legend=rev(rownames(module9_hub_matrix3)), 
       fill=rev(mycol), 
       border=rev(mycol) )
dev.off()







####module9正负边/度####

##方法一：根据从总网络的提取结果
#266
unique_labels <- unique(module9_gephi2$Label)
unique_labels
# Initialize empty matrix
module9_all_matrix <- matrix(0, nrow=2, ncol=length(unique_labels), dimnames=list(c("positive","negative"), unique_labels))

# Loop through each label
for (i in unique_labels) {
  
  # Subset edges with the current label
  module9_all_degree <- all_network_edge[all_network_edge$Source == i | all_network_edge$Target == i,] 
  
  # Subset edges with positive correlation
  module9_all_positive_degree <- module9_all_degree[module9_all_degree$e_cor >= 0,]
  
  
  # Subset edges with negative correlation
  module9_all_negative_degree <- module9_all_degree[module9_all_degree$e_cor <= 0,]
  
  # Create matrix with positive and negative degrees
  module9_all_i <- matrix(c(dim(module9_all_positive_degree)[1], dim(module9_all_negative_degree)[1]), nrow=2, ncol=1, byrow=TRUE, dimnames=list(c("positive","negative"), i))
  module9_all_i
  # Fill in matrix with values from current label
  module9_all_matrix[, i] <- module9_all_i[, 1]
}

print(module9_all_matrix)
colnames(module9_all_matrix) <- rownames(module9_gephi2)
module9_all_matrix

dim(module9_all_matrix)
write.table(module9_all_matrix,paste0(output,"module9_all_matrix.txt"),sep="\t",quote=F)


##方法二：基于网络重新计算得到的数据

module9_network_edge <- read.table("module9_edge_list2.txt", row.names=1,sep="\t", header=T, blank.lines.skip=F, check.names=F)
#读取数据的时候，header=T表示保留原始表格的表头.row.names=1表示保留数据表格的第一列，=2表示去除表格的第一列
module9_network_edge


# Initialize empty matrix
module9_all_matrix2 <- matrix(0, nrow=2, ncol=length(unique_labels), dimnames=list(c("positive","negative"), unique_labels))

# Loop through each label
for (i in unique_labels) {
  
  # Subset edges with the current label
  module9_all_degree <- module9_network_edge[module9_network_edge$Source == i | module9_network_edge$Target == i,] 
  
  # Subset edges with positive correlation
  module9_all_positive_degree <- module9_all_degree[module9_all_degree$e_cor >= 0,]
  
  
  # Subset edges with negative correlation
  module9_all_negative_degree <- module9_all_degree[module9_all_degree$e_cor <= 0,]
  
  # Create matrix with positive and negative degrees
  module9_all_i <- matrix(c(dim(module9_all_positive_degree)[1], dim(module9_all_negative_degree)[1]), nrow=2, ncol=1, byrow=TRUE, dimnames=list(c("positive","negative"), i))
  module9_all_i
  # Fill in matrix with values from current label
  module9_all_matrix2[, i] <- module9_all_i[, 1]
}

print(module9_all_matrix2)
colnames(module9_all_matrix2) <- rownames(module9_gephi2)
module9_all_matrix2

dim(module9_all_matrix2)

write.table(module9_all_matrix2,paste0(output,"module9_all_matrix2.txt"),sep="\t",quote=F)





####module10####
#提取数据
module10_hub <- module10_gephi2[which(module10_gephi2$degree >= 60 & module10_gephi2$closnesscentrality >= 0.3),] 
module10_hub
dim(module10_hub) #module10中没有核心微生物
write.table(module10_hub,paste0(output,"module10_hub.txt"),sep="\t",quote=F)

#从大网络边数据列表中提取出核心微生物所对应的边数据
all_network_edge <- read.table("all_network_edge_list-4-2.txt", row.names=1,sep="\t", header=T, blank.lines.skip=F, check.names=F)
#读取数据的时候，header=T表示保留原始表格的表头.row.names=1表示保留数据表格的第一列，=2表示去除表格的第一列
all_network_edge

####循环的编写####
# Get unique labels
unique_labels <- unique(module10_hub$Label)
unique_labels
# Initialize empty matrix
module10_hub_matrix <- matrix(0, nrow=2, ncol=length(unique_labels), dimnames=list(c("positive","negative"), unique_labels))

# Loop through each label
for (i in unique_labels) {
  
  # Subset edges with the current label
  module10_hub_degree <- all_network_edge[all_network_edge$Source == i | all_network_edge$Target == i,] 
  
  # Subset edges with positive correlation
  module10_hub_positive_degree <- module10_hub_degree[module10_hub_degree$e_cor >= 0,]
  
  
  # Subset edges with negative correlation
  module10_hub_negative_degree <- module10_hub_degree[module10_hub_degree$e_cor <= 0,]
  
  # Create matrix with positive and negative degrees
  module10_hub_i <- matrix(c(dim(module10_hub_positive_degree)[1], dim(module10_hub_negative_degree)[1]), nrow=2, ncol=1, byrow=TRUE, dimnames=list(c("positive","negative"), i))
  module10_hub_i
  # Fill in matrix with values from current label
  module10_hub_matrix[, i] <- module10_hub_i[, 1]
}

print(module10_hub_matrix)
colnames(module10_hub_matrix) <- rownames(module10_hub)
module10_hub_matrix
write.table(module10_hub_matrix,paste0(output,"module10_hub_matrix.txt"),sep="\t",quote=F)


#数据处理一下，以适合用barplot绘制柱状图
module10_hub_matrix <- as.data.frame(t(module10_hub_matrix))
module10_hub_matrix
sum <- apply(module10_hub_matrix,1,sum) 
module10_hub_matrix1 <- cbind(module10_hub_matrix,sum)
module10_hub_matrix1 <- as.data.frame(module10_hub_matrix1)
module10_hub_matrix1 <- module10_hub_matrix1[order(module10_hub_matrix1[,"sum"],decreasing = T),]
module10_hub_matrix1

module10_hub_matrix1 <- subset(module10_hub_matrix1, select = -sum)
module10_hub_matrix1 <- t(module10_hub_matrix1)
module10_hub_matrix1
#如果需要增减物种的数目，请调整此处的数值，由于核心微生物的数量为4，故无需进行top10
module10_hub_matrix2 <- module10_hub_matrix1[,1:10]
module10_hub_matrix2


##module10核心微生物出图
pdf(paste0(output,"module10_top10_hubtaxa.pdf"), width=7, height=7/3)
par(mfrow=c(1,3), mar=c(5,3,2,0.5))
mycol<- c("#D55E86","#63D07F") 

p1 <- barplot(module10_hub_matrix1,
              las=2, border=NA, axes=T, cex.names=.5,ylim=c(0,130),main="Module10",
              col=mycol[1:nrow(module10_hub_matrix1)],ylab="Degree ")

## legend 
plot.new()
legend("top", bty="n", cex=0.5, #x.intersp=0.1, #y.intersp=0.75,
       legend=rev(rownames(module10_hub_matrix1)), 
       fill=rev(mycol), 
       border=rev(mycol) )
dev.off()


####module10正负边/度####

##方法一：根据从总网络的提取结果
#106
unique_labels <- unique(module10_gephi2$Label)
unique_labels
# Initialize empty matrix
module10_all_matrix <- matrix(0, nrow=2, ncol=length(unique_labels), dimnames=list(c("positive","negative"), unique_labels))

# Loop through each label
for (i in unique_labels) {
  
  # Subset edges with the current label
  module10_all_degree <- all_network_edge[all_network_edge$Source == i | all_network_edge$Target == i,] 
  
  # Subset edges with positive correlation
  module10_all_positive_degree <- module10_all_degree[module10_all_degree$e_cor >= 0,]
  
  
  # Subset edges with negative correlation
  module10_all_negative_degree <- module10_all_degree[module10_all_degree$e_cor <= 0,]
  
  # Create matrix with positive and negative degrees
  module10_all_i <- matrix(c(dim(module10_all_positive_degree)[1], dim(module10_all_negative_degree)[1]), nrow=2, ncol=1, byrow=TRUE, dimnames=list(c("positive","negative"), i))
  module10_all_i
  # Fill in matrix with values from current label
  module10_all_matrix[, i] <- module10_all_i[, 1]
}

print(module10_all_matrix)
colnames(module10_all_matrix) <- rownames(module10_gephi2)
module10_all_matrix

dim(module10_all_matrix)
write.table(module10_all_matrix,paste0(output,"module10_all_matrix.txt"),sep="\t",quote=F)


##方法二：基于网络重新计算得到的数据

module10_network_edge <- read.table("module10_edge_list2.txt", row.names=1,sep="\t", header=T, blank.lines.skip=F, check.names=F)
#读取数据的时候，header=T表示保留原始表格的表头.row.names=1表示保留数据表格的第一列，=2表示去除表格的第一列
module10_network_edge


# Initialize empty matrix
module10_all_matrix2 <- matrix(0, nrow=2, ncol=length(unique_labels), dimnames=list(c("positive","negative"), unique_labels))

# Loop through each label
for (i in unique_labels) {
  
  # Subset edges with the current label
  module10_all_degree <- module10_network_edge[module10_network_edge$Source == i | module10_network_edge$Target == i,] 
  
  # Subset edges with positive correlation
  module10_all_positive_degree <- module10_all_degree[module10_all_degree$e_cor >= 0,]
  
  
  # Subset edges with negative correlation
  module10_all_negative_degree <- module10_all_degree[module10_all_degree$e_cor <= 0,]
  
  # Create matrix with positive and negative degrees
  module10_all_i <- matrix(c(dim(module10_all_positive_degree)[1], dim(module10_all_negative_degree)[1]), nrow=2, ncol=1, byrow=TRUE, dimnames=list(c("positive","negative"), i))
  module10_all_i
  # Fill in matrix with values from current label
  module10_all_matrix2[, i] <- module10_all_i[, 1]
}

print(module10_all_matrix2)
colnames(module10_all_matrix2) <- rownames(module10_gephi2)
module10_all_matrix2

dim(module10_all_matrix2)

write.table(module10_all_matrix2,paste0(output,"module10_all_matrix2.txt"),sep="\t",quote=F)





data<-read.table("degree_boxplot2.txt",sep="\t",header=TRUE)
data
#install.packages("ggplot2")
library("ggplot2")
p<-ggplot(data)+
  geom_boxplot(mapping =aes(x=group1, y=degree, color=group1),alpha=0.5,size=1.0,width=0.6)+
  geom_jitter(mapping =aes(x=group1, y=degree, color=group1), alpha=0.3,size=1.0)+
  scale_color_manual(limits=c("module3", "module5", "module6", "module9"),
                     values = c("#394A92","#497EB2","#D2352C","#D55E86"))
p
#分面横向排成一排
#p2<-p+ facet_grid(.~group2)
#p2
#分面排成2行2列
#p3<-p+ facet_wrap(.~group2)
#p3
#图形整体调整
mytheme<-theme_bw()+theme(axis.title=element_text(size=14),axis.text = element_text(size = 14),
                          panel.grid.major = element_line(color = "white"),
                          panel.grid.minor = element_line(color = "white"),
                          axis.text.x = element_text(size = 14, angle=30,vjust = 0.5,hjust = 0.5,color = "black"),
                          axis.text.y = element_text(size=14,color = "black"),
                          legend.position = "none",
                          plot.margin = unit(c(0.5, 0.5, 0.5, 0.5),"cm"))
#有关图例的调整代码
#图例字体大小legend.text = element_text(size=14),
#图例名称legend.title = element_blank(),
#图例有无legend.position=""
#p2+mytheme

p4<-p+mytheme
p4
#p5<-p+coord_cartesian(ylim = c(0,180),expand = TRUE)
#p5


#特别需要注意的是，R语言中进行统计分析，一定要保证不同分组变量的长度（个数）一致，否则会报错长度不一致。只能分析好差异显著性后进行标记

#install.packages("ggpubr")
#library("ggpubr")

#默认为"wilcox.test"这是两组的非参数检验，微生物生态方面我们常用非参数检验，但是这里我们是多个组的，因此，更改为Kruskal-Wallis检验
#compare_means(rootlength~group2,data=data, paired = TRUE,method = "kruskal.test")

#在这里我们就两两检验，看看那些组显著：
#df_summarise<-compare_means(rootlength~group1,data=data, paired = TRUE)
# compare_means的计算是以tbl格式返回的，我们可以转成dataframe
#df_summarise<-as.data.frame(df_summarise)

#write.table(df_summarise,"rootlength—wilcox.test.txt",quote= FALSE,row.names = T,
#            col.names = T,sep = "\t")

#手动添加显著性标记
#install.packages("doBy")
library(doBy)
rootlength_stat <- summaryBy(degree~group1+group2+variable, data, FUN = max)
names(rootlength_stat) <- c('group1', 'group2', 'degree')
rootlength_stat$sig <- c('a', 'b', 'a', 'c')

p6<-p4+geom_text(data=rootlength_stat, aes(x=group1,y=degree+15,label=sig),size=6,position=position_dodge(0.6))
p6

p7<-p6+xlab("")+ylab("Degree")
p7
ggsave("degree.pdf", width = 6, height = 8, units="cm")
ggsave("degree.tiff", width = 6, height = 8, units = "cm")
ggsave("degree.png", width = 6, height = 8, units = "cm")





module_edge <- read.table("module_edge.txt", header=T,row.names = 1,sep = "\t")
module_edge <- as.data.frame(module_edge)
module_edge 
module_edge3 <-module_edge[,1:ncol(module_edge)]
module_edge3<-as.matrix(module_edge3)
module_edge3
colnames(module_edge3)
##module3核心微生物出图
pdf(paste0(output,"module_edge.pdf"), width=7, height=7/6)
par(mfrow=c(1,5), mar=c(2,3,2,0.5))
mycol<- c("#D55E86","#63D07F") 


p1 <- barplot(module_edge3,
              las=2, border=NA, axes=T, cex.names=.5,
              col=mycol[1:nrow(module_edge)],ylab="Degree ")

## legend ,ylim=c(0,130),,main="Module3"
plot.new()
legend("top", bty="n", cex=0.5, #x.intersp=0.1, #y.intersp=0.75,
       legend=rev(rownames(module_edge)), 
       fill=rev(mycol), 
       border=rev(mycol) )
dev.off()






module_edge2 <- read.table("module_edge2.txt", header=T,row.names = 1,sep = "\t")
module_edge2 <- as.data.frame(module_edge2)
module_edge2
module_edge23 <-module_edge2[,1:ncol(module_edge2)]
module_edge23<-as.matrix(module_edge23)
module_edge23
colnames(module_edge23)
##module3核心微生物出图
pdf(paste0(output,"module_edge2.pdf"), width=7, height=8/4)
par(mfrow=c(1,4), mar=c(5,6,2,0.5))
mycol<- c("#D55E86","#63D07F") 


p1 <- barplot(module_edge23,
              las=2, border=NA, axes=T, cex.names=.5,ylim=c(0,5000),
              col=mycol[1:nrow(module_edge)],ylab="Degree ")

## legend ,ylim=c(0,130),,main="Module3"
plot.new()
legend("top", bty="n", cex=0.5, #x.intersp=0.1, #y.intersp=0.75,
       legend=rev(rownames(module_edge)), 
       fill=rev(mycol), 
       border=rev(mycol) )
dev.off()



####不同处理在同一模块中的贡献度####
###每个模块中各个处理间OTU的丰度差异，首先得到每个模块中OTU所对应的物种丰度，再做堆叠柱状图
### plotting average module response

## ROOT module 0
module0_rownames <- rownames(module0_gephi2)
module0_rownames #136

module0_rownames <-as.data.frame(module0_rownames)

#注意，根据两个矩阵各自指定的数据进行取交集
module0_abudance <- otu_norm_16s[rownames(otu_norm_16s) %in% module0_rownames$module0_rownames,]
module0_abudance

#module0_abudance <- sensitiv_otu[module0_rownames,]
#module0_abudance

#module0_abudance<-na.omit(module0_abudance)
#module0_abudance
#dim(module0_abudance)#29,4
#module0_abudance <-as.data.frame(module0_abudance)


#module0_abudance <- otu_norm_16s[rownames(module0_abudance),]
#module0_abudance

#dim(module0_abudance) #29,12

write.table(module0_abudance,paste0(output,"all_module0_abudance2.txt"),sep="\t",quote=F)

group <- read.table("mapping.txt",header=T,row.names = 1,sep = "\t")
group


pdf(paste0(output,"all_module0_abudance2.pdf"),width=7,height=7/3)
par(mfrow=c(1,3), mar=c(0.5,5,2,0))

CS_cols <- c("#32A852","#63D07F","#C9244C","#D55E86")
names(CS_cols) <- c("Ggt","S_135", "SA","H")

bargraph.CI(group$Group, colSums(module0_abudance)/1000, 
            las=2, ylab="cumulative relative abundance", cex.lab=.5, cex.axis=.7, cex.names=.7,
            err.width=.025, main="Module 0", col=CS_cols, border=F)
dev.off()



## ROOT module 1
module1_rownames <- rownames(module1_gephi2)
module1_rownames #162
module1_rownames <- as.data.frame(module1_rownames)
module1_rownames
#求otu_norm_16s和module1_rownames两个的交集
module1_abudance <- otu_norm_16s[rownames(otu_norm_16s) %in% module1_rownames$module1_rownames,]
module1_abudance

write.table(module1_abudance,paste0(output,"all_module1_abudance2.txt"),sep="\t",quote=F)

pdf(paste0(output,"all_module1_abudance2.pdf"),width=7,height=7/3)
par(mfrow=c(1,3), mar=c(0.5,5,2,0))

CS_cols <- c("#32A852","#63D07F","#C9244C","#D55E86")
names(CS_cols) <- c("Ggt","S_135", "SA","H")

bargraph.CI(group$Group, colSums(module1_abudance)/1000, 
            las=2, ylab="cumulative relative abundance", cex.lab=.5, cex.axis=.7, cex.names=.7,
            err.width=.025, main="Module 1", col=CS_cols, border=F)
dev.off()




## ROOT module 2
module2_rownames <- rownames(module2_gephi2)
module2_rownames

module2_rownames <- as.data.frame(module2_rownames)
module2_rownames

module2_abudance <- otu_norm_16s[rownames(otu_norm_16s) %in% module2_rownames$module2_rownames,]
module2_abudance

write.table(module2_abudance,paste0(output,"all_module2_abudance2.txt"),sep="\t",quote=F)

pdf(paste0(output,"all_module2_abudance2.pdf"),width=7,height=7/3)
par(mfrow=c(1,3), mar=c(0.5,5,2,0))

CS_cols <- c("#32A852","#63D07F","#C9244C","#D55E86")
names(CS_cols) <- c("Ggt","S_135", "SA","H")


bargraph.CI(group$Group, colSums(module2_abudance)/1000, 
            las=2, ylab="", cex.lab=.5, cex.axis=.7, cex.names=.7,
            err.width=.025, main="Module 2", col=CS_cols, border=F)

dev.off()


## ROOT module 3
module3_rownames <- rownames(module3_gephi2)
module3_rownames
module3_rownames <- as.data.frame(module3_rownames)
module3_rownames

module3_abudance <- otu_norm_16s[rownames(otu_norm_16s) %in% module3_rownames$module3_rownames,]
module3_abudance

write.table(module3_abudance,paste0(output,"all_module3_abudance2.txt"),sep="\t",quote=F)


pdf(paste0(output,"all_module3_abudance2.pdf"),width=7,height=7/3)
par(mfrow=c(1,3), mar=c(0.5,5,2,0))

CS_cols <- c("#32A852","#63D07F","#C9244C","#D55E86")
names(CS_cols) <- c("Ggt","S_135", "SA","H")


bargraph.CI(group$Group, colSums(module3_abudance)/1000, 
            las=2, ylab="", cex.lab=.5, cex.axis=.7, cex.names=.7,
            err.width=.025, main="Module 3", col=CS_cols, border=F)

dev.off()


## ROOT module 4
module4_rownames <- rownames(module4_gephi2)
module4_rownames
module4_rownames <- as.data.frame(module4_rownames)
module4_rownames

module4_abudance <- otu_norm_16s[rownames(module4_rownames) %in% module4_rownames$module4_rownames,]
module4_abudance
write.table(module4_abudance,paste0(output,"all_module4_abudance2.txt"),sep="\t",quote=F)


pdf(paste0(output,"all_module4_abudance2.pdf"),width=7,height=7/3)
par(mfrow=c(1,3), mar=c(0.5,5,2,0))

CS_cols <- c("#32A852","#63D07F","#C9244C","#D55E86")
names(CS_cols) <- c("Ggt","S_135", "SA","H")


bargraph.CI(group$Group, colSums(module4_abudance)/1000, 
            las=2, ylab="", cex.lab=.5, cex.axis=.7, cex.names=.7,
            err.width=.025, main="Module 4", col=CS_cols, border=F)

dev.off()




## ROOT module 5
module5_rownames <- rownames(module5_gephi2)
module5_rownames
module5_rownames <- as.data.frame(module5_rownames)
module5_rownames

module5_abudance <- otu_norm_16s[rownames(otu_norm_16s) %in% module5_rownames$module5_rownames,]
module5_abudance
write.table(module5_abudance,paste0(output,"all_module5_abudance2.txt"),sep="\t",quote=F)


pdf(paste0(output,"all_module5_abudance2.pdf"),width=7,height=7/3)
par(mfrow=c(1,3), mar=c(0.5,5,2,0))

CS_cols <- c("#32A852","#63D07F","#C9244C","#D55E86")
names(CS_cols) <- c("Ggt","S_135", "SA","H")


bargraph.CI(group$Group, colSums(module5_abudance)/1000, 
            las=2, ylab="", cex.lab=.5, cex.axis=.7, cex.names=.7,
            err.width=.025, main="Module 5", col=CS_cols, border=F)

dev.off()




## ROOT module 6
module6_rownames <- rownames(module6_gephi2)
module6_rownames
module6_rownames <- as.data.frame(module6_rownames)
module6_rownames

module6_abudance <- otu_norm_16s[rownames(otu_norm_16s) %in% module6_rownames$module6_rownames,]
module6_abudance
write.table(module6_abudance,paste0(output,"all_module6_abudance2.txt"),sep="\t",quote=F)


pdf(paste0(output,"all_module6_abudance2.pdf"),width=7,height=7/3)
par(mfrow=c(1,3), mar=c(0.5,5,2,0))

CS_cols <- c("#32A852","#63D07F","#C9244C","#D55E86")
names(CS_cols) <- c("Ggt","S_135", "SA","H")


bargraph.CI(group$Group, colSums(module6_abudance)/1000, 
            las=2, ylab="", cex.lab=.5, cex.axis=.7, cex.names=.7,
            err.width=.025, main="Module 6", col=CS_cols, border=F)

dev.off()




## ROOT module 7
module7_rownames <- rownames(module7_gephi2)
module7_rownames
module7_rownames <- as.data.frame(module7_rownames)
module7_rownames

module7_abudance <- otu_norm_16s[rownames(otu_norm_16s) %in% module7_rownames$module7_rownames,]
module7_abudance
write.table(module7_abudance,paste0(output,"all_module7_abudance2.txt"),sep="\t",quote=F)


pdf(paste0(output,"all_module7_abudance2.pdf"),width=7,height=7/3)
par(mfrow=c(1,3), mar=c(0.5,5,2,0))

CS_cols <- c("#32A852","#63D07F","#C9244C","#D55E86")
names(CS_cols) <- c("Ggt","S_135", "SA","H")


bargraph.CI(group$Group, colSums(module7_abudance)/1000, 
            las=2, ylab="", cex.lab=.5, cex.axis=.7, cex.names=.7,
            err.width=.025, main="Module 7", col=CS_cols, border=F)

dev.off()




## ROOT module 8
module8_rownames <- rownames(module8_gephi2)
module8_rownames
module8_rownames <- as.data.frame(module8_rownames)
module8_rownames

module8_abudance <- otu_norm_16s[rownames(otu_norm_16s) %in% module8_rownames$module8_rownames,]
module8_abudance
write.table(module8_abudance,paste0(output,"all_module8_abudance2.txt"),sep="\t",quote=F)


pdf(paste0(output,"all_module8_abudance2.pdf"),width=7,height=7/3)
par(mfrow=c(1,3), mar=c(0.5,5,2,0))

CS_cols <- c("#32A852","#63D07F","#C9244C","#D55E86")
names(CS_cols) <- c("Ggt","S_135", "SA","H")


bargraph.CI(group$Group, colSums(module8_abudance)/1000, 
            las=2, ylab="", cex.lab=.5, cex.axis=.7, cex.names=.7,
            err.width=.025, main="Module 8", col=CS_cols, border=F)

dev.off()




## ROOT module 9
module9_rownames <- rownames(module9_gephi2)
module9_rownames
module9_rownames <- as.data.frame(module9_rownames)
module9_rownames

module9_abudance <- otu_norm_16s[rownames(otu_norm_16s) %in% module9_rownames$module9_rownames,]
module9_abudance
write.table(module9_abudance,paste0(output,"all_module9_abudance2.txt"),sep="\t",quote=F)


pdf(paste0(output,"all_module9_abudance2.pdf"),width=7,height=7/3)
par(mfrow=c(1,3), mar=c(0.5,5,2,0))

CS_cols <- c("#32A852","#63D07F","#C9244C","#D55E86")
names(CS_cols) <- c("Ggt","S_135", "SA","H")


bargraph.CI(group$Group, colSums(module9_abudance)/1000, 
            las=2, ylab="", cex.lab=.5, cex.axis=.7, cex.names=.7,
            err.width=.025, main="Module 9", col=CS_cols, border=F)

dev.off()




## ROOT module 10
module10_rownames <- rownames(module10_gephi2)
module10_rownames
module10_rownames <- as.data.frame(module10_rownames)

module10_abudance <- otu_norm_16s[rownames(otu_norm_16s) %in% module10_rownames$module10_rownames,]
module10_abudance
write.table(module10_abudance,paste0(output,"all_module10_abudance2.txt"),sep="\t",quote=F)


pdf(paste0(output,"all_module10_abudance5.pdf"),width=7,height=7/3)
par(mfrow=c(1,3), mar=c(5,5,2,0))

CS_cols <- c("#32A852","#63D07F","#C9244C","#D55E86")
names(CS_cols) <- c("Ggt","S_135", "SA","H")


bargraph.CI(group$Group, colSums(module10_abudance)/1000, 
            las=2, ylab="", cex.lab=.5, cex.axis=.7, cex.names=.7,
            err.width=.025, main="Module 10", col=CS_cols, border=F)

dev.off()






## All

all_rownames <- rownames(network_feature)
all_rownames
all_rownames <- as.data.frame(all_rownames)
all_rownames
head(all_rownames)
all_abudance <- otu_norm_16s[rownames(otu_norm_16s) %in% all_rownames$all_rownames,]
all_abudance
write.table(all_abudance,paste0(output,"all_abudance2.txt"),sep="\t",quote=F)


pdf(paste0(output,"all_abudance2.pdf"),width=7,height=7/6)
par(mfrow=c(1,7), mar=c(0.5,3.5,2,0))

CS_cols <- c("dodgerblue4","orange","green","violetred4")
names(CS_cols) <- c("Ggt","S_135", "SA","H")


bargraph.CI(group$Group, colSums(all_abudance)/1000, 
            las=2, ylab="", cex.lab=.5, cex.axis=.7, cex.names=.7,
            err.width=.025, main="All", col=CS_cols, border=F)
dev.off()



## Legend
plot.new()
par(mar=c(0.5,0,2,0))
legend("left",  bty="n", cex=1, #x.intersp=0.1, y.intersp=1,
       legend=names(CS_cols), 
       fill=CS_cols, 
       border=CS_cols , xpd=T)

dev.off()




###不同模块OTU的物种分类单位
#### taxonomies of module OTUs 

## ROOT module 0
## defining bacteria of root modules

module0_rownames <- rownames(module0_gephi2)
module0_rownames

### counts of bacteria OTUs
module0_taxonomy <- as.data.frame(table(tax_filter_16s[module0_rownames, "labels"] ) )
module0_taxonomy
colnames(module0_taxonomy) <- c("phylum", "Module 0")
module0_taxonomy



## ROOT module 1
## defining bacteria of root modules

module1_rownames <- rownames(module1_gephi2)
module1_rownames

### counts of bacteria OTUs
module1_taxonomy <- as.data.frame(table(tax_filter_16s[module1_rownames, "labels"] ) )
module1_taxonomy
colnames(module1_taxonomy) <- c("phylum", "Module 1")
module1_taxonomy



## ROOT module 2
## defining bacteria of root modules

module2_rownames <- rownames(module2_gephi2)
module2_rownames

### counts of bacteria OTUs
module2_taxonomy <- as.data.frame(table(tax_filter_16s[module2_rownames, "labels"] ) )
module2_taxonomy
colnames(module2_taxonomy) <- c("phylum", "Module 2")
module2_taxonomy



## ROOT module 3
## defining bacteria of root modules

module3_rownames <- rownames(module3_gephi2)
module3_rownames

### counts of bacteria OTUs
module3_taxonomy <- as.data.frame(table(tax_filter_16s[module3_rownames, "labels"] ) )
module3_taxonomy
colnames(module3_taxonomy) <- c("phylum", "Module 3")
module3_taxonomy



## ROOT module 4
## defining bacteria of root modules

module4_rownames <- rownames(module4_gephi2)
module4_rownames

### counts of bacteria OTUs
module4_taxonomy <- as.data.frame(table(tax_filter_16s[module4_rownames, "labels"] ) )
module4_taxonomy
colnames(module4_taxonomy) <- c("phylum", "Module 4")
module4_taxonomy



## ROOT module 5
## defining bacteria of root modules

module5_rownames <- rownames(module5_gephi2)
module5_rownames

### counts of bacteria OTUs
module5_taxonomy <- as.data.frame(table(tax_filter_16s[module5_rownames, "labels"] ) )
module5_taxonomy
colnames(module5_taxonomy) <- c("phylum", "Module 5")
module5_taxonomy



## ROOT module 6
## defining bacteria of root modules

module6_rownames <- rownames(module6_gephi2)
module6_rownames

### counts of bacteria OTUs
module6_taxonomy <- as.data.frame(table(tax_filter_16s[module6_rownames, "labels"] ) )
module6_taxonomy
colnames(module6_taxonomy) <- c("phylum", "Module 6")
module6_taxonomy



## ROOT module 7
## defining bacteria of root modules

module7_rownames <- rownames(module7_gephi2)
module7_rownames

### counts of bacteria OTUs
module7_taxonomy <- as.data.frame(table(tax_filter_16s[module7_rownames, "labels"] ) )
module7_taxonomy
colnames(module7_taxonomy) <- c("phylum", "Module 7")
module7_taxonomy



## ROOT module 8
## defining bacteria of root modules

module8_rownames <- rownames(module8_gephi2)
module8_rownames

### counts of bacteria OTUs
module8_taxonomy <- as.data.frame(table(tax_filter_16s[module8_rownames, "labels"] ) )
module8_taxonomy
colnames(module8_taxonomy) <- c("phylum", "Module 8")
module8_taxonomy



## ROOT module 9
## defining bacteria of root modules

module9_rownames <- rownames(module9_gephi2)
module9_rownames

### counts of bacteria OTUs
module9_taxonomy <- as.data.frame(table(tax_filter_16s[module9_rownames, "labels"] ) )
module9_taxonomy
colnames(module9_taxonomy) <- c("phylum", "Module 9")
module9_taxonomy



## ROOT module 10
## defining bacteria of root modules

module10_rownames <- rownames(module10_gephi2)
module10_rownames

### counts of bacteria OTUs
module10_taxonomy <- as.data.frame(table(tax_filter_16s[module10_rownames, "labels"] ) )
module10_taxonomy
colnames(module10_taxonomy) <- c("phylum", "Module 10")
module10_taxonomy



## ROOT module 11
## defining bacteria of root modules

#module11_rownames <- rownames(module11_gephi2)
#module11_rownames

### counts of bacteria OTUs
#module11_taxonomy <- as.data.frame(table(tax_filter_16s[module11_rownames, "labels"] ) )
#module11_taxonomy
#colnames(module11_taxonomy) <- c("phylum", "Module 11")
#module11_taxonomy


module <- merge(module0_taxonomy, module1_taxonomy, all=T, by="phylum") 
module


module <- merge(module, module2_taxonomy, all=T, by="phylum") 
module

module <- merge(module, module3_taxonomy, all=T, by="phylum") 
module

module <- merge(module, module4_taxonomy, all=T, by="phylum") 
module

module <- merge(module, module5_taxonomy, all=T, by="phylum") 
module

module <- merge(module, module6_taxonomy, all=T, by="phylum") 
module

module <- merge(module, module7_taxonomy, all=T, by="phylum") 
module

module <- merge(module, module8_taxonomy, all=T, by="phylum") 
module

module <- merge(module, module9_taxonomy, all=T, by="phylum") 
module


module <- merge(module, module10_taxonomy, all=T, by="phylum") 
module

#module <- merge(module, module11_taxonomy, all=T, by="phylum") 
#module


all_taxonomy <- as.data.frame(table(tax_filter_16s[, "labels"] ) )
colnames(all_taxonomy) <- c("phylum", "all")

module <- merge(module, all_taxonomy, all=T, by="phylum") 
module

bacteria_modules <- module
bacteria_modules

bacteria_modules_mat <- bacteria_modules[2:13]
rownames(bacteria_modules_mat) <- bacteria_modules$phylum
bacteria_modules_mat[is.na(bacteria_modules_mat)] <- 0
colSums(bacteria_modules_mat)

bacteria_modules_prop <- t(t(bacteria_modules_mat)/colSums(bacteria_modules_mat) ) * 1
bacteria_modules_prop
colSums(bacteria_modules_prop)

write.table(bacteria_modules_prop,paste0(output,"bacteria_modules_prop.txt"),sep="\t",quote=F)


###手动添加颜色
phylum <- table(PHYLA_label_cols_16s_legend$labels)
phylum <- sort(phylum, decreasing = TRUE)
phylum
#这个示例中共计 12 类细菌门，因此手动分配 12 种颜色代指
#大家使用自己的数据时，视情况修改颜色数量和赋值
color <- c('#68AC57','#F4C28F','#394A92', '#497EB2', '#8E549E', '#394A92', '#E600E6', "lightgray",'#394A92', '#D2352C',"lightgray","lightgray","lightgray")
#, '#B3DE69', '#FDB462'
names(color) <- c(names(phylum[1]),names(phylum[2]),names(phylum[3]),names(phylum[4]),names(phylum[5]),names(phylum[6]),names(phylum[7]),names(phylum[8]),names(phylum[9]),names(phylum[10]),names(phylum[11]),names(phylum[12]),names(phylum[13]))

#names(color) <- names(c(phylum[1],phylum[2],phylum[3],phylum[4],phylum[5],phylum[6],phylum[7],phylum[8],phylum[9],phylum[10]))

color


##图4c
pdf(paste0(output,"all_module_taxonomy.pdf"), width=7, height=7/3)
par(mfrow=c(1,3), mar=c(2,4,1,0))

### bacteria
# PHYLA_label_cols_16s_legend
table(rownames(bacteria_modules_prop) %in% PHYLA_label_cols_16s_legend$labels) 
bp <- barplot(bacteria_modules_prop[,1:12],
              las=2, border=NA, axes=F, cex.names=.5,
              col = color)
              #col=PHYLA_label_cols_16s_legend[rownames(bacteria_modules_prop),]$cols )
text(bp, 1.1, labels=c(colSums(bacteria_modules_mat)[1:13], xpd=T, cex=.6, las=2))

## legend from Fig. S3
plot.new()
legend("left", bty="n", cex=0.5, #x.intersp=0.1, #y.intersp=0.75,
       legend=rev(PHYLA_label_cols_16s_legend$labels), 
       fill=rev(PHYLA_label_cols_16s_legend$cols), 
       border=rev(PHYLA_label_cols_16s_legend$cols) )
dev.off()

#cbind(bacteria_modules_prop[,1:16], bacteria_modules_prop[,4]),

####所有模块
modules_prop1 <- read.table("bacteria_modules_prop4.txt",header = T, row.names = 1,sep = "\t")
modules_prop1 

color <- c('#394A92','#68AC57','#8E549E','#F4C28F', '#497EB2', '#D2352C', '#E600E6', "lightgray")

color
# encoding="MacRoman",
#pdf(paste0(output,"all_module_taxonomy.pdf"), height=6, width=10, paper="a4r")
pdf(paste0(output,"all_module_tax.pdf"),width=7,height=7/3)
par(mfrow=c(1,3), mar=c(0.5,5,2,0))

#layout(matrix(c(1,2),1,2, byrow=F))

#par(oma=c(0,0,0,0), mar=c(6,6,1,1), xpd=NA)
phylum_bar_16s <- barplot(as.matrix(modules_prop1), col=color,
                          ylim=c(0,1), xaxt="n", border=NA, las=2)
axis(1, at=phylum_bar_16s, labels=colnames(modules_prop1), col.axis="black", las=2, cex.axis=0.6)
title(ylab="Relative abundance (%)")
#title(main="Bacteria Community")
#legend(43, 100, bty="n", cex=0.7, x.intersp=0.1, y.intersp=1,
#       legend=rev(PHYLA_label_cols_16s_legend$labels), 
#       fill=rev(PHYLA_label_cols_16s_legend$cols), 
#       border=rev(PHYLA_label_cols_16s_legend$cols) )
dev.off()



####高丰度4个模块
modules_prop2 <- read.table("bacteria_modules_prop5.txt",header = T, row.names = 1,sep = "\t")
modules_prop2 

color <- c('#394A92','#68AC57','#8E549E','#F4C28F', '#497EB2', '#D2352C', '#E600E6', "lightgray")

color

pdf(paste0(output,"Figure10.pdf"), encoding="MacRoman", height=3, width=5, paper="a4r")

layout(matrix(c(1,2),1,2, byrow=F))

par(oma=c(0,0,0,0), mar=c(6,6,1,1), xpd=NA)
phylum_bar_16s <- barplot(as.matrix(modules_prop2), col=color,
                          ylim=c(0,1), xaxt="n", border=NA, las=2)
axis(1, at=phylum_bar_16s, labels=colnames(modules_prop2), col.axis="black", las=2, cex.axis=0.6)
title(ylab="Relative abundance (%)")
#title(main="Bacteria Community")
#legend(43, 100, bty="n", cex=0.7, x.intersp=0.1, y.intersp=1,
#       legend=rev(PHYLA_label_cols_16s_legend$labels), 
#       fill=rev(PHYLA_label_cols_16s_legend$cols), 
#       border=rev(PHYLA_label_cols_16s_legend$cols) )
dev.off()







####R语言统计各网络节点的度频数分布####

####all_network####
all_network_node <- read.table("all_network_node_list-4.txt", row.names=1,sep="\t", header=T, blank.lines.skip=F, check.names=F)
#读取数据的时候，header=T表示保留原始表格的表头.row.names=1表示保留数据表格的第一列，=2表示去除表格的第一列
all_network_node


####循环的编写####
# Get unique labels
#unique_labels <- unique(all_network_node$degree)
#unique_labels

# Initialize empty matrix
#all_network_degree_matrix <- matrix(0, nrow=length(unique_labels), ncol=2, dimnames=list(unique_labels, c("degree","count")))

# Loop through each label
#for (i in unique_labels) {

# Subset edges with the current label
#  all_network_degree_count <- all_network_node[all_network_node$degree == i,] 
#  dim(all_network_degree_count)
# Create matrix with positive and negative degrees
#  all_network_degree_i <- matrix(c(i, dim(all_network_degree_count)[1]), nrow=1, ncol=2, byrow=TRUE, dimnames=list(i,c("degree","count")))
#  all_network_degree_i
# Fill in matrix with values from current label
#  all_network_degree_matrix[i, ] <- all_network_degree_i[, c(1,2)]
#}

#print(all_network_degree_matrix)

###注意上述循环报错如下:
#Error in `[<-`(`*tmp*`, i, , value = all_network_degree_i[, c(1, 2)]) :  subscript out of bounds

#经过chatGPT修改，发现all_network_degree_matrix[i, ] <- all_network_degree_i[, c(1,2)]应该是一个字符向量，而不是数字向量
##正确代码如下::


#已知网络节点列表时，对网络度分布进行统计，当未知时参见生信小白鱼的代码:200302-微生物共发生网络中节点度的幂律分布和在R中拟合
# Get unique labels
unique_labels <- unique(all_network_node$degree)
unique_labels

# Initialize empty matrix
all_network_degree_matrix <- matrix(0, nrow=length(unique_labels), ncol=2, dimnames=list(as.character(unique_labels), c("degree","count")))

# Loop through each label
for (i in unique_labels) {
  
  # Subset edges with the current label
  all_network_degree_count <- all_network_node[all_network_node$degree == i,] 
  dim(all_network_degree_count)
  
  # Create matrix with positive and negative degrees
  all_network_degree_i <- matrix(c(i, dim(all_network_degree_count)[1]), nrow=1, ncol=2, byrow=TRUE, dimnames=list(as.character(i),c("degree","count")))
  all_network_degree_i
  
  # Fill in matrix with values from current label
  all_network_degree_matrix[as.character(i), ] <- all_network_degree_i[, c(1,2)]
  #特别要注意，当unique函数的结果是字符时，i in unique()结果是字符，作为矩阵是不需要转换，但当时数字时，需要用as.character(i)来将数字i转化为字符i
}

print(all_network_degree_matrix)
all_network_degree_matrix <- as.data.frame(all_network_degree_matrix)
write.table(all_network_degree_matrix,paste0(output,"all_network_degree_matrix.txt"),sep="\t",quote=F)

#度分布统计
#degree_dist <- table(V(igraph)$degree)
degree_num <- as.numeric(all_network_degree_matrix$degree)
degree_count <- as.numeric(all_network_degree_matrix$count)

dat <- data.frame(degree = degree_num, count = degree_count)
head(dat)


#查看度分布
#可观察到微生物相关网络通常服从幂律分布
par(mfrow = c(1, 3))
hist(dat, xlab = 'Degree', ylab = 'Frequency', 
     main = 'Degree distribution')
plot(degree_count, degree_num, xlab = 'Degree', ylab = 'Count', 
     main = 'Degree distribution')
plot(degree_count, degree_num, log = 'xy', xlab = 'Log-degree', 
     ylab = 'Log-count', main = 'Log-log degree distribution')

##在R中，可通过非线性回归函数nls()拟合幂律分布。

#拟合，a 和 b 的初始值手动指定，数据不同需要多加尝试
mod <- nls(count ~ a*degree^b, data = dat, start = list(a = 2, b = 1.5))
summary(mod)

#提取关键值
a <- round(coef(mod)[1], 3)
b <- round(coef(mod)[2], 3)
a; b


#使用构建好的方程，通过预测变量拟合响应变量，再根据和观测值的差异，获得 R2
#SSre：the sum of the squares of the distances of the points from the fit
#SStot：the sum of the squares of the distances of the points from a horizontal line through the mean of all Y values
fit <- fitted(mod)
SSre <- sum((dat$count-fit)^2)
SStot <- sum((dat$count-mean(dat$count))^2)
R2 <- round(1 - SSre/SStot, 3)
R2


#p 值可以根据置换检验的原理获得
#将 count 的值分别随机置换 N 次（例如 999 次），通过随机置换数据后数据获取 R2（R2'）
#比较随机置换后值的 R2' 大于观测值的 R2 的频率，即为 p 值
p_num <- 1
dat_rand <- dat
for (i in 1:999) {
  dat_rand$count <- sample(dat_rand$count)
  SSre_rand <- sum((dat_rand$count-fit)^2)
  SStot_rand <- sum((dat_rand$count-mean(dat_rand$count))^2)
  R2_rand <- 1 - SSre_rand/SStot_rand
  if (R2_rand > R2) p_num <- p_num + 1
}
p_value <- p_num / (999+1)
p_value


#作图展示
library(ggplot2)

p <- ggplot(dat, aes(x = degree, y = count)) +
  geom_point(color = 'blue') +
  theme(panel.grid.major = element_line(color = 'gray'), panel.background = element_rect(color = 'black', fill = 'transparent')) +
  stat_smooth(method = 'nls', formula = y ~ a*x^b, method.args = list(start = list(a = 2, b = 1.5)), se = FALSE) +
  labs(x = 'Degree', y = 'Count')
p
#添加公式拟合的注释
label <- data.frame(
  formula = sprintf('italic(Y) == %.3f*italic(X)^%.3f', a, b),
  R2 = sprintf('italic(R^2) == %.3f', R2),
  p_value = sprintf('italic(P) < %.3f', p_value)
)
label


#注意：有关标签的x,y坐标需要根据第一次出图结果进行调整  
p2 <- p + geom_text(x = 60, y = 50, aes(label = formula), data = label, parse = TRUE, hjust = 0) +
  geom_text(x = 60, y = 40, aes(label = R2), data = label, parse = TRUE, hjust = 0) +
  geom_text(x = 60, y = 30, aes(label = p_value), data = label, parse = TRUE, hjust = 0)
p2

ggsave("all_network_degree.pdf", p2, width = 60, height = 50, units = "mm")
ggsave("all_network_degree.png", p2, width = 60, height = 50, units = "mm")
#ggsave("./H_zi_pi_p.tiff", H_zi_pi_p, width = 150, height = 100, units = "mm")








####Ggt_network####  
Ggt_network_node <- read.table("Ggt_network_node_list.txt", row.names=1,sep="\t", header=T, blank.lines.skip=F, check.names=F)
#读取数据的时候，header=T表示保留原始表格的表头.row.names=1表示保留数据表格的第一列，=2表示去除表格的第一列
Ggt_network_node

#已知网络节点列表时，对网络度分布进行统计，当未知时参见生信小白鱼的代码:200302-微生物共发生网络中节点度的幂律分布和在R中拟合
# Get unique labels
unique_labels <- unique(Ggt_network_node$degree)
unique_labels

# Initialize empty matrix
Ggt_network_degree_matrix <- matrix(0, nrow=length(unique_labels), ncol=2, dimnames=list(as.character(unique_labels), c("degree","count")))

# Loop through each label
for (i in unique_labels) {
  
  # Subset edges with the current label
  Ggt_network_degree_count <- Ggt_network_node[Ggt_network_node$degree == i,] 
  dim(Ggt_network_degree_count)
  
  # Create matrix with positive and negative degrees
  Ggt_network_degree_i <- matrix(c(i, dim(Ggt_network_degree_count)[1]), nrow=1, ncol=2, byrow=TRUE, dimnames=list(as.character(i),c("degree","count")))
  Ggt_network_degree_i
  
  # Fill in matrix with values from current label
  Ggt_network_degree_matrix[as.character(i), ] <- Ggt_network_degree_i[, c(1,2)]
  #特别要注意，当unique函数的结果是字符时，i in unique()结果是字符，作为矩阵是不需要转换，但当时数字时，需要用as.character(i)来将数字i转化为字符i
}

print(Ggt_network_degree_matrix)
Ggt_network_degree_matrix <- as.data.frame(Ggt_network_degree_matrix)
Ggt_network_degree_matrix

#结果中发现存在0值，需要剔除
Ggt_network_degree_matrix <- Ggt_network_degree_matrix[-8,]
Ggt_network_degree_matrix
write.table(Ggt_network_degree_matrix,paste0(output,"Ggt_network_degree_matrix.txt"),sep="\t",quote=F)


#度分布统计
#degree_dist <- table(V(igraph)$degree)
degree_num <- as.numeric(Ggt_network_degree_matrix$degree)
degree_count <- as.numeric(Ggt_network_degree_matrix$count)

dat <- data.frame(degree = degree_num, count = degree_count)
head(dat)


#查看度分布
#可观察到微生物相关网络通常服从幂律分布
par(mfrow = c(1, 3))
hist(dat, xlab = 'Degree', ylab = 'Frequency', 
     main = 'Degree distribution')
plot(degree_count, degree_num, xlab = 'Degree', ylab = 'Count', 
     main = 'Degree distribution')
plot(degree_count, degree_num, log = 'xy', xlab = 'Log-degree', 
     ylab = 'Log-count', main = 'Log-log degree distribution')

##在R中，可通过非线性回归函数nls()拟合幂律分布。

#拟合，a 和 b 的初始值手动指定，数据不同需要多加尝试
mod <- nls(count ~ a*degree^b, data = dat, start = list(a = 2, b = 1.5))
summary(mod)

#提取关键值
a <- round(coef(mod)[1], 3)
b <- round(coef(mod)[2], 3)
a; b


#使用构建好的方程，通过预测变量拟合响应变量，再根据和观测值的差异，获得 R2
#SSre：the sum of the squares of the distances of the points from the fit
#SStot：the sum of the squares of the distances of the points from a horizontal line through the mean of all Y values
fit <- fitted(mod)
SSre <- sum((dat$count-fit)^2)
SStot <- sum((dat$count-mean(dat$count))^2)
R2 <- round(1 - SSre/SStot, 3)
R2


#p 值可以根据置换检验的原理获得
#将 count 的值分别随机置换 N 次（例如 999 次），通过随机置换数据后数据获取 R2（R2'）
#比较随机置换后值的 R2' 大于观测值的 R2 的频率，即为 p 值
p_num <- 1
dat_rand <- dat
for (i in 1:999) {
  dat_rand$count <- sample(dat_rand$count)
  SSre_rand <- sum((dat_rand$count-fit)^2)
  SStot_rand <- sum((dat_rand$count-mean(dat_rand$count))^2)
  R2_rand <- 1 - SSre_rand/SStot_rand
  if (R2_rand > R2) p_num <- p_num + 1
}
p_value <- p_num / (999+1)
p_value


#作图展示
library(ggplot2)

Ggt_p <- ggplot(dat, aes(x = degree, y = count)) +
  geom_point(color = 'blue') +
  theme(panel.grid.major = element_line(color = 'gray'), panel.background = element_rect(color = 'black', fill = 'transparent')) +
  stat_smooth(method = 'nls', formula = y ~ a*x^b, method.args = list(start = list(a = 2, b = 1.5)), se = FALSE) +
  labs(x = 'Degree', y = 'Count')
Ggt_p
#添加公式拟合的注释
label <- data.frame(
  formula = sprintf('italic(Y) == %.3f*italic(X)^%.3f', a, b),
  R2 = sprintf('italic(R^2) == %.3f', R2),
  p_value = sprintf('italic(P) < %.3f', p_value)
)

Ggt_p2 <- Ggt_p + geom_text(x = 20, y = 13, aes(label = formula), data = label, parse = TRUE, hjust = 0) +
  geom_text(x = 20, y = 10, aes(label = R2), data = label, parse = TRUE, hjust = 0) +
  geom_text(x = 20, y = 7, aes(label = p_value), data = label, parse = TRUE, hjust = 0)

Ggt_p2

ggsave("Ggt_network_degree.pdf", Ggt_p2, width = 60, height = 50, units = "mm")
ggsave("Ggt_network_degree.png", Ggt_p2, width = 60, height = 50, units = "mm")
#ggsave("./H_zi_pi_p.tiff", H_zi_pi_p, width = 150, height = 100, units = "mm")




####S135_network####  
S135_network_node <- read.table("135_network_node_list.txt", row.names=1,sep="\t", header=T, blank.lines.skip=F, check.names=F)
#读取数据的时候，header=T表示保留原始表格的表头.row.names=1表示保留数据表格的第一列，=2表示去除表格的第一列
S135_network_node

#已知网络节点列表时，对网络度分布进行统计，当未知时参见生信小白鱼的代码:200302-微生物共发生网络中节点度的幂律分布和在R中拟合
# Get unique labels
unique_labels <- unique(S135_network_node$degree)
unique_labels

# Initialize empty matrix
S135_network_degree_matrix <- matrix(0, nrow=length(unique_labels), ncol=2, dimnames=list(as.character(unique_labels), c("degree","count")))

# Loop through each label
for (i in unique_labels) {
  
  # Subset edges with the current label
  S135_network_degree_count <- S135_network_node[S135_network_node$degree == i,] 
  dim(S135_network_degree_count)
  
  # Create matrix with positive and negative degrees
  S135_network_degree_i <- matrix(c(i, dim(S135_network_degree_count)[1]), nrow=1, ncol=2, byrow=TRUE, dimnames=list(as.character(i),c("degree","count")))
  S135_network_degree_i
  
  # Fill in matrix with values from current label
  S135_network_degree_matrix[as.character(i), ] <- S135_network_degree_i[, c(1,2)]
  #特别要注意，当unique函数的结果是字符时，i in unique()结果是字符，作为矩阵是不需要转换，但当时数字时，需要用as.character(i)来将数字i转化为字符i
}

print(S135_network_degree_matrix)
S135_network_degree_matrix <- as.data.frame(S135_network_degree_matrix)
rownames(S135_network_degree_matrix) <- NULL 
S135_network_degree_matrix

#结果中发现存在0值，需要剔除
S135_network_degree_matrix <- S135_network_degree_matrix[-40,]
S135_network_degree_matrix
write.table(S135_network_degree_matrix,paste0(output,"S135_network_degree_matrix.txt"),sep="\t",quote=F)


#度分布统计
#degree_dist <- table(V(igraph)$degree)
degree_num <- as.numeric(S135_network_degree_matrix$degree)
degree_count <- as.numeric(S135_network_degree_matrix$count)

dat <- data.frame(degree = degree_num, count = degree_count)
head(dat)


#查看度分布
#可观察到微生物相关网络通常服从幂律分布
par(mfrow = c(1, 3))
hist(dat, xlab = 'Degree', ylab = 'Frequency', 
     main = 'Degree distribution')
plot(degree_count, degree_num, xlab = 'Degree', ylab = 'Count', 
     main = 'Degree distribution')
plot(degree_count, degree_num, log = 'xy', xlab = 'Log-degree', 
     ylab = 'Log-count', main = 'Log-log degree distribution')

##在R中，可通过非线性回归函数nls()拟合幂律分布。

#拟合，a 和 b 的初始值手动指定，数据不同需要多加尝试
mod <- nls(count ~ a*degree^b, data = dat, start = list(a = 2, b = 1.5))
summary(mod)

#提取关键值
a <- round(coef(mod)[1], 3)
b <- round(coef(mod)[2], 3)
a; b


#使用构建好的方程，通过预测变量拟合响应变量，再根据和观测值的差异，获得 R2
#SSre：the sum of the squares of the distances of the points from the fit
#SStot：the sum of the squares of the distances of the points from a horizontal line through the mean of all Y values
fit <- fitted(mod)
SSre <- sum((dat$count-fit)^2)
SStot <- sum((dat$count-mean(dat$count))^2)
R2 <- round(1 - SSre/SStot, 3)
R2


#p 值可以根据置换检验的原理获得
#将 count 的值分别随机置换 N 次（例如 999 次），通过随机置换数据后数据获取 R2（R2'）
#比较随机置换后值的 R2' 大于观测值的 R2 的频率，即为 p 值
p_num <- 1
dat_rand <- dat
for (i in 1:999) {
  dat_rand$count <- sample(dat_rand$count)
  SSre_rand <- sum((dat_rand$count-fit)^2)
  SStot_rand <- sum((dat_rand$count-mean(dat_rand$count))^2)
  R2_rand <- 1 - SSre_rand/SStot_rand
  if (R2_rand > R2) p_num <- p_num + 1
}
p_value <- p_num / (999+1)
p_value


#作图展示
library(ggplot2)

S135_p <- ggplot(dat, aes(x = degree, y = count)) +
  geom_point(color = 'blue') +
  theme(panel.grid.major = element_line(color = 'gray'), panel.background = element_rect(color = 'black', fill = 'transparent')) +
  stat_smooth(method = 'nls', formula = y ~ a*x^b, method.args = list(start = list(a = 2, b = 1.5)), se = FALSE) +
  labs(x = 'Degree', y = 'Count')
S135_p
#添加公式拟合的注释
label <- data.frame(
  formula = sprintf('italic(Y) == %.3f*italic(X)^%.3f', a, b),
  R2 = sprintf('italic(R^2) == %.3f', R2),
  p_value = sprintf('italic(P) < %.3f', p_value)
)

S135_p2 <- S135_p + geom_text(x = 20, y = 10, aes(label = formula), data = label, parse = TRUE, hjust = 0) +
  geom_text(x = 20, y = 8, aes(label = R2), data = label, parse = TRUE, hjust = 0) +
  geom_text(x = 20, y = 6, aes(label = p_value), data = label, parse = TRUE, hjust = 0)

S135_p2

ggsave("S135_network_degree.pdf", S135_p2, width = 60, height = 50, units = "mm")
ggsave("S135_network_degree.png", S135_p2, width = 60, height = 50, units = "mm")
#ggsave("./H_zi_pi_p.tiff", H_zi_pi_p, width = 150, height = 100, units = "mm")





####H_network####  
H_network_node <- read.table("H_network_node_list.txt", row.names=1,sep="\t", header=T, blank.lines.skip=F, check.names=F)
#读取数据的时候，header=T表示保留原始表格的表头.row.names=1表示保留数据表格的第一列，=2表示去除表格的第一列
H_network_node

#已知网络节点列表时，对网络度分布进行统计，当未知时参见生信小白鱼的代码:200302-微生物共发生网络中节点度的幂律分布和在R中拟合
# Get unique labels
unique_labels <- unique(H_network_node$degree)
unique_labels

# Initialize empty matrix
H_network_degree_matrix <- matrix(0, nrow=length(unique_labels), ncol=2, dimnames=list(as.character(unique_labels), c("degree","count")))

# Loop through each label
for (i in unique_labels) {
  
  # Subset edges with the current label
  H_network_degree_count <- H_network_node[H_network_node$degree == i,] 
  dim(H_network_degree_count)
  
  # Create matrix with positive and negative degrees
  H_network_degree_i <- matrix(c(i, dim(H_network_degree_count)[1]), nrow=1, ncol=2, byrow=TRUE, dimnames=list(as.character(i),c("degree","count")))
  H_network_degree_i
  
  # Fill in matrix with values from current label
  H_network_degree_matrix[as.character(i), ] <- H_network_degree_i[, c(1,2)]
  #特别要注意，当unique函数的结果是字符时，i in unique()结果是字符，作为矩阵是不需要转换，但当时数字时，需要用as.character(i)来将数字i转化为字符i
}

print(H_network_degree_matrix)
H_network_degree_matrix <- as.data.frame(H_network_degree_matrix)
rownames(H_network_degree_matrix) <- NULL 
H_network_degree_matrix

#结果中发现存在0值，需要剔除
H_network_degree_matrix <- H_network_degree_matrix[-42,]
H_network_degree_matrix
write.table(H_network_degree_matrix,paste0(output,"H_network_degree_matrix.txt"),sep="\t",quote=F)


#度分布统计
#degree_dist <- table(V(igraph)$degree)
degree_num <- as.numeric(H_network_degree_matrix$degree)
degree_count <- as.numeric(H_network_degree_matrix$count)

dat <- data.frame(degree = degree_num, count = degree_count)
head(dat)


#查看度分布
#可观察到微生物相关网络通常服从幂律分布
par(mfrow = c(1, 3))
hist(dat, xlab = 'Degree', ylab = 'Frequency', 
     main = 'Degree distribution')
plot(degree_count, degree_num, xlab = 'Degree', ylab = 'Count', 
     main = 'Degree distribution')
plot(degree_count, degree_num, log = 'xy', xlab = 'Log-degree', 
     ylab = 'Log-count', main = 'Log-log degree distribution')

##在R中，可通过非线性回归函数nls()拟合幂律分布。

#拟合，a 和 b 的初始值手动指定，数据不同需要多加尝试
mod <- nls(count ~ a*degree^b, data = dat, start = list(a = 2, b = 1.5))
summary(mod)

#提取关键值
a <- round(coef(mod)[1], 3)
b <- round(coef(mod)[2], 3)
a; b


#使用构建好的方程，通过预测变量拟合响应变量，再根据和观测值的差异，获得 R2
#SSre：the sum of the squares of the distances of the points from the fit
#SStot：the sum of the squares of the distances of the points from a horizontal line through the mean of all Y values
fit <- fitted(mod)
SSre <- sum((dat$count-fit)^2)
SStot <- sum((dat$count-mean(dat$count))^2)
R2 <- round(1 - SSre/SStot, 3)
R2


#p 值可以根据置换检验的原理获得
#将 count 的值分别随机置换 N 次（例如 999 次），通过随机置换数据后数据获取 R2（R2'）
#比较随机置换后值的 R2' 大于观测值的 R2 的频率，即为 p 值
p_num <- 1
dat_rand <- dat
for (i in 1:999) {
  dat_rand$count <- sample(dat_rand$count)
  SSre_rand <- sum((dat_rand$count-fit)^2)
  SStot_rand <- sum((dat_rand$count-mean(dat_rand$count))^2)
  R2_rand <- 1 - SSre_rand/SStot_rand
  if (R2_rand > R2) p_num <- p_num + 1
}
p_value <- p_num / (999+1)
p_value


#作图展示
library(ggplot2)

H_p <- ggplot(dat, aes(x = degree, y = count)) +
  geom_point(color = 'blue') +
  theme(panel.grid.major = element_line(color = 'gray'), panel.background = element_rect(color = 'black', fill = 'transparent')) +
  stat_smooth(method = 'nls', formula = y ~ a*x^b, method.args = list(start = list(a = 2, b = 1.5)), se = FALSE) +
  labs(x = 'Degree', y = 'Count')
H_p
#添加公式拟合的注释
label <- data.frame(
  formula = sprintf('italic(Y) == %.3f*italic(X)^%.3f', a, b),
  R2 = sprintf('italic(R^2) == %.3f', R2),
  p_value = sprintf('italic(P) < %.3f', p_value)
)

H_p2 <- H_p + geom_text(x = 30, y = 12, aes(label = formula), data = label, parse = TRUE, hjust = 0) +
  geom_text(x = 30, y = 10, aes(label = R2), data = label, parse = TRUE, hjust = 0) +
  geom_text(x = 30, y = 8, aes(label = p_value), data = label, parse = TRUE, hjust = 0)

H_p2

ggsave("H_network_degree.pdf", H_p2, width = 60, height = 50, units = "mm")
ggsave("H_network_degree.png", H_p2, width = 60, height = 50, units = "mm")
#ggsave("./H_zi_pi_p.tiff", H_zi_pi_p, width = 150, height = 100, units = "mm")





####SA_network####  
SA_network_node <- read.table("SA_network_node_list.txt", row.names=1,sep="\t", header=T, blank.lines.skip=F, check.names=F)
#读取数据的时候，header=T表示保留原始表格的表头.row.names=1表示保留数据表格的第一列，=2表示去除表格的第一列
SA_network_node

#已知网络节点列表时，对网络度分布进行统计，当未知时参见生信小白鱼的代码:200302-微生物共发生网络中节点度的幂律分布和在R中拟合
# Get unique labels
unique_labels <- unique(SA_network_node$degree)
unique_labels

# Initialize empty matrix
SA_network_degree_matrix <- matrix(0, nrow=length(unique_labels), ncol=2, dimnames=list(as.character(unique_labels), c("degree","count")))

# Loop through each label
for (i in unique_labels) {
  
  # Subset edges with the current label
  SA_network_degree_count <- SA_network_node[SA_network_node$degree == i,] 
  dim(SA_network_degree_count)
  
  # Create matrix with positive and negative degrees
  SA_network_degree_i <- matrix(c(i, dim(SA_network_degree_count)[1]), nrow=1, ncol=2, byrow=TRUE, dimnames=list(as.character(i),c("degree","count")))
  SA_network_degree_i
  
  # Fill in matrix with values from current label
  SA_network_degree_matrix[as.character(i), ] <- SA_network_degree_i[, c(1,2)]
  #特别要注意，当unique函数的结果是字符时，i in unique()结果是字符，作为矩阵是不需要转换，但当时数字时，需要用as.character(i)来将数字i转化为字符i
}

print(SA_network_degree_matrix)
SA_network_degree_matrix <- as.data.frame(SA_network_degree_matrix)
rownames(SA_network_degree_matrix) <- NULL 
SA_network_degree_matrix

#结果中发现存在0值，需要剔除
SA_network_degree_matrix <- SA_network_degree_matrix[-78,]
SA_network_degree_matrix
write.table(SA_network_degree_matrix,paste0(output,"SA_network_degree_matrix.txt"),sep="\t",quote=F)


#度分布统计
#degree_dist <- table(V(igraph)$degree)
degree_num <- as.numeric(SA_network_degree_matrix$degree)
degree_count <- as.numeric(SA_network_degree_matrix$count)

dat <- data.frame(degree = degree_num, count = degree_count)
head(dat)


#查看度分布
#可观察到微生物相关网络通常服从幂律分布
par(mfrow = c(1, 3))
hist(dat, xlab = 'Degree', ylab = 'Frequency', 
     main = 'Degree distribution')
plot(degree_count, degree_num, xlab = 'Degree', ylab = 'Count', 
     main = 'Degree distribution')
plot(degree_count, degree_num, log = 'xy', xlab = 'Log-degree', 
     ylab = 'Log-count', main = 'Log-log degree distribution')

##在R中，可通过非线性回归函数nls()拟合幂律分布。

#拟合，a 和 b 的初始值手动指定，数据不同需要多加尝试
mod <- nls(count ~ a*degree^b, data = dat, start = list(a = 2, b = 1.5))
summary(mod)

#提取关键值
a <- round(coef(mod)[1], 3)
b <- round(coef(mod)[2], 3)
a; b


#使用构建好的方程，通过预测变量拟合响应变量，再根据和观测值的差异，获得 R2
#SSre：the sum of the squares of the distances of the points from the fit
#SStot：the sum of the squares of the distances of the points from a horizontal line through the mean of all Y values
fit <- fitted(mod)
SSre <- sum((dat$count-fit)^2)
SStot <- sum((dat$count-mean(dat$count))^2)
R2 <- round(1 - SSre/SStot, 3)
R2


#p 值可以根据置换检验的原理获得
#将 count 的值分别随机置换 N 次（例如 999 次），通过随机置换数据后数据获取 R2（R2'）
#比较随机置换后值的 R2' 大于观测值的 R2 的频率，即为 p 值
p_num <- 1
dat_rand <- dat
for (i in 1:999) {
  dat_rand$count <- sample(dat_rand$count)
  SSre_rand <- sum((dat_rand$count-fit)^2)
  SStot_rand <- sum((dat_rand$count-mean(dat_rand$count))^2)
  R2_rand <- 1 - SSre_rand/SStot_rand
  if (R2_rand > R2) p_num <- p_num + 1
}
p_value <- p_num / (999+1)
p_value


#作图展示
library(ggplot2)

SA_p <- ggplot(dat, aes(x = degree, y = count)) +
  geom_point(color = 'blue') +
  theme(panel.grid.major = element_line(color = 'gray'), panel.background = element_rect(color = 'black', fill = 'transparent')) +
  stat_smooth(method = 'nls', formula = y ~ a*x^b, method.args = list(start = list(a = 2, b = 1.5)), se = FALSE) +
  labs(x = 'Degree', y = 'Count')
SA_p
#添加公式拟合的注释
label <- data.frame(
  formula = sprintf('italic(Y) == %.3f*italic(X)^%.3f', a, b),
  R2 = sprintf('italic(R^2) == %.3f', R2),
  p_value = sprintf('italic(P) < %.3f', p_value)
)

SA_p2 <- H_p + geom_text(x = 30, y = 12, aes(label = formula), data = label, parse = TRUE, hjust = 0) +
  geom_text(x = 30, y = 10, aes(label = R2), data = label, parse = TRUE, hjust = 0) +
  geom_text(x = 30, y = 8, aes(label = p_value), data = label, parse = TRUE, hjust = 0)

SA_p2

ggsave("SA_network_degree.pdf", SA_p2, width = 60, height = 50, units = "mm")
ggsave("SA_network_degree.png", SA_p2, width = 60, height = 50, units = "mm")
#ggsave("./H_zi_pi_p.tiff", H_zi_pi_p, width = 150, height = 100, units = "mm")











####module0####  
module0_node <- read.table("All_module0_gephi2.txt", row.names=1,sep="\t", header=T, blank.lines.skip=F, check.names=F)
#读取数据的时候，header=T表示保留原始表格的表头.row.names=1表示保留数据表格的第一列，=2表示去除表格的第一列
module0_node

#已知网络节点列表时，对网络度分布进行统计，当未知时参见生信小白鱼的代码:200302-微生物共发生网络中节点度的幂律分布和在R中拟合
# Get unique labels
unique_labels <- unique(module0_node$degree)
unique_labels

# Initialize empty matrix
module0_degree_matrix <- matrix(0, nrow=length(unique_labels), ncol=2, dimnames=list(as.character(unique_labels), c("degree","count")))

# Loop through each label
for (i in unique_labels) {
  
  # Subset edges with the current label
  module0_degree_count <- module0_node[module0_node$degree == i,] 
  dim(module0_degree_count)
  
  # Create matrix with positive and negative degrees
  module0_degree_i <- matrix(c(i, dim(module0_degree_count)[1]), nrow=1, ncol=2, byrow=TRUE, dimnames=list(as.character(i),c("degree","count")))
  module0_degree_i
  
  # Fill in matrix with values from current label
  module0_degree_matrix[as.character(i), ] <- module0_degree_i[, c(1,2)]
  #特别要注意，当unique函数的结果是字符时，i in unique()结果是字符，作为矩阵是不需要转换，但当时数字时，需要用as.character(i)来将数字i转化为字符i
}

print(module0_degree_matrix)
module0_degree_matrix <- as.data.frame(module0_degree_matrix)
rownames(module0_degree_matrix) <- NULL 
module0_degree_matrix

#结果中发现存在0值，需要剔除
#module0_degree_matrix <- module0_degree_matrix[-42,]
#module0_degree_matrix
write.table(module0_degree_matrix,paste0(output,"module0_degree_matrix.txt"),sep="\t",quote=F)


#度分布统计
#degree_dist <- table(V(igraph)$degree)
degree_num <- as.numeric(module0_degree_matrix$degree)
degree_count <- as.numeric(module0_degree_matrix$count)

dat <- data.frame(degree = degree_num, count = degree_count)
head(dat)


#查看度分布
#可观察到微生物相关网络通常服从幂律分布
par(mfrow = c(1, 3))
hist(dat, xlab = 'Degree', ylab = 'Frequency', 
     main = 'Degree distribution')
plot(degree_count, degree_num, xlab = 'Degree', ylab = 'Count', 
     main = 'Degree distribution')
plot(degree_count, degree_num, log = 'xy', xlab = 'Log-degree', 
     ylab = 'Log-count', main = 'Log-log degree distribution')

##在R中，可通过非线性回归函数nls()拟合幂律分布。

#拟合，a 和 b 的初始值手动指定，数据不同需要多加尝试
mod <- nls(count ~ a*degree^b, data = dat, start = list(a = 2, b = 1.5))
summary(mod)

#提取关键值
a <- round(coef(mod)[1], 3)
b <- round(coef(mod)[2], 3)
a; b


#使用构建好的方程，通过预测变量拟合响应变量，再根据和观测值的差异，获得 R2
#SSre：the sum of the squares of the distances of the points from the fit
#SStot：the sum of the squares of the distances of the points from a horizontal line through the mean of all Y values
fit <- fitted(mod)
SSre <- sum((dat$count-fit)^2)
SStot <- sum((dat$count-mean(dat$count))^2)
R2 <- round(1 - SSre/SStot, 3)
R2


#p 值可以根据置换检验的原理获得
#将 count 的值分别随机置换 N 次（例如 999 次），通过随机置换数据后数据获取 R2（R2'）
#比较随机置换后值的 R2' 大于观测值的 R2 的频率，即为 p 值
p_num <- 1
dat_rand <- dat
for (i in 1:999) {
  dat_rand$count <- sample(dat_rand$count)
  SSre_rand <- sum((dat_rand$count-fit)^2)
  SStot_rand <- sum((dat_rand$count-mean(dat_rand$count))^2)
  R2_rand <- 1 - SSre_rand/SStot_rand
  if (R2_rand > R2) p_num <- p_num + 1
}
p_value <- p_num / (999+1)
p_value


#作图展示
library(ggplot2)

module0_p <- ggplot(dat, aes(x = degree, y = count)) +
  geom_point(color = 'blue') +
  theme(panel.grid.major = element_line(color = 'gray'), panel.background = element_rect(color = 'black', fill = 'transparent')) +
  stat_smooth(method = 'nls', formula = y ~ a*x^b, method.args = list(start = list(a = 2, b = 1.5)), se = FALSE) +
  labs(x = 'Degree', y = 'Count')
module0_p
#添加公式拟合的注释
label <- data.frame(
  formula = sprintf('italic(Y) == %.3f*italic(X)^%.3f', a, b),
  R2 = sprintf('italic(R^2) == %.3f', R2),
  p_value = sprintf('italic(P) < %.3f', p_value)
)

module0_p2 <- module0_p + geom_text(x = 30, y = 5, aes(label = formula), data = label, parse = TRUE, hjust = 0) +
  geom_text(x = 30, y = 4, aes(label = R2), data = label, parse = TRUE, hjust = 0) +
  geom_text(x = 30, y = 3, aes(label = p_value), data = label, parse = TRUE, hjust = 0)

module0_p2

ggsave("module0_degree.pdf", module0_p2, width = 60, height = 50, units = "mm")
ggsave("module0_degree.png", module0_p2, width = 60, height = 50, units = "mm")
#ggsave("./H_zi_pi_p.tiff", H_zi_pi_p, width = 150, height = 100, units = "mm")





####module1####  
module1_node <- read.table("All_module1_gephi2.txt", row.names=1,sep="\t", header=T, blank.lines.skip=F, check.names=F)
#读取数据的时候，header=T表示保留原始表格的表头.row.names=1表示保留数据表格的第一列，=2表示去除表格的第一列
module1_node

#已知网络节点列表时，对网络度分布进行统计，当未知时参见生信小白鱼的代码:200302-微生物共发生网络中节点度的幂律分布和在R中拟合
# Get unique labels
unique_labels <- unique(module1_node$degree)
unique_labels

# Initialize empty matrix
module1_degree_matrix <- matrix(0, nrow=length(unique_labels), ncol=2, dimnames=list(as.character(unique_labels), c("degree","count")))

# Loop through each label
for (i in unique_labels) {
  
  # Subset edges with the current label
  module1_degree_count <- module1_node[module1_node$degree == i,] 
  dim(module1_degree_count)
  
  # Create matrix with positive and negative degrees
  module1_degree_i <- matrix(c(i, dim(module1_degree_count)[1]), nrow=1, ncol=2, byrow=TRUE, dimnames=list(as.character(i),c("degree","count")))
  module1_degree_i
  
  # Fill in matrix with values from current label
  module1_degree_matrix[as.character(i), ] <- module1_degree_i[, c(1,2)]
  #特别要注意，当unique函数的结果是字符时，i in unique()结果是字符，作为矩阵是不需要转换，但当时数字时，需要用as.character(i)来将数字i转化为字符i
}

print(module1_degree_matrix)
module1_degree_matrix <- as.data.frame(module1_degree_matrix)
rownames(module1_degree_matrix) <- NULL 
module1_degree_matrix

#结果中发现存在0值，需要剔除
#module0_degree_matrix <- module0_degree_matrix[-42,]
#module0_degree_matrix
write.table(module1_degree_matrix,paste0(output,"module1_degree_matrix.txt"),sep="\t",quote=F)


#度分布统计
#degree_dist <- table(V(igraph)$degree)
degree_num <- as.numeric(module1_degree_matrix$degree)
degree_count <- as.numeric(module1_degree_matrix$count)

dat <- data.frame(degree = degree_num, count = degree_count)
head(dat)


#查看度分布
#可观察到微生物相关网络通常服从幂律分布
par(mfrow = c(1, 3))
hist(dat, xlab = 'Degree', ylab = 'Frequency', 
     main = 'Degree distribution')
plot(degree_count, degree_num, xlab = 'Degree', ylab = 'Count', 
     main = 'Degree distribution')
plot(degree_count, degree_num, log = 'xy', xlab = 'Log-degree', 
     ylab = 'Log-count', main = 'Log-log degree distribution')

##在R中，可通过非线性回归函数nls()拟合幂律分布。

#拟合，a 和 b 的初始值手动指定，数据不同需要多加尝试
mod <- nls(count ~ a*degree^b, data = dat, start = list(a = 2, b = 1.5))
summary(mod)

#提取关键值
a <- round(coef(mod)[1], 3)
b <- round(coef(mod)[2], 3)
a; b


#使用构建好的方程，通过预测变量拟合响应变量，再根据和观测值的差异，获得 R2
#SSre：the sum of the squares of the distances of the points from the fit
#SStot：the sum of the squares of the distances of the points from a horizontal line through the mean of all Y values
fit <- fitted(mod)
SSre <- sum((dat$count-fit)^2)
SStot <- sum((dat$count-mean(dat$count))^2)
R2 <- round(1 - SSre/SStot, 3)
R2


#p 值可以根据置换检验的原理获得
#将 count 的值分别随机置换 N 次（例如 999 次），通过随机置换数据后数据获取 R2（R2'）
#比较随机置换后值的 R2' 大于观测值的 R2 的频率，即为 p 值
p_num <- 1
dat_rand <- dat
for (i in 1:999) {
  dat_rand$count <- sample(dat_rand$count)
  SSre_rand <- sum((dat_rand$count-fit)^2)
  SStot_rand <- sum((dat_rand$count-mean(dat_rand$count))^2)
  R2_rand <- 1 - SSre_rand/SStot_rand
  if (R2_rand > R2) p_num <- p_num + 1
}
p_value <- p_num / (999+1)
p_value


#作图展示
library(ggplot2)

module1_p <- ggplot(dat, aes(x = degree, y = count)) +
  geom_point(color = 'blue') +
  theme(panel.grid.major = element_line(color = 'gray'), panel.background = element_rect(color = 'black', fill = 'transparent')) +
  stat_smooth(method = 'nls', formula = y ~ a*x^b, method.args = list(start = list(a = 2, b = 1.5)), se = FALSE) +
  labs(x = 'Degree', y = 'Count')
module1_p
#添加公式拟合的注释
label <- data.frame(
  formula = sprintf('italic(Y) == %.3f*italic(X)^%.3f', a, b),
  R2 = sprintf('italic(R^2) == %.3f', R2),
  p_value = sprintf('italic(P) < %.3f', p_value)
)

module1_p2 <- module1_p + geom_text(x = 50, y = 6, aes(label = formula), data = label, parse = TRUE, hjust = 0) +
  geom_text(x = 50, y = 5, aes(label = R2), data = label, parse = TRUE, hjust = 0) +
  geom_text(x = 50, y = 4, aes(label = p_value), data = label, parse = TRUE, hjust = 0)

module1_p2

ggsave("module1_degree.pdf", module1_p2, width = 60, height = 50, units = "mm")
ggsave("module1_degree.png", module1_p2, width = 60, height = 50, units = "mm")
#ggsave("./H_zi_pi_p.tiff", H_zi_pi_p, width = 150, height = 100, units = "mm")





####module2####  
module2_node <- read.table("All_module2_gephi2.txt", row.names=1,sep="\t", header=T, blank.lines.skip=F, check.names=F)
#读取数据的时候，header=T表示保留原始表格的表头.row.names=1表示保留数据表格的第一列，=2表示去除表格的第一列
module2_node

#已知网络节点列表时，对网络度分布进行统计，当未知时参见生信小白鱼的代码:200302-微生物共发生网络中节点度的幂律分布和在R中拟合
# Get unique labels
unique_labels <- unique(module2_node$degree)
unique_labels

# Initialize empty matrix
module2_degree_matrix <- matrix(0, nrow=length(unique_labels), ncol=2, dimnames=list(as.character(unique_labels), c("degree","count")))

# Loop through each label
for (i in unique_labels) {
  
  # Subset edges with the current label
  module2_degree_count <- module2_node[module2_node$degree == i,] 
  dim(module2_degree_count)
  
  # Create matrix with positive and negative degrees
  module2_degree_i <- matrix(c(i, dim(module2_degree_count)[1]), nrow=1, ncol=2, byrow=TRUE, dimnames=list(as.character(i),c("degree","count")))
  module2_degree_i
  
  # Fill in matrix with values from current label
  module2_degree_matrix[as.character(i), ] <- module2_degree_i[, c(1,2)]
  #特别要注意，当unique函数的结果是字符时，i in unique()结果是字符，作为矩阵是不需要转换，但当时数字时，需要用as.character(i)来将数字i转化为字符i
}

print(module2_degree_matrix)
module2_degree_matrix <- as.data.frame(module2_degree_matrix)
rownames(module2_degree_matrix) <- NULL 
module2_degree_matrix

#结果中发现存在0值，需要剔除
#module0_degree_matrix <- module0_degree_matrix[-42,]
#module0_degree_matrix
write.table(module2_degree_matrix,paste0(output,"module2_degree_matrix.txt"),sep="\t",quote=F)


#度分布统计
#degree_dist <- table(V(igraph)$degree)
degree_num <- as.numeric(module2_degree_matrix$degree)
degree_count <- as.numeric(module2_degree_matrix$count)

dat <- data.frame(degree = degree_num, count = degree_count)
head(dat)


#查看度分布
#可观察到微生物相关网络通常服从幂律分布
par(mfrow = c(1, 3))
hist(dat, xlab = 'Degree', ylab = 'Frequency', 
     main = 'Degree distribution')
plot(degree_count, degree_num, xlab = 'Degree', ylab = 'Count', 
     main = 'Degree distribution')
plot(degree_count, degree_num, log = 'xy', xlab = 'Log-degree', 
     ylab = 'Log-count', main = 'Log-log degree distribution')

##在R中，可通过非线性回归函数nls()拟合幂律分布。

#拟合，a 和 b 的初始值手动指定，数据不同需要多加尝试
mod <- nls(count ~ a*degree^b, data = dat, start = list(a = 2, b = 1.5))
summary(mod)

#提取关键值
a <- round(coef(mod)[1], 3)
b <- round(coef(mod)[2], 3)
a; b


#使用构建好的方程，通过预测变量拟合响应变量，再根据和观测值的差异，获得 R2
#SSre：the sum of the squares of the distances of the points from the fit
#SStot：the sum of the squares of the distances of the points from a horizontal line through the mean of all Y values
fit <- fitted(mod)
SSre <- sum((dat$count-fit)^2)
SStot <- sum((dat$count-mean(dat$count))^2)
R2 <- round(1 - SSre/SStot, 3)
R2


#p 值可以根据置换检验的原理获得
#将 count 的值分别随机置换 N 次（例如 999 次），通过随机置换数据后数据获取 R2（R2'）
#比较随机置换后值的 R2' 大于观测值的 R2 的频率，即为 p 值
p_num <- 1
dat_rand <- dat
for (i in 1:999) {
  dat_rand$count <- sample(dat_rand$count)
  SSre_rand <- sum((dat_rand$count-fit)^2)
  SStot_rand <- sum((dat_rand$count-mean(dat_rand$count))^2)
  R2_rand <- 1 - SSre_rand/SStot_rand
  if (R2_rand > R2) p_num <- p_num + 1
}
p_value <- p_num / (999+1)
p_value


#作图展示
library(ggplot2)

module2_p <- ggplot(dat, aes(x = degree, y = count)) +
  geom_point(color = 'blue') +
  theme(panel.grid.major = element_line(color = 'gray'), panel.background = element_rect(color = 'black', fill = 'transparent')) +
  stat_smooth(method = 'nls', formula = y ~ a*x^b, method.args = list(start = list(a = 2, b = 1.5)), se = FALSE) +
  labs(x = 'Degree', y = 'Count')
module2_p
#添加公式拟合的注释
label <- data.frame(
  formula = sprintf('italic(Y) == %.3f*italic(X)^%.3f', a, b),
  R2 = sprintf('italic(R^2) == %.3f', R2),
  p_value = sprintf('italic(P) < %.3f', p_value)
)

module2_p2 <- module2_p + geom_text(x = 30, y = 10, aes(label = formula), data = label, parse = TRUE, hjust = 0) +
  geom_text(x = 30, y = 8, aes(label = R2), data = label, parse = TRUE, hjust = 0) +
  geom_text(x = 30, y = 6, aes(label = p_value), data = label, parse = TRUE, hjust = 0)

module2_p2

ggsave("module2_degree.pdf", module2_p2, width = 60, height = 50, units = "mm")
ggsave("module2_degree.png", module2_p2, width = 60, height = 50, units = "mm")
#ggsave("./H_zi_pi_p.tiff", H_zi_pi_p, width = 150, height = 100, units = "mm")



####module3####  
module3_node <- read.table("All_module3_gephi2.txt", row.names=1,sep="\t", header=T, blank.lines.skip=F, check.names=F)
#读取数据的时候，header=T表示保留原始表格的表头.row.names=1表示保留数据表格的第一列，=2表示去除表格的第一列
module3_node

#已知网络节点列表时，对网络度分布进行统计，当未知时参见生信小白鱼的代码:200302-微生物共发生网络中节点度的幂律分布和在R中拟合
# Get unique labels
unique_labels <- unique(module3_node$degree)
unique_labels

# Initialize empty matrix
module3_degree_matrix <- matrix(0, nrow=length(unique_labels), ncol=2, dimnames=list(as.character(unique_labels), c("degree","count")))

# Loop through each label
for (i in unique_labels) {
  
  # Subset edges with the current label
  module3_degree_count <- module3_node[module3_node$degree == i,] 
  dim(module3_degree_count)
  
  # Create matrix with positive and negative degrees
  module3_degree_i <- matrix(c(i, dim(module3_degree_count)[1]), nrow=1, ncol=2, byrow=TRUE, dimnames=list(as.character(i),c("degree","count")))
  module3_degree_i
  
  # Fill in matrix with values from current label
  module3_degree_matrix[as.character(i), ] <- module3_degree_i[, c(1,2)]
  #特别要注意，当unique函数的结果是字符时，i in unique()结果是字符，作为矩阵是不需要转换，但当时数字时，需要用as.character(i)来将数字i转化为字符i
}

print(module3_degree_matrix)
module3_degree_matrix <- as.data.frame(module3_degree_matrix)
rownames(module3_degree_matrix) <- NULL 
module3_degree_matrix

#结果中发现存在0值，需要剔除
#module0_degree_matrix <- module0_degree_matrix[-42,]
#module0_degree_matrix
write.table(module3_degree_matrix,paste0(output,"module3_degree_matrix.txt"),sep="\t",quote=F)


#度分布统计
#degree_dist <- table(V(igraph)$degree)
degree_num <- as.numeric(module3_degree_matrix$degree)
degree_count <- as.numeric(module3_degree_matrix$count)

dat <- data.frame(degree = degree_num, count = degree_count)
head(dat)


#查看度分布
#可观察到微生物相关网络通常服从幂律分布
par(mfrow = c(1, 3))
hist(dat, xlab = 'Degree', ylab = 'Frequency', 
     main = 'Degree distribution')
plot(degree_count, degree_num, xlab = 'Degree', ylab = 'Count', 
     main = 'Degree distribution')
plot(degree_count, degree_num, log = 'xy', xlab = 'Log-degree', 
     ylab = 'Log-count', main = 'Log-log degree distribution')

##在R中，可通过非线性回归函数nls()拟合幂律分布。

#拟合，a 和 b 的初始值手动指定，数据不同需要多加尝试
mod <- nls(count ~ a*degree^b, data = dat, start = list(a = 2, b = 1.5))
summary(mod)

#提取关键值
a <- round(coef(mod)[1], 3)
b <- round(coef(mod)[2], 3)
a; b


#使用构建好的方程，通过预测变量拟合响应变量，再根据和观测值的差异，获得 R2
#SSre：the sum of the squares of the distances of the points from the fit
#SStot：the sum of the squares of the distances of the points from a horizontal line through the mean of all Y values
fit <- fitted(mod)
SSre <- sum((dat$count-fit)^2)
SStot <- sum((dat$count-mean(dat$count))^2)
R2 <- round(1 - SSre/SStot, 3)
R2


#p 值可以根据置换检验的原理获得
#将 count 的值分别随机置换 N 次（例如 999 次），通过随机置换数据后数据获取 R2（R2'）
#比较随机置换后值的 R2' 大于观测值的 R2 的频率，即为 p 值
p_num <- 1
dat_rand <- dat
for (i in 1:999) {
  dat_rand$count <- sample(dat_rand$count)
  SSre_rand <- sum((dat_rand$count-fit)^2)
  SStot_rand <- sum((dat_rand$count-mean(dat_rand$count))^2)
  R2_rand <- 1 - SSre_rand/SStot_rand
  if (R2_rand > R2) p_num <- p_num + 1
}
p_value <- p_num / (999+1)
p_value


#作图展示
library(ggplot2)

module3_p <- ggplot(dat, aes(x = degree, y = count)) +
  geom_point(color = 'blue') +
  theme(panel.grid.major = element_line(color = 'gray'), panel.background = element_rect(color = 'black', fill = 'transparent')) +
  stat_smooth(method = 'nls', formula = y ~ a*x^b, method.args = list(start = list(a = 2, b = 1.5)), se = FALSE) +
  labs(x = 'Degree', y = 'Count')
module3_p
#添加公式拟合的注释
label <- data.frame(
  formula = sprintf('italic(Y) == %.3f*italic(X)^%.3f', a, b),
  R2 = sprintf('italic(R^2) == %.3f', R2),
  p_value = sprintf('italic(P) < %.3f', p_value)
)

module3_p2 <- module3_p + geom_text(x = 30, y = 8, aes(label = formula), data = label, parse = TRUE, hjust = 0) +
  geom_text(x = 30, y = 6, aes(label = R2), data = label, parse = TRUE, hjust = 0) +
  geom_text(x = 30, y = 4, aes(label = p_value), data = label, parse = TRUE, hjust = 0)

module3_p2

ggsave("module3_degree.pdf", module3_p2, width = 60, height = 50, units = "mm")
ggsave("module3_degree.png", module3_p2, width = 60, height = 50, units = "mm")
#ggsave("./H_zi_pi_p.tiff", H_zi_pi_p, width = 150, height = 100, units = "mm")




####module4####  
module4_node <- read.table("All_module4_gephi2.txt", row.names=1,sep="\t", header=T, blank.lines.skip=F, check.names=F)
#读取数据的时候，header=T表示保留原始表格的表头.row.names=1表示保留数据表格的第一列，=2表示去除表格的第一列
module4_node

#已知网络节点列表时，对网络度分布进行统计，当未知时参见生信小白鱼的代码:200302-微生物共发生网络中节点度的幂律分布和在R中拟合
# Get unique labels
unique_labels <- unique(module4_node$degree)
unique_labels

# Initialize empty matrix
module4_degree_matrix <- matrix(0, nrow=length(unique_labels), ncol=2, dimnames=list(as.character(unique_labels), c("degree","count")))

# Loop through each label
for (i in unique_labels) {
  
  # Subset edges with the current label
  module4_degree_count <- module4_node[module4_node$degree == i,] 
  dim(module4_degree_count)
  
  # Create matrix with positive and negative degrees
  module4_degree_i <- matrix(c(i, dim(module4_degree_count)[1]), nrow=1, ncol=2, byrow=TRUE, dimnames=list(as.character(i),c("degree","count")))
  module4_degree_i
  
  # Fill in matrix with values from current label
  module4_degree_matrix[as.character(i), ] <- module4_degree_i[, c(1,2)]
  #特别要注意，当unique函数的结果是字符时，i in unique()结果是字符，作为矩阵是不需要转换，但当时数字时，需要用as.character(i)来将数字i转化为字符i
}

print(module4_degree_matrix)
module4_degree_matrix <- as.data.frame(module4_degree_matrix)
rownames(module4_degree_matrix) <- NULL 
module4_degree_matrix

#结果中发现存在0值，需要剔除
#module0_degree_matrix <- module0_degree_matrix[-42,]
#module0_degree_matrix
write.table(module4_degree_matrix,paste0(output,"module4_degree_matrix.txt"),sep="\t",quote=F)


#度分布统计
#degree_dist <- table(V(igraph)$degree)
degree_num <- as.numeric(module4_degree_matrix$degree)
degree_count <- as.numeric(module4_degree_matrix$count)

dat <- data.frame(degree = degree_num, count = degree_count)
head(dat)


#查看度分布
#可观察到微生物相关网络通常服从幂律分布
par(mfrow = c(1, 3))
hist(dat, xlab = 'Degree', ylab = 'Frequency', 
     main = 'Degree distribution')
plot(degree_count, degree_num, xlab = 'Degree', ylab = 'Count', 
     main = 'Degree distribution')
plot(degree_count, degree_num, log = 'xy', xlab = 'Log-degree', 
     ylab = 'Log-count', main = 'Log-log degree distribution')

##在R中，可通过非线性回归函数nls()拟合幂律分布。

#拟合，a 和 b 的初始值手动指定，数据不同需要多加尝试
mod <- nls(count ~ a*degree^b, data = dat, start = list(a = 2, b = 1.5))
summary(mod)

#提取关键值
a <- round(coef(mod)[1], 3)
b <- round(coef(mod)[2], 3)
a; b


#使用构建好的方程，通过预测变量拟合响应变量，再根据和观测值的差异，获得 R2
#SSre：the sum of the squares of the distances of the points from the fit
#SStot：the sum of the squares of the distances of the points from a horizontal line through the mean of all Y values
fit <- fitted(mod)
SSre <- sum((dat$count-fit)^2)
SStot <- sum((dat$count-mean(dat$count))^2)
R2 <- round(1 - SSre/SStot, 3)
R2


#p 值可以根据置换检验的原理获得
#将 count 的值分别随机置换 N 次（例如 999 次），通过随机置换数据后数据获取 R2（R2'）
#比较随机置换后值的 R2' 大于观测值的 R2 的频率，即为 p 值
p_num <- 1
dat_rand <- dat
for (i in 1:999) {
  dat_rand$count <- sample(dat_rand$count)
  SSre_rand <- sum((dat_rand$count-fit)^2)
  SStot_rand <- sum((dat_rand$count-mean(dat_rand$count))^2)
  R2_rand <- 1 - SSre_rand/SStot_rand
  if (R2_rand > R2) p_num <- p_num + 1
}
p_value <- p_num / (999+1)
p_value


#作图展示
library(ggplot2)

module4_p <- ggplot(dat, aes(x = degree, y = count)) +
  geom_point(color = 'blue') +
  theme(panel.grid.major = element_line(color = 'gray'), panel.background = element_rect(color = 'black', fill = 'transparent')) +
  stat_smooth(method = 'nls', formula = y ~ a*x^b, method.args = list(start = list(a = 2, b = 1.5)), se = FALSE) +
  labs(x = 'Degree', y = 'Count')
module4_p
#添加公式拟合的注释
label <- data.frame(
  formula = sprintf('italic(Y) == %.3f*italic(X)^%.3f', a, b),
  R2 = sprintf('italic(R^2) == %.3f', R2),
  p_value = sprintf('italic(P) < %.3f', p_value)
)

module4_p2 <- module4_p + geom_text(x = 20, y = 6.5, aes(label = formula), data = label, parse = TRUE, hjust = 0) +
  geom_text(x = 20, y = 5, aes(label = R2), data = label, parse = TRUE, hjust = 0) +
  geom_text(x = 20, y = 3, aes(label = p_value), data = label, parse = TRUE, hjust = 0)

module4_p2

ggsave("module4_degree.pdf", module4_p2, width = 60, height = 50, units = "mm")
ggsave("module4_degree.png", module4_p2, width = 60, height = 50, units = "mm")
#ggsave("./H_zi_pi_p.tiff", H_zi_pi_p, width = 150, height = 100, units = "mm")





####module5####  
module5_node <- read.table("All_module5_gephi2.txt", row.names=1,sep="\t", header=T, blank.lines.skip=F, check.names=F)
#读取数据的时候，header=T表示保留原始表格的表头.row.names=1表示保留数据表格的第一列，=2表示去除表格的第一列
module5_node

#已知网络节点列表时，对网络度分布进行统计，当未知时参见生信小白鱼的代码:200302-微生物共发生网络中节点度的幂律分布和在R中拟合
# Get unique labels
unique_labels <- unique(module5_node$degree)
unique_labels

# Initialize empty matrix
module5_degree_matrix <- matrix(0, nrow=length(unique_labels), ncol=2, dimnames=list(as.character(unique_labels), c("degree","count")))

# Loop through each label
for (i in unique_labels) {
  
  # Subset edges with the current label
  module5_degree_count <- module5_node[module5_node$degree == i,] 
  dim(module5_degree_count)
  
  # Create matrix with positive and negative degrees
  module5_degree_i <- matrix(c(i, dim(module5_degree_count)[1]), nrow=1, ncol=2, byrow=TRUE, dimnames=list(as.character(i),c("degree","count")))
  module5_degree_i
  
  # Fill in matrix with values from current label
  module5_degree_matrix[as.character(i), ] <- module5_degree_i[, c(1,2)]
  #特别要注意，当unique函数的结果是字符时，i in unique()结果是字符，作为矩阵是不需要转换，但当时数字时，需要用as.character(i)来将数字i转化为字符i
}

print(module5_degree_matrix)
module5_degree_matrix <- as.data.frame(module5_degree_matrix)
rownames(module5_degree_matrix) <- NULL 
module5_degree_matrix

#结果中发现存在0值，需要剔除
#module0_degree_matrix <- module0_degree_matrix[-42,]
#module0_degree_matrix
write.table(module5_degree_matrix,paste0(output,"module5_degree_matrix.txt"),sep="\t",quote=F)


#度分布统计
#degree_dist <- table(V(igraph)$degree)
degree_num <- as.numeric(module5_degree_matrix$degree)
degree_count <- as.numeric(module5_degree_matrix$count)

dat <- data.frame(degree = degree_num, count = degree_count)
head(dat)


#查看度分布
#可观察到微生物相关网络通常服从幂律分布
par(mfrow = c(1, 3))
hist(dat, xlab = 'Degree', ylab = 'Frequency', 
     main = 'Degree distribution')
plot(degree_count, degree_num, xlab = 'Degree', ylab = 'Count', 
     main = 'Degree distribution')
plot(degree_count, degree_num, log = 'xy', xlab = 'Log-degree', 
     ylab = 'Log-count', main = 'Log-log degree distribution')

##在R中，可通过非线性回归函数nls()拟合幂律分布。

#拟合，a 和 b 的初始值手动指定，数据不同需要多加尝试
mod <- nls(count ~ a*degree^b, data = dat, start = list(a = 2, b = 1.5))
summary(mod)

#提取关键值
a <- round(coef(mod)[1], 3)
b <- round(coef(mod)[2], 3)
a; b


#使用构建好的方程，通过预测变量拟合响应变量，再根据和观测值的差异，获得 R2
#SSre：the sum of the squares of the distances of the points from the fit
#SStot：the sum of the squares of the distances of the points from a horizontal line through the mean of all Y values
fit <- fitted(mod)
SSre <- sum((dat$count-fit)^2)
SStot <- sum((dat$count-mean(dat$count))^2)
R2 <- round(1 - SSre/SStot, 3)
R2


#p 值可以根据置换检验的原理获得
#将 count 的值分别随机置换 N 次（例如 999 次），通过随机置换数据后数据获取 R2（R2'）
#比较随机置换后值的 R2' 大于观测值的 R2 的频率，即为 p 值
p_num <- 1
dat_rand <- dat
for (i in 1:999) {
  dat_rand$count <- sample(dat_rand$count)
  SSre_rand <- sum((dat_rand$count-fit)^2)
  SStot_rand <- sum((dat_rand$count-mean(dat_rand$count))^2)
  R2_rand <- 1 - SSre_rand/SStot_rand
  if (R2_rand > R2) p_num <- p_num + 1
}
p_value <- p_num / (999+1)
p_value


#作图展示
library(ggplot2)

module5_p <- ggplot(dat, aes(x = degree, y = count)) +
  geom_point(color = 'blue') +
  theme(panel.grid.major = element_line(color = 'gray'), panel.background = element_rect(color = 'black', fill = 'transparent')) +
  stat_smooth(method = 'nls', formula = y ~ a*x^b, method.args = list(start = list(a = 2, b = 1.5)), se = FALSE) +
  labs(x = 'Degree', y = 'Count')
module5_p
#添加公式拟合的注释
label <- data.frame(
  formula = sprintf('italic(Y) == %.3f*italic(X)^%.3f', a, b),
  R2 = sprintf('italic(R^2) == %.3f', R2),
  p_value = sprintf('italic(P) < %.3f', p_value)
)

module5_p2 <- module5_p + geom_text(x = 35, y = 8, aes(label = formula), data = label, parse = TRUE, hjust = 0) +
  geom_text(x = 35, y = 6, aes(label = R2), data = label, parse = TRUE, hjust = 0) +
  geom_text(x = 35, y = 4, aes(label = p_value), data = label, parse = TRUE, hjust = 0)

module5_p2

ggsave("module5_degree.pdf", module5_p2, width = 60, height = 50, units = "mm")
ggsave("module5_degree.png", module5_p2, width = 60, height = 50, units = "mm")
#ggsave("./H_zi_pi_p.tiff", H_zi_pi_p, width = 150, height = 100, units = "mm")




####module6####  
module6_node <- read.table("All_module6_gephi2.txt", row.names=1,sep="\t", header=T, blank.lines.skip=F, check.names=F)
#读取数据的时候，header=T表示保留原始表格的表头.row.names=1表示保留数据表格的第一列，=2表示去除表格的第一列
module6_node

#已知网络节点列表时，对网络度分布进行统计，当未知时参见生信小白鱼的代码:200302-微生物共发生网络中节点度的幂律分布和在R中拟合
# Get unique labels
unique_labels <- unique(module6_node$degree)
unique_labels

# Initialize empty matrix
module6_degree_matrix <- matrix(0, nrow=length(unique_labels), ncol=2, dimnames=list(as.character(unique_labels), c("degree","count")))

# Loop through each label
for (i in unique_labels) {
  
  # Subset edges with the current label
  module6_degree_count <- module6_node[module6_node$degree == i,] 
  dim(module6_degree_count)
  
  # Create matrix with positive and negative degrees
  module6_degree_i <- matrix(c(i, dim(module6_degree_count)[1]), nrow=1, ncol=2, byrow=TRUE, dimnames=list(as.character(i),c("degree","count")))
  module6_degree_i
  
  # Fill in matrix with values from current label
  module6_degree_matrix[as.character(i), ] <- module6_degree_i[, c(1,2)]
  #特别要注意，当unique函数的结果是字符时，i in unique()结果是字符，作为矩阵是不需要转换，但当时数字时，需要用as.character(i)来将数字i转化为字符i
}

print(module6_degree_matrix)
module6_degree_matrix <- as.data.frame(module6_degree_matrix)
rownames(module6_degree_matrix) <- NULL 
module6_degree_matrix

#结果中发现存在0值，需要剔除
#module0_degree_matrix <- module0_degree_matrix[-42,]
#module0_degree_matrix
write.table(module6_degree_matrix,paste0(output,"module6_degree_matrix.txt"),sep="\t",quote=F)


#度分布统计
#degree_dist <- table(V(igraph)$degree)
degree_num <- as.numeric(module6_degree_matrix$degree)
degree_count <- as.numeric(module6_degree_matrix$count)

dat <- data.frame(degree = degree_num, count = degree_count)
head(dat)


#查看度分布
#可观察到微生物相关网络通常服从幂律分布
par(mfrow = c(1, 3))
hist(dat, xlab = 'Degree', ylab = 'Frequency', 
     main = 'Degree distribution')
plot(degree_count, degree_num, xlab = 'Degree', ylab = 'Count', 
     main = 'Degree distribution')
plot(degree_count, degree_num, log = 'xy', xlab = 'Log-degree', 
     ylab = 'Log-count', main = 'Log-log degree distribution')

##在R中，可通过非线性回归函数nls()拟合幂律分布。

#拟合，a 和 b 的初始值手动指定，数据不同需要多加尝试
mod <- nls(count ~ a*degree^b, data = dat, start = list(a = 2, b = 1.5))
summary(mod)

#提取关键值
a <- round(coef(mod)[1], 3)
b <- round(coef(mod)[2], 3)
a; b


#使用构建好的方程，通过预测变量拟合响应变量，再根据和观测值的差异，获得 R2
#SSre：the sum of the squares of the distances of the points from the fit
#SStot：the sum of the squares of the distances of the points from a horizontal line through the mean of all Y values
fit <- fitted(mod)
SSre <- sum((dat$count-fit)^2)
SStot <- sum((dat$count-mean(dat$count))^2)
R2 <- round(1 - SSre/SStot, 3)
R2


#p 值可以根据置换检验的原理获得
#将 count 的值分别随机置换 N 次（例如 999 次），通过随机置换数据后数据获取 R2（R2'）
#比较随机置换后值的 R2' 大于观测值的 R2 的频率，即为 p 值
p_num <- 1
dat_rand <- dat
for (i in 1:999) {
  dat_rand$count <- sample(dat_rand$count)
  SSre_rand <- sum((dat_rand$count-fit)^2)
  SStot_rand <- sum((dat_rand$count-mean(dat_rand$count))^2)
  R2_rand <- 1 - SSre_rand/SStot_rand
  if (R2_rand > R2) p_num <- p_num + 1
}
p_value <- p_num / (999+1)
p_value


#作图展示
library(ggplot2)

module6_p <- ggplot(dat, aes(x = degree, y = count)) +
  geom_point(color = 'blue') +
  theme(panel.grid.major = element_line(color = 'gray'), panel.background = element_rect(color = 'black', fill = 'transparent')) +
  stat_smooth(method = 'nls', formula = y ~ a*x^b, method.args = list(start = list(a = 2, b = 1.5)), se = FALSE) +
  labs(x = 'Degree', y = 'Count')
module6_p
#添加公式拟合的注释
label <- data.frame(
  formula = sprintf('italic(Y) == %.3f*italic(X)^%.3f', a, b),
  R2 = sprintf('italic(R^2) == %.3f', R2),
  p_value = sprintf('italic(P) < %.3f', p_value)
)

module6_p2 <- module6_p + geom_text(x = 35, y = 4.5, aes(label = formula), data = label, parse = TRUE, hjust = 0) +
  geom_text(x = 35, y = 4, aes(label = R2), data = label, parse = TRUE, hjust = 0) +
  geom_text(x = 35, y = 3.5, aes(label = p_value), data = label, parse = TRUE, hjust = 0)

module6_p2

ggsave("module6_degree.pdf", module6_p2, width = 60, height = 50, units = "mm")
ggsave("module6_degree.png", module6_p2, width = 60, height = 50, units = "mm")
#ggsave("./H_zi_pi_p.tiff", H_zi_pi_p, width = 150, height = 100, units = "mm")


####module7####  
module7_node <- read.table("All_module7_gephi2.txt", row.names=1,sep="\t", header=T, blank.lines.skip=F, check.names=F)
#读取数据的时候，header=T表示保留原始表格的表头.row.names=1表示保留数据表格的第一列，=2表示去除表格的第一列
module7_node

#已知网络节点列表时，对网络度分布进行统计，当未知时参见生信小白鱼的代码:200302-微生物共发生网络中节点度的幂律分布和在R中拟合
# Get unique labels
unique_labels <- unique(module7_node$degree)
unique_labels

# Initialize empty matrix
module7_degree_matrix <- matrix(0, nrow=length(unique_labels), ncol=2, dimnames=list(as.character(unique_labels), c("degree","count")))

# Loop through each label
for (i in unique_labels) {
  
  # Subset edges with the current label
  module7_degree_count <- module7_node[module7_node$degree == i,] 
  dim(module7_degree_count)
  
  # Create matrix with positive and negative degrees
  module7_degree_i <- matrix(c(i, dim(module7_degree_count)[1]), nrow=1, ncol=2, byrow=TRUE, dimnames=list(as.character(i),c("degree","count")))
  module7_degree_i
  
  # Fill in matrix with values from current label
  module7_degree_matrix[as.character(i), ] <- module7_degree_i[, c(1,2)]
  #特别要注意，当unique函数的结果是字符时，i in unique()结果是字符，作为矩阵是不需要转换，但当时数字时，需要用as.character(i)来将数字i转化为字符i
}

print(module7_degree_matrix)
module7_degree_matrix <- as.data.frame(module7_degree_matrix)
rownames(module7_degree_matrix) <- NULL 
module7_degree_matrix

#结果中发现存在0值，需要剔除
#module0_degree_matrix <- module0_degree_matrix[-42,]
#module0_degree_matrix
write.table(module7_degree_matrix,paste0(output,"module7_degree_matrix.txt"),sep="\t",quote=F)


#度分布统计
#degree_dist <- table(V(igraph)$degree)
degree_num <- as.numeric(module7_degree_matrix$degree)
degree_count <- as.numeric(module7_degree_matrix$count)

dat <- data.frame(degree = degree_num, count = degree_count)
head(dat)


#查看度分布
#可观察到微生物相关网络通常服从幂律分布
par(mfrow = c(1, 3))
hist(dat, xlab = 'Degree', ylab = 'Frequency', 
     main = 'Degree distribution')
plot(degree_count, degree_num, xlab = 'Degree', ylab = 'Count', 
     main = 'Degree distribution')
plot(degree_count, degree_num, log = 'xy', xlab = 'Log-degree', 
     ylab = 'Log-count', main = 'Log-log degree distribution')

##在R中，可通过非线性回归函数nls()拟合幂律分布。

#拟合，a 和 b 的初始值手动指定，数据不同需要多加尝试
mod <- nls(count ~ a*degree^b, data = dat, start = list(a = 2, b = 1.5))
summary(mod)

#提取关键值
a <- round(coef(mod)[1], 3)
b <- round(coef(mod)[2], 3)
a; b


#使用构建好的方程，通过预测变量拟合响应变量，再根据和观测值的差异，获得 R2
#SSre：the sum of the squares of the distances of the points from the fit
#SStot：the sum of the squares of the distances of the points from a horizontal line through the mean of all Y values
fit <- fitted(mod)
SSre <- sum((dat$count-fit)^2)
SStot <- sum((dat$count-mean(dat$count))^2)
R2 <- round(1 - SSre/SStot, 3)
R2


#p 值可以根据置换检验的原理获得
#将 count 的值分别随机置换 N 次（例如 999 次），通过随机置换数据后数据获取 R2（R2'）
#比较随机置换后值的 R2' 大于观测值的 R2 的频率，即为 p 值
p_num <- 1
dat_rand <- dat
for (i in 1:999) {
  dat_rand$count <- sample(dat_rand$count)
  SSre_rand <- sum((dat_rand$count-fit)^2)
  SStot_rand <- sum((dat_rand$count-mean(dat_rand$count))^2)
  R2_rand <- 1 - SSre_rand/SStot_rand
  if (R2_rand > R2) p_num <- p_num + 1
}
p_value <- p_num / (999+1)
p_value


#作图展示
library(ggplot2)

module7_p <- ggplot(dat, aes(x = degree, y = count)) +
  geom_point(color = 'blue') +
  theme(panel.grid.major = element_line(color = 'gray'), panel.background = element_rect(color = 'black', fill = 'transparent')) +
  stat_smooth(method = 'nls', formula = y ~ a*x^b, method.args = list(start = list(a = 2, b = 1.5)), se = FALSE) +
  labs(x = 'Degree', y = 'Count')
module7_p
#添加公式拟合的注释
label <- data.frame(
  formula = sprintf('italic(Y) == %.3f*italic(X)^%.3f', a, b),
  R2 = sprintf('italic(R^2) == %.3f', R2),
  p_value = sprintf('italic(P) < %.3f', p_value)
)

module7_p2 <- module7_p + geom_text(x = 20, y = 5.5, aes(label = formula), data = label, parse = TRUE, hjust = 0) +
  geom_text(x = 20, y = 5, aes(label = R2), data = label, parse = TRUE, hjust = 0) +
  geom_text(x = 20, y = 4, aes(label = p_value), data = label, parse = TRUE, hjust = 0)

module7_p2

ggsave("module7_degree.pdf", module7_p2, width = 60, height = 50, units = "mm")
ggsave("module7_degree.png", module7_p2, width = 60, height = 50, units = "mm")
#ggsave("./H_zi_pi_p.tiff", H_zi_pi_p, width = 150, height = 100, units = "mm")



####module8####  
module8_node <- read.table("All_module8_gephi2.txt", row.names=1,sep="\t", header=T, blank.lines.skip=F, check.names=F)
#读取数据的时候，header=T表示保留原始表格的表头.row.names=1表示保留数据表格的第一列，=2表示去除表格的第一列
module8_node

#已知网络节点列表时，对网络度分布进行统计，当未知时参见生信小白鱼的代码:200302-微生物共发生网络中节点度的幂律分布和在R中拟合
# Get unique labels
unique_labels <- unique(module8_node$degree)
unique_labels

# Initialize empty matrix
module8_degree_matrix <- matrix(0, nrow=length(unique_labels), ncol=2, dimnames=list(as.character(unique_labels), c("degree","count")))

# Loop through each label
for (i in unique_labels) {
  
  # Subset edges with the current label
  module8_degree_count <- module8_node[module8_node$degree == i,] 
  dim(module8_degree_count)
  
  # Create matrix with positive and negative degrees
  module8_degree_i <- matrix(c(i, dim(module8_degree_count)[1]), nrow=1, ncol=2, byrow=TRUE, dimnames=list(as.character(i),c("degree","count")))
  module8_degree_i
  
  # Fill in matrix with values from current label
  module8_degree_matrix[as.character(i), ] <- module8_degree_i[, c(1,2)]
  #特别要注意，当unique函数的结果是字符时，i in unique()结果是字符，作为矩阵是不需要转换，但当时数字时，需要用as.character(i)来将数字i转化为字符i
}

print(module8_degree_matrix)
module8_degree_matrix <- as.data.frame(module8_degree_matrix)
rownames(module8_degree_matrix) <- NULL 
module8_degree_matrix

#结果中发现存在0值，需要剔除
#module0_degree_matrix <- module0_degree_matrix[-42,]
#module0_degree_matrix
write.table(module8_degree_matrix,paste0(output,"module8_degree_matrix.txt"),sep="\t",quote=F)


#度分布统计
#degree_dist <- table(V(igraph)$degree)
degree_num <- as.numeric(module8_degree_matrix$degree)
degree_count <- as.numeric(module8_degree_matrix$count)

dat <- data.frame(degree = degree_num, count = degree_count)
head(dat)


#查看度分布
#可观察到微生物相关网络通常服从幂律分布
par(mfrow = c(1, 3))
hist(dat, xlab = 'Degree', ylab = 'Frequency', 
     main = 'Degree distribution')
plot(degree_count, degree_num, xlab = 'Degree', ylab = 'Count', 
     main = 'Degree distribution')
plot(degree_count, degree_num, log = 'xy', xlab = 'Log-degree', 
     ylab = 'Log-count', main = 'Log-log degree distribution')

##在R中，可通过非线性回归函数nls()拟合幂律分布。

#拟合，a 和 b 的初始值手动指定，数据不同需要多加尝试
mod <- nls(count ~ a*degree^b, data = dat, start = list(a = 2, b = 1.5))
summary(mod)

#提取关键值
a <- round(coef(mod)[1], 3)
b <- round(coef(mod)[2], 3)
a; b


#使用构建好的方程，通过预测变量拟合响应变量，再根据和观测值的差异，获得 R2
#SSre：the sum of the squares of the distances of the points from the fit
#SStot：the sum of the squares of the distances of the points from a horizontal line through the mean of all Y values
fit <- fitted(mod)
SSre <- sum((dat$count-fit)^2)
SStot <- sum((dat$count-mean(dat$count))^2)
R2 <- round(1 - SSre/SStot, 3)
R2


#p 值可以根据置换检验的原理获得
#将 count 的值分别随机置换 N 次（例如 999 次），通过随机置换数据后数据获取 R2（R2'）
#比较随机置换后值的 R2' 大于观测值的 R2 的频率，即为 p 值
p_num <- 1
dat_rand <- dat
for (i in 1:999) {
  dat_rand$count <- sample(dat_rand$count)
  SSre_rand <- sum((dat_rand$count-fit)^2)
  SStot_rand <- sum((dat_rand$count-mean(dat_rand$count))^2)
  R2_rand <- 1 - SSre_rand/SStot_rand
  if (R2_rand > R2) p_num <- p_num + 1
}
p_value <- p_num / (999+1)
p_value


#作图展示
library(ggplot2)

module8_p <- ggplot(dat, aes(x = degree, y = count)) +
  geom_point(color = 'blue') +
  theme(panel.grid.major = element_line(color = 'gray'), panel.background = element_rect(color = 'black', fill = 'transparent')) +
  stat_smooth(method = 'nls', formula = y ~ a*x^b, method.args = list(start = list(a = 2, b = 1.5)), se = FALSE) +
  labs(x = 'Degree', y = 'Count')
module8_p
#添加公式拟合的注释
label <- data.frame(
  formula = sprintf('italic(Y) == %.3f*italic(X)^%.3f', a, b),
  R2 = sprintf('italic(R^2) == %.3f', R2),
  p_value = sprintf('italic(P) < %.3f', p_value)
)

module8_p2 <- module8_p + geom_text(x = 15, y = 6, aes(label = formula), data = label, parse = TRUE, hjust = 0) +
  geom_text(x = 15, y = 5, aes(label = R2), data = label, parse = TRUE, hjust = 0) +
  geom_text(x = 15, y = 4, aes(label = p_value), data = label, parse = TRUE, hjust = 0)

module8_p2

ggsave("module8_degree.pdf", module8_p2, width = 60, height = 50, units = "mm")
ggsave("module8_degree.png", module8_p2, width = 60, height = 50, units = "mm")




####module9####  
module9_node <- read.table("All_module9_gephi2.txt", row.names=1,sep="\t", header=T, blank.lines.skip=F, check.names=F)
#读取数据的时候，header=T表示保留原始表格的表头.row.names=1表示保留数据表格的第一列，=2表示去除表格的第一列
module9_node

#已知网络节点列表时，对网络度分布进行统计，当未知时参见生信小白鱼的代码:200302-微生物共发生网络中节点度的幂律分布和在R中拟合
# Get unique labels
unique_labels <- unique(module9_node$degree)
unique_labels

# Initialize empty matrix
module9_degree_matrix <- matrix(0, nrow=length(unique_labels), ncol=2, dimnames=list(as.character(unique_labels), c("degree","count")))

# Loop through each label
for (i in unique_labels) {
  
  # Subset edges with the current label
  module9_degree_count <- module9_node[module9_node$degree == i,] 
  dim(module9_degree_count)
  
  # Create matrix with positive and negative degrees
  module9_degree_i <- matrix(c(i, dim(module9_degree_count)[1]), nrow=1, ncol=2, byrow=TRUE, dimnames=list(as.character(i),c("degree","count")))
  module9_degree_i
  
  # Fill in matrix with values from current label
  module9_degree_matrix[as.character(i), ] <- module9_degree_i[, c(1,2)]
  #特别要注意，当unique函数的结果是字符时，i in unique()结果是字符，作为矩阵是不需要转换，但当时数字时，需要用as.character(i)来将数字i转化为字符i
}

print(module9_degree_matrix)
module9_degree_matrix <- as.data.frame(module9_degree_matrix)
rownames(module9_degree_matrix) <- NULL 
module9_degree_matrix

#结果中发现存在0值，需要剔除
#module0_degree_matrix <- module0_degree_matrix[-42,]
#module0_degree_matrix
write.table(module9_degree_matrix,paste0(output,"module9_degree_matrix.txt"),sep="\t",quote=F)


#度分布统计
#degree_dist <- table(V(igraph)$degree)
degree_num <- as.numeric(module9_degree_matrix$degree)
degree_count <- as.numeric(module9_degree_matrix$count)

dat <- data.frame(degree = degree_num, count = degree_count)
head(dat)


#查看度分布
#可观察到微生物相关网络通常服从幂律分布
par(mfrow = c(1, 3))
hist(dat, xlab = 'Degree', ylab = 'Frequency', 
     main = 'Degree distribution')
plot(degree_count, degree_num, xlab = 'Degree', ylab = 'Count', 
     main = 'Degree distribution')
plot(degree_count, degree_num, log = 'xy', xlab = 'Log-degree', 
     ylab = 'Log-count', main = 'Log-log degree distribution')

##在R中，可通过非线性回归函数nls()拟合幂律分布。

#拟合，a 和 b 的初始值手动指定，数据不同需要多加尝试
mod <- nls(count ~ a*degree^b, data = dat, start = list(a = 2, b = 1.5))
summary(mod)

#提取关键值
a <- round(coef(mod)[1], 3)
b <- round(coef(mod)[2], 3)
a; b


#使用构建好的方程，通过预测变量拟合响应变量，再根据和观测值的差异，获得 R2
#SSre：the sum of the squares of the distances of the points from the fit
#SStot：the sum of the squares of the distances of the points from a horizontal line through the mean of all Y values
fit <- fitted(mod)
SSre <- sum((dat$count-fit)^2)
SStot <- sum((dat$count-mean(dat$count))^2)
R2 <- round(1 - SSre/SStot, 3)
R2


#p 值可以根据置换检验的原理获得
#将 count 的值分别随机置换 N 次（例如 999 次），通过随机置换数据后数据获取 R2（R2'）
#比较随机置换后值的 R2' 大于观测值的 R2 的频率，即为 p 值
p_num <- 1
dat_rand <- dat
for (i in 1:999) {
  dat_rand$count <- sample(dat_rand$count)
  SSre_rand <- sum((dat_rand$count-fit)^2)
  SStot_rand <- sum((dat_rand$count-mean(dat_rand$count))^2)
  R2_rand <- 1 - SSre_rand/SStot_rand
  if (R2_rand > R2) p_num <- p_num + 1
}
p_value <- p_num / (999+1)
p_value


#作图展示
library(ggplot2)

module9_p <- ggplot(dat, aes(x = degree, y = count)) +
  geom_point(color = 'blue') +
  theme(panel.grid.major = element_line(color = 'gray'), panel.background = element_rect(color = 'black', fill = 'transparent')) +
  stat_smooth(method = 'nls', formula = y ~ a*x^b, method.args = list(start = list(a = 2, b = 1.5)), se = FALSE) +
  labs(x = 'Degree', y = 'Count')
module9_p
#添加公式拟合的注释
label <- data.frame(
  formula = sprintf('italic(Y) == %.3f*italic(X)^%.3f', a, b),
  R2 = sprintf('italic(R^2) == %.3f', R2),
  p_value = sprintf('italic(P) < %.3f', p_value)
)

module9_p2 <- module9_p + geom_text(x = 30, y = 13, aes(label = formula), data = label, parse = TRUE, hjust = 0) +
  geom_text(x = 30, y = 11, aes(label = R2), data = label, parse = TRUE, hjust = 0) +
  geom_text(x = 30, y = 9, aes(label = p_value), data = label, parse = TRUE, hjust = 0)

module9_p2

ggsave("module9_degree.pdf", module9_p2, width = 60, height = 50, units = "mm")
ggsave("module9_degree.png", module9_p2, width = 60, height = 50, units = "mm")




####module10####  
module10_node <- read.table("All_module10_gephi2.txt", row.names=1,sep="\t", header=T, blank.lines.skip=F, check.names=F)
#读取数据的时候，header=T表示保留原始表格的表头.row.names=1表示保留数据表格的第一列，=2表示去除表格的第一列
module10_node

#已知网络节点列表时，对网络度分布进行统计，当未知时参见生信小白鱼的代码:200302-微生物共发生网络中节点度的幂律分布和在R中拟合
# Get unique labels
unique_labels <- unique(module10_node$degree)
unique_labels

# Initialize empty matrix
module10_degree_matrix <- matrix(0, nrow=length(unique_labels), ncol=2, dimnames=list(as.character(unique_labels), c("degree","count")))

# Loop through each label
for (i in unique_labels) {
  
  # Subset edges with the current label
  module10_degree_count <- module10_node[module10_node$degree == i,] 
  dim(module10_degree_count)
  
  # Create matrix with positive and negative degrees
  module10_degree_i <- matrix(c(i, dim(module10_degree_count)[1]), nrow=1, ncol=2, byrow=TRUE, dimnames=list(as.character(i),c("degree","count")))
  module10_degree_i
  
  # Fill in matrix with values from current label
  module10_degree_matrix[as.character(i), ] <- module10_degree_i[, c(1,2)]
  #特别要注意，当unique函数的结果是字符时，i in unique()结果是字符，作为矩阵是不需要转换，但当时数字时，需要用as.character(i)来将数字i转化为字符i
}

print(module10_degree_matrix)
module10_degree_matrix <- as.data.frame(module10_degree_matrix)
rownames(module10_degree_matrix) <- NULL 
module10_degree_matrix

#结果中发现存在0值，需要剔除
#module0_degree_matrix <- module0_degree_matrix[-42,]
#module0_degree_matrix
write.table(module10_degree_matrix,paste0(output,"module10_degree_matrix.txt"),sep="\t",quote=F)


#度分布统计
#degree_dist <- table(V(igraph)$degree)
degree_num <- as.numeric(module10_degree_matrix$degree)
degree_count <- as.numeric(module10_degree_matrix$count)

dat <- data.frame(degree = degree_num, count = degree_count)
head(dat)


#查看度分布
#可观察到微生物相关网络通常服从幂律分布
par(mfrow = c(1, 3))
hist(dat, xlab = 'Degree', ylab = 'Frequency', 
     main = 'Degree distribution')
plot(degree_count, degree_num, xlab = 'Degree', ylab = 'Count', 
     main = 'Degree distribution')
plot(degree_count, degree_num, log = 'xy', xlab = 'Log-degree', 
     ylab = 'Log-count', main = 'Log-log degree distribution')

##在R中，可通过非线性回归函数nls()拟合幂律分布。

#拟合，a 和 b 的初始值手动指定，数据不同需要多加尝试
mod <- nls(count ~ a*degree^b, data = dat, start = list(a = 2, b = 1.5))
summary(mod)

#提取关键值
a <- round(coef(mod)[1], 3)
b <- round(coef(mod)[2], 3)
a; b


#使用构建好的方程，通过预测变量拟合响应变量，再根据和观测值的差异，获得 R2
#SSre：the sum of the squares of the distances of the points from the fit
#SStot：the sum of the squares of the distances of the points from a horizontal line through the mean of all Y values
fit <- fitted(mod)
SSre <- sum((dat$count-fit)^2)
SStot <- sum((dat$count-mean(dat$count))^2)
R2 <- round(1 - SSre/SStot, 3)
R2


#p 值可以根据置换检验的原理获得
#将 count 的值分别随机置换 N 次（例如 999 次），通过随机置换数据后数据获取 R2（R2'）
#比较随机置换后值的 R2' 大于观测值的 R2 的频率，即为 p 值
p_num <- 1
dat_rand <- dat
for (i in 1:999) {
  dat_rand$count <- sample(dat_rand$count)
  SSre_rand <- sum((dat_rand$count-fit)^2)
  SStot_rand <- sum((dat_rand$count-mean(dat_rand$count))^2)
  R2_rand <- 1 - SSre_rand/SStot_rand
  if (R2_rand > R2) p_num <- p_num + 1
}
p_value <- p_num / (999+1)
p_value


#作图展示
library(ggplot2)

module10_p <- ggplot(dat, aes(x = degree, y = count)) +
  geom_point(color = 'blue') +
  theme(panel.grid.major = element_line(color = 'gray'), panel.background = element_rect(color = 'black', fill = 'transparent')) +
  stat_smooth(method = 'nls', formula = y ~ a*x^b, method.args = list(start = list(a = 2, b = 1.5)), se = FALSE) +
  labs(x = 'Degree', y = 'Count')
module10_p
#添加公式拟合的注释
label <- data.frame(
  formula = sprintf('italic(Y) == %.3f*italic(X)^%.3f', a, b),
  R2 = sprintf('italic(R^2) == %.3f', R2),
  p_value = sprintf('italic(P) < %.3f', p_value)
)

module10_p2 <- module10_p + geom_text(x = 25, y = 22, aes(label = formula), data = label, parse = TRUE, hjust = 0) +
  geom_text(x = 25, y = 18, aes(label = R2), data = label, parse = TRUE, hjust = 0) +
  geom_text(x = 25, y = 15, aes(label = p_value), data = label, parse = TRUE, hjust = 0)

module10_p2

ggsave("module10_degree.pdf", module10_p2, width = 60, height = 50, units = "mm")
ggsave("module10_degree.png", module10_p2, width = 60, height = 50, units = "mm")




####module10####  
module10_node <- read.table("All_module10_gephi2.txt", row.names=1,sep="\t", header=T, blank.lines.skip=F, check.names=F)
#读取数据的时候，header=T表示保留原始表格的表头.row.names=1表示保留数据表格的第一列，=2表示去除表格的第一列
module10_node

#已知网络节点列表时，对网络度分布进行统计，当未知时参见生信小白鱼的代码:200302-微生物共发生网络中节点度的幂律分布和在R中拟合
# Get unique labels
unique_labels <- unique(module10_node$degree)
unique_labels

# Initialize empty matrix
module10_degree_matrix <- matrix(0, nrow=length(unique_labels), ncol=2, dimnames=list(as.character(unique_labels), c("degree","count")))

# Loop through each label
for (i in unique_labels) {
  
  # Subset edges with the current label
  module10_degree_count <- module10_node[module10_node$degree == i,] 
  dim(module10_degree_count)
  
  # Create matrix with positive and negative degrees
  module10_degree_i <- matrix(c(i, dim(module10_degree_count)[1]), nrow=1, ncol=2, byrow=TRUE, dimnames=list(as.character(i),c("degree","count")))
  module10_degree_i
  
  # Fill in matrix with values from current label
  module10_degree_matrix[as.character(i), ] <- module10_degree_i[, c(1,2)]
  #特别要注意，当unique函数的结果是字符时，i in unique()结果是字符，作为矩阵是不需要转换，但当时数字时，需要用as.character(i)来将数字i转化为字符i
}

print(module10_degree_matrix)
module10_degree_matrix <- as.data.frame(module10_degree_matrix)
rownames(module10_degree_matrix) <- NULL 
module10_degree_matrix

#结果中发现存在0值，需要剔除
#module0_degree_matrix <- module0_degree_matrix[-42,]
#module0_degree_matrix
write.table(module10_degree_matrix,paste0(output,"module10_degree_matrix.txt"),sep="\t",quote=F)


#度分布统计
#degree_dist <- table(V(igraph)$degree)
degree_num <- as.numeric(module10_degree_matrix$degree)
degree_count <- as.numeric(module10_degree_matrix$count)

dat <- data.frame(degree = degree_num, count = degree_count)
head(dat)


#查看度分布
#可观察到微生物相关网络通常服从幂律分布
par(mfrow = c(1, 3))
hist(dat, xlab = 'Degree', ylab = 'Frequency', 
     main = 'Degree distribution')
plot(degree_count, degree_num, xlab = 'Degree', ylab = 'Count', 
     main = 'Degree distribution')
plot(degree_count, degree_num, log = 'xy', xlab = 'Log-degree', 
     ylab = 'Log-count', main = 'Log-log degree distribution')

##在R中，可通过非线性回归函数nls()拟合幂律分布。

#拟合，a 和 b 的初始值手动指定，数据不同需要多加尝试
mod <- nls(count ~ a*degree^b, data = dat, start = list(a = 2, b = 1.5))
summary(mod)

#提取关键值
a <- round(coef(mod)[1], 3)
b <- round(coef(mod)[2], 3)
a; b


#使用构建好的方程，通过预测变量拟合响应变量，再根据和观测值的差异，获得 R2
#SSre：the sum of the squares of the distances of the points from the fit
#SStot：the sum of the squares of the distances of the points from a horizontal line through the mean of all Y values
fit <- fitted(mod)
SSre <- sum((dat$count-fit)^2)
SStot <- sum((dat$count-mean(dat$count))^2)
R2 <- round(1 - SSre/SStot, 3)
R2


#p 值可以根据置换检验的原理获得
#将 count 的值分别随机置换 N 次（例如 999 次），通过随机置换数据后数据获取 R2（R2'）
#比较随机置换后值的 R2' 大于观测值的 R2 的频率，即为 p 值
p_num <- 1
dat_rand <- dat
for (i in 1:999) {
  dat_rand$count <- sample(dat_rand$count)
  SSre_rand <- sum((dat_rand$count-fit)^2)
  SStot_rand <- sum((dat_rand$count-mean(dat_rand$count))^2)
  R2_rand <- 1 - SSre_rand/SStot_rand
  if (R2_rand > R2) p_num <- p_num + 1
}
p_value <- p_num / (999+1)
p_value


#作图展示
library(ggplot2)

module10_p <- ggplot(dat, aes(x = degree, y = count)) +
  geom_point(color = 'blue') +
  theme(panel.grid.major = element_line(color = 'gray'), panel.background = element_rect(color = 'black', fill = 'transparent')) +
  stat_smooth(method = 'nls', formula = y ~ a*x^b, method.args = list(start = list(a = 2, b = 1.5)), se = FALSE) +
  labs(x = 'Degree', y = 'Count')
module10_p
#添加公式拟合的注释
label <- data.frame(
  formula = sprintf('italic(Y) == %.3f*italic(X)^%.3f', a, b),
  R2 = sprintf('italic(R^2) == %.3f', R2),
  p_value = sprintf('italic(P) < %.3f', p_value)
)

module10_p2 <- module10_p + geom_text(x = 12, y = 10, aes(label = formula), data = label, parse = TRUE, hjust = 0) +
  geom_text(x = 12, y = 8, aes(label = R2), data = label, parse = TRUE, hjust = 0) +
  geom_text(x = 12, y = 6, aes(label = p_value), data = label, parse = TRUE, hjust = 0)

module10_p2

ggsave("module10_degree.pdf", module10_p2, width = 60, height = 50, units = "mm")
ggsave("module10_degree.png", module10_p2, width = 60, height = 50, units = "mm")






####分离菌株鉴定及功能####
####遗传进化树构建####



rm(list=ls())
setwd("F:/360MoveData/Administrator/Desktop/功能菌株/")
##导入本次需要的R包##

library(ggtreeExtra)
library(ggtree)
library(treeio)
library(tidytree)
library(ggstar)
library(ggplot2)
library(ggnewscale)

####文涛老师公众号####
tree <- read.tree("./isolates_otu3.txt")
tree
# The attributes of tip point
IR3 <- read.table("./IR3.txt")
IR3
# The attributes of first ring
Phylum <- read.table("./Phylum.txt",header = TRUE, sep ="\t")
colnames(Phylum)[1]<-c("ID")
rownames(Phylum) <- Phylum$ID
Phylum <- subset(Phylum,select = -ID)
Phylum 

# The attributes of second ring
#dt3 <- read.csv("./secondring.csv")


#-提取otu表和tax表格#--------
#otu = as.data.frame(t(vegan_otu(ps )))
#tax = as.data.frame(vegan_tax(ps ))
#合并otu表格和tax表格#---------
IR3_Phylum = merge(IR3,Phylum,by = "row.names",all = F)
head(IR3_Phylum)



#--枝条分组信息先不做了--------统一绘图tax树的仪一起做。#---构造分组信息文件#--------
colnames(IR3_Phylum)[1] = c("id")
groupInfo <- split(row.names(IR3_Phylum), IR3_Phylum$Phylum) # OTU and phylum for group
groupInfo

tree <- groupOTU(tree, groupInfo)
unique(IR3_Phylum$Phylum)
p<- ggtree(tree,layout="circular", branch.length = "none", ladderize = FALSE,aes(color=Phylum)) %<+% IR3_Phylum +
#  geom_point(aes(fill = Phylum),pch = 21,size = 3)+
  theme(legend.position = "bottom")#+
#  geom_tiplab(size=2.5,hjust = 0)
#+  geom_tiplab(offset=1, color="black")

p = p + xlim(-10,NA)
p


##此处保存是为了获得物种注释的图例部分，由于目前关于该代码的图例部分还不完善
ggsave("./tree_legend.pdf",p,width = 20,height = 18,limitsize = FALSE)
ggsave("./tree_legend.png",p,width = 20,height = 18,limitsize = FALSE)
ggsave("./tree_legend.tiff",p,width = 20,height = 18,limitsize = FALSE)


#如果您发现在ggtree中添加物种注释后，与进化树中节点不匹配，请检查以下几点：

#1.检查物种注释的名称是否与进化树中节点的名称完全匹配。可能会存在小写字母和大写字母之间的区别，或者节点名称中存在下划线或短横线的情况。

#2.检查进化树的分支长度是否正确。分支长度错误可能会导致物种注释与节点不匹配。

#3.检查是否存在多个相同名称的节点或物种注释。这可能会导致注释错误匹配到节点或物种。

#如果以上步骤都没有解决问题，请尝试在创建进化树和添加注释时指定相同的树形对象和注释对象。例如：


####基础进化树####
#-提取可视化进化树中的数据#----
df <- p$data
# 去掉枝节点
df <- df[df$isTip, ]
head(df)
df$label
df

p1 = p + geom_tiplab2(data = df, aes(x +1, y =y,label= label ,color = Phylum,angle = angle ),size=3, offset=0.1)
p1




##自定义进化树物种分类颜色##
p11 <- p1 + scale_color_manual(
  limits=c("Bacillus tequilensis", "Pseudomonas geniculata", "Bacillus subtilis",
           "uncultured bacterium", "Paenibacillus polymyxa", "Brachybacterium avium",
           "Enterobacter hormaechei ", "Bacillus cereus","Bacillus velezensis",
           "Achromobacter mucicolens", "Streptomyces pratensis"),
  values=c("#9B3A4D","#E2AE79","#394A92","#68AC57","#497EB2","#8E549E","#D2352C",
          "#8CBDA7","#566CA5","#F4C28F","#70A0AC")
)
#,"#BAAFD1"
p11








####节点高亮####
#查看节点，为什么查看节点，因为我们要按照节点进行高亮处理
#p12<-p1 +geom_label2(aes(subset=!isTip,label=node,size=2))
#p12

#p13 <- p1 + 
#  #geom_highlight(node=33,fill = "#B2182B",alpha=0.5) #+
#geom_highlight(node=34,fill = "#E69F00",alpha=0.5) #+
#geom_highlight(node=35,fill = "#56B4E9",alpha=0.5) #+
#  geom_highlight(node=36,fill = "#009E73",alpha=0.5) +
#  geom_highlight(node=37,fill = "#F0E442",alpha=0.5) +
#geom_highlight(node=38,fill = "#0072B2",alpha=0.5) #+
#geom_highlight(node=39,fill = "#D55E00",alpha=0.5) #+
#  geom_highlight(node=40,fill = "#CC79A7",alpha=0.5) +
#geom_highlight(node=41,fill = "#CC6666",alpha=0.5) #+
#geom_highlight(node=42,fill = "#9999CC",alpha=0.5) #+
#  geom_highlight(node=43,fill = "#66CC99",alpha=0.5) +
#  geom_highlight(node=44,fill = "#999999",alpha=0.5) 

#p13



####节点标签高亮####
#最基本树
#ggtree(tree)
#去掉枝长信息
#ggtree(tree,branch.length = "none")
#添加文字标签
#ggtree(tree,branch.length = "none")+geom_tiplab(size=1.5)
#显示表示Bootstrap值的点
#ggtree(tree,branch.length = "none")+
#  geom_tiplab(size=1.5)+
#  geom_nodepoint(size="support",x=x-0.5)

#变成环形
#ggtree(tree,
#       branch.length = "none",layout="circular")+
#  geom_tiplab(size=3)+
#  geom_nodepoint(aes(size="support",x=x-0.5),color="#8f8fc3")

#添加色块
#ggtree(tree,branch.length = "none",layout="circular")+


#p11 <-  p+
#  geom_strip(taxa1="OTU80", taxa2="TSBH60", offset = 6.0, 
#             barsize =50 ,extend=0.5, color="#a9ddd4", label = "",
#             offset.text = 3)+
#  geom_strip(taxa1="TSBSH56_2", taxa2="YEMF10_3F25", offset = 6.0, 
#             barsize =50 ,extend=0.5, color="#9ec3db", label = "",
#             offset.text = 3)+
#  geom_strip(taxa1="TYGD10_1D11C", taxa2="TYGH10_2H29", offset = 6.0, 
#             barsize =50 ,extend=0.5, color="#cbc7de", label = "",
#             offset.text = 3)+
#  geom_strip(taxa1="YEMH10_2H25", taxa2="YEMD10_3SD18", offset = 6.0, 
#             barsize =50 ,extend=0.5, color="#BDACA4", label = "",
#             offset.text = 3)+
#  geom_tiplab(size=5,align = TRUE,offset = 0.8, color= "black")


#p11

#"#FFC125","#87CEFA","#7B68EE","#808080"灰色,"#800080"紫色,
#"#9ACD32",青绿"#D15FEE"浅紫,"#FFC0CB"淡粉,"#EE6A50"橙黄,"#8DEEEE",淡绿
#"#006400"深绿,"#800000"橙红,"#B0171F"深红,"#191970"深蓝

#"#E0C88D","#BAB18E","#844E1D","#882707","#B875F","#AE8966","#C19E67",
#"#AEA063","#D8C3AA","#BDB6AO","#E4C963","#BDB468","#986232","#97431B",
#"#815933","C86C25","#A96B37","#B8A594","#92734B","#776A4F","#D74F0E",
#"#AE440F","#BDACA4","#D0CDC6","#875E36","997D62","#7D7260","#9E6936",
#"#9F9686","#9B8E69","#9F794F","#B49A88","#4B2F14","#2C2213","#E8C01C",
#"#CD8C1D","#69552B","#8F6A38","#948031","B07E58"


#-----------为了更好的注释叶节点我需要构建phyloseq文件
#biofilm <- read.csv("./biofilm.csv")
mean = rowMeans(IR3)
df = cbind(df,mean)
colnames(df)
#----叶节点控制
#----x坐标即为叶节点的坐标，我们通过调整x值达到控制整个图形的目的
#p1 = p + geom_point(data = df, aes(x = x, y,size =mean,fill = Phylum ),
#                   inherit.aes = FALSE,pch = 21,color = "black")
#p1

##这里我直接将OTU表格映射到了外环上，所以密密麻麻的

#-------下面的工作为优化注释#---------
#-------模仿graphlan使用形状填充节点
#p$layers[[3]] = NULL
#p1 = p + geom_point(data = df, aes(x = x, y,size =mean,color = Phylum,shape = Phylum ),
#                    inherit.aes = FALSE)
#p1

#--------添加热图:方法一，将热图添加到内环
#mean = as.data.frame(mean)
#IR3
#p2 = gheatmap(p28,IR3, offset = 12.5, width=1.5, font.size=1, hjust=-.1,low = "white",high = "#f36e42")
#p2


##此处保存是为了获得抑菌活性热图的图例部分，由于目前关于该代码的图例部分还不完善
#ggsave("./tree_heatmap_legend.pdf",p2,width = 20,height = 18,limitsize = FALSE)
#ggsave("./tree_heatmap_legend.png",p2,width = 20,height = 18,limitsize = FALSE)
#ggsave("./tree_heatmap_legend.tiff",p2,width = 20,height = 18,limitsize = FALSE)



# #-按照分组添加热图-似乎热图只能使用同一套颜色银映射系统--所以后面谈添加热图只能使用geom_tile

# install.packages("ggnewscale")
library(ggnewscale)
##注：首先要知道树文件中tree$tip.label的顺序，然后将后续叠加图的所有行名都按照树的tip.label进行排序
####IAA####
#在excel中手动排序的IAA表格
dt_IAA <- read.csv("./IAA.csv")

library(ggnewscale)
p21 <- p11 + new_scale_fill()

p22 <- p21 +
  geom_fruit(
    data=dt_IAA,
    geom=geom_tile,
    mapping=aes(y=ID, alpha=IAA, fill=FALSE),
    color = "white",
    offset = 5.0,
    width=2,
    size = 0.1
  )+
  scale_alpha_continuous(
    range=c(0, 1),
    guide=guide_legend(keywidth = 1, keyheight = 0.3, order=5)
  ) +
  scale_fill_manual(
    values=c("#ae8ec2"),#,"#FFA500","#FF0000","#800000","#FFA500","#FF0000","#800000"
    guide=guide_legend(keywidth = 1, keyheight = 0.3, order=4)
  ) #+ theme(legend.position="none")##注意：注释完成才能产生图例
p22

##此处保存是为了获得IAA活性热图的图例部分，由于目前关于该代码的图例部分还不完善
ggsave("./tree_IAA_legend.pdf",p22,width = 20,height = 18,limitsize = FALSE)
ggsave("./tree_IAA_legend.png",p22,width = 20,height = 18,limitsize = FALSE)
ggsave("./tree_IAA_legend.tiff",p22,width = 20,height = 18,limitsize = FALSE)

####ACC####
dt_ACC <- read.csv("./ACC.csv")

library(ggnewscale)
p23 <- p21 + new_scale_fill()

p24 <- p23 +
  geom_fruit(
    data=dt_ACC,
    geom=geom_tile,
    mapping=aes(y=ID, alpha=ACC, fill=FALSE),
    color = "white",
    offset = 0.5,#特别需要注意的是此处的offset是指两层圈图之间的距离，并非是圈距离中心的距离
    width=2,
    size = 0.1
  )+
  scale_alpha_continuous(
    range=c(0, 1),
    guide=guide_legend(keywidth = 1, keyheight = 0.3, order=5)
  ) +
  scale_fill_manual(
    values=c("#69c4a9"),#,"#FFA500","#FF0000","#800000","#FFA500","#FF0000","#800000"
    guide=guide_legend(keywidth = 1, keyheight = 0.3, order=4)
  ) #+ theme(legend.position="none")
p24


##此处保存是为了获得ACC活性热图的图例部分，由于目前关于该代码的图例部分还不完善
ggsave("./tree_ACC_legend.pdf",p24,width = 20,height = 18,limitsize = FALSE)
ggsave("./tree_ACC_legend.png",p24,width = 20,height = 18,limitsize = FALSE)
ggsave("./tree_ACC_legend.tiff",p24,width = 20,height = 18,limitsize = FALSE)



####Sid####
dt_Sid <- read.csv("./Sid.csv")

library(ggnewscale)
p25 <- p21 + new_scale_fill()

p26 <- p25 +
  geom_fruit(
    data=dt_Sid,
    geom=geom_tile,
    mapping=aes(y=ID, alpha=Sid, fill=FALSE),
    color = "white",
    offset = 0.5,#特别需要注意的是此处的offset是指两层圈图之间的距离，并非是圈距离中心的距离
    width=2,
    size = 0.1
  )+
  scale_alpha_continuous(
    range=c(0, 1),
    guide=guide_legend(keywidth = 1, keyheight = 0.3, order=5)
  ) +
  scale_fill_manual(
    values=c("#2e8abf"),#,"#FFA500","#FF0000","#800000","#FFA500","#FF0000","#800000"
    guide=guide_legend(keywidth = 1, keyheight = 0.3, order=4)
  ) #+ theme(legend.position="none")
p26

##此处保存是为了获得Sid活性热图的图例部分，由于目前关于该代码的图例部分还不完善
#在p21的基础上获得图例
ggsave("./tree_Sid_legend.pdf",p26,width = 20,height = 18,limitsize = FALSE)
ggsave("./tree_Sid_legend.png",p26,width = 20,height = 18,limitsize = FALSE)
ggsave("./tree_Sid_legend.tiff",p26,width = 20,height = 18,limitsize = FALSE)




####P####
dt_P <- read.csv("./P.csv")

library(ggnewscale)
p27 <- p21 + new_scale_fill()

p28 <- p27 +
  geom_fruit(
    data=dt_P,
    geom=geom_tile,
    mapping=aes(y=ID, alpha=P, fill=FALSE),
    color = "white",
    offset = 0.5,#特别需要注意的是此处的offset是指两层圈图之间的距离，并非是圈距离中心的距离
    width=2,
    size = 0.1
  )+
  scale_alpha_continuous(
    range=c(0, 10),
    guide=guide_legend(keywidth = 1, keyheight = 0.3, order=5)
  ) +
  scale_fill_manual(
    values=c("#91d5da"),#,"#FFA500","#FF0000","#800000","#FFA500","#FF0000","#800000"
    guide=guide_legend(keywidth = 1, keyheight = 0.3, order=4)
  ) #+ theme(legend.position="none")
p28

##此处保存是为了获得P活性热图的图例部分，由于目前关于该代码的图例部分还不完善
#在p21的基础上获得图例
ggsave("./tree_P_legend.pdf",p28,width = 20,height = 18,limitsize = FALSE)
ggsave("./tree_P_legend.png",p28,width = 20,height = 18,limitsize = FALSE)
ggsave("./tree_P_legend.tiff",p28,width = 20,height = 18,limitsize = FALSE)



#---------外圈标记形状
p3 = p28 + new_scale_fill()+ geom_point(data = df, aes(x = x+29, y ),inherit.aes = FALSE,pch = 21,size = 2,fill = "red",color = "black")
p3

# library(ggcor) # + geom_star(data = df, aes(x = x+12, y ),inherit.aes = FALSE,fill = "blue",color = "black")
#p4 = p3 
#p4
# p4 + geom_point(data = df, aes(x = x+12, y,shape = Phylum,color=Phylum  ),inherit.aes = FALSE)
# 最外圈，添加矩形块，用于辨别丰度
#df$width= runif(length(df$x), 0, 4)
#p5 = p4 + geom_tile(data = df, aes(x = x+20, y ,width = width),color = "black",fill = "blue")
#p5

# 最外圈，添加矩形块，用于辨别丰度
#df$width= runif(length(df$x), 1, 5)
#p5 = p4 + geom_tile(data = df, aes(x = x+15, y ,width = width),color = "black",fill = "blue")
#p5

# 添加一条间隔线，可以设置宽度
#p6 = p5 + geom_bar(data = df, aes(x+8, y = y,color = "black"),inherit.aes = FALSE,stat = "identity",width = 0.3)
#p6

#geom_segment不能太小了，负责看上去怪异
#p7 = p6 + geom_segment(data = df, aes(x+18, y =y,xend = x+18 +width, yend = y),color = "red",size = 5)
#p7


# 添加标签
#p8 = p28 + geom_tiplab2(data = df, aes(x +13, y =y,label= label ,color = Phylum,angle = angle ),size=8, offset=0.1)
#p8
#-
ggsave("./tree.pdf",p28,width = 20,height = 18,limitsize = FALSE)
ggsave("./tree.png",p28,width = 20,height = 18,limitsize = FALSE)

##热图添加方式二：将热图添加至最外环
####IR####
library(reshape2)
IR3$ID <- rownames(IR3)
IR4 <- melt(IR3,id=c('ID'),
            measure.vars = c('Rs','Fg','Ss','Bb','Pc','Fon','Foc','Vm','Ggt'),
            variable.name = c('IR'),
            value.name = 'value')
IR4

p4 <- p28 +
  new_scale_fill() +
  geom_fruit(data = IR4,
             geom = geom_tile,
             mapping = aes(y=ID,x=IR,alpha= value,fill = IR),
             offset = 0.3,
             pwidth = 2.5,
             axis.params = list(
               axis="x",
               text.angle =-45,
               hjust=0
             )
           ) +
  scale_fill_manual(
    values = c("#f9d580","#ace4e4","#5ac8c8","#dfb0d6","#a5add3","#5698b9","#be64ac","#8c62aa","#3b4994"),
    guide=guide_legend(keywidth = 0.5,keyheight = 0.5,order=4)
  ) +
  scale_alpha_continuous(
    range=c(0,8),
    guide=guide_legend(keywidth = 0.5,keyheight = 0.5,order=5)
    
  )


p4

ggsave("./tree2.pdf",p4,width = 20,height = 18,limitsize = FALSE)
ggsave("./tree2.png",p4,width = 20,height = 18,limitsize = FALSE)
ggsave("./tree2.tiff",p4,width = 20,height = 18,limitsize = FALSE)



####Y叔公众号####

#为了重复
set.seed(1024)

#读取树文件
tr <- read.tree("./isolates_otu.nwk")
tr

#tip节点数
numtip <-length(tr$tip.label)
numtip



dat1 <-data.frame(ID=tr$tip.label,
                  Location=c(rep("HK",50),rep("TW",36),rep("SX",30),rep("GD",48),
                  rep("HN",20),rep("AH",20),rep("FJ",26)),








####复合菌群构建####

####响应面优化柱状图####
rm(list=ls())
setwd("E:/wheatrhizosphere/盆栽")


####小麦皿内####
##注：该数据是单个分组的数据
#读入数据：
syncom_wheat <- read.table('./syncom_wheat_rsm.txt',sep="\t", header=T, row.names=1)
head(syncom_wheat)
#View(syncom_wheat) #在Rstudio查看完整数据


##注意：如何给数据框指定的行或列进行求和

#方法一：当数据框较小时，指定的行货列数较少时，手动指定，
library(tidyverse)

syncom_wheat1 = syncom_wheat %>% #管道函数
  mutate(sum = r1+r2+r3+r4+r5+r6+r7+r8+r9+r10+r11+r12+r13+r14+r15+r16+r17+r18+r19+r20+r21+r22+r23+r24)#mutate为数据框添加新的行或列，此处是为了排序，对数据框按行求和，生成新的列sum
syncom_wheat1

syncom_wheat2 <- syncom_wheat1[order(syncom_wheat1[,"sum"],decreasing = T),]#根据数据框指定列对整个数据框进行排序，降序或升序
syncom_wheat2
syncom_wheat2 <- subset(syncom_wheat2, select = -sum)#求和排序后，去除求和列，方便出图
syncom_wheat2
#方法二：rowwise函数与c_across函数结合

#data53= data5 %>%
#  rowwise() %>%
#  mutate(total_length =sum(c_across(c(3:8))))
#data53


#方法二：rowwise函数与rowSums 函数结合
#rowSums 函数既可使用变量名，也可使用变量的位置索引对指定的列进行求和
#data54 = data5 %>%                                                

#  mutate(total_length = rowSums(.[,c(3:8)]))                       

#data54




#用到reshape2包里的melt函数：将宽数据转换为长数据
library(reshape2)

syncom_wheat5 <- melt(syncom_wheat2,id=c('Virus','Cell_lines'))
syncom_wheat5
#把后面三列按照Virus，Cell_lines合并,注意区别melt函数和gather函数

#查看转换结果：
head(syncom_wheat5)

####一、初步绘制柱状图：

ggplot(syncom_wheat5,aes(x=Cell_lines, y=value))+
  stat_summary(geom = 'bar',fun = 'mean',cex=1.3,width=.4)+
  geom_jitter(width = 0.3,size=1.2,color='black',pch=21,fill='white')+
  facet_grid(Virus~.,scales = 'free_y')+
  geom_blank(aes(y=value*1.2))+
  labs(x=NULL,y=NULL)+
  theme_test(base_size = 12)





####二、添加误差棒、修正顺序、设置颜色等：

#先设置两个用于添加误差棒的函数：
topbar <- function(x){      
  return(mean(x)+sd(x)/sqrt(length(x))) #误差采用了mean+-sem
}
bottombar <- function(x){
  return(mean(x)-sd(x)/sqrt(length(x)))
}

#根据数据表格指定的X轴顺序进行排列，比如：经过降序或者升序排列的数据
#本次是将不同复合菌株生物膜由大到小依次排列
celllist <- syncom_wheat5[1:22,]$Cell_lines
celllist <- factor(celllist,levels = celllist)
#按照不同分组来进行分面图的绘制
#Vlist <- c('Plant Height','Fresh Weight', "Dry Weight")
#Vlist <- factor(Vlist,levels = Vlist)


#绘图：
ggplot(syncom_wheat5,aes(factor(Cell_lines,levels = celllist),
                  value,
                  color=factor(Virus),#,levels = Vlist
                  fill=factor(Virus)))+#,levels = Vlist
  stat_summary(geom = 'bar',fun = 'mean',cex=1.3,width=.4)+
  stat_summary(geom = 'errorbar',
               fun.min = bottombar,fun.max = topbar,
               width=.2,cex=0.4,color='black')+ 
  geom_jitter(width = 0.3,size=0.6,color='black',
              pch=21,fill='white')+
  facet_grid(factor(Virus)~.,#,levels = Vlist
  scales = 'free_y')+
  labs(x=NULL,y=NULL)+
  geom_blank(aes(y=value*1.1))+
  theme_test(base_size = 12)+
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 60,vjust = 0.6),
        axis.text = element_text(color='black'),
        strip.background = element_rect(color = 'transparent',
                                        fill = 'transparent'))+
  scale_color_manual(values = "gray50")+
  scale_fill_manual(values = "gray50")
#'#4a8a53','#941319','#E77D13'



ggsave("syncom_wheat2.pdf", width = 24, height = 12, units="cm")
ggsave("syncom_wheat2.tiff", width = 24, height = 12, units = "cm")
ggsave("syncom_wheat2.png", width = 24, height = 12, units = "cm")







####黄瓜皿内####
##注：该数据是单个分组的数据
#读入数据：
syncom_cucumber <- read.table('./syncom_cucumber_rsm.txt',sep="\t", header=T, row.names=1)
head(syncom_cucumber)
#View(syncom_cucumber) #在Rstudio查看完整数据


##注意：如何给数据框指定的行或列进行求和

#方法一：当数据框较小时，指定的行货列数较少时，手动指定，
#library(tidyverse)

syncom_cucumber1 = syncom_cucumber %>% #管道函数
  mutate(sum = r1+r2+r3+r4+r5+r6+r7+r8+r9+r10+r11+r12+r13+r14+r15+r16+r17+r18+r19+r20+r21+r22+r23+r24)#mutate为数据框添加新的行或列，此处是为了排序，对数据框按行求和，生成新的列sum
syncom_cucumber1

syncom_cucumber2 <- syncom_cucumber1[order(syncom_cucumber1[,"sum"],decreasing = T),]#根据数据框指定列对整个数据框进行排序，降序或升序
syncom_cucumber2
syncom_cucumber2 <- subset(syncom_cucumber2, select = -sum)#求和排序后，去除求和列，方便出图
syncom_cucumber2
#方法二：rowwise函数与c_across函数结合

#data53= data5 %>%
#  rowwise() %>%
#  mutate(total_length =sum(c_across(c(3:8))))
#data53


#方法二：rowwise函数与rowSums 函数结合
#rowSums 函数既可使用变量名，也可使用变量的位置索引对指定的列进行求和
#data54 = data5 %>%                                                

#  mutate(total_length = rowSums(.[,c(3:8)]))                       

#data54




#用到reshape2包里的melt函数：将宽数据转换为长数据
library(reshape2)

syncom_cucumber5 <- melt(syncom_cucumber2,id=c('Virus','Cell_lines'))
syncom_cucumber5
#把后面三列按照Virus，Cell_lines合并,注意区别melt函数和gather函数

#查看转换结果：
head(syncom_cucumber5)

####一、初步绘制柱状图：

ggplot(syncom_cucumber5,aes(x=Cell_lines, y=value))+
  stat_summary(geom = 'bar',fun = 'mean',cex=1.3,width=.4)+
  geom_jitter(width = 0.3,size=1.2,color='black',pch=21,fill='white')+
  facet_grid(Virus~.,scales = 'free_y')+
  geom_blank(aes(y=value*1.2))+
  labs(x=NULL,y=NULL)+
  theme_test(base_size = 12)





####二、添加误差棒、修正顺序、设置颜色等：

#先设置两个用于添加误差棒的函数：
topbar <- function(x){      
  return(mean(x)+sd(x)/sqrt(length(x))) #误差采用了mean+-sem
}
bottombar <- function(x){
  return(mean(x)-sd(x)/sqrt(length(x)))
}

#根据数据表格指定的X轴顺序进行排列，比如：经过降序或者升序排列的数据
#本次是将不同复合菌株生物膜由大到小依次排列
celllist <- syncom_cucumber5[1:22,]$Cell_lines
celllist <- factor(celllist,levels = celllist)
#按照不同分组来进行分面图的绘制
#Vlist <- c('Plant Height','Fresh Weight', "Dry Weight")
#Vlist <- factor(Vlist,levels = Vlist)


#绘图：
ggplot(syncom_cucumber5,aes(factor(Cell_lines,levels = celllist),
                         value,
                         color=factor(Virus),#,levels = Vlist
                         fill=factor(Virus)))+#,levels = Vlist
  stat_summary(geom = 'bar',fun = 'mean',cex=1.3,width=.4)+
  stat_summary(geom = 'errorbar',
               fun.min = bottombar,fun.max = topbar,
               width=.2,cex=0.4,color='black')+ 
  geom_jitter(width = 0.3,size=0.6,color='black',
              pch=21,fill='white')+
  facet_grid(factor(Virus)~.,#,levels = Vlist
  scales = 'free_y')+
  labs(x=NULL,y=NULL)+
  geom_blank(aes(y=value*1.1))+
  theme_test(base_size = 12)+
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 60,vjust = 0.6),
        axis.text = element_text(color='black'),
        strip.background = element_rect(color = 'transparent',
                                        fill = 'transparent'))+
  scale_color_manual(values = "gray50")+
  scale_fill_manual(values = "gray50")
#'#4a8a53','#941319','#E77D13'



ggsave("syncom_cucumber.pdf", width = 12, height = 6, units="cm")
ggsave("syncom_cucumber.tiff", width = 12, height = 6, units = "cm")
ggsave("syncom_cucumber.png", width = 12, height = 6, units = "cm")














































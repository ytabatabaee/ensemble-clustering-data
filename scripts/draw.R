require(ggplot2);require(reshape2);require(scales);require(ggpubr);require(tidyr)

#res_limit_exps= read.csv('res_limit_exps.csv')
res_limit_exps= read.csv('res_limit_exps_leiden_cpm.csv')
res_limit_exps= read.csv('res_limit_exps_leiden_cpm_vary_res_final.csv')
#res_limit_exps= read.csv('res_limit_exps_leiden_mod_tree.csv')
#res_limit_exps$partition <- factor(res_limit_exps$partition, levels = c("Leiden-CPM(r=0.0001)", "SC(np=10)+Leiden-CPM(r=0.0001)", "SC(np=50)+Leiden-CPM(r=0.0001)", "SC(np=100)+Leiden-CPM(r=0.0001)"))

give.n <- function(x){
  return(c(y = median(x)*1.05, label = length(x))) 
}

res_limit_exps= read.csv('res_limit_exps_leiden_cpm_k.csv')
ggplot(aes(x= as.factor(n),y=cluster_size,fill=partition, color=partition), data=res_limit_exps)+
  #facet_wrap(~method,ncol=2)+
  geom_boxplot()+
  #geom_boxplot(outlier.size = 0)+
  #stat_summary(fun.data = mean_sdl,position = position_dodge(width=0.75))+
  stat_summary(fun.data = give.n, geom = "text", vjust=-1.5, position = position_dodge(width=0.75), col="black", size=2)+
  scale_y_continuous(name="Cluster size distribution")+
  scale_x_discrete(name="Size of clique")+
  #scale_x_discrete(name="Resolution value")+
  scale_fill_brewer(palette="Set2")+
  scale_color_brewer(palette = "Dark2")+
  theme_bw()+
  theme(legend.position='none', legend.direction = "vertical",legend.title=element_blank())
ggsave("res_limit_exps_leiden_cpm_k.pdf",width=7.8,height =3.3)

ggplot(aes(x= as.factor(n),y=cluster_size,fill=partition, color=partition), data=res_limit_exps)+
  #facet_wrap(~method,ncol=2)+
  geom_boxplot()+
  #geom_boxplot(outlier.size = 0)+
  #stat_summary(fun.data = mean_sdl,position = position_dodge(width=0.75))+
  stat_summary(fun.data = give.n, geom = "text", vjust=-1.5, position = position_dodge(width=0.75), col="black", size=2.3)+
  scale_y_continuous(name="Cluster size distribution")+
  scale_x_discrete(name="Number of cliques of size 10")+
  #scale_x_discrete(name="Resolution value")+
  scale_fill_brewer(palette="Set2")+
  scale_color_brewer(palette = "Dark2")+
  theme_bw()+
  guides(color = guide_legend(nrow = 2))+
  theme(legend.position='bottom', legend.direction = "vertical",legend.title=element_blank())
ggsave("res_limit_exps_leiden_cpm_k_legend.pdf",width=7.8,height =3.3)

res_limit_exps= read.csv('res_limit_exps_leiden_cpm_acc_k.csv')

ggplot(aes(x=as.factor(n), y=acc_value, fill=partition, color=partition, group=partition), data=res_limit_exps[(res_limit_exps$acc_measure=='NMI'),])+
  #facet_wrap(~acc_measure,ncol=3)+
  geom_point()+geom_line()+
  #stat_summary(fun.data = mean_sdl,position = position_dodge(width=0.75))+
  #stat_summary(fun.data = give.n, geom = "text", vjust=-1.5, position = position_dodge(width=0.75), col="black", size=3)+
  scale_y_continuous(name="NMI")+
  scale_x_discrete(name="Size of clique")+
  #scale_x_discrete(name="Resolution value")+
  #coord_cartesian(ylim=c(0,1))+
  scale_fill_brewer(palette="Set2")+
  scale_color_brewer(palette = "Dark2")+
  theme_bw()+
  theme(legend.position = "none", legend.direction = "horizontal",
        axis.text.x = element_text(angle=0))
ggsave("res_limit_exps_leiden_cpm_acc_k_nmi.pdf",width=2.5,height=2)

ggplot(aes(x=as.factor(n), y=acc_value, fill=partition, color=partition, group=partition), data=res_limit_exps[(res_limit_exps$acc_measure=='ARI'),])+
  #facet_wrap(~acc_measure,ncol=3)+
  geom_point()+geom_line()+
  #stat_summary(fun.data = mean_sdl,position = position_dodge(width=0.75))+
  #stat_summary(fun.data = give.n, geom = "text", vjust=-1.5, position = position_dodge(width=0.75), col="black", size=3)+
  scale_y_continuous(name="ARI")+
  scale_x_discrete(name="Size of clique")+
  coord_cartesian(ylim=c(0,1))+
  scale_fill_brewer(palette="Set2")+
  scale_color_brewer(palette = "Dark2")+
  theme_bw()+
  theme(legend.position = "none", legend.direction = "horizontal",
        axis.text.x = element_text(angle=0))
ggsave("res_limit_exps_leiden_cpm_acc_k_ari.pdf",width=2.5,height=2)

ggplot(aes(x=as.factor(k), y=acc_value, fill=partition, color=partition, group=partition), data=res_limit_exps[(res_limit_exps$acc_measure=='NMI'),])+
  #facet_wrap(~acc_measure,ncol=3)+
  geom_point()+geom_line()+
  #stat_summary(fun.data = mean_sdl,position = position_dodge(width=0.75))+
  #stat_summary(fun.data = give.n, geom = "text", vjust=-1.5, position = position_dodge(width=0.75), col="black", size=3)+
  scale_y_continuous(name="NMI")+
  scale_x_discrete(name="Number of cliques of size 10")+
  #scale_x_discrete(name="Resolution value")+
  #coord_cartesian(ylim=c(0,1))+
  scale_fill_brewer(palette="Set2")+
  scale_color_brewer(palette = "Dark2")+
  theme_bw()+
  theme(legend.position = "none", legend.direction = "horizontal",
        axis.text.x = element_text(angle=30))
ggsave("res_limit_exps_leiden_mod_nmi.pdf",width=2.5,height=2)

ggplot(aes(x=as.factor(k), y=acc_value, fill=partition, color=partition, group=partition), data=res_limit_exps[(res_limit_exps$acc_measure=='ARI'),])+
  #facet_wrap(~acc_measure,ncol=3)+
  geom_point()+geom_line()+
  #stat_summary(fun.data = mean_sdl,position = position_dodge(width=0.75))+
  #stat_summary(fun.data = give.n, geom = "text", vjust=-1.5, position = position_dodge(width=0.75), col="black", size=3)+
  scale_y_continuous(name="ARI")+
  scale_x_discrete(name="Number of cliques of size 10")+
  coord_cartesian(ylim=c(0,1))+
  scale_fill_brewer(palette="Set2")+
  scale_color_brewer(palette = "Dark2")+
  theme_bw()+
  theme(legend.position = "none", legend.direction = "horizontal",
        axis.text.x = element_text(angle=30))
ggsave("res_limit_exps_leiden_mod_ari.pdf",width=2.5,height=2)

res_limit_exps= read.csv('res_limit_exps_leiden_mod_tree_all.csv')
#res_limit_exps = res_limit_exps[res_limit_exps$partition!="ECG",]
#res_limit_exps = res_limit_exps[res_limit_exps$partition!="FastConsensus(Louvain)",]
#res_limit_exps = res_limit_exps[res_limit_exps$partition!="Strict(np=10,Leiden-mod)",]
#res_limit_exps = res_limit_exps[res_limit_exps$partition!="Strict(np=50,Leiden-mod)",]
res_limit_exps = res_limit_exps[res_limit_exps$partition!="FastEnsemble(Louvain)",]
res_limit_exps = res_limit_exps[res_limit_exps$partition!="Louvain",]
res_limit_exps$partition = factor(res_limit_exps$partition, levels=c('ECG','FastEnsemble(Leiden-mod)','Leiden-mod','FastConsensus(Louvain)','Strict(np=10,Leiden-mod)','Strict(np=50,Leiden-mod)'))

ggplot(aes(x= as.factor(k),y=cluster_size,fill=partition, color=partition), data=res_limit_exps)+
  #facet_wrap(~method,ncol=2)+
  geom_boxplot()+
  #geom_boxplot(outlier.size = 0)+
  #stat_summary(fun.data = mean_sdl,position = position_dodge(width=0.75))+
  #stat_summary(fun.data = give.n, geom = "text", vjust=-1.5, position = position_dodge(width=0.75), col="black", size=2.3)+
  scale_y_continuous(name="Cluster size distribution", breaks=c(seq(from=0,to=180,by=50), 10))+
  scale_x_discrete(name="Number of cliques of size 10")+
  geom_hline(yintercept=10, linetype="dotted", color = "grey40")+
  coord_cartesian(ylim=c(0,180))+
  #scale_x_discrete(name="Resolution value")+
  scale_fill_brewer(palette="Set2")+
  scale_color_brewer(palette = "Dark2")+
  theme_bw()+
  theme(legend.position=c(0.2,0.6), legend.direction = "vertical",legend.title=element_blank())
ggsave("res_limit_exps_leiden_mod_tree.pdf",width=8,height =3.3)

res_limit_exps= read.csv('res_limit_exps_leiden_mod_tree_acc_all.csv')
res_limit_exps = res_limit_exps[res_limit_exps$partition!="FastEnsemble(Louvain)",]
res_limit_exps = res_limit_exps[res_limit_exps$partition!="Louvain",]
#res_limit_exps = res_limit_exps[res_limit_exps$partition!="ECG",]
#res_limit_exps = res_limit_exps[res_limit_exps$partition!="FastConsensus(Louvain)",]
#res_limit_exps = res_limit_exps[res_limit_exps$partition!="Strict(np=10,Leiden-mod)",]
#res_limit_exps = res_limit_exps[res_limit_exps$partition!="Strict(np=50,Leiden-mod)",]
res_limit_exps$partition = factor(res_limit_exps$partition, levels=c('ECG','FastEnsemble(Leiden-mod)','Leiden-mod','FastConsensus(Louvain)','Strict(np=10,Leiden-mod)','Strict(np=50,Leiden-mod)'))

ggplot(aes(x=as.factor(k), y=acc_value, fill=partition, color=partition, group=partition), data=res_limit_exps[(res_limit_exps$acc_measure=='NMI'),])+
  #facet_wrap(~acc_measure,ncol=3)+
  geom_point()+geom_line()+
  #stat_summary(fun.data = mean_sdl,position = position_dodge(width=0.75))+
  #stat_summary(fun.data = give.n, geom = "text", vjust=-1.5, position = position_dodge(width=0.75), col="black", size=3)+
  scale_y_continuous(name="NMI")+
  scale_x_discrete(name="Number of cliques of size 10")+
  #scale_x_discrete(name="Resolution value")+
  coord_cartesian(ylim=c(0,1))+
  scale_fill_brewer(palette="Set2")+
  scale_color_brewer(palette = "Dark2")+
  theme_bw()+
  theme(legend.position = "none", legend.direction = "horizontal",
        axis.text.x = element_text(angle=30))
ggsave("res_limit_exps_leiden_mod_tree_nmi.pdf",width=2.5,height=2)

ggplot(aes(x=as.factor(k), y=acc_value, fill=partition, color=partition, group=partition), data=res_limit_exps[(res_limit_exps$acc_measure=='ARI'),])+
  #facet_wrap(~acc_measure,ncol=3)+
  geom_point()+geom_line()+
  #stat_summary(fun.data = mean_sdl,position = position_dodge(width=0.75))+
  #stat_summary(fun.data = give.n, geom = "text", vjust=-1.5, position = position_dodge(width=0.75), col="black", size=3)+
  scale_y_continuous(name="ARI")+
  scale_x_discrete(name="Number of cliques of size 10")+
  #scale_x_discrete(name="Resolution value")+
  coord_cartesian(ylim=c(0,1))+
  scale_fill_brewer(palette="Set2")+
  scale_color_brewer(palette = "Dark2")+
  theme_bw()+
  theme(legend.position = "none", legend.direction = "horizontal",
        axis.text.x = element_text(angle=30))
ggsave("res_limit_exps_leiden_mod_tree_ari.pdf",width=2.5,height=2)

  facet_wrap(~acc_measure,ncol=3)+
  geom_point()+geom_line()+
  #stat_summary(fun.data = mean_sdl,position = position_dodge(width=0.75))+
  #stat_summary(fun.data = give.n, geom = "text", vjust=-1.5, position = position_dodge(width=0.75), col="black", size=3)+
  scale_y_continuous(name="Accuracy")+
  #scale_x_discrete(name="Number of cliques of size 10")+
  #scale_x_discrete(name="Resolution value")+
  coord_cartesian(ylim=c(0,1))+
  scale_fill_brewer(palette="Set2")+
  scale_color_brewer(palette = "Dark2")+
  theme_bw()+
  theme(legend.position = "none", legend.direction = "horizontal",
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle=30))
ggsave("res_limit_exps_lieden_mod_tree_accuracy.pdf",width=5,height=2)


ggplot(aes(x=as.factor(k), y=acc_value, fill=partition, color=partition, group=partition), data=subset(res_limit_exps, acc_measure %in% c('FNR', 'FPR')))+
  facet_wrap(~acc_measure,ncol=3)+
  geom_point()+geom_line()+
  #stat_summary(fun.data = mean_sdl,position = position_dodge(width=0.75))+
  #stat_summary(fun.data = give.n, geom = "text", vjust=-1.5, position = position_dodge(width=0.75), col="black", size=3)+
  scale_y_continuous(name="Error")+
  scale_x_discrete(name="")+
  #scale_x_discrete(name="Resolution value")+
  #coord_cartesian(ylim=c(0,1))+
  scale_fill_brewer(palette="Set2")+
  scale_color_brewer(palette = "Dark2")+
  theme_bw()+
  theme(legend.position = "none", legend.direction = "horizontal",
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle=30))
ggsave("res_limit_exps_lieden_mod_tree_fpr_fnr.pdf",width=5,height=2)


res_limit_exps= read.csv('res_limit_exps_leiden_mod_ring.csv')
# making sure we get results with the default parameter (tr=0.8)
#res_limit_exps$partition[res_limit_exps$partition == 'FastEnsemble(Leiden-mod)'] <- 'FastEnsemble(Leiden-mod,tr=0.9)'
res_limit_exps = res_limit_exps[res_limit_exps$partition!="FastEnsemble(Leiden-mod)",]
res_limit_exps$partition[res_limit_exps$partition == 'FastEnsemble(Leiden-mod,tr=0.8)'] <- 'FastEnsemble(Leiden-mod)'
res_limit_exps$partition = factor(res_limit_exps$partition, levels=c('ECG','FastEnsemble(Leiden-mod)','Leiden-mod','FastConsensus(Louvain)','Strict(np=10,Leiden-mod)','Strict(np=50,Leiden-mod)'))

# for Min's talk at CNA
#res_limit_exps = res_limit_exps[res_limit_exps$k!=500,]
#res_limit_exps = res_limit_exps[res_limit_exps$k!=5000,]

ggplot(aes(x= as.factor(k),y=cluster_size,fill=partition, color=partition), data=subset(res_limit_exps, partition %in% c('ECG', 'FastEnsemble(Leiden-mod)', 'Leiden-mod', 'FastConsensus(Louvain)', 'Strict(np=10,Leiden-mod)','Strict(np=50,Leiden-mod)')))+
  #facet_wrap(~method,ncol=2)+
  geom_boxplot()+
  #geom_boxplot(outlier.size = 0)+
  #stat_summary(fun.data = mean_sdl,position = position_dodge(width=0.75))+
  #stat_summary(fun.data = give.n, geom = "text", vjust=-1.5, position = position_dodge(width=0.75), col="black", size=2.3)+
  scale_y_continuous(name="Cluster size distribution", breaks=c(seq(from=0,to=200,by=50), 10))+
  coord_cartesian(ylim=c(0,200))+
  scale_x_discrete(name="Number of cliques of size 10")+
  geom_hline(yintercept=10, linetype="dotted", color = "grey40")+
  #scale_x_discrete(name="Resolution value")+
  scale_fill_brewer(palette="Set2")+
  scale_color_brewer(palette = "Dark2")+
  theme_bw()+
  theme(legend.position=c(0.2,0.6), legend.direction = "vertical",legend.title=element_blank())
ggsave("res_limit_exps_leiden_mod_ring.pdf",width=8,height =3.3)

ggplot(aes(x= as.factor(k),y=cluster_size,fill=partition, color=partition), data=res_limit_exps)+
  #facet_wrap(~method,ncol=2)+
  geom_boxplot()+
  #stat_summary(fun.data = mean_sdl,position = position_dodge(width=0.75))+
  stat_summary(fun.data = give.n, geom = "text", vjust=-1.5, position = position_dodge(width=0.75), col="black", size=3)+
  scale_y_continuous(name="Cluster size distribution")+
  scale_x_discrete(name="Number of cliques of size 10")+
  #scale_x_discrete(name="Resolution value")+
  scale_fill_brewer(palette="Set2")+
  scale_color_brewer(palette = "Dark2")+
  theme_bw()+
  theme(legend.position='none', legend.direction = "horizontal", nrow=2)+
  guides(color=guide_legend(nrow=2, byrow=TRUE))
ggsave("res_limit_exps_leiden_cpm.pdf",width=7.4,height =3.3)

res_limit_exps= read.csv('res_limit_exps_leiden_cpm_vary_res_final.csv')
res_limit_exps$partition <- factor(res_limit_exps$partition, levels=c('Leiden-CPM', 'FastEnsemble(Leiden-CPM)', 'FastEnsemble(Leiden-CPM,tr=0.9)', 'Strict(np=10,Leiden-CPM)', 'Strict(np=50,Leiden-CPM)'))
res_limit_exps = res_limit_exps[res_limit_exps$partition!="FastEnsemble(Leiden-CPM,tr=0.9)",]


ggplot(aes(x= as.factor(res),y=cluster_size,fill=partition, color=partition), data=res_limit_exps)+
  #facet_wrap(~method,ncol=2)+
  geom_boxplot()+
  #stat_summary(fun.data = mean_sdl,position = position_dodge(width=0.75))+
  #stat_summary(fun.data = give.n, geom = "text", vjust=-1.5, position = position_dodge(width=0.75), col="black", size=2)+
  scale_y_continuous(name="Cluster size distribution", breaks=c(seq(from=0,to=570,by=200)))+
  geom_hline(yintercept=10, linetype="dotted", color = "grey40")+
  #geom_text(aes(color='black',y=0.5,x=0.5,label="1000 cliques of size 10"),size=2)+
  #scale_x_discrete(name="Number of cliques of size 10")+
  coord_cartesian(ylim=c(0,570))+
  scale_x_discrete(name="Resolution value")+
  #scale_fill_brewer(palette="Set2")+
  #scale_color_brewer(palette = "Dark2")+
  scale_fill_manual(values=c("#8da0cb","#fc8d62", "#a6d854", "#ffd92f"), name="")+
  scale_color_manual(values=c("#7570b3","#d95f02", "#66a61e", "#e6ab02"), name="")+
  theme_bw()+
  theme(legend.position=c(0.4,0.75), legend.direction = "vertical",legend.title=element_blank())
ggsave("res_limit_exps_leiden_cpm_vary_res_main.pdf",width=7.8,height =3.4)

res_limit_exps= read.csv('res_limit_exps_leiden_cpm_final.csv')
res_limit_exps$partition <- factor(res_limit_exps$partition, levels=c('Leiden-CPM', 'FastEnsemble(Leiden-CPM)', 'FastEnsemble(Leiden-CPM,tr=0.9)', 'Strict(np=10,Leiden-CPM)', 'Strict(np=50,Leiden-CPM)'))
res_limit_exps = res_limit_exps[res_limit_exps$partition!="FastEnsemble(Leiden-CPM,tr=0.9)",]

ggplot(aes(x= as.factor(k),y=cluster_size,fill=partition, color=partition), data=res_limit_exps)+
  #facet_wrap(~method,ncol=2)+
  geom_boxplot()+
  #stat_summary(fun.data = mean_sdl,position = position_dodge(width=0.75))+
  #stat_summary(fun.data = give.n, geom = "text", vjust=-1.5, position = position_dodge(width=0.75), col="black", size=2.1)+
  scale_y_continuous(name="Cluster size distribution", breaks=c(seq(from=0,to=60,by=20), 10))+
  #geom_text(aes(color='black',y=0.5,x=0.5,label="1000 cliques of size 10"),size=2)+
  scale_x_discrete(name="Number of cliques of size 10")+
  #scale_x_discrete(name="Resolution value")+
  geom_hline(yintercept=10, linetype="dotted", color = "grey40")+
  #scale_fill_brewer(palette="Set2")+
  scale_fill_manual(values=c("#8da0cb","#fc8d62", "#a6d854", "#ffd92f"), name="")+
  scale_color_manual(values=c("#7570b3","#d95f02", "#66a61e", "#e6ab02"), name="")+
  coord_cartesian(ylim=c(0,60))+
  #scale_color_brewer(palette = "Dark2")+
  theme_bw()+
  theme(legend.position='none', legend.direction = "vertical",legend.title=element_blank())
ggsave("res_limit_exps_leiden_cpm_main.pdf",width=7.8,height =3.4)

res_limit_exps= read.csv('res_limit_exps_leiden_cpm_vary_res_acc_final.csv')
res_limit_exps$partition <- factor(res_limit_exps$partition, levels=c('Leiden-CPM', 'FastEnsemble(Leiden-CPM)', 'FastEnsemble(Leiden-CPM,tr=0.9)', 'Strict(np=10,Leiden-CPM)', 'Strict(np=50,Leiden-CPM)'))
res_limit_exps = res_limit_exps[res_limit_exps$partition!="FastEnsemble(Leiden-CPM,tr=0.9)",]

ggplot(aes(x=as.factor(res), y=acc_value, fill=partition, color=partition, group=partition), data=res_limit_exps[(res_limit_exps$acc_measure=='NMI'),])+
  #facet_wrap(~acc_measure,ncol=3)+
  geom_point()+geom_line()+
  #stat_summary(fun.data = mean_sdl,position = position_dodge(width=0.75))+
  #stat_summary(fun.data = give.n, geom = "text", vjust=-1.5, position = position_dodge(width=0.75), col="black", size=3)+
  scale_y_continuous(name="NMI")+
  #scale_x_discrete(name="Number of cliques of size 10")+
  scale_x_discrete(name="Resolution value")+
  coord_cartesian(ylim=c(0,1))+
  #scale_fill_brewer(palette="Set2")+
  #scale_color_brewer(palette = "Dark2")+
  scale_fill_manual(values=c("#8da0cb","#fc8d62", "#a6d854", "#ffd92f"), name="")+
  scale_color_manual(values=c("#7570b3","#d95f02", "#66a61e", "#e6ab02"), name="")+
  theme_bw()+
  theme(legend.position = "none", legend.direction = "horizontal",
        axis.text.x = element_text(angle=30))
ggsave("res_limit_exps_leiden_cpm_nmi_all.pdf",width=2.5,height=2)

ggplot(aes(x=as.factor(res), y=acc_value, fill=partition, color=partition, group=partition), data=res_limit_exps[(res_limit_exps$acc_measure=='ARI'),])+
  #facet_wrap(~acc_measure,ncol=3)+
  geom_point()+geom_line()+
  #stat_summary(fun.data = mean_sdl,position = position_dodge(width=0.75))+
  #stat_summary(fun.data = give.n, geom = "text", vjust=-1.5, position = position_dodge(width=0.75), col="black", size=3)+
  scale_y_continuous(name="ARI")+
  #scale_x_discrete(name="Number of cliques of size 10")+
  scale_x_discrete(name="Resolution value")+
  coord_cartesian(ylim=c(0,1))+
  #scale_fill_brewer(palette="Set2")+
  #scale_color_brewer(palette = "Dark2")+
  scale_fill_manual(values=c("#8da0cb","#fc8d62", "#a6d854", "#ffd92f"), name="")+
  scale_color_manual(values=c("#7570b3","#d95f02", "#66a61e", "#e6ab02"), name="")+
  theme_bw()+
  theme(legend.position = "none", legend.direction = "horizontal",
        axis.text.x = element_text(angle=30))
ggsave("res_limit_exps_leiden_cpm_ari_all.pdf",width=2.5,height=2)

res_limit_exps= read.csv('res_limit_exps_leiden_cpm_acc_final.csv')
res_limit_exps$partition <- factor(res_limit_exps$partition, levels=c('Leiden-CPM', 'FastEnsemble(Leiden-CPM)', 'FastEnsemble(Leiden-CPM,tr=0.9)', 'Strict(np=10,Leiden-CPM)', 'Strict(np=50,Leiden-CPM)'))
res_limit_exps = res_limit_exps[res_limit_exps$partition!="FastEnsemble(Leiden-CPM,tr=0.9)",]

ggplot(aes(x=as.factor(k), y=acc_value, fill=partition, color=partition, group=partition), data=res_limit_exps[(res_limit_exps$acc_measure=='NMI'),])+
  #facet_wrap(~acc_measure,ncol=3)+
  geom_point()+geom_line()+
  #stat_summary(fun.data = mean_sdl,position = position_dodge(width=0.75))+
  #stat_summary(fun.data = give.n, geom = "text", vjust=-1.5, position = position_dodge(width=0.75), col="black", size=3)+
  scale_y_continuous(name="NMI")+
  scale_x_discrete(name="Number of cliques of size 10")+
  #scale_x_discrete(name="Resolution value")+
  coord_cartesian(ylim=c(0,1))+
  #scale_fill_brewer(palette="Set2")+
  #scale_color_brewer(palette = "Dark2")+
  scale_fill_manual(values=c("#8da0cb","#fc8d62", "#a6d854", "#ffd92f"), name="")+
  scale_color_manual(values=c("#7570b3","#d95f02", "#66a61e", "#e6ab02"), name="")+
  theme_bw()+
  theme(legend.position = "none", legend.direction = "horizontal",
        axis.text.x = element_text(angle=30))
ggsave("res_limit_exps_leiden_cpm_nmi_main2.pdf",width=2.5,height=2)

ggplot(aes(x=as.factor(k), y=acc_value, fill=partition, color=partition, group=partition), data=res_limit_exps[(res_limit_exps$acc_measure=='ARI'),])+
  #facet_wrap(~acc_measure,ncol=3)+
  geom_point()+geom_line()+
  #stat_summary(fun.data = mean_sdl,position = position_dodge(width=0.75))+
  #stat_summary(fun.data = give.n, geom = "text", vjust=-1.5, position = position_dodge(width=0.75), col="black", size=3)+
  scale_y_continuous(name="ARI")+
  scale_x_discrete(name="Number of cliques of size 10")+
  #scale_x_discrete(name="Resolution value")+
  coord_cartesian(ylim=c(0,1))+
  #scale_fill_brewer(palette="Set2")+
  #scale_color_brewer(palette = "Dark2")+
  scale_fill_manual(values=c("#8da0cb","#fc8d62", "#a6d854", "#ffd92f"), name="")+
  scale_color_manual(values=c("#7570b3","#d95f02", "#66a61e", "#e6ab02"), name="")+
  theme_bw()+
  theme(legend.position = "none", legend.direction = "horizontal",
        axis.text.x = element_text(angle=30))
ggsave("res_limit_exps_leiden_cpm_ari_main2.pdf",width=2.5,height=2)


res_limit_exps= read.csv('res_limit_exps_leiden_mod_ring_acc.csv')
res_limit_exps$partition[res_limit_exps$partition == 'FastEnsemble(Leiden-mod)'] <- 'FastEnsemble(Leiden-mod,tr=0.9)'
# making sure we get results with the default parameter (tr=0.8)
res_limit_exps = res_limit_exps[res_limit_exps$partition!="FastEnsemble(Leiden-mod,tr=0.9)",]
#res_limit_exps = res_limit_exps[res_limit_exps$partition!="Strict(np=10,Leiden-mod)",]
#res_limit_exps = res_limit_exps[res_limit_exps$partition!="Strict(np=50,Leiden-mod)",]
res_limit_exps$partition[res_limit_exps$partition == 'FastEnsemble(Leiden-mod,tr=0.8)'] <- 'FastEnsemble(Leiden-mod)'
res_limit_exps$partition = factor(res_limit_exps$partition, levels=c('ECG','FastEnsemble(Leiden-mod)','Leiden-mod','FastConsensus(Louvain)','Strict(np=10,Leiden-mod)','Strict(np=50,Leiden-mod)'))
# for Min's talk at CNA
#res_limit_exps = res_limit_exps[res_limit_exps$k!=500,]
#res_limit_exps = res_limit_exps[res_limit_exps$k!=5000,]

ggplot(aes(x=as.factor(k), y=acc_value, fill=partition, color=partition, group=partition), data=res_limit_exps[(res_limit_exps$acc_measure=='NMI'),])+
  #facet_wrap(~acc_measure,ncol=3)+
  geom_point()+geom_line()+
  #stat_summary(fun.data = mean_sdl,position = position_dodge(width=0.75))+
  #stat_summary(fun.data = give.n, geom = "text", vjust=-1.5, position = position_dodge(width=0.75), col="black", size=3)+
  scale_y_continuous(name="NMI")+
  scale_x_discrete(name="Number of cliques of size 10")+
  #scale_x_discrete(name="Resolution value")+
  coord_cartesian(ylim=c(0,1))+
  scale_fill_brewer(palette="Set2")+
  scale_color_brewer(palette = "Dark2")+
  theme_bw()+
  theme(legend.position = "none", legend.direction = "horizontal",
        axis.text.x = element_text(angle=30))
ggsave("res_limit_exps_leiden_mod_nmi.pdf",width=2.5,height=2)

ggplot(aes(x=as.factor(k), y=acc_value, fill=partition, color=partition, group=partition), data=res_limit_exps[(res_limit_exps$acc_measure=='ARI'),])+
  #facet_wrap(~acc_measure,ncol=3)+
  geom_point()+geom_line()+
  #stat_summary(fun.data = mean_sdl,position = position_dodge(width=0.75))+
  #stat_summary(fun.data = give.n, geom = "text", vjust=-1.5, position = position_dodge(width=0.75), col="black", size=3)+
  scale_y_continuous(name="ARI")+
  scale_x_discrete(name="Number of cliques of size 10")+
  coord_cartesian(ylim=c(0,1))+
  scale_fill_brewer(palette="Set2")+
  scale_color_brewer(palette = "Dark2")+
  theme_bw()+
  theme(legend.position = "none", legend.direction = "horizontal",
        axis.text.x = element_text(angle=30))
ggsave("res_limit_exps_leiden_mod_ari.pdf",width=2.5,height=2)

ggplot(aes(x=as.factor(k), y=acc_value, fill=partition, color=partition, group=partition), data=subset(res_limit_exps,acc_measure %in% c('NMI', 'ARI', 'F1-Score') & partition %in% c('Leiden-mod','FastEnsemble(Leiden-mod)','FastEnsemble(Leiden-mod,tr=0.9)','Louvain','FastEnsemble(Louvain)','Strict(np=10,Leiden-mod)','Strict(np=50,Leiden-mod)')))+ # 
  facet_wrap(~acc_measure,ncol=3)+
  geom_point()+geom_line()+
  #geom_bar(stat='identity', position = position_dodge2(preserve = "single"))+
  #stat_summary(fun.data = mean_sdl,position = position_dodge(width=0.75))+
  #stat_summary(fun.data = give.n, geom = "text", vjust=-1.5, position = position_dodge(width=0.75), col="black", size=3)+
  scale_y_continuous(name="Accuracy")+
  scale_x_discrete(name="Number of cliques of size 10")+
  #scale_x_discrete(name="Resolution value")+
  coord_cartesian(ylim=c(0,1))+
  scale_fill_brewer(palette="Set2")+
  scale_color_brewer(palette = "Dark2")+
  theme_bw()+
  guides(color=guide_legend(nrow=3, byrow=TRUE)) +
  theme(legend.position = "bottom", legend.direction = "horizontal", legend.title=element_blank(),
        axis.text.x = element_text(angle=30))
ggsave("res_limit_exps_leiden_mod_nmi_louvain.pdf",width=8.3,height=3.8)

ggplot(aes(x=as.factor(k), y=acc_value, fill=partition, color=partition, group=partition), data=res_limit_exps[(res_limit_exps$acc_measure=='ARI'),])+
  #facet_wrap(~acc_measure,ncol=3)+
  geom_point()+geom_line()+
  #stat_summary(fun.data = mean_sdl,position = position_dodge(width=0.75))+
  #stat_summary(fun.data = give.n, geom = "text", vjust=-1.5, position = position_dodge(width=0.75), col="black", size=3)+
  scale_y_continuous(name="ARI")+
  scale_x_discrete(name="Number of cliques of size 10")+
  coord_cartesian(ylim=c(0,1))+
  scale_fill_brewer(palette="Set2")+
  scale_color_brewer(palette = "Dark2")+
  theme_bw()+
  theme(legend.position = "none", legend.direction = "horizontal",
        axis.text.x = element_text(angle=30))
ggsave("res_limit_exps_leiden_mod_ari.pdf",width=2.5,height=2)

res_limit_exps= read.csv('res_limit_exps_leiden_cpm_vary_res_acc_final.csv')
res_limit_exps$partition <- factor(res_limit_exps$partition, levels=c('Leiden-CPM', 'FastEnsemble(Leiden-CPM)', 'FastEnsemble(Leiden-CPM,tr=0.9)', 'Strict(np=10,Leiden-CPM)', 'Strict(np=50,Leiden-CPM)'))
res_limit_exps= read.csv('res_limit_exps_leiden_cpm_acc_final.csv')

ggplot(aes(x=as.factor(res), y=acc_value, fill=partition, color=partition, group=partition), data=res_limit_exps[(res_limit_exps$acc_measure=='NMI'),])+
  #facet_wrap(~acc_measure,ncol=3)+
  geom_point()+geom_line()+
  #stat_summary(fun.data = mean_sdl,position = position_dodge(width=0.75))+
  #stat_summary(fun.data = give.n, geom = "text", vjust=-1.5, position = position_dodge(width=0.75), col="black", size=3)+
  scale_y_continuous(name="NMI")+
  scale_x_discrete(name="Number of Cliques")+
  #coord_cartesian(ylim=c(0,1))+
  scale_fill_brewer(palette="Set2")+
  scale_color_brewer(palette = "Dark2")+
  theme_bw()+
  theme(legend.position = "none", legend.direction = "horizontal",
        axis.text.x = element_text(angle=30))
ggsave("res_limit_exps_vary_res_leiden_cpm_nmi.pdf",width=2.5,height=2)

ggplot(aes(x=as.factor(res), y=acc_value, fill=partition, color=partition, group=partition), data=res_limit_exps[(res_limit_exps$acc_measure=='ARI'),])+
  #facet_wrap(~acc_measure,ncol=3)+
  geom_point()+geom_line()+
  #stat_summary(fun.data = mean_sdl,position = position_dodge(width=0.75))+
  #stat_summary(fun.data = give.n, geom = "text", vjust=-1.5, position = position_dodge(width=0.75), col="black", size=3)+
  scale_y_continuous(name="ARI")+
  scale_x_discrete(name="Resolution value")+
  #scale_x_discrete(name="Number of Cliques")+
  coord_cartesian(ylim=c(0,1))+
  scale_fill_brewer(palette="Set2")+
  scale_color_brewer(palette = "Dark2")+
  theme_bw()+
  theme(legend.position = "none", legend.direction = "horizontal",
        axis.text.x = element_text(angle=30))
ggsave("res_limit_exps_vary_res_leiden_cpm_ari.pdf",width=2.5,height=2)

ggplot(aes(x=as.factor(res), y=acc_value, fill=partition, color=partition, group=partition), data=res_limit_exps[(res_limit_exps$acc_measure=='NMI'),])+
  #facet_wrap(~acc_measure,ncol=3)+
  geom_point()+geom_line()+
  #stat_summary(fun.data = mean_sdl,position = position_dodge(width=0.75))+
  #stat_summary(fun.data = give.n, geom = "text", vjust=-1.5, position = position_dodge(width=0.75), col="black", size=3)+
  scale_y_continuous(name="NMI")+
  scale_x_discrete(name="Resolution value")+
  #coord_cartesian(ylim=c(0,1))+
  scale_fill_brewer(palette="Set2")+
  scale_color_brewer(palette = "Dark2")+
  theme_bw()+
  theme(legend.position = "none", legend.direction = "horizontal",
        axis.text.x = element_text(angle=30))
ggsave("res_limit_exps_leiden_cpm_vary_res_nmi.pdf",width=2.5,height=2)

ggplot(aes(x=as.factor(res), y=acc_value, fill=partition, color=partition, group=partition), data=res_limit_exps[(res_limit_exps$acc_measure=='ARI'),])+
  #facet_wrap(~acc_measure,ncol=3)+
  geom_point()+geom_line()+
  #stat_summary(fun.data = mean_sdl,position = position_dodge(width=0.75))+
  #stat_summary(fun.data = give.n, geom = "text", vjust=-1.5, position = position_dodge(width=0.75), col="black", size=3)+
  scale_y_continuous(name="ARI")+
  scale_x_discrete(name="Resolution value")+
  coord_cartesian(ylim=c(0,1))+
  scale_fill_brewer(palette="Set2")+
  scale_color_brewer(palette = "Dark2")+
  theme_bw()+
  theme(legend.position = "none", legend.direction = "horizontal",
        axis.text.x = element_text(angle=30))
ggsave("res_limit_exps_leiden_cpm_vary_res_ari.pdf",width=2.5,height=2)



ggplot(aes(x=as.factor(n), y=acc_value, fill=partition, color=partition), data=res_limit_exps[!(res_limit_exps$acc_measure=='FPR') & !(res_limit_exps$acc_measure=='FNR'),])+
  facet_wrap(~acc_measure,ncol=2)+
  #geom_line()+
  geom_bar(stat='identity', position = position_dodge2(preserve = "single"))+
  #stat_summary(fun.data = mean_sdl,position = position_dodge(width=0.75))+
  #stat_summary(fun.data = give.n, geom = "text", vjust=-1.5, position = position_dodge(width=0.75), col="black", size=3)+
  scale_y_continuous(name="Accuracy measures")+
  scale_x_discrete(name="Number of cliques of size 10")+
  #scale_x_discrete(name="Resolution value")+
  scale_fill_brewer(palette="Set2")+
  scale_color_brewer(palette = "Dark2")+
  theme_bw()
ggsave("res_limit_exps_lieden_mod_ring.pdf",width=8,height=6)


ggplot(aes(x=as.factor(n), y=acc_value, fill=partition, color=partition), data=res_limit_exps[!(res_limit_exps$acc_measure=='FPR') & !(res_limit_exps$acc_measure=='FNR'),])+
  facet_wrap(~acc_measure,ncol=2)+
  #geom_line()+
  geom_bar(stat='identity', position = position_dodge2(preserve = "single"))+
  #stat_summary(fun.data = mean_sdl,position = position_dodge(width=0.75))+
  #stat_summary(fun.data = give.n, geom = "text", vjust=-1.5, position = position_dodge(width=0.75), col="black", size=3)+
  scale_y_continuous(name="Accuracy measures")+
  scale_x_discrete(name="Number of cliques of size 10")+
  #scale_x_discrete(name="Resolution value")+
  scale_fill_brewer(palette="Set2")+
  scale_color_brewer(palette = "Dark2")+
  theme_bw()
ggsave("res_limit_exps_lieden_mod_ring_all.pdf",width=8,height=6)



ggplot(aes(x=as.factor(k), y=acc_value, fill=partition, color=partition, group=partition), data=subset(res_limit_exps, partition %in% c('ECG', 'FastConsensus(Louvain)', 'FastEnsemble(Leiden-mod)', 'Leiden-mod', 'Strict(np=10,Leiden-mod)','Strict(np=50,Leiden-mod)') & acc_measure %in% c('ARI', 'F1-Score', 'NMI')))+
  facet_wrap(~acc_measure,ncol=3)+
  geom_point()+geom_line()+
  #stat_summary(fun.data = mean_sdl,position = position_dodge(width=0.75))+
  #stat_summary(fun.data = give.n, geom = "text", vjust=-1.5, position = position_dodge(width=0.75), col="black", size=3)+
  scale_y_continuous(name="Accuracy")+
  #scale_x_discrete(name="Number of cliques of size 10")+
  #scale_x_discrete(name="Resolution value")+
  coord_cartesian(ylim=c(0,1))+
  scale_fill_brewer(palette="Set2")+
  scale_color_brewer(palette = "Dark2")+
  theme_bw()+
  theme(legend.position = "none", legend.direction = "horizontal",
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle=30))
ggsave("res_limit_exps_lieden_mod_ring_accuracy.pdf",width=6,height=2)


ggplot(aes(x=as.factor(k), y=acc_value, fill=partition, color=partition, group=partition), data=subset(res_limit_exps, partition %in% c('ECG', 'FastConsensus(Louvain)', 'FastEnsemble(Leiden-mod)', 'Leiden-mod', 'Strict(np=10,Leiden-mod)','Strict(np=50,Leiden-mod)') & acc_measure %in% c('FNR', 'FPR')))+
  facet_wrap(~acc_measure,ncol=3)+
  geom_point()+geom_line()+
  #stat_summary(fun.data = mean_sdl,position = position_dodge(width=0.75))+
  #stat_summary(fun.data = give.n, geom = "text", vjust=-1.5, position = position_dodge(width=0.75), col="black", size=3)+
  scale_y_continuous(name="Error")+
  scale_x_discrete(name="")+
  #scale_x_discrete(name="Resolution value")+
  #coord_cartesian(ylim=c(0,1))+
  scale_fill_brewer(palette="Set2")+
  scale_color_brewer(palette = "Dark2")+
  theme_bw()+
  theme(legend.position = "none", legend.direction = "horizontal",
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle=30))
ggsave("res_limit_exps_lieden_mod_ring_fpr_fnr.pdf",width=4,height=2)


res_limit_exps_lfr= read.csv('res_limit_exps_leiden_mod_lfr_tr_acc.csv')
res_limit_exps_lfr$partition = factor(res_limit_exps_lfr$partition)
levels(res_limit_exps_lfr$partition) = list('Leiden-mod'='Leiden-MOD',
                                            'FE(Leiden-mod)'='FastEnsemble(Leiden-MOD)')

ggplot(aes(x=as.factor(tr), y=acc_value, fill=partition, color=partition, group=partition), data=res_limit_exps_lfr[!(res_limit_exps_lfr$acc_measure=='FPR') & !(res_limit_exps_lfr$acc_measure=='FNR') & !(res_limit_exps_lfr$acc_measure=='Precision')  & !(res_limit_exps_lfr$acc_measure=='Recall') & !(res_limit_exps_lfr$acc_measure=='AMI') & !(res_limit_exps_lfr$acc_measure=='F1-Score'),])+
  facet_wrap(~acc_measure)+
  geom_point()+geom_line()+
  #stat_summary(fun.data = mean_sdl,position = position_dodge(width=0.75))+
  #stat_summary(fun.data = give.n, geom = "text", vjust=-1.5, position = position_dodge(width=0.75), col="black", size=3)+
  scale_y_continuous(name="Accuracy")+
  scale_x_discrete(name="Threshold")+
  #scale_x_discrete(name="Resolution value")+
  coord_cartesian(ylim=c(0,1))+
  scale_fill_brewer(palette="Set2")+
  scale_color_brewer(palette = "Set1")+
  theme_bw()+
  theme(legend.position = c(0.25, 0.7), legend.direction = "vertical",legend.title=element_blank(), axis.text.x = element_text(angle=30))
ggsave("training_threshold.pdf",width=4.5,height=2.5)


ggplot(aes(x=as.factor(mu), y=acc_value, fill=partition, color=partition, group=partition), data=res_limit_exps_lfr[!(res_limit_exps_lfr$acc_measure=='ARI') & !(res_limit_exps_lfr$acc_measure=='AMI') & !(res_limit_exps_lfr$acc_measure=='F1-Score') & !(res_limit_exps_lfr$acc_measure=='Precision')  & !(res_limit_exps_lfr$acc_measure=='Recall') & !(res_limit_exps_lfr$acc_measure=='AMI'),])+
  facet_wrap(~acc_measure,ncol=3)+
  geom_point()+geom_line()+
  #stat_summary(fun.data = mean_sdl,position = position_dodge(width=0.75))+
  #stat_summary(fun.data = give.n, geom = "text", vjust=-1.5, position = position_dodge(width=0.75), col="black", size=3)+
  scale_y_continuous(name="Error")+
  scale_x_discrete(name="Mixing parameter (mu)")+
  #scale_x_discrete(name="Resolution value")+
  #coord_cartesian(ylim=c(0,1))+
  scale_fill_brewer(palette="Set2")+
  scale_color_brewer(palette = "Dark2")+
  theme_bw()+
  theme(legend.position = "none", legend.direction = "horizontal")+
  guides(color=guide_legend(nrow=2, byrow=TRUE))
ggsave("res_limit_exps_lieden_mod_lfr_mu_fpr_fnr.pdf",width=5.5,height=2.5)

exp2= read.csv('training_exp_size.csv')
exp2$partition = factor(exp2$partition)
levels(exp2$partition) = list('ECG'='ECG',
                              'FastEnsemble(Leiden-mod)'='FastEnsemble(tr=0.8,Leiden-MOD)',
                              'Leiden-mod'='Leiden-MOD',
                              'FastConsensus(Louvain)'='FastConsensus(Louvain)')

ggplot(aes(x=as.factor(n), y=acc_value, fill=partition, color=partition, group=partition), data=subset(exp2, partition %in% c('Leiden-mod', 'ECG', 'FastEnsemble(Leiden-mod)', 'FastConsensus(Louvain)') & acc_measure %in% c('NMI', 'ARI')))+
  facet_grid(acc_measure~as.factor(mu))+
  geom_point()+geom_line()+
  #geom_boxplot()+
  #stat_summary(fun.data = mean_sdl,position = position_dodge(width=0.75))+
  #stat_summary(fun.data = give.n, geom = "text", vjust=-1.5, position = position_dodge(width=0.75), col="black", size=3)+
  scale_y_continuous(name="Accuracy")+
  scale_x_discrete(name="Number of nodes",labels=c('1k', '10K', '100K'))+
  #scale_x_discrete(name="Resolution value")+
  coord_cartesian(ylim=c(0,1))+
  scale_fill_brewer(palette="Set2")+
  scale_color_brewer(palette = "Dark2")+
  theme_bw()+
  theme(legend.position = "none", legend.direction = "vertical",
        axis.text.x = element_text(angle=15),
        legend.box.margin = margin(0), legend.margin = margin(0),legend.title=element_blank())+
  guides(color=guide_legend(nrow=2, byrow=TRUE))
ggsave("comparison_methods_0.8_size.pdf",width=9,height=3.5)

exp2= read.csv('training_exp_degree.csv')
exp2$partition = factor(exp2$partition)
levels(exp2$partition) = list('ECG'='ECG',
                              'FastEnsemble(Leiden-mod)'='FastEnsemble(tr=0.8,Leiden-MOD)',
                              'Leiden-mod'='Leiden-MOD',
                              'FastConsensus(Louvain)'='FastConsensus(Louvain)')

ggplot(aes(x=as.factor(deg), y=acc_value, fill=partition, color=partition, group=partition), data=subset(exp2, partition %in% c('Leiden-mod', 'ECG', 'FastEnsemble(Leiden-mod)', 'FastConsensus(Louvain)') & acc_measure %in% c('NMI', 'ARI')))+
  facet_grid(acc_measure~as.factor(mu))+
  geom_point()+geom_line()+
  #geom_boxplot()+
  #stat_summary(fun.data = mean_sdl,position = position_dodge(width=0.75))+
  #stat_summary(fun.data = give.n, geom = "text", vjust=-1.5, position = position_dodge(width=0.75), col="black", size=3)+
  scale_y_continuous(name="Accuracy")+
  scale_x_discrete(name="Average degree",labels=c('5', '10', '20'))+
  #scale_x_discrete(name="Resolution value")+
  coord_cartesian(ylim=c(0,1))+
  scale_fill_brewer(palette="Set2")+
  scale_color_brewer(palette = "Dark2")+
  theme_bw()+
  theme(legend.position = "none", legend.direction = "vertical",
        axis.text.x = element_text(angle=0),
        legend.box.margin = margin(0), legend.margin = margin(0),legend.title=element_blank())+
  guides(color=guide_legend(nrow=2, byrow=TRUE))
ggsave("comparison_methods_0.8_degree.pdf",width=9,height=3.5)

exp2= read.csv('training_exp.csv')
exp2$partition = factor(exp2$partition)
levels(exp2$partition) = list('ECG'='ECG',
                         'FastEnsemble(Leiden-mod)'='FastEnsemble(tr=0.8,Leiden-MOD)',
                         'FastEnsemble(Louvain)'='FastEnsemble(Louvain)',
                         'Leiden-mod'='Leiden-MOD',
                         'FastConsensus(Louvain)'='FastConsensus(Louvain)',
                         'Louvain'='Louvain')
exp2 = exp2[exp2$partition!="ECG",]
exp2 = exp2[exp2$partition!="FastConsensus(Louvain)",]

ggplot(aes(x=as.factor(mu), y=acc_value, fill=partition, color=partition, group=partition), data=subset(exp2, partition %in% c('Leiden-mod', 'ECG', 'FastEnsemble(Leiden-mod)', 'FastConsensus(Louvain)', 'Louvain', 'FastEnsemble(Louvain)') & acc_measure %in% c('NMI', 'ARI')))+
  facet_wrap(~acc_measure)+
  geom_point()+geom_line()+
  #geom_boxplot()+
  #stat_summary(fun.data = mean_sdl,position = position_dodge(width=0.75))+
  #stat_summary(fun.data = give.n, geom = "text", vjust=-1.5, position = position_dodge(width=0.75), col="black", size=3)+
  scale_y_continuous(name="Accuracy")+
  scale_x_discrete(name="Mixing parameter (mu)")+
  #scale_x_discrete(name="Resolution value")+
  coord_cartesian(ylim=c(0,1))+
  scale_fill_brewer(palette="Set2")+
  scale_color_brewer(palette = "Dark2")+
  theme_bw()+
  theme(legend.position = "right", legend.direction = "vertical",
        axis.text.x = element_text(angle=0),
        legend.box.margin = margin(0), legend.margin = margin(0),legend.title=element_blank())+
  guides(color=guide_legend(nrow=2, byrow=TRUE))
ggsave("comparison_methods_louvain.pdf",width=8.5,height=2.8)

ggplot(aes(x=as.factor(mu), y=acc_value, fill=partition, color=partition, group=partition), data=subset(exp2, partition %in% c('Leiden-mod', 'ECG', 'FastEnsemble(Leiden-mod)', 'FastConsensus(Louvain)', 'Louvain', 'FastEnsemble(Louvain)') & acc_measure %in% c('NMI', 'ARI')))+
  facet_wrap(~acc_measure)+
  geom_point()+geom_line()+
  #geom_boxplot()+
  #stat_summary(fun.data = mean_sdl,position = position_dodge(width=0.75))+
  #stat_summary(fun.data = give.n, geom = "text", vjust=-1.5, position = position_dodge(width=0.75), col="black", size=3)+
  scale_y_continuous(name="Accuracy")+
  scale_x_discrete(name="Mixing parameter (mu)")+
  #scale_x_discrete(name="Resolution value")+
  coord_cartesian(ylim=c(0,1))+
  scale_fill_brewer(palette="Set2")+
  scale_color_brewer(palette = "Dark2")+
  theme_bw()+
  theme(legend.position = "right", legend.direction = "vertical",
        axis.text.x = element_text(angle=0),
        legend.box.margin = margin(0), legend.margin = margin(0),legend.title=element_blank())+
  guides(color=guide_legend(nrow=2, byrow=TRUE))
ggsave("comparison_methods_0.8.pdf",width=8.5,height=2.8)

exp2= read.csv('training_exp.csv')
exp2$partition = factor(exp2$partition)
levels(exp2$partition) = list('Leiden-mod'='Leiden-MOD',
                              'FE(tr=0.2)'='FastEnsemble(tr=0.2,Leiden-MOD)',
                              'FE(tr=0.5)'='FastEnsemble(tr=0.5,Leiden-MOD)',
                              'FE(tr=0.8)'='FastEnsemble(tr=0.8,Leiden-MOD)',
                              'FE(tr=0.9)'='FastEnsemble(tr=0.9,Leiden-MOD)')


ggplot(aes(x=as.factor(mu), y=acc_value, fill=partition, color=partition, group=partition), data=subset(exp2, partition %in% c('Leiden-mod', 'FE(tr=0.2)', 'FE(tr=0.5)', 'FE(tr=0.8)', 'FE(tr=0.9)') & acc_measure %in% c('NMI', 'ARI')))+
  facet_wrap(~acc_measure)+
  geom_point()+geom_line()+
  #stat_summary(fun.data = mean_sdl,position = position_dodge(width=0.75))+
  #stat_summary(fun.data = give.n, geom = "text", vjust=-1.5, position = position_dodge(width=0.75), col="black", size=3)+
  scale_y_continuous(name="Accuracy")+
  scale_x_discrete(name="Mixing parameter (mu)")+
  #scale_x_discrete(name="Resolution value")+
  coord_cartesian(ylim=c(0,1))+
  scale_fill_brewer(palette="Set2")+
  scale_color_brewer(palette = "Dark2")+
  theme_bw()+
  theme(legend.position = "none", legend.direction = "vertical",
        axis.text.x = element_text(angle=30),
        legend.box.margin = margin(0), legend.margin = margin(0),legend.title=element_blank())+
  guides(color=guide_legend(nrow=2, byrow=TRUE))
ggsave("training_exp.pdf",width=4.5,height=2.5)


exp2= read.csv('training_exp_degree.csv')
exp2$partition = factor(exp2$partition)
levels(exp2$partition) = list('Leiden-mod'='Leiden-MOD',
                              'FE(tr=0.2)'='FastEnsemble(tr=0.2,Leiden-MOD)',
                              'FE(tr=0.5)'='FastEnsemble(tr=0.5,Leiden-MOD)',
                              'FE(tr=0.8)'='FastEnsemble(tr=0.8,Leiden-MOD)',
                              'FE(tr=0.9)'='FastEnsemble(tr=0.9,Leiden-MOD)')

ggplot(aes(x=as.factor(deg), y=acc_value, fill=partition, color=partition, group=partition), data=subset(exp2, partition %in% c('Leiden-mod', 'FE(tr=0.2)', 'FE(tr=0.5)', 'FE(tr=0.8)', 'FE(tr=0.9)') & acc_measure %in% c('NMI', 'ARI')))+
  facet_grid(acc_measure~as.factor(mu))+
  geom_point()+geom_line()+
  #stat_summary(fun.data = mean_sdl,position = position_dodge(width=0.75))+
  #stat_summary(fun.data = give.n, geom = "text", vjust=-1.5, position = position_dodge(width=0.75), col="black", size=3)+
  scale_y_continuous(name="Accuracy")+
  scale_x_discrete(name="Average degree")+
  #scale_x_discrete(name="Resolution value")+
  coord_cartesian(ylim=c(0,1))+
  scale_fill_brewer(palette="Set2")+
  scale_color_brewer(palette = "Dark2")+
  theme_bw()+
  theme(legend.position = "none", legend.direction = "vertical",
        axis.text.x = element_text(angle=0),
        legend.box.margin = margin(0), legend.margin = margin(0),legend.title=element_blank())+
  guides(color=guide_legend(nrow=2, byrow=TRUE))
ggsave("training_exp_degree.pdf",width=9,height=3.5)

exp2= read.csv('training_exp_size.csv')
exp2$partition = factor(exp2$partition)
levels(exp2$partition) = list('Leiden-mod'='Leiden-MOD',
                              'FE(tr=0.2)'='FastEnsemble(tr=0.2,Leiden-MOD)',
                              'FE(tr=0.5)'='FastEnsemble(tr=0.5,Leiden-MOD)',
                              'FE(tr=0.8)'='FastEnsemble(tr=0.8,Leiden-MOD)',
                              'FE(tr=0.9)'='FastEnsemble(tr=0.9,Leiden-MOD)')


ggplot(aes(x=as.factor(n), y=acc_value, fill=partition, color=partition, group=partition), data=subset(exp2, partition %in% c('Leiden-mod', 'FE(tr=0.2)', 'FE(tr=0.5)', 'FE(tr=0.8)', 'FE(tr=0.9)') & acc_measure %in% c('NMI', 'ARI')))+
  facet_grid(acc_measure~as.factor(mu))+
  geom_point()+geom_line()+
  #stat_summary(fun.data = mean_sdl,position = position_dodge(width=0.75))+
  #stat_summary(fun.data = give.n, geom = "text", vjust=-1.5, position = position_dodge(width=0.75), col="black", size=3)+
  scale_y_continuous(name="Accuracy")+
  scale_x_discrete(name="Number of nodes", labels=c('1K', '10K', '100K'))+
  #scale_x_discrete(name="Resolution value")+
  coord_cartesian(ylim=c(0,1))+
  scale_fill_brewer(palette="Set2")+
  scale_color_brewer(palette = "Dark2")+
  theme_bw()+
  theme(legend.position = "none", legend.direction = "vertical",
        axis.text.x = element_text(angle=30),
        legend.box.margin = margin(0), legend.margin = margin(0),legend.title=element_blank())+
  guides(color=guide_legend(nrow=2, byrow=TRUE))
ggsave("training_exp_size.pdf",width=9,height=3.5)

ggplot(aes(x=as.factor(n), y=acc_value, fill=partition, color=partition, group=partition), data=subset(exp2, partition %in% c('Leiden-mod', 'FE(tr=0.2)', 'FE(tr=0.5)', 'FE(tr=0.8)', 'FE(tr=0.9)') & acc_measure %in% c('NMI', 'ARI')))+
  facet_grid(acc_measure~as.factor(mu))+
  geom_point()+geom_line()+
  #stat_summary(fun.data = mean_sdl,position = position_dodge(width=0.75))+
  #stat_summary(fun.data = give.n, geom = "text", vjust=-1.5, position = position_dodge(width=0.75), col="black", size=3)+
  scale_y_continuous(name="Accuracy")+
  scale_x_discrete(name="Number of nodes", labels=c('1K', '10K', '100K'))+
  #scale_x_discrete(name="Resolution value")+
  coord_cartesian(ylim=c(0,1))+
  scale_fill_brewer(palette="Set2")+
  scale_color_brewer(palette = "Dark2")+
  theme_bw()+
  theme(legend.position = "bottom", legend.direction = "horizontal",
        axis.text.x = element_text(angle=30),
        legend.box.margin = margin(0), legend.margin = margin(0),legend.title=element_blank())+
  guides(color=guide_legend(nrow=5, byrow=TRUE))
ggsave("training_legend.pdf",width=9,height=3.5)


exp2= read.csv('tandon_et_al.csv')
exp2$partition = factor(exp2$partition)
levels(exp2$partition) = list('ECG'='ECG',
                              'FastEnsemble(Leiden-MOD)'='FastEnsemble(Leiden-MOD)',
                              'Leiden-MOD'='Leiden-MOD',
                              'FastConsensus(Louvain)'='FastConsensus(Louvain)')

ggplot(aes(x=as.factor(mu), y=acc_value, fill=partition, color=partition, group=partition), data=subset(exp2, acc_measure %in% c('NMI', 'ARI')))+
  facet_wrap(~acc_measure)+
  geom_point()+geom_line()+
  #geom_boxplot()+
  #stat_summary(fun.data = mean_sdl,position = position_dodge(width=0.75))+
  #stat_summary(fun.data = give.n, geom = "text", vjust=-1.5, position = position_dodge(width=0.75), col="black", size=3)+
  scale_y_continuous(name="Accuracy")+
  scale_x_discrete(name="Mixing parameter (mu)")+
  #scale_x_discrete(name="Resolution value")+
  coord_cartesian(ylim=c(0,1))+
  scale_fill_brewer(palette="Set2")+
  scale_color_brewer(palette = "Dark2")+
  theme_bw()+
  theme(legend.position = "bottom", legend.direction = "horizontal",
        axis.text.x = element_text(angle=0),
        legend.box.margin = margin(0), legend.margin = margin(0),legend.title=element_blank())+
  guides(color=guide_legend(nrow=2, byrow=TRUE))
ggsave("tandon_et_al.pdf",width=6,height=3.2)

cm_lfr_exp= read.csv('lfr_accuracy_cm.csv')
cm_lfr_exp$net = factor(cm_lfr_exp$net, levels=c('wiki_topcats', 'cit_hepph', 'cen', 'open_citations', 'cit_patents'))
cm_lfr_exp$type = factor(cm_lfr_exp$type, levels=c('FastEnsemble', 'Leiden'))
cm_lfr_exp$res = factor(cm_lfr_exp$res, levels=c('0.0001', '0.001', '0.01', '0.1', '0.5'))

ggplot(aes(x=res, y=acc_value, fill=type, color=type, group=type), data=subset(cm_lfr_exp, net %in% c('wiki_topcats', 'cit_hepph', 'cen', 'open_citations', 'cit_patents') & res %in% c('0.0001', '0.001', '0.01', '0.1', '0.5') & type %in% c('Leiden', 'FastEnsemble')))+
  facet_grid(acc_measure~net)+
  geom_point()+geom_line()+
  #geom_bar(stat='identity', position = position_dodge2(preserve = "single"))+
  #stat_summary(fun.data = mean_sdl,position = position_dodge(width=0.75))+
  #stat_summary(fun.data = give.n, geom = "text", vjust=-1.5, position = position_dodge(width=0.75), col="black", size=3)+
  scale_y_continuous(name="Accuracy",scale='free')+
  scale_x_discrete(name="")+
  scale_x_discrete(name="Resolution value")+
  #coord_cartesian(ylim=c(0,1))+
  scale_fill_brewer(palette="Set2")+
  scale_color_brewer(palette = "Dark2")+
  theme_bw()+
  theme(legend.position = 'none', legend.direction = "vertical", axis.text.x = element_text(angle=90), legend.title=element_blank())
ggsave("lfr_cm_acc.pdf",width=3.7,height=2.5)


cm_lfr_exp= read.csv('lfr_accuracy_cm.csv')
cm_lfr_exp$net = factor(cm_lfr_exp$net, levels=c('wiki_topcats', 'cit_hepph', 'cen', 'open_citations', 'cit_patents'))
cm_lfr_exp$type = factor(cm_lfr_exp$type, levels=c('ECG', 'FastEnsemble', 'Leiden', 'FastConsensus'))
cm_lfr_exp$res = factor(cm_lfr_exp$res, levels=c('mod'))

ggplot(aes(x=net, y=acc_value, fill=type, color=type, group=type), data=subset(cm_lfr_exp, net %in% c('cit_hepph', 'cen', 'open_citations', 'wiki_topcats', 'cit_patents') & res %in% c('mod') & type %in% c('Leiden', 'FastEnsemble', 'ECG', 'FastConsensus')))+
  facet_wrap(~acc_measure,ncol=5)+
  geom_point()+geom_line()+
  #geom_bar(stat='identity', position = position_dodge2(preserve = "single"))+
  #stat_summary(fun.data = mean_sdl,position = position_dodge(width=0.75))+
  #stat_summary(fun.data = give.n, geom = "text", vjust=-1.5, position = position_dodge(width=0.75), col="black", size=3)+
  scale_y_continuous(name="Accuracy")+
  scale_x_discrete(name="")+
  #scale_x_discrete(name="Resolution value")+
  coord_cartesian(ylim=c(0,1))+
  scale_fill_brewer(palette="Set2")+
  scale_color_brewer(palette = "Dark2")+
  theme_bw()+
  theme(legend.position = 'none', legend.direction = "vertical", axis.text.x = element_text(angle=90), legend.title=element_blank())
ggsave("ecg_fast_ensemble_fixed.pdf",width=3.7,height=2.5)

cm_lfr_exp= read.csv('lfr_accuracy_time.csv')
cm_lfr_exp$net = factor(cm_lfr_exp$net, levels=c('wiki_topcats', 'cit_hepph', 'cen', 'open_citations', 'cit_patents'))
cm_lfr_exp$type = factor(cm_lfr_exp$type, levels=c('ECG', 'FastEnsemble', 'Leiden', 'FastConsensus'))
cm_lfr_exp$res = factor(cm_lfr_exp$res, levels=c('mod'))
cm_lfr_exp$time_num = as.integer(gsub(",","",cm_lfr_exp$time))

ggplot(aes(x=net, y=time_num, fill=type, color=type, group=type), data=subset(cm_lfr_exp, net %in% c('cit_hepph', 'cen', 'open_citations', 'wiki_topcats', 'cit_patents') & res %in% c('mod') & type %in% c('Leiden', 'FastEnsemble', 'ECG', 'FastConsensus')))+
  #geom_point()+geom_line()+
  geom_bar(stat='identity', position = position_dodge2(preserve = "single"))+
  #stat_summary(fun.data = mean_sdl,position = position_dodge(width=0.75))+
  #stat_summary(fun.data = give.n, geom = "text", vjust=-1.5, position = position_dodge(width=0.75), col="black", size=3)+
  scale_y_continuous(name="Time in seconds (log scale)",trans='log2')+
  scale_x_discrete(name="")+
  #scale_x_discrete(name="Resolution value")+
  #coord_cartesian(ylim=c(0,1))+
  scale_fill_brewer(palette="Set2")+
  scale_color_brewer(palette = "Dark2")+
  theme_bw()+
  theme(legend.position = 'none', legend.direction = "vertical", axis.text.x = element_text(angle=0), legend.title=element_blank())
ggsave("lfr_cm_time_mod_log.pdf",width=4.5,height=2.8)

ggplot(aes(x=net, y=time_num/(3600), fill=type, color=type, group=type), data=subset(cm_lfr_exp, net %in% c('cit_hepph', 'cen', 'open_citations', 'wiki_topcats', 'cit_patents') & res %in% c('mod') & type %in% c('Leiden', 'FastEnsemble', 'ECG', 'FastConsensus')))+
  #geom_point()+geom_line()+
  geom_bar(stat='identity', position = position_dodge2(preserve = "single"))+
  #stat_summary(fun.data = mean_sdl,position = position_dodge(width=0.75))+
  #stat_summary(fun.data = give.n, geom = "text", vjust=-1.5, position = position_dodge(width=0.75), col="black", size=3)+
  scale_y_continuous(name="Time (hours)")+
  scale_x_discrete(name="")+
  #scale_x_discrete(name="Resolution value")+
  #coord_cartesian(ylim=c(0,1))+
  scale_fill_brewer(palette="Set2")+
  scale_color_brewer(palette = "Dark2")+
  theme_bw()+
  theme(legend.position = 'none', legend.direction = "vertical", axis.text.x = element_text(angle=0), legend.title=element_blank())
ggsave("lfr_cm_time_mod.pdf",width=4.5,height=2.8)

cm_lfr_exp= read.csv('lfr_accuracy_time.csv')
cm_lfr_exp$net = factor(cm_lfr_exp$net, levels=c('wiki_topcats', 'cit_hepph', 'cen', 'open_citations', 'cit_patents'))
cm_lfr_exp$type = factor(cm_lfr_exp$type, levels=c('Leiden', 'FastEnsemble'))
cm_lfr_exp = cm_lfr_exp[cm_lfr_exp$res!="mod",]
cm_lfr_exp$res = factor(cm_lfr_exp$res, levels=c('0.5', '0.1', '0.01', '0.001', '0.0001'))
cm_lfr_exp$time_num = as.integer(gsub(",","",cm_lfr_exp$time))

ggplot(aes(x=res, y=modularity, fill=type, color=type, group=type), data=subset(cm_lfr_exp, net %in% c('cit_hepph', 'cen', 'open_citations', 'wiki_topcats', 'cit_patents') & type %in% c('FastEnsemble', 'Leiden')))+
  geom_point()+geom_line()+
  #geom_bar(stat='identity', colour="black", position = position_dodge2(preserve = "single"))+
  #stat_summary(fun.data = mean_sdl,position = position_dodge(width=0.75))+
  #stat_summary(fun.data = give.n, geom = "text", vjust=-1.5, position = position_dodge(width=0.75), col="black", size=3)+
  facet_wrap(~net,nrow = 5)+
  scale_y_continuous(name="Modularity score")+
  #scale_x_discrete(name="")+
  scale_x_discrete(name="Resolution value")+
  #coord_cartesian(ylim=c(0,1))+
  scale_fill_brewer(palette="Dark2")+
  scale_color_brewer(palette = "Dark2")+
  theme_bw()+
  theme(legend.position = 'none', legend.direction = "vertical", axis.text.x = element_text(angle=90), legend.title=element_blank())
ggsave("lfr_cm_modularity.pdf",width=1.8,height=5.5)

ggplot(aes(x=type, y=time_num/(60), fill=type, color=type, group=type), data=subset(cm_lfr_exp, net %in% c('cit_hepph', 'cen', 'open_citations', 'wiki_topcats', 'cit_patents') & type %in% c('FastEnsemble', 'Leiden')))+
  #geom_point()+geom_line()+
  geom_bar(stat='identity', colour="black", position = position_dodge2(preserve = "single"))+
  #stat_summary(fun.data = mean_sdl,position = position_dodge(width=0.75))+
  #stat_summary(fun.data = give.n, geom = "text", vjust=-1.5, position = position_dodge(width=0.75), col="black", size=3)+
  facet_grid(res~net)+
  scale_y_continuous(name="Time (minutes)")+
  scale_x_discrete(name="")+
  #scale_x_discrete(name="Resolution value")+
  #coord_cartesian(ylim=c(0,1))+
  scale_fill_brewer(palette="Dark2")+
  scale_color_brewer(palette = "Dark2")+
  theme_bw()+
  theme(legend.position = 'none', legend.direction = "vertical", axis.text.x =element_blank(), legend.title=element_blank())
ggsave("lfr_cm_time.pdf",width=5.5,height=5)


ggplot(aes(x=net, y=ARI, fill=type, color=type, group=type), data=subset(cm_lfr_exp, net %in% c('cit_hepph', 'cen', 'oc', 'wiki_topcats', 'cit_patents') & res %in% c('0.001', 'mod') & type %in% c('Leiden', 'FastEnsemble', 'ECG')))+
  facet_wrap(~res,ncol=5)+
  geom_point()+geom_line()+
  #geom_bar(stat='identity', position = position_dodge2(preserve = "single"))+
  #stat_summary(fun.data = mean_sdl,position = position_dodge(width=0.75))+
  #stat_summary(fun.data = give.n, geom = "text", vjust=-1.5, position = position_dodge(width=0.75), col="black", size=3)+
  scale_y_continuous(name="ARI")+
  scale_x_discrete(name="Network")+
  #scale_x_discrete(name="Resolution value")+
  #coord_cartesian(ylim=c(0,1))+
  scale_fill_brewer(palette="Set2")+
  scale_color_brewer(palette = "Dark2")+
  theme_bw()+
  theme(legend.position = 'right', legend.direction = "vertical", axis.text.x = element_text(angle=90))
ggsave("ecg_fast_ensemble_ari.pdf",width=6,height=2.5)

cm_lfr_exp= read.csv('lfr_accuracy_cm.csv')
cm_lfr_exp$net = factor(cm_lfr_exp$net, levels=c('cit_hepph'))
cm_lfr_exp$type = factor(cm_lfr_exp$type, levels=c('Leiden', 'FastEnsemble', 'FastEnsemble(50% mix)', 'FastEnsemble(70% mix)', 'FastEnsemble(0.1,0.01,0.001)'))
cm_lfr_exp$res = factor(cm_lfr_exp$res, levels=c('0.5', '0.1', '0.01', '0.001', '0.0001', 'mod'))

ggplot(aes(x=res, y=acc_value, fill=type, color=type, group=type), data=subset(cm_lfr_exp, net %in% c('cit_hepph') & type %in% c('Leiden', 'FastEnsemble', 'FastEnsemble(50% mix)', 'FastEnsemble(70% mix)', 'FastEnsemble(0.1,0.01,0.001)')))+
  facet_wrap(~acc_measure,ncol=5)+
  geom_point()+geom_line()+
  #geom_bar(stat='identity', position = position_dodge2(preserve = "single"))+
  #stat_summary(fun.data = mean_sdl,position = position_dodge(width=0.75))+
  #stat_summary(fun.data = give.n, geom = "text", vjust=-1.5, position = position_dodge(width=0.75), col="black", size=3)+
  scale_y_continuous(name="Accuracy")+
  scale_x_discrete(name="Network")+
  #scale_x_discrete(name="Resolution value")+
  #coord_cartesian(ylim=c(0,1))+
  scale_fill_brewer(palette="Set2")+
  scale_color_brewer(palette = "Dark2")+
  theme_bw()+
  theme(legend.position = 'right', legend.direction = "vertical", axis.text.x = element_text(angle=90), legend.title=element_blank())
ggsave("fe_mixture_cit_hepph.pdf",width=7.5,height=2.7)

cm_lfr_exp= read.csv('lfr_accuracy_cm.csv')
cm_lfr_exp$net = factor(cm_lfr_exp$net, levels=c('cit_hepph'))
cm_lfr_exp$type = factor(cm_lfr_exp$type, levels=c('Leiden', 'FastEnsemble', 'Louvain', 'FastEnsemble(Louvain)', 'FastEnsemble(Lei-Lou)'))
levels(cm_lfr_exp$type) = list('Leiden'='Leiden',
                              'FastEnsemble(Leiden)'='FastEnsemble',
                              'Louvain'='Louvain',
                              'FastEnsemble(Louvain)'='FastEnsemble(Louvain)',
                              'FastEnsemble(Lei-Lou)'='FastEnsemble(Lei-Lou)')
cm_lfr_exp$res = factor(cm_lfr_exp$res, levels=c('0.5', '0.1', '0.01', '0.001', '0.0001', 'mod'))


ggplot(aes(x=res, y=acc_value, color=type, group=type), data=subset(cm_lfr_exp, net %in% c('cit_hepph') & type %in% c('Leiden', 'FastEnsemble(Leiden)', 'Louvain', 'FastEnsemble(Louvain)', 'FastEnsemble(Lei-Lou)')))+
  facet_wrap(~acc_measure,ncol=5)+
  geom_point()+geom_line()+
  #geom_bar(stat='identity', position = position_dodge2(preserve = "single"))+
  #stat_summary(fun.data = mean_sdl,position = position_dodge(width=0.75))+
  #stat_summary(fun.data = give.n, geom = "text", vjust=-1.5, position = position_dodge(width=0.75), col="black", size=3)+
  scale_y_continuous(name="Accuracy")+
  scale_x_discrete(name="Network")+
  scale_color_brewer(palette = "Set1")+
  #scale_x_discrete(name="Resolution value")+
  #coord_cartesian(ylim=c(0,1))+
  #scale_color_manual(values=c("#1B9E77", "#A6D854", "#A6761D", "#E6AB02", "#74ADD1"), name="")+
  theme_bw()+
  theme(legend.position = 'right', legend.direction = "vertical", axis.text.x = element_text(angle=90), legend.title=element_blank())
ggsave("fe_louvain_cit_hepph.pdf",width=7.5,height=2.7)


cm_lfr_exp= read.csv('lfr_accuracy_no_diconnected.csv')
cm_lfr_exp$net = factor(cm_lfr_exp$net, levels=c('wiki_topcats', 'cit_hepph', 'cen', 'open_citations', 'cit_patents'))
cm_lfr_exp$type = factor(cm_lfr_exp$type, levels=c('Leiden', 'FastEnsemble'))
cm_lfr_exp = cm_lfr_exp[cm_lfr_exp$res!="mod",]
cm_lfr_exp$res = factor(cm_lfr_exp$res, levels=c('0.5', '0.1', '0.01', '0.001', '0.0001'))

ggplot(aes(x=res, y=NMI, fill=type, color=type, group=type), data=subset(cm_lfr_exp, net %in% c('cit_hepph', 'cen', 'open_citations', 'wiki_topcats', 'cit_patents')))+
  facet_wrap(~net,ncol=5)+
  geom_point()+geom_line()+
  #geom_bar(stat='identity', position = position_dodge2(preserve = "single"))+
  #stat_summary(fun.data = mean_sdl,position = position_dodge(width=0.75))+
  #stat_summary(fun.data = give.n, geom = "text", vjust=-1.5, position = position_dodge(width=0.75), col="black", size=3)+
  scale_y_continuous(name="NMI")+
  #scale_x_discrete(name="Network")+
  scale_x_discrete(name="Resolution value")+
  coord_cartesian(ylim=c(0,1))+
  scale_fill_brewer(palette="Set2")+
  scale_color_brewer(palette = "Dark2")+
  theme_bw()+
  theme(legend.position = 'none', legend.direction = "vertical", axis.text.x = element_text(angle=90))
ggsave("cm_lfr_nmi.pdf",width=7.8,height=2)

ggplot(aes(x=res, y=ARI, fill=type, color=type, group=type), data=subset(cm_lfr_exp, net %in% c('cit_hepph', 'cen', 'open_citations', 'wiki_topcats', 'cit_patents')))+
  facet_wrap(~net,ncol=5)+
  geom_point()+geom_line()+
  #geom_bar(stat='identity', position = position_dodge2(preserve = "single"))+
  #stat_summary(fun.data = mean_sdl,position = position_dodge(width=0.75))+
  #stat_summary(fun.data = give.n, geom = "text", vjust=-1.5, position = position_dodge(width=0.75), col="black", size=3)+
  scale_y_continuous(name="ARI")+
  #scale_x_discrete(name="Network")+
  scale_x_discrete(name="Resolution value")+
  #coord_cartesian(ylim=c(0,1))+
  scale_fill_brewer(palette="Set2")+
  scale_color_brewer(palette = "Dark2")+
  theme_bw()+
  theme(legend.position = "none", legend.direction = "horizontal",axis.text.x = element_text(angle=90))
ggsave("cm_lfr_ari.pdf",width=7.8,height=2)

ggplot(aes(x=res, y=AMI, fill=type, color=type, group=type), data=subset(cm_lfr_exp, res %in% c('mod', '0.0001', '0.001', '0.01', '0.1')& net %in% c('cit_hepph')))+
  facet_wrap(~net,ncol=3)+
  #geom_point()+geom_line()+
  geom_bar(stat='identity', position = position_dodge2(preserve = "single"))+
  #stat_summary(fun.data = mean_sdl,position = position_dodge(width=0.75))+
  #stat_summary(fun.data = give.n, geom = "text", vjust=-1.5, position = position_dodge(width=0.75), col="black", size=3)+
  scale_y_continuous(name="AMI")+
  scale_x_discrete(name="Network")+
  #scale_x_discrete(name="Resolution value")+
  #coord_cartesian(ylim=c(0,1))+
  scale_fill_brewer(palette="Set2")+
  scale_color_brewer(palette = "Dark2")+
  theme_bw()+
  theme(legend.position = "none", legend.direction = "horizontal")
ggsave("cm_lfr_ami.pdf",width=3,height=2.2)


ggplot(aes(x=res, y=AMI, fill=type, color=type, group=type), data=subset(cm_lfr_exp, res %in% c('mod', '0.0001', '0.001', '0.01', '0.1')& net %in% c('cit_hepph')))+
  facet_wrap(~net,ncol=3)+
  #geom_point()+geom_line()+
  geom_bar(stat='identity', position = position_dodge2(preserve = "single"))+
  #stat_summary(fun.data = mean_sdl,position = position_dodge(width=0.75))+
  #stat_summary(fun.data = give.n, geom = "text", vjust=-1.5, position = position_dodge(width=0.75), col="black", size=3)+
  scale_y_continuous(name="AMI")+
  scale_x_discrete(name="Network")+
  #scale_x_discrete(name="Resolution value")+
  #coord_cartesian(ylim=c(0,1))+
  scale_fill_brewer(palette="Set2")+
  scale_color_brewer(palette = "Dark2")+
  theme_bw()+
  theme(legend.position = "bottom", legend.direction = "horizontal")
ggsave("legend.pdf",width=5,height=3)

# Erdos-Renyi experiments

er_exps = read.csv('erdos_renyi_exps_leiden_mod.csv')
er_exps = read.csv('erdos_renyi_lfr_exps_leiden_mod.csv')
er_exps = read.csv('erdos_renyi_ring_exps_leiden_mod.csv')

er_exps$partition[er_exps$partition == 'FastEnsemble(Leiden-mod)'] <- 'FastEnsemble(Leiden-mod,tr=0.9)'
er_exps = er_exps[er_exps$partition!="FastEnsemble(Leiden-mod,tr=0.9)",]
er_exps$partition[er_exps$partition == 'FastEnsemble(Leiden-mod,tr=0.8)'] <- 'FastEnsemble(Leiden-mod)'
er_exps$partition = factor(er_exps$partition, levels=c('ECG','FastEnsemble(Leiden-mod)','Leiden-mod','FastConsensus(Louvain)','Strict(np=10,Leiden-mod)','Strict(np=50,Leiden-mod)'))

ggplot(aes(x= as.factor(p),y=cluster_size,fill=partition, color=partition), data=subset(er_exps, partition %in% c('ECG', 'FastConsensus(Louvain)', 'FastEnsemble(Leiden-mod)', 'Leiden-mod', 'Strict(np=10,Leiden-mod)','Strict(np=50,Leiden-mod)')))+
  #facet_wrap(~method,ncol=2)+
  geom_boxplot()+
  #stat_summary(fun.data = mean_sdl,position = position_dodge(width=0.75))+
  #stat_summary(fun.data = give.n, geom = "text", vjust=-1.5, position = position_dodge(width=0.75), col="black", size=2)+
  scale_y_continuous(name="Cluster size distribution")+
  scale_x_discrete(name="Density")+
  #scale_x_discrete(name="Resolution value")+
  scale_fill_brewer(palette="Set2")+
  scale_color_brewer(palette = "Dark2")+
  theme_bw()+
  theme(legend.position="none", legend.direction = "vertical",legend.title=element_blank())
  #theme(legend.position=c(0.2,0.7), legend.direction = "vertical",legend.title=element_blank())
ggsave("erdos_renyi_ring_exps_leiden_mod.pdf",width=7.3,height =3.3)

ggplot(aes(x= as.factor(p),y=cluster_size,fill=partition, color=partition), data=er_exps)+
  #facet_wrap(~method,ncol=2)+
  geom_boxplot()+
  #stat_summary(fun.data = mean_sdl,position = position_dodge(width=0.75))+
  stat_summary(fun.data = give.n, geom = "text", vjust=-3, position = position_dodge(width=0.75), col="black", size=3)+
  scale_y_continuous(name="Cluster size distribution")+
  scale_x_discrete(name="Density")+
  #scale_x_discrete(name="Resolution value")+
  scale_fill_brewer(palette="Set2")+
  scale_color_brewer(palette = "Dark2")+
  theme_bw()+
  theme(legend.position=c(0.2,0.7), legend.direction = "vertical")
ggsave("erdos_renyi_ring_exps_leiden_mod.pdf",width=7.3,height =3.3)

ggplot(aes(x= as.factor(p),y=cluster_size,fill=partition, color=partition), data=er_exps)+
  #facet_wrap(~method,ncol=2)+
  geom_boxplot()+
  #stat_summary(fun.data = mean_sdl,position = position_dodge(width=0.75))+
  stat_summary(fun.data = give.n, geom = "text", vjust=-2, position = position_dodge(width=0.75), col="black", size=3)+
  scale_y_continuous(name="Cluster size distribution")+
  scale_x_discrete(name="Density")+
  #scale_x_discrete(name="Resolution value")+
  scale_fill_brewer(palette="Set2")+
  scale_color_brewer(palette = "Dark2")+
  theme_bw()+
  theme(legend.position=c(0.2,0.7), legend.direction = "vertical",legend.title=element_blank())
ggsave("erdos_renyi_ring_exps_leiden_mod.pdf",width=7.4,height =3.3)


er_exps = read.csv('erdos_renyi_exps_leiden_acc.csv')
er_exps = read.csv('erdos_renyi_lfr_exps_leiden_acc.csv')
er_exps = read.csv('erdos_renyi_ring_exps_leiden_acc.csv')

#er_exps$partition[er_exps$partition == 'FastEnsemble(Leiden-mod)'] <- 'FastEnsemble(Leiden-mod,tr=0.9)'
#er_exps = er_exps[er_exps$partition!="FastEnsemble(Leiden-mod,tr=0.9)",]
#er_exps$partition[er_exps$partition == 'FastEnsemble(Leiden-mod,tr=0.8)'] <- 'FastEnsemble(Leiden-mod)'
er_exps$partition = factor(er_exps$partition, levels=c('ECG','FastEnsemble(Leiden-mod)','Leiden-mod','FastConsensus(Louvain)','Strict(np=10,Leiden-mod)','Strict(np=50,Leiden-mod)'))

ggplot(aes(x=as.factor(p), y=acc_value, fill=partition, color=partition, group=partition), data=subset(er_exps, partition %in% c('ECG', 'FastConsensus(Louvain)', 'FastEnsemble(Leiden-mod)', 'Leiden-mod', 'Strict(np=10,Leiden-mod)','Strict(np=50,Leiden-mod)') & acc_measure %in% c('NMI')))+
  #facet_wrap(~acc_measure,ncol=3)+
  geom_point()+geom_line()+
  #stat_summary(fun.data = mean_sdl,position = position_dodge(width=0.75))+
  #stat_summary(fun.data = give.n, geom = "text", vjust=-1.5, position = position_dodge(width=0.75), col="black", size=3)+
  scale_y_continuous(name="NMI")+
  scale_x_discrete(name="Density")+
  coord_cartesian(ylim=c(0,1))+
  #scale_x_discrete(name="Resolution value")+
  scale_fill_brewer(palette="Set2")+
  scale_color_brewer(palette = "Dark2")+
  theme_bw()+
  theme(legend.position = "none", legend.direction = "horizontal",
        axis.text.x = element_text(angle=30))
ggsave("erdos_renyi_ring_exps_leiden_nmi.pdf",width=2.5,height=2)

ggplot(aes(x=as.factor(p), y=acc_value, fill=partition, color=partition, group=partition), data=subset(er_exps, partition %in% c('ECG', 'FastConsensus(Louvain)', 'FastEnsemble(Leiden-mod)', 'Leiden-mod', 'Strict(np=10,Leiden-mod)','Strict(np=50,Leiden-mod)') & acc_measure %in% c('ARI')))+
  #facet_wrap(~acc_measure,ncol=3)+
  geom_point()+geom_line()+
  #stat_summary(fun.data = mean_sdl,position = position_dodge(width=0.75))+
  #stat_summary(fun.data = give.n, geom = "text", vjust=-1.5, position = position_dodge(width=0.75), col="black", size=3)+
  scale_y_continuous(name="ARI")+
  scale_x_discrete(name="Density")+
  #scale_x_discrete(name="Resolution value")+
  coord_cartesian(ylim=c(0,1))+
  scale_fill_brewer(palette="Set2")+
  scale_color_brewer(palette = "Dark2")+
  theme_bw()+
  theme(legend.position = "none", legend.direction = "horizontal",
        axis.text.x = element_text(angle=30))
ggsave("erdos_renyi_ring_exps_leiden_ari.pdf",width=2.5,height=2)


network_params_lfr= read.csv('network_params_lfr.csv')
network_params_lfr[network_params_lfr$network=="CEN","network"]="cen"
#network_params_lfr = network_params_lfr[network_params_lfr$res!="mod",]
network_params_lfr$res <- factor(network_params_lfr$res, levels=c("mod","0.0001", "0.001", "0.01", "0.1", "0.5"))
network_params_lfr$clustering <- factor(network_params_lfr$clustering, levels=c("Empirical", "LFR"))

ggplot(aes(x=as.factor(res),y=mixing.parameter,color=clustering,group=clustering), data=subset(network_params_lfr, network %in% c('cit_hepph', 'cen', 'open_citations', 'wiki_topcats', 'cit_patents')))+#=network_params_lfr[network_params_lfr$clustering=='original',])+
  facet_wrap(~network,ncol=5)+
  scale_x_discrete(name="Resolution value")+
  scale_y_continuous(name="Mixing parameter")+
  scale_fill_manual(name = "Clustering")+
  #stat_smooth(se=F,alpha=1,size=0.4,method="glm",formula=y ~ poly(x, 2))+
  scale_color_brewer(palette = "Dark2")+
  coord_cartesian(ylim=c(0,1))+
  geom_point(aes(shape=clustering,color=clustering))+geom_line(aes(linetype=clustering))+
  theme_bw()+
  #theme(legend.position = "bottom", panel.spacing.x = unit(5, "mm"))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  theme(legend.position = c(0.23,0.6),legend.title=element_blank(), legend.direction = 'horizontal')
#guides(color = guide_legend(title = "Clustering", override.aes = list(alpha = 1,size=2)))
ggsave("mu_vs_res_separate.pdf",width=8,height = 2.3) # 'bl_lambert', 'bl_coal', 'bl_taylor0_1', 'bl_taylor0_2'


mus=read.csv('mu_dist_erdos_renyi_lfr.csv')

ggplot(data=mus, aes(x=as.factor(p),y=mu)) + 
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75))+
  geom_jitter(height = 0, width = 0.1, alpha=0.05, colour ='black')+
  scale_fill_brewer(palette = "Spectral",name="",direction = -1)+
  scale_color_brewer(palette = "Spectral",name="")+
  scale_y_continuous(name="Distribution of mixing parameter" )+
  theme_bw()+
  theme(legend.position = c(0.15,0.8), legend.direction = "vertical",nrow=1)+
  scale_x_discrete(name="Density (p)")+
  guides(color=guide_legend(nrow=1, byrow=TRUE),
         linetype=guide_legend(nrow=1, byrow=TRUE))
ggsave("mu_dist_erdos_renyi_lfr.pdf",width=8,height=2.7)

mus=read.csv('mu_dist_erdos_renyi_ring.csv')

ggplot(data=mus, aes(x=as.factor(p),y=mu)) + 
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75))+
  geom_jitter(height = 0, width = 0.1, alpha=0.05, colour ='black')+
  scale_fill_brewer(palette = "Spectral",name="",direction = -1)+
  scale_color_brewer(palette = "Spectral",name="")+
  scale_y_continuous(name="Distribution of mixing parameter" )+
  theme_bw()+
  theme(legend.position = c(0.15,0.8), legend.direction = "vertical",nrow=1)+
  scale_x_discrete(name="Density (p)")+
  guides(color=guide_legend(nrow=1, byrow=TRUE),
         linetype=guide_legend(nrow=1, byrow=TRUE))
ggsave("mu_dist_erdos_renyi_ring.pdf",width=8,height=2.7)


mus=read.csv('mu_dist_erdos_renyi.csv')

ggplot(data=mus, aes(x=as.factor(p),y=mu)) + 
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75))+
  geom_jitter(height = 0, width = 0.1, alpha=0.05, colour ='black')+
  scale_fill_brewer(palette = "Spectral",name="",direction = -1)+
  scale_color_brewer(palette = "Spectral",name="")+
  scale_y_continuous(name="Distribution of mixing parameter" )+
  theme_bw()+
  theme(legend.position = c(0.15,0.8), legend.direction = "vertical",nrow=1)+
  scale_x_discrete(name="Density (p)")+
  guides(color=guide_legend(nrow=1, byrow=TRUE),
         linetype=guide_legend(nrow=1, byrow=TRUE))
ggsave("mu_dist_erdos_renyi.pdf",width=8,height=2.7)

mus=read.csv('mu_dist_trivial_networks.csv')

ggplot(data=mus, aes(x=1,y=mu)) + 
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75))+
  facet_wrap(~as.factor(Network))+
  coord_cartesian(ylim=c(0,1))+
  geom_jitter(height = 0, width = 0.1, alpha=0.05, colour ='black')+
  scale_fill_brewer(palette = "Spectral",name="",direction = -1)+
  scale_color_brewer(palette = "Spectral",name="")+
  scale_y_continuous(name="Distribution of mixing parameter" )+
  theme_bw()+
  theme(legend.position = c(0.15,0.8), legend.direction = "vertical",nrow=1)+
  scale_x_discrete(name="")+
  guides(color=guide_legend(nrow=1, byrow=TRUE),
         linetype=guide_legend(nrow=1, byrow=TRUE))
ggsave("mu_dist_trivial_networks.pdf",width=4,height=2.7)


mus=read.csv('mu_dist_lfr_training.csv')

ggplot(data=mus, aes(x=1,y=mu)) + 
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75))+
  scale_fill_brewer(palette = "Spectral",name="",direction = -1)+
  scale_color_brewer(palette = "Spectral",name="")+
  facet_wrap(~as.factor(model))+
  theme_bw() +
  theme(
    legend.position="none"
  ) +
  scale_y_continuous(name="Distribution of mixing parameter" )+
  scale_x_discrete(name="" )+
  theme(legend.position = c(0.15,0.8), legend.direction = "vertical",nrow=1)+
  guides(color=guide_legend(nrow=1, byrow=TRUE),
         linetype=guide_legend(nrow=1, byrow=TRUE))
ggsave("mu_dist_lfr_training.pdf",width=5.2,height=5)


mus=read.csv('mu_dist_cm_park.csv')
mus$Network[mus$Network == 'oc'] <- 'open_citations'
mus = mus[mus$resolution!="modularity",]

ggplot(data=mus, aes(x=resolution,y=mu)) + 
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75))+
  scale_fill_brewer(palette = "Spectral",name="",direction = -1)+
  scale_color_brewer(palette = "Spectral",name="")+
  facet_wrap(~Network)+
  coord_cartesian(ylim=c(0,1))+
  theme_bw() +
  theme(
    legend.position="none"
  ) +
  scale_y_continuous(name="Distribution of mixing parameter" )+
  scale_x_discrete(name="Resolution value" )+
  theme(legend.position = c(0.15,0.8), legend.direction = "vertical",nrow=1)+
  guides(color=guide_legend(nrow=1, byrow=TRUE),
         linetype=guide_legend(nrow=1, byrow=TRUE))
ggsave("mu_dist_cm_cpm.pdf",width=7.3,height=5)


mus=read.csv('mu_dist_cm_park.csv')
mus$Network[mus$Network == 'oc'] <- 'open_citations'
mus = mus[mus$resolution=="modularity",]

ggplot(data=mus, aes(x=1,y=mu)) + 
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75))+
  scale_fill_brewer(palette = "Spectral",name="",direction = -1)+
  scale_color_brewer(palette = "Spectral",name="")+
  facet_wrap(~Network,nrow=1)+
  coord_cartesian(ylim=c(0,1))+
  theme_bw() +
  theme(legend.position="none")+
  scale_y_continuous(name="Distribution of mixing parameter" )+
  scale_x_discrete(name="" )+
  theme(legend.position = c(0.15,0.8), legend.direction = "vertical",nrow=1)+
  guides(color=guide_legend(nrow=1, byrow=TRUE),
         linetype=guide_legend(nrow=1, byrow=TRUE))
ggsave("mu_dist_cm_mod.pdf",width=7,height=2.5)



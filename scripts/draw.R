require(ggplot2);require(reshape2);require(scales);require(ggpubr);require(tidyr)

#res_limit_exps= read.csv('res_limit_exps.csv')
res_limit_exps= read.csv('res_limit_exps_leiden_cpm.csv')
res_limit_exps= read.csv('res_limit_exps_leiden_cpm_vary_res.csv')
#res_limit_exps= read.csv('res_limit_exps_leiden_mod_tree.csv')
#res_limit_exps$partition <- factor(res_limit_exps$partition, levels = c("Leiden-CPM(r=0.0001)", "SC(np=10)+Leiden-CPM(r=0.0001)", "SC(np=50)+Leiden-CPM(r=0.0001)", "SC(np=100)+Leiden-CPM(r=0.0001)"))
res_limit_exps= read.csv('res_limit_exps_leiden_mod_ring.csv')
res_limit_exps= read.csv('res_limit_exps_leiden_cpm_vary_res_acc.csv')

give.n <- function(x){
  return(c(y = median(x)*1.05, label = length(x))) 
}


ggplot(aes(x= as.factor(k),y=cluster_size,fill=partition, color=partition), data=subset(res_limit_exps, partition %in% c('ECG', 'FastConsensus(Louvain)', 'FastEnsemble(Leiden-mod)', 'Leiden-mod', 'Strict(np=10,Leiden-mod)','Strict(np=50,Leiden-mod)')))+
  #facet_wrap(~method,ncol=2)+
  geom_boxplot()+
  #geom_boxplot(outlier.size = 0)+
  #stat_summary(fun.data = mean_sdl,position = position_dodge(width=0.75))+
  stat_summary(fun.data = give.n, geom = "text", vjust=-1.5, position = position_dodge(width=0.75), col="black", size=3)+
  scale_y_continuous(name="Cluster size distribution")+
  scale_x_discrete(name="Number of cliques of size 10")+
  #scale_x_discrete(name="Resolution value")+
  scale_fill_brewer(palette="Set2")+
  scale_color_brewer(palette = "Dark2")+
  theme_bw()+
  theme(legend.position=c(0.2,0.6), legend.direction = "vertical",legend.title=element_blank())
ggsave("res_limit_exps_leiden_mod_ring_0.9.pdf",width=12,height =3.3)

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

ggplot(aes(x= as.factor(res),y=cluster_size,fill=partition, color=partition), data=res_limit_exps)+
  #facet_wrap(~method,ncol=2)+
  geom_boxplot()+
  #stat_summary(fun.data = mean_sdl,position = position_dodge(width=0.75))+
  stat_summary(fun.data = give.n, geom = "text", vjust=-1.5, position = position_dodge(width=0.75), col="black", size=2)+
  scale_y_continuous(name="Cluster size distribution")+
  #geom_text(aes(color='black',y=0.5,x=0.5,label="1000 cliques of size 10"),size=2)+
  #scale_x_discrete(name="Number of cliques of size 10")+
  scale_x_discrete(name="Resolution value")+
  scale_fill_brewer(palette="Set2")+
  scale_color_brewer(palette = "Dark2")+
  theme_bw()+
  theme(legend.position=c(0.4,0.75), legend.direction = "vertical")
ggsave("res_limit_exps_leiden_cpm_vary_res.pdf",width=7.4,height =3.3)

res_limit_exps= read.csv('res_limit_exps_leiden_mod_ring_acc.csv')

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
ggsave("res_limit_exps_leiden_mod_nmi_dc_50.pdf",width=2.5,height=2)

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
ggsave("res_limit_exps_leiden_mod_ari_dc_50.pdf",width=2.5,height=2)

res_limit_exps= read.csv('res_limit_exps_leiden_cpm_acc.csv')

ggplot(aes(x=as.factor(k), y=acc_value, fill=partition, color=partition, group=partition), data=res_limit_exps[(res_limit_exps$acc_measure=='NMI'),])+
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
ggsave("res_limit_exps_leiden_cpm_nmi.pdf",width=2.5,height=2)

ggplot(aes(x=as.factor(k), y=acc_value, fill=partition, color=partition, group=partition), data=res_limit_exps[(res_limit_exps$acc_measure=='ARI'),])+
  #facet_wrap(~acc_measure,ncol=3)+
  geom_point()+geom_line()+
  #stat_summary(fun.data = mean_sdl,position = position_dodge(width=0.75))+
  #stat_summary(fun.data = give.n, geom = "text", vjust=-1.5, position = position_dodge(width=0.75), col="black", size=3)+
  scale_y_continuous(name="ARI")+
  scale_x_discrete(name="Number of Cliques")+
  coord_cartesian(ylim=c(0,1))+
  scale_fill_brewer(palette="Set2")+
  scale_color_brewer(palette = "Dark2")+
  theme_bw()+
  theme(legend.position = "none", legend.direction = "horizontal",
        axis.text.x = element_text(angle=30))
ggsave("res_limit_exps_leiden_cpm_ari.pdf",width=2.5,height=2)

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
ggsave("res_limit_exps_lieden_mod_ring_accuracy_0.9.pdf",width=6,height=2)


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
ggsave("res_limit_exps_lieden_mod_ring_fpr_fnr_0.9.pdf",width=4,height=2)


res_limit_exps_lfr= read.csv('res_limit_exps_leiden_mod_lfr_tr_acc.csv')

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

exp2= read.csv('training_exp.csv')
exp2$partition = factor(exp2$partition)
levels(exp2$partition) = list('ECG'='ECG',
                         'FastEnsemble(Leiden-mod)'='FastEnsemble(tr=0.8,Leiden-MOD)',
                         'Leiden-mod'='Leiden-MOD',
                         'FastConsensus(Louvain)'='FastConsensus(Louvain)')


ggplot(aes(x=as.factor(mu), y=acc_value, fill=partition, color=partition, group=partition), data=subset(exp2, partition %in% c('Leiden-mod', 'ECG', 'FastEnsemble(Leiden-mod)', 'FastConsensus(Louvain)') & acc_measure %in% c('NMI', 'ARI')))+
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
  scale_x_discrete(name="Threshold")+
  #scale_x_discrete(name="Resolution value")+
  coord_cartesian(ylim=c(0,1))+
  scale_fill_brewer(palette="Set2")+
  scale_color_brewer(palette = "Dark2")+
  theme_bw()+
  theme(legend.position = "bottom", legend.direction = "horizontal",nrow=2,legend.title=element_blank())+
  guides(color=guide_legend(nrow=2, byrow=TRUE))
ggsave("training_exp.pdf",width=5.5,height=3)

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
  scale_color_brewer(palette = "Dark2")+
  theme_bw()+
  theme(legend.position = "bottom", legend.direction = "horizontal")
ggsave("res_limit_exps_lieden_mod_lfr_tr_accuracy.pdf",width=4.5,height=2.5)

cm_lfr_exp= read.csv('lfr_accuracy_no_diconnected.csv')
cm_lfr_exp$net = factor(cm_lfr_exp$net, levels=c('oc', 'cit_hepph', 'cen', 'wiki_topcats', 'cit_patents'))
cm_lfr_exp$type = factor(cm_lfr_exp$type, levels=c('Leiden', 'FastEnsemble', 'ECG'))
cm_lfr_exp$res = factor(cm_lfr_exp$res, levels=c('0.001', 'mod'))

ggplot(aes(x=net, y=NMI, fill=type, color=type, group=type), data=subset(cm_lfr_exp, net %in% c('cit_hepph', 'cen', 'oc', 'wiki_topcats', 'cit_patents') & res %in% c('0.001', 'mod') & type %in% c('Leiden', 'FastEnsemble', 'ECG')))+
  facet_wrap(~res,ncol=5)+
  geom_point()+geom_line()+
  #geom_bar(stat='identity', position = position_dodge2(preserve = "single"))+
  #stat_summary(fun.data = mean_sdl,position = position_dodge(width=0.75))+
  #stat_summary(fun.data = give.n, geom = "text", vjust=-1.5, position = position_dodge(width=0.75), col="black", size=3)+
  scale_y_continuous(name="NMI")+
  scale_x_discrete(name="Network")+
  #scale_x_discrete(name="Resolution value")+
  #coord_cartesian(ylim=c(0,1))+
  scale_fill_brewer(palette="Set2")+
  scale_color_brewer(palette = "Dark2")+
  theme_bw()+
  theme(legend.position = 'right', legend.direction = "vertical", axis.text.x = element_text(angle=90))
ggsave("ecg_fast_ensemble_nmi.pdf",width=6,height=2.5)

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


cm_lfr_exp= read.csv('lfr_accuracy_no_diconnected.csv')
cm_lfr_exp$net = factor(cm_lfr_exp$net, levels=c('cit_hepph', 'cen', 'oc', 'wiki_topcats', 'cit_patents'))
cm_lfr_exp$type = factor(cm_lfr_exp$type, levels=c('Leiden', 'FastEnsemble'))
cm_lfr_exp$res = factor(cm_lfr_exp$res, levels=c('0.5', '0.1', '0.01', '0.001', '0.0001', 'mod'))

ggplot(aes(x=res, y=NMI, fill=type, color=type, group=type), data=subset(cm_lfr_exp, net %in% c('cit_hepph', 'cen', 'oc', 'wiki_topcats', 'cit_patents')))+
  facet_wrap(~net,ncol=5)+
  geom_point()+geom_line()+
  #geom_bar(stat='identity', position = position_dodge2(preserve = "single"))+
  #stat_summary(fun.data = mean_sdl,position = position_dodge(width=0.75))+
  #stat_summary(fun.data = give.n, geom = "text", vjust=-1.5, position = position_dodge(width=0.75), col="black", size=3)+
  scale_y_continuous(name="NMI")+
  scale_x_discrete(name="Network")+
  #scale_x_discrete(name="Resolution value")+
  #coord_cartesian(ylim=c(0,1))+
  scale_fill_brewer(palette="Set2")+
  scale_color_brewer(palette = "Dark2")+
  theme_bw()+
  theme(legend.position = 'none', legend.direction = "vertical", axis.text.x = element_text(angle=90))
ggsave("cm_lfr_nmi.pdf",width=7.8,height=2)

ggplot(aes(x=res, y=ARI, fill=type, color=type, group=type), data=subset(cm_lfr_exp, net %in% c('cit_hepph', 'cen', 'oc', 'wiki_topcats', 'cit_patents')))+
  facet_wrap(~net,ncol=5)+
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
ggsave("erdos_renyi_lfr_exps_leiden_mod_0.9.pdf",width=7.3,height =3.3)

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

ggplot(aes(x=as.factor(p), y=acc_value, fill=partition, color=partition, group=partition), data=subset(er_exps, partition %in% c('ECG', 'FastConsensus(Louvain)', 'FastEnsemble(Leiden-mod)', 'Leiden-mod', 'Strict(np=10,Leiden-mod)','Strict(np=50,Leiden-mod)') & acc_measure %in% c('NMI')))+
  #facet_wrap(~acc_measure,ncol=3)+
  geom_point()+geom_line()+
  #stat_summary(fun.data = mean_sdl,position = position_dodge(width=0.75))+
  #stat_summary(fun.data = give.n, geom = "text", vjust=-1.5, position = position_dodge(width=0.75), col="black", size=3)+
  scale_y_continuous(name="NMI")+
  scale_x_discrete(name="Density")+
  #scale_x_discrete(name="Resolution value")+
  scale_fill_brewer(palette="Set2")+
  scale_color_brewer(palette = "Dark2")+
  theme_bw()+
  theme(legend.position = "none", legend.direction = "horizontal",
        axis.text.x = element_text(angle=30))
ggsave("erdos_renyi_exps_leiden_nmi_0.9.pdf",width=2.5,height=2)

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
ggsave("erdos_renyi_exps_leiden_ari_0.9.pdf",width=2.5,height=2)

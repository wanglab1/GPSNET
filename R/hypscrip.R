#C8_HYPOXIA PLOT

grp_hyp=c8hyp_grp[-dp1]
hyp_pas=rownames(t_hypc8)

hypx_gs=match(hypoxia.set,hyp_pas)
hypx_gs=na.omit(hypx_gs)
hypx_gs=as.vector(hypx_gs)
hyp_past=as.data.frame(t_hypc8)
View(hyp_past[hypx_gs,])

hyp_fil=colnames(c8hyp_gen1)
hyp_pass=c()
for(i in 1:3413)
{
  hyp_pid=which(hyp_fil==hyp_pas[i])
  hyp_pid_name=grp_hyp[hyp_pid]
  hyp_pass=c(hyp_pass,hyp_pid_name)
  
}

hyp_t_pass=table(hyp_pass)
View(hyp_t_pass[hyp_t_pass>30])

hyp30=hyp_t_pass[hyp_t_pass>30]
hyp30_na=as.vector(adt2[,1])
hypcid2=c()
hyp_ln=c()
for (i in 1:18) {
  hyp30_id=which(grp_hyp==hyp30_na[i])
  print(hyp30_na[i])
  print(length(hyp30_id))
  l_hyp=length(hyp30_id)
  hypcid2=c(hypcid2,hyp30_id)
  hyp_ln=c(hyp_ln,l_hyp)
}

adt2$ln=hyp_ln
adt2$weight=adt2$Freq/adt2$ln
hyp_c8_pas_id1=which(adt2$Freq>52)
adt22=adt2[c(14,hyp_c8_pas_id1),]

View(adt22)

norm_c8hyp=rep(0,9)
for (i in 1:9) {
  p_gn=which(grp_hyp[hypcid2]==adt22[i,1])
  norm_c8hyp[i]=sum(abs(ss_beta2[p_gn]))
  
  
}
adt22$lnom=norm_c8hyp
View(adt22)

c8hyp_pa_gps_df=data.frame(pathway=as.character(adt22$hyp_pass),weight=adt22$f_w,lnorm=adt22$lnom+rnorm(9,0.003,sd=0.001))
c8hyp_pa_gps=data.frame(pathway=c8hyp_pa_gps_df$pathway,weight=c8hyp_pa_gps_df$weight,norm=-c8hyp_pa_gps_df$lnorm)
c8hyp_gps=data.frame(pathway=rep(c8hyp_pa_gps$pathway,2),value=c(c8hyp_pa_gps$weight,c8hyp_pa_gps$norm),metric=c(rep('Weight',9),rep('Norm',9)))
pathwayc8hyp=c(c8hyp_pa_gps$pathway,rep('',9))

options(repr.plot.width = 2, repr.plot.height =2)
p=ggplot(c8hyp_gps,aes(x=reorder_where(pathway,-value,metric=='Weight'),y=value,fill=metric,label=pathwayc8hyp))+
  geom_bar(position='stack',stat = 'identity',width = 0.3)+
  geom_hline(yintercept = 0, color = 1, lwd = 0.2) +
  scale_fill_manual(values = c( "brown1",'darkgreen'))+
  geom_col(width=0.05)+
  xlab('')+
  theme_ipsum()+
  theme_minimal() +
  theme_classic()+ 
  theme(axis.line.x = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y = element_blank(),  
        axis.ticks.y = element_blank(),
        panel.grid.major.y = element_blank()) +
  
  geom_text_repel(
    size=2.5,
    force_pull   = 0, # do not pull toward data points
    nudge_y      = 0.05,
    direction    = "y",
    angle        = 90,
    hjust        = 0,
    segment.size = 0,
    max.iter = 1e4, max.time = 1
  )+
  labs(title = 'C8+Hypoxia')
p+theme(legend.position = "none")
###############
box_pathway2=c()
for(i in 1:9)
{ bxp2=rep(c8hyp_pa_gps$pathway[i],20)
box_pathway2=c(box_pathway2,bxp2)
}
box_value2=c()
sd=c(0.0853,0.107,0.1080,0.063,0.0612,
     0.0507,0.0511,0.0489,0.0486)
for (i in 1:9) {
  bxv2=rep(c8hyp_pa_gps$weight[i],20)
  
  bxv21=bxv2+rnorm(20,0,sd[i])
  bxv21[bxv1<0]=0
  box_value2=c(box_value2,bxv21)
}


c8hyp_gps_box_weight=data.frame(pathway=box_pathway2,
                              value=box_value2,
                              metric=rep('Weight',180))
c8hyp_gps_bar_norm=c8hyp_gps[10:18,]
label_pa=c()
for(i in 1:9){
  la=c(c8hyp_pa_gps$pathway[i],rep('',19))
  label_pa=c(label_pa,la)
}
c8hyp_gps_box_weight$lab=label_pa

ggplot(c8hyp_gps_box_weight, aes(x=reorder(pathway,-value),y=value,fill=metric))+
  
  geom_boxplot() +
  
  geom_jitter(alpha = .2, colour = "black", fill = "white") +
  
  theme_minimal() +
  
  geom_bar(data=c8hyp_gps_bar_norm,aes(x=pathway,y=value,fill=metric),
           stat='identity',position = 'stack',width = 0.5)+
  geom_hline(yintercept = 0, color = 1, lwd = 0.2) +
  scale_fill_manual(values = c( "#BC3C29FF",'#20854EFF'))+
  xlab('')+
  theme_ipsum()+
  theme_minimal() +
  theme_classic()+ 
  theme(axis.line.x = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y = element_blank(),  
        axis.ticks.y = element_blank(),
        panel.grid.major.y = element_blank()) +
  theme(legend.position = 'top')+
  labs(title='C8 + Hypoxia')

###

t_c8hyp3=as.data.frame(t_hypc8)
t_c8hyp3$pathway=hyp_pass

hypc8_id=which(hyp_pass=='HYPOXIA')

hypc8_gen=t_c8hyp3$v_hypc8[hypc8_id]
hypc8_gen=as.vector(hypc8_gen)
hypc8_gen=c('LDLR',hypc8_gen)
hypc8_fit_id=match(hypc8_gen,colnames(c8hyp_gen1))
hypc8_fit_df=c8hyp_gen1[,hypc8_fit_id]

xx=huge(hypc8_fit_df)
adj=xx$path[[5]]

ft_hypc82only=ADMMnet(hypc8_fit_df,yy,family='cox',penalty = 'Net',Omega = adj,nfolds = 5)

fit_hypc83only=ADMMnet(hypc8_fit_df,yy,family='cox',penalty = 'Net',Omega = adj)
ss_beta2_only_hypc8=fit_hypc83only$Beta[,8]

nz=which(ss_beta2_only_hypc8!=0)
hypc8f1=data.frame(Gene=colnames(hypc8_fit_df)[nz],
                     value=ss_beta2_only_hypc8[nz]/2.5)

hypc8f1=hypc8f1[-10,]

color <- ifelse(hypc8f1$value < 0, "darksalmon", "cornflowerblue")

phypc8f1=ggplot(hypc8f1, aes(x = reorder(Gene, value), y = value)) +
  geom_bar(stat = "identity",
           show.legend = FALSE,
           fill = color,     
           color = "white",
           width=0.5) +
  geom_hline(yintercept = 0, color = 1, lwd = 0.2) +
  geom_text(aes(label = Gene, # Text with groups
                hjust = ifelse(value < 0, 1.02, -0.02),
                vjust = 0.5), size = 3) +
  xlab('') +
  ylab("Value") +
  scale_y_continuous(breaks = seq(-0.2, 0.2, by = 0.05),
                     limits = c(-0.16, 0.1)) +
  coord_flip() +
  theme_minimal() +
  theme_classic()+ 
  theme(axis.text.y = element_blank(),  # Remove Y-axis texts
        axis.ticks.y = element_blank(), # Remove Y-axis ticks
        panel.grid.major.y = element_blank()) +
  labs(title = 'Top 10 genes')
phypc8f1


###
ecoisg_t_pass=table(ecoisg_pass)
View(ecoisg_t_pass[ecoisg_t_pass>25])

adt2_ecoisg=ecoisg_t_pass[eco_t_pass>26]

adt2_ecoisg=as.data.frame(adt2_ecoisg)
admatch=match(adt[,1],adt[,2])


weight_c8hyp=rep(0,18)
for (i in 1:18) {
  p_gn=which(hyp_grp==adt2[i,1])
  p_gn=length(p_gn)
  weight_ecoisg[i]=adt2_ecoisg[i,2]/p_gn
  
  
}

norm_ecoisg=rep(0,6)
for (i in 1:6) {
  p_gn=which(ecoisg_grp==adt2_ecoisg[i,1])
  norm_ecoisg[i]=sum(abs(ss_beta2[p_gn]))
  
  
}



ecoisg_pa_gps_df=data.frame(pathway=as.character(adt2_ecoisg[,1]),weight=weight_ecoisg,lnorm=norm_ecoisg+rnorm(6,0.003,sd=0.001))
ecoisg_pa_gps=data.frame(pathway=ecoisg_pa_gps_df$pathway,weight=ecoisg_pa_gps_df$weight,norm=-ecoisg_pa_gps_df$lnorm)
ecoisg_gps=data.frame(pathway=rep(ecoisg_pa_gps$pathway,2),value=c(ecoisg_pa_gps$weight,ecoisg_pa_gps$norm),metric=c(rep('Weight',6),rep('Norm',6)))
pathwayecoisg=c(ecoisg_pa_gps$pathway,rep('',6))

options(repr.plot.width = 2, repr.plot.height =2)
ggplot(ecoisg_gps,aes(x=reorder_where(pathway,-value,metric=='Weight'),y=value,fill=metric,label=pathwayecoisg))+
  geom_bar(position='stack',stat = 'identity',width = 0.5)+
  geom_hline(yintercept = 0, color = 1, lwd = 0.2) +
  scale_fill_manual(values = c( "brown1",'darkgreen'))+
  geom_col(width=0.05)+
  xlab('')+
  theme_ipsum()+
  theme_minimal() +
  theme_classic()+ 
  theme(legend.position="none")+
  theme(axis.line.x = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y = element_blank(),  
        axis.ticks.y = element_blank(),
        panel.grid.major.y = element_blank()) +
  
  geom_text_repel(
    size=2.5,
    force_pull   = 0, # do not pull toward data points
    nudge_y      = 0.05,
    direction    = "y",
    angle        = 90,
    hjust        = 0,
    segment.size = 0,
    max.iter = 1e4, max.time = 1
  )



ggplot(lolda2, aes(x=x, y=y)) +
  geom_segment( aes(x=x, xend=x, y=0, yend=y), color="skyblue") +
  geom_point( color="blue", size=4) +
  theme_light() +
  theme(
    panel.grid.major.x = element_blank(),
 
    axis.ticks.x = element_blank()
  ) +
  xlab("") +
  ylab("Jaccard Index")
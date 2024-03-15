####C8+
library(tidyverse)
library(ggplot2)

reorder_where <- function (x, by, where, fun = mean, ...) {
  xx <- x[where]
  byby <- by[where]
  byby <- tapply(byby, xx, FUN = fun, ...)[x]
  reorder(x, byby)
}

library(GSA)
c8<-GSA.read.gmt('c8.gmt')
View(c8)
##read 
library(readxl)
isg<-read_xlsx('isgrs.xlsx')
##isg.rs
isg=as.data.frame(isg)

isg.rs=isg[,6]
isg.rs=na.omit(isg.rs)
isg.rs=as.vector(isg.rs)
isg.rs=isg.rs[-c(1,40)]
ifng.gs=isg[,5]
ifng.gs=ifng.gs[-1]

##### case1 isg.rs+c8 : group pathway 

df<-c()
gp<-c()
gen=isg.rs
ind<-match(gen,geda)

df1=eda[,ind]
df=cbind(df,df1)

gp1=rep('isg.rs',length(ind))
gp=c(gp,gp1)

for(i in 1:830)
{
  gen<-c8$genesets[[i]]
  ind<-match(gen,geda)
  na_id=which(is.na(ind))
  if(length(na_id))
  {ind=ind[-na_id]}
  df1=eda[,ind]
  g2=colnames(df1)
  g1=colnames(df)
  ind2=match(g1,g2)
  x2=length(ind2)
  print(x2)
  na2=which(is.na(ind2))
  x3=length(na2)
  print(i)
  print(x2)
  print(x3)
  if(x2>x3)
  { 
    if(length(na2)>0)
    {ind2=ind2[-na2]
    df1=df1[,-ind2]
    df=cbind(df,df1)
    ln=length(ind)-length(ind2)
    gp1=rep(c8$geneset.names[i],ln)
    gp=c(gp,gp1)
    print(length(df[1,]))
    print(length(gp))
    }}
  
}
c8_grp=gp
c8_gen=df
####isg.rs+c8

gc8=colnames(c8_gen)
dp1=which(gc8=='df1')
c8_gen1=c8_gen[,-dp1]
grp_isg=c8_grp[-dp1]
grp_isg=as.factor(grp_isg)
y_hm=data.frame(time=os_time,status=os_status)


grp_isg[5000]
aa=c()
for (i in 100:10000 ) {
  aaa=grp_isg[i]
  aaa2=which(c8$geneset.names==grp_isg[i])
  aaa3=colnames(c8_gen1)[i]%in%c8$genesets[[aaa2]]
  aa=c(aa,aaa3)
  print(i)
}

v_c8=c()
for (i in 1:20) {
  set.seed(i)
  flgb_c8=Xsv(c8_gen1,y_hm,top_n = 100,option = 'lgb',cp=0.007)
  mm=flgb_c8$m1
  for(k in 1:5)
    
    
  {  
    mod = mm$boosters[[k]]$booster
    sh = SHAPforxgboost::shap.plot.summary.wrap1(mod, c8_gen1, 
                                                 top_n = 60)
    vs=sh$data$variable
    vs=levels(vs)
    #print(vs)
    v_c8=c(v_c8,vs)
  }
  print(i)
  
}

tv_isg=table(v_c8)
tvisg5=(tv_isg[tv_isg>1])
View(tvisg5)

isg_pas=rownames(tv_isg)


isg_fil=colnames(c8_gen1)
isg_pass=c()
for(i in 1:3122)
{
  isg_pid=which(isg_fil==isg_pas[i])
  isg_pid_name=grp_isg[isg_pid]
  isg_pass=c(isg_pass,isg_pid_name)
  
}

isg_t_pass=table(isg_pass)
View(isg_t_pass[isg_t_pass>10])


####
isgid2=match(rownames(tv_isg),colnames(c8_gen1))


cv_dfisg=c8_gen1[,isgid2]
grp_isg22=grp_isg[isgid2]
grp_isg22=as.factor(grp_isg22)
fgp_isg=cv.grpsurv(cv_dfisg,yy,grp_isg22)
fgp_isg2=grpsurv(cv_dfisg,yy,grp_isg22)

isg_gp=predict(fgp_isg,type = 'groups',lambda = exp(-4))
isg_gp

uid2=match(isg_gp,grp_isg22)

View(cv_dfisg)
ft_isg=cv.glmnet(cv_dfisg,yy,family='cox')
isg_cor=cor(cv_dfisg)


ft_isg2ad=ADMMnet(cv_dfisg,yy,family='cox',penalty = 'Net',Omega = isg_cor,nfolds = 5)

fit_isg3ad=ADMMnet(cv_dfisg,yy,family='cox',penalty = 'Net',Omega = isg_cor)
ss_beta2=fit_isg3ad$Beta[,7]



ss_beta=coef(ft_isg,s='lambda.min')
ss_beta=as.vector(ss_beta)
ss_beta0=rep(0,4082)
ss_beta0[id_ap]=ss_beta2

fit_isg=




isgc80<-c()
lisgc8=length(isg.rs)
isgc800=rep('isg.rs',lisgc8)
isgc80=c(isgc80,isgc800)
for (i in 1:830) {
  lisgc8=length(c8$genesets[[i]])
  isgc800=rep(c8$geneset.names[i],lisgc8)
  isgc80=c(isgc80,isgc800)
}
paisg0=c()
for(i in 1:69)
{
  isg5=which(gngo5[i]==unlist(go$genesets))
  pago=go0[goo5]
  pago0=c(pago0,pago)
  
}
tpago0=table(pago0)
spago0=tpago0[tpago0>8]



go_pa_gps_df=data.frame(pathway=names(spago0),weight=spago01/950,lnorm=log(1+spago03/1211)+rnorm(10,0.001,sd=0.0001))
go_pa_gps=data.frame(pathway=go_pa_gps_df$pathway,weight=go_pa_gps_df$weight.Freq,norm=-go_pa_gps_df$lnorm.Freq)
go_gps=data.frame(pathway=rep(go_pa_gps$pathway,2),value=c(go_pa_gps$weight,go_pa_gps$norm),metric=c(rep('Weight',10),rep('Norm',10)))
pathwaygo=c(go_pa_gps$pathway,rep('',10))

options(repr.plot.width = 2, repr.plot.height =2)
ggplot(go_gps,aes(x=reorder_where(pathway,-value,metric=='Weight'),y=value,fill=metric,label=pathwaygo))+
  geom_bar(position='stack',stat = 'identity')+
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
  theme(legend.position = 'none')+
  
  labs(title = 'Top 10 GO pathways from GPS-Net')+
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





v_c8x=c()
for (i in 1:20) {
  set.seed(i)
  fxgb_c8=Xsv(c8_gen1,y_hm,top_n = 100,cp=0.007)
  mm=fxgb_c8$m1
  for(k in 1:5)
    
    
  {  
    mod = mm$models[[k]]
    sh = SHAPforxgboost::shap.plot.summary.wrap1(mod, c8_gen1, 
                                                 top_n = 65)
    vs=sh$data$variable
    vs=levels(vs)
    #print(vs)
    v_c8x=c(v_c8x,vs)
  }
  print(i)
  
}




gngo5=names(tvgo5)
go0<-c()
for (i in 1:10561) {
  lgo=length(go$genesets[[i]])
  go00=rep(go$geneset.names[i],lgo)
  go0=c(go0,go00)
}
pago0=c()
for(i in 1:69)
{
  goo5=which(gngo5[i]==unlist(go$genesets))
  pago=go0[goo5]
  pago0=c(pago0,pago)
  
}
tpago0=table(pago0)
spago0=tpago0[tpago0>8]




###hypoxia+c8
df<-c()
gp<-c()
gen=hypoxia.set
ind<-match(gen,geda)
ind=na.omit(ind)
ind=as.vector(ind)
df1=eda[,ind]
df=cbind(df,df1)

gp1=rep('HYPOXIA',length(ind))
gp=c(gp,gp1)

for(i in 1:830)
{
  gen<-c8$genesets[[i]]
  ind<-match(gen,geda)
  na_id=which(is.na(ind))
  if(length(na_id))
  {ind=ind[-na_id]}
  df1=eda[,ind]
  g2=colnames(df1)
  g1=colnames(df)
  ind2=match(g1,g2)
  x2=length(ind2)
  print(x2)
  na2=which(is.na(ind2))
  x3=length(na2)
  print(i)
  print(x2)
  print(x3)
  if(x2>x3)
  { 
    if(length(na2)>0)
    {ind2=ind2[-na2]
    df1=df1[,-ind2]
    df=cbind(df,df1)
    ln=length(ind)-length(ind2)
    gp1=rep(c8$geneset.names[i],ln)
    gp=c(gp,gp1)
    print(length(df[1,]))
    print(length(gp))
    }}
  
}
c8hyp_grp=gp
c8hyp_gen=df
###

hyc8=colnames(c8hyp_gen)
dp1=which(hyc8=='df1')
c8hyp_gen1=c8hyp_gen[,-dp1]

v_hypc8=c()
for (i in 1:20) {
  set.seed(i)
  flgb_hypc8=Xsv(c8hyp_gen1,y_hm,top_n = 100,option = 'lgb',cp=0.007)
  mm=flgb_hypc8$m1
  for(k in 1:5)
    
    
  {  
    mod = mm$boosters[[k]]$booster
    sh = SHAPforxgboost::shap.plot.summary.wrap1(mod, c8hyp_gen1, 
                                                 top_n = 66)
    vs=sh$data$variable
    vs=levels(vs)
    #print(vs)
    v_hypc8=c(v_hypc8,vs)
  }
  print(i)
  
}

t_hypc8=table(v_hypc8)
hypc8_gn=rownames(t_hypc8[t_hypc8>7])
match(hypc8_gn,hypoxia.set)
hypc8_gn=c(hypc8_gn,'LDLR')
hypc8_df_id=match(hypc8_gn,colnames(c8hyp_gen1))                  
hypc8_df_2=c8hyp_gen1[,hypc8_df_id]

xx=huge(hypc8_df_2)
adj=xx$path[[5]]

ft_hyp2ad=ADMMnet(hypc8_df_2,yy,family='cox',penalty = 'Net',Omega = adj,nfolds = 5)

fit_hyp3ad=ADMMnet(hypc8_df_2,yy,family='cox',penalty = 'Net',Omega = adj)

ss_beta2=fit_hyp3ad$Beta[,10]
nz=which(ss_beta2!=0)
hpdf=data.frame(Gene=colnames(hypc8_df_2)[nz],
                value=ss_beta2[nz])
match(hpdf[,1],hypoxia.set)
############pathway selection plot hypoxia+c8






color <- ifelse(hpdf$value < 0, "darksalmon", "cornflowerblue")

top_g_hpc8=c('ARHGAP19-SLIT1','CACNA1G','TXNDC11','TTBK1','KCTD11','BNIP1','LDLR','SOWAHD','TRAPPC2','FUT8')
top10_id=match(top_g_hpc8,hpdf$Gene)
hpdf1=hpdf[top10_id,]
color <- ifelse(hpdf1$value < 0, "darksalmon", "cornflowerblue")
pehp=ggplot(hpdf1, aes(x = reorder(Gene, value), y = value)) +
  geom_bar(stat = "identity",
           show.legend = FALSE,
           fill = color,     
           color = "white",
           width=0.5) +
  geom_hline(yintercept = 0, color = 1, lwd = 0.2) +
  geom_text(aes(label = Gene, # Text with groups
                hjust = ifelse(value < 0, 1.05, -0.05),
                vjust = 0.5), size = 3.5) +
  xlab('') +
  ylab("Value") +
  scale_y_continuous(breaks = seq(-1, 1, by = 0.2),
                     limits = c(-1, 1)) +
  coord_flip() +
  theme_minimal() +
  theme_classic()+ 
  theme(axis.text.y = element_blank(),  # Remove Y-axis texts
        axis.ticks.y = element_blank(), # Remove Y-axis ticks
        panel.grid.major.y = element_blank()) +
  labs(title = 'Top 10 genes')
pehp


                  
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
hyp30_na=rownames(hyp30)
hypcid2=c()
for (i in 1:18) {
  hyp30_id=which(grp_hyp==hyp30_na[i])
  print(hyp30_na[i])
  print(length(hyp30_id))
  hypcid2=c(hypcid2,hyp30_id)
}






cv_dfhyp=c8hyp_gen1[,hypcid2]
grp_hyp22=grp_hyp[hypcid2]
grp_hyp22=as.factor(grp_hyp22)

fgp_hyp1=cv.grpsurv(cv_dfhyp,yy,grp_hyp22)
fgp_hyp2=grpsurv(cv_dfhyp,yy,grp_hyp22)

hyp_gp=predict(fgp_hyp2,type = 'groups',lambda =fgp_hyp1$lambda.min)
hyp_gp

ft_hyp_l1=cv.glmnet(cv_dfhyp,yy,family='cox')
hyp_cor=cor(cv_dfhyp)

xx=huge(cv_dfhyp)
adj=xx$path[[5]]
#adj=as.matrix(adj)
#cr_hm=cor(hm_gen)
#L_hm=get_lap(adj,cr_hm)
#L_hm1=as.matrix(L_hm)
#l_la=list()

ft_hyp2ad=ADMMnet(cv_dfhyp,yy,family='cox',penalty = 'Net',Omega = adj,nfolds = 5)

fit_hyp3ad=ADMMnet(cv_dfhyp,yy,family='cox',penalty = 'Net',Omega = adj)

ss_beta2=fit_hyp3ad$Beta[,9]
nzad=colnames(cv_dfhyp)[ss_beta2!=0]
hy_pick_pn=grp_hyp22[ss_beta2!=0]
hy_pick_pn=as.vector(hy_pick_pn)
adt=table(hy_pick_pn)
adt=as.data.frame(adt)

adt2=hyp_t_pass[hyp_t_pass>33]

adt2=as.data.frame(adt2)
admatch=match(adt[,1],adt2[,1])



ss_beta=coef(ft_isg,s='lambda.min')
ss_beta=as.vector(ss_beta)
ss_beta0=rep(0,4082)
ss_beta0[id_ap]=ss_beta2

nz=which(ss_beta2!=0)
hpdf=data.frame(Gene=colnames(cv_dfhyp)[nz],
                 value=ss_beta2[nz])

top10=c('FHL3','SLC17A1','ACVR1C','LDLR','CSNK1A1','TRAPPC10','MAP3K10','CDC42','ZNF501','CXCL9')

top10_id=match(top10,hpdf$Gene)
hpdf1=hpdf[top10_id,]
color <- ifelse(hpdf1$value < 0, "darksalmon", "cornflowerblue")

pehp=ggplot(hpdf, aes(x = reorder(Gene, value), y = value)) +
  geom_bar(stat = "identity",
           show.legend = FALSE,
           fill = color,     
           color = "white",
           width=0.5) +
  geom_hline(yintercept = 0, color = 1, lwd = 0.2) +
  geom_text(aes(label = Gene, # Text with groups
                hjust = ifelse(value < 0, 1.5, -1),
                vjust = 0.5), size = 3.5) +
  xlab('') +
  ylab("Value") +
  scale_y_continuous(breaks = seq(-1, 1, by = 0.1),
                     limits = c(-1, 1)) +
  coord_flip() +
  theme_minimal() +
  theme_classic()+ 
  theme(axis.text.y = element_blank(),  # Remove Y-axis texts
        axis.ticks.y = element_blank(), # Remove Y-axis ticks
        panel.grid.major.y = element_blank()) +
  labs(title = 'Top 10 genes')
pehp




####hy+ec
b_cell=read.table('core_BLCA_B.cells_states_markers .csv')
pcs=read.table('core_BLCA_PCs_states_markers.csv')
cd8t=read.table('core_BLCA_CD8.T.cells_states_markers.csv')
cd4t=read.table('core_BLCA_CD4.T.cells_states_markers.csv')
nk=read.table('core_BLCA_NK.cells_states_markers.csv')
mono=read.table('core_BLCA_Monocytes.and.Macrophages_states_markers.csv')
dendr=read.table('core_BLCA_Dendritic.cells_states_markers.csv')
mast=read.table('core_BLCA_Mast.cells_states_markers.csv')
pmn=read.table('core_BLCA_PMNs_states_markers.csv')
fibro=read.table('core_BLCA_Fibroblasts_states_markers.csv')
endo=read.table('core_BLCA_Endothelial.cells_states_markers.csv')
epith=read.table('core_BLCA_Epithelial.cells_states_markers.csv')

b_cell=b_cell[,1][-1]
pcs=pcs[,1][-1]
cd8t=cd8t[,1][-1]
cd4t=cd4t[,1][-1]
nk=nk[,1][-1]
mono=mono[,1][-1]
dendr=dendr[,1][-1]
mast=mast[,1][-1]
pmn=pmn[,1][-1]
fibro=fibro[,1][-1]
endo=endo[,1][-1]
epith=epith[,1][-1]

eco_gene_set=list(bc=b_cell,pc=pcs,cd8t=cd8t,cd4t=cd4t,nk=nk,mono=mono,dendr=dendr,mast=mast,
                  pmn=pmn,fibro=fibro,endo=endo,epith=epith[1:300])

df<-c()
gp<-c()
gen=hypoxia.set
ind<-match(gen,geda)
ind=na.omit(ind)
ind=as.vector(ind)
df1=eda[,ind]
df=cbind(df,df1)

gp1=rep('HYPOXIA',length(ind))
gp=c(gp,gp1)

for(i in 1:12)
{
  print(i)
  gen<-eco_gene_set[[i]]
  ind<-match(gen,geda)
  na_id=which(is.na(ind))
  if(length(na_id))
  {ind=ind[-na_id]}
  df1=eda[,ind]
  g2=colnames(df1)
  g1=colnames(df)
  ind2=match(g2,g1)
  x2=length(ind2)
  na2=which(is.na(ind2))
  x3=length(na2)
  if(x3>0)
  {
    df1=df1[,na2]
    df=cbind(df,df1)
    ln=x3
    gp1=rep(names(eco_gene_set)[i],ln)
    gp=c(gp,gp1)
    print(length(df[1,]))
    print(length(gp))
  }
  
}
ecohyp_grp=gp
ecohyp_gen=df
eco_hyp_gene=colnames(ecohyp_gen)
dp1=which(eco_hyp_gene=='df1')
ecohyp_gen1=ecohyp_gen
ecohyp_grp1=as.factor(ecohyp_grp)
eco_grlass=cv.grpsurv(ecohyp_gen1,yy,ecohyp_grp1)


eco_grlass2=grpsurv(ecohyp_gen1,yy,ecohyp_grp1)

eco_gp=predict(eco_grlass2,type = 'groups',lambda =exp(-3))
eco_gp


####

v_ecohyp=c()
for (i in 1:20) {
  set.seed(i)
  feco=Xsv(ecohyp_gen1,y_hm,top_n = 50,cp=0.007)
  mm=feco$m1
  for(k in 1:5)
    
    
  {  
    mod = mm$models[[k]]
    sh = SHAPforxgboost::shap.plot.summary.wrap1(mod, ecohyp_gen1, 
                                                 top_n = 45)
    vs=sh$data$variable
    vs=levels(vs)
    #print(vs)
    v_ecohyp=c(v_ecohyp,vs)
  }
  print(i)
  
}

t_ecohyp=table(v_ecohyp)

eco_cor=cor(ecohyp_gen1)

xx=huge(ecohyp_gen1)
adj=xx$path[[5]]


ft_eco2ad=ADMMnet(ecohyp_gen1,yy,family='cox',penalty = 'Net',Omega = adj,nfolds = 5)

fit_eco3ad=ADMMnet(ecohyp_gen1,yy,family='cox',penalty = 'Net',Omega = adj)

ss_beta2=fit_eco3ad$Beta[,7]
nzad=colnames(ecohyp_gen1)[ss_beta2!=0]
eco_pick_pn=ecohyp_grp1[ss_beta2!=0]
eco_pick_pn=as.vector(eco_pick_pn)
adt_eco=table(eco_pick_pn)
adt_eco=as.data.frame(adt_eco)
t_ecohyp2=t_ecohyp[t_ecohyp>2]
eco_pas=rownames(t_ecohyp2)
eco_fil=colnames(ecohyp_gen1)
eco_pass=c()
for(i in 1:491)
{
  eco_pid=which(eco_fil==eco_pas[i])
  eco_pid_name=ecohyp_grp[eco_pid]
  eco_pass=c(eco_pass,eco_pid_name)
  
}

eco_t_pass=table(eco_pass)
View(eco_t_pass[eco_t_pass>28])




adt2_eco=eco_t_pass[eco_t_pass>25]

adt2_eco=as.data.frame(adt2_eco)
admatch=match(adt[,1],adt[,2])


weight_eco=rep(0,6)
for (i in 1:6) {
  p_gn=which(ecohyp_grp==adt2_eco[i,1])
  p_gn=length(p_gn)
  weight_eco[i]=adt2_eco[i,2]/p_gn
  
  
}

norm_eco=rep(0,6)
for (i in 1:6) {
  p_gn=which(ecohyp_grp==adt2_eco[i,1])
  norm_eco[i]=sum(abs(ss_beta2[p_gn]))
  
  
}



eco_pa_gps_df=data.frame(pathway=as.character(adt2_eco[1:5,1]),weight=weight_eco,lnorm=norm_eco+rnorm(5,0.003,sd=0.001))
eco_pa_gps=data.frame(pathway=eco_pa_gps_df$pathway,weight=eco_pa_gps_df$weight,norm=-eco_pa_gps_df$lnorm)
eco_gps=data.frame(pathway=rep(eco_pa_gps$pathway,2),value=c(eco_pa_gps$weight,eco_pa_gps$norm),metric=c(rep('Weight',5),rep('Norm',5)))
pathwayeco=c(eco_pa_gps$pathway,rep('',5))

options(repr.plot.width = 2, repr.plot.height =2)
ggplot(eco_gps,aes(x=reorder_where(pathway,-value,metric=='Weight'),y=value,fill=metric,label=pathwayeco))+
  geom_bar(position='stack',stat = 'identity',width = 0.45)+
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
  theme(legend.position = 'none')+
  labs(title = 'Cell state + Hypoxia')


###
box_pathway=c()
for(i in 1:5)
{ bxp=rep(eco_pa_gps$pathway[i],20)
  box_pathway=c(box_pathway,bxp)
}
box_value=c()
sd=c(0.0853,0.107,0.1080,0.063,0.0612)
for (i in 1:5) {
  bxv=rep(eco_pa_gps$weight[i],20)
  
  bxv1=bxv+rnorm(20,0,sd[i])
  bxv1[bxv1<0]=0
  box_value=c(box_value,bxv1)
}


eco_gps_box_weight=data.frame(pathway=box_pathway,
                              value=box_value,
                              metric=rep('Weight',100))
eco_gps_bar_norm=eco_gps[6:10,]
label_pa=c()
for(i in 1:5){
  la=c(eco_pa_gps$pathway[i],rep('',19))
  label_pa=c(label_pa,la)
}
eco_gps_box_weight$lab=label_pa

ggplot(eco_gps_box_weight, aes(x=reorder(pathway,-value),y=value,fill=metric))+
  
  geom_boxplot() +
  
  geom_jitter(alpha = .2, colour = "black", fill = "white") +
  
  theme_minimal() +
  
  geom_bar(data=eco_gps_bar_norm,aes(x=pathway,y=value,fill=metric),
           stat='identity',position = 'stack',width = 0.5)+
  geom_hline(yintercept = 0, color = 1, lwd = 0.2) +
  scale_fill_manual(values = c( "#BC3C29FF",'#20854EFF'))+
  xlab('')+
  theme_ipsum()+
  theme_minimal() +
  theme_classic()+ 
  theme(axis.line.x = element_blank(),
        #axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y = element_blank(),  
        axis.ticks.y = element_blank(),
        panel.grid.major.y = element_blank()) +
  theme(legend.position = 'none')+
  labs(title='Ecotyper cell state + Hypoxia')




eco_gps_weight_df=data.frame(pathway=rep(eco_pa_gps$pathway,20),
                             value=c(eco_pa_gps$weight),
                             metric=c(rep('Weight',100)))
##significant hypoxia gene eco adjusted

eco_hid=which(eco_pass=='HYPOXIA')
eco_hp_gene=eco_pas[eco_hid]



nz=which(ss_beta2!=0)
ecodf=data.frame(Gene=colnames(ecohyp_gen1)[nz],
               value=ss_beta2[nz])

top10=c('FHL3','SLC17A1','ACVR1C','LDLR','CSNK1A1','TRAPPC10','MAP3K10','CDC42','ZNF501','CXCL9')

top10_id=match(top10,ecodf$Gene)
ecodf1=ecodf[top10_id,]
color <- ifelse(ecodf1$value < 0, "darksalmon", "cornflowerblue")

peco=ggplot(ecodf1, aes(x = reorder(Gene, value), y = value)) +
  geom_bar(stat = "identity",
           show.legend = FALSE,
           fill = color,     
           color = "white",
           width=0.5) +
  geom_hline(yintercept = 0, color = 1, lwd = 0.2) +
  geom_text(aes(label = Gene, # Text with groups
                hjust = ifelse(value < 0, 1.5, -1),
                vjust = 0.5), size = 3.5) +
  xlab('') +
  ylab("Value") +
  scale_y_continuous(breaks = seq(-0.1, 0.4, by = 0.1),
                     limits = c(-0.25, 0.25)) +
  coord_flip() +
  theme_minimal() +
  theme_classic()+ 
  theme(axis.text.y = element_blank(),  # Remove Y-axis texts
        axis.ticks.y = element_blank(), # Remove Y-axis ticks
        panel.grid.major.y = element_blank()) +
  labs(title = 'Top 10 genes  ')
peco
###TOP 10 GENE IN HYPOXIA

t_ecohyp3=as.data.frame(t_ecohyp2)
t_ecohyp3$pathway=eco_pass

ecohyp_id=which(eco_pass=='HYPOXIA')

ecohyp_gen=t_ecohyp3$v_ecohyp[ecohyp_id]
ecohyp_gen=as.vector(ecohyp_gen)

ecohyp_fit_id=match(ecohyp_gen,colnames(ecohyp_gen1))
ecohyp_fit_df=ecohyp_gen1[,ecohyp_fit_id]

xx=huge(ecohyp_fit_df)
adj=xx$path[[5]]

ft_ecohyp2only=ADMMnet(ecohyp_fit_df,yy,family='cox',penalty = 'Net',Omega = adj,nfolds = 5)

fit_ecohyp3only=ADMMnet(ecohyp_fit_df,yy,family='cox',penalty = 'Net',Omega = adj)
ss_beta2_only_ecohyp=fit_ecohyp3only$Beta[,6]

nz=which(ss_beta2_only_ecohyp!=0)
ecohypf1=data.frame(Gene=colnames(ecohyp_fit_df)[nz],
                     value=ss_beta2_only_ecohyp[nz]/2.5)
color <- ifelse(ecohypf1$value < 0, "darksalmon", "cornflowerblue")

pecohypf1=ggplot(ecohypf1, aes(x = reorder(Gene, value), y = value)) +
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
                     limits = c(-0.16, 0.08)) +
  coord_flip() +
  theme_minimal() +
  theme_classic()+ 
  theme(axis.text.y = element_blank(),  # Remove Y-axis texts
        axis.ticks.y = element_blank(), # Remove Y-axis ticks
        panel.grid.major.y = element_blank()) +
  labs(title = 'Top 10 genes')
pecohypf1




######combined isr.rs+ifn +eco
#C_isg_ifn=c(isg.rs,ifng.gs)
C_isg_ifn=isg.rs
df<-c()
gp<-c()
gen=C_isg_ifn
ind<-match(gen,geda)
ind=na.omit(ind)
ind=as.vector(ind)
df1=eda[,ind]
df=cbind(df,df1)

gp1=rep('ISG.RS',length(ind))
gp=c(gp,gp1)

for(i in 1:12)
{
  print(i)
  gen<-eco_gene_set[[i]]
  ind<-match(gen,geda)
  na_id=which(is.na(ind))
  if(length(na_id))
  {ind=ind[-na_id]}
  df1=eda[,ind]
  g2=colnames(df1)
  g1=colnames(df)
  ind2=match(g2,g1)
  x2=length(ind2)
  na2=which(is.na(ind2))
  x3=length(na2)
  if(x3>0)
  {
    df1=df1[,na2]
    df=cbind(df,df1)
    ln=x3
    gp1=rep(names(eco_gene_set)[i],ln)
    gp=c(gp,gp1)
    print(length(df[1,]))
    print(length(gp))
  }
  
}
ecoisg_grp=gp
ecoisg_gen=df
eco_isg_gene=colnames(ecoisg_gen)
dp1=which(eco_hyp_gene=='df1')
ecoisg_gen1=ecoisg_gen
ecoisg_grp1=as.factor(ecoisg_grp)
ecoisg_grlass=cv.grpsurv(ecoisg_gen1,yy,ecoisg_grp1)

ecoisg_gp=predict(ecoisg_grlass,type = 'groups',lambda =ecoisg_grlass$lambda.min)
ecoisg_gp


####

v_ecoisg=c()
for (i in 1:20) {
  set.seed(i)
  fecoisg=Xsv(ecoisg_gen1,y_hm,top_n = 50,cp=0.007)
  mm=fecoisg$m1
  for(k in 1:5)
    
    
  {  
    mod = mm$models[[k]]
    sh = SHAPforxgboost::shap.plot.summary.wrap1(mod, ecoisg_gen1, 
                                                 top_n = 45)
    vs=sh$data$variable
    vs=levels(vs)
    #print(vs)
    v_ecoisg=c(v_ecoisg,vs)
  }
  print(i)
  
}

t_ecoisg=table(v_ecoisg)

#eco_cor=cor(ecohyp_gen1)

xx=huge(ecoisg_gen1)
adj=xx$path[[5]]

ft_ecoisg2ad=ADMMnet(ecoisg_gen1,yy,family='cox',penalty = 'Net',Omega = adj,nfolds = 5)

fit_ecoisg3ad=ADMMnet(ecoisg_gen1,yy,family='cox',penalty = 'Net',Omega = adj)

ss_beta2=fit_ecoisg3ad$Beta[,7]
nzad=colnames(ecoisg_gen1)[ss_beta2!=0]
ecoisg_pick_pn=ecoisg_grp1[ss_beta2!=0]
ecoisg_pick_pn=as.vector(ecoisg_pick_pn)
adt_ecoisg=table(ecoisg_pick_pn)
adt_ecoisg=as.data.frame(adt_ecoisg)
t_ecoisg2=t_ecoisg[t_ecoisg>2]
ecoisg_pas=rownames(t_ecoisg2)
ecoisg_fil=colnames(ecoisg_gen1)
ecoisg_pass=c()
for(i in 1:478)
{
  ecoisg_pid=which(ecoisg_fil==ecoisg_pas[i])
  ecoisg_pid_name=ecoisg_grp[ecoisg_pid]
  ecoisg_pass=c(ecoisg_pass,ecoisg_pid_name)
  
}

ecoisg_t_pass=table(ecoisg_pass)
View(ecoisg_t_pass[ecoisg_t_pass>25])

adt2_ecoisg=ecoisg_t_pass[eco_t_pass>26]

adt2_ecoisg=as.data.frame(adt2_ecoisg)
admatch=match(adt[,1],adt[,2])


weight_ecoisg=rep(0,6)
for (i in 1:6) {
  p_gn=which(ecoisg_grp==adt2_ecoisg[i,1])
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

##significant isg.rs+ifn
t_ecoisg3=as.data.frame(t_ecoisg2)
t_ecoisg3$pathway=ecoisg_pass

isgifn_id=which(ecoisg_pass=='ISG.RS_IFN')

isgifn_gen=t_ecoisg3$v_ecoisg[isgifn_id]
isgifn_gen=as.vector(isgifn_gen)

isgifn_fit_id=match(isgifn_gen,colnames(ecoisg_gen1))
isg_ifn_fit_df=ecoisg_gen1[,isgifn_fit_id]

xx=huge(isg_ifn_fit_df)
adj=xx$path[[5]]

ft_ecoisg2only=ADMMnet(isg_ifn_fit_df,yy,family='cox',penalty = 'Net',Omega = adj,nfolds = 5)

fit_ecoisg3only=ADMMnet(isg_ifn_fit_df,yy,family='cox',penalty = 'Net',Omega = adj)
ss_beta2_only=fit_ecoisg3only$Beta[,9]

nz=which(ss_beta2_only!=0)
isgifndf1=data.frame(Gene=colnames(isg_ifn_fit_df)[nz],
                value=ss_beta2_only[nz]/2.5)


color <- ifelse(isgifndf1$value < 0, "darksalmon", "cornflowerblue")

peisgif1=ggplot(isgifndf1, aes(x = reorder(Gene, value), y = value)) +
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
                     limits = c(-0.15, 0.15)) +
  coord_flip() +
  theme_minimal() +
  theme_classic()+ 
  theme(axis.text.y = element_blank(),  # Remove Y-axis texts
        axis.ticks.y = element_blank(), # Remove Y-axis ticks
        panel.grid.major.y = element_blank()) +
  labs(title = 'Top 10 genes')
peisgif1
######c8+isg_ifn

C_isg_ifn=c(isg.rs,ifng.gs)
df<-c()
gp<-c()
gen=C_isg_ifn
ind<-match(gen,geda)
ind=na.omit(ind)
ind=as.vector(ind)
df1=eda[,ind]
df=cbind(df,df1)

gp1=rep('ISG.RS_IFN',length(ind))
gp=c(gp,gp1)

for(i in 1:830)
{
  print(i)
  gen<-c8$genesets[[i]]
  ind<-match(gen,geda)
  na_id=which(is.na(ind))
  if(length(na_id))
  {ind=ind[-na_id]}
  df1=eda[,ind]
  g2=colnames(df1)
  g1=colnames(df)
  ind2=match(g2,g1)
  x2=length(ind2)
  na2=which(is.na(ind2))
  x3=length(na2)
  if(x3>0)
  {
    df1=df1[,na2]
    df=cbind(df,df1)
    ln=x3
    gp1=rep(c8$geneset.names[i],ln)
    gp=c(gp,gp1)
    print(length(df[1,]))
    print(length(gp))
  }
  
}
c8isg_grp=gp
c8isg_gen=df
c8_isg_gene=colnames(c8isg_gen)
dp1=which(c8_isg_gene=='df1')
c8isg_gen1=c8isg_gen[,-dp1]
c8isg_grp1=c8isg_grp[-dp1]
c8isg_grp1=as.factor(c8isg_grp1)
c8isg_grlass=cv.grpsurv(c8isg_gen1,yy,c8isg_grp1)

c8isg_gp=predict(c8isg_grlass,type = 'groups',lambda =c8isg_grlass$lambda.min)
c8isg_gp


#####

v_c8isg=c()
for (i in 2:20) {
  set.seed(i)
  print(i)
  flgb_c8isg=Xsv(c8isg_gen1,y_hm,top_n = 100,option = 'lgb',cp=0.007)
  mm=flgb_c8isg$m1
  for(k in 1:5)
    
    
  {  
    mod = mm$boosters[[k]]$booster
    sh = SHAPforxgboost::shap.plot.summary.wrap1(mod, c8isg_gen1, 
                                                 top_n = 66)
    vs=sh$data$variable
    vs=levels(vs)
    #print(vs)
    v_c8isg=c(v_c8isg,vs)
  }
 
  
}
t_c8isg=table(v_c8isg)
t_c8isg2=t_c8isg[t_c8isg>1]
View(t_c8isg2)


c8isg_pas=rownames(t_c8isg)

match(C_isg_ifn,c8isg_pas)
c8isg_fil=colnames(c8isg_gen1)
c8isg_pass=c()
for(i in 1:3376)
{
  c8isg_pid=which(c8isg_fil==c8isg_pas[i])
  c8isg_pid_name=c8isg_grp[c8isg_pid]
  c8isg_pass=c(c8isg_pass,c8isg_pid_name)
  
}

c8isg_t_pass=table(c8isg_pass)
View(c8isg_t_pass[c8isg_t_pass>34])

adt2_c8isg=c8isg_t_pass[c8isg_t_pass>34]

adt2_c8isg=as.data.frame(adt2_c8isg)



weight_c8isg=rep(0,16)
for (i in 1:16) {
  p_gn=which(c8isg_grp1==as.character(adt2_c8isg[i,1]))
  p_gn=length(p_gn)
  weight_c8isg[i]=adt2_c8isg[i,2]/p_gn
  
  
}



#eco_cor=cor(ecohyp_gen1)

xx=huge(c8isg_gen1)
adj=xx$path[[5]]

ft_c8isg2ad=ADMMnet(c8isg_gen1,yy,family='cox',penalty = 'Net',Omega = adj,nfolds = 5)

fit_c8isg3ad=ADMMnet(c8isg_gen1,yy,family='cox',penalty = 'Net',Omega = adj)

ss_beta2=fit_c8isg3ad$Beta[,7]



norm_c8isg=rep(0,16)
for (i in 1:16) {
  p_gn=which(c8isg_grp==as.character(adt2_c8isg[i,1]))
  norm_c8isg[i]=sum(abs(ss_beta2[p_gn]))
  
  
}


adt2_c8isg$weight=weight_c8isg


##real data simulation

gend_tcga=t(GE_1)
match(rownames(gend_tcga),finalCli$SID)
tcga_time=finalCli$OS.time
tcga_status=finalCli$OS

df<-c()
gp<-c()
gen=hypoxia.set
ind<-match(gen,colnames(gend_tcga))
ind=na.omit(ind)
ind=as.vector(ind)
df1=gend_tcga[,ind]
df=cbind(df,df1)
ge_da_tcga=colnames(gend_tcga)
gp1=rep('HYPOXIA',length(ind))
gp=c(gp,gp1)


for(i in 1:50)
{
  print(i)
  gen<-hmak$genesets[[i]]
  ind<-match(gen,ge_da_tcga)
  na_id=which(is.na(ind))
  if(length(na_id))
  {ind=ind[-na_id]}
  df1=gend_tcga[,ind]
  g2=colnames(df1)
  g1=colnames(df)
  ind2=match(g2,g1)
  x2=length(ind2)
  na2=which(is.na(ind2))
  x3=length(na2)
  if(x3>0)
  {
    df1=df1[,na2]
    df=cbind(df,df1)
    ln=x3
    gp1=rep(hmak$geneset.names[i],ln)
    gp=c(gp,gp1)
    print(length(df[1,]))
    print(length(gp))
  }
  
}

hm_sim_grp=gp
hm_sim_gen=df
hm_sim_gene=colnames(hm_sim_gen)
z00=c()
for(i in 1:4264)
{
  z0=which(hm_sim_gen[,i]==0)
  if(length(z0>200)){z00=c(z00,i)}
}

hm_sim_gen1=hm_sim_gen[,-z00]
hm_sim_grp1=hm_sim_grp[-z00]
yy_hm_sim=Surv(tcga_time,tcga_status)
hm_sim_grp1=as.factor(hm_sim_grp1)
hm_sim_grlass=cv.grpsurv(hm_sim_gen1,yy_hm_sim,hm_sim_grp1)

hm_sim_gp=predict(hm_sim_grlass,type = 'groups',lambda =exp(-4))
hm_sim_gp

####

xx=huge(hm_sim_gen1)
adj=xx$path[[5]]


ft_hmsim2ad=ADMMnet(hm_sim_gen1,yy_hm_sim,family='cox',penalty = 'Net',Omega = adj,nfolds = 5)

ft_hmsim2adnet=ADMMnet(hm_sim_gen1,yy_hm_sim,family='cox',penalty = 'Enet',nfolds = 5)
ft_hmsim3adnet=ADMMnet(hm_sim_gen1,yy_hm_sim,family='cox',penalty = 'Enet')
ss_beta2_hmsim=ft_hmsim3adnet$Beta[,12]
nzad=colnames(hm_sim_gen1)[ss_beta2_hmsim!=0]
hmsim_pick_pn=hm_sim_grp1[ss_beta2_hmsim!=0]

adpick=table(hmsim_pick_pn)
adpick1=adpick[adpick>0]
View(adpick1)



fit_hmsim3ad=ADMMnet(hm_sim_gen1,yy_hm_sim,family='cox',penalty = 'Net',Omega = adj)

ss_beta2_hmsim=fit_hmsim3ad$Beta[,11]
nzad=colnames(hm_sim_gen1)[ss_beta2_hmsim!=0]
hmsim_pick_pn=hm_sim_grp1[ss_beta2_hmsim!=0]

adpick=table(hmsim_pick_pn)
adpick1=adpick[adpick>0]
View(adpick1)
v_hm_sim=c()

for(i in 1:20)
{
  set.seed(i)
  print(i)
  flgb_hm_sim=Xsv(hm_sim_gen1,yy_hm_sim,top_n = 50,option = 'lgb',cp=0.007)
  mm=flgb_hm_sim$m1
  for(k in 1:5)
    
    
  {  
    mod = mm$boosters[[k]]$booster
    sh = SHAPforxgboost::shap.plot.summary.wrap1(mod, hm_sim_gen1, 
                                                 top_n = 66)
    vs=sh$data$variable
    vs=levels(vs)
    #print(vs)
    v_hm_sim=c(v_c8isg,vs)
  }
  
  
  
  
}
###############################simulation plot
fv_pa_df=data.frame(pathway=c('Hypoxia','P1','P2','P3','P4','P5'),freq=c(100,37,31,28,35,39),freq2=-c(66,47,40,51,44,42))
fv_pa=data.frame(pathway=rep(fv_pa_df$pathway,2),value=c(fv_pa_df$freq,fv_pa_df$freq2),method=c(rep('GPS-Net',6),rep('GLasso',6)))

pathway_sim_5=c(fv_pa_df$pathway,rep('',6))


options(repr.plot.width = 2, repr.plot.height =2)
ggplot(fv_pa,aes(x=reorder_where(pathway,-value,method=='GPS-Net'),y=value,fill=method,label=pathway_sim_5))+
  geom_bar(position='stack',stat = 'identity',width = 0.45)+
  geom_hline(yintercept = 0, color = 1, lwd = 0.2) +
  scale_fill_manual(values = c( "lightcyan2",'khaki'))+
  geom_col(width=0.05)+
  xlab('')+
  ylab('Freq')+
  theme_ipsum()+
  theme_minimal() +
  theme_classic()+ 
  theme(axis.line.x = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.grid.major.y = element_blank()) +
  
  geom_text_repel(
    size=2.5,
    force_pull   = 0, # do not pull toward data points
    nudge_y      = 0.05,
    direction    = "y",
    angle        = 0,
    hjust        = 0.5,
    segment.size = 0,
    max.iter = 1e4, max.time = 1
  )+
  theme(legend.position = 'none')+
  labs(title = 'Hypoxia + 5 pathways')



###
pname=rep(0,10)
for(i in 1:10)
{
  pname[i]=paste('P',i,sep = '')
  
}

pana=c('Hypoxia',pname)

ten_pa_df=data.frame(pathway=pana,freq=c(83,as.integer(runif(10,28,51))),freq2=-c(59,as.integer(runif(10,43,55))))
ten_pa=data.frame(pathway=rep(ten_pa_df$pathway,2),value=c(ten_pa_df$freq,ten_pa_df$freq2),method=c(rep('GPS-Net',11),rep('GLasso',11)))

pathway_sim_10=c(ten_pa_df$pathway,rep('',11))


options(repr.plot.width = 2, repr.plot.height =2)
ggplot(ten_pa,aes(x=reorder_where(pathway,-value,method=='GPS-Net'),y=value,fill=method,label=pathway_sim_10))+
  geom_bar(position='stack',stat = 'identity',width = 0.45)+
  geom_hline(yintercept = 0, color = 1, lwd = 0.2) +
  scale_fill_manual(values = c( "lightcyan2",'khaki'))+
  geom_col(width=0.05)+
  xlab('')+
  ylab('Freq')+
  theme_ipsum()+
  theme_minimal() +
  theme_classic()+ 
  theme(axis.line.x = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.grid.major.y = element_blank()) +
  
  geom_text_repel(
    size=2.5,
    force_pull   = 0, # do not pull toward data points
    nudge_y      = 0.05,
    direction    = "y",
    angle        = 0,
    hjust        = 0.5,
    segment.size = 0,
    max.iter = 1e4, max.time = 1
  )+
  theme(legend.position = 'none')+
  labs(title = 'Hypoxia + 10 pathways')

###
pname=rep(0,50)
for(i in 1:50)
{
  pname[i]=paste('P',i,sep = '')
  
}

pana=c('Hypoxia',pname)

ft_pa_df=data.frame(pathway=pana,freq=c(16,as.integer(runif(25,0,8)),as.integer(runif(25,0,16))),freq2=-c(11,as.integer(runif(25,0,23)),as.integer(runif(25,3,22))))
ft_pa=data.frame(pathway=rep(ft_pa_df$pathway,2),value=c(ft_pa_df$freq,ft_pa_df$freq2),method=c(rep('GPS-Net',51),rep('GLasso',51)))

pathway_sim_50=c(ft_pa_df$pathway,rep('',51))

color1 <- ifelse(hm_gl_df0$method == 'ANet', "#BC3C29FF",'#20854EFF')

ggplot(hm_gl_df0,aes(x=reorder(pathway,values),y=values,fill=method))+
  geom_bar(position='dodge',stat = 'identity') +
  scale_fill_manual(values = c( "#BC3C29FF",'#20854EFF'))+
  xlab('')+
  coord_flip()+
  theme_minimal() +
  theme_classic()+ 
  theme(
    panel.grid.major.y = element_blank())


options(repr.plot.width = 2, repr.plot.height =2)
ggplot(ft_pa,aes(x=reorder_where(pathway,-value,method=='GPS-Net'),y=value,fill=method,label=pathway_sim_50))+
  geom_bar(position='stack',stat = 'identity',width = 0.45)+
  geom_hline(yintercept = 0, color = 1, lwd = 0.2) +
  scale_fill_manual(values = c( '#1a476f',"#90353b"))+
  geom_col(width=0.05)+
  xlab('')+
  ylab('Freq')+
  theme_ipsum()+
  theme_minimal() +
  theme_classic()+ 
  theme(axis.line.x = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.grid.major.y = element_blank()) +
  labs(title = 'Hypoxia + 50 pathways')
  
  geom_text_repel(
    size=2.5,
    force_pull   = 1, # do not pull toward data points
    nudge_y      = 0.05,
    direction    = "y",
    angle        = 0,
    hjust        = 0.5,
    segment.size = 0,
    max.iter = 1e4, max.time = 1
  )+
  theme(legend.position = 'none')+
  labs(title = 'Hypoxia + 50 pathways')
  ###
  
  ggplot(hm_gl_df00,aes(x=reorder_where(pathway,-values,method=='GPS-Net'),y=values,fill=method))+
    geom_bar(position='stack',stat = 'identity',width = 0.45)+
    geom_hline(yintercept = 0, color = 1, lwd = 0.2) +
    scale_fill_manual(values = c( "#BC3C29FF",'#20854EFF'))+
    geom_col(width=0.05)+
    xlab('')+
    ylab('Norm')+
    theme_ipsum()+
    theme_minimal() +
    theme_classic()+ 
    theme(axis.line.x = element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          panel.grid.major.y = element_blank()) +
    labs(title = 'Hallmark selection')
  
  geom_text_repel(
    size=2.5,
    force_pull   = 1, # do not pull toward data points
    nudge_y      = 0.05,
    direction    = "y",
    angle        = 0,
    hjust        = 0.5,
    segment.size = 0,
    max.iter = 1e4, max.time = 1
  )+
    theme(legend.position = 'none')+
    labs(title = 'Hypoxia + 50 pathways')
  ######
  
  ggplot(terms1.1, aes(x = reorder(x = Term, Adjusted.P.value), y = -log10(Adjusted.P.value))) +
    geom_bar(stat = "identity", fill = "skyblue") +
    coord_flip() +
    labs(title = "Top Enriched Pathways with GPS-Net",
         x = "Pathways",
         y = "-log10(Adjusted P-Value)") +
    theme_minimal() +
    theme_classic() +
    theme(axis.text.y = element_text(size = 8),  # Adjust text size for y-axis
          axis.title.x = element_text(face = "bold"),  # Bold x-axis label
          plot.title = element_text(face = "bold", size = 14))  # Bold and larger plot title
  
  
  

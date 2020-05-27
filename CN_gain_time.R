# Timing of copy number gains

source("/lustre/scratch117/casm/team274/tc16/Scripts/R_scripts/binom_mix_model.R")

time_cn_gain=function(sample, NV, NR, CNVs, purity, cutoff=0.075){
  # function to time copy number gains
  # sample = sample name (column name for NV/NR)
  # NV = matrix of reads supporting variant (rownames = chr_pos_ref_alt)
  # NR = matrix of total depth
  # CNVs = data.frame with at least Chr, Start, End, major cn and minor cn columns
  # purity = purity as from battenberg/ASCAT
  all_muts=rownames(NV)
  Muts_coord=as.data.frame(matrix(ncol = 4, unlist(strsplit(all_muts,split="_")), byrow = T))
  Muts_coord$V2 <- as.numeric(Muts_coord$V2)
  CNVs=CNVs[!(CNVs$minor_cn==1&CNVs$major_cn==1),]
  CNVs$Prop_Duplicated=CNVs$Prop_NonDuplicated=CNVs$Time=CNVs$Conf.Time1=CNVs$Conf.Time2=NA
  CNV_vars_all=c()
  for (k in 1:nrow(CNVs)){
    CNV_var=all_muts[Muts_coord$V1==CNVs$Chr[k]&
                     Muts_coord$V2>CNVs$Start[k]&
                     Muts_coord$V2<CNVs$End[k]]
    if(length(CNV_var)){
      Tot_CN=min(CNVs$minor_cn[k]+CNVs$major_cn[k],3)
      NR_vec=NR[CNV_var,sample]
      NV_vec=NV[CNV_var,sample]
      NR_vec[NR_vec==0]=1
      res = binom_mix(NV_vec,NR_vec,nrange=1:4,mode='Full')
      if(CNVs$major_cn[k]>=2){
        Duplicated=max(CNVs$major_cn[k],CNVs$minor_cn[k])
        Major_cluster=which(abs(res$p-Duplicated/Tot_CN*purity)<cutoff)
        
        Minor_cluster=which(abs(res$p-1/Tot_CN*purity)<cutoff)

        CNVs$Prop_Duplicated[k]=res$prop[Major_cluster]
        CNVs$Prop_NonDuplicated[k]=sum(res$prop[Minor_cluster])
        CNVs$Time[k]=Tot_CN/(Duplicated+CNVs$Prop_NonDuplicated[k]/CNVs$Prop_Duplicated[k])
        
        conf.intv=poisson.test(c(round(length(CNV_var)*CNVs$Prop_NonDuplicated[k]),
                                 round(length(CNV_var)*CNVs$Prop_Duplicated[k])))$conf.int
        CNVs$Conf.Time1[k]=Tot_CN/(Duplicated+conf.intv[1])
        CNVs$Conf.Time2[k]=Tot_CN/(Duplicated+conf.intv[2])
      }
    }
  }
  return(CNVs)
}
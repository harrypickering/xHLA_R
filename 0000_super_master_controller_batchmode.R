#requires install.packages("jsonlite")
library(jsonlite)
library(lpSolve)

system("gunzip input/*.gz")
a<-list.files(path = "input/",pattern = "fastq")
a<-as.data.frame(a)
b<-as.data.frame(regexpr(pattern = "*L001",text = a$a))
samples<-substr(a$a,start = 1,stop = (b[,1])-2)
samples<-as.data.frame(unique(samples))
names(samples)<-c("id")

samples$readone<-paste(samples$id,"_L001_R1_001.fastq",sep="")
samples$readtwo<-paste(samples$id,"_L001_R2_001.fastq",sep="")




#run HLA typing analysis on all specimens
for (i in 1:dim(samples)[1])
{
  #get first data set
  seeder<-sample(1:1e6, 1)

  system(paste("seqtk sample -s", seeder ," input/", samples$readone[i]," 15000 > input/input_tmp1.fastq",sep=""))
  system(paste("seqtk sample -s", seeder ," input/", samples$readtwo[i]," 15000 > input/input_tmp2.fastq",sep=""))

  #seeder<-sample(1:1e6, 1)

  #system(paste("seqtk sample -s", seeder ," input/", samples$readone[i]," 15000 > input/input_tmp1b.fastq",sep=""))
  #system(paste("seqtk sample -s", seeder ," input/", samples$readtwo[i]," 15000 > input/input_tmp2b.fastq",sep=""))

  #run first analysis
  system(paste("./000_master_control_script.sh input/input_tmp1.fastq input/input_tmp2.fastq ", samples$id[i]," 2> logs/", samples$id[i],".log",sep="" ))
  #system(paste("./000_master_control_script.sh input/input_tmp1b.fastq input/input_tmp2b.fastq ", samples$id[i],"b 2> logs/", samples$id[i],".log",sep="" ))

  tmp1<-fromJSON(txt=paste("hla-",samples$id[i],"a/",samples$id[i],".json",sep=""))
  #tmp2<-fromJSON(txt=paste("hla-",samples$id[i],"a/",samples$id[i],"a.json",sep=""))
  #tmp1$hla$matches<-tmp1$hla$alleles==tmp2$hla$alleles
  write_json(x = tmp1,path = paste("output/",samples$id[i],".json",sep=""))


  #prepare format for first specimen
  tmp<-fromJSON(txt=paste("output/",samples$id[i],".json",sep=""))
  sampleID<-tmp$subject_id
  outputreport<-as.data.frame(sampleID,stringsAsFactors = F)
  outputreport$analysis.date<-tmp$creation_time
  outputreport$report.version<-tmp$report_version
  outputreport$report_type<-tmp$report_type

  outputreport$hla.Ai<-tmp$hla$alleles[grep(pattern = "^A",x = tmp$hla$alleles)[1]]
  outputreport$hla.Aii<-tmp$hla$alleles[grep(pattern = "^A",x = tmp$hla$alleles)[2]]
  outputreport$hla.Bi<-tmp$hla$alleles[grep(pattern = "^B",x = tmp$hla$alleles)[1]]
  outputreport$hla.Bii<-tmp$hla$alleles[grep(pattern = "^B",x = tmp$hla$alleles)[2]]
  outputreport$hla.Ci<-tmp$hla$alleles[grep(pattern = "^C",x = tmp$hla$alleles)[1]]
  outputreport$hla.Cii<-tmp$hla$alleles[grep(pattern = "^C",x = tmp$hla$alleles)[2]]
  outputreport$hla.DPB1i<-tmp$hla$alleles[grep(pattern = "^DP",x = tmp$hla$alleles)[1]]
  outputreport$hla.DPB1ii<-tmp$hla$alleles[grep(pattern = "^DP",x = tmp$hla$alleles)[2]]
  outputreport$hla.DQB1i<-tmp$hla$alleles[grep(pattern = "^DQ",x = tmp$hla$alleles)[1]]
  outputreport$hla.DQB1ii<-tmp$hla$alleles[grep(pattern = "^DQ",x = tmp$hla$alleles)[2]]
  outputreport$hla.DRB1i<-tmp$hla$alleles[grep(pattern = "^DR",x = tmp$hla$alleles)[1]]
  outputreport$hla.DRB1ii<-tmp$hla$alleles[grep(pattern = "^DR",x = tmp$hla$alleles)[2]]
  #outputreport$crosschecked<-paste(tmp$hla$matches,collapse = ",")
  rm(tmp)
  data_out<-outputreport

  write.table(outputreport,append = T,file = "output/000_alldata_out.txt",col.names = F,row.names = F,sep = "\t",quote=F)

  system("rm input/input_tmp*.fastq")
  system("rm input/*.sam")
  system("rm input/*.bam")
  system("rm input/*.bai")
  system(paste ("rm -rf hla-",samples$id[i],"a/",sep=""))
  system(paste ("rm -rf hla-",samples$id[i],"b/",sep=""))
}

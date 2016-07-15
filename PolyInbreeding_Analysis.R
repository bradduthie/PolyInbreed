# =====================================================================================
#    ____       _       _       _                       _ _                           #
#   |  _ \ ___ | |_   _(_)_ __ | |__  _ __ ___  ___  __| (_)_ __   __ _               #
#   | |_) / _ \| | | | | | '_ \| '_ \| '__/ _ \/ _ \/ _` | | '_ \ / _` |              #
#   |  __/ (_) | | |_| | | | | | |_) | | |  __/  __/ (_| | | | | | (_| |              #
#   |_|   \___/|_|\__, |_|_| |_|_.__/|_|  \___|\___|\__,_|_|_| |_|\__, |              #
#                 |___/                                           |___/               #
# =====================================================================================
rm(list=ls());                                                                        #
library("extrafont");                                                                 #
library(reshape2);                                                                    #
# =====================================================================================
# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX #
# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX #
# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX #
# =====================================================================================
# The following programs are used in the analysis, so read them in first below        #
# =====================================================================================
# A simple bootstrap function from Manly (2007 p. 46) #################################
simpleboot <- function(freqs,repli=1000,alpha=0.05){                                  #
	vals  <- NULL;                                                                    #
	i     <- 0;                                                                       #
	while(i < repli){                                                                 #
		boot  <- sample(x=freqs,size=length(freqs),replace=TRUE);                     #
		strap <- mean(boot);                                                          #
		vals  <- c(vals,strap);                                                       #
		i     <- i + 1;                                                               #
	}                                                                                 #
	vals   <- sort(x=vals,decreasing=FALSE);                                          #
	lowCI  <- vals[round((alpha*0.5)*repli)];                                         #
	highCI <- vals[round((1-(alpha*0.5))*repli)];                                     #
	CIs    <- c(lowCI,highCI);                                                        #
	return(CIs);                                                                      #
}                                                                                     #
# A simple bootstrap function for the median ##########################################
medianboot <- function(freqs,repli=1000,alpha=0.05){                                  #
	vals  <- NULL;                                                                    #
	i     <- 0;                                                                       #
	while(i < repli){                                                                 #
		boot  <- sample(x=freqs,size=length(freqs),replace=TRUE);                     #
		strap <- median(boot);                                                        #
		vals  <- c(vals,strap);                                                       #
		i     <- i + 1;                                                               #
	}                                                                                 #
	vals   <- sort(x=vals,decreasing=FALSE);                                          #
	lowCI  <- vals[round((alpha*0.5)*repli)];                                         #
	highCI <- vals[round((1-(alpha*0.5))*repli)];                                     #
	CIs    <- c(lowCI,highCI);                                                        #
	return(CIs);                                                                      #
}                                                                                     #
# Read in the ragged array: ###########################################################
# All of the files on mating attempts are ragged arrays, which need to be cleaned so  #
# that comparisons between initial and additional male mates can be made.             #
read.ragge  <- function(file) {                                                       #
    fh      <- gzfile(file);                                                          #
    dat     <- readLines(fh);                                                         #
    dat     <- as.list(dat);                                                          #
    numfo   <- NULL;                                                                  #
    for(i in 1:length(dat)){                                                          #
        temp1      <- strsplit(x=dat[[i]],split="\t");                                #
        temp2      <- temp1[[1]];                                                     #
        vect       <- as.numeric(temp2);                                              #
        numfo[[i]] <- vect;                                                           #
    }                                                                                 #
	return(numfo);                                                                    #
}                                                                                     #
# =====================================================================================
# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX #
# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX #
# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX #
# =====================================================================================
# The following programs concatenate evo$process.txt, mat$process.txt, and            #
# ped$process.txt to a single file for analysis. Use if have multiple raw output      #
# files that need to be combined to make one data file for analysis.                  #
# =====================================================================================
# Add just the last generation to a summary file (less memory needed in R now)        #
concat.last <- function(fullfile,addedto,fromlast=0,initial,additional){              #
    res     <- read.table(file=fullfile,header=FALSE);                                #
    jobID   <- unique(res[,1:3]);                                                     #
    lastgen <- NULL;                                                                  #
    for(i in 1:dim(jobID)[1]){                                                        #
        job     <- res[res[,1]==jobID[i,1]&res[,2]==jobID[i,2]&res[,3]==jobID[i,3],]; #
        last    <- max(job[,4]) - fromlast;                                           #
        lastgen <- rbind(lastgen,job[job[,4]==last,]);                                #
        print(as.numeric(as.vector(job[1,1:5])));                                     #
    }                                                                                 #
    for(i in 1:dim(lastgen)[1]){                                                      #
        line <- c(initial,additional,as.numeric(as.vector(lastgen[i,])));             #
        cat(line,file=addedto,sep="\t",append=TRUE);                                  #
        cat("\n",file=addedto,append=TRUE);                                           #
        print(c(i,line));                                                             #
    }                                                                                 #
    rm(res);                                                                          #
}                                                                                     #
# concat.last(fullfile="evo$process.txt", addedto="evo.txt",                          #
#             fromlast=0,initial=NA,additional=NA);                                   #
# Below is for use of the ped files ###################################################
clean.last <- function(fullfile,addedto,initial,additional){                          #
    res     <- read.table(file=fullfile,header=FALSE);                                #
    lfm     <- res[res[,3]==0 & res[,4]==1 & res[,8] >= 0,];                          #
    rm(res);                                                                          #
    dat     <- cbind(lfm[,13:14],lfm[,1],lfm[,7],lfm[,9],lfm[,10:12]);                #
    nsim    <- -1;                                                                    #
    for(i in 1:dim(dat)[1]){                                                          #
        line <- c(initial,additional,as.numeric(as.vector(dat[i,])));                 #
        cat(line,file=addedto,sep="\t",append=TRUE);                                  #
        cat("\n",file=addedto,append=TRUE);                                           #
        if(line[5] != nsim){                                                          #
            print(c(i,line));                                                         #
            nsim <- line[5];                                                          #
        }                                                                             #
    }                                                                                 #
}                                                                                     #
#clean.last(fullfile="ped$process.txt",addedto="ped.txt",initial=NA,additional=NA);   #
# Below is for use of the mates files #################################################
clean.mate <- function(fullfile,addedto,initial,additional){                          #
    res     <- read.ragge(file=fullfile);                                             #
    for(i in 1:length(res)){                                                          #
        Pr <- res[[i]][1];                                                            #
        Bi <- res[[i]][2];                                                            #
        Ci <- res[[i]][3];                                                            #
        Py <- res[[i]][5];                                                            #
        Nm <- res[[i]][6];                                                            #
        IM <- res[[i]][8];                                                            #
        if(length(res[[i]]) > 8){                                                     #
            if(res[[i]][7] == res[[i]][9]){ # If a female rejected additional males   # 
                Nm <- 1;                    # and instead remated with her initial    #
                AM <- NA;                                                             #
            }else{                                                                    #
                sumk <- 0;                                                            #
                nmbk <- 0;                                                            #
                EPMs <- (length(res[[i]]) - 8) / 2;                                   #
                for(j in 1:EPMs){                                                     #
                    sumk <- sumk + res[[i]][(2*j)+8];                                 #
                    nmbk <- nmbk + 1;                                                 #
                }                                                                     #
                AM <- sumk / nmbk;                                                    #
            }                                                                         #
        }else{                                                                        #
            AM <- NA;                                                                 #
        }                                                                             #
        line <- c(initial,additional,Pr,Bi,Ci,Py,Nm,IM,AM);                           #
        cat(line,file=addedto,sep="\t",append=TRUE);                                  #
        cat("\n",file=addedto,append=TRUE);                                           #
        if(i %% 10000 == 0){                                                          #
            print(c(i,line));                                                         #
        }                                                                             #
    }                                                                                 #
}                                                                                     #
#clean.mate(fullfile="mat$process.txt",addedto="mat.txt",initial=NA,additional=NA);   #
# =====================================================================================
# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX #
# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX #
# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX #
# =====================================================================================
# NOTE: From here down, the analysis will run using three file: 1) evo.txt, 2) ped.txt, and 3) mat.txt ; each of these files is a concatination of the three correponding files (see above) printed off after Polyinbreeding simulations:
# 1) evo.txt : Key summary statistics of mean allele values, f-coefficients, and correlation between traits from the last generation of a simulation (see above).
# 2) ped.txt : The complete pedigree of individuals from the all simulations in the last two generations (see above).
# 3) mat.xt : All of a females mate choices, and her kinship with each, in the last two generations (see above).
# These files will be used to generate all of the figures in the publication, and in the supporting information. The might not work for *different* parameter values that users might wish to explore, but the code below can be modified to adapt figures to new paramter values.
# =====================================================================================
# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX #
# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX #
# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX #
# =====================================================================================
# Read in the key information, and process it in a way that can be used easily for figures
last.gen <- read.table("evo.txt",header=FALSE);
evo      <- cbind(last.gen[,1:2],last.gen[,4:5],last.gen[,8],last.gen[,7]);
colnames(evo) <- c("Ini","Add","Beta","Cost","P","W");

ped      <- read.table("ped.txt",header=FALSE);
colnames(ped) <- c("Ini","Add","Beta","Cost","Process","Fcoeff","AdOff","gW","gP","Neut");
sims     <- unique(ped[,1:4]);
summ <- NULL;
for(i in 1:dim(sims)[1]){
    use  <- which(ped[,1]==sims[i,1] & 
                  ped[,2]==sims[i,2] & 
                  ped[,3]==sims[i,3] & 
                  ped[,4]==sims[i,4]);
    work   <- ped[use,];
    gPyt   <- work[,9];
    gPyt[gPyt < 0] <- 0;
    mns.P  <- tapply(X=gPyt,    INDEX=work[,5],FUN=mean);
    mns.W  <- tapply(X=work[,8],INDEX=work[,5],FUN=mean);
    mns.k  <- tapply(X=work[,6],INDEX=work[,5],FUN=mean);
    simblk <- matrix(dat=as.numeric(rep(sims[i,],length(mns.P))),ncol=4,byrow=TRUE);
    block  <- cbind(simblk,mns.P,mns.W,mns.k);
    summ   <- rbind(summ,block);
    print(as.numeric(sims[i,]));
}
colnames(summ) <- c("Ini","Add","Beta","Cost","mns.P","mns.W","mns.k");

alle <- NULL;
for(i in 1:dim(sims)[1]){
    use  <- which(evo[,1]==sims[i,1] & 
                  evo[,2]==sims[i,2] & 
                  evo[,3]==sims[i,3] & 
                  evo[,4]==sims[i,4]);
    work   <- evo[use,];
    if(sum(is.na(work[,5])) > 0){ # Should not happen, but removes NA values just in case.
        rmov   <- which(is.na(work[,5]));
        work   <- work[-rmov,];
    }
    CIs.P  <- simpleboot(freqs=work[,5]);
    CIs.W  <- simpleboot(freqs=work[,6]);
    Pvals  <- c(mean(work[,5]),CIs.P);
    Wvals  <- c(mean(work[,6]),CIs.W);
    line   <- c(as.numeric(sims[i,]),Pvals,Wvals);
    alle   <- rbind(alle,line);
    print(as.numeric(line));
}
rownames(alle) <- NULL;
colnames(alle) <- c("Ini","Add","Beta","Cost","Pmn","Pdn","Pup","Wmn","Wdn","Wup");
# The function below sorts the first two columns of a table (needed for later)
bcsort <- function(tab,col1=1,col2=2){
    ord <- order(tab[,col1]);
    tab <- tab[ord,];
    ntb <- NULL;
    for(i in unique(tab[,col1])){
        tmp <- tab[tab[,col1]==i,];
        tor <- order(tmp[,col2]);
        ntb <- rbind(ntb,tmp[tor,]);
    }
    return(ntb);
}
# =====================================================================================
# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX #
# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX #
# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX #
# =====================================================================================
# TODO: BUILD FIGURE 1                                                                #
# ====================================================================================#
setEPS(); # postscript below for final publication?
cairo_ps("gI_vals.eps",family="Arial",height=7,width=7);
par(mar=c(0.2,0.2,0.2,0.2),oma=c(6,6,1,1),lwd=2);
use <- which(summ[,1] == 100 & summ[,2] == 100);
aus <- which(alle[,1] == 100 & alle[,2] == 100);
dat <- as.data.frame(cbind(summ[use,3],summ[use,4],summ[use,6]));
colnames(dat) <- c("bs","cs","vl");
sda <- melt(dat, id = c('bs','cs'));
sda <- sda[,-3];
ale <- alle[aus,];
fourblock <- c(1,1.8,2.6,3.4,4.2);
fullblock <- c(fourblock,fourblock+6,fourblock+12,fourblock+18);
boxplots.dat <- boxplot(as.numeric(value)~bs+cs,data=sda,xlim=c(0,fullblock[20]+2),
                at=fullblock,lwd=2.25,xaxt="n",yaxt="n",
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5,ylim=c(-520,230));
axis(side=2,at=c(-400,-200,0,200),cex.axis=2,cex.lab=2);
for(i in 1:dim(ale)[1]){
    if(ale[i,6] > 0 & ale[i,7] > 0){
        lsh <- fullblock[i] - 0.25;
        rsh <- fullblock[i] + 0.25;
        polygon(x=c(lsh,rsh,rsh,lsh),y=c(-10,-10,350,350),border=NA,col="grey70");
    }
}
boxplots.dat <- boxplot(as.numeric(value)~bs+cs,data=sda,xlim=c(0,fullblock[20]+2),
                at=fullblock,lwd=2.25,add=TRUE,xaxt="n",yaxt="n",
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5,ylim=c(0,240));
abline(h=0,lty="dotted",lwd=0.8);
axis(side=1,at=c(fourblock[3],fourblock[3]+6,fourblock[3]+12,fourblock[3]+18),
     cex.axis=2,cex.lab=2,labels=c("0.00","0.01","0.02","0.03"));
polygon(x=c(fullblock[11]-1,fullblock[20]+2.5,fullblock[20]+2.5,fullblock[11]-1),
        y=c(240-50,240-50,240,240),border=NA,col="black");
polygon(x=c(fullblock[11]-1,fullblock[20]+2.5,fullblock[20]+2.5,fullblock[11]-1),
        y=c(-538,-538,-538+50,-538+50),border=NA,col="black");
text(x=fullblock[11]-1,y=210,cex=1.5,labels="Inbreeding preference",pos=4,col="white");
text(x=fullblock[11]-1,y=-518,cex=1.5,labels="Inbreeding avoidance",pos=4,col="white");
segments(x0=fullblock[20]+2.5, y0=-538, x1 = fullblock[20]+2.5, y1 = 240, lwd=3);
segments(x0=fullblock[11]-1, y0=240, x1 = fullblock[20]+2.5, y1 = 240, lwd=3);
segments(x0=fullblock[11]-1, y0=240-50, x1 = fullblock[20]+2.5, y1 = 240-50, lwd=3);
segments(x0=fullblock[11]-1, y0=240-50, x1 = fullblock[11]-1, y1 = 240, lwd=3);
segments(x0=fullblock[11]-1, y0=-538, x1 = fullblock[20]+2.5, y1 = -538, lwd=3);
segments(x0=fullblock[11]-1, y0=-538+50, x1 = fullblock[20]+2.5, y1 = -538+50, lwd=3);
segments(x0=fullblock[11]-1, y0=-538+50, x1 = fullblock[11]-1, y1 = -538, lwd=3);

mtext(expression(paste("Cost of polyandry (",c[P],")")),
	outer=TRUE,side=1,line=3.75,cex=2);
mtext(expression(paste("Mean inbreeding strategy phenotype value (",I[p],")")),
	outer=TRUE,side=2,line=3.25,cex=1.5);

dev.off();
# ====================================================================================#

# =====================================================================================
# TODO: BUILD FIGURE 2                                                                #
# ====================================================================================#
# The chunk of code below just separates Beta values evenly so the histograms are clear
even.Beta <- function(vector){
    length.vector <- length(vector);
    for(i in 1:length.vector){
        check <- vector[i];
        if(check == 0){
            vector[i] <- 1;
        }
        if(check == 0.2){
            vector[i] <- 2;
        }
        if(check == 1){
            vector[i] <- 3; 
        }
        if(check == 2){
            vector[i] <- 4;
        }
        if(check == 5){
            vector[i] <- 5;
        }
    }
    return(vector);
}
# -------------------------------------------------------------------------------------
em       <- last.gen; # Using the original evo.txt file now.
if(sum(is.na(em[,7])) > 0){ # Should not happen, but if NAs 
    rm       <- which(is.na(em[,7]));
    em       <- em[-rm,];
}
# -------------------------------------------------------------------------------------
setEPS(); # postscript below for final publication?
cairo_ps("Pvalues.eps",family="Arial",height=7,width=7);
par(mfrow=c(3,3),mar=c(0.2,0.2,0.2,0.2),oma=c(6,6,1,1));
# ----------------------- Panel `a':
use <- which(em[,1]==100 & em[,2]==100);
dat <- em[use,];
nBn <- even.Beta(dat[,4]);
fourblock <- c(1,1.8,2.6,3.4,4.2);
fullblock <- c(fourblock,fourblock+6,fourblock+12,fourblock+18);

sdd <- as.data.frame(cbind(dat[,4],dat[,5],dat[,8]));
colnames(sdd) <- c("bs","cs","vl");
sda <- melt(sdd, id = c('bs','cs'));
sda <- sda[,-3];
boxplots.dat <- boxplot(as.numeric(value)~bs+cs,data=sda,type="n",
                xlim=c(1,max(fullblock)),yaxt="n",xaxt="n",
                at=fullblock,lwd=1.5,cex.lab=1.75,ylim=c(-18,16),
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5);
TST <- NULL;
for(j in unique(dat[,5])){
    for(i in unique(dat[,4])){
        sub <- which(dat[,4]==i & dat[,5]==j);
        cis <- simpleboot(dat[sub,8]);
        ttt <- cis[1] > 0 & cis[2] > 0;
        TST <- rbind(TST,c(j,i,cis,ttt));
    }
}
TST <- bcsort(TST);
for(i in 1:dim(TST)[1]){
    if(TST[i,5] == 1){
        lsh <- fullblock[i] - 0.25;
        rsh <- fullblock[i] + 0.25; 
        polygon(x=c(lsh,rsh,rsh,lsh),y=c(-40,-40,20,20),border=NA,col="grey70");
    }
}
boxplots.dat <- boxplot(as.numeric(value)~bs+cs,data=sda,
                xlim=c(1,max(fullblock)),yaxt="n",xaxt="n",
                at=fullblock,lwd=1.5,cex.lab=1.75,ylim=c(-18,16),
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5,add=TRUE);
abline(h=0,lwd=0.8,lty="dotted");
text(x=0,y=14,cex=1.75,labels=expression(paste("(A) ",S['100,100'])),pos=4);
axis(side=2,at=c(-20,-10,0,10),cex.axis=1.5);
# ----------------------- Panel `b':
use <- which(em[,1]==100 & em[,2]==10);
dat <- em[use,];
nBn <- even.Beta(dat[,4]);
fourblock <- c(1,1.8,2.6,3.4,4.2);
fullblock <- c(fourblock,fourblock+6,fourblock+12,fourblock+18);

sdd <- as.data.frame(cbind(dat[,4],dat[,5],dat[,8]));
colnames(sdd) <- c("bs","cs","vl");
sda <- melt(sdd, id = c('bs','cs'));
sda <- sda[,-3];
boxplots.dat <- boxplot(as.numeric(value)~bs+cs,data=sda,type="n",
                xlim=c(1,max(fullblock)),yaxt="n",xaxt="n",
                at=fullblock,lwd=1.5,cex.lab=1.75,ylim=c(-18,16),
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5);
TST <- NULL;
for(j in unique(dat[,5])){
    for(i in unique(dat[,4])){
        sub <- which(dat[,4]==i & dat[,5]==j);
        cis <- simpleboot(dat[sub,8]);
        ttt <- cis[1] > 0 & cis[2] > 0;
        TST <- rbind(TST,c(j,i,cis,ttt));
    }
}
TST <- bcsort(TST);
for(i in 1:dim(TST)[1]){
    if(TST[i,5] == 1){
        lsh <- fullblock[i] - 0.25;
        rsh <- fullblock[i] + 0.25; 
        polygon(x=c(lsh,rsh,rsh,lsh),y=c(-40,-40,20,20),border=NA,col="grey70");
    }
}
boxplots.dat <- boxplot(as.numeric(value)~bs+cs,data=sda,
                xlim=c(1,max(fullblock)),yaxt="n",xaxt="n",
                at=fullblock,lwd=1.5,cex.lab=1.75,ylim=c(-18,16),
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5,add=TRUE);
abline(h=0,lwd=0.8,lty="dotted");
text(x=0,y=14,cex=1.75,labels=expression(paste("(B) ",S['100,10'])),pos=4);
# ----------------------- Panel `c':
use <- which(em[,1]==100 & em[,2]==2);
dat <- em[use,];
nBn <- even.Beta(dat[,4]);
fourblock <- c(1,1.8,2.6,3.4,4.2);
fullblock <- c(fourblock,fourblock+6,fourblock+12,fourblock+18);

sdd <- as.data.frame(cbind(dat[,4],dat[,5],dat[,8]));
colnames(sdd) <- c("bs","cs","vl");
sda <- melt(sdd, id = c('bs','cs'));
sda <- sda[,-3];
boxplots.dat <- boxplot(as.numeric(value)~bs+cs,data=sda,type="n",
                xlim=c(1,max(fullblock)),yaxt="n",xaxt="n",
                at=fullblock,lwd=1.5,cex.lab=1.75,ylim=c(-18,16),
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5);
TST <- NULL;
for(j in unique(dat[,5])){
    for(i in unique(dat[,4])){
        sub <- which(dat[,4]==i & dat[,5]==j);
        cis <- simpleboot(dat[sub,8]);
        ttt <- cis[1] > 0 & cis[2] > 0;
        TST <- rbind(TST,c(j,i,cis,ttt));
    }
}
TST <- bcsort(TST);
for(i in 1:dim(TST)[1]){
    if(TST[i,5] == 1){
        lsh <- fullblock[i] - 0.25;
        rsh <- fullblock[i] + 0.25; 
        polygon(x=c(lsh,rsh,rsh,lsh),y=c(-40,-40,20,20),border=NA,col="grey70");
    }
}
boxplots.dat <- boxplot(as.numeric(value)~bs+cs,data=sda,
                xlim=c(1,max(fullblock)),yaxt="n",xaxt="n",
                at=fullblock,lwd=1.5,cex.lab=1.75,ylim=c(-18,16),
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5,add=TRUE);
abline(h=0,lwd=0.8,lty="dotted");
text(x=0,y=14,cex=1.75,labels=expression(paste("(C) ",S['100,2'])),pos=4);
# ----------------------- Panel `d':
use <- which(em[,1]==10 & em[,2]==100);
dat <- em[use,];
nBn <- even.Beta(dat[,4]);
fourblock <- c(1,1.8,2.6,3.4,4.2);
fullblock <- c(fourblock,fourblock+6,fourblock+12,fourblock+18);

sdd <- as.data.frame(cbind(dat[,4],dat[,5],dat[,8]));
colnames(sdd) <- c("bs","cs","vl");
sda <- melt(sdd, id = c('bs','cs'));
sda <- sda[,-3];
boxplots.dat <- boxplot(as.numeric(value)~bs+cs,data=sda,type="n",
                xlim=c(1,max(fullblock)),yaxt="n",xaxt="n",
                at=fullblock,lwd=1.5,cex.lab=1.75,ylim=c(-18,16),
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5);
TST <- NULL;
for(j in unique(dat[,5])){
    for(i in unique(dat[,4])){
        sub <- which(dat[,4]==i & dat[,5]==j);
        cis <- simpleboot(dat[sub,8]);
        ttt <- cis[1] > 0 & cis[2] > 0;
        TST <- rbind(TST,c(j,i,cis,ttt));
    }
}
TST <- bcsort(TST);
for(i in 1:dim(TST)[1]){
    if(TST[i,5] == 1){
        lsh <- fullblock[i] - 0.25;
        rsh <- fullblock[i] + 0.25; 
        polygon(x=c(lsh,rsh,rsh,lsh),y=c(-40,-40,20,20),border=NA,col="grey70");
    }
}
boxplots.dat <- boxplot(as.numeric(value)~bs+cs,data=sda,
                xlim=c(1,max(fullblock)),yaxt="n",xaxt="n",
                at=fullblock,lwd=1.5,cex.lab=1.75,ylim=c(-18,16),
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5,add=TRUE);
abline(h=0,lwd=0.8,lty="dotted");
text(x=0,y=14,cex=1.75,labels=expression(paste("(D) ",S['10,100'])),pos=4);
axis(side=2,at=c(-20,-10,0,10),cex.axis=1.5);
# ----------------------- Panel `e':
use <- which(em[,1]==10 & em[,2]==10);
dat <- em[use,];
nBn <- even.Beta(dat[,4]);
fourblock <- c(1,1.8,2.6,3.4,4.2);
fullblock <- c(fourblock,fourblock+6,fourblock+12,fourblock+18);

sdd <- as.data.frame(cbind(dat[,4],dat[,5],dat[,8]));
colnames(sdd) <- c("bs","cs","vl");
sda <- melt(sdd, id = c('bs','cs'));
sda <- sda[,-3];
boxplots.dat <- boxplot(as.numeric(value)~bs+cs,data=sda,type="n",
                xlim=c(1,max(fullblock)),yaxt="n",xaxt="n",
                at=fullblock,lwd=1.5,cex.lab=1.75,ylim=c(-18,16),
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5);
TST <- NULL;
for(j in unique(dat[,5])){
    for(i in unique(dat[,4])){
        sub <- which(dat[,4]==i & dat[,5]==j);
        cis <- simpleboot(dat[sub,8]);
        ttt <- cis[1] > 0 & cis[2] > 0;
        TST <- rbind(TST,c(j,i,cis,ttt));
    }
}
TST <- bcsort(TST);
for(i in 1:dim(TST)[1]){
    if(TST[i,5] == 1){
        lsh <- fullblock[i] - 0.25;
        rsh <- fullblock[i] + 0.25; 
        polygon(x=c(lsh,rsh,rsh,lsh),y=c(-40,-40,20,20),border=NA,col="grey70");
    }
}
boxplots.dat <- boxplot(as.numeric(value)~bs+cs,data=sda,
                xlim=c(1,max(fullblock)),yaxt="n",xaxt="n",
                at=fullblock,lwd=1.5,cex.lab=1.75,ylim=c(-18,16),
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5,add=TRUE);
abline(h=0,lwd=0.8,lty="dotted");
text(x=0,y=14,cex=1.75,labels=expression(paste("(E) ",S['10,10'])),pos=4);
# ----------------------- Panel `f':
use <- which(em[,1]==10 & em[,2]==2);
dat <- em[use,];
nBn <- even.Beta(dat[,4]);
fourblock <- c(1,1.8,2.6,3.4,4.2);
fullblock <- c(fourblock,fourblock+6,fourblock+12,fourblock+18);

sdd <- as.data.frame(cbind(dat[,4],dat[,5],dat[,8]));
colnames(sdd) <- c("bs","cs","vl");
sda <- melt(sdd, id = c('bs','cs'));
sda <- sda[,-3];
boxplots.dat <- boxplot(as.numeric(value)~bs+cs,data=sda,type="n",
                xlim=c(1,max(fullblock)),yaxt="n",xaxt="n",
                at=fullblock,lwd=1.5,cex.lab=1.75,ylim=c(-18,16),
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5);
TST <- NULL;
for(j in unique(dat[,5])){
    for(i in unique(dat[,4])){
        sub <- which(dat[,4]==i & dat[,5]==j);
        cis <- simpleboot(dat[sub,8]);
        ttt <- cis[1] > 0 & cis[2] > 0;
        TST <- rbind(TST,c(j,i,cis,ttt));
    }
}
TST <- bcsort(TST);
for(i in 1:dim(TST)[1]){
    if(TST[i,5] == 1){
        lsh <- fullblock[i] - 0.25;
        rsh <- fullblock[i] + 0.25; 
        polygon(x=c(lsh,rsh,rsh,lsh),y=c(-40,-40,20,20),border=NA,col="grey70");
    }
}
boxplots.dat <- boxplot(as.numeric(value)~bs+cs,data=sda,
                xlim=c(1,max(fullblock)),yaxt="n",xaxt="n",
                at=fullblock,lwd=1.5,cex.lab=1.75,ylim=c(-18,16),
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5,add=TRUE);
abline(h=0,lwd=0.8,lty="dotted");
text(x=0,y=14,cex=1.75,labels=expression(paste("(F) ",S['10,2'])),pos=4);
# ----------------------- Panel `g':
use <- which(em[,1]==2 & em[,2]==100);
dat <- em[use,];
nBn <- even.Beta(dat[,4]);
fourblock <- c(1,1.8,2.6,3.4,4.2);
fullblock <- c(fourblock,fourblock+6,fourblock+12,fourblock+18);

sdd <- as.data.frame(cbind(dat[,4],dat[,5],dat[,8]));
colnames(sdd) <- c("bs","cs","vl");
sda <- melt(sdd, id = c('bs','cs'));
sda <- sda[,-3];
boxplots.dat <- boxplot(as.numeric(value)~bs+cs,data=sda,type="n",
                xlim=c(1,max(fullblock)),yaxt="n",xaxt="n",
                at=fullblock,lwd=1.5,cex.lab=1.75,ylim=c(-18,16),
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5);
TST <- NULL;
for(j in unique(dat[,5])){
    for(i in unique(dat[,4])){
        sub <- which(dat[,4]==i & dat[,5]==j);
        cis <- simpleboot(dat[sub,8]);
        ttt <- cis[1] > 0 & cis[2] > 0;
        TST <- rbind(TST,c(j,i,cis,ttt));
    }
}
TST <- bcsort(TST);
for(i in 1:dim(TST)[1]){
    if(TST[i,5] == 1){
        lsh <- fullblock[i] - 0.25;
        rsh <- fullblock[i] + 0.25; 
        polygon(x=c(lsh,rsh,rsh,lsh),y=c(-40,-40,20,20),border=NA,col="grey70");
    }
}
boxplots.dat <- boxplot(as.numeric(value)~bs+cs,data=sda,
                xlim=c(1,max(fullblock)),yaxt="n",xaxt="n",
                at=fullblock,lwd=1.5,cex.lab=1.75,ylim=c(-18,16),
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5,add=TRUE);
abline(h=0,lwd=0.8,lty="dotted");
text(x=0,y=14,cex=1.75,labels=expression(paste("(G) ",S['2,100'])),pos=4);
axis(side=2,at=c(-20,-10,0,10),cex.axis=1.5);
axis(side=1,at=c(fourblock[3],fourblock[3]+6,fourblock[3]+12,fourblock[3]+18),
     cex.axis=1.5,cex.lab=1.5,labels=c("0.00","0.01","0.02","0.03"));
# ----------------------- Panel `h':
use <- which(em[,1]==2 & em[,2]==10);
dat <- em[use,];
nBn <- even.Beta(dat[,4]);
fourblock <- c(1,1.8,2.6,3.4,4.2);
fullblock <- c(fourblock,fourblock+6,fourblock+12,fourblock+18);

sdd <- as.data.frame(cbind(dat[,4],dat[,5],dat[,8]));
colnames(sdd) <- c("bs","cs","vl");
sda <- melt(sdd, id = c('bs','cs'));
sda <- sda[,-3];
boxplots.dat <- boxplot(as.numeric(value)~bs+cs,data=sda,type="n",
                xlim=c(1,max(fullblock)),yaxt="n",xaxt="n",
                at=fullblock,lwd=1.5,cex.lab=1.75,ylim=c(-18,16),
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5);
TST <- NULL;
for(j in unique(dat[,5])){
    for(i in unique(dat[,4])){
        sub <- which(dat[,4]==i & dat[,5]==j);
        cis <- simpleboot(dat[sub,8]);
        ttt <- cis[1] > 0 & cis[2] > 0;
        TST <- rbind(TST,c(j,i,cis,ttt));
    }
}
TST <- bcsort(TST);
for(i in 1:dim(TST)[1]){
    if(TST[i,5] == 1){
        lsh <- fullblock[i] - 0.25;
        rsh <- fullblock[i] + 0.25; 
        polygon(x=c(lsh,rsh,rsh,lsh),y=c(-40,-40,20,20),border=NA,col="grey70");
    }
}
boxplots.dat <- boxplot(as.numeric(value)~bs+cs,data=sda,
                xlim=c(1,max(fullblock)),yaxt="n",xaxt="n",
                at=fullblock,lwd=1.5,cex.lab=1.75,ylim=c(-18,16),
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5,add=TRUE);
abline(h=0,lwd=0.8,lty="dotted");
text(x=0,y=14,cex=1.75,labels=expression(paste("(H) ",S['2,10'])),pos=4);
axis(side=1,at=c(fourblock[3],fourblock[3]+6,fourblock[3]+12,fourblock[3]+18),
     cex.axis=1.5,cex.lab=1.5,labels=c("0.00","0.01","0.02","0.03"));
# ----------------------- Panel `i':
use <- which(em[,1]==2 & em[,2]==2);
dat <- em[use,];
nBn <- even.Beta(dat[,4]);
fourblock <- c(1,1.8,2.6,3.4,4.2);
fullblock <- c(fourblock,fourblock+6,fourblock+12,fourblock+18);

sdd <- as.data.frame(cbind(dat[,4],dat[,5],dat[,8]));
colnames(sdd) <- c("bs","cs","vl");
sda <- melt(sdd, id = c('bs','cs'));
sda <- sda[,-3];
boxplots.dat <- boxplot(as.numeric(value)~bs+cs,data=sda,type="n",
                xlim=c(1,max(fullblock)),yaxt="n",xaxt="n",
                at=fullblock,lwd=1.5,cex.lab=1.75,ylim=c(-18,16),
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5);
TST <- NULL;
for(j in unique(dat[,5])){
    for(i in unique(dat[,4])){
        sub <- which(dat[,4]==i & dat[,5]==j);
        cis <- simpleboot(dat[sub,8]);
        ttt <- cis[1] > 0 & cis[2] > 0;
        TST <- rbind(TST,c(j,i,cis,ttt));
    }
}
TST <- bcsort(TST);
for(i in 1:dim(TST)[1]){
    if(TST[i,5] == 1){
        lsh <- fullblock[i] - 0.25;
        rsh <- fullblock[i] + 0.25; 
        #polygon(x=c(lsh,rsh,rsh,lsh),y=c(-40,-40,20,20),border=NA,col="grey70");
    }
}
boxplots.dat <- boxplot(as.numeric(value)~bs+cs,data=sda,
                xlim=c(1,max(fullblock)),yaxt="n",xaxt="n",
                at=fullblock,lwd=1.5,cex.lab=1.75,ylim=c(-18,16),
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5,add=TRUE);
abline(h=0,lwd=0.8,lty="dotted");
text(x=0,y=14,cex=1.75,labels=expression(paste("(I) ",S['2,2'])),pos=4);
axis(side=1,at=c(fourblock[3],fourblock[3]+6,fourblock[3]+12,fourblock[3]+18),
     cex.axis=1.5,cex.lab=1.5,labels=c("0.00","0.01","0.02","0.03"));

mtext(expression(paste("Cost of polyandry (",c[P],")")),
	outer=TRUE,side=1,line=3.75,cex=1.5);

mtext(expression(paste("Mean polyandry allele value (",P[a],")")),
	outer=TRUE,side=2,line=3.25,cex=1.5);

dev.off();
# ====================================================================================#

# =====================================================================================
# TODO: BUILD FIGURE 3                                                                #
# ====================================================================================#
setEPS(); # postscript below for final publication?
cairo_ps("FullFact_gP.eps",family="Arial",height=7,width=7);
par(mfrow=c(3,3),mar=c(0.2,0.2,0.2,0.2),oma=c(6,6,1,1));
# ----------------------- Panel `a':
use <- which(summ[,1] == 100 & summ[,2] == 100);
aus <- which(alle[,1] == 100 & alle[,2] == 100);
dat <- as.data.frame(cbind(summ[use,3],summ[use,4],summ[use,5]));
dat <- bcsort(dat,col1=2,col2=1);
colnames(dat) <- c("bs","cs","vl");
sda <- melt(dat, id = c('bs','cs'));
sda <- sda[,-3];
ale <- alle[aus,];
ale <- bcsort(ale,col1=4,col2=3);
fourblock <- c(1,1.8,2.6,3.4,4.2);
fullblock <- c(fourblock,fourblock+6,fourblock+12,fourblock+18);
boxplots.dat <- boxplot(as.numeric(value)~bs+cs,data=sda,
                at=fullblock,lwd=1.5,xaxt="n",yaxt="n",
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5,ylim=c(0,240));
axis(side=2,at=c(0,100,200),cex.axis=2,cex.lab=2);
for(i in 1:dim(ale)[1]){
    if(ale[i,6] > 0 & ale[i,7] > 0){
        lsh <- fullblock[i] - 0.25;
        rsh <- fullblock[i] + 0.25;
        polygon(x=c(lsh,rsh,rsh,lsh),y=c(-10,-10,350,350),border=NA,col="grey70");
    }
}
boxplots.dat <- boxplot(as.numeric(value)~bs+cs,data=sda,
                at=fullblock,lwd=1.5,add=TRUE,xaxt="n",yaxt="n",
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5,ylim=c(0,240));
text(x=9,y=220,cex=2,labels=expression(paste("(A) ",S['100,100'])),pos=4);
# ----------------------- Panel `b':
use <- which(summ[,1] == 100 & summ[,2] == 10);
aus <- which(alle[,1] == 100 & alle[,2] == 10);
dat <- as.data.frame(cbind(summ[use,3],summ[use,4],summ[use,5]));
dat <- bcsort(dat,col1=2,col2=1);
colnames(dat) <- c("bs","cs","vl");
sda <- melt(dat, id = c('bs','cs'));
sda <- sda[,-3];
ale <- alle[aus,];
ale <- bcsort(ale,col1=4,col2=3);
fourblock <- c(1,1.8,2.6,3.4,4.2);
fullblock <- c(fourblock,fourblock+6,fourblock+12,fourblock+18);
boxplots.dat <- boxplot(as.numeric(value)~bs+cs,data=sda,
                at=fullblock,lwd=1.5,xaxt="n",yaxt="n",
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5,ylim=c(0,240));
for(i in 1:dim(ale)[1]){
    if(ale[i,6] > 0 & ale[i,7] > 0){
        lsh <- fullblock[i] - 0.25;
        rsh <- fullblock[i] + 0.25;
        polygon(x=c(lsh,rsh,rsh,lsh),y=c(-10,-10,350,350),border=NA,col="grey70");
    }
}
boxplots.dat <- boxplot(as.numeric(value)~bs+cs,data=sda,
                at=fullblock,lwd=1.5,add=TRUE,xaxt="n",yaxt="n",
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5,ylim=c(0,240));
text(x=9,y=220,cex=2,labels=expression(paste("(B) ",S['100,10'])),pos=4);
# ----------------------- Panel `c':
use <- which(summ[,1] == 100 & summ[,2] == 2);
aus <- which(alle[,1] == 100 & alle[,2] == 2);
dat <- as.data.frame(cbind(summ[use,3],summ[use,4],summ[use,5]));
dat <- bcsort(dat,col1=2,col2=1);
colnames(dat) <- c("bs","cs","vl");
sda <- melt(dat, id = c('bs','cs'));
sda <- sda[,-3];
ale <- alle[aus,];
ale <- bcsort(ale,col1=4,col2=3);
fourblock <- c(1,1.8,2.6,3.4,4.2);
fullblock <- c(fourblock,fourblock+6,fourblock+12,fourblock+18);
boxplots.dat <- boxplot(as.numeric(value)~bs+cs,data=sda,
                at=fullblock,lwd=1.5,xaxt="n",yaxt="n",
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5,ylim=c(0,240));
for(i in 1:dim(ale)[1]){
    if(ale[i,6] > 0 & ale[i,7] > 0){
        lsh <- fullblock[i] - 0.25;
        rsh <- fullblock[i] + 0.25;
        polygon(x=c(lsh,rsh,rsh,lsh),y=c(-10,-10,350,350),border=NA,col="grey70");
    }
}
boxplots.dat <- boxplot(as.numeric(value)~bs+cs,data=sda,
                at=fullblock,lwd=1.5,add=TRUE,xaxt="n",yaxt="n",
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5,ylim=c(0,240));
text(x=9,y=220,cex=2,labels=expression(paste("(C) ",S['100,2'])),pos=4);
# ----------------------- Panel `d':
use <- which(summ[,1] == 10 & summ[,2] == 100);
aus <- which(alle[,1] == 10 & alle[,2] == 100);
dat <- as.data.frame(cbind(summ[use,3],summ[use,4],summ[use,5]));
dat <- bcsort(dat,col1=2,col2=1);
colnames(dat) <- c("bs","cs","vl");
sda <- melt(dat, id = c('bs','cs'));
sda <- sda[,-3];
ale <- alle[aus,];
ale <- bcsort(ale,col1=4,col2=3);
fourblock <- c(1,1.8,2.6,3.4,4.2);
fullblock <- c(fourblock,fourblock+6,fourblock+12,fourblock+18);
boxplots.dat <- boxplot(as.numeric(value)~bs+cs,data=sda,
                at=fullblock,lwd=1.5,xaxt="n",yaxt="n",
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5,ylim=c(0,240));
axis(side=2,at=c(0,100,200),cex.axis=2,cex.lab=2);
for(i in 1:dim(ale)[1]){
    if(ale[i,6] > 0 & ale[i,7] > 0){
        lsh <- fullblock[i] - 0.25;
        rsh <- fullblock[i] + 0.25;
        polygon(x=c(lsh,rsh,rsh,lsh),y=c(-10,-10,350,350),border=NA,col="grey70");
    }
}
boxplots.dat <- boxplot(as.numeric(value)~bs+cs,data=sda,
                at=fullblock,lwd=1.5,add=TRUE,xaxt="n",yaxt="n",
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5,ylim=c(0,240));
text(x=9,y=220,cex=2,labels=expression(paste("(D) ",S['10,100'])),pos=4);
# ----------------------- Panel `e':
use <- which(summ[,1] == 10 & summ[,2] == 10);
aus <- which(alle[,1] == 10 & alle[,2] == 10);
dat <- as.data.frame(cbind(summ[use,3],summ[use,4],summ[use,5]));
dat <- bcsort(dat,col1=2,col2=1);
colnames(dat) <- c("bs","cs","vl");
sda <- melt(dat, id = c('bs','cs'));
sda <- sda[,-3];
ale <- alle[aus,];
ale <- bcsort(ale,col1=4,col2=3);
fourblock <- c(1,1.8,2.6,3.4,4.2);
fullblock <- c(fourblock,fourblock+6,fourblock+12,fourblock+18);
boxplots.dat <- boxplot(as.numeric(value)~bs+cs,data=sda,
                at=fullblock,lwd=1.5,xaxt="n",yaxt="n",
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5,ylim=c(0,240));
for(i in 1:dim(ale)[1]){
    if(ale[i,6] > 0 & ale[i,7] > 0){
        lsh <- fullblock[i] - 0.25;
        rsh <- fullblock[i] + 0.25;
        polygon(x=c(lsh,rsh,rsh,lsh),y=c(-10,-10,350,350),border=NA,col="grey70");
    }
}
boxplots.dat <- boxplot(as.numeric(value)~bs+cs,data=sda,
                at=fullblock,lwd=1.5,add=TRUE,xaxt="n",yaxt="n",
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5,ylim=c(0,240));
text(x=9,y=220,cex=2,labels=expression(paste("(E) ",S['10,10'])),pos=4);
# ----------------------- Panel `f':
use <- which(summ[,1] == 10 & summ[,2] == 2);
aus <- which(alle[,1] == 10 & alle[,2] == 2);
dat <- as.data.frame(cbind(summ[use,3],summ[use,4],summ[use,5]));
dat <- bcsort(dat,col1=2,col2=1);
colnames(dat) <- c("bs","cs","vl");
sda <- melt(dat, id = c('bs','cs'));
sda <- sda[,-3];
ale <- alle[aus,];
ale <- bcsort(ale,col1=4,col2=3);
fourblock <- c(1,1.8,2.6,3.4,4.2);
fullblock <- c(fourblock,fourblock+6,fourblock+12,fourblock+18);
boxplots.dat <- boxplot(as.numeric(value)~bs+cs,data=sda,
                at=fullblock,lwd=1.5,xaxt="n",yaxt="n",
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5,ylim=c(0,240));
for(i in 1:dim(ale)[1]){
    if(ale[i,6] > 0 & ale[i,7] > 0){
        lsh <- fullblock[i] - 0.25;
        rsh <- fullblock[i] + 0.25;
        polygon(x=c(lsh,rsh,rsh,lsh),y=c(-10,-10,350,350),border=NA,col="grey70");
    }
}
boxplots.dat <- boxplot(as.numeric(value)~bs+cs,data=sda,
                at=fullblock,lwd=1.5,add=TRUE,xaxt="n",yaxt="n",
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5,ylim=c(0,240));
text(x=9,y=220,cex=2,labels=expression(paste("(F) ",S['10,2'])),pos=4);
# ----------------------- Panel `g':
use <- which(summ[,1] == 2 & summ[,2] == 100);
aus <- which(alle[,1] == 2 & alle[,2] == 100);
dat <- as.data.frame(cbind(summ[use,3],summ[use,4],summ[use,5]));
dat <- bcsort(dat,col1=2,col2=1);
colnames(dat) <- c("bs","cs","vl");
sda <- melt(dat, id = c('bs','cs'));
sda <- sda[,-3];
ale <- alle[aus,];
ale <- bcsort(ale,col1=4,col2=3);
fourblock <- c(1,1.8,2.6,3.4,4.2);
fullblock <- c(fourblock,fourblock+6,fourblock+12,fourblock+18);
boxplots.dat <- boxplot(as.numeric(value)~bs+cs,data=sda,
                at=fullblock,lwd=1.5,xaxt="n",yaxt="n",
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5,ylim=c(0,240));
axis(side=2,at=c(0,100,200),cex.axis=2,cex.lab=2);
axis(side=1,at=c(fourblock[3],fourblock[3]+6,fourblock[3]+12,fourblock[3]+18),
     cex.axis=1.5,cex.lab=1.5,labels=c("0.00","0.01","0.02","0.03"));
for(i in 1:dim(ale)[1]){
    if(ale[i,6] > 0 & ale[i,7] > 0){
        lsh <- fullblock[i] - 0.25;
        rsh <- fullblock[i] + 0.25;
        polygon(x=c(lsh,rsh,rsh,lsh),y=c(-10,-10,350,350),border=NA,col="grey70");
    }
}
boxplots.dat <- boxplot(as.numeric(value)~bs+cs,data=sda,
                at=fullblock,lwd=1.5,add=TRUE,xaxt="n",yaxt="n",
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5,ylim=c(0,240));
text(x=9,y=220,cex=2,labels=expression(paste("(G) ",S['2,100'])),pos=4);
# ----------------------- Panel `h':
use <- which(summ[,1] == 2 & summ[,2] == 10);
aus <- which(alle[,1] == 2 & alle[,2] == 10);
dat <- as.data.frame(cbind(summ[use,3],summ[use,4],summ[use,5]));
dat <- bcsort(dat,col1=2,col2=1);
colnames(dat) <- c("bs","cs","vl");
sda <- melt(dat, id = c('bs','cs'));
sda <- sda[,-3];
ale <- alle[aus,];
ale <- bcsort(ale,col1=4,col2=3);
fourblock <- c(1,1.8,2.6,3.4,4.2);
fullblock <- c(fourblock,fourblock+6,fourblock+12,fourblock+18);
boxplots.dat <- boxplot(as.numeric(value)~bs+cs,data=sda,
                at=fullblock,lwd=1.5,xaxt="n",yaxt="n",
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5,ylim=c(0,240));
axis(side=1,at=c(fourblock[3],fourblock[3]+6,fourblock[3]+12,fourblock[3]+18),
     cex.axis=1.5,cex.lab=1.5,labels=c("0.00","0.01","0.02","0.03"));
for(i in 1:dim(ale)[1]){
    if(ale[i,6] > 0 & ale[i,7] > 0){
        lsh <- fullblock[i] - 0.25;
        rsh <- fullblock[i] + 0.25;
        polygon(x=c(lsh,rsh,rsh,lsh),y=c(-10,-10,350,350),border=NA,col="grey70");
    }
}
boxplots.dat <- boxplot(as.numeric(value)~bs+cs,data=sda,
                at=fullblock,lwd=1.5,add=TRUE,xaxt="n",yaxt="n",
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5,ylim=c(0,240),pos=4);
text(x=9,y=220,cex=2,labels=expression(paste("(H) ",S['2,10'])),pos=4);
# ----------------------- Panel `i':
use <- which(summ[,1] == 2 & summ[,2] == 2);
aus <- which(alle[,1] == 2 & alle[,2] == 2);
dat <- as.data.frame(cbind(summ[use,3],summ[use,4],summ[use,5]));
dat <- bcsort(dat,col1=2,col2=1);
colnames(dat) <- c("bs","cs","vl");
sda <- melt(dat, id = c('bs','cs'));
sda <- sda[,-3];
ale <- alle[aus,];
ale <- bcsort(ale,col1=4,col2=3);
fourblock <- c(1,1.8,2.6,3.4,4.2);
fullblock <- c(fourblock,fourblock+6,fourblock+12,fourblock+18);
boxplots.dat <- boxplot(as.numeric(value)~bs+cs,data=sda,
                at=fullblock,lwd=1.5,xaxt="n",yaxt="n",
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5,ylim=c(0,240));
axis(side=1,at=c(fourblock[3],fourblock[3]+6,fourblock[3]+12,fourblock[3]+18),
     cex.axis=1.5,cex.lab=1.5,labels=c("0.00","0.01","0.02","0.03"));
for(i in 1:dim(ale)[1]){
    if(ale[i,6] > 0 & ale[i,7] > 0){
        lsh <- fullblock[i] - 0.25;
        rsh <- fullblock[i] + 0.25;
        #polygon(x=c(lsh,rsh,rsh,lsh),y=c(-10,-10,350,350),border=NA,col="grey70");
    }
}
boxplots.dat <- boxplot(as.numeric(value)~bs+cs,data=sda,
                at=fullblock,lwd=1.5,add=TRUE,xaxt="n",yaxt="n",
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5,ylim=c(0,240));
text(x=9,y=220,cex=2,labels=expression(paste("(I) ",S['2,2'])),pos=4);

mtext(expression(paste("Cost of polyandry (",c[P],")")),
	outer=TRUE,side=1,line=3.75,cex=1.5);

mtext(expression(paste("Mean polyandry phenotype value (",P[p],")")),
	outer=TRUE,side=2,line=3.25,cex=1.5);

dev.off();
# ====================================================================================#

# =====================================================================================
# TODO: BUILD FIGURE 4                                                                #
# ====================================================================================#
setEPS(); # postscript below for final publication?
cairo_ps("Pvalues_Nmin1.eps",family="Arial",height=10,width=7);
par(mfrow=c(2,1),mar=c(1,5,1,1),lwd=2);
# ----------------------- Panel `A':
use <- which(em[,1]==-1 & em[,2]==100);
dat <- em[use,];
nBn <- even.Beta(dat[,4]);
fourblock <- c(1,1.8,2.6,3.4,4.2);
fullblock <- c(fourblock,fourblock+6,fourblock+12,fourblock+18);

sdd <- as.data.frame(cbind(dat[,4],dat[,5],dat[,8]));
colnames(sdd) <- c("bs","cs","vl");
sda <- melt(sdd, id = c('bs','cs'));
sda <- sda[,-3];
boxplots.dat <- boxplot(as.numeric(value)~bs+cs,data=sda,type="n",xaxt="n",
                xlim=c(1,max(fullblock)),yaxt="n",xaxt="n",
                at=fullblock,lwd=1.5,cex.lab=1.35,ylim=c(-16,14),
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5,
                xlab="",
                ylab=expression(paste("Mean polyandry allele value (",P[a],")")));
TST <- NULL;
for(j in unique(dat[,5])){
    for(i in unique(dat[,4])){
        sub <- which(dat[,4]==i & dat[,5]==j);
        cis <- simpleboot(dat[sub,8]);
        ttt <- cis[1] > 0 & cis[2] > 0;
        TST <- rbind(TST,c(j,i,cis,ttt));
    }
}
for(i in 1:dim(TST)[1]){
    if(TST[i,5] == 1){
        lsh <- fullblock[i] - 0.25;
        rsh <- fullblock[i] + 0.25; 
        polygon(x=c(lsh,rsh,rsh,lsh),y=c(-40,-40,20,20),border=NA,col="grey70");
    }
}
boxplots.dat <- boxplot(as.numeric(value)~bs+cs,data=sda,xaxt="n",
                xlim=c(1,max(fullblock)),yaxt="n",xaxt="n",
                at=fullblock,lwd=1.5,cex.lab=1.75,ylim=c(-16,14),
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5,add=TRUE);
abline(h=0,lwd=0.8,lty="dotted");
axis(side=2,at=c(-20,-10,0,10),cex.axis=1.75,cex.lab=2);
text(x=fullblock[18],y=13,cex=2.5,labels="A",pos=4);

par(mar=c(5,5,0.1,1),lwd=2);
# ----------------------- Panel `B':
use <- which(summ[,1] == -1 & summ[,2] == 100);
aus <- which(alle[,1] == -1 & alle[,2] == 100);
dat <- as.data.frame(cbind(summ[use,3],summ[use,4],summ[use,5]));
colnames(dat) <- c("bs","cs","vl");
sda <- melt(dat, id = c('bs','cs'));
sda <- sda[,-3];
ale <- alle[aus,];
fourblock <- c(1,1.8,2.6,3.4,4.2);
fullblock <- c(fourblock,fourblock+6,fourblock+12,fourblock+18);
boxplots.dat <- boxplot(as.numeric(value)~bs+cs,data=sda,xlim=c(1,max(fullblock)),
                xlab=expression(paste("Cost of polyandry (",c[P],")")),
                ylab=expression(paste("Mean polyandry phenotype value (",P[p],")")),
                at=fullblock,lwd=1.5,xaxt="n",yaxt="n",cex.lab=1.35,
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5,ylim=c(0,210));
axis(side=2,at=c(0,100,200),cex.axis=1.75,cex.lab=2);
axis(side=1,at=c(fourblock[3],fourblock[3]+6,fourblock[3]+12,fourblock[3]+18),
     cex.axis=1.5,cex.lab=1.55,labels=c("0.00","0.01","0.02","0.03"));
for(i in 1:dim(ale)[1]){
    if(TST[i,5] == 1){
        lsh <- fullblock[i] - 0.25;
        rsh <- fullblock[i] + 0.25;
        polygon(x=c(lsh,rsh,rsh,lsh),y=c(-10,-10,350,350),border=NA,col="grey70");
    }
}
boxplots.dat <- boxplot(as.numeric(value)~bs+cs,data=sda,xlim=c(1,max(fullblock)),
                at=fullblock,lwd=1.5,add=TRUE,xaxt="n",yaxt="n",
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5,ylim=c(0,210));
text(x=fullblock[18],y=200,cex=2.5,labels="B",pos=4);

dev.off()
# ====================================================================================#
# =====================================================================================
# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX #
# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX #
# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX #
# =====================================================================================
# Some new information from the mat.txt file is needed to make Figures 5 & 6
# The code below prepares this information before the figures are built.
# Four new functions are used to tie everything together properly
mat      <- read.table("mat.txt",header=FALSE);

kbeta    <- which(mat[,4]==0|mat[,4]==0.2|mat[,4]==1|mat[,4]==2|mat[,4]==5);
pol      <- mat[kbeta,];

kcost    <- which(pol[,5]==0|pol[,5]==0.01|pol[,5]==0.02|pol[,5]==0.03);
pol      <- pol[kcost,];
sims     <- unique(cbind(pol[,1:5]));

f.bypr <- function(pol,initial,additional){
    use      <- which(pol[,1] == initial & pol[,2] == additional);
    tus      <- pol[use,];
    pan      <- tus[tus[,7]>1,];
    rm       <- which(pan[,8] > 1 | pan[,9] > 1); #Errors cleaning the file.
    if(length(rm) > 0){
        pan      <- pan[-rm,];
    }
    bypr <- NULL;
    for(i in unique(pan[,3])){
        for(j in unique(pan[,4])){
            for(k in unique(pan[,5])){
               calc <- which(pan[,3]==i & pan[,4]==j & pan[,5]==k);
               if(length(calc) > 0){
                   if(length(calc) > 4){
                       ppp  <- t.test(pan[calc,9]-pan[calc,8])$p.value;
                   }else{
                       ppp  <- NA;
                   }
                   bypr <- rbind(bypr,c(i,j,k,mean(pan[calc,9]-pan[calc,8]),ppp));
               } 
            }
        }
    }
    return(bypr);
}

f.pomn <- function(pol,initial,additional){
    use      <- which(pol[,1] == initial & pol[,2] == additional);
    tus      <- pol[use,];
    ret      <- NULL;
    for(i in unique(tus[,4])){
        for(j in unique(tus[,5])){
            pmat <- tus[tus[,4]==i & tus[,5]==j,];
            pmns <- tapply(X=pmat[,7],INDEX=pmat[,3],FUN=mean);
            pprp <- sum(pmns > 1) / length(pmns);
            pmns <- mean(pmns[pmns > 1]);
            ret  <- rbind(ret,c(i,j,pmns,pprp));
        }
    }
    return(ret);
}

f.CIva <- function(bypr){
    CIva <- NULL;
    for(i in unique(bypr[,2])){
        for(j in unique(bypr[,3])){
            samp <- bypr[bypr[,2]==i & bypr[,3]==j,];
            difs <- samp[,4];
            pvls <- samp[!is.na(samp[,5]),5]<0.05;
            prps <- sum(pvls) / length(pvls);
            difs <- difs[!is.na(difs)];
            CIsm <- simpleboot(difs);
            CIva <- rbind(CIva,c(i,j,mean(difs),CIsm,prps));
        }
    }
    return(CIva);
}

bporder <- function(bp){
    ret <- bp[order(bp[,2]),];
    new <- NULL;
    for(i in unique(ret[,2])){
        tmp <- ret[ret[,2]==i,];
        odr <- tmp[order(tmp[,1]),];
        new <- rbind(new,odr);
    }
    return(new);
}
# We can now use the above functions to actually make Figures 5 and 6 below.
# =====================================================================================
# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX #
# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX #
# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX #
# =====================================================================================
# =====================================================================================
# TODO: BUILD FIGURE 5                                                                #
# ====================================================================================#
# -------------- Build the figure. ============================================ XXX #
setEPS(); # postscript below for final publication?
cairo_ps("FullFact_mn.eps",family="Arial",height=7,width=7);
par(mfrow=c(3,3),mar=c(0.2,0.2,0.2,0.2),oma=c(6,6,1,6));
# ----------------------- Panel `a':
bp1 <- f.bypr(pol=pol,initial=100,additional=100);
bp2 <- f.CIva(bypr=bp1);
bp2 <- bporder(bp2);
tbp <- NULL;
for(i in unique(bp2[,2])){
    tb  <- bp2[bp2[,2]==i,];
    tb  <- tb[order(tb[,1]),];
    tbp <- rbind(tbp,tb);
}
bp2 <- tbp;
bp3 <- f.pomn(pol=pol,initial=100,additional=100);
bp3 <- bporder(bp3);
fourblock <- c(1,1.8,2.6,3.4,4.2);
fullblock <- c(fourblock,fourblock+6,fourblock+12,fourblock+18);
belowpts  <- -0.0045;
plot(x=fullblock[1],y=bp2[1,3],type="n",xlim=c(0,max(fullblock)+1),
     ylim=c(-0.012 + belowpts,0.015),xaxt="n",yaxt="n",xlab="",ylab="");
abline(h=0,lwd=0.9,lty="dotted");
poin <- 21;
for(i in 1:dim(bp2)[1]){
    if(poin == (21 + 5)){
        poin <- 21;
    }
    for(j in unique(bp2[,1])){
        for(k in unique(bp2[,2])){
            if(bp2[i,1] == j & bp2[i,2] == k){
                points(x=fullblock[i],y=bp2[i,3],pch=poin,cex=0.8,bg="black");
                arrows(x0=fullblock[i],x1=fullblock[i],y0=bp2[i,3],y1=bp2[i,4],
                       angle=90,length=0.015,lwd=1.00);
                arrows(x0=fullblock[i],x1=fullblock[i],y0=bp2[i,3],y1=bp2[i,5],
                       angle=90,length=0.015,lwd=1.00);
            }
        }
    }
    poin <- poin + 1;
}
axis(side=2,at=c(-0.01,0,0.01),cex.axis=1.75);
par(new = TRUE);
plot(x=fullblock[1],y=bp2[1,3],type="n",xlim=c(0,max(fullblock)+1),ylim=c(0,40),yaxt="n",
     xaxt="n",xlab="",ylab="");
polygon(x=c(0,fullblock[20]+1,fullblock[20]+1,0),y=c(0,0,8,8),border="black",col="grey80");
segments(x0=0.1,y0=4,x1=fullblock[20]+1,y1=4,lwd=1,lty="dotted");
for(i in 1:dim(bp3)[1]){
    for(j in unique(bp3[,1])){
        for(k in unique(bp3[,2])){
            if(bp3[i,1] == j & bp3[i,2] == k){
                lsh <- fullblock[i] - 0.25;
                rsh <- fullblock[i] + 0.25; 
                polygon(x=c(lsh,rsh,rsh,lsh),y=c(8,8,bp3[i,3]+8,
                        bp3[i,3]+8),border=NA,col="grey70");
                polygon(x=c(lsh,rsh,rsh,lsh),y=c(0,0,bp3[i,4]*8,
                        bp3[i,4]*8),border=NA,col="black");
            }
        }
    }
}
polygon(x=c(0,fullblock[20]+1,fullblock[20]+1,0),y=c(0,0,8,8),border="black",col=NA);
par(new = TRUE);
plot(x=fullblock[1],y=bp2[1,3],type="n",xlim=c(0,max(fullblock)+1),
     ylim=c(-0.012 + belowpts,0.015),yaxt="n",xaxt="n",xlab="",ylab="");
for(i in 1:dim(bp2)[1]){
    if(poin == (21 + 5)){
        poin <- 21;
    }
    for(j in unique(bp2[,1])){
        for(k in unique(bp2[,2])){
            if(bp2[i,1] == j & bp2[i,2] == k){
                points(x=fullblock[i],y=bp2[i,3],pch=poin,cex=0.8,bg="black");
                arrows(x0=fullblock[i],x1=fullblock[i],y0=bp2[i,3],y1=bp2[i,4],
                       angle=90,length=0.015,lwd=1.00);
                arrows(x0=fullblock[i],x1=fullblock[i],y0=bp2[i,3],y1=bp2[i,5],
                       angle=90,length=0.015,lwd=1.00);
            }
        }
    }
    poin <- poin + 1;
}
text(x=-1,y=0.0125,cex=2,labels=expression(paste("(A) ",S['100,100'])),pos=4);
box();
# ----------------------- Panel `b':
bp1 <- f.bypr(pol=pol,initial=100,additional=10);
bp2 <- f.CIva(bypr=bp1);
bp2 <- bporder(bp2);
tbp <- NULL;
for(i in unique(bp2[,2])){
    tb  <- bp2[bp2[,2]==i,];
    tb  <- tb[order(tb[,1]),];
    tbp <- rbind(tbp,tb);
}
bp2 <- tbp;
bp3 <- f.pomn(pol=pol,initial=100,additional=10);
bp3 <- bporder(bp3);
fourblock <- c(1,1.8,2.6,3.4,4.2);
fullblock <- c(fourblock,fourblock+6,fourblock+12,fourblock+18);
belowpts  <- -0.0045;
plot(x=fullblock[1],y=bp2[1,3],type="n",xlim=c(0,max(fullblock)+1),
     ylim=c(-0.012 + belowpts,0.015),xaxt="n",yaxt="n",xlab="",ylab="");
abline(h=0,lwd=0.9,lty="dotted");
poin <- 21;
for(i in 1:dim(bp2)[1]){
    if(poin == (21 + 5)){
        poin <- 21;
    }
    for(j in unique(bp2[,1])){
        for(k in unique(bp2[,2])){
            if(bp2[i,1] == j & bp2[i,2] == k){
                points(x=fullblock[i],y=bp2[i,3],pch=poin,cex=0.8,bg="black");
                arrows(x0=fullblock[i],x1=fullblock[i],y0=bp2[i,3],y1=bp2[i,4],
                       angle=90,length=0.015,lwd=1.00);
                arrows(x0=fullblock[i],x1=fullblock[i],y0=bp2[i,3],y1=bp2[i,5],
                       angle=90,length=0.015,lwd=1.00);
            }
        }
    }
    poin <- poin + 1;
}
par(new = TRUE);
plot(x=fullblock[1],y=bp2[1,3],type="n",xlim=c(0,max(fullblock)+1),ylim=c(0,40),yaxt="n",
     xaxt="n",xlab="",ylab="");
polygon(x=c(0,fullblock[20]+1,fullblock[20]+1,0),y=c(0,0,8,8),border="black",col="grey80");
segments(x0=0.1,y0=4,x1=fullblock[20]+1,y1=4,lwd=1,lty="dotted");
for(i in 1:dim(bp3)[1]){
    for(j in unique(bp3[,1])){
        for(k in unique(bp3[,2])){
            if(bp3[i,1] == j & bp3[i,2] == k){
                lsh <- fullblock[i] - 0.25;
                rsh <- fullblock[i] + 0.25; 
                polygon(x=c(lsh,rsh,rsh,lsh),y=c(8,8,bp3[i,3]+8,
                        bp3[i,3]+8),border=NA,col="grey70");
                polygon(x=c(lsh,rsh,rsh,lsh),y=c(0,0,bp3[i,4]*8,
                        bp3[i,4]*8),border=NA,col="black");
            }
        }
    }
}
polygon(x=c(0,fullblock[20]+1,fullblock[20]+1,0),y=c(0,0,8,8),border="black",col=NA);
par(new = TRUE);
plot(x=fullblock[1],y=bp2[1,3],type="n",xlim=c(0,max(fullblock)+1),
     ylim=c(-0.012 + belowpts,0.015),yaxt="n",xaxt="n",xlab="",ylab="");
for(i in 1:dim(bp2)[1]){
    if(poin == (21 + 5)){
        poin <- 21;
    }
    for(j in unique(bp2[,1])){
        for(k in unique(bp2[,2])){
            if(bp2[i,1] == j & bp2[i,2] == k){
                points(x=fullblock[i],y=bp2[i,3],pch=poin,cex=0.8,bg="black");
                arrows(x0=fullblock[i],x1=fullblock[i],y0=bp2[i,3],y1=bp2[i,4],
                       angle=90,length=0.015,lwd=1.00);
                arrows(x0=fullblock[i],x1=fullblock[i],y0=bp2[i,3],y1=bp2[i,5],
                       angle=90,length=0.015,lwd=1.00);
            }
        }
    }
    poin <- poin + 1;
}
text(x=-1,y=0.0125,cex=2,labels=expression(paste("(B) ",S['100,10'])),pos=4);
box();
# ----------------------- Panel `c':
bp1 <- f.bypr(pol=pol,initial=100,additional=2);
bp2 <- f.CIva(bypr=bp1);
bp2 <- bporder(bp2);
tbp <- NULL;
for(i in unique(bp2[,2])){
    tb  <- bp2[bp2[,2]==i,];
    tb  <- tb[order(tb[,1]),];
    tbp <- rbind(tbp,tb);
}
bp2 <- tbp;
bp3 <- f.pomn(pol=pol,initial=100,additional=2);
bp3 <- bporder(bp3);
fourblock <- c(1,1.8,2.6,3.4,4.2);
fullblock <- c(fourblock,fourblock+6,fourblock+12,fourblock+18);
belowpts  <- -0.0045;
plot(x=fullblock[1],y=bp2[1,3],type="n",xlim=c(0,max(fullblock)+1),
     ylim=c(-0.012 + belowpts,0.015),xaxt="n",yaxt="n",xlab="",ylab="");
abline(h=0,lwd=0.9,lty="dotted");
poin <- 21;
for(i in 1:dim(bp2)[1]){
    if(poin == (21 + 5)){
        poin <- 21;
    }
    for(j in unique(bp2[,1])){
        for(k in unique(bp2[,2])){
            if(bp2[i,1] == j & bp2[i,2] == k){
                points(x=fullblock[i],y=bp2[i,3],pch=poin,cex=0.8,bg="black");
                arrows(x0=fullblock[i],x1=fullblock[i],y0=bp2[i,3],y1=bp2[i,4],
                       angle=90,length=0.015,lwd=1.00);
                arrows(x0=fullblock[i],x1=fullblock[i],y0=bp2[i,3],y1=bp2[i,5],
                       angle=90,length=0.015,lwd=1.00);
            }
        }
    }
    poin <- poin + 1;
}
par(new = TRUE);
plot(x=fullblock[1],y=bp2[1,3],type="n",xlim=c(0,max(fullblock)+1),ylim=c(0,40),yaxt="n",
     xaxt="n",xlab="",ylab="");
polygon(x=c(0,fullblock[20]+1,fullblock[20]+1,0),y=c(0,0,8,8),border="black",col="grey80");
segments(x0=0.1,y0=4,x1=fullblock[20]+1,y1=4,lwd=1,lty="dotted");
for(i in 1:dim(bp3)[1]){
    for(j in unique(bp3[,1])){
        for(k in unique(bp3[,2])){
            if(bp3[i,1] == j & bp3[i,2] == k){
                lsh <- fullblock[i] - 0.25;
                rsh <- fullblock[i] + 0.25; 
                polygon(x=c(lsh,rsh,rsh,lsh),y=c(8,8,bp3[i,3]+8,
                        bp3[i,3]+8),border=NA,col="grey70");
                polygon(x=c(lsh,rsh,rsh,lsh),y=c(0,0,bp3[i,4]*8,
                        bp3[i,4]*8),border=NA,col="black");
            }
        }
    }
}
polygon(x=c(0,fullblock[20]+1,fullblock[20]+1,0),y=c(0,0,8,8),border="black",col=NA);
axis(side=4,at=c(8,18,28,38),labels=c("0","10","20","30"),cex.axis=1.75);
par(new = TRUE);
plot(x=fullblock[1],y=bp2[1,3],type="n",xlim=c(0,max(fullblock)+1),
     ylim=c(-0.012 + belowpts,0.015),yaxt="n",xaxt="n",xlab="",ylab="");
for(i in 1:dim(bp2)[1]){
    if(poin == (21 + 5)){
        poin <- 21;
    }
    for(j in unique(bp2[,1])){
        for(k in unique(bp2[,2])){
            if(bp2[i,1] == j & bp2[i,2] == k){
                points(x=fullblock[i],y=bp2[i,3],pch=poin,cex=0.8,bg="black");
                arrows(x0=fullblock[i],x1=fullblock[i],y0=bp2[i,3],y1=bp2[i,4],
                       angle=90,length=0.015,lwd=1.00);
                arrows(x0=fullblock[i],x1=fullblock[i],y0=bp2[i,3],y1=bp2[i,5],
                       angle=90,length=0.015,lwd=1.00);
            }
        }
    }
    poin <- poin + 1;
}
text(x=-1,y=0.0125,cex=2,labels=expression(paste("(C) ",S['100,2'])),pos=4);
par(new = TRUE);
plot(x=fullblock[1],y=bp2[1,3],type="n",xlim=c(0,max(fullblock)+1),ylim=c(0,40),yaxt="n",
     xaxt="n",xlab="",ylab="");
polygon(x=c(0,fullblock[20]+1,fullblock[20]+1,0),y=c(0,0,8,8),border="black",col="grey80");
segments(x0=0.1,y0=4,x1=fullblock[20]+1,y1=4,lwd=1,lty="dotted");
for(i in 1:dim(bp3)[1]){
    for(j in unique(bp3[,1])){
        for(k in unique(bp3[,2])){
            if(bp3[i,1] == j & bp3[i,2] == k){
                lsh <- fullblock[i] - 0.25;
                rsh <- fullblock[i] + 0.25; 
                polygon(x=c(lsh,rsh,rsh,lsh),y=c(0,0,bp3[i,4]*8,
                        bp3[i,4]*8),border=NA,col="black");
            }
        }
    }
}
box();
# ----------------------- Panel `d':
bp1 <- f.bypr(pol=pol,initial=10,additional=100);
bp2 <- f.CIva(bypr=bp1);
bp2 <- bporder(bp2);
tbp <- NULL;
for(i in unique(bp2[,2])){
    tb  <- bp2[bp2[,2]==i,];
    tb  <- tb[order(tb[,1]),];
    tbp <- rbind(tbp,tb);
}
bp2 <- tbp;
bp3 <- f.pomn(pol=pol,initial=10,additional=100);
bp3 <- bporder(bp3);
fourblock <- c(1,1.8,2.6,3.4,4.2);
fullblock <- c(fourblock,fourblock+6,fourblock+12,fourblock+18);
belowpts  <- -0.0045;
plot(x=fullblock[1],y=bp2[1,3],type="n",xlim=c(0,max(fullblock)+1),
     ylim=c(-0.012 + belowpts,0.015),xaxt="n",yaxt="n",xlab="",ylab="");
abline(h=0,lwd=0.9,lty="dotted");
poin <- 21;
for(i in 1:dim(bp2)[1]){
    if(poin == (21 + 5)){
        poin <- 21;
    }
    for(j in unique(bp2[,1])){
        for(k in unique(bp2[,2])){
            if(bp2[i,1] == j & bp2[i,2] == k){
                points(x=fullblock[i],y=bp2[i,3],pch=poin,cex=0.8,bg="black");
                arrows(x0=fullblock[i],x1=fullblock[i],y0=bp2[i,3],y1=bp2[i,4],
                       angle=90,length=0.015,lwd=1.00);
                arrows(x0=fullblock[i],x1=fullblock[i],y0=bp2[i,3],y1=bp2[i,5],
                       angle=90,length=0.015,lwd=1.00);
            }
        }
    }
    poin <- poin + 1;
}
axis(side=2,at=c(-0.01,0,0.01),cex.axis=1.75);
par(new = TRUE);
plot(x=fullblock[1],y=bp2[1,3],type="n",xlim=c(0,max(fullblock)+1),ylim=c(0,40),yaxt="n",
     xaxt="n",xlab="",ylab="");
polygon(x=c(0,fullblock[20]+1,fullblock[20]+1,0),y=c(0,0,8,8),border="black",col="grey80");
segments(x0=0.1,y0=4,x1=fullblock[20]+1,y1=4,lwd=1,lty="dotted");
for(i in 1:dim(bp3)[1]){
    for(j in unique(bp3[,1])){
        for(k in unique(bp3[,2])){
            if(bp3[i,1] == j & bp3[i,2] == k){
                lsh <- fullblock[i] - 0.25;
                rsh <- fullblock[i] + 0.25; 
                polygon(x=c(lsh,rsh,rsh,lsh),y=c(8,8,bp3[i,3]+8,
                        bp3[i,3]+8),border=NA,col="grey70");
                polygon(x=c(lsh,rsh,rsh,lsh),y=c(0,0,bp3[i,4]*8,
                        bp3[i,4]*8),border=NA,col="black");
            }
        }
    }
}
polygon(x=c(0,fullblock[20]+1,fullblock[20]+1,0),y=c(0,0,8,8),border="black",col=NA);
par(new = TRUE);
plot(x=fullblock[1],y=bp2[1,3],type="n",xlim=c(0,max(fullblock)+1),
     ylim=c(-0.012 + belowpts,0.015),yaxt="n",xaxt="n",xlab="",ylab="");
for(i in 1:dim(bp2)[1]){
    if(poin == (21 + 5)){
        poin <- 21;
    }
    for(j in unique(bp2[,1])){
        for(k in unique(bp2[,2])){
            if(bp2[i,1] == j & bp2[i,2] == k){
                points(x=fullblock[i],y=bp2[i,3],pch=poin,cex=0.8,bg="black");
                arrows(x0=fullblock[i],x1=fullblock[i],y0=bp2[i,3],y1=bp2[i,4],
                       angle=90,length=0.015,lwd=1.00);
                arrows(x0=fullblock[i],x1=fullblock[i],y0=bp2[i,3],y1=bp2[i,5],
                       angle=90,length=0.015,lwd=1.00);
            }
        }
    }
    poin <- poin + 1;
}
text(x=-1,y=0.0125,cex=2,labels=expression(paste("(D) ",S['10,100'])),pos=4);
par(new = TRUE);
plot(x=fullblock[1],y=bp2[1,3],type="n",xlim=c(0,max(fullblock)+1),ylim=c(0,40),yaxt="n",
     xaxt="n",xlab="",ylab="");
polygon(x=c(0,fullblock[20]+1,fullblock[20]+1,0),y=c(0,0,8,8),border="black",col="grey80");
segments(x0=0.1,y0=4,x1=fullblock[20]+1,y1=4,lwd=1,lty="dotted");
for(i in 1:dim(bp3)[1]){
    for(j in unique(bp3[,1])){
        for(k in unique(bp3[,2])){
            if(bp3[i,1] == j & bp3[i,2] == k){
                lsh <- fullblock[i] - 0.25;
                rsh <- fullblock[i] + 0.25; 
                polygon(x=c(lsh,rsh,rsh,lsh),y=c(0,0,bp3[i,4]*8,
                        bp3[i,4]*8),border=NA,col="black");
            }
        }
    }
}
box();
# ----------------------- Panel `e':
bp1 <- f.bypr(pol=pol,initial=10,additional=10);
bp2 <- f.CIva(bypr=bp1);
bp2 <- bporder(bp2);
tbp <- NULL;
for(i in unique(bp2[,2])){
    tb  <- bp2[bp2[,2]==i,];
    tb  <- tb[order(tb[,1]),];
    tbp <- rbind(tbp,tb);
}
bp2 <- tbp;
bp3 <- f.pomn(pol=pol,initial=10,additional=10);
bp3 <- bporder(bp3);
fourblock <- c(1,1.8,2.6,3.4,4.2);
fullblock <- c(fourblock,fourblock+6,fourblock+12,fourblock+18);
belowpts  <- -0.0045;
plot(x=fullblock[1],y=bp2[1,3],type="n",xlim=c(0,max(fullblock)+1),
     ylim=c(-0.012 + belowpts,0.015),xaxt="n",yaxt="n",xlab="",ylab="");
abline(h=0,lwd=0.9,lty="dotted");
poin <- 21;
for(i in 1:dim(bp2)[1]){
    if(poin == (21 + 5)){
        poin <- 21;
    }
    for(j in unique(bp2[,1])){
        for(k in unique(bp2[,2])){
            if(bp2[i,1] == j & bp2[i,2] == k){
                points(x=fullblock[i],y=bp2[i,3],pch=poin,cex=0.8,bg="black");
                arrows(x0=fullblock[i],x1=fullblock[i],y0=bp2[i,3],y1=bp2[i,4],
                       angle=90,length=0.015,lwd=1.00);
                arrows(x0=fullblock[i],x1=fullblock[i],y0=bp2[i,3],y1=bp2[i,5],
                       angle=90,length=0.015,lwd=1.00);
            }
        }
    }
    poin <- poin + 1;
}
par(new = TRUE);
plot(x=fullblock[1],y=bp2[1,3],type="n",xlim=c(0,max(fullblock)+1),ylim=c(0,40),yaxt="n",
     xaxt="n",xlab="",ylab="");
polygon(x=c(0,fullblock[20]+1,fullblock[20]+1,0),y=c(0,0,8,8),border="black",col="grey80");
segments(x0=0.1,y0=4,x1=fullblock[20]+1,y1=4,lwd=1,lty="dotted");
for(i in 1:dim(bp3)[1]){
    for(j in unique(bp3[,1])){
        for(k in unique(bp3[,2])){
            if(bp3[i,1] == j & bp3[i,2] == k){
                lsh <- fullblock[i] - 0.25;
                rsh <- fullblock[i] + 0.25; 
                polygon(x=c(lsh,rsh,rsh,lsh),y=c(8,8,bp3[i,3]+8,
                        bp3[i,3]+8),border=NA,col="grey70");
                polygon(x=c(lsh,rsh,rsh,lsh),y=c(0,0,bp3[i,4]*8,
                        bp3[i,4]*8),border=NA,col="black");
            }
        }
    }
}
polygon(x=c(0,fullblock[20]+1,fullblock[20]+1,0),y=c(0,0,8,8),border="black",col=NA);
par(new = TRUE);
plot(x=fullblock[1],y=bp2[1,3],type="n",xlim=c(0,max(fullblock)+1),
     ylim=c(-0.012 + belowpts,0.015),yaxt="n",xaxt="n",xlab="",ylab="");
for(i in 1:dim(bp2)[1]){
    if(poin == (21 + 5)){
        poin <- 21;
    }
    for(j in unique(bp2[,1])){
        for(k in unique(bp2[,2])){
            if(bp2[i,1] == j & bp2[i,2] == k){
                points(x=fullblock[i],y=bp2[i,3],pch=poin,cex=0.8,bg="black");
                arrows(x0=fullblock[i],x1=fullblock[i],y0=bp2[i,3],y1=bp2[i,4],
                       angle=90,length=0.015,lwd=1.00);
                arrows(x0=fullblock[i],x1=fullblock[i],y0=bp2[i,3],y1=bp2[i,5],
                       angle=90,length=0.015,lwd=1.00);
            }
        }
    }
    poin <- poin + 1;
}
text(x=-1,y=0.0125,cex=2,labels=expression(paste("(E) ",S['10,10'])),pos=4);
par(new = TRUE);
plot(x=fullblock[1],y=bp2[1,3],type="n",xlim=c(0,max(fullblock)+1),ylim=c(0,40),yaxt="n",
     xaxt="n",xlab="",ylab="");
polygon(x=c(0,fullblock[20]+1,fullblock[20]+1,0),y=c(0,0,8,8),border="black",col="grey80");
segments(x0=0.1,y0=4,x1=fullblock[20]+1,y1=4,lwd=1,lty="dotted");
for(i in 1:dim(bp3)[1]){
    for(j in unique(bp3[,1])){
        for(k in unique(bp3[,2])){
            if(bp3[i,1] == j & bp3[i,2] == k){
                lsh <- fullblock[i] - 0.25;
                rsh <- fullblock[i] + 0.25; 
                polygon(x=c(lsh,rsh,rsh,lsh),y=c(8,8,bp3[i,3]+8,
                        bp3[i,3]+8),border=NA,col="grey70");
                polygon(x=c(lsh,rsh,rsh,lsh),y=c(0,0,bp3[i,4]*8,
                        bp3[i,4]*8),border=NA,col="black");
            }
        }
    }
}
box();
# ----------------------- Panel `f':
bp1 <- f.bypr(pol=pol,initial=10,additional=2);
bp2 <- f.CIva(bypr=bp1);
bp2 <- bporder(bp2);
tbp <- NULL;
for(i in unique(bp2[,2])){
    tb  <- bp2[bp2[,2]==i,];
    tb  <- tb[order(tb[,1]),];
    tbp <- rbind(tbp,tb);
}
bp2 <- tbp;
bp3 <- f.pomn(pol=pol,initial=10,additional=2);
bp3 <- bporder(bp3);
fourblock <- c(1,1.8,2.6,3.4,4.2);
fullblock <- c(fourblock,fourblock+6,fourblock+12,fourblock+18);
belowpts  <- -0.0045;
plot(x=fullblock[1],y=bp2[1,3],type="n",xlim=c(0,max(fullblock)+1),
     ylim=c(-0.012 + belowpts,0.015),xaxt="n",yaxt="n",xlab="",ylab="");
abline(h=0,lwd=0.9,lty="dotted");
poin <- 21;
for(i in 1:dim(bp2)[1]){
    if(poin == (21 + 5)){
        poin <- 21;
    }
    for(j in unique(bp2[,1])){
        for(k in unique(bp2[,2])){
            if(bp2[i,1] == j & bp2[i,2] == k){
                points(x=fullblock[i],y=bp2[i,3],pch=poin,cex=0.8,bg="black");
                arrows(x0=fullblock[i],x1=fullblock[i],y0=bp2[i,3],y1=bp2[i,4],
                       angle=90,length=0.015,lwd=1.00);
                arrows(x0=fullblock[i],x1=fullblock[i],y0=bp2[i,3],y1=bp2[i,5],
                       angle=90,length=0.015,lwd=1.00);
            }
        }
    }
    poin <- poin + 1;
}
par(new = TRUE);
plot(x=fullblock[1],y=bp2[1,3],type="n",xlim=c(0,max(fullblock)+1),ylim=c(0,40),yaxt="n",
     xaxt="n",xlab="",ylab="");
polygon(x=c(0,fullblock[20]+1,fullblock[20]+1,0),y=c(0,0,8,8),border="black",col="grey80");
segments(x0=0.1,y0=4,x1=fullblock[20]+1,y1=4,lwd=1,lty="dotted");
for(i in 1:dim(bp3)[1]){
    for(j in unique(bp3[,1])){
        for(k in unique(bp3[,2])){
            if(bp3[i,1] == j & bp3[i,2] == k){
                lsh <- fullblock[i] - 0.25;
                rsh <- fullblock[i] + 0.25; 
                polygon(x=c(lsh,rsh,rsh,lsh),y=c(8,8,bp3[i,3]+8,
                        bp3[i,3]+8),border=NA,col="grey70");
                polygon(x=c(lsh,rsh,rsh,lsh),y=c(0,0,bp3[i,4]*8,
                        bp3[i,4]*8),border=NA,col="black");
            }
        }
    }
}
polygon(x=c(0,fullblock[20]+1,fullblock[20]+1,0),y=c(0,0,8,8),border="black",col=NA);
axis(side=4,at=c(8,18,28,38),labels=c("0","10","20","30"),cex.axis=1.75);
par(new = TRUE);
plot(x=fullblock[1],y=bp2[1,3],type="n",xlim=c(0,max(fullblock)+1),
     ylim=c(-0.012 + belowpts,0.015),yaxt="n",xaxt="n",xlab="",ylab="");
for(i in 1:dim(bp2)[1]){
    if(poin == (21 + 5)){
        poin <- 21;
    }
    for(j in unique(bp2[,1])){
        for(k in unique(bp2[,2])){
            if(bp2[i,1] == j & bp2[i,2] == k){
                points(x=fullblock[i],y=bp2[i,3],pch=poin,cex=0.8,bg="black");
                arrows(x0=fullblock[i],x1=fullblock[i],y0=bp2[i,3],y1=bp2[i,4],
                       angle=90,length=0.015,lwd=1.00);
                arrows(x0=fullblock[i],x1=fullblock[i],y0=bp2[i,3],y1=bp2[i,5],
                       angle=90,length=0.015,lwd=1.00);
            }
        }
    }
    poin <- poin + 1;
}
text(x=-1,y=0.0125,cex=2,labels=expression(paste("(F) ",S['10,2'])),pos=4);
par(new = TRUE);
plot(x=fullblock[1],y=bp2[1,3],type="n",xlim=c(0,max(fullblock)+1),ylim=c(0,40),yaxt="n",
     xaxt="n",xlab="",ylab="");
polygon(x=c(0,fullblock[20]+1,fullblock[20]+1,0),y=c(0,0,8,8),border="black",col="grey80");
segments(x0=0.1,y0=4,x1=fullblock[20]+1,y1=4,lwd=1,lty="dotted");
for(i in 1:dim(bp3)[1]){
    for(j in unique(bp3[,1])){
        for(k in unique(bp3[,2])){
            if(bp3[i,1] == j & bp3[i,2] == k){
                lsh <- fullblock[i] - 0.25;
                rsh <- fullblock[i] + 0.25; 
                polygon(x=c(lsh,rsh,rsh,lsh),y=c(0,0,bp3[i,4]*8,
                        bp3[i,4]*8),border=NA,col="black");
            }
        }
    }
}
box();
# ----------------------- Panel `g':
bp1 <- f.bypr(pol=pol,initial=2,additional=100);
bp2 <- f.CIva(bypr=bp1);
bp2 <- bporder(bp2);
tbp <- NULL;
for(i in unique(bp2[,2])){
    tb  <- bp2[bp2[,2]==i,];
    tb  <- tb[order(tb[,1]),];
    tbp <- rbind(tbp,tb);
}
bp2 <- tbp;
bp3 <- f.pomn(pol=pol,initial=2,additional=100);
bp3 <- bporder(bp3);
fourblock <- c(1,1.8,2.6,3.4,4.2);
fullblock <- c(fourblock,fourblock+6,fourblock+12,fourblock+18);
belowpts  <- -0.0045;
plot(x=fullblock[1],y=bp2[1,3],type="n",xlim=c(0,max(fullblock)+1),
     ylim=c(-0.012 + belowpts,0.015),xaxt="n",yaxt="n",xlab="",ylab="");
abline(h=0,lwd=0.9,lty="dotted");
poin <- 21;
for(i in 1:dim(bp2)[1]){
    if(poin == (21 + 5)){
        poin <- 21;
    }
    for(j in unique(bp2[,1])){
        for(k in unique(bp2[,2])){
            if(bp2[i,1] == j & bp2[i,2] == k){
                points(x=fullblock[i],y=bp2[i,3],pch=poin,cex=0.8,bg="black");
                arrows(x0=fullblock[i],x1=fullblock[i],y0=bp2[i,3],y1=bp2[i,4],
                       angle=90,length=0.015,lwd=1.00);
                arrows(x0=fullblock[i],x1=fullblock[i],y0=bp2[i,3],y1=bp2[i,5],
                       angle=90,length=0.015,lwd=1.00);
            }
        }
    }
    poin <- poin + 1;
}
axis(side=2,at=c(-0.01,0,0.01),cex.axis=1.75);
par(new = TRUE);
plot(x=fullblock[1],y=bp2[1,3],type="n",xlim=c(0,max(fullblock)+1),ylim=c(0,40),yaxt="n",
     xaxt="n",xlab="",ylab="");
polygon(x=c(0,fullblock[20]+1,fullblock[20]+1,0),y=c(0,0,8,8),border="black",col="grey80");
segments(x0=0.1,y0=4,x1=fullblock[20]+1,y1=4,lwd=1,lty="dotted");
for(i in 1:dim(bp3)[1]){
    for(j in unique(bp3[,1])){
        for(k in unique(bp3[,2])){
            if(bp3[i,1] == j & bp3[i,2] == k){
                lsh <- fullblock[i] - 0.25;
                rsh <- fullblock[i] + 0.25; 
                polygon(x=c(lsh,rsh,rsh,lsh),y=c(8,8,bp3[i,3]+8,
                        bp3[i,3]+8),border=NA,col="grey70");
                polygon(x=c(lsh,rsh,rsh,lsh),y=c(0,0,bp3[i,4]*8,
                        bp3[i,4]*8),border=NA,col="black");
            }
        }
    }
}
polygon(x=c(0,fullblock[20]+1,fullblock[20]+1,0),y=c(0,0,8,8),border="black",col=NA);
par(new = TRUE);
plot(x=fullblock[1],y=bp2[1,3],type="n",xlim=c(0,max(fullblock)+1),
     ylim=c(-0.012 + belowpts,0.015),yaxt="n",xaxt="n",xlab="",ylab="");
for(i in 1:dim(bp2)[1]){
    if(poin == (21 + 5)){
        poin <- 21;
    }
    for(j in unique(bp2[,1])){
        for(k in unique(bp2[,2])){
            if(bp2[i,1] == j & bp2[i,2] == k){
                points(x=fullblock[i],y=bp2[i,3],pch=poin,cex=0.8,bg="black");
                arrows(x0=fullblock[i],x1=fullblock[i],y0=bp2[i,3],y1=bp2[i,4],
                       angle=90,length=0.015,lwd=1.00);
                arrows(x0=fullblock[i],x1=fullblock[i],y0=bp2[i,3],y1=bp2[i,5],
                       angle=90,length=0.015,lwd=1.00);
            }
        }
    }
    poin <- poin + 1;
}
text(x=-1,y=0.0125,cex=2,labels=expression(paste("(G) ",S['2,100'])),pos=4);
axis(side=1,at=c(fourblock[3],fourblock[3]+6,fourblock[3]+12,fourblock[3]+18),
     cex.axis=1.5,cex.lab=1.35,labels=c("0.00","0.01","0.02","0.03"));
box();
# ----------------------- Panel `h':
bp1 <- f.bypr(pol=pol,initial=2,additional=10);
bp2 <- f.CIva(bypr=bp1);
bp2 <- bporder(bp2);
tbp <- NULL;
for(i in unique(bp2[,2])){
    tb  <- bp2[bp2[,2]==i,];
    tb  <- tb[order(tb[,1]),];
    tbp <- rbind(tbp,tb);
}
bp2 <- tbp;
bp3 <- f.pomn(pol=pol,initial=2,additional=10);
bp3 <- bporder(bp3);
fourblock <- c(1,1.8,2.6,3.4,4.2);
fullblock <- c(fourblock,fourblock+6,fourblock+12,fourblock+18);
belowpts  <- -0.0045;
plot(x=fullblock[1],y=bp2[1,3],type="n",xlim=c(0,max(fullblock)+1),
     ylim=c(-0.012 + belowpts,0.015),xaxt="n",yaxt="n",xlab="",ylab="");
abline(h=0,lwd=0.9,lty="dotted");
poin <- 21;
for(i in 1:dim(bp2)[1]){
    if(poin == (21 + 5)){
        poin <- 21;
    }
    for(j in unique(bp2[,1])){
        for(k in unique(bp2[,2])){
            if(bp2[i,1] == j & bp2[i,2] == k){
                points(x=fullblock[i],y=bp2[i,3],pch=poin,cex=0.8,bg="black");
                arrows(x0=fullblock[i],x1=fullblock[i],y0=bp2[i,3],y1=bp2[i,4],
                       angle=90,length=0.015,lwd=1.00);
                arrows(x0=fullblock[i],x1=fullblock[i],y0=bp2[i,3],y1=bp2[i,5],
                       angle=90,length=0.015,lwd=1.00);
            }
        }
    }
    poin <- poin + 1;
}
par(new = TRUE);
plot(x=fullblock[1],y=bp2[1,3],type="n",xlim=c(0,max(fullblock)+1),ylim=c(0,40),yaxt="n",
     xaxt="n",xlab="",ylab="");
polygon(x=c(0,fullblock[20]+1,fullblock[20]+1,0),y=c(0,0,8,8),border="black",col="grey80");
segments(x0=0.1,y0=4,x1=fullblock[20]+1,y1=4,lwd=1,lty="dotted");
for(i in 1:dim(bp3)[1]){
    for(j in unique(bp3[,1])){
        for(k in unique(bp3[,2])){
            if(bp3[i,1] == j & bp3[i,2] == k){
                lsh <- fullblock[i] - 0.25;
                rsh <- fullblock[i] + 0.25; 
                polygon(x=c(lsh,rsh,rsh,lsh),y=c(8,8,bp3[i,3]+8,
                        bp3[i,3]+8),border=NA,col="grey70");
                polygon(x=c(lsh,rsh,rsh,lsh),y=c(0,0,bp3[i,4]*8,
                        bp3[i,4]*8),border=NA,col="black");
            }
        }
    }
}
polygon(x=c(0,fullblock[20]+1,fullblock[20]+1,0),y=c(0,0,8,8),border="black",col=NA);
par(new = TRUE);
plot(x=fullblock[1],y=bp2[1,3],type="n",xlim=c(0,max(fullblock)+1),
     ylim=c(-0.012 + belowpts,0.015),yaxt="n",xaxt="n",xlab="",ylab="");
for(i in 1:dim(bp2)[1]){
    if(poin == (21 + 5)){
        poin <- 21;
    }
    for(j in unique(bp2[,1])){
        for(k in unique(bp2[,2])){
            if(bp2[i,1] == j & bp2[i,2] == k){
                points(x=fullblock[i],y=bp2[i,3],pch=poin,cex=0.8,bg="black");
                arrows(x0=fullblock[i],x1=fullblock[i],y0=bp2[i,3],y1=bp2[i,4],
                       angle=90,length=0.015,lwd=1.00);
                arrows(x0=fullblock[i],x1=fullblock[i],y0=bp2[i,3],y1=bp2[i,5],
                       angle=90,length=0.015,lwd=1.00);
            }
        }
    }
    poin <- poin + 1;
}
text(x=-1,y=0.0125,cex=2,labels=expression(paste("(H) ",S['2,10'])),pos=4);
axis(side=1,at=c(fourblock[3],fourblock[3]+6,fourblock[3]+12,fourblock[3]+18),
     cex.axis=1.5,cex.lab=1.35,labels=c("0.00","0.01","0.02","0.03"));
par(new = TRUE);
plot(x=fullblock[1],y=bp2[1,3],type="n",xlim=c(0,max(fullblock)+1),ylim=c(0,40),yaxt="n",
     xaxt="n",xlab="",ylab="");
polygon(x=c(0,fullblock[20]+1,fullblock[20]+1,0),y=c(0,0,8,8),border="black",col="grey80");
segments(x0=0.1,y0=4,x1=fullblock[20]+1,y1=4,lwd=1,lty="dotted");
for(i in 1:dim(bp3)[1]){
    for(j in unique(bp3[,1])){
        for(k in unique(bp3[,2])){
            if(bp3[i,1] == j & bp3[i,2] == k){
                lsh <- fullblock[i] - 0.25;
                rsh <- fullblock[i] + 0.25; 
                polygon(x=c(lsh,rsh,rsh,lsh),y=c(0,0,bp3[i,4]*8,
                        bp3[i,4]*8),border=NA,col="black");
            }
        }
    }
}
box();
# ----------------------- Panel `i':
bp1 <- f.bypr(pol=pol,initial=2,additional=2);
bp2 <- f.CIva(bypr=bp1);
bp2 <- bporder(bp2);
tbp <- NULL;
for(i in unique(bp2[,2])){
    tb  <- bp2[bp2[,2]==i,];
    tb  <- tb[order(tb[,1]),];
    tbp <- rbind(tbp,tb);
}
bp2 <- tbp;
bp3 <- f.pomn(pol=pol,initial=2,additional=2);
bp3 <- bporder(bp3);

fourblock <- c(1,1.8,2.6,3.4,4.2);
fullblock <- c(fourblock,fourblock+6,fourblock+12,fourblock+18);
belowpts  <- -0.0045;
bp2[bp2[,4] <= -1*0.015,4] <- -0.015; # Avoids image run-off into the bottom.
plot(x=fullblock[1],y=bp2[1,3],type="n",xlim=c(0,max(fullblock)+1),
     ylim=c(-0.012 + belowpts,0.015),xaxt="n",yaxt="n",xlab="",ylab="");
abline(h=0,lwd=0.9,lty="dotted");
poin <- 21;
for(i in 1:dim(bp2)[1]){
    if(poin == (21 + 5)){
        poin <- 21;
    }
    for(j in unique(bp2[,1])){
        for(k in unique(bp2[,2])){
            if(bp2[i,1] == j & bp2[i,2] == k){
                points(x=fullblock[i],y=bp2[i,3],pch=poin,cex=0.8,bg="black");
                arrows(x0=fullblock[i],x1=fullblock[i],y0=bp2[i,3],y1=bp2[i,4],
                       angle=90,length=0.015,lwd=1.00);
                arrows(x0=fullblock[i],x1=fullblock[i],y0=bp2[i,3],y1=bp2[i,5],
                       angle=90,length=0.015,lwd=1.00);
            }
        }
    }
    poin <- poin + 1;
}
par(new = TRUE);
plot(x=fullblock[1],y=bp2[1,3],type="n",xlim=c(0,max(fullblock)+1),ylim=c(0,40),yaxt="n",
     xaxt="n",xlab="",ylab="");
polygon(x=c(0,fullblock[20]+1,fullblock[20]+1,0),y=c(0,0,8,8),border="black",col="grey80");
segments(x0=0.1,y0=4,x1=fullblock[20]+1,y1=4,lwd=1,lty="dotted");
for(i in 1:dim(bp3)[1]){
    for(j in unique(bp3[,1])){
        for(k in unique(bp3[,2])){
            if(bp3[i,1] == j & bp3[i,2] == k){
                lsh <- fullblock[i] - 0.25;
                rsh <- fullblock[i] + 0.25; 
                polygon(x=c(lsh,rsh,rsh,lsh),y=c(8,8,bp3[i,3]+8,
                        bp3[i,3]+8),border=NA,col="grey70");
                polygon(x=c(lsh,rsh,rsh,lsh),y=c(0,0,bp3[i,4]*8,
                        bp3[i,4]*8),border=NA,col="black");
            }
        }
    }
}
polygon(x=c(0,fullblock[20]+1,fullblock[20]+1,0),y=c(0,0,8,8),border="black",col=NA);
axis(side=4,at=c(8,18,28,38),labels=c("0","10","20","30"),cex.axis=1.75);
par(new = TRUE);
plot(x=fullblock[1],y=bp2[1,3],type="n",xlim=c(0,max(fullblock)+1),
     ylim=c(-0.012 + belowpts,0.015),yaxt="n",xaxt="n",xlab="",ylab="");
for(i in 1:dim(bp2)[1]){
    if(poin == (21 + 5)){
        poin <- 21;
    }
    for(j in unique(bp2[,1])){
        for(k in unique(bp2[,2])){
            if(bp2[i,1] == j & bp2[i,2] == k){
                points(x=fullblock[i],y=bp2[i,3],pch=poin,cex=0.8,bg="black");
                arrows(x0=fullblock[i],x1=fullblock[i],y0=bp2[i,3],y1=bp2[i,4],
                       angle=90,length=0.015,lwd=1.00);
                arrows(x0=fullblock[i],x1=fullblock[i],y0=bp2[i,3],y1=bp2[i,5],
                       angle=90,length=0.015,lwd=1.00);
            }
        }
    }
    poin <- poin + 1;
}
text(x=-1,y=0.0125,cex=2,labels=expression(paste("(I) ",S['2,2'])),pos=4);
par(new = TRUE);
plot(x=fullblock[1],y=bp2[1,3],type="n",xlim=c(0,max(fullblock)+1),ylim=c(0,40),yaxt="n",
     xaxt="n",xlab="",ylab="");
polygon(x=c(0,fullblock[20]+1,fullblock[20]+1,0),y=c(0,0,8,8),border="black",col="grey80");
segments(x0=0.1,y0=4,x1=fullblock[20]+1,y1=4,lwd=1,lty="dotted");
for(i in 1:dim(bp3)[1]){
    for(j in unique(bp3[,1])){
        for(k in unique(bp3[,2])){
            if(bp3[i,1] == j & bp3[i,2] == k){
                lsh <- fullblock[i] - 0.25;
                rsh <- fullblock[i] + 0.25; 
                polygon(x=c(lsh,rsh,rsh,lsh),y=c(0,0,bp3[i,4]*8,
                        bp3[i,4]*8),border=NA,col="black");
            }
        }
    }
}
axis(side=1,at=c(fourblock[3],fourblock[3]+6,fourblock[3]+12,fourblock[3]+18),
     cex.axis=1.5,cex.lab=1.35,labels=c("0.00","0.01","0.02","0.03"));
box();
# --------------------------------------------
mtext(expression(paste("Cost of polyandry (",c[P],")")),
	outer=TRUE,side=1,line=3.75,cex=1.5);

mtext(expression(paste("Mean inbreeding adjustment (",k[adj],")")),
	outer=TRUE,side=2,line=3.25,cex=1.5);

mtext(expression(paste("Mean number of mates")),
	outer=TRUE,side=4,line=3.75,cex=1.5);

dev.off();
# =====================================================================================

# =====================================================================================
# TODO: BUILD FIGURE 6                                                                #
# ====================================================================================#
setEPS(); # postscript below for final publication?
cairo_ps("Nmin1_mn.eps",family="Arial",height=6,width=6);
par(mar=c(1,1,1,1),oma=c(3.5,3.5,0.1,3.5));
# ----------------------- Panel `a':
bp1 <- f.bypr(pol=pol,initial=-1,additional=100);
bp2 <- f.CIva(bypr=bp1);
bp2 <- bporder(bp2);
bp3 <- f.pomn(pol=pol,initial=-1,additional=100);
bp3 <- bporder(bp3);
fourblock <- c(1,1.8,2.6,3.4,4.2);
fullblock <- c(fourblock,fourblock+6,fourblock+12,fourblock+18);
belowpts  <- -0.0045;
plot(x=fullblock[1],y=bp2[1,3],type="n",xlim=c(0,max(fullblock)+1),
     ylim=c(-0.011 + belowpts,0.014),xaxt="n",yaxt="n",xlab="",ylab="");
abline(h=0,lwd=0.9,lty="dotted");
poin <- 21;
for(i in 1:dim(bp2)[1]){
    if(poin == (21 + 5)){
        poin <- 21;
    }
    for(j in unique(bp2[,1])){
        for(k in unique(bp2[,2])){
            if(bp2[i,1] == j & bp2[i,2] == k){
                points(x=fullblock[i],y=bp2[i,3],pch=poin,cex=0.8,bg="black");
                arrows(x0=fullblock[i],x1=fullblock[i],y0=bp2[i,3],y1=bp2[i,4],
                       angle=90,length=0.025,lwd=1.50);
                arrows(x0=fullblock[i],x1=fullblock[i],y0=bp2[i,3],y1=bp2[i,5],
                       angle=90,length=0.025,lwd=1.50);
            }
        }
    }
    poin <- poin + 1;
}
axis(side=2,at=c(-0.01,0,0.01),cex.axis=1.75);
par(new = TRUE);
plot(x=fullblock[1],y=bp2[1,3],type="n",xlim=c(0,max(fullblock)+1),ylim=c(0,45),yaxt="n",
     xaxt="n",xlab="",ylab="");
polygon(x=c(0,fullblock[20]+1,fullblock[20]+1,0),y=c(0,0,8,8),border="black",col="grey80");
segments(x0=0.1,y0=4,x1=fullblock[20]+1,y1=4,lwd=1,lty="dotted");
for(i in 1:dim(bp3)[1]){
    for(j in unique(bp3[,1])){
        for(k in unique(bp3[,2])){
            if(bp3[i,1] == j & bp3[i,2] == k){
                lsh <- fullblock[i] - 0.25;
                rsh <- fullblock[i] + 0.25; 
                polygon(x=c(lsh,rsh,rsh,lsh),y=c(8,8,bp3[i,3]+8,
                        bp3[i,3]+8),border=NA,col="grey70");
                polygon(x=c(lsh,rsh,rsh,lsh),y=c(0,0,bp3[i,4]*8,
                        bp3[i,4]*8),border=NA,col="black");
            }
        }
    }
}
polygon(x=c(0,fullblock[20]+1,fullblock[20]+1,0),y=c(0,0,8,8),border="black",col=NA);
axis(side=4,at=c(8,18,28,38),labels=c("0","10","20","30"),cex.axis=1.5);
par(new = TRUE);
plot(x=fullblock[1],y=bp2[1,3],type="n",xlim=c(0,max(fullblock)+1),
     ylim=c(-0.011 + belowpts,0.014),yaxt="n",xaxt="n",xlab="",ylab="");
for(i in 1:dim(bp2)[1]){
    if(poin == (21 + 5)){
        poin <- 21;
    }
    for(j in unique(bp2[,1])){
        for(k in unique(bp2[,2])){
            if(bp2[i,1] == j & bp2[i,2] == k){
                points(x=fullblock[i],y=bp2[i,3],pch=poin,cex=0.8,bg="black");
                arrows(x0=fullblock[i],x1=fullblock[i],y0=bp2[i,3],y1=bp2[i,4],
                       angle=90,length=0.025,lwd=1.50);
                arrows(x0=fullblock[i],x1=fullblock[i],y0=bp2[i,3],y1=bp2[i,5],
                       angle=90,length=0.025,lwd=1.50);
            }
        }
    }
    poin <- poin + 1;
}
par(new = TRUE);
plot(x=fullblock[1],y=bp2[1,3],type="n",xlim=c(0,max(fullblock)+1),ylim=c(0,45),yaxt="n",
     xaxt="n",xlab="",ylab="");
polygon(x=c(0,fullblock[20]+1,fullblock[20]+1,0),y=c(0,0,8,8),border="black",col="grey80");
segments(x0=0.1,y0=4,x1=fullblock[20]+1,y1=4,lwd=1,lty="dotted");
for(i in 1:dim(bp3)[1]){
    for(j in unique(bp3[,1])){
        for(k in unique(bp3[,2])){
            if(bp3[i,1] == j & bp3[i,2] == k){
                lsh <- fullblock[i] - 0.25;
                rsh <- fullblock[i] + 0.25; 
                polygon(x=c(lsh,rsh,rsh,lsh),y=c(0,0,bp3[i,4]*8,
                        bp3[i,4]*8),border=NA,col="black");
            }
        }
    }
}
axis(side=1,at=c(fourblock[3],fourblock[3]+6,fourblock[3]+12,fourblock[3]+18),
     cex.axis=1.5,cex.lab=1.35,labels=c("0.00","0.01","0.02","0.03"));
box();

mtext(expression(paste("Cost of polyandry (",c[P],")")),
	outer=TRUE,side=1,line=1.75,cex=1.5);

mtext(expression(paste("Mean inbreeding adjustment (",k[adj],")")),
	outer=TRUE,side=2,line=1.75,cex=1.5);

mtext(expression(paste("Mean number of mates")),
	outer=TRUE,side=4,line=1.75,cex=1.5);

dev.off();

# =====================================================================================

# =====================================================================================
# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX #
# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX #
# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX #
# =====================================================================================
# TODO: All of the code below will print the supplemental material figures            #
# =====================================================================================
# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX #
# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX #
# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX #
# =====================================================================================
# =====================================================================================
# TODO: BUILD FIGURE S1-1                                                             #
# ====================================================================================#
setEPS(); # postscript below for final publication?
cairo_ps("Ivalues.eps",family="Arial",height=7,width=7);
par(mfrow=c(3,3),mar=c(0.2,0.2,0.2,0.2),oma=c(6,6,1,1));
# ----------------------- Panel `a':
use <- which(em[,1]==100 & em[,2]==100);
dat <- em[use,];
nBn <- even.Beta(dat[,4]);
fourblock <- c(1,1.8,2.6,3.4,4.2);
fullblock <- c(fourblock,fourblock+6,fourblock+12,fourblock+18);

sdd <- as.data.frame(cbind(dat[,4],dat[,5],dat[,7]));
colnames(sdd) <- c("bs","cs","vl");
sda <- melt(sdd, id = c('bs','cs'));
sda <- sda[,-3];
boxplots.dat <- boxplot(as.numeric(value)~bs+cs,data=sda,type="n",
                xlim=c(1,max(fullblock)),yaxt="n",xaxt="n",
                at=fullblock,lwd=1.5,cex.lab=1.75,ylim=c(-29,14),
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5);
TST <- NULL;
for(j in unique(dat[,5])){
    for(i in unique(dat[,4])){
        sub <- which(dat[,4]==i & dat[,5]==j);
        cis <- simpleboot(dat[sub,7]);
        ttt <- cis[1] * cis[2];
        ttt <- max(c(ttt,0));
        ttt[ttt > 0] <-1;
        TST <- rbind(TST,c(j,i,cis,ttt));
    }
}
TST <- bcsort(TST);
for(i in 1:dim(TST)[1]){
    if(TST[i,5] == 1){
        lsh <- fullblock[i] - 0.25;
        rsh <- fullblock[i] + 0.25; 
        polygon(x=c(lsh,rsh,rsh,lsh),y=c(-40,-40,15,15),border=NA,col="grey70");
    }
}
boxplots.dat <- boxplot(as.numeric(value)~bs+cs,data=sda,
                xlim=c(1,max(fullblock)),yaxt="n",xaxt="n",
                at=fullblock,lwd=1.5,cex.lab=1.75,ylim=c(-29,14),
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5,add=TRUE);
abline(h=0,lwd=0.8,lty="dotted");
polygon(x=c(-1,12,12,-1),y=c(10,10,18,18),border=NA,col="black");
text(x=-0.5,y=11.75,cex=1.75,labels=expression(paste("(A) ",S['100,100'])),pos=4,col="white");
axis(side=2,at=c(-20,-10,0,10),cex.axis=1.5);
# ----------------------- Panel `b':
use <- which(em[,1]==100 & em[,2]==10);
dat <- em[use,];
nBn <- even.Beta(dat[,4]);
fourblock <- c(1,1.8,2.6,3.4,4.2);
fullblock <- c(fourblock,fourblock+6,fourblock+12,fourblock+18);

sdd <- as.data.frame(cbind(dat[,4],dat[,5],dat[,7]));
colnames(sdd) <- c("bs","cs","vl");
sda <- melt(sdd, id = c('bs','cs'));
sda <- sda[,-3];
boxplots.dat <- boxplot(as.numeric(value)~bs+cs,data=sda,type="n",
                xlim=c(1,max(fullblock)),yaxt="n",xaxt="n",
                at=fullblock,lwd=1.5,cex.lab=1.75,ylim=c(-29,14),
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5);
TST <- NULL;
for(j in unique(dat[,5])){
    for(i in unique(dat[,4])){
        sub <- which(dat[,4]==i & dat[,5]==j);
        cis <- simpleboot(dat[sub,7]);
        ttt <- cis[1] * cis[2];
        ttt <- max(c(ttt,0));
        ttt[ttt > 0] <-1;
        TST <- rbind(TST,c(j,i,cis,ttt));
    }
}
TST <- bcsort(TST);
for(i in 1:dim(TST)[1]){
    if(TST[i,5] == 1){
        lsh <- fullblock[i] - 0.25;
        rsh <- fullblock[i] + 0.25; 
        polygon(x=c(lsh,rsh,rsh,lsh),y=c(-40,-40,15,15),border=NA,col="grey70");
    }
}
boxplots.dat <- boxplot(as.numeric(value)~bs+cs,data=sda,
                xlim=c(1,max(fullblock)),yaxt="n",xaxt="n",
                at=fullblock,lwd=1.5,cex.lab=1.75,ylim=c(-29,14),
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5,add=TRUE);
abline(h=0,lwd=0.8,lty="dotted");
polygon(x=c(-1,12,12,-1),y=c(10,10,18,18),border=NA,col="black");
text(x=-0.5,y=11.75,cex=1.75,labels=expression(paste("(B) ",S['100,10'])),pos=4,col="white");
# ----------------------- Panel `c':
use <- which(em[,1]==100 & em[,2]==2);
dat <- em[use,];
nBn <- even.Beta(dat[,4]);
fourblock <- c(1,1.8,2.6,3.4,4.2);
fullblock <- c(fourblock,fourblock+6,fourblock+12,fourblock+18);

sdd <- as.data.frame(cbind(dat[,4],dat[,5],dat[,7]));
colnames(sdd) <- c("bs","cs","vl");
sda <- melt(sdd, id = c('bs','cs'));
sda <- sda[,-3];
boxplots.dat <- boxplot(as.numeric(value)~bs+cs,data=sda,type="n",
                xlim=c(1,max(fullblock)),yaxt="n",xaxt="n",
                at=fullblock,lwd=1.5,cex.lab=1.75,ylim=c(-29,14),
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5);
TST <- NULL;
for(j in unique(dat[,5])){
    for(i in unique(dat[,4])){
        sub <- which(dat[,4]==i & dat[,5]==j);
        cis <- simpleboot(dat[sub,7]);
        ttt <- cis[1] * cis[2];
        ttt <- max(c(ttt,0));
        ttt[ttt > 0] <-1;
        TST <- rbind(TST,c(j,i,cis,ttt));
    }
}
TST <- bcsort(TST);
for(i in 1:dim(TST)[1]){
    if(TST[i,5] == 1){
        lsh <- fullblock[i] - 0.25;
        rsh <- fullblock[i] + 0.25; 
        polygon(x=c(lsh,rsh,rsh,lsh),y=c(-40,-40,15,15),border=NA,col="grey70");
    }
}
boxplots.dat <- boxplot(as.numeric(value)~bs+cs,data=sda,
                xlim=c(1,max(fullblock)),yaxt="n",xaxt="n",
                at=fullblock,lwd=1.5,cex.lab=1.75,ylim=c(-29,14),
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5,add=TRUE);
abline(h=0,lwd=0.8,lty="dotted");
polygon(x=c(-1,12,12,-1),y=c(10,10,18,18),border=NA,col="black");
text(x=-0.5,y=11.75,cex=1.75,labels=expression(paste("(C) ",S['100,2'])),pos=4,col="white");
# ----------------------- Panel `d':
use <- which(em[,1]==10 & em[,2]==100);
dat <- em[use,];
nBn <- even.Beta(dat[,4]);
fourblock <- c(1,1.8,2.6,3.4,4.2);
fullblock <- c(fourblock,fourblock+6,fourblock+12,fourblock+18);

sdd <- as.data.frame(cbind(dat[,4],dat[,5],dat[,7]));
colnames(sdd) <- c("bs","cs","vl");
sda <- melt(sdd, id = c('bs','cs'));
sda <- sda[,-3];
boxplots.dat <- boxplot(as.numeric(value)~bs+cs,data=sda,type="n",
                xlim=c(1,max(fullblock)),yaxt="n",xaxt="n",
                at=fullblock,lwd=1.5,cex.lab=1.75,ylim=c(-29,14),
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5);
TST <- NULL;
for(j in unique(dat[,5])){
    for(i in unique(dat[,4])){
        sub <- which(dat[,4]==i & dat[,5]==j);
        cis <- simpleboot(dat[sub,7]);
        ttt <- cis[1] * cis[2];
        ttt <- max(c(ttt,0));
        ttt[ttt > 0] <-1;
        TST <- rbind(TST,c(j,i,cis,ttt));
    }
}
TST <- bcsort(TST);
for(i in 1:dim(TST)[1]){
    if(TST[i,5] == 1){
        lsh <- fullblock[i] - 0.25;
        rsh <- fullblock[i] + 0.25; 
        polygon(x=c(lsh,rsh,rsh,lsh),y=c(-40,-40,15,15),border=NA,col="grey70");
    }
}
boxplots.dat <- boxplot(as.numeric(value)~bs+cs,data=sda,
                xlim=c(1,max(fullblock)),yaxt="n",xaxt="n",
                at=fullblock,lwd=1.5,cex.lab=1.75,ylim=c(-29,14),
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5,add=TRUE);
abline(h=0,lwd=0.8,lty="dotted");
polygon(x=c(-1,12,12,-1),y=c(10,10,18,18),border=NA,col="black");
text(x=-0.5,y=11.75,cex=1.75,labels=expression(paste("(D) ",S['10,100'])),pos=4,col="white");
axis(side=2,at=c(-20,-10,0,10),cex.axis=1.5);
# ----------------------- Panel `e':
use <- which(em[,1]==10 & em[,2]==10);
dat <- em[use,];
nBn <- even.Beta(dat[,4]);
fourblock <- c(1,1.8,2.6,3.4,4.2);
fullblock <- c(fourblock,fourblock+6,fourblock+12,fourblock+18);

sdd <- as.data.frame(cbind(dat[,4],dat[,5],dat[,7]));
colnames(sdd) <- c("bs","cs","vl");
sda <- melt(sdd, id = c('bs','cs'));
sda <- sda[,-3];
boxplots.dat <- boxplot(as.numeric(value)~bs+cs,data=sda,type="n",
                xlim=c(1,max(fullblock)),yaxt="n",xaxt="n",
                at=fullblock,lwd=1.5,cex.lab=1.75,ylim=c(-29,14),
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5);
TST <- NULL;
for(j in unique(dat[,5])){
    for(i in unique(dat[,4])){
        sub <- which(dat[,4]==i & dat[,5]==j);
        cis <- simpleboot(dat[sub,7]);
        ttt <- cis[1] * cis[2];
        ttt <- max(c(ttt,0));
        ttt[ttt > 0] <-1;
        TST <- rbind(TST,c(j,i,cis,ttt));
    }
}
TST <- bcsort(TST);
for(i in 1:dim(TST)[1]){
    if(TST[i,5] == 1){
        lsh <- fullblock[i] - 0.25;
        rsh <- fullblock[i] + 0.25; 
        polygon(x=c(lsh,rsh,rsh,lsh),y=c(-40,-40,15,15),border=NA,col="grey70");
    }
}
boxplots.dat <- boxplot(as.numeric(value)~bs+cs,data=sda,
                xlim=c(1,max(fullblock)),yaxt="n",xaxt="n",
                at=fullblock,lwd=1.5,cex.lab=1.75,ylim=c(-29,14),
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5,add=TRUE);
abline(h=0,lwd=0.8,lty="dotted");
polygon(x=c(-1,12,12,-1),y=c(10,10,18,18),border=NA,col="black");
text(x=-0.5,y=11.75,cex=1.75,labels=expression(paste("(E) ",S['10,10'])),pos=4,col="white");
# ----------------------- Panel `f':
use <- which(em[,1]==10 & em[,2]==2);
dat <- em[use,];
nBn <- even.Beta(dat[,4]);
fourblock <- c(1,1.8,2.6,3.4,4.2);
fullblock <- c(fourblock,fourblock+6,fourblock+12,fourblock+18);

sdd <- as.data.frame(cbind(dat[,4],dat[,5],dat[,7]));
colnames(sdd) <- c("bs","cs","vl");
sda <- melt(sdd, id = c('bs','cs'));
sda <- sda[,-3];
boxplots.dat <- boxplot(as.numeric(value)~bs+cs,data=sda,type="n",
                xlim=c(1,max(fullblock)),yaxt="n",xaxt="n",
                at=fullblock,lwd=1.5,cex.lab=1.75,ylim=c(-29,14),
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5);
TST <- NULL;
for(j in unique(dat[,5])){
    for(i in unique(dat[,4])){
        sub <- which(dat[,4]==i & dat[,5]==j);
        cis <- simpleboot(dat[sub,7]);
        ttt <- cis[1] * cis[2];
        ttt <- max(c(ttt,0));
        ttt[ttt > 0] <-1;
        TST <- rbind(TST,c(j,i,cis,ttt));
    }
}
TST <- bcsort(TST);
for(i in 1:dim(TST)[1]){
    if(TST[i,5] == 1){
        lsh <- fullblock[i] - 0.25;
        rsh <- fullblock[i] + 0.25; 
        polygon(x=c(lsh,rsh,rsh,lsh),y=c(-40,-40,15,15),border=NA,col="grey70");
    }
}
boxplots.dat <- boxplot(as.numeric(value)~bs+cs,data=sda,
                xlim=c(1,max(fullblock)),yaxt="n",xaxt="n",
                at=fullblock,lwd=1.5,cex.lab=1.75,ylim=c(-29,14),
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5,add=TRUE);
abline(h=0,lwd=0.8,lty="dotted");
polygon(x=c(-1,12,12,-1),y=c(10,10,18,18),border=NA,col="black");
text(x=-0.5,y=11.75,cex=1.75,labels=expression(paste("(F) ",S['10,2'])),pos=4,col="white");
# ----------------------- Panel `g':
use <- which(em[,1]==2 & em[,2]==100);
dat <- em[use,];
nBn <- even.Beta(dat[,4]);
fourblock <- c(1,1.8,2.6,3.4,4.2);
fullblock <- c(fourblock,fourblock+6,fourblock+12,fourblock+18);

sdd <- as.data.frame(cbind(dat[,4],dat[,5],dat[,7]));
colnames(sdd) <- c("bs","cs","vl");
sda <- melt(sdd, id = c('bs','cs'));
sda <- sda[,-3];
boxplots.dat <- boxplot(as.numeric(value)~bs+cs,data=sda,type="n",
                xlim=c(1,max(fullblock)),yaxt="n",xaxt="n",
                at=fullblock,lwd=1.5,cex.lab=1.75,ylim=c(-29,14),
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5);
TST <- NULL;
for(j in unique(dat[,5])){
    for(i in unique(dat[,4])){
        sub <- which(dat[,4]==i & dat[,5]==j);
        cis <- simpleboot(dat[sub,7]);
        ttt <- cis[1] * cis[2];
        ttt <- max(c(ttt,0));
        ttt[ttt > 0] <-1;
        TST <- rbind(TST,c(j,i,cis,ttt));
    }
}
TST <- bcsort(TST);
for(i in 1:dim(TST)[1]){
    if(TST[i,5] == 1){
        lsh <- fullblock[i] - 0.25;
        rsh <- fullblock[i] + 0.25; 
        polygon(x=c(lsh,rsh,rsh,lsh),y=c(-40,-40,15,15),border=NA,col="grey70");
    }
}
boxplots.dat <- boxplot(as.numeric(value)~bs+cs,data=sda,
                xlim=c(1,max(fullblock)),yaxt="n",xaxt="n",
                at=fullblock,lwd=1.5,cex.lab=1.75,ylim=c(-29,14),
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5,add=TRUE);
abline(h=0,lwd=0.8,lty="dotted");
polygon(x=c(-1,12,12,-1),y=c(10,10,18,18),border=NA,col="black");
text(x=-0.5,y=11.75,cex=1.75,labels=expression(paste("(G) ",S['2,100'])),pos=4,col="white");
axis(side=2,at=c(-20,-10,0,10),cex.axis=1.5);
axis(side=1,at=c(fourblock[3],fourblock[3]+6,fourblock[3]+12,fourblock[3]+18),
     cex.axis=1.5,cex.lab=1.5,labels=c("0.00","0.01","0.02","0.03"));
# ----------------------- Panel `h':
use <- which(em[,1]==2 & em[,2]==10);
dat <- em[use,];
nBn <- even.Beta(dat[,4]);
fourblock <- c(1,1.8,2.6,3.4,4.2);
fullblock <- c(fourblock,fourblock+6,fourblock+12,fourblock+18);

sdd <- as.data.frame(cbind(dat[,4],dat[,5],dat[,7]));
colnames(sdd) <- c("bs","cs","vl");
sda <- melt(sdd, id = c('bs','cs'));
sda <- sda[,-3];
boxplots.dat <- boxplot(as.numeric(value)~bs+cs,data=sda,type="n",
                xlim=c(1,max(fullblock)),yaxt="n",xaxt="n",
                at=fullblock,lwd=1.5,cex.lab=1.75,ylim=c(-29,14),
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5);
TST <- NULL;
for(j in unique(dat[,5])){
    for(i in unique(dat[,4])){
        sub <- which(dat[,4]==i & dat[,5]==j);
        cis <- simpleboot(dat[sub,7]);
        ttt <- cis[1] * cis[2];
        ttt <- max(c(ttt,0));
        ttt[ttt > 0] <-1;
        TST <- rbind(TST,c(j,i,cis,ttt));
    }
}
TST <- bcsort(TST);
for(i in 1:dim(TST)[1]){
    if(TST[i,5] == 1){
        lsh <- fullblock[i] - 0.25;
        rsh <- fullblock[i] + 0.25; 
        polygon(x=c(lsh,rsh,rsh,lsh),y=c(-40,-40,15,15),border=NA,col="grey70");
    }
}
boxplots.dat <- boxplot(as.numeric(value)~bs+cs,data=sda,
                xlim=c(1,max(fullblock)),yaxt="n",xaxt="n",
                at=fullblock,lwd=1.5,cex.lab=1.75,ylim=c(-29,14),
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5,add=TRUE);
abline(h=0,lwd=0.8,lty="dotted");
polygon(x=c(-1,12,12,-1),y=c(10,10,18,18),border=NA,col="black");
text(x=-0.5,y=11.75,cex=1.75,labels=expression(paste("(H) ",S['2,10'])),pos=4,col="white");
axis(side=1,at=c(fourblock[3],fourblock[3]+6,fourblock[3]+12,fourblock[3]+18),
     cex.axis=1.5,cex.lab=1.5,labels=c("0.00","0.01","0.02","0.03"));
# ----------------------- Panel `i':
use <- which(em[,1]==2 & em[,2]==2);
dat <- em[use,];
nBn <- even.Beta(dat[,4]);
fourblock <- c(1,1.8,2.6,3.4,4.2);
fullblock <- c(fourblock,fourblock+6,fourblock+12,fourblock+18);

sdd <- as.data.frame(cbind(dat[,4],dat[,5],dat[,7]));
colnames(sdd) <- c("bs","cs","vl");
sda <- melt(sdd, id = c('bs','cs'));
sda <- sda[,-3];
boxplots.dat <- boxplot(as.numeric(value)~bs+cs,data=sda,type="n",
                xlim=c(1,max(fullblock)),yaxt="n",xaxt="n",
                at=fullblock,lwd=1.5,cex.lab=1.75,ylim=c(-29,14),
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5);
TST <- NULL;
for(j in unique(dat[,5])){
    for(i in unique(dat[,4])){
        sub <- which(dat[,4]==i & dat[,5]==j);
        cis <- simpleboot(dat[sub,7]);
        ttt <- cis[1] * cis[2];
        ttt <- max(c(ttt,0));
        ttt[ttt > 0] <-1;
        TST <- rbind(TST,c(j,i,cis,ttt));
    }
}
TST <- bcsort(TST);
for(i in 1:dim(TST)[1]){
    if(TST[i,5] == 1){
        lsh <- fullblock[i] - 0.25;
        rsh <- fullblock[i] + 0.25; 
        polygon(x=c(lsh,rsh,rsh,lsh),y=c(-40,-40,15,15),border=NA,col="grey70");
    }
}
boxplots.dat <- boxplot(as.numeric(value)~bs+cs,data=sda,
                xlim=c(1,max(fullblock)),yaxt="n",xaxt="n",
                at=fullblock,lwd=1.5,cex.lab=1.75,ylim=c(-29,14),
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5,add=TRUE);
abline(h=0,lwd=0.8,lty="dotted");
polygon(x=c(-1,12,12,-1),y=c(10,10,18,18),border=NA,col="black");
text(x=-0.5,y=11.75,cex=1.75,labels=expression(paste("(I) ",S['2,2'])),pos=4,col="white");
axis(side=1,at=c(fourblock[3],fourblock[3]+6,fourblock[3]+12,fourblock[3]+18),
     cex.axis=1.5,cex.lab=1.5,labels=c("0.00","0.01","0.02","0.03"));

mtext(expression(paste("Cost of polyandry (",c[P],")")),
	outer=TRUE,side=1,line=3.75,cex=1.5);

mtext(expression(paste("Inbreeding mean allele value (",I[a],")")),
	outer=TRUE,side=2,line=3.25,cex=1.5);

dev.off();

# =====================================================================================

# =====================================================================================
# TODO: BUILD FIGURE S1-2                                                             #
# ====================================================================================#
setEPS(); # postscript below for final publication?
cairo_ps("Ivalues_Nmin1.eps",family="Arial",height=7,width=7);
par(mar=c(5,5,2,2));
use <- which(em[,1]==-1 & em[,2]==100);
dat <- em[use,];
nBn <- even.Beta(dat[,4]);
fourblock <- c(1,1.8,2.6,3.4,4.2);
fullblock <- c(fourblock,fourblock+6,fourblock+12,fourblock+18);

sdd <- as.data.frame(cbind(dat[,4],dat[,5],dat[,7]));
colnames(sdd) <- c("bs","cs","vl");
sda <- melt(sdd, id = c('bs','cs'));
sda <- sda[,-3];
boxplots.dat <- boxplot(as.numeric(value)~bs+cs,data=sda,type="n",
                xlim=c(1,max(fullblock)),yaxt="n",xaxt="n",
                at=fullblock,lwd=1.5,cex.lab=1.75,ylim=c(-29,14),
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5,
                xlab=expression(paste("Cost of polyandry (",c[P],")")),
                ylab=expression(paste("Inbreeding mean allele value (",I[a],")")));
TST <- NULL;
for(j in unique(dat[,5])){
    for(i in unique(dat[,4])){
        sub <- which(dat[,4]==i & dat[,5]==j);
        cis <- simpleboot(dat[sub,7]);
        ttt <- cis[1] * cis[2];
        ttt <- max(c(ttt,0));
        ttt[ttt > 0] <-1;
        TST <- rbind(TST,c(j,i,cis,ttt));
    }
}
TST <- bcsort(TST);
for(i in 1:dim(TST)[1]){
    if(TST[i,5] == 1){
        lsh <- fullblock[i] - 0.25;
        rsh <- fullblock[i] + 0.25; 
        polygon(x=c(lsh,rsh,rsh,lsh),y=c(-40,-40,15,15),border=NA,col="grey70");
    }
}
boxplots.dat <- boxplot(as.numeric(value)~bs+cs,data=sda,
                xlim=c(1,max(fullblock)),yaxt="n",xaxt="n",
                at=fullblock,lwd=1.5,cex.lab=1.75,ylim=c(-29,14),
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5,add=TRUE);
abline(h=0,lwd=0.8,lty="dotted");
axis(side=2,at=c(-20,-10,0,10),cex.axis=1.5);
axis(side=1,at=c(fourblock[3],fourblock[3]+6,fourblock[3]+12,fourblock[3]+18),
     cex.axis=1.5,cex.lab=1.5,labels=c("0.00","0.01","0.02","0.03"));

dev.off();

# =====================================================================================

# =====================================================================================
# TODO: BUILD FIGURE S1-3                                                             #
# ====================================================================================#
setEPS(); # postscript below for final publication?
cairo_ps("FullFact_gI.eps",family="Arial",height=7,width=7);
par(mfrow=c(3,3),mar=c(0.2,0.2,0.2,0.2),oma=c(6,6,1,1));
# ----------------------- Panel `a':
use <- which(summ[,1] == 100 & summ[,2] == 100);
aus <- which(alle[,1] == 100 & alle[,2] == 100);
dat <- as.data.frame(cbind(summ[use,3],summ[use,4],summ[use,6]));
colnames(dat) <- c("bs","cs","vl");
sda <- melt(dat, id = c('bs','cs'));
sda <- sda[,-3];
ale <- alle[aus,];
fourblock <- c(1,1.8,2.6,3.4,4.2);
fullblock <- c(fourblock,fourblock+6,fourblock+12,fourblock+18);
boxplots.dat <- boxplot(as.numeric(value)~bs+cs,data=sda,
                at=fullblock,lwd=1.5,xaxt="n",yaxt="n",
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5,ylim=c(-520,250));
axis(side=2,at=c(-400,-200,0,200),cex.axis=2,cex.lab=2);
for(i in 1:dim(ale)[1]){
    if(ale[i,6] > 0 & ale[i,7] > 0){
        lsh <- fullblock[i] - 0.25;
        rsh <- fullblock[i] + 0.25;
        #polygon(x=c(lsh,rsh,rsh,lsh),y=c(-10,-10,350,350),border=NA,col="grey70");
    }
}
boxplots.dat <- boxplot(as.numeric(value)~bs+cs,data=sda,
                at=fullblock,lwd=1.5,add=TRUE,xaxt="n",yaxt="n",
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5,ylim=c(-520,250));
abline(h=0,lty="dotted",lwd=0.8);
text(x=9.5,y=200,cex=2,labels=expression(paste("(A) ",S['100,100'])),pos=4);
# ----------------------- Panel `b':
use <- which(summ[,1] == 100 & summ[,2] == 10);
aus <- which(alle[,1] == 100 & alle[,2] == 10);
dat <- as.data.frame(cbind(summ[use,3],summ[use,4],summ[use,6]));
colnames(dat) <- c("bs","cs","vl");
sda <- melt(dat, id = c('bs','cs'));
sda <- sda[,-3];
ale <- alle[aus,];
fourblock <- c(1,1.8,2.6,3.4,4.2);
fullblock <- c(fourblock,fourblock+6,fourblock+12,fourblock+18);
boxplots.dat <- boxplot(as.numeric(value)~bs+cs,data=sda,
                at=fullblock,lwd=1.5,xaxt="n",yaxt="n",
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5,ylim=c(-520,250));
for(i in 1:dim(ale)[1]){
    if(ale[i,6] > 0 & ale[i,7] > 0){
        lsh <- fullblock[i] - 0.25;
        rsh <- fullblock[i] + 0.25;
        #polygon(x=c(lsh,rsh,rsh,lsh),y=c(-10,-10,350,350),border=NA,col="grey70");
    }
}
boxplots.dat <- boxplot(as.numeric(value)~bs+cs,data=sda,
                at=fullblock,lwd=1.5,add=TRUE,xaxt="n",yaxt="n",
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5,ylim=c(-520,250));
abline(h=0,lty="dotted",lwd=0.8);
text(x=10,y=200,cex=2,labels=expression(paste("(B) ",S['100,10'])),pos=4);
# ----------------------- Panel `c':
use <- which(summ[,1] == 100 & summ[,2] == 2);
aus <- which(alle[,1] == 100 & alle[,2] == 2);
dat <- as.data.frame(cbind(summ[use,3],summ[use,4],summ[use,6]));
colnames(dat) <- c("bs","cs","vl");
sda <- melt(dat, id = c('bs','cs'));
sda <- sda[,-3];
ale <- alle[aus,];
fourblock <- c(1,1.8,2.6,3.4,4.2);
fullblock <- c(fourblock,fourblock+6,fourblock+12,fourblock+18);
boxplots.dat <- boxplot(as.numeric(value)~bs+cs,data=sda,
                at=fullblock,lwd=1.5,xaxt="n",yaxt="n",
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5,ylim=c(-520,250));
for(i in 1:dim(ale)[1]){
    if(ale[i,6] > 0 & ale[i,7] > 0){
        lsh <- fullblock[i] - 0.25;
        rsh <- fullblock[i] + 0.25;
        #polygon(x=c(lsh,rsh,rsh,lsh),y=c(-10,-10,350,350),border=NA,col="grey70");
    }
}
boxplots.dat <- boxplot(as.numeric(value)~bs+cs,data=sda,
                at=fullblock,lwd=1.5,add=TRUE,xaxt="n",yaxt="n",
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5,ylim=c(-520,250));
abline(h=0,lty="dotted",lwd=0.8);
text(x=10,y=200,cex=2,labels=expression(paste("(C) ",S['100,2'])),pos=4);
# ----------------------- Panel `d':
use <- which(summ[,1] == 10 & summ[,2] == 100);
aus <- which(alle[,1] == 10 & alle[,2] == 100);
dat <- as.data.frame(cbind(summ[use,3],summ[use,4],summ[use,6]));
colnames(dat) <- c("bs","cs","vl");
sda <- melt(dat, id = c('bs','cs'));
sda <- sda[,-3];
ale <- alle[aus,];
fourblock <- c(1,1.8,2.6,3.4,4.2);
fullblock <- c(fourblock,fourblock+6,fourblock+12,fourblock+18);
boxplots.dat <- boxplot(as.numeric(value)~bs+cs,data=sda,
                at=fullblock,lwd=1.5,xaxt="n",yaxt="n",
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5,ylim=c(-520,250));
axis(side=2,at=c(-400,-200,0,200),cex.axis=2,cex.lab=2);
for(i in 1:dim(ale)[1]){
    if(ale[i,6] > 0 & ale[i,7] > 0){
        lsh <- fullblock[i] - 0.25;
        rsh <- fullblock[i] + 0.25;
        #polygon(x=c(lsh,rsh,rsh,lsh),y=c(-10,-10,350,350),border=NA,col="grey70");
    }
}
boxplots.dat <- boxplot(as.numeric(value)~bs+cs,data=sda,
                at=fullblock,lwd=1.5,add=TRUE,xaxt="n",yaxt="n",
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5,ylim=c(-520,250));
abline(h=0,lty="dotted",lwd=0.8);
text(x=10,y=200,cex=2,labels=expression(paste("(D) ",S['10,100'])),pos=4);
# ----------------------- Panel `e':
use <- which(summ[,1] == 10 & summ[,2] == 10);
aus <- which(alle[,1] == 10 & alle[,2] == 10);
dat <- as.data.frame(cbind(summ[use,3],summ[use,4],summ[use,6]));
colnames(dat) <- c("bs","cs","vl");
sda <- melt(dat, id = c('bs','cs'));
sda <- sda[,-3];
ale <- alle[aus,];
fourblock <- c(1,1.8,2.6,3.4,4.2);
fullblock <- c(fourblock,fourblock+6,fourblock+12,fourblock+18);
boxplots.dat <- boxplot(as.numeric(value)~bs+cs,data=sda,
                at=fullblock,lwd=1.5,xaxt="n",yaxt="n",
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5,ylim=c(-520,250));
for(i in 1:dim(ale)[1]){
    if(ale[i,6] > 0 & ale[i,7] > 0){
        lsh <- fullblock[i] - 0.25;
        rsh <- fullblock[i] + 0.25;
        #polygon(x=c(lsh,rsh,rsh,lsh),y=c(-10,-10,350,350),border=NA,col="grey70");
    }
}
boxplots.dat <- boxplot(as.numeric(value)~bs+cs,data=sda,
                at=fullblock,lwd=1.5,add=TRUE,xaxt="n",yaxt="n",
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5,ylim=c(-520,250));
abline(h=0,lty="dotted",lwd=0.8);
text(x=10,y=200,cex=2,labels=expression(paste("(E) ",S['10,10'])),pos=4);
# ----------------------- Panel `f':
use <- which(summ[,1] == 10 & summ[,2] == 2);
aus <- which(alle[,1] == 10 & alle[,2] == 2);
dat <- as.data.frame(cbind(summ[use,3],summ[use,4],summ[use,6]));
colnames(dat) <- c("bs","cs","vl");
sda <- melt(dat, id = c('bs','cs'));
sda <- sda[,-3];
ale <- alle[aus,];
fourblock <- c(1,1.8,2.6,3.4,4.2);
fullblock <- c(fourblock,fourblock+6,fourblock+12,fourblock+18);
boxplots.dat <- boxplot(as.numeric(value)~bs+cs,data=sda,
                at=fullblock,lwd=1.5,xaxt="n",yaxt="n",
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5,ylim=c(-520,250));
for(i in 1:dim(ale)[1]){
    if(ale[i,6] > 0 & ale[i,7] > 0){
        lsh <- fullblock[i] - 0.25;
        rsh <- fullblock[i] + 0.25;
        #polygon(x=c(lsh,rsh,rsh,lsh),y=c(-10,-10,350,350),border=NA,col="grey70");
    }
}
boxplots.dat <- boxplot(as.numeric(value)~bs+cs,data=sda,
                at=fullblock,lwd=1.5,add=TRUE,xaxt="n",yaxt="n",
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5,ylim=c(-520,250));
abline(h=0,lty="dotted",lwd=0.8);
text(x=10,y=200,cex=2,labels=expression(paste("(F) ",S['10,2'])),pos=4);
# ----------------------- Panel `g':
use <- which(summ[,1] == 2 & summ[,2] == 100);
aus <- which(alle[,1] == 2 & alle[,2] == 100);
dat <- as.data.frame(cbind(summ[use,3],summ[use,4],summ[use,6]));
colnames(dat) <- c("bs","cs","vl");
sda <- melt(dat, id = c('bs','cs'));
sda <- sda[,-3];
ale <- alle[aus,];
fourblock <- c(1,1.8,2.6,3.4,4.2);
fullblock <- c(fourblock,fourblock+6,fourblock+12,fourblock+18);
boxplots.dat <- boxplot(as.numeric(value)~bs+cs,data=sda,
                at=fullblock,lwd=1.5,xaxt="n",yaxt="n",
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5,ylim=c(-520,250));
axis(side=2,at=c(-400,-200,0,200),cex.axis=2,cex.lab=2);
axis(side=1,at=c(fourblock[3],fourblock[3]+6,fourblock[3]+12,fourblock[3]+18),
     cex.axis=1.5,cex.lab=1.5,labels=c("0.00","0.01","0.02","0.03"));
for(i in 1:dim(ale)[1]){
    if(ale[i,6] > 0 & ale[i,7] > 0){
        lsh <- fullblock[i] - 0.25;
        rsh <- fullblock[i] + 0.25;
        #polygon(x=c(lsh,rsh,rsh,lsh),y=c(-10,-10,350,350),border=NA,col="grey70");
    }
}
boxplots.dat <- boxplot(as.numeric(value)~bs+cs,data=sda,
                at=fullblock,lwd=1.5,add=TRUE,xaxt="n",yaxt="n",
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5,ylim=c(-520,250));
abline(h=0,lty="dotted",lwd=0.8);
text(x=10,y=200,cex=2,labels=expression(paste("(G) ",S['2,100'])),pos=4);
# ----------------------- Panel `h':
use <- which(summ[,1] == 2 & summ[,2] == 10);
aus <- which(alle[,1] == 2 & alle[,2] == 10);
dat <- as.data.frame(cbind(summ[use,3],summ[use,4],summ[use,6]));
colnames(dat) <- c("bs","cs","vl");
sda <- melt(dat, id = c('bs','cs'));
sda <- sda[,-3];
ale <- alle[aus,];
fourblock <- c(1,1.8,2.6,3.4,4.2);
fullblock <- c(fourblock,fourblock+6,fourblock+12,fourblock+18);
boxplots.dat <- boxplot(as.numeric(value)~bs+cs,data=sda,
                at=fullblock,lwd=1.5,xaxt="n",yaxt="n",
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5,ylim=c(-520,250));
axis(side=1,at=c(fourblock[3],fourblock[3]+6,fourblock[3]+12,fourblock[3]+18),
     cex.axis=1.5,cex.lab=1.5,labels=c("0.00","0.01","0.02","0.03"));
for(i in 1:dim(ale)[1]){
    if(ale[i,6] > 0 & ale[i,7] > 0){
        lsh <- fullblock[i] - 0.25;
        rsh <- fullblock[i] + 0.25;
        #polygon(x=c(lsh,rsh,rsh,lsh),y=c(-10,-10,350,350),border=NA,col="grey70");
    }
}
boxplots.dat <- boxplot(as.numeric(value)~bs+cs,data=sda,
                at=fullblock,lwd=1.5,add=TRUE,xaxt="n",yaxt="n",
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5,ylim=c(-520,250),pos=4);
abline(h=0,lty="dotted",lwd=0.8);
text(x=10,y=200,cex=2,labels=expression(paste("(H) ",S['2,10'])),pos=4);
# ----------------------- Panel `i':
use <- which(summ[,1] == 2 & summ[,2] == 2);
aus <- which(alle[,1] == 2 & alle[,2] == 2);
dat <- as.data.frame(cbind(summ[use,3],summ[use,4],summ[use,6]));
colnames(dat) <- c("bs","cs","vl");
sda <- melt(dat, id = c('bs','cs'));
sda <- sda[,-3];
ale <- alle[aus,];
fourblock <- c(1,1.8,2.6,3.4,4.2);
fullblock <- c(fourblock,fourblock+6,fourblock+12,fourblock+18);
boxplots.dat <- boxplot(as.numeric(value)~bs+cs,data=sda,
                at=fullblock,lwd=1.5,xaxt="n",yaxt="n",
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5,ylim=c(-520,250));
axis(side=1,at=c(fourblock[3],fourblock[3]+6,fourblock[3]+12,fourblock[3]+18),
     cex.axis=1.5,cex.lab=1.5,labels=c("0.00","0.01","0.02","0.03"));
for(i in 1:dim(ale)[1]){
    if(ale[i,6] > 0 & ale[i,7] > 0){
        lsh <- fullblock[i] - 0.25;
        rsh <- fullblock[i] + 0.25;
        #polygon(x=c(lsh,rsh,rsh,lsh),y=c(-10,-10,350,350),border=NA,col="grey70");
    }
}
boxplots.dat <- boxplot(as.numeric(value)~bs+cs,data=sda,
                at=fullblock,lwd=1.5,add=TRUE,xaxt="n",yaxt="n",
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5,ylim=c(-520,250));
abline(h=0,lty="dotted",lwd=0.8);
text(x=10,y=200,cex=2,labels=expression(paste("(I) ",S['2,2'])),pos=4);

mtext(expression(paste("Cost of polyandry (",c[P],")")),
	outer=TRUE,side=1,line=3.75,cex=2);

mtext(expression(paste("Mean inbreeding strategy phenotype value (",I[p],")")),
	outer=TRUE,side=2,line=3.25,cex=1.5);

dev.off();

# =====================================================================================

# =====================================================================================
# TODO: BUILD FIGURE S1-4                                                             #
# ====================================================================================#
setEPS(); # postscript below for final publication?
cairo_ps("gI_vals_Nmin1.eps",family="Arial",height=7,width=7);
par(mar=c(0.2,0.2,0.2,0.2),oma=c(6,6,1,1),lwd=2);
# ----------------------- Panel `a':
use <- which(summ[,1] == -1 & summ[,2] == 100);
aus <- which(alle[,1] == -1 & alle[,2] == 100);
dat <- as.data.frame(cbind(summ[use,3],summ[use,4],summ[use,6]));
colnames(dat) <- c("bs","cs","vl");
sda <- melt(dat, id = c('bs','cs'));
sda <- sda[,-3];
ale <- alle[aus,];
fourblock <- c(1,1.8,2.6,3.4,4.2);
fullblock <- c(fourblock,fourblock+6,fourblock+12,fourblock+18);
boxplots.dat <- boxplot(as.numeric(value)~bs+cs,data=sda,
                at=fullblock,lwd=2.25,xaxt="n",yaxt="n",
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5,ylim=c(-450,230));
axis(side=2,at=c(-400,-200,0,200),cex.axis=2,cex.lab=2);
for(i in 1:dim(ale)[1]){
    if(ale[i,6] > 0 & ale[i,7] > 0){
        lsh <- fullblock[i] - 0.25;
        rsh <- fullblock[i] + 0.25;
        #polygon(x=c(lsh,rsh,rsh,lsh),y=c(-500,-500,350,350),border=NA,col="grey70");
    }
}
boxplots.dat <- boxplot(as.numeric(value)~bs+cs,data=sda,
                at=fullblock,lwd=2.25,add=TRUE,xaxt="n",yaxt="n",
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5,ylim=c(0,240));

abline(h=0,lty="dotted",lwd=0.8);
axis(side=1,at=c(fourblock[3],fourblock[3]+6,fourblock[3]+12,fourblock[3]+18),
     cex.axis=2,cex.lab=2,labels=c("0.00","0.01","0.02","0.03"));

mtext(expression(paste("Cost of polyandry (",c[P],")")),
	outer=TRUE,side=1,line=3.75,cex=2);

mtext(expression(paste("Mean inbreeding strategy phenotype value (",I[p],")")),
	outer=TRUE,side=2,line=3.25,cex=1.5);

dev.off();

# =====================================================================================

# =====================================================================================
# TODO: BUILD FIGURE S1-5                                                             #
# ====================================================================================#

setEPS(); # postscript below for final publication?
cairo_ps("fvalues.eps",family="Arial",height=7,width=7);
par(mfrow=c(3,3),mar=c(0.2,0.2,0.2,0.2),oma=c(6,6,1,1));
# ----------------------- Panel `a':
use <- which(em[,1]==100 & em[,2]==100);
dat <- em[use,];
nBn <- even.Beta(dat[,4]);
fourblock <- c(1,1.8,2.6,3.4,4.2);
fullblock <- c(fourblock,fourblock+6,fourblock+12,fourblock+18);

sdd <- as.data.frame(cbind(dat[,4],dat[,5],dat[,13]));
colnames(sdd) <- c("bs","cs","vl");
sda <- melt(sdd, id = c('bs','cs'));
sda <- sda[,-3];
boxplots.dat <- boxplot(as.numeric(value)~bs+cs,data=sda,type="n",
                xlim=c(1,max(fullblock)),yaxt="n",xaxt="n",
                at=fullblock,lwd=1.5,cex.lab=1.75,ylim=c(0,1.2),
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5);
axis(side=2,at=c(0,0.25,0.5,0.75,1.0),cex.axis=1.5);
text(x=0,y=1.1,cex=1.75,labels=expression(paste("(A) ",S['100,100'])),pos=4);
# ----------------------- Panel `b':
use <- which(em[,1]==100 & em[,2]==10);
dat <- em[use,];
nBn <- even.Beta(dat[,4]);
fourblock <- c(1,1.8,2.6,3.4,4.2);
fullblock <- c(fourblock,fourblock+6,fourblock+12,fourblock+18);

sdd <- as.data.frame(cbind(dat[,4],dat[,5],dat[,13]));
colnames(sdd) <- c("bs","cs","vl");
sda <- melt(sdd, id = c('bs','cs'));
sda <- sda[,-3];
boxplots.dat <- boxplot(as.numeric(value)~bs+cs,data=sda,type="n",
                xlim=c(1,max(fullblock)),yaxt="n",xaxt="n",
                at=fullblock,lwd=1.5,cex.lab=1.75,ylim=c(0,1.2),
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5);
text(x=0,y=1.1,cex=1.75,labels=expression(paste("(B) ",S['100,10'])),pos=4);
# ----------------------- Panel `c':
use <- which(em[,1]==100 & em[,2]==2);
dat <- em[use,];
nBn <- even.Beta(dat[,4]);
fourblock <- c(1,1.8,2.6,3.4,4.2);
fullblock <- c(fourblock,fourblock+6,fourblock+12,fourblock+18);

sdd <- as.data.frame(cbind(dat[,4],dat[,5],dat[,13]));
colnames(sdd) <- c("bs","cs","vl");
sda <- melt(sdd, id = c('bs','cs'));
sda <- sda[,-3];
boxplots.dat <- boxplot(as.numeric(value)~bs+cs,data=sda,type="n",
                xlim=c(1,max(fullblock)),yaxt="n",xaxt="n",
                at=fullblock,lwd=1.5,cex.lab=1.75,ylim=c(0,1.2),
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5);
text(x=0,y=1.1,cex=1.75,labels=expression(paste("(C) ",S['100,2'])),pos=4);
# ----------------------- Panel `d':
use <- which(em[,1]==10 & em[,2]==100);
dat <- em[use,];
nBn <- even.Beta(dat[,4]);
fourblock <- c(1,1.8,2.6,3.4,4.2);
fullblock <- c(fourblock,fourblock+6,fourblock+12,fourblock+18);

sdd <- as.data.frame(cbind(dat[,4],dat[,5],dat[,13]));
colnames(sdd) <- c("bs","cs","vl");
sda <- melt(sdd, id = c('bs','cs'));
sda <- sda[,-3];
boxplots.dat <- boxplot(as.numeric(value)~bs+cs,data=sda,type="n",
                xlim=c(1,max(fullblock)),yaxt="n",xaxt="n",
                at=fullblock,lwd=1.5,cex.lab=1.75,ylim=c(0,1.2),
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5);
text(x=0,y=1.1,cex=1.75,labels=expression(paste("(D) ",S['10,100'])),pos=4);
axis(side=2,at=c(0,0.25,0.5,0.75,1.0),cex.axis=1.5);
# ----------------------- Panel `e':
use <- which(em[,1]==10 & em[,2]==10);
dat <- em[use,];
nBn <- even.Beta(dat[,4]);
fourblock <- c(1,1.8,2.6,3.4,4.2);
fullblock <- c(fourblock,fourblock+6,fourblock+12,fourblock+18);

sdd <- as.data.frame(cbind(dat[,4],dat[,5],dat[,13]));
colnames(sdd) <- c("bs","cs","vl");
sda <- melt(sdd, id = c('bs','cs'));
sda <- sda[,-3];
boxplots.dat <- boxplot(as.numeric(value)~bs+cs,data=sda,type="n",
                xlim=c(1,max(fullblock)),yaxt="n",xaxt="n",
                at=fullblock,lwd=1.5,cex.lab=1.75,ylim=c(0,1.2),
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5);
text(x=0,y=1.1,cex=1.75,labels=expression(paste("(E) ",S['10,10'])),pos=4);
# ----------------------- Panel `f':
use <- which(em[,1]==10 & em[,2]==2);
dat <- em[use,];
nBn <- even.Beta(dat[,4]);
fourblock <- c(1,1.8,2.6,3.4,4.2);
fullblock <- c(fourblock,fourblock+6,fourblock+12,fourblock+18);

sdd <- as.data.frame(cbind(dat[,4],dat[,5],dat[,13]));
colnames(sdd) <- c("bs","cs","vl");
sda <- melt(sdd, id = c('bs','cs'));
sda <- sda[,-3];
boxplots.dat <- boxplot(as.numeric(value)~bs+cs,data=sda,type="n",
                xlim=c(1,max(fullblock)),yaxt="n",xaxt="n",
                at=fullblock,lwd=1.5,cex.lab=1.75,ylim=c(0,1.2),
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5);
text(x=0,y=1.1,cex=1.75,labels=expression(paste("(F) ",S['10,2'])),pos=4);
# ----------------------- Panel `g':
use <- which(em[,1]==2 & em[,2]==100);
dat <- em[use,];
nBn <- even.Beta(dat[,4]);
fourblock <- c(1,1.8,2.6,3.4,4.2);
fullblock <- c(fourblock,fourblock+6,fourblock+12,fourblock+18);

sdd <- as.data.frame(cbind(dat[,4],dat[,5],dat[,13]));
colnames(sdd) <- c("bs","cs","vl");
sda <- melt(sdd, id = c('bs','cs'));
sda <- sda[,-3];
boxplots.dat <- boxplot(as.numeric(value)~bs+cs,data=sda,type="n",
                xlim=c(1,max(fullblock)),yaxt="n",xaxt="n",
                at=fullblock,lwd=1.5,cex.lab=1.75,ylim=c(0,1.2),
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5);
text(x=0,y=1.1,cex=1.75,labels=expression(paste("(G) ",S['2,100'])),pos=4);
axis(side=2,at=c(0,0.25,0.5,0.75,1.0),cex.axis=1.5);
axis(side=1,at=c(fourblock[3],fourblock[3]+6,fourblock[3]+12,fourblock[3]+18),
     cex.axis=1.5,cex.lab=1.5,labels=c("0.00","0.01","0.02","0.03"));
# ----------------------- Panel `h':
use <- which(em[,1]==2 & em[,2]==10);
dat <- em[use,];
nBn <- even.Beta(dat[,4]);
fourblock <- c(1,1.8,2.6,3.4,4.2);
fullblock <- c(fourblock,fourblock+6,fourblock+12,fourblock+18);

sdd <- as.data.frame(cbind(dat[,4],dat[,5],dat[,13]));
colnames(sdd) <- c("bs","cs","vl");
sda <- melt(sdd, id = c('bs','cs'));
sda <- sda[,-3];
boxplots.dat <- boxplot(as.numeric(value)~bs+cs,data=sda,type="n",
                xlim=c(1,max(fullblock)),yaxt="n",xaxt="n",
                at=fullblock,lwd=1.5,cex.lab=1.75,ylim=c(0,1.2),
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5);
text(x=0,y=1.1,cex=1.75,labels=expression(paste("(H) ",S['2,10'])),pos=4);
axis(side=1,at=c(fourblock[3],fourblock[3]+6,fourblock[3]+12,fourblock[3]+18),
     cex.axis=1.5,cex.lab=1.5,labels=c("0.00","0.01","0.02","0.03"));
# ----------------------- Panel `i':
use <- which(em[,1]==2 & em[,2]==2);
dat <- em[use,];
nBn <- even.Beta(dat[,4]);
fourblock <- c(1,1.8,2.6,3.4,4.2);
fullblock <- c(fourblock,fourblock+6,fourblock+12,fourblock+18);

sdd <- as.data.frame(cbind(dat[,4],dat[,5],dat[,13]));
colnames(sdd) <- c("bs","cs","vl");
sda <- melt(sdd, id = c('bs','cs'));
sda <- sda[,-3];
boxplots.dat <- boxplot(as.numeric(value)~bs+cs,data=sda,type="n",
                xlim=c(1,max(fullblock)),yaxt="n",xaxt="n",
                at=fullblock,lwd=1.5,cex.lab=1.75,ylim=c(0,1.2),
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5);

text(x=0,y=1.1,cex=1.75,labels=expression(paste("(I) ",S['2,2'])),pos=4);
axis(side=1,at=c(fourblock[3],fourblock[3]+6,fourblock[3]+12,fourblock[3]+18),
     cex.axis=1.5,cex.lab=1.5,labels=c("0.00","0.01","0.02","0.03"));

mtext(expression(paste("Cost of polyandry (",c[P],")")),
	outer=TRUE,side=1,line=3.75,cex=1.5);

mtext(expression(paste("Mean inbreeding coefficient value (f)")),
	outer=TRUE,side=2,line=3.25,cex=1.5);

dev.off();

# =====================================================================================

# =====================================================================================
# TODO: BUILD FIGURE S1-6                                                             #
# ====================================================================================#

setEPS(); # postscript below for final publication?
cairo_ps("fvalues_Nmin1.eps",family="Arial",height=7,width=7);
par(mar=c(5,5,1,1));

use <- which(em[,1]==100 & em[,2]==100);
dat <- em[use,];
nBn <- even.Beta(dat[,4]);
fourblock <- c(1,1.8,2.6,3.4,4.2);
fullblock <- c(fourblock,fourblock+6,fourblock+12,fourblock+18);

sdd <- as.data.frame(cbind(dat[,4],dat[,5],dat[,13]));
colnames(sdd) <- c("bs","cs","vl");
sda <- melt(sdd, id = c('bs','cs'));
sda <- sda[,-3];
boxplots.dat <- boxplot(as.numeric(value)~bs+cs,data=sda,type="n",
                xlim=c(1,max(fullblock)),yaxt="n",xaxt="n",
                at=fullblock,lwd=1.5,cex.lab=1.75,ylim=c(0,1.2),
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5,
                xlab=expression(paste("Cost of polyandry (",c[P],")")),
                ylab=expression(paste("Mean inbreeding coefficient value (f)")));
axis(side=2,at=c(0,0.25,0.5,0.75,1.0),cex.axis=1.5);
axis(side=1,at=c(fourblock[3],fourblock[3]+6,fourblock[3]+12,fourblock[3]+18),
     cex.axis=1.5,cex.lab=1.5,labels=c("0.00","0.01","0.02","0.03"));
dev.off();

# =====================================================================================

# =====================================================================================
# TODO: BUILD FIGURE S1-9                                                             #
# ====================================================================================#

setEPS(); # postscript below for final publication?
cairo_ps("Nvalues.eps",family="Arial",height=7,width=7);
par(mfrow=c(3,3),mar=c(0.2,0.2,0.2,0.2),oma=c(6,6,1,1));
# ----------------------- Panel `a':
use <- which(em[,1]==100 & em[,2]==100);
dat <- em[use,];
nBn <- even.Beta(dat[,4]);
fourblock <- c(1,1.8,2.6,3.4,4.2);
fullblock <- c(fourblock,fourblock+6,fourblock+12,fourblock+18);

sdd <- as.data.frame(cbind(dat[,4],dat[,5],dat[,9]));
colnames(sdd) <- c("bs","cs","vl");
sda <- melt(sdd, id = c('bs','cs'));
sda <- sda[,-3];
boxplots.dat <- boxplot(as.numeric(value)~bs+cs,data=sda,type="n",
                xlim=c(1,max(fullblock)),yaxt="n",xaxt="n",
                at=fullblock,lwd=1.5,cex.lab=1.75,ylim=c(-20,20),
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5);
TST <- NULL;
for(j in unique(dat[,5])){
    for(i in unique(dat[,4])){
        sub <- which(dat[,4]==i & dat[,5]==j);
        cis <- simpleboot(dat[sub,9]);
        ttt <- cis[1] * cis[2];
        ttt <- max(c(ttt,0));
        ttt[ttt > 0] <-1;
        TST <- rbind(TST,c(j,i,cis,ttt));
    }
}
for(i in 1:dim(TST)[1]){
    if(TST[i,5] == 1){
        lsh <- fullblock[i] - 0.25;
        rsh <- fullblock[i] + 0.25; 
    }
}
boxplots.dat <- boxplot(as.numeric(value)~bs+cs,data=sda,
                xlim=c(1,max(fullblock)),yaxt="n",xaxt="n",
                at=fullblock,lwd=1.5,cex.lab=1.75,ylim=c(-20,20),
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5,add=TRUE);
abline(h=0,lwd=0.8,lty="dotted");
text(x=0,y=18,cex=1.75,labels=expression(paste("(A) ",S['100,100'])),pos=4);
axis(side=2,at=c(-20,-10,0,10),cex.axis=1.5);
# ----------------------- Panel `b':
use <- which(em[,1]==100 & em[,2]==10);
dat <- em[use,];
nBn <- even.Beta(dat[,4]);
fourblock <- c(1,1.8,2.6,3.4,4.2);
fullblock <- c(fourblock,fourblock+6,fourblock+12,fourblock+18);

sdd <- as.data.frame(cbind(dat[,4],dat[,5],dat[,9]));
colnames(sdd) <- c("bs","cs","vl");
sda <- melt(sdd, id = c('bs','cs'));
sda <- sda[,-3];
boxplots.dat <- boxplot(as.numeric(value)~bs+cs,data=sda,type="n",
                xlim=c(1,max(fullblock)),yaxt="n",xaxt="n",
                at=fullblock,lwd=1.5,cex.lab=1.75,ylim=c(-20,20),
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5);
TST <- NULL;
for(j in unique(dat[,5])){
    for(i in unique(dat[,4])){
        sub <- which(dat[,4]==i & dat[,5]==j);
        cis <- simpleboot(dat[sub,9]);
        ttt <- cis[1] * cis[2];
        ttt <- max(c(ttt,0));
        ttt[ttt > 0] <-1;
        TST <- rbind(TST,c(j,i,cis,ttt));
    }
}
for(i in 1:dim(TST)[1]){
    if(TST[i,5] == 1){
        lsh <- fullblock[i] - 0.25;
        rsh <- fullblock[i] + 0.25; 
    }
}
boxplots.dat <- boxplot(as.numeric(value)~bs+cs,data=sda,
                xlim=c(1,max(fullblock)),yaxt="n",xaxt="n",
                at=fullblock,lwd=1.5,cex.lab=1.75,ylim=c(-20,20),
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5,add=TRUE);
abline(h=0,lwd=0.8,lty="dotted");
text(x=0,y=18,cex=1.75,labels=expression(paste("(B) ",S['100,10'])),pos=4);
# ----------------------- Panel `c':
use <- which(em[,1]==100 & em[,2]==2);
dat <- em[use,];
nBn <- even.Beta(dat[,4]);
fourblock <- c(1,1.8,2.6,3.4,4.2);
fullblock <- c(fourblock,fourblock+6,fourblock+12,fourblock+18);

sdd <- as.data.frame(cbind(dat[,4],dat[,5],dat[,9]));
colnames(sdd) <- c("bs","cs","vl");
sda <- melt(sdd, id = c('bs','cs'));
sda <- sda[,-3];
boxplots.dat <- boxplot(as.numeric(value)~bs+cs,data=sda,type="n",
                xlim=c(1,max(fullblock)),yaxt="n",xaxt="n",
                at=fullblock,lwd=1.5,cex.lab=1.75,ylim=c(-20,20),
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5);
TST <- NULL;
for(j in unique(dat[,5])){
    for(i in unique(dat[,4])){
        sub <- which(dat[,4]==i & dat[,5]==j);
        cis <- simpleboot(dat[sub,9]);
        ttt <- cis[1] * cis[2];
        ttt <- max(c(ttt,0));
        ttt[ttt > 0] <-1;
        TST <- rbind(TST,c(j,i,cis,ttt));
    }
}
for(i in 1:dim(TST)[1]){
    if(TST[i,5] == 1){
        lsh <- fullblock[i] - 0.25;
        rsh <- fullblock[i] + 0.25; 
    }
}
boxplots.dat <- boxplot(as.numeric(value)~bs+cs,data=sda,
                xlim=c(1,max(fullblock)),yaxt="n",xaxt="n",
                at=fullblock,lwd=1.5,cex.lab=1.75,ylim=c(-20,20),
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5,add=TRUE);
abline(h=0,lwd=0.8,lty="dotted");
text(x=0,y=18,cex=1.75,labels=expression(paste("(C) ",S['100,2'])),pos=4);
# ----------------------- Panel `d':
use <- which(em[,1]==10 & em[,2]==100);
dat <- em[use,];
nBn <- even.Beta(dat[,4]);
fourblock <- c(1,1.8,2.6,3.4,4.2);
fullblock <- c(fourblock,fourblock+6,fourblock+12,fourblock+18);

sdd <- as.data.frame(cbind(dat[,4],dat[,5],dat[,9]));
colnames(sdd) <- c("bs","cs","vl");
sda <- melt(sdd, id = c('bs','cs'));
sda <- sda[,-3];
boxplots.dat <- boxplot(as.numeric(value)~bs+cs,data=sda,type="n",
                xlim=c(1,max(fullblock)),yaxt="n",xaxt="n",
                at=fullblock,lwd=1.5,cex.lab=1.75,ylim=c(-20,20),
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5);
TST <- NULL;
for(j in unique(dat[,5])){
    for(i in unique(dat[,4])){
        sub <- which(dat[,4]==i & dat[,5]==j);
        cis <- simpleboot(dat[sub,9]);
        ttt <- cis[1] * cis[2];
        ttt <- max(c(ttt,0));
        ttt[ttt > 0] <-1;
        TST <- rbind(TST,c(j,i,cis,ttt));
    }
}
for(i in 1:dim(TST)[1]){
    if(TST[i,5] == 1){
        lsh <- fullblock[i] - 0.25;
        rsh <- fullblock[i] + 0.25; 
    }
}
boxplots.dat <- boxplot(as.numeric(value)~bs+cs,data=sda,
                xlim=c(1,max(fullblock)),yaxt="n",xaxt="n",
                at=fullblock,lwd=1.5,cex.lab=1.75,ylim=c(-20,20),
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5,add=TRUE);
abline(h=0,lwd=0.8,lty="dotted");
text(x=0,y=18,cex=1.75,labels=expression(paste("(D) ",S['10,100'])),pos=4);
axis(side=2,at=c(-20,-10,0,10),cex.axis=1.5);
# ----------------------- Panel `e':
use <- which(em[,1]==10 & em[,2]==10);
dat <- em[use,];
nBn <- even.Beta(dat[,4]);
fourblock <- c(1,1.8,2.6,3.4,4.2);
fullblock <- c(fourblock,fourblock+6,fourblock+12,fourblock+18);

sdd <- as.data.frame(cbind(dat[,4],dat[,5],dat[,9]));
colnames(sdd) <- c("bs","cs","vl");
sda <- melt(sdd, id = c('bs','cs'));
sda <- sda[,-3];
boxplots.dat <- boxplot(as.numeric(value)~bs+cs,data=sda,type="n",
                xlim=c(1,max(fullblock)),yaxt="n",xaxt="n",
                at=fullblock,lwd=1.5,cex.lab=1.75,ylim=c(-20,20),
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5);
TST <- NULL;
for(j in unique(dat[,5])){
    for(i in unique(dat[,4])){
        sub <- which(dat[,4]==i & dat[,5]==j);
        cis <- simpleboot(dat[sub,9]);
        ttt <- cis[1] * cis[2];
        ttt <- max(c(ttt,0));
        ttt[ttt > 0] <-1;
        TST <- rbind(TST,c(j,i,cis,ttt));
    }
}
for(i in 1:dim(TST)[1]){
    if(TST[i,5] == 1){
        lsh <- fullblock[i] - 0.25;
        rsh <- fullblock[i] + 0.25; 
    }
}
boxplots.dat <- boxplot(as.numeric(value)~bs+cs,data=sda,
                xlim=c(1,max(fullblock)),yaxt="n",xaxt="n",
                at=fullblock,lwd=1.5,cex.lab=1.75,ylim=c(-20,20),
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5,add=TRUE);
abline(h=0,lwd=0.8,lty="dotted");
text(x=0,y=18,cex=1.75,labels=expression(paste("(E) ",S['10,10'])),pos=4);
# ----------------------- Panel `f':
use <- which(em[,1]==10 & em[,2]==2);
dat <- em[use,];
nBn <- even.Beta(dat[,4]);
fourblock <- c(1,1.8,2.6,3.4,4.2);
fullblock <- c(fourblock,fourblock+6,fourblock+12,fourblock+18);

sdd <- as.data.frame(cbind(dat[,4],dat[,5],dat[,9]));
colnames(sdd) <- c("bs","cs","vl");
sda <- melt(sdd, id = c('bs','cs'));
sda <- sda[,-3];
boxplots.dat <- boxplot(as.numeric(value)~bs+cs,data=sda,type="n",
                xlim=c(1,max(fullblock)),yaxt="n",xaxt="n",
                at=fullblock,lwd=1.5,cex.lab=1.75,ylim=c(-20,20),
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5);
TST <- NULL;
for(j in unique(dat[,5])){
    for(i in unique(dat[,4])){
        sub <- which(dat[,4]==i & dat[,5]==j);
        cis <- simpleboot(dat[sub,9]);
        ttt <- cis[1] * cis[2];
        ttt <- max(c(ttt,0));
        ttt[ttt > 0] <-1;
        TST <- rbind(TST,c(j,i,cis,ttt));
    }
}
for(i in 1:dim(TST)[1]){
    if(TST[i,5] == 1){
        lsh <- fullblock[i] - 0.25;
        rsh <- fullblock[i] + 0.25; 
    }
}
boxplots.dat <- boxplot(as.numeric(value)~bs+cs,data=sda,
                xlim=c(1,max(fullblock)),yaxt="n",xaxt="n",
                at=fullblock,lwd=1.5,cex.lab=1.75,ylim=c(-20,20),
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5,add=TRUE);
abline(h=0,lwd=0.8,lty="dotted");
text(x=0,y=18,cex=1.75,labels=expression(paste("(F) ",S['10,2'])),pos=4);
# ----------------------- Panel `g':
use <- which(em[,1]==2 & em[,2]==100);
dat <- em[use,];
nBn <- even.Beta(dat[,4]);
fourblock <- c(1,1.8,2.6,3.4,4.2);
fullblock <- c(fourblock,fourblock+6,fourblock+12,fourblock+18);

sdd <- as.data.frame(cbind(dat[,4],dat[,5],dat[,9]));
colnames(sdd) <- c("bs","cs","vl");
sda <- melt(sdd, id = c('bs','cs'));
sda <- sda[,-3];
boxplots.dat <- boxplot(as.numeric(value)~bs+cs,data=sda,type="n",
                xlim=c(1,max(fullblock)),yaxt="n",xaxt="n",
                at=fullblock,lwd=1.5,cex.lab=1.75,ylim=c(-20,20),
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5);
TST <- NULL;
for(j in unique(dat[,5])){
    for(i in unique(dat[,4])){
        sub <- which(dat[,4]==i & dat[,5]==j);
        cis <- simpleboot(dat[sub,9]);
        ttt <- cis[1] * cis[2];
        ttt <- max(c(ttt,0));
        ttt[ttt > 0] <-1;
        TST <- rbind(TST,c(j,i,cis,ttt));
    }
}
for(i in 1:dim(TST)[1]){
    if(TST[i,5] == 1){
        lsh <- fullblock[i] - 0.25;
        rsh <- fullblock[i] + 0.25; 
    }
}
boxplots.dat <- boxplot(as.numeric(value)~bs+cs,data=sda,
                xlim=c(1,max(fullblock)),yaxt="n",xaxt="n",
                at=fullblock,lwd=1.5,cex.lab=1.75,ylim=c(-20,20),
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5,add=TRUE);
abline(h=0,lwd=0.8,lty="dotted");
text(x=0,y=18,cex=1.75,labels=expression(paste("(G) ",S['2,100'])),pos=4);
axis(side=2,at=c(-20,-10,0,10),cex.axis=1.5);
axis(side=1,at=c(fourblock[3],fourblock[3]+6,fourblock[3]+12,fourblock[3]+18),
     cex.axis=1.5,cex.lab=1.5,labels=c("0.00","0.01","0.02","0.03"));
# ----------------------- Panel `h':
use <- which(em[,1]==2 & em[,2]==10);
dat <- em[use,];
nBn <- even.Beta(dat[,4]);
fourblock <- c(1,1.8,2.6,3.4,4.2);
fullblock <- c(fourblock,fourblock+6,fourblock+12,fourblock+18);

sdd <- as.data.frame(cbind(dat[,4],dat[,5],dat[,9]));
colnames(sdd) <- c("bs","cs","vl");
sda <- melt(sdd, id = c('bs','cs'));
sda <- sda[,-3];
boxplots.dat <- boxplot(as.numeric(value)~bs+cs,data=sda,type="n",
                xlim=c(1,max(fullblock)),yaxt="n",xaxt="n",
                at=fullblock,lwd=1.5,cex.lab=1.75,ylim=c(-20,20),
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5);
TST <- NULL;
for(j in unique(dat[,5])){
    for(i in unique(dat[,4])){
        sub <- which(dat[,4]==i & dat[,5]==j);
        cis <- simpleboot(dat[sub,9]);
        ttt <- cis[1] * cis[2];
        ttt <- max(c(ttt,0));
        ttt[ttt > 0] <-1;
        TST <- rbind(TST,c(j,i,cis,ttt));
    }
}
for(i in 1:dim(TST)[1]){
    if(TST[i,5] == 1){
        lsh <- fullblock[i] - 0.25;
        rsh <- fullblock[i] + 0.25; 
    }
}
boxplots.dat <- boxplot(as.numeric(value)~bs+cs,data=sda,
                xlim=c(1,max(fullblock)),yaxt="n",xaxt="n",
                at=fullblock,lwd=1.5,cex.lab=1.75,ylim=c(-20,20),
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5,add=TRUE);
abline(h=0,lwd=0.8,lty="dotted");
text(x=0,y=18,cex=1.75,labels=expression(paste("(H) ",S['2,10'])),pos=4);
axis(side=1,at=c(fourblock[3],fourblock[3]+6,fourblock[3]+12,fourblock[3]+18),
     cex.axis=1.5,cex.lab=1.5,labels=c("0.00","0.01","0.02","0.03"));
# ----------------------- Panel `i':
use <- which(em[,1]==2 & em[,2]==2);
dat <- em[use,];
nBn <- even.Beta(dat[,4]);
fourblock <- c(1,1.8,2.6,3.4,4.2);
fullblock <- c(fourblock,fourblock+6,fourblock+12,fourblock+18);

sdd <- as.data.frame(cbind(dat[,4],dat[,5],dat[,9]));
colnames(sdd) <- c("bs","cs","vl");
sda <- melt(sdd, id = c('bs','cs'));
sda <- sda[,-3];
boxplots.dat <- boxplot(as.numeric(value)~bs+cs,data=sda,type="n",
                xlim=c(1,max(fullblock)),yaxt="n",xaxt="n",
                at=fullblock,lwd=1.5,cex.lab=1.75,ylim=c(-20,20),
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5);
TST <- NULL;
for(j in unique(dat[,5])){
    for(i in unique(dat[,4])){
        sub <- which(dat[,4]==i & dat[,5]==j);
        cis <- simpleboot(dat[sub,9]);
        ttt <- cis[1] * cis[2];
        ttt <- max(c(ttt,0));
        ttt[ttt > 0] <-1;
        TST <- rbind(TST,c(j,i,cis,ttt));
    }
}
for(i in 1:dim(TST)[1]){
    if(TST[i,5] == 1){
        lsh <- fullblock[i] - 0.25;
        rsh <- fullblock[i] + 0.25; 
    }
}
boxplots.dat <- boxplot(as.numeric(value)~bs+cs,data=sda,
                xlim=c(1,max(fullblock)),yaxt="n",xaxt="n",
                at=fullblock,lwd=1.5,cex.lab=1.75,ylim=c(-20,20),
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5,add=TRUE);
abline(h=0,lwd=0.8,lty="dotted");
text(x=0,y=18,cex=1.75,labels=expression(paste("(I) ",S['2,2'])),pos=4);
axis(side=1,at=c(fourblock[3],fourblock[3]+6,fourblock[3]+12,fourblock[3]+18),
     cex.axis=1.5,cex.lab=1.5,labels=c("0.00","0.01","0.02","0.03"));

mtext(expression(paste("Cost of polyandry (",c[P],")")),
	outer=TRUE,side=1,line=3.75,cex=1.5);

mtext(expression(paste("Neutral mean allele value (",eta[a],")")),
	outer=TRUE,side=2,line=3.25,cex=1.5);

dev.off();

# =====================================================================================

# =====================================================================================
# TODO: BUILD FIGURE S1-10                                                            #
# ====================================================================================#

setEPS(); # postscript below for final publication
cairo_ps("Nvalues_Nmin1.eps",family="Arial",height=7,width=7);
par(mar=c(5,5,2,2));
use <- which(em[,1]==-1 & em[,2]==100);
dat <- em[use,];
nBn <- even.Beta(dat[,4]);
fourblock <- c(1,1.8,2.6,3.4,4.2);
fullblock <- c(fourblock,fourblock+6,fourblock+12,fourblock+18);

sdd <- as.data.frame(cbind(dat[,4],dat[,5],dat[,9]));
colnames(sdd) <- c("bs","cs","vl");
sda <- melt(sdd, id = c('bs','cs'));
sda <- sda[,-3];
boxplots.dat <- boxplot(as.numeric(value)~bs+cs,data=sda,type="n",
                xlim=c(1,max(fullblock)),yaxt="n",xaxt="n",
                at=fullblock,lwd=1.5,cex.lab=1.75,ylim=c(-20,20),
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5,
                xlab=expression(paste("Cost of polyandry (",c[P],")")),
                ylab=expression(paste("Neutral mean allele value (",eta[a],")")));
TST <- NULL;
for(j in unique(dat[,5])){
    for(i in unique(dat[,4])){
        sub <- which(dat[,4]==i & dat[,5]==j);
        cis <- simpleboot(dat[sub,9]);
        ttt <- cis[1] * cis[2];
        ttt <- max(c(ttt,0));
        ttt[ttt > 0] <-1;
        TST <- rbind(TST,c(j,i,cis,ttt));
    }
}
for(i in 1:dim(TST)[1]){
    if(TST[i,5] == 1){
        lsh <- fullblock[i] - 0.25;
        rsh <- fullblock[i] + 0.25; 
    }
}
boxplots.dat <- boxplot(as.numeric(value)~bs+cs,data=sda,
                xlim=c(1,max(fullblock)),yaxt="n",xaxt="n",
                at=fullblock,lwd=1.5,cex.lab=1.75,ylim=c(-29,14),
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5,add=TRUE);
abline(h=0,lwd=0.8,lty="dotted");
axis(side=2,at=c(-20,-10,0,10),cex.axis=1.5);
axis(side=1,at=c(fourblock[3],fourblock[3]+6,fourblock[3]+12,fourblock[3]+18),
     cex.axis=1.5,cex.lab=1.5,labels=c("0.00","0.01","0.02","0.03"));

dev.off();

# =====================================================================================

# =====================================================================================
# TODO: BUILD FIGURE S1-11                                                            #
# ====================================================================================#

setEPS(); # postscript below for final publication
cairo_ps("Ndistr.eps",family="Arial",height=7,width=7);
par(mar=c(5,5,1,1),lwd=2);
hist(em[,9],breaks=100,xlim=c(-20,20),col="grey40",main="",freq=FALSE,
     xlab=expression(paste("Mean neutral allele value (",eta[a],")")),
     ylab="Probability density",cex.lab=1.5,cex.axis=1.5,lwd=2);

dev.off();

# ====================================================================================#


# =====================================================================================
# TODO: BUILD FIGURE S1-12                                                            #
# ====================================================================================#

# -------------- Build the figure. ============================================ XXX #
setEPS(); # postscript below for final publication?
cairo_ps("adjhist.eps",family="Arial",height=7,width=7);
par(mfrow=c(3,3),mar=c(0.2,0.2,0.2,0.2),oma=c(6,6,1,1));
# ----------------------- Panel `a':
bp1 <- f.bypr(pol=pol,initial=100,additional=100);
bp2 <- f.CIva(bypr=bp1);
bp2 <- bporder(bp2);
bp3 <- f.pomn(pol=pol,initial=100,additional=100);
bp3 <- bporder(bp3);
dat <- as.data.frame(cbind(bp1[,4],bp1[,2:3]));
colnames(dat) <- c("Corr","Beta","cost");
sda <- melt(dat, id = c('Beta','cost'));
sda <- sda[,-3];
fourblock <- c(1,1.8,2.6,3.4,4.2);
fullblock <- c(fourblock,fourblock+6,fourblock+12,fourblock+18);
boxplots.dat <- boxplot(as.numeric(value)~Beta+cost,data=sda,xlim=c(0,fullblock[20]+2),
                at=fullblock,lwd=2.25,xaxt="n",yaxt="n",
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5,ylim=c(-0.05,0.05));
axis(side=2,at=c(-0.04,-0.02,0,0.02,0.04),cex.axis=1.7,cex.lab=2);
boxplots.dat <- boxplot(as.numeric(value)~Beta+cost,data=sda,xlim=c(0,fullblock[20]+2),
                at=fullblock,lwd=2.25,add=TRUE,xaxt="n",yaxt="n",
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5,ylim=c(-0.05,0.05));
abline(h=0,lty="dotted",lwd=0.8);
text(x=-0.0025,y=0.0425,cex=1.75,labels=expression(paste("(A) ",S['100,100'])),pos=4);
# ----------------------- Panel `b':
bp1 <- f.bypr(pol=pol,initial=100,additional=10);
bp2 <- f.CIva(bypr=bp1);
bp2 <- bporder(bp2);
bp3 <- f.pomn(pol=pol,initial=100,additional=10);
bp3 <- bporder(bp3);
dat <- as.data.frame(cbind(bp1[,4],bp1[,2:3]));
colnames(dat) <- c("Corr","Beta","cost");
sda <- melt(dat, id = c('Beta','cost'));
sda <- sda[,-3];
fourblock <- c(1,1.8,2.6,3.4,4.2);
fullblock <- c(fourblock,fourblock+6,fourblock+12,fourblock+18);
boxplots.dat <- boxplot(as.numeric(value)~Beta+cost,data=sda,xlim=c(0,fullblock[20]+2),
                at=fullblock,lwd=2.25,xaxt="n",yaxt="n",
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5,ylim=c(-0.05,0.05));
boxplots.dat <- boxplot(as.numeric(value)~Beta+cost,data=sda,xlim=c(0,fullblock[20]+2),
                at=fullblock,lwd=2.25,add=TRUE,xaxt="n",yaxt="n",
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5,ylim=c(-0.05,0.05));
abline(h=0,lty="dotted",lwd=0.8);
text(x=-0.0025,y=0.0425,cex=1.75,labels=expression(paste("(B) ",S['100,10'])),pos=4);
# ----------------------- Panel `c':
bp1 <- f.bypr(pol=pol,initial=100,additional=2);
bp2 <- f.CIva(bypr=bp1);
bp2 <- bporder(bp2);
bp3 <- f.pomn(pol=pol,initial=100,additional=2);
bp3 <- bporder(bp3);
dat <- as.data.frame(cbind(bp1[,4],bp1[,2:3]));
colnames(dat) <- c("Corr","Beta","cost");
sda <- melt(dat, id = c('Beta','cost'));
sda <- sda[,-3];
fourblock <- c(1,1.8,2.6,3.4,4.2);
fullblock <- c(fourblock,fourblock+6,fourblock+12,fourblock+18);
boxplots.dat <- boxplot(as.numeric(value)~Beta+cost,data=sda,xlim=c(0,fullblock[20]+2),
                at=fullblock,lwd=2.25,xaxt="n",yaxt="n",
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5,ylim=c(-0.05,0.05));
boxplots.dat <- boxplot(as.numeric(value)~Beta+cost,data=sda,xlim=c(0,fullblock[20]+2),
                at=fullblock,lwd=2.25,add=TRUE,xaxt="n",yaxt="n",
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5,ylim=c(-0.05,0.05));
abline(h=0,lty="dotted",lwd=0.8);
text(x=-0.0025,y=0.0425,cex=1.75,labels=expression(paste("(C) ",S['100,2'])),pos=4);
# ----------------------- Panel `d':
bp1 <- f.bypr(pol=pol,initial=10,additional=100);
bp2 <- f.CIva(bypr=bp1);
bp2 <- bporder(bp2);
bp3 <- f.pomn(pol=pol,initial=10,additional=100);
bp3 <- bporder(bp3);
dat <- as.data.frame(cbind(bp1[,4],bp1[,2:3]));
colnames(dat) <- c("Corr","Beta","cost");
sda <- melt(dat, id = c('Beta','cost'));
sda <- sda[,-3];
fourblock <- c(1,1.8,2.6,3.4,4.2);
fullblock <- c(fourblock,fourblock+6,fourblock+12,fourblock+18);
boxplots.dat <- boxplot(as.numeric(value)~Beta+cost,data=sda,xlim=c(0,fullblock[20]+2),
                at=fullblock,lwd=2.25,xaxt="n",yaxt="n",
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5,ylim=c(-0.05,0.05));
axis(side=2,at=c(-0.04,-0.02,0,0.02,0.04),cex.axis=1.7,cex.lab=2);
boxplots.dat <- boxplot(as.numeric(value)~Beta+cost,data=sda,xlim=c(0,fullblock[20]+2),
                at=fullblock,lwd=2.25,add=TRUE,xaxt="n",yaxt="n",
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5,ylim=c(-0.05,0.05));
abline(h=0,lty="dotted",lwd=0.8);
text(x=-0.0025,y=0.0425,cex=1.75,labels=expression(paste("(D) ",S['10,100'])),pos=4);
# ----------------------- Panel `e':
bp1 <- f.bypr(pol=pol,initial=10,additional=10);
bp2 <- f.CIva(bypr=bp1);
bp2 <- bporder(bp2);
bp3 <- f.pomn(pol=pol,initial=10,additional=10);
bp3 <- bporder(bp3);
dat <- as.data.frame(cbind(bp1[,4],bp1[,2:3]));
colnames(dat) <- c("Corr","Beta","cost");
sda <- melt(dat, id = c('Beta','cost'));
sda <- sda[,-3];
fourblock <- c(1,1.8,2.6,3.4,4.2);
fullblock <- c(fourblock,fourblock+6,fourblock+12,fourblock+18);
boxplots.dat <- boxplot(as.numeric(value)~Beta+cost,data=sda,xlim=c(0,fullblock[20]+2),
                at=fullblock,lwd=2.25,xaxt="n",yaxt="n",
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5,ylim=c(-0.05,0.05));
boxplots.dat <- boxplot(as.numeric(value)~Beta+cost,data=sda,xlim=c(0,fullblock[20]+2),
                at=fullblock,lwd=2.25,add=TRUE,xaxt="n",yaxt="n",
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5,ylim=c(-0.05,0.05));
abline(h=0,lty="dotted",lwd=0.8);
text(x=-0.0025,y=0.0425,cex=1.75,labels=expression(paste("(E) ",S['10,10'])),pos=4);
# ----------------------- Panel `f':
bp1 <- f.bypr(pol=pol,initial=10,additional=2);
bp2 <- f.CIva(bypr=bp1);
bp2 <- bporder(bp2);
bp3 <- f.pomn(pol=pol,initial=10,additional=2);
bp3 <- bporder(bp3);
dat <- as.data.frame(cbind(bp1[,4],bp1[,2:3]));
colnames(dat) <- c("Corr","Beta","cost");
sda <- melt(dat, id = c('Beta','cost'));
sda <- sda[,-3];
fourblock <- c(1,1.8,2.6,3.4,4.2);
fullblock <- c(fourblock,fourblock+6,fourblock+12,fourblock+18);
boxplots.dat <- boxplot(as.numeric(value)~Beta+cost,data=sda,xlim=c(0,fullblock[20]+2),
                at=fullblock,lwd=2.25,xaxt="n",yaxt="n",
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5,ylim=c(-0.05,0.05));
boxplots.dat <- boxplot(as.numeric(value)~Beta+cost,data=sda,xlim=c(0,fullblock[20]+2),
                at=fullblock,lwd=2.25,add=TRUE,xaxt="n",yaxt="n",
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5,ylim=c(-0.05,0.05));
abline(h=0,lty="dotted",lwd=0.8);
text(x=-0.0025,y=0.0425,cex=1.75,labels=expression(paste("(F) ",S['10,2'])),pos=4);
# ----------------------- Panel `g':
bp1 <- f.bypr(pol=pol,initial=2,additional=100);
bp2 <- f.CIva(bypr=bp1);
bp2 <- bporder(bp2);
bp3 <- f.pomn(pol=pol,initial=2,additional=100);
bp3 <- bporder(bp3);
dat <- as.data.frame(cbind(bp1[,4],bp1[,2:3]));
colnames(dat) <- c("Corr","Beta","cost");
sda <- melt(dat, id = c('Beta','cost'));
sda <- sda[,-3];
fourblock <- c(1,1.8,2.6,3.4,4.2);
fullblock <- c(fourblock,fourblock+6,fourblock+12,fourblock+18);
boxplots.dat <- boxplot(as.numeric(value)~Beta+cost,data=sda,xlim=c(0,fullblock[20]+2),
                at=fullblock,lwd=2.25,xaxt="n",yaxt="n",
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5,ylim=c(-0.05,0.05));
axis(side=2,at=c(-0.04,-0.02,0,0.02,0.04),cex.axis=1.7,cex.lab=2);
boxplots.dat <- boxplot(as.numeric(value)~Beta+cost,data=sda,xlim=c(0,fullblock[20]+2),
                at=fullblock,lwd=2.25,add=TRUE,xaxt="n",yaxt="n",
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5,ylim=c(-0.05,0.05));
abline(h=0,lty="dotted",lwd=0.8);
axis(side=1,at=c(fourblock[3],fourblock[3]+6,fourblock[3]+12,fourblock[3]+18),
     cex.axis=2,cex.lab=2,labels=c("0.00","0.01","0.02","0.03"));
text(x=-0.0025,y=0.0425,cex=1.75,labels=expression(paste("(G) ",S['2,100'])),pos=4);
# ----------------------- Panel `h':
bp1 <- f.bypr(pol=pol,initial=2,additional=10);
bp2 <- f.CIva(bypr=bp1);
bp2 <- bporder(bp2);
bp3 <- f.pomn(pol=pol,initial=2,additional=10);
bp3 <- bporder(bp3);
dat <- as.data.frame(cbind(bp1[,4],bp1[,2:3]));
colnames(dat) <- c("Corr","Beta","cost");
sda <- melt(dat, id = c('Beta','cost'));
sda <- sda[,-3];
fourblock <- c(1,1.8,2.6,3.4,4.2);
fullblock <- c(fourblock,fourblock+6,fourblock+12,fourblock+18);
boxplots.dat <- boxplot(as.numeric(value)~Beta+cost,data=sda,xlim=c(0,fullblock[20]+2),
                at=fullblock,lwd=2.25,xaxt="n",yaxt="n",
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5,ylim=c(-0.05,0.05));
boxplots.dat <- boxplot(as.numeric(value)~Beta+cost,data=sda,xlim=c(0,fullblock[20]+2),
                at=fullblock,lwd=2.25,add=TRUE,xaxt="n",yaxt="n",
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5,ylim=c(-0.05,0.05));
abline(h=0,lty="dotted",lwd=0.8);
axis(side=1,at=c(fourblock[3],fourblock[3]+6,fourblock[3]+12,fourblock[3]+18),
     cex.axis=2,cex.lab=2,labels=c("0.00","0.01","0.02","0.03"));
text(x=-0.0025,y=0.0425,cex=1.75,labels=expression(paste("(H) ",S['2,10'])),pos=4);
# ----------------------- Panel `i':
bp1 <- f.bypr(pol=pol,initial=2,additional=2);
bp2 <- f.CIva(bypr=bp1);
bp2 <- bporder(bp2);
bp3 <- f.pomn(pol=pol,initial=2,additional=2);
bp3 <- bporder(bp3);
dat <- as.data.frame(cbind(bp1[,4],bp1[,2:3]));
colnames(dat) <- c("Corr","Beta","cost");
sda <- melt(dat, id = c('Beta','cost'));
sda <- sda[,-3];
fourblock <- c(1,1.8,2.6,3.4,4.2);
fullblock <- c(fourblock,fourblock+6,fourblock+12,fourblock+18);
boxplots.dat <- boxplot(as.numeric(value)~Beta+cost,data=sda,xlim=c(0,fullblock[20]+2),
                at=fullblock,lwd=2.25,xaxt="n",yaxt="n",
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5,ylim=c(-0.05,0.05));
boxplots.dat <- boxplot(as.numeric(value)~Beta+cost,data=sda,xlim=c(0,fullblock[20]+2),
                at=fullblock,lwd=2.25,add=TRUE,xaxt="n",yaxt="n",
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5,ylim=c(-0.05,0.05));
abline(h=0,lty="dotted",lwd=0.8);
axis(side=1,at=c(fourblock[3],fourblock[3]+6,fourblock[3]+12,fourblock[3]+18),
     cex.axis=2,cex.lab=2,labels=c("0.00","0.01","0.02","0.03"));
text(x=-0.0025,y=0.0425,cex=1.75,labels=expression(paste("(I) ",S['2,2'])),pos=4);
# --------------------------------------------
mtext(expression(paste("Cost of polyandry (",c[P],")")),
	outer=TRUE,side=1,line=3.75,cex=1.5);

mtext(expression(paste("Mean inbreeding adjustment (",k[adj],")")),
	outer=TRUE,side=2,line=3.25,cex=1.5);

dev.off();
# =====================================================================================


# =====================================================================================
# TODO: BUILD FIGURE S1-13                                                            #
# ====================================================================================#
setEPS(); # postscript below for final publication?
cairo_ps("adjhist_Nmin1.eps",family="Arial",height=7,width=7);
par(mar=c(0.2,0.2,0.2,0.2),oma=c(6,6,1,1),lwd=2);
bp1 <- f.bypr(pol=pol,initial=-1,additional=100);
bp2 <- f.CIva(bypr=bp1);
bp2 <- bporder(bp2);
bp3 <- f.pomn(pol=pol,initial=-1,additional=100);
bp3 <- bporder(bp3);
dat <- as.data.frame(cbind(bp1[,4],bp1[,2:3]));
colnames(dat) <- c("Corr","Beta","cost");
sda <- melt(dat, id = c('Beta','cost'));
sda <- sda[,-3];
fourblock <- c(1,1.8,2.6,3.4,4.2);
fullblock <- c(fourblock,fourblock+6,fourblock+12,fourblock+18);
boxplots.dat <- boxplot(as.numeric(value)~Beta+cost,data=sda,xlim=c(0,fullblock[20]+2),
                at=fullblock,lwd=2.25,xaxt="n",yaxt="n",
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5,ylim=c(-0.05,0.05));
axis(side=2,at=c(-0.04,-0.03,-0.02,-0.01,0,0.01,0.02,0.03,0.04),cex.axis=1.7,cex.lab=2);
boxplots.dat <- boxplot(as.numeric(value)~Beta+cost,data=sda,xlim=c(0,fullblock[20]+2),
                at=fullblock,lwd=2.25,add=TRUE,xaxt="n",yaxt="n",
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5,ylim=c(-0.05,0.05));
abline(h=0,lty="dotted",lwd=0.8);
axis(side=1,at=c(fourblock[3],fourblock[3]+6,fourblock[3]+12,fourblock[3]+18),
     cex.axis=2,cex.lab=2,labels=c("0.00","0.01","0.02","0.03"));
mtext(expression(paste("Cost of polyandry (",c[P],")")),
	outer=TRUE,side=1,line=3.75,cex=2);
mtext(expression(paste("Mean inbreeding adjustment (",k[adj],") in population")),
	outer=TRUE,side=2,line=3.25,cex=1.5);
dev.off()
# ====================================================================================#


# =====================================================================================
# TODO: BUILD FIGURE S1-14                                                            #
# ====================================================================================#
setEPS(); # postscript below for final publication?
cairo_ps("trait_cors.eps",family="Arial",height=7,width=7);
par(mfrow=c(3,3),mar=c(0.2,0.2,0.2,0.2),oma=c(6,6,1,1));
# ----------------------- Panel `a':
use <- which(last.gen[,1] == 100 & last.gen[,2] == 100);
crr <- cbind(last.gen[use,1:2],last.gen[use,4:5],last.gen[use,14]);
dat <- as.data.frame(cbind(crr[,5],crr[,3:4]));
colnames(dat) <- c("Corr","bs","cs");
sda <- melt(dat, id = c('bs','cs'));
sda <- sda[,-3];
fourblock <- c(1,1.8,2.6,3.4,4.2);
fullblock <- c(fourblock,fourblock+6,fourblock+12,fourblock+18);
boxplots.dat <- boxplot(as.numeric(value)~bs+cs,data=sda,type="n",
                xlim=c(1,max(fullblock)),yaxt="n",xaxt="n",
                at=fullblock,lwd=1.5,cex.lab=1.75,ylim=c(-0.3,0.3),
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5);
TST <- NULL;
for(j in unique(dat[,3])){
    for(i in unique(dat[,2])){
        sub <- which(dat[,2]==i & dat[,3]==j);
        cis <- simpleboot(dat[sub,1]);
        ttt <- cis[1] * cis[2];
        ttt <- max(c(ttt,0));
        ttt[ttt > 0] <-1;
        TST <- rbind(TST,c(j,i,cis,ttt));
    }
}
TST <- bcsort(TST);
for(i in 1:dim(TST)[1]){
    if(TST[i,5] == 1){
        lsh <- fullblock[i] - 0.25;
        rsh <- fullblock[i] + 0.25; 
        #polygon(x=c(lsh,rsh,rsh,lsh),y=c(-40,-40,15,15),border=NA,col="grey70");
    }
}
boxplots.dat <- boxplot(as.numeric(value)~bs+cs,data=sda,
                xlim=c(1,max(fullblock)),yaxt="n",xaxt="n",
                at=fullblock,lwd=1.5,cex.lab=1.75,ylim=c(-0.3,0.3),
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5,add=TRUE);
abline(h=0,lwd=0.8,lty="dotted");
axis(side=2,at=c(-0.2,0,0.2),cex.axis=1.5);
text(x=0,y=0.26,cex=1.75,labels=expression(paste("(A) ",S['100,100'])),pos=4);
# ----------------------- Panel `b':
use <- which(last.gen[,1] == 100 & last.gen[,2] == 10);
crr <- cbind(last.gen[use,1:2],last.gen[use,4:5],last.gen[use,14]);
dat <- as.data.frame(cbind(crr[,5],crr[,3:4]));
colnames(dat) <- c("Corr","bs","cs");
sda <- melt(dat, id = c('bs','cs'));
sda <- sda[,-3];
fourblock <- c(1,1.8,2.6,3.4,4.2);
fullblock <- c(fourblock,fourblock+6,fourblock+12,fourblock+18);
boxplots.dat <- boxplot(as.numeric(value)~bs+cs,data=sda,type="n",
                xlim=c(1,max(fullblock)),yaxt="n",xaxt="n",
                at=fullblock,lwd=1.5,cex.lab=1.75,ylim=c(-0.3,0.3),
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5);
TST <- NULL;
for(j in unique(dat[,3])){
    for(i in unique(dat[,2])){
        sub <- which(dat[,2]==i & dat[,3]==j);
        cis <- simpleboot(dat[sub,1]);
        ttt <- cis[1] * cis[2];
        ttt <- max(c(ttt,0));
        ttt[ttt > 0] <-1;
        TST <- rbind(TST,c(j,i,cis,ttt));
    }
}
TST <- bcsort(TST);
for(i in 1:dim(TST)[1]){
    if(TST[i,5] == 1){
        lsh <- fullblock[i] - 0.25;
        rsh <- fullblock[i] + 0.25; 
        #polygon(x=c(lsh,rsh,rsh,lsh),y=c(-40,-40,15,15),border=NA,col="grey70");
    }
}
boxplots.dat <- boxplot(as.numeric(value)~bs+cs,data=sda,
                xlim=c(1,max(fullblock)),yaxt="n",xaxt="n",
                at=fullblock,lwd=1.5,cex.lab=1.75,ylim=c(-0.3,0.3),
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5,add=TRUE);
abline(h=0,lwd=0.8,lty="dotted");
text(x=0,y=0.26,cex=1.75,labels=expression(paste("(B) ",S['100,10'])),pos=4);
# ----------------------- Panel `c':
use <- which(last.gen[,1] == 100 & last.gen[,2] == 2);
crr <- cbind(last.gen[use,1:2],last.gen[use,4:5],last.gen[use,14]);
dat <- as.data.frame(cbind(crr[,5],crr[,3:4]));
colnames(dat) <- c("Corr","bs","cs");
sda <- melt(dat, id = c('bs','cs'));
sda <- sda[,-3];
fourblock <- c(1,1.8,2.6,3.4,4.2);
fullblock <- c(fourblock,fourblock+6,fourblock+12,fourblock+18);
boxplots.dat <- boxplot(as.numeric(value)~bs+cs,data=sda,type="n",
                xlim=c(1,max(fullblock)),yaxt="n",xaxt="n",
                at=fullblock,lwd=1.5,cex.lab=1.75,ylim=c(-0.3,0.3),
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5);
TST <- NULL;
for(j in unique(dat[,3])){
    for(i in unique(dat[,2])){
        sub <- which(dat[,2]==i & dat[,3]==j);
        cis <- simpleboot(dat[sub,1]);
        ttt <- cis[1] * cis[2];
        ttt <- max(c(ttt,0));
        ttt[ttt > 0] <-1;
        TST <- rbind(TST,c(j,i,cis,ttt));
    }
}
TST <- bcsort(TST);
for(i in 1:dim(TST)[1]){
    if(TST[i,5] == 1){
        lsh <- fullblock[i] - 0.25;
        rsh <- fullblock[i] + 0.25; 
        #polygon(x=c(lsh,rsh,rsh,lsh),y=c(-40,-40,15,15),border=NA,col="grey70");
    }
}
boxplots.dat <- boxplot(as.numeric(value)~bs+cs,data=sda,
                xlim=c(1,max(fullblock)),yaxt="n",xaxt="n",
                at=fullblock,lwd=1.5,cex.lab=1.75,ylim=c(-0.3,0.3),
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5,add=TRUE);
abline(h=0,lwd=0.8,lty="dotted");
text(x=0,y=0.26,cex=1.75,labels=expression(paste("(C) ",S['100,2'])),pos=4);
# ----------------------- Panel `d':
use <- which(last.gen[,1] == 10 & last.gen[,2] == 100);
crr <- cbind(last.gen[use,1:2],last.gen[use,4:5],last.gen[use,14]);
dat <- as.data.frame(cbind(crr[,5],crr[,3:4]));
colnames(dat) <- c("Corr","bs","cs");
sda <- melt(dat, id = c('bs','cs'));
sda <- sda[,-3];
fourblock <- c(1,1.8,2.6,3.4,4.2);
fullblock <- c(fourblock,fourblock+6,fourblock+12,fourblock+18);
boxplots.dat <- boxplot(as.numeric(value)~bs+cs,data=sda,type="n",
                xlim=c(1,max(fullblock)),yaxt="n",xaxt="n",
                at=fullblock,lwd=1.5,cex.lab=1.75,ylim=c(-0.3,0.3),
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5);
TST <- NULL;
for(j in unique(dat[,3])){
    for(i in unique(dat[,2])){
        sub <- which(dat[,2]==i & dat[,3]==j);
        cis <- simpleboot(dat[sub,1]);
        ttt <- cis[1] * cis[2];
        ttt <- max(c(ttt,0));
        ttt[ttt > 0] <-1;
        TST <- rbind(TST,c(j,i,cis,ttt));
    }
}
TST <- bcsort(TST);
for(i in 1:dim(TST)[1]){
    if(TST[i,5] == 1){
        lsh <- fullblock[i] - 0.25;
        rsh <- fullblock[i] + 0.25; 
        #polygon(x=c(lsh,rsh,rsh,lsh),y=c(-40,-40,15,15),border=NA,col="grey70");
    }
}
boxplots.dat <- boxplot(as.numeric(value)~bs+cs,data=sda,
                xlim=c(1,max(fullblock)),yaxt="n",xaxt="n",
                at=fullblock,lwd=1.5,cex.lab=1.75,ylim=c(-0.3,0.3),
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5,add=TRUE);
abline(h=0,lwd=0.8,lty="dotted");
axis(side=2,at=c(-0.2,0,0.2),cex.axis=1.5);
text(x=0,y=0.26,cex=1.75,labels=expression(paste("(D) ",S['10,100'])),pos=4);
# ----------------------- Panel `e':
use <- which(last.gen[,1] == 10 & last.gen[,2] == 10);
crr <- cbind(last.gen[use,1:2],last.gen[use,4:5],last.gen[use,14]);
dat <- as.data.frame(cbind(crr[,5],crr[,3:4]));
colnames(dat) <- c("Corr","bs","cs");
sda <- melt(dat, id = c('bs','cs'));
sda <- sda[,-3];
fourblock <- c(1,1.8,2.6,3.4,4.2);
fullblock <- c(fourblock,fourblock+6,fourblock+12,fourblock+18);
boxplots.dat <- boxplot(as.numeric(value)~bs+cs,data=sda,type="n",
                xlim=c(1,max(fullblock)),yaxt="n",xaxt="n",
                at=fullblock,lwd=1.5,cex.lab=1.75,ylim=c(-0.3,0.3),
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5);
TST <- NULL;
for(j in unique(dat[,3])){
    for(i in unique(dat[,2])){
        sub <- which(dat[,2]==i & dat[,3]==j);
        cis <- simpleboot(dat[sub,1]);
        ttt <- cis[1] * cis[2];
        ttt <- max(c(ttt,0));
        ttt[ttt > 0] <-1;
        TST <- rbind(TST,c(j,i,cis,ttt));
    }
}
TST <- bcsort(TST);
for(i in 1:dim(TST)[1]){
    if(TST[i,5] == 1){
        lsh <- fullblock[i] - 0.25;
        rsh <- fullblock[i] + 0.25; 
        #polygon(x=c(lsh,rsh,rsh,lsh),y=c(-40,-40,15,15),border=NA,col="grey70");
    }
}
boxplots.dat <- boxplot(as.numeric(value)~bs+cs,data=sda,
                xlim=c(1,max(fullblock)),yaxt="n",xaxt="n",
                at=fullblock,lwd=1.5,cex.lab=1.75,ylim=c(-0.3,0.3),
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5,add=TRUE);
abline(h=0,lwd=0.8,lty="dotted");
text(x=0,y=0.26,cex=1.75,labels=expression(paste("(E) ",S['10,10'])),pos=4);
# ----------------------- Panel `f':
use <- which(last.gen[,1] == 10 & last.gen[,2] == 2);
crr <- cbind(last.gen[use,1:2],last.gen[use,4:5],last.gen[use,14]);
dat <- as.data.frame(cbind(crr[,5],crr[,3:4]));
colnames(dat) <- c("Corr","bs","cs");
sda <- melt(dat, id = c('bs','cs'));
sda <- sda[,-3];
fourblock <- c(1,1.8,2.6,3.4,4.2);
fullblock <- c(fourblock,fourblock+6,fourblock+12,fourblock+18);
boxplots.dat <- boxplot(as.numeric(value)~bs+cs,data=sda,type="n",
                xlim=c(1,max(fullblock)),yaxt="n",xaxt="n",
                at=fullblock,lwd=1.5,cex.lab=1.75,ylim=c(-0.3,0.3),
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5);
TST <- NULL;
for(j in unique(dat[,3])){
    for(i in unique(dat[,2])){
        sub <- which(dat[,2]==i & dat[,3]==j);
        cis <- simpleboot(dat[sub,1]);
        ttt <- cis[1] * cis[2];
        ttt <- max(c(ttt,0));
        ttt[ttt > 0] <-1;
        TST <- rbind(TST,c(j,i,cis,ttt));
    }
}
TST <- bcsort(TST);
for(i in 1:dim(TST)[1]){
    if(TST[i,5] == 1){
        lsh <- fullblock[i] - 0.25;
        rsh <- fullblock[i] + 0.25; 
        #polygon(x=c(lsh,rsh,rsh,lsh),y=c(-40,-40,15,15),border=NA,col="grey70");
    }
}
boxplots.dat <- boxplot(as.numeric(value)~bs+cs,data=sda,
                xlim=c(1,max(fullblock)),yaxt="n",xaxt="n",
                at=fullblock,lwd=1.5,cex.lab=1.75,ylim=c(-0.3,0.3),
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5,add=TRUE);
abline(h=0,lwd=0.8,lty="dotted");
text(x=0,y=0.26,cex=1.75,labels=expression(paste("(F) ",S['10,2'])),pos=4);
# ----------------------- Panel `g':
use <- which(last.gen[,1] == 2 & last.gen[,2] == 100);
crr <- cbind(last.gen[use,1:2],last.gen[use,4:5],last.gen[use,14]);
dat <- as.data.frame(cbind(crr[,5],crr[,3:4]));
colnames(dat) <- c("Corr","bs","cs");
sda <- melt(dat, id = c('bs','cs'));
sda <- sda[,-3];
fourblock <- c(1,1.8,2.6,3.4,4.2);
fullblock <- c(fourblock,fourblock+6,fourblock+12,fourblock+18);
boxplots.dat <- boxplot(as.numeric(value)~bs+cs,data=sda,type="n",
                xlim=c(1,max(fullblock)),yaxt="n",xaxt="n",
                at=fullblock,lwd=1.5,cex.lab=1.75,ylim=c(-0.3,0.3),
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5);
TST <- NULL;
for(j in unique(dat[,3])){
    for(i in unique(dat[,2])){
        sub <- which(dat[,2]==i & dat[,3]==j);
        cis <- simpleboot(dat[sub,1]);
        ttt <- cis[1] * cis[2];
        ttt <- max(c(ttt,0));
        ttt[ttt > 0] <-1;
        TST <- rbind(TST,c(j,i,cis,ttt));
    }
}
TST <- bcsort(TST);
for(i in 1:dim(TST)[1]){
    if(TST[i,5] == 1){
        lsh <- fullblock[i] - 0.25;
        rsh <- fullblock[i] + 0.25; 
        #polygon(x=c(lsh,rsh,rsh,lsh),y=c(-40,-40,15,15),border=NA,col="grey70");
    }
}
boxplots.dat <- boxplot(as.numeric(value)~bs+cs,data=sda,
                xlim=c(1,max(fullblock)),yaxt="n",xaxt="n",
                at=fullblock,lwd=1.5,cex.lab=1.75,ylim=c(-0.3,0.3),
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5,add=TRUE);
abline(h=0,lwd=0.8,lty="dotted");
axis(side=2,at=c(-0.2,0,0.2),cex.axis=1.5);
axis(side=1,at=c(fourblock[3],fourblock[3]+6,fourblock[3]+12,fourblock[3]+18),
     cex.axis=2,cex.lab=2,labels=c("0.00","0.01","0.02","0.03"));
text(x=0,y=0.26,cex=1.75,labels=expression(paste("(G) ",S['2,100'])),pos=4);
# ----------------------- Panel `h':
use <- which(last.gen[,1] == 2 & last.gen[,2] == 10);
crr <- cbind(last.gen[use,1:2],last.gen[use,4:5],last.gen[use,14]);
dat <- as.data.frame(cbind(crr[,5],crr[,3:4]));
colnames(dat) <- c("Corr","bs","cs");
sda <- melt(dat, id = c('bs','cs'));
sda <- sda[,-3];
fourblock <- c(1,1.8,2.6,3.4,4.2);
fullblock <- c(fourblock,fourblock+6,fourblock+12,fourblock+18);
boxplots.dat <- boxplot(as.numeric(value)~bs+cs,data=sda,type="n",
                xlim=c(1,max(fullblock)),yaxt="n",xaxt="n",
                at=fullblock,lwd=1.5,cex.lab=1.75,ylim=c(-0.3,0.3),
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5);
TST <- NULL;
for(j in unique(dat[,3])){
    for(i in unique(dat[,2])){
        sub <- which(dat[,2]==i & dat[,3]==j);
        cis <- simpleboot(dat[sub,1]);
        ttt <- cis[1] * cis[2];
        ttt <- max(c(ttt,0));
        ttt[ttt > 0] <-1;
        TST <- rbind(TST,c(j,i,cis,ttt));
    }
}
TST <- bcsort(TST);
for(i in 1:dim(TST)[1]){
    if(TST[i,5] == 1){
        lsh <- fullblock[i] - 0.25;
        rsh <- fullblock[i] + 0.25; 
        #polygon(x=c(lsh,rsh,rsh,lsh),y=c(-40,-40,15,15),border=NA,col="grey70");
    }
}
boxplots.dat <- boxplot(as.numeric(value)~bs+cs,data=sda,
                xlim=c(1,max(fullblock)),yaxt="n",xaxt="n",
                at=fullblock,lwd=1.5,cex.lab=1.75,ylim=c(-0.3,0.3),
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5,add=TRUE);
abline(h=0,lwd=0.8,lty="dotted");
axis(side=1,at=c(fourblock[3],fourblock[3]+6,fourblock[3]+12,fourblock[3]+18),
     cex.axis=2,cex.lab=2,labels=c("0.00","0.01","0.02","0.03"));
text(x=0,y=0.26,cex=1.75,labels=expression(paste("(H) ",S['2,10'])),pos=4);
# ----------------------- Panel `i':
use <- which(last.gen[,1] == 2 & last.gen[,2] == 2);
crr <- cbind(last.gen[use,1:2],last.gen[use,4:5],last.gen[use,14]);
dat <- as.data.frame(cbind(crr[,5],crr[,3:4]));
colnames(dat) <- c("Corr","bs","cs");
sda <- melt(dat, id = c('bs','cs'));
sda <- sda[,-3];
fourblock <- c(1,1.8,2.6,3.4,4.2);
fullblock <- c(fourblock,fourblock+6,fourblock+12,fourblock+18);
boxplots.dat <- boxplot(as.numeric(value)~bs+cs,data=sda,type="n",
                xlim=c(1,max(fullblock)),yaxt="n",xaxt="n",
                at=fullblock,lwd=1.5,cex.lab=1.75,ylim=c(-0.3,0.3),
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5);
TST <- NULL;
for(j in unique(dat[,3])){
    for(i in unique(dat[,2])){
        sub <- which(dat[,2]==i & dat[,3]==j);
        cis <- simpleboot(dat[sub,1]);
        ttt <- cis[1] * cis[2];
        ttt <- max(c(ttt,0));
        ttt[ttt > 0] <-1;
        TST <- rbind(TST,c(j,i,cis,ttt));
    }
}
TST <- bcsort(TST);
for(i in 1:dim(TST)[1]){
    if(TST[i,5] == 1){
        lsh <- fullblock[i] - 0.25;
        rsh <- fullblock[i] + 0.25; 
        #polygon(x=c(lsh,rsh,rsh,lsh),y=c(-40,-40,15,15),border=NA,col="grey70");
    }
}
boxplots.dat <- boxplot(as.numeric(value)~bs+cs,data=sda,
                xlim=c(1,max(fullblock)),yaxt="n",xaxt="n",
                at=fullblock,lwd=1.5,cex.lab=1.75,ylim=c(-0.3,0.3),
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5,add=TRUE);
abline(h=0,lwd=0.8,lty="dotted");
axis(side=1,at=c(fourblock[3],fourblock[3]+6,fourblock[3]+12,fourblock[3]+18),
     cex.axis=2,cex.lab=2,labels=c("0.00","0.01","0.02","0.03"));
text(x=0,y=0.26,cex=1.75,labels=expression(paste("(I) ",S['2,2'])),pos=4);

mtext(expression(paste("Cost of polyandry (",c[P],")")),
	outer=TRUE,side=1,line=3.75,cex=1.5);

mtext(expression(paste("Correlation between ",I[g]," and ",P[g]," genotypes")),
	outer=TRUE,side=2,line=3.25,cex=1.5);

dev.off();

# =====================================================================================

# =====================================================================================
# TODO: BUILD FIGURE S1-15                                                            #
# ====================================================================================#
setEPS(); # postscript below for final publication?
cairo_ps("trait_cors_soc.eps",family="Arial",height=7,width=7);
par(mar=c(0.2,0.2,0.2,0.2),oma=c(6,6,1,1),lwd=2);
use <- which(last.gen[,1] == -1 & last.gen[,2] == 100);
crr <- cbind(last.gen[use,1:2],last.gen[use,4:5],last.gen[use,14]);
dat <- as.data.frame(cbind(crr[,5],crr[,3:4]));
colnames(dat) <- c("Corr","Beta","cost");
sda <- melt(dat, id = c('Beta','cost'));
sda <- sda[,-3];
fourblock <- c(1,1.8,2.6,3.4,4.2);
fullblock <- c(fourblock,fourblock+6,fourblock+12,fourblock+18);
boxplots.dat <- boxplot(as.numeric(value)~Beta+cost,data=sda,xlim=c(0,fullblock[20]+2),
                at=fullblock,lwd=2.25,xaxt="n",yaxt="n",
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5,ylim=c(-0.2,0.2));
axis(side=2,at=c(-0.2,-0.1,0,0.1,0.2),cex.axis=2,cex.lab=2);
boxplots.dat <- boxplot(as.numeric(value)~Beta+cost,data=sda,xlim=c(0,fullblock[20]+2),
                at=fullblock,lwd=2.25,add=TRUE,xaxt="n",yaxt="n",
                col=c("white","grey70","grey60","grey45","grey30"),
                tick=c(2.2,6.4,10.6,14.8,19),cex.axis=1.5,ylim=c(-0.2,0.2));
abline(h=0,lty="dotted",lwd=0.8);
axis(side=1,at=c(fourblock[3],fourblock[3]+6,fourblock[3]+12,fourblock[3]+18),
     cex.axis=2,cex.lab=2,labels=c("0.00","0.01","0.02","0.03"));
mtext(expression(paste("Cost of polyandry (",c[P],")")),
	outer=TRUE,side=1,line=3.75,cex=2);
mtext(expression(paste("Correlation between ",I[g]," and ",P[g]," genotypes")),
	outer=TRUE,side=2,line=3.25,cex=1.5);

dev.off();
# ====================================================================================#

# =====================================================================================
# TODO: BUILD FIGURE S1-16                                                            #
# ====================================================================================#

setEPS(); # postscript below for final publication
cairo_ps("trait_cors_distr.eps",family="Arial",height=7,width=7);
par(mar=c(5,5,1,1),lwd=2);
hist(last.gen[,14],breaks=100,xlim=c(-0.35,0.35),col="grey40",main="",freq=FALSE,
     xlab=expression(paste("Correlation between ",I[g]," and ",P[g]," genotypes")),
     ylab="Probability density",cex.lab=1.5,cex.axis=1.5,lwd=2);

dev.off();

# ====================================================================================#

# =====================================================================================
# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX #
# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX #
# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX #
# =====================================================================================
# Bonus below -- not even intersting enough for SI, but keeping here just in case.    #
# =====================================================================================
# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX #
# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX #
# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX #
# =====================================================================================

# =====================================================================================
# TODO: LOOK AT INBREEDING ADJUSTMENT WHEN POLYANDRY TO ADJUST INBREEDING OCCURS      #
# ====================================================================================#
f.bypr2 <- function(pol,initial,additional){
    use      <- which(pol[,1] == initial & pol[,2] == additional);
    tus      <- pol[use,];
    negIn  <- NULL;
    evoada <- evo[evo[,1]==initial & evo[,2]==additional,];
    pedada <- ped[ped[,1]==initial & ped[,2]==additional,];
    matada <- mat[mat[,1]==initial & mat[,2]==additional,];
    count  <- 0;
    avoidn <- NULL;
    for(i in unique(pedada[,5])){
        tempped <- pedada[pedada[,5]==i,];
        saveped <- NULL;
        for(j in unique(tempped[,3])){
            for(k in unique(tempped[,4])){
                checkped <- tempped[tempped[,3]==j & tempped[,4]==k,];
                checkgW  <- mean(checkped[,8]);
                if(!is.na(checkgW) & checkgW < 0 & j > 0.2){
                    saveped <- rbind(saveped,checkped[1,1:5]);
                }
                if(!is.na(checkgW) & checkgW > 0 & j <= 0.2){
                    saveped <- rbind(saveped,checkped[1,1:5]);
                }
            }
        }
        avoidn <- rbind(avoidn,saveped);
        count  <- count + 1;
    }
    print("Finished finding inbreeding preference");
    negval <- NULL
    for(i in 1:dim(avoidn)[1]){
        tempchk <- which(tus[,4]==avoidn[i,3] & tus[,5]==avoidn[i,4] & tus[,3]==avoidn[i,5]);
        negval  <- rbind(negval,tus[tempchk,]);
    }
    print("Finished weeding out inbreeding preference");
    tus      <- negval;
    pan      <- tus[tus[,7]>1,];
    rm       <- which(pan[,8] > 1 | pan[,9] > 1); #Errors cleaning the file.
    if(length(rm) > 0){
        pan      <- pan[-rm,];
    }
    bypr <- NULL;
    for(i in unique(pan[,3])){
        for(j in unique(pan[,4])){
            for(k in unique(pan[,5])){
               calc <- which(pan[,3]==i & pan[,4]==j & pan[,5]==k);
               if(length(calc) > 0){
                   if(length(calc) > 4){
                       ppp  <- t.test(pan[calc,9]-pan[calc,8])$p.value;
                   }else{
                       ppp  <- NA;
                   }
                   bypr <- rbind(bypr,c(i,j,k,mean(pan[calc,9]-pan[calc,8]),ppp));
               } 
            }
        }
    }
    return(bypr);
}

setEPS(); # postscript below for final publication?
cairo_ps("Nmin1_mn2.eps",family="Arial",height=6,width=6);
par(mar=c(1,1,1,1),oma=c(3.5,3.5,0.1,3.5));
# ----------------------- Panel `a':
bp1 <- f.bypr(pol=pol,initial=-1,additional=100);
bp2 <- f.CIva(bypr=bp1);
bp2 <- bporder(bp2);
bp3 <- f.pomn(pol=pol,initial=-1,additional=100);
bp3 <- bporder(bp3);
bp4 <- f.bypr2(pol=pol,initial=-1,additional=100);
bp5 <- f.CIva(bypr=bp4);
bp5 <- bporder(bp5);
fourblock <- c(1,1.8,2.6,3.4,4.2);
fullblock <- c(fourblock,fourblock+6,fourblock+12,fourblock+18);
belowpts  <- -0.0045;
plot(x=fullblock[1],y=bp2[1,3],type="n",xlim=c(0,max(fullblock)+1),
     ylim=c(-0.011 + belowpts,0.014),xaxt="n",yaxt="n",xlab="",ylab="");
abline(h=0,lwd=0.9,lty="dotted");
poin <- 21;
for(i in 1:dim(bp2)[1]){
    if(poin == (21 + 5)){
        poin <- 21;
    }
    for(j in unique(bp2[,1])){
        for(k in unique(bp2[,2])){
            if(bp2[i,1] == j & bp2[i,2] == k){
                points(x=fullblock[i],y=bp2[i,3],pch=poin,cex=0.8,bg="black");
                arrows(x0=fullblock[i],x1=fullblock[i],y0=bp2[i,3],y1=bp2[i,4],
                       angle=90,length=0.025,lwd=1.50);
                arrows(x0=fullblock[i],x1=fullblock[i],y0=bp2[i,3],y1=bp2[i,5],
                       angle=90,length=0.025,lwd=1.50);
            }
        }
    }
    poin <- poin + 1;
}
axis(side=2,at=c(-0.01,0,0.01),cex.axis=1.75);
par(new = TRUE);
plot(x=fullblock[1],y=bp2[1,3],type="n",xlim=c(0,max(fullblock)+1),ylim=c(0,45),yaxt="n",
     xaxt="n",xlab="",ylab="");
polygon(x=c(0,fullblock[20]+1,fullblock[20]+1,0),y=c(0,0,8,8),border="black",col="grey80");
segments(x0=0.1,y0=4,x1=fullblock[20]+1,y1=4,lwd=1,lty="dotted");
for(i in 1:dim(bp3)[1]){
    for(j in unique(bp3[,1])){
        for(k in unique(bp3[,2])){
            if(bp3[i,1] == j & bp3[i,2] == k){
                lsh <- fullblock[i] - 0.25;
                rsh <- fullblock[i] + 0.25; 
                polygon(x=c(lsh,rsh,rsh,lsh),y=c(8,8,bp3[i,3]+8,
                        bp3[i,3]+8),border=NA,col="grey70");
                polygon(x=c(lsh,rsh,rsh,lsh),y=c(0,0,bp3[i,4]*8,
                        bp3[i,4]*8),border=NA,col="black");
            }
        }
    }
}
polygon(x=c(0,fullblock[20]+1,fullblock[20]+1,0),y=c(0,0,8,8),border="black",col=NA);
axis(side=4,at=c(8,18,28,38),labels=c("0","10","20","30"),cex.axis=1.5);
par(new = TRUE);
plot(x=fullblock[1],y=bp2[1,3],type="n",xlim=c(0,max(fullblock)+1),
     ylim=c(-0.011 + belowpts,0.014),yaxt="n",xaxt="n",xlab="",ylab="");
for(i in 1:dim(bp2)[1]){
    if(poin == (21 + 5)){
        poin <- 21;
    }
    for(j in unique(bp2[,1])){
        for(k in unique(bp2[,2])){
            if(bp2[i,1] == j & bp2[i,2] == k){
                points(x=fullblock[i],y=bp2[i,3],pch=poin,cex=0.8,bg="black");
                arrows(x0=fullblock[i],x1=fullblock[i],y0=bp2[i,3],y1=bp2[i,4],
                       angle=90,length=0.025,lwd=1.50);
                arrows(x0=fullblock[i],x1=fullblock[i],y0=bp2[i,3],y1=bp2[i,5],
                       angle=90,length=0.025,lwd=1.50);
                points(x=fullblock[i],y=bp5[i,3],pch=poin,cex=0.8,bg="white");
                points(x=fullblock[i],y=bp2[i,3],pch=poin,cex=0.8,bg="black");
            }
        }
    }
    poin <- poin + 1;
}
par(new = TRUE);
plot(x=fullblock[1],y=bp2[1,3],type="n",xlim=c(0,max(fullblock)+1),ylim=c(0,45),yaxt="n",
     xaxt="n",xlab="",ylab="");
polygon(x=c(0,fullblock[20]+1,fullblock[20]+1,0),y=c(0,0,8,8),border="black",col="grey80");
segments(x0=0.1,y0=4,x1=fullblock[20]+1,y1=4,lwd=1,lty="dotted");
for(i in 1:dim(bp3)[1]){
    for(j in unique(bp3[,1])){
        for(k in unique(bp3[,2])){
            if(bp3[i,1] == j & bp3[i,2] == k){
                lsh <- fullblock[i] - 0.25;
                rsh <- fullblock[i] + 0.25; 
                polygon(x=c(lsh,rsh,rsh,lsh),y=c(0,0,bp3[i,4]*8,
                        bp3[i,4]*8),border=NA,col="black");
            }
        }
    }
}
axis(side=1,at=c(fourblock[3],fourblock[3]+6,fourblock[3]+12,fourblock[3]+18),
     cex.axis=1.5,cex.lab=1.35,labels=c("0.00","0.01","0.02","0.03"));
box();

mtext(expression(paste("Cost of polyandry (",c[P],")")),
	outer=TRUE,side=1,line=1.75,cex=1.5);

mtext(expression(paste("Mean inbreeding adjustment (",k[adj],")")),
	outer=TRUE,side=2,line=1.75,cex=1.5);

mtext(expression(paste("Mean number of mates")),
	outer=TRUE,side=4,line=1.75,cex=1.5);

dev.off();

# ====================================================================================#

# =====================================================================================
# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX #
# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX #
# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX #
# =====================================================================================
#                  _____ _       _     _              _ 
#                 |  ___(_)_ __ (_)___| |__   ___  __| |
#                 | |_  | | '_ \| / __| '_ \ / _ \/ _` |
#                 |  _| | | | | | \__ \ | | |  __/ (_| |
#                 |_|   |_|_| |_|_|___/_| |_|\___|\__,_|
#                                                       
# =====================================================================================
# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX #
# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX #
# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX #
# =====================================================================================






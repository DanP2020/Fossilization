#This program will simulate speciation and extinction based on four variables. User will determine
#the number of lineages to begin with, the speciation rate, the extinction rate, and number of years
#over which the program will run.

evolution <- function(species.count=100,run.count=100,generation.count=30,extinction.rate=.01,speciation.rate=.99,overrun=1,spec.increment=.01,exti.increment=.01){
  
    holder.spec<-array(0,c(run.count,generation.count+1,overrun))
    holder.exti<-array(0,c(run.count,generation.count+1,overrun))
    holder.species<-array(0,c(run.count,generation.count+1,overrun))
    holder.counts<-array(0,c(run.count+1,7,overrun))
    dimnames(holder.counts)[[2]]<-c("Net.Spec","Net.Ext","Abs.Spec","Abs.Ext","Turn.Pulse","Turn.Pulse+1","Turn.Pulse+2")
    
    finallist<-list(0)            #Initiate list for overall matrix storage
    if(overrun > 1){
      pb<-txtProgressBar(min = 1, max = overrun, style = 3) #setup the progress bar
    }
    speciation.mat <- matrix(0,run.count,generation.count+1) #setup matrix for storage of speciation counts (col = generations, row = run #)
    extinction.mat <- matrix(0,run.count,generation.count+1) #setup matrix for storage of Extinction counts (col = generations, row = run #)
    species.mat <- matrix(0,run.count,generation.count+1) #setup matrix for storage of Total Species counts (col = generations, row = run #)
    
    spec.rate<-speciation.rate
    exti.rate<-extinction.rate
    
    for(d in 1:overrun){
    
      spec.rate<-speciation.rate - (d-1)*spec.increment
      exti.rate<-extinction.rate + (d-1)*exti.increment
      
    for(c in 1:run.count){
      
      loop.count<-species.count
      
      lineage <- matrix(0,loop.count,generation.count+1)   #Generates empty matrix for the addition of lineage splits
      lineage[,1]<-1:nrow(lineage)                            #Makes first column the lineage names  
      
      for(a in 1:generation.count){   #Sets up loop for cycling through evolutionary time
        for(b in 1:nrow(lineage)){    #Sets up loop for speciation, extinction, continuation determination
          chance<-runif(1,0,1)        #Generates a different random number for each species during each cycle
          
            if(lineage[b,a] != 0 & lineage[b,a] != 3){    #Sets up conditional - species that are already extinct do not go through the other conditional
              if(chance <= exti.rate){       #Conditional setup for whether a species will: Go extinct
                lineage[b,a+1] <- 0  
              } else if (chance >= spec.rate){   #: Speciate
                lineage[b,a+1] <- 2
                lineage<-rbind(lineage,rep(3,generation.count+1))
                lineage[nrow(lineage),1]<-lineage[b,1]
                lineage[nrow(lineage),a+1]<-1
              } else {                   #: Continue unimpeded to the next cycle
                lineage[b,a+1] <- 1
              }
          
          }
        }
      }
      
      xx<<-lineage
      
      species.mat[,1]<-species.count
      speciation.mat[,1]<-0
      extinction.mat[,1]<-0
      
     
      for(a in 1:generation.count){
        speciation.mat[c,a+1]<- sum(lineage[,a+1] == 2)    #counts number of speciations in each generation and inserts into speciation.mat
        species.mat[c,a+1]<- sum(lineage[,a+1] > 0 & lineage[,a+1] <3)       #counts number of total species in each generation and inserts into species.mat
        
      for(a in 1:generation.count){
        exti.counter<-0
        for(b in 1:length(lineage[,1])){
          if(a < length(lineage[1,])){
              if((lineage[b,a+1] == 0 & lineage[b,a] > 0)){
                exti.counter<-exti.counter+1
              }
            }
          }
        extinction.mat[c,a+1]<-exti.counter
        }
      }
    }
    
      holder.spec[,,d]<- speciation.mat
      holder.exti[,,d] <- extinction.mat
      holder.species[,,d] <- species.mat
      holder.counts[1,1,d]<-spec.rate
      holder.counts[1,2,d]<-exti.rate
      holder.counts[1,3,d]<-spec.rate
      holder.counts[1,4,d]<-exti.rate
      holder.counts[1,5,d]<-0
      
      
      for(a in 1:run.count){
        tempcount.spec<-0
        tempcount.exti<-0
        tempcount.abs.spec<-0
        tempcount.abs.exti<-0
        tempcount.turnover<-0
        tempcount.turnover1<-0
        tempcount.turnover2<-0
        
        for(b in 1:(generation.count)){
         
          if((species.mat[a,b+1]-species.mat[a,b])>=10){
            tempcount.spec<-tempcount.spec+1
          }
          
          if((species.mat[a,b]-species.mat[a,b+1])>=6){
            tempcount.exti<-tempcount.exti+1
          }
          if(speciation.mat[a,b+1]>=10){
            tempcount.abs.spec<-tempcount.abs.spec+1
          }
          
          if(extinction.mat[a,b+1]>=6){
            tempcount.abs.exti<-tempcount.abs.exti+1
          }
          
          if(extinction.mat[a,b+1]>=6 & speciation.mat[a,b+1]>=10){
            tempcount.turnover<-tempcount.turnover+1
          }
          
          if(b > 1 & b < generation.count){
            if(extinction.mat[a,b+1]>=6 & speciation.mat[a,b+1]>=10){
              tempcount.turnover1<-tempcount.turnover1+1
            }
            if(extinction.mat[a,b+2]>=6 & speciation.mat[a,b+1]>=10){
              tempcount.turnover1<-tempcount.turnover1+1
            }
            if(extinction.mat[a,b]>=6 & speciation.mat[a,b+1]>=10){
              tempcount.turnover1<-tempcount.turnover1+1
            }
          } else {
            if(extinction.mat[a,b+1]>=6 & speciation.mat[a,b+1]>=10){
              tempcount.turnover1<-tempcount.turnover1+1
            }
          }
          
          if(b > 2 & b < (generation.count-1)){
            if(extinction.mat[a,b+1]>=6 & speciation.mat[a,b+1]>=10){
              tempcount.turnover2<-tempcount.turnover2+1
            }
            if(extinction.mat[a,b+2]>=6 & speciation.mat[a,b+1]>=10){
              tempcount.turnover2<-tempcount.turnover2+1
            }
            if(extinction.mat[a,b]>=6 & speciation.mat[a,b+1]>=10){
              tempcount.turnover2<-tempcount.turnover2+1
            }
            if(extinction.mat[a,b+3]>=6 & speciation.mat[a,b+1]>=10){
              tempcount.turnover2<-tempcount.turnover2+1
            }
            if(extinction.mat[a,b-1]>=6 & speciation.mat[a,b+1]>=10){
              tempcount.turnover2<-tempcount.turnover2+1
            }
          } else {
            if(extinction.mat[a,b+1]>=6 & speciation.mat[a,b+1]>=10){
              tempcount.turnover2<-tempcount.turnover2+1
            }
          }
          
          holder.counts[a+1,1,d]<-tempcount.spec
          holder.counts[a+1,2,d]<-tempcount.exti
          holder.counts[a+1,3,d]<-tempcount.abs.spec
          holder.counts[a+1,4,d]<-tempcount.abs.exti
          holder.counts[a+1,5,d]<-tempcount.turnover
          holder.counts[a+1,6,d]<-tempcount.turnover1
          holder.counts[a+1,7,d]<-tempcount.turnover2
          
        }
      }
      holder.counts[1,5,d]<-(length(which(holder.counts[2:run.count+1,5,d] > 1)))/run.count
      holder.counts[1,6,d]<-(length(which(holder.counts[2:run.count+1,6,d] > 1)))/run.count
      holder.counts[1,7,d]<-(length(which(holder.counts[2:run.count+1,7,d] > 1)))/run.count
                  
  if(overrun > 1){
      setTxtProgressBar(pb,d)          #progress bar ticker
  }
    }
    
  finallist<-list(speciation = holder.spec,extinction = holder.exti,species = holder.species,counts = holder.counts)
  
  
if(overrun > 1){    
  close(pb)
}
  
  graph.mat<-matrix(0,overrun,5)
  
  for(a in 1:overrun){
    graph.mat[a,]<-holder.counts[1,1:5,a]
  }
  
  graph.mat<-graph.mat*100
  
  plot.new()
  barplot(graph.mat[,5],names.arg = graph.mat[,4],main="Turnover Pulses in Simulated Data",xlab="Speciation Rate",ylab="Percent of runs with two or more Turnovers",cex.names = .8,col = rainbow(100), ylim=c(0,100),axes=FALSE,las=2)
  box()
  axis(2,at=seq(0,100,10))
  abline(h=95,lty = 2,lwd=2)
  
  return(finallist)
    
}
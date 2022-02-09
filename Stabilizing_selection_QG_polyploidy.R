# This is the simulation code and functions used for the paper "The evolution of the additive variance of a trait under stabilizing selection after autopolyploidization"
# You can change the variance of environmental effects manually at Line 205 (set to 0.05 by default)

# A function to introduce mutations within newly formed gametes

Mutation<-function(haplotype, U, L, var.add.eff)
{
  
  Nb_mut<-rpois(1,U)
  
  position_mut<-c(sample(1:L, Nb_mut, replace = F))
  
  if(Nb_mut != 0){
    
    for(i in 1:length(position_mut)){
      
      haplotype[position_mut[i]]= haplotype[position_mut[i]] + rnorm(n = 1, mean = 0, sd = var.add.eff^0.5)
      
    }
    
  }
  
  return(haplotype)
  
}

# A function to select the reproducers 

Selection <-function(Npop, phenotype,om_2,opt)
{
  
  selected.parent = 0
  list.parent = c(NULL)
  
  repeat{
    
    fit.max = max(exp(-((phenotype-opt)^2)/(2*om_2)))
    
    rand.repro1 = sample(1:Npop, 1)
    
    fit.samp = exp(-((phenotype[rand.repro1]-opt)^2)/(2*om_2))
    
    if( (fit.samp/fit.max) >= runif(1, 0, 1) ){
      
      selected.parent = selected.parent + 1
      
      list.parent[selected.parent]=rand.repro1
      
      # Selfing or outcrossing
      
      repeat{
        
        rand.repro2 = sample(1:Npop, 1)
        fit.samp = exp(-((phenotype[rand.repro2]-opt)^2)/(2*om_2))
        
        if( rand.repro1 != rand.repro2 && (fit.samp/fit.max) >= runif(1, 0, 1) ){
          
          selected.parent = selected.parent + 1
          
          list.parent[selected.parent]=rand.repro2
          
          break
        }
        
      }
      
    }
    
    if(selected.parent == 2*Npop) break
    
  }
  
  return(list.parent)
  
}

# A function to generate the genome of offsprings 

Genome.offspring <-function(list.parent, L, Npop, U, genome, step, var.add.eff,gen)
{
  
  if(step == 0){
    
    for(n in 1:Npop){
      
      haplo.temp1 = c(NULL)
      haplo.temp2 = c(NULL)
      
      for(i in 1:L){
        
        add_loc1 = c(genome[2*list.parent[2*n], i], genome[2*list.parent[2*n] - 1, i])
        haplo.temp1[i] = sample(add_loc1, 1)
        
        add_loc2 = c(genome[2*list.parent[2*n - 1], i], genome[2*list.parent[2*n - 1] - 1, i])
        haplo.temp2[i] = sample(add_loc2, 1)
        
      }
      
      haplo.temp.mut1 = Mutation(haplo.temp1, U, L, var.add.eff)  
      haplo.temp.mut2 = Mutation(haplo.temp2, U, L, var.add.eff)
      
      if(n == 1){genotype.off = haplo.temp.mut1
      
      genotype.off = rbind(genotype.off, haplo.temp.mut2)
      
      }
      
      else{genotype.off = rbind(genotype.off, haplo.temp.mut1, haplo.temp.mut2)}
    }
    
  }
  
  else{
    
    if(gen == 1){
      
      for(n in 1:Npop){
        
        haplo.temp1 = c(NULL)
        haplo.temp2 = c(NULL)
        haplo.temp3 = c(NULL)
        haplo.temp4 = c(NULL)
        
        for(i in 1:L){
          
          add_loc1 = c(genome[2*list.parent[2*n], i], genome[2*list.parent[2*n] - 1, i])
          rand_1 = sample(add_loc1, 2, replace = F)
          haplo.temp1[i] = rand_1[1]
          haplo.temp2[i] = rand_1[2]
          
          add_loc2 = c(genome[2*list.parent[2*n - 1], i], genome[2*list.parent[2*n - 1] - 1, i])
          rand_2 = sample(add_loc2, 2, replace = F)
          haplo.temp3[i] = rand_2[1]
          haplo.temp4[i] = rand_2[2]
          
        }
        
        haplo.temp.mut1 = Mutation(haplo.temp1, U, L, var.add.eff)  
        haplo.temp.mut2 = Mutation(haplo.temp2, U, L, var.add.eff)
        haplo.temp.mut3 = Mutation(haplo.temp3, U, L, var.add.eff)  
        haplo.temp.mut4 = Mutation(haplo.temp4, U, L, var.add.eff)
        
        if(n == 1){genotype.off = haplo.temp.mut1
        
        genotype.off = rbind(genotype.off, haplo.temp.mut2, haplo.temp.mut3, haplo.temp.mut4)
        
        }
        
        else{genotype.off = rbind(genotype.off, haplo.temp.mut1, haplo.temp.mut2, haplo.temp.mut3, haplo.temp.mut4)}
        
      }
    }
    
    else{
      for(n in 1:Npop){
        
        haplo.temp1 = c(NULL)
        haplo.temp2 = c(NULL)
        haplo.temp3 = c(NULL)
        haplo.temp4 = c(NULL)
        
        for(i in 1:L){
          
          add_loc1 = c(genome[4*list.parent[2*n], i], genome[4*list.parent[2*n] - 1, i], genome[4*list.parent[2*n] -2, i], genome[4*list.parent[2*n] - 3, i])
          rand_1 = sample(add_loc1, 2, replace = F)
          haplo.temp1[i] = rand_1[1]
          haplo.temp2[i] = rand_1[2]
          
          add_loc2 = c(genome[4*list.parent[2*n - 1], i], genome[4*list.parent[2*n - 1] - 1, i], genome[4*list.parent[2*n - 1] - 2, i], genome[4*list.parent[2*n - 1] - 3, i])
          rand_2 = sample(add_loc2, 2, replace = F)
          haplo.temp3[i] = rand_2[1]
          haplo.temp4[i] = rand_2[2]
          
        }
        
        haplo.temp.mut1 = Mutation(haplo.temp1, U, L, var.add.eff)  
        haplo.temp.mut2 = Mutation(haplo.temp2, U, L, var.add.eff)
        haplo.temp.mut3 = Mutation(haplo.temp3, U, L, var.add.eff)  
        haplo.temp.mut4 = Mutation(haplo.temp4, U, L, var.add.eff)
        
        if(n == 1){genotype.off = haplo.temp.mut1
        
        genotype.off = rbind(genotype.off, haplo.temp.mut2, haplo.temp.mut3, haplo.temp.mut4)
        
        }
        
        else{genotype.off = rbind(genotype.off, haplo.temp.mut1, haplo.temp.mut2, haplo.temp.mut3, haplo.temp.mut4)}
      }
    }
  }
  
  return(genotype.off)
  
}

# A function to compute the phenotypic values of individuals

Phenotype.offspring <-function(Genotype, Npop)
{
  
  pheno.off = c(NULL)
  
  for(i in 1:Npop){
    
    pheno.off[i] = Genotype[i] + rnorm(1, 0, 1) # genotype + environmental effects
    
  }
  
  return(pheno.off)
  
}

# A function to compute the genotypic values of individuals

genotype.offspring <-function(Genome, Npop, step, dosage){
  
  if(step == 0){
    
    geno.off = c(NULL)
    
    for(i in 1:Npop){
      
      geno.off[i] = sum(Genome[2*i,]) + sum(Genome[(2*i) - 1,])
      
    }
    
  }
  
  else{
    
    geno.off = c(NULL)
    
    for(i in 1:Npop){
      
      geno.off[i] = (1+dosage)*(sum(Genome[4*i,]) + sum(Genome[(4*i) - 1,])+ sum(Genome[(4*i) - 2,])+ sum(Genome[(4*i) - 3,]))
      
    }
  }
  
  return(geno.off)
  
}

# A function to compute check is population is at equilibrium

is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol

# A function to decompose the genetic variance of the quantitative trait under study

variances <- function(genome, genotype,Npop, L, step, dosage){
  
  mean.genotype = mean(genotype)
  genetic.variance = (sum((genotype - mean.genotype)^2))/Npop
  
  par_locus<-0
  
  if(step == 0){
    for(i in 1:L){
      
      phenotype_loc = c(NULL)
      
      for(j in 1:Npop){
        
        phenotype_loc[j] = genome[2*j, i] + genome[(2*j) - 1, i]
        
      }
      
      mean.phenotype = mean(phenotype_loc)
      par_locus = par_locus + (sum((phenotype_loc - mean.phenotype)^2))/Npop
      
    }
  }
  
  else{
    for(i in 1:L){
      
      phenotype_loc = c(NULL)
      
      for(j in 1:Npop){
        
        phenotype_loc[j] = (1+dosage)*(genome[4*j, i] + genome[(4*j) - 1, i] + genome[(4*j) - 2, i] + genome[(4*j) - 3, i])
        
      }
      
      mean.phenotype = mean(phenotype_loc)
      par_locus = par_locus + (sum((phenotype_loc - mean.phenotype)^2))/Npop
      
    }
  }
  
  cov.genetic = genetic.variance - par_locus
  
  summary.variance = c(genetic.variance, par_locus, cov.genetic)
  
  return(summary.variance)
  
}

# A function to compute the mean frequency of the ancestral beneficial allele

frequence <- function(genome,L,Npop, step){ 
  
  freq.all0=c(NULL)
  
  if(step == 0){
    
    for(i in 1:L){
      
      nb.all1 = 0
      
      for(j in 1:(2*Npop)){
        if(genome[j,i] == 0){nb.all1 = nb.all1 + 1}
      }
      
      freq.all1 = nb.all1/(2*Npop)
      
      freq.all0[i]=freq.all1
    }
  }
  
  else{
    
    for(i in 1:L){
      
      nb.all1 = 0
      
      for(j in 1:(2*Npop)){
        if(genome[j,i] == 0){nb.all1 = nb.all1 + 1}
      }
      
      freq.all1 = nb.all1/(4*Npop)
      
      freq.all0[i]=freq.all1
    }
  }
  
  return(freq.all0)
}

# A function to compute the mean effect of alleles in the population

mean.all.effect<- function(genome){
  
  add_eff = mean(as.numeric(row.names(table(genome))))
  
  return(add_eff)
  
}

# The main simulation function

Simulation.programme <- function(Number.sim, L, var.add.eff, dosage, Npop, U, om_2){
  
  Nsimul = 0
  
  repeat{
    
    Nsimul = Nsimul + 1  
    
    mean.fitness = c(NULL)
    
    opt = 0
    
    om.2 = om_2
    
    ### Initialisation
    
    geno.init = c(rep(0, times=L))
    
    for(i in 1:(2*Npop)){
      
      if(i == 1)(genome = geno.init)
      else{genome = rbind(genome, geno.init)}
      
    }
    
    phenotype = c(rnorm(Npop, 0, 1))
    
    gen = 0
    
    repeat{
      
      gen = gen + 1
      
      fitness = exp(-(phenotype^2)/(2*om_2))
      
      mean.fitness[gen] = mean(fitness)
      
      Parent.sel = Selection(Npop, phenotype, om_2,opt)
      genome = Genome.offspring(Parent.sel, L, Npop, U, genome, step = 0, var.add.eff,gen)
      genotype = genotype.offspring(genome, Npop, step = 0, dosage)
      phenotype = Phenotype.offspring(genotype, Npop)
      
      if( is.wholenumber(gen/1000) == TRUE && (gen/1000) > 1 ){
        
        mean.1<-mean(mean.fitness[(gen-999):gen])
        mean.2<-mean(mean.fitness[(gen-1999):(gen-1000)])
        
        if( abs(1 - (mean.1/mean.2)) <= 0.01){break}
        
      }
      
    }
    
    ## Step 2
    
    gen = 0
    
    dyn_fit = c(NULL)
    dyn_varg = c(NULL)
    dyn_var_add = c(NULL)
    dyn_cov = c(NULL)
    dyn_freq0 = c(NULL)
    dyn_all_effect = c(NULL)
    
    for(g in 1:2001){
      
      gen = gen + 1  
      
      if(gen == 1){
        
        fitness = exp(-(phenotype^2)/(2*om_2))
        
        dyn_fit[gen] = mean(fitness)
        
        var_all = variances(genome, genotype, Npop, L, step = 0, dosage)
        freq = frequence(genome, L, Npop, step = 0)
        
        dyn_varg[gen] = var_all[1]
        dyn_var_add[gen] = var_all[2]
        dyn_cov[gen] = var_all[3]
        dyn_freq0[gen] = mean(freq, na.rm = T)
        dyn_all_effect[gen] = mean.all.effect(genome)
        
        Parent.sel = Selection(Npop, phenotype, om_2,opt)
        genome = Genome.offspring(Parent.sel, L, Npop, U, genome, step = 1, var.add.eff,gen)
        genotype = genotype.offspring(genome, Npop, step = 1, dosage)
        phenotype = Phenotype.offspring(genotype, Npop)
        
      }
      
      else{
        
        fitness = exp(-(phenotype^2)/(2*om_2))
        
        dyn_fit[gen] = mean(fitness)
        
        var_all = variances(genome, genotype, Npop, L, step = 1, dosage)
        freq = frequence(genome, L, Npop, step = 1)
        
        dyn_varg[gen] = var_all[1]
        dyn_var_add[gen] = var_all[2]
        dyn_cov[gen] = var_all[3]
        dyn_freq0[gen] = mean(freq, na.rm = T)
        dyn_all_effect[gen] = mean.all.effect(genome)
        
        Parent.sel = Selection(Npop, phenotype, om_2,opt)
        genome = Genome.offspring(Parent.sel, L, Npop, U, genome, step = 1, var.add.eff,gen)
        genotype = genotype.offspring(genome, Npop, step = 1, dosage)
        phenotype = Phenotype.offspring(genotype, Npop)
        
      }
      
    }
    
    if(Nsimul == 1){
      
      recap.fit = dyn_fit
      recap.vg = dyn_varg
      recap.va = dyn_var_add
      recap.cov = dyn_cov
      recap.freq = dyn_freq0
      recap.all = dyn_all_effect
      
    }
    
    else{
      
      recap.fit = rbind(recap.fit, dyn_fit)
      recap.vg = rbind(recap.vg, dyn_varg)
      recap.va = rbind(recap.va, dyn_var_add)
      recap.cov = rbind(recap.cov, dyn_cov)
      recap.freq = rbind(recap.freq, dyn_freq0)
      recap.all = rbind(recap.all, dyn_all_effect)
      
    }
    
    if(Nsimul == Number.sim) break
    
  }
  
  write.table(recap.fit, file=paste("fit_d",dosage,"_L",L,"_om",om_2,"_N",Npop,".txt",sep=""),row.names = F, dec = ".")
  
  write.table(recap.vg, file=paste("vg_d",dosage,"_L",L,"_om",om_2,"_N",Npop,".txt",sep=""),row.names = F, dec = ".")
  
  write.table(recap.va, file=paste("va_d",dosage,"_L",L,"_om",om_2,"_N",Npop,".txt",sep=""),row.names = F, dec = ".")
  
  write.table(recap.cov, file=paste("cov_d",dosage,"_L",L,"_om",om_2,"_N",Npop,".txt",sep=""),row.names = F, dec = ".")
  
  write.table(recap.freq, file=paste("freq_d",dosage,"_L",L,"_om",om_2,"_N",Npop,".txt",sep=""),row.names = F, dec = ".")
  
  write.table(recap.all, file=paste("all_d",dosage,"_L",L,"_om",om_2,"_N",Npop,".txt",sep=""),row.names = F, dec = ".")
  
}

# Example of simulation parameters set

Simulation.programme(10, 50, 0.05, 0, 1000, 0.01, 1)

# This is the simulation code and functions used for the paper "The evolution of the additive variance of a trait under stabilizing selection after autopolyploidization"
# You can change the variance of environmental effects manually at Line 205 (set to 0.05 by default)

# A function to introduce mutations within newly formed gametes

Mutation<-function(haplotype, U, L, var.add.eff)
{
  
  Nb_mut<-rpois(1,U) ## The number of mutations to be done 
  
  position_mut<-c(sample(1:L, Nb_mut, replace = F)) ## The positions of the mutations on the haploid genome
  
  if(Nb_mut != 0){
    
    ## If you have at least 1 mutation, do the mutations
    
    for(i in 1:length(position_mut)){
      
      haplotype[position_mut[i]]= haplotype[position_mut[i]] + rnorm(n = 1, mean = 0, sd = var.add.eff^0.5)
      
    }
    
  }
  
  return(haplotype)
  
}

# A function to select the reproducers 

Selection <-function(Npop, phenotype,om_2,opt)
{
  
  selected.parent = 0 ## A counter to see how many parents have been selected
  list.parent = c(NULL) ## A vector to store the the pedigree of parents
  
  repeat{
    
    fit.max = max(exp(-((phenotype-opt)^2)/(2*om_2))) ## Inferring the highest fitness value in the population at generation t
    
    rand.repro1 = sample(1:Npop, 1) ## Select one random individual
    
    fit.samp = exp(-((phenotype[rand.repro1]-opt)^2)/(2*om_2)) ## Inferring the fitness of the selected indidivual
    
    if( (fit.samp/fit.max) >= runif(1, 0, 1) ){
      
      ## Selection process, the ration of fit.samp/fit.max need to be higher than the runif number for the parent to be selected
      
      selected.parent = selected.parent + 1 ## Add one parent to the counter
      
      list.parent[selected.parent]=rand.repro1 ## Add the coordinate of the parent to the list 
      
      # Outcrossing
      
      repeat{
        
        rand.repro2 = sample(1:Npop, 1) ## Sample a second random individual
        fit.samp = exp(-((phenotype[rand.repro2]-opt)^2)/(2*om_2))  ## Infer the fitness of the second random individual
        
        if( rand.repro1 != rand.repro2 && (fit.samp/fit.max) >= runif(1, 0, 1) ){
          
          ## If they are two different individuals and fitness is high enough, reproduction is allowed
          
          selected.parent = selected.parent + 1 ## Add +1 to the counter of selected parents
          
          list.parent[selected.parent]=rand.repro2 ## Store the identity of the parent
          
          break
        }
        
      }
      
    }
    
    if(selected.parent == 2*Npop) break ## Break when you have enough parents
    
  }
  
  return(list.parent) ## Return the list of selected parents
  
}

# A function to generate the genome of offsprings 

Genome.offspring <-function(list.parent, L, Npop, U, genome, step, var.add.eff,gen)
{
  
  if(step == 0){
    
    ## Step 0 is a fully diploid population
    
    for(n in 1:Npop){
      
      ## Empty vectors for the haploid gametes of parents
      
      haplo.temp1 = c(NULL)
      haplo.temp2 = c(NULL)
      
      for(i in 1:L){
        
        ## add_loc1 is the alleles of the first parent at a given loci and haplo.temp1 is sampling randomly one of the two alleles
        
        add_loc1 = c(genome[2*list.parent[2*n], i], genome[2*list.parent[2*n] - 1, i])
        haplo.temp1[i] = sample(add_loc1, 1)
        
        ## the same for add_loc2
        
        add_loc2 = c(genome[2*list.parent[2*n - 1], i], genome[2*list.parent[2*n - 1] - 1, i])
        haplo.temp2[i] = sample(add_loc2, 1)
        
      }
      
      ## Two haploid genomes from each parent in wich we are introducing mutations
      
      haplo.temp.mut1 = Mutation(haplo.temp1, U, L, var.add.eff)  
      haplo.temp.mut2 = Mutation(haplo.temp2, U, L, var.add.eff)
      
      ## Building the new genome file with the genomes of the offsprings
      
      if(n == 1){genotype.off = haplo.temp.mut1
      
      genotype.off = rbind(genotype.off, haplo.temp.mut2)
      
      }
      
      else{genotype.off = rbind(genotype.off, haplo.temp.mut1, haplo.temp.mut2)}
    }
    
  }
  
  else{
    
    ## Just the same as before but for tetraploid individuals
    
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
  
  # An empty vector of phenotypic value
  
  for(i in 1:Npop){
    
    pheno.off[i] = Genotype[i] + rnorm(1, 0, 1) # genotype + environmental effects
    
    # Adding an environmental effect to the genotypic value to have the phenotypic value
    
  }
  
  return(pheno.off)
  
}

# A function to compute the genotypic values of individuals

genotype.offspring <-function(Genome, Npop, step, dosage){
  
  if(step == 0){
    
    ## Step 0 is for a fully diploid population
    
    geno.off = c(NULL)
    
    ## An empty vector of genotypic values
    
    for(i in 1:Npop){
      
      ## With additivity, genotype is just the sum of all values stored in the genome fill
      
      geno.off[i] = sum(Genome[2*i,]) + sum(Genome[(2*i) - 1,])
      
    }
    
  }
  
  else{
    
    ## The same for a fully tetraploid population
    
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
  
  Nsimul = 0  ## Is counting the number of simulations performed per parameter set
  
  repeat{
    
    Nsimul = Nsimul + 1  ## Add +1 to the count at each new simulation 
    
    mean.fitness = c(NULL)  ## Vector to stock population fitness
    
    opt = 0 ## The phenotypic optimum before the environmental change
    
    om.2 = om_2 ## The strength of stabilizing selection
    
    ### Initialisation
    
    geno.init = c(rep(0, times=L)) ## Generaiting the null haplotypes with only zero values at all loci
    
    for(i in 1:(2*Npop)){
      
      ## Npop is population size, and merge the null haplotypes to generation N individuals genetically similar at the beginning of simulations
      
      if(i == 1)(genome = geno.init)
      else{genome = rbind(genome, geno.init)}
      
    }
    
    phenotype = c(rnorm(Npop, 0, 1)) ## A first vector of phenotypic values (genotype = 0, and a random environmental effect with rnorm)
    
    gen = 0 ## Initiating the counting of generations 
    
    repeat{
      
      gen = gen + 1 ## Adding one generation to the counter
      
      fitness = exp(-(phenotype^2)/(2*om_2)) ## Computing the fitness of individuals based on phenotype-fitness function
      
      mean.fitness[gen] = mean(fitness) ## Inferring the mean fitness value of the population at generation t
      
      ## To understand the functions, have a look to the annotations in each function 
      
      ## To prepare the new generation we are making 1) Selection of the parents
      ## 2) Generating the genome of offsprings
      ## 3) Inferring the genotypes of offsprings
      ## 4) Inferring the phenotypes of offsprings
      
      Parent.sel = Selection(Npop, phenotype, om_2,opt) # 1)
      genome = Genome.offspring(Parent.sel, L, Npop, U, genome, step = 0, var.add.eff,gen) # 2)
      genotype = genotype.offspring(genome, Npop, step = 0, dosage) # 3)
      phenotype = Phenotype.offspring(genotype, Npop) # 4)
      
      ## Check if you have run enough generation
      ## And check if we reached the equilibirum for genetic diversity
      
      if( is.wholenumber(gen/1000) == TRUE && (gen/1000) > 1 ){
        
        mean.1<-mean(mean.fitness[(gen-999):gen]) ## Variance of the first 1000 generations
        mean.2<-mean(mean.fitness[(gen-1999):(gen-1000)]) ## Variance of the second 1000 generations
        
        if( abs(1 - (mean.1/mean.2)) <= 0.01){break} ## Have a look if it's deviate to much or not
        
      }
      
    }
    
    ## Step 2
    
    # We initiate the evolution of polyploidy
    
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
    
    # If we did enough simulation, break the repaet function
  }
  
  # Write the results in the working directory as .txt files
  
  write.table(recap.fit, file=paste("fit_d",dosage,"_L",L,"_om",om_2,"_N",Npop,".txt",sep=""),row.names = F, dec = ".")
  
  write.table(recap.vg, file=paste("vg_d",dosage,"_L",L,"_om",om_2,"_N",Npop,".txt",sep=""),row.names = F, dec = ".")
  
  write.table(recap.va, file=paste("va_d",dosage,"_L",L,"_om",om_2,"_N",Npop,".txt",sep=""),row.names = F, dec = ".")
  
  write.table(recap.cov, file=paste("cov_d",dosage,"_L",L,"_om",om_2,"_N",Npop,".txt",sep=""),row.names = F, dec = ".")
  
  write.table(recap.freq, file=paste("freq_d",dosage,"_L",L,"_om",om_2,"_N",Npop,".txt",sep=""),row.names = F, dec = ".")
  
  write.table(recap.all, file=paste("all_d",dosage,"_L",L,"_om",om_2,"_N",Npop,".txt",sep=""),row.names = F, dec = ".")
  
}

# Example of simulation parameters set

Simulation.programme(10, 50, 0.05, 0, 1000, 0.01, 1)

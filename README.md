# CoExpansionValidation
A R package for validating ABC co-expansion inference 

Tests of synchronous community expansion are commonly conducted using hierarchical Approximate Bayesian Computation (hABC), a statistical framework for inferring the degree of concordance across species. However, this framework is often used without demonstrating adequate performance.

Below we will show how to use this __CoExpansionValidation__ R package to assess the performance of genetic datasets for detecting co-expansion events across species. It will go through  an example step-by-step, showing how to generate one pseudo-observed dataset, and use ABC to "infer" the number of co-expansion events. Obviously, to apply this type of performance assessment in a empirical setting, many replicates of pseudo-observed datasets needs to be generated for inference.


## __Installation__
  
external programs that need to be installed before running:
  
 *  [BayeSSC - Serial Simcoal](http://www.stanford.edu/group/hadlylab/ssc/) 
 *  [msReject module in msBayes](http://msbayes.sourceforge.net/) 

The BayeSSC is used for simulation, and msReject is used for inference. While functions in this package does not directly use [hBayeSSC](https://github.com/UH-Bioinformatics/hBayeSSC), steps of hyperstat calculation and running msReject are set to be identical to [hBayeSSC](https://github.com/UH-Bioinformatics/hBayeSSC). Future development will expand on other types of summary statistics.

To install this package using devtools:

```r
devtools::install_github("huatengh/CoExpansionValidation",build_vignettes = TRUE)
```


```{r setup}
library(CoExpansionValidation)
options(stringsAsFactors = F)
```


 
## __1.Generating a pseudo-observed dataset__  
  
### i) Determine the number of co-expansion events  
  
Two options:  

* user specifies a list of species and assigns species to events randomly or evenly:  

```{r specify the number of co-expansion events}
species<-1:10 # or a vector of species names, can be characters
coevents<-2
#assigning species randomly
species.assignment<-assign_species_to_events(species = species,nco.events = coevents,even = F)
species.assignment
#assigning species evenly
species.assignment<-assign_species_to_events(species = species,nco.events = coevents,even = T)
species.assignment
```
  
* user specifies the species' names and the _alpha_ for Dirichlet process, and simulate the number of co-expansion events  

```{r simulate the number of co-expansion events}
species<-1:10
alpha<-1.1
species.assignment<-generate_coevent_number_dirichlet(species = species,alpha = alpha,maxevent=4)
species.assignment
coevents<-length(species.assignment)
```

### ii) Sampling the time for co-expansion events

User can provide the time for each event (time from present), or randomly sample the time from a range. Buffer can be added around co-expansion events, such that they are at least certain distance away from each other. No requirement on the unit of time yet--below user can provide species' generation time to turn them into number of generations.

```{r generate the expansion time}
# if user provides the time, make sure its length equals 
# the coevents specified above
exp.time<-c(30000,50000) 

#randomly draw time from a range
time.range<-c(30000,50000)
buffer<-5000
exp.time<-generate_cotime_with_buffer(time.range = time.range,nco.events =coevents,buffer =  buffer)
exp.time
```
User can get a vector of the sampled time for each species, easier to write out to a file

```{r get the expansion time for each species}
x<-species_exp_time(species.assignment,exp.time)
x
```

### iii) Simulate the "observed" summary statistics with BayeSSC 

Here, we provide a wrapper function for running BayeSCC from R. BayeSSC itself has very detail instructions on [their webpage](https://web.stanford.edu/group/hadlylab/ssc/index.html#IO).   User need to give:

* the path to the BayeSCC executable (absolute or relative to current working directory; including the executable name)
  
And then, there are two options:   
  
a) a configuration file: a table (TAB delimited file) consists of columns with following header names

__Column name__     |	__Description__
------------------- | ---------------
species(1)	        | Name of species
locinum             | Number of loci
nsam 	              | Number of haplotypes to be simulated
nsites    	        | Number of base pairs
tstv    	          | % transitions
gamma   	          | Gamma distribution for substitution rate heterogeneity, two numbers separated by space
mutation rate(2) 	  | The locus mutation rate per generation
gen(3)   	          | The generation length
Ne                  | The Effective population size
popratio            | The ratio between ancestral and current population size, <1 for expansion and >1 for shrinkage 
eventtime.generation(4) | The time of the historical event in unit of generation(optional)

  (1): species can appear in multiple rows if having multiple types of loci; using one row for one locus at a time works as well. However, all loci in one species share the same generation time, population size and expansion/shrinkage history. For these columns, function below will only read  from the first row of a species.  
  (2): mutation rate (and population size, and other number columns) can take BayeSCC-style prior distributions. For example, {U:1,299} for uniform distribution between 1 and 299. see [BayeSCC webpage](https://web.stanford.edu/group/hadlylab/ssc/index.html#IO) for details.  
  (3): generation length needs to be in the same unit as the time range and buffer provided before. That is, if MYA was used, then here would need the generation time in unit of million years.  
  (4): optional, if not provided, will be calculated automatically with the generation time and previously simulated events' time.

```{r run BayeSCC from a configuration file}
path_to_bayessc<-"BayeSSC.exe" #windows version, linux or mac starts with ./
#check the conf data frame --an example of the configuration table-- included in the package
head(conf)
#read your own configuration file:
#configurefile<-"test.conf"
#conf<-read.table(configurefile,header=T,sep="\t",stringsAsFactors = F)

#code below will simulate summary statistics with BayeSSC,see the function's help for details 
simulatedobs<-runbayeSSC_with_conf(BayeSSCallocation = path_to_bayessc,conf = conf,prefix = 'temp',species.assignment = species.assignment,exp.time = exp.time,intern=F)

#the package comes with a simulatedobs data frame, user can check how it looks like
colnames(simulatedobs)
head(simulatedobs[,1:5])
```



b) if all loci in all species can use the same .par file for BayeSCC (all loci in a species are of the same type) 

* .par file for BayeSCC. The format is fixed for only modelling one expansion event for each species.User can write out the templatepar file included in this package and edit it. Do not add or delete lines from the template. For the line under "1 historical event", the time of the event can set to be any integer; the wrapper function will swap in the generated (or set) expansion time from above. Note that "the number of loci" is actually the length in bp for DNA. Users provide other parameters-- population size, sample size and etc.  
* the number of loci need to be generated per species 
* generation time (if not set, default to 1)  

Below is a template file, already scanned in as part of the package data, templateparfile 

> BayeSSC//Number of population samples
1  
//Population sizes  
{U:500000,500000}  
//Sample sizes  
40  
//Growth rates  
0  
//Number of migration matrices : If 0 : No migration between demes  
0  
//Historical event format:  
1 historical event  
399137 0 0 1 0.01 0 0  
//mutation rate  
{U:0.00008,0.00008}  
//Number of independent loci  
800  
//Data type, tvts  
DNA 0.33  
//Gamma distribution for mutation  
0 0  
>   

```{r run BayeSCC from par file}
path_to_bayessc<-"BayeSSC.exe" #windows version, linux or mac starts with ./
#show the templateparfile included in this package
templateparfile
# you can write this template out to a file and edit it with any text editor
cat(paste(templateparfile,sep = '',collapse = '\n'),"\n",sep='',file="test.par")
bayessc_par_file<-'test.par'
nloci<-10
gen<-1 # or can be a vector with generation time for each species
#extract information from the par file to a configuration table
conf<-par_to_config(bayessc_par_file,species,nloci,gen)
#run BayeSSC with the configuration table
simulatedobs<-runbayeSSC_with_conf(BayeSSCallocation = path_to_bayessc,conf = conf,prefix = 'partemp',species.assignment = species.assignment,exp.time = exp.time,intern=TRUE)
```
### iv) calculate hyperstats with hBaySSC 
     
We calculate hyperstats from the "observed" statistics. In our paper, this step was performed with [hBayeSSC](https://github.com/UH-Bioinformatics/hBayeSSC), but it depends on python 2.X which is no longer supported. This package provides a conversion of the code to R. Here, we need


* the simulated summary statistics from previous step

Following code will calcualte hyperstats and write it out to a file that can be used for the third step of running msReject.


```{r hyperstats from hBayeSSC}
#if you write the observed summary stats into a file previously
obsfile<-"temp_obs_file"
simulatedobs<-read.table(obsfile,sep="\t",header=T,stringsAsFactors = F)

#calculate hyperstat across all loci
hyperstat<-calculate_hyperstat(simulatedobs)

#make a vector containing these hyperstats together the 'real' number of events 
#and expansion time for each species 
a<-rep(0,3)
a[1]<-'temp' #or any name you want to give to this pseudo-observed dataset
a[2]<-length(exp.time) #true number of events
a[3]<-length(species) #total number of species
names(a)<-c("uid","nevent","nspecies")
species.time<-species_exp_time(species.assignment = species.assignment,exp.time = exp.time)
a<-c(a,species.time)
hyperstat<-c(a,hyperstat)

#write the hyperstat to a file
outfile<-"temp_hyperstat_file"
cat(paste(hyperstat,sep='',collapse = "\t"),"\n",sep='',file = outfile)
    
```

## __2.ABC simulation with BayeSSC__  
ABC simulation is very similar to generating test data. The biggest difference is that the event time should be autogenerated according to a prior rather than directly specified by user.  

For the prior of the number expansion event, hBayeSSC used a flat distribution-- the number of co-expansion species are uniformly distributed from 1 to n (the total number of species) . Here, we adopted the [PyMsbayes-style](http://joaks1.github.io/PyMsBayes/) prior, using a gamma distribution for the alpha parameter of the dirichelete process. User need to decide on the two parameters for the gamma distribution: concentrationShape and consentrationScale. Check [prior selection in PyMsbayes](http://joaks1.github.io/PyMsBayes/tutorials/selecting-priors.html) for how to select these two parameters.  

Similar to step 1, this step needs:

* path to bayeSSC  
* the time range for possible expansion events  
* buffer between any two events  
* number of replicates  
* Either     
  + configuration file  
  or if all species can use the same par file for BayeSSC  
  + list of species  
  + number of loci per species
  + generation time
  + .par file 
  
```{r ABC simulation with configuration file}
path_to_bayessc="BayeSSC.exe"
concentrationShape=20.0
concentrationscale=0.15
time.range<-c(30000,50000)
buffer<-500
npod<-10
conf<-conf
#running ABC simulation
reference.table<-ABC_simulation_with_conf(npod=npod,conf=conf,time.range=time.range,buffer=buffer,concentrationscale=concentrationscale,concentrationShape=concentrationShape,BayeSSCallocation=path_to_bayessc,prefix='temp',do.parallel=2,write.reference.file = F)

#the package comes with a toy example of reference.table
colnames(reference.table)
head(reference.table[,1:5])
```

```{r ABC simulation from par file}
path_to_bayessc="BayeSSC.exe"
concentrationShape=20.0
concentrationscale=0.15
time.range<-c(30000,50000)
buffer<-500
npod<-100
# you can write the template out to a file and edit it with any text editor
cat(paste(templateparfile,sep = '',collapse = '\n'),"\n",sep='',file="test.par")
bayessc_par_file<-'test.par'
nloci<-10
species<-1:10
gen<-1 # or can be a vector with one generation time for each species
conf<-par_to_config(bayessc_par_file,species,nloci,gen)
reference.table<-ABC_simulation_with_conf(npod=npod,conf=conf,time.range=time.range,buffer=buffer,concentrationscale=concentrationscale,concentrationShape=concentrationShape,BayeSSCallocation=path_to_bayessc,prefix='temp',do.parallel=2)

```

## __3.Generating Inference for Test Data__  

### __i) Running msReject to get posterior samples__  

For installing/compiling msReject see [msbayes webpage](http://msbayes.sourceforge.net/). User need to provide
* the "observed" hyperstats  
* the simulated reference data table  
* the path to msReject  
* sampling tolerance  
* prefix to the posterior file  

```{r run MsReject}
path_to_msreject<-"./msReject"
samplingtolerance<-0.5 #just a toy example, should be much smaller for real runs
prefix<-'temp'
posterior<-run_msreject(hyperstat=hyperstat,reference.table=reference.table,MsRejectallocation=path_to_msreject,samplingtolerance=samplingtolerance,prefix=prefix)
```
### __ii) Infer the number of co-expansion events__  

With the posterior file after running msReject, we can use the VGAM, abc and locfit package to do the final acceptance and parameter estimation, some scripts are on the [hBayeSSC webpage](https://github.com/UH-Bioinformatics/hBayeSSC#msreject-module). Here are some codes for just looking at the mode of event number in the posterior samples.

```{r the mode of number of events}
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}
getmode(posterior$nevent)
```

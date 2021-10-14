# SexCoal
This is a coalescent simulation for pseudo-autosomal region of a sex chromosome evolving with a single potential locus under sexually antagonistic selection.

The code is first detailed in Kirkpatrick and Guerrero (2014), while this version adds demographic change to a single population and is detailed in Dagilis et al (in prep.)

The code can be compiled using a standar c++ compiler, and requires the boost library.


# Compiling the code

To compile the code, make sure you have the Boost library available. Then simply run:

```
c++ *.cpp -o SaSimSnps
```

# Parameters necessary to run the code

The compiled code takes a single parameter - an integer $n$ which determines the parameter file called $inputn$.

The parameter file consists of the following lines

```
Number of run
Population size
Mutation rate
Recombination rate
Migration rate between populations (currently non-functional)
Number of epochs - each epoch rescales effective population sizes; n of these
Epoch breakpoints (in generations); n-1 of these
Epoch scaling factors (the relative size of the population at each epoch, compared to the modern population size); n of these
Ages of the sexually antagonistic site and the Y, separated by a space. set to 0 assumes infinite age, no sweep
Location of the SDR site in rho, >= 0. (recommend keeping this at 0)
Location of the sexually antagonistic site, in rho
Male to female recombination rate ratio
Frequency of the sexually antagonistic site on the Y
Frequency of the sexually antagonistic site on the X
Frequency of the sexually antagonistic site on the Y previous to Y sweep
Buffer, in rho, to ignore around the SDR and sexually antagonistic site (avoids infinite coalescent times under some parameters, a value of 0.001 or so is sufficient)
Sampling scheme: 0- Only use the number of chromosomes from the list below - set populations, X or Y status and sexually antagonistic locus at random; 1- Only use the number of populations and  ; 2- Use chromosomes listed below, but set status at antagonistic site at random using the frequencies on X and Y specified
Location of neutral sites - only used to initiate the size of the chromosome simulated, current version places mutations at random sites
Chromosomes sampled, each listed as "a b c", where a=0 or 1, for population 0 or 1, b=0 or 1, for the X or Y, and c= 0 or 1 for status at sexually antagonistic locus 
```

For example, the below is a parameter file to simulate 14 X and 14 Y chromosomes of the Japan Sea stickleback, with a sweep on the Y 306000 generations ago.

```

```


If the above file is named input1, you can then run the code simply by running:

```
SaSimSnps 1
```
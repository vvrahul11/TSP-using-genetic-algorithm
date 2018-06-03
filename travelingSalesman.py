import random as rnd
import math
from copy import deepcopy
import matplotlib.pyplot as plt

class individual:
    def __init__(self):
        """this creates an individual"""
        self.my_genes = list()                    ## create an empty list to store the gene 
        self.my_fitness = 10000000000000

    def addNewGenes(self,genes):
        for gene in genes:                        ## begin a loop to create a random number 'r' in each step ##
            r = rnd.random()
            if r < 0.5:                           ## if r<0.5 
                self.my_genes.append(gene)        ## Add an item to the end of the list ##
            else:
                self.my_genes = [gene] + self.my_genes[0:] ## Add an item to the beginnning of the list ##
    def hasAllGenes(self,genes):
        for gene in genes:
            if gene not in self.genes:
                return False
        return True

    def replaceGenes(self,genes):
        self.my_genes = genes

    def showMe(self):
        print self.my_genes,self.my_fitness       ## print the genes and its corresponding distance ##

    def returMySequence(self):
        return '-'.join([str(an_integer) for an_integer in self.my_genes])

    def myFitness(self,city_distances):
        my_fitness = list()                       ## Define an empty array for the distance ##
        for i in range(0,len(self.my_genes)-1):   ## Begin a for loop ##
            path_address = tuple([self.my_genes[i],self.my_genes[i+1]]) ## store the path address for city1 and city2 ##
            my_fitness.append(city_distances[path_address])             ## append it to the empty list defined 1st ##
        self.my_fitness = sum(my_fitness)                               ## find the sum of distance for each chromosome##

    def rtnFitness(self):                         ## to check the fitness ##
        return self.my_fitness
            


    def returnMutatedGenes(self,probability_of_mutation):
        for i in range(0,len(self.my_genes)-1):   ## Begin a for loop ##            
            r = rnd.random()
            genes = deepcopy(self.my_genes)
            if r <= probability_of_mutation:               ## create a random number and check it is > or < than the probability of mutation ##o
                gene_1_pos = rnd.randint(0,len(genes)-1)   ## pick one positions within the chromosome for mutation ##
                gene_2_pos = rnd.randint(0,len(genes)-1)   ## pick another positions within the chromosome for mutation ##
                gene_1 = genes[gene_1_pos]
                gene_2 = genes[gene_2_pos]

                genes[gene_1_pos] = gene_2## assign the second positon value to the second positon ##
                genes[gene_2_pos] = gene_1## assign the first position value to the first position ##
        return genes

class individuals:                             ## Individual class ##
    def __init__(self):
        self.my_individuals = list()
        self.best_individual = list()#something... will get replaced anyway
    def returnNumberOfIndividuals(self):               ## for returning the number of individual from my_individual ##
        return len(self.my_individuals)
    def __str__(self):
        return self.best_individual.returMySequence()

    def returnRandomIndividual(self):
        i = int(rnd.random() * len(self.my_individuals))
        return self.my_individuals[i]

    def generate(self,genes,number_of_individuals):
        for i in range(0,number_of_individuals):
            an_individual = individual()
            an_individual.addNewGenes(genes)
            self.my_individuals.append(an_individual) ## call the individual class ##
        self.best_individual = self.my_individuals[0]     ## Assume first individual chromosome is the best ##
    def findBestIndividual(self):                       ## find the best individual ##
        for i in self.my_individuals:               ## Run a loop through all the chromosomes ##
            if i.rtnFitness() < self.best_individual.rtnFitness():   ## if the best chromosome assumed is less in fitness to the next chromosome in the population, ##
                self.best_individual = i
    def addIndividual(self,individual):
        self.my_individuals.append(individual)
        self.best_individual = self.my_individuals[0]     ## Assume first individual chromosome is the best ##
    def evaluateFitnesses(self,city_distances):           ## to evaluate the fitness ie. to find the total distance travelled for all gene.
        for i in range(0,len(self.my_individuals)):
            an_individual = self.my_individuals[i]
            an_individual.myFitness(city_distances)       ## Call the myFitness function and pass city_distances
    def geneSequenceCount(self):    
        sequence_library = dict()
        for individual in self.my_individuals:
            sequence = individual.returMySequence()
            if sequence not in sequence_library.keys():
                sequence_library[sequence] = 0
            sequence_library[sequence] += 1
        for key,value in sequence_library.items():
            print key,value
    def crossover(self,indiv_1,indiv_2,probability_of_crossover):  ## for implementing the crossing over operator ##
        genes_1 = indiv_1.my_genes                                  
        genes_2 = indiv_2.my_genes
        genes_1 = list()                                           ## two empty lists are created
        genes_2 = list()
        for i in range(0,len(genes_1)-1):                          ## begin a loop for the entire length of genes_1 ##
            r = rnd.random()
            if r >= probability_of_crossover:                      ## if ranom number r is less than p.c.o ##
                gene1 = genes_1[i]                                 ## do cross over ##
                gene2 = genes_2[i]
                genes_1.append(gene1)
                genes_2.append(gene2)

            if r < probability_of_crossover:
                for j in range(i,len(genes_1)-1):
                    genes_1 = genes_2[j]
                    genes_2 = genes_1[j]

        ok_genepool = genes_1
        final_individual_1 = individual();final_individual_1.replaceGenes(genes_1)
        final_individual_2 = individual();final_individual_2.replaceGenes(genes_2)
        if final_individual_1.hasAllGenes(ok_genepool) and final_individual_2.hasAllGenes(ok_genepool):
            return [final_individual_1,final_individual_2]
        else:
            return [indiv_1,indiv_2]

class Generations:
    def __init__(self):
        self.populations = list()
    def generateNewGeneration(self,old_generation,disruptive_probabilities,selection_option="best"):
        new_generation = individuals()                                                               ## create an object and call the individual class ##
        number_of_individuals_in_old_generation = old_generation.returnNumberOfIndividuals()
        if selection_option == "best":
            best_individual = old_generation.best_individual
            number_of_individuals = old_generation.returnNumberOfIndividuals()
            half_number_of_individuals = number_of_individuals/2
            ## Crossover
            fraction_to_cross = 0.25
            for i in range(0,old_generation.returnNumberOfIndividuals()):
                if (i/1.0 % fraction_to_cross) == 0.0:
                    random_individual = old_generation.returnRandomIndividual()
                    new_individuals = old_generation.crossover(random_individual,best_individual,disruptive_probabilities['crossover'])

            ## Mutation
            fraction_to_mutate = 0.25
            for i in range(0,old_generation.returnNumberOfIndividuals()):
                if (i/1.0 % fraction_to_mutate) == 0.0:
                    r = rnd.random()
                    if r < disruptive_probabilities['mutation']:
                        new_individual = individual()
                        new_individual.replaceGenes(best_individual.returnMutatedGenes(disruptive_probabilities['mutation']))
                    else:
                        new_individual = best_individual
                    new_generation.addIndividual(new_individual)
        while len(new_generation.my_individuals) > number_of_individuals_in_old_generation:
            print 'removed an individual'
            throw_away_element = new_generation.my_individuals.pop()
        while len(new_generation.my_individuals) < number_of_individuals_in_old_generation:
            print 'added an individual'
            new_generation.addIndividual(best_individual)
        return new_generation

    def run(self,number_of_generations,number_of_individuals,disruptive_probabilities,cities,distances):
        new_generation = individuals()                                  ## call the function individual ##
        new_generation.generate(cities,number_of_individuals)           ## generate the new genes ##
        new_generation.evaluateFitnesses(distances)                     ## evaluate the genes ##
        new_generation.findBestIndividual()                             ## find the best individual ##
        print 'best sequence in generation',new_generation
        print 'sequence frequencies',new_generation.geneSequenceCount()
        print 20*'-'
        self.populations.append(new_generation)
        for i in range(1,number_of_generations):
            the_old_generation = self.populations[len(self.populations)-1]     ## in each step take the last generation created ##
            the_new_generation = self.generateNewGeneration(                   ## pass the arguments for creating a new generation ##
                the_old_generation,
                disruptive_probabilities,
                "best")
            the_new_generation.evaluateFitnesses(distances)
            the_new_generation.findBestIndividual()
            print 'best sequence in generation',the_new_generation
            print 'sequence frequencies',
            the_new_generation.geneSequenceCount()
            print 20*'-'
            self.populations.append(the_new_generation)
    def showMe(self):                                                       ## this function plot the final best individual plot
        x = list()
        y = list()
        for i in range(0,len(self.populations)):
            x.append(i)                                                      ## x-axis the length of the population ##
            y.append(self.populations[i].best_individual.rtnFitness())       ## y-axis the distance taken to travel ##
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(x,y)
        plt.show()


#execution
cities = [i for i in range(0,10)]                   # Define the cities ##
city_positions = dict()                          # Create an empty dictionary ##
xymax = 1000.0                                   ## xy coordinate max= 1000 ##
for city in cities:
    city_positions[city] = [rnd.random()*xymax,rnd.random()*xymax]    ## create a random number and multiply it with xymax to generate the city coordinate position##

city_distances = dict()                         ## create an empty dictionary ##
for city1 in city_positions.keys():             ## Loop through the geneated coordinates for city1 ##
    for city2 in city_positions.keys():         ## Loop through the geneated coordinates for the city2 ##
        key = tuple([city1,city2])              ## create a key for city1 and city 2 ##
        coordinates1 = city_positions[city1]    ## get the value of the city1 position ##
        coordinates2 = city_positions[city2]    ## get the value of the city2 position ##
        city_distances[key] = math.sqrt(
            (coordinates1[0] - coordinates2[0])**2 + (coordinates1[1] - coordinates2[1])**2)   ## find the distance between the city1 and city2 ##

number_of_individuals = 10                     ## this is for creating 100 initial chromosomes ##
number_of_generations=500
myGenerations = Generations()
myGenerations.run(
    number_of_generations,
    number_of_individuals,
    {'mutation':0.5,
    'crossover':0.8},
    cities,
    city_distances)
myGenerations.showMe()

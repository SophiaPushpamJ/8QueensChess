# -*- coding: utf-8 -*-
"""
Created on Fri Jun 26 15:33:28 2020

@author: Sophia
"""

import random
import numpy as np
#from pprint import pprint
import math

class Chess:
    def __init__(self, N, maxFitness, initialPopulationCount):
        self.N = N
        self.maxFitness = maxFitness
        self.initialPopulationCount = initialPopulationCount
        self.initalPopList = []
        
    def createChromosome(self):
        '''
        This function creates a chromosome of length N with no repeated position for queen

        Returns
        -------
        chromosome.

        '''
        self.randomChromosome = random.sample(range(0,self.N),self.N)        
        return self.randomChromosome
    
    #createChromosome()
    
    def transposeChromosome(self,chromosome):
        '''
        This function creates a 2D array of the chromosome

        Returns
        -------
        2D array.

        '''
        visual = []
        
        for ind in range(self.N):
            twoDRep = []
            for j in range(self.N):
                if j == chromosome[ind]:
                    twoDRep.append("X")
                else:
                    twoDRep.append("0")
            visual.append(twoDRep)
        
        transposeList = []
        for value in range(0,self.N):
            newList = []
            for j in range(0,self.N):
                newList.append(visual[j][value])
            transposeList.append(newList)
        return transposeList
    
    #transposeChromosome([0,1,2,3,4,5,6,7])
    
    def visualizeChromosome(self,chromosome):
        '''
        This function gives a visual representation of the input chromosome.

        Parameters
        ----------
        chromosome : TYPE
            DESCRIPTION.

        Returns
        -------
        prints the  string
        None
            DESCRIPTION.

        '''
        visual = self.transposeChromosome(chromosome)
        for row in visual:
            for element in row:
                print(element,end=" | ")
            print("")
        return "" 
    
    #visualizeChromosome([0,1,2,3,4,5,6,7])
            
    def calcFitnessScore(self,chromosome):
        '''
        Calculate the fitness value of each chromosomes by reducing the row clash and diagonal clash

        Returns
        -------
        fitnessScore

        '''
        rowClashCounter = 0
        for gene in chromosome:
            if chromosome.count(gene) > 1:
                rowClashCounter += 1
        
        diagonalClashCounter = 0
        transposeList =  self.transposeChromosome(chromosome)
        arr = np.array(transposeList)
        diags = [arr[::-1,:].diagonal(i) for i in range(-arr.shape[0]+1,arr.shape[1])] #upper diagnol value
        diags.extend(arr.diagonal(i) for i in range(arr.shape[1]-1,-arr.shape[0],-1)) #lower diagnol value
        possibleDiagList = [n.tolist() for n in diags]
        
        for combination in possibleDiagList:
            if combination.count("X") > 1 :
                diagonalClashCounter += combination.count("X")
     
        fitnessScore = maxFitness - (rowClashCounter + diagonalClashCounter)
        return fitnessScore

    #calcFitnessScore([0,1,2,3,4,5,6,7])    
    
    def genInitalPopulation(self):
        '''
        This function creates the initial population

        Returns
        -------
        the entire initial population

        '''
        
        for i in range(0,self.initialPopulationCount):
            chromosome = self.createChromosome()
            fitness = self.calcFitnessScore(chromosome)
            self.initalPopList.append({"chromosome":chromosome, "fitScore":fitness})
        self.initalPopList = sorted(self.initalPopList, key = lambda x:x["fitScore"], reverse = True)
        return self.initalPopList
    
    #genInitalPopulation()
    
    def chooseParent(self,population):
        '''
        This function chooses 2 parents randomly from the entire population

        Parameters
        ----------
        population : TYPE
            DESCRIPTION.

        Returns
        -------
        parents.

        '''
        #parentList = []
        indexValList = random.sample(range(0,len(population)),2)
        #print(population[2])
        parentList = [population[element] for element in indexValList]
           
        return parentList
    
    #chooseParent(genInitalPopulation())
   
    def performCrossover(self, parentList):
        '''
        This function performs crossover between parents with the maximum fitness score and generates 2 children

        Returns
        -------
        children
        

        '''
        parent1 = parentList[0]
        parent2 = parentList[1]
        indexLower = math.floor(self.N/4)
        indexUpper = math.floor(self.N/2)
        coValChild1 = parent1[indexLower:indexUpper+2] #crossOver value for Child1
        rmValPar2 = [element for element in parent2 if element not in coValChild1]
        child1 = rmValPar2[0:indexLower] + coValChild1 + rmValPar2[indexLower:]#major parent1 minor parent2
        coValChild2 = parent2[indexLower:indexUpper+2] #crossOver value for Child2
        rmValPar1 = [element for element in parent1 if element not in coValChild2]
        child2 = rmValPar1[0:indexLower] + coValChild2 + rmValPar1[indexLower:]#major parent2 minor parent1
        childList = [child1,child2]
        return childList
    
    #performCrossover([[1,2,3,4,5,6,7,0],[0,5,4,8,7,9,3,2]])
    
    def performMutation(self,childChromosome):
        '''
        This function will mutate the a child chromosome

        Returns
        -------
        mutated child sequence.

        '''
        neChild = random.sample(childChromosome, len(childChromosome))
        return neChild
    
    #performMutation([1,2,3,4,5,6,7,0])
 
    

if __name__ == '__main__':
    QueenCount = 8 # setting value as 8 since it is a 8*8 chess queen problem
    InitialPopCount = 100 #intial population count is set as 1000
    maxFitness = int((QueenCount * (QueenCount-1))/2)
    chess = Chess(QueenCount,maxFitness,InitialPopCount) #8 queens hence passing value as 8 and passing the initial population count as 1000
    population = chess.genInitalPopulation()
    maxFitValue = population[0]
    Generation = 0
    finalSequence = []
    
    
    def Check():
        globals()
        maxFitValue = population[0]
        if maxFitValue["fitScore"] == maxFitness:
            finalSequence.append(maxFitValue["chromosome"])
            print("finalSequence:{} generated in Generation: {}".format(finalSequence[0],Generation))
            chess.visualizeChromosome(finalSequence[0])
            return False
        else:
            #print("MaxSeq:{} generated in Generation: {}".format(population[0]["chromosome"],Generation))
            return True  
    
    while Check():
        parentList = []
        childList = []
        newChildList = []
        parentList = chess.chooseParent(population)
        parent = [parentList[i]["chromosome"] for i in range(0,len(parentList))]
        childList = chess.performCrossover(parent)
        childList.append(chess.performMutation(childList[random.randint(0,1)]))
        del population[-3:]
        for i in range(len(childList)):
            population.append({"chromosome":childList[i],"fitScore":chess.calcFitnessScore(childList[i])})
            population = sorted(population, key = lambda x:x["fitScore"], reverse = True)
        Generation +=1


from Chapter1 import *
import math
from random import randint, uniform


def MotifEnumeration(Dna, k, d):
	Patterns = []

	for line in Dna:
		Neighborhood = []
		for i in xrange(len(line) - k + 1):
			Neighborhood.append(Neighbors(line [i: i + k], d))
		Patterns.append([item for sublist in Neighborhood for item in sublist])
	
	Suburb = set(Patterns[0]) 
	for s in Patterns[1:]:
		Suburb.intersection_update(s)

	Purge = [s for s in Suburb] # Make list out of set fo printing function

	return Purge

# MotifEnumeration(DNA, 5, 2)

def Entropy(Varone,Vartwo):
	entropy1 = Varone * math.log(Varone,2)
	entropy2 = Vartwo * math.log(Vartwo,2) 
	entropy = entropy1 + entropy2 
	print -round(entropy, 4)

def PatternMotifs(k): # Create all possible patterns given a k size
	nucletoids = ["C", "T", "G", "A"]

	if k == 1:
		return nucletoids

	Patterns = []
	SuffixPattern = PatternMotifs(k - 1) # Recursive function.

	for Text in SuffixPattern:
		for x in nucletoids:
			Patterns.append(x + Text) 

	return Patterns

# print PatternMotifs(3)

def DistanceBetweenPatternAndString(Pattern, Dna):
	k = len (Pattern)
	distance = 0

	for text in Dna:
		HamDistance = k + 1

		for i in xrange(len(text) - k):
			H = HammingDistance(Pattern, text[i: i + k])
			if HamDistance > H:
				HamDistance = H

		distance += HamDistance

	return distance

# print DistanceBetweenPatternAndString(Pat, Dn)

def MedianString(Dna, k):
	distance = 2 ** 32 # Simulate infinity. Max int value in 32 bits.
	Median = ""

	for pattern in PatternMotifs(k):
		d = DistanceBetweenPatternAndString(pattern, Dna)
		if distance > d:
			distance = d
			Median = pattern

	return Median

def StructureProfile(Grid): # Orders information from Sceptic into a list of lists. Creates a Grid.
	for i in xrange(len(Grid)):
		Grid[i] = [float(x) for x in Grid[i].split(" ")]
	return Grid

def LineProbability(text, Grid):
	Prob = 1
	nucletoids = {nuc:index for index,nuc in enumerate("ACGT")}
	for index,value in enumerate(text):
		Prob *= Grid [nucletoids [value] ] [index]
	return Prob

def ProfileProbableKmer(Text, k, Grid): #Returns the Pattern that best suits a Probability Grid.
	LineProb = [-1, None]
	for i in xrange(len(Text) - k + 1):
		P = LineProbability(Text[i: i + k], Grid)
		if LineProb[0] < P:
			LineProb = [P, Text [i: i + k] ]

	return LineProb[1]

def SingleGrid(text): #Creates a grid with the probability for that single kmer.
	Grid = []
	for i in range(4):
		Grid.append([1] * len(text))

	PositionInGrid = 0

	for i in text:
		if i == "A" or  i == "a":
			Grid[0][PositionInGrid] += 1
			PositionInGrid += 1
			
		elif i == "C" or i =="c":
			Grid[1][PositionInGrid] += 1
			PositionInGrid += 1
			
		elif i == "G" or i =="g":
			Grid[2][PositionInGrid] += 1
			PositionInGrid += 1
		else:
			Grid[3][PositionInGrid] += 1
			PositionInGrid += 1
	
	return Grid

def ProbabilityProfile(list): # Creates a Probability Grid according to a list of kmers.
	if len(list) == 1:
		SingleGrid(list[0])

	else:

		Grid = SingleGrid(list[0])

		for text in list[1:]:
			SecondGrid = SingleGrid(text)
			for i in xrange(len(Grid)):
				Grid[i] = [x + y for x,y in zip(Grid[i],SecondGrid[i])]
			
	for i in xrange(len(Grid[0])):
		average = Grid[0][i] + Grid[1][i] + Grid[2][i] + Grid[3][i]
		if average != 0:
			Grid[0][i] /= float(average)
			Grid[1][i] /= float(average)
			Grid[2][i] /= float(average)
			Grid[3][i] /= float(average)

	return Grid

def ProfileCount(list): # Creates a Probability Grid according to a list of kmers.
	if len(list) == 1:
		return SingleGrid(list[0])

	else:

		Grid = SingleGrid(list[0])

		for text in list[1:]:
			SecondGrid = SingleGrid(text)
			for i in xrange(len(Grid)):
				Grid[i] = [x + y for x,y in zip(Grid[i],SecondGrid[i])]

	return Grid

def profile_with_pseudocounts(motifs):
	'''Returns the profile of the dna list motifs.'''
	prof = []
	for i in xrange(len(motifs[0])):
		col = ''.join([motifs[j][i] for j in xrange(len(motifs))])
		prof.append([float(col.count(nuc)+1)/float(len(col)+4) for nuc in 'ACGT'])
	return prof

def score(motifs):
	'''Returns the score of the dna list motifs.'''
	score = 0
	for i in xrange(len(motifs[0])):
		motif = ''.join([motifs[j][i] for j in xrange(len(motifs))])
		score += min([HammingDistance(motif, homogeneous*len(motif)) for homogeneous in 'ACGT'])
	return score


def ScoreMotif(list): #Score the similarity of all kmers
	
	Grid = SingleGrid(list[0])

	for text in list[1:]:
		SecondGrid = SingleGrid(text)
		for i in xrange(len(Grid)):
			Grid[i] = [x + y for x,y in zip(Grid[i],SecondGrid[i])]

	Score = 0


	Grid2 = []
	for i in range(len(list[0])):
		Grid2.append([0] * 4)
	

	for j in xrange(len(list[0])):
		for n in xrange(4):
			Grid2[j][n] = Grid[n][j]


	for i in xrange(len(Grid2)):
		Grid2[i].remove(max(Grid2[i]))
		Score += sum(Grid2[i])
	
	return Score

def GreedyMotifSearch(Dna, k, t):
	BestMotifs = [x[0: k] for x in Dna] # List of first kmers in each DNA line

	for i in xrange(len(Dna[0]) - k + 1): # All kmers in the first line of DNA
		Motif = [Dna[0][i: i + k]]

		for j in xrange(1, t):
			Grid = ProbabilityProfile(Motif)	
			Motif.append(ProfileProbableKmer(Dna[j], k, Grid))

		if ScoreMotif(Motif) < ScoreMotif(BestMotifs):
			BestMotifs = Motif
			

	return BestMotifs

def RandomizedMotifSearch(Dna, k, t):
	Motifs = []
	length = len(Dna[0]) - k

	for text in Dna:
		r = randint(0, length)
		Motifs.append(text[r: r + k])
	
	BestMotifs = Motifs

	for i in xrange(2 * 64):
		Profile = ProbabilityProfile(Motifs)
		Motifs = [ProfileProbableKmer(x, k, Profile) for x in Dna]

		if ScoreMotif(Motifs) < ScoreMotif(BestMotifs):
			BestMotifs = Motifs
		else:
			 return BestMotifs

def NumRuns(Dna, k, N):
	Motifs = []
	for i in xrange(N):
		Motifs.append(RandomizedMotifSearch(Dna, k, 1))

	BestMotifs = Motifs[0]
	BestScore = ScoreMotif(Motifs[0])

	for j in Motifs:
		Score = ScoreMotif(j)

		if Score < BestScore:
			BestMotifs = j
			BestScore = Score

	return BestMotifs

def Weighted_choice(choices):
   total = sum(choices.values())
   r = uniform(0,total)
   upto = 0
   for c, w in choices.iteritems():
      if upto + w >= r:
        return c
      upto += w
   assert False, "Shouldn't get here"

def GibbsSampler(Dna, k, t, N): # N is the number of randomstarts
	length = len(Dna[0]) - k
	BestMotifs = []

	for s in xrange(N):

		#Randomly select kmer motifs
		Motifs = []
		RandomDnaNumber = {}
		counter = 0

		for strand in Dna:

			r = randint(0, length)
			RandomDnaNumber[counter] = r
			counter += 1

			Motifs.append(strand[r: r + k])

		
		if BestMotifs == []:
			BestMotifs = Motifs

		#Take out one random line of Dna, and construct Profile from the rest	
		LuckyNumber = randint(0,t - 1)
		Ejected = Dna[LuckyNumber] 
		del Motifs[LuckyNumber]

		Grid = ProbabilityProfile(Motifs)


		#Find the probability of each kmer in the Dna line that was taken out.

		EjectedMotifs = {}
		for i in xrange(len(Ejected) - k + 1):
			E = Ejected[i: i + k]
			EjectedMotifs[i] = LineProbability(E, Grid)


		#Draw a kmer from the probabilities of the kmers in Ejected, and replace original choice with it.
				
		WeightedDraw = Weighted_choice(EjectedMotifs)
				
		Section =  RandomDnaNumber[LuckyNumber] #Get the random number that was used to construct motifs
		Dna[LuckyNumber] = Ejected[:Section] + Ejected[WeightedDraw : WeightedDraw + k] + Ejected[Section + k:]
		Motifs.append(Ejected[WeightedDraw : WeightedDraw + k])

		if ScoreMotif(Motifs) < ScoreMotif(BestMotifs):
			BestMotifs = Motifs

	return BestMotifs

def GibbsNumRun(Dna, k, t, NStarts, NLoops):
	Adn = Dna[:]
	BestMotifs = GibbsSampler(Adn, k, t, NStarts)

	for z in xrange(2, NLoops):

		Adn = Dna[:]
		MotifsRun = GibbsSampler(Adn, k, t, NStarts)

		if ScoreMotif(MotifsRun) < ScoreMotif(BestMotifs):
			BestMotifs = MotifsRun

	return BestMotifs

def printingMotifs(list):
	for text in list:
		print text


# Dna = NewFile[1:]

with open("/Users/kanos/Downloads/dataset_163_4.txt") as input_data:
		k,t,N = map(int, input_data.readline().split())
		Dna = [line.strip() for line in input_data.readlines()]

printingMotifs(GibbsNumRun(Dna,k,t,2000,35))



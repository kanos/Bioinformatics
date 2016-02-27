

def CountingPattern(Text,Pattern):
	count = 0
	for i in xrange(0, len(Text) - (len(Pattern) - 1)):
		if Text[i:i + len(Pattern)] == Pattern:
			count += 1
	return count

# CountingWords(text,pattern)

def FrequentWords(Text, k):
	FrequentPatterns = []
	FinalPatterns = []
	Count = []
	for i in xrange(0, len(Text) - k + 1):
		Pattern = Text[i: i+k]
		Count.insert(i, CountingPattern(Text, Pattern))
		MaxCount = max(Count)
	for i in xrange(0, len(Text) - k + 1):
		if Count[i] == MaxCount:
			FrequentPatterns.append(Text[i: i + k])
	for i in FrequentPatterns:
		if i not in FinalPatterns:
			FinalPatterns.append(i)
	print FinalPatterns

# FrequentWords(text,k)


def SymbolToNumber(Number):
	if Number == 0:
		return "A"
	elif Number == 1:
		return "C"
	elif Number == 2:
		return "G"
	else:
		return "T"


def NumbertoPattern(Number, kmer):
	Sequence = []
	letters = ""
	for i in range(kmer):
		Sequence.append(Number%4)
		Number = Number / 4
	for i in Sequence:
		letters += SymbolToNumber(i)
	letters = letters[::-1] #Reverse order of letters
	return letters

def PatternToNumber(Pattern):
	position_in_string = len(Pattern)-1
	Number_value = 0
	for i in Pattern:
		if i == "A" or i =="a":
			#Number_value += 0 # Operation not actually neded.
			position_in_string -= 1
		elif i == "C" or i =="c":
			Number_value += (1 * 4**position_in_string)
			position_in_string -= 1
		elif i == "G" or i =="g":
			Number_value += (2 * 4**position_in_string)
			position_in_string -= 1
		else:
			Number_value += (3 * 4**position_in_string)
			position_in_string -= 1
	return Number_value

def ComputingFrequencies(Text,k):
	FrequencyArray = [0] * (4**k)
	for i in xrange(0, ( len(Text) - k)):
		Pattern = Text[i: i + k]
		j = PatternToNumber(Pattern)
		FrequencyArray[j] += 1
	return FrequencyArray

def ReverseComplement(Text):
	Reverse = ""
	for i in Text:
		if i == 'A':
			Reverse += "T"		
		elif i == 'C':		
			Reverse += "G"
		elif i == 'G':
			Reverse += "C"
		if i == 'T':
			Reverse += "A"
	Reverse = Reverse[::-1]
	return Reverse

def printing (text):
	print ' '.join(map(str, text)) 

def FindPattern(Pattern, Genome):
	places = []
	lengthPat = len(Pattern)
	for i in xrange(0, len(Genome) - lengthPat + 1):
		if Genome[i:i + lengthPat] == Pattern:
			places.append(i)
	return printing(places)


def ClumbFinding1(Genome, k, t, L):
	FrequentPatterns = []
	for i in xrange(0, 4**k - 1):
		Clumb[i] = 0

	for i in xrange(0, len(Genome) - L):
		Text = Genome[i : i + L]
		FrequencyArray = ComputingFrequencies(Text, k)
		for i in xrange(0, 4**k - 1):
			if FrequencyArray(i) >= t: 
				Clumb[i] = 1

	for i in xrange(0, 4**k - 1):
		if Clumb[i] == 1:
			Pattern = PatternToNumber(i,k)
			FrequentPatterns.append(Pattern)

	return FrequentPatterns

def ClumbFinding(Genome, k, L, t):
	FrequentPatterns = []
	Clump = [0] * (4**k)
	Text = Genome[0:L]
	FrequencyArray = ComputingFrequencies(Text, k)

	for i in xrange(0, 4**k - 1):
		if FrequencyArray[i] >= t:
			Clump[i] = 1

	for i in xrange(1, len(Genome) - L):
		FirstPattern = Genome[i - 1: i + k - 1]
		index = PatternToNumber(FirstPattern)
		FrequencyArray[index] -= 1
		LastPattern = Genome[i + L - k: i + L ]
		index = PatternToNumber(LastPattern)
		FrequencyArray[index] += 1
		if FrequencyArray[index] >= t:
			Clump[index] = 1

	for i in xrange(0, 4**k - 1):
		if Clump[i] == 1:
			Pattern = NumbertoPattern(i,k)
			FrequentPatterns.append(Pattern)	

	return len(FrequentPatterns)

def Skew(Genome):
	skewList = [0]
	skew = 0
	for i in Genome:
		if i == "C":
			skew -= 1
			skewList.append(skew)
		elif i == "G":
			skew += 1
			skewList.append(skew)
		else:
			skewList.append(skew)
	return skewList

# Skew("GAGCCACCGCGATA")

def MinimumSkew(Genome):
	skewList = Skew(Genome)
	Minimum = skewList[0]

	for i in skewList:
		if i < Minimum:
			Minimum = i

	indexes = [i for i, x in enumerate(skewList) if x == Minimum]

	return indexes

# MinimumSkew(text)

def HammingDistance(text1, text2):
	distance = 0
	for i in xrange(0,len(text1)):
		if text1[i] != text2[i]:
			distance += 1
	return distance

def AproximatePatternMatching(Pattern, Text, d):
	count = []
	for i in xrange(0, len(Text) - (len(Pattern) - 1)):
		if HammingDistance(Pattern, Text[i: i + len(Pattern)]) <= d:
			count.append(i)
	return count 

def ApproximatePatternCount(Pattern, Text, d):
	count = 0
	for i in xrange(0, len(Text) - (len(Pattern) - 1)):
		if HammingDistance(Pattern, Text[i:i + len(Pattern)]) <= d:
			count += 1
	return count


def Neighbors(Pattern, d): # Take pattern and give all posible mutations with the confine of d Hamming Distance.
	if d == 0:
		return Pattern

	nucletoids = ["C", "T", "G", "A"]

	if len(Pattern) == 1:
		return nucletoids
	

	Neighborhood = []
	SuffixNeighbors = Neighbors(Pattern[1:], d) # Recursive function.

	for Text in SuffixNeighbors:
		if HammingDistance(Pattern[1:], Text) < d:
			for x in nucletoids:
				Neighborhood.append(x + Text) 
		else: 
			Neighborhood.append(Pattern[0] + Text) # If Hamming Distance is bigger than D, add Pattern letter to list.

	return Neighborhood


def FrequentWordsWithMismatches(Text, k , d):
	FinalPatterns = []
	FrequencyArray = [0] * (4**k)
	Close = [0] * (4**k)

	for i in xrange(0, len(Text) - k):
		Neighborhood = Neighbors(Text[i : i + k], d)
		for pattern in Neighborhood:
			index = PatternToNumber(pattern)
			Close[index] = 1

	for i in xrange(0, 4**k - 1):
		if Close[i] == 1:
			Pattern = NumbertoPattern(i, k)
			FrequencyArray[i] = ApproximatePatternCount(Pattern, Text, d)	

	MaxCount = max(FrequencyArray)			

	for i in xrange(0, 4**k - 1):
		if FrequencyArray[i] == MaxCount:
			Pattern = NumbertoPattern(i,k)
			FinalPatterns.append(Pattern)

	print FinalPatterns

def FrequentWordsWithMismatchesDict(Text, k ,d):
	NeighborhoodDict = {}

	for i in xrange(0, len(Text) - k):
		Neighborhood = Neighbors(Text[i: i + k],d)
		for i in Neighborhood:
			NeighborhoodDict[i] = NeighborhoodDict.get(i,0) + 1
		
	MaxCount = max(NeighborhoodDict.values())

	FinalPatterns = [x for x,y in NeighborhoodDict.items() if y==MaxCount]

	return FinalPatterns

def FrequentWordsWithMismatchesAndReverse(Text, k ,d):
	NeighborhoodDict = {}

	for i in xrange(0, len(Text) - k):

		Neighborhood = Neighbors(Text[i: i + k],d)
		for text in Neighborhood:
			NeighborhoodDict[text] = NeighborhoodDict.get(text,0) + 1

		Reverse = ReverseComplement(Text[i: i + k])
		ReverseNeighborhood = Neighbors(Reverse,d)
		for text2 in ReverseNeighborhood:
			NeighborhoodDict[text2] = NeighborhoodDict.get(text2,0) + 1		

	MaxCount = max(NeighborhoodDict.values())

	FinalPatterns = [x for x,y in NeighborhoodDict.items() if y==MaxCount]

	return FinalPatterns



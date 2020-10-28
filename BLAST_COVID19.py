import os

SARS_N_INPUT_FILENAME = "SARS_N.txt"

def fileToString(filename):
	inputFile = open(filename, "r")
	inputlines = inputFile.readlines()
	CompleteString = ""
	for line in inputlines:
		workingString = line.replace(" ", "")
		workingString = workingString.replace("\n", "")
		CompleteString = CompleteString + workingString
	return CompleteString

print(fileToString(SARS_N_INPUT_FILENAME))
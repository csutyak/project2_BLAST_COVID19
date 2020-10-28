import os

SARS_N_FILENAME = "SARS_N.txt"
COVID19_N_FILENAME = "Covid19_N.txt"


def fileToString(filename):
	inputFile = open(filename, "r")
	inputlines = inputFile.readlines()
	CompleteString = ""
	for line in inputlines:
		workingString = line.replace(" ", "")
		workingString = workingString.replace("\n", "")
		CompleteString = CompleteString + workingString
	return CompleteString

sars_n_string = fileToString(SARS_N_FILENAME)
covid_n_string = fileToString(COVID19_N_FILENAME)


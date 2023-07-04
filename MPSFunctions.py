


def preProcessing1(df, threshold = 10):
    import pandas as pd
    import numpy as np
    import pandas as pd
    
           
    # df = pd.DataFrame(df)
    # df = df.drop(index = 0)
    # df.set_index(df["Hybridization REF"], inplace = True)
    # df.drop(columns = ["Hybridization REF"], inplace = True)                                                   
    # df = df.replace("0.0000", "0")
    # df = df.astype('float')
    
    # # df = df.set_index(df["Entrez_Gene_Id"])
    # # df.drop(columns = ["Hugo_Symbol", "Entrez_Gene_Id"], inplace = True)
    # # df = df.dropna(axis = "rows") 

    # # setting the index of df
  
    # df = df.dropna(axis = "rows") 
    
    df = np.log2(df + 1)
    # Creating threshold for mean
    
    mean = df.mean(axis = 1)
    mean.sort_values(inplace = True)

    divideBy = int(100 / threshold)
    percentile = int(len(df.index) /  divideBy)
    meanPercentile = mean[percentile:]    

    #Creating threshold for standard deviation
    std = df.std(axis = 1)
    std.sort_values(inplace = True)

    stdPercentile = std[percentile:]
    

    df = df[df.index.isin(meanPercentile.index & stdPercentile.index)]

    #applying z score
    from scipy.stats import zscore
    df = df.apply(zscore, axis = 1)
    

    return df


"""
Creates Empty Bins for Bins of Gene Expressions and fills them.
Number of bins is default at 10. , inplace = True
df: df must only have the columns of patient gene expression
"""
def binCreation(df):
    row_length = len(df.index)
    col_length = len(df.columns)
    tempMatrix = np.zeros([row_length,col_length])
    tempDF = pd.DataFrame(tempMatrix)
    df = df.reset_index()
    df = pd.concat([df, tempDF], axis = 1)
    df.set_index(df.iloc[:,0], inplace = True)
    df.drop(df.columns[[0]],axis = 1, inplace = True)

    # names the columns
    newColNames = []
    lenColNames = int(len(df.columns) / 2)
    colNames = df.columns[0:lenColNames]
    for i in colNames:
        newColNames.append(i)
    for i in colNames:
        tempName = i + " bin"
        newColNames.append(tempName)
    df.columns = newColNames

    binDivisor = int(len(df.index) / 10)
    NUMBER_OF_PATIENTS = int(len(df.columns) / 2)

    for idx, columnName in enumerate(list(df)):
        if idx >= NUMBER_OF_PATIENTS : 
            break

        df = df.sort_values(by = [columnName])
        for i in range(1,10):
            start = int((i-1) * binDivisor)
            end = int(i*binDivisor)
            df.iloc[start:end,idx+NUMBER_OF_PATIENTS] = i
        for i in range(10,11):
            start = int((i-1) * binDivisor+1)
            end = int(len(df.index)) # number of genes
            df.iloc[start:end,idx+NUMBER_OF_PATIENTS] = i

    return df


"""
Creating modules based off level 5 data z-score
df: the dataframe of experiments
threshold: the z-score of threshold. for instance 2.58 would imply p < 0.01 (values on right of gaussian curve)
sign
sign: either '+' or '-' to determine what you want the threshold to be
numModule: number of modules that are created from dataset
returns a list of lists of modules that fit the criteria (need to be like module[0] -> convert to df to do anything with the data)
"""
def moduleCreator(df, threshold, sign, numModules):
    cols = df.columns.tolist()
    moduleList = []
    count = 0
    for i in range(len(df.columns)):
        if count == numModules:
            break
        column = df.iloc[:, i:i+1]
        if sign == "-":
            moduleList.append(column[column[cols[i]] < threshold])
            
        if sign == '+':
            moduleList.append(column[column[cols[i]] > threshold])
        else:
            print("Pleae enter either '+' or '-' for sign")
            break
        count+=1
    return moduleList


"""
Creates a new column on dataframe for module. Populates that column with zero if 
gene isn't in the module and a 1 if the gene is in the module. Uses REF/Entrez ID

df: the dataframe that's used: specialized for the firebrowse data
module: module whose genes that you are looking for
String columnName: name of column in the dataframe with the module information
returns the same dataframe but with a column with gene inclusion of the module
REF = entrez id
"""

def moduleColumnInclusion(df, module, columnName):
    # initializing column in the dataframe
    df[columnName] = 0

    # splitting up dataframe to make sure that REF and gene symbol are separated
    df.reset_index(inplace = True)
    df[['Gene_Symbol', 'REF']] = df["Hybridization REF"].apply(lambda x: pd.Series(str(x).split("|")))
    df.drop(columns = ["Hybridization REF", 'REF'], inplace = True)
    df.replace("?",np.nan, inplace = True)
    df.dropna(axis = "rows", inplace = True)
    df.set_index(df['Gene_Symbol'], inplace = True)
    df.drop(columns = ['Gene_Symbol'], inplace = True)

    # making sure index is of same type as the module
    # df.index = df.index.astype("str")
    # module = pd.DataFrame(module[0])

    common_genes = np.intersect1d(df.index, module.index)
    for i in df.index:
        if i in common_genes:
            # print(1)
            df.loc[i, columnName] = 1
    # dropping colums
    # df.drop(columns = ['level_0', 'index'], inplace = True)
    
    df.sort_index(ascending = True, inplace = True)

    return df



import math

"""
Calculates a single iteration (or value) of mutualInformation. Needs to be summed up for full calculation.
int jointDistribution: probability of event a and event b
int marginalDistribution1: probability of event a | b
int marginalDistribution2: probability of event b | a
@return int: a single value of mutualInformation calculation
"""
def mutualInformation(jointDistribution, marginalDistribution1, marginalDistribution2):
    return jointDistribution * math.log(jointDistribution /((marginalDistribution1 * marginalDistribution2) +0.000000000000000000001)+0.00000000000000000000000001)

possibleBins = [1,2,3,4,5,6,7,8,9,10]
possibleInclusion = [0,1]

"""
Calculates total mutual information in a dataframe.
df dataframe: dataframe of information
int[] bins: list of bin numbersfirst_module
int[] moduleInclusion: the list of numbers where if a gene is in a module or not
int binPatientColumn: the integer location of the column where the bin of last patient is located
int moduleColumn: the integer location of the module to check inclusion
@return int[]:  a complete mutual information calculation for all patients and provided modules.
"""

def finalMutualInformation(dataframe,bins,moduleInclusion,binPatientColumn,moduleColumn):
    finalSum = 0
    for i in bins:
        a = dataframe.iloc[:,binPatientColumn] == i # first patient
        for j in moduleInclusion: # first module that appears
            b = (dataframe.iloc[:,moduleColumn] == j)
            c = (a&b)
            length = len(dataframe.index)
            sum_c = c.sum()
            sum_a = a.sum()
            sum_b = b.sum()

            d = sum_a / length
            e = sum_b / length
            f = sum_c / length
            finalSum += mutualInformation(f,d,e)


    return finalSum

"""
Function to create a dataframe of randomly shuffled modules that have the same length as genes in module.
# idea: create list of 1s whose length matches length of number of genes in the modules. then, create list of 0s
# that is the length of number of genes minus number of genes in module
# append lists to each other
# shuffle
# make shuffled thing a series
# append to a dataframe
int numRandomization: the number of randomizations desired (1000, 10000)
module= column of dataframe that has the column we are lodf.drop(columns = ['level_0', 'index'], inplace = True)oking for
df: dataframe we are using: need it to reset index and also for number of genes in the dataframe (length of index of desired dataframe)
returns: df -- a dat# idea: create list of 1s whose length matches length of number of genes in the modules. then, create list of 0s
d
"""

def zeroListMaker(n):
    listofzeroes = [0] * n 
    return listofzeroes

def oneListMaker(n):
    listofones = [1] * n
    return listofones

def randomModules(df, module, numRandomizations = 1000):
    # choosing a random seed
    np.random.seed(101)
    df.reset_index(inplace = True, drop = True)
    
    # df.drop(columns = ['REF'])
    index = len(df.index)
    finalDataFrame = pd.DataFrame()
    lengthModule = module.value_counts()[1]
    for i in range(numRandomizations):
        includedGeneLength = oneListMaker(lengthModule)
        otherGene = zeroListMaker(index - lengthModule)
        allGene = pd.Series(includedGeneLength + otherGene)
        pd.Series(np.random.shuffle(allGene))
        finalDataFrame[str(i) + " iteration"] = allGene
         # finalDataFrame = pd.concat([finalDataFrame, allGene],axis = 1)

    return finalDataFrame
    

# calculating MI for all patients
# need to keep patient IDs somewhere
# remember that MI is calculated between bins and moduleInclusion
from sklearn import metrics
# metrics.mutual_info_score(labels_pred = new_z_Score_Values_3.iloc[:,1212], labels_true = new_z_Score_Values_3.iloc[:,2424])
# allMI = []
# for i in range(1212, 2164):
#     allMI.append(metrics.mutual_info_score( new_z_Score_Values_3.iloc[:,i], new_z_Score_Values_3.iloc[:,2424]))

"""
Function for calculating mutual information for multiple modules and multiple patients/bins.
Need bins dataframe, module dataframe, possibleBins, possibleInclusion, number of bins to traverse, number of columns to traverse
For information on calculating mutual information, please refer to finalMutualInformation_1 documentation
returns a list of mutual information z_Score_Values_2
"""
#patients and bins are in the same df

# binDF = new_z_Score_Values_3.iloc[:,1062:2124]
# moduleDF = new_z_Score_Values_3.iloc[:,2124:2125]
def fullFinalMutualInformation(binDF, moduleDF, numModules = 1):
    
    finalAllMI = []
    for i in range(len(binDF.columns)): # for every single value in bin
        binCol = binDF.iloc[:,i]
        for j in range(numModules):
            # print(binDF.iloc[:,i])
            # print(moduleDF.iloc[:,j])
            finalAllMI.append(metrics.mutual_info_score(binCol, moduleDF.iloc[:,j]))
    
    return finalAllMI
# print(allMI)

def zScorewithNull(calcMI, nullMIScores): 
    val = ( (calcMI - nullMIScores.mean()) / nullMIScores.std() ) 
    return val

"""
Function that calculates MPS Score
Uses zScorewithNull function
calcMI: list calculated mutual information scores
nullMIScores: list of null MI scores
df: dataframe of gene expression
moduleDF: df with module inclusion
moduleNum: which module in the moduleDF is selected
"""

def calcMPS(allMI, nullMIScores, df, moduleDF, moduleNum = 0):
    #z scoring mutual information based off null distribution
    val = zScorewithNull(allMI, nullMIScores)

    # makes every value less than zero, zero
    val = pd.DataFrame(val.clip(lower = 0))

    for idx in range(len(val)):
        matrix = np.corrcoef(moduleDF.iloc[:,moduleNum],df.iloc[:,idx])
        if matrix[0][1] < 0:
            val.iloc[idx,0] = val.iloc[idx,0] * -1
            
    
    return val

def write_dataframe_to_csv(dataframe, filename):
    import csv

    dataframe.to_csv(filename, index=False)


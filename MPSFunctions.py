

"""
Process RNA-seq data with thresholding, log2, and z-score across axis = 1 (genes).
df: the initial dataframe with genes for rows and patients for columns
threshold: the percentile you want to cut off for std and mean. (10 percentile means top 90% genes in std and mean are selected)   
"""


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
    return jointDistribution * math.log(jointDistribution /((marginalDistribution1 * marginalDistribution2) +0.00000000000000000000000001)+0.00000000000000000000000001)

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

"""
Function that writes a dataframe to a csv file.
dataframe: dataframe to write to file
filename: name of the file  
"""

def write_dataframe_to_csv(dataframe, filename):
    import csv

    dataframe.to_csv(filename, index=False)


""" 
Function to plot MPS positive and MPS negative curve.
df: dataframe with survival data and/or MPS data
MPS: MPS column of information in the dataframe
event: whether the patient died or not, represented with 1, 0
time: the time the patient was involved in the study
patientID: the id of the patients from the MPS differentiation
"""

def MPSStratification(df, MPS, event, time, patientID):
    
    MPS_Scores_Positive = df[df[MPS].apply(lambda x: x > 0)]
    MPS_Scores_Negative = df[df[MPS].apply(lambda x: x < 0)]

    print(len(MPS_Scores_Positive))
    print(len(MPS_Scores_Negative))

    patient_IDs_positive = MPS_Scores_Positive[patientID].tolist()
    patient_IDs_negative = MPS_Scores_Negative[patientID].tolist()

    survival_data_2_pos = df[df[patientID].isin(patient_IDs_positive)]
    survival_data_2_neg = df[df[patientID].isin(patient_IDs_negative)]

    kmf_pos = KaplanMeierFitter()
    kmf_neg = KaplanMeierFitter()
    kmf_pos.fit(durations=survival_data_2_pos[time], event_observed=survival_data_2_pos[event], label="MPS Positive")
    kmf_neg.fit(durations=survival_data_2_neg[time], event_observed=survival_data_2_neg[event], label='MPS Negative')
    ax = kmf_pos.plot(show_censors=True)
    ax = kmf_neg.plot(ax=ax, show_censors=True)

    # Perform Cox Hazard Ratio test
    cph = CoxPHFitter()
    editDF = pd.concat([MPS_Scores_Positive, MPS_Scores_Negative])
    
    print(editDF)

    cph.fit(editDF[[time, event, MPS]], duration_col=time, event_col=event) # takes into account all of the MPS, want only the MPS_Scores_Positive and MPS_Scores_Negative
    hazard_ratio = cph.hazard_ratios_[MPS]

    # Perform log-rank test
    results = logrank_test(survival_data_2_pos[time], survival_data_2_neg[time], survival_data_2_pos[event],
                           survival_data_2_neg[event])
    p_value = results.p_value

    # Print the p-value and hazard ratio
    ax.annotate('p-value: {:.4f}'.format(p_value), xy=(0.1, 0.3), xycoords='axes fraction', fontsize=12)
    ax.annotate('Hazard Ratio: {:.2f}'.format(hazard_ratio), xy=(0, 0.1), xycoords='axes fraction', fontsize=12)

    # Set the plot labels and title
    ax.set_xlabel('Months')
    ax.set_ylabel('Survival Probability')
    ax.set_title('Survival Curves')


# def MPSStratification(df, MPS, event, time, patientID):
#     MPS_Scores_Positive = df[df[MPS].apply(lambda x : x > 0)]
#     MPS_Scores_Negative = df[df[MPS].apply(lambda x : x < 0)]

#     print(len(MPS_Scores_Positive))
#     print(len(MPS_Scores_Negative))

#     patient_IDs_positive = MPS_Scores_Positive[patientID].tolist()
#     patient_IDs_negative = MPS_Scores_Negative[patientID].tolist()

#     survival_data_2_pos = df[df[patientID].isin(patient_IDs_positive)]
#     survival_data_2_neg = df[df[patientID].isin(patient_IDs_negative)]

#     print(survival_data_2_pos)
#     print(survival_data_2_neg)

#     kmf_pos = KaplanMeierFitter()
#     kmf_neg = KaplanMeierFitter()
#     kmf_pos.fit(durations = survival_data_2_pos[time], event_observed = survival_data_2_pos[event], label = "MPS Positive")
#     kmf_neg.fit(durations = survival_data_2_neg[time], event_observed = survival_data_2_neg[event], label = 'MPS Negative')
#     ax = kmf_pos.plot(show_censors = True)
#     ax = kmf_neg.plot(ax = ax, show_censors = True)

#     #Perform Cox Hazard Ratio test
    
#     cph = CoxPHFitter()
#     print(survival_data_2_neg)
#     print(survival_data_2_pos)
#     cph.fit(pd.concat([MPS_Scores_Positive, MPS_Scores_Negative]), duration_col=time, event_col=event)
#     hazard_ratio = cph.hazard_ratios_[MPS]


    

#     # Perform log-rank test
#     results = logrank_test(survival_data_2_pos[time], survival_data_2_neg[time], survival_data_2_pos[event], survival_data_2_neg[event])
#     p_value = results.p_value

#     # Print the p-value and hazard ratio
#     ax.annotate('p-value: {:.4f}'.format(p_value), xy=(0.1, 0.3), xycoords='axes fraction', fontsize=12)
#     ax.annotate('Hazard Ratio: {:.2f}'.format(hazard_ratio), xy=(0, 0.1), xycoords='axes fraction', fontsize=12)

#     # Set the plot labels and title
#     ax.set_xlabel('Time')
#     ax.set_ylabel('Survival Probability')
#     ax.set_title('Survival Curves')

#     plt.show()
"""  
Function to calculate differentialExpression between level 3 and level 4 data of LINCS
df: the dataframe to apply the operation to
returns calculated level 4 data
"""

def differentialExpression(df):
    mad_func = lambda x: np.median(np.abs(x - np.median(x)))

    mad_values = df.apply(mad_func, axis = 1)

    denom = 1.4826 * mad_values

    denom = list(denom)

    # for i in range(len(denom)):
    #     if denom[i] ==0:
    #         denom[i] = 0.1


    median_values = df.median(axis = "columns")

    numerator = df.subtract(median_values, axis = "index")

    next_result = numerator.divide(denom, axis = "index")

    return next_result


"""  
way to check level 3 to level 4 data
"""

corrList = []
anotherCorrList = []
for idx in range(len(combine.columns)):
    x = combine.iloc[idx,:]
    y = combine.iloc[idx,:]

    combine = pd.concat([x,y], axis = 1)
    combine.replace([np.inf, -np.inf], np.nan, inplace=True)
    combine = combine.dropna(axis = 0)
    x= combine.iloc[:,0]
    y = combine.iloc[:,1]
    
    matrix = np.corrcoef(x,y) #ABY001_SUDHL4_XH_X1_B15
    corrList.append(matrix[0][1])
    r, p = stats.pearsonr(x, y)
    anotherCorrList.append(r)

    


def dictCreation(key, value):
    idToSymbol = zip(key, value) # zip creates a dictionary between gene id and the gene symbol
    idToSymbol = list(idToSymbol)
    idToSymbol = dict(idToSymbol)

    return idToSymbol

def iteratingThroughDict(key, dictionary):

    map_list = []
    # int_first_module = list(map(int, first_module)) # fastest way to convert string to integer
    for i in key: # the column with entrez id in raw_gene_exp file
        for key in dictionary:
            if i == key:
                map_list.append(dictionary[key])

    return map_list

"""  
Scatterplot grapher that does a lot of useful operations such as correlatios, line of best fit
x = first column of info
y = second column of info (must match type and length of x)
title: the name of the graph, if inputted
"""

def graphScatter(x, y, title = ""):
    
    plt.axhline(0, color='black', linestyle='-')
    plt.axvline(0, color='black', linestyle='-')
    plt.scatter(x, y, s= 10, alpha = 0.05)

    corr = np.corrcoef(x, y)[0, 1]
    plt.text(0.75, 0.95, f'Correlation: {corr:.2f}', ha='center', va='top', transform=plt.gca().transAxes)
    m, b = np.polyfit(x, y, deg=1)
    plt.plot(x, m*x + b, color='red')

    plt.text(0.75, 0.85, f'Slope of Line of Best Fit: {m:.2f}', ha='center', va='top', transform=plt.gca().transAxes)
    plt.text(0.75, 0.80, f'y-intercept: {b:.2f}', ha='center', va='top', transform=plt.gca().transAxes)

    plt.scatter(x, y, s= 10, alpha = 0.2, color = 'darkblue')

    x = pd.DataFrame(x)
    y = pd.DataFrame(y)
    plt.xlabel(f"{x.columns}")
    plt.ylabel(f"{y.columns}")

    plt.title(f"{title}")
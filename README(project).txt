This is the README file of the LINCS project and will be updated on a weekly basis.



Week 1:
Initial foray into LINCS. Learning about the data. Learned how to download the data.
Created first module (biologically coherent set of genes) based off Z-score of level 5 data.
Calculated its Mutual Information.
General Tips I've learned:
Slice a few columns of the dataframe to manipulate before using the entire dataframe.
Use lambda x to run operations.

Week 2:
Began looking into MPS scores and understanding KM curves. Mutual information in itself does not carry a lot of information as it is from 0 to 
infinity. Thus, to see if it's mutual information actually matters, we randomize module inclusion (each randomized module inclusion has the
same length of the first module inclusion). Default randomization is 1000. Then, use this randomized module inclusion to calculate MI with 
existing bins. Using this "null" mutual information, calculated null mean and null std which was then applied to each MI previously calculated
from patients. Using this z-score, I "signed" it by taking the correlation sign between gene expression and module inclusion, creating MPS positive
and MPS negative.
From there, I learned KM curves and graphed the MPS positive and MPS negative. However, I used the threshold of 0.01 for my MPS, so I only had 
25 patients meet that statistical measure.

Week 3:
Broadened scope of MPS scores, creating KM curve that kind of looks like Balaji's. Now, trying to use np.randomshuffle and clean up my code.
Need to have everything be a function, tired of taking hours going back and doing things.

Tuesday:
created KM curve that matches Balaji's but my MPS scores are weird. weakly positively correlated with Balaji's when they should be strongly correlated. started creating lots of functions and trying to get rid of things.
got stuck on calculating MI for randomization values. got proper tcga dataset from balaji and need to explore it. 

Wednesday:
Brief explanation of my project so far:
[Describe difference between RNA-seq and LINCS]

Essentially, 720,000 experiments, we have 12,000 genes whose expression has been measured. There are five categories of data, using the most processed form. This means that the gene expression values have been z-scored across their replicate experiments.
I create modules, which are biologically coherent set of genes, from each experiment, with z-score greater than 2.58 (p < 0.01). I then use this module from the dataset and stop using LINCS.
Then, I use a dataset of gene expression of BRCA patients (1082 patients [some with primary tumors, others with secondary tumors, need to figure that out]). I am trying to see if these modules I create are "informative" of the gene expression of the patients. To do this,
I sort each gene expression and split the values ten different ways into different discrete "bins". I then see if these bins are mutually informative of the module inclusion. I calculate a mutual information score for each patient and each module. But what is mutual information?
[briefly go into math and its differences from correlation]

Mutual information ranges from 0 to infinity, so how do I know if my mutual information values are actually good? To remedy this, I create a null distribution based on randomizing module inclusion and then calculate MI. If my MI for a patient is in top 1 % then, that is significant.
To sign my mutual information, I find the Pearson correlation between gene expression and gene module inclusion and sign my values accordingly. Finally, I use TCGA survival data to plot KM curves and see if my modules stratify patient survival.

[explain KM curves]

today, created functions up to MPS with new data from firebrowse. runs pretty smoothly. now need to read in survival data

Thursday: 
Spent a long time creating functions for the past two days just having to change a lot due to transitioning to gene symbol to check work with Balaji. should have just started with this.
made a new file with just my functions
going to start importing my functions from another file just to clean up my jupyter notebook since it's kind of hellish having to scroll past everything over and over again   
Action Items: Need to figure out the gene differences
            : figure out the gene symbol mapping
            : need to create a km curve tomorrow that matches balaji's

Friday: should make my functions general, so I can just preprocess data myself. did that. 
        had journal club today, so didn't really start work until 2 pm.  simran was able to recreate balaji's results just using his data
        i'm still struggling. the genes don't match up but balaji said it shouldn't matter too much since those genes are RNA anti-sense proteins and most of the stuff that we deal with is protein coding


Saturday: unsighed MI matches
Sunday: just got off the call with Balaji. Need to think about a few things: need to save every coding file I create to make sure I can go back and see what I've been doing. learned about another way to randomize modules. what do I want to do in the lab by beginning of September.

Week 3
Monday: narrwoed it down to MI. my MI calculation is off. I don't know why. has to do with module inclusion, I believe. should be 109 genes instead of 112?
Figured it out. Needed to use common genes that Balaji provided to filter 109 genes out of 112.

Tuesday: not a lot of programming since I came in late. researched drugs I was interested in.

Wednesday: met with balaji. there are genes in common genes that aren't in the raw rsem data he provided, which is strange. he says he will look into it more himself. calculated common_genes from like 5 cancers that's how he got common_genes0.1. need to move forward with my LINCS stuff

Thursday: created figure for journal club. calculated p-value for KM curve. converted some level 3 data to level 4 but some has really low correlation; idk why, maybe it's because of my zeroes. started working on calculating hazard randomization

Friday: trying to calculate hazard ratio. will try to figure out level 3 to level 4 conversion problems

Saturday/Sunday: break

Week 4
So far this week, finished verifying how the level 3 to level 4 and level 4 to level 5 data is transformed. Now that I understand the LINCS data, I can use it to create modules which I can input into my MPS functions and see if the KM curves are significant.

Monday: converting between level 4 and level 5 data. some values in the metadata like "ABY001_A375_XH:ADO-TRASTUZUMAB_EMTANSINE:5:24" just doesn't exist in the level 5 or level 4 data.

Tuesday: figuring out how to clone repositories and how to constantly save my files on GitHub
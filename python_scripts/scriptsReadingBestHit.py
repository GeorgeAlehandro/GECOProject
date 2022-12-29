# Simple Way to Read TSV Files in Python using pandas
# importing pandas library
import pandas as pd
def BlastpAlignmentLengthFilter(filename):
    best_hit_df = pd.read_csv(filename, sep='\t')
    filt = best_hit_df[(best_hit_df['identpercentage'] > 30) & (best_hit_df[['qlen','slen']].min(axis=1)/best_hit_df['alignmentlength']*100 > 80)]
    # printing data
    print((filt["alignmentlength"]*filt["identpercentage"]).sum()/filt["qlen"].sum()/100)
BlastpAlignmentLengthFilter('~/GECO/for_blastp/TOP1_interestVSproche1')
BlastpAlignmentLengthFilter('~/GECO/for_blastp/TOP1_interestVSproche2')
BlastpAlignmentLengthFilter('~/GECO/for_blastp/TOP1_interestVSloin1')
BlastpAlignmentLengthFilter('~/GECO/for_blastp/TOP1_interestVSloin2')

def BlastnAlignmentLengthFilter(filename):
    best_hit_df = pd.read_csv(filename, sep='\t')
    filt = best_hit_df[(best_hit_df['identpercentage'] >= 70) & (best_hit_df[['qlen','slen']].min(axis=1)/best_hit_df['alignmentlength']*100 >= 70)]
    return filt

# Calcule ANI
interestVSproche1 = BlastnAlignmentLengthFilter('~/GECO/for_blastn/TOP1_interestVSproche1')
interestVSproche2 = BlastnAlignmentLengthFilter('~/GECO/for_blastn/TOP1_interestVSproche2')
interestVSloin1 = BlastnAlignmentLengthFilter('~/GECO/for_blastn/TOP1_interestVSloin1')
interestVSloin2 = BlastnAlignmentLengthFilter('~/GECO/for_blastn/TOP1_interestVSloin2')

#ANI blastn
def CalculateANI(table):
    print((table["alignmentlength"]*table["identpercentage"]).sum()/table["qlen"].sum()/100)
print("ANI")
CalculateANI(interestVSproche1)
CalculateANI(interestVSproche2)

CalculateANI(interestVSloin1)
CalculateANI(interestVSloin2)

#AF blastn
def CalculateAF(table):
    print(table["qlen"].sum()/1988613)
print("AF")
CalculateAF(interestVSproche1)
CalculateAF(interestVSproche2)

CalculateAF(interestVSloin1)
CalculateAF(interestVSloin2)
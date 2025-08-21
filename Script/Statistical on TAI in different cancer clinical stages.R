library(myTAI)
library(readxl)

data_cancer_8stage <- read_excel("merge_rna_stage.xlsx")

# Hourglass model analysis
PlotSignature( ExpressionSet = data_cancer_8stage,
               measure       = "TAI",  # "TAI" for Transcriptome Age Index; "TDI" for Transcriptome Divergence Index; "TPI" for Transcriptome Polymorphism Index
               TestStatistic = "FlatLineTest",  # Options: FlatLineTest, ReductiveHourglassTest(HLH), EarlyConservationTest(LHH), LateConservationTest(HHL), ReverseHourglassTest(LHL)
               # modules = list(early = 1, mid = 2:3, late = 4),
               xlab          = "Ontogeny", 
               ylab          = "TAI" )